#!/usr/bin/env python
"""
Important Notice:
--references has been removed, please use --reference (multiple uses are ok)

align FASTQ, SAM, or BAM file (gzip and bzip2 supported)
with bowtie, produces BAM files (sorted and indexed)
output will be in ./align
In align, there will be a folder for each reference genome
    e.g. align/ref1 , align/ref2, align/ref1, align/ref2
In reference folder, there will be a folder for unique or random alignments
    (unique means only maps to one spot, random is maq-like behavior)
    e.g. align/ref1/random, align/ref1/unique, etc.
"""
import pysam
import logging
import os
from os import getcwd, environ, curdir, makedirs, devnull
from os.path import join, exists, splitext
from textwrap import dedent
import sys
from tempfile import mkdtemp
from subprocess import Popen, PIPE
import scripter
from scripter import assert_path, path_to_executable, Usage, \
                     exit_on_Usage, InvalidFileException, get_logger, \
                     Environment, critical, debug
from seriesoftubes.tubes.polledpipe import PolledPipe
from seriesoftubes.tubes import wait_for_job
from seriesoftubes.fnparsers import BowtieFilenameParser
    
from pkg_resources import get_distribution
__version__ = get_distribution('seriesoftubes').version
VERSION = __version__

def main():
    e = Environment(version=VERSION, doc=__doc__)
    e.set_filename_parser(BowtieFilenameParser)
    # let bowtie do the multiprocessing
    e.override_num_cpus(1)
    parser = e.argument_parser
    parser.add_argument('--path-to-bowtie', nargs='?',
                        default=path_to_executable('bowtie', '/usr/local/bowtie-*',
                                                   environ='SOT_PATH_TO_BOWTIE'),
                        help='The path to the bowtie executable')
    parser.add_argument('--path-to-samtools', nargs='?',
                        default=path_to_executable('samtools', '/usr/local/samtools*',
                                                   environ='SOT_PATH_TO_SAMTOOLS'),
                        help='The path to the samtools executable')
    # fix aliases, should be --ref too
    parser.add_argument('--reference', dest='references', action='append',
                        help=dedent('''\
    Reference genome to align against (either a bowtie index name or file, or a
    fasta file). This flag may be called multiple times (which will cause each 
    reference to be aligned to separately). If no references are specified,
    we'll look the for environment variable SOT_DEFAULT_REFERENCES, which 
    should be given as a list, e.g. "foo foo2 foo3"'''),
                        )
    parser.add_argument('--no-unique', dest='unique', action='store_false',
                        help='do not produce unique/ alignment folder') 
    parser.add_argument('--no-random', dest='random', action='store_false',
                        help='do not produce random/ alignment folder') 
    parser.add_argument('--ignore-quality', dest='use_quality',
                        action='store_false',
                        help=dedent('''\
    Use -v mode with bowtie, allows only n mismatches total. Also applies to
    counter-references if any are called'''))
    parser.add_argument('--mismatches', default='2',
                        help=dedent('''\
    allow n mismatches, in the seed (default) or total if --ignore-quality (-v mode)'''))
    parser.add_argument('--quals-type', default='solex1.3',
                        choices=['solexa', 'solexa1.3', 'phred64', 'phred33',
                                 'integer'],
                        help='Valid options are integer, solexa1.3, solexa, phred33, or phred64 (see bowtie for more info)')
    parser.add_argument('--max-quality', default='70',
                        help=dedent('''\
    specify maximum quality scores of all mismatched positions (default is 70),
    ignored in --ignore-quality (-v) mode'''))
    parser.add_argument('--seed-length', dest='seed_len', default='28',
                        help='use seed length of m (default is 28)')
    cparser = parser.add_argument_group('counter-alignments',
                                        description=dedent('''\
    specify counter-reference genome(s)/sequence(s) to use for filtering out
    unwanted reads.'''))
    cparser.add_argument('--counter-reference', dest='counter_references',
                         action='append',
                        help=dedent('''\
    Optional counter-reference genome/sequences to align against (either a bowtie
    index name or file, or a fasta file). This flag may be called multiple times.
    All counter-references will be concatenated into one index, and reads will
    be aligned in --no-unique (-M 1) mode. Any reads which align will be saved
    in a separate directory called 'bad_reads' and not aligned against the
    reference genomes/sequences. If no counter-references are specified, we'll 
    look the for environment variable SOT_DEFAULT_COUTNER_REFERENCES,
    which should be given as a list, e.g. "foo foo2 foo3"'''),
                        )
    cparser.add_argument('--counter-mismatches', default=None,
                        help=dedent('''\
    allow n mismatches to counter-reference(s), in the seed (default) or total
    if --ignore-quality (-v mode). Default: same as references'''))
    cparser.add_argument('--counter-max-quality', default='70',
                        help=dedent('''\
    specify maximum quality scores of all mismatched positions when aligning to
    counter-reference(s) (default is 70), ignored in --ignore-quality (-v) mode'''))
    context = e.get_context()
    new_references = validate_references(**context)
    new_counter_references = cat_counter_references(**context)
    e.update_context({'references': new_references,
                      'counter_references': new_counter_references})
    e.do_action(align)

def fastq_to_bowtie(fasta_file, target_dir=curdir, path_to_bowtie='bowtie'):
    """given a filename, makes a bowtie index
    if that file is a FASTA file
    """
    if exists(fasta_file):
        f = open(fasta_file, 'rU')
        for line in f:
            if line.startswith('#'): continue
            elif line.startswith('>'):
                args = [path_to_bowtie + '-build', fasta_file,
                        join(target_dir, fasta_file)]
                debug(' '.join(args))
                P = Popen(args, stdout=open(devnull, 'w'), stderr=PIPE)
                stderr = P.communicate()[1]
                if stderr.splitlines()[0].startswith('Error'): return None
                else: return join(getcwd(), target_dir, fasta_file)
    return None

def find_bowtie_index(r, path_to_bowtie='bowtie'):
    """check for bowtie index as given. return True if found, else return False"""
    args = [path_to_bowtie + '-inspect', '-v', '-s', r]
    debug(' '.join(args))
    P = Popen(args, stdout=open(devnull, 'w'), stderr=PIPE, cwd=mkdtemp())
    stderr = P.communicate()[1].splitlines()
    if not stderr[0].startswith('Could not locate'):
        for line in stderr:
            if line.startswith('Opening'):
                index_ebwt1 = line[(1+line.find('"')):line.rfind('"')]
                index_basename = index_ebwt1[0:index_ebwt1.find('.1.ebwt')]
                return index_basename
    rprime = join(getcwd(), r)
    args = [path_to_bowtie + '-inspect', '-v', '-s', rprime]
    debug(' '.join(args))
    P = Popen(args, stdout=open(devnull, 'w'), stderr=PIPE, cwd=mkdtemp())
    stderr = P.communicate()[1].splitlines()
    if not stderr[0].startswith('Could not locate'):
        for line in stderr:
            if line.startswith('Opening'):
                index_ebwt1 = line[(1+line.find('"')):line.rfind('"')]
                index_basename = index_ebwt1[0:index_ebwt1.find('.1.ebwt')]
                return index_basename
    return None

def cat_counter_references(counter_references=None, target_dir=curdir,
                            path_to_bowtie='bowtie',
                            logger=None, **kwargs):
    if counter_references is None: return
    try: makedirs(target_dir, mode=0755)
    except OSError: pass
    debug('Validating counter-references and building counter-reference index')
    valid_references = validate_references(references=counter_references,
                                           target_dir=target_dir,
                                           path_to_bowtie=path_to_bowtie,
                                           logger=logger,
                                           environ_key='SOT_DEFAULT_COUNTER_REFERENCES')
    crefs_fa = open(join(target_dir, 'counter_references.fa'), 'w')
    for ref in counter_references:
        Popen([path_to_bowtie + '-inspect', ref], stdout=crefs_fa).wait()
    crefs_index = join(target_dir, counter_references)
    args = [path_to_bowtie + '-build', crefs_fa, crefs_index]
    P = Popen(args, stderr=PIPE)
    stderr = P.communicate()[1]
    if stderr.startswith('Error'):
        critical(stderr)
        critical('No counter-references will be used.')
    return crefs_index

def validate_references(references=None, path_to_bowtie='bowtie',
                        logger=None, environ_key='SOT_DEFAULT_REFERENCES',
                        target_dir=curdir,
                        **kwargs):
    makedirs(target_dir, mode=0755)
    debug('Validating references')
    new_references = []
    if references is None:
        if environ.has_key(environ_key):
            references = environ[environ_key].split()
        else:
            critical('no reference genomes specified')
            return []
    
    for r in references:
        bowtie_index = find_bowtie_index(r, path_to_bowtie=path_to_bowtie)
        if bowtie_index is None:
            if exists(r):
                debug('Attempting to build bowtie index from %s' % r)
                new_index = fastq_to_bowtie(r, target_dir=target_dir,
                                            path_to_bowtie=path_to_bowtie)
                if new_index is not None:
                    new_references.append(new_index)
                    continue
                else:
                    critical('Failed to build bowtie index.')
            critical('bowtie could not find the index for %s', r)
            critical('we will not align to %s', r)
        else:
            new_references.append(bowtie_index)
    return new_references

@exit_on_Usage
def align(fp_obj, references=[], counter_references=None,
          random=True, unique=True, max_quality='70',
          quals_type='solexa1.3',  mismatches='2', seed_len='28',
          counter_mismatches=None,
          use_quality=True, logging_level=10, num_cpus=1, **kwargs):
    if not unique and not random: raise Usage('Nothing to do')
    common_flags = ['-y', '-a', '--time', '--best', '--chunkmbs', '1024',
                    '--strata', '--sam']
    common_flags.extend(['-p', str(num_cpus)])
    logger = get_logger(logging_level)
    uniqueness = {}
    if unique: uniqueness.update({'unique': ['-m','1']})
    if random: uniqueness.update({'random': ['-M','1']})
    stdout_buffer = []
    common_flags.extend(['-l', seed_len])
    
    new_sources = []
    if counter_references is not None:
        #counter align first
        flags = [item for item in common_flags]
        flags.extend(uniqueness['random'])
        if counter_mismatches is None: counter_mismatches = mismatches
        if use_quality:
            flags.extend(['-e', max_quality])
            flags.extend(['-n', counter_mismatches])
        else:
            flags.extend(['-v', counter_mismatches])
            
        if fp_obj.paired_end: flags.extend(['-X','600'])
        new_filenames = counteralign_once(fp_obj, flags, counter_references,
                              logger=logger, **kwargs)
        
        # after alignment
        fp_obj.input_file = new_filenames[0]
        fp_obj.second_file = new_filenames[1]
        fp_obj.use_pysam = True
        fp_obj.format = 'BAM'
    
    for match_type, match_flag in uniqueness.items():
        flags = [item for item in common_flags]
        flags.extend(match_flag)
        
        source = fp_obj.fastq_source
        # In case we're from that bad day...
##        if source=='081124_HWI-EAS355_0001_Meghan' and match_type=='unique':
##            flags.extend(['-v', '3'])
#        elif use_quality:
        if use_quality:
            flags.extend(['-e', max_quality])
            flags.extend(['-n', mismatches])
        else:
            flags.extend(['-v', mismatches])

        if fp_obj.paired_end: flags.extend(['-X','600'])
        
        for ref in references:
            s = align_once(fp_obj, flags, ref, match_type, logger=logger, **kwargs)
            stdout_buffer.append(s)
    return '\n'.join([s for s in stdout_buffer if s is not None])

def counteralign_once(fp_obj, flags, ref, **kwargs):
    """Produce counter-alignements"""
    refname = os.path.basename(ref)
    output_dir, output_file = os.path.split(fp_obj.tmp_filename(refname))
    fp_obj.check_output_dir(output_dir)
    filename1 = os.path.abspath(fp_obj.input_file)
    second_file = fp_obj.second_file
    if second_file is not None: filename2 = os.path.abspath(second_file)
    else: filename2 = None
    if use_quality:
        if fp_obj.use_pysam: flags.append('--phred33-quals')
        else: flags.append(''.join(['--', quals_type, '-quals']))
    if fp_obj.paired_end:
        file_args = [ref, '--12', '-']
        logger.info('Automagically interpreting %s files', fp_obj.format)
    else:
        logger.info('Automagically interpreting %s file', fp_obj.format)
        file_args = [ref, '-']
    bowtie_args = [path_to_bowtie] + flags + file_args
    
    # finish parsing input here
    input_stderr = PolledPipe(logger=logger, level=logging.ERROR)
    bowtie_stderr = PolledPipe(logger=logger, level=logging.ERROR)
    if fp_obj.use_pysam:
        if fp_obj.paired_end:
            in_args = [sys.executable, '-m', 'seriesoftubes.converters.bamtotab',
                       filename1]
        else:
            in_args = [sys.executable, '-m', 'seriesoftubes.converters.bamtofastq',
                       filename1]
    elif fp_obj.paired_end and fp_obj.format =='FASTQ':
        in_args = [sys.executable, '-m', 'seriesoftubes.converters.fastqtotab',
                   filename1, filename2]
    elif fp_obj.format == 'FASTQ':
        in_args = [sys.executable, '-m', 'seriesoftubes.converters.cat',
                   filename1]
    else:
        logger.critical("Couldn't figure out what to do with file %s of format %s",
                        fp_obj.input_file, fp_obj.format)
    logger.info(' '.join(in_args))
    input_reader = Popen(in_args, stdout=PIPE, stderr=input_stderr.w,
                         bufsize=-1)
    logger.info('Launching bowtie (output will be piped to samtools)')
    logger.info(' '.join(bowtie_args))
    bowtie_aligner = Popen(bowtie_args, stdin=input_reader.stdout,
                           stdout=PIPE, stderr=bowtie_stderr.w,
                           bufsize=-1)
    logger.info('Only unaligned reads will be saved.')
    samtools_args = [path_to_samtools, 'view', '-b', '-S', '-o',
                     '-f', '0x4',# ONLY SAVE UNALIGNED READS
                     join(output_dir, output_file), '-']
    logger.info('Launching samtools to encode bowtie output as BAM')
    logger.info(' '.join(samtools_args))
    samtools_stdout = PolledPipe(logger=logger, level=logging.WARN)
    samtools_stderr = PolledPipe(logger=logger, level=logging.ERROR)
    samtools_viewer = Popen(samtools_args, stdin=bowtie_aligner.stdout,
                            stdout=samtools_stdout.w,
                            stderr=samtools_stderr.w, bufsize=-1)
    
    logger.debug('Waiting for bowtie to finish')
    pollables = [input_stderr, bowtie_stderr, samtools_stdout, samtools_stderr]
    wait_for_job(bowtie_aligner, pollables, logger)
    
    if not bowtie_aligner.returncode == 0:
        logger.critical("bowtie did not run properly [%d]",
                        bowtie_aligner.returncode)
        samtools_viewer.terminate()
        samtools_viewer.poll()
        logger.critical("samtools terminated")
        return
    
    logger.debug('Alignment successfully completed')
    logger.debug('Waiting for samtools to finish')
    wait_for_job(samtools_viewer, [samtools_stdout, samtools_stderr], logger)
    if not samtools_viewer.returncode == 0:
        logger.critical("samtools view did not run properly [%d]",
                        samtools_viewer.returncode)
        return
        
    logger.debug('Unsorted BAM file successfully written')
    
    return (join(output_dir, output_file), None)

def align_once(fp_obj, flags, ref, match_type, use_quality=False,
               quals_type='solexa1.3',
               path_to_bowtie= None, path_to_samtools=None, logger=None,
               **kwargs):
    refname = os.path.basename(ref)
    path_to_unsorted = fp_obj.tmp_filename(refname, match_type)
    output_dir = os.path.split(path_to_unsorted)[0]
    fp_obj.check_output_dir(output_dir)
    filename1 = os.path.abspath(fp_obj.input_file)
    second_file = fp_obj.second_file
    if second_file is not None: filename2 = os.path.abspath(second_file)
    else: filename2 = None
    if use_quality:
        if fp_obj.use_pysam: flags.append('--phred33-quals')
        else: flags.append(''.join(['--', quals_type, '-quals']))
    if fp_obj.paired_end:
        file_args = [ref, '--12', '-']
        logger.info('Automagically interpreting %s files', fp_obj.format)
    else:
        logger.info('Automagically interpreting %s file', fp_obj.format)
        file_args = [ref, '-']
    bowtie_args = [path_to_bowtie] + flags + file_args
    
    # finish parsing input here
    input_stderr = PolledPipe(logger=logger, level=logging.ERROR)
    bowtie_stderr = PolledPipe(logger=logger, level=logging.ERROR)
    if fp_obj.use_pysam:
        if fp_obj.paired_end:
            in_args = [sys.executable, '-m', 'seriesoftubes.converters.bamtotab',
                       filename1]
        else:
            in_args = [sys.executable, '-m', 'seriesoftubes.converters.bamtofastq',
                       filename1]
    elif fp_obj.paired_end and fp_obj.format =='FASTQ':
        in_args = [sys.executable, '-m', 'seriesoftubes.converters.fastqtotab',
                   filename1, filename2]
    elif fp_obj.format == 'FASTQ':
        in_args = [sys.executable, '-m', 'seriesoftubes.converters.cat',
                   filename1]
    else:
        logger.critical("Couldn't figure out what to do with file %s of format %s",
                        fp_obj.input_file, fp_obj.format)
    logger.info(' '.join(in_args))
    input_reader = Popen(in_args, stdout=PIPE, stderr=input_stderr.w,
                         bufsize=-1)
    logger.info('Launching bowtie (output will be piped to samtools)')
    logger.info(' '.join(bowtie_args))
    bowtie_aligner = Popen(bowtie_args, stdin=input_reader.stdout,
                           stdout=PIPE, stderr=bowtie_stderr.w,
                           bufsize=-1)
    
    samtools_args = [path_to_samtools, 'view', '-b', '-S', '-o',
                     path_to_unsorted, '-']
    logger.info('Launching samtools to encode bowtie output as BAM')
    logger.info(' '.join(samtools_args))
    samtools_stdout = PolledPipe(logger=logger, level=logging.WARN)
    samtools_stderr = PolledPipe(logger=logger, level=logging.ERROR)
    samtools_viewer = Popen(samtools_args, stdin=bowtie_aligner.stdout,
                            stdout=samtools_stdout.w,
                            stderr=samtools_stderr.w, bufsize=-1)
    
    logger.debug('Waiting for bowtie to finish')
    pollables = [input_stderr, bowtie_stderr, samtools_stdout, samtools_stderr]
    wait_for_job(bowtie_aligner, pollables, logger)
    
    if not bowtie_aligner.returncode == 0:
        logger.critical("bowtie did not run properly [%d]",
                        bowtie_aligner.returncode)
        samtools_viewer.terminate()
        samtools_viewer.poll()
        logger.critical("samtools terminated")
        return
    
    logger.debug('Alignment successfully completed')
    logger.debug('Waiting for samtools to finish')
    wait_for_job(samtools_viewer, [samtools_stdout, samtools_stderr], logger)
    if not samtools_viewer.returncode == 0:
        logger.critical("samtools view did not run properly [%d]",
                        samtools_viewer.returncode)
        return
        
    logger.debug('Unsorted BAM file successfully written')
    
    logger.info('Launching samtools again to sort BAM output')
    output_dir, output_file = os.path.split(path_to_unsorted)
    bam_file = os.path.splitext(output_file)[0]
    sorter_args = [path_to_samtools, 'sort', output_file, bam_file]
    logger.info(' '.join(sorter_args))
    samtools_stdout = PolledPipe(logger=logger, level=logging.WARN)
    samtools_stderr = PolledPipe(logger=logger, level=logging.ERROR)
    samtools_sorter = Popen(sorter_args, stdout=samtools_stdout.w,
                            stderr=samtools_stderr.w, cwd=output_dir)
    wait_for_job(samtools_sorter, [samtools_stdout, samtools_stderr], logger)
    if not samtools_sorter.returncode == 0:
        logger.critical("samtools sort did not run properly [%d]",
                        samtools_sorter.returncode)
        return
    
    # don't destroy the files until we're sure we succeeded!
    assert_path(os.path.join(output_dir, bam_file + '.bam'))
    logger.debug('Removing unsorted file %s', path_to_unsorted)
    os.remove(path_to_unsorted)
    
    logger.debug('Launching samtools again to index sorted BAM output') 
    samtools_stdout = PolledPipe(logger=logger, level=logging.WARN)
    samtools_stderr = PolledPipe(logger=logger, level=logging.ERROR)
    index_args = [path_to_samtools, 'index', bam_file + '.bam']
    samtools_indexer = Popen(index_args, stdout=samtools_stdout.w,
                            stderr=samtools_stderr.w, cwd=output_dir)
    wait_for_job(samtools_indexer, [samtools_stdout, samtools_stderr], logger)
    if not samtools_indexer.returncode == 0:
        logger.critical("samtools index did not run properly [%d]",
                        samtools_indexer.returncode)
        return
    
    # Make sure indexing succeeds
    assert_path(os.path.join(output_dir, bam_file + '.bam.bai'))
    return

if __name__=="__main__": main()
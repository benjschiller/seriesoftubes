#!/usr/bin/env python
"""
align FASTQ, SAM, or BAM file (gzip and bzip2 supported)
with bowtie2, produces BAM files (sorted and indexed)
output will be in ./align2
In align2, there will be a folder for each reference genome
    e.g. align2/ref1 , align2/ref2, align2/ref1, align2/ref2
"""
import pysam
import logging
import os
from os import getcwd, environ, curdir, makedirs, devnull
from os.path import join, exists, splitext
from textwrap import dedent
import sys
from tempfile import mkdtemp
from subprocess import Popen, PIPE, STDOUT
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
    # let bowtie2 do the multiprocessing
    e.override_num_cpus(1)
    parser = e.argument_parser
    parser.add_argument('--path-to-bowtie2', nargs='?',
                        default=path_to_executable('bowtie2', '/usr/local/bowtie2-*',
                                                   environ='SOT_PATH_TO_BOWTIE2'),
                        help='The path to the bowtie2 executable')
    parser.add_argument('--path-to-samtools', nargs='?',
                        default=path_to_executable('samtools', '/usr/local/samtools*',
                                                   environ='SOT_PATH_TO_SAMTOOLS'),
                        help='The path to the samtools executable')
    # fix aliases, should be --ref too
    parser.add_argument('--reference', dest='references', action='append',
                        help=dedent('''\
    Reference genome to align against (either a bowtie2 index name or file, or a
    fasta file). This flag may be called multiple times (which will cause each 
    reference to be aligned to separately). If no references are specified,
    we'll look the for environment variable SOT_DEFAULT_REFERENCES, which 
    should be given as a list, e.g. "foo foo2 foo3"'''),
                        )
    parser.add_argument('--ignore-quality', dest='use_quality',
                        action='store_false',
                        help=dedent('''\
    Ignore quality scores if available. Also applies to
    counter-references if any are called'''))
    cparser = parser.add_argument_group('counter-alignments',
                                        description=dedent('''\
    specify counter-reference genome(s)/sequence(s) to use for filtering out
    unwanted reads.'''))
    cparser.add_argument('--counter-reference', dest='counter_references',
                         action='append',
                        help=dedent('''\
    Optional counter-reference genome/sequences to align against (either a bowtie2
    index name or file, or a fasta file). This flag may be called multiple times.
    All counter-references will be concatenated into one index, and reads will
    be aligned in --fast mode. Any reads which align will be saved
    in a separate directory called 'counteraligned' and not aligned against the
    reference genomes/sequences. If no counter-references are specified, we'll 
    look the for environment variable SOT_DEFAULT_COUTNER_REFERENCES,
    which should be given as a list, e.g. "foo foo2 foo3"'''),
                        )
    context = e.get_context()
    new_references = validate_references(**context)
    new_counter_references = cat_counter_references(**context)
    e.update_context({'references': new_references,
                      'counter_references': new_counter_references})
    sequence = e.get_sequence(**context)
    e._sequence = merge_pairs(sequence)
    e.do_action(align2)

def merge_pairs(old_sequence):
    """pair sequence files that differ only by 1/2 in one position
    """
    new_sequence = []
    sequence = []
    for fp_obj in old_sequence:
        if fp_obj.paired_end: new_sequence.append(fp_obj)
        else: sequence.append(fp_obj)
    input_files = [fp_obj.input_file for fp_obj in sequence]
    mates = []
    for i in xrange(len(input_files)):
        if i in mates: continue
        these_mates = []
        for j in xrange(i, len(input_files)):
            identity = [(not x==y) for x,y in zip(input_files[i], input_files[j])]
            if sum(identity) == 1:
                index = identity.index(True)
                if (input_files[i][index] == '1' and
                    input_files[j][index] == '2') or \
                   (input_files[i][index] == '2' and
                    input_files[j][index] == '1'):
                    these_mates.append(j)
                    print i,j,input_files[i],input_files[j]
        if len(these_mates) == 0: continue
        elif len(these_mates) == 1:
            j = these_mates[0]
        else:
            # User input required
            print "Ambiguous filename pairing"
            print "Please select the correct mate pair for %s:" % input_files[i]
            while True:
                for j in these_mates:
                    print "[%d] %s" % (j, input_files[j])
                rawinput = raw_input("Enter your choice: ")
                try: choice = int(rawinput.strip())
                except ValueError: continue
                if not choice in these_mates:
                    print "%d is not a valid choice" % choice
                else:
                    j = choice
                    break
        identity = [(not x==y) for x,y in zip(input_files[i], input_files[j])]
        index = identity.index(True)
        if input_files[i][index] == '1' and input_files[j][index] == '2':
            sequence[i].second_file = sequence[j].input_file
            sequence[i].paired_end = True
            mates.append(j)
        elif input_files[j][index] == '1' and input_files[i][index] == '2':
            sequence[j].second_file = sequence[i].input_file
            sequence[j].paired_end = True
            mates.append(i)

    new_sequence.extend([x for i, x in enumerate(sequence) if not i in mates])
    return new_sequence

def fasta_to_bowtie2(fasta_file, target_dir=curdir, path_to_bowtie2='bowtie2'):
    """given a filename, makes a bowtie2 index
    if that file is a FASTA file
    """
    if exists(fasta_file):
        f = open(fasta_file, 'rU')
        for line in f:
            if line.startswith('#'): continue
            elif line.startswith('>'):
                args = [path_to_bowtie2 + '-build', fasta_file,
                        join(target_dir, fasta_file)]
                debug(' '.join(args))
                P = Popen(args, stdout=open(devnull, 'w'), stderr=PIPE)
                stderr = P.communicate()[1]
                if stderr.splitlines()[0].startswith('Error'): return None
                else: return join(getcwd(), target_dir, fasta_file)
    return None

def find_bowtie2_index(r, path_to_bowtie2='bowtie2'):
    """check for bowtie2 index as given. return True if found, else return False"""
    args = [path_to_bowtie2 + '-inspect', '-v', '-s', r]
    debug(' '.join(args))
    P = Popen(args, stdout=open(devnull, 'w'), stderr=PIPE, cwd=mkdtemp())
    stderr = P.communicate()[1].splitlines()
    if not stderr[0].startswith('Could not locate'):
        for line in stderr:
            if line.startswith('Opening'):
                index_bt2 = line[(1 + line.find('"')):line.rfind('"')]
                index_basename = index_bt2[0:index_bt2.find('.1.bt2')]
                return index_basename
    for d in [getcwd(), os.path.split(path_to_bowtie2)[0],
              join(os.path.split(path_to_bowtie2)[0], 'indexes')]:
        rprime = join(d, r)
        args = [path_to_bowtie2 + '-inspect', '-v', '-s', rprime]
        debug(' '.join(args))
        P = Popen(args, stdout=open(devnull, 'w'), stderr=PIPE, cwd=mkdtemp())
        stderr = P.communicate()[1].splitlines()
        if not stderr[0].startswith('Could not locate'):
            for line in stderr:
                if line.startswith('Opening'):
                    index_bt2 = line[(1 + line.find('"')):line.rfind('"')]
                    index_basename = index_bt2[0:index_bt2.find('.1.bt2')]
                    return index_basename
    return None

def cat_counter_references(counter_references=None, target_dir=curdir,
                            path_to_bowtie2='bowtie2',
                            logger=None, **kwargs):
    if counter_references is None: return
    try: makedirs(target_dir, mode=0755)
    except OSError: pass
    debug('Validating counter-references and building counter-reference index')
    valid_references = validate_references(references=counter_references,
                                           target_dir=target_dir,
                                           path_to_bowtie2=path_to_bowtie2,
                                           logger=logger,
                                           environ_key='SOT_DEFAULT_COUNTER_REFERENCES')
    crefs_fa = open(join(target_dir, 'counter_references.fa'), 'w')
    for ref in counter_references:
        Popen([path_to_bowtie2 + '-inspect', ref], stdout=crefs_fa).wait()
    crefs_index = join(target_dir, counter_references)
    args = [path_to_bowtie2 + '-build', crefs_fa, crefs_index]
    P = Popen(args, stderr=PIPE)
    stderr = P.communicate()[1]
    if stderr.startswith('Error'):
        critical(stderr)
        critical('No counter-references will be used.')
    return crefs_index

def validate_references(references=None, path_to_bowtie2='bowtie2',
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
        bowtie2_index = find_bowtie2_index(r, path_to_bowtie2=path_to_bowtie2)
        if bowtie2_index is None:
            if exists(r):
                debug('Attempting to build bowtie2 index from %s' % r)
                new_index = fasta_to_bowtie2(r, target_dir=target_dir,
                                            path_to_bowtie2=path_to_bowtie2)
                if new_index is not None:
                    new_references.append(new_index)
                    continue
                else:
                    critical('Failed to build bowtie2 index.')
            critical('bowtie2 could not find the index for %s', r)
            critical('we will not align to %s', r)
        else:
            new_references.append(bowtie2_index)
    return new_references

def convert_to_fastq(fp_obj, logger=None):
    """Convert a SAM or BAM file to FASTQ file(s) for alignment
    """
    input_file = fp_obj.input_file
    output_dir = fp_obj.output_dir
    fastq_dir = join(output_dir, 'fastq_input')
    fp_obj.check_output_dir(fastq_dir)
    protoname = fp_obj.protoname
    if fp_obj.paired_end:
        fastq_filenames = (join(fastq_dir, '%s.1.txt.gz' % protoname),
                           join(fastq_dir, '%s.2.txt.gz' % protoname))
        logger.info('Converting file %s to FASTQ files %s, %s',
                    input_file, fastq_filenames[0], fastq_filenames[1])
        in_args = [sys.executable, '-m', 'seriesoftubes.converters.bamtofastq2',
                   '--gzip', input_file, fastq_filenames[0], fastq_filenames[1]]
    else:
        fastq_filename = join(fastq_dir, '%s.txt.gz' % protoname)
        logger.info('Converting file %s to FASTQ file %s',
                    input_file, fastq_filename)
        in_args = [sys.executable, '-m', 'seriesoftubes.converters.bamtofastq2',
                   input_file, fastq_filename]
    logger.debug('Launching %s', ' '.join(in_args))
    polledpipe = PolledPipe(logger=logger, level=logging.ERROR)
    job = Popen(in_args, stdout=polledpipe.w, stderr=STDOUT)
    wait_for_job(job, [polledpipe], logger)
    if fp_obj.paired_end:
        logger.debug('Settings input_file to %s', fastq_filenames[0])
        fp_obj.input_file = fastq_filenames[0]
        logger.debug('Settings second_file to %s', fastq_filenames[1])
        fp_obj.second_file = fastq_filenames[1]
    else:
        logger.debug('Settings input_file to %s', fastq_filename)
        fp_obj.input_file = fastq_filenames[0]
    logger.debug('Setting use_pysam to False')
    fp_obj.use_pysam = False
    logger.debug('Setting format to FASTQ')
    fp_obj.format = 'FASTQ'
    logger.debug('Ignoring open_func, it will not be used')
    if not job.returncode == 0:
        logger.critical('Conversion FAILED!')
    else:
        logger.info('Conversion successful')
    return
    
@exit_on_Usage
def align2(fp_obj, references=[], counter_references=None,
          unique=True, seed_len='28',
          use_quality=True, logging_level=10, num_cpus=1, **kwargs):
    common_flags = ['--time', '-p', str(num_cpus), '-L', seed_len]
    if fp_obj.paired_end: common_flags.extend(['-X', '600'])
    
    logger = get_logger(logging_level)
    if references is None or references == []:
        logger.critical('Nothing to do')
        return 
    stdout_buffer = []
    
    if fp_obj.format == 'BAM' or fp_obj.format == 'SAM':
        convert_to_fastq(fp_obj, logger=logger)
    
    if not fp_obj.format == 'FASTQ':
        logger.critical('%s only supports FASTQ files' % __file__)
        return
    
    if counter_references is not None:
        #counter align first
        flags = [item for item in common_flags]
        
        flags.append('--fast')
        new_filenames = counteralign_once(fp_obj, flags, counter_references,
                                          logger=logger, **kwargs)
        # after alignment
        fp_obj.input_file = new_filenames[0]
        fp_obj.second_file = new_filenames[1]
    
    flags = [item for item in common_flags]
    
    for ref in references:
        s = align_once(fp_obj, flags, ref, logger=logger, **kwargs)
        stdout_buffer.append(s)
    return '\n'.join([s for s in stdout_buffer if s is not None])

def make_paired_name(name1, name2):
    identity = [(not x==y) for x,y in zip(name1,name2)]
    if not sum(identity) == 1: raise ValueError('Not valid names')
    else:
        index = identity.index(True)
        return identity[0:index] + '%' + identity[(index + 1):]

def counteralign_once(**kwargs):
    """Produce counter-alignements"""
    if use_quality:
        if fp_obj.use_pysam: flags.append('--phred33')
        else: flags.append('--phred64')
        
    refname = os.path.basename(ref)
    output_dir, output_file = os.path.split(fp_obj.tmp_filename(refname))
    fp_obj.check_output_dir(output_dir)
    fp_obj.check_output_dir(join(output_dir, 'counteraligned'))
    filename1 = os.path.abspath(fp_obj.input_file)
    second_file = fp_obj.second_file
    if second_file is not None: filename2 = os.path.abspath(second_file)
    else: filename2 = None
    
    if fp_obj.paired_end:
        try:
            paired_file = make_paired_name(input_file, second_file)
            counteraligned = os.path.abspath(join(output_dir, 'counteraligned',
                                                  paired_file))
        except ValueError:
            counteraligned = os.path.abspath(join(output_dir, 'counteraligned',
                                                  input_file))
        file_args = ['-x', ref, '-1', filename1, '-2', filename2,
                     '--al-conc-gz', counteraligned,
                     '--un-conc-gz', join(output_dir, paired_file)]
        new_filenames = (join(output_dir, input_file),
                         join(output_dir, second_file))
    else:
        file_args = ['-x', ref, '-U', filename1,
                     '--al-gz', join(output_dir, 'counteraligned', input_file),
                     '--un-gz', join(output_dir, input_file)]
        new_filenames = (join(output_dir, input_file), None)
        
    bowtie2_args = [path_to_bowtie2] + flags + file_args
    
    # finish parsing input here
    bowtie2_stdout = PolledPipe(logger=logger, level=logging.ERROR)
    bowtie2_stderr = PolledPipe(logger=logger, level=logging.ERROR)
    logger.info('Launching bowtie2 (output will be piped to samtools for BAM encoding)')
    logger.info(' '.join(bowtie2_args))
    bowtie2_aligner = Popen(bowtie2_args, stdout=open(devnull, 'w'), stderr=bowtie2_stderr.w,
                           bufsize=-1)
    logger.info(' '.join(in_args))
    logger.info('counteraligned reads will be saved as GZIPed FASTQ files in counteraligned/')
    
    logger.debug('Waiting for bowtie2 to finish')
    pollables = [bowtie2_stderr]
    wait_for_job(bowtie2_aligner, pollables, logger)
    
    if not bowtie2_aligner.returncode == 0:
        logger.critical("bowtie2 did not run properly [%d]",
                        bowtie2_aligner.returncode)
        return
    
    logger.debug('Alignment successfully completed')
        
    return new_filenames

def align_once(fp_obj, flags, ref, use_quality=False,
               path_to_bowtie2=None, path_to_samtools=None, logger=None,
               **kwargs):
    if use_quality:
        if fp_obj.use_pysam: flags.append('--phred33')
        else: flags.append('--phred64')
        
    refname = os.path.basename(ref)
    path_to_unsorted = fp_obj.tmp_filename(refname)
    output_dir = os.path.split(path_to_unsorted)[0]
    fp_obj.check_output_dir(output_dir)
    filename1 = os.path.abspath(fp_obj.input_file)
    second_file = fp_obj.second_file
    if second_file is not None: filename2 = os.path.abspath(second_file)
    else: filename2 = None
    
    if fp_obj.paired_end:
        file_args = ['-x', ref, '-1', filename1, '-2', filename2]  
    else:
        file_args = ['-x', ref, '-U', filename1]
  
    bowtie2_args = [path_to_bowtie2] + flags + file_args
    
    # finish parsing input here
    bowtie2_stderr = PolledPipe(logger=logger, level=logging.ERROR)
    logger.info('Launching bowtie2 (output will be piped to samtools for BAM encoding)')
    logger.info(' '.join(bowtie2_args))
    bowtie2_aligner = Popen(bowtie2_args, stdout=PIPE, stderr=bowtie2_stderr.w,
                           bufsize=-1)
    
    samtools_args = [path_to_samtools, 'view', '-b', '-S', '-o',
                     path_to_unsorted, '-']
    logger.info('Launching samtools to encode bowtie2 output as BAM')
    logger.info(' '.join(samtools_args))
    samtools_stdout = PolledPipe(logger=logger, level=logging.WARN)
    samtools_stderr = PolledPipe(logger=logger, level=logging.ERROR)
    samtools_viewer = Popen(samtools_args, stdin=bowtie2_aligner.stdout,
                            stdout=samtools_stdout.w,
                            stderr=samtools_stderr.w, bufsize= -1)
    
    logger.debug('Waiting for bowtie2 to finish')
    pollables = [bowtie2_stderr, samtools_stdout, samtools_stderr]
    wait_for_job(bowtie2_aligner, pollables, logger)
    
    if not bowtie2_aligner.returncode == 0:
        logger.critical("bowtie2 did not run properly [%d]",
                        bowtie2_aligner.returncode)
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

if __name__ == "__main__": main()

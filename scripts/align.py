#!/usr/bin/env python
"""
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
from os import getcwd
from os.path import join
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
                        default=path_to_executable('bowtie', '/usr/local/bowtie-*'),
                        help='The path to the bowtie executable')
    parser.add_argument('--path-to-samtools', nargs='?',
                        default=path_to_executable('samtools', '/usr/local/samtools*'),
                        help='The path to the samtools executable')
    # fix aliases, should be --ref too
    parser.add_argument('--references', nargs='*',
                        default=['hg19', 'hg18'],
                        help='Reference genomes to align against (requires the appropriate bowtie index)',
                        )
    parser.add_argument('--no-unique', dest='unique', action='store_false',
                        help='do not produce unique/ alignment folder') 
    parser.add_argument('--no-random', dest='random', action='store_false',
                        help='do not produce random/ alignment folder') 
    parser.add_argument('--ignore-quality', dest='use_quality',
                        action='store_false',
                        help='Use -v mode with bowtie, allows only n mismatches total')
    parser.add_argument('--mismatches', default='2',
                        help="allow n mismatches, in the seed (default) or total if --ignore-quality")
    parser.add_argument('--quals-type', default='solex1.3',
                        choices=['solexa', 'solexa1.3', 'phred64', 'phred33',
                                 'integer'],
                        help='Valid options are integer, solexa1.3, solexa, phred33, or phred64 (see bowtie for more info)')
    parser.add_argument('--max-quality', default='70',
                        help='specify maximum quality scores of all mismatched positions (default is 70), ignored in -v mode')
    parser.add_argument('--seed-length', dest='seed_len', default='28',
                        help='use seed length of m (default is 28)')
    context = e.get_context()
    new_references = validate_references(**context)
    e.update_context({'references': new_references})
    e.do_action(align)

def validate_references(references=None, path_to_bowtie='bowtie',
                        logger=None, **kwargs):
    debug('Validating references')
    new_references = []
    for r in references:
        args = [path_to_bowtie, r, 'foo_does_not_exist.fakefile']
        P = Popen(args, stderr=PIPE, cwd=mkdtemp())
        err_msg = P.communicate()[1].splitlines()[0]
        if err_msg.find('Bowtie index') != -1:
            # try current working directory, complain on failure
            rprime = join(getcwd(), r) 
            args = [path_to_bowtie, rprime, 'foo_does_not_exist.fakefile']
            P = Popen(args, stderr=PIPE, cwd=mkdtemp())
            err_msg = P.communicate()[1].splitlines()[0]
            if err_msg.find('Bowtie index') != -1:
                critical('bowtie could not find the index for %s', r)
                critical('we will not align to %s', r)
            else:
                new_references.append(rprime)
        else:
            new_references.append(r)
    return new_references

@exit_on_Usage
def align(fp_obj, references=[], random=True, unique=True, max_quality='70',
          quals_type='solexa1.3',  mismatches='2', seed_len='28',
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
    for match_type, match_flag in uniqueness.items():
        flags = [item for item in common_flags]
        flags.extend(match_flag)
        
        source = fp_obj.fastq_source
        # In case we're from that bad day...
        if source=='081124_HWI-EAS355_0001_Meghan' and match_type=='unique':
            flags.extend(['-v', '3'])
        elif use_quality:
            flags.extend(['-e', max_quality])
            flags.extend(['-n', mismatches])
        else:
            flags.extend(['-v', mismatches])

        if fp_obj.paired_end: flags.extend(['-X','600'])
        
        for ref in references:
            s = align_once(fp_obj, flags, ref, match_type, logger=logger, **kwargs)
            stdout_buffer.append(s)
    return '\n'.join([s for s in stdout_buffer if s is not None])

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
    logger.info('Launching botwie (output will be piped to samtools)')
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
                            stdout=sys.stdout,
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
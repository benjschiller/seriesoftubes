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
import gzip
import bz2
import pysam
import logging
import os
import platform
import select
import sys
import time
from subprocess import Popen, PIPE
import scripter
from scripter import assert_path, path_to_executable, Usage, \
                     exit_on_Usage, InvalidFileException, get_logger, \
                     Environment, FilenameParser

if platform.system() == 'Windows':
    raise RuntimeError("""Microsoft Windows is not and will not be supported.
This is due to an underlying design problem in the OS.
Use cygwin if you need to run this on Windows.""")
    
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
    # fix aliases elsewhere
    parser.add_argument('--references', aliases=['--ref'], nargs='*',
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
    parser.add_argument('--seed-length', default='28',
                        help='use seed length of m (default is 28)')
    e.do_action(align)

        
class PolledPipe(object):
    """
    A PolledPipe object has two attributes
    r - pipe read file descriptor
    w - pipe write file descriptor
    and three methods
    poll() -- poll the read end of the pipe
    readlines() -- read availble lines if the poll says the pipe is ready
    log(level=logging.error) -- emit the result of readline (if not None)
    
    a select.poll() object has r registered to it
    """
    def __init__(self, logger=None, level=None):
        self._logger = logger
        self._level = level
        r, w = os.pipe()
        self.r = r
        self._r_file = os.fdopen(r, 'r', 0)
        self.w = w
        self._poll = select.poll()
        self._poll.register(r)
        return
    
    def poll(self, timeout=0):
        self._poll.poll(timeout)
        
    def readlines(self, timeout=0):
        results = self.poll(timeout)
        if results is None: raise StopIteration
        for result in results:
            event = result[1]
            if (event & select.POLLERR) == select.POLLERR:
                raise IOError('Something went wrong trying to read from a pipe')
            elif (event & select.POLLNVAL) == select.POLLNVAL:
                raise IOError('Invalid request: descriptor for pipe not open')
            elif (event & select.POLLHUP) == select.POLLHUP:
                raise IOError('Pipe file descriptor already hung up')
            elif (event & select.POLLIN) == select.POLLIN or \
                 (event & select.POLLPRI) == select.POLLPRI:
                yield self._r_file.readline()
        raise StopIteration
    
    def log(self, level=None):
        """
        emit the results of readlines
        """
        if self._logger is None:
            raise RuntimeWarning('PolledPipe did not have a logger attached')
        if level is None: level = self._level
        if level is None: level = logging.ERROR
        for line in self.readlines():
            self._logger.log(level, line)
        return

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
    path_to_unsorted = fp_obj.tmp_filename(ref, match_type)
    output_dir = os.path.split(path_to_unsorted)[0]
    fp_obj.check_output_dir(output_dir)
    filename1 = os.path.abspath(fp_obj.input_file)
    second_file = fp_obj.second_file
    if second_file is not None: filename2 = os.path.abspath(second_file)
    else: filename2 = None
    
    # finish parsing input here
    input_stderr = PolledPipe(logger=logger, level=logging.ERROR)
    if fp_obj.use_pysam:
        logger.info('Automagically interpreting SAM/BAM file')
        in_args = [sys.executable, '-m', 'seriesoftubes.converters.bamtotab',
                   filename1]
        logger.info(' '.join(in_args))
        input_stream = Popen(in_args, stdout=PIPE, stderr=input_stderr.w,
                             bufsize=-1)
    else:
        compression = fp_obj.compression
        logger.debug('Automagically interpreting FASTQ file (compression is %s)',
                     compression)
        in_args = [sys.executable, '-m',
                   'seriesoftubes.converters.fastqtotab',
                   '--compression', compression, filename1]
        if filename2 is not None:
            in_args.append(filename2)
        logger.info(' '.join(in_args))
        input_stream = Popen(in_args, stdout=PIPE, stderr=input_stderr.w,
                             bufsize=-1)
    
    if use_quality:
        if fp_obj.use_pysam: flags.append('--phred33-quals')
        else: flags.append(''.join(['--', quals_type, '-quals']))
    file_args = [ref, '--12', '-']
    bowtie_args = [path_to_bowtie] + flags + file_args
    logger.info('Launching botwie (output will be sent to PIPE)')
    logger.info(' '.join(bowtie_args))
    bowtie_stderr = PolledPipe(logger=logger, level=logging.ERROR)
    bowtie_aligner = Popen(bowtie_args, stdin=input_stream.stdout,
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
    wait_for_job(bowtie_aligner,
                 [input_stderr, bowtie_stderr,
                  samtools_stdout, samtools_stderr], logger)
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

def wait_for_job(job, logs=[], logger=None):
    otime = time.time()
    while True:
        atime = time.time()
        if otime - atime > 60 and logger is not None:
            logger.debug('Still running')
        otime = atime
        time.sleep(3)
        if job.poll() is not None: break
        for log in logs: log.log()
    for log in logs: log.log()

def get_pair_info(illumina_name):
    """
    take a filename from GERALD output and figure out whether it's the first
    or second read (or single-end read) and return a tuple
    (pair_index, second_file_name, new_output_name)
    """
    name_parts = illumina_name.split('.')
    for i in range(len(name_parts)):
        part = name_parts[i]
        subparts = part.split('_') 
        if len(subparts) > 0:
            if subparts[0] == 's':
                if len(subparts) == 1: continue
                pair_index = subparts[2]
                # lane = subparts[1]
                if pair_index == '1':
                    join_part = '_'.join(subparts[0:2] + subparts[3:])
                    new_output_name = '.'.join(name_parts[0:i] + [join_part] +\
                                                name_parts[i+1:])
                    second_part = '_'.join(subparts[0:2] + ['2'] + subparts[3:])
                    second_name = '.'.join(name_parts[0:i] + [second_part] + 
                                           name_parts[i+1:])
                    return (pair_index, second_name, new_output_name)
                elif pair_index == '2':
                    raise InvalidFileException
                    return
    return None

class BowtieFilenameParser(FilenameParser):
    def __init__(self, filename, *args, **kwargs):
        super(BowtieFilenameParser, self).__init__(filename, *args, **kwargs)
        compression, format = self.discover_file_format2(filename)
        if format == 'SAM' or format == 'BAM':
            self.use_pysam = True
            # try to open the file so we're sure it works
            f = pysam.Samfile(filename)
            del f
        elif format == 'FASTQ':
            self.use_pysam = False
            self.compression = compression
        else:
            raise RuntimeError('Dubious file format')
        self.second_file = None
        self.check_paired_end()
        if len(self.protoname.split('.')) > 6:
            self.fastq_source = self.protoname.split('.')[6]
        else:
            self.fastq_source = 'Unknown'
    
    @staticmethod
    def discover_file_format2(filename):
        f = open(filename, 'rb')
        head = f.read(3)
        f.close()
        # check magic words for compression
        if head == '\x1f\x8b\x08':
            compression = 'gz'
            open_func = gzip.GzipFile
        elif head=='\x42\x5a\x68':
            open_func = bz2.BZ2File
            compression = 'bz2'
        else:
            open_func = open
            compression = 'none'
        uncompressed = open_func(filename)
        head2 = uncompressed.read(4)
        # check for BAM
        if head2 == 'BAM\x01': return (compression, 'BAM')
        # check for SAM 
        head2b = uncompressed.readline()
        if head2 == '@HD\t':
            if head2b[0:3] == 'VN:': return (compression, 'SAM')
        # check FASTQ
        title = head2 + head2b
        seq = uncompressed.readline()
        title2 = uncompressed.readline()
        qual = uncompressed.readline()
        if len(seq) == len(qual) and \
           (title.startswith('@') or title2.startswith('+')) and \
           (title2[1:].strip() == '' or title[1:] == title2[1:]):
            return (compression, 'FASTQ')
        # otherwise give up
        else: return (None, None)
            
    def check_paired_end(self):
        # check if this is a paired-end file
        # if so, grab its partner
        seqfile_name = os.path.basename(self.input_file)
        pair_info = get_pair_info(seqfile_name)
        if pair_info is not None:
            pair_index = pair_info[0]
            second_name = pair_info[1]
            new_name = pair_info[2]
            scripter.debug('NOTICE: Detected paired read file.')
            if pair_index == '1':
                scripter.debug('Attempting to find second file.')

                self.second_file = os.sep.join([self.input_dir, second_name])
                self.protoname = os.path.splitext(new_name)[0]
                scripter.debug('Found %s', self.second_file)
                try:
                    assert_path(self.second_file)
                    self.paired_end = True
                except IOError:
                    scripter.debug('Failed to find paired end file')
                    self.paired_end = False
            elif pair_index == '2':
                scripter.debug('This is the second file, ignoring it.')
                raise InvalidFileException
            else:
                scripter.debug('Failed to find paired end')
                self.paired_end = False
        else: self.paired_end = False

    def tmp_filename(self, ref, match_type):
        return os.path.join(self.output_dir, ref, match_type,
                                       self.with_extension('tmp'))

if __name__=="__main__": main()
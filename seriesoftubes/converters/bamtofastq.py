"""Convert BAM/SAM to FASTQ format*
@name
sequence
+name
quality score (phred33)

files may be SAM or BAM (autodetected)
If the file(s) contain paired-end sequences, we will write to two files
    (in the current working directory)
If the files contain single end sequences, we will write to stdout by default
Output is always to stdout (err goes to stderr, redirect it if you need to)
"""
import pysam
import os
import select
from scripter import path_to_executable 
from argparse import ArgumentParser
from copy import copy
from subprocess import Popen, PIPE
from os import mkfifo, getcwd, devnull
from os.path import join, exists, abspath
import subprocess
from sys import argv, stdin, stdout, stderr, exit, executable
from gzip import GzipFile
from .discover import PATH_TO_GZIP, gzip_class_factory

class UnpairedBAMToFastqConverter(object):
    """Works with unpaired SAM/BAM file"""
    def __init__(self, file_, wd=None, stderr=None, logger=None):
        self.require_bam(file_)
        if wd is None: wd = getcwd()
        fifofile = join(wd, '._sot_fifo')
        if exists(fifofile):
            raise RuntimeError('%s already exists' % fifofile)
        mkfifo(fifofile)
        self.fifo = abspath(fifofile)
        self._args = [executable, '-m', 'seriesoftubes.converters.bamtofastq',
                   file_, '--single-stdout', fifofile]
        self._logger = logger
        self._stderr = stderr
    
    def launch(self):
        if self._logger is not None:
            self._logger.info('Launching %s', ' '.join(self._args))
        self.subprocess = Popen(self._args, stdout=open(devnull, 'w'),
                                stderr=self._stderr, bufsize=-1)
        
    def get_fifo_readers(self):
        return (open(self.fifos[0], 'r'), open(self.fifos[1], 'r'))
    
    def require_bam(self, filename):
        with open(filename, 'rb') as f:
            head = f.read(3)
        # check magic words for compression
        if head == '\x1f\x8b\x08':
            open_func = GzipFile
        else:
            open_func = open
        uncompressed = open_func(filename)
        head2 = uncompressed.read(4)
        # check for BAM
        if head2 == 'BAM\x01': return
        # check for SAM 
        if head2 == '@HD\t': return
        else: raise ValueError('Not a SAM/BAM file') 

def main():
    """
    what to do if we execute the module as a script
    
    bamtofastq can only convert files (not stdin) because of the paired-end problem
    """
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('files', nargs='+', help='List of input files')
    parser.add_argument('--no-gzip',
                        help='Do not compress output')
    parser.add_argument('--single-stdout',
                        help='Save single-end reads to here (default: stdout)',
                        default=stdout)
    args = parser.parse_args()
    context = vars(args)
    read_files(**context)
    
def pair_writer(out1, out2):
    def writer(read1, read2):
        out1.write('@%s\n%s\n+\n%s\n' % (read1.qname, read1.seq, read1.qual))
        out2.write('@%s\n%s\n+\n%s\n' % (read2.qname, read2.seq, read2.qual))
    return writer

def read_files(files=None, no_gzip=False, single_stdout=stdout):
    """
    actually reads the SAM/BAM files
    """
    for file in files:
        if file is None: continue
        f = pysam.Samfile(file)
        #check if first read is paired
        aread = f.next()
        f.close()
        f = pysam.Samfile(file)
        if aread.is_paired:
            if no_gzip:
                print 'Detected paired-end reads, redirecting output to text files'
                file1 = file + '_1.txt'
                file2 = file + '_2.txt'
                fh1 = open(file1, 'w')
                fh2 = open(file2, 'w')
            elif PATH_TO_GZIP is not None:
                print 'Detected paired-end reads, redirecting output to .gz text files (using system gzip)'
                file1 = file + '_1.txt.gz'
                file2 = file + '_2.txt.gz'
                open_func = gzip_class_factory(PATH_TO_GZIP)
                fh1 = open_func(file1, 'wb')
                fh2 = open_func(file2, 'wb')
            else:
                print 'Detected paired-end reads, redirecting output to .gz text files'
                file1 = file + '_1.txt.gz'
                file2 = file + '_2.txt.gz'
                fh1 = GzipFile(file1, 'wb')
                fh2 = GzipFile(file2, 'wb')
            is_paired = False
            write = pair_writer(fh1, fh2)
            incomplete_pairs = []
            for aread in f:
                is_paired = False
                qname = aread.qname
                for i in xrange(len(incomplete_pairs)):
                    if incomplete_pairs[i].qname == qname:
                        mate_read = incomplete_pairs.pop(i)
                        # figure out order
                        if aread.flag & 0x4 == 0x4:
                            write(aread, mate_read)
                        else:
                            write(mate_read, aread)
                        is_paired = True
                        break
                if not is_paired: incomplete_pairs.append(aread)
            unpaired = len(incomplete_pairs)
            if not unpaired == 0:
                raise RuntimeError('%d unpaired reads remaining' % unpaired) 
        else:
            if no_gzip:
                open_func = open
            elif PATH_TO_GZIP is not None:
                open_func = gzip_class_factory(PATH_TO_GZIP)
            else:
                open_func = GzipFile
            fh1 = open_func(file1, 'wb')
            for aread in f:
                qname = aread.qname or ''
                seq = aread.seq or ''
                qual = aread.qual or ''
                rec = '@%s\n%s\n+\n%s\n' % (qname, seq, qual)
                fh1.write(rec)
            fh1.close()
    exit(0)        
        
if __name__ == '__main__': main()
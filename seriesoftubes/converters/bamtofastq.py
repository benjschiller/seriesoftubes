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
import subprocess
from sys import argv, stdin, stdout, stderr, exit
from gzip import GzipFile
from .discover import PATH_TO_GZIP

def main():
    """
    what to do if we execute the module as a script
    
    bamtotab can only convert files (not stdin) because of the paired-end problem
    """
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('files', nargs='+', help='List of input files')
    parser.add_argument('--no-gzip',
                        help='Do not compress paired-end output files')
    parser.add_argument('--no-stdout',
                        help='Save single-end reads to text files too')
    args = parser.parse_args()
    context = vars(args)
    read_files(**context)

def read_files(files=None, no_gzip=False, no_stdout=False):
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
                file2 = files + '_2.txt'
                fh1 = open(file1, 'w')
                fh2 = open(file2, 'w')
            elif PATH_TO_GZIP is not None:
                print 'Detected paired-end reads, redirecting output to .gz text files (using system gzip)'
                file1 = file + '_1.txt.gz'
                file2 = files + '_2.txt.gz'
                open_func = gzip_class_factory(PATH_TO_GZIP)
                fh1 = open_func(file1, 'wb')
                fh2 = open_func(file2, 'wb')
            else:
                print 'Detected paired-end reads, redirecting output to .gz text files'
                file1 = file + '_1.txt.gz'
                file2 = files + '_2.txt.gz'
                fh1 = GzipFile(file1, 'wb')
                fh2 = GzipFile(file2, 'wb')
            print file1
            print file2
            for aread in f:
                qname = aread.qname or ''
                seq = aread.seq or ''
                qual = aread.qual or ''
                rec = '@{0!s}\n{1!s}\n+{0!s}\n{2!s}\n'.format(qname, seq, qual)
                if aread.is_read1: fh1.write(rec)
                elif aread.is_read2: fh2.write(rec)
                else: raise ValueError("This shouldn't happen")
        else:
            fh2 = None
            if no_stdout:
                print 'Redirecting output to text files'
                file1 = file + '.txt'
                fh1 = open(file1, 'w')
            else:
                fh1 = stdout
            for aread in f:
                qname = aread.qname or ''
                seq = aread.seq or ''
                qual = aread.qual or ''
                rec = '@{0!s}\n{1!s}\n+{0!s}\n{2!s}\n'.format(qname, seq, qual)
                fh1.write(rec)
        if fh1 is not None: fh1.close()
        if fh2 is not None: fh2.close()
    exit(0)        
        
if __name__ == '__main__': main()
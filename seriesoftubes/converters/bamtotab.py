import pysam
import os
import select
from scripter import path_to_executable 
from argparse import ArgumentParser
from copy import copy
from subprocess import Popen, PIPE
import subprocess
from sys import argv, stdin, stdout, stderr, exit

def main():
    """
    what to do if we execute the module as a script
    bamtotab can only convert files because of the paired-end problem
    """
    parser = ArgumentParser(description="""Convert BAM/SAM to tab-format
NAME    SEQ1    QUAL1    SEQ2    QUAL2

files may be SAM or BAM (autodetected)
Output is always to stdout (err goes to stderr, redirect it if you need to)""")
    parser.add_argument('files', nargs='+', help='List of input files')
    args = parser.parse_args()
    context = vars(args)
    read_files(context['files'])

def read_files(files):
    for file in files:
        if file is None: continue
        f = pysam.Samfile(file)
        for aread in f:
            qname = aread.qname or ''
            seq = aread.seq or ''
            qual = aread.qual or ''
            if aread.is_paired:
                if aread.is_read2: continue
                whence = f.tell()
                aread2 = f.mate(aread)
                seq2 = aread2.seq or ''
                qual2 = aread2.qual or ''
                f.seek(whence)
                print '%s\t%s\t%s\t%s\t%s' % (qname, seq, qual, seq2, qual2)
            else:
                print '%s\t%s\t%s' % (qname, seq, qual)
    exit(0)        
        
if __name__ == '__main__': main()
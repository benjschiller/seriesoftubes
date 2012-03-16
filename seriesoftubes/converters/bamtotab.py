"""Convert BAM/SAM to tab-format

NAME    SEQ1    QUAL1    SEQ2    QUAL2

files may be SAM or BAM (autodetected)

Output is always to stdout (err goes to stderr, redirect it if you need to)"""
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
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('files', nargs='+', help='List of input files')
    args = parser.parse_args()
    context = vars(args)
    read_files(context['files'])

def read_files(files):
    for file in files:
        if file is None: continue
        f = pysam.Samfile(file)
        last = None
        for aread in f:
            qname = aread.qname or ''
            seq = aread.seq or ''
            qual = aread.qual or ''
            if last is not None:
                name, seq1, qual1 = last
                print '%s\t%s\t%s\t%s\t%s' % (name, seq1, qual1, seq, qual)
                last = None
            elif aread.is_paired:
                last  = (qname, seq, qual) 
            else:
                print '%s\t%s\t%s' % (qname, seq, qual)
    exit(0)        
        
if __name__ == '__main__': main()
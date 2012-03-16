"""Convert single-end or paired-end sequencing files to tab-format
NAME    SEQ1    QUAL1    (SEQ2    QUAL2)

Supports FASTQ (plaintext, gzip and bz2), SAM, BAM
Output is always to stdout (err goes to stderr, redirect it if you need to)
"""
from argparse import ArgumentParser
import bz2
import gzip
import bamtotab, fastqtotab
from .discover import discover_file_format

def main():
    """
    what to do if we execute the module as a script
    """
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('file_read1', help='File with either single or both paired reads or the first of paired reads')
    parser.add_argument('file_read2', nargs='?', help='(optional) File with second of paired reads')
    args = parser.parse_args()
    context = vars(args)
    file1 = context['file_read1']
    file2 = context['file_read2']
    open_func, format = discover_file_format(file1)
    if format in ['BAM', 'SAM']: bamtotab.read_files([file1])
    elif format == 'FASTQ': fastqtotab.read_files(file1, file2, open_func)
    else: raise RuntimeWarning('Dubious file format') 

if __name__ == '__main__': main()
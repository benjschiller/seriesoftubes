"""Convert single or pairs of FASTQ files to tab-format

NAME    SEQ1    QUAL1    (SEQ2    QUAL2)

Output is always to stdout (err goes to stderr, redirect it if you need to)
"""
from argparse import ArgumentParser
from sys import exit
from .discover import discover_file_format
from itertools import izip

def main():
    """
    what to do if we execute the module as a script
    
    
    fastqtotab can only convert files because of the paired-end problem
    """
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('fastq_file_read1', help='File with either single reads or the first of paired reads')
    parser.add_argument('fastq_file_read2', nargs='?', help='(optional) File with second of paired reads')
    args = parser.parse_args()
    context = vars(args)
    file1 = context['fastq_file_read1']
    file2 = context['fastq_file_read2']
    open_func = discover_file_format(file1)[0]
    read_files(file1, file2, open_func)

def read_files(file1, file2, open_func):
    """
    We do not validate records here
    """
    fh = open_func(file1)
    i = 0
    if file2 is not None:
        fh2= open_func(file2)
        # rec = [title, seq, qual, seq2, qual]
        rec = [None, None, None, None, None]
        line_fmt = '{0!s}\t{1!s}\t{2!s}\t{3!s}\t{4!s}'
        for line, line2 in izip(fh, fh2):
            if i == 0: rec[0] = line[1:].rstrip()
            elif i == 1:
                rec[1] = line.rstrip()
                rec[3] = line2.rstrip()
            elif i == 2: pass
            elif i == 3:
                rec[2] = line.rstrip()
                rec[4] = line2.rstrip()
                print line_fmt.format(*rec)
                i = -1
            i += 1
    else:
        # rec = [title, seq, qual]
        rec = [None, None, None]
        line_fmt = '{0!s}\t{1!s}\t{2!s}'
        for line in fh:
            if i == 0: rec[0] = line[1:].rstrip()
            elif i == 1: rec[1] = line.rstrip()
            elif i == 2: pass
            elif i == 3:
                rec[2] = line.rstrip()
                print line_fmt.format(*rec)
                i = -1
            i += 1
    exit(0)        
        
if __name__ == '__main__': main()
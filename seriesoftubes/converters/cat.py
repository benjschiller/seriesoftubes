"""Writes file(s) to stdout, decompress them if needed
supports gzip, bzip2

Output is always to stdout (err goes to stderr, redirect it if you need to)
"""
from argparse import ArgumentParser
from copy import copy
from sys import argv, stdin, stdout, stderr, exit
from .discover import discover_file_format
from itertools import izip

def main():
    """
    what to do if we execute the module as a script
    
    cat
    """
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('files', nargs='+', help='files')
    args = parser.parse_args()
    context = vars(args)
    files = context['files']
    for file in files:
        open_func = discover_file_format(file)[0]
        for line in open_func(file):
            print line
    exit(0)        
        
if __name__ == '__main__': main()
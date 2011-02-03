#!/usr/bin/env python
'''
Finds and sequences errors (Ns) at each base position in SAM files

DEPRECATED : will be replaced with detect_bam_errors.py

Default options: --target=errors.SAM
'''
raise UserWarning('This script is broken right now')
import os
import collections
from collections import defaultdict
import scripter
try:
    import pysam
except ImportError:
    detect_errors = detect_errors_nopysam

VERSION = "2.4"

def main():
    e = scripter.Environment(version=VERSION, doc=__doc__)
    e.set_source_dir = 'alignments.SAM'
    e.set_target_dir = 'errors.SAM'
    e.set_filename_parser(FilenameParser)
    e.do_action(detect_errors)

def parse_MD_value(aligned_read):
    #windows hack
    s = None
    s = dict(aligned_read.tags())['MD']
    if s is None: return None
    start = 0  
    end = 0
    mpos = 0
    mismatched_bases = []
    slen = len(s)
    while True:
        if start is slen:
            break
        # if we're at a letter
        if s[start].isalpha():
            mismatched_bases.append( mpos )
            start += 1
            end += 1
            mpos += 1
            continue
        # if we're at a number
        while True:
            if end is slen: break
            elif s[end].isdigit(): end += 1
            else: break
        number = s[start:end]
        start = end
        mpos += int(number)
    if is_reverse_direction(line):
        return [seq_length - x + 1 for x in mismatched_bases]
    else:
        return mismatched_bases

def is_reverse_direction(line):
    '''True if the read is mapped to the minus strand
    requires the file follows SAM specification
    '''
    return line.split()[1] is '16'

def detect_errors(parsed_filename, verbose=False, **kwargs):
    if verbose:
        stdout_buffer = ''.join(['Determining basewise error rates in ',
                                 parsed_filename.input_file])
    else:
        stdout_buffer = ''
    
    pysam_file = pysam.Samfile(parsed_filename.input_file)
    errors = defaultdict(lambda: 0)
    first_good_line = True
    max_seq_length = 0
    i = 0

    for line in pysam_file:
        i += 1
        seq_length = len(line.split()[9])
        if seq_length > max_seq_length: max_seq_length = seq_length
        mismatched_bases = parse_MD_value(line)
        for base in mismatched_bases: errors[base] += 1

    if verbose:
        stdout_buffer = os.linesep.join([stdout_buffer,
                                 ' '.join(['Found', str(i), 'sequences']),
	                             ' '.join(['Writing error counts to',
                                 parsed_filename.output_filename]), ''])
    newfile = open(parsed_filename.output_filename,'w')
    for pos in range(max_seq_length):
        count = errors[pos]
        newfile.write(''.join(['\t'.join([str(pos), str(count)]), os.linesep]))
    newfile.close()

    return stdout_buffer

def detect_errors_nopysam(parsed_filename, verbose=False, **kwargs):
    raise DeprecatedWarning('This function is deprecated. Install pysam' +
                            'to use the newer version')
    if verbose:
        stdout_buffer = ''.join(['Determining basewise error rates in ',
                                 parsed_filename.input_file])
    else:
        stdout_buffer = ''

    
    f = open(parsed_filename.input_file)
    errors = defaultdict(lambda: 0)
    first_good_line = True
    max_seq_length = 0
    valid_lines = (line for line in f if line[0] is not '@')
    i = 0

    for line in valid_lines:
        i += 1
        seq_length = len(line.split()[9])
        if seq_length > max_seq_length: max_seq_length = seq_length
        mismatched_bases = parse_MD_value(line)
        for base in mismatched_bases: errors[base] += 1

    if verbose:
        stdout_buffer = os.linesep.join([stdout_buffer,
                                 ' '.join(['Found', str(i), 'sequences']),
	                             ' '.join(['Writing error counts to',
                                 parsed_filename.output_filename]), ''])
    newfile = open(parsed_filename.output_filename,'w')
    for pos in range(max_seq_length):
        count = errors[pos]
        newfile.write(''.join(['\t'.join([str(pos), str(count)]), os.linesep]))
    newfile.close()

    return stdout_buffer

class FilenameParser(scripter.FilenameParser):
    def __init__(self, filename, *args, **kwargs):
        super(FilenameParser, self).__init__(filename, *args, **kwargs)
        self.output_filename = os.path.join(self.output_dir,
                         self.with_extension(os.extsep.join(['errors','txt'])))

if __name__== "__main__": main()

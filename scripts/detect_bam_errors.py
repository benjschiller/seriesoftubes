#!/usr/bin/env python
'''
Finds and sequences errors (Ns) at each base position in BAM files

Default options: --target=errors.BAM
'''
import os
import collections
from collections import defaultdict
import scripter
from scripter import Usage, exit_on_Usage
try: import pysam
except ImportError: raise Usage('This script requires pysam')

VERSION = "2.4"

def main():
    e = scripter.Environment(version=VERSION, doc=__doc__)
    e.set_source_dir('alignments.BAM')
    e.set_target_dir('errors.BAM')
    e.set_filename_parser(FilenameParser)
    e.do_action(detect_errors)

def parse_MD_value(aligned_read):
    try: s = dict(aligned_read.tags)['MD']
    except KeyError: return []
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
    if aligned_read.is_reverse:
        return map(lambda x: aligned_read.qlen - x + 1, mismatched_bases)
    else:
        return mismatched_bases

@exit_on_Usage
def detect_errors(parsed_filename, verbose=False, **kwargs):
    stdout_buffer = ''
    stdout_buffer += 'Determining basewise error rates in {!s}\n'.format( 
                                             parsed_filename.input_file)
    seq_length = 0
    errors = defaultdict(lambda: 0)
    num_reads = 0

    with pysam.Samfile(parsed_filename.input_file) as pysam_file:
        for aligned_read in pysam_file:
            num_reads += 1
            seq_length = max(aligned_read.qlen, seq_length)
            mismatched_bases = parse_MD_value(aligned_read)
            for base in mismatched_bases: errors[base] += 1

    if verbose:
        stdout_buffer += 'Found {!s} sequences\n'.format(num_reads)  
        stdout_buffer += 'Writing error counts to{!s}\n'.format(
                                        parsed_filename.output_filename)

    with open(parsed_filename.output_filename, 'w') as newfile:
        for pos in xrange(seq_length):
            newfile.write('{!s}\t{!s}\n'.format(pos, errors[pos]))

    return stdout_buffer

class FilenameParser(scripter.FilenameParser):
    def __init__(self, filename, *args, **kwargs):
        super(FilenameParser, self).__init__(filename, *args, **kwargs)
        self.output_filename = os.path.join(self.output_dir,
                         self.with_extension(os.extsep.join(['errors','txt'])))

if __name__== "__main__": main()
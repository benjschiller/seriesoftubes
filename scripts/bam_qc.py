#!/usr/bin/env python
'''
acts on bam files
Finds and sequences errors (Ns) at each base position in BAM files

Default options: --target=errors.BAM
'''
import os.path
import collections
from collections import defaultdict
import scripter
from scripter import print_debug, exit_on_Usage
import pysam
VERSION = "2.4"

def main():
    e = scripter.Environment(version=VERSION, doc=__doc__)
    e.set_source_dir('alignments.BAM')
    e.set_target_dir(os.path.join('bam_qc', 'duplicate_counts'))
    e.do_action(calculate_duplicates, stay_open=True)
    e.set_target_dir(os.path.join('bam_qc', 'error_counts'))
    e.do_action(detect_errors)

def calculate_duplicates(parsed_filename, debug=False, **kwargs):
    '''
    note: you must be looking at a sorted file, or this won't work
    and you should not use random alignments, because this will not be accurate
    '''
    bam_file = pysam.Samfile(parsed_filename.input_file, 'rb')

    current_frag = ('', 0, '', 0)
    d = collections.defaultdict(lambda: 0)
    dup_counter = 1
    unaligned = 0

    if debug: print_debug('Calculating duplicate rate for',
                          parsed_filename.input_file)

    for read in bam_file:
        # skip mate if it exists
        if read.is_paired and read.is_read2:
            continue
        # break if we've reached unaligned reads (end of file)
        if read.rname is -1:
            d[dup_counter] += 1
            unaligned = 1
            break

        last_frag = current_frag
        current_frag = (read.rname, read.pos, read.mrnm, read.mpos)

        # check if the FULL fragment is identical
        # if so, increment the duplicate counter
        if current_frag == last_frag:
            dup_counter += 1
        # if not, increment the dictionary and reset the counter
        else:
            d[dup_counter] += 1
            dup_counter = 1

    for read in bam_file: unaligned += 1

    unique_tags = sum(d.values())
    total_tags = sum([k * v for k, v in d.items()])

    output_filename = os.path.join(parsed_filename.output_dir,
                                   parsed_filename.with_extension('txt'))
    
    with open(output_filename, 'w') as output_file:
        output_file.write('#{!s}\n'.format(output_filename))
        for k in sorted(d.keys()):
            values = (k, d[k], float(d[k])/float(total_tags))
            output_file.write('{:3G}\t{!s}\t{:.2E}\n'.format(*values))
        output_file.write('Unique\t{!s}\t{:.2E}\n'.format(unique_tags,
                                     float(unique_tags)/float(total_tags)))
        output_file.write('Total\t{!s}\t1.0\n'.format(total_tags))
        output_file.write('Unaligned\t{!s}\t-\n'.format(unaligned))
    
    return

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

    output_filename = os.path.join(parsed_filename.output_dir,
                                   parsed_filename.with_extension('errors.txt'))
    with open(output_filename, 'w') as newfile:
        for pos in xrange(seq_length):
            newfile.write('{!s}\t{!s}\n'.format(pos, errors[pos]))

    return stdout_buffer

if __name__== "__main__": main()
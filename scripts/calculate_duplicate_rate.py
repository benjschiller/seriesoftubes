#!/usr/bin/env python
'''
acts on bam files
'''
import os.path
import subprocess
import collections
import scripter
import pysam
VERSION = "2.4"

def main():
    e = scripter.Environment(version=VERSION, doc=__doc__)
    e.set_source_dir('alignments.BAM')
    e.set_target_dir('duplicates.rates')
    e.do_action(calculate_duplicates)

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

if __name__=="__main__": main()
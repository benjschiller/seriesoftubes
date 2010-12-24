#!/usr/bin/env python
'''
acts on bam files
'''
import os
import subprocess
import collections
import scripter
import pysam
scripter.SCRIPT_DOC = __doc__
scripter.SCRIPT_VERSION = "2.2"
scripter.SOURCE_DIR = 'alignments.BAM'
scripter.TARGET_DIR = 'duplicate.rates'

def calculate_duplicates(parsed_filename, debug=False, **kwargs):
    '''
    note: you must be looking at a sorted file, or this won't work
    and you should not use random alignments, because this will not be accurate
    '''
    try: bam_file = pysam.Samfile(parsed_filename.input_file, 'rb')
    except: return

    current_frag = ('', 0, '', 0)
    d = collections.defaultdict(lambda: 0)
    dup_counter = 1
    unaligned = 0

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
    output_file = open(output_filename, 'w')
    output_file.write('#')
    output_file.write(output_filename)
    output_file.write(os.linesep)

    for k in sorted(d.keys()):
        output_file.write('\t'.join([str(x) for x in 
                                    [k,d[k],float(d[k])/float(total_tags)]]))
        output_file.write(os.linesep)

    output_file.write('\t'.join(['Unique', str(unique_tags),
                                 str(float(unique_tags)/float(total_tags))]))
    output_file.write(os.linesep)
    output_file.write('\t'.join(['Total', str(total_tags), '1.0']))
    output_file.write(os.linesep)
    output_file.write('\t'.join(['Unaligned', str(unaligned), '-']))
    output_file.write(os.linesep)

    output_file.close()
    stdout_buffer = ' '.join(['Calculated duplicate rate for',
                             parsed_filename.input_file])

    if debug: return stdout_buffer
    else: return ''

if __name__=="__main__":
	scripter.perform(calculate_duplicates)

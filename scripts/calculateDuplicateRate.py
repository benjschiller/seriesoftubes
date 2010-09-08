#!/usr/bin/env python
'''
acts on bam files
'''
import os
import subprocess
import collections
import scripter
scripter.SCRIPT_DOC = __doc__
scripter.SCRIPT_VERSION = "2.1"
scripter.SOURCE_DIR = 'alignments.BAM'
scripter.TARGET_DIR = 'duplicate.rates'

PATH_TO_SAMTOOLS='/usr/local/bin/samtools'

def calculate_duplicates(parsed_filename, **kwargs):
    '''
    note: you must be looking at a sorted file, or this won't work
    and you should not use random alignments, because this will not be accurate
    '''
    if not parsed_filename.file_extension == 'bam': return ''
    contents = subprocess.Popen([PATH_TO_SAMTOOLS, 'view',
                                 parsed_filename.input_file], 
                                 stdout=subprocess.PIPE)
    content_stream = contents.stdout
    return actually_calculate_duplicates(parsed_filename,
                                         content_stream, **kwargs)

def actually_calculate_duplicates(parsed_filename, content_stream, 
                                  debug=False, **kwargs):
    these_coords = ('', '')
    d = collections.defaultdict(lambda: 0)
    dup_counter = 0

    for line in content_stream:
        previous_coords = these_coords
        ref_name = line.split()[2]
        ref_coord = line.split()[3]
        these_coords = (ref_name, ref_coord)
        if ref_name == '*':
            dup_counter = 0
            continue

        if previous_coords == these_coords:
            # increment the duplicate counter
            dup_counter += 1
        else:
            # Check if we were looking at duplicates
            if dup_counter > 0: d[dup_counter] += 1
            # Reset the duplicate counter
            dup_counter = 0
            # Increment the number of 1-time occurrences
            d[1] += 1

    if dup_counter > 0: d[dup_counter] += 1

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

    output_file.close()
    stdout_buffer = ' '.join(['Calculated duplicate rate for',
                             parsed_filename.input_file])

    if debug: return stdout_buffer
    else: return ''

if __name__=="__main__":
	scripter.perform(calculate_duplicates)

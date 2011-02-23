#!/usr/bin/env python
'''
removes PCR duplicates from BAM files using samtools
'''
import os.path
import subprocess
import scripter
import pysam

VERSION = "2.4"

# check for sam vs bam later
def main():
    e = scripter.Environment(version=VERSION, doc=__doc__)
    e.set_source_dir('alignments.BAM')
    e.set_target_dir('unique_tags.BAM')
    e.do_action(remove_duplicates)

def remove_duplicates(parsed_filename, **kwargs):
    if not parsed_filename.file_extension=='bam': return
    stdout_buffer = ''

    path_to_output_file = os.path.join(parsed_filename.output_dir, 
                                       parsed_filename.with_extension('bam'))
    
    # fix this later
    try:
        x = pysam.rmdup('-S', parsed_filename.input_file, path_to_output_file)
        stdout_buffer += x
    except: pass
    try:
        x += pysam.index(path_to_output_file)
        stdout_buffer += x
    except: pass
    
    return stdout_buffer

if __name__=="__main__": main()
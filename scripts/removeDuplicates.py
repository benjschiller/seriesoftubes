#!/usr/bin/env python
'''
removes PCR duplicates from BAM files using samtools
'''
import os
import subprocess
import scripter
scripter.SCRIPT_DOC = __doc__
scripter.SCRIPT_VERSION = "2.1"
scripter.SOURCE_DIR = 'alignments.BAM'
scripter.TARGET_DIR = 'unique_tags.BAM'

PATH_TO_SAMTOOLS = '/usr/local/bin/samtools'

def remove_duplicates(parsed_filename, **kwargs):
    if not parsed_filename.file_extension=='bam': return
    stdout_buffer = ''

    steps = []
    path_to_output_file = os.path.join(parsed_filename.output_dir, 
                                       parsed_filename.with_extension('bam'))
    steps.append([PATH_TO_SAMTOOLS, 'rmdup', '-S', parsed_filename.input_file,
                  path_to_output_file])
    steps.append([PATH_TO_SAMTOOLS, 'index', path_to_output_file])

    for step in steps:
        job = subprocess.Popen(step, stdout=subprocess.PIPE, 
                               stderr=subprocess.STDOUT)
        (stdout_data, stderr_data) = job.communicate()

        if stdout_data.strip()=='':
            stdout_buffer = os.linesep.join([stdout_buffer, '', '', 
                                         ' '.join(step)])
        else:
            stdout_buffer = os.linesep.join([stdout_buffer, '', '', 
                                         ' '.join(step), '', stdout_data])
    return stdout_buffer

if __name__=="__main__":
	scripter.perform(remove_duplicates)

#!/usr/bin/env python
'''
Converts SAM to BAM files
and also builds an index for the BAM files

Default options: --target=alignments.BAM
'''
import os
import subprocess
import scripter
scripter.SCRIPT_DOC = __doc__
scripter.SCRIPT_VERSION = "2.1"
scripter.SOURCE_DIR = 'alignments.SAM'
scripter.TARGET_DIR = 'alignments.BAM'

PATH_TO_SAMTOOLS='/usr/local/bin/samtools'

class FilenameParser(scripter.FilenameParser):
    def __init__(self, filename, *args, **kwargs):
        super(FilenameParser, self).__init__(filename, *args,**kwargs)

        f = lambda x: os.path.join(self.output_dir, self.with_extension(x))
        self.output_file = f('bam')
        self.temp_file = f('tmp')

def convert_sam_to_bam(parsed_filename, verbose=False, **kwargs):

    stdout_buffer = ''

    steps = []
    steps.append([PATH_TO_SAMTOOLS, 'view', '-S', '-b', '-o', 
                  parsed_filename.temp_file, parsed_filename.input_file])
    steps.append([PATH_TO_SAMTOOLS, 'sort', parsed_filename.temp_file, 
                  os.path.splitext(parsed_filename.output_file)[0]])
    steps.append([PATH_TO_SAMTOOLS,'index',parsed_filename.output_file])

    for step in steps:
        job = subprocess.Popen(step, stdout=subprocess.PIPE, 
                               stderr=subprocess.STDOUT)
        (stdout_data,stderr_data) = job.communicate()
        stdout_buffer = os.linesep.join([stdout_buffer, '', '', 
                                         ' '.join(step), '', stdout_data,
                                         '', ''])
    if verbose:
        stdout_buffer = os.linesep.join([stdout_buffer, 
                           ' '.join(['Removing temporary file',
                            os.path.basename(parsed_filename.temp_file)])])

    os.remove(parsed_filename.temp_file)

    return stdout_buffer

if __name__=="__main__":
	scripter.perform(convert_sam_to_bam, FilenameParser)

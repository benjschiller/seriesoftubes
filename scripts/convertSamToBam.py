#!/usr/bin/env python
'''
Converts SAM to BAM files
and also builds an index for the BAM files

Default options: --target=alignments.BAM
'''
import os
import subprocess
import scripter
from scripter import print_debug
import pysam
scripter.SCRIPT_DOC = __doc__
scripter.SCRIPT_VERSION = "2.2"
scripter.SOURCE_DIR = 'alignments.SAM'
scripter.TARGET_DIR = 'alignments.BAM'

class FilenameParser(scripter.FilenameParser):
    def __init__(self, filename, *args, **kwargs):
        super(FilenameParser, self).__init__(filename, *args,**kwargs)

        f = lambda x: os.path.join(self.output_dir, self.with_extension(x))
        self.output_file = f('bam')
        self.temp_file = f('tmp')

def convert_sam_to_bam(parsed_filename, verbose=False, debug=False, **kwargs):

    stdout_buffer = ''

    if debug: print_debug('Converting', parsed_filename.input_file, 
                          'to temporary BAM file', parsed_filename.temp_file)
    # this methods raises an error, we need to catch it
    try:
        pysam.view('-b', '-S', '-o' + parsed_filename.temp_file, 
                   parsed_filename.input_file)
    except pysam.SamtoolsError, converting_error:
        converting = str(converting_error)
        if not converting.startswith("'[samopen]"):
            raise converting_error
    if debug: print_debug(converting)

    if debug: print_debug('Sorting BAM file...')
    sorting = pysam.sort(parsed_filename.temp_file,
                         os.path.splitext(parsed_filename.output_file)[0])
    sorting = os.linesep.join(sorting)
    if debug: print_debug(sorting)

    if debug: print_debug('Indexing sorted BAM file...')
    indexing = pysam.index(parsed_filename.output_file)
    indexing = os.linesep.join(indexing)
    if debug: print_debug(indexing)

    if debug: print_debug('Removing temporary file', parsed_filename.temp_file)
    elif verbose: stdout_buffer = os.linesep.join([stdout_buffer, 
                                                  'Converting...', converting,
                                                  'Sorting...', sorting,
                                                  'Indexing', indexing])

    os.remove(parsed_filename.temp_file)

    return stdout_buffer

if __name__=="__main__":
	scripter.perform(convert_sam_to_bam, FilenameParser)

#!/usr/bin/env python
'''
Converts SAM to BAM files
and also builds an index for the BAM files

--no-unmapped   do not produce BAM files with only unmapped reads
--no-mapped     do not produce BAM files with only mapped reads

Default options: --source=alignments.SAM --target=alignments.BAM
'''
import os
import subprocess
import scripter
from scripter import print_debug
import pysam
scripter.SCRIPT_DOC = __doc__
scripter.SCRIPT_VERSION = "2.3"
scripter.SOURCE_DIR = 'alignments.SAM'
scripter.TARGET_DIR = 'alignments.BAM'
scripter.SCRIPT_LONG_OPTS = ["extra-samtools-options="]

def check_script_options(options):
     specific_options = {}
     if options.has_key('extra-samtools-options'):
         specific_options['extra_samtools_options'
                          ] = options['extra-samtools-options'].strip().split()
     else:
         specific_options['extra_samtools_options'] = []

     return specific_options

class FilenameParser(scripter.FilenameParser):
    def __init__(self, filename, *args, **kwargs):
        super(FilenameParser, self).__init__(filename, *args,**kwargs)

        f = lambda x: os.path.join(self.output_dir, self.with_extension(x))
        self.output_file = f('bam')
        self.temp_file = f('tmp')

def pysam_wrapper(pysam_args, pysam_func=None, expected_error='',
                  except_merging=True):
    '''wraps a pysam function, captures errors that match an expected error
       string (if not provided, captures all SamtoolsError's)
       and then returns the error string

       or if there was no error, it returns the output as a string

       this function guarantees you a string, altho it may throw an exception
    '''
    try:
        execution = pysam_func(*pysam_args)
    except pysam.SamtoolsError, execution_error:
        execution = str(execution_error)
        if execution.lstrip("'").startswith(expected_error):
            pass
        elif except_merging and 'merging from' in execution:
            pass
        else:
            raise pysam.SamtoolsError(execution)

    if type(execution) is list:
        if len(execution) is 0: return ''
        else: return os.linesep.join(execution)
    elif type(execution) is str:
        return execution
    else:
        return ''

def convert_sam_to_bam(parsed_filename, verbose=False, debug=False,
                       extra_samtools_options=[], **kwargs):

    stdout_buffer = ''

    if debug: print_debug('Converting', parsed_filename.input_file, 
                          'to temporary BAM file', parsed_filename.temp_file)
    # this methods raises an error, we need to catch it
    converting = pysam_wrapper(extra_samtools_options + 
                               ['-b', '-S', '-o' + parsed_filename.temp_file, 
                               parsed_filename.input_file],
                               pysam_func=pysam.view,
                               expected_error='[samopen]')
    if debug: print_debug(converting)

    if debug: print_debug('Sorting BAM file...')
    sorting = pysam_wrapper((parsed_filename.temp_file,
                            os.path.splitext(parsed_filename.output_file)[0]),
                            pysam_func=pysam.sort,
                            expected_error='[')
    if debug: print_debug(sorting)

    if debug: print_debug('Indexing sorted BAM file...')
    indexing = pysam_wrapper((parsed_filename.output_file,),
                             pysam_func=pysam.index)
    if debug: print_debug(indexing)

    if debug: print_debug('Removing temporary file', parsed_filename.temp_file)
    elif verbose: stdout_buffer = os.linesep.join([stdout_buffer, 
                                                  'Converting...', converting,
                                                  'Sorting...', sorting,
                                                  'Indexing', indexing])

    os.remove(parsed_filename.temp_file)

    return stdout_buffer

if __name__=="__main__":
    scripter.check_script_options = check_script_options
    scripter.perform(convert_sam_to_bam, FilenameParser)

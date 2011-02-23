#!/usr/bin/env python
#Copyright (c) 2010 Benjamin Schiller
#
#This file is part of seriesoftubes.
#
#seriesoftubes is free software; you can redistribute it and/or modify it under
#the terms of the Artistic Free License 2.0. All rights not covered by the 
#Artistic Free License 2.0 are reserved. The Artistic Free License 2.0 will take
#precedence in the event of a conflict with this copyright notice.
#
# See <http://www.opensource.org/licenses/artistic-license-2.0.php> for more information.
'''
Converts SAM to BAM files
and also builds an index for the BAM files

--no-unfiltered   do not produce BAM files with only unmapped reads
--no-filtered     do not produce BAM files with only mapped reads

Default options: --source=alignments.SAM --target=alignments.BAM
'''
import os
import subprocess
import scripter
from scripter import print_debug
import pysam
VERSION = "2.4"

def main():
    long_opts = ["no-filtered", "no-unfiltered"]
    e = scripter.Environment(long_opts=long_opts, version=VERSION, doc=__doc__)
    filtered, unfiltered = check_script_options(e.get_options())
    e.set_source_dir('alignments.SAM')
    e.set_target_dir('alignments.BAM')
    e.set_filename_parser(FilenameParser)
    if unfiltered:
        e.do_action(convert_sam_to_bam)
    if filtered:
        e.update_script_kwargs({'extra-samtools-options': ['-F', '0x4']})
        e.set_target_dir('alignments_filtered.BAM')
        e.do_action(convert_sam_to_bam)

def check_script_options(options):
#     specific_options = {}
#     if options.has_key('extra-samtools-options'):
#         specific_options['extra_samtools_options'
#                          ] = options['extra-samtools-options'].strip().split()
#     else:
#         specific_options['extra_samtools_options'] = []
#
#     return specific_options
    filtered = True
    unfiltered = True
    if options.has_key("no-filtered"): filtered = False
    if options.has_key("no-unfiltered"): unfiltered = False
    return filtered, unfiltered

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
        else: return '\n'.join(execution)
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
    elif verbose: stdout_buffer = '\n'.join([stdout_buffer, 
                                             'Converting...', converting,
                                             'Sorting...', sorting,
                                             'Indexing', indexing])

    os.remove(parsed_filename.temp_file)

    return stdout_buffer

if __name__=="__main__": main()

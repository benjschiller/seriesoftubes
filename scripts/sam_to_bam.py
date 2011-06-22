#!/usr/bin/env python
'''
This script is deprecated
Converts SAM to BAM files
and also builds an index for the BAM files
Output is two folders
    all.BAM/        sorted, indexed BAM files with mapped + unmapped reads
    aligned.BAM/    sorted, indexed BAM files with mapped reads only

--no-filtering    do not produce aligned.BAM/ folder
--no-indexing     do not produce .bam.bai index files
--no-sorting      do not sort BAM files, just convert

Default options: --source=alignments.SAM
'''
raise DeprecationWarning
import os
import re
import pysam
import scripter
from scripter import print_debug
VERSION = "2.4"

def main():
    long_opts = ["no-filtering", "no-indexing", "no-sorting"]
    e = scripter.Environment(long_opts=long_opts, version=VERSION, doc=__doc__)
    new_kwargs = check_script_options(e.get_options())
    e.update_script_kwargs(new_kwargs)
    e.set_source_dir('alignments.SAM')
    e.set_target_dir('alignments.BAM')
    e.set_filename_parser(FilenameParser)
    e.do_action(convert_sam_to_bam)

def check_script_options(options):
    filtering = not options.has_key("no-filtering")
    sorting = not options.has_key("no-sorting")
    indexing = not options.has_key("no-indexing")
    return {'filtering': filtering, 'indexing': indexing, 'sorting': sorting}

class FilenameParser(scripter.FilenameParser):
    def __init__(self, filename, *args, **kwargs):
        super(FilenameParser, self).__init__(filename, *args,**kwargs)

        f = lambda x: os.path.join(self.output_dir, self.with_extension(x))
        self.output_file = f('bam')
        self.temp_file = f('tmp')

def convert_sam_to_bam(parsed_filename, verbose=False, debug=False,
                       sorting=True, indexing=True, filtering=True, **kwargs):
    temp_file = parsed_filename.temp_file
    input_file = parsed_filename.input_file
    output_file = parsed_filename.output_file

    if sorting:
        if debug: print_debug('Converting', input_file, 
                              'to temporary BAM file', temp_file)
        pysam.view('-b', '-S', '-o' + parsed_filename.temp_file, 
                   parsed_filename.input_file)
        if debug: print_debug('Sorting', input_file, 
                              'to temporary BAM file', temp_file)
        pysam.sort(temp_file,
                   os.path.splitext(output_file)[0])
        if debug: print_debug('Removing temporary file', temp_file)
        os.remove(parsed_filename.temp_file)
    else:
        if debug: print_debug('Converting', input_file, 
                              'to temporary BAM file', output_file)
        pysam.view('-b', '-S', '-o' + parsed_filename.output_file, 
                   parsed_filename.input_file)
    if indexing:
        if debug: print_debug('Indexing', output_file) 
        pysam.index(output_file)
    
    if filtering:
        odir, ofile = os.path.split(output_file)
        fildir = re.sub('alignments.BAM', 'alignments_filtered.BAM', odir, 1)
        parsed_filename.check_output_dir(fildir)
        filtered_file = os.path.join(fildir, ofile)
        output_files = ' and '.join([output_file, filtered_file])
        
        if debug: print_debug('Copying aligned reads from {!s} to {!s}'.\
                              format(output_file, filtered_file))
        pysam.view('-F 0x4', '-b', '-o' + filtered_file, output_file)
    else:
        output_files = output_file
    
    if sorting and indexing: attrs='(sorted and indexed)'
    elif indexing: attrs='(indexed)'
    elif sorting: attrs='(sorted)'
    else: attrs=''
    return 'Converted {!s} to {!s} {!s}'.format(input_file,
                                                  output_files,
                                                  attrs)

if __name__=="__main__": main()

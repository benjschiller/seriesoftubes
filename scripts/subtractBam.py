#!/usr/bin/env python
'''
--sam-in            Input files are SAM (not BAM)
--sam-out           Output SAM files (not BAM)
--remove-all        Remove all reads in mapped BAM/SAM file, not just mapped ones

subtracts bam files based on sequence names
(for example, use to remove reads that also map rRNA)

outputs to alignments_filtered.BAM

by default, looks in alignments.BAM for .bam files
and expects corresponding files in mapped/alignments.BAM
ignores the species folder name i.e. mapped/alignments.BAM(/foo)

e.g
alignments.BAM/crypto_h99_grubii/random/crypto_small_RNA.dcr.2010.s_5.no_linker.bam
mapped/alignments.BAM/crypto_h99_grubii_rRNA/random/crypto_small_RNA.dcr.2010.s_5.no_linker.bam
'''
import glob
import os
import subprocess
import collections
import scripter
import pysam
scripter.SCRIPT_DOC = __doc__
scripter.SCRIPT_VERSION = "2.2"
scripter.SOURCE_DIR = 'alignments.BAM'
scripter.TARGET_DIR = 'alignments_filtered.BAM'
BOOLEAN_OPTS = ["sam-in", "sam-out", "remove-all"]
scripter.SCRIPT_LONG_OPTS = BOOLEAN_OPTS

def check_script_options(options):
    specific_options = {}

    for option in BOOLEAN_OPTS:
        pyoption = "_".join(option.split("-"))
        specific_options[pyoption] = options.has_key(option)

    return specific_options

def get_removable_reads(parsed_filename, sam_in=False, remove_all=False,
                        **kwargs):
    '''
returns a list of read names that are mapped from a BAM file

looks at parsed_filename.
'''
    removeable_reads = []

    if sam_in: open_opts = 'rb'
    else: open_opts = 'r'
    bam_file = pysam.Samfile(parsed_filename.mapped_filename, open_opts)

    if remove_all:
        return [read.qname for read in bam_file]

    for read in bam_file:
        if is_mapped(read):
            removeable_reads.append(read.qname)

    return removable_reads

def is_mapped(read):
    '''
true if a pysam.AlignedRead is mapped
also checks if read is paired end; if so, returns if the pair is mapped
'''
    if read.is_paired:
        if read.is_proper_pair: return True
        else: return False
    else:
        if read.is_unmapped: return False
        else: return True

def remove_reads(parsed_filename, reads_to_remove, sam_in=False, sam_out=False,
                 debug=False, **kwargs):
    '''
note: you must be looking at a sorted file, or this won't work
'''
    if sam_in: open_opts = 'rb'
    else: open_opts = 'r'
    bam_file = pysam.Samfile(parsed_filename.mapped_file, open_opts)

    if sam_out: write_opts = 'wb'
    else: write_opts = 'w'
    out_bam_file = pysam.Samfile(parsed_filename.output_file, write_opts,
                                 template = bam_file)
    for read in bam_file:
        if read.qname in reads_to_remove: continue
        else: out_bam_file.write(read)

    if debug: return stdout_buffer
    else: return ''

class SubtractBamFilenameParser(scripter.FilenameParser):
    def __init__(self, filename, verbose=False, sam_out=False, *args, **kwargs):
        super(BowtieFilenameParser, self).__init__(filename, *args, **kwargs)

        # check for the mapped_file
        potential_filenames = glob.glob(os.path.join(
                                            ['mapped',
                                            self.input_dir[0],
                                            '*',
                                            self.input_dir[2:],
                                            os.path.basename(self.input_file)]))

        if len(potential_filenames) is 1:
            self.mapped_file = potential_filenames[0]
        elif len(potential_filenames) is 0:
            raise scripter.Usage('Could not find mapped file')
        else:
            raise scripter.Usage('Ambiguous mapped file', *potential_filenames)

        if sam_out:
            self.output_file = self.output_dir + self.with_extension('sam')
        else:
            self.output_file = self.output_dir + self.with_extension('bam')

if __name__=="__main__":
    scripter.check_script_options = check_script_options
	scripter.perform(calculate_duplicates, SubtractBamFilenameParser)

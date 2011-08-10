#!/usr/bin/env python
'''
subtracts BAM (or SAM) files based on sequence names
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
import scripter
from scripter import InvalidFileException
import pysam
from pkg_resources import get_distribution
__version__ = get_distribution('seriesoftubes').version
VERSION = __version__

def main():
    e = scripter.Environment(version=VERSION, doc=__doc__)
    parser = e.argument_parser
    parser.add_argument('--sam-out', action='store_true',
                        help='Output SAM files (not BAM)')
    parser.add_argument('--remove-all', action='store_true',
                        help='Remove all reads in mapped BAM/SAM file, not just mapped ones')
    parser.set_defaults(**{'target': 'filtered'})
    e.set_filename_parser(SubtractBamFilenameParser) 
    e.do_action(remove_reads)
    
def is_mapped(read):
    '''
    true if a pysam.AlignedRead is properly mapped, otherwise false
    (checks if read is paired end; if so, returns if the pair is mapped)
    '''
    if read.is_paired: return read.is_proper_pair
    else: return not read.is_unmapped

def remove_reads(parsed_filename, remove_all=False, sam_out=False,
                 debug=False, **kwargs):
    '''
note: you must be looking at a sorted file, or this won't work
'''

    if sam_out: write_opts = 'w'
    else: write_opts = 'wb'

    with pysam.Samfile(parsed_filename.mapped_file) as mapped:
        if remove_all:
            reads_to_remove = set([read.qname for read in mapped])
        else:
            reads_to_remove = set([read.qname for read in mapped
                                   if is_mapped(read)])

    if debug:
        scripter.debug('Found {!s} reads in {!s}'.format(len(reads_to_remove),
                                                   parsed_filename.mapped_file))
    
    with pysam.Samfile(parsed_filename.input_file) as bam_file:
        with pysam.Samfile(parsed_filename.output_file, write_opts,
                           template = bam_file) as out_bam_file:
            for read in bam_file:
                if read.qname not in reads_to_remove:
                    out_bam_file.write(read)

    return

class SubtractBamFilenameParser(scripter.FilenameParser):
    def __init__(self, filename, sam_out=False, *args, **kwargs):
        super(SubtractBamFilenameParser, self).__init__(filename,
                                                        sam_out=sam_out,
                                                        *args,
                                                        **kwargs)
        fext = os.path.splitext(filename)[1].rstrip(os.extsep)
        if not (fext == 'sam' or fext =='bam'): raise InvalidFileException
        if not self.is_dummy_file:
            # check for the mapped_file
            input_dir_parts = self.input_dir.split(os.sep)
            glob_path = ['mapped', input_dir_parts[0], '*'] + \
                        input_dir_parts[2:] + \
                        [os.path.basename(self.input_file)]
            potential_filenames = glob.glob(os.sep.join(glob_path))

            if len(potential_filenames) is 1:
                self.mapped_file = potential_filenames[0]
            elif len(potential_filenames) is 0:
                raise scripter.Usage('Could not find mapped file')
            else:
                raise scripter.Usage('Ambiguous mapped file', *potential_filenames)
            scripter.debug('Mapped file will be', self.mapped_file)

            if sam_out:
                self.output_file = os.sep.join([self.output_dir,
                                               self.with_extension('sam')])
            else:
                self.output_file = os.sep.join([self.output_dir,
                                               self.with_extension('bam')])
            scripter.debug('Output file will be', self.output_file)

if __name__=="__main__": main()

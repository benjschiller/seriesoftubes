#!/usr/bin/env python
"""
find_sites searches a pairs of BED (or MACS xls) and FASTA files
for matches to a TFBS motif
"""
import os
import re
import scripter
import Bio.Motif
from bioplus.sitefinder import find_sites, USE_MOODS
from scripter import Usage, InvalidFileException, get_logger
if not USE_MOODS:
    raise ImportError("""WARNING: MOODS is not installed. You may obtain it from
http://www.cs.helsinki.fi/group/pssmfind/""")
from pkg_resources import get_distribution
__version__ = get_distribution('seriesoftubes').version
VERSION = __version__

def main():
    e = scripter.Environment(doc=__doc__, version=VERSION)
    parser = e.argument_parser
    parser.add_argument('--motif', required=True,
                        help='Path to file containing motif')
    parser.add_argument('--motif-number',
                        help='Motif number within file (e.g. 1, 2, 3) [Default is to assume only one motif in file]')
    parser.add_argument('--motif-type', default='MEME',
                        help='motif type (see Bio.Motif for more info)')
    e.set_filename_parser(FilenameParser)
    e.do_action(action)
   
def get_motif(foo, motif_number, motif_type):
    """
    get motif_number motif from file foo
    """
    if not os.path.exists(foo):
        raise Usage(foo, "not found")
    motifs = Bio.Motif.parse(open(foo), motif_type)
    for i in range(motif_number):
        try:
            motif = motifs.next()
        except StopIteration:
            raise Usage("%s only contains %d motifs") % (foo ,i)
    return motif

def action(fp_obj, motif_file=None, motif_type=None, motif_number=1,
           debug=False, **kwargs):
    logger = get_logger()
    logger.debug("trying to find sites for",
                          fp_obj.input_file)
    motif = get_motif(motif_file, motif_type, motif_type)
    stdout_buffer = find_sites(fp_obj.input_file,
                               fp_obj.fasta_file,
                               motif, bed=fp_obj.is_bed, xls=fp_obj.is_xls,
                               output_dir=fp_obj.output_dir,
                               src_fnc = "find_sites.py",
                               debug=debug)
    return stdout_buffer

class FilenameParser(scripter.FilenameParser):
    def __init__(self, filename, include_width_in_name=False,
                 target_dir=None, motif_file='unknown_motif',
                 *args, **kwargs):
        fext = os.path.splitext(filename)[1].lstrip(os.extsep)
        if not fext == 'bed':
            self.is_bed = True
            self.is_xls = False
        elif fext == 'xls':
            self.is_bed = False
            self.is_xls = True
        else: raise InvalidFileException
        motif_name = re.sub('\W', '_', os.path.abspath(motif_file))
        target_dir = target_dir + os.sep + motif_name
        super(FilenameParser, self).__init__(filename,
                                             target_dir = target_dir,
                                             *args, **kwargs)
        self.fasta_file = None
        for file_extension in ['fa', 'fasta', 'FA', 'FASTA']:
            fasta_file = os.path.join(self.input_dir,
                            os.extsep.join([self.protoname, file_extension]))
            scripter.debug("Trying", fasta_file)
            if os.path.exists(fasta_file):
                self.fasta_file = fasta_file
                scripter.debug("Using", fasta_file)
                break
        if self.fasta_file is None:
            raise Usage("Could not find the FASTA file for ", self.input_file)
            # or in the future, fetch the sequence?

if __name__=="__main__": main()
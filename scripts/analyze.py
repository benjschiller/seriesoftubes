#!/usr/bin/env python
"""
find_sites searches a pairs of BED (or MACS xls) and FASTA files
for matches to a TFBS motif
"""
from os.path import exists
import Bio.Motif
from bioplus.sitefinder import find_sites, USE_MOODS
from seriesoftubes.fnparsers import PeaksFilenameParser
from scripter import Environment, get_logger, Usage
if not USE_MOODS:
    raise ImportError("""WARNING: MOODS is not installed. You may obtain it from
http://www.cs.helsinki.fi/group/pssmfind/""")
from pkg_resources import get_distribution
__version__ = get_distribution('seriesoftubes').version
VERSION = __version__

def main():
    e = Environment(doc=__doc__, version=VERSION)
    e.set_filename_parser(PeaksFilenameParser)
    parser = e.argument_parser
    parser.add_argument('--motif', required=True, dest='motif_file',
                        help='Path to file containing motif')
    parser.add_argument('--motif-number', type=int, default=1,
                        help='Motif number within file (e.g. 1, 2, 3) [Default is 1st]')
    parser.add_argument('--motif-type', default='MEME',
                        help='motif type (see Bio.Motif for more info)')
    parser.add_argument('--genome',
                        help='Reference genome (path to 2bit file)')
    parser.set_defaults(**{'target': 'analysis'})
    e.do_action(action)
   
def get_motif(foo, motif_number, motif_type):
    """
    get motif_number motif from file foo
    """
    if not exists(foo):
        raise Usage(foo, "not found")
    motifs = Bio.Motif.parse(open(foo), motif_type)
    for i in xrange(motif_number):
        try:
            motif = motifs.next()
        except StopIteration:
            raise Usage("%s only contains %d motifs") % (foo ,i)
    return motif

def action(fp_obj, motif_file=None, motif_type=None, motif_number=1,
           **kwargs):
    logger = get_logger()
    logger.debug("trying to find sites for %s", fp_obj.input_file)
    motif = get_motif(motif_file, motif_number, motif_type)
    stdout_buffer = find_sites(fp_obj.input_file,
                               fp_obj.fasta_file,
                               motif, bed=fp_obj.is_bed, xls=fp_obj.is_xls,
                               output_dir=fp_obj.output_dir,
                               src_fnc = __file__)
    return stdout_buffer
                
if __name__=="__main__": main()
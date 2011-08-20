#!/usr/bin/env python
"""
fetch the FASTA sequences for the peaks
"""
import os
import operator
import platform
import argparse
import twobitreader
import scripter
from scripter import Usage, InvalidFileException, get_logger
from pkg_resources import get_distribution
__version__ = get_distribution('seriesoftubes').version
VERSION = __version__

def positive_int(i):
    """
    takes a string specifying an integer
    returns a positive integer or raises an ArgumentTypeError
    """
    try:
        I = int(i)
    except ValueError, TypeError:
        msg = '%d is not an integer' 
        raise argparse.ArgumentTypeError(msg)
    if not I > 0:
        msg = '%d is not a positive integer' 
        raise argparse.ArgumentTypeError(msg)
    if I%2 == 0: return I + 1
    else: return I
    
def main():
    e = scripter.Environment(doc=__doc__, version=VERSION)
    parser = e.argument_parser
    parser.add_argument('--include-width-in-name', action='store_true',
                        help='Include width in directory name for output (e.g. peaks_150.FASTA)')
    parser.add_argument('--bed', action='store_true',
                        help='Make a BED file too saying where sequences come from')
    parser.add_argument('--sort', action='store_true',
                        help="""Sort peaks if possible by
    number of reads under peak (MACS-xls)
    number of reads under summit (MACS-subpeaks)
    p-value (MACS-peaks)""")
    sgroup = parser.add_mutually_exclusive_group(required=True)
#                                             'source',
#                                        'Source input files (use only one)')
    sgroup.add_argument('--from-MACS-subpeaks', action='store_true',
                        help='Fetch the peaks from MACS _subpeaks.bed '
                        '(uses the summit position in column 4)')
    sgroup.add_argument('--from-MACS-xls', action='store_true',
                        help='Fetch peaks from MACS _peaks.xls')
    sgroup.add_argument('--from-MACS-bed', action='store_true',
                        help='Fetch peaks from MACS _peaks.bed')
    sgroup.add_argument('--peaks', action='store_true',
                        help='Assume the BED file contains the positions of peaks')
    sgroup.add_argument('--summits', action='store_true',
                        help='Assume the BED file contains the positions of peak summits')
    parser.add_argument('--path-to-gbdb',
help='Location of "gdbdb" or 2bit files. If gbdb is not in /gbdb or C:\gbdb, specify the path here')
    ggroup = parser.add_mutually_exclusive_group(required=True)
    ggroup.add_argument('--ref',
               help='Use 2bit file foo as reference genome (Looks also for {path-to-gbdb}/foo/foo.2bit))')
    ggroup.add_argument('--hg18', const='hg18', action='store_const',
               help='Shortcut for --ref hg18') 
    ggroup.add_argument('--hg19', const='hg19', action='store_const',
               help='Shortcut for --ref hg19') 
    ggroup.add_argument('--mm9', const='mm9', action='store_const',
               help='Shortcut for --ref mm9') 
    parser.add_argument('--width', type=positive_int,
                        help="""Get sequence x distance from the peak center/summit
    if width is even, we will use width + 1. Default is to use interval.""")
    parser.add_argument('--npeaks', type=positive_int,
                        help='Include only top # peaks (if not sorted, use first peaks)')
    parser.set_defaults(**{'target': 'sequences'})
    e.set_filename_parser(FilenameParser)
    e.do_action(get_sequences)

def find_2bit_file(ref, path_to_gbdb=None):
    if ref is None: raise Usage("No reference genome specified")
    if os.path.exists(ref): return ref
    fname = "{!s}{!s}2bit".format(ref, os.extsep)
    if path_to_gbdb is not None:
        # check in its own dir
        specified_path = os.path.join(path_to_gbdb, ref, fname)
        if os.path.exists(specified_path): return specified_path
        # check in the main dir
        specified_path = os.path.join(path_to_gbdb, fname)
        if os.path.exists(specified_path): return specified_path
    # try absolute paths
    new_path = os.path.join("gbdb", ref, fname)
    if platform.system() == 'Windows':
        home_drive = os.environ['HOMEDRIVE']
        windows_path = "{!s}{!s}{!s}".format(home_drive, os.sep, new_path)
        if os.path.exists(windows_path): return windows_path
    else: # assume *nix
        unix_path = os.sep + new_path
        if os.path.exists(unix_path): return unix_path
    # give up
    raise Usage("Could not find 2bit file ", fname)


def detect_reference(parsed_filename, path_to_gbdb=None):
    # finish
    return find_2bit_file("hg19", path_to_gbdb=path_to_gbdb)

def write_to_fasta(file_handle, sequence, name=''):
    # write first line
    file_handle.write('>{!s}\n'.format(name))
    # write sequence
    num_lines = (len(sequence)+ 59 ) / 60
    for j in xrange(num_lines):
        start = 0 + j*60
        end = 60 + j*60
        file_handle.write('{!s}\n'.format(sequence[start:end]))
    return

def get_sequences(parsed_filename, from_MACS_subpeaks=True, from_MACS_xls=False,
           from_MACS_bed=False, peaks=False, summits=False,
           ref=None, width=None, path_to_gbdb=None, sort=False,
           npeaks=None, bed=False, logger=None,
           debug = False,
           **kwargs):
    ref = find_2bit_file(ref, path_to_gbdb)
    logger = get_logger()
    if from_MACS_subpeaks:
        sort_item = lambda x: int(x[3]) 
        if width is not None:
            center_coord = lambda x: int(x[4]) # the 5th col
            start_coord = lambda x: center_coord(x) - width/2
            end_coord = lambda x: center_coord(x) + width/2 + 1
        else:
            start_coord = lambda x: int(x[1])
            end_coord = lambda x: int(x[2])
    if from_MACS_xls:
        sort_item = lambda x: int(x[5])
        if width is not None:
            center_coord = lambda x: int(x[1]) + int(x[4]) # start + 5th col
            start_coord = lambda x: center_coord(x) - width/2
            end_coord = lambda x: center_coord(x) + width/2 
        else:
            start_coord = lambda x: int(x[1])
            end_coord = lambda x: int(x[2]) - 1
    if from_MACS_bed or peaks:
        sort_item = lambda x: float(x[4])
        if width is not None:
            center_coord = lambda x: (int(x[1]) + int(x[2])) / 2
            start_coord = lambda x: center_coord(x) - width/2
            end_coord = lambda x: center_coord(x) + width/2 + 1
        else:
            start_coord = lambda x: int(x[1])
            end_coord = lambda x: int(x[2])
    if summits:
        if width is None: raise Usage('--width required for --summits')
        sort_item = lambda x: float(x[4]) 
        center_coord = lambda x: int(x[1])
        start_coord = lambda x: center_coord(x) - width/2
        end_coord = lambda x: center_coord(x) + width/2 + 1

    # output file
    output_file = os.path.join(parsed_filename.output_dir,
                               os.extsep.join([parsed_filename.protoname,
                                               'fa']))
    if bed: bed_filename = os.path.join(parsed_filename.output_dir,
                                os.extsep.join([parsed_filename.protoname,
                                                'bed']))

    if ref is None:
        ref = detect_reference(parsed_filename, path_to_gbdb)
    logger.debug("Using", ref, "to extract sequences", "from", 
                    parsed_filename.input_file, "to", output_file)
        
    ref_genome = twobitreader.TwoBitFile(ref)
    output_handle = open(output_file, 'w')
    input_file = open(parsed_filename.input_file, 'rU')
    input_lines = (line.rstrip() for line in input_file)

    logger.debug("Determining chromsome lengths for", ref)
    chrom_lengths = ref_genome.sequence_sizes()
    valid_chroms = chrom_lengths.keys()

    if bed: bed_file = open(bed_filename, 'w')
    logger.debug("Now converting {!s}".format(bed_filename)) 
    # if we're going to sort, we'll have to hold the list in memory
    num_peaks = 0
    bed_template = "{!s}\t{!s}\t{!s}\tpeak{!s}\t{!s}\t+\n"
    seq_list = [] 
    for line in input_lines:
        words = line.split()
        if len(words) == 0: continue
        # determine chrom/contig
        chrom = words[0].lower()
        # skip if not a valid contig/chrom
        if not chrom in valid_chroms: continue
        # determine start and end
        start = start_coord(words)
        end = end_coord(words)
        # skip if we've already fallen off the chromosome
        if int(start) > chrom_lengths[chrom]: continue
        # truncate chrom start for MACS
        num_peaks += 1
        start = max(int(start), 0)
        end = end # don't need to truncate end, twobitreader does that
        name = "peak_{!s}_{!s}_{!s}_{!s}".format(num_peaks, chrom, start, end)
        sitem = (chrom, start, end, name, sort_item(words))
        if sort:
            seq_list.append(sitem)
        else:
            sequence = ref_genome[chrom][start:end]
            write_to_fasta(output_handle, sequence, name=name)
            if bed:
                bed_file.write(bed_template.format(*sitem))

    if sort:
        logger.debug("Sorting {!s} peaks".format(num_peaks))
        sorted_list = sorted(seq_list, key=operator.itemgetter(4), reverse=True)
        if npeaks is None or num_peaks < npeaks:
            npeaks = num_peaks
            logger.debug("Writing results for \
all {!s} peaks".format(npeaks))
        else:
            logger.debug("Writing results for \
top {!s} peaks".format(npeaks))
                
        for x in xrange(npeaks):
            sitem = sorted_list[x]
            chrom, start, end = sitem[0:3]
            sequence = ref_genome[chrom][start:end]
            write_to_fasta(output_handle, sequence, name=sitem[3])
            if bed: bed_file.write(bed_template.format(*sitem))
            
    output_handle.close()
    if bed: bed_file.close()
    return


class FilenameParser(scripter.FilenameParser):
    def __init__(self, filename, include_width_in_name=False, target=None,
                 width = None,
                 from_MACS_subpeaks=False, from_MACS_xls=False,
                 from_MACS_bed=False, peaks=False, summits=False,
                 *args, **kwargs):
        fext = os.path.splitext(filename)[1].lstrip(os.extsep)
        if not (fext == 'bed' or fext == 'xls'): raise InvalidFileException
        if target is None:
            scripter.error('Cannot use width in name, width not specified')
            if include_width_in_name: target = "peaks_%dbp" % width
            else: target= 'peaks'
        super(FilenameParser, self).__init__(filename,
                                             target=target,
                                             *args, **kwargs)
        if from_MACS_subpeaks:
            if not self.input_file.endswith("_subpeaks.bed"):
                raise InvalidFileException
        elif from_MACS_xls:
            if not self.input_file.endswith("_peaks.xls"):
                raise InvalidFileException
        elif from_MACS_bed:
            if not self.protoname.endswidth("_peaks.bed"):
                raise InvalidFileException
        elif peaks or summits:
            if not self.file_extension=="bed":
                raise InvalidFileException

if __name__=="__main__": main()

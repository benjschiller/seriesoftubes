#!/usr/bin/env python
"""
fetch the FASTA sequences for the peaks

requires 2bit file and the twoBitToFa, twoBitInfo tools
from UCSC genome browser project

Source input files (use only one):
--from-MACS-subpeaks        Fetch the peaks from MACS _subpeaks.bed
                              (uses the summit position in column 4)
--from-MACS-xls             Fetch peaks from MACS _peaks.xls
--from-MACS-bed             Fetch peaks from MACS _peaks.bed
--peaks                     Assume the BED file contains the positions
                                of peaks
--summits                   Assume the BED file contains the positions
                                of peak summits

Genome options (use only one):
--ref=foo                   Use 2bit file foo as reference genome
                                (Looks also for /gbdb/foo/foo.2bit)
--hg18, --hg19, --mm9       aliases for --ref hg18, etc.
--path-to-gbdb=foo          If gbdb is not in /gbdb or C:\gbdb,
                                specify the path here
                                 
Options for --peaks:
--from-center               Grab some interval of sequence from center of peak
                                (for use with --width)

Common options:
--width=x                   Get sequence x distance from the peak center/summit
                                (implies --from-center when applicable)
                                (default: 150)
                                (note: if width is even, we will use width + 1)
--include-width-in-name     Include width in directory name for output
                                (e.g. peaks_150.FASTA)
--sort                      Sort peaks if possible by
                                number of reads under peak (MACS-xls)
                                number of reads under summit (MACS-subpeaks)
                                p-value (MACS-peaks)
--npeaks=#                  Include only the top # of peaks. If not --sort,
                                we'll just take the first #.
--bed                       Make a BED file too saying where sequences come from
"""
import os
import operator
import platform
import twobitreader
import scripter
from scripter import Usage, print_debug, InvalidFileException
VERSION = "2.4"
SOURCE_DIR = 'fromBAM.MACS'
TARGET_DIR = 'peaks.FASTA'

def main():
#    if platform.system() == 'Windows': windows_exit()
    
    source_opts = ["from-MACS-subpeaks", "from-MACS-xls",
                   "from-MACS-bed", "peaks", "summits"]
    boolean_opts = ["include-width-in-name", "from-center", "sort",
                    "bed"] + source_opts
    long_opts = ["ref=", "hg19", "hg18", "mm9", "path-to-gbdb=",
                 "width=","npeaks="] + boolean_opts
    e = scripter.Environment(long_opts=long_opts, doc=__doc__, version=VERSION)
    e.parse_boolean_opts(boolean_opts)
    all_opts = check_script_options(e.get_options(), sources=source_opts,
                                    debug=e.is_debug())
    e.update_script_kwargs(all_opts)
    e.set_source_dir(SOURCE_DIR)
    e.set_target_dir(TARGET_DIR)
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

def check_script_options(options, sources=[], debug=False):
    sopts = {}
   
    # Find all reference genomes and their paths
    options_count = sum([options.has_key("hg18"), options.has_key("hg19"),
                        options.has_key("mm9"), options.has_key("ref")])
    
    if options.has_key('path-to-gbdb'):
        path_to_gbdb = options['path-to-gbdb']
        sopts["path_to_gbdb"] = path_to_gbdb
    else:
        path_to_gbdb = None
    ref = None
    if options_count == 0:
        ref = None
    elif options_count > 1:
        raise Usage("More than one reference genome specified")
    else:
        if options.has_key("hg18"): ref = "hg18"
        elif options.has_key("hg19"): ref = "hg19"
        elif options.has_key("mm9"): ref = "mm9"
        elif options.has_key("ref"): ref = options["ref"]
        sopts["ref"] = find_2bit_file(ref, path_to_gbdb)

    source_count = sum([options.has_key(x) for x in sources])
    if not source_count == 1:
        raise Usage("Must specify exactly one of",
                    " ".join(["".join(["--", x]) for x in sources]))    

    if options.has_key("npeaks"): sopts["npeaks"] = int(options["npeaks"])

    if options.has_key("width"):
        try: width = int(options["width"])
        except ValueError: raise Usage("width must be an integer")
        if not width > 0: raise Usage("width must be a positive integer")
        if width%2 is 0: width += 1
        # even intervals are not symmetric, so make them odd by +1
        sopts["width"] = width
        sopts["from_center"] = True

    return sopts

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
           ref=None, width=151, from_center=False, path_to_gbdb=None,
           sort=False, npeaks=None, bed=False,
           debug = False,
           **kwargs):
    if from_MACS_subpeaks:
        sort_item = lambda x: int(x[3]) 
        if from_center:
            center_coord = lambda x: int(x[4]) # the 5th col
            start_coord = lambda x: center_coord(x) - width/2
            end_coord = lambda x: center_coord(x) + width/2 + 1
        else:
            start_coord = lambda x: int(x[1])
            end_coord = lambda x: int(x[2])
    if from_MACS_xls:
        sort_item = lambda x: int(x[5])
        if from_center:
            center_coord = lambda x: int(x[1]) + int(x[4]) # start + 5th col
            start_coord = lambda x: center_coord(x) - width/2
            end_coord = lambda x: center_coord(x) + width/2 
        else:
            start_coord = lambda x: int(x[1])
            end_coord = lambda x: int(x[2]) - 1
    if from_MACS_bed or peaks:
        sort_item = lambda x: float(x[4])
        if from_center:
            center_coord = lambda x: (int(x[1]) + int(x[2])) / 2
            start_coord = lambda x: center_coord(x) - width/2
            end_coord = lambda x: center_coord(x) + width/2 + 1
        else:
            start_coord = lambda x: int(x[1])
            end_coord = lambda x: int(x[2])
    if summits:
        sort_item = lambda x: float(x[4]) 
        from_center = True # implied logically
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
    if debug:
        print_debug("Using", ref, "to extract sequences", "from", 
                    parsed_filename.input_file, "to", output_file)
        
    ref_genome = twobitreader.TwoBitFile(ref)
    output_handle = open(output_file, 'w')
    input_file = open(parsed_filename.input_file, 'rU')
    input_lines = (line.rstrip() for line in input_file)

    if debug: print_debug("Determining chromsome lengths for", ref)
    chrom_lengths = ref_genome.sequence_sizes()
    valid_chroms = chrom_lengths.keys()

    if bed: bed_file = open(bed_filename, 'w')
    if debug: print_debug("Now converting {!s}".format(bed_filename)) 
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
        if debug: print_debug("Sorting {!s} peaks".format(num_peaks))
        sorted_list = sorted(seq_list, key=operator.itemgetter(4), reverse=True)
        if npeaks is None or num_peaks < npeaks:
            npeaks = num_peaks
            if debug: print_debug("Writing results for \
all {!s} peaks".format(npeaks))
        else:
            if debug: print_debug("Writing results for \
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
    def __init__(self, filename, include_width_in_name=False, target_dir=None,
                 debug=False, *args, **kwargs):
        fext = os.path.splitext(filename)[1].lstrip(os.extsep)
        if not (fext == 'bed' or fext == 'xls'): raise InvalidFileException
        if include_width_in_name:
            target_dir = "peaks_{!s}bp.FASTA".format(kwargs['width'])
        super(FilenameParser, self).__init__(filename, debug=debug,
                                             target_dir=target_dir,
                                             *args, **kwargs)
        if kwargs['from_MACS_subpeaks']:
            if not self.input_file.endswith("_subpeaks.bed"):
                raise InvalidFileException
        elif kwargs['from_MACS_xls']:
            if not self.input_file.endswith("_peaks.xls"):
                raise InvalidFileException
        elif kwargs['from_MACS_bed']:
            if not self.protoname.endswidth("_peaks.bed"):
                raise InvalidFileException
        elif kwargs['peaks'] or kwargs['summits']:
            if not self.file_extension=="bed":
                raise InvalidFileException

if __name__=="__main__": main()

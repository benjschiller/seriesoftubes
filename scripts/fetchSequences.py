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
--detect-ref                autodetect the reference genome for each file
                                only works if you have a directory named
                                from.MACS with (this is the default)
--path-to-gbdb=foo          If gbdb is not in /gbdb or C:/gbdb,
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
import subprocess
import StringIO
import scripter
from scripter import Usage, print_debug, extend_buffer
scripter.SCRIPT_DOC = __doc__
scripter.SCRIPT_VERSION = "2.2"
scripter.SOURCE_DIR = 'fromBAM.MACS'
scripter.TARGET_DIR = 'peaks.FASTA'
scripter.ALLOWED_EXTENSIONS = ['bed', 'xls']
SOURCE_OPTS = ["from-MACS-subpeaks", "from-MACS-xls",
               "from-MACS-bed", "peaks", "summits"]
BOOLEAN_OPTS = ["include-width-in-name", "from-center", "sort", "bed",
                "detect-ref"] + SOURCE_OPTS
scripter.SCRIPT_LONG_OPTS = ["ref=", "hg19", "hg18", "mm9", "detect-ref",
                             "width=","npeaks="] + BOOLEAN_OPTS

global PATH_TO_UCSC_TOOLS
global PATH_TO_GBDB
PATH_TO_UCSC_TOOLS =  "/usr/local/ucsc-bin"
PATH_TO_GBDB = None

def path_to_executable(name, directory=PATH_TO_UCSC_TOOLS):
    return scripter.path_to_executable(name, directory=directory)

def find_2bit_file(ref, path_to_gbdb=PATH_TO_GBDB):
    if ref is None:
        raise Usage("No reference genome specified")
    if os.path.exists(ref):
        return ref
    else:
        new_path = os.path.join("gbdb", ref, os.extsep.join([ref, "2bit"]))
        if path_to_gbdb is not None:
            specified_path = os.path.join(path_to_gbdb, ref,
                                          os.extsep.join([ref, "2bit"]))
            if os.path.exists(specified_path):
                return specified_path
        if os.path.exists(os.sep + new_path):
            return os.sep + new_path
        elif os.path.exists("C:" + os.sep + new_path):
            return "C" + os.sep + new_path
        else:
            raise Usage("Could not find 2bit file", ref)

def check_script_options(options, debug=False):
    sopts = {}
   
    # Find all reference genomes and their paths
    options_count = sum([options.has_key("hg18"), options.has_key("hg19"),
                        options.has_key("mm9"), options.has_key("ref")])
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
        sopts["ref"] = find_2bit_file(ref)
    
    for option in BOOLEAN_OPTS:
        pyoption = "_".join(option.split("-"))
        sopts[pyoption] = options.has_key(option)

    source_count = sum([options.has_key(x) for x in SOURCE_OPTS])
    if not source_count == 1:
        raise Usage("Must specify exactly one of",
                    " ".join(["".join(["--", x]) for x in SOURCE_OPTS]))

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

def detect_reference(parsed_filename):
    # finish
    return find_2bit_file("hg19")

def action(parsed_filename, from_MACS_subpeaks=True, from_MACS_xls=False,
           from_MACS_bed=False, peaks=False, summits=False,
           ref=None, detect_ref=True, width=151, from_center=False,
           path_to_twobittofa=path_to_executable("twoBitToFa"),
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


    stdout_buffer =""
    # output file
    output_file = os.path.join(parsed_filename.output_dir,
                               os.extsep.join([parsed_filename.protoname,
                                               'fa']))
    if bed: bed_filename = os.path.join(parsed_filename.output_dir,
                                os.extsep.join([parsed_filename.protoname,
                                                'bed']))

    if ref is None:
        ref = detect_reference(parsed_filename)
    if debug:
        print_debug("Using", ref, "to extract sequences", "from", 
                    parsed_filename.input_file, "to", output_file)
    input_file = open(parsed_filename.input_file)
    input_lines = (line.rstrip() for line in input_file)

    if debug: print_debug("Determining chromsome lengths for", ref)
    chrom_lengths = get_chrom_lengths(ref, debug=debug)
    valid_chroms = chrom_lengths.keys()


    if bed: bed_file = open(bed_filename, 'w')
    if debug: print_debug("Now converting") 
    step = [path_to_twobittofa, "-seqList=/dev/stdin", ref, "/dev/stdout"]
    job = subprocess.Popen(step,
                           stdin=subprocess.PIPE,
                           stdout=open(output_file, 'w'),
                           stderr=subprocess.PIPE)
    # if we're going to sort, we'll have to hold the list in memory
    if sort:
        seq_list = []
        sort_list = []
        bed_list = []
    else:
        seq_list = job.stdin
    num_peaks = 0
    for line in input_lines:
        words = line.split()
        # determine chrom/contig
        chrom = words[0].lower()
        # skip if not a valid contig/chrom
        if not chrom in valid_chroms: continue
        # determine start and end
        start = str(start_coord(words))
        end = str(end_coord(words))
        # skip if we've already fallen off the chromosome
        if int(start) > chrom_lengths[chrom]: continue
        # truncate chrom ends
        start = str(max(int(start), 0))
        end = str(min(int(end), chrom_lengths[chrom]))
        num_peaks += 1
        if debug:
            if num_peaks%1000 == 0: print_debug(str(num_peaks))

        if bed:
            bed_line = "\t".join([chrom, start, end, "peak" + str(num_peaks),
                                  str(sort_item(words)), "+"]) + os.linesep
        if sort:
            seq_list.append("".join([chrom, ":",
                                         start, "-", end, os.linesep]))
            sort_list.append(sort_item(words))
            if bed: bed_list.append(bed_line)
        else:
            seq_list.writelines("".join([chrom, ":",
                                         start, "-", end, os.linesep]))
            if bed: bed_file.write(bed_line)

    if sort:
        if debug: print_debug("Sorting")
        sorterator = (sitem for sitem in sorted(zip(sort_list, seq_list)))
        if bed:
            sorterator = (sitem for sitem in 
                                sorted(zip(sort_list, seq_list, bed_list)))
        if npeaks is not None and npeaks < num_peaks:
            if debug: print_debug("Taking top", str(npeaks), "peaks")
            sorted_list = []
            for x in range(npeaks):
                sitem = sorterator.next()
                sorted_list.append(sitem[1])
                if bed: bed_file.write(sitem[2])
        else:
            if debug: print_debug("Using all", str(num_peaks), "peaks")
            sorted_list = [sitem[1] for sitem in sorterator]
            if bed: bed_file.writelines([sitem[2] for sitem in sorterator])
        (stdout_data, stderr_data) = job.communicate(
                                        os.linesep.join(sorted_list))
    else:
        (stdout_data, stderr_data) = job.communicate()

    stdout_buffer = extend_buffer(stdout_buffer, stderr_data)
    return stdout_buffer

def get_chrom_lengths(ref, path_to_twoBitInfo=path_to_executable("twoBitInfo"),
                      debug=False, **kwargs):
    chrom_lengths = {}
    job = subprocess.Popen([path_to_twoBitInfo, ref, '/dev/stdout'],
                            stdout=subprocess.PIPE)
    (stdout_data, stderr_data) = job.communicate()
    for chrom in stdout_data.split(os.linesep):
        if chrom.strip()=='': continue
        (chrom_name, chrom_size) = chrom.split()
        chrom_lengths[chrom_name] = int(chrom_size)
        if debug: print_debug(chrom_name, "has length", chrom_size)
    return chrom_lengths

class FilenameParser(scripter.FilenameParser):
    def __init__(self, filename, include_width_in_name=False,
                 debug=False, *args, **kwargs):
        if include_width_in_name:
            target_dir = "".join(["peaks", str(kwargs['width']), ".FASTA"])
        else:
            target_dir = None
        super(FilenameParser, self).__init__(filename, debug=debug,
                                             target_dir=target_dir,
                                             *args, **kwargs)

        if not self.is_dummy_file:
            if kwargs['from_MACS_subpeaks']:
                if not self.input_file.endswith("_subpeaks.bed"):
                    self.is_invalid = True
            elif kwargs['from_MACS_xls']:
                if not self.input_file.endswith("_peaks.xls"):
                    self.is_invalid = True
            elif kwargs['from_MACS_bed']:
                if not self.protoname.endswidth("_peaks.bed"):
                    self.is_invalid = True
            elif kwargs['peaks'] or kwargs['summits']:
                if not self.file_extension=="bed": self.is_invalid = True

if __name__=="__main__":
    scripter.check_script_options = check_script_options
    scripter.perform(action, FilenameParser)

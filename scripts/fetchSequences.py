#!/usr/bin/env python
"""
fetch the FASTA sequences for the peaks

requires 2bit file and the twoBitToFa tool
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
                                 
Options for --peaks:
--from-center               Grab some interval of sequence from center of peak
                                (for use with --width)

Common options:
--width=x                   Get sequence x distance from the peak center/summit
                                (implies --from-center when applicable)
                                (default: 150)
--include-width-in-name     Include width in directory name for output
                                (e.g. peaks_150.FASTA)
"""
import os
import subprocess
import scripter
from scripter import Usage
scripter.SCRIPT_DOC = __doc__
scripter.SCRIPT_VERSION = "2.2"
scripter.SOURCE_DIR = 'fromBAM.MACS'
scripter.TARGET_DIR = 'peaks.FASTA'
scripter.ALLOWED_EXTENSIONS = ['bed', 'xls']
SOURCE_OPTS = ["from-MACS-subpeaks", "from-MACS-xls",
               "from-MACS-bed", "peaks", "summits"]
BOOLEAN_LONG_OPTS = SOURCE_OPTS + ["include-width-in-name", "from-center"]
scripter.SCRIPT_LONG_OPTS = ["ref=", "hg19", "hg18", "mm9", "detect-ref"
                             "width="] + BOOLEAN_OPTS

PATH_TO_TWOBITTOFA =  "/usr/local/ucsc-bin/twoBitToFa"

def find_2bit_file(ref):
    if ref == None:
        raise Usage("No reference genome specified")
    if os.path.exists(ref):
        return ref
    else:
        new_path = os.path.join("gbdb", ref, os.extsep.join([ref, "2bit"]))
        if os.path.exists(new_path):
            return new_path
        else:
            raise Usage("Could not find 2bit file", ref)

def check_script_options(options, debug=False):
    sopts = {}
   
    # Find all reference genomes and their paths
    option_count = sum([options.has_key("hg18"), options.has_key("hg19"),
                        options.has_key("mm9"), options.has_key("ref")])
    ref = None
    if options_count == 0:
        # ref is None
    elif options_count > 1:
        raise Usage("More than one reference genome specified")
    else:
        elif options.has_key("hg18"): ref = "hg18"
        elif options.has_key("hg19"): ref = "hg19"
        elif options.has_key("mm9"): ref = "mm9"
        elif options.has_key("ref"): ref = options["ref"]
        sopts["ref"] = find_2bit_file(ref)
    
    for option in BOOLEAN_LONG_OPTS:
        pyoption = "_".join(option.split("-"))
        sopts[pyoption] = options.has_key(option)

    source_count = sum([options.has_key(x)] for x in SOURCE_OPTS)
    if not source_count == 1:
        raise Usage("Must specify exactly one of"
                    " ".join(["".join(["--", x]) for x in SOURCE_OPTS]))

    if options.has_key("width"):
        try: width = int(options["width"]
        except ValueError: raise Usage("width must be an integer")
        if not width > 0: raise Usage("width must be a positive integer")
        sopts["width"] = width
        sopts["from_center"] = True

    return specific_options

def action(parsed_filename, from_MACS_subpeaks=True, from_MACS_xls=False,
           from_MACS_bed=False, peaks=False, summits=False,
           ref=None, width=150, from_center=False,
           **kwargs):
    stdout_buffer =""
    stdout_buffer = parsed_filename.input_file
    return stdout_buffer


class FilenameParser(scripter.FilenameParser):
    def __init__(self, filename, include_width_in_name=False,
                 debug=False,
                 *args, **kwargs):
        if include_width_in_name:
            target_dir = "".join(["peaks", str(kwargs['width']), ".FASTA"])
        else
            target_dir = None
        super(FilenameParser, self).__init__(filename, debug=debug,
                                             target_dir=target_dir
                                             *args, **kwargs)
        if not self.is_dummy_file:
            if kwargs['from_MACS_subpeaks']:
                if not self.protoname.endswidth("_subpeaks.bed"):
                    self.is_invalid = True
            elif kwargs['from_MACS_xls']
                if not self.input_file.endswith("_peaks.xls")
                    self.is_invalid = True
            elif kwargs['from_MACS_bed']
                if not self.protoname.endswidth("_peaks.bed")
                    self.is_invalid = True
            elif kwargs['peaks'] or kwargs['summits']:
                if not self.file_extension=="bed": self.is_invalid = True

if __name__=="__main__":
    scripter.check_script_options = check_script_options
    scripter.perform(action, FilenameParser)

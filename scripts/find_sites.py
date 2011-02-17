#!/usr/bin/env python
"""
findSites searches a pairs of BED (or MACS xls) and FASTA files
for matches to a TFBS motif
*** findSites absolutely requires bioplus ***

Required flags
--motif=foo             Fetch the motif from file FOO

Optional flags:
--bed                   Assume we have a bed file (default)
--xls                   Assume we have an xls file
--motif-type=[MEME]     Type of motif file
                            Default: MEME
--motif-number=n        Motif number within file (e.g. 1, 2, 3, ...)
                            Default: Assume only one motif.
"""
import os
import scripter
from scripter import Usage, print_debug, extend_buffer, InvalidFileException
import Bio
import Bio.Motif
import bioplus
import bioplus.sitefinder
from bioplus.sitefinder import find_sites
try: import MOODS
except ImportError: print_debug("WARNING: MOODS is not installed")

SCRIPT_VERSION = "2.4"

def main():
    boolean_opts = ["bed", "xls"]
    long_opts = ["motif=", "motif-type=", "motif-number="] + boolean_opts
    e = scripter.Environment(long_opts=long_opts, doc=__doc__, version=VERSION)
    e.parse_boolean_opts(boolean_opts)
    all_opts = check_script_options(e.get_options(), debug=e.is_debug())
    e.update_script_kwargs(all_opts)
    e.set_source_dir('peaks.FASTA')
    e.set_target_dir('sites.analysis')
    e.set_filename_parser(FilenameParser)
    e.do_action(action)
    
def check_script_options(options, debug=False):
    sopts = {}
  
    # motif flag is REQUIRED
    if not options.has_key("motif"):
        raise Usage("--motif=[location of motif file] is required")
    else:
        motif_file = options['motif']
        if not os.path.exists(motif_file):
            raise Usage(motif_file, "not found")
        if options.has_key("motif-type"):
            motif_type = options['motif-type']
        else:
            motif_type = "MEME"
        if options.has_key("motif_number"):
            motif_number = options["motif_number"]
            try:
                motifs = Bio.Motif.parse(open(motif_file), motif_type)
            except ValueError:
                raise Usage(motif_file, "does not contain a valid motif")
            for i in range(1, motif_number):
                try:
                    motif = motifs.next()
                except StopIteration:
                    raise Usage(motif_file, "only contains", str(i), "motifs")
        else:
            try:
                motif = Bio.Motif.read(open(motif_file), motif_type)
            except ValueError:
                raise Usage(motif_file, "does not contain a valid motif")
        sopts['motif'] = motif
   
    if sopts['bed'] and sopts['xls']:
        raise Usage("can only specify one of --bed, --xls")
    elif not sopts['bed'] and not sopts['xls']:
        sopts['bed'] = True

    return sopts

def action(parsed_filename, motif, bed=True, xls=False,
           silent=False, debug=False, **kwargs):
    if debug: print_debug("trying to find sites for",
                          parsed_filename.input_file)
    stdout_buffer = find_sites(parsed_filename.input_file,
                               parsed_filename.fasta_file,
                               motif, bed=bed, xls=xls,
                               output_dir=parsed_filename.output_dir,
                               src_fnc = "findSites.py",
                               debug=debug)
    if not silent: return stdout_buffer

class FilenameParser(scripter.FilenameParser):
    def __init__(self, filename, include_width_in_name=False,
                 debug=False, *args, **kwargs):
        fext = os.path.splitext(filename)[0].lstrip(os.extsep)
        if not (fext == 'bed' or fext == 'xls'): raise InvalidFileException
        super(FilenameParser, self).__init__(filename, debug=debug,
                                             *args, **kwargs)
        if self.is_dummy_file: return
        self.fasta_file = None
        for file_extension in ['fa', 'fasta', 'FA', 'FASTA']:
            fasta_file = os.path.join(self.input_dir,
                            os.extsep.join([self.protoname, file_extension]))
            if debug: print_debug("Trying", fasta_file)
            if os.path.exists(fasta_file):
                self.fasta_file = fasta_file
                if debug: print_debug("Using", fasta_file)
                break
        if self.fasta_file is None:
            raise Usage("Could not find the FASTA file for ", self.input_file)
            # or in the future, fetch the sequence?

if __name__=="__main__": main()
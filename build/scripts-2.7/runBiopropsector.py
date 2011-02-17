#!C:\Python27\python.exe
'''
finds a motif using Bioprospector
default source: peaks.FASTA
default target: motifs.Bioprospector

--GR                        use settings for the glucocorticoid receptor
                            i.e., -W 17, -a 1
--width=#                   assume motif width is #
                            (default: 10)
--second-width=#            assume second half of motif has width #
                            (default: 0)
--gap=                      assume the gap between two blocks is eactly #
--min-gap=                  assume the gap between two blocks is at least #
--max-gap=                  assume the gap between two blocks is at most #
--palindrome                assume motif is a palindrome of two blocks,
                            same width
--background-file=foo       use FASTA file foo too compute background model
--pre-background-file=foo   use precomputed background model in file foo
--background-only           only compute the background model,
                            do not find motifs
--background                compute and save the background model
                                (otherwise I think bioprospector does it
                                 on the fly)
--oops                      assume exactly one motif per sequence
--num-iters=#               try to find each motif # times
--num-motifs=#              report the # best motifs
                            num-motifs can be at most num-iters
--degenerate                assume the motif is very degenerate
--no-reverse-complement     do not search reverse complement of sequence

-----
from Bioprospector (these flags will be based along)
'''
import os
import subprocess
import StringIO
import scripter
import Bio.SeqIO
from scripter import Usage, valid_int, print_debug, extend_buffer
scripter.SCRIPT_DOC = __doc__
scripter.SCRIPT_VERSION = "2.2"
scripter.SOURCE_DIR = 'peaks.FASTA'
scripter.TARGET_DIR = 'motifs.Bioprospector'
scripter.ALLOWED_EXTENSIONS = ['fa','fasta','FA','FASTA']
BOOLEAN_OPTS = ["GR", "degenerate", "palindrome", "oops", "background-only",
                "background", "no-reverse-complement"]
scripter.SCRIPT_LONG_OPTS = ["width=", "second-width=",
                             "background-file=", "pre-background-file=",
                             "gap=", "min-gap=", "max-gap="
                             "num-iters=", "num-motifs="] + BOOLEAN_OPTS

global PATH_TO_BIN
global BIOPROSPECTOR_NAME
global GENOMEBG_NAME
PATH_TO_BIN = "/usr/local/bin"
BIOPROSPECTOR_NAME = "Bioprospector"
GENOMEBG_NAME = "genomebg"

def path_to_executable(name, directory=PATH_TO_BIN):
    return scripter.path_to_executable(name, directory=directory)

def valid_gap(gap, gap_type=""):
    """
    checks if a gap, min-gap, or max-gap is valid
    returns the gap as an int
    """
    msg = ' '.join([gap_type, " must be an integer 0 - 50"])
 
    return valid_int(gap, msg, 0, 50)

def valid_width(width, width_type=""):
    """
    checks if a width, second-width is valid
    returns the width as an integer
    """
    msg = ' '.join([width_type, "width must be an integer 4 - 50"])
    return valid_int(width, msg, 4, 50)

def valid_num_iters(num_iters):
    msg = "num-iters must be an integer 1 - 200"
    return valid_int(num_iters, msg, 1, 200)

def valid_num_motifs(num_motifs, num_iters):
    msg = "num-motifs must be an integer 1 - " + num_iters
    return valid_int(num_motifs, msg, 1, num_iters)

def check_script_options(options, debug=False):
    bpopts = {}

    for option in BOOLEAN_OPTS:
        pyoption = "_".join(option.split("-"))
        bpopts[pyoption] = options.has_key(option)

    if bpopts["background_only"]: return bpopts

    # parse width, second_width
    if options.has_key("width"):
        width = valid_width(options["width"])
        bpopts["width"] = str(width)
    else:
        bpopts["width"] = "10"
    if options.has_key("second-width"):
        second_width = valid_width(options["second-width"])
        bpopts["width"] = str(second_width)
    else:
        bpopts["width"] = "0"

    # parse gap, min_gap, max_gap 
    valid_gap_diff_msg = "max-gap - min-gap must be smaller than 20"
    if (options.has_key("min-gap") or
      options.has_key("max-gap")) and options.has_key("gap"):
        raise Usage("Do not use --gap with --min-gap or --max-gap")
    elif options.has_key("gap"):
        gap = valid_gap(options["gap"], "gap")
        bpopts["min_gap"] = str(gap)
        bpopts["max_gap"] = str(gap)
    elif options.has_key("min-gap") and options.has_key("max-gap"):
        min_gap = valid_gap(options["min-gap"], "min-gap")
        max_gap = valid_gap(options["max-gap"], "max-gap")
        if max_gap - min_gap >= 20:
            raise Usage("max-gap - min-gap must be smaller than 20")
        bpopts["min_gap"] = str(min_gap)
        bpopts["max_gap"] = str(max_gap)
    elif options.has_key("max-gap"):
        max_gap = valid_gap(options["max-gap"], "max-gap")
        min_gap = max(max_gap - 19, 0)
        bpopts["min_gap"] = str(min_gap)
        bpopts["max_gap"] = str(max_gap)
    elif options.has_key("min-gap"):
        min_gap = valid_gap(options["min-gap"], "min-gap")
        max_gap = min(min_gap + 19, 50)
        bpopts["min_gap"] = str(min_gap, 50)
        bpopts["max_gap"] = str(max_gap)
    else:
        bpopts["min_gap"] = "0" 
        bpopts["max_gap"] = "0"

    # check background file
    if options.has_key("background-file") and \
        options.has_key("pre-background-file"):
        raise Usage("Do not use both background-file and pre-backround-file")
    elif options.has_key("background-file"):
        bg_file = options["background-file"]
        if not os.path.exists(bg_file):
            raise Usage(bg_file, "does not exist or could not be read")
        bpopts["background_file"] = bg_file
    elif options.has_key("pre-background-file"):
        pre_bg_file = options["pre-background-file"]
        if not os.path.exists(pre_bg_file):
            raise Usage(pre_bg_file, "does not exist or could not be read")
        bpopts["pre_background_file"] = pre_bg_file

    # check how many iterations
    if options.has_key("num-iters"):
        num_iters = valid_iters(options["num-iters"])
        bpopts["num_iters"] = str(num_iters)
    else:
        num_iters = 40
        bpopts["num_iters"] = str(num_iters)

    # check how many motifs to report
    if options.has_key("num-motifs"):
        num_motifs = valid_num_motifs(options["num-motifs"], num_iters)
        bpopts["num_motifs"] = str(num_motifs)
    elif num_iters < 5:
        bpopts["num_motifs"] = str(num_motifs)

    return bpopts

def _oj(a,b):
    """"
    option join
    a,b -> -a b
    """
    return "".join(["-", a, " ", b])

def action(parsed_filename, background=False, background_only="False",
           path_to_bioprospector=path_to_executable(BIOPROSPECTOR_NAME),
           pre_background_file=None, background_file=None,
           debug=False, width="10", second_width="0",
           min_gap="0", max_gap="0",
           num_iters="40", num_motifs="5", oops=False, 
           degenerate = False, no_reverse_complement=False, palindrome=False,
           **kwargs):
    """
    """
    stdout_buffer =  ""
    if parsed_filename.protoname.endswith("_subpeaks"):
        pass
        # do something?
    if background_only: return compute_background(parsed_filename)
    if background:
        stdout_buffer = extend_buffer(stdout_buffer,
                                      compute_background(parsed_filename), 2)
    if GR:
        oops = True
        width = "17"

    bp_options = [_oj("W",width), _oj("w", second_width),
                  _oj("n", num_iters), _oj("r", num_motifs),
                  _oj("G", max_gap), _oj("g", min_gap)]
    if pre_background_file is not None:
        bp_options.append(_oj("f", pre_background_file))
    if _background_file is not None:
        bp_options.append(_oj("b", background_file))
    if palindrome: bp_options.append(_oj("p", 1))
    if no_rc: bp_options.append(_oj("d", 1))
    if oops: bp_options.append(_oj("a", 1))
    if degenerate: bp_options.append(_oj("h", 1))
    return stdout_buffer

def compute_background(parsed_filename,
                       path_to_genomebg=path_to_executable(GENOMEBG_NAME)):
    """compute the background model with genomebg"""
    stdout_buffer =  ""

    bg_output_file = os.path.join(parsed_filename.output_dir,
                                  os.extsep.join([parsed_filename.protoname,
                                                  "background"]))

    input_ = StringIO.StringIO()
    fasta_records = Bio.SeqIO.parse(open(parsed_filename.input_file), "fasta")
    for r in fasta_records:
        input_.writelines([">", r.description, os.linesep,
                           str(r.seq), os.linesep])
    input_.seek(0)
    step = path_to_genomebg + ["-i", "/dev/stdin"]
    job = subprocess.Popen(step,
                           stdin=input_,
                           stdout=open(bg_output_file, 'w'),
                           stderr=subprocess.PIPE)
    (stdout_data, stderr_data) = job.communicate()
    stdout_buffer = extend_buffer(stdout_buffer, ' '.join(step))
    stdout_buffer = extend_buffer(stdout_buffer, stderr_data, 1)

    return stdout_buffer

if __name__=="__main__":
    scripter.check_script_options = check_script_options
    scripter.perform(action)

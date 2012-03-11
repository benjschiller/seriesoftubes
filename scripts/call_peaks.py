#!/usr/bin/env python
'''
calls peaks using MACS v2
use --config to specify file with matched sample/controls
'''
from ConfigParser import ConfigParser
import sys
from errno import ENOENT, EACCES
from os import access, extsep, strerror, R_OK, getcwd
from os.path import exists, join
import platform
from subprocess import Popen, STDOUT, PIPE
import pysam
from scripter import Environment, path_to_executable, get_logger, Usage, exit_on_Usage
from seriesoftubes.fnparsers import BAMFilenameParser
from bioplus.genometools import guess_bam_genome, genome, NoMatchFoundError
from pkg_resources import get_distribution, VersionConflict
__version__ = get_distribution('seriesoftubes').version
VERSION = __version__
TARGET_DIR = 'peaks'
MACS_VERSION = get_distribution('MACS>=2.0.10').version

def main():
    e = Environment(doc=__doc__, version=VERSION)
    parser = e.argument_parser
    parser.add_argument('-g', '--genome-size', dest='user_gsize', default=None,
                        help='Optional user-specified genome size')
    parser.add_argument('--path-to-macs',
                        default=path_to_executable("macs2"),
                        help="optional path to macs2 executable")
    parser.add_argument('--path-to-R',
                        default=path_to_executable(["R64", "R"]),
                        help="optional path to R executable")
    parser.add_argument('--no-subpeaks', dest='subpeaks', action='store_false',
                        default=True,
                        help='do not call subpeaks with --call-summits')
    parser.set_defaults(**{'target': 'peaks'})
    e.set_filename_parser(BAMFilenameParser)
    e.set_config_reader(read_setup_file)
    e.set_config_writer(write_setup_file)
    e.do_action(run_macs, stay_open=True)
    sys.exit()
    
def write_setup_file(controls=None, *args, **kwargs):
    setup_file = open(join(TARGET_DIR, 'setup.txt'), 'w')
    for sample, value in controls.items():
        name, control = value
        if control == None: control = 'None'
        record = '[%s]\nsample = %s\ncontrol = %s\n' % (name, sample, control)
        setup_file.write(record)
    setup_file.close()

def read_setup_file(setup_file):
        if not exists(setup_file):
            raise IOError(ENOENT, strerror(ENOENT), setup_file)
        if not access(setup_file, R_OK):
            raise IOError(EACCES, strerror(EACCES), setup_file)
        config_parser = ConfigParser()
        config_parser.readfp(open(setup_file, 'rU'))
        controls = {} # {sample: (name, control)}
        for section in config_parser.sections():
            name = section
            sample = config_parser.get(section, 'sample')
            if config_parser.has_option(section, 'control'):
                control = config_parser.get(section, 'control')
                if control.strip() == 'None': control = None
            else:
                control = None
            controls[sample] = (name, control)
        return {'controls': controls}
            

def decide_format(input_file, control_file, logger=None): 
    # See if we have paired-end files
    s = pysam.Samfile(input_file)
    is_paired = [s.next().is_paired for i in xrange(100000)]
    if control_file is not None:
        t = pysam.Samfile(control_file)
        is_paired_control = [s.next().is_paired for i in xrange(100000)]
    else:
        is_paired_control = [True]
    if all(is_paired) and all(is_paired_control):
        if logger is not None:
            logger.warn('Detected paired end files')
            logger.warn('Using new BAMPEParser instead of BAMParser')
        return 'BAMPE'
    else: return 'BAM'


@exit_on_Usage
def run_macs(fp_obj, subpeaks=True, path_to_macs=None, logging_level=10,
             user_gsize=None,
             **kwargs):
    """Run MACS on a BAM file
    """
    logger = get_logger(logging_level)
    if path_to_macs is None:
        path_to_macs = path_to_executable(["macs14", "macs14.py", "macs"])

    input_file = fp_obj.input_file
    control_file = fp_obj.control_file    
    logger.debug('Processing', input_file)
    if control_file is not None: logger.debug('with control', control_file)

    # determine genome name and size
    if user_gsize: genome_size = user_gsize
    else:
        try: genome_build = guess_bam_genome(input_file)
        except NoMatchFoundError:
            raise Usage('Could not determine genome size for %s' % bam_file)
        gname = ''.join([x for x in genome_build if x.isalpha()])
        if gname == 'hg': return 'hs'
        elif gname in ['mm', 'ce', 'dm']: return gname
        else: return '%.1e' % sum(genome(genome_build).itervalues())
    
    fmt = decide_format(input_file, control_file, logger)
    macs_options = ['-f %s' % fmt, # correct file format BAM or BAMPE
                    '-B', #bedgraph
                    '-S', #whole genome
                    '-g %s' % genome_size,
                    '-n %s' % fb_obj.run_name.replace(' ', '_'), # run name
                    '-t %s' % input_file] # treatment
    if control_file is not None:
        macs_options.append('-c %s' % control_file)
    if subpeaks: macs_options.append('--call-summits')

    step = [path_to_macs] + macs_options
    if platform.system() is 'Windows': step.insert(sys.executable, 0)
    
    logger.debug('Launching', ' '.join(step))
    cwd = join(getcwd(), fp_obj.output_dir)
    job = Popen(step, stdout=PIPE, stderr=STDOUT, cwd=cwd)

    stdout_buffer = '%s\n\n%s\n' % (' '.join(step), job.communicate()[0])
    logger.debug(stdout_buffer)
    
    return stdout_buffer

if __name__=="__main__": main()

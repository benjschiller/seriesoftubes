#!/usr/bin/env python
'''
calls peaks using MACS
use --config to specify file with matched sample/controls
'''
from ConfigParser import ConfigParser
import sys
import os
from errno import ENOENT, EACCES
from os import access, extsep, strerror, R_OK
from os.path import exists
import platform
import subprocess
import scripter
import pysam
from scripter import Environment, path_to_executable, get_logger
from seriesoftubes.fnparsers import BAMFilenameParser
from pkg_resources import get_distribution, VersionConflict
__version__ = get_distribution('seriesoftubes').version
VERSION = __version__
TARGET_DIR = 'peaks'

def main():
    e = Environment(doc=__doc__, version=VERSION)
    parser = e.argument_parser
    parser.add_argument('--path-to-macs',
                        default=path_to_executable(["macs2",
                                                    "macs14",
                                                    "macs14.py",
                                                    "macs"]),
                        help="optional path to macs executable (macs2 will be \
used if found)")
    parser.add_argument('--path-to-R',
                        default=path_to_executable(["R64", "R"]),
                        help="optional path to R executable")
    parser.add_argument('--no-diag', dest='diag', action='store_false', default=True,
                        help='disable generation of diagnostic xls file')
    parser.add_argument('--no-subpeaks', dest='subpeaks', action='store_false',
                        default=True,
                        help='do not call subpeaks with PeakSplitter')
    g = parser.add_mutually_exclusive_group()
    g.add_argument('--no-fix', dest='fix', action='store_false', default=True,
                   help='do not auto-fix subpeaks filename')
    g.add_argument('--fix-only', action='store_true', default=True,
                   help='only fix subpeaks filename, do nothing else')
    parser.set_defaults(**{'target': 'peaks'})
    e.set_filename_parser(BAMFilenameParser)
    e.set_config_reader(read_setup_file)
    e.set_config_writer(write_setup_file)
    e.do_action(run_macs, stay_open=True)
    sys.exit()
    
def write_setup_file(controls=None, *args, **kwargs):
    setup_file = open(os.path.join(TARGET_DIR, 'setup.txt'), 'w')
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
            
def _prepend_cwd(d):
    '''prepends the current working directory to a directory d'''
    return os.path.join(os.getcwd(), d)

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

def run_macs(fp_obj,
           diag=True, subpeaks=True,
           fix_only=False, fix=True, 
           path_to_macs=None, path_to_R=None, logging_level=10,
           **kwargs):
    """Run MACS on a BAM file and produce the pdf from the .R model"""
    is_v2 = os.path.basename(kwargs['path_to_macs'])=='macs2'
    logger = get_logger(logging_level)
    if fix_only:
        return fix_subpeaks_filename(fp_obj, logger=logger)
    if path_to_macs is None:
        path_to_macs = path_to_executable(["macs14", "macs14.py", "macs"])

    input_file = fp_obj.input_file
    control_file = fp_obj.control_file
    if is_v2:
        fmt = decide_format(input_file, control_file, logger)
    else:
        fmt = 'BAM'
    logger.debug('Processing', input_file)
    if control_file is not None:
        logger.debug('with control', control_file)

    macs_options = ['--format %s' % fmt]
    macs_options.append('--bdg')
    macs_options.append('--single-profile')
    if subpeaks: macs_options.append('--call-subpeaks')
    if diag: macs_options.append('--diag')
    run_name = '_'.join(fp_obj.run_name.split())
    macs_options.append('--name %s' % run_name)
    macs_options.append('--treatment %s' % _prepend_cwd(input_file)) 
    if control_file is not None:
        macs_options.append('--control %s' % _prepend_cwd(control_file))

    step = [path_to_macs] + macs_options
    if platform.system() is 'Windows': step.insert(sys.executable, 0)
    
    logger.debug('Launching', ' '.join(step))
    job = subprocess.Popen(step,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT,
                           cwd=_prepend_cwd(fp_obj.output_dir))

    stdout_buffer = '%s\n\n%s\n' % (' '.join(step), job.communicate()[0])
    logger.debug(stdout_buffer)
    if subpeaks and fix: fix_subpeaks_filename(fp_obj)
    
    R_file = os.path.join(os.getcwd(), fp_obj.output_dir,
                          ''.join([fp_obj.protoname, '_model.r']))
    R_output = Rmodel_to_pdf(R_file, path_to_R=path_to_R, logger=logger)
    return '%s\n%s' % (stdout_buffer, R_output)

def Rmodel_to_pdf(R_file, path_to_R=None, logger=logger):
    if path_to_R is None:
        try:
            path_to_R = path_to_executable(["R64", "R"], max_depth=3)
        except scripter.Usage:
            msg = 'R is not installed. Could not process model pdf'
            logger.error(msg)
            return '%s\n' % msg
    if not os.path.exists(R_file):
        logger.warn('Warning: %s does not exist', R_file)
        return 'Nothing to do\n'
    
    logger.debug('Processing %s', R_file)
    job = subprocess.Popen([path_to_R, '--vanilla'], 
                           stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT,
                           stdin=open(R_file),
                           cwd=_prepend_cwd(parsed_filename.output_dir))
    stdout_buffer = '%s\n\n%s\n' % (' '.join(step), job.communicate()[0])
    logger.debug(stdout_buffer)
    return stdout_buffer

def fix_subpeaks(parsed_filename, logger=None):
    """
    there's a few problems with the subpeaks file that we fix here
    
    1) Breaks MACS naming convention, rename it
    2) First line is not compatible with BED format, delete it
    """
    if logger is None: logger = get_logger()
    fn_parts = parsed_filename.protoname.split(extsep)
    path_to_subpeaks = os.path.join(parsed_filename.output_dir,
                        extsep.join([fn_parts[0] + '_peaks'] + \
                                        fn_parts[1:] + ['subpeaks', 'bed']))
    new_path_to_subpeaks = os.path.join(parsed_filename.output_dir,
                            extsep.join([fn_parts[0]] + fn_parts[1:] +
                                            ['_subpeaks', 'bed']))
    # check if the file exists
    if os.path.exists(new_path_to_subpeaks):
        logger.error('Cannot move %s to %s. Latter file already exists',
                    path_to_subpeaks, new_path_to_subpeaks)
        return
    try: old_file = open(path_to_subpeaks)
    except EnvironmentError:
        logger.error('Could not open %s for reading', path_to_subpeaks)
        return  
    try: new_file = open(path_to_subpeaks, 'w')
    except EnvironmentError:
        logger.error('Could not open %s for writing', new_path_to_subpeaks)
        return
    try:
        old_file.readline()
        new_file.writelines(old_file.readlines())
    except EnvironmentError:
        logger.error('Something went wrong copying %s to %s',
                     path_to_subpeaks, new_path_to_subpeaks)
        return
        old_file.close()
        new_file.close()
    try: os.remove(path_to_subpeaks)
    except EnvironmentError:
        logger.warn('Copied %s to %s but could not remove original',
                    path_to_subpeaks, new_path_to_subpeaks)
        return
    logger.info('Moved %s to %s successfully',
                path_to_subpeaks, new_path_to_subpeaks)
    return

if __name__=="__main__": main()

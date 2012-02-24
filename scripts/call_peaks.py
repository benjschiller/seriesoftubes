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
from seriesoftubes.fnparsers import AlignmentsFilenameParser
from pkg_resources import get_distribution, VersionConflict
__version__ = get_distribution('seriesoftubes').version
VERSION = __version__
TARGET_DIR = 'peaks'

def enforce_MACS_version():
    """
    Make sure we MACS > v1.4.1 and/or >v2.0
    """
    v1 = None
    v2 = None
    try:
        v1 = get_distribution('MACS>=1.4.1,<2').version
    except VersionConflict:
        pass
    try:
        v2 = get_distribution('MACS>2').version
    except VersionConflict:
        pass
    if v1 is None and v2 is None:
        raise VersionConflict()
    return v1, v2

#MACS_VERSION1, MACS_VERSION2 = check_MACS_version()  

def main():
    e = Environment(doc=__doc__, version=VERSION)
    parser = e.argument_parser
    parser.add_argument('--path-to-macs',
                        default=path_to_executable(["macs2",
                                                    "macs14",
                                                    "macs14.py",
                                                    "macs"]),
                        help="optional path to macs executable (macs2 will be \
used if found")
    parser.add_argument('--path-to-R',
                        default=path_to_executable(["R64", "R"]),
                        help="optional path to R executable")
    parser.add_argument('--no-wig', dest='wig', action='store_false',
                        default=True,
                        help='do not generate wiggle files (also disable subpeaks')
    parser.add_argument('--chrom-wigs', dest='single_wig', action='store_false',
                        default=True,
                        help='generate a wig file per chromosome instead of a single wig')
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
    e.set_filename_parser(AlignmentsFilenameParser)
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
        return controls
            
def _prepend_cwd(d):
    '''prepends the current working directory to a directory d'''
    return os.path.join(os.getcwd(), d)

def _join_opt(a,b):
    '''joins an option with its value
    _join_opt('--option,'value') -> '--option=value'
    '''
    return '='.join([a,b])

def run_macs(*args, **kwargs):
    if os.path.basename(kwargs['path_to_macs'])=='macs2':
        return run_macs2(*args, **kwargs)
    else:
        return run_macs14(*args, **kwargs)
    
def run_macs2(parsed_filename, logging_level=10, **kwargs):
    raise NotImplementError
    logger = get_logger(logging_level)
    input_file = parsed_filename.input_file
    control_file = parsed_filename.control_file

    # See if we have paired-end files
    s = pysam.Samfile(input_file)
    is_paired = [s.next().is_paired for i in xrange(100000)]
    if control_file is not None:
        t = pysam.Samfile(control_file)
        is_paired_control = [s.next().is_paired for i in xrange(100000)]
    else:
        is_paired_control = [True]
    if all(is_paired + is_paired_control):
        logger.warn('Detected paired end files')
        logger.debug('Checking for pybedtools')
        try:
            import pybedtools
            run_macs2_with_bedgraph(parsed_filename, **kwargs)
        except ImportError:
            logger.error('pybedtools is not installed. Falling back on MACS fixed-width track features')
            run_macs2_alone(parsed_filename, **kwargs)
    else:
        run_macs2_alone(parsed_filename, **kwargs)
    return

def run_macs14(parsed_filename,
           wig=True, single_wig=True, diag=True, subpeaks=True,
           fix_only=False, fix=True, 
           path_to_macs=None, path_to_R=None, logging_level=10,
           **kwargs):
    """Run MACS on a BAM file and produce the pdf from the .R model"""
    logger = get_logger(logging_level)
    if fix_only:
        return fix_subpeaks_filename(parsed_filename)

    input_file = parsed_filename.input_file
    control_file = parsed_filename.control_file
    logger.debug('Processing', input_file)
    if control_file is not None:
        logger.debug('with control', control_file)

    macs_options = ['--format=BAM']
    if wig:
        macs_options.append('--wig')
        if single_wig: macs_options.append('--single-profile')
        if subpeaks: macs_options.append('--call-subpeaks')
    if diag: macs_options.append('--diag')
    run_name = '_'.join(parsed_filename.run_name.split())
    macs_options.append(_join_opt('--name', run_name))
    macs_options.append(_join_opt('--treatment',
                                  _prepend_cwd(input_file)))
    if parsed_filename.control_file is not None:
        macs_options.append(_join_opt('--control',
                                _prepend_cwd(control_file)))

    if path_to_macs is None:
        path_to_macs = path_to_executable(["macs14", "macs14.py", "macs"])
    if platform.system() is 'Windows':
        step = [sys.executable, path_to_macs] + macs_options
    else:
        step = [path_to_macs] + macs_options
    
    logger.debug('Launching', ' '.join(step))
    job = subprocess.Popen(step,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT,
                           cwd=_prepend_cwd(parsed_filename.output_dir))

    stdout_data = job.communicate()[0]
    stdout_buffer = '{0!s}\n\n{1!s}'.format(' '.join(step), stdout_data)
   
    if subpeaks and fix:
        output = fix_subpeaks_filename(parsed_filename)
        logger.debug('{0!s}\n{1!s}'.format(stdout_buffer, output))

    # now process R file
    if make_pdf:
        R_file = os.path.join(os.getcwd(), parsed_filename.output_dir,
                              ''.join([parsed_filename.protoname, '_model.r']))

        logger.debug('Processing', R_file)

        if not os.path.exists(R_file):
            logger.debug('Warning:', R_file, 'does not exist')
        else:
            R_pointer = open(R_file)
            if path_to_R is None:
                path_to_R = path_to_executable(["R64", "R"], max_depth=3)
            step =[path_to_R, '--vanilla'] 
            job = subprocess.Popen(step,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT,
                                   stdin=R_pointer,
                                   cwd=_prepend_cwd(parsed_filename.output_dir))
            stdout_data = job.communicate()[0]
            R_pointer.close()
            stdout_buffer = '{0!s}\n\n{1!s}\n{2!s}\n'.format(stdout_buffer,
                                                             ' '.join(step),
                                                             stdout_data)

    return stdout_buffer

def fix_subpeaks_filename(parsed_filename, debug=False):
    """
    there's a few problems with the subpeaks file that we fix here
    
    1) Breaks MACS naming convention, rename it
    2) First line is not compatible with BED format, delete it
    """
    fn_parts = parsed_filename.protoname.split(extsep)
    path_to_subpeaks = os.path.join(parsed_filename.output_dir,
                        extsep.join([fn_parts[0] + '_peaks'] + \
                                        fn_parts[1:] + ['subpeaks', 'bed']))
    new_path_to_subpeaks = os.path.join(parsed_filename.output_dir,
                            extsep.join([fn_parts[0]] + fn_parts[1:] +
                                            ['_subpeaks', 'bed']))
    # check if the file exists
    if os.path.exists(new_path_to_subpeaks):
        return 'Cannot move {0!s} to {1!s}. Latter file already exists'.format(
                                        path_to_subpeaks, new_path_to_subpeaks)
    try:
        old_file = open(path_to_subpeaks)
    except EnvironmentError:
        return 'Could not open {0!s} for reading'.format(path_to_subpeaks)
    try:
        new_file = open(path_to_subpeaks, 'w')
    except EnvironmentError:
        return 'Could not open {0!s} for writing'.format(new_path_to_subpeaks)
    try:
        old_file.readline()
        new_file.writelines(old_file.readlines())
    except EnvironmentError:
        return 'Something went wrong copying {0!s} to {1!s}'.format(
                                        path_to_subpeaks, new_path_to_subpeaks)
        old_file.close()
        new_file.close()
    try:
        os.remove(path_to_subpeaks)
    except EnvironmentError:
        return 'Copied {0!s} to {1!s} but could not remove original'.format(path_to_subpeaks,
                                                                          new_path_to_subpeaks)
    return 'Moved {0!s} to {1!s} successfully'

if __name__=="__main__": main()

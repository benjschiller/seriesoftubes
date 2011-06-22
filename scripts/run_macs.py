#!/usr/bin/env python
'''
calls peaks using MACS

--setup=setup.txt   specify file with matched sample/controls

--no-wig            do not generate wiggle files (also disable subpeaks)
--chrom-wigs        generate a wig file per chromosome instead of a single wig
--no-diag           disable generation of diagnostic xls file
--no-subpeaks       do not call subpeaks with PeakSplitter
--no-pdf            do not generate a pdf of the model using R
--no-fix            do not auto-fix subpeaks filename
--fix-only          only fix subpeaks filename, do nothing else
'''
from ConfigParser import ConfigParser
import sys
import os
from errno import ENOENT, EACCES
from os import access, extsep, strerror, R_OK
from os.path import exists
import platform
import subprocess
from scripter import Environment, FilenameParser, InvalidFileException, Usage, \
                     path_to_executable,  print_debug

VERSION = "2.4"

def main():
    long_opts = ["no-wig", "no-diag", "chrom-wigs", "no-subpeaks",
                 "no-pdf", "fix-only", "no-fix", "setup="]
    e = Environment(long_opts=long_opts, doc=__doc__, version=VERSION)
    path_to_macs = path_to_executable(["macs14", "macs14.py", "macs"])
    if e.is_debug(): print_debug('Found MACS at', path_to_macs)
    path_to_R = path_to_executable(["R64", "R"])
    if e.is_debug(): print_debug('Found R at', path_to_R)
    e.set_source_dir('aligned.BAM')
    e.set_target_dir('fromBAM.macs')
    e.update_script_kwargs(check_script_options(e.get_options()))
    e.update_script_kwargs({'path_to_R': path_to_R,
                            'path_to_macs': path_to_macs})
    e.set_filename_parser(MacsFilenameParser)
    controls = e.get_script_kwargs()['controls']
    write_setup_file(controls)
    e.do_action(run_macs, stay_open=True)
    
def write_setup_file(controls):
    setup_file = open(os.path.join('fromBAM.macs', 'setup.txt'), 'w')
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
            
def check_script_options(options):
    specific_options = {}

    specific_options['wig'] = not options.has_key('no-wig')
    specific_options['single_wig'] = not options.has_key('chrom-wigs')
    specific_options['diag'] = not options.has_key('no-diag')
    specific_options['subpeaks'] = not options.has_key('no-subpeaks')
    specific_options['make_pdf'] = not options.has_key('no-pdf')
    specific_options['fix_only'] = options.has_key('fix-only')
    specific_options['fix'] = not options.has_key('fix')

    if options.has_key("setup"):
        specific_options['controls'] = read_setup_file(options['setup'])
    
    if not specific_options['fix'] and specific_options['fix_only']:
        raise Usage('Cannot specify both --no-fix and --fix-only')

    return specific_options

class MacsFilenameParser(FilenameParser):
    def __init__(self, filename, controls = {}, debug=False, *args, **kwargs):
        if not os.path.splitext(filename)[1] == '.bam':
            raise InvalidFileException(filename)
        super(MacsFilenameParser, self).__init__(filename, *args, **kwargs)
        
        sample = self.protoname
        # check controls
        if controls.has_key(sample):
            run_name, control = controls[sample]
            if debug: print_debug(sample, 'has control', control)
            if control is None:
                self.control_file = None
            else:
                self.control_file = os.path.join(self.input_dir,
                                                 control + '.bam')
        elif sample in [v[1] for v in controls.values()]:
            if debug: print_debug(sample, 'is a control, aborting')
            raise InvalidFileException
        else:
            if debug: print_debug(sample, 'has no control indicated, continuing anyway')
            # not in setup.txt, make an entry in controls
            self.control_file = None
            run_name = sample
            controls[sample] = (sample, None)
            
        self.run_name = run_name  
        self.output_dir = os.path.join(self.output_dir, run_name)
        self.check_output_dir(self.output_dir)

def _prepend_cwd(d):
    '''prepends the current working directory to a directory d'''
    return os.path.join(os.getcwd(), d)

def _join_opt(a,b):
    '''joins an option with its value
    _join_opt('--option,'value') -> '--option=value'
    '''
    return '='.join([a,b])

def run_macs(parsed_filename, debug=False, silent=False,
           wig=True, single_wig=True, diag=True, subpeaks=True,
           make_pdf=True, fix_only=False, fix=True, 
           path_to_macs=None, path_to_R=None,
           **kwargs):
    """Run MACS on a BAM file and produce the pdf from the .R model"""
    if fix_only:
        return fix_subpeaks_filename(parsed_filename)

    input_file = parsed_filename.input_file
    control_file = parsed_filename.control_file
    if debug:
        print_debug('Processing', input_file)
        if control_file is not None:
            print_debug('with control', control_file)

    macs_options = ['--format=BAM']
    if wig:
        macs_options.append('--wig')
        if single_wig: macs_options.append('--single-wig')
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
    
    if debug: print_debug('Launching', ' '.join(step))
    job = subprocess.Popen(step,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT,
                           cwd=_prepend_cwd(parsed_filename.output_dir))

    stdout_data = job.communicate()[0]
    stdout_buffer = '{0!s}\n\n{1!s}'.format(' '.join(step), stdout_data)
   
    if subpeaks and fix:
        output = fix_subpeaks_filename(parsed_filename)
        if debug: print_debug('{0!s}\n{1!s}'.format(stdout_buffer, output))

    # now process R file
    if make_pdf:
        R_file = os.path.join(os.getcwd(), parsed_filename.output_dir,
                              ''.join([parsed_filename.protoname, '_model.r']))

        if debug: print_debug('Processing', R_file)

        if not os.path.exists(R_file):
            if not silent:
                print_debug('Warning:', R_file, 'does not exist')
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

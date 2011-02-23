#!/usr/bin/env python
'''
calls peaks using MACS

--no-wig            do not generate wiggle files (also disable subpeaks)
--chrom-wigs        generate a wig file per chromosome instead of a single wig
--no-diag           disable generation of diagnostic xls file
--no-subpeaks       do not call subpeaks with PeakSplitter
--no-pdf            do not generate a pdf of the model using R
--no-fix            do not auto-fix subpeaks filename
--fix-only          only fix subpeaks filename, do nothing else
'''
import sys
import os
import platform
import subprocess
import scripter
from scripter import path_to_executable, print_debug
VERSION = "2.4"

# HARDCODED CONTROLS #
CONTROLS = {
    'A549.GR.1uM.Dex.90m.samantha_cooper.081124_HWI-EAS355_0001_Meghan.s_1': 'A549.Input.1uM.Dex.90m.samantha_cooper.090326_HWI-EAS355_0006_Meghan_040109run.s_8',
    'Nalm6.GR.1uM.Dex.90m.samantha_cooper.080821_HWI-EAS355_0003_sam_meghan.s_2': 'Nalm6.Input.1uM.Dex.90m.samantha_cooper.080821_HWI-EAS355_0003_sam_meghan.s_1.all',
    'U2os-GRalpha.GR.1uM.Dex.90m.samantha_cooper.080821_HWI-EAS355_0003_sam_meghan.s_3': 'U2OS-GRalpha.Input.1uM.Dex.90m.samantha_cooper.090326_HWI-EAS355_0006_Meghan_040109run.s_7',
    'U2os-GRalpha.GR.1uM.Dex.90m.samantha_cooper.081124_HWI-EAS355_0001_Meghan.s_4': 'U2os-GRalpha.GR.0uM.EtOH.90m.samantha_cooper.081124_HWI-EAS355_0001_Meghan.s_2',
    'U2os-GRgamma.GR.1uM.Dex.90m.samantha_cooper.081124_HWI-EAS355_0001_Meghan.s_3': 'U2os-GRalpha.GR.0uM.EtOH.90m.samantha_cooper.081124_HWI-EAS355_0001_Meghan.s_2.all',
}
# HARDCODED CONTROLS #

def main():
    long_opts = ["no-wig", "no-diag", "chrom-wigs", "no-subpeaks",
                 "no-pdf", "fix-only", "no-fix"]
    e = scripter.Environment(long_opts=long_opts, doc=__doc__, version=VERSION)
    path_to_macs = path_to_executable(["macs14", "macs14.py", "macs"])
    if e.is_debug(): print_debug('Found MACS at', path_to_macs)
    path_to_R = path_to_executable(["R64", "R"])
    if e.is_debug(): print_debug('Found R at', path_to_R)
    e.set_source_dir('alignments.BAM')
    e.set_target_dir('fromBAM.macs')
    e.update_script_kwargs(check_script_options(e.get_options()))
    e.update_script_kwargs({'path_to_R': path_to_R,
                            'path_to_macs': path_to_macs})
    e.set_filename_parser(FilenameParser)
    e.do_action(run_macs)

def check_script_options(options):
    specific_options = {}

    specific_options['wig'] = not options.has_key('no-wig')
    specific_options['single_wig'] = not options.has_key('chrom-wigs')
    specific_options['diag'] = not options.has_key('no-diag')
    specific_options['subpeaks'] = not options.has_key('no-subpeaks')
    specific_options['make_pdf'] = not options.has_key('no-pdf')
    specific_options['fix_only'] = options.has_key('fix-only')
    specific_options['fix'] = not options.has_key('fix')
    
    if not specific_options['fix'] and specific_options['fix_only']:
        raise scripter.Usage('Cannot specify both --no-fix and --fix-only')

    return specific_options

class FilenameParser(scripter.FilenameParser):
    def __init__(self, filename, *args, **kwargs):
        if not os.path.splitext(filename)[1] == '.bam':
            raise scripter.InvalidFileException(filename)
        super(FilenameParser, self).__init__(filename, *args, **kwargs)
        self.output_dir = os.path.join(self.output_dir, self.protoname)
        self.check_output_dir(self.output_dir)

        self.has_control = False
        name, ext  = os.path.splitext(self.protoname)
        ext = ext.lstrip(os.extsep)

        if not ext.startswith('barcode') and not ext == 'all':
            name = self.protoname
            ext = None

        if CONTROLS.has_key(name) and (ext=='all' or ext is None):
            self.has_control = True
            self.is_control = False
            self.control_file = os.path.join(self.input_dir,
                                    os.extsep.join([CONTROLS[name],
                                                    selff.file_extension]))
        elif name in CONTROLS.values() or self.protoname in CONTROLS.values():
            self.has_control = False
            self.is_control = True
        else:
            self.has_control = False
            self.is_control = False

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

    if parsed_filename.is_control: return

    if debug:
        print_debug('Processing', parsed_filename.input_file)
        if parsed_filename.has_control:
            print_debug('with control', parsed_filename.control_file)

    macs_options = ['--format=BAM']
    if wig:
        macs_options.append('--wig')
        if single_wig: macs_options.append('--single-wig')
        if subpeaks: macs_options.append('--call-subpeaks')
    if diag: macs_options.append('--diag')
    macs_options.append(_join_opt('--name', parsed_filename.protoname))
    macs_options.append(_join_opt('--treatment',
                                  _prepend_cwd(parsed_filename.input_file)))
    if parsed_filename.has_control:
        macs_options.append(_join_opt('--control',
                                _prepend_cwd(parsed_filename.control_file)))

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

    (stdout_data, stderr_data) = job.communicate()
    stdout_buffer = '{!s}\n\n{!s}'.format(join(step), stdout_data)
   
    if subpeaks and fix:
        output = fix_subpeaks_file(parsed_filename)
        if debug: print_debug('{!s}\n{!s}'.format(stdout_buffer, output))

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
            (stdout_data, stderr_data) = job.communicate()
            R_pointer.close()
            stdout_buffer = '{!s}\n\n{!s}\n{!s}\n'.format(stdout_buffer,
                                                          ' '.join(step),
                                                          stdout_data)

    return stdout_buffer

def fix_subpeaks_file(parsed_filename, debug=False):
    """
    there's a few problems with the subpeaks file that we fix here
    
    1) Breaks MACS naming convention, rename it
    2) First line is not compatible with BED format, delete it
    """
    fn_parts = parsed_filename.protoname.split(os.extsep)
    path_to_subpeaks = os.path.join(parsed_filename.output_dir,
                        os.extsep.join([fn_parts[0], 'subpeaks'] + \
                                        fn_parts[1:-1] + [fn_parts[-1] +
                                                         '_peaks', 'bed']))
    new_path_to_subpeaks = os.path.join(parsed_filename.output_dir,
                            os.extsep.join([fn_parts[0]] + fn_parts[1:-1] + \
                                           [fn_parts[-1] + '_subpeaks', 'bed']))
    # check if the file exists
    if os.path.exists(new_path_to_subpeaks):
        return 'Cannot move {!s} to {!s}. Latter file already exists'.format(
                                        path_to_subpeaks, new_path_to_subpeaks)
    msg = ' '.join(['Renamed', path_to_subpeaks, 'to', new_path_to_subpeaks])
    try:
        old_file = open(path_to_subpeaks)
    except OSError:
        return 'Could not open {!s} for reading'.format(path_to_subpeaks)
    try:
        new_file = open(path_to_subpeaks, 'w')
    except OSError:
        return 'Could not open {!s} for writing'.format(new_path_to_subpeaks)
    try:
        discard_first_line = old_file.readline()
        new_file.writelines(old_file.readlines())
    except OSError:
        return 'Something went wrong copying {!s} to {!s}'.format(
                                        path_to_subpeaks, new_path_to_subpeaks)
        old_file.close()
        new_file.close()
    try:
        os.remove(path_to_subpeaks)
    except OSError:
        return 'Copied {!s} to {!s} but could not remove original'.format(path_to_subpeaks,
                                                                          new_path_to_subpeaks)
    return 'Moved {!s} to {!s} successfully'

if __name__=="__main__": main()

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
import subprocess
import time
import pkg_resources
import scripter
scripter.SCRIPT_DOC = __doc__
scripter.SCRIPT_VERSION = "2.4"
scripter.SOURCE_DIR = 'alignments.BAM'
scripter.TARGET_DIR = 'fromBAM.macs'
scripter.ALLOWED_EXTENSIONS = ['bam']
scripter.SCRIPT_LONG_OPTS = ["no-wig", "no-diag", "chrom-wigs", "subpeaks",
                             "no-pdf", "fix-only", "no-fix"]

PATH_TO_MACS = '/Library/Frameworks/Python.framework/Versions/2.7/bin/macs14'
PATH_TO_R = '/usr/bin/R64'

# HARDCODED CONTROLS #
CONTROLS = {
    'A549.GR.1uM.Dex.90m.samantha_cooper.081124_HWI-EAS355_0001_Meghan.s_1': 'A549.Input.1uM.Dex.90m.samantha_cooper.090326_HWI-EAS355_0006_Meghan_040109run.s_8',
    'Nalm6.GR.1uM.Dex.90m.samantha_cooper.080821_HWI-EAS355_0003_sam_meghan.s_2': 'Nalm6.Input.1uM.Dex.90m.samantha_cooper.080821_HWI-EAS355_0003_sam_meghan.s_1.all',
    'U2os-GRalpha.GR.1uM.Dex.90m.samantha_cooper.080821_HWI-EAS355_0003_sam_meghan.s_3': 'U2OS-GRalpha.Input.1uM.Dex.90m.samantha_cooper.090326_HWI-EAS355_0006_Meghan_040109run.s_7',
    'U2os-GRalpha.GR.1uM.Dex.90m.samantha_cooper.081124_HWI-EAS355_0001_Meghan.s_4': 'U2os-GRalpha.GR.0uM.EtOH.90m.samantha_cooper.081124_HWI-EAS355_0001_Meghan.s_2',
    'U2os-GRgamma.GR.1uM.Dex.90m.samantha_cooper.081124_HWI-EAS355_0001_Meghan.s_3': 'U2os-GRalpha.GR.0uM.EtOH.90m.samantha_cooper.081124_HWI-EAS355_0001_Meghan.s_2.all',
}
# HARDCODED CONTROLS #

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
        super(FilenameParser, self).__init__(filename, *args, **kwargs)
        self.output_dir = os.path.join(self.output_dir, self.protoname)
        if not self.is_dummy_file: self.check_output_dir(self.output_dir)

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
                                                    self.file_extension]))
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

def _print_debug(*args):
    statement = ' '.join(args)
    print >>sys.stderr, statement

def action(parsed_filename, debug=False, silent=False,
           wig=True, single_wig=True, diag=True, subpeaks=True,
           make_pdf=True, fix_only=False, fix=True,**kwargs):
    """Run MACS on a BAM file and produce the pdf from the .R model"""
    if fix_only:
        return fix_subpeaks_filename(parsed_filename)

    if parsed_filename.is_control: return

    if debug:
        _print_debug('Processing', parsed_filename.input_file)
        if parsed_filename.has_control:
            _print_debug('with control', parsed_filename.control_file)

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

    step = [PATH_TO_MACS] + macs_options

    if debug: _print_debug('Launching', ' '.join(step))
    job = subprocess.Popen(step,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT,
                           cwd=_prepend_cwd(parsed_filename.output_dir))

    (stdout_data, stderr_data) = job.communicate()
    stdout_buffer = os.linesep.join([' '.join(step), '', stdout_data])
   
    if subpeaks and fix:
        output = fix_subpeaks_filename(parsed_filename)
        if debug: _print_debug(os.linesep.join([stdout_buffer, output]))

    # now process R file
    if make_pdf:
        R_file = os.path.join(os.getcwd(), parsed_filename.output_dir,
                              ''.join([parsed_filename.protoname, '_model.r']))

        if debug: _print_debug('Processing', R_file)

        if not os.path.exists(R_file):
            if not silent:
                _print_debug('Warning:', R_file, 'does not exist')
        else:
            R_pointer = open(R_file)
            step =[PATH_TO_R, '--vanilla'] 
            job = subprocess.Popen(step,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT,
                                   stdin=R_pointer,
                                   cwd=_prepend_cwd(parsed_filename.output_dir))
            (stdout_data, stderr_data) = job.communicate()
            R_pointer.close()
            stdout_buffer = os.linesep.join([stdout_buffer, '', '',
                                            ' '.join(step), stdout_data, ''])

    return stdout_buffer

def fix_subpeaks_filename(parsed_filename, debug=False):
    fn_parts = parsed_filename.protoname.split(os.extsep)
    path_to_subpeaks = os.path.join(parsed_filename.output_dir,
                        os.extsep.join([fn_parts[0], 'subpeaks'] + \
                                        fn_parts[1:-1] + [fn_parts[-1] +
                                                         '_peaks', 'bed']))
    new_path_to_subpeaks = os.path.join(parsed_filename.output_dir,
                            os.extsep.join([fn_parts[0]] + fn_parts[1:-1] + \
                                           [fn_parts[-1] + '_subpeaks', 'bed']))
    msg = ' '.join(['Renamed', path_to_subpeaks, 'to', new_path_to_subpeaks])
    try: os.rename(path_to_subpeaks, new_path_to_subpeaks)
    except OSError:
        msg = ' '.join(['Could not move', path_to_subpeaks,
                        'to', new_path_to_subpeaks])
    return msg

if __name__=="__main__":
    scripter.check_script_options = check_script_options
    scripter.perform(action, FilenameParser)

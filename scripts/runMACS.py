#!/usr/bin/env python
'''
calls peaks using MACS
'''
import os
import subprocess
import time
import scripter
scripter.SCRIPT_DOC = __doc__
scripter.SCRIPT_VERSION = "2.0"
scripter.SOURCE_DIR = 'alignments.BAM'
scripter.TARGET_DIR = 'fromBAM.macs'

PATH_TO_MACS = '/Library/Frameworks/Python.framework/Versions/2.7/bin/macs14'
PATH_TO_R = '/usr/bin/R64'

# HARDCODED CONTROLS #
CONTROLS = {
    'A549.GR.1uM.Dex.90m.samantha_cooper.081124_HWI-EAS355_0001_Meghan.s_1.all': 'A549.Input.1uM.Dex.90m.samantha_cooper.090326_HWI-EAS355_0006_Meghan_040109run.s_8',
    'A549.GR.1uM.Dex.90m.samantha_cooper.081124_HWI-EAS355_0001_Meghan.s_1.barcodeGACG': 'A549.Input.1uM.Dex.90m.samantha_cooper.090326_HWI-EAS355_0006_Meghan_040109run.s_8',
    'A549.GR.1uM.Dex.90m.samantha_cooper.081124_HWI-EAS355_0001_Meghan.s_1.barcodeTCAT': 'A549.Input.1uM.Dex.90m.samantha_cooper.090326_HWI-EAS355_0006_Meghan_040109run.s_8',
    'Nalm6.GR.1uM.Dex.90m.samantha_cooper.080821_HWI-EAS355_0003_sam_meghan.s_2.all': 'Nalm6.Input.1uM.Dex.90m.samantha_cooper.080821_HWI-EAS355_0003_sam_meghan.s_1.all',
    'Nalm6.GR.1uM.Dex.90m.samantha_cooper.080821_HWI-EAS355_0003_sam_meghan.s_2.barcodeGACG': 'Nalm6.Input.1uM.Dex.90m.samantha_cooper.080821_HWI-EAS355_0003_sam_meghan.s_1.barcodeGACG',
    'Nalm6.GR.1uM.Dex.90m.samantha_cooper.080821_HWI-EAS355_0003_sam_meghan.s_2.TCAT': 'Nalm6.Input.1uM.Dex.90m.samantha_cooper.080821_HWI-EAS355_0003_sam_meghan.s_1.TCAT',
    'U2os-GRalpha.GR.1uM.Dex.90m.samantha_cooper.080821_HWI-EAS355_0003_sam_meghan.s_3': 'U2OS-GRalpha.Input.1uM.Dex.90m.samantha_cooper.090326_HWI-EAS355_0006_Meghan_040109run.s_7',
    'U2os-GRalpha.GR.1uM.Dex.90m.samantha_cooper.081124_HWI-EAS355_0001_Meghan.s_4': 'U2os-GRalpha.GR.0uM.EtOH.90m.samantha_cooper.081124_HWI-EAS355_0001_Meghan.s_2',
    'U2os-GRgamma.GR.1uM.Dex.90m.samantha_cooper.081124_HWI-EAS355_0001_Meghan.s_3.all': 'U2os-GRalpha.GR.0uM.EtOH.90m.samantha_cooper.081124_HWI-EAS355_0001_Meghan.s_2.all',
    'U2os-GRgamma.GR.1uM.Dex.90m.samantha_cooper.081124_HWI-EAS355_0001_Meghan.s_3.barcodeGACG': 'U2os-GRalpha.GR.0uM.EtOH.90m.samantha_cooper.081124_HWI-EAS355_0001_Meghan.s_2.barcodeGACG',
    'U2os-GRgamma.GR.1uM.Dex.90m.samantha_cooper.081124_HWI-EAS355_0001_Meghan.s_3.barcodeTCAT': 'U2os-GRalpha.GR.0uM.EtOH.90m.samantha_cooper.081124_HWI-EAS355_0001_Meghan.s_2.barcodeTCAT'
}
# HARDCODED CONTROLS #

class FilenameParser(scripter.FilenameParser):
    def __init__(self, filename, *args, **kwargs):
        super(FilenameParser, self).__init__(filename, *args, **kwargs)
        self.output_dir = os.path.join(self.output_dir, self.protoname)
        self.check_output_dir(self.output_dir)

        if CONTROLS.has_key(self.protoname):
            self.has_control = True
            self.is_control = False
            self.control_file = os.path.join(self.input_dir,
                                    os.extsep.join([CONTROLS[self.protoname],
                                                    self.file_extension]))
        elif self.protoname in CONTROLS.values():
            self.is_control = True
            self.has_control = False
        else:
            self.has_control = False
            self.is_control = False

def run_macs(parsed_filename, **kwargs):
    macs_options_1_3 = ['--tsize=30', '--format=BAM', '--bw=300', '--wig',
                        '--diag', '--pvalue=1e-5', '--mfold=32',
                        '--name='+parsed_filename.protoname,
                        '--treatment='+parsed_filename.input_file]
    macs_options = ['--format=BAM', '--call-subpeaks', '--wig', '--diag',
                    '--name='+parsed_filename.protoname, 
                    '--treatment='+parsed_filename.input_file]

    if parsed_filename.has_control:
        macs_options.extend('--control=' + parsed_filename.control_file)

    step = [PATH_TO_MACS] + macs_options
    job = subprocess.Popen(step,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           cwd=parsed_filename.output_dir)
    (stdout_data, stderr_data) = job.communicate()
    stdout_buffer = os.linesep.join([' '.join(step), '', stdout_data, ''])
    return stdout_buffer

def process_R_file(parsed_filename, **kwargs):
    R_file = os.path.join(parsed_filename.output_dir,
                          ''.join([parsed_filename.protoname, '_model.r']))
# wait for R file to exist
    while True:
        if os.path.exists(R_file): break
        else: time.sleep(10)

    R_pointer = open(R_file)
    step =[PATH_TO_R, '--vanilla'] 
    job = subprocess.Popen(step, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT,
                            stdin=R_pointer, cwd=parsed_filename.output_dir)
    (stdout_data, stderr_data) = job.communicate()
    R_pointer.close()
    stdout_buffer = os.linesep.join([' '.join(step), stdout_data, ''])
    return stdout_buffer

def action(parsed_filename, **kwargs):
    if parsed_filename.file_extension == 'bam' and not \
       parsed_filename.is_control:
        stdout1 = run_macs(parsed_filename, **kwargs)
        stdout2 = process_R_file(parsed_filename, **kwargs)
        return os.linesep.join([stdout1, stdout2])

if __name__=="__main__":
	scripter.perform(action, FilenameParser)

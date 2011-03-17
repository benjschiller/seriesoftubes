#!/usr/bin/env python
'''
Removes 3' linkers from FASTQ files

--linker=NNNNN  specifies linker

--no-target,--same-dir  keep new files in the same directory as the input files

'''
import os
import Bio.SeqIO.QualityIO
import scripter
from scripter import print_debug, assert_path, exit_on_Usage

VERSION = "2.4"
TARGET_DIR = 'fastq.no_linker'

@exit_on_Usage
def main():
    long_opts = ["no-target", "same-dir", "linker="] 
    
    e = scripter.Environment(long_opts=long_opts, doc=__doc__, version=VERSION)
    options = e.get_options()
    if options.has_key('same-dir') or options.has_key('no-target'):
        e.target_dir = os.curdir
    else:
        e.target_dir = TARGET_DIR
    e.set_filename_parser(BarcodeFilenameParser)
    require_linker(e)
    e.do_action(action)
    return

def require_linker(environment):
    """
    enforces that the environment has specified a linker
    and adds it to the script kwargs
    
    for now, that means it was included as a command-line option
    """
    options = environment.get_options()
    if not options.has_key('linker'): raise scripter.Usage('linker required')
    linker = options['linker']
    if not linker.isalpha() and len(linker)>0:
        raise scripter.Usage('invalid linker')
    environment.add_script_kwarg('linker', linker)
    return

def action(parsed_filename, linker='', verbose=False, **kwargs):

    if linker == '': raise scripter.Usage('no linker provided')
    f = open(parsed_filename.input_file, "r")
    record_generator = Bio.SeqIO.QualityIO.FastqGeneralIterator(f)

    i = 0
    linker_only = 0
    too_short = 0

    output_file = open(parsed_filename.output_filename, 'w')
    for title, seq, qual in record_generator:
        i += 1
        linker_index = seq.find(linker)
        
        if linker_index == -1:
            output_file.write("@%s\n%s\n+%s\n%s\n" % (title, seq, title, qual))
        
        if linker_index==0:
            linker_only += 1
            continue
        elif linker_index==1 or linker_index==2 or linker_index==3:
            too_short += 1
            continue
        else:
            output_file.write("@%s\n%s\n+%s\n%s\n" % (title, seq[0:linker_index],
                                                  title, qual[0:linker_index]))
    # close and exit #
    f.close()
    output_file.close()
    return '\n'.join([parsed_filename.output_filename,
                            'Processed ' + str(i) + ' records',
                            str(linker_only) + ' linker only dimers',
                            str(too_short) + ' sequences too short (1-3 bp)'])


class BarcodeFilenameParser(scripter.FilenameParser):
    def __init__(self, filename, verbose=False, target_dir=TARGET_DIR,
                 no_target = False, *args, **kwargs):
        super(BarcodeFilenameParser, self).__init__(filename,
                                                    target_dir=target_dir,
                                                    *args, **kwargs)
        # check for old-style
        if self.protoname.split('.')[-1]=='all':
            self.protoname = self.protoname[0:-4]
        if 'no_target': self.output_dir = self.input_dir

        # check if this is a paired-end file
        # if so, grab its partner
        illumina_name = os.path.basename(self.input_file)
        if illumina_name.count('_') >= 3:
            if verbose: print_debug('NOTICE: Detected paired read file.')
            iln_parts = illumina_name.split('_')
            self.paired_end = True
            if iln_parts[2] == '1':
                if verbose: print_debug('Attempting to find second file.')

                self.second_file = os.sep.join([self.input_dir,
                                                '_'.join(iln_parts[0:2] + ['2']
                                                         + iln_parts[3:])])
                self.protoname2 = os.path.splitext(
                                        os.path.basename(self.second_file))[0]
                if verbose: print_debug('Found', self.second_file)
                try: assert_path(self.second_file)
                except IOError:
                    if verbose: print_debug('Failed to find paired end file')
                    self.paired_end = False
            elif iln_parts[2] == '2':
                if verbose: print_debug('This is the second file, ignoring it.')
                self.is_invalid = True
            else:
                if verbose: print_debug('Failed to find paired end')
                self.paired_end = False
        else: self.paired_end = False

        self.output_filename = os.path.join(self.output_dir,
                            os.extsep.join([self.protoname,
                                            'no_linker',
                                             self.file_extension]))

if __name__== "__main__": main()
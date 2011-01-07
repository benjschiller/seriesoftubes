#!/usr/bin/env python
'''
Removes 3' linkers from FASTQ files

--linker=NNNNN  specifies linker

--no-target,--same-dir  keep new files in the same directory as the input files

'''
import os
import scripter
from scripter import print_debug, assert_path
import Bio.SeqIO
scripter.SCRIPT_DOC = __doc__
scripter.SCRIPT_VERSION = "2.3"
BOOLEAN_OPTS = ["pool-rc", "no-pool-rc", "no-target", "same-dir"]
scripter.SCRIPT_LONG_OPTS = ["linker="] + BOOLEAN_OPTS
#scripter.SOURCE_DIR = None
scripter.TARGET_DIR = os.curdir

def check_script_options(options):
    sopts = {}

    for option in BOOLEAN_OPTS:
        pyoption = "_".join(option.split("-"))
        sopts[pyoption] = options.has_key(option)

    if not options.has_key('linker'):
        raise scripter.Usage('linker required')
    else:
        linker = options['linker']
        if not linker.isalpha() and len(linker)>0:
            raise scripter.Usage('invalid linker')
        sopts['linker'] = linker

    if sopts['same_dir']: sopts['no_target'] = True

    return sopts

class IlluminaID:
    def __init__(self, name):
        name_parts = name.split(':')

        self.instrument = name_parts[0]
        self.title = name_parts[1]
        self.x = name_parts[2]

        last_part = name_parts[-1]
        p = last_part.find('#')
        s = last_part.find('/')

        self.y = last_part[0:p]
        self.barcode = last_part[p+1:s]
        self.member_of_pair = last_part[s+1:]

    def __str__(self):
        last_part = ''.join([self.y, '#', self.barcode,
                             '/', self.member_of_pair])
        return ':'.join([self.instrument, self.title, self.x, last_part])

def action(parsed_filename, linker='', verbose=False, **kwargs):
    stdout_buffer = ''

    if linker == '': raise scripter.Usage('no linker provided')
    f = open(parsed_filename.input_file, "rU")
    record_generator = Bio.SeqIO.parse(f, "fastq-illumina")

    linker_only = 0
    too_short = 0

    # and make a unmatched file
    output_file = open(parsed_filename.output_filename, 'w')
    # and if we're doing paired end, make a mismatched file
    for record in record_generator:
        linker_index = str(record.seq).find(linker)
        if linker_index is not -1:
            record = record[0:linker_index]
        
        if linker_index==0:
            linker_only += 1
            continue
        elif linker_index==1 or linker_index==2 or linker_index==3:
            too_short += 1
            continue
        else:
            Bio.SeqIO.write(record, output_file, "fastq-illumina")
    # close and exit #
    f.close()
    output_file.close()
    return os.linesep.join([parsed_filename.output_filename,
                            str(linker_only) + ' linker only dimers',
                            str(too_short) + ' sequences too short (1-3 bp)'])


def _complement(nucleotide):
    '''returns the complement of a single nucleotide'''
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                  'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    return complement.setdefault(nucleotide,'N')

def _rc(seq):
    '''returns the reverse complement of a sequence'''
    L = [_complement(s) for s in seq]
    L.reverse()
    return ''.join(L)

class BarcodeFilenameParser(scripter.FilenameParser):
    def __init__(self, filename, verbose=False, *args, **kwargs):
        super(BarcodeFilenameParser, self).__init__(filename,
                                                    *args, **kwargs)
        # check for old-style
        if self.protoname.split('.')[-1]=='all':
            self.protoname = self.protoname[0:-4]
        if kwargs['no_target']: self.output_dir = self.input_dir

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

if __name__== "__main__":
    scripter.check_script_options = check_script_options

    FilenameParser = BarcodeFilenameParser
    scripter.perform(action, FilenameParser)

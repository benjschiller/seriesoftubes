#!/usr/bin/env python
'''
Separates FASTQ files by barcodes

--barcodes=SEQ1,SEQ2    (default: TCAT, GACG, AGTC, CTGA)
--pool-rc               pool reverse complements
--no-pool-rc            do not pool reverse complements (default)
--no-target,--same-dir  keep new files in the same directory as the input files
'''
import os
import scripter
scripter.SCRIPT_DOC = __doc__
scripter.SCRIPT_VERSION = "2.1"
scripter.SCRIPT_LONG_OPTS = ["barcodes=", "pool-rc", "no-pool-rc",
                             "no-target", "same-dir"]
scripter.SOURCE_DIR = 'sequences.FASTQ'
scripter.TARGET_DIR = os.curdir

DEFAULT_BARCODES = ['TCAT', 'GACG', 'AGTC', 'CTGA']

def check_script_options(options):
    specific_options = {}
    if options.has_key('barcodes'):
        specific_options['barcodes'] = options['barcodes'].split(',')
    else: # default to TCAT,GACG
        specific_options['barcodes'] = DEFAULT_BARCODES

    if options.has_key('pool-rc') and options.has_key('no-pool-rc'):
        raise scripter.Usage('conflicting options: --pool-rc and',
                             '--no-pool-rc')
    elif options.has_key('pool-rc'):
        specific_options['pool_rc'] = True
    else:
        specific_options['pool_rc'] = False

    if options.has_key('no-target') or options.has_key('same-dir'):
        specific_options['no_target'] = True
    else:
        specific_options['no_target'] = False

    return specific_options

def action(parsed_filename, **kwargs):
    stdout_buffer = ''
    verbose = kwargs['verbose']
    debug = kwargs['debug']
    split_names = split_file(parsed_filename, **kwargs)
    if verbose:
        stdout_buffer = ''.join([stdout_buffer, 'Split ', filename,
                               'as ', ", ".join(fixed_names), os.linesep])
    return stdout_buffer

def get_barcode(x):
    return x.strip()[-6:-2]

def match_barcode(current_barcode, allowed_barcodes, mismatches=1,
                  pool_rc=False):
    accepted = []
    for barcode in allowed_barcodes:
        if pool_rc: valid_barcodes = [barcode, _rc(barcode)]
        else: valid_barcodes = [barcode]
        for barcode_ in valid_barcodes:
            barcode_length = len(barcode_)
            score = 0
            for i in range(barcode_length):
                if barcode_[i] == current_barcode[i]: score += 1
            if score >= (barcode_length - 1):
                accepted.append(barcode_)
                continue
    if len(accepted) == 1: return accepted[0]
    else: return ''

# read write four lines at a time #
def split_file(parsed_filename, pool_rc=False, barcodes=DEFAULT_BARCODES,
               **kwargs):
    f = open(parsed_filename.input_file)
    i = 1
    current_barcode = ''
    matched_barcode = ''
    files = {}
    filenames = []
    for barcode in barcodes:
        fname = parsed_filename.output_filename(barcode)
        filenames.append( fname )
        files[barcode] = open( fname , 'w' )
    for line in f:
        if i == 1:
            current_barcode = get_barcode(line)
            matched_barcode = match_barcode(current_barcode, barcodes) 
        if len(matched_barcode) is not 0: files[matched_barcode].write(line)
        i = i%4 + 1
    # close and exit #
    f.close()
    for f_ in files.values(): f_.close()
    return filenames

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
    def __init__(self, filename, *args, **kwargs):
        super(BarcodeFilenameParser, self).__init__(filename,
                                                    *args, **kwargs)
        if self.protoname.split('.')[-1]=='all':
            self.protoname = self.protoname[0:-4]
        if kwargs['no_target']: self.output_dir = self.input_dir

    def output_filename(self, barcode):
        return os.path.join(self.output_dir,
                            os.extsep.join([self.protoname,
                                            ''.join(['barcode', barcode]),
                                             self.file_extension]))

if __name__== "__main__":
    scripter.check_script_options = check_script_options
    
    FilenameParser = BarcodeFilenameParser
    scripter.perform(action, FilenameParser)

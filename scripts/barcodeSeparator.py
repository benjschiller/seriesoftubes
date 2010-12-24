#!/usr/bin/env python
'''
Separates FASTQ files by barcodes

--barcodes=SEQ1,SEQ2    (default: TCAT, GACG, AGTC, CTGA)
--pool-rc               pool reverse complements
--no-pool-rc            do not pool reverse complements (default)
--no-target,--same-dir  keep new files in the same directory as the input files
--strip-after-barcode=n strip n bases after the barcode is removed (5' end)
                        (by default this 1 now, and is ignored if GERALD
                         handled the barcoding)
'''
import os
import scripter
from scripter import print_debug, assert_path
import Bio.SeqIO
scripter.SCRIPT_DOC = __doc__
scripter.SCRIPT_VERSION = "2.3"
BOOLEAN_OPTS = ["pool-rc", "no-pool-rc", "no-target", "same-dir"]
scripter.SCRIPT_LONG_OPTS = ["barcodes=",
                             "strip-after-barcode="] + BOOLEAN_OPTS
#scripter.SOURCE_DIR = None
scripter.TARGET_DIR = os.curdir

DEFAULT_BARCODES = ['TCAT', 'GACG', 'AGTC', 'CTGA']

def check_script_options(options):
    sopts = {}
    if options.has_key('barcodes'):
        sopts['barcodes'] = options['barcodes'].split(',')
    else: # default to TCAT,GACG
        sopts['barcodes'] = DEFAULT_BARCODES

    for option in BOOLEAN_OPTS:
        pyoption = "_".join(option.split("-"))
        sopts[pyoption] = options.has_key(option)

    if options.has_key('pool-rc') and options.has_key('no-pool-rc'):
        raise scripter.Usage('conflicting options: --pool-rc and',
                             '--no-pool-rc')

    if sopts['same_dir']: sopts['no_target'] = True

    return sopts

def action(parsed_filename, **kwargs):
    stdout_buffer = ''
    verbose = kwargs['verbose']
    debug = kwargs['debug']
    split_names = split_file(parsed_filename, **kwargs)
    if verbose:
        stdout_buffer = ''.join([stdout_buffer, 'Split ',
                                 parsed_filename.input_file,
                                 'as ', ", ".join(fixed_names), os.linesep])
    return stdout_buffer

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
               strip_length = 1, **kwargs):
    f = open(parsed_filename.input_file, "rU")
    record_generator = Bio.SeqIO.parse(f, "fastq-illumina")
    if parsed_filename.paired_end:
        paired_end = True
        f2 = open(parsed_filename.second_file, "rU")
        record_generator2 = Bio.SeqIO.parse(f2, "fastq-illumina")
    else: paired_end = False

    matched_barcode = ''
    files = {}
    filenames = []
    for barcode in barcodes:
        fname = parsed_filename.output_filename(barcode)
        filenames.append(fname)
        files[barcode] = open(fname, 'w')
    # and make a unmatched file
    unmatched_filename = parsed_filename.output_filename("unmatched")
    filenames.append(unmatched_filename)
    unmatched_file = open(unmatched_filename, 'w')
    # and if we're doing paired end, make a mismatched file
    if paired_end:
        matched_barcode2 = ''
        files2 = {}
        filenames2 = []
        for barcode in barcodes:
            fname2 = parsed_filename.output_filename2(barcode)
            filenames2.append(fname2)
            files2[barcode] = open(fname2, 'w')

        unmatched_filename2 = parsed_filename.output_filename2("unmatched")
        filenames2.append(unmatched_filename2)
        unmatched_file2 = open(unmatched_filename2, 'w')

        mismatched_filename = parsed_filename.output_filename("mismatched")
        filenames.append(mismatched_filename)
        mismatched_file = open(mismatched_filename, 'w')

        mismatched_filename2 = parsed_filename.output_filename2("mismatched")
        filenames2.append(mismatched_filename2)
        mismatched_file2 = open(mismatched_filename2, 'w')
    for record in record_generator:
        # check to see if barcode is in name, as with old-style GERALD
        il_ID = IlluminaID(record.id)
        barcode = il_ID.barcode
        if barcode=='0':
            matched_barcode = match_barcode(str(record.seq), barcodes) 
            if len(matched_barcode) is not 0:
                record = record[len(matched_barcode) + strip_length:]
                il_ID.barcode = matched_barcode
                record.id = str(il_ID)
                record.name = str(il_ID)
                record.description = str(il_ID)
        else:
            matched_barcode = match_barcode(barcode, barcodes)

        # check if we have a matching paired read
        if paired_end:
            record2 = record_generator2.next()
            il_ID2 = IlluminaID(record2.id)
            barcode2 = il_ID2.barcode
            if barcode2=='0':
                matched_barcode2 = match_barcode(str(record2.seq), barcodes) 
                if len(matched_barcode2) is not 0:
                    record2 = record2[len(matched_barcode2) + strip_length:]
                    il_ID2.barcode = matched_barcode2
                    record2.id = str(il_ID2)
                    record2.name = str(il_ID2)
                    record2.description = str(il_ID2)
            else:
                matched_barcode2 = match_barcode(barcode2, barcodes)

            if len(matched_barcode) is not 0 and \
              matched_barcode == matched_barcode2:
                Bio.SeqIO.write(record,
                                files[matched_barcode], "fastq-illumina")
                Bio.SeqIO.write(record2,
                                files2[matched_barcode2], "fastq-illumina")
            elif len(matched_barcode) is 0 or len(matched_barcode2) is 0:
                Bio.SeqIO.write(record,
                                unmatched_file, "fastq-illumina")
                Bio.SeqIO.write(record2,
                                unmatched_file2, "fastq-illumina")
            else:
                Bio.SeqIO.write(record,
                                mismatched_file, "fastq-illumina")
                Bio.SeqIO.write(record2,
                                mismatched_file2, "fastq-illumina")
        else:
            if len(matched_barcode) is not 0:
                Bio.SeqIO.write(record,
                                files[matched_barcode], "fastq-illumina")
            else:
                Bio.SeqIO.write(record,
                                unmatched_file, "fastq-illumina")
    # close and exit #
    f.close()
    for f_ in files.values(): f_.close()
    unmatched_file.close()
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

    def output_filename(self, barcode):
        return os.path.join(self.output_dir,
                            os.extsep.join([self.protoname,
                                            '_'.join(['barcode', barcode]),
                                             self.file_extension]))
    def output_filename2(self, barcode):
        return os.path.join(self.output_dir,
                            os.extsep.join([self.protoname2,
                                            '_'.join(['barcode', barcode]),
                                             self.file_extension]))

if __name__== "__main__":
    scripter.check_script_options = check_script_options
    
    FilenameParser = BarcodeFilenameParser
    scripter.perform(action, FilenameParser)

#!/usr/bin/env python
'''
Separates FASTQ files by barcodes

--barcodes=SEQ1,SEQ2    (default: TCAT, GACG, AGTC, CTGA)
--no-target,--same-dir  keep new files in the same directory as the input files
--strip-after-barcode=n strip n bases after the barcode is removed (5' end)
                        (by default this 1 now, and is ignored if GERALD
                         handled the barcoding)
'''
import os
import operator
import scripter
from scripter import print_debug, assert_path, InvalidFileException
import Bio.SeqIO.QualityIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
VERSION = "2.4"
DEFAULT_BARCODES = ['TCAT', 'GACG', 'AGTC', 'CTGA']

def main():
    boolean_opts = ["no-target", "same-dir"]
    long_opts = ["barcodes=", "strip-after-barcode="] + boolean_opts
    e = scripter.Environment(long_opts=long_opts, version=VERSION, doc=__doc__)
    e.set_target_dir(os.curdir)
    e.parse_boolean_opts(boolean_opts)
    e.update_script_kwargs(check_script_options(e.get_options()))
    e.set_filename_parser(BarcodeFilenameParser)
    e.do_action(split_file)

def check_script_options(options):
    sopts = {}
    if options.has_key('barcodes'):
        sopts['barcodes'] = options['barcodes'].split(',')
    else: # default to TCAT,GACG
        sopts['barcodes'] = DEFAULT_BARCODES

    if options.has_key('same-dir') or options.has_key('no-target'):
        sopts['no_target'] = True

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

def match_barcode(current_barcode, allowed_barcodes, mismatches=1):
    accepted = []
    for barcode in allowed_barcodes:
        barcode_length = len(barcode)
        m = 0
        i = 0
        while True:
            i += 1
            if m > mismatches: break
            if i == barcode_length:
                accepted.append(barcode)
                break
            if not barcode[i] == current_barcode[i]: m+=1
#        score = sum(map(operator.eq, current_barcode, barcode))
#        if not score + mismatches < barcode_length:
#            accepted.append(barcode)
        if len(accepted) == 1: return accepted[0]
    else: return ''

# read write four lines at a time #
def split_file(parsed_filename, barcodes=DEFAULT_BARCODES,
               strip_length = 1, **kwargs):
    f = open(parsed_filename.input_file, "rU")
    record_generator = FastqGeneralIterator(f)
    if parsed_filename.paired_end:
        paired_end = True
        f2 = open(parsed_filename.second_file, "rU")
        record_generator2 = FastqGeneralIterator(f2)
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
    for title, seq, qual in record_generator:
        # check to see if barcode is in name, as with old-style GERALD
        il_ID = IlluminaID(title)
        barcode = il_ID.barcode
        if barcode=='0':
            matched_barcode = match_barcode(seq, barcodes) 
            if not len(matched_barcode) == 0:
                seq = seq[len(matched_barcode) + strip_length:]
                qual = qual[len(matched_barcode) + strip_length:]
                il_ID.barcode = matched_barcode
                title = str(il_ID)
        else: matched_barcode = match_barcode(barcode, barcodes)

        # check if we have a matching paired read
        if paired_end:
            title2, seq2, qual2 = record_generator2.next()
            il_ID2 = IlluminaID(title2)
            barcode2 = il_ID2.barcode
            if barcode2=='0':
                matched_barcode2 = match_barcode(seq2, barcodes) 
                if not len(matched_barcode2) == 0:
                    seq2 = seq2[len(matched_barcode) + strip_length:]
                    qual2 = qual2[len(matched_barcode) + strip_length:]
                    il_ID2.barcode = matched_barcode2
                    title2 = str(il_ID)
            else: matched_barcode2 = match_barcode(barcode2, barcodes)

            if len(matched_barcode) is not 0 and \
              matched_barcode == matched_barcode2:
                files[matched_barcode].write("@%s\n%s\n+%s\n%s\n" % 
                                             (title, seq, title, qual))
                files2[matched_barcode2].write("@%s\n%s\n+%s\n%s\n" % 
                                             (title2, seq2, title2, qual2))
            elif len(matched_barcode) == 0 or len(matched_barcode2) == 0:
                unmatched_file.write("@%s\n%s\n+%s\n%s\n" % 
                                             (title, seq, title, qual))
                unmatched_file2.write("@%s\n%s\n+%s\n%s\n" % 
                                             (title2, seq2, title2, qual2))
            else:
                mismatched_file.write("@%s\n%s\n+%s\n%s\n" % 
                                             (title, seq, title, qual))
                mismatched_file2.write("@%s\n%s\n+%s\n%s\n" % 
                                             (title2, seq2, title2, qual2))
        else:
            if len(matched_barcode) is not 0:
                files[matched_barcode].write("@%s\n%s\n+%s\n%s\n" % 
                                             (title, seq, title, qual))
            else:
                unmatched_file.write("@%s\n%s\n+%s\n%s\n" % 
                                             (title, seq, title, qual))
    # close and exit #
    f.close()
    for f_ in files.values(): f_.close()
    unmatched_file.close()
    
    if verbose:
        return 'Split {!s} as {!s}'.format(parsed_filename.input_file,
                                           join(filenames))
    else: return

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
                raise InvalidFileException(self.input_file)
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

if __name__== "__main__": main()
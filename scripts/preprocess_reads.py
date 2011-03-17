#!/usr/bin/env python
"""
preprocess_reads.py is meant to be run on GERALD output
                    although it will run on any FASTQ files*

We can
+ Separate reads by a set of variable-length barcodes
+ Cleave linker/adaptor sequence from the 3' ends of reads
+ Cleave adaptor sequence from the 5' end
    - before barcode
    + after barcode
? Check for positional bias / reading errors
    - before removing adaptors/barcodes (sequencer bias)
    - after removing adaptors/barcodes (library bias)
+ Removing trailing Ns from sequences
+ Discard sequences that are less than 4 nucleotides in length
- Produce FASTQ sequence files ready for immediate alignment

Command-line flags
NOT IMPLEMENTED YET:
--config=foo            Use configuration in file foo
--save-config=foo       Save configuration for this run to file foo 
IMPLEMENTED:
--barcodes=SEQ1,SEQ2    (default: TCAT,GACG,AGTC,CTGA)
--no-target,--same-dir  keep new files in the same directory as the input files
--strip-before-barcode=n strip n bases after the barcode is removed (5' end)
                        (by default this 1 now, and is ignored if GERALD
                         handled the barcoding)
--strip-after-barcode=n strip n bases after the barcode is removed (5' end)
                        (by default this 1 now, and is ignored if GERALD
                         handled the barcoding)

*it expects that files are named s_?_sequence.* (single-end reads) or
                                 s_?_[12]_sequence.* (paired-end reads)
"""
import os
from itertools import imap, izip, chain
from functools import partial
import scripter
from scripter import print_debug, assert_path, InvalidFileException
from Bio.SeqIO.QualityIO import FastqGeneralIterator
VERSION = "2.4"
DEFAULT_BARCODES = ['TCAT', 'GACG', 'AGTC', 'CTGA']

def main():
    boolean_opts = ["no-target", "same-dir"]
    long_opts = ["barcodes=", "strip-after-barcode=", "linker=",
                 "strip-before-barcode="] + boolean_opts
    e = scripter.Environment(long_opts=long_opts, version=VERSION, doc=__doc__)
    e.set_target_dir(os.curdir)
    e.parse_boolean_opts(boolean_opts)
    e.update_script_kwargs(check_script_options(e.get_options()))
    e.set_filename_parser(BarcodeFilenameParser)
    e.do_action(splitter)

def check_script_options(options):
    sopts = {}
    if not options.has_key('linker'):
        linker = None
    else:
        linker = options['linker']
        if not linker.isalpha() and len(linker)>0:
            raise scripter.Usage('invalid linker')
    sopts['linker'] = linker
    if options.has_key('barcodes'):
        sopts['barcodes'] = options['barcodes'].split(',')
    else: # default to TCAT,GACG
        sopts['barcodes'] = DEFAULT_BARCODES

    if options.has_key('strip-after-barcode'):
        sopts['strip_after_barcode'] = options['strip-after-barcode']
    else: # default to TCAT,GACG
        sopts['strip_after_barcode'] = 1
        
    if options.has_key('strip-before-barcode'):
        sopts['strip_before-barcode'] = options['strip-before-barcode']
    else: # default to TCAT,GACG
        sopts['strip_before_barcode'] = 0
        
    if options.has_key('same-dir') or options.has_key('no-target'):
        sopts['no_target'] = True

    return sopts

#class IlluminaID:
#    def __init__(self, name):
#        name_parts = name.split(':')
#
#        self.instrument = name_parts[0]
#        self.title = name_parts[1]
#        self.x = name_parts[2]
#
#        last_part = name_parts[-1]
#        p = last_part.find('#')
#        s = last_part.find('/')
#
#        self.y = last_part[0:p]
#        self.barcode = last_part[p+1:s]
#        self.member_of_pair = last_part[s+1:]
#
#    def __str__(self):
#        last_part = ''.join([self.y, '#', self.barcode,
#                             '/', self.member_of_pair])
#        return ':'.join([self.instrument, self.title, self.x, last_part])

def match_barcode(seq, barcodes, mismatches=1):
    """
    try to match seq to a list barcodes
    allow mismatches (default 1)
    returns the match (not seq) or None
    """
    accepted = None
    ne = str.__ne__
    barcode_lengths = map(len, barcodes)
    max_barcode_length = max(barcode_lengths)
    if max_barcode_length > len(seq):
        print 'slow', len(seq)
        # slow implementation
        for barcode in barcodes:
            barcode_length = len(barcode)
            comparable = seq[0:barcode_length]
            if not len(comparable) == barcode_length: continue
            hamming = sum(imap(ne, barcode, seq[0:barcode_length]))
            if mismatches > hamming:
                if accepted is None: accepted = barcode
                else: return None
        return accepted
    else:
        accepted = None
        for barcode in barcodes:
            if mismatches > sum(imap(ne, barcode, seq)):
                if accepted is None: accepted = barcode
                else: return None
        return accepted
#        matching = lambda x: mismatches > sum(imap(ne, x, seq))
#        good = map(matching, barcodes)
#        if good.count(True) == 1: return barcodes[good.index(True)]
#        else: return None

def pretrim_record_5prime(record, trim_length=0):
    '''
    trim an assigned record from the 5' end
    expects (barcode, (title, seq, qual))
    '''
    if trim_length == 0: return record
    title, seq, qual = record
    return (title, seq[trim_length:], qual[trim_length:])

def trim_record_5prime(record, trim_length=0):
    '''
    trim an assigned record from the 5' end
    expects (barcode, (title, seq, qual))
    '''
    barcode, (title, seq, qual) = record
    return (barcode, (title, seq[trim_length:], qual[trim_length:]))
    
def trim_record_3prime(record, trim_length=0):
    '''
    trim a record from the 3' end
    expects (barcode, (title, seq, qual))
    '''
    barcode, (title, seq, qual) = record
    return (barcode, (title, seq[0:-trim_length], qual[0:-trim_length]))

def trim_trailing_Ns(record):
    '''
    trim any trailing 3' 'N's
    expects record is (barcode, (title, seq, qual))
    returns truncated (barcode, (title, seq, qual))
    '''
    barcode, (title, seq, qual) = record
    trimmed_seq = seq.rstrip('N')
    i = len(trimmed_seq)
    return (barcode, (title, seq[0:i], qual[0:i]))

def splitter(pf, **kwargs):
    if pf.paired_end: return split_paired_files(pf, **kwargs)
    else: return split_file(pf, **kwargs)

def split_file(parsed_filename,
               barcodes=DEFAULT_BARCODES, linker=None, min_length=4,
               verbose=False, **kwargs):
    f = open(parsed_filename.input_file, "rU")
    records = FastqGeneralIterator(f)
    
    barcoded_files = {}
    filenames = []
    output_filename = parsed_filename.output_filename
    for barcode in barcodes:
        fname = output_filename(barcode)
        filenames.append(fname)
        barcoded_files[barcode] = open(fname, 'w')
        
    # and make a unmatched file
    unmatched_filename = output_filename("unmatched")
    filenames.append(unmatched_filename)
    unmatched_file = open(unmatched_filename, 'w')
    
    # assigned_records is a list of (barcode, record) items
    writer = partial(write_record, barcoded_files=barcoded_files,
                     unmatched_file=unmatched_file, linker=linker,
                     min_length=min_length)
    final_records = apply_plan(records, barcodes=barcodes, linker=linker,
                               **kwargs)
    results = map(writer, final_records)
    linker_only = results.count('linker')
    too_short = results.count('short')
    record_count = len(results)
    
    # close and exit #
    f.close()
    for f_ in barcoded_files.values(): f_.close()
    unmatched_file.close()
    
    if verbose:
        stdout_buffer = 'Split {!s} as {!s}'.format(parsed_filename.input_file,
                                           ', '.join(filenames))
        stdout_buffer += '\n'.join([parsed_filename.output_filename,
                    'Processed {!s} records'.format(record_count),
                    '{!s} linker only dimers'.format(linker_only),
                    '{!s} sequences too short (1-3 bp)'.format(too_short)])
        return stdout_buffer
    else: return
    
def assign_record(record, barcodes=DEFAULT_BARCODES):
    """
    Assign a record to a barcode
    returns (barcode, new_record)
    if unmatched, returns (None, record)
    """
    title, seq, qual = record
    title_head, last_part = title.rsplit(':', 1)
    pound_loc = last_part.find('#')
    slash_loc = pound_loc + last_part[pound_loc:].find('/')
    barcode = last_part[pound_loc+1:slash_loc]
    if barcode=='0':
        assigned_barcode = match_barcode(seq, barcodes) 
        if assigned_barcode is not None:
            barcode_len = len(assigned_barcode)
            seq = seq[barcode_len:]
            qual = qual[barcode_len:]
            barcode = assigned_barcode
            last_part = ''.join([last_part[0:pound_loc+1], barcode,
                                 last_part[slash_loc:]])
            title = ':'.join((title_head, last_part))
            record = (title, seq, qual)
    elif barcode.isdigit():
        # then we have a numbered index from Illumina, just use it as-is
        assigned_barcode = barcode
    elif barcode.isalpha():
        # then we already extracted the barcode at some point, try to match it
        assigned_barcode = match_barcode(barcode, barcodes)
    else:
        # couldn't do anything, give up
        assigned_barcode = None
    return (assigned_barcode, record)

def write_record(assigned_record, barcoded_files=None, unmatched_file=None,
                 linker=None, min_length=4):
    """
    write the assigned_record to the correct barcode file
    an assigned_record is a (barcode, (title, seq, qual)) (all strs)
    
    will not work unless you provide a dictionary of barcoded_files and an
    unmatched_file object
    """
    barcode, (title, seq, qual) = assigned_record
    # check if record 1 is valid
    if len(seq) == 0 and barcode == linker: problem = 'linker'
    elif len(seq) < min_length: problem = 'short'
    else: problem = None
    # produce lines
    line_fmt = "@{!s}\n{!s}\n+{!s}\n{!s}\n"
    line = line_fmt.format(title, seq, title, qual)
    if problem is None:
        if barcode is None: unmatched_file.write(line)
        else: barcoded_files[barcode].write(line)
    return problem
        
def write_record_pair(assigned_record_pair,
                       barcoded_file_pairs=None,
                       unmatched_files=None, mismatched_files=None,
                       linker=None, min_length=4):
    """
    write the assigned_record to the correct barcode file
    an assigned_record is a (barcode, (title, seq, qual)) (all strs)
    
    will not work unless you provide a dictionary of barcoded_file_pairss and
    unmatched_files and mismatched_files (pairs of file object)
    """
    barcode, (title, seq, qual) = assigned_record_pair[0]
    barcode2, (title2, seq2, qual2) = assigned_record_pair[1]
    # check if record 1 is valid
    if len(seq) == 0 and barcode == linker: problem = 'linker'
    elif len(seq) < min_length: problem = 'short'
    else: problem = None
    # check if record 2 is valid
    if len(seq2) == 0 and barcode2 == linker: problem2 = 'linker'
    elif len(seq2) < min_length: problem2 = 'short'
    else: problem2 = None
    # produce lines
    if problem is not None and problem2 is not None:
        return (problem, problem2)
    line_fmt = "@{!s}\n{!s}\n+{!s}\n{!s}\n"
    line1 = line_fmt.format(title, seq, title, qual)
    line2 = line_fmt.format(title2, seq2, title2, qual2)
    if barcode == barcode2:
        if problem is None and problem2 is None:
            if barcode is None:
                unmatched_files[0].write(line1)
                unmatched_files[1].write(line2)
            else:
                barcoded_file_pairs[barcode][0].write(line1)
                barcoded_file_pairs[barcode][1].write(line2)
        elif problem is None:
            mismatched_files[0].write(line1)
        else: #problem2 is None
            mismatched_files[1].write(line2)
    else:
        if problem is None: mismatched_files[0].write(line1)
        if problem2 is None: mismatched_files[1].write(line2)
    return (problem, problem2)
        
def is_empty_record(record):
#    barcode, (title, seq, qual) = record
    return len(record[1][1]) == 0

def pair_has_empty_record(record_pair):
    return is_empty_record(record_pair[0]) or is_empty_record(record_pair[1])
        
def apply_plan(records, barcodes=DEFAULT_BARCODES,
               min_length=4,
               strip_after_barcode = 1, strip_before_barcode=0,
               **kwargs):
    assigner = partial(assign_record, barcodes=barcodes)
    if strip_before_barcode > 0:
        pretrimmer = partial(pretrim_record_5prime,
                          trim_length=strip_before_barcode)
        pretrimmed_records = imap(pretrimmer, records)
        assigned_records = imap(assigner, pretrimmed_records)
    else:
        assigned_records = imap(assigner, records)
    if strip_after_barcode > 0:
        trimmer = partial(trim_record_5prime, trim_length=strip_after_barcode)
        trimmed_records = imap(trimmer, assigned_records)
        further_trimmed_records = imap(trim_trailing_Ns, trimmed_records)
        return further_trimmed_records
    else:
        further_trimmed_records = imap(trim_trailing_Ns, assigned_records)
        return further_trimmed_records

def split_paired_files(parsed_filename,
                       barcodes=DEFAULT_BARCODES, linker=None,
                min_length=4,
                verbose=False, **kwargs):
    f = open(parsed_filename.input_file, "rU")
    f2 = open(parsed_filename.second_file, "rU")
    records = FastqGeneralIterator(f)
    records2 = FastqGeneralIterator(f2)
        
    barcoded_file_pairs = {}
    filenames = []
    output_filename = parsed_filename.output_filename
    output_filename2 = parsed_filename.output_filename2
    for barcode in barcodes:
        fname = output_filename(barcode)
        fname2 = output_filename2(barcode)
        filenames.extend((fname, fname2))
        barcoded_file_pairs[barcode] = (open(fname, 'w'), open(fname2, 'w'))
        
    # and make a unmatched file
    unmatched_filename = output_filename("unmatched")
    unmatched_filename2 = output_filename2("unmatched")
    unmatched_files = (open(unmatched_filename, 'w'),
                       open(unmatched_filename, 'w'))

    mismatched_filename = output_filename("mismatched")
    mismatched_filename2 = output_filename2("mismatched")
    mismatched_files = (open(mismatched_filename, 'w'),
                        open(mismatched_filename, 'w'))
    filenames.extend((unmatched_filename, unmatched_filename2,
                      mismatched_filename, mismatched_filename2))

    writer = partial(write_record_pair,
                     barcoded_file_pairs=barcoded_file_pairs,
                     unmatched_files=unmatched_files,
                     mismatched_files=mismatched_files,
                     linker=linker, min_length=min_length)
    final_records = apply_plan(records, barcodes=barcodes, linker=linker,
                               **kwargs)
    final_records2 = apply_plan(records2, barcodes=barcodes, linker=linker,
                                **kwargs)
    final_record_pairs = izip(final_records, final_records2)
    results = list(chain.from_iterable(map(writer, final_record_pairs)))
    linker_only = results.count('linker')
    too_short = results.count('short')
    record_count = len(results)

    # close and exit #
    f.close()
    for f_, f2_ in barcoded_file_pairs.values():
        f_.close()
        f2_.close()
    unmatched_files[0].close()
    unmatched_files[1].close()
    mismatched_files[0].close()
    mismatched_files[1].close()
    
    if verbose:
        stdout_buffer = 'Split {!s}, {!s} as {!s}'.format(parsed_filename.input_file,
                                                parsed_filename.second_file,
                                           ', '.join(filenames))
        stdout_buffer += '\n'.join([parsed_filename.output_filename,
                            'Processed {!s} records'.format(record_count),
                            '{!s} linker only dimers'.format(linker_only),
                            '{!s} sequences too short (1-3 bp)'.format(too_short)])
        return stdout_buffer
    else: return


def cleave_3prime_linker(record, linker='', verbose=False, **kwargs):
    """
    cleave_linker takes a record in the format
    (barcode, (title, seq, qual))
    and cleaves the linker sequence from the 3' end, if it is present
    """
    if linker == '': raise scripter.Usage('no linker provided')
    barcode, (title, seq, qual) = record

    linker_index = seq.find(linker)
    
    if linker_index == -1: return record
    elif linker_index==0: return (linker, (title, seq, qual))
    else:
        return (barcode, (title, seq[0:linker_index], qual[0:linker_index]))

#def reduce_counts():
#    if debug:
#        stdout_buffer = 'Determining basewise error rates in {!s}\n'.format(
#                                                    parsed_filename.input_file)
#    else: stdout_buffer = ''
#    counter = lambda x: (x == 'A', x =='T', x == 'G', x == 'C', x == 'N')
#    seqs = imap(str.upper, imap(itemgetter(1), records))
#    reducto = reduce(add, izip(imap(counter, seqs)))
#            
#    f = open(output_filename, 'w')
#    first_row = ("cycle", "pyindex", "A", "T", "G", "C", "N")
#    f.write('\t'.join(["{!s}"]*6).format(*first_row))
#    counts = list(counts)
#    for i in range(count_len):
#        row = (i+2, i, counts[i][0], counts[i][1], counts[i][2], counts[i][3],
#               counts[i][4])
#        f.write('\t'.join(["{!s}"]*6).format(*row))
#    f.close()
#    return stdout_buffer


class BarcodeFilenameParser(scripter.FilenameParser):
    def __init__(self, filename, verbose=False, *args, **kwargs):
        super(BarcodeFilenameParser, self).__init__(filename,
                                                    *args, **kwargs)
        protoname = self.protoname
        # check for old-style
        if os.path.splitext(protoname)[-3:] == 'all':
            protoname = protoname[0:-4]
        if kwargs['no_target']: self.output_dir = self.input_dir
        
        # check if this is a paired-end file
        # if so, grab its partner
        input_file = self.input_file
        illumina_name = os.path.basename(input_file)
        if illumina_name.count('_') >= 3:
            if verbose: print_debug('NOTICE: Detected paired read file.')
            iln_parts = illumina_name.split('_')
            if iln_parts[2] == '1':
                if verbose: print_debug('Attempting to find second file.')

                second_file = os.sep.join([self.input_dir,
                                           '_'.join(iln_parts[0:2] + ['2'] 
                                                    + iln_parts[3:])])
                self.protoname2 = os.path.splitext(
                                        os.path.basename(second_file))[0]
                try:
                    assert_path(second_file)
                    if verbose: print_debug('Found', second_file)
                    self.second_file = second_file
                    paired_end = True
                except IOError:
                    if verbose: print_debug('Failed to find paired end file')
                    paired_end = False
            elif iln_parts[2] == '2':
                if verbose: print_debug('This is the second file, ignoring it.')
                raise InvalidFileException(input_file)
            else:
                if verbose: print_debug('Failed to find paired end')
                paired_end = False
        else: paired_end = False
        self.paired_end = paired_end

    def output_filename(self, barcode):
        return os.path.join(self.output_dir,
                            os.extsep.join([self.protoname,
                                            'barcode_' + barcode,
                                            self.file_extension]))
    def output_filename2(self, barcode):
        return os.path.join(self.output_dir,
                            os.extsep.join([self.protoname2,
                                            'barcode_' + barcode,
                                            self.file_extension]))

if __name__== "__main__": main()
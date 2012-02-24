#!/usr/bin/env python
"""
preprocess_reads.py is meant to be run on GERALD output
                    although it will run on any FASTQ files*

We can
+ Separate reads by a set of variable-length barcodes
+ Cleave linker/adaptor sequence from the 3' ends of reads
+ Cleave adaptor sequence from the 5' end
    + before barcode
    + after barcode
+ Removing trailing Ns from sequences
+ Discard sequences that are less than 4 nucleotides in length
+ Produce gzipped FASTQ sequence files ready for immediate alignment

A configuration file 'preprocess_reads.cfg' is saved in target
directory (unless one is provided by the user).

*it expects that files are named s_?_sequence.* (single-end reads) or
                                 s_?_[12]_sequence.* (paired-end reads)
"""
import os
from itertools import imap, izip, chain
from functools import partial
import scripter
from seriesoftubes.converters.discover import discover_file_format
from seriesoftubes.fnparsers import BarcodeFilenameParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from ConfigParser import SafeConfigParser
from errno import ENOENT, EACCES
from os import access, strerror, R_OK
from os.path import exists
import argparse
from gzip import GzipFile
from pkg_resources import get_distribution
__version__ = get_distribution('seriesoftubes').version
VERSION = __version__

def valid_seq(s):
        if not s.isalpha() and len(s)>0:
            msg = '%s is not a valid sequence'
            raise argparse.ArgumentTypeError(msg)
        else:
            return s.upper()
        
def main():
    e = scripter.Environment(version=VERSION, doc=__doc__)
    parser = e.argument_parser
    parser.add_argument('--strip-after-barcode', default=1, type=int,
                        help="""strip n bases after the barcode is removed (5' end)
(by default this 1 now, and is ignored if GERALD handled the barcoding)""")
    parser.add_argument('--strip-before-barcode', default=0, type=int,
                        help="""strip n bases before the barcode is removed (5' end)
(by default this 0 now, and is ignored if GERALD handled the barcoding)""")
    parser.add_argument('--max-length', type=int,
                        help='truncate final sequences to n bases (default: ignore)')
    bgroup = parser.add_argument_group('barcodes',
                                       'Specify sequence barcodes in the sample(s)')
    bgroup.add_argument('-b', '--barcodes', action='append', type=valid_seq,
                        help="Specify a barcode sequence. May be invoked multiple times")
    bgroup.add_argument('--kry-barcodes', dest='barcodes', action='store_const',
                        help='Alias for -bTCAT -bGACG -bAGTC -bCTGA',
                        const=['TCAT', 'GACG', 'AGTC', 'CTGA'])
    parser.add_argument('--linker', default=None, type=valid_seq,
                        help="Specify a 3' adaptor/linker sequence that we should clip off of each read")
    parser.set_defaults(**{'target': 'processed'})
    e.set_filename_parser(BarcodeFilenameParser)
    e.set_config_reader(read_config)
    e.set_config_writer(write_config)
    e.do_action(splitter)

def write_config(barcodes=None, max_length=None, strip_after_barcode=0,
                 strip_before_barcode=0, linker=None,
                 *args, **kwargs):
    config = SafeConfigParser()
    config.add_section('main')
    if barcodes is not None:
        config.set('main', 'barcodes', ','.join(barcodes))
    if linker is not None:
        config.set('main', 'linker', linker)
    if max_length is not None:
        config.set('main', 'max-length', str(max_length))
    config.set('main', 'strip-after-barcode', str(strip_after_barcode))
    config.set('main', 'strip-before-barcode', str(strip_before_barcode))
    with open('preprocess_reads.cfg', 'wb') as configfile:
        config.write(configfile)

def read_config(setup_file):
    if not exists(setup_file):
        raise IOError(ENOENT, strerror(ENOENT), setup_file)
    if not access(setup_file, R_OK):
        raise IOError(EACCES, strerror(EACCES), setup_file)
    config = SafeConfigParser()
    config.readfp(open(setup_file, 'rU'))
    context = {}
    if config.has_option('main', 'barcodes'):
        context['barcodes'] = config.get('main', 'barcodes').split(',')
    else: context['barcodes'] = None
    if config.has_option('main', 'max-length'):
        context['max_length'] = int(config.get('main', 'max-length'))
    else: context['max_length'] = None
    if config.has_option('main', 'linker'):
        context['linker'] = config.get('main', 'linker')
    else: context['linker'] = None
    context['strip_after_barcode'] = int(config.get('main', 'strip-after-barcode'))
    context['strip_before_barcode'] = int(config.get('main', 'strip-before-barcode'))
    return context

def match_barcode(seq, barcodes, mismatches=1):
    """
    try to match seq to a list of barcodes
    allow mismatches (default 1)
    returns the match (from barcodes) or None
    """
    accepted = None
    ne = str.__ne__
    barcode_lengths = map(len, barcodes)
    max_barcode_length = max(barcode_lengths)
    if max_barcode_length > len(seq):
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
        # fast implementation requires this condition
        accepted = None
        for barcode in barcodes:
            if mismatches > sum(imap(ne, barcode, seq)):
                if accepted is None: accepted = barcode
                else: return None
        return accepted

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
    
def truncate_record(record, max_length=0):
    '''
    truncate a record so that it is at most max_length
    starting at the 5' end 
    expects (barcode, (title, seq, qual))
    '''
    barcode, (title, seq, qual) = record
    return (barcode, (title, seq[0:max_length], qual[0:max_length]))

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
    logger = scripter.get_logger()
    if pf.paired_end: return split_paired_files(pf, logger=logger, **kwargs)
    else: return split_file(pf, logger=logger, **kwargs)

def split_file(fp_obj, no_gzip=False,
               barcodes=None, linker=None, min_length=4,
               max_length=None, logger=None,
               verbose=False, **kwargs):
    filename = fp_obj.input_file
    open_func, format = discover_file_format(filename)
    f = open_func(filename, "rU")
    records = FastqGeneralIterator(f)
    
    barcoded_files = {}
    filenames = []
    output_filename = partial(fp_obj.output_filename, no_gzip=no_gzip)
    if no_gzip: open_func = open
    else: open_func = GzipFile
    if barcodes is not None:
        processed_file = None
        for barcode in barcodes:
            fname = output_filename(barcode)
            filenames.append(fname)
            barcoded_files[barcode] = open_func(fname, 'w')
        
        # and make a unmatched file
        unmatched_filename = output_filename("unmatched")
        filenames.append(unmatched_filename)
        unmatched_file = open_func(unmatched_filename, 'w')
    else:
        barcoded_files = None
        unmatched_file = None
        processed_filename = output_filename("processed", is_barcode=False)
        filenames.append(processed_filename)
        processed_file = open_func(processed_filename, 'w')

    # assigned_records is a list of (barcode, record) items
    writer = partial(write_record, barcoded_files=barcoded_files,
                     unmatched_file=unmatched_file,
                     processed_file=processed_file,
                     linker=linker, min_length=min_length)
    final_records = apply_plan(records, barcodes=barcodes, linker=linker,
                               max_length=max_length,
                               **kwargs)
        
    results = map(writer, final_records)
    linker_only = results.count('linker')
    too_short = results.count('short')
    record_count = len(results)
    
    # close and exit #
    f.close()
    if barcoded_files is not None:
        for f_ in barcoded_files.values(): f_.close()
    if unmatched_file is not None: unmatched_file.close()
    if processed_file is not None: processed_file.close()
    
    logger.info('Split %s as %s ', fp_obj.input_file, ', '.join(filenames))
    logger.info('Processed %s records', record_count)
    logger.info('%s linker only dimers', linker_only)
    logger.info('%s sequences too short (1-3 bp)', too_short)
    
def assign_record(record, barcodes=None):
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
                 processed_file=None,
                 linker=None, min_length=4):
    """
    write the assigned_record to the correct barcode file
    an assigned_record is a (barcode, (title, seq, qual)) (all strs)
    
    will not work unless you provide a dictionary of barcoded_files and an
    unmatched_file object
    """
    barcode, (title, seq, qual) = assigned_record
    # check if record 1 is valid
    if len(seq) == 0 and barcode == linker and barcode is not None:
        problem = 'linker'
    elif len(seq) < min_length:
        problem = 'short'
    else:
        problem = None
        
    # produce lines
    line_fmt = "@{0!s}\n{1!s}\n+{0!s}\n{2!s}\n"
    line = line_fmt.format(title, seq, qual)
    is_processed = barcoded_files is None
    # abort write if there's a problem with the read
    if problem is not None:
        return problem
    if is_processed: processed_file.write(line) 
    elif barcode is None: unmatched_file.write(line)
    else: barcoded_files[barcode].write(line)
        
def write_record_pair(assigned_record_pair,
                      barcoded_file_pairs=None,
                      unmatched_files=None,
                      processed_files=None,
                      orphaned_read_files=None,
                      mismatched_files=None,
                      linker=None, min_length=4):
    """
    write the assigned_record to the correct barcode file
    an assigned_record is a (barcode, (title, seq, qual)) (all strs)
    
    will not work unless you provide a dictionary of barcoded_file_pairss and
    unmatched_files and mismatched_files (2-tuples of file object)
    """
    barcode, (title, seq, qual) = assigned_record_pair[0]
    # check if record 1 is valid
    if len(seq) == 0 and barcode == linker: problem = 'linker'
    elif len(seq) < min_length: problem = 'short'
    else: problem = None
    
    barcode2, (title2, seq2, qual2) = assigned_record_pair[1]
    # check if record 2 is valid
    if len(seq2) == 0 and barcode2 == linker: problem2 = 'linker'
    elif len(seq2) < min_length: problem2 = 'short'
    else: problem2 = None

    # produce lines
    line_fmt = "@{0!s}\n{1!s}\n+{0!s}\n{2!s}\n"
    line = line_fmt.format(title, seq, qual)
    line2 = line_fmt.format(title2, seq2, title2, qual2)
    
    is_processed = barcoded_file_pairs is None
    
    # abort write if both reads have problems
    if problem is not None or problem2 is not None:
    # if only one read has a problem, use write_record instead
        if problem2 is not None:
            if is_processed: orphaned_read_files[0].write(line)
            else: mismatched_files[0].write(line)
        else:
            if is_processed: orphaned_read_files[1].write(line2)
            else: mismatched_files[1].write(line2)
        return (problem, problem2)
        
    # select output files
    if barcoded_file_pairs is None:
        output_file = processed_files[0]
        output_file2 = processed_files[1]
    elif barcode is None and barcode2 is None:
        output_file = unmatched_files[0]
        output_file2 = unmatched_files[1]
    elif not barcode == barcode2:
        output_file = mismatched_files[0]
        output_file2 = mismatched_files[1]
    else:
        output_file = barcoded_file_pairs[barcode][0]
        output_file2 = barcoded_file_pairs[barcode][1]

    # write and return
    output_file.write(line)
    output_file2.write(line2)
    return (problem, problem2)
        
def is_empty_record(record):
#    barcode, (title, seq, qual) = record
    return len(record[1][1]) == 0

def pair_has_empty_record(record_pair):
    return is_empty_record(record_pair[0]) or is_empty_record(record_pair[1])
        
def apply_plan(records, barcodes=None,
               min_length=4, max_length=None,
               strip_after_barcode = 1, strip_before_barcode=0,
               **kwargs):
    if barcodes is None:
        assigner = partial(lambda r: (None, r), barcodes=barcodes)
    else:
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
    else:
        further_trimmed_records = imap(trim_trailing_Ns, assigned_records)
    if max_length is not None:
        truncator = partial(truncate_record, max_length=max_length)
        neat_records = imap(truncator, further_trimmed_records)
    else: neat_records = further_trimmed_records
    return neat_records
        

def split_paired_files(fp_obj, no_gzip=False,
                       barcodes=None, linker=None,
                       min_length=4, logger=None,
                       verbose=False, **kwargs):
    filename = fp_obj.input_file
    filename2 = fp_obj.second_file
    open_func = discover_file_format(filename)[0]
    open_func2 = discover_file_format(filename2)[0]
    f = open_func(filename, "rU")
    f2 = open_func2(filename2, "rU")
    records = FastqGeneralIterator(f)
    records2 = FastqGeneralIterator(f2)
        
    barcoded_file_pairs = {}
    filenames = []
    if no_gzip: open_func = open
    else: open_func = GzipFile
    output_filename = partial(fp_obj.output_filename, no_gzip=no_gzip)
    output_filename2 = partial(fp_obj.output_filename2, no_gzip=no_gzip)
    if barcodes is not None:
        processed_files = None
        orphaned_read_files = None
        for barcode in barcodes:
            fname = output_filename(barcode)
            fname2 = output_filename2(barcode)
            filenames.extend((fname, fname2))
            barcoded_file_pairs[barcode] = (open_func(fname, 'w'), 
                                            open_func(fname2, 'w'))
        
        # and make a unmatched file
        unmatched_filename = output_filename("unmatched")
        unmatched_filename2 = output_filename2("unmatched")
        unmatched_files = (open_func(unmatched_filename, 'w'),
                           open_func(unmatched_filename2, 'w'))
    
        mismatched_filename = output_filename("mismatched")
        mismatched_filename2 = output_filename2("mismatched")
        mismatched_files = (open_func(mismatched_filename, 'w'),
                            open_func(mismatched_filename2, 'w'))
        filenames.extend((unmatched_filename, unmatched_filename2,
                          mismatched_filename, mismatched_filename2))
    else:
        barcoded_file_pairs = None
        unmatched_files = None
        mismatched_files = None
        orphaned_read_filename = output_filename("orphaned",
                                                 is_barcoded=False)
        orphaned_read_filename2 = output_filename2("mismatched",
                                                   is_barcoded=False)
        orphaned_read_files = (open_func(orphaned_read_filename, 'w'),
                               open_func(orphaned_read_filename2, 'w'))
        processed_filename = output_filename("processed", is_barcoded=False)
        processed_filename2 = output_filename2("processed", is_barcoded=False)
        processed_files = (open_func(processed_filename, 'w'),
                           open_func(processed_filename2, 'w'))
        filenames.extend((orphaned_read_filename, orphaned_read_filename2),
                         (processed_filename, processed_filename2))
        
    writer = partial(write_record_pair,
                     barcoded_file_pairs=barcoded_file_pairs,
                     unmatched_files=unmatched_files,
                     mismatched_files=mismatched_files,
                     processed_files=processed_files,
                     orphaned_read_files=orphaned_read_files,
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
    
    logger.info('Split %s, %s as %s', fp_obj.input_file, fp_obj.second_file,
                ', '.join(filenames))
    logger.info('Processed %s records', record_count)
    logger.info('%s linker only dimers', linker_only)
    logger.info('%s sequences too short (1-3 bp)', too_short)


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

if __name__== "__main__": main()
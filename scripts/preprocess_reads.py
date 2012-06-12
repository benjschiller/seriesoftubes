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
from scripter import Environment, get_logger
from seriesoftubes.converters.discover import discover_file_format, \
         PATH_TO_GZIP, gzip_class_factory
from seriesoftubes.fnparsers import BarcodeFilenameParser
from seriesoftubes.cPreprocess import *
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from ConfigParser import SafeConfigParser
from errno import ENOENT, EACCES
from os import access, strerror, R_OK
from os.path import exists
from argparse import ArgumentTypeError
from gzip import GzipFile
from pkg_resources import get_distribution
__version__ = get_distribution('seriesoftubes').version
VERSION = __version__

def valid_seq(s):
        if not s.isalpha() and len(s)>0:
            msg = '%s is not a valid sequence'
            raise ArgumentTypeError(msg)
        else:
            return s.upper()
        
def main():
    e = Environment(version=VERSION, doc=__doc__)
    parser = e.argument_parser
    parser.add_argument('--strip-after-barcode', default=1, type=int,
                        help="""strip n bases after the barcode is removed (5' end)
(by default this 1 now, and is ignored if GERALD handled the barcoding)""")
    parser.add_argument('--strip-before-barcode', default=0, type=int,
                        help="""strip n bases before the barcode is removed (5' end)
(by default this 0 now, and is ignored if GERALD handled the barcoding)""")
    parser.add_argument('--min-length', type=int, default=4,
                        help='require sequences to be at least n total bases of non-N sequence (default: ignore)')
    parser.add_argument('--max-length', type=int, default=-1,
                        help='truncate final sequences to n bases (default: ignore)')
    parser.add_argument('--no-gzip', default=False, action='store_true',
                        help = 'Do not gzip output files')
    bgroup = parser.add_argument_group('barcodes',
                                       'Specify sequence barcodes in the sample(s)')
    bgroup.add_argument('-b', '--barcodes', action='append', type=valid_seq,
                        help="Specify a barcode sequence. May be invoked multiple times")
    bgroup.add_argument('--kry-barcodes', dest='barcodes', action='store_const',
                        help='Alias for -bTCAT -bGACG -bAGTC -bCTGA',
                        const=['TCAT', 'GACG', 'AGTC', 'CTGA'])
    parser.add_argument('--linker', default='', type=valid_seq,
                        help="Specify a 3' adaptor/linker sequence that we should clip off of each read")
    parser.set_defaults(**{'target': 'processed'})
    e.set_filename_parser(BarcodeFilenameParser)
    e.set_config_reader(read_config)
    e.set_config_writer(write_config)
    e.do_action(splitter)

def write_config(barcodes=None, min_length=None, max_length=None, strip_after_barcode=0,
                 strip_before_barcode=0, linker=None, target_dir = os.curdir,
                 target=None,
                 *args, **kwargs):
    config = SafeConfigParser()
    config.add_section('main')
    if barcodes is not None:
        config.set('main', 'barcodes', ','.join(barcodes))
    config.set('main', 'linker', linker)
    config.set('main', 'min-length', str(min_length))
    config.set('main', 'max-length', str(max_length))
    config.set('main', 'strip-after-barcode', str(strip_after_barcode))
    config.set('main', 'strip-before-barcode', str(strip_before_barcode))
    config_f = os.path.join(target_dir, 'preprocess_reads.cfg')
    if not os.path.exists(target_dir):
        os.makedirs(target_dir, mode=0755)
    with open(config_f, 'wb') as configfile:
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
    if config.has_option('main', 'max-length'):
        context['max_length'] = int(config.get('main', 'max-length'))
    if config.has_option('main', 'min-length'):
        context['min_length'] = int(config.get('main', 'min-length'))
    if config.has_option('main', 'linker'):
        context['linker'] = config.get('main', 'linker')
    context['strip_after_barcode'] = int(config.get('main', 'strip-after-barcode'))
    context['strip_before_barcode'] = int(config.get('main', 'strip-before-barcode'))
    return context

def splitter(pf, **kwargs):
    logger = get_logger()
    if pf.paired_end:
        return split_paired_files(pf, logger=logger,
                                  **kwargs)
    else:
        return split_file(pf, logger=logger,
                          **kwargs)

def split_file(fp_obj, no_gzip=False,
               barcodes=[], linker='', min_length=4,
               max_length=-1, logger=None,
               strip_before_barcode=0,
               strip_after_barcode=0, 
               **kwargs):
    if logger is None: logger = get_logger()
    filename = fp_obj.input_file
    open_func, format_ = discover_file_format(filename)
    if not format_ == 'FASTQ':
        logger.error('Only FASTQ files are supported at this time')
        return
    f = open_func(filename, "rU")
    records = FasterFastqIterator(f)
    
    barcoded_files = {}
    filenames = []
    output_filename = partial(fp_obj.output_filename, no_gzip=no_gzip)
    if no_gzip: open_func = open
    elif PATH_TO_GZIP is not None:
        open_func = gzip_class_factory(PATH_TO_GZIP)
    else: open_func = GzipFile
    if barcodes is None: barcodes = []
    if len(barcodes) > 0:
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

    writer_args = {'barcoded_files': barcoded_files,
              'unmatched_file': unmatched_file,
              'processed_file': processed_file}
    results = apply_plan(records, writer_args, barcodes=barcodes, linker=linker,
                               min_length=min_length, max_length=max_length,
                               strip_after_barcode=strip_after_barcode,
                               strip_before_barcode=strip_before_barcode,
                               logger=logger)
    linker_only = results['linker']
    too_short = results['short']
    record_count = results['all']
    # close and exit #
    f.close()
    if barcoded_files is not None:
        logger.debug('closing barcoded files')
        for f_ in barcoded_files.values(): f_.close()
    if unmatched_file is not None:
        logger.debug('closing unmatched file')
        unmatched_file.close()
    if processed_file is not None:
        logger.debug('closing output file')
        processed_file.close()
    
    logger.info('Split %s as %s ', fp_obj.input_file, ', '.join(filenames))
    logger.info('Processed %s records', record_count)
    logger.info('%s linker only dimers', linker_only)
    logger.info('%s sequences too short (1-3 bp)', too_short)
    
def split_paired_files(fp_obj, no_gzip=False,
                       barcodes=None, linker='',
                       min_length=4,
                       max_length = -1,
                       strip_before_barcode=0,
                       strip_after_barcode=0,
                       logger=None,
                       **kwargs):
    filename = fp_obj.input_file
    filename2 = fp_obj.second_file
    open_func, format1 = discover_file_format(filename)
    open_func2, format2 = discover_file_format(filename2)
    if not (format1 == 'FASTQ' and format2 == 'FASTQ'):
        logger.error('Only FASTQ files are supported at this time')
        return
    f = open_func(filename, "rU")
    f2 = open_func2(filename2, "rU")
    records = FasterFastqIterator(f)
    records2 = FasterFastqIterator(f2)
        
    barcoded_file_pairs = {}
    filenames = []
    if no_gzip:
        open_func = open
    elif PATH_TO_GZIP is not None:
        open_func = gzip_class_factory(PATH_TO_GZIP)
    else:
        open_func = GzipFile
    output_filename = partial(fp_obj.output_filename, no_gzip=no_gzip)
    output_filename2 = partial(fp_obj.output_filename2, no_gzip=no_gzip)
    if barcodes is None: barcodes = []
    if len(barcodes) > 0:
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
        
    writer_args = {'barcoded_file_pairs': barcoded_file_pairs,
                   'unmatched_files': unmatched_files,
                   'mismatched_files': mismatched_files,
                   'processed_files': processed_files,
                   'orphaned_read_files': orphaned_read_files,
                   'linker': linker,
                   'min_length': min_length
                   }

    results = apply_plan_pe(records, records2, writer_args,
                            barcodes=barcodes, linker=linker,
                            min_length=min_length,
                            max_length=max_length,
                            strip_after_barcode=strip_after_barcode,
                            strip_before_barcode=strip_before_barcode,
                            logger=logger)
    linker_only = results['linker']
    too_short = results['short']
    record_count = results['all']

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


if __name__== "__main__": main()

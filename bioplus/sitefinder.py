"""
tools for dealing with binding sites (instances of sequence motifs)
"""
import os.path
from seqtools import rc
from tabfile import TabFile, MacsFile, BedFile
import operator
import Bio.Motif
import Bio.SeqIO
try:
    import MOODS
    USE_MOODS = True
except ImportError:
    USE_MOODS = False

def _clean_sequence(seq):
    '''just in case we are fed masked sequence, this returns
    the longest contiguous sequence without any N/n/X/x'''
    if 'N' in seq: seq = max(seq.split('N'),key=len)
    if 'n' in seq: seq = max(seq.split('n'),key=len)
    if 'X' in seq: seq = max(seq.split('X'),key=len)
    if 'x' in seq: seq = max(seq.split('x'),key=len)
    return seq

def _reverse_complement(forward_seq):
    '''just make sure we can reverse complement DNA sequences'''
    try:
        site_seq = forward_seq.reverse_complement()
    except AttributeError:
        site_seq = rc(forward_seq)
    return site_seq

def MOODS_search(seq, motif, thresholds=0):
    '''an equivalent to Motif.search_pwm()'''
    if not USE_MOODS: raise RuntimeError('MOODS could not be imported')
    sequence = seq
    matrix_ = MOODS.transpose([ map(lambda x: x[1], sorted(x.items()))
                for x in motif.log_odds()])
    # Note: algorithm = 'lf' fails due to segmentation fault
    results_per_matrix = MOODS.search(sequence, [matrix_], thresholds,
            bg=None, algorithm='pla', q=7, absolute_threshold=True,
            both_strands=True, combine=True)

    # format as Motif.search_pwm results
    search_results = results_per_matrix[0]
    # figure out direction of reverse results
    # do we need to reverse it?
    results_sorted_like_Bio_Motif = sorted(search_results, 
            key=operator.itemgetter(0),
            cmp=lambda x,y: cmp(abs(x), abs(y)))
    return results_sorted_like_Bio_Motif

def search_peak(peak_ID, peak, peakseq, motif):
    '''provide information about matches to a motif in a
    peak region, and about the region
    
    peak MUST provide EITHER
    (1) the following public methods
    chrom = reference (e.g. chr1, chrX)
    chromStart = start coordinate, 0-based
    chromEnd = end coordinate, open
    or (2) the following public method
    coordinates = a tuple containing (chrom, chomStart, chromEnd)

    peak may optionally provide the following methods
    tags (if not found, we will replace with 'NA')
    summit (if not found, we will use the peak center)
    misc (a list of anything else)

    For each peak, the best motif hit is returned where
    best is defined as the motif hit with the most information
    and closest to the center (in the case of ties)

    Note site position is 0-based, in contrast with
    earlier versions of biotools

    returns a tuple of four things:
    peak info
    peak BED row
    a list of info about sites (motif matches)
    a list of BED rows for sites
    '''
   
    clean_peak_seq = _clean_sequence(peakseq)
    clean_peak_length = len(clean_peak_seq)
    clean_offset = peakseq.find(clean_peak_seq)

    site_count = 0
    total_Ri = 0
    best_Ri = 0
    best_seq = 'NA'
    best_position = 'NA'
    best_strand = 'NA'
    sites_info_rows = []
    sites_bed_rows = []

    motif_length = len(motif)

    try:
        peak_coordinates = peak.coordinates()
    except AttributeError:
        peak_coordinates = [peak.chrom(),
                    peak.chrom_start(),
                    peak.chrom_end()]

    try:
        peak_center = peak.summit()
        if int(peak_center) > int(peak_coordinates[1]):
            peak_center -= int(peak_coordinates[1])
            
    except AttributeError:
        peak_center = clean_peak_length / 2 + clean_offset

    try:
        peak_intensity = peak.tags()
    except (ValueError, AttributeError):
        try:
            peak_intensity = peak.tagsv2()
        except AttributeError:
            peak_intensity = 'NA'

    try: peak_misc = peak.misc()
    except AttributeError: peak_misc = []

    if USE_MOODS: search_results = MOODS_search(clean_peak_seq, motif)
    else: search_results = motif.search_pwm(clean_peak_seq)

    for position, Ri in search_results:
        site_count += 1
        total_Ri += Ri
        
        if position > 0:
            offset = position + clean_offset
            site_seq = clean_peak_seq[offset:offset + motif_length]
            strand = '+'
        else:
            offset = clean_offset + clean_peak_length + position - motif_length
            site_seq = _reverse_complement(clean_peak_seq[offset:
                            offset + motif_length])
            strand = '-'

        site_ID = '_'.join([str(x) for x in
                    [peak_ID, 'motif', site_count]])

        chrom = peak_coordinates[0]
        chrom_start = int(peak_coordinates[1])
        site_BED_row = [chrom, chrom_start + offset,
                    chrom_start + offset + motif_length,
                    site_ID, Ri, strand]

        site_info_row = site_BED_row + [offset, site_seq,
                        peak_ID, clean_peak_length,
                        len(peak), peak_center]

        sites_bed_rows.append(site_BED_row)
        sites_info_rows.append(site_info_row)

        # check if this is the best site
        if Ri > best_Ri or \
            Ri == best_Ri and (abs(peak_center - offset) - motif_length/2) < \
                    best_position:
            best_Ri = Ri
            best_seq = str(site_seq)
            best_position = offset
            best_strand = strand
    

    peak_info_row = peak_coordinates + [peak_ID, peak_intensity,
        site_count, total_Ri,
        best_Ri, best_seq, best_position, best_strand, 
        clean_peak_length, peak_center, ','.join(peak_misc)]

    peaks_BED_row = peak_coordinates + [peak_ID, 1000, '+',
        int(peak_coordinates[1]) + peak_center,
        int(peak_coordinates[1]) + peak_center + 1]

    return (peak_info_row, peaks_BED_row, sites_info_rows, sites_bed_rows)

def find_sites(peaks_file, fasta_file, motif, bed=True, xls=False,
               output_dir = None, motif_type = 'MEME',
               src_fnc="find_sites", **kwargs):
    '''
findSites(peaks_file,FASTAfile,motif) takes the NAME_peaks.xls file outputed
by MACS, as well as a FASTAfile, and finds instances of the motif specified by
motif (a Bio.Motif object). It will output two new files for peaks and sites
called NAME.peaks.info and NAME.sites.info. It will also create files called
NAMES.peaks.bed and NAME.sites.bed which are proper BED files (scores are tag
density, and information content, respectively). All files are 0-based,
half-open in line with the BED convention. MACS coordinates are corrected
accordingly.

f.peaks.info contains
Peak (1) chr, (2) start (3) end
(4) Peak ID
(5) Relative summit
(6) Number of unique tags in peak region
(7) -10*log10(pvalue)
(8) fold_enrichment
(9) FDR
(10) # motif instances found
(11) Total Ri for discovered motif instances
(12) Greatest Ri of any motif in peak region
(13) Sequence of that motif instance
(14) Position (offset) of that motif (left-end)

f.peaks.bed contains
Peak (1) chr, (2) start (3) end
(4) Peak ID
(5) Number of unique tags in peak region
(6) Strand .
(7) Summit position (absolute)
(8) Summit position + 1

f.sites.info contains
Site (1) chr (2) start (3) end
(4) Unique Site ID (internally generated)
(5) The motif information content Ri, in bits
(6) motif orientation, best score (+) or (-)
---- BED file ends here ----
(7) the motif sequence (e.g., ACAACA)
(8) Position (offset) of that motif (left-end)
(9) peak ID, fetched from MACS
(10) used peak length
(11) true peak length
(11) peak summit offset
    '''

    if type(motif) is str:
        motif = Bio.Motif.read(open(motif), motif_type)


    # start the output file
    prefix = os.path.splitext(os.path.basename(peaks_file))[0]
    if output_dir is not None: prefix = os.path.join(output_dir, prefix)
    sites_info = TabFile(os.extsep.join([prefix, 'sites', 'info']), 'w')
    sites_bed = TabFile(os.extsep.join([prefix, 'sites', 'bed']), 'w')
    peaks_info = TabFile(os.extsep.join([prefix, 'peaks', 'info']), 'w')
    peaks_bed = TabFile(os.extsep.join([prefix, 'peaks', 'bed']), 'w')

    peaks_cols = ['chr', 'start', 'end', 'peak_ID', 'peak_intensity',
                  'site_count', 'total_Ri', 'best_Ri', 'best_seq',
                  'best_offset', 'best_strand', 'clean_peak_length',
                  'peak_summit', 'peak_misc']
    peaks_msg = os.linesep.join(['# This file was generated by ' + src_fnc,
                                 '# comments are retained from original file',
                                 "\t".join(peaks_cols), ''])
    peaks_info.write(peaks_msg)
    sites_cols = ['chr', 'start', 'end', 'site_ID', 'Ri', 'strand', 'offset',
                  'motif_seq', 'peak_ID', 'peak_length', 'reported_peak_length',
                  'peak_summit']
    sites_msg = os.linesep.join(['# This file was generated by ' + src_fnc,
                                 "\t".join(sites_cols), ''])
    sites_info.write(sites_msg)

    if bed: peak_generator = BedFile(peaks_file)
    elif xls: peak_generator = MacsFile(peaks_file)
    else: raise ValueError('Neither bed nor xls')
    # peakSeqs is a generator
    peak_seqs = (r.seq for r in Bio.SeqIO.parse(open(fasta_file, 'rU'),
                                                'fasta'))
    nosites = 0
    peaknumber = 0
    
    for peak in peak_generator:
#            if peaknumber%10000 is 0: print peaknumber
        peaknumber += 1
        seq = peak_seqs.next()
        # Generate a peak ID
        peak_ID = '{!s}_{!s}'.format(prefix, peaknumber) 
        (peak_info, peak_bed,
         sites_info_rows, sites_bed_rows) = search_peak(peak_ID, peak, seq,
                                                        motif)
        peaks_info.write_row(peak_info)
        peaks_bed.write_row(peak_bed)
        sites_info.write_rows(sites_info_rows)
        sites_bed.write_rows(sites_bed_rows)
        if len(sites_info_rows) is 0: nosites += 1

    sites_info.close()
    sites_bed.close()
    peaks_info.close()
    peaks_bed.close()
    message = "There were {!s} of {!s} peaks with no identifiable \
sites in {!s} using a cutoff of 0 bits".format(nosites, peaknumber, fasta_file)
    stdout_buffer = message
    # get the motif
    motif_str = ''
    try:
        motif_str = os.linesep.join([
                        ", ".join([ "".join([str(base), ": ", str(odds)]) 
                                for base, odds in position.items()])
                        for position in motif.log_odds()])
    except AttributeError:
        motif_str = str(motif)
    message = os.linesep.join([message, 'The following motif was used',
                               motif_str])
    # print message and write it to a log
    g = open(prefix + '.log', 'w')
    g.write(message)
    g.close()
    return stdout_buffer
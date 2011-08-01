"""
a pure python implementation of the SAM/BAM file formats
"""
import time
from itertools import imap
from struct import calcsize, pack

ID1 = 31
ID2 = 139
CM = 8
FLG = 4
XFL = 0
OS = 255
XLEN = 6
SI1 = 66
SI2 = 67
SLEN = 2

if not calcsize('B') == 1 or not calcsize('H') == 2 or not calcsize('I') == 4:
    raise ImportError('The system you are using does not have ' \
                      'properly sized types char, short and int')

def create_BGZF_header(file_obj, mtime=None, BSIZE=0):
    """
    write a valid BGZF header to file_obj at its current position
    returns the file position afterward
    """
    if mtime is None: mtime = int(time.time())
    h = pack('<BBBBIBBHBBHH', ID1, ID2, CM, FLG, mtime, XFL, OS, XLEN,
             SI1, SI2, SLEN, BSIZE)
    return h

from codecs import make_identity_dict

def make_seq_encoding_map():
    m = dict((i, 15) for i in range(256))
    for i, c in enumerate('=ACMGRSVTWYHKDBN'):
        m[ord(c)] = i
    return m

SEQ_ENCODING_MAP = make_seq_encoding_map()                

def bam_record(sam_record):
    """
    convert a SAM alignment record to a BAM alignment record
    """
    r = sam_record.split('\t')
    # sam fields
    QNAME = r[0]
    FLAG = r[1]
    RNAME= r[2]
    POS = r[3]
    MAPQ = r[4]
    CIGAR = r[5]
    RNEXT = r[6]
    PNEXT = r[7]
    TLEN = r[8]
    SEQ = r[9]
    QUAL = r[10]
    if len(r) > 11: opt_fields = r[11:]
    else: opt_fields = []
    # bam fields
    pos = POS - 1
    read_name = QNAME + '\0'
    l_read_name = len(read_name)
    next_pos = PNEXT - 1
    tlen = TLEN
    bam_record_header = pack('iiiIIiiii', block_size, refID, pos, bin_mq_nl,
                             flag_nc, l_seq, next_refID, next_pos, tlen)
    cigar = 
    seq_nybbles, l_seq = charmap_encode(SEQ, strict, SEQ_ENCODING_MAP)
    high_nybbles = seq_nybbles[::2]
    low_nybbles = seq_nybbles[::2]
    

def convert_to_bam_iterator(iter):
    """
    converts an iterator iter of (name, seq, qual) records
    into an interator of BAM alignment records
    note: these records are NOT compressed
    
    intended to be used with GeneralFastqIterator from Biopython
    """
    return starmap(bam_alignment_record, iter)
import gzip
import bz2

def discover_file_format(filename):
    f = open(filename, 'rb')
    head = f.read(3)
    f.close()
    # check magic words for compression
    if head == '\x1f\x8b\x08':
        print 'gzip'
        open_func = gzip.GzipFile
    elif head=='\x42\x5a\x68':
        open_func = bz2.BZ2File
    else:
        open_func = open
    uncompressed = open_func(filename)
    head2 = uncompressed.read(4)
    # check for BAM
    if head2 == 'BAM\x01': return (open_func, 'BAM')
    # check for SAM 
    head2b = uncompressed.readline()
    if head2 == '@HD\t':
        if head2b[0:3] == 'VN:': return (open_func, 'SAM')
    # check fastq
    title = head2 + head2b
    seq = uncompressed.readline()
    title2 = uncompressed.readline()
    qual = uncompressed.readline()
    if len(seq) == len(qual) and \
       (title.startswith('@') or title2.startswith('+')) and \
       (title2[1:].strip() == '' or title[1:] == title2[1:]):
        return (open_func, 'FASTQ')
    # otherwise give up
    else: return (None, None)

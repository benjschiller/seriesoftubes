from itertools import izip
import numpy
def MACS_track_to_iter(track):
    i = 0
    for end_pos, value in izip(*track):
        while i < end_pos:
            i += 1
            yield value
            
def array_to_bedgraph(a, chrom, write_zero_values=False,
                      precision=0):
    """ generator that takes a numpy array and yields tuples bedGraph file
    note: the array must be compatible with enumerate()
    chrom should be specified too
    
    precision specifies how many places past the decimal to retain (default 1)
    
    To coerce into a string use something like
    for line in array_to_bedgraph(a, chrom):
        '%s\t%d\t%d\t%.1f\n' % line
    """
    last_pos = 0
    last_value = round(a[0], precision)
    for pos, value in enumerate(numpy.around(a, decimals=precision)):
        if value == last_value: continue
        else:
            if last_value != 0 or write_zero_values:
                yield (chrom, last_pos, pos, last_value) 
            last_pos = pos
            last_value = value
    if last_value != 0:
        yield (chrom, last_pos, pos, last_value)

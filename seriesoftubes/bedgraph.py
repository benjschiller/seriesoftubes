from itertools import repeat

def in_window_factory(chrom, start, end):
    """
    function factory
    
    returns a function that checks if a bedgraph line overlaps a given
    interval. along with UCSC convention, end coordinate is open and start
    coordinate is 0-based
    (also compatible with BED files)
    """
    def in_window(line):
        parts = line.rstrip('\r\n').split()
        if len(parts) < 3: return False
        this_chrom, this_start, this_end = parts[0:3]
        if not this_chrom == chrom: return False
        elif not int(this_start) < end: return False
        elif not int(this_end) > start: return False
        else: return '\t'.join(parts)     
    return in_window

def in_windows_factory(intervals):
    """
    function factory
    
    returns a function that checks if a bedgraph line overlaps any given
    interval. assumes we already filtered by chrom
    
    intervals must be given as tuples (int start, int end)
    
    along with UCSC convention, end coordinate is open and start
    coordinate is 0-based
    (also compatible with BED files)
    """
    for rec in intervals:
        rec = 0 
    return in_window

def window_to_iter(chrom, start, end, filename, my_type=int):
    """
    generator that returns the values in a region [start,end) on chrom
    
    from a bedgraph file called filename
    
    note: bedgraph values are assumed to be ints unless my_type is set to
    something  else
    """
    in_window = in_window_factory(chrom, start, end)
    output = filter(in_window, open(filename, 'rU'))
    if len(output) == 0:
        for i in repeat(0, end - start):
            yield i
    first_range = output[0].rstrip('\r\n').split('\t')
    last_range = output[0].rstrip('\r\n').split('\t')
    
    previous_end = start
    for line in output:
        this_range = line.rstrip('\r\n').split()
        this_start = int(this_range[1])
        for i in xrange(previous_end, this_start): yield 0
        this_end = int(this_range[2])
        this_value = my_type(this_range[3])
        for i in xrange(this_start, this_end): yield this_value
        previous_end = this_end
        
    for i in xrange(last_end, end): yield 0
    raise StopIteration
        
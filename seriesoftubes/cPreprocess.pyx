# cython: profile=True
from cpython cimport bool

DEF ERR_TOO_SHORT = 1
DEF ERR_LINKER = 2

cdef struct int_pair:
    int left
    int right

def FasterFastqIterator(fastq_handle):
    """Cython-ized version of FastqGeneralIterator
    """
    cdef:
        str line, title_line, seq_string, second_title, quality_string
        int seq_len
    #We need to call handle.readline() at least four times per record,
    #so we'll save a property look up each time:
    
    #Skip any text before the first record (e.g. blank lines, comments?)
    handle_readline = fastq_handle.readline
    while True:
        line = handle_readline()
#        if line == "" : return #Premature end of file, or just empty?
        if line[0] == "@":
            break

    while True:
        if line[0] != "@":
            raise ValueError("Records in Fastq files should start with '@' character")
        title_line = line[1:].rstrip()
        #Will now be at least one line of quality data - in most FASTQ files
        #just one line! We therefore use string concatenation (if needed)
        #rather using than the "".join(...) trick just in case it is multiline:
        seq_string = handle_readline().rstrip()
        #There may now be more sequence lines, or the "+" quality marker line:
        while True:
            line = handle_readline()
            if line == '':
                raise ValueError("End of file without quality information.")
            ## disable second title checking
            if line[0] == "+":
#                #The title here is optional, but if present must match!
                second_title = line[1:].rstrip()
#                if second_title and second_title != title_line:
#                    raise ValueError("Sequence and quality captions differ.")
                break
            seq_string += line.rstrip() #removes trailing newlines
        #This is going to slow things down a little, but assuming
        #this isn't allowed we should try and catch it here:
        ## disable character checking
#        if " " in seq_string or "\t" in seq_string:
#            raise ValueError("Whitespace is not allowed in the sequence.")
        seq_len = len(seq_string)

        #Will now be at least one line of quality data...
        quality_string = handle_readline().rstrip()
        #There may now be more quality data, or another sequence, or EOF
        while True:
            line = handle_readline()
            if line == '': break #end of file
            if line[0] == "@":
                #This COULD be the start of a new sequence. However, it MAY just
                #be a line of quality data which starts with a "@" character.  We
                #should be able to check this by looking at the sequence length
                #and the amount of quality data found so far.
                if len(quality_string) >= seq_len:
                    #We expect it to be equal if this is the start of a new record.
                    #If the quality data is longer, we'll raise an error below.
                    break
                #Continue - its just some (more) quality data.
            quality_string += line.rstrip()
        
        if seq_len != len(quality_string):
            raise ValueError("Lengths of sequence and quality values differs "
                             " for %s (%i and %i)." \
                             % (title_line, seq_len, len(quality_string)))

        #Return the record and then continue...
        yield (title_line, seq_string, quality_string)
        if line == '' : return #StopIteration at end of file
#    assert False, "Should not reach this line"

cdef struct Read:
    char *title
    char *seq
    char *qual

cdef Read asRead(str title, str seq, str qual):
    cdef Read read
    read.title = <bytes>title
    read.seq = <bytes>seq
    read.qual = <bytes>qual
    return read

cdef struct Record:
    char *barcode
    Read read

cdef Record asRecord(str barcode, str title, str seq, str qual):
    cdef Record record
    record.barcode = <bytes>barcode
    record.read = asRead(title, seq, qual)
    return record

cdef struct RecordPair:
    Record first
    Record second

cdef RecordPair pair(Record record1, Record record2):
    cdef RecordPair rp
    rp.first = record1
    rp.second = record2
    return rp

cdef char *match_barcode(char *seq, list barcodes, int mismatches=1):
    """
    try to match seq to a list of barcodes
    allow mismatches (default 1)
    returns the match (from barcodes) or None
    """
    cdef:
        list barcode_lengths = map(len, barcodes)
        int n_barcodes = len(barcodes)
        int max_barcode_length = max(barcode_lengths)
        int seq_len = len(seq)
        int barcode_length
        char *accepted = ''
        char *barcode, comparable
        int hamming, i, j
    for i in range(n_barcodes):
        barcode = barcodes[i]
        barcode_length = barcode_lengths[i]
        if seq_len < barcode_length: continue
        hamming = 0
        for j in range(barcode_length):
            if barcode[j] != seq[j]: hamming += 1
        if mismatches > hamming:
            if accepted == b'': accepted = barcode
            else: return b''
    return accepted

cdef Read pretrim_record_5prime(Read read, int trim_length):
    '''
    trim an assigned record from the 5' end
    expects (barcode, (title, seq, qual))
    '''
    read.seq += trim_length
    read.qual += trim_length 
    return read

cdef Record trim_record_5prime(Record record, int trim_length):
    '''
    trim an assigned record from the 5' end
    expects (barcode, (title, seq, qual))
    '''
    record.read.seq += trim_length
    record.read.qual += trim_length
    return record
    
cdef Record truncate_record(Record record, int max_length):
    '''
    truncate a record so that it is at most max_length
    starting at the 5' end 
    expects (barcode, (title, seq, qual))
    '''
    record.read.seq[max_length] = b'\0'
    record.read.qual[max_length] = b'\0'
    return record

cdef Record trim_record_3prime(Record record, int trim_length):
    '''
    trim a record from the 3' end
    expects (barcode, (title, seq, qual))
    '''
    if trim_length == 0: return record
    record.read.seq[-trim_length] = b'\0'
    record.read.qual[-trim_length] = b'\0'
    return record

cdef Record trim_trailing_Ns(Record record):
    '''
    trim any trailing 3' 'N's
    expects record is (barcode, (title, seq, qual))
    returns truncated (barcode, (title, seq, qual))
    '''
    cdef int end = len(str.rstrip(record.read.seq, 'N'))
    record.read.seq[end] = b'\0'
    record.read.qual[end] = b'\0'
    return record

cdef Record cleave_linker(Record record, char *linker):
    cdef:
        int i = str.find(record.read.seq, linker)
    if i == -1: return record
    else:
        record.barcode = linker
        record.read.seq[i] = b'\0'
        record.read.qual[i] = b'\0'
        return record

cdef Read as_read(bytes title, bytes seq, bytes qual):
    cdef Read read
    read.title = title
    read.seq = seq
    read.qual = qual
    return read

cdef Record assign_read(Read read, list barcodes):
    """
    Assign a record to a barcode
    returns (barcode, new_record)
    if unmatched, returns (None, record)
    """
#    title, seq, qual = record/
#    title_head, last_part = title.rsplit(':', 1)
    cdef:
        bytes title_head, last_part
        int pound_loc, slash_loc
        char *barcode
        int barcode_len
        int n_barcodes = len(barcodes)
        Record record
    record.read = read

    if n_barcodes == 0:
        record.barcode = b''
        return record
    title_head, last_part = read.title.rsplit(':', 1)
    if last_part.isalpha():
        # CASAVA 1.8 file
        record.barcode = match_barcode(<bytes>last_part.rstrip(), barcodes)
        return record 
    pound_loc = last_part.find('#')
    slash_loc = last_part[pound_loc:].find('/')
    if slash_loc != -1:
        barcode = <bytes>last_part[(pound_loc+1):(pound_loc+slash_loc)]
    else:
        barcode = <bytes>last_part[(pound_loc+1):]
    if barcode==b'0':
        record.barcode = match_barcode(read.seq, barcodes)
        barcode_len = len(record.barcode)
        record.read.seq += barcode_len
        record.read.qual += barcode_len
        if not record.barcode == <bytes>'':
            last_part = last_part[0:(pound_loc+1)] + barcode + last_part[(pound_loc+slash_loc):]
            record.read.title = <bytes>('%s:%s' %( title_head , last_part ))
    elif barcode.isdigit():
        # then we have a numbered index from Illumina, just use it as-is
        record.barcode = barcode
    elif barcode.isalpha():
        # then we already extracted the barcode at some point, try to match it
        record.barcode = match_barcode(barcode, barcodes)
    else:
        record.barcode = b''
    return record

cdef bool too_short(Read read, int min_length):
    cdef char c
    cdef int i = 0, L = 0
    for c in read.seq:
        L += 1
        if c == b'N': i += 1
    return L - i < min_length

cdef int write_record(Record record, barcoded_files=None,
                        unmatched_file=None,
                        processed_file=None,
                        str linker='', int min_length=4):
    """
    write the assigned_record to the correct barcode file
    an assigned_record is a (barcode, (title, seq, qual)) (all strs)
    
    will not work unless you provide a dictionary of barcoded_files and an
    unmatched_file object
    """
    cdef:
        str barcode = record.barcode
        str title = record.read.title
        str seq = record.read.seq
        str qual = record.read.qual
        str line = "@%s\n%s\n+%s\n%s\n" % (title, seq, title, qual)
        bool is_linker = (len(seq) == 0 and
                          barcode == linker
                          and not barcode == '')
        bool is_too_short = too_short(record.read, min_length)
        bool is_processed = barcoded_files is None
        bool no_barcode = (barcode == '')
    # check if record 1 is valid
    if is_linker: return ERR_LINKER
    elif is_too_short: return ERR_TOO_SHORT
        
    # produce lines
    if is_processed: processed_file.write(line) 
    elif no_barcode: unmatched_file.write(line)
    else: barcoded_files[barcode].write(line)
    return 0
        
cdef int_pair write_record_pair(RecordPair assigned_record_pair,
                      barcoded_file_pairs=None,
                      unmatched_files=None,
                      processed_files=None,
                      orphaned_read_files=None,
                      mismatched_files=None,
                      linker='', int min_length=4):
    """
    write the assigned_record to the correct barcode file
    an assigned_record is a (barcode, (title, seq, qual)) (all strs)
    
    will not work unless you provide a dictionary of barcoded_file_pairss and
    unmatched_files and mismatched_files (2-tuples of file object)
    """
    cdef:
        int_pair ret
        Record record = assigned_record_pair.first
        str barcode = record.barcode
        str title = record.read.title
        str seq = record.read.seq
        str qual = record.read.qual
        int problem = 0
        str line = "@%s\n%s\n+%s\n%s\n" % (title, seq, title, qual)
        bool is_linker1 = (len(seq) == 0
                           and not (linker == '')
                           and barcode == linker)
        bool is_too_short1 = too_short(record.read, min_length)
                
        Record record2 = assigned_record_pair.second
        str barcode2 = record2.barcode
        str title2 = record2.read.title
        str seq2 = record2.read.seq
        str qual2 = record2.read.qual
        int problem2 = 0
        str line2 = "@%s\n%s\n+%s\n%s\n" % (title2, seq2, title2, qual2)
        bool is_linker2 = (len(seq2) == 0
                           and not (linker == '')
                           and barcode2 == linker)
        bool is_too_short2 = too_short(record.read, min_length)
        
        str line_fmt
        bool is_processed = (barcoded_file_pairs is None)

    # check if record 1 is valid
    if is_linker1: problem = ERR_LINKER
    elif is_too_short1: problem = ERR_TOO_SHORT
    
    # check if record 2 is valid
    if is_linker2: problem2 = ERR_LINKER
    elif is_too_short2: problem2 = ERR_TOO_SHORT
    ret.left = problem
    ret.right = problem2

    # abort write if both reads have problems
    if (not problem=='') or (not problem2 == ''):
    # if only one read has a problem, use write_record instead
        if not problem2 == '':
            if is_processed: orphaned_read_files[0].write(line)
            else: mismatched_files[0].write(line)
        else:
            if is_processed: orphaned_read_files[1].write(line2)
            else: mismatched_files[1].write(line2)
        return ret
        
    # select output files
    if barcoded_file_pairs is None:
        output_file = processed_files[0]
        output_file2 = processed_files[1]
    elif barcode == '' and barcode2 == '':
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
    return ret


cpdef dict apply_plan(reads, writer_args, list barcodes=[], str linker='',
                 int min_length=4, int max_length = -1,
                 int strip_after_barcode = 1, int strip_before_barcode=0,
                 logger=None):
    cdef:
        int i, n_short=0, n_linker=0
        char *llinker = <bytes>linker
        int result
        Record record
        tuple read, read2
    for read in reads:
        i += 1
        record = apply_plan_to_read(read, barcodes, llinker,
                 min_length, max_length,
                 strip_after_barcode, strip_before_barcode)
        result = write_record(record,
                              writer_args['barcoded_files'],
                              writer_args['unmatched_file'],
                              writer_args['processed_file'],
                              linker,
                              min_length)
        if result == ERR_TOO_SHORT: n_short += 1
        elif result == ERR_LINKER: n_linker += 1
    return {'all': i, 'short': n_short, 'linker': n_linker}

cpdef dict apply_plan_pe(reads, reads2, writer_args, list barcodes=[],
                         str linker='',
                 int min_length=4, int max_length = -1,
                 int strip_after_barcode = 1, int strip_before_barcode=0,
                 logger=None):
    cdef:
        int i, n_short = 0, n_linker = 0
        char *llinker = <bytes>linker
        int_pair result
        Record record
        tuple read, read2
    for read in reads:
        i += 1
        read2 = read2.next()
        record = apply_plan_to_read(read, barcodes, llinker,
                 min_length, max_length,
                 strip_after_barcode, strip_before_barcode)
        record2 = apply_plan_to_read(read2, barcodes, llinker,
                 min_length, max_length,
                 strip_after_barcode, strip_before_barcode)
        record_pair = pair(record, record2)
        result = write_record_pair(record_pair,
                              writer_args['barcoded_files'],
                              writer_args['unmatched_file'],
                              writer_args['processed_file'],
                              linker,
                              min_length)
        if result.left == ERR_LINKER or result.right == ERR_LINKER: n_linker += 1
        else:
            if result.left == ERR_TOO_SHORT: n_short += 1
            if result.right == ERR_TOO_SHORT: n_short += 1
    return {'all': i, 'short': n_short, 'linker': n_linker}

cdef Record apply_plan_to_read(tuple t, list barcodes, char *linker,
                 int min_length, int max_length,
                 int strip_after_barcode, int strip_before_barcode):
    cdef:
        Record record
        Read read
    read = as_read(t[0], t[1], t[2])
    if strip_before_barcode > 0:
        read = pretrim_record_5prime(read, strip_before_barcode)
    record = assign_read(read, barcodes)
    if strip_after_barcode > 0:
        record = trim_record_5prime(record, strip_after_barcode)
    if not linker == b'':
        record = cleave_linker(record, linker)
    record = trim_trailing_Ns(record)
    if max_length >= 0:
        record = truncate_record(record, max_length)
    return record
        
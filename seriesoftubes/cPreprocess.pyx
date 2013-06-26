# cython: profile=True, c_string_type=str, c_string_encoding=ascii
from cpython cimport bool
from libc.stdlib cimport malloc, free
from cpython.version cimport PY_MAJOR_VERSION

DEF ERR_TOO_SHORT = 1
DEF ERR_LINKER = 2
DEF MAX_FASTQ_SIZE = 131072


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

#cdef Read asRead(str title, str seq, str qual):
#    cdef Read read
#    read.title = <bytes>title
#    read.seq = <bytes>seq
#    read.qual = <bytes>qual
#    return read

cdef struct Record:
    char *barcode
    Read *read

#cdef Record asRecord(str barcode, str title, str seq, str qual):
#    cdef Record record
#    record.barcode = <bytes>barcode
#    record.read = asRead(title, seq, qual)
#    return record

#cdef struct s_RecordPair:
#    Record *first
#    Record *second
#
#ctypedef s_RecordPair RecordPair
#cdef RecordPair pair(Record record1, Record record2):
#    cdef RecordPair rp
#    rp.first = record1
#    rp.second = record2
#    return rp

cdef bytes match_barcode(bytes seq, list barcodes, int mismatches=1):
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
        bytes accepted = b''
        bytes barcode
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
            else:
                accepted = b''
                break
    return accepted

cdef void pretrim_read_5prime(Read *read, int trim_length):
    '''
    trim an assigned record from the 5' end
    expects (barcode, (title, seq, qual))
    '''
    read[0].seq += trim_length
    read[0].qual += trim_length

cdef void trim_read_5prime(Read *read, int trim_length):
    '''
    trim an assigned record from the 5' end
    expects (barcode, (title, seq, qual))
    '''
    read[0].seq += trim_length
    read[0].qual += trim_length

cdef void truncate_read(Read *read, int max_length):
    '''
    truncate a record so that it is at most max_length
    starting at the 5' end
    expects (barcode, (title, seq, qual))
    '''
    read[0].seq[max_length] = b'\0'
    read[0].qual[max_length] = b'\0'

cdef void trim_read_3prime(Read *read, int trim_length):
    '''
    trim a record from the 3' end
    expects (barcode, (title, seq, qual))
    '''
    if trim_length != 0:
        read[0].seq[-trim_length] = b'\0'
        read[0].qual[-trim_length] = b'\0'

# Violates NCBI SRA requirements# Violates NCBI SRA requirements
#cdef Record trim_trailing_Ns(Record record):
#    '''
#    trim any trailing 3' 'N's
#    expects record is (barcode, (title, seq, qual))
#    returns truncated (barcode, (title, seq, qual))
#    '''
#    cdef int end = len(str.rstrip(record.read.seq, 'N'))
#    record.read.seq[end] = b'\0'
#    record.read.qual[end] = b'\0'
#    return record

cdef void cleave_linker(Record *record, char *linker):
    cdef:
        Read *read = record[0].read
        int i
    i = str.find(read.seq, linker)
    if i != -1:
        record[0].barcode = linker
        read[0].seq[i] = b'\0'
        read[0].qual[i] = b'\0'

cdef Read as_read(bytes title, bytes seq, bytes qual):
    cdef Read read
    read.title = title
    read.seq = seq
    read.qual = qual
    return read

cdef Read *as_read2(tuple t):
    cdef Read read
    read.title = <bytes>t[0]
    read.seq = <bytes>t[1]
    read.qual = <bytes>t[2]

cdef bytes assign_read(Read *read, list barcodes):
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
        bytes barcode
        int barcode_len
        int n_barcodes = len(barcodes)

    if n_barcodes == 0:
        barcode = b''
    else:
        title_head, last_part = read.title.rsplit(':', 1)
        if last_part.isalpha():
            # CASAVA 1.8 file
            barcode = match_barcode(<bytes>last_part.rstrip(), barcodes)
            return barcode
        # MiSeq files
        if read.title.count(' ') == 1:
            title_head2, last_part2 = read.title.split(' ', 1)
            if last_part2.count(':') == 3:
                # regular barcoded read, extract from sequence
                barcode = match_barcode(<bytes>read.seq, barcodes)
                barcode_len = len(barcode)
                read.seq += barcode_len
                read.qual += barcode_len
                return barcode
        # old files
        pound_loc = last_part.find('#')
        slash_loc = last_part[pound_loc:].find('/')
        if slash_loc != -1:
            barcode = <bytes>last_part[(pound_loc+1):(pound_loc+slash_loc)]
        else:
            barcode = <bytes>last_part[(pound_loc+1):]
        if barcode==b'0':
            barcode = match_barcode(<bytes>read.seq, barcodes)
            barcode_len = len(barcode)
            read.seq += barcode_len
            read.qual += barcode_len
        elif barcode.isdigit():
            # then we have a numbered index from Illumina, just use it as-is
            #record.barcode = barcode
            pass
        elif barcode.isalpha():
            # then we already extracted the barcode at some point, try to match it
            barcode = match_barcode(barcode, barcodes)
        else:
            barcode = b''
    return barcode

cdef bytes assign_read_no_clip(Read *read, list barcodes, int offset):
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
        bytes barcode = b''
        int barcode_len
        int n_barcodes = len(barcodes)

    if n_barcodes == 0:
        barcode = b''
    else:
        title_head, last_part = read[0].title.rsplit(':', 1)
        if last_part.isalpha():
            # CASAVA 1.8 file
            barcode = match_barcode(<bytes>last_part.rstrip(), barcodes)
        # MiSeq files
        if read.title.count(' ') == 1:
            title_head2, last_part2 = read.title.split(' ', 1)
            if last_part2.count(':') == 3:
                # regular barcoded read, extract from sequence
                barcode = match_barcode(<bytes>read.seq, barcodes)
                barcode_len = len(barcode)
                return barcode
        else:
            pound_loc = last_part.find('#')
            slash_loc = last_part[pound_loc:].find('/')
            if slash_loc != -1:
                barcode = <bytes>last_part[(pound_loc+1):(pound_loc+slash_loc)]
            else:
                barcode = <bytes>last_part[(pound_loc+1):]
            if barcode == b'0':
                barcode = match_barcode(read.seq + offset, barcodes)
                barcode_len = len(barcode)
            elif barcode.isdigit():
                # then we have a numbered index from Illumina, just use it as-is
                pass
            elif barcode.isalpha():
                # then we already extracted the barcode at some point, try to match it
                barcode = match_barcode(barcode, barcodes)
            else:
                barcode = b''
    return barcode

cdef bool too_short(Read *read, int min_length):
    cdef char c
    cdef int i = 0, L = 0
    for c in bytes(read.seq):
        L += 1
        if c == b'N':
            i += 1
    return L - i < min_length

cdef int write_record(Record *record, barcoded_files=None,
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

cdef int_pair write_record_pair(Record *record, Record *record2,
                                object barcoded_file_pairs=None,
                                object unmatched_files=None,
                                object processed_files=None,
                                object orphaned_read_files=None,
                                object mismatched_files=None,
                                str linker='', int min_length=4):
    """
    write the assigned_record to the correct barcode file
    an assigned_record is a (barcode, (title, seq, qual)) (all strs)

    will not work unless you provide a dictionary of barcoded_file_pairs and
    unmatched_files and mismatched_files (2-tuples of file object)
    """
    cdef:
        int_pair ret

        object barcode, title, seq, qual, line
        int problem = 0, problem2 = 0
        bool is_linker1, is_too_short1

        object barcode2, title2, seq2, qual2, line2
        bool is_linker2, is_too_short2
        bool is_processed = (barcoded_file_pairs is None)

    ret.left = 0
    ret.right = 0
    barcode = record.barcode
    title = record.read.title
    seq = record.read.seq
    qual = record.read.qual
    line = "@%s\n%s\n+%s\n%s\n" % (title, seq, title, qual)
    is_linker1 = len(seq) == 0 and (not linker == '') and barcode == linker
    is_too_short1 = too_short(record.read, min_length)

    barcode2 = record2.barcode
    title2 = record2.read.title
    seq2 = record2.read.seq
    qual2 = record2.read.qual
    line2 = "@%s\n%s\n+%s\n%s\n" % (title2, seq2, title2, qual2)
    is_linker2 = len(seq2) == 0 and (not linker == '') and barcode2 == linker
    is_too_short2 = too_short(record2.read, min_length)

    # check if record 1 is valid
    if is_linker1: problem = ERR_LINKER
    elif is_too_short1: problem = ERR_TOO_SHORT

    # check if record 2 is valid
    if is_linker2: problem2 = ERR_LINKER
    elif is_too_short2: problem2 = ERR_TOO_SHORT
    ret.left = problem
    ret.right = problem2

    # abort write if both reads have problems
    if (not problem == 0) or (not problem2 == 0):
        if not problem2 == 0:
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
                 bool no_clipping=False,
                 logger=None):
    cdef:
        int i=0, n_short=0, n_linker=0
        char *llinker = <bytes>linker
        int result
        Record *record = <Record *>malloc(sizeof(Record))
        Read *read = <Read *>malloc(sizeof(Read))
    record.read = read
    while True:
        i += 1
        try: t = reads.next()
        except StopIteration: break
        read.title = <char *>malloc(len(t[0]))
        read.seq = <char *>malloc(len(t[1]))
        record.barcode = <char *>malloc(len(t[1]))
        read.qual = <char *>malloc(len(t[2]))
        read.title = <bytes>(t[0])
        read.seq = <bytes>(t[1])
        read.qual = <bytes>(t[2])
        apply_plan_to_read(record,
                 read, barcodes, llinker,
                 min_length, max_length,
                 strip_after_barcode, strip_before_barcode,
                 no_clipping)
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
                 bool no_clipping=False, logger=None):
    cdef:
        int i = 0, n_short = 0, n_linker = 0
        char *llinker = <bytes>linker
        char *dummy = ''
        int_pair result
        Record *record = <Record *>malloc(sizeof(Record))
        Record *record2 = <Record *>malloc(sizeof(Record))
        Read *read = <Read *>malloc(sizeof(Read))
        Read *read2 = <Read *>malloc(sizeof(Read))
        tuple t, t2
    record.read = read
    record2.read = read2
    # Allocate ~1 MB of memory each
    # Should be enough for any sequence???
    read.title = <char *>malloc(MAX_FASTQ_SIZE)
    read.seq = <char *>malloc(MAX_FASTQ_SIZE)
    record.barcode = <char *>malloc(MAX_FASTQ_SIZE)
    read.qual = <char *>malloc(MAX_FASTQ_SIZE)
    read2.title = <char *>malloc(MAX_FASTQ_SIZE)
    read2.seq = <char *>malloc(MAX_FASTQ_SIZE)
    record2.barcode = <char *>malloc(MAX_FASTQ_SIZE)
    read2.qual = <char *>malloc(MAX_FASTQ_SIZE)
    while True:
        i += 1
#        read = as_read(t[0], t[1], t[2])
        try:
            t = reads.next()
        except StopIteration:
            try: reads2.next()
            except StopIteration: pass
            else: raise SyntaxWarning('More reads left in second file')
            finally: break

        t2 = reads2.next()

        read.title = <bytes>(t[0])
        read.seq = <bytes>(t[1])
        read.qual = <bytes>(t[2])
        read2.title = <bytes>(t2[0])
        read2.seq = <bytes>(t2[1])
        read2.qual = <bytes>(t2[2])

        apply_plan_to_read(record,
                 read, barcodes, llinker,
                 min_length, max_length,
                 strip_after_barcode, strip_before_barcode,
                 no_clipping)
        apply_plan_to_read(record2,
                 read2, barcodes, llinker,
                 min_length, max_length,
                 strip_after_barcode, strip_before_barcode,
                 no_clipping)

        result = write_record_pair(record, record2,
                              writer_args['barcoded_file_pairs'],
                              writer_args['unmatched_files'],
                              writer_args['processed_files'],
                              writer_args['orphaned_read_files'],
                              writer_args['mismatched_files'],
                              linker,
                              min_length)
        if result.left == ERR_LINKER or result.right == ERR_LINKER:
            n_linker += 1
        else:
            if result.left == ERR_TOO_SHORT: n_short += 1
            if result.right == ERR_TOO_SHORT: n_short += 1
    return {'all': i, 'short': n_short, 'linker': n_linker}

cdef void apply_plan_to_read(Record *record,
                 Read *read, list barcodes, char *linker,
                 int min_length, int max_length,
                 int strip_after_barcode, int strip_before_barcode,
                 bool no_clipping):
    cdef bytes barcode
    if no_clipping:
        barcode = assign_read_no_clip(read, barcodes,
                                             strip_before_barcode)
    else:
        if strip_before_barcode > 0:
            pretrim_read_5prime(read, strip_before_barcode)
        barcode = assign_read(read, barcodes)
        if strip_after_barcode > 0: trim_read_5prime(read, strip_after_barcode)
        if not linker == b'': cleave_linker(record, linker)
#    record = trim_trailing_Ns(record)
        if max_length >= 0: truncate_read(read, max_length)
    record.barcode = barcode


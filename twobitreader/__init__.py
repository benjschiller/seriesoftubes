import array
import bisect
import ctypes
import os
import os.path
import struct
import sys
from itertools import chain, islice, imap, izip

def true_long_type():
    """
    OS X uses an 8-byte long, so make sure this long is the right size
    and switch to Int if needed
    """
    long_size = ctypes.sizeof(ctypes.c_long)
    if long_size == 4:
        return 'L'
    elif long_size == 8:
        if ctypes.sizeof(ctypes.c_int)==4: return 'I'
    raise ImportError("Couldn't find a valid 4-byte long type equivalent to LONG")
         

LONG = true_long_type()

def byte_to_bases(x):
    c = (x >> 4) & 0xf
    f = x & 0xf
    cc = (c >> 2) & 0x3
    cf = c & 0x3
    fc = (f >> 2) & 0x3
    ff = f & 0x3
    return map(bits_to_base, (cc, cf, fc, ff))

def bits_to_base(x):
    if x is 0: return 'T'
    if x is 1: return 'C'
    if x is 2: return 'A'
    if x is 3: return 'G'

def base_to_bin(x):
    if x == 'T': return '00'
    if x == 'C': return '01'
    if x == 'A': return '10'
    if x == 'G': return '11'

def create_byte_table():
    """create BYTE_TABLE"""
    d = {}
    for x in xrange(256):
        d[x] = byte_to_bases(x)
    return d

def split16(x):
    """
    split a 16-bit number into integer representation
    of its course and fine parts in binary representation
    """
    c = (x >> 8) & 0xff
    f = x & 0xff
    return c, f

def create_twobyte_table():
    """create TWOBYTE_TABLE"""
    d = {}
    for x in xrange(66536):
        c, f = split16(x)
        d[x] = byte_to_bases(c) + byte_to_bases(f)
    return d

BYTE_TABLE = create_byte_table()
TWOBYTE_TABLE = create_twobyte_table()

def longs_to_char_array(longs, first_base_offset, last_base_offset, array_size):
    """
    takes in a iterable of longs and converts them to bases in a char array
    returns a ctypes string buffer
    """
    longs_len = len(longs)
    dna = ctypes.create_string_buffer(array_size)
    # translate from 32-bit blocks to bytes
    # this method ensures correct endianess
    bytes = array.array('B')
    bytes.fromstring(longs.tostring())
    # first block
    first_block = ''.join([''.join(BYTE_TABLE[bytes[x]]) for x in range(4)])
    i = 16 - first_base_offset
    if array_size < i: i = array_size
    dna[0:i] = first_block[first_base_offset:first_base_offset + i]
    if longs_len == 1: return dna
    # middle blocks (implicitly skipped if they don't exist)
    for byte in bytes[4:-4]:
        dna[i:i + 4] = BYTE_TABLE[byte]
        i += 4
    # last block
    last_block = ''.join([''.join(BYTE_TABLE[bytes[x]]) for x in range(-4,0)])
    dna[i:i + last_base_offset] = last_block[0:last_base_offset]
    return dna

class TwoBitFile(dict):
    """
python-level reader for .2bit files (i.e., from UCSC genome browser)
(note: no writing support)


TwoBitFile inherits from dict
You may access sequences by name, e.g.
>>> genome = TwoBitFile('hg18.2bit')
>>> chr20 = genome['chr20']

Sequences are returned as TwoBitSequence objects
You may access intervals by slicing or using str() to dump the entire entry
e.g.
>>> chr20[100100:100200]
'ttttcctctaagataatttttgccttaaatactattttgttcaatactaagaagtaagataacttccttttgttggtat
ttgcatgttaagtttttttcc'
>>> whole_chr20 = str(chr20)

Fair warning: dumping the entire chromosome requires a lot of memory

See TwoBitSequence for more info
    """
    
    def __init__(self, foo):
        super(dict, self).__init__(self)
        if not os.path.exists(foo):
            raise TwoBitReadError('{!s} not found'.format(foo))
        if not os.access(foo, os.R_OK):
            raise TwoBitReadError('Cannot open {!s} for reading'.format(foo))
        self._filename = foo
        self._file_handle = open(foo, 'rb')
        self._load_header()
        self._load_index()
        for name, offset in self._offset_dict.iteritems():
            self[name] = TwoBitSequence(self._file_handle, offset,
                                        self._byteswapped)
        return        
        
    def _load_header(self):
        file_handle = self._file_handle
        #
        # load header
        #
#        try: raw_header = file_handle.read(16)
#        except IOError: raise TwoBitReadError('Premature file end')
        header = array.array(LONG)
        try: header.fromfile(file_handle, 4)
        except IOError: raise TwoBitReadError('Premature file end')
        # check signature -- must be 0x1A412743
        # if not, swap bytes
        byteswapped = False
#        endianess = '='
#        (signature, version, sequence_count, reserved) = struct.unpack("=IIII", raw_header)
        (signature, version, sequence_count, reserved) = header
        if not signature == 0x1A412743:
            byteswapped = True
#            if sys.byteorder=='little': endianess = '>'
#            else: endianess = '<'
#            swapped_header = array.array(LONG, [signature, version,
#                                               sequence_count, reserved])
#            swapped_header.byteswap()
            header.byteswap()
#            (signature, version, sequence_count, reserved) = swapped_header
            (signature, version, sequence_count, reserved) = header
            if not signature == 0x1A412743:
                raise TwoBitReadError('Invalid 2-bit file signature in header')
        if not version == 0: 
            raise TwoBitReadError('Invalid 2-bit file version in header')
        if not reserved == 0: 
            raise TwoBitReadError('Invalid 2-bit file reserved header field')
#        self._endianess = endianess
        self._byteswapped = byteswapped
        self._sequence_count = sequence_count
        
    def _load_index(self):
        #
        # load index
        #
        file_handle = self._file_handle
#        endianess = self._endianess
        byteswapped = self._byteswapped
        remaining = self._sequence_count
        sequence_offsets = []
        file_handle.seek(16)
        while True:
            if remaining == 0: break
            name_size = array.array('B')
            try: name_size.fromfile(file_handle, 1)
            except IOError: raise TwoBitReadError('Premature file end')
            if byteswapped: name_size.byteswap()
            name = array.array('c')
            if byteswapped: name.byteswap()
            try: name.fromfile(file_handle, name_size[0])
            except IOError: raise TwoBitReadError('Premature file end')
#            try: raw_index = name.fromfile(file_handle, name_size[0])
#            except IOError: raise TwoBitReadError('Premature file end')
            offset = array.array(LONG)
            try: offset.fromfile(file_handle, 1)
            except IOError: raise TwoBitReadError('Premature file end')
            if byteswapped: offset.byteswap()
            sequence_offsets.append((name.tostring(), offset[0]))
#            try: raw_index = file_handle.read(name_size + 4)
#            except IOError: raise TwoBitReadError('Premature file end')
            # index is (name, offset)
#            index = struct.unpack(endianess + 'c'*name_size + LONG, raw_index)
#            sequence_offset = (''.join(index[0:-1]), index[-1])
#            sequence_offsets.append(sequence_offset)
            remaining -= 1
        self._sequence_offsets = sequence_offsets
        self._offset_dict = dict(sequence_offsets)

    def sequence_sizes(self):
        """returns a dictionary with the sizes of each sequence"""
        d = {}
        file_handle = self._file_handle
        endianess = self._endianess
        for name, offset in self._offset_dict.iteritems():
            file_handle.seek(offset)
            raw_record_header = file_handle.read(4)
            dna_size = struct.unpack(endianess + 'I', raw_record_header)[0]
            d[name] = dna_size
        return d

class TwoBitSequence(object):
    """
A TwoBitSequence object refers to an entry in a TwoBitFile

You may access intervals by slicing or using str() to dump the entire entry
e.g.
>>> genome = TwoBitFile('hg18.2bit')
>>> chr20 = genome['chr20']
>>> chr20[100100:100200] # slicing returns a string
'ttttcctctaagataatttttgccttaaatactattttgttcaatactaagaagtaagataacttccttttgttggtat
ttgcatgttaagtttttttcc'
>>> whole_chr20 = str(chr20) # get whole chr as string

Fair warning: dumping the entire chromosome requires a lot of memory
 
Note that we follow python/UCSC conventions:
Coordinates are 0-based, end-open
If you attempt to access a slice past the end of the sequence,
it will be truncated at the end.

Your computer probably doesn't have enough memory to load a whole genome
but if you want to string-ize your TwoBitFile, here's a recipe:

x = TwoBitFile('my.2bit')
d = x.dict()
for k,v in d.iteritems(): d[k] = str(v)
    """
    def __init__(self, file_handle, offset, byteswapped=False):
        self._file_handle = file_handle
        self._original_offset = offset
        self._byteswapped = byteswapped
#        if byteswapped:
#            if sys.byteorder=='little': endianess = '>'
#            else: endianess = '<'
#        else: endianess = '='
#        self._endianess = endianess
        file_handle.seek(offset)
        header = array.array(LONG)
        header.fromfile(file_handle, 2)
        if byteswapped: header.byteswap()
        dna_size, n_block_count = header
#        dna_size, n_block_count = struct.unpack(endianess + 'LL',
#                                                raw_record_header)
        self._dna_size = dna_size
        self._packed_dna_size = (dna_size + 15) / 16 # this is 32-bit fragments
#        raw_n_info = file_handle.read(8*n_block_count)
        n_block_starts = array.array(LONG)
        n_block_sizes = array.array(LONG)
        n_block_starts.fromfile(file_handle, n_block_count)
        if byteswapped: n_block_starts.byteswap()
        n_block_sizes.fromfile(file_handle, n_block_count)
        if byteswapped: n_block_sizes.byteswap()
#        n_info = struct.unpack(endianess + 'LL'*n_block_count, raw_n_info)
        self._n_block_starts = n_block_starts
        self._n_block_sizes= n_block_sizes
        mask_rawc = array.array(LONG)
        mask_rawc.fromfile(file_handle, 1)
        if byteswapped: mask_rawc.byteswap()
        mask_block_count = mask_rawc[0]
#        mask_block_count = struct.unpack(endianess + LONG,
#                                         file_handle.read(4))[0]
        mask_block_starts = array.array(LONG)
        mask_block_starts.fromfile(file_handle, mask_block_count)
        if byteswapped: mask_block_starts.byteswap()
        mask_block_sizes = array.array(LONG)
        mask_block_sizes.fromfile(file_handle, mask_block_count)
        if byteswapped: mask_block_sizes.byteswap()
        self._mask_block_starts = mask_block_starts
        self._mask_block_sizes = mask_block_sizes
        file_handle.read(4)
        self._offset = file_handle.tell()

    def __len__(self):
        return self._dna_size

    def __getslice__(self, min, max=None):
        return self.get_slice(min, max)

    def get_slice(self, min, max=None):
        """
        get_slice returns only a sub-sequence
        """
        # handle negative coordinates
        dna_size = self._dna_size
        if max < 0:
            if max < -dna_size: raise IndexError('index out of range')
            max = dna_size + 1 + max
        if min < 0:
            if max < -dna_size: raise IndexError('index out of range')
            min = dna_size + 1 + min
        # make sure there's a proper range
        if min > max and max is not None: return ''
        if max == 0: return ''
        # load all the data
        if max > dna_size: max = dna_size
        file_handle = self._file_handle
        byteswapped = self._byteswapped
        n_block_starts = self._n_block_starts
        n_block_sizes = self._n_block_sizes
        mask_block_starts = self._mask_block_starts
        mask_block_sizes = self._mask_block_sizes
        offset = self._offset
        packed_dna_size = self._packed_dna_size

        # region_size is how many bases the region is       
        if max is None: region_size = dna_size - min
        else: region_size = max - min
        
        # first_block, last_block are the first/last 32-bit blocks we need
        # blocks start at 0
        first_block = min / 16
        last_block = max / 16
        # don't read past seq end
        if last_block >= packed_dna_size: last_block = packed_dna_size -1
        # +1 we still need to read block
        blocks_to_read = last_block - first_block + 1
        
        # jump directly to desired file location
        local_offset = offset + first_block * 4
        file_handle.seek(local_offset)
        
        # note we won't actually read the last base
        # this is a python slice first_base_offset:16*blocks+last_base_offset
        first_base_offset = min % 16
        last_base_offset = max % 16
        
        fourbyte_dna = array.array(LONG)
        fourbyte_dna.fromfile(file_handle, blocks_to_read)
        if byteswapped: fourbyte_dna.byteswap()
        string_as_array = longs_to_char_array(fourbyte_dna, first_base_offset,
                                              last_base_offset, region_size)
        for start, size in izip(n_block_starts, n_block_sizes):
            end = start + size
            if end <= min: continue
            if start > max: break
            if start < min: start = min
            if end > max: end = max 
            start -= min
            end -= min
            string_as_array[start:end] = 'N'*(end-start)
        lower = str.lower
        first_useful = bisect.bisect_right(mask_block_starts, min) - 1
        if first_useful == -1: first_useful = 0
        last_useful = 1 + bisect.bisect_right(mask_block_starts, max,
                                             lo=first_useful)
        for start, size in izip(mask_block_starts[first_useful:last_useful],
                                mask_block_sizes[first_useful:last_useful]):
            end = start + size
            if end <= min: continue
            if start > max: break
            if start < min: start = min
            if end > max: end = max 
            start -= min
            end -= min
            string_as_array[start:end] = lower(string_as_array[start:end])
        return string_as_array[0:]

    def __str__(self):
        """
        returns the entire chromosome
        """
        return self.__getslice__(0, None)
    
class TwoBitReadError(IOError):
    """
    Base exception for TwoBit module
    """
    def __init__(self, *args):
        self.msg = ''.join(args)
    def __str__(self):
        return self.msg

def print_specification():
    """
    Prints the twoBit file format specification I got from the Internet.
    This is only here for reference
    """
    return """
From http://www.its.caltech.edu/~alok/reviews/blatSpecs.html

.2bit files

A .2bit file can store multiple DNA sequence (up to 4 gig total) in a compact randomly accessible format. The two bit files contain masking information as well as the DNA itself. The file begins with a 16 byte header containing the following fields:

signature - the number 0x1A412743 in the architecture of the machine that created the file.
version - zero for now. Readers should abort if they see a version number higher than 0.
sequenceCount - the number of sequences in the file
reserved - always zero for now.
All fields are 32 bits unless noted. If the signature value is not as given, the reader program should byte swap the signature and see if the swapped version matches. If so all multiple-byte entities in the file will need to be byte-swapped. This enables these binary files to be used unchanged on different architectures.

The header is followed by a file index. There is one entry in the index for each sequence. Each index entry contains three fields:

nameSize - a byte containing the length of the name field
name - this contains the sequence name itself, and is variable length depending on nameSize.
offset - 32 bit offset of the sequence data relative to the start of the file
The index is followed by the sequence records. These contain 9 fields:

dnaSize - number of bases of DNA in the sequence.
nBlockCount - the number of blocks of N's in the file (representing unknown sequence).
nBlockStarts - a starting position for each block of N's
nBlockSizes - the size of each block of N's
maskBlockCount - the number of masked (lower case) blocks
maskBlockStarts - starting position for each masked block
maskBlockSizes - the size of each masked block
packedDna - the dna packed to two bits per base as so: 00 - T, 01 - C, 10 - A, 11 - G. The first base is in the most significant 2 bits byte, and the last base in the least significant 2 bits, so that the sequence TCAG would be represented as 00011011. The packedDna field will be padded with 0 bits as necessary so that it takes an even multiple of 32 bit in the file, as this improves i/o performance on some machines.
.nib files
"""
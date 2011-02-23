import struct
import os
import os.path
from itertools import chain, islice, imap

def create_table():
    trans = {'00': 'T', '01': 'C', '10': 'A', '11': 'G'}
    four_trans = {}
    for x in range(256):
        y= bin(x)[2:]
        y = '0'*(8-len(y)) + y
        one, two, three, four = y[0:2], y[2:4], y[4:6], y[6:8]
        bases = trans[one] + trans[two] + trans[three] + trans[four]
        four_trans[x] = bases
    return four_trans
    
FOUR_TABLE = create_table()

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
        self._load_header_and_index()
        for k,v in self._offset_dict.iteritems():
            self[k] = TwoBitSequence(self._file_handle, v, self._endianess)
        return        
        
    def _load_header_and_index(self):
        too_short_error = TwoBitReadError('Premature file end')
        file_handle = self._file_handle
        #
        # load header
        #
        try: raw_header = file_handle.read(16)
        except IOError: raise too_short_error
        # check signature -- must be 0x1A412743
        # try little endian first
        endianess = '<'
        signature = struct.unpack(endianess + "I", raw_header[0:4])[0]
        if not signature == 0x1A412743:
            # try big endian next
            self._endianness = '>'
            signature = struct.unpack(endianess + "I", raw_header[0:4])[0]
            if not signature == 0x1A412743:
                raise TwoBitReadError('Invalid 2-bit file signature in header')
        self._endianess = endianess
        #else: pass
        # load the remaining header entries
        (version, sequence_count,
         reserved) = struct.unpack(self._endianess + "III", raw_header[4:])
        if not version == 0: 
            raise TwoBitReadError('Invalid 2-bit file version in header')
        if not reserved == 0: 
            raise TwoBitReadError('Invalid 2-bit file reserved header field')
        self._sequence_count = sequence_count
        #
        # load index
        #
        remaining = sequence_count
        sequence_offsets = []
        while True:
            if remaining == 0: break
            try: raw_name_size = file_handle.read(1)
            except IOError: raise too_short_error
            name_size = struct.unpack(endianess + "B", raw_name_size)[0]
            try: raw_index = file_handle.read(name_size + 4)
            except IOError: raise too_short_error
            # index is (name, offset)
            index = struct.unpack(endianess + 'c'*name_size + 'I', raw_index)
            sequence_offset = (''.join(index[0:-1]), index[-1])
            sequence_offsets.append(sequence_offset)
            remaining -= 1
        self._sequence_offsets = sequence_offsets
        self._offset_dict = dict(sequence_offsets)
        self._list = [item[0] for item in sequence_offsets] 
        # index loaded, done now   

    def sequence_sizes(self):
        """returns a dictionary with the sizes of each sequence"""
        d = {}
        for name, offset in self._offset_dict.iteritems():
            file_handle = self._file_handle
            endianess = self._endianess
            file_handle.seek(offset)
            raw_record_header =  file_handle.read(4)
            dna_size = struct.unpack(endianess + 'I',
                                     raw_record_header)[0]
            d[name] = dna_size
        return d

#class _TwoBitIter(list):

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
    def __init__(self, file_handle, offset, endianess):
        self._file_handle = file_handle
        self._original_offset = offset
        self._endianess = endianess
        file_handle.seek(offset)
        raw_record_header =  file_handle.read(8)
        dna_size, n_block_count = struct.unpack(endianess + 'II',
                                                raw_record_header)
        self._dna_size = dna_size
        self._padding = (dna_size%4) / 2
        self._packed_dna_size = (dna_size + 2) / 4
        self._n_block_count = n_block_count
        raw_n_info = file_handle.read(8*n_block_count)
        n_info = struct.unpack(endianess + 'II'*n_block_count, raw_n_info)
        self._n_block_starts = n_info[0:n_block_count]
        self._n_block_sizes= n_info[n_block_count:]
        mask_block_count = struct.unpack(endianess + 'I',
                                         file_handle.read(4))[0]
        self._mask_block_count = mask_block_count
        raw_mask_info = file_handle.read(8 * mask_block_count)
        mask_info = struct.unpack(endianess + 'II'*mask_block_count,
                                  raw_mask_info)
        self._mask_block_starts = mask_info[0:mask_block_count]
        self._mask_block_sizes= mask_info[mask_block_count:]
        # skip reserved 32-bit word
        # note sure why this is three instead of 4
        file_handle.read(4)
        self._offset = file_handle.tell()
        # continue

    def __len__(self):
        return self._dna_size

    def __getslice__(self, min, max=None):
        """
        get_range returns only a sub-sequence
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
        packed_dna_size = self._packed_dna_size
        padding = self._padding
        file_handle = self._file_handle
        endianess = self._endianess
        n_block_count = self._n_block_count
        n_block_starts = self._n_block_starts
        n_block_sizes = self._n_block_sizes
        mask_block_count = self._mask_block_count
        mask_block_starts = self._mask_block_starts
        mask_block_sizes = self._mask_block_sizes
        offset = self._offset
       
        if max is None: region_size = dna_size - min
        else: region_size = max - min
        read_size = (region_size + 2) / 4
        read_offset = offset + min / 4 + (packed_dna_size%4) 
        if read_size + padding <= packed_dna_size:
            read_size += padding
        left_padding = min % 4
        if max is not None: right_padding = max % 4
        else: right_padding = 0
        if left_padding > 0:
            if read_size < packed_dna_size:
                read_size += 1
        if right_padding > 0:
            if read_size < packed_dna_size:
                read_size += 1
        
        # do work
        file_handle.seek(read_offset)
        packed_dna = file_handle.read(read_size)
        byte_dna = struct.unpack(endianess + 'B'*read_size, packed_dna)
        if max is None: stop = None
        else: stop = max - min + left_padding + padding
        string_as_list = list(islice(chain(*[FOUR_TABLE[byte] for byte in byte_dna]),
                                     left_padding + padding, stop))
        for i in xrange(n_block_count):
            start = n_block_starts[i]
            end = start + n_block_sizes[i]
            if start > max or end == max: break
            if start < min: start = min
            start -= min
            if end > max: end = max 
            end -= min
            if end < 1: continue
            string_as_list[start:end] = ['N']*(end-start)
        lower = str.lower
        for i in xrange(mask_block_count):
            start = mask_block_starts[i]
            end = start + mask_block_sizes[i]
            if start > max or end == max: break
            if start < min: start = min
            start -= min
            if end > max: end = max 
            end -= min
            if end < 1: continue
            string_as_list[start:end] = [lower(x) for x in \
                                                    string_as_list[start:end]]
        dna_seq = ''.join(string_as_list)
        return dna_seq

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
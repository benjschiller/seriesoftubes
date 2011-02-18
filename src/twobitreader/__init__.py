import struct
import os
import os.path

class TwoBitReader(object):
    """
    python-level reader for .2bit files (i.e., from UCSC genome browser)
    
    TwoBitReader.dict() returns a dictionary for access to sequences by name
    TwoBitReader.sequences() returns a generator that yields sequences
    """
    def __init__(self, foo):
        too_short_error = TwoBitReaderError('Premature file end')
        
        if not os.path.exists(foo):
            raise TwoBitReaderError('{!s} not found'.format(foo))
        if not os.access(foo, os.R_OK):
            raise TwoBitReaderError('Cannot open {!s} for reading'.format(foo))
        self._filename = foo
        file_handle = open(foo, 'rb')
        self._file_handle = file_handle
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
                raise TwoBitReaderError('Invalid 2-bit file signature in header')
        self._endianess = endianess
        #else: pass
        # load the remaining header entries
        (version, sequence_count,
         reserved) = struct.unpack(self._endianess + "III", raw_header[4:])
        if not version == 0: 
            raise TwoBitReaderError('Invalid 2-bit file version in header')
        if not reserved == 0: 
            raise TwoBitReaderError('Invalid 2-bit file reserved header field')
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
            print sequence_offset
            sequence_offsets.append(sequence_offset)
            remaining -= 1
        self._sequence_offsets = sequence_offsets
        self._dict = dict(sequence_offsets)
        self._list = [item[0] for item in sequence_offsets] 
        # index loaded, done now   

    def dict(self):
        return _TwoBitDict(self._file_handle, self._dict)

class _TwoBitDict(dict):
    def __init__(self, file_handle, offset_dict):
        self._offset_dict = offset_dict
        
    def __getitem__(self, i):
        self._file_handle.seek(self._offset_dict, i)
    
    def __setitem__(self, i, value):
        raise NotImplementedError

class TwoBitReaderError(IOError):
    """
    Base exception for TwoBitReader class
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
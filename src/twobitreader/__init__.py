import struct
import os
import os.path

class TwoBitReader(object):
    def __init__(self, foo):
        if not os.path.exists(foo):
            raise IOError('{!s} not found'.format(foo))
        if not os.access(foo, os.R_OK):
            raise IOError('Cannot open {!s} for reading'.format(foo))
        self._filename = foo
        self._file = open(foo, 'rb')
        # load header
        raw_header = self._file.read(16)
        self.header = struct.unpack_from("iiii", raw_header)
        print self.header
        (signature, version, sequence_count, reserved) = self.header
        # load index

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
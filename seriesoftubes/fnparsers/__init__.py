"""
implements FilenameParsers for scripter
"""

from bowtie import BowtieFilenameParser
#from bwa import bwaFilenameParser
from illumina import BarcodeFilenameParser
from aligned import BAMFilenameParser
from peaks import PeaksFilenameParser
__all__ = [
    'BowtieFilenameParser',
#    'bwaFilenameParser',
    'BarcodeFilenameParser',
    'BAMFilenameParser',
    'PeaksFilenameParser'
]

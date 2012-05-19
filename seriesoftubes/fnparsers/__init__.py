"""
implements FilenameParsers for scripter
"""

from bowtie import BowtieFilenameParser
from illumina import BarcodeFilenameParser
from aligned import BAMFilenameParser
from peaks import PeaksFilenameParser
__all__ = ['BowtieFilenameParser', 'BarcodeFilenameParser',
           'BAMFilenameParser', 'PeaksFilenameParser']
"""
implements FilenameParsers for scripter
"""

from bowtie import BowtieFilenameParser
from illumina import BarcodeFilenameParser
from aligned import AlignmentsFilenameParser
__all__ = ['BowtieFilenameParser', 'BarcodeFilenameParser',
           'AlignmentsFilenameParser']
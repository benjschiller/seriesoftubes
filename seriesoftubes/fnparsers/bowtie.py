import pysam
import os
import os.path
import gzip
import bz2
import scripter

class BowtieFilenameParser(scripter.FilenameParser):
    def __init__(self, filename, *args, **kwargs):
        super(BowtieFilenameParser, self).__init__(filename, *args, **kwargs)
        compression, format = self.discover_file_format2(filename)
        if format == 'SAM' or format == 'BAM':
            self.use_pysam = True
            # try to open the file so we're sure it works
            f = pysam.Samfile(filename)
            del f
        elif format == 'FASTQ':
            self.use_pysam = False
            self.compression = compression
        else:
            raise RuntimeError('Dubious file format')
        self.second_file = None
        self.check_paired_end()
        if len(self.protoname.split('.')) > 6:
            self.fastq_source = self.protoname.split('.')[6]
        else:
            self.fastq_source = 'Unknown'
    
    @staticmethod
    def discover_file_format2(filename):
        f = open(filename, 'rb')
        head = f.read(3)
        f.close()
        # check magic words for compression
        if head == '\x1f\x8b\x08':
            compression = 'gz'
            open_func = gzip.GzipFile
        elif head=='\x42\x5a\x68':
            open_func = bz2.BZ2File
            compression = 'bz2'
        else:
            open_func = open
            compression = 'none'
        uncompressed = open_func(filename)
        head2 = uncompressed.read(4)
        # check for BAM
        if head2 == 'BAM\x01': return (compression, 'BAM')
        # check for SAM 
        head2b = uncompressed.readline()
        if head2 == '@HD\t':
            if head2b[0:3] == 'VN:': return (compression, 'SAM')
        # check FASTQ
        title = head2 + head2b
        seq = uncompressed.readline()
        title2 = uncompressed.readline()
        qual = uncompressed.readline()
        if len(seq) == len(qual) and \
           (title.startswith('@') or title2.startswith('+')) and \
           (title2[1:].strip() == '' or title[1:] == title2[1:]):
            return (compression, 'FASTQ')
        # otherwise give up
        else: return (None, None)
            
    def check_paired_end(self):
        # check if this is a paired-end file
        # if so, grab its partner
        seqfile_name = os.path.basename(self.input_file)
        pair_info = get_pair_info(seqfile_name)
        if pair_info is not None:
            pair_index = pair_info[0]
            second_name = pair_info[1]
            new_name = pair_info[2]
            scripter.debug('NOTICE: Detected paired read file.')
            if pair_index == '1':
                scripter.debug('Attempting to find second file.')

                self.second_file = os.sep.join([self.input_dir, second_name])
                self.protoname = os.path.splitext(new_name)[0]
                scripter.debug('Found %s', self.second_file)
                try:
                    scripter.assert_path(self.second_file)
                    self.paired_end = True
                except IOError:
                    scripter.debug('Failed to find paired end file')
                    self.paired_end = False
            elif pair_index == '2':
                scripter.debug('This is the second file, ignoring it.')
                raise scripter.InvalidFileException
            else:
                scripter.debug('Failed to find paired end')
                self.paired_end = False
        else:
            scripter.debug('This file contains single-end reads.')
            self.paired_end = False

    def tmp_filename(self, ref, match_type):
        return os.path.join(self.output_dir, ref, match_type,
                                       self.with_extension('tmp'))

def get_pair_info(illumina_name):
    """
    take a filename from GERALD output and figure out whether it's the first
    or second read (or single-end read) and return a tuple
    (pair_index, second_file_name, new_output_name)
    """
    name_parts = illumina_name.split('.')
    for i in range(len(name_parts)):
        part = name_parts[i]
        subparts = part.split('_') 
        if len(subparts) > 0:
            if subparts[0] == 's':
                if len(subparts) == 1: continue
                pair_index = subparts[2]
                # lane = subparts[1]
                if pair_index == '1':
                    join_part = '_'.join(subparts[0:2] + subparts[3:])
                    new_output_name = '.'.join(name_parts[0:i] + [join_part] +\
                                                name_parts[i+1:])
                    second_part = '_'.join(subparts[0:2] + ['2'] + subparts[3:])
                    second_name = '.'.join(name_parts[0:i] + [second_part] + 
                                           name_parts[i+1:])
                    return (pair_index, second_name, new_output_name)
                elif pair_index == '2':
                    raise scripter.InvalidFileException
                    return
    return None
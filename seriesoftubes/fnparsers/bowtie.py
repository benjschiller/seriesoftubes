import pysam
import os
import os.path
import gzip
import bz2
import scripter
from ..converters.discover import discover_file_format

class BowtieFilenameParser(scripter.FilenameParser):
    def __init__(self, filename, *args, **kwargs):
        super(BowtieFilenameParser, self).__init__(filename, *args, **kwargs)
        open_func, format = discover_file_format(filename)
#        self.split_file = False
        self.format = format
        self.open_func = open_func
        self.second_file = None
        if format == 'SAM' or format == 'BAM':
            self.use_pysam = True
            # try to open the file so we're sure it works
            f = pysam.Samfile(filename)
            aread = f.next()
            self.paired_end = aread.is_paired
            del f, aread
            self.fastq_source = 'Unknown'
        elif format == 'FASTQ':
            self.use_pysam = False
            self.check_paired_end()
            if len(self.protoname.split('.')) > 6:
                self.fastq_source = self.protoname.split('.')[6]
            else:
                self.fastq_source = 'Unknown'
        else:
            raise RuntimeError('Dubious file format')
    
    def check_paired_end(self):
        # check if this is a paired-end file
        # if so, grab its partner
        seqfile_name = os.path.basename(self.input_file)
        pair_info = get_pair_info(seqfile_name)
        if pair_info is None:
            pair_info = get_new_pair_info(seqfile_name)
#            if pair_info is not None: self.split_file = True
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

    def tmp_filename(self, ref, match_type=None):
        if match_type is None:
            return os.path.join(self.output_dir, ref,
                                self.with_extension('tmp'))
        else:
            return os.path.join(self.output_dir, ref, match_type,
                                self.with_extension('tmp'))

def get_new_pair_info(illumina_name):
    """
    take a filename from CASAVA 1.8 output and figure out whether it's the
    first or second read (or single-end read) and return a tuple
    (pair_index, second_file_name, new_output_name).
    """
    end = illumina_name.find('.fastq')
    if end == -1: return None
    name = illumina_name[0:end]
    parts = name.split('_')
    if len(parts) < 3: return None
    # align split files one at a time
#    if not parts[-1].split('.')[0] == '001':
#        scripter.debug('Not the first of a split file, ignoring it')
#        raise scripter.InvalidFileException
    read = parts[-2]
    if read == 'R1' or read =='R2':
        num = parts[-1]
        second_file = '_'.join(parts[0:-2]) + '_R2_' + num + illumina_name[end:]
        output_name = '_'.join(parts[0:-2]) + '_' + num + '.fastq'
        return (read[1], second_file, output_name)
    return None

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
                join_part = '_'.join(subparts[0:2] + subparts[3:])
                new_output_name = '.'.join(name_parts[0:i] + [join_part] +\
                                            name_parts[i+1:])
                second_part = '_'.join(subparts[0:2] + ['2'] + subparts[3:])
                second_name = '.'.join(name_parts[0:i] + [second_part] + 
                                       name_parts[i+1:])
                return (pair_index, second_name, new_output_name)
    return None
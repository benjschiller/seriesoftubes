import scripter
import os


def get_new_pair_info(illumina_name):
    """
    take a filename from CASAVA 1.8 output and figure out whether it's the
    first or second read (or single-end read) and return a tuple
    (pair_index, second_file_name, new_output_name).
    """
    end = illumina_name.find('.fastq')
    if end == -1:
        return None
    name = illumina_name[0:end]
    parts = name.split('_')
    if len(parts) < 3:
        return None
    read = parts[-2]
    if read == 'R1' or read == 'R2':
        num = parts[-1]
        second_file = '_'.join(parts[0:-2]) + '_R2_' + num + \
                      illumina_name[end:]
        return (read, second_file)
    return None


class BarcodeFilenameParser(scripter.FilenameParser):
    def __init__(self, filename, verbose=False, *args, **kwargs):
        super(BarcodeFilenameParser, self).__init__(filename,
                                                    *args, **kwargs)
        protoname = self.protoname
        # check for old-style
        if os.path.splitext(protoname)[-3:] == 'all':
            protoname = protoname[0:-4]

        # check if this is a paired-end file
        # if so, grab its partner
        input_file = self.input_file
        illumina_name = os.path.basename(input_file)

        # try new style first
        new_info = get_new_pair_info(illumina_name)
        if new_info is not None:
            scripter.debug('NOTICE: Detected new-style paired read file.')
            read = new_info[0]
            if read == 'R2':
                scripter.debug('This is the second file, ignoring it.')
                raise scripter.InvalidFileException(input_file)
            elif read == 'R1':
                second_file = os.path.join(self.input_dir, new_info[1])
                try:
                    scripter.assert_path(second_file)
                    scripter.debug('Found %s', second_file)
                    self.second_file = second_file
                    self.protoname2 = os.path.splitext(
                        os.path.basename(second_file))[0]
                    paired_end = True
                except IOError:
                    scripter.debug('Failed to find paired end file')
                    paired_end = False
            else:
                scripter.debug('Failed to find paired end')
                paired_end = False
        elif illumina_name.count('_') >= 3:
            scripter.debug('NOTICE: Detected paired read file.')
            iln_parts = illumina_name.split('_')
            if iln_parts[2] == '1':
                scripter.debug('Attempting to find second file.')

                second_file = os.sep.join([self.input_dir,
                                           '_'.join(iln_parts[0:2] + ['2']
                                                    + iln_parts[3:])])
                try:
                    scripter.assert_path(second_file)
                    scripter.debug('Found %s', second_file)
                    self.second_file = second_file
                    self.protoname2 = os.path.splitext(
                        os.path.basename(second_file))[0]
                    paired_end = True
                except IOError:
                    scripter.debug('Failed to find paired end file')
                    paired_end = False
            elif iln_parts[2] == '2':
                scripter.debug('This is the second file, ignoring it.')
                raise scripter.InvalidFileException(input_file)
            else:
                scripter.debug('Failed to find paired end')
                paired_end = False
        else:
            paired_end = False
        self.paired_end = paired_end

    def output_filename(self, barcode, is_barcode=True, no_gzip=False):
        file_ext = self.file_extension
        if is_barcode:
            ret = os.path.join(self.output_dir,
                               os.extsep.join([self.protoname,
                                               'barcode_' + barcode,
                                               file_ext]))
        else:
            ret = os.path.join(self.output_dir,
                               os.extsep.join([self.protoname,
                                               barcode, file_ext]))
        if ret.endswith("gz") or no_gzip:
            return ret
        else:
            return os.extsep.join([ret, "gz"])

    def output_filename2(self, barcode, is_barcode=True, no_gzip=False):
        file_ext = self.file_extension
        if is_barcode:
            ret = os.path.join(self.output_dir,
                               os.extsep.join([self.protoname2,
                                               'barcode_' + barcode,
                                               file_ext]))
        else:
            ret = os.path.join(self.output_dir,
                               os.extsep.join([self.protoname2,
                                               barcode,
                                               file_ext]))
        if ret.endswith("gz") or no_gzip:
            return ret
        else:
            return os.extsep.join([ret, "gz"])

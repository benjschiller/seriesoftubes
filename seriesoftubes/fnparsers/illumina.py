import scripter
import os

class BarcodeFilenameParser(scripter.FilenameParser):
    def __init__(self, filename, verbose=False, *args, **kwargs):
        super(BarcodeFilenameParser, self).__init__(filename,
                                                    *args, **kwargs)
        protoname = self.protoname
        # check for old-style
        if os.path.splitext(protoname)[-3:] == 'all':
            protoname = protoname[0:-4]
        if kwargs['no_target']: self.output_dir = self.input_dir
        
        # check if this is a paired-end file
        # if so, grab its partner
        input_file = self.input_file
        illumina_name = os.path.basename(input_file)
        if illumina_name.count('_') >= 3:
            scripter.debug('NOTICE: Detected paired read file.')
            iln_parts = illumina_name.split('_')
            if iln_parts[2] == '1':
                scripter.debug('Attempting to find second file.')

                second_file = os.sep.join([self.input_dir,
                                           '_'.join(iln_parts[0:2] + ['2'] 
                                                    + iln_parts[3:])])
                self.protoname2 = os.path.splitext(
                                        os.path.basename(second_file))[0]
                try:
                    scripter.assert_path(second_file)
                    scripter.debug('Found %s', second_file)
                    self.second_file = second_file
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
        else: paired_end = False
        self.paired_end = paired_end

    def output_filename(self, barcode, is_barcode=True, no_gzip=False):
        file_ext = self.file_extension
        if no_gzip: end = '' 
        else: end = '.gz'
        if is_barcode:
            return os.path.join(self.output_dir,
                                os.extsep.join([self.protoname,
                                                'barcode_' + barcode,
                                                file_ext + end]))
        else:
            return os.path.join(self.output_dir,
                                os.extsep.join([self.protoname,
                                                barcode,
                                                file_ext + end]))
            
    def output_filename2(self, barcode, is_barcode=True, no_gzip=False):
        file_ext = self.file_extension
        if no_gzip: end = '' 
        else: end = '.gz'
        if is_barcode:
            return os.path.join(self.output_dir,
                                os.extsep.join([self.protoname2,
                                                'barcode_' + barcode,
                                                file_ext + end]))
        else:
            return os.path.join(self.output_dir,
                                os.extsep.join([self.protoname2,
                                                barcode,
                                                file_ext + end]))

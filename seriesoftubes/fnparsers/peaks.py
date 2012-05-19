from twobitreader import TwoBitFile, TwoBitFileError, twobit_reader
from twobitreader.download import save_genome
from urllib2 import URLError
from scripter import Usage, InvalidFileException, FilenameParser, \
                     debug, warning, error
from os.path import abspath, exists, extsep, join, splitext
from re import sub
from os import sep

def try_to_find_genome(genome):
    """
    takes a string that may refer to a file or genome hosted on UCSC
    and tries to open the 2bit file or download it if needed
    """
    if exists(genome):
        try:
            t = TwoBitFile(genome)
            return t
        except TwoBitFileError:
            warning('%s is not a valid .2bit file', genome)
        return None
    else:
        try:
            save_genome(genome)
            try:
                t = TwoBitFile("%s.2bit" % genome)
                return t
            except TwoBitFileError:
                error('Downloaded %s.2bit but %s.2bit is not a valid .2bit file', genome)
        except URLError:
            warning('Download of %s.2bit failed', genome)
        return None
        

class PeaksFilenameParser(FilenameParser):
    def __init__(self, filename, include_width_in_name=False,
                 target=None, motif_file='unknown_motif', genome=None,
                 *args, **kwargs):
        fext = splitext(filename)[1].lstrip(extsep)
        if fext == 'bed':
            self.is_bed = True
            self.is_xls = False
        elif fext == 'xls':
            self.is_bed = False
            self.is_xls = True
        else: raise InvalidFileException
        motif_name = sub('\W', '_', abspath(motif_file))
        target = target + sep + motif_name
        super(PeaksFilenameParser, self).__init__(filename,
                                             target = target,
                                             *args, **kwargs)
        self.fasta_file = None
        for file_extension in ['fa', 'fasta', 'FA', 'FASTA']:
            fasta_file = join(self.input_dir,
                            extsep.join([self.protoname, file_extension]))
            debug("Trying", fasta_file)
            if exists(fasta_file):
                self.fasta_file = fasta_file
                debug("Using", fasta_file)
                break
        if self.fasta_file is None:
            warning('Could not find the FASTA file for %s',
                          self.input_file)
            if genome is None:
                raise Usage("Could not find the FASTA file for ", self.input_file,
                            " and no genome was specified")
            else:
                t = try_to_find_genome(genome)
                if t is None:
                    raise Usage("Could not find the FASTA file for ", self.input_file,
                                " and failed to use %s" % genome)
                else:
                    fasta_file = join(self.input_dir, '%s.fa' % self.protoname)
                    debug('Creating FASTA file %s for %s using %s',
                                   fasta_file, self.input_file, genome)
                    input_fhd = open(self.input_file, 'rU')
                    fasta_fhd = open(fasta_file, 'w')
                    twobit_reader(t, input_stream=input_fhd,
                                  write=fasta_fhd.write)
                    fasta_fhd.close()
                    self.fasta_file = fasta_file
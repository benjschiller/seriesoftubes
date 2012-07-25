from gzip import GzipFile
from bz2 import BZ2File
try:
    from scripter import path_to_executable, Usage
    from subprocess import Popen, PIPE
    from os import pipe
    from io import TextIOBase
    try: PATH_TO_GZIP = path_to_executable('gzip')
    except Usage: pass
except ImportError: pass
from sys import stderr

# slow for reading, fast for writing
def gzip_class_factory(path_to_gzip='gzip'):
    class gzip_open_func(object):
        """gzip open func
        modes:
        (r) read using gzip.GzipFile
        (w) write using system gzip
        (P) PIPE from `gzip -d` for 
        """
        def __init__(self, filename, mode='r'):
            self._mode = mode[0]
            if mode[0] == 'w':
                args = [path_to_gzip, '-f', '-c', '-']
                stdout = open(filename, 'wb')
                self._stdout = stdout
                self.proc = Popen(args, stdin = PIPE, bufsize=0,
                                  stderr=stderr,
                                  stdout=stdout)
                self.write = self.proc.stdin.write
#                self.write = self._fakewrite
                self.close = self._close_w
                self.filename = filename
            elif mode[0] == 'r':
                self._file = GzipFile(filename, 'rb')
                self.read = self._file.read
                self.readline = self._file.readline
                self.close = self._file.close
            else:
                raise NotImplementedError
        
#        def _fakewrite(self, x):
#            print x
#            if x.count('\x1a') > 0:
#                raise ValueError(x)
#            self.proc.stdin.write(x)
        def __iter__(self):
            if not self._mode == 'r':
                raise IOError('cannot iterate over file in write mode')
            return self._file.__iter__()
        
        def _close_w(self):
            """returns returncode
            """
            self.proc.stdin.flush()
            self._stdout.flush()
#            self.proc.stdout.flush()
            self._stdout.close()
#            self.proc.communicate('\x1a')
            return
        
    return gzip_open_func

def discover_file_format(filename):
    """discover the format of a file
    returns a tuple (open_function, 'FORMAT')
    
    open_function will either be open, gzip.GzipFile, bz2.BZ2File, or None
    FORMAT can be 'BAM', 'SAM', 'FASTQ', 'FASTA', or None
    """
    f = open(filename, 'rb')
    head = f.read(3)
    f.close()
    # check magic words for compression
    if head == '\x1f\x8b\x08':
        if PATH_TO_GZIP is None:
            open_func = GzipFile
        else:
            open_func = gzip_class_factory(PATH_TO_GZIP)
    elif head=='\x42\x5a\x68':
        open_func = BZ2File
    else:
        open_func = open
    uncompressed = open_func(filename)
    head2 = uncompressed.read(4)
    head2b = uncompressed.readline()
    # check for BAM
    if head2 == 'BAM\x01': return (open_func, 'BAM')
    # check for SAM 
    if head2 == '@HD\t':
        if head2b[0:3] == 'VN:': return (open_func, 'SAM')
    # check fastq
    title = head2 + head2b
    seq = uncompressed.readline()
    title2 = uncompressed.readline()
    qual = uncompressed.readline()
    # illumina broke FASTQ convention, check for @
    if len(seq) == len(qual) and \
       (title.startswith('@') or title2.startswith('+')) and \
       (title2[1:].strip() == '' or title[1:] == title2[1:]):
        return (open_func, 'FASTQ')
    # otherwise give up
    else: return (None, None)

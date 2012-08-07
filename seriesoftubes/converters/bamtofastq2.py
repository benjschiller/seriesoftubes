"""
Another module for converting bamtofastq, defaults to writing gzfiles
Accepts only one input file

Includes class for converting paired-end sequencing files to two fastq pipes

Supports SAM, BAM
"""
from argparse import ArgumentParser
from gzip import GzipFile
from os import mkfifo, getcwd, devnull
import os
from os.path import join, exists, abspath
from pysam import Samfile
from sys import executable, stdin
from subprocess import Popen
import select
from select import PIPE_BUF
from .discover import PATH_TO_GZIP, gzip_class_factory

class PairedBAMToFastqConverter(object):
    """Works with any SAM/BAM file"""
    def __init__(self, file_, wd=None, stderr=None, logger=None):
        self.require_bam(file_)
        if wd is None: wd = getcwd()
        fifofile1 = join(wd, '._sot_fifo1')
        fifofile2 = join(wd, '._sot_fifo2')
        if exists(fifofile1):
            raise RuntimeError('%s already exists' % fifofile1)
        if exists(fifofile2):
            raise RuntimeError('%s already exists' % fifofile2)
        mkfifo(fifofile1)
        mkfifo(fifofile2)
        self.fifos = (abspath(fifofile1), abspath(fifofile2))
        self._args = [executable, '-m', __module__,
                      file_, fifofile1, fifofile2]
        self._logger = logger
        self._stderr = stderr
    
    def launch(self):
        if self._logger is not None:
            self._logger.info('Launching %s', ' '.join(self._args))
        self.subprocess = Popen(self._args, stdout=open(devnull, 'w'),
                                stderr=self._stderr, bufsize=-1)
    
    def get_fifo_readers(self):
        return (open(self.fifos[0], 'r'), open(self.fifos[1], 'r'))
    
    def require_bam(self, filename):
        with open(filename, 'rb') as f:
            head = f.read(3)
        # check magic words for compression
        if head == '\x1f\x8b\x08':
            open_func = GzipFile
        else:
            open_func = open
        uncompressed = open_func(filename)
        head2 = uncompressed.read(4)
        # check for BAM
        if head2 == 'BAM\x01': return
        # check for SAM 
        if head2 == '@HD\t': return
        else: raise ValueError('Not a SAM/BAM file') 

def pair_writer(out1, out2):
    def writer(read1, read2):
        POLLOUT = select.POLLOUT
        rec = '@%s\n%s\n+\n%s\n' % (read1.qname, read1.seq, read1.qual)
        N = len(rec) / PIPE_BUF + 1 - (len(rec) % PIPE_BUF == 0)
        i = 0
        while True:
            if i == N: break
            while True:
                x = select.select([], [out1], [], 2)
                if x[1] == [out1]: break
            y = len(rec[(i*PIPE_BUF):((i+1)*PIPE_BUF)])
            x = os.write(out1, rec[(i*PIPE_BUF):((i+1)*PIPE_BUF)])
            if not x==y: raise IOError('parital write')
            i += 1
            
        rec = '@%s\n%s\n+\n%s\n' % (read2.qname, read2.seq, read2.qual)
        N = len(rec) / PIPE_BUF + 1 - (len(rec) == PIPE_BUF)
        i = 0
        while True:
            if i == N: break
            while True:
                x = select.select([], [out2], [], 2)
                if x[1] == [out2]: break
            y = len(rec[(i*PIPE_BUF):((i+1)*PIPE_BUF)])
            x = os.write(out2, rec[(i*PIPE_BUF):((i+1)*PIPE_BUF)])
            if not x==y: raise IOError('parital write')
            i += 1
    return writer

def gzwriter(out1, out2):
    def writer(read1, read2):
        out1.write('@%s\n%s\n+\n%s\n' % (read1.qname, read1.seq, read1.qual))
        out2.write('@%s\n%s\n+\n%s\n' % (read2.qname, read2.seq, read2.qual))
    return writer

def main():
    """
    what to do if we execute the module as a script
    (not intended for user by user)
    """
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('infile', default=stdin,
                        help='BAM/SAM input file (default: stdin)')
    parser.add_argument('--gzip', action='store_true', default=False, 
                        help='Compress paired-end output files')
    parser.add_argument('outfile1', help='Output file for first mate / single reads (default: stdout)')
    parser.add_argument('outfile2', help='(required for paired-end files) output filename for second mate')
    args = parser.parse_args()
    context = vars(args)
    outfile1 = context['outfile1']
    outfile2 = context['outfile2']
    f = Samfile(context['infile'])
    incomplete_pairs = []
    if context['gzip']:
        if PATH_TO_GZIP is not None:
            open_func = gzip_class_factory(PATH_TO_GZIP)
            fh1 = open_func(outfile1, 'w')
            fh2 = open_func(outfile2, 'w')
        else:
            fh1 = GzipFile(outfile1, 'wb')
            fh2 = GzipFile(outfile2, 'wb')
        is_paired = False
        gzwrite = gzwriter(fh1, fh2)
        for aread in f:
            is_paired = False
            qname = aread.qname
            for i in xrange(len(incomplete_pairs)):
                if incomplete_pairs[i].qname == qname:
                    mate_read = incomplete_pairs.pop(i)
                    # figure out order
                    if aread.flag & 0x4 == 0x4:
                        gzwrite(aread, mate_read)
                    else:
                        gzwrite(mate_read, aread)
                    is_paired = True
                    break
            if not is_paired: incomplete_pairs.append(aread)
        unpaired = len(incomplete_pairs)
        out1.close()
        out2.close()
        f.close()
    else:
        if not exists(outfile1): os.mknod(outfile1)
        if outfile2 is not None:
            if not exists(outfile2): os.mknod(outfile2)
        out1 = os.open(outfile1, os.O_WRONLY|os.O_NONBLOCK)
        out2 = os.open(outfile2, os.O_WRONLY|os.O_NONBLOCK)
        is_paired = False
        write = pair_writer(out1, out2)
        for aread in f:
            is_paired = False
            qname = aread.qname
            for i in xrange(len(incomplete_pairs)):
                if incomplete_pairs[i].qname == qname:
                    mate_read = incomplete_pairs.pop(i)
                    # figure out order
                    if aread.flag & 0x4 == 0x4:
                        write(aread, mate_read)
                    else:
                        write(mate_read, aread)
                    is_paired = True
                    break
            if not is_paired: incomplete_pairs.append(aread)
        unpaired = len(incomplete_pairs)
        os.close(out1)
        os.close(out2)
        f.close()
    if not unpaired == 0:
        raise RuntimeError('%d unpaired reads remaining' % unpaired) 

if __name__ == '__main__': main()

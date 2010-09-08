#!/usr/bin/env python
'''
Finds and sequences errors (Ns) at each base position in FASTQ files

Default options: --target=errors.FASTQ
'''
import bioplus.positionmatrix
from bioplus.positionmatrix import positionMatrix
import Bio.SeqIO
import os
import scripter
scripter.SCRIPT_DOC = __doc__
scripter.SCRIPT_VERSION = "2.1"
scripter.SOURCE_DIR = 'sequences.FASTQ'
scripter.TARGET_DIR = 'errors.FASTQ'

def FASTQ_generator(foo):
    seqRecords = Bio.SeqIO.parse(open(foo, "rU"), "fastq-illumina")
    return (str(x.seq) for x in seqRecords)

def detect_errors(parsed_filename, verbose=False, debug=False, **kwargs):
    if verbose:
        stdout_buffer = ''.join(['Determining basewise error rates in',
                                 parsed_filename.input_file, os.linesep])
    else:
        stdout_buffer = ''
    seqs = FASTQ_generator(parsed_filename.input_file)
    pms = positionMatrix()
    def row(*args):	return "".join(["\t".join([str(x) for x in args]),"\n"])
    firstRow = row("cycle", "pyindex", "A", "T", "G", "C", "N")
    pms.countSeqs( seqs )
    output_filename = os.path.join(parsed_filename.output_dir,
                                   parsed_filename.with_extension(
                                        os.extsep.join(['errors','txt'])))
    f = open(output_filename,'w')
    f.write(firstRow)
    for i in range(len(pms)):
        p = pms[i]
        f.write(row(i+2, i, p['A'], p['T'], p['G'] ,p['C'], p['N']))
    f.close()
    return stdout_buffer

if __name__=="__main__":
    scripter.perform(detect_errors)

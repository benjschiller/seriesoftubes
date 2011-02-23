#!/usr/bin/env python
'''
Finds and sequences errors (Ns) at each base position in FASTQ files

Default options: --target=errors.FASTQ
'''
#import bioplus.positionmatrix
#from bioplus.positionmatrix import positionMatrix
import Bio.SeqIO.QualityIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import os
import operator
import scripter

VERSION = "2.4"

def main():
    e = scripter.Environment(version=VERSION, doc=__doc__)
    e.source_dir = 'sequences.FASTQ'
    e.target_dir = 'errors.FASTQ'
    e.do_action(detect_errors)

def detect_errors(parsed_filename, verbose=False, debug=False, **kwargs):
    if debug:
        stdout_buffer = 'Determining basewise error rates in {!s}\n'.format(
                                                    parsed_filename.input_file)
    else: stdout_buffer = ''
    
    handle = open(parsed_filename.input_file, 'rU')
    record_generator = FastqGeneralIterator(handle)
    
    title, seq, qual = record_generator.next()
    counter = lambda x: [x == 'A', x =='T', x == 'G', x == 'C', x == 'N']
    counts = map(counter, seq.upper())
    count_len = len(counts)
    j=0
    for title, seq, qual in record_generator:
        j+=1
        if debug and j%(5*10**5)==0: print j
        more_counts = map(counter, seq.upper())
        counts = (map(operator.add, counts[i], more_counts[i])
                  for i in range(count_len))
            
    output_filename = os.path.join(parsed_filename.output_dir,
                                   parsed_filename.with_extension(
                                        os.extsep.join(['errors','txt'])))
    f = open(output_filename, 'w')
    first_row = ("cycle", "pyindex", "A", "T", "G", "C", "N")
    f.write('\t'.join(["{!s}"]*6).format(*first_row))
    counts = list(counts)
    for i in range(count_len):
        row = (i+2, i, counts[i][0], counts[i][1], counts[i][2], counts[i][3],
               counts[i][4])
        f.write('\t'.join(["{!s}"]*6).format(*row))
    f.close()
    return stdout_buffer

if __name__=="__main__": main()
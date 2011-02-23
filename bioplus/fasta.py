# this module contains meta-tools for large FASTA files
# and requires python
from Bio.SeqIO import parse, write
import random
import itertools
import os.path

def count(foo):
    '''takes a file named foo returns the number lines'''
    f = parse(open(foo,'rU'),'fasta')
    n = 0
    for dummyX in f:
        n += 1
    return n
    
def random(foo,n):
    '''takes a file foo and returns n random sequences from it'''
    max_n = count(foo)
    record_numbers = itertools.repeat(random.randint(1,max_n),times=n)
    seq_recs = parse(open(foo,'rU'),'fasta')
    i = 0
    seqs = []
    for rec in seq_recs:
        i += 1
        if i in record_numbers: seqs.append(rec)
    return seqs

def reader(foo):
    '''
    '''
    list_of_seqs = [record.seq for record in parse(open(foo),'fasta')]
    return list_of_seqs

def writer(foo, L):
    '''
    write to FASTA file foo the list of SeqRecord objects L.
    Warning: overwrites foo, does not append
    '''
    write(L, open(foo, 'w'), 'fasta')

def random_files(foo,n,R):
    '''
    takes a FASTA file foo and creates (in the current directory) R random
    files each containing n random sequences from foo, named foo_random[0-R].fa
    '''
    for filenumber in R:
        prefix = foo.split(os.path.sep)[-1]
        boo = prefix + '_random' + pad(filenumber,R) + '.fa'
        writer(boo,random(foo,n))

def pad(x,y):
    '''takes two integers x and y, and returns str(x) with enough 0s to match the length of str(y)'''
    strX = str(x)
    strY = str(y)
    return strX.zfill(len(strY))

def truncate_lines(f,n):
    '''
    truncate_lines(f,n) truncates lines in a file to at most n letters.
    See also truncate_seqs
    '''
    seqs = open(f,'r')
    # USE FILENAME CORRECTION SCHEME
    tseqs = open(f + '.' + str(n),'w')
    for seq in seqs: tseqs.writeline(seq[0:n]+os.linesep)
    seqs.close()
    tseqs.close()

def truncate_seqs(f, n):
    '''
    truncate FASTA seqs truncates sequences in a FASTA file named f to at most
    n bases, writing a new FASTA file named f.n.fa
    '''
    seq_recs = parse(open(f,'rU'),'fasta')
    # USE FILENAME CORRECTION SCHEME
    foo = f+str(n)+'.fa'
    writer(foo, [rec[0:n] for rec in seq_recs])

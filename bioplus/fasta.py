'''meta-tools for large FASTA files'''
import random
import itertools
import os.path
from Bio.SeqIO import parse, write
import Bio.Seq

def count(foo):
    '''takes a file named foo returns the number lines'''
    f = parse(open(foo,'rU'),'fasta')
    n = 0
    for dummyX in f:
        n += 1
    return n
    
def random_seq(foo, n):
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
    generator yielding Bio.Seq.Seq objects from a FASTA file
    '''
    for record in parse(open(foo, 'rU'), 'fasta'):
        yield record.seq

def writer(foo, iterable):
    '''
    writes SeqRecord objects from iterable to FASTA file foo.
    Warning: overwrites foo, does not append
    '''
    write(iterable, open(foo, 'w'), 'fasta')

def random_files(foo, n, R):
    '''
    takes a FASTA file foo and creates (in the current directory) R random
    files each containing n random sequences from foo, named foo_random[0-R].fa
    '''
    for filenumber in R:
        prefix = foo.split(os.path.sep)[-1]
        boo = prefix + '_random' + pad(filenumber, R) + '.fa'
        writer(boo,random(foo,n))

def pad(x, y):
    '''
    takes two integers x and y, and returns str(x) with enough 0s to match
    the length of str(y)
    '''
    return str(x).zfill(len(str(y)))

def truncate_lines(f, n):
    '''
    truncate_lines(f,n) truncates lines in a file to at most n characters
    See truncate_seqs to truncate sequences instead of lines
    '''
    with open(f,'r') as seqs, open('{!s}.{!s}'.format(f, n), 'w') as tseqs:
        tseqs.writelines((seq[0:n]+'\n' for seq in seqs))

def truncate_seqs(f, n):
    '''
    truncate FASTA seqs truncates sequences in a FASTA file named f to at most
    n bases, writing a new FASTA file named f.n.fa
    '''
    seq_recs = parse(open(f,'rU'),'fasta')
    # USE FILENAME CORRECTION SCHEME
    foo = f+str(n)+'.fa'
    writer(foo, (rec[0:n] for rec in seq_recs))

def permute_fasta(f):
    '''
    takes a FASTA file and returns a new FASTA file with each sequence randomly
    permuted (separately, such that its % A,T,G,C doesn't change)
    '''
    mute = Bio.Seq.Seq.tomutable
    shuffle = random.shuffle
    with open(f + '_permuted.fa', 'w') as output:
        with open(f, 'rU') as fobj:
            for seq_rec in parse(fobj, 'fasta'):
                seq_rec.seq = mute(seq_rec.seq)
                shuffle(seq_rec.seq)
                write(seq_rec, output, 'fasta')
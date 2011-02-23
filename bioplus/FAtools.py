# this module contains meta-tools for large FASTA files
# and requires python
import Bio.SeqIO
import random
import itertools
import os.path

def FAcount(foo):
    '''takes a file named foo returns the number lines'''
    f = Bio.SeqIO.parse(open(foo,'rU'),'fasta')
    n = 0
    for dummyX in f:
        n += 1
    return n
    
def FArandom(foo,n):
    '''takes a file foo and returns n random sequences from it'''
    maxN = FAcount(foo)
    recordNumbers = itertools.repeat(random.randint(1,maxN),times=n)
    seqRecs = Bio.SeqIO.parse(open(foo,'rU'),'fasta')
    i = 0
    seqs = []
    for rec in seqRecs:
        i += 1
        if i in recordNumbers: seqs.append(rec)
    return seqs

def FAreader(foo):
    '''
    '''
    list_of_seqs = [record.seq for record in Bio.SeqIO.parse(open(foo),'fasta')]
    return list_of_seqs

def FAwriter(foo,L):
    '''write to FASTA file foo the list of SeqRecord objects L. Warning: overwrites foo, does not append'''
    Bio.SeqIO.write(L,open(foo,'w'),'fasta')

def FArandomfiles(foo,n,R):
    '''takes a FASTA file foo and creates (in the current directory) R random files each containing n random sequences from foo, named foo_random[0-R].fa'''
    for filenumber in R:
        prefix = foo.split(os.path.sep)[-1]
        boo = prefix + '_random' + pad(filenumber,R) + '.fa'
        FAwriter(boo,FArandom(foo,n))

def pad(x,y):
    '''takes two integers x and y, and returns str(x) with enough 0s to match the length of str(y)'''
    strX = str(x)
    strY = str(y)
    return strX.zfill(len(strY))

def truncateLines(f,n):
    '''truncateLines(f,n) truncates lines in a file to at most n letters. See also truncateFASTA'''
    seqs = open(f,'r')
    # USE FILENAME CORRECTION SCHEME
    tseqs = open(f + '.' + str(n),'w')
    for seq in seqs: tseqs.writeline(seq[0:n]+os.linesep)
    seqs.close()
    tseqs.close()

def truncateSeqs(f,n):
    '''truncate FASTA seqs truncates sequences in a FASTA file named f to at most n bases, writing a new FASTA file named f.n.fa'''
    '''takes a file foo and returns n random sequences from it'''
    seqRecs = Bio.SeqIO.parse(open(f,'rU'),'fasta')
    # USE FILENAME CORRECTION SCHEME
    foo = f+str(n)+'.fa'
    Bio.SeqIO.write([ rec[0:n] for rec in seqRecs  ],open(foo,'w'),'fasta')

# this module wraps the built-in module random
# and provides some useful random sequence generators, etc.

import random
import itertools

def random_seq(n=1,GC=0.5):
	'''random_seq provides a random nucleotide (A, T, G, or C). You may optionally provide n, a positive integer, which will cause random_seq to return a string of n nucleotides. You may also optionally provide GC, the probability of encountering a G or C, which must be on the closed interval [0,1]. The probability of encountering an A or T is calculated as 1 - GC.
	'''
	return ''.join( list(random_seq_generator(n,GC)) )

def random_seq_generator(n=1,GC=0.5):
	'''
	random_seq_generator acts like random_seq, but returns a generator that
	returns the nucleotides one by one
	'''
	myError = ValueError('randomN requires a positive integer n (default = 1) \
and a probability GC (float 0.0 to 1.0)')
	#AT = 1 - GC
	if not type(GC)==float or GC == 0 or GC == 1: raise myError 
	elif not GC >= 0 and GC <= 1: raise myError
	elif not type(n)==int:	raise myError
	elif n < 1: raise myError
	else:
		randomGenerator = itertools.repeat(random.random())
		gc_or_at = lambda x: random.choice(['G','C']) \
							 if x < GC else random.choice(['A','T'])
		return itertools.imap(gc_or_at, randomGenerator)

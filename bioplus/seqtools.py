import collections
import operator
from itertools import izip, imap
from numpy.core.fromnumeric import std 

def get_first(seq1, seq2, ignore_case=True):
	"""
	returns the alphabetically prior of seq1 and seq2
	if they are the same or we fail to order, return seq1
	
	if ignore_case=True, we apply str.upper to each letter
	"""
	if ignore_case:
		upper = str.upper
		for a, b in izip(seq1, seq2):
			A = upper(a)
			B = upper(b)
			if A == B: continue
			elif A > B: return seq1
			elif B > A: return seq2
	else:
		for a, b in izip(seq1, seq2):
			if a == b: continue
			elif a > b: return seq1
			elif b > a: return seq2
	return seq1
	
def count_seqs(iterable, reverse_complement=True, ignore_case=True):
	'''
count_seqs takes as input an iterable that yields sequences, and then
returns the counts of each sequence as a defaultdict (similar to dict, can
be recast as dict type).

Sorting: Use dictSort to sort, if needed.

Reverse complements: if rc=True, count_seqs will count reverse complements
as the same sequence and use the alphabetically prior sequence as the key.
if rc=False, reverse complements will be treated as different sequences.

	'''
	seq_counts = collections.defaultdict(lambda: 0)
	if reverse_complement:
		orient_seq = lambda s: get_first(s, rc(s), ignore_case=ignore_case)
		oriented_seqs = imap(orient_seq, iterable)
		for seq in oriented_seqs: seq_counts[seq]+=1
	else: # if not rc
		for seq in iterable: seq_counts[seq]+=1
	return seq_counts

def count_compare(a, b, rc=False):
	"""
count_compare takes two default dictionaries (defaultdict objects),
a and b, each with default factory "lambda: 0" and key-value pairs that specify
the counts (value) of each sequence (key), and returns a new dictionary whose
values are tuples (value_a,value_b), the values from a and b.

Reverse complements: if rc=True, count_seqs will count reverse complements as
the same sequence and use the alphabetically prior sequence as the key.
if rc=False, reverse complements will be treated as different sequences.
If you have already ensured that the keys meet this condition, you should
use rc=False, but rc=True is also safe.

note: defaultdict is in the collections module of the standard library
	"""
	if not (type(a)==collections.defaultdict and type(b)==collections.defaultdict):
		raise ValueError('a and b must be default dictionaries (collections.defaultdict)')
	# proceed
	comparedCounts = {}
	for k in set( a.keys(), b.keys() ): comparedCounts[k] = (a[k],b[k])
	return comparedCounts

def analyze_sites(some_list):
	'''
	requires the following input format:
	seq score cdist cons
	cdist = distance to center, cons = conservation value
	'''
	# column definitions for reference / calling
	x0 = 0 # x0 is sequence
	x1 = 1 # x1 is score (Ri)
	x2 = 2 # x2 is peaklength
	x3 = 3 # x3 is cdist
	x4 = 4 # x4 is conservation score

	# columns are seq, score, peaklength, cdist, cons
	# build an array scalar with the correct column types
	for x in some_list:
		# forces sequences to uppercase
		x[x0] = x[x0].upper()
		# forces length to int
		x[x2] = int(x[x2])
		# forces cdist to int
		x[x3] = int(x[x3])
		# forces cons score to floating point
		# and correct nan
		if x[x4]=='nan': x[x4] = 0.0
		else: x[x4] = float(x[x4])
	# sort the list by sequence (arbitrarily)
	sorted_list = sorted(some_list, key=operator.itemgetter(x0))
	# Part 3: continue with analysis
	analyzedSeqs = []
	counter = len(sorted_list)
	while counter > 0:
		site = sorted_list.pop()
		counter -= 1
		seq = site[x0]
		score = site[x1]
		# start a list of dists, absdists, cscores, we'll sum later
		i = 1
		plengths = [ site[x2] ]
		dists = [ site[x3] ]
		absdists = [ abs(site[x3]) ]
		cscores = [ site[x4] ]
		while True: # grab sequences that match
			try: y = sorted_list.pop()
			except IndexError: break
			counter -= 1
			# look if we match the last entry
			if y[x0] == seq:
				i += 1
				plengths.append(y[x2])
				dists.append(y[x3])
				absdists.append( abs( y[x3] ) )
				cscores.append(y[x4])
			# if we don't, escape to previous loop
			else:
				# append back
				sorted_list.append(y)
				break
		# append seq, sum of scores, # of instances
		avgPlength = float(sum(plengths))/i
		avgDist = float(sum(dists))/i
		avgAbsdist = float(sum(absdists))/i
		distErr = std(absdists)
		analyzedSeqs.append([seq, score, i, sum(cscores), avgPlength, avgDist,
							avgAbsdist,distErr])
	analyzedSeqs.reverse()
	return sorted(analyzedSeqs, key=operator.itemgetter(x2))

def complement(nucleotide):
	'''returns the complement of a single nucleotide'''
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
	return complement.setdefault(nucleotide,'N')

def rc(seq):
	'''returns the reverse complement of a sequence'''
	L = [complement(s) for s in seq]
	L.reverse()
	return ''.join(L)

reversecomplement = rc
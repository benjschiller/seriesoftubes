'''tools for dealing with motifs that have spacers, i.e. stretches of unbound DNA within the motif'''

from tabfile import TabFile
from positionmatrix import genDict
import collections
import operator

def spacerGC(L, spacerOffset=6, spacerLength=3):
	'''spacerGC takes as input a list of [15 bp GBSs (strings)] and returns the number of sequences that have 0,1,2,3 G/Cs in the 3 bp spacer as an array in that order'''
	# gc counts of 0, 1, 2, 3
	GCcounts = [0, 0, 0, 0]
	for s in L:
		spacer = s[0][spacerOffset:spacerOffset+spacerLength].upper()
		GCs = spacer.count('C') + spacer.count('G')
		GCcounts[GCs] += 1
	return GCcounts

def center_region(f, max_dist=75, motif_length=17):
	'''returns a function that specifies whether a given motif is in +/-
	x bp from the peak_center

	requires the tabFile object f to determine the indices properly
	'''
	column_dict = f.column_dict()
	peak_summit = column_dict['peak_summit']
	for offset_name in ('motif_offset', 'offset', 'max_offset'):
		if column_dict.has_key(offset_name):
			site_offset = column_dict[offset_name]
			break
	return lambda x: int(x[peak_summit]) > (int(x[site_offset]) -
						max_dist) and \
			int(x[peak_summit]) < (int(x[site_offset]) +
						max_dist - motif_length)

def count_spacers_from_info(foo, cutoff=None, region_rule=None,
	region_width=None, spacer_offset=8, spacer_length=3, output_file=None):
	'''count spacers from a .sites.info or .peaks.info file
	
	optionally you may supply
	cutoff, a minimum cutoff (float or int)
	region_rule, a function that selects the column 
	'''
	input_file = TabFile(foo, colNames=True)
	rows = (x for x in input_file)
	conditions = [lambda x: x[7].isalpha, # col 7 is a sequence
				  lambda x: x[7] is not '-', # col 7 is not -
				  lambda x: x[7] is not 'NA', # col 7 is not NA
				  lambda x: x[7].strip() is not '-'] # col 7 is not missing

	if cutoff is not None:
		conditions.append(lambda x: float(x[4])>cutoff)

	if region_rule is 'center_region':
		if region_width is not None:
			conditions.append(center_region(input_file,
						max_dist=region_width/2))
		else:
			conditions.append(center_region(input_file,
							max_dist=75))
	elif region_rule is not None:
		conditions.append(lambda x: region_rule(x))

	selected_rows = (x[7].upper() for x in rows if
				all([f(x) for f in conditions]))

	returnZero = lambda: 0
	spacers = collections.defaultdict(returnZero)
	for s in selected_rows:
		if not s== '-' and not s == 'NA':
			spacer = s[spacer_offset:spacer_offset + spacer_length].upper()
			spacers[spacer] += 1
	if output_file is None:
		output_file = raw_input('Output file name: ')
	with TabFile(output_file, 'w') as f:
		f.write_table(sorted(spacers.iteritems(),
						     key=operator.itemgetter(1), reverse=True))
	return spacers

def count_letters(K):
	n=range(len(K[0]))
	counts = []
	for j in n:
		counts.append(genDict({}))
	for x in K:
		for j in n:
			counts[j][x[j]]+=1
	return counts

def countLetters(d):
	counts = []
	n=range(len(d.keys()[0]))
	for j in n:
		counts.append(collections.defaultdict(lambda: 0))
	for x in d.keys():
		for j in n:
			counts[j][x[j]]+=d[x]
	return counts

#def runSpacer():
#	'''runs the full analysis on all cell types'''
#	# home
#	home = '/work/solexa/'
#	peakdir = home + 'peaks/'
#	cells = ['A549','Nalm6', 'U2OS']
#	shifts = ['','_shifted']
#	
#	for cellType in cells:
#		for s in shifts:
#			f = tabFile(peakdir+cellType+s+'.peaks.bed')
#			seqs = f.readCol(-2)
#			spacerCounts(seqs, spacerOffset = 7, fname = cellType+s+'.spacers')


import math
import operator
from collections import defaultdict
from numpy.core.multiarray import zeros
from tabfile import TabFile

class CharDict(defaultdict):
    """
    CharDict is a generalized dictionary that inherits defaultdict.
    default values are 0. Accepts any keys.
    """

    def __init__(self, dict={}):
        def returnZero():
            return 0
        super(CharDict,self).__init__(returnZero, dict)

    def __add__(self, x, y):
        z = CharDict()
        for k in set(x.keys() + y.keys()):
            z[k] = x[k] + y[k]
        return z

    def uncertainty(self):
        """
        calculates the uncertainty H in this position
        (specified by CharDict). 
        reats 0*log(0) as 0
        """
        H = 0
        for pN in self.itervalues():
            if pN==0: H = 0
            else: H = -pN * math.log(pN, 2)
        return H

class PositionWeightMatrix(object):
    """
    Stores counts of nucleotide bases at each position. objects are immutable.
    sequences may be added to the counts, but the object may not be modified
    in situ
    """

    def __init__(self,n=None):
        raise DeprecationWarning("Use Bio.Motif. PositionWeightMatrix will be \
removed soon")
        self._errN = 'n must be a positive integer'
        if not n==None: self._initialize_positions(n)
        self._is_probs=False
        self._is_Ri=False

    def _initialize_positions(self,n):
        self._L = []
        if type(n)==int:
            if n > 0:
                self._n = n
                for dummyVariable in range(self._n):
                    self._L.append( CharDict() )
            else: raise ValueError(self._errN)
        else: raise TypeError(self._errN)

    def __add__(self, x,y):
        n = x._n
        if n == y._n:
            z = PositionWeightMatrix(n)
            for i in range(n):
                z[i] = x[i] + y[i]
        else: raise ValueError('PositionWeightMatrix objects are not the same \
length (number of positions)')
        return z

    def __getitem__(self, y):
        return self._L[y]
    
    def __len__(self):
        return len(self._L)

    def count_file(self, seqsFile, n=0):
        """uses a tabFile with a list of sequences, in column n (by default n=0, the first column) and extracts counts"""
        if self._is_probs: raise UserWarning('Already converted to probabilities')
        # open file
        seqsFile.open()
        # read first sequence and set siteLength
        rows = (row for row in seqsFile)
        row = rows.next()
        site = row[n]
        siteLength = len(site)
        self._initialize_positions(siteLength)
        # initialize the object
        for i in range(self._n): self._L[i][ site[i].upper() ] += 1
        # read remaining sequences
        while True:
            try:
                row = rows.next()
                site = row[n]
            except StopIteration: break
            if len(site)==siteLength:
                for i in range(self._n): self._L[i][ site[i].upper() ] += 1
            else:
                # clean up
                del self._L
                del self._n
                seqsFile.close()
                raise ValueError('One of the sequences you are trying to add is not the correct length ('+str(self._n)+'): '+site)
        self._n = siteLength

    def count_seqs(self, L, debug=False):
        """adds a list of sequences to the counts"""
        if self._is_probs: raise UserWarning('Already converted to probabilities')
        firstSite = True
        n = 0
        m = 0
        for site in L:
            n += 1
            if n%(10**6)==0:
                m += 1
                if debug: print str(n)
            if firstSite:
                siteLength=len(site)
                self._initialize_positions(siteLength)
                firstSite = False
            for i in range(self._n):
                if len(site)==siteLength: self._L[i][ site[i].upper() ] += 1
                else:
                    # clean up
                    del self._L
                    del self._n
                    raise ValueError('One of the sequences you are trying to add is not the correct length ('+str(self._n)+'): '+site)

    def import_from_MEME(self,filename,n=1,mode='biotools'):
        """imports a motif from the output of MEME (meme.txt)

        if there are multiple motifs in the output, we will use motif n (the first is n=1, which is also the default)
        """
        import Bio.Motif.Parsers.MEME as MEME
        f = open(filename)
        MEME_object = MEME.read(f)
        motif_name = 'Motif ' + str(n)
        biopython_motif = MEME_object.get_motif_by_name(motif_name)
        if mode=='biopython': return biopython_motif
        if mode=='biotools':
            internal_n = len(biopython_motif)
            # this next line is instead of initializePositions
            biotools_motif = [CharDict(biopython_motif[i]) for i in range(internal_n)]
            self._L = biotools_motif
            self._n = internal_n    
        else: raise UserWarning('Not a valid mode.')

    def rc(self):
        """returns the reverse complement of this object"""
        new = PositionWeightMatrix(self._n)
        # complement the object
        for i in range(self._n):
            for base in self._L[i].keys():
                complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
                if complement.has_key(base): new[i][complement[base]] = self._L[i][base]
                else: new[i][base] = self._L[i][base]
        new._L.reverse() # reverse the object
        return new

    def __repr__(self):
        return self._L.__repr__()

    def make_probs(self,trueProbs=False):
        """normalizes everything to 1"""
        if self._is_Ri:
            for x in self._L:
                for k in x.keys():
                    x[k] = 2**(x[k]-2)
            self._is_Ri = False
        else:
            for x in self._L:
                total = sum( x.values() )
                zeros = x.values().count(0)
                if trueProbs or zeros==0:
                    for k in x.keys():
                        x[k] = float(x[k]) / total
                else:
                    # fake one occurrence
                    total += 1
                    for k in x.keys():
                        if x[k]==0: x[k]=1
                        x[k] = float(x[k]) / total
        self._is_probs = True

    def make_Ri(self):
        """changes from counts or probabilities to Ri, information content"""
        if not self._is_Ri:
            if not self._is_probs: self.make_probs()
            for p in self._L:
                for k in p.keys():
                    p[k] = 2+math.log(p[k],2)
            self._is_probs = False
            self._is_Ri = True
        else:
            print 'Already Ri'

    def seq_Ri(self,s):
        """seqRi returns the information content Ri in bits of a sequences, as measured with the given positionmatrix"""
        if not self._is_Ri: self.makeRi()
        Ri = 0
        if len(s) != self._n:
            raise UserWarning('Cannot evaluate a sequence which is not the exact length of the position matrix')
        for x in range(self._n): Ri += self[x][s[x].upper()]
        return Ri

    def uncertainty(self):
        """returns the uncertainty H(l) of the matrix as a list. Use sum() for the total uncertainty.

Note: this function calls uncertainty() from the baseDict instance, and as such it can be overwritten implicitly. baseDict.uncertainty() treats 0*log(0) as 0"""
        if not self._is_probs: self.make_probs()
        return [position.uncertainty() for position in self]

    def Rs(self):
        """returns the Schneider Rs value, which is the expectation of Ri over all possible sequences, calculated as the sum of 2-uncertainty."""    
        if not self._is_probs: self.make_probs()
        return sum([2 - position.uncertainty() for position in self])

def KL(p,q):
    """returns a list of the KL divergence (relative entropy) at each position from positionmatrix p to positionmatrix q. use sum() for the sum"""
    if not len(p)==len(q): raise SyntaxError('Length of p and q must be the same, instead length of p is ' + len(p) + ' and length of q is ' + len(q))
    else: n = len(p)
    for i in xrange(n):
        KLi = 0
        for j in ['A','G','T','C']:
            KLi += p[i][j] * math.log(p[i][j] / q[i][j], 2)
        KL.append(KLi)
    return KL

# requires numpy, maybe relax requierment?
# needs to return an mmMatrix object
# needs to be able to save to file
# needs to be able to make and save image
def joint_matrix(sites):
    """takes as input a filename and returns the joint Rate matrix
    for the list of sequences contained in that file

    Joint rates R(X;Y_ are defined as 
    R(X;Y) = - sum over X,Y p(x,y) * I(X;Y)
    I(X;y) = - sum over X,Y p(x,y) * log2[p(x,y)/(p(x)p(y))]
    """
    bases = ['A','C','G','T']
    indexDictionary = {} # the index dictionary
    for i in range(4):
        for j in range(4):
            ssPair = bases[i] + bases[j]
            indexDictionary[ssPair]=i,j

    site_length = len(sites[0])
# initialize the matrix
    di_counts = zeros([site_length,site_length],dtype='(4,4)int')

    def add_seq(m,s,n,b):
        """adds the dinucleotide counts from one sequence to the mm_matrix (an array, passed by refence). requires the length n"""
        for i in range(n):
            for j in range(n):
                m[i,j][ b[s[i]+s[j]] ] += 1

# count pairs over every sequence
    for site in sites:
        add_seq(di_counts,site.upper(),site_length,indexDictionary)

# convert to probabilities
    di_probs = zeros([site_length,site_length],dtype='(4,4)float')
    total_seqs = di_counts[0,0].sum()
    for i in range(site_length):
        for j in range(site_length):
                for ii in range(4):
                    for jj in range(4):
                        di_probs[i,j][ii,jj] = di_counts[i,j][ii,jj] / float(total_seqs)

    mm_matrix = zeros([site_length,site_length],dtype='float')
    for i in range(site_length):
        for j in range(site_length):
            # sum over all dinucleotide combinations
            pM = di_probs[i,j]

            # Determine Iij
            Iij = 0.0
            for x in range(4):
                for y in range(4):            
                    px = pM[x,:].sum()
                    py = pM[:,y].sum()
                    pxy = pM[x,y]
                    if any([pxy==0, py==0, px==0]): continue
                    Iij += pxy * math.log(pxy/px/py, 2)

            # Determine Rij
            Rij = 0.0
            for x in range(4):
                for y in range(4):            
                    pxy = pM[x, y]
                    Rij -= pxy * Iij

            mm_matrix[i][j] = Rij
    return (di_counts, di_probs, mm_matrix)

def spacerGC(L, spacerOffset=6, spacerLength=3):
    """
    spacerGC takes as input a list of [15 bp GBSs (strings)] and
    returns the number of sequences that have 0,1,2,3 G/Cs in the 3 bp spacer 
    as an array in that order
    """
    # gc counts of 0, 1, 2, 3
    GCcounts = [0, 0, 0, 0]
    for s in L:
        spacer = s[0][spacerOffset:spacerOffset+spacerLength].upper()
        GCs = spacer.count('C') + spacer.count('G')
        GCcounts[GCs] += 1
    return GCcounts

def center_region(f, max_dist=75, motif_length=17):
    """returns a function that specifies whether a given motif is in +/-
    x bp from the peak_center

    requires the tabFile object f to determine the indices properly
    """
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
    """
    count spacers from a .sites.info or .peaks.info file
    
    optionally you may supply
    cutoff, a minimum cutoff (float or int)
    region_rule, a function that selects the column 
    """
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

    spacers = CharDict(dict)
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

def count_letters(L):
    n=xrange(len(L[0]))
    counts = []
    for j in n:
        counts.append(CharDict())
    for x in L:
        for j in n:
            counts[j][x[j]]+=1
    return counts
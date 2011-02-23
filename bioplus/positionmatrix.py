import math
from numpy.core.multiarray import zeros
from collections import defaultdict

class genDict(defaultdict):
    '''genDict is a generalized dictionary that inerherits defaultdict. default values are 0. Accepts any keys.
    '''

    def __init__(self,dict={}):
        def returnZero():
            return 0
        super(genDict,self).__init__(returnZero,dict)

    def __add__(self, x, y):
        z = genDict()
        for k in set(x.keys()+y.keys()):
            z[k] = x[k] + y[k]
        return z
    
    def uncertainty(self):
        '''calculates the uncertainty H in this position (specified by genDict). Treats 0*log(0) as 0'''
        H = 0
        for pN in self.values():
            if pN==0: H = 0
            else: H = -pN * math.log(pN, 2)
        return H

class baseDict(dict):
    '''baseDict is deprecated. Use genDict.

    baseDict is a dictionary that always and only contains the keys 'A','T','G','C' and corresponding values. default values are 0
    '''
    def __init__(self):
        self['A'] = 0
        self['T'] = 0
        self['G'] = 0
        self['C'] = 0
    
    def __add__(self, x, y):
        z = baseDict()
        z['A'] = x['A'] + y['A']
        z['T'] = x['T'] + y['T']
        z['G'] = x['G'] + y['G']
        z['C'] = x['C'] + y['C']
        return z
    
    def uncertainty(self):
        '''calculates the uncertainty H in this position. Treats 0*log(0) as 0'''
        H = 0
        for pN in self.values():
            if pN==0: H = 0
            else: H = -pN * math.log(pN, 2)
        return H

class positionMatrix:
    '''Stores counts of nucleotide bases at each position. objects are immutable. sequences may be added to the counts, but the object may not be modified in situ'''

    def __init__(self,n=None):
        self._errN = 'n must be a positive integer'
        if not n==None: self._initializePositions(n)
        self._isProbs=False
        self._isRi=False

    def _initializePositions(self,n):
        self._L = []
        if type(n)==int:
            if n > 0:
                self._n = n
                for dummyVariable in range(self._n):
                    self._L.append( genDict() )
            else: raise SyntaxError(self._errN)
        else: raise TypeError(self._errN)

    def __add__(self, x,y):
        n = x._n
        if n == y._n:
            z = positionMatrix(n)
            for i in range(n):
                z[i] = x[i] + y[i]
        else: raise ValueError('baseCounts objects are not the same length (number of positions)')
        return z

    def __getitem__(self, y):
        return self._L[y]
    
    def __len__(self):
        return len(self._L)

    def countFile(self, seqsFile, n=0):
        '''uses a tabFile with a list of sequences, in column n (by default n=0, the first column) and extracts counts'''
        if self._isProbs: raise UserWarning('Already converted to probabilities')
        # open file
        seqsFile.open()
        # read first sequence and set siteLength
        rows = ( row for row in seqsFile )
        row = rows.next()
        site = row[n]
        siteLength = len(site)
        self._initializePositions(siteLength)
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

    def countSeqs(self, L, debug=False):
        '''adds a list of sequences to the counts'''
        if self._isProbs: raise UserWarning('Already converted to probabilities')
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
                self._initializePositions(siteLength)
                firstSite = False
            for i in range(self._n):
                if len(site)==siteLength: self._L[i][ site[i].upper() ] += 1
                else:
                    # clean up
                    del self._L
                    del self._n
                    raise ValueError('One of the sequences you are trying to add is not the correct length ('+str(self._n)+'): '+site)

    def import_from_MEME(self,filename,n=1,mode='biotools'):
        '''imports a motif from the output of MEME (meme.txt)

        if there are multiple motifs in the output, we will use motif n (the first is n=1, which is also the default)
        '''
        import Bio.Motif.Parsers.MEME as MEME
        f = open(filename)
        MEME_object = MEME.read(f)
        motif_name = 'Motif ' + str(n)
        biopython_motif = MEME_object.get_motif_by_name(motif_name)
        if mode=='biopython': return biopython_motif
        if mode=='biotools':
            internal_n = len(biopython_motif)
            # this next line is instead of initializePositions
            biotools_motif = [genDict(biopython_motif[i]) for i in range(internal_n)]
            self._L = biotools_motif
            self._n = internal_n    
        else: raise UserWarning('Not a valid mode.')

    def rc(self):
        '''returns the reverse complement of this object'''
        new = positionMatrix(self._n)
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

    def makeProbs(self,trueProbs=False):
        '''normalizes everything to 1'''
        if self._isRi:
            for x in self._L:
                for k in x.keys():
                    x[k] = 2**(x[k]-2)
            self._isRi = False
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
        self._isProbs = True

    def makeRi(self):
        '''changes from counts or probabilities to Ri, information content'''
        if not self._isRi:
            if not self._isProbs: self.makeProbs()
            for p in self._L:
                for k in p.keys():
                    p[k] = 2+math.log(p[k],2)
            self._isProbs = False
            self._isRi = True
        else:
            print 'Already Ri'

    def seqRi(self,s):
        '''seqRi returns the information content Ri in bits of a sequences, as measured with the given positionmatrix'''
        if not self._isRi: self.makeRi()
        Ri = 0
        if len(s) != self._n:
            raise UserWarning('Cannot evaluate a sequence which is not the exact length of the position matrix')
        for x in range(self._n): Ri += self[x][s[x].upper()]
        return Ri

    def uncertainty(self):
        '''returns the uncertainty H(l) of the matrix as a list. Use sum() for the total uncertainty.

Note: this function calls uncertainty() from the baseDict instance, and as such it can be overwritten implicitly. baseDict.uncertainty() treats 0*log(0) as 0'''
        if not self._isProbs: self.makeProbs()
        return [position.uncertainty() for position in self]

    def Rs(self):
        '''returns the Schneider Rs value, which is the expectation of Ri over all possible sequences, calculated as the sum of 2-uncertainty.'''    
        if not self._isProbs: self.makeProbs()
        return sum([2 - position.uncertainty() for position in self])

def KL(p,q):
    '''returns a list of the KL divergence (relative entropy) at each position from positionmatrix p to positionmatrix q. use sum() for the sum'''
    if not len(p)==len(q): raise SyntaxError('Length of p and q must be the same, instead length of p is ' + len(p) + ' and length of q is ' + len(q))
    else: n = len(p)
    for i in range(n):
        KLi = 0
        for j in ['A','G','T','C']:
            KLi += p[i][j] * math.log( p[i][j] / q[i][j] , 2)
        KL.append(KLi)
    return KL

def addRiToIterable(iterable, p, n=0):
    '''takes an iterable that yields lists (or something else with the + operator defined) and returns a generator that adds (using +) the Ri (information content) of the sequence which must be contained in position n (i.e. x[n] where x is the item in iterable)

the Ri (information content) is calculated using p, which must define the Ri method, such that p.Ri(seq) returns the Ri for the sequence seq.
'''
    return ( x + p.Ri(x[n]) for x in iterable )

# requires numpy, maybe relax requierment?
# needs to return an mmMatrix object
# needs to be able to save to file
# needs to be able to make and save image
def joint_matrix(sites):
    '''takes as input a filename and returns the joint Rate matrix
    for the list of sequences contained in that file

    Joint rates R(X;Y_ are defined as 
    R(X;Y) = - sum over X,Y p(x,y) * I(X;Y)
    I(X;y) = - sum over X,Y p(x,y) * log2[p(x,y)/(p(x)p(y))]
    '''
    bases = ['A','C','G','T']
    indexDictionary = {} # the index dictionary
    for i in range(4):
        for j in range(4):
            ssPair = bases[i] + bases[j]
            indexDictionary[ssPair]=i,j

    siteLength = len(sites[0])
# initialize the matrix
    diCounts = zeros([siteLength,siteLength],dtype='(4,4)int')

    def addSeq(m,s,n,b):
        '''adds the dinucleotide counts from one sequence to the mmMatrix (an array, passed by refence). requires the length n'''
        for i in range(n):
            for j in range(n):
                m[i,j][ b[s[i]+s[j]] ] += 1

# count pairs over every sequence
    for site in sites:
        addSeq(diCounts,site.upper(),siteLength,indexDictionary)

# convert to probabilities
    diProbs = zeros([siteLength,siteLength],dtype='(4,4)float')
    totalSeqs = diCounts[0,0].sum()
    for i in range(siteLength):
        for j in range(siteLength):
                for ii in range(4):
                    for jj in range(4):
                        diProbs[i,j][ii,jj] = diCounts[i,j][ii,jj] / float(totalSeqs)

    mmMatrix = zeros([siteLength,siteLength],dtype='float')
    for i in range(siteLength):
        for j in range(siteLength):
            # sum over all dinucleotide combinations
            pM = diProbs[i,j]

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

            mmMatrix[i][j] = Rij
    return (diCounts, diProbs, mmMatrix)

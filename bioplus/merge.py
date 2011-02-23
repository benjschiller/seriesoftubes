'''tools for merging and sorting'''

def merge(left, right, n=None, order='descending'):
    '''merge takes two sorted lists and merges them so that elements are ordered. if the third parameter n is specified, then the lists are assumed to have list elements and are sorted on the col n (first col. = 0)

order may be ascending or descending (default: descending)
'''
    if order=='descending': return mergedesc(left,right,n)
    elif order=='ascending': return mergeasc(left,right,n)
    else: raise ValueError('merge requires that you specify order=\'descending\' or order=\'ascending\'')

def mergeasc(left,right,n=None):
    '''merge in ascending order, see merge'''
    result = []
    i , j = 0, 0
    if n==None:
        while(i < len(left) and j < len(right)):
                if (left[i] <= right[j]):
                    result.append(left[i])
                    i = i + 1
                else:
                    result.append(right[j])
                    j = j + 1
    else:   
        while(i < len(left) and j < len(right)):
            if (left[i][n] <= right[j][n]):
                result.append(left[i])
                i = i + 1
            else:
                result.append(right[j])
                j = j + 1
    result += left[i:]
    result += right[j:]
    return result

def mergedesc(left,right,n=None):
    '''merge in descending order, see merge'''
    result = []
    i , j = 0, 0
    if n==None:
        while(i < len(left) and j < len(right)):
                if (left[i] >= right[j]):
                    result.append(left[i])
                    i = i + 1
                else:
                    result.append(right[j])
                    j = j + 1
    else:   
        while(i < len(left) and j < len(right)):
            if (left[i][n] >= right[j][n]):
                result.append(left[i])
                i = i + 1
            else:
                result.append(right[j])
                j = j + 1
    result += left[i:]
    result += right[j:]
    return result

def mergesort(list,n=None, order='descending'):
    '''mergesort sorts a single list, by its elements (n=None), or by the nth elelement in each element (n=0,1,2...). The returned order may be ascending (order='ascending') or descending (order='decending', default)'''
    if order=='descending': return mergesortdesc(list,n)
    elif order=='ascending': return mergesortasc(list,n)
    else: raise ValueError('merge requires that you specify order=\'descending\' or order=\'ascending\'')

def mergesortasc(list,n=None):
    '''mergesort in ascending order, see mergesort'''
    if len(list) < 2:
        return list
    else:
        middle = len(list) / 2
        left = mergesortasc(list[:middle],n)
        right = mergesortasc(list[middle:],n)
        return mergeasc(left, right,n)

def mergesortdesc(list,n=None):
    '''mergesort in descending order, see mergesort'''
    if len(list) < 2:
        return list
    else:
        middle = len(list) / 2
        left = mergesortdesc(list[:middle],n)
        right = mergesortdesc(list[middle:],n)
        return mergedesc(left, right,n)

def argmergesort(L):
    '''argmergesort sorts a list of depth 0 and returns the arguments that sort it'''
    K = mergesort( [[i,L[i]] for i in range(len(L))] , 1 )
    return [x[0] for x in K]

def mergeLists(a,b):
    '''mergeList(a,b) appends each item in b to a, only if that item is not already present in a'''
    for x in b:
        if a.count(x)==0: a.append(x)
    return a

def dictSort(d):
    '''sorts a dictionary according to its values. returns a list (dictionaries are sorted on keys, not values)'''
    k = d.keys()
    v = d.values()
    s = argmergesort(v)
    print s
    e  = []
    for i in s:
        e.append( [ k[i] , v[i] ] )
    return e

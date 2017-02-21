import numpy as np

from scipy.sparse import csr_matrix
from scipy.spatial import Delaunay
from scipy.stats import binom
from scipy.stats import binom_test
from scipy.stats import gamma, norm



def ppa_agree(locs,y):
    tri = Delaunay(locs)
    indptr,indices = tri.vertex_neighbor_vertices
    adj = csr_matrix((np.ones(len(indices)),indices,indptr), shape=(len(locs), len(locs)))

    # convert to one-hot matrix
    if len(y.shape)==1:
        ymat = np.zeros((len(y),max(y)+1))
        ymat[np.arange(len(y)),y] = 1
    else:
        ymat = y
    
    tot_agree = np.dot((ymat.T*adj),ymat).trace()
    py = ymat.mean(axis=0)
    base_rate = np.trace(np.outer(py,py))        
    return tot_agree,base_rate,adj


def ppa(locs,y,sparse=True,bootstrap_N=200):
    '''
    y can be a 1-D array of binary observations
    or it can be a matrix of categorical observations (one-hot)
    '''
    tot_agree, base_rate, adj = ppa_agree(locs,y)
    tot_edge = adj.sum()

  
    pval = None
    
    boots = []
    Yrand = np.copy(y)
    for _ in xrange(bootstrap_N):
        np.random.shuffle(Yrand) # should work on matrix or list
        boots.append(ppa_agree(locs,Yrand)[0])

    boots = np.array(boots)
    pval = (sum(b >= tot_agree for b in boots) + 1)/float(len(boots) + 1)
        
    return pval, base_rate, adj.sum()

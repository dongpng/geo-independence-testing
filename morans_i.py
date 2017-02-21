import math

import numpy as np
import scipy
import scipy.stats

np.seterr(divide='ignore', invalid='ignore') 

# Cutoffs to use for distance matrices
CUTOFFS_MUNI = [x for x in range(30000,80000,5000)]
CUTOFFS_PROV = [x for x in range(60000,110000,5000)]


def getDistanceMatrixPolygons(polygons):
    """ Return a matrix with distances between polygons 
    polygons: an array of polygon objects
    """
    # Based on centroids. 
    centroids = [poly.centroid for poly in polygons]
    N = len(centroids)
    distance_matrix = np.zeros((N,N))
    for i in xrange(N):
        for j in xrange(i + 1, N):
            distance_matrix[i,j] = centroids[i].distance(centroids[j])
            distance_matrix[j,i] = distance_matrix[i,j]
    return distance_matrix


def getWeightMatrix(distance_matrix, cutoff, normalize=False):
    """
    Applies a cutoff to generate a weight matrix
    Row-normalized and diagonals are set to zero.
    
    distance matrix: matrix of distances between polygons
    cutoff: the cutoff to apply
    """
    result = 1 * (distance_matrix < cutoff)
    np.fill_diagonal(result, 0)
    if normalize:
        row_sums = result.sum(axis=1, dtype='float')
        result = result / row_sums[:, np.newaxis]
    return result


def moransI_fast(freq_bins, weights_matrix, mu=None, tot_weight=None, bin_var=None):
    """
    A faster version of Moran's I.

    freq_bins: array of frequencies
    weights_matrix: matrix of weights between bins
    """
    if mu is None:
        mu = freq_bins.mean()
    if tot_weight is None:
        tot_weight = weights_matrix.sum()
    if bin_var is None:
        bin_var = freq_bins.var()

    numerator = 0.
    N = freq_bins.size
    
    for i in range(N):
        numerator += np.dot(weights_matrix[i,:],(freq_bins-mu)) * (freq_bins[i] - mu)       
    return (1. / tot_weight) * (numerator / bin_var)


def moransI(freq_bins, weights_matrix):
    """
    Calculates Moran's I.
    In Grieve et al, the spatial weighting function is a simple indicator function,
    with value=1 if the distance is less than the cutoff, and zero otherwise.

    freq_bins: array of frequencies
    weights_matrix: matrix of weights between bins
    """
    mu = freq_bins.mean()
    tot_weight = 0.
    numerator = 0.
    N = freq_bins.size
    bin_var = freq_bins.var()
    for i in range(N):
        for j in range(N):
            val = weights_matrix[i][j] * (freq_bins[i] - mu) * (freq_bins[j] - mu)
            tot_weight += weights_matrix[i][j]
            numerator += val
            
    morans_I = (1. / tot_weight) * (numerator / bin_var)
    return morans_I


def morans_I_pval_boots(bins, weight_matrix, N_samp=500):
    """
    Return Morans'I, Z-score and p-value
    bins: array of frequencies
    weight_matrix: matrix of weight between the bins
    N_samp: Number of samples
    """
    mu = bins.mean()
    bin_var = bins.var()
    tot_weight = weight_matrix.sum()

        
    mI = moransI_fast(bins, weight_matrix, mu, tot_weight, bin_var)
    
    boots = []
    fb_rand = bins.copy()
    for _ in xrange(N_samp):
        np.random.shuffle(fb_rand)
        boots.append(moransI_fast(fb_rand, weight_matrix, mu, tot_weight, bin_var))
    boots = np.array(boots)
    
    pval = (sum(b >= mI for b in boots) + 1)/float(len(boots) + 1)
    
    return mI,pval



def morans_I_pval_cutoffsweep(bins, distance_matrix, cutoffs, N_samp=500):
    """ Calculate Moran's I for a range of cutoff values
    
    Return the Moran's I values, the p-values, and values for best p-value 
    
    bins: arra of frequencies
    distance_matrix: the distance matrix
    N_samp: Number of samples
    """
    # Init
    MoransI_vals = []
    MoransI_pvals = []
    MoransI_zScores = []
    best_Morans_I_val = None
    best_pval = float("inf")
    best_cutoff = None
    
    for cutoff in cutoffs:
        weight_matrix_c = getWeightMatrix(distance_matrix, cutoff)
        mI,pval = morans_I_pval_boots(bins, weight_matrix_c, N_samp)
    
        MoransI_vals.append(mI)
        MoransI_pvals.append(pval)
        MoransI_zScores.append(z)
        if pval < best_pval:
            best_pval = pval
            best_Morans_I_val = mI
            best_cutoff = cutoff

    #None are z-values, not used anymore
    return best_Morans_I_val, best_pval, None, best_cutoff, MoransI_vals, MoransI_pvals, MoransI_zScores



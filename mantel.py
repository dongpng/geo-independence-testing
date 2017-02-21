import math
import time

import matplotlib.pyplot as plt
import numpy as np
import scipy

from numpy.random import RandomState
from scipy.spatial.distance import cdist, pdist, squareform
from scipy.stats import gamma

import HSIC

SIGMAS_MANTEL = [x for x in range(35000,85000,5000)]


def matrix_mantel(X, method="Delta", sigma=None):
    """ Returns an array of the upper triangle of the distance/similarity matrix between labels
    
    Y: Array with labels
    method: Delta, Euclidian
    """

    if method == "Delta":
        return pdist(X, 'hamming')
    elif method == "Euclidian":
        return pdist(X, 'euclidean')
    elif method == "Gaussian":
        return HSIC.kernelMatrixGaussian(X,X, sigma)[np.triu_indices(X.shape[0], 1)]
    elif method == "DeltaSim":
        return 1- pdist(X, 'hamming')
    elif method == "AbsDiffFreq":
        return pdist(X, 'chebyshev')
    elif method == "Euclidian_threshold":
        result = pdist(X, 'euclidean')
        mdn = np.median(result)
        result[result < mdn] = 0
        result[result >= mdn] = 1
        return result


def mantel_pval(X, Y, kX="Euclidian", kY="Delta", sigmaX=None, sigmaY=None, 
                N_samp=500, random_seed=None, return_boots=False):
    """ Calculates the Mantel test 
    A faster method for calculating the p-values
    X: Array of locations
    Y: Array of observations 
    
    kX: method for calculating the X matrix
    kY: method for calculating the Y matrix
    sigmaX: param for the Gaussian kernel
    sigmaY: param for the Gaussian kernel
    N_samp: Number of samples for bootstrap 
    random_seed: for calculating p values
    return_boots: return bootstrap vaules?
    """
    prng = RandomState(random_seed)
    
    # Calculate two distance matrices
    if not sigmaX and kX == "Gaussian":
        sigmaX = HSIC.getSigmaGaussian(X,X,200)
    if not sigmaY and kY == "Gaussian":
        sigmaY = HSIC.getSigmaGaussian(Y,Y,200)
    
    # Calculate two distance matrices
    A = matrix_mantel(X, kX, sigmaX)
    B = matrix_mantel(Y, kY, sigmaY)
    
    # Just pearson correlation
    A_minus_mean = A-A.mean()
    B_mean = B.mean()
    B_minus_mean = B-B_mean
    lright = math.sqrt((A_minus_mean**2).sum())
    lleft = math.sqrt(((B_minus_mean)**2).sum())
    top = A_minus_mean.dot(B_minus_mean)
    mantel_val = top/(lright * lleft)
    
    # Calculating the p-value
    pval = 1.0
    Yrand = np.copy(Y)
    boots = []
    for i in xrange(N_samp):
        prng.shuffle(Yrand)
        B_tmp = matrix_mantel(Yrand, kY, sigmaY)
        
        if return_boots:
            boots.append(A_minus_mean.dot(B_tmp-B_mean))
        if A_minus_mean.dot(B_tmp-B_mean) >= top: #can use B_mean instead of B_tmp.mean()
            pval += 1
            
    pval /= N_samp + 1
    
    if return_boots:
        return mantel_val, pval, boots
    else:
        return mantel_val, pval
    


def mantel_pval_old(X, Y, kX="Euclidian", kY="Delta", N_samp=500, random_seed=None):
    """ Calculates the Mantel test 
    
    X: Array of locations
    Y: Array of observations 
    N_samp: Number of samples for bootstrap 
    """
    prng = RandomState(random_seed)
    
    A = matrix_mantel(X, kX)
    B = matrix_mantel(Y, kY)
    
    # Just pearson correlation
    mantel_val = scipy.stats.pearsonr(A,B)[0]
    
    # Calculating the p-value
    boots = []
    Yrand = np.copy(Y)
    for _ in xrange(N_samp):
        prng.shuffle(Yrand)
        boots.append(scipy.stats.pearsonr(A, matrix_mantel(Yrand, kY))[0])
    
    boots = np.array(boots)
    pval = (sum(b >= mantel_val for b in boots) + 1)/float(len(boots) + 1)
    return mantel_val, pval


def mantel_pval_bandwidth_sweep(X, Y, kernelX="Gaussian", kernelY="Delta", N_samp=500, random_seed=None):
    """" Calculate mantel by sweeping over bandwidth values """
    mantel_vals = []
    mantel_pvals = []
    best_mantel_val = None
    best_pval = float("inf")
    for s in SIGMAS_MANTEL:
        val, pval = mantel_pval(X, Y, kX=kernelX, kY=kernelY, sigmaX=s, 
                                N_samp=N_samp, random_seed=random_seed)
        mantel_vals.append(val)
        mantel_pvals.append(pval)
        
        if pval < best_pval:
            best_pval = pval
            best_mantel_val = val

    return best_mantel_val, best_pval, mantel_vals, mantel_pvals


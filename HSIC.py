import math
import random
import time

import numpy as np
import scipy

from numpy.random import RandomState
from scipy.spatial.distance import cdist, pdist, squareform
from scipy.stats import gamma, norm


SIGMAS_HSIC = [x for x in range(35000,85000,5000)]


def kernelMatrixGaussian(m, m2, sigma=None):
    """
    Calculates kernel matrix with a Gaussian kernel
    
    m: rows are data points
    m2: rows are data points
    sigma: the bandwidth of the Gaussian kernel. 
           If not provided, the median distance between points will be used.
         
    """
    
    pairwise_distances = cdist(m, m2, 'sqeuclidean')

    # If sigma is not provided, set sigma based on median distance heuristic.
    if sigma is None:
        sigma = math.sqrt(0.5 * np.median(pairwise_distances[pairwise_distances>0])) 
    gamma = -1.0/(2 * sigma**2)
    return np.exp(gamma * pairwise_distances)

def columnDistanceGaussian(col1, col2, sigma):
    gamma = -1.0/(2 * sigma**2)
    
    result = np.array([scipy.spatial.distance.sqeuclidean(x,y) for x,y in zip(col1, col2)])
    return np.exp(gamma * np.array(result))

def columnDistanceLinear(col1, col2):  
    return np.array([np.dot(x,y) for x,y in zip(col1, col2)])


def getSigmaGaussian(m, m2, sample_size = 200, sigma_multiply=0): 
    """ Calculate sigma for a gaussian kernel based on observations.
    
    m: rows are data points
    m2: rows are data points
    sample_size: maximum number of observations to take into account. 
                 If number of observation is larger, take random sample
    """
    if m.shape[0] > sample_size:
        prng = RandomState(m.shape[0]) # To have the same bandwidth for the same samples
        ind = prng.choice(m.shape[0], sample_size)
        m = m[ind]
        m2 = m2[ind]
    pairwise_distances = cdist(m, m2, 'sqeuclidean')
    
    distance_result = np.median(pairwise_distances[pairwise_distances>0])
    if sigma_multiply != 0:
        distance_result += sigma_multiply * np.std(pairwise_distances[pairwise_distances>0])
    return math.sqrt(0.5 * distance_result) 

def kernelMatrixLinear(m, m2):
    """ Calculates kernel matrix with a linear kernel
    
    m: rows are data points
    m2: rows are data points   
    """
    return np.dot(m, m2.T)

def kernelMatrixDelta(m, m2):
    """ 1 if items are the same. 0 otherwise """
    return 1 - cdist(m, m2, 'hamming')


def columnDistanceDelta(col1, col2):  
    return np.array([1 if x==y else 0 for x,y in zip(col1, col2)])

def HSIC_pval_old(X, Y, N_samp=100, kernelX="Gaussian", kernelY="Gaussian", sigmaX=None, sigmaY=None):
    """ Calculates HSIC and p-value 
    
    Old implementation that calculates complete Gramm matrices
    
    X: Data. Each row is a datapoint.
    Y: Data. Each row is a datapoint.
    N_samp: Number of samples
    kernelX: Kernel to use (Gaussian or Linear)
    kernelY: Kernel to use (Gaussian or Linear)
    """
    timeA = time.time()
    m,_ = X.shape
    
    # Calculate Gram matrices
    sigmaX = getSigmaGaussian(X,X,200) if sigmaX is None else sigmaX
    sigmaY = getSigmaGaussian(Y,Y,200) if sigmaY is None else sigmaY
    K = kernelMatrixGaussian(X,X,sigmaX) if kernelX == "Gaussian" else kernelMatrixLinear(X,X)
    L = kernelMatrixGaussian(Y,Y,sigmaY) if kernelY == "Gaussian" else kernelMatrixLinear(Y,Y)
    
    # Centering matrix
    H = np.identity(m) - 1.0/m * np.ones((m,m)) 
    Kc = np.mat(H) * np.mat(K)*np.mat(H)
    
    # Dividing by m here, although some papers use m-1
    HSIC = np.trace(np.dot(np.transpose(Kc),L))/m**2 
    
    boots = []
    Yrand = np.copy(Y)
    for _ in xrange(N_samp):
        np.random.shuffle(Yrand)
        L = kernelMatrixGaussian(Yrand,Yrand) if kernelY == "Gaussian" else kernelMatrixLinear(Yrand,Yrand)
        boots.append(np.trace(np.dot(np.transpose(Kc),L))/m**2)
    
    boots = np.array(boots)
    pval = (sum(b >= HSIC for b in boots) + 1)/float(len(boots) + 1)
    return HSIC, pval



def HSIC_pval(X, Y, N_samp=500, kernelX="Gaussian", kernelY="Gaussian", eta = 0.001, 
              sigmaX=None, sigmaY=None,
              p_method="boots", return_boots=False):
    """ Calculates HSIC and p-value 
    
    Gram matrices are approximated using incomplete Cholesky decomposition.
    
    X: Data. Each row is a datapoint.
    Y: Data. Each row is a datapoint.
    N_samp: Number of samples
    kernelX: Kernel to use (Gaussian, Linear, Delta)
    kernelY: Kernel to use (Gaussian, Linear, Delta)
    eta: Threshold for incomplete Cholesky decomposition
    sigmaX: sigma for X when using Gaussian kernel
    sigmaY: sigma for Y when using Gaussian kernel
    """
    timeA = time.time()
    m,_ = X.shape
    
    sigmaX = getSigmaGaussian(X,X,200) if sigmaX is None else sigmaX
    sigmaY = getSigmaGaussian(Y,Y,200) if sigmaY is None else sigmaY
    
    A,max_rankA = incompleteCholeskyKernel(X, m, kernelX, sigmaX, eta)
    B,max_rankB = incompleteCholeskyKernel(Y, m, kernelY, sigmaY, eta)
    
    centered_A = A.T - A.T.mean(axis=0)
    tmp = B * np.mat(centered_A)
    HSIC = np.trace(tmp * tmp.T)/m**2
    
    boots = []
    Yrand = np.copy(Y)
    for _ in xrange(N_samp):
        np.random.shuffle(Yrand)
        
        B, max_rankB = incompleteCholeskyKernel(Yrand, m, kernelY, sigmaY, eta)
        
        tmp = np.mat(B) * np.mat(centered_A)
        boots.append(np.trace(tmp * tmp.T)/m**2)

    boots = np.array(boots)
    
    if p_method == "boots":
        pval = (sum(b >= HSIC for b in boots) + 1)/float(len(boots) + 1)
    else: #gamma
        fit_alpha, fit_loc, fit_beta= gamma.fit(boots)
        pval = 1 - gamma.cdf(HSIC, fit_alpha, scale=fit_beta, loc=fit_loc)
                
    if return_boots:
        return HSIC, pval, boots
    else:
        return HSIC, pval



def HSIC_pval_bandwidth_sweep(locs, has_word, N_samp=500, kernelX="Gaussian", kernelY="Gaussian", eta=0.001):
    """" Calculate HSIC by sweeping over bandwidth values """
    HSIC_vals = []
    HSIC_pvals = []
    best_HSIC_val = None
    best_pval = float("inf")
    for s in SIGMAS_HSIC:
        val,pval = HSIC_pval(locs, has_word, N_samp, kernelX, kernelY, eta, sigmaX=s)
        HSIC_vals.append(val)
        HSIC_pvals.append(pval)
        
        if pval < best_pval:
            best_pval = pval
            best_HSIC_val = val

    return best_HSIC_val, best_pval, HSIC_vals, HSIC_pvals



def incompleteCholesky(K, k, eta = 0.01):
    """ Incomplete Cholesky decomposition 
    Based on algorithm in Kernel Methods for Pattern Analysis, chapter
    Elementary algorithms in feature space, fragment 5.4 
    
    K: the matrix
    k: numbers of rows for new matrix
    eta: threshold
    """
    ell,_ = K.shape
    I = []
    R = np.zeros((ell,ell))
    d = np.diagonal(K).copy() 
    a = max(d)
    I.append(np.argmax(d))
    j = 0
    while a > eta and j < k:
        nu_j = math.sqrt(a)
        for i in xrange(ell):
            R[j,i] = (K[I[j],i] - np.dot(R[:,i].T,R[:,I[j]]))/nu_j
        d = d - R[j,:]**2
        a = max(d)
        I.append(np.argmax(d))
        j += 1
        
    return R[:j,], j


def incompleteCholeskyKernel(X, maxrank, kernel, sigma = None, eta = 0.001):
    """ Incomplete Cholesky decomposition 
    Based on algorithm in Kernel Methods for Pattern Analysis, chapter
    Elementary algorithms in feature space, fragment 5.4.
    Doesn't need to compute Gram matrix beforehand.
    
    K: the matrix
    k: numbers of rows for new matrix
    kernel: kernel to use
    sigma: in case of Gaussian kernel
    eta: threshold
    """
    maxrank = min(maxrank, 100)
    ell,_ = X.shape
    I = []
    R = np.zeros((maxrank,ell))
    
    d = None
    if kernel == "Gaussian":
        d = columnDistanceGaussian(X, X, sigma)
    elif kernel == "Linear":
        d = columnDistanceLinear(X,X)
    elif kernel == "Delta":
        d = columnDistanceDelta(X,X)
        
    a = max(d)
    I.append(np.argmax(d))
    j = 0
    while j < maxrank and a > eta:
        nu_j = math.sqrt(a)
        x_elem = np.atleast_2d(X[I[j]])
        
        K_tmp = None
        if kernel == "Gaussian":
            K_tmp = kernelMatrixGaussian(x_elem, X, sigma)
        elif kernel == "Linear":
            K_tmp = kernelMatrixLinear(x_elem, X)
        elif kernel == "Delta":
            K_tmp = kernelMatrixDelta(x_elem, X)
            

        for i in xrange(ell):
            R[j,i] = (K_tmp[0][i] - np.dot(R[:,i].T,R[:,I[j]]))/nu_j
        d = d - R[j,:]**2
        a = max(d)
        I.append(np.argmax(d))
        j += 1
        
    return R[:j,], j




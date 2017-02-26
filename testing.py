import math
import numpy as np
import random
import unittest

import HSIC
import mantel
import morans_i

class MyTest(unittest.TestCase):


    ##### Moran's I ######
    def test_morans_i(self):
        weights_matrix = [[0,1,1,1,0,0], [1,0,0,1,1,0], [1,0,0,1,0,0], [1,1,1,0,1,1], [0,1,0,1,0,1], [0,0,0,1,1,0]]
        freqs = np.asarray ([14,12,11,25,28,30])
        
        self.assertAlmostEqual(0.1676, morans_i.moransI_fast(np.asarray(freqs), np.asarray(weights_matrix)), places=4)


    ##### Mantel ######
    def test_mantel_matrices(self):
        # Testing the label matrix (distance matrices)
        Y2 = np.array([1,-1,1])[np.newaxis, :].T
        label_matrix =  mantel.matrix_mantel(Y2, "Delta")
        expected_upper = [1.,0,1.]
        np.testing.assert_array_equal(label_matrix, expected_upper)

        Y3 = np.array([2,1,3,2])[np.newaxis, :].T
        label_matrix =  mantel.matrix_mantel(Y3, "Delta")
        expected_upper = [1.,1.,0.,1.,1.,1.]
        np.testing.assert_array_equal(label_matrix, expected_upper)


    def test_mantel(self):
        X = np.array([[1,0],[1,1],[1,0],[1,1],[5,4],[5,5]]*1)
        Y = np.array([1,1,1,1,-1,-1]*1)[np.newaxis, :].T
        self.assertEqual(0.9, math.floor(mantel.mantel_pval(X,Y)[0]*10)/10)
        self.assertEqual(0.0, math.floor(mantel.mantel_pval(X,Y)[1]*10)/10)


    def test_mantel_pval(self):
        X = np.random.randn(100,2)
        Y = np.random.choice([0,1], size=(100,1))
        mantel1,pval1 = mantel.mantel_pval(X, Y, random_seed=235325)
        mantel2,pval2 = mantel.mantel_pval_old(X, Y, random_seed=235325)
        self.assertAlmostEqual(mantel1, mantel2, places=7)
        self.assertEqual(pval1, pval2)

        X = np.random.randn(200,2)
        Y = np.random.choice([0,1], size=(200,1))
        mantel1,pval1 = mantel.mantel_pval(X, Y, random_seed=85842)
        mantel2,pval2 = mantel.mantel_pval_old(X, Y, random_seed=85842)
        self.assertAlmostEqual(mantel1, mantel2, places=7)
        self.assertEqual(pval1, pval2)



    def test_mantel2(self):
        # Case with one point in the center having a different value
        X = np.array([[0,0],[0,1],[0,2],[1,0],[1,1],[1,2],[2,0],[2,1],[2,2]])
        Y = np.array([[0],[0],[0],[0],[1],[0],[0],[0],[0]])
        
        self.assertAlmostEqual(-0.4000349, mantel.mantel_pval(X, Y,)[0], places=7)

        # Triangle
        X = np.array([[0,0],[0,1],[0,2],[1,0],[1,1],[1,2],[2,0],[2,1],[2,2]])
        Y = np.array([[0],[0],[0],[1],[0],[0],[1],[1],[0]])
        self.assertAlmostEqual(0.2971499, mantel.mantel_pval(X, Y)[0], places=7)

        X = np.array([[1],[0.2],[0.5],[-1]])
        Y = np.array([[0],[1],[0],[1]])
        self.assertAlmostEqual(0.241648837, mantel.mantel_pval(X, Y, kX="Euclidian", kY="AbsDiffFreq")[0], places=7)


    ##### HSIC ######
    def test_HSIC(self):
        # Comparing the output values with Arthur Gretton's implementation
        X1 = np.array([1,2,3,4]).reshape((4,1))
        Y1 = np.array([-1,0,-1,0]).reshape((4,1))
        X2 = np.array([[1,0],[2,1],[3,5],[4,10]])
        Y2 = np.array([[-1,3],[0,6],[-1,3],[0,-1]])
        X3 = np.array([[1,0,8],[2,1,2],[3,5,0],[4,10,5]])
        Y3 = np.array([[-1,3,6],[0,6,9],[-1,3,3],[0,-1,-1]])
        self.assertAlmostEqual(0.10755, HSIC.HSIC_pval(X1,X1,5,"Gaussian","Gaussian")[0], places=5)
        self.assertAlmostEqual(0.014440, HSIC.HSIC_pval(X1,Y1,5,"Gaussian","Gaussian")[0], places=5)
        self.assertAlmostEqual(0.072315, HSIC.HSIC_pval(X2,Y2,5,"Gaussian","Gaussian")[0], places=5)
        self.assertAlmostEqual(0.069226, HSIC.HSIC_pval(X3,Y3,5,"Gaussian","Gaussian")[0], places=5)
        self.assertAlmostEqual(0.10755, HSIC.HSIC_pval_full_gram(X1,X1,5,"Gaussian","Gaussian")[0], places=5)
        self.assertAlmostEqual(0.014440, HSIC.HSIC_pval_full_gram(X1,Y1,5,"Gaussian","Gaussian")[0], places=5)
        self.assertAlmostEqual(0.072315, HSIC.HSIC_pval_full_gram(X2,Y2,5,"Gaussian","Gaussian")[0], places=5)
        self.assertAlmostEqual(0.069226, HSIC.HSIC_pval_full_gram(X3,Y3,5,"Gaussian","Gaussian")[0], places=5)

    def test_Cholesky(self):
        A = np.array([[8,2,3], [2,1,5], [3,5,4]], dtype=np.float)
        expected1 = np.array([[2.8284e+00,7.0711e-01,1.0607e+00],
                              [2.6191e-16,2.5065e+00,1.6956e+00]])
        np.testing.assert_array_almost_equal(expected1, HSIC.incompleteCholesky(A,2)[0], decimal=4)

        B = np.array([[5,3,4,1], [1,9,4,2] ,[2,7,4,4] , [1,0,0,1]], dtype=np.float)
        expected2 = np.array([[0.33333, 3.00000, 1.33333, 0.66667],
                              [2.21108, 0.90453, 1.60806, 0.35176]])
        np.testing.assert_array_almost_equal(expected2, HSIC.incompleteCholesky(B,2)[0], decimal=4)

        # Test cholesky without precomputing kernel
        random.seed(2523523)
        X = np.random.rand(100,5)

        # With computing the Gram matrix first and a Gaussian kernel
        K = HSIC.kernelMatrixGaussian(X,X)
        A,_ = HSIC.incompleteCholesky(K, 10)

        # Without
        A2,_ = HSIC.incompleteCholeskyKernel(X, 10, "Gaussian", HSIC.getSigmaGaussian(X,X))

        # Should be equal
        np.testing.assert_array_almost_equal(A, A2, decimal=6)

        Y = np.random.choice([-1,1], size=(100,1))
        K = HSIC.kernelMatrixLinear(Y,Y)
        A,_ = HSIC.incompleteCholesky(K, 10)
        A2,_ = HSIC.incompleteCholeskyKernel(Y, 10, "Linear")
        np.testing.assert_array_almost_equal(A, A2, decimal=6)

        Y = np.random.choice([0,1,2], size=(100,1))
        K = HSIC.kernelMatrixDelta(Y,Y)
        A,_ = HSIC.incompleteCholesky(K, 10)
        A2,_ = HSIC.incompleteCholeskyKernel(Y, 10, "Delta")
        np.testing.assert_array_almost_equal(A, A2, decimal=6)


    def test_kernel_matrix(self):
        Y1 = np.array([1,0,1,0]).reshape((4,1))
        expected1 = np.array([[1,0,1,0],[0,1,0,1],[1,0,1,0],[0,1,0,1]])
        np.testing.assert_array_equal(expected1, HSIC.kernelMatrixDelta(Y1,Y1))

        Y2 = np.array([1,2,3,4]).reshape((4,1))
        expected2 = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
        np.testing.assert_array_equal(expected2, HSIC.kernelMatrixDelta(Y2,Y2))

        np.testing.assert_array_equal(np.array([1,0,0,0]), HSIC.columnDistanceDelta(Y1,Y2))


if __name__ == '__main__':
    unittest.main()
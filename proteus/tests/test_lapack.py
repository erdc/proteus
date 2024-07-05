import math
import pytest
import numpy as np
from proteus import clapack

@pytest.fixture(autouse=True)
def simple_dense_matrix():
    dense_factor = clapack.DenseFactor(3)
    array = np.array([[1.,0.,1.],[2.,1.,2.],[1.,0.,-1.]])
    mat = np.ndarray(shape=(3,3),
                     dtype=float,
                     buffer=array)
    yield dense_factor, mat

def test_denseFactor_init():
    dense_factor = clapack.DenseFactor(2)
    assert dense_factor.n == 2
    assert dense_factor.lu[0][0] == 0.
    assert dense_factor.pivots[0] == 0

def test_blasCopy(simple_dense_matrix):
    dense_factor, mat = simple_dense_matrix
    clapack.blasCopy(3, mat, dense_factor)
    assert np.array_equal(dense_factor.lu, mat)
    
def test_denseFactorPrepare(simple_dense_matrix):
    dense_factor, mat = simple_dense_matrix
    
    clapack.denseFactorPrepare(3,
                                   mat,
                                   dense_factor)

    lu_expected = np.array([[1.,0.,1.],[2.,1.,0.],[1.,0.,-2.]])
    assert np.array_equal(dense_factor.lu,lu_expected)
    pivots_expected = np.array([1,2,3])
    assert np.array_equal(dense_factor.pivots,pivots_expected)

def test_denseFactorSolve(simple_dense_matrix):
    dense_factor, mat = simple_dense_matrix

    clapack.denseFactorPrepare(3,
                                   mat,
                                   dense_factor)

    b = np.array([0.,0.,0.], dtype='float64')
    clapack.denseFactorSolve(3,
                                 mat,
                                 dense_factor,
                                 b)
    assert np.array_equal(b, np.array([0.,0.,0.]))

    b = np.array([1.,1.,1.], dtype='float64')
    clapack.denseFactorSolve(3,
                                 mat,
                                 dense_factor,
                                 b)
    assert np.array_equal(b, np.array([1.,-1.,0.]))

def test_denseCalculateEigenvalues(simple_dense_matrix):
    dense_factor, mat = simple_dense_matrix

    work = np.zeros((3*5),'d')
    eigvals_r = np.zeros((3,),'d')
    eigvals_i = np.zeros((3,),'d')
    rightEigenvectors = np.zeros((3,3),'d')
    leftEigenvectors = np.zeros((3,3),'d')
    
    clapack.denseCalculateEigenvalues(b'V',
                                      b'V',
                                      3,
                                      mat,
                                      3,
                                      eigvals_r,
                                      eigvals_i,
                                      leftEigenvectors,
                                      3,
                                      rightEigenvectors,
                                      3,
                                      work,
                                      15)

    real_eig_vals = np.array([math.sqrt(2), -math.sqrt(2), 1.])
    imag_eig_vals = np.zeros(3)

    left_eig_vectors = np.array([[ 0.14464074,  0.98766876,  0.05991216],
                                 [-0.34919364, -0.4091058 ,  0.84302802],
                                 [ 0.        ,  1.        ,  0.        ]])

    right_eig_vectors = np.array([[ 0.92387953,  0.        ,  0.38268343],
                                  [-0.38268343,  0.        ,  0.92387953],
                                  [-0.93704257,  0.15617376, -0.31234752]])

    assert np.allclose(eigvals_r, real_eig_vals)
    assert np.array_equal(eigvals_i, imag_eig_vals)
    assert np.allclose(leftEigenvectors, left_eig_vectors)
    assert np.allclose(rightEigenvectors, right_eig_vectors)
import math
import pytest
import numpy as np
from proteus import (clapackTest)

@pytest.fixture(autouse=True)
def simple_dense_matrix():
    dense_factor = clapackTest.DenseFactor(3)
    array = np.array([[1.,0.,1.],[2.,1.,2.],[1.,0.,-1.]])
    mat = np.ndarray(shape=(3,3),
                     dtype=float,
                     buffer=array)
    yield dense_factor, mat

def test_denseFactor_init():
    dense_factor = clapackTest.DenseFactor(2)
    assert dense_factor.n == 2
    assert dense_factor.lu[0][0] == 0.
    assert dense_factor.pivots[0] == 0

def test_blasCopy(simple_dense_matrix):
    dense_factor, mat = simple_dense_matrix
    clapackTest.blasCopy(3, mat, dense_factor)
    assert np.array_equal(dense_factor.lu, mat)
    
def test_denseFactorPrepare(simple_dense_matrix):
    dense_factor, mat = simple_dense_matrix
    
    clapackTest.denseFactorPrepare(3,
                                   mat,
                                   dense_factor)

    lu_expected = np.array([[1.,0.,1.],[2.,1.,0.],[1.,0.,-2.]])
    assert np.array_equal(dense_factor.lu,lu_expected)
    pivots_expected = np.array([1,2,3])
    assert np.array_equal(dense_factor.pivots,pivots_expected)

def test_denseFactorSolve(simple_dense_matrix):
    dense_factor, mat = simple_dense_matrix

    clapackTest.denseFactorPrepare(3,
                                   mat,
                                   dense_factor)

    b = np.array([0.,0.,0.], dtype='float64')
    clapackTest.denseFactorSolve(3,
                                 mat,
                                 dense_factor,
                                 b)
    assert np.array_equal(b, np.array([0.,0.,0.]))

    b = np.array([1.,1.,1.], dtype='float64')
    clapackTest.denseFactorSolve(3,
                                 mat,
                                 dense_factor,
                                 b)
    assert np.array_equal(b, np.array([1.,-1.,0.]))

def test_denseCalculateEigenvalues(simple_dense_matrix):
    assert 1 == 0

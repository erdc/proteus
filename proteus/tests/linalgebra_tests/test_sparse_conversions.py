#!/usr/bin/env python
"""
Test module for some linear algebra tools functions.
"""
from proteus.iproteus import *
from proteus import superluWrappers
from proteus import LinearAlgebraTools
from petsc4py import PETSc as p4pyPETSc
import numpy
import os, sys
import pytest

@pytest.mark.LinearAlgebraTools
def test_superlu_2_dense():
    """ Tests the superlu wrapper to dense matrix function. """

    vals    = numpy.array([10.,-2.,3.,9.,3.,7.,8.,7.,3.,
                           8.,7.,5.,8.,9.,9.,13.,4.,2.,-1.])
    col_idx = numpy.array([0, 4 ,0 ,1 ,5 ,1 ,2 ,3 ,0 ,2 ,
                           3 ,4 ,1 ,3 ,4 ,5  ,1 ,4 , 5 ],
                          dtype='int32')
    row_idx = numpy.array([0, 2, 5, 8, 12, 16, 19],
                          dtype='int32')
    size_n  = len(row_idx)-1

    superlu_mat = LinearAlgebraTools.SparseMat(size_n,
                                               size_n,
                                               len(vals),
                                               vals,
                                               col_idx,
                                               row_idx)

    dense_mat = LinearAlgebraTools.superlu_sparse_2_dense(superlu_mat)
    comparison_mat = numpy.loadtxt(os.path.join(os.path.dirname(__file__),
                                                'sparse_mat_1.txt'))
    assert numpy.allclose(dense_mat,comparison_mat)

@pytest.mark.LinearAlgebraTools
def test_sparse_2_dense():
    '''
    This function tests the petsc4py_sparse_2_dense function in the
    LinearAlgebraTools module.
    '''
    from proteus import LinearAlgebraTools
    vals    = [10.,-2.,3.,9.,3.,7.,8.,7.,3.,8.,7.,5.,8.,9.,9.,13.,4.,2.,-1.]
    col_idx = [0  , 4 ,0 ,1 ,5 ,1 ,2 ,3 ,0 ,2 ,3 ,4 ,1 ,3 ,4 ,5  ,1 ,4 , 5 ]
    row_idx = [0, 2, 5, 8, 12, 16, 19]
    size_n  = len(row_idx)-1
    petsc_mat = p4pyPETSc.Mat().createAIJ(size = (size_n,size_n), \
                                          csr  = (row_idx, col_idx, vals))
    dense_mat = LinearAlgebraTools.petsc4py_sparse_2_dense \
                (petsc_mat)
    comparison_mat = numpy.loadtxt(os.path.join(os.path.dirname(__file__),
                                                'sparse_mat_1.txt'))
    assert numpy.allclose(dense_mat,comparison_mat)

@pytest.mark.LinearAlgebraTools
def test_superlu_2_petsc():
    """ Tests the function superlu_2_petsc4py """
    vals    = numpy.array([10.,-2.,3.,9.,3.,7.,8.,7.,3.,
                           8.,7.,5.,8.,9.,9.,13.,4.,2.,-1.])
    col_idx = numpy.array([0, 4 ,0 ,1 ,5 ,1 ,2 ,3 ,0 ,2 ,
                           3 ,4 ,1 ,3 ,4 ,5  ,1 ,4 , 5 ],
                          dtype='int32')
    row_idx = numpy.array([0, 2, 5, 8, 12, 16, 19],
                          dtype='int32')
    size_n  = len(row_idx)-1

    superlu_mat = LinearAlgebraTools.SparseMat(size_n,
                                               size_n,
                                               len(vals),
                                               vals,
                                               col_idx,
                                               row_idx)
    petsc4py_mat = LinearAlgebraTools.superlu_2_petsc4py(superlu_mat)

    dense_mat = LinearAlgebraTools.petsc4py_sparse_2_dense(petsc4py_mat)
    comparison_mat = numpy.loadtxt(os.path.join(os.path.dirname(__file__),
                                                'sparse_mat_1.txt'))
    assert numpy.allclose(dense_mat,comparison_mat)


if __name__ == '__main__':
    test_sparse_2_dense()
    test_superlu_2_dense()
    test_superlu_2_petsc()

#!/usr/bin/env python
"""
Test module for some linear algebra tools functions.
"""
from proteus.iproteus import *
from proteus import superluWrappers
from proteus import LinearAlgebraTools
from petsc4py import PETSc as p4py
import proteus.test_utils.TestTools
import numpy
import os, sys

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
    petsc_mat = p4py.Mat().createAIJ(size = (size_n,size_n), \
                                          csr  = (row_idx, col_idx, vals))
    dense_mat = LinearAlgebraTools.petsc4py_sparse_2_dense \
                (petsc_mat)
    comparison_mat = numpy.loadtxt(os.path.join(os.path.dirname(__file__),
                                                'sparse_mat_1.txt'))
    assert numpy.allclose(dense_mat,comparison_mat)

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

class TestPetscMat(proteus.test_utils.TestTools.BasicTest):

    def test_01(self):
        Qp = numpy.array([[1, 0, 3],
                          [1.5, -2.1, 0],
                          [0, 4., 3]])
        vals_Qp = [1., 3., 1.5, -2.1, 4., 3.]
        col_idx_Qp = [0, 2, 0, 1, 1, 2]
        row_idx_Qp = [0, 2, 4, 6]
        size_n = len(row_idx_Qp)-1
        petsc_Qp = p4py.Mat().createAIJ(size=(size_n,size_n),
                                        csr = (row_idx_Qp,
                                               col_idx_Qp,
                                               vals_Qp))
        proteus_Qp = LinearAlgebraTools.petsc4py_sparse_2_dense(petsc_Qp)
        assert (Qp==proteus_Qp).all()
        
    def test_02(self):
        Fp = numpy.array([[3.2, 1.1, 0.],
                           [6.3, 0. , 1.],
                           [-5.1, 0., 0.]])
        vals_Fp = [3.2, 1.1, 6.3, 1., -5.1]
        col_idx_Fp = [0, 1, 0, 2, 0]
        row_idx_Fp = [0, 2, 4, 5]
        size_n = len(row_idx_Fp)-1
        petsc_Fp = p4py.Mat().createAIJ(size=(size_n,size_n),
                                        csr = (row_idx_Fp,
                                               col_idx_Fp,
                                               vals_Fp))
        proteus_Fp = LinearAlgebraTools.petsc4py_sparse_2_dense(petsc_Fp)
        assert (Fp==proteus_Fp).all()

    def test_03(self):
        Ap = numpy.array([[3., 1., 2.],
                          [7., 4., 0.],
                          [3., 8., 1.1]])
        vals_Ap = [3., 1., 2., 7., 4., 3., 8., 1.1]
        col_idx_Ap = [0, 1, 2, 0, 1, 0, 1, 2]
        row_idx_Ap = [0, 3, 5, 8]
        size_n = len(row_idx_Ap) - 1
        petsc_Ap = p4py.Mat().createAIJ(size=(size_n,size_n),
                                        csr = (row_idx_Ap,
                                               col_idx_Ap,
                                               vals_Ap))
        proteus_Ap = LinearAlgebraTools.petsc4py_sparse_2_dense(petsc_Ap)
        assert (Ap==proteus_Ap).all()
    
if __name__ == '__main__':
    pass

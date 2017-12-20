#!/usr/bin/env python
"""
Test module for PETSc wrapper functions and classes.
"""
from proteus.iproteus import *
from proteus import LinearAlgebraTools as LAT
from petsc4py import PETSc as p4pyPETSc
import numpy
import pytest

@pytest.mark.LinearAlgebraTools
def test_PETScWrapper_Matrix_getSubMatrix():
    '''Tests PETScWrapper_Matrix getSubMatrix method. '''
    vals    = [10.,-2.,3.,9.,3.,7.,8.,7.,3.,8.,7.,5.,8.,9.,9.,13.,4.,2.,-1.]
    col_idx = [0  , 4 ,0 ,1 ,5 ,1 ,2 ,3 ,0 ,2 ,3 ,4 ,1 ,3 ,4 ,5  ,1 ,4 , 5 ]
    row_idx = [0, 2, 5, 8, 12, 16, 19]
    size_n  = len(row_idx)-1
    petsc_mat = LAT.PETScWrapper_Matrix(p4pyPETSc.Mat().createAIJ(size = (size_n,size_n), \
                                                                  csr  = (row_idx, col_idx, vals)))
    idx_set_row = p4pyPETSc.IS().createGeneral([1,2,5])
    idx_set_col = p4pyPETSc.IS().createGeneral([1,2,3])
    A = petsc_mat.getSubMatrix(idx_set_row, idx_set_col)

    vals_sub_mat = [9., 7., 8., 7., 4.]
    col_idx_sub_mat = [0, 0, 1, 2, 0]
    row_idx_sub_mat = [0, 1, 4, 5]

    assert numpy.array_equal(row_idx_sub_mat, A.PETSc_Mat.getValuesCSR()[0])    
    assert numpy.array_equal(col_idx_sub_mat, A.PETSc_Mat.getValuesCSR()[1])
    assert numpy.array_equal(vals_sub_mat, A.PETSc_Mat.getValuesCSR()[2])

if __name__ == '__main__':
    pass

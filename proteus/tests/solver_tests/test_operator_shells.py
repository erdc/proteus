#!/usr/bin/env python
"""

Test module for the convection-diffusion operator.

"""
import proteus.test_utils.TestTools
from proteus.iproteus import *
from proteus import Comm
from proteus import LinearAlgebraTools

import os
import sys
import inspect
import petsc4py
import numpy as np
import pytest

from petsc4py import PETSc as p4pyPETSc
from scipy.sparse import csr_matrix

proteus.test_utils.TestTools.addSubFolders( inspect.currentframe() )

class TestOperatorShells(proteus.test_utils.TestTools.BasicTest):

    def test_pcd_shell(self):
        '''  Tests the pcd_shell operators produce correct output. '''
        vals_Qp = [1., 3., 1.5, -2.1, 4., 3.]
        col_idx_Qp = [0, 2, 0, 1, 1, 2]
        row_idx_Qp = [0, 2, 4, 6]
        vals_Fp = [3.2, 1.1, 6.3, 1., -5.1]
        col_idx_Fp = [0, 1, 0, 2, 0]
        row_idx_Fp = [0, 2, 4, 5]
        vals_Ap = [3., 1., 2., 7., 4., 3., 8., 1.1]
        col_idx_Ap = [0, 1, 2, 0, 1, 0, 1, 2]
        row_idx_Ap = [0, 3, 5, 8]
        size_n = len(row_idx_Qp) - 1
        petsc_matQp = p4pyPETSc.Mat().createAIJ(size = (size_n,size_n), 
                                                csr = (row_idx_Qp,col_idx_Qp,vals_Qp))
        petsc_matFp = p4pyPETSc.Mat().createAIJ(size = (size_n,size_n), 
                                                csr = (row_idx_Fp,col_idx_Fp,vals_Fp))
        petsc_matAp = p4pyPETSc.Mat().createAIJ(size = (size_n,size_n), 
                                                csr = (row_idx_Ap,col_idx_Ap,vals_Ap))
        PCD_shell = LinearAlgebraTools.PCDInv_shell(petsc_matQp,petsc_matFp,petsc_matAp)
        x_vec = np.ones(size_n)
        y_vec = np.zeros(size_n)
        x_PETSc_vec = p4pyPETSc.Vec().createWithArray(x_vec)
        y_PETSc_vec = p4pyPETSc.Vec().createWithArray(y_vec)
        A = None
        PCD_shell.apply(A,x_PETSc_vec,y_PETSc_vec)
        true_vec = np.mat('[0.52511998  -0.13008364  -0.03200969]')
        assert np.allclose(y_vec,true_vec)

    def test_BAinvBt_shell(self):
        ''' Test the B_Ainv_Bt_shell '''
        vals_A =    [5.5,3.6,1.8,7.1,7.9,2.1,1.0]
        col_idx_A = [0  ,1  ,2  ,1  ,2  ,0  ,2 ]
        row_idx_A = [0, 3, 5, 7]
        vals_B    =  [1.1, 6.3, 7.3, 3.6, 6.3, 2.4, 6.2, 1.2, 1.0]
        col_idx_B =  [0  , 2  , 0  , 1  , 2  , 1  , 2  , 0  , 2  ]
        row_idx_B =  [0, 2, 5, 7, 9]
        num_B_rows = len(row_idx_B)-1
        num_B_cols = len(row_idx_A)-1
        petsc_matA = p4pyPETSc.Mat().createAIJ(size = (num_B_cols,num_B_cols) ,
                                               csr = (row_idx_A, col_idx_A, vals_A))
        petsc_matB = p4pyPETSc.Mat().createAIJ(size = (num_B_rows,num_B_cols) ,
                                               csr = (row_idx_B, col_idx_B, vals_B))
        BinvABt_shell = LinearAlgebraTools.B_Ainv_Bt_shell(petsc_matA,petsc_matB)
        x_vec = np.ones(num_B_rows)
        y_vec = np.zeros(num_B_rows)
        x_PETSc_vec = p4pyPETSc.Vec().createWithArray(x_vec)
        y_PETSc_vec = p4pyPETSc.Vec().createWithArray(y_vec)
        A = None
        BinvABt_shell.mult(A,x_PETSc_vec,y_PETSc_vec)
        true_solu = np.mat('[64.60463468 60.7745635 35.16771019 15.33818394]')
        assert np.allclose(y_vec,true_solu)

    def test_lsc_shell(self):
        ''' Test for the lsc operator shell '''
        vals_A =     [5.5,7.1,1.0]
        col_idx_A =  [0 , 1 , 2  ]
        row_idx_A =  [0, 1, 2, 3]
        vals_B    =  [1.1, 6.3, 7.3, 3.6, 6.3]
        col_idx_B =  [0  , 2  , 0  , 1  , 2  ]
        row_idx_B =  [0, 2, 5]
        vals_F  =    [3.2, 1.1, 6.3, 1., -5.1]
        col_idx_F  = [0, 1, 0, 2, 0]
        row_idx_F  = [0, 2, 4, 5]
        num_B_rows = len(row_idx_B) - 1
        num_B_cols = len(row_idx_A) - 1
        petsc_matA = p4pyPETSc.Mat().createAIJ(size = (num_B_cols,num_B_cols), 
                                                csr = (row_idx_A,col_idx_A,vals_A))
        petsc_matB = p4pyPETSc.Mat().createAIJ(size = (num_B_rows,num_B_cols), 
                                                csr = (row_idx_B,col_idx_B,vals_B))
        petsc_matF = p4pyPETSc.Mat().createAIJ(size = (num_B_cols,num_B_cols), 
                                                csr = (row_idx_F,col_idx_F,vals_F))
        LSC_shell = LinearAlgebraTools.LSCInv_shell(petsc_matA,petsc_matB,petsc_matF)
        x_vec = np.ones(num_B_rows)
        y_vec = np.zeros(num_B_rows)
        x_PETSc_vec = p4pyPETSc.Vec().createWithArray(x_vec)
        y_PETSc_vec = p4pyPETSc.Vec().createWithArray(y_vec)
        A = None
        LSC_shell.apply(A,x_PETSc_vec,y_PETSc_vec)
        true_solu = np.mat('[-0.01096996,0.00983216]')
        assert np.allclose(y_vec,true_solu)

if __name__ == '__main__':
    pass

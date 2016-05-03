#!/usr/bin/env python
"""

Test module for the convection-diffusion operator.

"""
from proteus.iproteus import *
from proteus import Comm
from proteus import LinearAlgebraTools
comm = Comm.get()
Profiling.logLevel=7
Profiling.verbose=True
import numpy.testing as npt
from nose.tools import ok_ as ok
from nose.tools import eq_ as eq
from nose.tools import set_trace
from petsc4py import PETSc as p4pyPETSc

from scipy.sparse import csr_matrix
import numpy as np


def test_pcd_shell():
    '''
    Tests the pcd_shell operators produce correct output.
    '''
    from proteus import LinearAlgebraTools
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
    petsc_matQp = p4pyPETSc.Mat().createAIJ(size = (size_n,size_n), \
                                            csr = (row_idx_Qp,col_idx_Qp,vals_Qp))
    petsc_matFp = p4pyPETSc.Mat().createAIJ(size = (size_n,size_n), \
                                            csr = (row_idx_Fp,col_idx_Fp,vals_Fp))
    petsc_matAp = p4pyPETSc.Mat().createAIJ(size = (size_n,size_n), \
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

def test_pcd_matrix():
    '''
    Test construction of the PCD operators.
    '''

if __name__ == '__main__':
    from proteus import Comm
    comm = Comm.init()
    test_pcd_shell()

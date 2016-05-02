#!/usr/bin/env python
"""

Test module for 2D Quadrilateral Meshes

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
                (petsc_mat.getValuesCSR())
    comparison_mat = np.loadtxt('sparse_mat_1.txt')
    assert np.allclose(dense_mat,comparison_mat)

def test_Qp_mat():
    ''' Verify that Qp returns the correct pressure mass matrix  '''
    import stokes_2d_p
    import stokes_2d_n
    pList = [stokes_2d_p]
    nList = [stokes_2d_n]    
    so = default_so
    so.tnList = [0.,1.]
    so.name = pList[0].name
    so.sList=[default_s]
    opts.verbose=True
    opts.profile=True
    ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
    # *** 1 *** : test the pressure mass matrix is being created properly.
    smoother = LinearSolvers.NavierStokes3D_Qp(L=ns.modelList[0].par_jacobianList[1])
    operator_constructor = LinearSolvers.schurOperatorConstructor(smoother, 'stokes')
    Qp_raw = operator_constructor.getQp()
    Qp_dense = LinearAlgebraTools.petsc4py_sparse_2_dense(Qp_raw.getValuesCSR())
    pressure_mass_matrix = np.loadtxt('pressure_mass_matrix.txt')
    assert np.allclose(pressure_mass_matrix,Qp_dense)
    # *** 2 *** : test solver does not generate an error
    ns.calculateSolution('test_Qp_mat')

if __name__ == '__main__':
    from proteus import Comm
    comm = Comm.init()
    test_sparse_2_dense()
    test_Qp_mat()

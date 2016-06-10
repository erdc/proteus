#!/usr/bin/env python
"""

Test module for the pressure mass matrix schur complement

"""
import os
import tables
import numpy as np
from proteus.iproteus import *
from proteus import Comm
from proteus import LinearAlgebraTools
comm = Comm.get()
#Profiling.logLevel=7
#Profiling.verbose=True
import numpy.testing as npt
from nose.tools import ok_ as ok
from nose.tools import eq_ as eq
from nose.tools import set_trace
from petsc4py import PETSc as p4pyPETSc

from scipy.sparse import csr_matrix

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
    comparison_mat = np.loadtxt('sparse_mat_1.txt')
    assert np.allclose(dense_mat,comparison_mat)

def test_pressure_mass_matrix_shell():
    ''' This function tests pressure mass matrix shells '''
    from proteus import LinearAlgebraTools
    vals = [1.,3.,1.5,-2.1,4.,3.]
    col_idx = [0, 2, 0, 1, 1, 2]
    row_idx = [0, 2, 4, 6]
    size_n = len(row_idx) - 1
    petsc_mat = p4pyPETSc.Mat().createAIJ(size = (size_n,size_n), \
                                          csr  = (row_idx, col_idx, vals))
    Qp_shell = LinearAlgebraTools.MatrixShell(petsc_mat)
    Qp_inv_shell = LinearAlgebraTools.MatrixInvShell(petsc_mat)
    x_vec = np.ones(size_n)
    y_vec = np.zeros(size_n)
    x_PETSc_vec = p4pyPETSc.Vec().createWithArray(x_vec)
    y_PETSc_vec = p4pyPETSc.Vec().createWithArray(y_vec)
    A = None
    Qp_shell.mult(A, x_PETSc_vec, y_PETSc_vec)
    comparison_vec = np.array([4., -0.6, 7.])
    assert np.allclose(y_vec,comparison_vec)
    x_vec = np.ones(size_n)
    y_vec = np.zeros(size_n)
    x_PETSc_vec = p4pyPETSc.Vec().createWithArray(x_vec)
    y_PETSc_vec = p4pyPETSc.Vec().createWithArray(y_vec)
    Qp_inv_shell.apply(A, x_PETSc_vec, y_PETSc_vec)
    comparison_vec = np.array([1.02564103, 0.25641026, -0.00854701])
    assert np.allclose(y_vec,comparison_vec)

class TestStokesOperatorConstruction():
    """ This class tests the operator construction for a 2D Stokes Poiseulle Flow problem """
    
    def setUp(self):
        import stokes_2d_p
        import stokes_2d_n
        pList = [stokes_2d_p]
        nList = [stokes_2d_n]    
        so = default_so
        so.tnList = [0.,1.]
        so.name = pList[0].name
        so.sList=[default_s]
        opts.verbose=False
        opts.profile=True
        self.ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
        smoother = LinearSolvers.NavierStokes3D_Qp(L=self.ns.modelList[0].par_jacobianList[1])
        self.operator_constructor = LinearSolvers.SchurOperatorConstructor(smoother, 'stokes')
        self.ns.modelList[0].par_jacobianList[1].pde.scale_dt = False

    def tearDown(self):
        if os.path.exists('proteus.log'):
            os.remove('proteus.log')
        else:
            pass
        os.remove('poiseulleFlow.h5')
        os.remove('poiseulleFlow.xmf')

    def test_01_MassMatrix(self):
        ''' Verify that Q returns the correct pressure and velocity mass matrix  '''
        Qp_petsc = self.operator_constructor.getQp()
        Qp_dense = LinearAlgebraTools.petsc4py_sparse_2_dense(Qp_petsc)
        Qv_petsc = self.operator_constructor.getQv()
        Qv_dense = LinearAlgebraTools.petsc4py_sparse_2_dense(Qv_petsc)

        pressure_mass_matrix = np.loadtxt('pressure_mass_matrix.txt')
        velocity_mass_matrix = np.loadtxt('velocity_mass_matrix.txt')

        assert np.allclose(pressure_mass_matrix,Qp_dense)
        assert np.allclose(velocity_mass_matrix,Qv_dense)

    def test_02_PressureLaplace(self):
        ''' Test the pressure Lapacian matrix '''
        # ARB : I don't remember why I need to do this!
        self.ns.modelList[0].levelModelList[1].dphi[(0,0)].dof[:] = 1.0
        self.ns.modelList[0].levelModelList[1].dphi[(1,1)].dof[:] = 1.0
        self.ns.modelList[0].levelModelList[1].dphi[(2,2)].dof[:] = 1.0
        Ap_raw = self.operator_constructor.getAp()
        Ap_dense = LinearAlgebraTools.petsc4py_sparse_2_dense(Ap_raw)
        pressure_laplace_matrix = np.loadtxt('pressure_laplace_matrix.txt')
        assert np.allclose(pressure_laplace_matrix,Ap_dense)

    def test_03_B(self):
        ''' Test the B operator '''
        # TODO (ARB) : Need a specific test for B

    def test_04_BBt(self):
        ''' Test that -(B^{t})' = B '''
        B_petsc = self.operator_constructor.getB()
        B_dense = LinearAlgebraTools.petsc4py_sparse_2_dense(B_petsc)
        Bt_petsc = self.operator_constructor.getBt()
        Bt_dense = LinearAlgebraTools.petsc4py_sparse_2_dense(Bt_petsc)
        assert np.allclose(-1.0*np.transpose(Bt_dense),B_dense)

    def test_05_FullRun(self):
        """ Test the Poiseulle Flow runs and produces expected h5 output """
        self.ns.calculateSolution('test_Qp_mat')
        test_path = os.path.dirname(os.path.abspath(__file__))
        expected = tables.openFile(os.path.join(test_path,
                                   'poiseulleFlow_expected.h5'),
                                   'r')
        actual = tables.openFile('poiseulleFlow.h5','r')
        assert np.allclose(expected.root.velocity_t1,
                           actual.root.velocity_t1)
        expected.close()
        actual.close()
        


if __name__ == '__main__':
    from proteus import Comm
    comm = Comm.init()   
    test = TestStokesOperatorConstruction()
    test.setUp()
#    test.test_02_PressureLaplace()
    test.test_05_FullRun()

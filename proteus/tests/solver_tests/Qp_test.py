#!/usr/bin/env python
"""

Test module for the pressure mass matrix schur complement

"""
from proteus.iproteus import *
from proteus import Comm
from proteus import LinearAlgebraTools
import proteus.tests.TestTools

import os,sys,inspect
import tables
import numpy as np
from petsc4py import PETSc as p4pyPETSc
from scipy.sparse import csr_matrix

proteus.tests.TestTools.addSubFolders(inspect.currentframe())
import stokes_2d_p
import stokes_2d_n

def test_pressure_mass_matrix_shell():
    ''' This function tests pressure mass matrix shells '''
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

class TestStokesOperatorConstruction(proteus.tests.TestTools.SimulationTest):
    """ This class tests the operator construction on a 2D Stokes Poiseulle Flow problem """

    def setup_method(self,method):
        reload(stokes_2d_p)
        reload(stokes_2d_n)
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
        self._setRelativePath()

    def teardown_method(self,method):
        FileList = ['proteus.log',
                    'poiseulleFlow.h5',
                    'poiseulleFlow.xmf',                    
                    'rdomain.edge',
                    'rdomain.ele',
                    'rdomain.neig',
                    'rdomain.node',
                    'rdomain.poly']
        self.remove_files(FileList)

    def _setRelativePath(self):
        self.scriptdir = os.path.dirname(__file__)

    def test_01_MassMatrix(self):
        ''' Verify that Q returns the correct pressure and velocity mass matrix  '''
        Qp_petsc = self.operator_constructor.getQp()
        Qp_dense = LinearAlgebraTools.petsc4py_sparse_2_dense(Qp_petsc)
        Qv_petsc = self.operator_constructor.getQv()
        Qv_dense = LinearAlgebraTools.petsc4py_sparse_2_dense(Qv_petsc)

        rel_path_1 = "comparison_files/pressure_mass_matrix.txt"
        rel_path_2 = "comparison_files/velocity_mass_matrix.txt"

        pressure_mass_matrix = np.loadtxt(os.path.join(self.scriptdir,rel_path_1))
        velocity_mass_matrix = np.loadtxt(os.path.join(self.scriptdir,rel_path_2))

        assert np.allclose(pressure_mass_matrix,Qp_dense)
        assert np.allclose(velocity_mass_matrix,Qv_dense)

    def test_02_PressureLaplace(self):
        ''' Test the pressure Lapacian matrix '''
        Ap_raw = self.operator_constructor.getAp()
        Ap_dense = LinearAlgebraTools.petsc4py_sparse_2_dense(Ap_raw)
        rel_path = "comparison_files/pressure_laplace_matrix.txt"
        pressure_laplace_matrix = np.loadtxt(os.path.join(self.scriptdir,rel_path))
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
    
if __name__ == '__main__':
    pass

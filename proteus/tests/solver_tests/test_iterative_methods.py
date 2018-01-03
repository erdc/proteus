#!/usr/bin/env python
"""

Test module for solver shell operator.

"""
import proteus.test_utils.TestTools
from proteus.iproteus import *
from proteus import Comm
from proteus import LinearAlgebraTools as LAT
from proteus import LinearSolvers as LS

import os
import sys
import inspect
import petsc4py as p4pyPETSc
import numpy as np
import pytest
import pickle

from petsc4py import PETSc as p4pyPETSc

proteus.test_utils.TestTools.addSubFolders( inspect.currentframe() )

def create_petsc_vecs(matrix_A):
    """
    Creates a right-hand-side and solution PETSc4Py vector for
    testing ksp solves.

    Parameters
    ----------
    matrix_A: :class:`p4pyPETSc.Mat`
        Global matrix object

    Returns
    -------
    vec_lst: tuple
        This is a list of :class:`pypyPETSc.Vec` where the first is
        a vector of ones (usually to act as a RHS-vector) while the
        second vector is a vector of zeros (usually to act as a
        storage vector for the solution).
    """
    b = p4pyPETSc.Vec().create()
    x = p4pyPETSc.Vec().create()
    b.createWithArray(np.ones(matrix_A.getSizes()[0][0]))
    x.createWithArray(np.zeros(matrix_A.getSizes()[0][0]))
    return (b, x)

def initialize_schur_ksp_obj(matrix_A, schur_approx):
    """
    Creates a right-hand-side and solution PETSc4Py vector for
    testing ksp solves.

    Parameters
    ----------
    matrix_A: :class:`p4pyPETSc.Mat`
        Global matrix object.
    schur_approx: :class:`LS.SchurPrecon`

    Returns
    -------
    ksp_obj: :class:`p4pyPETSc.KSP`
    """
    ksp_obj = p4pyPETSc.KSP().create()
    ksp_obj.setOperators(matrix_A,matrix_A)
    pc = schur_approx.pc
    ksp_obj.setPC(pc)
    ksp_obj.setFromOptions()
    pc.setFromOptions()
    pc.setOperators(matrix_A,matrix_A)
    pc.setUp()
    schur_approx.setUp(ksp_obj)
    ksp_obj.setUp()
    ksp_obj.pc.setUp()
    return ksp_obj

@pytest.fixture()
def initialize_petsc_options(request):
    """Initializes schur complement petsc options. """
    petsc_options = p4pyPETSc.Options()
    petsc_options.setValue('ksp_type','gmres')
    petsc_options.setValue('ksp_gmres_restart',500)
    petsc_options.setValue('ksp_atol',1e-20)
    petsc_options.setValue('ksp_gmres_modifiedgramschmidt','')
    petsc_options.setValue('pc_fieldsplit_type','schur')
    petsc_options.setValue('pc_fieldsplit_schur_fact_type','upper')
    petsc_options.setValue('pc_fieldsplit_schur_precondition','user')
    petsc_options.setValue('fieldsplit_velocity_ksp_type','preonly')
    petsc_options.setValue('fieldsplit_velocity_pc_type', 'lu')
    petsc_options.setValue('fieldsplit_pressure_ksp_type','preonly')

@pytest.fixture()
def load_nse_cavity_matrix(request):
    """Loads a Navier-Stokes matrix drawn from the MPRANS module. """
    A = LAT.petsc_load_matrix(os.path.join
                              (os.path.dirname(__file__),
                               'import_modules/NSE_cavity_matrix'))
    yield A

@pytest.fixture()
def load_nse_step_matrix(request):
    """
    Loads a Navier-Stokes matrix for the backwards step problem from
    the MPRANS module.  This matrix is constructed using no-slip
    boundary conditions, and weakly enforced Dirichlet conditions.
    """
    A = LAT.petsc_load_matrix(os.path.join
                              (os.path.dirname(__file__),
                               'import_modules/NSE_step_no_slip'))
    yield A

@pytest.mark.LinearSolvers
@pytest.mark.skip(reason='this test is completed in a different PR')
def test_Schur_Sp_solve_global_null_space(load_nse_cavity_matrix,
                                          initialize_petsc_options):
    """Tests a KSP solve using the Sp Schur complement approximation.
    For this test, the global matrix has a null space because the
    boundary conditions are pure Dirichlet. """
    mat_A = load_nse_cavity_matrix
    b, x = create_petsc_vecs(mat_A)
    petsc_options = initialize_petsc_options

    solver_info = LS.ModelInfo(3,'interlaced')
    schur_approx = LS.Schur_Sp(mat_A,
                               '',
                               True,
                               solver_info=solver_info)
    ksp_obj = initialize_schur_ksp_obj(mat_A,schur_approx)
    ksp_obj.solve(b,x)

    assert ksp_obj.converged == True
    assert ksp_obj.its == 35
    assert np.allclose(ksp_obj.norm, 0.0007464632)
    assert ksp_obj.reason == 2

@pytest.mark.LinearSolvers
@pytest.mark.skip(reason='this test is completed in a different PR')
def test_Schur_Sp_solve(load_nse_step_matrix,
                        initialize_petsc_options):
    """Tests a KSP solve using the Sp Schur complement approximation.
       For this test, the global matrix does not have a null space."""
    mat_A = load_nse_step_matrix
    b, x = create_petsc_vecs(mat_A)

    solver_info = LS.ModelInfo(3, 'interlaced')
    schur_approx = LS.Schur_Sp(mat_A,
                               '',
                               solver_info=solver_info)
    ksp_obj = initialize_schur_ksp_obj(mat_A, schur_approx)
    ksp_obj.solve(b,x)

    assert ksp_obj.converged == True
    assert ksp_obj.its == 45
    assert np.allclose(ksp_obj.norm, 394.7036050627)
    assert ksp_obj.reason == 2

class TestIterativeMethods(proteus.test_utils.TestTools.BasicTest):

    def setup_method(self,method):
        self.petsc_options = p4pyPETSc.Options()
        self._scriptdir = os.path.dirname(__file__)
        self.quad_mass_matrix = np.load(os.path.join(self._scriptdir,
                                        'import_modules/quad_mass_matrix.npy'))

    def teardown_method(self,method):
        pass

    @pytest.mark.LinearAlgebraTools
    def test_dense_numpy_2_petsc4py(self):
        A_petsc = LAT.dense_numpy_2_petsc4py(self.quad_mass_matrix)
        A_new = LAT.petsc4py_sparse_2_dense(A_petsc)
        assert np.linalg.norm(self.quad_mass_matrix - A_new) == 0

    @pytest.mark.LinearSolvers
    def test_chebyshev_iteration_1(self):
        '''  Tests the pcd_shell operators produce correct output. '''
        A = self.quad_mass_matrix
        n = self.quad_mass_matrix.shape[0]
        alpha = 1./4
        beta = 9./4
        x0 = np.zeros(n)
        b1 = np.ones(n)
        for i in range(0,n,2):
            b1[i] = 0.
        A_petsc = LAT.dense_numpy_2_petsc4py(A)
        x0_petsc = p4pyPETSc.Vec().createWithArray(x0)
        b1_petsc = p4pyPETSc.Vec().createWithArray(b1)
        solver = LS.ChebyshevSemiIteration(A_petsc,
                                           alpha,
                                           beta,
                                           True)
        solver.apply(b1_petsc, x0_petsc, 20)
        expected = np.load(os.path.join(self._scriptdir,'import_modules/sol_10.npy'))
        actual = x0_petsc
        assert np.allclose(expected,actual.getArray())

    @pytest.mark.LinearSolvers
    def test_chebyshev_iteration_2(self):
        '''  Tests the pcd_shell operators produce correct output. '''
        A = np.diag(1./np.diag(self.quad_mass_matrix)).dot(self.quad_mass_matrix)
        n = self.quad_mass_matrix.shape[0]
        alpha = 1./4
        beta = 9./4
        x0 = np.zeros(n)
        b1 = np.zeros(n)
        for i in range(0,n):
            b1[i] = i
        A_petsc = LAT.dense_numpy_2_petsc4py(A)
        x0_petsc = p4pyPETSc.Vec().createWithArray(x0)
        b1_petsc = p4pyPETSc.Vec().createWithArray(b1)
        solver = LS.ChebyshevSemiIteration(A_petsc,
                                           alpha,
                                           beta,
                                           save_iterations=True)
        solver.apply(b1_petsc, x0_petsc, 20)
        expected = np.load(os.path.join(self._scriptdir,'import_modules/sol_20_lst.npy'))
        for i,item in enumerate(expected):
            assert np.allclose(item,solver.iteration_results[i],1e-12)

if __name__ == '__main__':
    pass

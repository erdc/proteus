#!/usr/bin/env python
"""

Test module for solver shell operator.

"""
from __future__ import division
from builtins import range
from past.utils import old_div
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

def initialize_asm_ksp_obj(matrix_A):
    """
    Creates a right-hand-side and solution PETSc4Py vector for
    testing ksp solves.

    Parameters
    ----------
    matrix_A: :class:`p4pyPETSc.Mat`
        Global matrix object.

    Returns
    -------
    ksp_obj: :class:`p4pyPETSc.KSP`
    """
    ksp_obj = p4pyPETSc.KSP().create()
    ksp_obj.setOperators(matrix_A,matrix_A)
    ksp_obj.setFromOptions()
    ksp_obj.setUp()
    return ksp_obj


def build_amg_index_sets(L_sizes):
    """
    Create PETSc index sets for the velocity components of a saddle
    point matrix

    Parameters
    ----------
    L_sizes : 
        Sizes of original saddle-point system

    Returns:
    --------
    Index_Sets : lst
        List of velocity index sets
    """
    neqns = L_sizes[0][0]
    velocityDOF=[]
    for start in range(1,3):
        velocityDOF.append(np.arange(start=start,
                                     stop=1+neqns,
                                     step=3,
                                     dtype='i'))
        velocityDOF_full=np.vstack(velocityDOF).transpose().flatten()
    velocity_u_DOF = []
    velocity_u_DOF.append(np.arange(start=1,
                                    stop=1+neqns,
                                    step=3,
                                    dtype='i'))
    velocity_u_DOF_full = np.vstack(velocity_u_DOF).transpose().flatten()
    velocity_v_DOF = []
    velocity_v_DOF.append(np.arange(start=2,
                                    stop=2+neqns,
                                    step=3,
                                    dtype='i'))
    velocity_v_DOF_full = np.vstack(velocity_v_DOF).transpose().flatten()
    isvelocity = p4pyPETSc.IS()
    isvelocity.createGeneral(velocityDOF_full)
    isu = p4pyPETSc.IS()
    isu.createGeneral(velocity_u_DOF_full)
    isv = p4pyPETSc.IS()
    isv.createGeneral(velocity_v_DOF_full)
    return [isvelocity, isu, isv]


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
def initialize_velocity_block_petsc_options(request):
    petsc_options = p4pyPETSc.Options()
    petsc_options.setValue('ksp_type','gmres')
    petsc_options.setValue('ksp_gmres_restart',100)
    petsc_options.setValue('ksp_pc_side','right')
    petsc_options.setValue('ksp_atol',1e-8)
    petsc_options.setValue('ksp_gmres_modifiedgramschmidt','')
    petsc_options.setValue('pc_type','hypre')
    petsc_options.setValue('pc_type_hypre_type','boomeramg')

def load_matrix(mat_file_name):
    """Load a matrix """
    A = LAT.petsc_load_matrix(os.path.join
                              (os.path.dirname(__file__),
                               'import_modules/'+mat_file_name))
    return A
    
@pytest.fixture()
def load_nse_cavity_matrix(request):
    """Loads a Navier-Stokes matrix drawn from the MPRANS module. """
    A = LAT.petsc_load_matrix(os.path.join
                              (os.path.dirname(__file__),
                               'import_modules/NSE_cavity_matrix.bin'))
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
                               'import_modules/NSE_step_no_slip.bin'))
    yield A

@pytest.fixture()
def load_small_step_matrix(request):
    """
    Loads a small example of a backwards facing step matrix for
    testing purposes. (Note: this matrix does not have advection)
    """
    A = LAT.petsc_load_matrix(os.path.join
                              (os.path.dirname(__file__),
                               'import_modules/saddle_point_small.bin'))
    yield A

@pytest.fixture()
def load_medium_step_matrix(request):
    """
    Loads a medium sized backwards facing step matrix for studying
    different AMG preconditioners. (Note: this matrix does not have
    advection)
    """
    A = LAT.petsc_load_matrix(os.path.join
                              (os.path.dirname(__file__),
                               'import_modules/saddle_point_matrix.bin'))
    yield A

@pytest.fixture()
def load_rans2p_step_newton_1(request):
    A = LAT.petsc_load_matrix(os.path.join
                              (os.path.dirname(__file__),
                               'import_modules/rans2p_step_newton_1.bin'))
    yield A

@pytest.fixture()
def load_rans2p_step_newton_5(request):
    A = LAT.petsc_load_matrix(os.path.join
                              (os.path.dirname(__file__),
                               'import_modules/rans2p_step_newton_5.bin'))
    yield A

@pytest.mark.LinearSolvers
def test_Schur_Sp_solve_global_null_space(load_nse_cavity_matrix,
                                          initialize_petsc_options):
    """Tests a KSP solve using the Sp Schur complement approximation.
    For this test, the global matrix has a null space because the
    boundary conditions are pure Dirichlet. """
    mat_A = load_nse_cavity_matrix
    b, x = create_petsc_vecs(mat_A)
    petsc_options = initialize_petsc_options

    solver_info = LS.ModelInfo('interlaced',
                               3,
                               bdy_null_space=True)
    schur_approx = LS.Schur_Sp(L=mat_A,
                               prefix='',
                               solver_info=solver_info)
    ksp_obj = initialize_schur_ksp_obj(mat_A,schur_approx)
    ksp_obj.solve(b,x)

    assert ksp_obj.converged == True
    assert ksp_obj.its == 35
    assert np.allclose(ksp_obj.norm, 0.0007464632)
    assert ksp_obj.reason == 2

@pytest.mark.LinearSolvers
def test_Schur_Sp_solve(load_nse_step_matrix,
                        initialize_petsc_options):
    """Tests a KSP solve using the Sp Schur complement approximation.
       For this test, the global matrix does not have a null space."""
    mat_A = load_nse_step_matrix
    b, x = create_petsc_vecs(mat_A)

    solver_info = LS.ModelInfo('interlaced', 3)
    schur_approx = LS.Schur_Sp(mat_A,
                               '',
                               solver_info=solver_info)
    ksp_obj = initialize_schur_ksp_obj(mat_A, schur_approx)
    ksp_obj.solve(b,x)

    assert ksp_obj.converged == True
    assert ksp_obj.its == 45
    assert np.allclose(ksp_obj.norm, 394.7036050627)
    assert ksp_obj.reason == 2
    
@pytest.mark.amg
def test_amg_basic(load_small_step_matrix,
                   initialize_velocity_block_petsc_options):
    mat_A = load_small_step_matrix

    petsc_options = initialize_velocity_block_petsc_options
    L_sizes = mat_A.getSizes()
    index_sets = build_amg_index_sets(L_sizes)

    #Initialize ksp object
    F_ksp = initialize_asm_ksp_obj(mat_A.createSubMatrix(index_sets[0],
                                                      index_sets[0]))
    b, x = create_petsc_vecs(mat_A.createSubMatrix(index_sets[0],
                                                index_sets[0]))   
    F_ksp.solve(b,x)
    assert F_ksp.its == 9

@pytest.mark.amg
def test_amg_iteration_performance(load_medium_step_matrix,
                                   initialize_velocity_block_petsc_options):
    mat_A = load_medium_step_matrix
    petsc_options = initialize_velocity_block_petsc_options
    L_sizes = mat_A.getSizes()
    index_sets = build_amg_index_sets(L_sizes)

    F_ksp = initialize_asm_ksp_obj(mat_A.createSubMatrix(index_sets[0],
                                                      index_sets[0]))
    b, x = create_petsc_vecs(mat_A.createSubMatrix(index_sets[0],
                                                index_sets[0]))

    F_ksp.solve(b,x)
    assert F_ksp.its == 41

def test_amg_step_problem_01(load_rans2p_step_newton_1,
                             initialize_velocity_block_petsc_options):
    mat_A = load_rans2p_step_newton_1
    petsc_options = initialize_velocity_block_petsc_options
    L_sizes = mat_A.getSizes()
    index_sets = build_amg_index_sets(L_sizes)

    F_ksp = initialize_asm_ksp_obj(mat_A.createSubMatrix(index_sets[0],
                                                      index_sets[0]))
    b, x = create_petsc_vecs(mat_A.createSubMatrix(index_sets[0],
                                                index_sets[0]))
    F_ksp.solve(b,x)
    assert F_ksp.its == 59

def test_amg_step_problem_02(load_rans2p_step_newton_5,
                             initialize_velocity_block_petsc_options):
    mat_A = load_rans2p_step_newton_5
    petsc_options = initialize_velocity_block_petsc_options
    L_sizes = mat_A.getSizes()
    index_sets = build_amg_index_sets(L_sizes)

    F_ksp = initialize_asm_ksp_obj(mat_A.createSubMatrix(index_sets[0],
                                                      index_sets[0]))
    b, x = create_petsc_vecs(mat_A.createSubMatrix(index_sets[0],
                                                index_sets[0]))
    F_ksp.solve(b,x)
    assert F_ksp.its == 60

class TestSmoothingAlgorithms(proteus.test_utils.TestTools.BasicTest):

    def setup_method(self,method):
        self._scriptdir = os.path.dirname(__file__)
        self.saddle_point_matrix=LAT.petsc_load_matrix(os.path.join(self._scriptdir,
                                                                    'import_modules/saddle_point_small'))
        # self.saddle_point_matrix = LAT.petsc_load_matrix(os.path.join(self._scriptdir,
        #                                                               'import_modules/saddle_point_matrix'))
    def test_matrix_splitting_1(self):
        vals_F  =    [3.2, 1.1, 5.4, 6.3, 1., -5.1, 1.2]
        col_idx_F  = [0, 1, 2, 0, 2, 0, 1]
        row_idx_F  = [0, 3, 5, 7]

        num_v_unkwn = len(row_idx_F) - 1

        petsc_matF = LAT.csr_2_petsc(size = (num_v_unkwn,num_v_unkwn),
                                     csr = (row_idx_F,col_idx_F,vals_F))

        A = LAT.split_PETSc_Mat(petsc_matF)
        A[0].axpy(1.0,A[1])
        assert np.allclose(A[0].getValuesCSR()[2], petsc_matF.getValuesCSR()[2])

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
        alpha = old_div(1.,4)
        beta = old_div(9.,4)
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
        A = np.diag(old_div(1.,np.diag(self.quad_mass_matrix))).dot(self.quad_mass_matrix)
        n = self.quad_mass_matrix.shape[0]
        alpha = old_div(1.,4)
        beta = old_div(9.,4)
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

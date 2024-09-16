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
import petsc4py
import numpy as np
import pytest
import pickle

from petsc4py import PETSc

proteus.test_utils.TestTools.addSubFolders( inspect.currentframe() )

def create_petsc_vecs(matrix_A):
    """
    Creates a right-hand-side and solution PETSc4Py vector for
    testing ksp solves.

    Parameters
    ----------
    matrix_A: :class:`PETSc.Mat`
        Global matrix object

    Returns
    -------
    vec_lst: tuple
        This is a list of :class:`pypyPETSc.Vec` where the first is
        a vector of ones (usually to act as a RHS-vector) while the
        second vector is a vector of zeros (usually to act as a
        storage vector for the solution).
    """
    b = PETSc.Vec().create()
    x = PETSc.Vec().create()
    b.createWithArray(np.ones(matrix_A.getSizes()[0][0],'d'))
    x.createWithArray(np.zeros(matrix_A.getSizes()[0][0],'d'))
    return (b, x)


def initialize_asm_ksp_obj(matrix_A):
    """
    Creates a right-hand-side and solution PETSc4Py vector for
    testing ksp solves.

    Parameters
    ----------
    matrix_A: :class:`PETSc.Mat`
        Global matrix object.

    Returns
    -------
    ksp_obj: :class:`PETSc.KSP`
    """
    ksp_obj = PETSc.KSP().create()
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
    isvelocity = PETSc.IS()
    isvelocity.createGeneral(velocityDOF_full)
    isu = PETSc.IS()
    isu.createGeneral(velocity_u_DOF_full)
    isv = PETSc.IS()
    isv.createGeneral(velocity_v_DOF_full)
    return [isvelocity, isu, isv]



def load_matrix(mat_file_name):
    """Load a matrix """
    A = LAT.petsc_load_matrix(os.path.join
                              (os.path.dirname(__file__),
                               'import_modules/'+mat_file_name))
    return A
    

class TestSmoothingAlgorithms(proteus.test_utils.TestTools.BasicTest):

    def setup_method(self,method):
        self._scriptdir = os.path.dirname(__file__)
        self.saddle_point_matrix=LAT.petsc_load_matrix(os.path.join(self._scriptdir,
                                                                    'import_modules/saddle_point_small'))
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
        self.petsc_options = PETSc.Options()
        self.petsc_options.clear()
        for k in self.petsc_options.getAll(): self.petsc_options.delValue(k)
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
        x0_petsc = PETSc.Vec().createWithArray(x0)
        b1_petsc = PETSc.Vec().createWithArray(b1)
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
        x0_petsc = PETSc.Vec().createWithArray(x0)
        b1_petsc = PETSc.Vec().createWithArray(b1)
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

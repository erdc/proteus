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

from petsc4py import PETSc as p4pyPETSc

proteus.test_utils.TestTools.addSubFolders( inspect.currentframe() )


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

#!/usr/bin/env python
""" Test modules for Driven Cavity Stokes preconditioners. """

import proteus.test_utils.TestTools as TestTools
import proteus.LinearAlgebraTools as LAT
import proteus.LinearSolvers as LS
from proteus.iproteus import *


import os
import sys
import inspect
import numpy as np
import tables
import pickle
import petsc4py
from petsc4py import PETSc
import pytest

TestTools.addSubFolders( inspect.currentframe() )
import stokesDrivenCavity_2d_p
import stokesDrivenCavity_2d_n       
def create_petsc_vecs(matrix_A):
    """
    Creates a right-hand-side and solution PETSc vector for
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
    b.createWithArray(np.ones(matrix_A.getSizes()[0][0]))
    x.createWithArray(np.zeros(matrix_A.getSizes()[0][0]))
    return (b, x)

@pytest.mark.LinearSolvers
@pytest.mark.modelTest
class TestStokes(proteus.test_utils.TestTools.SimulationTest):
    """Run a Stokes test with mumps LU factorization """

    def setup_method(self):
        reload(stokesDrivenCavity_2d_p)
        reload(stokesDrivenCavity_2d_n)
        self.pList = [stokesDrivenCavity_2d_p]
        self.nList = [stokesDrivenCavity_2d_n]
        self.so = default_so
        self.so.tnList = [0.,1.]
        self.so.name = self.pList[0].name
        self.so.sList = self.pList[0].name
        self.so.sList = [default_s]

    def teardown_method(self):
        """Tear down function. """
        Profiling.closeLog()
        FileList = ['proteus_default.log',
                    'proteus.log',
                    #'rdomain.ele',
                    #'rdomain.edge',
                    #'rdomain.neig',
                    #'rdomain.node',
                    #'rdomain.poly',
                    'drivenCavityStokesTrial.h5',
                    'drivenCavityStokesTrial.xmf']
        self.remove_files(FileList)

    def _setPETSc(self):
        petsc4py.PETSc.Options().setValue("ksp_type","fgmres")
        petsc4py.PETSc.Options().setValue("ksp_atol",1e-20)
        petsc4py.PETSc.Options().setValue("ksp_atol",1e-12)
        petsc4py.PETSc.Options().setValue("pc_fieldsplit_type","schur")
        petsc4py.PETSc.Options().setValue("pc_fieldsplit_schur_fact_type","upper")
        petsc4py.PETSc.Options().setValue("fieldsplit_velocity_ksp_type","preonly")
        petsc4py.PETSc.Options().setValue("fieldsplit_velocity_pc_type","lu")
        petsc4py.PETSc.Options().setValue("fieldsplit_pressure_ksp_type","preonly")

    def _setPETSc_LU(self):
        petsc4py.PETSc.Options().setValue("ksp_type","preonly")
        petsc4py.PETSc.Options().setValue("pc_type","lu")
        petsc4py.PETSc.Options().setValue("pc_factor_mat_solver_package","mumps")

    def _runTest(self):
        Profiling.openLog('proteus.log',11)
        Profiling.verbose = True
        self._scriptdir = os.path.dirname(__file__)
        self.ns = NumericalSolution.NS_base(self.so,
                                            self.pList,
                                            self.nList,
                                            self.so.sList,
                                            opts)
        self.ns.calculateSolution('stokes')
        actual = tables.open_file('drivenCavityStokesTrial.h5','r')
        expected_path = 'comparison_files/' + 'comparison_' + 'drivenCavityStokes' + '_velocity_t1.csv'
        #write comparison file
        #np.array(actual.root.velocity_t1).tofile(os.path.join(self._scriptdir, expected_path),sep=",")
        np.testing.assert_almost_equal(np.fromfile(os.path.join(self._scriptdir, expected_path),sep=","),np.array(actual.root.velocity_t1).flatten(),decimal=2)
        actual.close()

    @pytest.mark.slowTest
    def test_01_FullRun(self):
        stokesDrivenCavity_2d_n.linearSmoother = proteus.LinearSolvers.Schur_Qp
        self._setPETSc()
        self._runTest()
        relpath = 'comparison_files/Qp_expected.log'
        actual_log = TestTools.NumericResults.build_from_proteus_log('proteus.log')
        expected_log = TestTools.NumericResults.build_from_proteus_log(os.path.join(self._scriptdir,
                                                                                    relpath))
        plot_lst = [(1.0,0,0),(1.0,1,0),(1.0,2,0)]
        L1 = expected_log.get_ksp_resid_it_info(plot_lst)
        L2 = actual_log.get_ksp_resid_it_info(plot_lst)
        assert L1 == L2

    @pytest.mark.slowTest
    def test_02_FullRun(self):
        self._setPETSc_LU()
        self._runTest()

def initialize_schur_ksp_obj(matrix_A, schur_approx):
    """
    Creates a right-hand-side and solution PETSc4Py vector for
    testing ksp solves.

    Parameters
    ----------
    matrix_A: :class:`PETSc.Mat`
        Global matrix object.
    schur_approx: :class:`LS.SchurPrecon`

    Returns
    -------
    ksp_obj: :class:`PETSc.KSP`
    """
    ksp_obj = PETSc.KSP().create()
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
def load_nse_cavity_matrix(request):
    """Loads a Navier-Stokes matrix drawn from the MPRANS module. """
    A = LAT.petsc_load_matrix('dump_stokes_drivenCavityStokesTrial_1.0par_j_1')
    yield A

@pytest.fixture()
def initialize_petsc_options(request):
    """Initializes schur complement petsc options. """
    petsc_options = PETSc.Options()
    petsc_options.setValue('ksp_type','gmres')
    petsc_options.setValue('ksp_gmres_restart',500)
    petsc_options.setValue('ksp_atol',1e-16)
    petsc_options.setValue('ksp_rtol',1.0e-12)
    petsc_options.setValue('ksp_gmres_modifiedgramschmidt','')
    petsc_options.setValue('pc_fieldsplit_type','schur')
    petsc_options.setValue('pc_fieldsplit_schur_fact_type','upper')
    petsc_options.setValue('pc_fieldsplit_schur_precondition','user')
    petsc_options.setValue('fieldsplit_velocity_ksp_type','preonly')
    petsc_options.setValue('fieldsplit_velocity_pc_type', 'lu')
    petsc_options.setValue('fieldsplit_pressure_ksp_type','preonly')

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
    assert ksp_obj.its == 89
    assert ksp_obj.norm < np.linalg.norm(b)*1.0e-9 + 1.0e-16
    assert ksp_obj.reason == 2

if __name__ == '__main__':
    pass

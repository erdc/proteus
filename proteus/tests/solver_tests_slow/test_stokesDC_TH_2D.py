#!/usr/bin/env python
""" Test modules for Driven Cavity Stokes preconditioners. """

import proteus.test_utils.TestTools as TestTools
from proteus.iproteus import *

Profiling.logLevel = 7
Profiling.verbose = True

import os
import sys
import inspect
import numpy
import tables
import pickle
import petsc4py
import pytest

TestTools.addSubFolders( inspect.currentframe() )
import stokesDrivenCavity_2d_p
import stokesDrivenCavity_2d_n        

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
        FileList = ['proteus_default.log',
                    'proteus.log',
                    'rdomain.ele',
                    'rdomain.edge',
                    'rdomain.neig',
                    'rdomain.node',
                    'rdomain.poly',
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
        Profiling.openLog('proteus.log',7)
        self._scriptdir = os.path.dirname(__file__)
        self.ns = NumericalSolution.NS_base(self.so,
                                            self.pList,
                                            self.nList,
                                            self.so.sList,
                                            opts)
        self.ns.calculateSolution('stokes')
        relpath = 'comparison_files/drivenCavityStokes_expected.h5'
        expected = tables.open_file(os.path.join(self._scriptdir,relpath))
        actual = tables.open_file('drivenCavityStokesTrial.h5','r')
        assert numpy.allclose(expected.root.velocity_t1,
                              actual.root.velocity_t1,
                              atol=1e-2)
        expected.close()
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

if __name__ == '__main__':
    pass

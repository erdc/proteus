#!/usr/bin/env python
""" Test modules for Driven Cavity Stokes preconditioners. """

import proteus.tests.TestTools
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

proteus.tests.TestTools.addSubFolders( inspect.currentframe() )
import stokesDrivenCavity_2d_p
import stokesDrivenCavity_2d_n        

class TestStokes(proteus.tests.TestTools.SimulationTest):
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
        petsc4py.PETSc.Options().setValue("ksp_atol",1e-6)
        petsc4py.PETSc.Options().setValue("pc_fieldsplit_type","schur")
        petsc4py.PETSc.Options().setValue("pc_fieldsplit_schur_fact_type","full")
        petsc4py.PETSc.Options().setValue("fieldsplit_velocity_ksp_type","preonly")
        petsc4py.PETSc.Options().setValue("fieldsplit_velocity_pc_type","lu")
        petsc4py.PETSc.Options().setValue("fieldsplit_pressure_ksp_type","fgmres")
        petsc4py.PETSc.Options().setValue("fieldsplit_pressure_ksp_atol",1e-8)
        petsc4py.PETSc.Options().setValue("fieldsplit_pressure_ksp_rtol",1e-8)

    def _setPETSc_LU(self):
        petsc4py.PETSc.Options().setValue("ksp_type","preonly")
        petsc4py.PETSc.Options().setValue("pc_type","lu")
        petsc4py.PETSc.Options().setValue("pc_factor_mat_solver_package","mumps")

    def _runTest(self):
        self._scriptdir = os.path.dirname(__file__)
        self.ns = NumericalSolution.NS_base(self.so,
                                            self.pList,
                                            self.nList,
                                            self.so.sList,
                                            opts)
        self.ns.calculateSolution('stokes')
        relpath = 'comparison_files/drivenCavityStokes_expected.h5'
        expected = tables.openFile(os.path.join(self._scriptdir,relpath))
        actual = tables.openFile('./drivenCavityStokesTrial.h5','r')
        assert numpy.allclose(expected.root.velocity_t1,
                              actual.root.velocity_t1,
                              atol=1e-2)
        expected.close()
        actual.close()        
        
    @pytest.mark.skip(reason="in development")
    def test_01_FullRun(self):
        self._setPETSc()
        self._runTest()

    def test_02_FullRun(self):
        self._setPETSc_LU()
        self._runTest()

if __name__ == '__main__':
    pass

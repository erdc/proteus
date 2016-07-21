#!/usr/bin/env python
""" Test module for the pressure mass matrix preconditioner for a Stokes problem. """

import os
import sys
import inspect
import numpy
import tables
import pickle
import petsc4py

import proteus.tests.TestTools
from proteus.iproteus import *


Profiling.logLevel = 7
Profiling.verbose = True

class TestStokesQp(proteus.tests.TestTools.SimulationTest):
    """Run a stokes test with the pressure mass matrix. """
    # def __init__(self):
    #     proteus.tests.TestTools.addSubFolders( inspect.currentframe() )
    #     import stokesDrivenCavity_2d_p
    #     import stokesDrivenCavity_2d_n
    #     pList = [stokesDrivenCavity_2d_p]
    #     nList = [stokesDrivenCavity_2d_n]
    #     so = default_so
    #     so.tnList = [0.,1.]
    #     so.name = pList[0].name
    #     so.sList = [default_s]
    #     opts.verbose = True
    #     opts.profile = True
    #     self._setPETSc()
    #     self._scriptdir = os.path.dirname(__file__)
    #     self.ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)        

    @classmethod
    def setup_class(self):
        proteus.tests.TestTools.addSubFolders( inspect.currentframe() )
        import stokesDrivenCavity_2d_p
        import stokesDrivenCavity_2d_n
        pList = [stokesDrivenCavity_2d_p]
        nList = [stokesDrivenCavity_2d_n]
        so = default_so
        so.tnList = [0.,1.]
        so.name = pList[0].name
        so.sList = pList[0].name
        so.sList = [default_s]
        self._setPETSc()
        self._scriptdir = os.path.dirname(__file__)
        self.ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)

    @classmethod
    def teardown_class(self):
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
        for file in FileList:
            if os.path.isfile(file):
                os.remove(file)

    @classmethod
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

    def test_01_FullRun(self):
        self.ns.calculateSolution('stokes')
        relpath = "comparison_files/drivenCavityStokes_expected.h5"
        expected = tables.openFile(os.path.join(self._scriptdir,relpath))
        actual = tables.openFile('./drivenCavityStokesTrial.h5','r')
        assert numpy.allclose(expected.root.velocity_t1,
                              actual.root.velocity_t1)        
        NA = proteus.tests.TestTools.NumericResults('proteus.log')
        relpath = "comparison_files/expected_NA"
        expected_NA_file = open(os.path.join(self._scriptdir,relpath),'r')
        expected_NA = pickle.load(expected_NA_file)
        assert (expected_NA == NA.data_dictionary)

if __name__ == '__main__':
    A = TestStokesQp()
    A.test_01_FullRun()
#    A.teardown_class()

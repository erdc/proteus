#!/usr/bin/env python
""" Test modules for Driven Cavity Stokes preconditioners. """

import proteus.test_utils.TestTools as TestTools
from proteus.iproteus import *

Profiling.logLevel = 7
Profiling.verbose = False

import os
import sys
import inspect
import numpy as np
import h5py
import pickle
import petsc4py
import pytest

proteus.test_utils.TestTools.addSubFolders( inspect.currentframe() )
from proteus import defaults
from NavierStokes_ST_LS_SO_VV import NavierStokes_ST_LS_SO_VV
import_modules = os.path.join(os.path.dirname(os.path.realpath(__file__)),'import_modules')
@pytest.mark.LinearSolvers
@pytest.mark.modelTest
@pytest.mark.navierstokesTest
class Test_NSE_Driven_Cavity(proteus.test_utils.TestTools.SimulationTest):

    def setup_method(self):
        nseDrivenCavity_2d_p = defaults.load_physics('nseDrivenCavity_2d_p',import_modules)
        nseDrivenCavity_2d_n = defaults.load_numerics('nseDrivenCavity_2d_n',import_modules)
        self.pList = [nseDrivenCavity_2d_p]
        self.nList = [nseDrivenCavity_2d_n]
        self.so = default_so
        self.so.name = self.pList[0].name
        self.so.sList = self.pList[0].name
        self.so.sList = [default_s]
        self._scriptdir = os.path.dirname(__file__)
        Profiling.openLog("proteus.log",10)
        Profiling.logAllProcesses = True

    def teardown_method(self):
        Profiling.closeLog()

    def _runTest(self):
        self.ns = NumericalSolution.NS_base(self.so,
                                            self.pList,
                                            self.nList,
                                            self.so.sList,
                                            opts)
        self.ns.calculateSolution('stokes')

        # The produced output has diverged from the old comparison
        # output. It needs to be confirmed that the new ouput is
        # in fact correct and drivenCavityNSE_LSC_expected.h5 should
        # be updated accordingly.
        actual = h5py.File('drivenCavityNSETrial.h5','r')
        relpath = 'comparison_files/drivenCavityNSE_LSC_expected.csv'
        #np.savetxt(os.path.join(self._scriptdir,relpath),actual.root.velocity_t7.read(),delimiter=',')
        expected = np.loadtxt(os.path.join(self._scriptdir,relpath),delimiter=',')
        assert np.allclose(expected,
                           actual['velocity_t7'][:],
                           rtol=1e-8, atol=1e-8) 
        actual.close()

    @pytest.mark.slowTest
    def test_01_FullRun(self):
        """Runs with nsedriven cavity with the following settings.
        - RE enforced through viscosity
        - THQuads
        - LSC
        """
        self.pList[0].coefficients = NavierStokes_ST_LS_SO_VV(epsFact=0.0,
                                                              sigma=0.0,
                                                              rho_0=1.0,
                                                              nu_0=1.0,
                                                              rho_1=1.0,
                                                              nu_1=1.0,
                                                              g=[0.0,0.0],
                                                              nd=2,
                                                              LS_model=None,
                                                              KN_model=None,
                                                              epsFact_density=None,
                                                              stokes=False);

        self.pList[0].coefficients.variableNames = ['p','u','v']

        self.pList[0].dirichletConditions = {0:self.pList[0].getDBCp,
                                             1:self.pList[0].getDBCu,
                                             2:self.pList[0].getDBCv}
        self.nList[0].nLevels = 1
        self.nList[0].linearSmoother = self.nList[0].Schur_LSC
        self.so.tnList = self.nList[0].tnList
        self._setPETSc(petsc_file = os.path.join(self._scriptdir,'import_modules/petsc.options.schur_lsc'))
        self._runTest()

if __name__ == '__main__':
    pass
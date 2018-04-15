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

proteus.test_utils.TestTools.addSubFolders( inspect.currentframe() )
from import_modules import step2d_so
from import_modules import step2d
from import_modules import twp_navier_stokes_step2d_p
from import_modules import twp_navier_stokes_step2d_n

@pytest.mark.LinearSolvers
@pytest.mark.modelTest
@pytest.mark.navierstokesTest
class Test_NSE_Driven_Cavity(proteus.test_utils.TestTools.SimulationTest):

    def setup_method(self):
        reload(step2d_so)
        reload(step2d)
        reload(twp_navier_stokes_step2d_p)
        reload(twp_navier_stokes_step2d_n)
        self.pList = [twp_navier_stokes_step2d_p]
        self.nList = [twp_navier_stokes_step2d_n]
        self.pList[0].name = 'test_1'        
        self.so = step2d_so
        self.so.name = self.pList[0].name
        self.so.sList = self.pList[0].name
        self.so.sList = [default_s]
        self._scriptdir = os.path.dirname(__file__)
        Profiling.openLog("proteus.log",10)
        Profiling.logAllProcesses = True

    def teardown_method(self):
        Profiling.closeLog()

    def _runTest(self):
        assert 1 == 1
        self.ns = NumericalSolution.NS_base(self.so,
                                            self.pList,
                                            self.nList,
                                            self.so.sList,
                                            opts)
        self.ns.calculateSolution('stokes')

        # relpath = 'comparison_files/drivenCavityNSE_LSC_expected.h5'
        # expected = tables.open_file(os.path.join(self._scriptdir,relpath))
        # actual = tables.open_file('drivenCavityNSETrial.h5','r')
        # assert numpy.allclose(expected.root.velocity_t7,
        #                       actual.root.velocity_t7,
        #                       atol=1e-2)
        # expected.close()
        # actual.close()
        # relpath = 'comparison_files/drivenCavityNSE_LSC_expected.log'
        actual_log = TestTools.NumericResults.build_from_proteus_log('proteus.log')
        # expected_log = TestTools.NumericResults.build_from_proteus_log(os.path.join(self._scriptdir,
        #                                                                             relpath))
        L1 = actual_log.get_ksp_resid_it_info([(' test_1 ',1e+18,0,0)])
        L2 = actual_log.get_ksp_resid_it_info([(' test_1 ',1e+18,0,1)])
        L3 = actual_log.get_ksp_resid_it_info([(' test_1 ',1e+18,0,2)])
        L4 = actual_log.get_ksp_resid_it_info([(' test_1 ',1e+18,0,3)])
        L5 = actual_log.get_ksp_resid_it_info([(' test_1 ',1e+18,0,4)])
        L6 = actual_log.get_ksp_resid_it_info([(' test_1 ',1e+18,0,5)])        

        assert L1[0][1]==25
        assert L2[0][1]==28
        assert L3[0][1]==41
        assert L4[0][1]==37
        assert L5[0][1]==38
        assert L6[0][1]==38
        # plot_lst = [(3.7,0,3),(3.2,0,2),(2.7,0,2),(2.2,0,1),(1.7,0,1)]
        # L1 = expected_log.get_ksp_resid_it_info(plot_lst)
        # L2 = actual_log.get_ksp_resid_it_info(plot_lst)
        # assert L1 == L2

    @pytest.mark.slowTest
    def test_01_FullRun(self):
        """Runs with nsedriven cavity with the following settings.
        """
        self._setPETSc(petsc_file = os.path.join(self._scriptdir,'import_modules/petsc.options.schur'))
        self._runTest()

if __name__ == '__main__':
    pass

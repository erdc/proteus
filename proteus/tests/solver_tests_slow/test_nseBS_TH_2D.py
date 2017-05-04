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
import nseBackwardsFacingStep_2d_n
import nseBackwardsFacingStep_2d_p
from NavierStokes_ST_LS_SO_VV import NavierStokes_ST_LS_SO_VV

class Test_NSE_BackwardsFacingStep(proteus.test_utils.TestTools.SimulationTest):

    def setup_method(self):
        reload(nseBackwardsFacingStep_2d_p)
        reload(nseBackwardsFacingStep_2d_n)
        self.pList = [nseBackwardsFacingStep_2d_p]
        self.nList = [nseBackwardsFacingStep_2d_n]
        self.so = default_so
        self.so.name = self.pList[0].name
        self.so.sList = self.pList[0].name
        self.so.sList = [default_s]
        self._scriptdir = os.path.dirname(__file__)
        Profiling.openLog("proteus.log",10)
        Profiling.logAllProcesses = True

    def _runTest(self):
        self.ns = NumericalSolution.NS_base(self.so,
                                            self.pList,
                                            self.nList,
                                            self.so.sList,
                                            opts)
        self.ns.calculateSolution('stokes')
        
        relpath = 'comparison_files/TrialTriPCD_2D.h5'
        expected = tables.openFile(os.path.join(self._scriptdir,relpath))
        actual = tables.openFile('backwardsFacingStepTrial.h5','r')
        assert numpy.allclose(expected.root.velocity_t5,
                              actual.root.velocity_t5,
                              atol=1e-2)
        expected.close()
        actual.close()
        
        relpath = 'comparison_files/nseBackwardsFacingStep_2D_TriPCD.log'
        actual_log = TestTools.NumericResults.build_from_proteus_log('proteus.log')
        expected_log = TestTools.NumericResults.build_from_proteus_log(os.path.join(self._scriptdir,
                                                                                    relpath))
        plot_lst = [(1.0,0,1),(1.7,0,1),(2.2,0,1),(2.7,0,2),(3.0,0,2)]
        L1 = expected_log.get_ksp_resid_it_info(plot_lst)
        L2 = actual_log.get_ksp_resid_it_info(plot_lst)
        assert L1 == L2        
        
    def test_01_FullRun(self):
        """Runs with nsedriven cavity with the following settings.
        - RE enforced through viscosity
        - THQuads
        - LSC
        """
        nseBackwardsFacingStep_2d_p.dirichletConditions = {0:nseBackwardsFacingStep_2d_p.getDBCp,
                                                           1:nseBackwardsFacingStep_2d_p.getDBCu,
                                                           2:nseBackwardsFacingStep_2d_p.getDBCv}
        
        nseBackwardsFacingStep_2d_p.Advectivefluxboundaryconditions =  {0:nseBackwardsFacingStep_2d_p.getAdvFluxBCp,
                                                                        1:nseBackwardsFacingStep_2d_p.getAdvFluxBCu,
                                                                        2:nseBackwardsFacingStep_2d_p.getAdvFluxBCv}
        
        nseBackwardsFacingStep_2d_n.nLevels = 1
        nseBackwardsFacingStep_2d_n.linearSmoother = nseBackwardsFacingStep_2d_n.NavierStokes3D_PCD
        self.so.tnList = self.nList[0].tnList
        self._setPETSc(petsc_file = os.path.join(self._scriptdir,'import_modules/petsc.options.schur_pcd'))
        self._runTest()

if __name__ == '__main__':
    pass

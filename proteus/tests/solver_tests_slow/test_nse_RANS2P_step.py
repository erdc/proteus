#!/usr/bin/env python
""" Test modules for Driven Cavity Stokes preconditioners. """
from __future__ import absolute_import

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
try:
    from .import_modules import step2d_so
    from .import_modules import step2d
    from .import_modules import twp_navier_stokes_step2d_p
    from .import_modules import twp_navier_stokes_step2d_n
except:
    from import_modules import step2d_so
    from import_modules import step2d
    from import_modules import twp_navier_stokes_step2d_p
    from import_modules import twp_navier_stokes_step2d_n

def load_simulation(context_options_str=None):
    """
    Loads a two-phase step problem with settings

    Parameters
    ----------
    settings:

    Returns:
    --------

    """
    from proteus import Context
    Context.contextOptionsString=context_options_str

    reload(step2d_so)
    reload(step2d)
    reload(twp_navier_stokes_step2d_p)
    reload(twp_navier_stokes_step2d_n)

    pList = [twp_navier_stokes_step2d_p]
    nList = [twp_navier_stokes_step2d_n]
    pList[0].name = 'test_1'        
    so = step2d_so
    so.name = pList[0].name
    so.sList = pList[0].name
    so.sList = [default_s]
    _scriptdir = os.path.dirname(__file__)
    Profiling.openLog("proteus.log",10)
    Profiling.logAllProcesses = True
    ns = NumericalSolution.NS_base(so,
                                   pList,
                                   nList,
                                   so.sList,
                                   opts)
    return ns

def runTest(ns, name):
    ns.calculateSolution(name)
    actual_log = TestTools.NumericResults.build_from_proteus_log('proteus.log')
    return actual_log

@pytest.mark.skip
def test_01_FullRun():
    """ Runs two-dimensional step problem with the settings:
        * Strongly enforced Free-slip BC.
        * Pressure Projection Stablization.
        * he = 0.05
    """
    TestTools.SimulationTest._setPETSc(petsc_file = os.path.join(os.path.dirname(__file__),
                                                                 'import_modules/petsc.options.schur'))
    context_options_str='he=0.05'
    ns = load_simulation(context_options_str)
    actual_log = runTest(ns,'test_1')

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

@pytest.mark.skip
def test_02_FullRun():
    """ Runs two-dimensional step problem with the settings:
        * Strongly enforced no-slip BC.
        * Pressure Projection Stablization.
        * he = 0.05
    """
    TestTools.SimulationTest._setPETSc(petsc_file = os.path.join(os.path.dirname(__file__),
                                                                 'import_modules/petsc.options.schur'))
    context_options_str="boundary_condition_type='ns'"
    ns = load_simulation(context_options_str)
    actual_log = runTest(ns,'test_1')

    L1 = actual_log.get_ksp_resid_it_info([(' test_1 ',1e+18,0,0)])
    L2 = actual_log.get_ksp_resid_it_info([(' test_1 ',1e+18,0,1)])
    L3 = actual_log.get_ksp_resid_it_info([(' test_1 ',1e+18,0,2)])
    L4 = actual_log.get_ksp_resid_it_info([(' test_1 ',1e+18,0,3)])
    L5 = actual_log.get_ksp_resid_it_info([(' test_1 ',1e+18,0,4)])
    L6 = actual_log.get_ksp_resid_it_info([(' test_1 ',1e+18,0,5)])

    assert L1[0][1]==33
    assert L2[0][1]==35
    assert L3[0][1]==46
    assert L4[0][1]==44
    assert L5[0][1]==44
    assert L6[0][1]==44

if __name__ == '__main__':
    pass

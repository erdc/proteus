#!/usr/bin/env python
""" 
Test module for periodic boundary conditions and null space class.
"""
import proteus.test_utils.TestTools as TestTools
from proteus import defaults, default_s, Profiling, Context, Comm, NumericalSolution
from proteus.iproteus import opts

Profiling.logLevel = 7
Profiling.verbose = True

import os
import sys
import inspect
import numpy as np
import h5py
import pickle
import petsc4py
import pytest

from . import duct

@pytest.fixture()
def load_periodic_duct(request):
    nList = []
    pList = []
    sList = []
    script_dir = os.path.dirname(__file__)
    duct = defaults.load_source("duct",os.path.join(script_dir,'duct.py'))
    so = defaults.load_system('duct',path = script_dir)
    print(duct.opts)
    from petsc4py import PETSc
    OptDB=PETSc.Options()
    print(OptDB)
    OptDB.clear()
    for (pModule,nModule) in so.pnList:
        if not isinstance(pModule, defaults.Physics_base):
            pList.append(defaults.load_physics(pModule))
            if pList[-1].name == None:
                pList[-1].name = pModule
            nList.append(defaults.load_numerics(nModule))
        else:
            pList.append(pModule)
            nList.append(nModule)
    if so.sList == []:
        for i in range(len(so.pnList)):
            s = default_s
            sList.append(s)
    else:
        sList = so.sList
    yield so, pList, nList, sList

@pytest.fixture()
def load_periodic_opts_2D(request):
    Context.contextOptionsString="periodic=True grid=True nd=2 nnx=42 triangles=False spaceOrder=1 weak=True coord=True pc_type='selfp_petsc'" 
    script_dir = os.path.dirname(__file__)
    relpath = 'petsc/petsc.options.schur.selfp_petsc.superlu'
    petscOptionsFile = os.path.join(script_dir,relpath)
    Comm.argv = TestTools.fixture_set_petsc_options_from_file(petscOptionsFile)
    comm = Comm.init()

@pytest.fixture()
def load_periodic_opts_3D(request):
    Context.contextOptionsString="periodic=True nd=3 coord=True pc_type='selfp_petsc'"
    script_dir = os.path.dirname(__file__)
    relpath = 'petsc/petsc.options.schur.selfp_petsc.superlu'
    petscOptionsFile = os.path.join(script_dir,relpath)
    Comm.argv = TestTools.fixture_set_petsc_options_from_file(petscOptionsFile)
    comm = Comm.init()

@pytest.fixture()
def load_periodic_opts_3D_T2(request):
    Context.contextOptionsString="periodic=True nd=3 coord=True pc_type='selfp_petsc' A_block_AMG=True" 
    script_dir = os.path.dirname(__file__)
    relpath = 'petsc/petsc.options.schur.selfp_petsc.gamg.superlu'
    petscOptionsFile = os.path.join(script_dir,relpath)
    Comm.argv = TestTools.fixture_set_petsc_options_from_file(petscOptionsFile)
    comm = Comm.init()

def test_periodic_2D(load_periodic_opts_2D,
                     load_periodic_duct):
    ns = NumericalSolution.NS_base(load_periodic_duct[0],
                                   load_periodic_duct[1],
                                   load_periodic_duct[2],
                                   load_periodic_duct[3],
                                   opts)
    ns.calculateSolution('test_run')
    script_dir = os.path.dirname(__file__)
    actual = h5py.File('ductq1t12dpghe0.0975609756097561.h5')

    expected_path = 'comparison_files/' + 'comparison_' + 'ductq1t12dpghe0.0975609756097561' + '_velocity_t25.csv'
    #write comparison file
    #np.array(actual['velocity_t25']).tofile(os.path.join(script_dir, expected_path),sep=",")
    np.testing.assert_almost_equal(np.fromfile(os.path.join(script_dir, expected_path),sep=","),np.array(actual['velocity_t25']).flatten(),decimal=7)
    actual.close()



def test_periodic_3D(load_periodic_opts_3D,
                     load_periodic_duct):
    ns = NumericalSolution.NS_base(load_periodic_duct[0],
                                   load_periodic_duct[1],
                                   load_periodic_duct[2],
                                   load_periodic_duct[3],
                                   opts)
    ns.calculateSolution('test_run')

    script_dir = os.path.dirname(__file__)
    actual = h5py.File('ductp1t13dpghe0.4.h5')

    expected_path = 'comparison_files/' + 'comparison_' + 'ductp1t13dpghe0.4' + '_velocity_t25.csv'
    #write comparison file
    #np.array(actual['velocity_t25']).tofile(os.path.join(script_dir, expected_path),sep=",")
    np.testing.assert_almost_equal(np.fromfile(os.path.join(script_dir, expected_path),sep=","),np.array(actual['velocity_t25']).flatten(),decimal=10)
    actual.close()

def test_periodic_3D_amg(load_periodic_opts_3D_T2,
                         load_periodic_duct):
    ns = NumericalSolution.NS_base(load_periodic_duct[0],
                                   load_periodic_duct[1],
                                   load_periodic_duct[2],
                                   load_periodic_duct[3],
                                   opts)
    ns.calculateSolution('test_run')

    script_dir = os.path.dirname(__file__)
    actual = h5py.File('ductp1t13dpghe0.4.h5')
    expected_path = 'comparison_files/' + 'comparison_' + 'ductp1t13dpghe0.4' + '_velocity_t25.csv'
    np.testing.assert_almost_equal(np.fromfile(os.path.join(script_dir, expected_path),sep=","),np.array(actual['velocity_t25']).flatten(),decimal=6)
    actual.close()

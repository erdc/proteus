#!/usr/bin/env python
""" 
Test module for periodic boundary conditions and null space class.
"""
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

import duct

@pytest.fixture()
def load_periodic_duct(request):
    nList = []
    pList = []
    sList = []
    reload(duct)
    script_dir = os.path.dirname(__file__)
    so = proteus.defaults.load_system('duct',
                                      path = script_dir)
    for (pModule,nModule) in so.pnList:
        if not isinstance(pModule, proteus.defaults.Physics_base):
            pList.append(proteus.defaults.load_physics(pModule))
            if pList[-1].name == None:
                pList[-1].name = pModule
            nList.append(proteus.defaults.load_numerics(nModule))
        else:
            pList.append(pModule)
            nList.append(nModule)
    #
    if so.sList == []:
        for i in range(len(so.pnList)):
            s = default_s
            sList.append(s)
    else:
        sList = so.sList
    yield so, pList, nList, sList

@pytest.fixture()
def load_periodic_opts_2D(request):
    opts.contextOptions = "periodic=True grid=True nd=2 nnx=42 triangles=False spaceOrder=1 weak=True coord=True pc_type='selfp_petsc'" 
    proteus.Context.contextOptionsString=opts.contextOptions
    script_dir = os.path.dirname(__file__)
    relpath = 'petsc/petsc.options.schur.selfp_petsc.superlu'
    opts.petscOptionsFile = os.path.join(script_dir,relpath)
    proteus.Comm.argv = TestTools.fixture_set_petsc_options_from_file(opts.petscOptionsFile)
    comm = Comm.init()

@pytest.fixture()
def load_periodic_opts_3D(request):
    opts.contextOptions = "periodic=True nd=3 coord=True pc_type='selfp_petsc'"
    proteus.Context.contextOptionsString=opts.contextOptions
    script_dir = os.path.dirname(__file__)
    relpath = 'petsc/petsc.options.schur.selfp_petsc.superlu'
    opts.petscOptionsFile = os.path.join(script_dir,relpath)
    proteus.Comm.argv = TestTools.fixture_set_petsc_options_from_file(opts.petscOptionsFile)
    comm = Comm.init()

@pytest.fixture()
def load_periodic_opts_3D_T2(request):
    opts.contextOptions = "periodic=True nd=3 coord=True pc_type='selfp_petsc' A_block_AMG=True" 
    proteus.Context.contextOptionsString=opts.contextOptions
    script_dir = os.path.dirname(__file__)
    relpath = 'petsc/petsc.options.schur.selfp_petsc.amg'
    opts.petscOptionsFile = os.path.join(script_dir,relpath)
    proteus.Comm.argv = TestTools.fixture_set_petsc_options_from_file(opts.petscOptionsFile)
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
    relpath = 'comparison_files/basic_2d_test.h5'
    expected = tables.open_file(os.path.join(script_dir,relpath))
    actual = tables.open_file('ductq1t12dpghe0.0975609756098.h5')

    assert numpy.allclose(expected.root.velocity_t25.read(),
                          actual.root.velocity_t25.read(),
                          atol=1e-8)
    expected.close()
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
    relpath = 'comparison_files/basic_3d_test.h5'
    expected = tables.open_file(os.path.join(script_dir,relpath))
    actual = tables.open_file('ductp1t13dpghe0.2.h5')

    assert numpy.allclose(expected.root.velocity_t25.read(),
                          actual.root.velocity_t25.read(),
                          atol=1e-10)
    expected.close()
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
    relpath = 'comparison_files/basic_3d_test.h5'
    expected = tables.open_file(os.path.join(script_dir,relpath))
    actual = tables.open_file('ductp1t13dpghe0.2.h5')

    assert numpy.allclose(expected.root.velocity_t25.read(),
                          actual.root.velocity_t25.read(),
                          atol=1e-6)
    expected.close()
    actual.close()

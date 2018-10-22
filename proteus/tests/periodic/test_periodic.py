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
    so = proteus.defaults.load_system('duct')
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
    opts.petscOptionsFile = './petsc/petsc.options.schur.selfp_petsc.superlu'
    petsc_argv = [sys.argv[0]]
    with open(opts.petscOptionsFile) as petsc_file:
        data = petsc_file.readlines()
    def strip_comments(line):
        if '#' in line:
            line = line[:line.index('#')]
        return line
    stripped_data = [strip_comments(line) for line in data]
    petsc_argv += '\n'.join(stripped_data).split()
    proteus.Comm.argv = petsc_argv
    comm = Comm.init()

@pytest.fixture()
def load_periodic_opts_3D(request):
    opts.contextOptions = "periodic=True nd=3 coord=True pc_type='selfp_petsc'" 
    proteus.Context.contextOptionsString=opts.contextOptions
    opts.petscOptionsFile = './petsc/petsc.options.schur.selfp_petsc.superlu'
    petsc_argv = [sys.argv[0]]
    with open(opts.petscOptionsFile) as petsc_file:
        data = petsc_file.readlines()
    def strip_comments(line):
        if '#' in line:
            line = line[:line.index('#')]
        return line
    stripped_data = [strip_comments(line) for line in data]
    petsc_argv += '\n'.join(stripped_data).split()
    proteus.Comm.argv = petsc_argv
    comm = Comm.init()

def test_periodic_2D(load_periodic_opts_2D,
                     load_periodic_duct):
    ns = NumericalSolution.NS_base(load_periodic_duct[0],
                                   load_periodic_duct[1],
                                   load_periodic_duct[2],
                                   load_periodic_duct[3],
                                   opts)
    ns.calculateSolution('test_run')

def test_periodic_3D(load_periodic_opts_3D,
                     load_periodic_duct):
    ns = NumericalSolution.NS_base(load_periodic_duct[0],
                                   load_periodic_duct[1],
                                   load_periodic_duct[2],
                                   load_periodic_duct[3],
                                   opts)
    ns.calculateSolution('test_run')


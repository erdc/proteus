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
import pytest

from petsc4py import PETSc as p4pyPETSc

TestTools.addSubFolders( inspect.currentframe() )
import cavity2d
import twp_navier_stokes_cavity_2d_so
import twp_navier_stokes_cavity_2d_p
import twp_navier_stokes_cavity_2d_n

def clean_up_directory():
    FileList = ['forceHistory_p',
                'forceHistory_v',
                'cavity2D',
                'mesh',
                'momentHistory',
                'proteus',
                'wettedAreaHistory',
                'twp_navier_stokes_cavity_2d']
    mesh_ext = ['asy', 'edge', 'ele', 'neig', 'node',
                'ply', 'poly', 'txt', 'xmf', 'h5', 'log']
    TestTools.removeFiles(prefix_ext_tuple=(FileList,mesh_ext))

@pytest.fixture()
def load_cavity_problem(request):
    reload(cavity2d)
    reload(twp_navier_stokes_cavity_2d_so)
    reload(twp_navier_stokes_cavity_2d_p)
    reload(twp_navier_stokes_cavity_2d_n)
    pList = [twp_navier_stokes_cavity_2d_p]
    nList = [twp_navier_stokes_cavity_2d_n]
    so = twp_navier_stokes_cavity_2d_so
    so.sList = [default_s]
    yield pList, nList, so

@pytest.fixture()
def initialize_tp_pcd_options(request):
    petsc_options = p4pyPETSc.Options()
    petsc_options.setValue('rans2p_ksp_type','gmres')
    petsc_options.setValue('rans2p_ksp_gmres_restart',500)
    petsc_options.setValue('rans2p_ksp_atol',1e-20)
    petsc_options.setValue('rans2p_ksp_gmres_modifiedgramschmidt','')
    petsc_options.setValue('rans2p_pc_fieldsplit_type','schur')
    petsc_options.setValue('rans2p_pc_fieldsplit_schur_fact_type','upper')
    petsc_options.setValue('rans2p_pc_fieldsplit_schur_precondition','user')
    petsc_options.setValue('rans2p_fieldsplit_velocity_ksp_type','preonly')
    petsc_options.setValue('rans2p_fieldsplit_velocity_pc_type','lu')
    petsc_options.setValue('rans2p_fieldsplit_pressure_ksp_type','preonly')
    petsc_options.setValue('innerTPPCDsolver_Qp_visc_ksp_type','preonly')
    petsc_options.setValue('innerTPPCDsolver_Qp_visc_pc_type','lu')
    petsc_options.setValue('innerTPPCDsolver_Qp_dens_ksp_type','preonly')
    petsc_options.setValue('innerTPPCDsolver_Qp_dens_pc_type','lu')
    petsc_options.setValue('innerTPPCDsolver_Ap_rho_ksp_type','preonly')
    petsc_options.setValue('innerTPPCDsolver_Ap_rho_ksp_constant_null_space','')
    petsc_options.setValue('innerTPPCDsolver_Ap_rho_pc_type','hypre')
    petsc_options.setValue('innerTPPCDsolver_Ap_rho_pc_hypre_type','boomeramg')

@pytest.mark.Stabilization
@pytest.mark.LinearSolvers
def test_bochev_pressure_cavity(load_cavity_problem,
                                initialize_tp_pcd_options):
    initialize_tp_pcd_options
    ns =NumericalSolution.NS_base(load_cavity_problem[2],
                                  load_cavity_problem[0],
                                  load_cavity_problem[1],
                                  load_cavity_problem[2].sList,
                                  opts)
    ns.calculateSolution('bochev_pressure')
    script_dir = os.path.dirname(__file__)
    relpath = 'comparison_files/twp_navier_stokes_cavity_2d.h5'
    expected = tables.open_file(os.path.join(script_dir,relpath))
    actual = tables.open_file('twp_navier_stokes_cavity_2d.h5','r')
    numpy.allclose(expected.root.p_t1,
                   actual.root.p_t1)
    expected.close()
    actual.close()
    clean_up_directory()

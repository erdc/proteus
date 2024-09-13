#!/usr/bin/env python
""" Test modules for Driven Cavity Stokes preconditioners. """
import proteus.test_utils.TestTools as TestTools
import proteus
from proteus import LinearAlgebraTools as LAT
from proteus import LinearSolvers as LS
from proteus import Profiling
from proteus import NumericalSolution
from proteus import iproteus
#from proteus.iproteus import *

import os
import sys
import inspect
import numpy as np
import h5py
import pickle
import petsc4py
from petsc4py import PETSc
import pytest

proteus.test_utils.TestTools.addSubFolders( inspect.currentframe() )
from proteus import iproteus
from proteus import (defaults, default_p, default_n)
from importlib import reload
opts=iproteus.opts
import_modules = os.path.join(os.path.dirname(__file__),'import_modules')
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
    from proteus import default_s
    reload(PETSc)
    reload(iproteus)
    reload(default_p)
    reload(default_n)
    reload(default_s)
    Profiling.openLog("proteus.log",11)
    Profiling.verbose=True
    Context.contextOptionsString=context_options_str

    step2d_so = defaults.load_system('step2d_so',import_modules)
    twp_navier_stokes_step2d_p = defaults.load_physics('twp_navier_stokes_step2d_p',import_modules)
    twp_navier_stokes_step2d_n = defaults.load_numerics('twp_navier_stokes_step2d_n',import_modules)
    
    pList = [twp_navier_stokes_step2d_p]
    nList = [twp_navier_stokes_step2d_n]
    pList[0].name = 'step2d'        
    so = step2d_so
    so.name = pList[0].name
    so.sList = pList[0].name
    so.sList = [default_s]
    _scriptdir = os.path.dirname(__file__)
    Profiling.logAllProcesses = False
    ns = NumericalSolution.NS_base(so,
                                   pList,
                                   nList,
                                   so.sList,
                                   opts)
    return ns

def initialize_petsc_options():
    """Initializes schur complement petsc options. """
    petsc_options = PETSc.Options()
    petsc_options.clear()
    for k in petsc_options.getAll(): petsc_options.delValue(k)
    petsc_options.setValue('ksp_type', 'fgmres')
    petsc_options.setValue('ksp_rtol', 1.0e-8)
    petsc_options.setValue('ksp_atol', 1.0e-8)
    petsc_options.setValue('ksp_gmres_restart', 300)
    petsc_options.setValue('ksp_gmres_modifiedgramschmidt', 1)
    petsc_options.setValue('ksp_pc_side','right')
    petsc_options.setValue('pc_type', 'fieldsplit')
    petsc_options.setValue('pc_fieldsplit_type', 'schur')
    petsc_options.setValue('pc_fieldsplit_schur_fact_type', 'upper')
    petsc_options.setValue('pc_fieldsplit_schur_precondition', 'user')
    # Velocity block options
    petsc_options.setValue('fieldsplit_velocity_ksp_type', 'gmres')
    petsc_options.setValue('fieldsplit_velocity_ksp_gmres_modifiedgramschmidt', 1)
    petsc_options.setValue('fieldsplit_velocity_ksp_atol', 1e-5)
    petsc_options.setValue('fieldsplit_velocity_ksp_rtol', 1e-5)
    petsc_options.setValue('fieldsplit_velocity_ksp_pc_side', 'right')
    petsc_options.setValue('fieldsplit_velocity_fieldsplit_u_ksp_type', 'preonly')
    petsc_options.setValue('fieldsplit_velocity_fieldsplit_u_pc_type', 'hypre')
    petsc_options.setValue('fieldsplit_velocity_fieldsplit_u_pc_hypre_type', 'boomeramg')
    petsc_options.setValue('fieldsplit_velocity_fieldsplit_u_pc_hypre_boomeramg_coarsen_type', 'HMIS')
    petsc_options.setValue('fieldsplit_velocity_fieldsplit_v_ksp_type', 'preonly')
    petsc_options.setValue('fieldsplit_velocity_fieldsplit_v_pc_type', 'hypre')
    petsc_options.setValue('fieldsplit_velocity_fieldsplit_v_pc_hypre_type', 'boomeramg')
    petsc_options.setValue('fieldsplit_velocity_fieldsplit_v_pc_hypre_boomeramg_coarsen_type', 'HMIS')
    petsc_options.setValue('fieldsplit_velocity_fieldsplit_w_ksp_type', 'preonly')
    petsc_options.setValue('fieldsplit_velocity_fieldsplit_w_pc_type', 'hypre')
    petsc_options.setValue('fieldsplit_velocity_fieldsplit_w_pc_hypre_type', 'boomeramg')
    petsc_options.setValue('fieldsplit_velocity_fieldsplit_w_pc_hypre_boomeramg_coarsen_type', 'HMIS')
    #PCD Schur Complement options
    petsc_options.setValue('fieldsplit_pressure_ksp_type', 'preonly')
    petsc_options.setValue('innerTPPCDsolver_Qp_visc_ksp_type', 'preonly')
    petsc_options.setValue('innerTPPCDsolver_Qp_visc_pc_type', 'lu')
    petsc_options.setValue('innerTPPCDsolver_Qp_visc_pc_factor_mat_solver_type', 'superlu_dist')
    petsc_options.setValue('innerTPPCDsolver_Qp_dens_ksp_type', 'preonly')
    petsc_options.setValue('innerTPPCDsolver_Qp_dens_pc_type', 'lu')
    petsc_options.setValue('innerTPPCDsolver_Qp_dens_pc_factor_mat_solver_type', 'superlu_dist')
    petsc_options.setValue('innerTPPCDsolver_Ap_rho_ksp_type', 'richardson')
    petsc_options.setValue('innerTPPCDsolver_Ap_rho_ksp_max_it', 1)
    #petsc_options.setValue('innerTPPCDsolver_Ap_rho_ksp_constant_null_space',1)
    petsc_options.setValue('innerTPPCDsolver_Ap_rho_pc_type', 'hypre')
    petsc_options.setValue('innerTPPCDsolver_Ap_rho_pc_hypre_type', 'boomeramg')
    petsc_options.setValue('innerTPPCDsolver_Ap_rho_pc_hypre_boomeramg_strong_threshold', 0.5)
    petsc_options.setValue('innerTPPCDsolver_Ap_rho_pc_hypre_boomeramg_interp_type', 'ext+i-cc')
    petsc_options.setValue('innerTPPCDsolver_Ap_rho_pc_hypre_boomeramg_coarsen_type', 'HMIS')
    petsc_options.setValue('innerTPPCDsolver_Ap_rho_pc_hypre_boomeramg_agg_nl', 2)
    return petsc_options

def initialize_schur_ksp_obj(matrix_A, schur_approx):
    """
    Creates a right-hand-side and solution PETSc4Py vector for
    testing ksp solves.

    Parameters
    ----------
    matrix_A: :class:`PETSc.Mat`
        Global matrix object.
    schur_approx: :class:`LS.SchurPrecon`

    Returns
    -------
    ksp_obj: :class:`PETSc.KSP`
    """
    ksp_obj = PETSc.KSP().create()
    ksp_obj.setOperators(matrix_A,matrix_A)
    pc = schur_approx.pc
    ksp_obj.setPC(pc)
    ksp_obj.setFromOptions()
    pc.setFromOptions()
    pc.setOperators(matrix_A,matrix_A)
    pc.setUp()
    schur_approx.setUp(ksp_obj)
    ksp_obj.setUp()
    ksp_obj.pc.setUp()
    return ksp_obj

def runTest(ns, name):
    ns.calculateSolution(name)
    actual_log = TestTools.NumericResults.build_from_proteus_log('proteus.log')
    return actual_log

@pytest.mark.LinearSolvers
def test_step_slip_FullRun():
    """ Runs two-dimensional step problem with the settings:
        * Strongly enforced Free-slip BC.
        * Pressure Projection Stablization.
        * he = 0.05
    """
    petsc_options = initialize_petsc_options()
    context_options_str='he=0.05'
    ns = load_simulation(context_options_str)
    actual_log = runTest(ns,'test_1')

    L1 = actual_log.get_ksp_resid_it_info([(' step2d ',1.0,0,0)])
    L2 = actual_log.get_ksp_resid_it_info([(' step2d ',1.0,0,1)])
    L3 = actual_log.get_ksp_resid_it_info([(' step2d ',1.0,0,2)])
    print(L1,L2,L3)
    assert L1[0][1]==2
    assert L2[0][1]==11
    assert L3[0][1]<=10

@pytest.mark.LinearSolvers
def test_step_noslip_FullRun():
    """ Runs two-dimensional step problem with the settings:
        * Strongly enforced no-slip BC.
        * Pressure Projection Stablization.
        * he = 0.05
    """
    petsc_options = initialize_petsc_options
    context_options_str="boundary_condition_type='ns'"
    ns = load_simulation(context_options_str)
    actual_log = runTest(ns,'test_2')

    L1 = actual_log.get_ksp_resid_it_info([(' step2d ',1.0,0,0)])
    L2 = actual_log.get_ksp_resid_it_info([(' step2d ',1.0,0,1)])
    L3 = actual_log.get_ksp_resid_it_info([(' step2d ',1.0,0,2)])
    print(L1,L2,L3)
    assert L1[0][1]==2
    assert L2[0][1]==20
    assert L3[0][1]==23

@pytest.mark.LinearSolvers
def test_Schur_Sp_solve():
    """Tests a KSP solve using the Sp Schur complement approximation.
       For this test, the global matrix does not have a null space."""
    mat_A = load_matrix_step_noslip()
    petsc_options = initialize_petsc_options()
    b, x = create_petsc_vecs(mat_A)

    solver_info = LS.ModelInfo('interlaced', 3)
    schur_approx = LS.Schur_Sp(mat_A,
                               '',
                               solver_info=solver_info)
    ksp_obj = initialize_schur_ksp_obj(mat_A, schur_approx)
    ksp_obj.solve(b,x)

    assert ksp_obj.is_converged == True
    assert ksp_obj.reason == 2
    assert float(ksp_obj.norm) < 1.0e-5
    assert ksp_obj.its == 63
    
def create_petsc_vecs(matrix_A):
    """
    Creates a right-hand-side and solution PETSc vector for
    testing ksp solves.

    Parameters
    ----------
    matrix_A: :class:`PETSc.Mat`
        Global matrix object

    Returns
    -------
    vec_lst: tuple
        This is a list of :class:`pypyPETSc.Vec` where the first is
        a vector of ones (usually to act as a RHS-vector) while the
        second vector is a vector of zeros (usually to act as a
        storage vector for the solution).
    """
    b = PETSc.Vec().create()
    x = PETSc.Vec().create()
    b.createWithArray(np.ones(matrix_A.getSizes()[0][0],'d'))
    x.createWithArray(np.zeros(matrix_A.getSizes()[0][0],'d'))
    return (b, x)

def initialize_asm_ksp_obj(matrix_A):
    """
    Creates a right-hand-side and solution PETSc vector for
    testing ksp solves.

    Parameters
    ----------
    matrix_A: :class:`PETSc.Mat`
        Global matrix object.

    Returns
    -------
    ksp_obj: :class:`PETSc.KSP`
    """
    ksp_obj = PETSc.KSP().create()
    ksp_obj.setOperators(matrix_A,matrix_A)
    ksp_obj.setFromOptions()
    ksp_obj.setUp()
    return ksp_obj

def build_amg_index_sets(L_sizes):
    """
    Create PETSc index sets for the velocity components of a saddle
    point matrix

    Parameters
    ----------
    L_sizes :
        Sizes of original saddle-point system

    Returns:
    --------
    Index_Sets : lst
        List of velocity index sets
    """
    neqns = L_sizes[0][0]
    velocityDOF=[]
    for start in range(1,3):
        velocityDOF.append(np.arange(start=start,
                                     stop=1+neqns,
                                     step=3,
                                     dtype='i'))
    velocityDOF_full=np.vstack(velocityDOF).transpose().flatten()
    velocity_u_DOF = []
    velocity_u_DOF.append(np.arange(start=0,
                                    stop=2*neqns//3,
                                    step=2,
                                    dtype='i'))
    velocity_u_DOF_full = np.vstack(velocity_u_DOF).transpose().flatten()
    velocity_v_DOF = []
    velocity_v_DOF.append(np.arange(start=1,
                                    stop=1+2*neqns//3,
                                    step=2,
                                    dtype='i'))
    velocity_v_DOF_full = np.vstack(velocity_v_DOF).transpose().flatten()
    isvelocity = PETSc.IS()
    isvelocity.createGeneral(velocityDOF_full)
    isu = PETSc.IS()
    isu.createGeneral(velocity_u_DOF_full)
    isv = PETSc.IS()
    isv.createGeneral(velocity_v_DOF_full)
    return [isvelocity, isu, isv]

def clear_petsc_options():
    petsc_options = PETSc.Options()
    petsc_options.clear()
    for k in petsc_options.getAll(): petsc_options.delValue(k)

def initialize_velocity_block_petsc_options():
    petsc_options = PETSc.Options()
    petsc_options.clear()
    for k in petsc_options.getAll(): petsc_options.delValue(k)
    petsc_options.setValue('ksp_type','gmres')
    petsc_options.setValue('ksp_gmres_restart',100)
#    petsc_options.setValue('ksp_pc_side','right')
    petsc_options.setValue('ksp_atol',1e-8)
    petsc_options.setValue('ksp_gmres_modifiedgramschmidt','')
    petsc_options.setValue('pc_type','hypre')
    petsc_options.setValue('pc_type_hypre_type','boomeramg')
    return petsc_options

def initialize_velocity_block_petsc_options_2():
    petsc_options = PETSc.Options()
    petsc_options.clear()
    for k in petsc_options.getAll(): petsc_options.delValue(k)
    petsc_options.setValue('ksp_type','gmres')
    petsc_options.setValue('ksp_gmres_restart',100)
    petsc_options.setValue('ksp_atol',1e-8)
    petsc_options.setValue('ksp_gmres_modifiedgramschmidt','')
    return petsc_options

def initialize_velocity_block_petsc_options_3():
    petsc_options = PETSc.Options()
    petsc_options.clear()
    for k in petsc_options.getAll(): petsc_options.delValue(k)
    petsc_options.setValue('ksp_type','gmres')
    petsc_options.setValue('ksp_gmres_restart',100)
    petsc_options.setValue('ksp_pc_side','right')
    petsc_options.setValue('ksp_atol',1e-8)
    petsc_options.setValue('ksp_gmres_modifiedgramschmidt','')
    petsc_options.setValue('pc_type','hypre')
    petsc_options.setValue('pc_type_hypre_type','boomeramg')
    return petsc_options

def load_matrix_step_slip():
    """
    Loads a medium sized backwards facing step matrix for studying
    different AMG preconditioners.
    """
    A = LAT.petsc_load_matrix('dump_test_1_step2d_1.0par_j_0')
    return A

def load_matrix_step_noslip():
    """
    Loads a medium sized backwards facing step matrix for studying
    different AMG preconditioners.
    """
    A = LAT.petsc_load_matrix('dump_test_2_step2d_1.0par_j_0')
    return A

def test_amg_iteration_matrix_noslip():
    mat_A = load_matrix_step_noslip()
    petsc_options = initialize_velocity_block_petsc_options()
    L_sizes = mat_A.getSizes()
    index_sets = build_amg_index_sets(L_sizes)
    F_ksp = initialize_asm_ksp_obj(mat_A.createSubMatrix(index_sets[0],
                                                         index_sets[0]))
    b, x = create_petsc_vecs(mat_A.createSubMatrix(index_sets[0],
                                                   index_sets[0]))
    F_ksp.solve(b,x)
    assert F_ksp.its == 51

    PETSc.Options().setValue('pc_hypre_boomeramg_relax_type_all','sequential-Gauss-Seidel')
    F_ksp = initialize_asm_ksp_obj(mat_A.createSubMatrix(index_sets[0],
                                                         index_sets[0]))
    b, x = create_petsc_vecs(mat_A.createSubMatrix(index_sets[0],
                                                   index_sets[0]))

    F_ksp.solve(b,x)
    assert F_ksp.its == 55

    clear_petsc_options()
    initialize_velocity_block_petsc_options()

    PETSc.Options().setValue('pc_hypre_boomeramg_coarsen_type','PMIS')
    F_ksp = initialize_asm_ksp_obj(mat_A.createSubMatrix(index_sets[0],
                                                         index_sets[0]))
    b, x = create_petsc_vecs(mat_A.createSubMatrix(index_sets[0],
                                                   index_sets[0]))

    F_ksp.solve(b,x)
    assert F_ksp.its == 77

    clear_petsc_options()
    initialize_velocity_block_petsc_options()

    PETSc.Options().setValue('pc_hypre_boomeramg_relax_type_all','sequential-Gauss-Seidel')
    PETSc.Options().setValue('pc_hypre_boomeramg_coarsen_type','PMIS')
    F_ksp = initialize_asm_ksp_obj(mat_A.createSubMatrix(index_sets[0],
                                                         index_sets[0]))
    b, x = create_petsc_vecs(mat_A.createSubMatrix(index_sets[0],
                                                   index_sets[0]))

    F_ksp.solve(b,x)
    assert F_ksp.its == 88

def test_amg_iteration_matrix_slip():
    mat_A = load_matrix_step_slip()
    petsc_options = initialize_velocity_block_petsc_options_2()
    L_sizes = mat_A.getSizes()
    index_sets = build_amg_index_sets(L_sizes)

    F_ksp = initialize_asm_ksp_obj(mat_A.createSubMatrix(index_sets[0],
                                                         index_sets[0]))

    F_ksp.pc.setType('fieldsplit')
    F_ksp.pc.setFieldSplitIS(('v1',index_sets[1]),('v2',index_sets[2]))

    F_ksp.pc.getFieldSplitSubKSP()[0].setType('richardson')
    F_ksp.pc.getFieldSplitSubKSP()[1].setType('richardson')
    F_ksp.pc.getFieldSplitSubKSP()[0].pc.setType('hypre')
    F_ksp.pc.getFieldSplitSubKSP()[0].pc.setHYPREType('boomeramg')
    F_ksp.pc.getFieldSplitSubKSP()[1].pc.setType('hypre')
    F_ksp.pc.getFieldSplitSubKSP()[1].pc.setHYPREType('boomeramg')

    b, x = create_petsc_vecs(mat_A.createSubMatrix(index_sets[0],
                                                   index_sets[0]))
    F_ksp.solve(b,x)
    assert F_ksp.its == 7

@pytest.mark.amg
def test_amg_basic():
    mat_A = load_matrix_step_noslip()

    petsc_options = initialize_velocity_block_petsc_options()
    L_sizes = mat_A.getSizes()
    index_sets = build_amg_index_sets(L_sizes)

    #Initialize ksp object
    F_ksp = initialize_asm_ksp_obj(mat_A.createSubMatrix(index_sets[0],
                                                      index_sets[0]))
    b, x = create_petsc_vecs(mat_A.createSubMatrix(index_sets[0],
                                                index_sets[0]))   
    F_ksp.solve(b,x)
    assert F_ksp.its == 51

@pytest.mark.amg
def test_amg_iteration_performance():
    mat_A = load_matrix_step_noslip()
    petsc_options = initialize_velocity_block_petsc_options_2()
    L_sizes = mat_A.getSizes()
    index_sets = build_amg_index_sets(L_sizes)

    F_ksp = initialize_asm_ksp_obj(mat_A.createSubMatrix(index_sets[0],
                                                      index_sets[0]))
    b, x = create_petsc_vecs(mat_A.createSubMatrix(index_sets[0],
                                                index_sets[0]))

    F_ksp.solve(b,x)
    petsc_options.view()
    assert F_ksp.its == 174

@pytest.mark.amg
def test_amg_step_problem_noslip():
    mat_A = load_matrix_step_noslip()
    petsc_options = initialize_velocity_block_petsc_options_3()
    L_sizes = mat_A.getSizes()
    index_sets = build_amg_index_sets(L_sizes)

    F_ksp = initialize_asm_ksp_obj(mat_A.createSubMatrix(index_sets[0],
                                                      index_sets[0]))
    b, x = create_petsc_vecs(mat_A.createSubMatrix(index_sets[0],
                                                index_sets[0]))
    F_ksp.solve(b,x)
    assert F_ksp.its == 80

@pytest.mark.amg
def test_amg_step_problem_slip():
    mat_A = load_matrix_step_slip()
    petsc_options = initialize_velocity_block_petsc_options()
    L_sizes = mat_A.getSizes()
    index_sets = build_amg_index_sets(L_sizes)

    F_ksp = initialize_asm_ksp_obj(mat_A.createSubMatrix(index_sets[0],
                                                      index_sets[0]))
    b, x = create_petsc_vecs(mat_A.createSubMatrix(index_sets[0],
                                                index_sets[0]))
    F_ksp.solve(b,x)
    assert F_ksp.its == 68


if __name__ == '__main__':
    pass

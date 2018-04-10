#!/usr/bin/env python
"""

Test module for solver shell operator.

"""
import proteus.test_utils.TestTools
from proteus.iproteus import *
from proteus import Comm
from proteus import LinearAlgebraTools as LAT
from proteus import LinearSolvers as LS

import os
import sys
import inspect
import petsc4py as p4pyPETSc
import numpy as np
import pytest
import pickle
import time

from petsc4py import PETSc as p4pyPETSc

proteus.test_utils.TestTools.addSubFolders( inspect.currentframe() )

def create_petsc_vecs(matrix_A):
    """
    Creates a right-hand-side and solution PETSc4Py vector for
    testing ksp solves.

    Parameters
    ----------
    matrix_A: :class:`p4pyPETSc.Mat`
        Global matrix object

    Returns
    -------
    vec_lst: tuple
        This is a list of :class:`pypyPETSc.Vec` where the first is
        a vector of ones (usually to act as a RHS-vector) while the
        second vector is a vector of zeros (usually to act as a
        storage vector for the solution).
    """
    b = p4pyPETSc.Vec().create()
    x = p4pyPETSc.Vec().create()
    b.createWithArray(np.ones(matrix_A.getSizes()[0][0]))
    x.createWithArray(np.zeros(matrix_A.getSizes()[0][0]))
    return (b, x)

def initialize_asm_ksp_obj(matrix_A):
    """
    Creates a right-hand-side and solution PETSc4Py vector for
    testing ksp solves.

    Parameters
    ----------
    matrix_A: :class:`p4pyPETSc.Mat`
        Global matrix object.

    Returns
    -------
    ksp_obj: :class:`p4pyPETSc.KSP`
    """
    ksp_obj = p4pyPETSc.KSP().create()
    ksp_obj.setOperators(matrix_A,matrix_A)
    ksp_obj.setFromOptions()
    ksp_obj.setUp()
    return ksp_obj

def initialize_asm_ksp_obj_diagonal_pc(mat_A,
                                       index_sets):
    ksp_obj = p4pyPETSc.KSP().create()
    ksp_obj.setOperators(mat_A.getSubMatrix(index_sets[0],index_sets[0]),
                         mat_A.getSubMatrix(index_sets[0],index_sets[0]))
    ksp_obj.setFromOptions()
    ksp_obj.pc.setFieldSplitFields(2,['u_block',0],['v_block',1])
    ksp_obj.pc.getFieldSplitSubKSP()[0].pc.setType('lu')
    ksp_obj.pc.getFieldSplitSubKSP()[1].pc.setType('lu')    
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
    velocity_u_DOF.append(np.arange(start=1,
                                    stop=1+neqns,
                                    step=3,
                                    dtype='i'))
    velocity_u_DOF_full = np.vstack(velocity_u_DOF).transpose().flatten()
    velocity_v_DOF = []
    velocity_v_DOF.append(np.arange(start=2,
                                    stop=2+neqns,
                                    step=3,
                                    dtype='i'))
    velocity_v_DOF_full = np.vstack(velocity_v_DOF).transpose().flatten()
    isvelocity = p4pyPETSc.IS()
    isvelocity.createGeneral(velocityDOF_full)
    isu = p4pyPETSc.IS()
    isu.createGeneral(velocity_u_DOF_full)
    isv = p4pyPETSc.IS()
    isv.createGeneral(velocity_v_DOF_full)
    return [isvelocity, isu, isv]

def clear_petsc_options():
    for key in p4pyPETSc.Options().getAll():
        p4pyPETSc.Options().delValue(key)

def initialize_velocity_block_petsc_options():
    petsc_options = p4pyPETSc.Options()
    petsc_options.setValue('ksp_type','gmres')
    petsc_options.setValue('ksp_gmres_restart',100)
    petsc_options.setValue('ksp_pc_side','right')
    petsc_options.setValue('ksp_atol',1e-8)
    petsc_options.setValue('ksp_gmres_modifiedgramschmidt','')
    petsc_options.setValue('pc_type','hypre')
    petsc_options.setValue('pc_type_hypre_type','boomeramg')

def initialize_velocity_diagonal_petsc_options():
    petsc_options = p4pyPETSc.Options()
    petsc_options.setValue('ksp_type','fgmres')
    petsc_options.setValue('ksp_gmres_restart',100)
    petsc_options.setValue('ksp_pc_side','right')
    petsc_options.setValue('ksp_atol',1e-8)
    petsc_options.setValue('ksp_gmres_modifiedgramschmidt','')
    petsc_options.setValue('pc_type','fieldsplit')
#    petsc_options.setValue('pc_fieldsplit_detect_saddle_point','')
#    petsc_options.setValue('pc_fieldsplit_type','additive')
    
@pytest.fixture()
def load_saddle_point_matrix_1(request):
    """
    Loads a small example of a backwards facing step matrix for
    testing purposes. (Note: this matrix does not have advection)
    """
    A = LAT.petsc_load_matrix(os.path.join
                              (os.path.dirname(__file__),
                               'import_modules/saddle_point_mat_1'))
    yield A

@pytest.fixture()
def load_broad_crested_weir_matrix_1(request):
    """
    Loads a saddle-point matrix pulled from the broad-crested weir
    simulation.  See bcw_matrix_1_simulation_info.txt for simulation
    information.
    """
    A = LAT.petsc_load_matrix(os.path.join
                              (os.path.dirname(__file__),
                               'import_modules/bcw_matrix_1'))
    yield A
    
@pytest.mark.amg
def test_amg_bcw_matrix_1(load_broad_crested_weir_matrix_1):
    mat_A = load_broad_crested_weir_matrix_1
    petsc_options = initialize_velocity_block_petsc_options
    L_sizes = mat_A.getSizes()
    index_sets = build_amg_index_sets(L_sizes)
    trial_times = []

    # Trial 1
    F_ksp = initialize_asm_ksp_obj(mat_A.getSubMatrix(index_sets[0],
                                                      index_sets[0]))
    b, x = create_petsc_vecs(mat_A.getSubMatrix(index_sets[0],
                                                index_sets[0]))
    start = time.time()
    F_ksp.solve(b,x)
    end = time.time()
    trial_times.append(end-start)

    assert F_ksp.its==8

    # Trial 2
    clear_petsc_options()
    initialize_velocity_block_petsc_options()
    
    p4pyPETSc.Options().setValue('pc_hypre_boomeramg_relax_type_all','sequential-Gauss-Seidel')
    F_ksp = initialize_asm_ksp_obj(mat_A.getSubMatrix(index_sets[0],
                                                      index_sets[0]))
    b, x = create_petsc_vecs(mat_A.getSubMatrix(index_sets[0],
                                                index_sets[0]))

    start = time.time()
    F_ksp.solve(b,x)
    end = time.time()
    trial_times.append(end-start)

    assert F_ksp.its == 7

    # Trial 3
    clear_petsc_options()
    initialize_velocity_block_petsc_options()

    p4pyPETSc.Options().setValue('pc_hypre_boomeramg_coarsen_type','PMIS')
    F_ksp = initialize_asm_ksp_obj(mat_A.getSubMatrix(index_sets[0],
                                                      index_sets[0]))
    b, x = create_petsc_vecs(mat_A.getSubMatrix(index_sets[0],
                                                index_sets[0]))

    start = time.time()
    F_ksp.solve(b,x)
    end = time.time()
    trial_times.append(end-start)
    
    assert F_ksp.its == 6

    # Trial 4
    clear_petsc_options()
    initialize_velocity_block_petsc_options()

    p4pyPETSc.Options().setValue('pc_hypre_boomeramg_relax_type_all','sequential-Gauss-Seidel')
    p4pyPETSc.Options().setValue('pc_hypre_boomeramg_coarsen_type','PMIS')
    F_ksp = initialize_asm_ksp_obj(mat_A.getSubMatrix(index_sets[0],
                                                      index_sets[0]))
    b, x = create_petsc_vecs(mat_A.getSubMatrix(index_sets[0],
                                                index_sets[0]))

    start = time.time()
    F_ksp.solve(b,x)
    end = time.time()
    trial_times.append(end-start)

    assert F_ksp.its ==8

@pytest.mark.amg
@pytest.mark.skip
def test_amg_iteration_matrix_1(load_saddle_point_matrix_1):
    mat_A = load_saddle_point_matrix_1
    petsc_options = initialize_velocity_block_petsc_options
    L_sizes = mat_A.getSizes()
    index_sets = build_amg_index_sets(L_sizes)

    F_ksp = initialize_asm_ksp_obj(mat_A.getSubMatrix(index_sets[0],
                                                      index_sets[0]))
    b, x = create_petsc_vecs(mat_A.getSubMatrix(index_sets[0],
                                                index_sets[0]))
#    F_ksp.solve(b,x)
#    assert F_ksp.its == 76

    p4pyPETSc.Options().setValue('pc_hypre_boomeramg_relax_type_all','sequential-Gauss-Seidel')
    F_ksp = initialize_asm_ksp_obj(mat_A.getSubMatrix(index_sets[0],
                                                      index_sets[0]))
    b, x = create_petsc_vecs(mat_A.getSubMatrix(index_sets[0],
                                                index_sets[0]))

#    F_ksp.solve(b,x)
#    assert F_ksp.its == 80

    clear_petsc_options()
    initialize_velocity_block_petsc_options()

    p4pyPETSc.Options().setValue('pc_hypre_boomeramg_coarsen_type','PMIS')
    F_ksp = initialize_asm_ksp_obj(mat_A.getSubMatrix(index_sets[0],
                                                      index_sets[0]))
    b, x = create_petsc_vecs(mat_A.getSubMatrix(index_sets[0],
                                                index_sets[0]))

#    F_ksp.solve(b,x)
#    assert F_ksp.its = 197

    clear_petsc_options()
    initialize_velocity_block_petsc_options()

    p4pyPETSc.Options().setValue('pc_hypre_boomeramg_relax_type_all','sequential-Gauss-Seidel')
    p4pyPETSc.Options().setValue('pc_hypre_boomeramg_coarsen_type','PMIS')
    F_ksp = initialize_asm_ksp_obj(mat_A.getSubMatrix(index_sets[0],
                                                      index_sets[0]))
    b, x = create_petsc_vecs(mat_A.getSubMatrix(index_sets[0],
                                                index_sets[0]))

#    F_ksp.solve(b,x)
#    assert F_ksp.its = 231

@pytest.mark.amg
@pytest.mark.skip
def test_amg_iteration_matrix_1_block(load_saddle_point_matrix_1):
    mat_A = load_saddle_point_matrix_1
    L_sizes = mat_A.getSizes()
    index_sets = build_amg_index_sets(L_sizes)
    trial_times = []

    # Trial 1
    clear_petsc_options()
    initialize_velocity_diagonal_petsc_options()

    p4pyPETSc.Options().setValue('fieldsplit_u_block_ksp_type','preonly')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_ksp_type','preonly')
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_type', 'lu')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_type', 'lu')    
    
    F_ksp = initialize_asm_ksp_obj_diagonal_pc(mat_A,
                                               index_sets)
    b, x = create_petsc_vecs(mat_A.getSubMatrix(index_sets[0],
                                                index_sets[0]))
    start = time.time()
    F_ksp.solve(b,x)
    end = time.time()
    trial_times.append(end-start)
    
    assert F_ksp.its == 9

    # Trial 2
    clear_petsc_options()
    initialize_velocity_diagonal_petsc_options()

    p4pyPETSc.Options().setValue('fieldsplit_u_block_ksp_type','gmres')
    p4pyPETSc.Options().setValue('fieldsplit_u_block_ksp_gmres_modifiedgramschmidt','')
#    p4pyPETSc.Options().setValue('fieldsplit_u_block_ksp_pc_side','right')
    
    p4pyPETSc.Options().setValue('fieldsplit_v_block_ksp_type','gmres')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_ksp_gmres_modifiedgramschmidt','')
#    p4pyPETSc.Options().setValue('fieldsplit_v_block_ksp_pc_side','right')
    
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_type', 'hypre')    
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_type_hypre_type','boomeramg')
    
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_type', 'hypre')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_type_hypre_type','boomeramg')
    
    F_ksp = initialize_asm_ksp_obj_diagonal_pc(mat_A,
                                               index_sets)
    b, x = create_petsc_vecs(mat_A.getSubMatrix(index_sets[0],
                                                index_sets[0]))
    start = time.time()
    F_ksp.solve(b,x)
    end = time.time()
    trial_times.append(end - start)
    
    assert F_ksp.its == 9
    assert F_ksp.pc.getFieldSplitSubKSP()[0].its == 6
    assert F_ksp.pc.getFieldSplitSubKSP()[1].its == 5

    # Trial 3
    clear_petsc_options()
    initialize_velocity_diagonal_petsc_options()

    p4pyPETSc.Options().setValue('fieldsplit_u_block_ksp_type','gmres')
    p4pyPETSc.Options().setValue('fieldsplit_u_block_ksp_gmres_modifiedgramschmidt','')
    
    p4pyPETSc.Options().setValue('fieldsplit_v_block_ksp_type','gmres')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_ksp_gmres_modifiedgramschmidt','')
    
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_type', 'hypre')    
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_type_hypre_type','boomeramg')
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_hypre_boomeramg_relax_type_all','sequential-Gauss-Seidel')
    
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_type', 'hypre')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_type_hypre_type','boomeramg')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_hypre_boomeramg_relax_type_all','sequential-Gauss-Seidel')    
    
    F_ksp = initialize_asm_ksp_obj_diagonal_pc(mat_A,
                                               index_sets)
    b, x = create_petsc_vecs(mat_A.getSubMatrix(index_sets[0],
                                                index_sets[0]))
    start = time.time()
    F_ksp.solve(b,x)
    end = time.time()
    trial_times.append(end - start)

    assert F_ksp.its == 9
    assert F_ksp.pc.getFieldSplitSubKSP()[0].its == 6
    assert F_ksp.pc.getFieldSplitSubKSP()[1].its == 5

    # Trial 4
    clear_petsc_options()
    initialize_velocity_diagonal_petsc_options()

    p4pyPETSc.Options().setValue('fieldsplit_u_block_ksp_type','gmres')
    p4pyPETSc.Options().setValue('fieldsplit_u_block_ksp_gmres_modifiedgramschmidt','')
    p4pyPETSc.Options().setValue('fieldsplit_u_block_ksp_pc_side','right')
    
    p4pyPETSc.Options().setValue('fieldsplit_v_block_ksp_type','gmres')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_ksp_gmres_modifiedgramschmidt','')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_ksp_pc_side','right')
    
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_type', 'hypre')    
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_type_hypre_type','boomeramg')
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_hypre_boomeramg_relax_type_all','sequential-Gauss-Seidel')
    
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_type', 'hypre')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_type_hypre_type','boomeramg')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_hypre_boomeramg_relax_type_all','sequential-Gauss-Seidel')    
    
    F_ksp = initialize_asm_ksp_obj_diagonal_pc(mat_A,
                                               index_sets)
    b, x = create_petsc_vecs(mat_A.getSubMatrix(index_sets[0],
                                                index_sets[0]))
    start = time.time()
    F_ksp.solve(b,x)
    end = time.time()
    trial_times.append(end - start)

    assert F_ksp.its == 9
    assert F_ksp.pc.getFieldSplitSubKSP()[0].its == 6
    assert F_ksp.pc.getFieldSplitSubKSP()[1].its == 6
    
    # Trial 5
    clear_petsc_options()
    initialize_velocity_diagonal_petsc_options()

    p4pyPETSc.Options().setValue('fieldsplit_u_block_ksp_type','gmres')
    p4pyPETSc.Options().setValue('fieldsplit_u_block_ksp_gmres_modifiedgramschmidt','')
#    p4pyPETSc.Options().setValue('fieldsplit_u_block_ksp_pc_side','right')
    
    p4pyPETSc.Options().setValue('fieldsplit_v_block_ksp_type','gmres')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_ksp_gmres_modifiedgramschmidt','')
#    p4pyPETSc.Options().setValue('fieldsplit_v_block_ksp_pc_side','right')
    
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_type', 'hypre')    
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_type_hypre_type','boomeramg')
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_hypre_boomeramg_relax_type_all','sequential-Gauss-Seidel')
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_hypre_boomeramg_coarsen_type','PMIS')
    
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_type', 'hypre')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_type_hypre_type','boomeramg')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_hypre_boomeramg_relax_type_all','sequential-Gauss-Seidel')    
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_hypre_boomeramg_coarsen_type','PMIS')    
    
    F_ksp = initialize_asm_ksp_obj_diagonal_pc(mat_A,
                                               index_sets)
    b, x = create_petsc_vecs(mat_A.getSubMatrix(index_sets[0],
                                                index_sets[0]))
    start = time.time()
    F_ksp.solve(b,x)
    end = time.time()
    trial_times.append(end - start)

    assert F_ksp.its == 9
    assert F_ksp.pc.getFieldSplitSubKSP()[0].its == 17
    assert F_ksp.pc.getFieldSplitSubKSP()[1].its == 13

    # Trial 6
    clear_petsc_options()
    initialize_velocity_diagonal_petsc_options()

    p4pyPETSc.Options().setValue('fieldsplit_u_block_ksp_type','gmres')
    p4pyPETSc.Options().setValue('fieldsplit_u_block_ksp_gmres_modifiedgramschmidt','')
    p4pyPETSc.Options().setValue('fieldsplit_u_block_ksp_pc_side','right')
    
    p4pyPETSc.Options().setValue('fieldsplit_v_block_ksp_type','gmres')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_ksp_gmres_modifiedgramschmidt','')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_ksp_pc_side','right')
    
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_type', 'hypre')    
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_type_hypre_type','boomeramg')
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_hypre_boomeramg_relax_type_all','sequential-Gauss-Seidel')
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_hypre_boomeramg_coarsen_type','PMIS')
    
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_type', 'hypre')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_type_hypre_type','boomeramg')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_hypre_boomeramg_relax_type_all','sequential-Gauss-Seidel')    
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_hypre_boomeramg_coarsen_type','PMIS')    
    
    F_ksp = initialize_asm_ksp_obj_diagonal_pc(mat_A,
                                               index_sets)
    b, x = create_petsc_vecs(mat_A.getSubMatrix(index_sets[0],
                                                index_sets[0]))
    start = time.time()
    F_ksp.solve(b,x)
    end = time.time()
    trial_times.append(end - start)

    assert F_ksp.its == 9
    assert F_ksp.pc.getFieldSplitSubKSP()[0].its == 17
    assert F_ksp.pc.getFieldSplitSubKSP()[1].its == 14

    # Trial 7
    clear_petsc_options()
    initialize_velocity_diagonal_petsc_options()

    p4pyPETSc.Options().setValue('fieldsplit_u_block_ksp_type','gmres')
    p4pyPETSc.Options().setValue('fieldsplit_u_block_ksp_gmres_modifiedgramschmidt','')
    p4pyPETSc.Options().setValue('fieldsplit_u_block_ksp_pc_side','right')
    
    p4pyPETSc.Options().setValue('fieldsplit_v_block_ksp_type','gmres')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_ksp_gmres_modifiedgramschmidt','')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_ksp_pc_side','right')
    
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_type', 'hypre')    
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_type_hypre_type','boomeramg')
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_hypre_boomeramg_relax_type_all','l1scaled-SOR/Jacobi')
    
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_type', 'hypre')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_type_hypre_type','boomeramg')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_hypre_boomeramg_relax_type_all','l1scaled-SOR/Jacobi')    
    
    F_ksp = initialize_asm_ksp_obj_diagonal_pc(mat_A,
                                               index_sets)
    b, x = create_petsc_vecs(mat_A.getSubMatrix(index_sets[0],
                                                index_sets[0]))
    start = time.time()
    F_ksp.solve(b,x)
    end = time.time()
    trial_times.append(end - start)

    assert F_ksp.its == 9
    assert F_ksp.pc.getFieldSplitSubKSP()[0].its == 6
    assert F_ksp.pc.getFieldSplitSubKSP()[1].its == 5

@pytest.mark.amg
def test_amg_iteration_matrix_1_block(load_broad_crested_weir_matrix_1):
    mat_A = load_broad_crested_weir_matrix_1
    L_sizes = mat_A.getSizes()
    index_sets = build_amg_index_sets(L_sizes)
    trial_times = []

    # Trial 1
    clear_petsc_options()
    initialize_velocity_diagonal_petsc_options()

    p4pyPETSc.Options().setValue('fieldsplit_u_block_ksp_type','preonly')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_ksp_type','preonly')
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_type', 'lu')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_type', 'lu')    
    
    F_ksp = initialize_asm_ksp_obj_diagonal_pc(mat_A,
                                               index_sets)
    b, x = create_petsc_vecs(mat_A.getSubMatrix(index_sets[0],
                                                index_sets[0]))
    start = time.time()
    F_ksp.solve(b,x)
    end = time.time()
    trial_times.append(end-start)

    assert F_ksp.its == 11
    
    # Trial 2
    clear_petsc_options()
    initialize_velocity_diagonal_petsc_options()

    p4pyPETSc.Options().setValue('fieldsplit_u_block_ksp_type','gmres')
    p4pyPETSc.Options().setValue('fieldsplit_u_block_ksp_gmres_modifiedgramschmidt','')
#    p4pyPETSc.Options().setValue('fieldsplit_u_block_ksp_pc_side','right')
    
    p4pyPETSc.Options().setValue('fieldsplit_v_block_ksp_type','gmres')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_ksp_gmres_modifiedgramschmidt','')
#    p4pyPETSc.Options().setValue('fieldsplit_v_block_ksp_pc_side','right')
    
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_type', 'hypre')    
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_type_hypre_type','boomeramg')
    
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_type', 'hypre')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_type_hypre_type','boomeramg')
    
    F_ksp = initialize_asm_ksp_obj_diagonal_pc(mat_A,
                                               index_sets)
    b, x = create_petsc_vecs(mat_A.getSubMatrix(index_sets[0],
                                                index_sets[0]))
    start = time.time()
    F_ksp.solve(b,x)
    end = time.time()
    trial_times.append(end - start)
    
    assert F_ksp.its == 11
    assert F_ksp.pc.getFieldSplitSubKSP()[0].its == 3
    assert F_ksp.pc.getFieldSplitSubKSP()[1].its == 3

    # Trial 3
    clear_petsc_options()
    initialize_velocity_diagonal_petsc_options()

    p4pyPETSc.Options().setValue('fieldsplit_u_block_ksp_type','gmres')
    p4pyPETSc.Options().setValue('fieldsplit_u_block_ksp_gmres_modifiedgramschmidt','')
    
    p4pyPETSc.Options().setValue('fieldsplit_v_block_ksp_type','gmres')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_ksp_gmres_modifiedgramschmidt','')
    
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_type', 'hypre')    
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_type_hypre_type','boomeramg')
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_hypre_boomeramg_relax_type_all','sequential-Gauss-Seidel')
    
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_type', 'hypre')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_type_hypre_type','boomeramg')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_hypre_boomeramg_relax_type_all','sequential-Gauss-Seidel')    
    
    F_ksp = initialize_asm_ksp_obj_diagonal_pc(mat_A,
                                               index_sets)
    b, x = create_petsc_vecs(mat_A.getSubMatrix(index_sets[0],
                                                index_sets[0]))
    start = time.time()
    F_ksp.solve(b,x)
    end = time.time()
    trial_times.append(end - start)

    assert F_ksp.its == 11
    assert F_ksp.pc.getFieldSplitSubKSP()[0].its == 4
    assert F_ksp.pc.getFieldSplitSubKSP()[1].its == 3
    
    # Trial 4
    clear_petsc_options()
    initialize_velocity_diagonal_petsc_options()

    p4pyPETSc.Options().setValue('fieldsplit_u_block_ksp_type','gmres')
    p4pyPETSc.Options().setValue('fieldsplit_u_block_ksp_gmres_modifiedgramschmidt','')
    p4pyPETSc.Options().setValue('fieldsplit_u_block_ksp_pc_side','right')
    
    p4pyPETSc.Options().setValue('fieldsplit_v_block_ksp_type','gmres')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_ksp_gmres_modifiedgramschmidt','')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_ksp_pc_side','right')
    
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_type', 'hypre')    
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_type_hypre_type','boomeramg')
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_hypre_boomeramg_relax_type_all','sequential-Gauss-Seidel')
    
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_type', 'hypre')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_type_hypre_type','boomeramg')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_hypre_boomeramg_relax_type_all','sequential-Gauss-Seidel')    
    
    F_ksp = initialize_asm_ksp_obj_diagonal_pc(mat_A,
                                               index_sets)
    b, x = create_petsc_vecs(mat_A.getSubMatrix(index_sets[0],
                                                index_sets[0]))
    start = time.time()
    F_ksp.solve(b,x)
    end = time.time()
    trial_times.append(end - start)

    assert F_ksp.its == 11
    assert F_ksp.pc.getFieldSplitSubKSP()[0].its == 4
    assert F_ksp.pc.getFieldSplitSubKSP()[1].its == 3
    
    # Trial 5
    clear_petsc_options()
    initialize_velocity_diagonal_petsc_options()

    p4pyPETSc.Options().setValue('fieldsplit_u_block_ksp_type','gmres')
    p4pyPETSc.Options().setValue('fieldsplit_u_block_ksp_gmres_modifiedgramschmidt','')
#    p4pyPETSc.Options().setValue('fieldsplit_u_block_ksp_pc_side','right')
    
    p4pyPETSc.Options().setValue('fieldsplit_v_block_ksp_type','gmres')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_ksp_gmres_modifiedgramschmidt','')
#    p4pyPETSc.Options().setValue('fieldsplit_v_block_ksp_pc_side','right')
    
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_type', 'hypre')    
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_type_hypre_type','boomeramg')
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_hypre_boomeramg_relax_type_all','sequential-Gauss-Seidel')
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_hypre_boomeramg_coarsen_type','PMIS')
    
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_type', 'hypre')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_type_hypre_type','boomeramg')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_hypre_boomeramg_relax_type_all','sequential-Gauss-Seidel')    
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_hypre_boomeramg_coarsen_type','PMIS')    
    
    F_ksp = initialize_asm_ksp_obj_diagonal_pc(mat_A,
                                               index_sets)
    b, x = create_petsc_vecs(mat_A.getSubMatrix(index_sets[0],
                                                index_sets[0]))
    start = time.time()
    F_ksp.solve(b,x)
    end = time.time()
    trial_times.append(end - start)

    assert F_ksp.its == 11
    assert F_ksp.pc.getFieldSplitSubKSP()[0].its == 4
    assert F_ksp.pc.getFieldSplitSubKSP()[1].its == 4

    # Trial 6
    clear_petsc_options()
    initialize_velocity_diagonal_petsc_options()

    p4pyPETSc.Options().setValue('fieldsplit_u_block_ksp_type','gmres')
    p4pyPETSc.Options().setValue('fieldsplit_u_block_ksp_gmres_modifiedgramschmidt','')
    p4pyPETSc.Options().setValue('fieldsplit_u_block_ksp_pc_side','right')
    
    p4pyPETSc.Options().setValue('fieldsplit_v_block_ksp_type','gmres')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_ksp_gmres_modifiedgramschmidt','')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_ksp_pc_side','right')
    
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_type', 'hypre')    
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_type_hypre_type','boomeramg')
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_hypre_boomeramg_relax_type_all','sequential-Gauss-Seidel')
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_hypre_boomeramg_coarsen_type','PMIS')
    
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_type', 'hypre')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_type_hypre_type','boomeramg')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_hypre_boomeramg_relax_type_all','sequential-Gauss-Seidel')    
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_hypre_boomeramg_coarsen_type','PMIS')    
    
    F_ksp = initialize_asm_ksp_obj_diagonal_pc(mat_A,
                                               index_sets)
    b, x = create_petsc_vecs(mat_A.getSubMatrix(index_sets[0],
                                                index_sets[0]))
    start = time.time()
    F_ksp.solve(b,x)
    end = time.time()
    trial_times.append(end - start)

    assert F_ksp.its == 11
    assert F_ksp.pc.getFieldSplitSubKSP()[0].its == 4
    assert F_ksp.pc.getFieldSplitSubKSP()[1].its == 4

    # Trial 7
    clear_petsc_options()
    initialize_velocity_diagonal_petsc_options()

    p4pyPETSc.Options().setValue('fieldsplit_u_block_ksp_type','gmres')
    p4pyPETSc.Options().setValue('fieldsplit_u_block_ksp_gmres_modifiedgramschmidt','')
    p4pyPETSc.Options().setValue('fieldsplit_u_block_ksp_pc_side','right')
    
    p4pyPETSc.Options().setValue('fieldsplit_v_block_ksp_type','gmres')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_ksp_gmres_modifiedgramschmidt','')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_ksp_pc_side','right')
    
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_type', 'hypre')    
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_type_hypre_type','boomeramg')
    p4pyPETSc.Options().setValue('fieldsplit_u_block_pc_hypre_boomeramg_relax_type_all','l1scaled-SOR/Jacobi')
    
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_type', 'hypre')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_type_hypre_type','boomeramg')
    p4pyPETSc.Options().setValue('fieldsplit_v_block_pc_hypre_boomeramg_relax_type_all','l1scaled-SOR/Jacobi')    
    
    F_ksp = initialize_asm_ksp_obj_diagonal_pc(mat_A,
                                               index_sets)
    b, x = create_petsc_vecs(mat_A.getSubMatrix(index_sets[0],
                                                index_sets[0]))
    start = time.time()
    F_ksp.solve(b,x)
    end = time.time()
    trial_times.append(end - start)

    assert F_ksp.its == 11
    assert F_ksp.pc.getFieldSplitSubKSP()[0].its == 4
    assert F_ksp.pc.getFieldSplitSubKSP()[1].its == 3
    import pdb ; pdb.set_trace()
     
if __name__ == '__main__':
    pass

#!/usr/bin/env python
"""
Test module for linear boundary value problems in 2d

This module solves equations of the form

.. _math::

  \nabla \cdot \left( a(x) \nabla u \right) = f(x)

"""
from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
Profiling.logLevel=7
Profiling.verbose=True
def test_c0p1_serial():
    ###
    #import p and n files as modules 
    import poisson_het_2d_p
    import poisson_het_2d_c0pk_n
    #lists of p and n files (number of models long)
    pList = [poisson_het_2d_p]
    nList = [poisson_het_2d_c0pk_n]
    #split operator
    so = default_so
    so.name = pList[0].name = "poisson_2d_c0p1_serial"+"pe"+`comm.size()`
    so.sList=[default_s]
    #parun options
    opts.logLevel=7
    opts.verbose=True
    opts.profile=True
    ###update numerics for test
    nList[0].femSpaces = {0:default_n.C0_AffineLinearOnSimplexWithNodalBasis}
    nList[0].parallel = False
    #rectangular grid [nn x nn] with nLevels of refinement
    nList[0].nn = 21
    nList[0].nLevels = 2
    #Sparse LU
    nList[0].linearSolver=default_n.LU
    nList[0].multilevelLinearSolver=default_n.LU
    #make sure partitioning set correctly even in serial 
    nList[0].parallelPartitioningType = default_n.MeshParallelPartitioningTypes.element
    ###
    #create a numerical solution (manages solution process)
    ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
    #solve the problem
    failed = ns.calculateSolution('poisson_2d_c0p1_serial')
    assert(not failed)
def test_c0p1():
    ###
    #import p and n files as modules 
    import poisson_het_2d_p
    import poisson_het_2d_c0pk_n
    #lists of p and n files (number of models long)
    pList = [poisson_het_2d_p]
    nList = [poisson_het_2d_c0pk_n]
    #split operator
    so = default_so
    so.name = pList[0].name = "poisson_2d_c0p1"+"pe"+`comm.size()`
    so.sList=[default_s]
    #parun options
    opts.logLevel=7
    opts.verbose=True
    opts.profile=True
    ###update numerics for test
    #
    nList[0].femSpaces = {0:default_n.C0_AffineLinearOnSimplexWithNodalBasis}
    nList[0].parallel = True
    #rectangular grid [nn x nn] with nLevels of refinement. Must be 1 for parallel
    nList[0].nn = 21
    nList[0].nLevels = 1
    #PETSc solvers from petsc4py
    nList[0].linearSolver=default_n.KSP_petsc4py
    nList[0].multilevelLinearSolver=default_n.KSP_petsc4py
    #nodal partitioning
    nList[0].parallelPartitioningType = default_n.MeshParallelPartitioningTypes.node
    #numerical flux needed for boundary conditions in parallel
    nList[0].numericalFluxType = default_n.Advection_DiagonalUpwind_Diffusion_SIPG_exterior
    #velocity postprocessing
    nList[0].conservativeFlux = {0:'pwl'}
    ###create a numerical solution (manages solution process)
    ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
    #update the linear solver object with options not able to set from n file
    linearSolver = ns.nlsList[-1].solverList[-1].linearSolver
    assert 'ksp' in dir(linearSolver), "couldn't find ksp in {} with dir {}".format(linearSolver,dir(linearSolver))
    linearSolver.ksp.setType('cg')
    linearSolver.ksp.max_it = 5000
    ###
    failed = ns.calculateSolution('poisson_2d_c0p1')
    assert(not failed)
def test_c0p2_serial():
    ###
    #import p and n files as modules 
    import poisson_het_2d_p
    import poisson_het_2d_c0pk_n
    #lists of p and n files (number of models long)
    pList = [poisson_het_2d_p]
    nList = [poisson_het_2d_c0pk_n]
    #split operator 
    so = default_so
    so.name = pList[0].name = "poisson_2d_c0p2_serial"+"pe"+`comm.size()`
    so.sList=[default_s]
    #parun options
    opts.logLevel=7
    opts.verbose=True
    opts.profile=True
    ###update numerics for test
    nList[0].femSpaces = {0:default_n.C0_AffineQuadraticOnSimplexWithNodalBasis}
    nList[0].parallel = False
    #rectangular grid [nn x nn] with nLevels of refinement
    nList[0].nn = 21
    nList[0].nLevels = 2
    #Sparse LU
    nList[0].linearSolver=default_n.LU
    nList[0].multilevelLinearSolver=default_n.LU
    #make sure partitioning set correctly even in serial 
    nList[0].parallelPartitioningType = default_n.MeshParallelPartitioningTypes.element
    ###
    #create a numerical solution (manages solution process)
    ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
    #solve the problem
    failed = ns.calculateSolution('poisson_2d_c0p1_serial')
    assert(not failed)

def test_c0p2():
    ###
    #import p and n files as modules 
    import poisson_het_2d_p
    import poisson_het_2d_c0pk_n
    #lists of p and n files (number of models long)
    pList = [poisson_het_2d_p]
    nList = [poisson_het_2d_c0pk_n]
    #split operator
    so = default_so
    so.name = pList[0].name = "poisson_2d_c0p2"+"pe"+`comm.size()`
    so.sList=[default_s]
    #parun options
    opts.logLevel=7
    opts.verbose=True
    opts.profile=True
    ###update numerics for test
    #
    nList[0].femSpaces = {0:default_n.C0_AffineQuadraticOnSimplexWithNodalBasis}
    nList[0].parallel = True
    #rectangular grid [nn x nn] with nLevels of refinement. Must be 1 for parallel
    nList[0].nn = 21
    nList[0].nLevels = 1
    #PETSc solvers from petsc4py
    nList[0].linearSolver=default_n.KSP_petsc4py
    nList[0].multilevelLinearSolver=default_n.KSP_petsc4py
    #nodal partitioning
    nList[0].parallelPartitioningType = default_n.MeshParallelPartitioningTypes.node
    #numerical flux needed for boundary conditions in parallel
    nList[0].numericalFluxType = default_n.Advection_DiagonalUpwind_Diffusion_SIPG_exterior
    #velocity postprocessing
    nList[0].conservativeFlux = {0:'pwl-bdm'}
    #
    ###create a numerical solution (manages solution process)
    ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
    #update the linear solver object with options not able to set from n file 
    linearSolver = ns.nlsList[-1].solverList[-1].linearSolver
    assert 'ksp' in dir(linearSolver), "couldn't find ksp in {} with dir {}".format(linearSolver,dir(linearSolver))
    linearSolver.ksp.setType('fgmres')
    linearSolver.ksp.max_it = 10000
    #preconditioner (via petsc4py)
    linearSolver.ksp.getPC().setType('asm')
    linearSolver.ksp.getPC().setASMOverlap(2)
    linearSolver.ksp.getPC().setASMType(1) #3--basic, 1 -- restrict, 2 -- interpolate
    #
    failed = ns.calculateSolution('poisson_2d_c0p2')
    assert(not failed)

def test_ncp1_serial():
    ###
    #import p and n files as modules 
    import poisson_het_2d_p
    import poisson_het_2d_ncp1_n
    #lists of p and n files (number of models long)
    pList = [poisson_het_2d_p]
    nList = [poisson_het_2d_ncp1_n]
    #split operator
    so = default_so
    so.name = pList[0].name = "poisson_2d_ncp1_serial"+"pe"+`comm.size()`
    so.sList=[default_s]
    #parun options
    opts.logLevel=7
    opts.verbose=True
    opts.profile=True
    ###update numerics for test
    nList[0].femSpaces = {0:default_n.NC_AffineLinearOnSimplexWithNodalBasis}
    nList[0].parallel = False
    #rectangular grid [nn x nn] with nLevels of refinement
    nList[0].nn = 21
    nList[0].nLevels = 2
    #Sparse LU
    nList[0].linearSolver=default_n.LU
    nList[0].multilevelLinearSolver=default_n.LU
    #make sure partitioning set correctly even in serial 
    nList[0].parallelPartitioningType = default_n.MeshParallelPartitioningTypes.element
    ###
    #create a numerical solution (manages solution process)
    ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
    #solve the problem
    failed = ns.calculateSolution('poisson_2d_ncp1_serial')
    assert(not failed)
def test_ncp1():
    ###
    #import p and n files as modules 
    import poisson_het_2d_p
    import poisson_het_2d_ncp1_n
    #lists of p and n files (number of models long)
    pList = [poisson_het_2d_p]
    nList = [poisson_het_2d_ncp1_n]
    #split operator
    so = default_so
    so.name = pList[0].name = "poisson_2d_ncp1"+"pe"+`comm.size()`
    so.sList=[default_s]
    #parun options
    opts.logLevel=7
    opts.verbose=True
    opts.profile=True
    ###update numerics for test
    #
    nList[0].femSpaces = {0:default_n.NC_AffineLinearOnSimplexWithNodalBasis}
    nList[0].parallel = True
    #rectangular grid [nn x nn] with nLevels of refinement. Must be 1 for parallel
    nList[0].nn = 21
    nList[0].nLevels = 1
    #PETSc solvers from petsc4py
    nList[0].linearSolver=default_n.KSP_petsc4py
    nList[0].multilevelLinearSolver=default_n.KSP_petsc4py
    #nodal partitioning
    nList[0].parallelPartitioningType = default_n.MeshParallelPartitioningTypes.node
    #numerical flux needed for boundary conditions in parallel
    nList[0].numericalFluxType = default_n.Advection_DiagonalUpwind_Diffusion_SIPG_exterior
    #
    ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
    #update the linear solver object with options not able to set from n file
    linearSolver = ns.nlsList[-1].solverList[-1].linearSolver
    assert 'ksp' in dir(linearSolver), "couldn't find ksp in {} with dir {}".format(linearSolver,dir(linearSolver))
    linearSolver.ksp.setType('fgmres')
    linearSolver.ksp.max_it = 10000
    #preconditioner (via petsc4py)
    linearSolver.ksp.getPC().setType('asm')
    linearSolver.ksp.getPC().setASMOverlap(2)
    linearSolver.ksp.getPC().setASMType(1) #3--basic, 1 -- restrict, 2 -- interpolate
    #
    failed = ns.calculateSolution('poisson_2d_ncp1')
    assert(not failed)

def test_dgp1_serial(numerical_flux_flag):
    ###
    #import p and n files as modules 
    import poisson_het_2d_p
    import poisson_het_2d_dgpk_n
    #lists of p and n files (number of models long)
    pList = [poisson_het_2d_p]
    nList = [poisson_het_2d_dgpk_n]
    #split operator 
    so = default_so
    so.name = pList[0].name = "poisson_2d_dgp1_serial"+"pe"+`comm.size()`
    so.sList=[default_s]
    #parun options
    opts.logLevel=7
    opts.verbose=True
    opts.profile=True
    ###update numerics for test
    #pick the type of DG flux
    if numerical_flux_flag == 'SIPG':
        nList[0].numericalFluxType = default_n.Advection_DiagonalUpwind_Diffusion_SIPG
    elif numerical_flux_flag == 'IIPG':
        nList[0].numericalFluxType = default_n.Advection_DiagonalUpwind_Diffusion_IIPG
    elif numerical_flux_flag == 'LDG':
        nList[0].numericalFluxType = default_n.Diffusion_LDG
    else:
        nList[0].numericalFluxType = default_n.Advection_DiagonalUpwind_Diffusion_NIPG
    #
    nList[0].femSpaces = {0:default_n.DG_AffineLinearOnSimplexWithNodalBasis}
    nList[0].parallel = False
    #rectangular grid [nn x nn] with nLevels of refinement
    nList[0].nn = 21
    nList[0].nLevels = 2
    #Sparse LU
    nList[0].linearSolver=default_n.LU
    nList[0].multilevelLinearSolver=default_n.LU
    #make sure partitioning set correctly even in serial 
    nList[0].parallelPartitioningType = default_n.MeshParallelPartitioningTypes.element
    ###
    #create a numerical solution (manages solution process)
    ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
    #solve the problem
    failed = ns.calculateSolution('poisson_2d_dgp1_serial')
    assert(not failed)

def test_dgp1(numerical_flux_flag):
    ###
    #import p and n files as modules 
    import poisson_het_2d_p
    import poisson_het_2d_dgpk_n
    #lists of p and n files (number of models long)
    pList = [poisson_het_2d_p]
    nList = [poisson_het_2d_dgpk_n]
    #split operator 
    so = default_so
    so.name = pList[0].name = "poisson_2d_dgp1"+"pe"+`comm.size()`
    so.sList=[default_s]
    #parun options
    opts.logLevel=7
    opts.verbose=True
    opts.profile=True
    ###update numerics for test
    nList[0].femSpaces = {0:default_n.DG_AffineLinearOnSimplexWithNodalBasis}
    nList[0].parallel = True
    nList[0].nLayersOfOverlapForParallel = 1
    nList[0].parallelPartitioningType = default_n.MeshParallelPartitioningTypes.element
    #pick type of dg flux
    if numerical_flux_flag == 'SIPG':
        nList[0].numericalFluxType = default_n.Advection_DiagonalUpwind_Diffusion_SIPG
    elif numerical_flux_flag == 'IIPG':
        nList[0].numericalFluxType = default_n.Advection_DiagonalUpwind_Diffusion_IIPG
    elif numerical_flux_flag == 'LDG':
        nList[0].numericalFluxType = default_n.Diffusion_LDG
        nList[0].nLayersOfOverlapForParallel = 2
    else:
        nList[0].numericalFluxType = default_n.Advection_DiagonalUpwind_Diffusion_NIPG
    #rectangular grid [nn x nn] with nLevels of refinement (has to be 1 for parallel)
    nList[0].nn = 21
    nList[0].nLevels = 1
    #PETSc solvers from petsc4py
    nList[0].linearSolver=default_n.KSP_petsc4py
    nList[0].multilevelLinearSolver=default_n.KSP_petsc4py
    #
    ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
    #update the linear solver object with options not able to set from n file
    linearSolver = ns.nlsList[-1].solverList[-1].linearSolver
    assert 'ksp' in dir(linearSolver), "couldn't find ksp in {} with dir {}".format(linearSolver,dir(linearSolver))
    linearSolver.ksp.setType('fgmres')
    linearSolver.ksp.max_it = 10000
    linearSolver.ksp.atol = 1.0e-10; linearSolver.ksp.rtol=0.0
    #preconditioner (via petsc4py)
    linearSolver.ksp.getPC().setType('asm')
    linearSolver.ksp.getPC().setASMOverlap(2)
    linearSolver.ksp.getPC().setASMType(1) #3--basic, 1 -- restrict, 2 -- interpolate
    #
    failed = ns.calculateSolution('poisson_2d_dgp1')
    assert(not failed)

def test_dgp2_serial(numerical_flux_flag):
    ###
    #import p and n files as modules 
    import poisson_het_2d_p
    import poisson_het_2d_dgpk_n
    #lists of p and n files (number of models long)
    pList = [poisson_het_2d_p]
    nList = [poisson_het_2d_dgpk_n]
    #split operator 
    so = default_so
    so.name = pList[0].name = "poisson_2d_dgp2_serial"+"pe"+`comm.size()`
    so.sList=[default_s]
    #parun options
    opts.logLevel=7
    opts.verbose=True
    opts.profile=True
    ###update numerics for test
    #pick the type of DG flux
    if numerical_flux_flag == 'SIPG':
        nList[0].numericalFluxType = default_n.Advection_DiagonalUpwind_Diffusion_SIPG
    elif numerical_flux_flag == 'IIPG':
        nList[0].numericalFluxType = default_n.Advection_DiagonalUpwind_Diffusion_IIPG
    elif numerical_flux_flag == 'LDG':
        nList[0].numericalFluxType = default_n.Diffusion_LDG
    else:
        nList[0].numericalFluxType = default_n.Advection_DiagonalUpwind_Diffusion_NIPG
    #
    nList[0].femSpaces = {0:default_n.DG_AffineQuadraticOnSimplexWithNodalBasis}
    nList[0].parallel = False
    #rectangular grid [nn x nn] with nLevels of refinement
    nList[0].nn = 21
    nList[0].nLevels = 2
    #Sparse LU
    nList[0].linearSolver=default_n.LU
    nList[0].multilevelLinearSolver=default_n.LU
    nList[0].parallelPartitioningType = default_n.MeshParallelPartitioningTypes.element
    #velocity postprocessor
    nList[0].conservativeFlux = {0:'dg-point-eval'}
    ###
    #create a numerical solution (manages solution process)
    ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
    #solve the problem
    failed = ns.calculateSolution('poisson_2d_dgp2_serial')
    assert(not failed)

def test_dgp2(numerical_flux_flag):
    ###
    #import p and n files as modules 
    import poisson_het_2d_p
    import poisson_het_2d_dgpk_n
    #lists of p and n files (number of models long)
    pList = [poisson_het_2d_p]
    nList = [poisson_het_2d_dgpk_n]
    #split operator 
    so = default_so
    so.name = pList[0].name = "poisson_2d_dgp2"+"pe"+`comm.size()`
    so.sList=[default_s]
    #parun options
    opts.logLevel=7
    opts.verbose=True
    opts.profile=True
    ###update numerics for test
    nList[0].femSpaces = {0:default_n.DG_AffineQuadraticOnSimplexWithNodalBasis}
    nList[0].parallel = True
    nList[0].nLayersOfOverlapForParallel = 1
    nList[0].parallelPartitioningType = default_n.MeshParallelPartitioningTypes.element
    #pick type of dg flux
    if numerical_flux_flag == 'SIPG':
        nList[0].numericalFluxType = default_n.Advection_DiagonalUpwind_Diffusion_SIPG
    elif numerical_flux_flag == 'IIPG':
        nList[0].numericalFluxType = default_n.Advection_DiagonalUpwind_Diffusion_IIPG
    elif numerical_flux_flag == 'LDG':
        nList[0].numericalFluxType = default_n.Diffusion_LDG
        nList[0].nLayersOfOverlapForParallel = 2
    else:
        nList[0].numericalFluxType = default_n.Advection_DiagonalUpwind_Diffusion_NIPG
    #rectangular grid [nn x nn] with nLevels of refinement (has to be 1 for parallel)
    nList[0].nn = 21
    nList[0].nLevels = 1
    #PETSc solvers from petsc4py
    nList[0].linearSolver=default_n.KSP_petsc4py
    nList[0].multilevelLinearSolver=default_n.KSP_petsc4py
    #velocity postprocessor
    nList[0].conservativeFlux = {0:'dg-point-eval'}
    ###
    #create a numerical solution (manages solution process)
    ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
    #update the linear solver object with options not able to set from n file
    linearSolver = ns.nlsList[-1].solverList[-1].linearSolver
    assert 'ksp' in dir(linearSolver), "couldn't find ksp in {} with dir {}".format(linearSolver,dir(linearSolver))
    linearSolver.ksp.setType('fgmres')
    linearSolver.ksp.max_it = 10000
    linearSolver.ksp.atol = 1.0e-10; linearSolver.ksp.rtol=0.0
    #preconditioner (via petsc4py)
    linearSolver.ksp.getPC().setType('asm')
    linearSolver.ksp.getPC().setASMOverlap(2)
    linearSolver.ksp.getPC().setASMType(1) #3--basic, 1 -- restrict, 2 -- interpolate
    #
    failed = ns.calculateSolution('poisson_2d_dgp2')
    assert(not failed)

def test_dgp1_superlu_dist(numerical_flux_flag):
    ###
    #import p and n files as modules 
    import poisson_het_2d_p
    import poisson_het_2d_dgpk_n
    #lists of p and n files (number of models long)
    pList = [poisson_het_2d_p]
    nList = [poisson_het_2d_dgpk_n]
    #split operator 
    so = default_so
    so.name = pList[0].name = "poisson_2d_dgp1"+"pe"+`comm.size()`
    so.sList=[default_s]
    #parun options
    opts.logLevel=7
    opts.verbose=True
    opts.profile=True
    ###update numerics for test
    nList[0].femSpaces = {0:default_n.DG_AffineLinearOnSimplexWithNodalBasis}
    nList[0].parallel = True
    nList[0].nLayersOfOverlapForParallel = 1
    nList[0].parallelPartitioningType = default_n.MeshParallelPartitioningTypes.element
    #pick type of dg flux
    if numerical_flux_flag == 'SIPG':
        nList[0].numericalFluxType = default_n.Advection_DiagonalUpwind_Diffusion_SIPG
    elif numerical_flux_flag == 'IIPG':
        nList[0].numericalFluxType = default_n.Advection_DiagonalUpwind_Diffusion_IIPG
    elif numerical_flux_flag == 'LDG':
        nList[0].numericalFluxType = default_n.Diffusion_LDG
        nList[0].nLayersOfOverlapForParallel = 2
    else:
        nList[0].numericalFluxType = default_n.Advection_DiagonalUpwind_Diffusion_NIPG
    #rectangular grid [nn x nn] with nLevels of refinement (has to be 1 for parallel)
    nList[0].nn = 21
    nList[0].nLevels = 1
    #PETSc solvers from petsc4py
    nList[0].linearSolver=default_n.KSP_petsc4py
    nList[0].multilevelLinearSolver=default_n.KSP_petsc4py
    #
    ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
    #update the linear solver object with options not able to set from n file
    linearSolver = ns.nlsList[-1].solverList[-1].linearSolver
    assert 'ksp' in dir(linearSolver), "couldn't find ksp in {} with dir {}".format(linearSolver,dir(linearSolver))
    #solve with superlu_dist (for petsc this is preconditioner only)
    linearSolver.ksp.setType('preonly')
    linearSolver.ksp.atol = 1.0e-10; linearSolver.ksp.rtol=0.0
    #preconditioner (via petsc4py) is superlu_dist
    linearSolver.ksp.getPC().setType('lu')
    linearSolver.ksp.getPC().setFactorSolverPackage('superlu_dist')
    #
    failed = ns.calculateSolution('poisson_2d_dgp1_superlu_dist')
    assert(not failed)

def test_dgp2_superlu_dist(numerical_flux_flag):
    ###
    #import p and n files as modules 
    import poisson_het_2d_p
    import poisson_het_2d_dgpk_n
    #lists of p and n files (number of models long)
    pList = [poisson_het_2d_p]
    nList = [poisson_het_2d_dgpk_n]
    #split operator 
    so = default_so
    so.name = pList[0].name = "poisson_2d_dgp2"+"pe"+`comm.size()`
    so.sList=[default_s]
    #parun options
    opts.logLevel=7
    opts.verbose=True
    opts.profile=True
    ###update numerics for test
    nList[0].femSpaces = {0:default_n.DG_AffineQuadraticOnSimplexWithNodalBasis}
    nList[0].parallel = True
    nList[0].nLayersOfOverlapForParallel = 1
    nList[0].parallelPartitioningType = default_n.MeshParallelPartitioningTypes.element
    #pick type of dg flux
    if numerical_flux_flag == 'SIPG':
        nList[0].numericalFluxType = default_n.Advection_DiagonalUpwind_Diffusion_SIPG
    elif numerical_flux_flag == 'IIPG':
        nList[0].numericalFluxType = default_n.Advection_DiagonalUpwind_Diffusion_IIPG
    elif numerical_flux_flag == 'LDG':
        nList[0].numericalFluxType = default_n.Diffusion_LDG
        nList[0].nLayersOfOverlapForParallel = 2
    else:
        nList[0].numericalFluxType = default_n.Advection_DiagonalUpwind_Diffusion_NIPG
    #rectangular grid [nn x nn] with nLevels of refinement (has to be 1 for parallel)
    nList[0].nn = 21
    nList[0].nLevels = 1
    #PETSc solvers from petsc4py
    nList[0].linearSolver=default_n.KSP_petsc4py
    nList[0].multilevelLinearSolver=default_n.KSP_petsc4py
    #velocity postprocessor
    nList[0].conservativeFlux = {0:'dg-point-eval'}
    ###
    #create a numerical solution (manages solution process)
    ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
    #update the linear solver object with options not able to set from n file
    linearSolver = ns.nlsList[-1].solverList[-1].linearSolver
    assert 'ksp' in dir(linearSolver), "couldn't find ksp in {} with dir {}".format(linearSolver,dir(linearSolver))
    #solve with superlu_dist (for petsc this is preconditioner only)
    linearSolver.ksp.setType('preonly')
    linearSolver.ksp.atol = 1.0e-10; linearSolver.ksp.rtol=0.0
    #preconditioner (via petsc4py) is superlu_dist
    linearSolver.ksp.getPC().setType('lu')
    linearSolver.ksp.getPC().setFactorSolverPackage('superlu_dist')
    #solve the problem
    failed = ns.calculateSolution('poisson_2d_dgp2_superlu_dist')
    assert(not failed)


if __name__ == '__main__':
    import optparse

    usage = "usage: %prog [options] "
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("--parallel",
                      help="testing in parallel? ",
                      action="store_true",
                      dest="parallel",
                      default=False)
    parser.add_option("-F","--fem_space",
                      help="which finite element space to test: c0p1, c0p2, ncp1, dgp1, dgp2, all",
                      action="store",
                      default='c0p1')
    parser.add_option("-N","--dg_flux",
                      help="which dg flux to use: SIPG, IIPG, LDG, NIPG",
                      action="store",
                      default='NIPG')
    
    (run_opts,args) = parser.parse_args()

    if run_opts.parallel:
        if run_opts.fem_space in ['c0p2','all']:
            test_c0p2()
        if run_opts.fem_space in ['ncp1','all']:
            test_ncp1()
        if run_opts.fem_space in ['c0p1','all']:
            test_c0p1()
        if run_opts.fem_space in ['dgp1','all']:
            if run_opts.dg_flux == 'SIPG':
                test_dgp1_superlu_dist(run_opts.dg_flux)
            #currently fails if SIPG is chosen
            test_dgp1(run_opts.dg_flux)
        if run_opts.fem_space in ['dgp2','all']:
            if run_opts.dg_flux == 'SIPG':
                test_dgp2_superlu_dist(run_opts.dg_flux)
            test_dgp2(run_opts.dg_flux)
        
    else:
        if run_opts.fem_space in ['c0p2','all']:
            test_c0p2_serial()
        if run_opts.fem_space in ['ncp1','all']:
            test_ncp1_serial()
        if run_opts.fem_space in ['c0p1','all']:
            test_c0p1_serial()
        if run_opts.fem_space in ['dgp1','all']:
            test_dgp1_serial(run_opts.dg_flux)
        if run_opts.fem_space in ['dgp2','all']:
            test_dgp2_serial(run_opts.dg_flux)

    Profiling.logEvent("Closing Log")
    try:
        Profiling.closeLog()
    except:
        pass

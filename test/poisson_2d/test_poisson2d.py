#!/usr/bin/env python
"""
Test module for linear boundary value problems (serial)

This module solves equations of the form

.. _math::

  \nabla \cdot \left( a(x) \nabla u \right) = f(x)

"""
from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
Profiling.logLevel=7
Profiling.verbose=True
def test_c0p1():
    import poisson_het_2d_p
    import poisson_het_2d_c0pk_n
    pList = [poisson_het_2d_p]
    nList = [poisson_het_2d_c0pk_n]    
    so = default_so
    so.name = pList[0].name = "poisson_2d_c0p1"+"pe"+`comm.size()`
    so.sList=[default_s]
    opts.logLevel=7
    opts.verbose=True
    opts.profile=True
    opts.gatherArchive=True
    nList[0].femSpaces[0]  = default_n.C0_AffineLinearOnSimplexWithNodalBasis
    nList[0].linearSolver=default_n.KSP_petsc4py
    nList[0].multilevelLinearSolver=default_n.KSP_petsc4py
    nList[0].numericalFluxType = default_n.Advection_DiagonalUpwind_Diffusion_SIPG_exterior
    #nList[0].linearSolver=default_n.LU
    #nList[0].multilevelLinearSolver=default_n.LU
    ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
    ns.calculateSolution('poisson_2d_c0p1')

def test_c0p2():
    import poisson_het_2d_p
    import poisson_het_2d_c0pk_n
    pList = [poisson_het_2d_p]
    nList = [poisson_het_2d_c0pk_n]    
    so = default_so
    so.name = pList[0].name = "poisson_2d_c0p2"+"pe"+`comm.size()`
    so.sList=[default_s]
    opts.logLevel=7
    opts.verbose=True
    opts.profile=True
    opts.gatherArchive=True
    nList[0].femSpaces[0]  = default_n.C0_AffineQuadraticOnSimplexWithNodalBasis
    nList[0].linearSolver=default_n.KSP_petsc4py
    nList[0].multilevelLinearSolver=default_n.KSP_petsc4py
    #set ksp options
    from petsc4py import PETSc
    OptDB = PETSc.Options()
    OptDB.setValue('ksp_type','preonly')
    OptDB.setValue('pc_type','lu')
    OptDB.setValue('pc_factor_mat_solver_package','superlu_dist')
    nList[0].numericalFluxType = default_n.Advection_DiagonalUpwind_Diffusion_SIPG_exterior
    #nList[0].linearSolver=default_n.LU
    #nList[0].multilevelLinearSolver=default_n.LU
    ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
    ns.calculateSolution('poisson_2d_c0p2')

def test_load_vector():
    import poisson_het_2d_p
    import poisson_het_2d_c0pk_n
    pList = [poisson_het_2d_p]
    nList = [poisson_het_2d_c0pk_n]    
    so = default_so
    so.name = pList[0].name = "poisson_2d_c0p1"+"pe"+`comm.size()`
    so.sList=[default_s]
    opts.logLevel=7
    opts.verbose=True
    opts.profile=True
    opts.gatherArchive=True
    nList[0].femSpaces[0]  = default_n.C0_AffineLinearOnSimplexWithNodalBasis
    nList[0].linearSolver=default_n.LU
    nList[0].multilevelLinearSolver=default_n.LU
    ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
    ns.calculateSolution('poisson_2d_c0p1')
    
    #test load vector calculation in a crude way
    #see if sum over test functions of residual is same as the sum over the test functions of the load vector
    #should be if have strong bc's and no other terms in residual besides stiffness and load 
    #transport model on the final grid
    finest_model = ns.modelList[0].levelModelList[-1]
    #cheating a little bit to get the global test space dimension from the global trial space
    nr = finest_model.u[0].femSpace.dim
    import numpy as np
    import pdb
    pdb.set_trace()
    r = np.zeros((nr,),'d')
    f = np.zeros((nr,),'d')
    utmp = np.zeros((nr,),'d')
    finest_model.getResidual(utmp,r)
    finest_model.getLoadVector(f)
    r_sum = r.sum()
    f_sum = f.sum()
    assert abs(r_sum-f_sum) < 1.0e-6, "expected near equality got r_sum=%s and f_sum= %s " % (r_sum,f_sum)
if __name__ == '__main__':
    #test_c0p1()
    #test_c0p2()
    test_load_vector()
    Profiling.logEvent("Closing Log")
    try:
        Profiling.closeLog()
    except:
        pass

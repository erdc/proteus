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
    import sp_gw_p
    import sp_gw_c0p1_n
    pList = [sp_gw_p]
    nList = [sp_gw_c0p1_n]    
    so = default_so
    so.name = pList[0].name = "single_phase_gw_c0p1"+"pe"+`comm.size()`
    so.sList=[default_s]
    opts.logLevel=7
    opts.verbose=True
    opts.profile=True
    opts.gatherArchive=True
    nList[0].femSpaces[0]  = default_n.C0_AffineLinearOnSimplexWithNodalBasis
    nList[0].linearSolver=default_n.KSP_petsc4py
    nList[0].multilevelLinearSolver=default_n.KSP_petsc4py
    nList[0].numericalFluxType = default_n.Advection_DiagonalUpwind_Diffusion_SIPG_exterior
    #set ksp options
    from petsc4py import PETSc
    OptDB = PETSc.Options()
    OptDB.setValue('ksp_type','bcgsl')
    OptDB.setValue('pc_type','asm')
    OptDB.setValue('pc_asm_type','restrict')
    OptDB.setValue('pc_asm_overlap','2')
    OptDB.setValue('ksp_atol','1.0e-10')
    OptDB.setValue('ksp_rtol','0.0')
    #nList[0].linearSolver=default_n.LU
    #nList[0].multilevelLinearSolver=default_n.LU
    ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
    ns.calculateSolution('single_phase_gw')
    assert(True)
def test_ncp1():
    import sp_gw_p
    import sp_gw_ncp1_n
    pList = [sp_gw_p]
    nList = [sp_gw_ncp1_n]    
    so = default_so
    so.name = pList[0].name = "single_phase_gw_c0p1"+"pe"+`comm.size()`
    so.sList=[default_s]
    opts.logLevel=7
    opts.verbose=True
    opts.profile=True
    opts.gatherArchive=True
    nList[0].linearSolver=default_n.KSP_petsc4py
    nList[0].multilevelLinearSolver=default_n.KSP_petsc4py
    nList[0].numericalFluxType = default_n.Advection_DiagonalUpwind_Diffusion_SIPG_exterior 
    #set ksp options
    from petsc4py import PETSc
    OptDB = PETSc.Options()
    OptDB.setValue('ksp_type','preonly')
    OptDB.setValue('pc_type','lu')
    OptDB.setValue('pc_factor_mat_solver_package','superlu_dist')
    #nList[0].linearSolver=default_n.LU
    #nList[0].multilevelLinearSolver=default_n.LU
    ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
    ns.calculateSolution('single_phase_gw')
    assert(True)



if __name__ == '__main__':
    test_c0p1()
    test_ncp1()
    Profiling.logEvent("Closing Log")
    try:
        Profiling.closeLog()
    except:
        pass
#    test_c0q1()
#    test_c0q2()

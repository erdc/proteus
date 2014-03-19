#!/usr/bin/env python
r"""
Test module for linear boundary value problems (serial)

This module solves equations of the form

.. math::

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

def test_mass_and_stiff_jacobians():
    ###setup fine grid problem
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
    #for simplicity go ahead and compute the fine grid solutions
    ns.calculateSolution('single_phase_gw')
    ###get the mass matrix
    #transport model on the final grid
    finest_model = ns.modelList[0].levelModelList[-1]
    #make a copy of the transport model's jacobian
    mass_jacobian,mass_vals= LinearAlgebraTools.SparseMatCopy(finest_model.jacobian)
    stiff_jacobian,stiff_vals= LinearAlgebraTools.SparseMatCopy(finest_model.jacobian)
    #now try to assemble the mass jacobian
    finest_model.getMassJacobian(mass_jacobian)
    finest_model.getSpatialJacobian(stiff_jacobian)
    mass_jacobian.fwrite('mass_jacobian.out')
    stiff_jacobian.fwrite('stiff_jacobian.out')
    finest_model.jacobian.fwrite('full_jacobian.out')
    #just test to see that (mass/dt + stiff)*x = J*x
    dt = finest_model.timeIntegration.dt
    nr,nc=mass_jacobian.shape
    assert stiff_jacobian.shape[0] == nr
    assert stiff_jacobian.shape[1] == nc
    import numpy as np
    x=np.random.rand(nc); 
    m_split=np.zeros(nr,'d'); b_split=np.zeros(nr,'d')
    b_full=np.zeros(nr,'d')
    #J.x
    finest_model.jacobian.matvec(x,b_full)
    #M.x
    mass_jacobian.matvec(x,m_split)
    #M.x/dt
    m_split /= dt
    #K.x
    stiff_jacobian.matvec(x,b_split)
    b_split += m_split 
    diff = b_full-b_split
    err = np.dot(diff,diff)
    tol = 1.0e-6
    print 'err= ',err
    assert err < tol
    return mass_jacobian,stiff_jacobian,finest_model.jacobian


if __name__ == '__main__':
    test_c0p1()
    test_ncp1()
    Profiling.logEvent("Closing Log")
    try:
        Profiling.closeLog()
    except:
        pass

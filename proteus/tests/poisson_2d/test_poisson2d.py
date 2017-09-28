#!/usr/bin/env python
"""
Test module for linear boundary value problems (serial)

This module solves equations of the form

.. math::

  \nabla \cdot \left( a(x) \nabla u \right) = f(x)

"""
from proteus.iproteus import *
import os
import numpy as np
import poisson_het_2d_p
import poisson_het_2d_c0pk_n

class TestPoisson2D():

    @classmethod
    def setup_class(cls):
        pass

    @classmethod
    def teardown_class(cls):
        pass

    def setup_method(self,method):
        self.aux_names = []
        reload(poisson_het_2d_p)
        self.meshdir = os.path.dirname(os.path.abspath(__file__))
        
    def teardown_method(self,method):
        filenames = []
        for aux_name in self.aux_names:
            filenames.extend([aux_name+'.'+ext for ext in ['h5','xmf']])
        filenames.append('proteus.log')
        for f in filenames:
            if os.path.exists(f):
                try:
                    os.remove(f)
                except OSError, e:
                    print ("Error: %s - %s" %(e.filename,e.strerror))
            else:
                pass
            
    def test_c0p1(self,use_strong_constraints=False):
        reload(poisson_het_2d_c0pk_n)
        pList = [poisson_het_2d_p]
        nList = [poisson_het_2d_c0pk_n]
        reload(default_so)
        so = default_so
        so.name = pList[0].name = "poisson_2d_c0p1"+"pe"+`comm.size()`
        reload(default_s)
        so.sList=[default_s]
        opts.logLevel=7
        opts.verbose=True
        opts.profile=True
        opts.gatherArchive=True
        nList[0].femSpaces[0]  = default_n.C0_AffineLinearOnSimplexWithNodalBasis
        nList[0].linearSolver=default_n.KSP_petsc4py
        nList[0].multilevelLinearSolver=default_n.KSP_petsc4py
        nList[0].numericalFluxType = default_n.Advection_DiagonalUpwind_Diffusion_SIPG_exterior
        soln_name = 'poisson_2d_c0p1'
        if use_strong_constraints == True:
            nList[0].numericalFluxType = nList[0].Exterior_StrongFlux
            soln_name = 'poisson_2d_c0p1_strong_dirichlet'
        ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
        ns.calculateSolution(soln_name)
        self.aux_names.append(pList[0].name)
        del ns

    def test_c0p1_strong(self):
        return self.test_c0p1(use_strong_constraints=True)

    def test_c0p2(self,
                  use_strong_constraints=False,
                  test_superlu=False):
        reload(poisson_het_2d_c0pk_n)
        pList = [poisson_het_2d_p]
        nList = [poisson_het_2d_c0pk_n]
        reload(default_so)
        so = default_so
        so.name = pList[0].name = "poisson_2d_c0p2"+"pe"+`comm.size()`
        reload(default_s)
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
        if test_superlu == True:
            OptDB.setValue('ksp_type','preonly')
            OptDB.setValue('pc_type','lu')
            OptDB.setValue('pc_factor_mat_solver_package','superlu_dist')
        nList[0].numericalFluxType = default_n.Advection_DiagonalUpwind_Diffusion_SIPG_exterior
        soln_name = 'poisson_2d_c0p2'
        if use_strong_constraints == True:
            nList[0].numericalFluxType = nList[0].Exterior_StrongFlux
            soln_name = 'poisson_2d_c0p2_strong_dirichlet'
        #nList[0].linearSolver=default_n.LU
        #nList[0].multilevelLinearSolver=default_n.LU
        ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
        ns.calculateSolution(soln_name)
        self.aux_names.append(pList[0].name)
        del ns

    def test_c0p2_strong(self):
        return self.test_c0p2(use_strong_constraints=True)
    def test_2dm(self,use_strong_constraints=False,nLevels=3):
        reload(poisson_het_2d_p)
        reload(poisson_het_2d_c0pk_n)
        poisson_het_2d_p.meshfile=os.path.join(self.meshdir,'square4x4')
        assert os.path.isfile(poisson_het_2d_p.meshfile+'.3dm'), "could not find meshfile {0}".format(poisson_het_2d_p.meshfile+'.3dm')
        poisson_het_2d_p.x0 = (0.0,0.0,0.)
        poisson_het_2d_p.L  = (1.,1.,1.)
        poisson_het_2d_c0pk_n.nLevels = nLevels
        pList = [poisson_het_2d_p]
        nList = [poisson_het_2d_c0pk_n]
        reload(default_so)
        so = default_so
        so.name = pList[0].name = "poisson_2d_c0p1_2dm"+"pe"+`comm.size()`
        reload(default_s)
        so.sList=[default_s]
        opts.logLevel=7
        opts.verbose=True
        opts.profile=True
        opts.gatherArchive=True
        nList[0].femSpaces[0]  = default_n.C0_AffineLinearOnSimplexWithNodalBasis
        nList[0].linearSolver=default_n.KSP_petsc4py
        nList[0].multilevelLinearSolver=default_n.KSP_petsc4py
        nList[0].numericalFluxType = default_n.Advection_DiagonalUpwind_Diffusion_SIPG_exterior
        soln_name = 'poisson_2d_c0p1_2dm'
        if use_strong_constraints == True:
            nList[0].numericalFluxType = nList[0].Exterior_StrongFlux
        soln_name = 'poisson_2d_c0p1_strong_dirichlet'
        #nList[0].linearSolver=default_n.LU
        #nList[0].multilevelLinearSolver=default_n.LU
        ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
        ns.calculateSolution(soln_name)
        self.aux_names.append(pList[0].name)
        del ns

def compute_load_vector(use_weak_dirichlet=False):
    import poisson_het_2d_p
    reload(poisson_het_2d_p)
    import poisson_het_2d_c0pk_n
    pList = [poisson_het_2d_p]
    nList = [poisson_het_2d_c0pk_n]
    reload(default_so)
    so = default_so
    so.name = pList[0].name = "poisson_2d_c0p1"+"pe"+`comm.size()`
    reload(default_s)
    so.sList=[default_s]
    opts.logLevel=7
    opts.verbose=True
    opts.profile=True
    opts.gatherArchive=True
    nList[0].femSpaces[0]  = default_n.C0_AffineLinearOnSimplexWithNodalBasis
    nList[0].linearSolver=default_n.LU
    nList[0].multilevelLinearSolver=default_n.LU
    if use_weak_dirichlet:
        nList[0].linearSolver=default_n.KSP_petsc4py
        nList[0].multilevelLinearSolver=default_n.KSP_petsc4py
        nList[0].numericalFluxType = default_n.Advection_DiagonalUpwind_Diffusion_SIPG_exterior

    ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
    ns.calculateSolution('poisson_2d_c0p1')
    #test load vector calculation in a crude way
    #see if sum over test functions of residual is same as the sum over the test functions of the load vector
    #should be if have strong bc's and no other terms in residual besides stiffness and load
    #transport model on the final grid
    finest_model = ns.modelList[0].levelModelList[-1]
    #cheating a little bit to get the global test space dimension from the global trial space
    nr = finest_model.u[0].femSpace.dim
    r = np.zeros((nr,),'d')
    f = np.zeros((nr,),'d')
    utmp = np.zeros((nr,),'d')
    finest_model.getResidual(utmp,r)
    finest_model.getLoadVector(f)
    del ns
    return r,f

def notest_load_vector():
    for name,use_num_flux in zip(['Strong_Dir','Weak_Dir'],[False,True]):
        r,f = compute_load_vector(use_num_flux)
        test = npt.assert_almost_equal
        test.descrption = 'test_load_vector_{}'.format(name)
        yield test,r,f
        
if __name__ == '__main__':
    pass

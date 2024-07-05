#!/usr/bin/env python
"""
Test module for linear boundary value problems (serial)

This module solves equations of the form

.. _math::

  \nabla \cdot \left( a(x) \nabla u \right) = f(x)

"""
import os
import numpy as np
import numpy.testing as npt
import pytest

from proteus.iproteus import *

from . import sp_gw_p
from . import sp_gw_c0p1_n
from . import sp_gw_ncp1_n

class TestSinglePhaseGW(object):

    @classmethod
    def setup_class(cls):
        pass

    @classmethod
    def teardown_class(cls):
        pass

    def setup_method(self,method):
        self.aux_names = []

    def teardown_method(self,method):
        filenames = []
        for aux_name in self.aux_names:
            filenames.extend([aux_name+'.'+ext for ext in ['h5','xmf']])
        filenames.append('proteus.log')
        for f in filenames:
            if os.path.exists(f):
                try:
                    os.remove(f)
                except OSError as e:
                    print ("Error: %s - %s" %(e.filename,e.strerror))
            else:
                pass

    def test_c0p1(self):
        reload(sp_gw_p)
        reload(sp_gw_c0p1_n)
        pList = [sp_gw_p]
        nList = [sp_gw_c0p1_n]
        so = default_so
        sim_name = "single_phase_gw_c0p1"
        so.name = pList[0].name = sim_name +"pe"+repr(comm.size())
        so.sList=[default_s]
        opts.logLevel=7
        opts.verbose=True
        opts.profile=True
        opts.gatherArchive=True
        nList[0].femSpaces[0]  = default_n.C0_AffineLinearOnSimplexWithNodalBasis
        nList[0].numericalFluxType = default_n.Advection_DiagonalUpwind_Diffusion_SIPG_exterior
        nList[0].linearSolver=default_n.LU
        nList[0].multilevelLinearSolver=default_n.LU
        ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
        failed = ns.calculateSolution(sim_name)
        self.aux_names.append(pList[0].name)
        assert (not failed)

    @pytest.mark.skipif(os.getenv('TEST_PROFILE')=="proteus-conda-osx", reason="need to fix locally on osx")
    def test_ncp1(self):
        reload(sp_gw_p)
        reload(sp_gw_ncp1_n)
        pList = [sp_gw_p]
        nList = [sp_gw_ncp1_n]
        so = default_so
        sim_name = "single_phase_gw_ncp1"
        so.name = pList[0].name = sim_name +"pe"+repr(comm.size())
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
        failed = ns.calculateSolution(sim_name)
        self.aux_names.append(pList[0].name)
        assert (not failed)

    def test_mass_and_stiff_jacobians(self):
        ###setup fine grid problem
        reload(sp_gw_p)
        reload(sp_gw_c0p1_n)
        pList = [sp_gw_p]
        nList = [sp_gw_c0p1_n]
        pList[0].ndays=500
        pList[0].nDTout=1
        so = default_so
        sim_name = "single_phase_gw_c0p1"
        so.name = pList[0].name = sim_name +"pe"+repr(comm.size())
        self.aux_names.append(pList[0].name)
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
        mass_jacobian = finest_model.initializeMassJacobian()
        stiff_jacobian= finest_model.initializeSpatialJacobian()
        #now try to assemble the mass jacobian
        finest_model.getMassJacobian(mass_jacobian)
        finest_model.getSpatialJacobian(stiff_jacobian)

        #just test to see that (mass/dt + stiff)*x = J*x
        dt = finest_model.timeIntegration.dt
        nr,nc=mass_jacobian.shape
        assert stiff_jacobian.shape[0] == nr
        assert stiff_jacobian.shape[1] == nc
        import numpy as np
        #get constant portion of the load vector
        f = np.zeros((nr,),'d')
        finest_model.getLoadVector(f)

        ntries=5
        for i in range(ntries):
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
            npt.assert_almost_equal(b_split,b_full)

if __name__ == '__main__':
    pass

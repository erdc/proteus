#!/usr/bin/env python
"""
Test module for surface tension
"""
from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
Profiling.logLevel=1
Profiling.verbose=True
import os
import numpy as np
import tables
import pytest
from proteus import default_so
from . import (parameters,
               risingBubble_so, risingBubble,
               vof_p,
               vof_n,
               ls_p,
               ls_n,
               redist_p,
               redist_n,               
               ls_consrv_p,
               ls_consrv_n,
               twp_navier_stokes_p,
               twp_navier_stokes_n,
               pressureincrement_p,
               pressureincrement_n,
               pressure_p,
               pressure_n,
               pressureInitial_p,
               pressureInitial_n)


class TestSurfaceTension():

    @classmethod
    def setup_class(cls):
        pass
    
    @classmethod
    def teardown_class(cls):
        pass
    
    def reload_modules(self):
        reload(default_so)
        reload(risingBubble)
        reload(risingBubble_so)
        reload(vof_p)
        reload(vof_n)
        reload(ls_p)
        reload(ls_n)
        reload(redist_p)
        reload(redist_n)               
        reload(ls_consrv_p)
        reload(ls_consrv_n)
        reload(twp_navier_stokes_p)
        reload(twp_navier_stokes_n)
        reload(pressureincrement_p)
        reload(pressureincrement_n)
        reload(pressure_p)
        reload(pressure_n)
        reload(pressureInitial_p)
        reload(pressureInitial_n)
        
    def setup_method(self,method):
        self._scriptdir = os.path.dirname(__file__)
        
    def teardown_method(self,method):
        pass
    
    def test_2D_with_supg(self):
        # Set parameters for this test
        parameters.ct.USE_SUPG_NS=1
        parameters.ct.ARTIFICIAL_VISCOSITY_NS=1
        parameters.ct.nd=2
        # RELOAD MODULES
        self.reload_modules()
        pnList = [(vof_p,               vof_n),
                  (ls_p,                ls_n),
                  (redist_p,            redist_n),
                  (ls_consrv_p,         ls_consrv_n),
                  (twp_navier_stokes_p, twp_navier_stokes_n),
                  (pressureincrement_p, pressureincrement_n),
                  (pressure_p,          pressure_n),
                  (pressureInitial_p,   pressureInitial_n)]
        self.so = risingBubble_so
        pList=[]
        nList=[]
        sList=[]
        for (pModule,nModule) in pnList:
            pList.append(pModule)
            if pList[-1].name == None:
                pList[-1].name = pModule.__name__
            nList.append(nModule)
        for i in range(len(pnList)):
            sList.append(default_s)
        self.so.name += "_2D_supg"
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               pList,
                                               nList,
                                               sList,
                                               opts)
        ns.calculateSolution('2D_supg')
        # COMPARE VS SAVED FILES #
        expected_path = 'comparison_files/risingBubble_2D_supg.h5'
        expected = tables.open_file(os.path.join(self._scriptdir,expected_path))
        actual = tables.open_file('risingBubble_2D_supg.h5','r')
        assert np.allclose(expected.root.phi_t2,actual.root.phi_t2,atol=1e-10)
        assert np.allclose(expected.root.p_t2,actual.root.p_t2,atol=1e-10)
        assert np.allclose(expected.root.velocity_t2,actual.root.velocity_t2,atol=1e-7)        
        expected.close()
        actual.close()

    def test_2D_with_EV(self):
        # Set parameters for this test
        parameters.ct.USE_SUPG_NS=0
        parameters.ct.ARTIFICIAL_VISCOSITY_NS=2
        parameters.ct.nd=2 
        # RELOAD MODULES
        self.reload_modules()
        pnList = [(vof_p,               vof_n),
                  (ls_p,                ls_n),
                  (redist_p,            redist_n),
                  (ls_consrv_p,         ls_consrv_n),
                  (twp_navier_stokes_p, twp_navier_stokes_n),
                  (pressureincrement_p, pressureincrement_n),
                  (pressure_p,          pressure_n),
                  (pressureInitial_p,   pressureInitial_n)]
        self.so = risingBubble_so
        pList=[]
        nList=[]
        sList=[]
        for (pModule,nModule) in pnList:
            pList.append(pModule)
            if pList[-1].name == None:
                pList[-1].name = pModule.__name__
            nList.append(nModule)
        for i in range(len(pnList)):            
            sList.append(default_s)                           
        self.so.name += "_2D_ev"
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               pList,
                                               nList,
                                               sList,
                                               opts)
        ns.calculateSolution('2D_ev')
        # COMPARE VS SAVED FILES #
        expected_path = 'comparison_files/risingBubble_2D_ev.h5'
        expected = tables.open_file(os.path.join(self._scriptdir,expected_path))
        actual = tables.open_file('risingBubble_2D_ev.h5','r')
        assert np.allclose(expected.root.phi_t2,actual.root.phi_t2,atol=1e-10)
        assert np.allclose(expected.root.p_t2,actual.root.p_t2,atol=1e-10)
        assert np.allclose(expected.root.velocity_t2,actual.root.velocity_t2,atol=1e-7)        
        expected.close()
        actual.close()

    def test_3D_with_SUPG(self):
        # Set parameters for this test
        parameters.ct.USE_SUPG_NS=1
        parameters.ct.ARTIFICIAL_VISCOSITY_NS=1
        parameters.ct.nd=3
        # RELOAD MODULES
        self.reload_modules()
        pnList = [(vof_p,               vof_n),
                  (ls_p,                ls_n),
                  (redist_p,            redist_n),
                  (ls_consrv_p,         ls_consrv_n),
                  (twp_navier_stokes_p, twp_navier_stokes_n),
                  (pressureincrement_p, pressureincrement_n),
                  (pressure_p,          pressure_n),
                  (pressureInitial_p,   pressureInitial_n)]
        self.so = risingBubble_so
        pList=[]
        nList=[]
        sList=[]
        for (pModule,nModule) in pnList:
            pList.append(pModule)
            if pList[-1].name == None:
                pList[-1].name = pModule.__name__
            nList.append(nModule)
        for i in range(len(pnList)):            
            sList.append(default_s)                            
        self.so.name += "_3D_supg" 
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               pList,
                                               nList,
                                               sList,
                                               opts)
        ns.calculateSolution('3D_supg')
        # COMPARE VS SAVED FILES #
        expected_path = 'comparison_files/risingBubble_3D_supg.h5'
        expected = tables.open_file(os.path.join(self._scriptdir,expected_path))
        actual = tables.open_file('risingBubble_3D_supg.h5','r')
        assert np.allclose(expected.root.phi_t2,actual.root.phi_t2,atol=1e-10)
        assert np.allclose(expected.root.p_t2,actual.root.p_t2,atol=1e-10)
        assert np.allclose(expected.root.velocity_t2,actual.root.velocity_t2,atol=1e-10)        
        expected.close()
        actual.close()

    def test_3D_with_EV(self):
        # Set parameters for this test
        parameters.ct.USE_SUPG_NS=0
        parameters.ct.ARTIFICIAL_VISCOSITY_NS=2
        parameters.ct.nd=3
        # RELOAD MODULES
        self.reload_modules()
        pnList = [(vof_p,               vof_n),
                  (ls_p,                ls_n),
                  (redist_p,            redist_n),
                  (ls_consrv_p,         ls_consrv_n),
                  (twp_navier_stokes_p, twp_navier_stokes_n),
                  (pressureincrement_p, pressureincrement_n),
                  (pressure_p,          pressure_n),
                  (pressureInitial_p,   pressureInitial_n)]
        self.so = risingBubble_so
        pList=[]
        nList=[]
        sList=[]
        for (pModule,nModule) in pnList:
            pList.append(pModule)
            if pList[-1].name == None:
                pList[-1].name = pModule.__name__
            nList.append(nModule)
        for i in range(len(pnList)):            
            sList.append(default_s)                            
        self.so.name += "_3D_ev" 
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               pList,
                                               nList,
                                               sList,
                                               opts)
        ns.calculateSolution('3D_ev')
        # COMPARE VS SAVED FILES #
        expected_path = 'comparison_files/risingBubble_3D_ev.h5'
        expected = tables.open_file(os.path.join(self._scriptdir,expected_path))
        actual = tables.open_file('risingBubble_3D_ev.h5','r')
        assert np.allclose(expected.root.phi_t2,actual.root.phi_t2,atol=1e-10)
        assert np.allclose(expected.root.p_t2,actual.root.p_t2,atol=1e-10)
        assert np.allclose(expected.root.velocity_t2,actual.root.velocity_t2,atol=1e-10)        
        expected.close()
        actual.close()                        


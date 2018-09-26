#!/usr/bin/env python
"""
Test module for surface tension
"""
from builtins import range
from builtins import object
from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
Profiling.logLevel=1
Profiling.verbose=True
import os
import numpy as np
import tables
from . import (parameters,
               NS_convergence_so, NS_convergence,
               twp_navier_stokes_p,
               twp_navier_stokes_n,
               pressureincrement_p,
               pressureincrement_n,
               pressure_p,
               pressure_n,
               pressureInitial_p,
               pressureInitial_n)


class TestSurfaceTension(object):

    @classmethod
    def setup_class(cls):
        pass
    
    @classmethod
    def teardown_class(cls):
        pass
    
    def reload_modules(self):
        reload(NS_convergence)
        reload(NS_convergence_so)
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
    
    def test_no_stab(self):
        # Set parameters for this test
        parameters.ct.USE_SUPG_NS=0
        parameters.ct.ARTIFICIAL_VISCOSITY_NS=0
        # RELOAD MODULES
        self.reload_modules()
        pnList = [(twp_navier_stokes_p, twp_navier_stokes_n),
                  (pressureincrement_p, pressureincrement_n),
                  (pressure_p,          pressure_n),
                  (pressureInitial_p,   pressureInitial_n)]
        self.so = NS_convergence_so
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
        self.so.name += "_no_stab"
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               pList,
                                               nList,
                                               sList,
                                               opts)
        ns.calculateSolution('no_stab')
        # COMPARE VS REFERENCE #
        actual = tables.open_file('NS_convergence_no_stab.h5','r')
        assert np.isclose(np.amax(actual.root.p_t11) - np.amin(actual.root.p_t11),
                          0.57387646058,
                          atol=1e-10)
        assert np.isclose(np.amax(actual.root.u_t11) - np.amin(actual.root.u_t11),
                          1.68689338205,
                          atol=1e-10)
        assert np.isclose(np.amax(actual.root.v_t11) - np.amin(actual.root.v_t11),
                          1.68689335359,
                          atol=1e-10)        
        actual.close()
        
    def test_supg_with_shock_capturing(self):
        # Set parameters for this test
        parameters.ct.USE_SUPG_NS=1
        parameters.ct.ARTIFICIAL_VISCOSITY_NS=1
        # RELOAD MODULES
        self.reload_modules()
        pnList = [(twp_navier_stokes_p, twp_navier_stokes_n),
                  (pressureincrement_p, pressureincrement_n),
                  (pressure_p,          pressure_n),
                  (pressureInitial_p,   pressureInitial_n)]
        self.so = NS_convergence_so
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
        self.so.name += "_supg"
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               pList,
                                               nList,
                                               sList,
                                               opts)
        ns.calculateSolution('supg')
        # COMPARE VS REFERENCE #
        actual = tables.open_file('NS_convergence_supg.h5','r')
        assert np.isclose(np.amax(actual.root.p_t11) - np.amin(actual.root.p_t11),
                          0.691428927609,
                          atol=1e-10)
        assert np.isclose(np.amax(actual.root.u_t11) - np.amin(actual.root.u_t11),
                          1.68689322528,
                          atol=1e-10)
        assert np.isclose(np.amax(actual.root.v_t11) - np.amin(actual.root.v_t11),
                          1.68689322528,
                          atol=1e-10)        
        actual.close()
        
    def test_ev(self):
        # Set parameters for this test
        parameters.ct.USE_SUPG_NS=0
        parameters.ct.ARTIFICIAL_VISCOSITY_NS=2
        # RELOAD MODULES
        self.reload_modules()
        pnList = [(twp_navier_stokes_p, twp_navier_stokes_n),
                  (pressureincrement_p, pressureincrement_n),
                  (pressure_p,          pressure_n),
                  (pressureInitial_p,   pressureInitial_n)]
        self.so = NS_convergence_so
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
        self.so.name += "_ev" 
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               pList,
                                               nList,
                                               sList,
                                               opts)
        ns.calculateSolution('ev')
        # COMPARE VS REFERENCE #
        actual = tables.open_file('NS_convergence_ev.h5','r')
        assert np.isclose(np.amax(actual.root.p_t11) - np.amin(actual.root.p_t11),
                          0.594072307827,
                          atol=1e-10)
        assert np.isclose(np.amax(actual.root.u_t11) - np.amin(actual.root.u_t11),
                          1.68689325714,
                          atol=1e-10)
        assert np.isclose(np.amax(actual.root.v_t11) - np.amin(actual.root.v_t11),
                          1.68689318939,
                          atol=1e-10)        
        actual.close()



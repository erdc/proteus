#!/usr/bin/env python
"""
Test module for Elliptic Re-distancing
"""
from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
Profiling.logLevel=1
Profiling.verbose=True
import os
import numpy as np
import tables
from . import (vortex2D, vortex2D_so,
               ncls_p, ncls_n,
               rdls_p, rdls_n)

class TestEllipticRedistancing():

    @classmethod
    def setup_class(cls):
        pass
    
    @classmethod
    def teardown_class(cls):
        pass
    
    def setup_method(self,method):
        """Initialize the test problem. """
        self._scriptdir = os.path.dirname(__file__)                
               
    def teardown_method(self,method):
        pass
    
    def test_ELLIPTIC_REDISTANCING_0(self):
        # Set parameters for test #
        vortex2D.ct.ELLIPTIC_REDISTANCING = False
        reload(rdls_p)
        reload(vortex2D_so)
        pnList = [(ncls_p,ncls_n),
                  (rdls_p,rdls_n)]
        self.so = vortex2D_so
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
        self.so.name += "_ELLIPTIC_REDIST_0"
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               pList,
                                               nList,
                                               sList,
                                               opts)
        ns.calculateSolution('rdls')
        actual = tables.open_file('vortex_c0p1_level_1_ELLIPTIC_REDIST_0.h5','r')
        assert np.isclose(np.amax(actual.root.u_t2),0.137239911295,atol=1e-10)
        actual.close()
        
    def test_ELLIPTIC_REDISTANCING_1(self):
        # Set parameters for test #
        vortex2D.ct.ELLIPTIC_REDISTANCING = True
        vortex2D.ct.ELLIPTIC_REDISTANCING_TYPE = 1
        reload(rdls_p)
        reload(vortex2D_so)
        pnList = [(ncls_p,ncls_n),
                  (rdls_p,rdls_n)]
        self.so = vortex2D_so
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
        self.so.name += "_ELLIPTIC_REDIST_1"
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               pList,
                                               nList,
                                               sList,
                                               opts)
        ns.calculateSolution('rdls')
        actual = tables.open_file('vortex_c0p1_level_1_ELLIPTIC_REDIST_1.h5','r')
        assert np.isclose(np.amax(actual.root.u_t2),0.144279865323,atol=1e-10)
        actual.close()
        
    def test_ELLIPTIC_REDISTANCING_2(self):
        # Set parameters for test #
        vortex2D.ct.ELLIPTIC_REDISTANCING = True
        vortex2D.ct.ELLIPTIC_REDISTANCING_TYPE = 2
        reload(rdls_p)
        reload(vortex2D_so)
        pnList = [(ncls_p,ncls_n),
                  (rdls_p,rdls_n)]
        self.so = vortex2D_so
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
        self.so.name += "_ELLIPTIC_REDIST_2"
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               pList,
                                               nList,
                                               sList,
                                               opts)
        ns.calculateSolution('rdls')
        actual = tables.open_file('vortex_c0p1_level_1_ELLIPTIC_REDIST_2.h5','r')

        print np.amax(actual.root.u_t2)
        assert np.isclose(np.amax(actual.root.u_t2),0.112863525457,atol=1e-10)
        actual.close()
        
    def test_ELLIPTIC_REDISTANCING_3(self):
        # Set parameters for test #
        vortex2D.ct.ELLIPTIC_REDISTANCING = True
        vortex2D.ct.ELLIPTIC_REDISTANCING_TYPE = 3
        reload(rdls_p)
        reload(vortex2D_so)
        pnList = [(ncls_p,ncls_n),
                  (rdls_p,rdls_n)]
        self.so = vortex2D_so
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
        self.so.name += "_ELLIPTIC_REDIST_3"
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               pList,
                                               nList,
                                               sList,
                                               opts)
        ns.calculateSolution('rdls')
        actual = tables.open_file('vortex_c0p1_level_1_ELLIPTIC_REDIST_3.h5','r')
        assert np.isclose(np.amax(actual.root.u_t2),0.112863525462,atol=1e-10)        
        actual.close()

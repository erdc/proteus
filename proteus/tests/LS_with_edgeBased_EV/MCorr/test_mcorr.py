#!/usr/bin/env python
"""
Test module for the conservative LS with EV
"""

from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
Profiling.logLevel=2
Profiling.verbose=True
import numpy as np
import tables
import thelper_cons_ls
import thelper_cons_ls_so
import thelper_vof_p
import thelper_vof_n
import thelper_ncls_p
import thelper_ncls_n
import thelper_redist_p
import thelper_redist_n
import thelper_MCorr_p
import thelper_MCorr_n

class TestMCorr():

    @classmethod
    def setup_class(cls):
        pass
    
    @classmethod
    def teardown_class(cls):
        pass
    
    def setup_method(self,method):        
        """Initialize the test problem. """
        reload(thelper_cons_ls)
        self._scriptdir = os.path.dirname(__file__)
        
    def teardown_method(self,method):
        pass

    def test_supg(self):
        thelper_cons_ls.ct.STABILIZATION_TYPE_ncls=0
        thelper_cons_ls.ct.DO_REDISTANCING=False
        thelper_cons_ls.ct.STABILIZATION_TYPE_vof=0
        reload(thelper_cons_ls_so)
        reload(thelper_ncls_p)
        reload(thelper_ncls_n)
        reload(thelper_redist_p)
        reload(thelper_redist_n)
        reload(thelper_vof_p)
        reload(thelper_vof_n)
        reload(thelper_MCorr_p)
        reload(thelper_MCorr_n)

        self.so = __import__("thelper_cons_ls_so")
        self.pList=[]
        self.nList=[]
        self.sList=[]
        for (pModule,nModule) in self.so.pnList:
            self.pList.append(__import__(pModule))
            if self.pList[-1].name == None:
                self.pList[-1].name = pModule
            self.nList.append(__import__(nModule))
        for i in range(len(self.so.pnList)):
            s = default_s
            self.sList.append(s)                               
        self.so.name += "_supg"
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               self.pList,
                                               self.nList,
                                               self.sList,
                                               opts)
        ns.calculateSolution('vof')
        # COMPARE VS SAVED FILES #
        expected_path = 'comparison_files/cons_ls_level_3_supg.h5'
        expected = tables.openFile(os.path.join(self._scriptdir,expected_path))
        actual = tables.openFile('cons_ls_level_3_supg.h5','r')
        assert np.allclose(expected.root.vof_t2,
                           actual.root.vof_t2,
                           atol=1e-10)
        expected.close()
        actual.close()

    def test_edge_based_EV(self):
        thelper_cons_ls.ct.STABILIZATION_TYPE_ncls=1
        thelper_cons_ls.ct.DO_REDISTANCING=True
        thelper_cons_ls.ct.STABILIZATION_TYPE_vof=1
        reload(thelper_cons_ls_so)
        reload(thelper_ncls_p)
        reload(thelper_ncls_n)
        reload(thelper_vof_p)
        reload(thelper_vof_n)
        reload(thelper_MCorr_p)
        reload(thelper_MCorr_n)

        self.so = __import__("thelper_cons_ls_so")
        self.pList=[]
        self.nList=[]
        self.sList=[]
        for (pModule,nModule) in self.so.pnList:
            self.pList.append(__import__(pModule))
            if self.pList[-1].name == None:
                self.pList[-1].name = pModule
            self.nList.append(__import__(nModule))
        for i in range(len(self.so.pnList)):
            s = default_s
            self.sList.append(s)
                       
        self.so.name += "_edge_based_EV"
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               self.pList,
                                               self.nList,
                                               self.sList,
                                               opts)
        ns.calculateSolution('vof')
        # COMPARE VS SAVED FILES #
        expected_path = 'comparison_files/cons_ls_level_3_edge_based_EV.h5'
        expected = tables.openFile(os.path.join(self._scriptdir,expected_path))
        actual = tables.openFile('cons_ls_level_3_edge_based_EV.h5','r')
        assert np.allclose(expected.root.vof_t2,
                           actual.root.vof_t2,
                           atol=1e-10)
        expected.close()
        actual.close()        
        

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
import sys
sys.path.append('import_files')
import cons_ls
import cons_ls_so
import vof_p
import vof_n
import ncls_p
import ncls_n
import redist_p
import redist_n
import MCorr_p
import MCorr_n

class TestMCorr():

    @classmethod
    def setup_class(cls):
        pass
    
    @classmethod
    def teardown_class(cls):
        pass
    
    def setup_method(self,method):        
        """Initialize the test problem. """
        reload(cons_ls)
        self._scriptdir = os.path.dirname(__file__)
        
    def teardown_method(self,method):
        pass

    def test_supg(self):
        cons_ls.ct.STABILIZATION_TYPE_ncls=0
        cons_ls.ct.DO_REDISTANCING=False
        cons_ls.ct.STABILIZATION_TYPE_vof=0
        reload(cons_ls_so)
        reload(ncls_p)
        reload(ncls_n)
        reload(redist_p)
        reload(redist_n)
        reload(vof_p)
        reload(vof_n)
        reload(MCorr_p)
        reload(MCorr_n)

        self.so = __import__("cons_ls_so")
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
        cons_ls.ct.STABILIZATION_TYPE_ncls=1
        cons_ls.ct.DO_REDISTANCING=True
        cons_ls.ct.STABILIZATION_TYPE_vof=1
        reload(cons_ls_so)
        reload(ncls_p)
        reload(ncls_n)
        reload(vof_p)
        reload(vof_n)
        reload(MCorr_p)
        reload(MCorr_n)

        self.so = __import__("cons_ls_so")
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
        

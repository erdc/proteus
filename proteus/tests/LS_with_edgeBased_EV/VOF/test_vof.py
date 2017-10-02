#!/usr/bin/env python
"""
Test module for VOF with EV
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
import vof
import vof_p
import vof_n

class TestVOF():

    @classmethod
    def setup_class(cls):
        pass
    
    @classmethod
    def teardown_class(cls):
        pass
    
    def setup_method(self,method):
        """Initialize the test problem. """
        reload(vof)
        self.pList = [vof_p]
        self.nList = [vof_n]        
        self.sList = [default_s]
        self.so = default_so
        self.so.tnList = self.nList[0].tnList
        self._scriptdir = os.path.dirname(__file__)                
        self.sim_names = []
        self.aux_names = []
               
    def teardown_method(self,method):
        pass
    
    def test_vof(self):
        ########
        # SUPG #
        ########
        vof.ct.STABILIZATION_TYPE = 0 # SUPG
        vof.ct.FCT = False
        reload(vof_n)
        reload(vof_p)
        self.so.name = self.pList[0].name+"_SUPG"        
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               self.pList,
                                               self.nList,
                                               self.sList,
                                               opts)
        self.sim_names.append(ns.modelList[0].name)
        ns.calculateSolution('vof')
        # COMPARE VS SAVED FILES #
        expected_path = 'comparison_files/vof_level_3_SUPG.h5'
        expected = tables.openFile(os.path.join(self._scriptdir,expected_path))
        actual = tables.openFile('vof_level_3_SUPG.h5','r')
        assert np.allclose(expected.root.u_t2,
                           actual.root.u_t2,
                           atol=1e-10)
        expected.close()
        actual.close()
        
    def test_EV1(self):
        #######################
        # ENTROPY VISCOSITY 1 # Polynomial entropy
        #######################
        vof.ct.STABILIZATION_TYPE = 1 # EV
        vof.ct.ENTROPY_TYPE = 1 #polynomial
        vof.ct.cE = 1.0
        vof.ct.FCT = True
        reload(vof_n)
        reload(vof_p)
        self.so.name = self.pList[0].name+"_EV1"        
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               self.pList,
                                               self.nList,
                                               self.sList,
                                               opts)
        self.sim_names.append(ns.modelList[0].name)
        ns.calculateSolution('vof')
        # COMPARE VS SAVED FILES #
        expected_path = 'comparison_files/vof_level_3_EV1.h5'
        expected = tables.openFile(os.path.join(self._scriptdir,expected_path))
        actual = tables.openFile('vof_level_3_EV1.h5','r')
        assert np.allclose(expected.root.u_t2,
                           actual.root.u_t2,
                           atol=1e-10)
        expected.close()
        actual.close()
        
    def test_EV2(self): 
        vof.ct.STABILIZATION_TYPE = 1 # EV
        vof.ct.ENTROPY_TYPE = 2 #logarithmic
        vof.ct.cE = 0.1
        vof.ct.FCT = True
        reload(vof_n)
        reload(vof_p)
        self.so.name = self.pList[0].name+"_EV2"
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               self.pList,
                                               self.nList,
                                               self.sList,
                                               opts)
        self.sim_names.append(ns.modelList[0].name)
        ns.calculateSolution('vof')
        # COMPARE VS SAVED FILES #
        expected_path = 'comparison_files/vof_level_3_EV2.h5'
        expected = tables.openFile(os.path.join(self._scriptdir,expected_path))
        actual = tables.openFile('vof_level_3_EV2.h5','r')
        assert np.allclose(expected.root.u_t2,
                           actual.root.u_t2,
                           atol=1e-10)
        expected.close()
        actual.close()
        
    def test_SmoothnessBased(self):
        vof.ct.STABILIZATION_TYPE = 2 # EV
        vof.ct.FCT = True
        reload(vof_n)
        reload(vof_p)
        self.so.name = self.pList[0].name+"_SmoothnessBased"
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               self.pList,
                                               self.nList,
                                               self.sList,
                                               opts)
        self.sim_names.append(ns.modelList[0].name)
        ns.calculateSolution('vof')
        # COMPARE VS SAVED FILES #
        expected_path = 'comparison_files/vof_level_3_SmoothnessBased.h5'
        expected = tables.openFile(os.path.join(self._scriptdir,expected_path))
        actual = tables.openFile('vof_level_3_SmoothnessBased.h5','r')
        assert np.allclose(expected.root.u_t2,
                           actual.root.u_t2,
                           atol=1e-10)
        expected.close()
        actual.close()

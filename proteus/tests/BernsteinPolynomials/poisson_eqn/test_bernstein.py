#!/usr/bin/env python
"""
Test module for Berstein basis FE
"""
from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
Profiling.logLevel=1
Profiling.verbose=True
import os
import numpy as np
import tables
from . import (poisson_p,
               poisson_n)

class TestVOF():

    @classmethod
    def setup_class(cls):
        pass
    
    @classmethod
    def teardown_class(cls):
        pass
    
    def setup_method(self,method):
        """Initialize the test problem. """
        reload(poisson_p)
        self.pList = [poisson_p]
        self.nList = [poisson_n]        
        self.sList = [default_s]
        self.so = default_so
        #self.so.tnList = self.nList[0].tnList
        self._scriptdir = os.path.dirname(__file__)                
        self.sim_names = []
        self.aux_names = []
               
    def teardown_method(self,method):
        pass
    
    def test_2D_hex(self):
        poisson_p.ct.nd = 2
        poisson_p.ct.useHex = True

        reload(poisson_p)
        reload(poisson_n)
        self.so.name = "2D_"+self.pList[0].name+"_hex_degree2"
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               self.pList,
                                               self.nList,
                                               self.sList,
                                               opts)
        self.sim_names.append(ns.modelList[0].name)
        ns.calculateSolution('poisson')
        
        # COMPARE VS SAVED FILES #
        #expected_path = 'comparison_files/vof_level_3_SUPG.h5'
        #expected = tables.open_file(os.path.join(self._scriptdir,expected_path))
        #actual = tables.open_file('vof_level_3_SUPG.h5','r')
        #assert np.allclose(expected.root.u_t2,
        #                   actual.root.u_t2,
        #                   atol=1e-10)
        #expected.close()
        #actual.close()

    def test_2D_simplex(self):
        poisson_p.ct.nd = 2
        poisson_p.ct.useHex = False

        reload(poisson_p)
        reload(poisson_n)
        self.so.name = "2D_"+self.pList[0].name+"_simplex_degree2"
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               self.pList,
                                               self.nList,
                                               self.sList,
                                               opts)
        self.sim_names.append(ns.modelList[0].name)
        ns.calculateSolution('poisson')
        
        # COMPARE VS SAVED FILES #
        #expected_path = 'comparison_files/vof_level_3_SUPG.h5'
        #expected = tables.open_file(os.path.join(self._scriptdir,expected_path))
        #actual = tables.open_file('vof_level_3_SUPG.h5','r')
        #assert np.allclose(expected.root.u_t2,
        #                   actual.root.u_t2,
        #                   atol=1e-10)
        #expected.close()
        #actual.close()        
        

    def test_3D_hex(self):
        poisson_p.ct.nd = 3
        poisson_p.ct.useHex = True

        reload(poisson_p)
        reload(poisson_n)
        self.so.name = "3D_"+self.pList[0].name+"_hex_degree2"
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               self.pList,
                                               self.nList,
                                               self.sList,
                                               opts)
        self.sim_names.append(ns.modelList[0].name)
        ns.calculateSolution('poisson')

    def test_3D_simplex(self):
        poisson_p.ct.nd = 3
        poisson_p.ct.useHex = False

        reload(poisson_p)
        reload(poisson_n)
        self.so.name = "3D_"+self.pList[0].name+"_simplex_degree2"
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               self.pList,
                                               self.nList,
                                               self.sList,
                                               opts)
        self.sim_names.append(ns.modelList[0].name)
        ns.calculateSolution('poisson')        

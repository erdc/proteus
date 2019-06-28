#!/usr/bin/env python
"""
Test module for SWEs
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
import pytest
from proteus import default_so


class TestSWEs(object):

    @classmethod
    def setup_class(cls):
        pass

    @classmethod
    def teardown_class(cls):
        pass

    def setup_method(self,method):
        self._scriptdir = os.path.dirname(__file__)

    def teardown_method(self,method):
        pass

    def run_test(self,sw_p,sw_n,name):
        pnList = [(sw_p, sw_n)]
        self.so = default_so
        self.so.tnList = sw_n.tnList
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
        self.so.name = name
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               pList,
                                               nList,
                                               sList,
                                               opts)
        ns.calculateSolution(name)

    def compare_files(self,path,name):
        # COMPARE VS SAVED FILES #
        expected_path = path+'/'+name+'.h5'
        expected = tables.open_file(os.path.join(self._scriptdir,expected_path))
        actual = tables.open_file(name+'.h5','r')
        assert np.allclose(expected.root.h_t2,actual.root.h_t2,atol=1e-10)
        assert np.allclose(expected.root.velocity_t2,actual.root.velocity_t2,atol=1e-10)
        expected.close()
        actual.close()

    def test_dam_over_bumps(self):
        from .dam_over_bumps import (sw_p,sw_n)
        name = "SWEs_dam_over_bumps"
        # run problem
        self.run_test(sw_p,sw_n,name)
        # compare output files
        self.compare_files('dam_over_bumps/comparison_files',name)

    def test_oneD_dambreak_flat_bottom(self):
        from .oneD_dambreak_flat_bottom import (sw_p,sw_n)
        name = "SWEs_oneD_dambreak_flat_bottom"
        # run problem
        self.run_test(sw_p,sw_n,name)
        # compare output files
        self.compare_files('oneD_dambreak_flat_bottom/comparison_files',name)

    def test_oneD_paraboloid_with_friction(self):
        from .paraboloid_with_friction.oneD import (sw_p,sw_n)
        name = "SWEs_oneD_paraboloid_with_friction"
        # run problem
        self.run_test(sw_p,sw_n,name)
        # compare output files
        self.compare_files('paraboloid_with_friction/oneD/comparison_files',name)

    def test_twoD_paraboloid_with_friction(self):
        from .paraboloid_with_friction.twoD import (sw_p,sw_n)
        name = "SWEs_twoD_paraboloid_with_friction"
        # run problem
        self.run_test(sw_p,sw_n,name)
        # compare output files
        self.compare_files('paraboloid_with_friction/twoD/comparison_files',name)

    def test_gauges(self):
        from .test_gauges import (sw_p,sw_n)
        name = "SWEs_test_gauges"
        # run problem
        self.run_test(sw_p,sw_n,name)
        # compare output files
        self.compare_files('test_gauges/comparison_files',name)

    def test_reflecting_BCs(self):
        from .test_reflecting_BCs import (sw_p,sw_n)
        name = "SWEs_test_reflecting_BCs"
        # run problem
        self.run_test(sw_p,sw_n,name)
        # compare output files
        self.compare_files('test_reflecting_BCs/comparison_files',name)

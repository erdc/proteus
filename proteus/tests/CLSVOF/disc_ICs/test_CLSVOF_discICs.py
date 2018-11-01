#!/usr/bin/env python
"""
Test module for testing disc ICs for CLSVOF
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
from . import (parameters,
               clsvof,
               clsvof_p,
               clsvof_n)

class TestSWEs(object):

    @classmethod
    def setup_class(cls):
        pass

    @classmethod
    def teardown_class(cls):
        pass

    def reload_modules(self):
        reload(default_so)
        reload(clsvof)
        reload(clsvof_p)
        reload(clsvof_n)
                        
    def setup_method(self,method):
        self._scriptdir = os.path.dirname(__file__)

    def teardown_method(self,method):
        pass

    def run_test(self,clsvof_p,clsvof_n,name):
        pnList = [(clsvof_p, clsvof_n)]
        self.so = default_so
        self.so.tnList = clsvof_n.tnList
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
        assert np.allclose(expected.root.u_t1,actual.root.u_t1,atol=1e-10)
        expected.close()
        actual.close()

    def test_case_1(self):
        from . import (clsvof_p,clsvof_n)
        name = "test_case_1"
        # run problem
        self.run_test(clsvof_p,clsvof_n,name)
        # compare output files
        self.compare_files('comparison_files',name)

    def test_case_2(self):
        parameters.ct.test_case=2
        self.reload_modules()        
        name = "test_case_2"
        # run problem
        self.run_test(clsvof_p,clsvof_n,name)
        # compare output files
        self.compare_files('comparison_files',name)        

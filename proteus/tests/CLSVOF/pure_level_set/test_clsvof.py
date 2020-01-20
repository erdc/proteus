#!/usr/bin/env python
"""
Test module for CLSVOF
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

class TestCLSVOF(object):

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

    def compare_files(self,path,name):
        # COMPARE VS SAVED FILES #
        #expected_path = path+'/'+name+'.h5'
        #expected = tables.open_file(os.path.join(self._scriptdir,expected_path))
        #actual = tables.open_file(name+'.h5','r')
        #assert np.allclose(expected.root.u_t1,actual.root.u_t1,atol=1e-10)
        #expected.close()

        actual = tables.open_file(name+'.h5','r')
        expected_path = 'comparison_files/' + 'comparison_' + name + '_u_t2.csv'
        #write comparison file
        #np.array(actual.root.u_t2).tofile(os.path.join(self._scriptdir, expected_path),sep=",")
        np.testing.assert_almost_equal(np.fromfile(os.path.join(self._scriptdir, expected_path),sep=","),np.array(actual.root.u_t2).flatten(),decimal=10)
        actual.close()

    def test_case_1(self):
        # Set parameters for this test
        parameters.ct.test_case=1
        # RELOAD MODULES
        self.reload_modules()
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
        self.so.name = "clsvof_test_case_1"
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               pList,
                                               nList,
                                               sList,
                                               opts)
        ns.calculateSolution('test_case_1')
        # COMPARE VS SAVED FILES #

        self.compare_files('comparison_files',self.so.name)
        #expected_path = 'comparison_files/clsvof_test_case_1.h5'
        #expected = tables.open_file(os.path.join(self._scriptdir,expected_path))
        #actual = tables.open_file('clsvof_test_case_1.h5','r')
        #assert np.allclose(expected.root.u_t2,actual.root.u_t2,atol=1e-10)
        #expected.close()
        #actual.close()

    def test_case_2(self):
        # Set parameters for this test
        parameters.ct.test_case=2
        # RELOAD MODULES
        self.reload_modules()
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
        self.so.name = "clsvof_test_case_2"
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               pList,
                                               nList,
                                               sList,
                                               opts)
        ns.calculateSolution('test_case_2')
        # COMPARE VS SAVED FILES #
        self.compare_files('comparison_files',self.so.name)
        #expected_path = 'comparison_files/clsvof_test_case_2.h5'
        #expected = tables.open_file(os.path.join(self._scriptdir,expected_path))
        #actual = tables.open_file('clsvof_test_case_2.h5','r')
        #assert np.allclose(expected.root.u_t2,actual.root.u_t2,atol=1e-10)
        #expected.close()
        #actual.close()

    def test_case_3(self):
        # Set parameters for this test
        parameters.ct.test_case=3
        # RELOAD MODULES
        self.reload_modules()
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
        self.so.name = "clsvof_test_case_3"
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               pList,
                                               nList,
                                               sList,
                                               opts)
        ns.calculateSolution('test_case_3')
        # COMPARE VS SAVED FILES #
        self.compare_files('comparison_files',self.so.name)
        #expected_path = 'comparison_files/clsvof_test_case_3.h5'
        #expected = tables.open_file(os.path.join(self._scriptdir,expected_path))
        #actual = tables.open_file('clsvof_test_case_3.h5','r')
        #assert np.allclose(expected.root.u_t2,actual.root.u_t2,atol=1e-10)
        #expected.close()
        #actual.close()

    def test_case_4(self):
        # Set parameters for this test
        parameters.ct.test_case=4
        # RELOAD MODULES
        self.reload_modules()
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
        self.so.name = "clsvof_test_case_4"
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               pList,
                                               nList,
                                               sList,
                                               opts)
        ns.calculateSolution('test_case_4')
        # COMPARE VS SAVED FILES #
        self.compare_files('comparison_files',self.so.name)
        #expected_path = 'comparison_files/clsvof_test_case_4.h5'
        #expected = tables.open_file(os.path.join(self._scriptdir,expected_path))
        #actual = tables.open_file('clsvof_test_case_4.h5','r')
        #assert np.allclose(expected.root.u_t2,actual.root.u_t2,atol=1e-10)
        #expected.close()
        #actual.close()

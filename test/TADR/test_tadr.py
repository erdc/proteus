#!/usr/bin/env python
"""
Test module for TADR with EV
"""
from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
Profiling.logLevel=2
Profiling.verbose=True
import numpy as np
import h5py
from . import thelper_tadr
from . import thelper_tadr_p
from . import thelper_tadr_n

class TestTADR(object):

    @classmethod
    def setup_class(cls):
        pass

    @classmethod
    def teardown_class(cls):
        pass

    def setup_method(self,method):
        """Initialize the test problem. """
        reload(default_n)
        reload(thelper_tadr)
        self.pList = [thelper_tadr_p]
        self.nList = [thelper_tadr_n]
        self.sList = [default_s]
        self.so = default_so
        self.so.tnList = self.nList[0].tnList
        self._scriptdir = os.path.dirname(__file__)
        self.sim_names = []
        self.aux_names = []

    def teardown_method(self,method):
        pass

    def test_supg(self):
        ########
        # SUPG #
        ########
        thelper_tadr.ct.STABILIZATION_TYPE = 0 # SUPG
        thelper_tadr.ct.FCT = False
        reload(default_n)
        reload(thelper_tadr_p)
        reload(thelper_tadr_n)
        tname = "_SUPG"
        self.so.name = self.pList[0].name+tname
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               self.pList,
                                               self.nList,
                                               self.sList,
                                               opts)
        self.sim_names.append(ns.modelList[0].name)
        ns.calculateSolution('tadr')
        # COMPARE VS SAVED FILES #
        
        with h5py.File("tadr_level_0"+tname+".h5",'r') as actual:
            expected_path = os.path.join(self._scriptdir,"comparison_files","tadr_level_0"+tname+ "_u_t2.csv")
            #write comparison file
            #(actual['u_t2'][:]).tofile(expected_path,sep=",")
            np.testing.assert_almost_equal(np.fromfile(expected_path,sep=","),actual['u_t2'][:],decimal=10)

    def test_TaylorGalerkin(self):
        ##################
        # TaylorGalerkin #
        ##################
        thelper_tadr.ct.STABILIZATION_TYPE = 1 # Taylor Galerkin
        thelper_tadr.ct.FCT = False
        reload(default_n)
        reload(thelper_tadr_p)
        reload(thelper_tadr_n)
        tname = "_TaylorGalerkin"
        self.so.name = self.pList[0].name+tname
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               self.pList,
                                               self.nList,
                                               self.sList,
                                               opts)
        self.sim_names.append(ns.modelList[0].name)
        ns.calculateSolution('tadr')
        # COMPARE VS SAVED FILES #
        with h5py.File("tadr_level_0"+tname+".h5",'r') as actual:
            expected_path = os.path.join(self._scriptdir,"comparison_files","tadr_level_0"+tname+ "_u_t2.csv")
            #write comparison file
            #(actual['u_t2'][:]).tofile(expected_path,sep=",")
            np.testing.assert_almost_equal(np.fromfile(expected_path,sep=","),actual['u_t2'][:],decimal=10)

    def test_EV1(self):
        #######################
        # ENTROPY VISCOSITY 1 # Polynomial entropy
        #######################
        thelper_tadr.ct.STABILIZATION_TYPE = 2 # EV
        thelper_tadr.ct.ENTROPY_TYPE = 1 #polynomial
        thelper_tadr.ct.cE = 1.0
        thelper_tadr.ct.FCT = True
        reload(default_n)
        reload(thelper_tadr_p)
        reload(thelper_tadr_n)
        tname = "_EV1"
        self.so.name = self.pList[0].name+tname
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               self.pList,
                                               self.nList,
                                               self.sList,
                                               opts)
        self.sim_names.append(ns.modelList[0].name)
        ns.calculateSolution('tadr')
        # COMPARE VS SAVED FILES #
        with h5py.File("tadr_level_0"+tname+".h5",'r') as actual:
            expected_path = os.path.join(self._scriptdir,"comparison_files", "tadr_level_0"+tname+ "_u_t2.csv")
            #write comparison file
            #(actual['u_t2'][:]).tofile(expected_path,sep=",")
            np.testing.assert_almost_equal(np.fromfile(expected_path,sep=","),actual['u_t2'][:],decimal=10)

    def test_EV2(self):
        thelper_tadr.ct.STABILIZATION_TYPE = 2 # EV
        thelper_tadr.ct.ENTROPY_TYPE = 1 #logarithmic
        thelper_tadr.ct.cE = 0.1
        thelper_tadr.ct.FCT = True
        reload(default_n)
        reload(thelper_tadr_p)
        reload(thelper_tadr_n)
        tname = "_EV2"
        self.so.name = self.pList[0].name+tname
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               self.pList,
                                               self.nList,
                                               self.sList,
                                               opts)
        self.sim_names.append(ns.modelList[0].name)
        ns.calculateSolution('tadr')
        # COMPARE VS SAVED FILES #
        with h5py.File("tadr_level_0"+tname+".h5",'r') as actual:
            expected_path = os.path.join(self._scriptdir,"comparison_files","tadr_level_0"+tname+ "_u_t2.csv")
            #write comparison file
            #(actual['u_t2'][:]).tofile(expected_path,sep=",")
            np.testing.assert_almost_equal(np.fromfile(expected_path,sep=","),actual['u_t2'][:],decimal=10)

    def test_SmoothnessBased(self):
        thelper_tadr.ct.STABILIZATION_TYPE = 3 # Smoothness based
        thelper_tadr.ct.FCT = True
        reload(default_n)
        reload(thelper_tadr_p)
        reload(thelper_tadr_n)
        tname = "_SmoothnessBased"
        self.so.name = self.pList[0].name+tname
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               self.pList,
                                               self.nList,
                                               self.sList,
                                               opts)
        self.sim_names.append(ns.modelList[0].name)
        ns.calculateSolution('tadr')
        # COMPARE VS SAVED FILES #
        with h5py.File("tadr_level_0"+tname+".h5",'r') as actual:
            expected_path = os.path.join(self._scriptdir,"comparison_files","tadr_level_0"+tname+ "_u_t2.csv")
            #write comparison file
            #(actual['u_t2'][:]).tofile(expected_path,sep=",")
            np.testing.assert_almost_equal(np.fromfile(expected_path,sep=","),actual['u_t2'][:],decimal=10)

    def test_stab4(self):
        thelper_tadr.ct.STABILIZATION_TYPE = 4 # Proposed by D.Kuzmin
        thelper_tadr.ct.FCT = True
        reload(default_n)
        reload(thelper_tadr_p)
        reload(thelper_tadr_n)
        tname="_stab4"
        self.so.name = self.pList[0].name+tname
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               self.pList,
                                               self.nList,
                                               self.sList,
                                               opts)
        self.sim_names.append(ns.modelList[0].name)
        ns.calculateSolution('tadr')
        # COMPARE VS SAVED FILES #
        with h5py.File("tadr_level_0"+tname+".h5",'r') as actual:
            expected_path = os.path.join(self._scriptdir,"comparison_files","tadr_level_0"+tname+ "_u_t2.csv")
            #write comparison file
            #(actual['u_t2'][:]).tofile(expected_path,sep=",")
            np.testing.assert_almost_equal(np.fromfile(expected_path,sep=","),actual['u_t2'][:],decimal=10)

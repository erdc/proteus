#!/usr/bin/env python
"""
Test module for CLSVOF with RANS3PF
"""
from builtins import range
from builtins import object
from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
Profiling.logLevel=7
Profiling.verbose=True
import numpy as np
import tables
import pytest
from proteus import default_so
from . import (parameters,
               multiphase,
               multiphase_so,
               clsvof_p,
               clsvof_n,
               twp_navier_stokes_p,
               twp_navier_stokes_n,
               pressureincrement_p,
               pressureincrement_n,
               pressure_p,
               pressure_n,
               pressureInitial_p,
               pressureInitial_n)

class TestCLSVOF_with_RANS3PF(object):

    @classmethod
    def setup_class(cls):
        pass

    @classmethod
    def teardown_class(cls):
        pass

    def reload_modules(self):
        reload(default_so)
        reload(multiphase)
        reload(multiphase_so)
        reload(clsvof_p)
        reload(clsvof_n)
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
        FileList = ['multiphase_2D_falling_bubble.h5','multiphase_2D_falling_bubble.xmf',
                    'multiphase_3D_falling_bubble.h5','multiphase_3D_falling_bubble.xmf',
                    ]
        for file in FileList:
            if os.path.isfile(file):
                os.remove(file)
            else:
                pass

    def test_2D_falling_bubble(self):
        # Set parameters for this test
        parameters.ct.test_case=1
        # RELOAD MODULES
        self.reload_modules()
        pnList = [(clsvof_p,               clsvof_n),
                  (twp_navier_stokes_p, twp_navier_stokes_n),
                  (pressureincrement_p, pressureincrement_n),
                  (pressure_p,          pressure_n),
                  (pressureInitial_p,   pressureInitial_n)]
        self.so = multiphase_so
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
        self.so.name += "_2D_falling_bubble"
        # NUMERICAL SOLUTION #
        try:
            ns = proteus.NumericalSolution.NS_base(self.so,
                                               pList,
                                               nList,
                                               sList,
                                               opts)
        except:
            assert 0, "Numerical setup failed"
        try:
            ns.calculateSolution('2D_falling_bubble')
        except:
            assert 0, "Calculate solution failed" 

        # COMPARE VS SAVED FILES #
        actual = tables.open_file('multiphase_2D_falling_bubble.h5','r')
        expected_path = 'comparison_files/' + 'comparison_2D_phi_t2.csv'
        #write comparison file
        #np.array(actual.root.phi_t2).tofile(os.path.join(self._scriptdir, expected_path),sep=",")
        np.testing.assert_almost_equal(np.fromfile(os.path.join(self._scriptdir, expected_path),sep=","),np.array(actual.root.phi_t2),decimal=10)
        actual.close()

    def test_3D_falling_bubble(self):
        # Set parameters for this test
        parameters.ct.test_case=2
        # RELOAD MODULES
        self.reload_modules()
        pnList = [(clsvof_p,               clsvof_n),
                  (twp_navier_stokes_p, twp_navier_stokes_n),
                  (pressureincrement_p, pressureincrement_n),
                  (pressure_p,          pressure_n),
                  (pressureInitial_p,   pressureInitial_n)]
        self.so = multiphase_so
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
        self.so.name += "_3D_falling_bubble"
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               pList,
                                               nList,
                                               sList,
                                               opts)
        ns.calculateSolution('3D_falling_bubble')
        # COMPARE VS SAVED FILES #
        actual = tables.open_file('multiphase_3D_falling_bubble.h5','r')
        expected_path = 'comparison_files/' + 'comparison_3D_phi_t2.csv'
        #write comparison file
        #np.array(actual.root.phi_t2).tofile(os.path.join(self._scriptdir, expected_path),sep=",")
        np.testing.assert_almost_equal(np.fromfile(os.path.join(self._scriptdir, expected_path),sep=","),np.array(actual.root.phi_t2),decimal=10)

        actual.close()        


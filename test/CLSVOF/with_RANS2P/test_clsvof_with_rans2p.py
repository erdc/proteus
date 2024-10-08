#!/usr/bin/env python
"""
Test module for clsvof with rans2p
"""
from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
Profiling.logLevel=1
Profiling.verbose=True
import os
import numpy as np
import h5py
import pytest
from proteus import default_so
from . import (multiphase_so, multiphase,
               clsvof_p,
               clsvof_n,
               twp_navier_stokes_p,
               twp_navier_stokes_n)

class TestCLSVOFWithRans2p(object):

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
        
    def setup_method(self,method):
        self._scriptdir = os.path.dirname(__file__)
        
    def teardown_method(self,method):
        FileList = ['mesh.ele',
                    'mesh.edge',
                    'mesh.node',
                    'mesh.neigh',
                    'mesh.face',
                    'mesh.poly',
                    ]
        for file in FileList:
            if os.path.isfile(file):
                os.remove(file)
            else:
                pass

    @pytest.mark.skip(reason="Not reproducible on both python2 and python3")
    def test_2D_multiphase(self):
        # RELOAD MODULES
        self.reload_modules()
        pnList = [(twp_navier_stokes_p, twp_navier_stokes_n),
                  (clsvof_p,               clsvof_n)]
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
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               pList,
                                               nList,
                                               sList,
                                               opts)
        ns.calculateSolution('2D_falling_bubble')
        # COMPARE VS SAVED FILES #
        expected_path = 'comparison_files/multiphase_2D_falling_bubble.h5'
        expected= h5py.File(os.path.join(self._scriptdir,expected_path))
        actual = h5py.File('multiphase_2D_falling_bubble.h5','r')
        assert np.allclose(expected['phi_t2'],actual['phi_t2'],atol=1e-10)       
        expected.close()
        actual.close()

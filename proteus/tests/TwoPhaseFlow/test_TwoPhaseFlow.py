#!/usr/bin/env python
"""
Test module for TwoPhaseFlow
"""
import pytest
import tables
import numpy as np
import proteus.defaults
from proteus import Context
from proteus import default_so
from proteus.iproteus import *
import os
import sys
Profiling.logLevel=1
Profiling.verbose=True
from proteus.tests import Norms


L2_norm_risingBubble_phi = 0.9350911663163574 
L2_norm_risingBubble_velocity = np.array([0.02598942, 0.0315501,  0.        ])
L2_norm_marin_phi = 1.8423547464722558 
L2_norm_marin_velocity=np.array([0.45341437, 0.03755329, 0.82602671])
L2_norm_moses_phi = 0.8269543135413808 
L2_norm_moses_velocity = np.array([0.72064355, 0.12947117, 0.99516898])
L2_norm_dambreak_phi=2.797989916796235
L2_norm_dambreak_velocity=np.array([0.48407022, 0.95495096, 0.        ])
L2_norm_dambreakSolverOptions_phi = 2.7959767243349942 
L2_norm_dambreakSolverOptions_velocity = np.array([0.7568211,  0.70052942, 0.        ])
L2_norm_buckling_phi=0.5681116923610815
L2_norm_buckling_velocity=np.array( [0.19366006, 0.38872359, 0.        ])
L2_norm_filling_phi=0.109810825700778
L2_norm_filling_velocity = np.array([0.020761056644454,0.016814557793586,0.000000000000000])


class TestTwoPhaseFlow(object):

    def teardown_method(self, method):
        """ Tear down function """
        FileList = ['marin.h5',
                    'moses.h5',
                    'damBreak.h5',
                    'damBreak_solver_options.h5',
                    'TwoDimBucklingFlow.h5',
                    'filling.h5'
                    ]
        for file in FileList:
            if os.path.isfile(file):
                os.remove(file)
            else:
                pass


    def setup_method(self,method):
        self._scriptdir = os.path.dirname(__file__)
        self.path = proteus.__path__[0]+"/tests/TwoPhaseFlow/"

    def compare_vs_saved_files(self,name,comparePhi,compareVelocity):
        actual = tables.open_file(name+'.h5','r')
        L2_norm_phi = Norms.get_L2_norm(actual,actual.root.phi_t2)
        L2_norm_velocity = Norms.get_L2_vectorNorm(actual,actual.root.velocity_t2)
        actual.close()
        #print("%.15f,%.15f,%.15f,%.15f" % (L2_norm_phi,L2_norm_velocity[0],L2_norm_velocity[1],L2_norm_velocity[2]))
        #import sys; sys.exit()
        np.testing.assert_almost_equal(L2_norm_phi,comparePhi)
        np.testing.assert_almost_equal(L2_norm_velocity,compareVelocity)

    # *** 2D tests *** #
    def test_risingBubble(self): #uses structured triangle mesh
        os.system("parun --TwoPhaseFlow --path " + self.path + " "
                  "risingBubble.py -l5 -v -C 'final_time=0.1 dt_output=0.1 refinement=1'")
        self.compare_vs_saved_files("risingBubble",L2_norm_risingBubble_phi,L2_norm_risingBubble_velocity)
        self.teardown_method(self)

    def test_damBreak(self):
        os.system("parun --TwoPhaseFlow --path " + self.path + " "
                  "damBreak.py -l5 -v -C 'final_time=0.1 dt_output=0.1 he=0.1'")
        self.compare_vs_saved_files("damBreak",L2_norm_dambreak_phi,L2_norm_dambreak_velocity)
        self.teardown_method(self)

    def test_damBreak_solver_options(self):
        os.system("parun --TwoPhaseFlow --path " + self.path + " "
                  "damBreak_solver_options.py -l5 -v -C 'final_time=0.1 dt_output=0.1 he=0.5'")
        self.compare_vs_saved_files("damBreak_solver_options",L2_norm_dambreakSolverOptions_phi,L2_norm_dambreakSolverOptions_velocity)
        self.teardown_method(self)

#    @pytest.mark.skip(reason="long test")
    def test_TwoDimBucklingFlow(self):
        os.system("parun --TwoPhaseFlow --path " + self.path + " "
                  "TwoDimBucklingFlow.py -l5 -v -C 'final_time=0.1 dt_output=0.1 he=0.09'")
        self.compare_vs_saved_files("TwoDimBucklingFlow",L2_norm_buckling_phi,L2_norm_buckling_velocity)
        self.teardown_method(self)

#    @pytest.mark.skip(reason="long test")
    def test_fillingTank(self):
        os.system("parun --TwoPhaseFlow --path " + self.path + " "
                  "fillingTank.py -l5 -v -C 'final_time=0.02 dt_output=0.02 he=0.01'")
        self.compare_vs_saved_files("fillingTank",L2_norm_filling_phi,L2_norm_filling_velocity)
        self.teardown_method(self)

    # *** 3D tests *** #
    def test_marin(self):
        os.system("parun --TwoPhaseFlow --path " + self.path + " "
                  "marin.py -l5 -v -C 'final_time=0.1 dt_output=0.1 he=0.5'")
        self.compare_vs_saved_files("marin",L2_norm_marin_phi,L2_norm_marin_velocity)
        self.teardown_method(self)

    def test_moses(self):
        os.system("parun --TwoPhaseFlow --path " + self.path + " "
                  "moses.py -l5 -v -C 'final_time=0.1 dt_output=0.1 he=0.5'")
        self.compare_vs_saved_files("moses",L2_norm_moses_phi,L2_norm_moses_velocity)
        self.teardown_method(self)

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

class TestTwoPhaseFlow(object):

    def teardown_method(self, method):
        """ Tear down function """
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


    def setup_method(self,method):
        self._scriptdir = os.path.dirname(__file__)
        self.path = proteus.__path__[0]+"/tests/TwoPhaseFlow/"

    def compare_vs_saved_files(self,name):
        expected_path = 'comparison_files/' + name + '.h5'
        expected = tables.open_file(os.path.join(self._scriptdir,expected_path))
        actual = tables.open_file(name+'.h5','r')
        assert np.allclose(expected.root.phi_t2,actual.root.phi_t2,atol=1e-8)
        assert np.allclose(expected.root.velocity_t2,actual.root.velocity_t2,atol=1e-8)
        expected.close()
        actual.close()

    # *** 2D tests *** #
    def test_risingBubble(self):
        os.system("parun --TwoPhaseFlow --path " + self.path + " "
                  "risingBubble.py -l5 -v -C 'final_time=0.1 dt_output=0.1 refinement=1'")
        self.compare_vs_saved_files("risingBubble")

    def test_damBreak(self):
        os.system("parun --TwoPhaseFlow --path " + self.path + " "
                  "damBreak.py -l5 -v -C 'final_time=0.1 dt_output=0.1 he=0.1'")
        self.compare_vs_saved_files("damBreak")

    def test_damBreak_solver_options(self):
        os.system("parun --TwoPhaseFlow --path " + self.path + " "
                  "damBreak_solver_options.py -l5 -v -C 'final_time=0.1 dt_output=0.1 he=0.05'")
        self.compare_vs_saved_files("damBreak_solver_options")

#    @pytest.mark.skip(reason="long test")
    def test_TwoDimBucklingFlow(self):
        os.system("parun --TwoPhaseFlow --path " + self.path + " "
                  "TwoDimBucklingFlow.py -l5 -v -C 'final_time=0.1 dt_output=0.1 he=0.05'")
        self.compare_vs_saved_files("TwoDimBucklingFlow")

#    @pytest.mark.skip(reason="long test")
    def test_fillingTank(self):
        os.system("parun --TwoPhaseFlow --path " + self.path + " "
                  "fillingTank.py -l5 -v -C 'final_time=0.1 dt_output=0.1 he=0.01'")
        self.compare_vs_saved_files("fillingTank")

    # *** 3D tests *** #
    def test_marin(self):
        os.system("parun --TwoPhaseFlow --path " + self.path + " "
                  "marin.py -l5 -v -C 'final_time=0.1 dt_output=0.1 he=0.5'")
        self.compare_vs_saved_files("marin")
        self.teardown_method(self)

    def test_moses(self):
        os.system("parun --TwoPhaseFlow --path " + self.path + " "
                  "moses.py -l5 -v -C 'final_time=0.1 dt_output=0.1 he=0.5'")
        self.compare_vs_saved_files("moses")
        self.teardown_method(self)

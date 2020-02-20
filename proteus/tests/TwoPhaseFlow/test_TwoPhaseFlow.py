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

    def setup_method(self,method):
        self._scriptdir = os.path.dirname(__file__)
        self.path = proteus.__path__[0]+"/tests/TwoPhaseFlow/"

    def teardown_method(self, method):
        """ Tear down function """
        FileList = ['marin.h5','marin.xmf'
                    'moses.h5','moses.xmf'
                    'damBreak.h5','damBreak.xmf'
                    'damBreak_solver_options.h5','damBreak_solver_options.xmf'
                    'TwoDimBucklingFlow.h5','TwoDimBucklingFlow.xmf'
                    'filling.h5','filling.xmf'
                    ]
        for file in FileList:
            if os.path.isfile(file):
                os.remove(file)
            else:
                pass

    def compare_vs_saved_files(self,name):
        actual = tables.open_file(name+'.h5','r')

        expected_path = 'comparison_files/' + 'comparison_' + name + '_phi_t2.csv'
        #write comparison file
        #np.array(actual.root.phi_t2).tofile(os.path.join(self._scriptdir, expected_path),sep=",")
        np.testing.assert_almost_equal(np.fromfile(os.path.join(self._scriptdir, expected_path),sep=","),np.array(actual.root.phi_t2).flatten(),decimal=10)

        expected_path = 'comparison_files/' + 'comparison_' + name + '_velocity_t2.csv'
        #write comparison file
        #np.array(actual.root.velocity_t2).tofile(os.path.join(self._scriptdir, expected_path),sep=",")
        np.testing.assert_almost_equal(np.fromfile(os.path.join(self._scriptdir, expected_path),sep=","),np.array(actual.root.velocity_t2).flatten(),decimal=10)

        actual.close()

    # *** 2D tests *** #
    def test_risingBubble(self): #uses structured triangle mesh
        os.system("parun --TwoPhaseFlow --path " + self.path + " "
                  "risingBubble.py -l5 -v -C 'final_time=0.1 dt_output=0.1 refinement=1'")
        self.compare_vs_saved_files("risingBubble")

    def test_damBreak(self):
        os.system("parun --TwoPhaseFlow --path " + self.path + " "
                  "damBreak.py -l5 -v -C 'final_time=0.1 dt_output=0.1 he=0.1'")
        self.compare_vs_saved_files("damBreak")

    @pytest.mark.skip(reason="numerics are very sensitive, hashdist build doesn't pass but conda does")
    def test_damBreak_solver_options(self):
        os.system("parun --TwoPhaseFlow --path " + self.path + " "
                  "damBreak_solver_options.py -l5 -v -C 'final_time=0.1 dt_output=0.1 he=0.1'")
        self.compare_vs_saved_files("damBreak_solver_options")

#    @pytest.mark.skip(reason="long test")
    def test_TwoDimBucklingFlow(self):
        os.system("parun --TwoPhaseFlow --path " + self.path + " "
                  "TwoDimBucklingFlow.py -l5 -v -C 'final_time=0.1 dt_output=0.1 he=0.09'")
        self.compare_vs_saved_files("TwoDimBucklingFlow")

#    @pytest.mark.skip(reason="long test")
    @pytest.mark.skip(reason="need to redo after history revision")
    def test_fillingTank(self):
        os.system("parun --TwoPhaseFlow --path " + self.path + " "
                  "fillingTank.py -l5 -v -C 'final_time=0.02 dt_output=0.02 he=0.01'")
        self.compare_vs_saved_files("fillingTank")

    # *** 3D tests *** #
    def test_marin(self):
        os.system("parun --TwoPhaseFlow --path " + self.path + " "
                  "marin.py -l5 -v -C 'final_time=0.1 dt_output=0.1 he=0.5'")
        self.compare_vs_saved_files("marin")

    def test_moses(self):
        os.system("parun --TwoPhaseFlow --path " + self.path + " "
                  "moses.py -l5 -v -C 'final_time=0.1 dt_output=0.1 he=0.5'")
        self.compare_vs_saved_files("moses")

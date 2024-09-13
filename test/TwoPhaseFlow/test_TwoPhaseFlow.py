#!/usr/bin/env python
"""
Test module for TwoPhaseFlow
"""
import h5py
import pytest
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
        self.path = self._scriptdir

    def teardown_method(self, method):
        """ Tear down function """
        FileList = ['marin.h5','marin.xmf'
                    'moses.h5','moses.xmf'
                    'damBreak.h5','damBreak.xmf'
                    'TwoDimBucklingFlow.h5','TwoDimBucklingFlow.xmf'
                    'filling.h5','filling.xmf'
                    ]
        for file in FileList:
            if os.path.isfile(file):
                os.remove(file)
            else:
                pass

    def compare_vs_saved_files(self,name,write=False):
        actual = h5py.File(name+'.h5','r')

        expected_path = 'comparison_files/' + 'comparison_' + name + '_phi_t2.csv'
        #write comparison file
        if(write):
            np.array(actual['phi_t2']).tofile(os.path.join(self._scriptdir, expected_path),sep=",")
        np.testing.assert_almost_equal(np.fromfile(os.path.join(self._scriptdir, expected_path),sep=","),np.array(actual['phi_t2'][:]).flatten(),decimal=6)

        expected_path = 'comparison_files/' + 'comparison_' + name + '_velocity_t2.csv'
        #write comparison file
        if(write):
            np.array(actual['velocity_t2']).tofile(os.path.join(self._scriptdir, expected_path),sep=",")
        np.testing.assert_almost_equal(np.fromfile(os.path.join(self._scriptdir, expected_path),sep=","),np.array(actual['velocity_t2']).flatten(),decimal=6)

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

    def test_damBreak_hotstart(self):
        os.system("parun --TwoPhaseFlow --path " + self.path + " "
                  "damBreak.py -l5 -v -H -C 'final_time=0.1 dt_output=0.1 he=0.1 hotstart=True'")

    def test_TwoDimBucklingFlow(self):
        os.system("parun --TwoPhaseFlow --path " + self.path + " "
                  "TwoDimBucklingFlow.py -l5 -v -C 'final_time=0.1 dt_output=0.1 he=0.09'")
        self.compare_vs_saved_files("TwoDimBucklingFlow")

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

    @pytest.mark.skip(reason="PUMI is broken")
    def test_damBreak_genPUMI(self):
        os.system("parun --TwoPhaseFlow --genPUMI --path " + self.path + " "
                  "damBreak.py -l5 -v -C 'final_time=0.1 dt_output=0.1 he=0.1'")

    @pytest.mark.skip(reason="PUMI is broken")
    def test_damBreak_runPUMI(self):
        os.system("parun --TwoPhaseFlow --path " + self.path + " "
                  "damBreak_PUMI.py -l5 -v -C 'final_time=0.1 dt_output=0.1 he=0.1 adapt=0'")
        self.compare_vs_saved_files("damBreak_PUMI")

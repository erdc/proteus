#!/usr/bin/env python
"""
Test module for level set transport
"""
from __future__ import print_function
from builtins import range
from builtins import object
from proteus.iproteus import *
import os
import numpy as np
import tables

class TestRotation2D(object):

    @classmethod
    def setup_class(cls):
        pass

    @classmethod
    def teardown_class(cls):
        pass

    def setup_method(self,method):
        self.aux_names = []
        self.meshdir = os.path.dirname(os.path.abspath(__file__))
        self._scriptdir = os.path.dirname(os.path.abspath(__file__))
        
    def teardown_method(self,method):
        filenames = []
        for aux_name in self.aux_names:
            filenames.extend([aux_name+'.'+ext for ext in ['h5','xmf']])
        filenames.append('proteus.log')
        for f in filenames:
            if os.path.exists(f):
                try:
                    os.remove(f)
                except OSError as e:
                    print ("Error: %s - %s" %(e.filename,e.strerror))
            else:
                pass
            
    def test_rotation2D(self,use_strong_constraints=False):
        from proteus import default_s, default_so, default_p, default_n
        reload(default_s)
        reload(default_so)
        reload(default_p)
        reload(default_n)
        from . import (ls_rotation_2d_p,
                       redist_rotation_2d_p,
                       vof_rotation_2d_p,
                       ls_consrv_rotation_2d_p,
                       ls_rotation_2d_n,
                       redist_rotation_2d_n,
                       vof_rotation_2d_n,
                       ls_consrv_rotation_2d_n,
                       ls_rotation_2d_so)
        opts.logLevel=7
        opts.verbose=True
        opts.profile=True
        opts.gatherArchive=True
        sList=[]
        if ls_rotation_2d_so.sList == []:
            for i in range(len(ls_rotation_2d_so.pnList)):
                s = default_s
                sList.append(s)
        else:
            sList = ls_rotation_2d_so.sList
        ns = NumericalSolution.NS_base(ls_rotation_2d_so,
                                       [ls_rotation_2d_p,
                                        redist_rotation_2d_p,
                                        vof_rotation_2d_p,
                                        ls_consrv_rotation_2d_p],
                                       [ls_rotation_2d_n,
                                        redist_rotation_2d_n,
                                        vof_rotation_2d_n,
                                        ls_consrv_rotation_2d_n],
                                       sList,
                                       opts)
        ns.calculateSolution(ls_rotation_2d_so.name)
        self.aux_names.append(ls_rotation_2d_so.name)
        # COMPARE VS SAVED FILES #
        actual = tables.open_file(ls_rotation_2d_so.name+'.h5','r')
        expected_path = 'comparison_files/' + 'comparison_' + ls_rotation_2d_so.name + '_u_t10.csv'
        #write comparison file
        #np.array(actual.root.u_t10).tofile(os.path.join(self._scriptdir, expected_path),sep=",")
        np.testing.assert_almost_equal(np.fromfile(os.path.join(self._scriptdir, expected_path),sep=","),np.array(actual.root.u_t10).flatten(),decimal=10)

        expected_path = 'comparison_files/' + 'comparison_' + ls_rotation_2d_so.name + '_phid_t10.csv'
        #write comparison file
        #np.array(actual.root.phid_t10).tofile(os.path.join(self._scriptdir, expected_path),sep=",")
        np.testing.assert_almost_equal(np.fromfile(os.path.join(self._scriptdir, expected_path),sep=","),np.array(actual.root.phid_t10).flatten(),decimal=10)

        expected_path = 'comparison_files/' + 'comparison_' + ls_rotation_2d_so.name + '_vof_t10.csv'
        #write comparison file
        #np.array(actual.root.vof_t10).tofile(os.path.join(self._scriptdir, expected_path),sep=",")
        np.testing.assert_almost_equal(np.fromfile(os.path.join(self._scriptdir, expected_path),sep=","),np.array(actual.root.vof_t10).flatten(),decimal=10)

        actual.close()
        del ns
        
if __name__ == '__main__':
    pass

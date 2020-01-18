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
from . import (ls_vortex_3d_p,
               redist_vortex_3d_p,
               vof_vortex_3d_p,
               ls_consrv_vortex_3d_p,
               ls_vortex_3d_n,
               redist_vortex_3d_n,
               vof_vortex_3d_n,
               ls_consrv_vortex_3d_n,
               ls_vortex_3d_so)

class TestVortex3D(object):

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
            
    def test_vortex3D(self,use_strong_constraints=False):
        from proteus import default_s
        reload(default_s)
        opts.logLevel=7
        opts.verbose=True
        opts.profile=True
        opts.gatherArchive=True
        sList=[]
        if ls_vortex_3d_so.sList == []:
            for i in range(len(ls_vortex_3d_so.pnList)):
                s = default_s
                sList.append(s)
        else:
            sList = ls_vortex_3d_so.sList

        ns = NumericalSolution.NS_base(ls_vortex_3d_so,
                                       [ls_vortex_3d_p,
                                        redist_vortex_3d_p,
                                        vof_vortex_3d_p,
                                        ls_consrv_vortex_3d_p],
                                       [ls_vortex_3d_n,
                                        redist_vortex_3d_n,
                                        vof_vortex_3d_n,
                                        ls_consrv_vortex_3d_n],
                                       sList,
                                       opts)
        try:
            ns.calculateSolution(ls_vortex_3d_so.name)
        except:
            assert 0, "Calculate solution failed"
        self.aux_names.append(ls_vortex_3d_so.name)
        # COMPARE VS SAVED FILES #
        actual = tables.open_file(ls_vortex_3d_so.name+'.h5','r')
        expected_path = 'comparison_files/' + 'comparison_3D_u_t80.csv'
        #write comparison file
        #np.array(actual.root.u_t80).tofile(os.path.join(self._scriptdir, expected_path),sep=",")
        np.testing.assert_almost_equal(np.fromfile(os.path.join(self._scriptdir, expected_path),sep=","),np.array(actual.root.u_t80),decimal=10)

        expected_path = 'comparison_files/' + 'comparison_3D_phid_t80.csv'
        #write comparison file
        #np.array(actual.root.phid_t80).tofile(os.path.join(self._scriptdir, expected_path),sep=",")
        np.testing.assert_almost_equal(np.fromfile(os.path.join(self._scriptdir, expected_path),sep=","),np.array(actual.root.phid_t80),decimal=10)

        expected_path = 'comparison_files/' + 'comparison_3D_vof_t80.csv'
        #write comparison file
        #np.array(actual.root.vof_t80).tofile(os.path.join(self._scriptdir, expected_path),sep=",")
        np.testing.assert_almost_equal(np.fromfile(os.path.join(self._scriptdir, expected_path),sep=","),np.array(actual.root.vof_t80),decimal=10)

        actual.close()
        del ns
        
if __name__ == '__main__':
    pass

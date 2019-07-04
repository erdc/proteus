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
from . import (ls_vortex_2d_p,
               redist_vortex_2d_p,
               vof_vortex_2d_p,
               ls_consrv_vortex_2d_p,
               ls_vortex_2d_n,
               redist_vortex_2d_n,
               vof_vortex_2d_n,
               ls_consrv_vortex_2d_n,
               ls_vortex_2d_so)

class TestVortex2D(object):

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
            
    def test_vortex2D(self,use_strong_constraints=False):
        from proteus import default_s
        reload(default_s)
        opts.logLevel=7
        opts.verbose=True
        opts.profile=True
        opts.gatherArchive=True
        sList=[]
        if ls_vortex_2d_so.sList == []:
            for i in range(len(ls_vortex_2d_so.pnList)):
                s = default_s
                sList.append(s)
        else:
            sList = ls_vortex_2d_so.sList

        ns = NumericalSolution.NS_base(ls_vortex_2d_so,
                                       [ls_vortex_2d_p,
                                        redist_vortex_2d_p,
                                        vof_vortex_2d_p,
                                        ls_consrv_vortex_2d_p],
                                       [ls_vortex_2d_n,
                                        redist_vortex_2d_n,
                                        vof_vortex_2d_n,
                                        ls_consrv_vortex_2d_n],
                                       sList,
                                       opts)
        ns.calculateSolution(ls_vortex_2d_so.name)
        self.aux_names.append(ls_vortex_2d_so.name)
        # COMPARE VS SAVED FILES #
        expected_path = ls_vortex_2d_so.name+'_expected.h5'
        expected = tables.open_file(os.path.join(self._scriptdir,expected_path))
        actual = tables.open_file(ls_vortex_2d_so.name+'.h5','r')
        assert np.allclose(expected.root.u_t80,
                           actual.root.u_t80,
                           atol=1e-10)
        assert np.allclose(expected.root.phid_t80,
                           actual.root.phid_t80,
                           atol=1e-10)
        assert np.allclose(expected.root.vof_t80,
                           actual.root.vof_t80,
                           atol=1e-10)
        expected.close()
        actual.close()
        del ns
        
    def test_vortex2D_exactHeaviside(self,use_strong_constraints=False):
        from proteus import default_s
        import vortex2D
        reload(default_s)
        reload(vortex2D)
        vortex2D.useExact=True
        reload(ls_consrv_vortex_2d_p)
        opts.logLevel=7
        opts.verbose=True
        opts.profile=True
        opts.gatherArchive=True
        sList=[]
        if ls_vortex_2d_so.sList == []:
            for i in range(len(ls_vortex_2d_so.pnList)):
                s = default_s
                sList.append(s)
        else:
            sList = ls_vortex_2d_so.sList
        ls_vortex_2d_so.name +='_exactHeaviside'
        ns = NumericalSolution.NS_base(ls_vortex_2d_so,
                                       [ls_vortex_2d_p,
                                        redist_vortex_2d_p,
                                        vof_vortex_2d_p,
                                        ls_consrv_vortex_2d_p],
                                       [ls_vortex_2d_n,
                                        redist_vortex_2d_n,
                                        vof_vortex_2d_n,
                                        ls_consrv_vortex_2d_n],
                                       sList,
                                       opts)
        ns.calculateSolution(ls_vortex_2d_so.name)
        self.aux_names.append(ls_vortex_2d_so.name)
        # COMPARE VS SAVED FILES #
        expected_path = ls_vortex_2d_so.name+'_expected.h5'
        expected = tables.open_file(os.path.join(self._scriptdir,expected_path))
        actual = tables.open_file(ls_vortex_2d_so.name+'.h5','r')
        assert np.allclose(expected.root.u_t80,
                           actual.root.u_t80,
                           atol=1e-10)
        assert np.allclose(expected.root.phid_t80,
                           actual.root.phid_t80,
                           atol=1e-10)
        assert np.allclose(expected.root.vof_t80,
                           actual.root.vof_t80,
                           atol=1e-10)
        expected.close()
        actual.close()
        del ns
        
if __name__ == '__main__':
    pass

#!/usr/bin/env python
"""
Test module for linear boundary value problems (serial)

This module solves equations of the form

.. math::

  \nabla \cdot \left( a(x) \nabla u \right) = f(x)

"""
from proteus.iproteus import *
import os
import numpy as np
import tables
from . import re_gl_6_3d_p
from . import re_gl_6_3d_n
from . import sm_gl_6_3d_p
from . import sm_gl_6_3d_n

class TestRichards():

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
                except OSError, e:
                    print ("Error: %s - %s" %(e.filename,e.strerror))
            else:
                pass
            
    def test_richards(self,use_strong_constraints=False):
        pList = [re_gl_6_3d_p]
        nList = [re_gl_6_3d_n]
        reload(default_so)
        so = default_so
        so.name = pList[0].name = "richards"
        reload(default_s)
        so.sList=[default_s]
        opts.logLevel=7
        opts.verbose=True
        opts.profile=True
        opts.gatherArchive=True
        ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
        ns.calculateSolution(so.name)
        self.aux_names.append(so.name)
        # COMPARE VS SAVED FILES #
        expected_path = so.name+'_expected.h5'
        expected = tables.open_file(os.path.join(self._scriptdir,expected_path))
        actual = tables.open_file(so.name+'.h5','r')
        assert np.allclose(expected.root.pressure_head_t1,
                           actual.root.pressure_head_t1,
                           atol=1e-10)
        expected.close()
        actual.close()
        del ns
        
    def test_elastoplastic(self,use_strong_constraints=False):
        pList = [sm_gl_6_3d_p]
        nList = [sm_gl_6_3d_n]
        reload(default_so)
        so = default_so
        so.name = pList[0].name = "elastoplastic"
        reload(default_s)
        so.sList=[default_s]
        opts.logLevel=7
        opts.verbose=True
        opts.profile=True
        opts.gatherArchive=True
        ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
        ns.calculateSolution(so.name)
        self.aux_names.append(so.name)
        # COMPARE VS SAVED FILES #
        expected_path = so.name+'_expected.h5'
        expected = tables.open_file(os.path.join(self._scriptdir,expected_path))
        actual = tables.open_file(so.name+'.h5','r')
        assert np.allclose(expected.root.displacement_t1,
                           actual.root.displacement_t1,
                           atol=1e-10)
        expected.close()
        actual.close()
        del ns
        
if __name__ == '__main__':
    pass

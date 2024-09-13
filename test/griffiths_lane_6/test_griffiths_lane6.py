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
import h5py
from . import re_gl_6_3d_p
from . import re_gl_6_3d_n
from . import sm_gl_6_3d_p
from . import sm_gl_6_3d_n

from proteus import Quadrature
from proteus import MeshTools
from proteus import FemTools

class TestRichards(object):

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
        actual = h5py.File(so.name+'.h5','r')
        expected_path = 'comparison_files/' + 'comparison_3D_pressure_t1.csv'
        #write comparison file
        #np.array(actual['pressure_head_t1]).tofile(os.path.join(self._scriptdir, expected_path),sep=",")
        np.testing.assert_almost_equal(np.fromfile(os.path.join(self._scriptdir, expected_path),sep=","),np.array(actual['pressure_head_t1']),decimal=10)

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
        try:
            ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
        except:
            assert 0, "Failed at NS_base"
        ns.calculateSolution(so.name)
        self.aux_names.append(so.name)
        # COMPARE VS SAVED FILES #
        actual = h5py.File(so.name+'.h5','r')
        expected_path = 'comparison_files/' + 'comparison_3D_displacement_t1.csv'
        #write comparison file
        #np.array(actual['displacement_t1']).tofile(os.path.join(self._scriptdir, expected_path),sep=",")
        np.testing.assert_almost_equal(np.fromfile(os.path.join(self._scriptdir, expected_path),sep=","),np.array(actual['displacement_t1']).flatten(),decimal=10)

        actual.close()
        del ns
        
if __name__ == '__main__':
    pass

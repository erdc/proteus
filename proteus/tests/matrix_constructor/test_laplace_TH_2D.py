#!/usr/bin/env python
"""
Test module for 2D Quadrilateral Meshes
"""

import os,sys,inspect

cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
if cmd_folder not in sys.path:
    sys.path.insert(0,cmd_folder)

cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],
                                                              "import_modules")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0,cmd_subfolder)

from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
Profiling.logLevel=7
Profiling.verbose=True
import numpy
import numpy.testing as npt
import laplace_template_TH_2D_8 as L_2d
from nose.tools import ok_ as ok
from nose.tools import eq_ as eq
from nose.tools import set_trace

class TestLaplaceConstruction2D():
    """ Verify construction of 2D laplace operator with Taylor Hood
  
    Test mesh include 8 elements.
    """
    
    @classmethod
    def setup_class(self):
        """ Initialize the test problem """
        self.Laplace_object = L_2d.ns
        self._setRelativePath()
    @classmethod
    def teardown_class(self):
        """ Tear down function """
        FileList = ['Laplace_matrix_test.xmf',
                    'Laplace_matrix_test.h5',
                    'proteus.log',
                    'rdomain.poly',
                    "reference_triangle_2d.ele",
                    "reference_triangle_2d.node",
                    "reference_triangle_2d.face",
                    "reference_triangle_2d.poly"]
        for file in FileList:
            if os.path.isfile(file):
                os.remove(file)
    @classmethod
    def _setRelativePath(self):
        self.scriptdir = os.path.dirname(__file__)

    def test_1(self):
        """ Tests the attachMassOperator function in one-level-transport """
        self.Laplace_object.calculateSolution('tt')
        mm = self.Laplace_object.modelList[0].levelModelList[0]
        op_constructor = LinearSolvers.OperatorConstructor(mm)
        op_constructor.attachLaplaceOperator()
        laplace = LinearAlgebraTools.superlu_sparse_2_dense(op_constructor.LaplaceOperator)
        rel_path = "comparison_files/laplace_TH_mesh.npy"
        comparison_mat_p = numpy.load(os.path.join(self.scriptdir,rel_path))
        assert numpy.allclose(laplace,comparison_mat_p,rtol=1e-04,atol=1e-04)

if __name__ == '__main__':
    pass

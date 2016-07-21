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
import BOperator_setup_template_2D as B_2d
from nose.tools import ok_ as ok
from nose.tools import eq_ as eq
from nose.tools import set_trace

class TestMassConstruction2D():
    """ Verify construction of 2D Mass Matrix using transport coefficients """
    @classmethod
    def setup_class(self):
        """ Initialize the test problem """
        self.Boperator_object = B_2d.ns
        self._setRelativePath()

    @classmethod
    def teardown_class(self):
        """ Tear down function """
        FileList = ['Laplace_matrix_test.xmf',
                    'Laplace_matrix_test.h5',
                    'proteus.log',
                    'reference_triangle_2d.ele',
                    'reference_triangle_2d.node',
                    'reference_triangle_2d.poly']
        for file in FileList:
            if os.path.isfile(file):
                os.remove(file)
    @classmethod
    def _setRelativePath(self):
        self.scriptdir = os.path.dirname(__file__)

    def test_1(self):
        """ An initial test of the coefficient class. """
        self.Boperator_object.modelList[0].levelModelList[0].calculateCoefficients()
        rowptr, colind, nzval = self.Boperator_object.modelList[0].levelModelList[0].jacobian.getCSRrepresentation()
        self.Boperator_object.modelList[0].levelModelList[0].scale_dt = False
        self.Asys_rowptr = rowptr.copy()
        self.Asys_colptr = colind.copy()
        self.Asys_nzval = nzval.copy()
        nn = len(self.Asys_rowptr)-1
        self.Asys = LinearAlgebraTools.SparseMatrix(nn,nn,
                                                    self.Asys_nzval.shape[0],
                                                    self.Asys_nzval,
                                                    self.Asys_colptr,
                                                    self.Asys_rowptr)
        self.petsc4py_A = self.Boperator_object.modelList[0].levelModelList[0].getJacobian(self.Asys)
        B_mat = LinearAlgebraTools.superlu_sparse_2_dense(self.petsc4py_A)
        rel_path = "comparison_files/B_mat_reference_element_1.npy"
        comparison_mat = numpy.load(os.path.join(self.scriptdir, rel_path))
        assert numpy.allclose(B_mat,comparison_mat)

    def test_2(self):
        """ Tests the attachMassOperator function in one-level-transport """
        mm = self.Boperator_object.modelList[0].levelModelList[0]
        op = LinearSolvers.OperatorConstructor(mm)
        op.attachBOperator()
        B_mat = LinearAlgebraTools.superlu_sparse_2_dense(op.BOperator)
        rel_path = "comparison_files/B_mat_reference_element_1.npy"
        comparison_mat = numpy.load(os.path.join(self.scriptdir, rel_path))
        assert numpy.allclose(B_mat,comparison_mat)

if __name__ == '__main__':
    tt = TestMassConstruction2D()
    tt.setUp()
    tt.test_2()
    tt.tearDown()

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
import advection_template_TH_2D as A_2d
from nose.tools import ok_ as ok
from nose.tools import eq_ as eq
from nose.tools import set_trace

class TestMassConstruction2D():
    """ Verify construction of 2D Mass Matrix using transport coefficients """
    def __init__(self):
        pass

    def setUp(self):
        """ Initialize the test problem """
        self.Advection_object = A_2d.ns
        
    def tearDown(self):
        """ Tear down function """
        FileList = ['Advection_matrix_test.xmf',
                    'Advection_matrix_test.h5',
                    'proteus.log',
                    "reference_triangle_2d.ele",
                    "reference_triangle_2d.node",
                    "reference_triangle_2d.face",
                    "reference_triangle_2d.poly"]
        for file in FileList:
            if os.path.isfile(file):
                os.remove(file)

    def test_1(self):
        """ An initial test of the advection operator coefficient class """
        self.Advection_object.modelList[0].levelModelList[0].calculateCoefficients()
        rowptr,colind,nzval = self.Advection_object.modelList[0].levelModelList[0].jacobian.getCSRrepresentation()
        self.Advection_object.modelList[0].levelModelList[0].scale_dt = False
        self.Asys_rowptr = rowptr.copy()
        self.Asys_colptr = colind.copy()
        self.Asys_nzval = nzval.copy()
        nn = len(self.Asys_rowptr)-1
        self.Asys = LinearAlgebraTools.SparseMatrix(nn,nn,
                                                    self.Asys_nzval.shape[0],
                                                    self.Asys_nzval,
                                                    self.Asys_colptr,
                                                    self.Asys_rowptr)
        self.petsc4py_A = self.Advection_object.modelList[0].levelModelList[0].getJacobian(self.Asys)
        A = LinearAlgebraTools.superlu_sparse_2_dense(self.petsc4py_A)
        comparison_A = numpy.load('./comparison_files/advection_reference_triangle_2d.npy')
        assert numpy.allclose(A,comparison_A)

    def test_2(self):
        """ Tests the attachAdvectionOperator function """
        mm = self.Advection_object.modelList[0].levelModelList[0]
        u_field = numpy.ones((1,6,2))
        mm.attachAdvectionOperator(u_field)
        comparison_A = numpy.load('./comparison_files/advection_reference_triangle_2d.npy')
        A = LinearAlgebraTools.superlu_sparse_2_dense(mm.AdvectionOperator)
        assert numpy.allclose(A,comparison_A)

if __name__ == '__main__':
    tt = TestMassConstruction2D()
    tt.setUp()
    tt.test_2()
    tt.tearDown()

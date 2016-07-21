#!/usr/bin/env python
"""
Test module for Laplace Matrix builder
"""
import os,sys,inspect

cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
if cmd_folder not in sys.path:
    sys.path.insert(0,cmd_folder)

cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"import_modules")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0,cmd_subfolder)

import os.path
import numpy
from proteus.iproteus import *
from proteus import Comm
from proteus import LinearAlgebraTools
import laplace_template_C0P1_2D as tp_2d
comm = Comm.get()
Profiling.logLevel = 7
Profiling.verbose = True

class TestLaplaceConstruction2D():
    """ Run tests to verify the construction of the 2D Laplace operator """
    @classmethod
    def setup_class(self):
        """ Initialize the test problem """
        self.laplace_object = tp_2d.ns

    @classmethod
    def teardown_class(self):
        """ Tear down function """
        FileList = ["Laplace_matrix_test.h5",
                    "Laplace_matrix_test.xmf",
                    "reference_triangle_2d.ele",
                    "reference_triangle_2d.node",
                    "reference_triangle_2d.poly",
                    "proteus.log",
                    "proteus_default.log"]
        for file in FileList:
            if os.path.isfile(file):
                os.remove(file)

    def test_1(self):
        """ Initial test to check whether this is working """
        self.laplace_object.modelList[0].levelModelList[0].calculateCoefficients()
        rowptr, colind, nzval = self.laplace_object.modelList[0].levelModelList[0].jacobian.getCSRrepresentation()
        self.laplace_object.modelList[0].levelModelList[0].scale_dt = False
        self.Asys_rowptr = rowptr.copy()
        self.Asys_colptr = colind.copy()
        self.Asys_nzval = nzval.copy()
        nn = len(self.Asys_rowptr)-1
        self.Asys = LinearAlgebraTools.SparseMatrix(nn,nn,
                                                    self.Asys_nzval.shape[0],
                                                    self.Asys_nzval,
                                                    self.Asys_colptr,
                                                    self.Asys_rowptr)
        self.petsc4py_A = self.laplace_object.modelList[0].levelModelList[0].getSpatialJacobian(self.Asys)
        laplace_mat = LinearAlgebraTools.superlu_sparse_2_dense(self.petsc4py_A)
        comparison_mat = numpy.loadtxt('./comparison_files/laplace_reference_c0p1_2D.txt')
        assert numpy.allclose(laplace_mat,comparison_mat)

if __name__ == "__main__":
    tlc2d = TestLaplaceConstruction2D()
    tlc2d.setUp()
    tlc2d.test_1()
    tlc2d.tearDown()

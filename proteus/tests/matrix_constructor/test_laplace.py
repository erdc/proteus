#!/usr/bin/env python
"""
Test module for Laplace Matrix builder
"""
import pdb
import os
import os.path
import numpy
from proteus.iproteus import *
from proteus import Comm
from proteus import LinearAlgebraTools
import laplace_setup_template as tp
comm = Comm.get()
Profiling.logLevel = 7
Profiling.verbose = True

class TestLaplaceConstruction():
    """ Run tests to verify the construction of the Laplace operator """
    def __init__(self):
        pass

    def setUp(self):
        """ Initialize the test problem """
        self.laplace_object = tp.ns

    def tearDown(self):
        """ Tear down function """
        FileList = ["Laplace_matrix_test.h5",
                    "Laplace_matrix_test.xmf",
                    "reference_triangle.ele",
                    "reference_triangle.node",
                    "reference_triangle.poly",
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
        comparison_mat = numpy.loadtxt('laplace_reference_c0p1.txt')
        assert numpy.allclose(laplace_mat,comparison_mat)

if __name__ == "__main__":
    tlc = TestLaplaceConstruction()
    tlc.setUp()
    tlc.test_1()
    tlc.tearDown()
    

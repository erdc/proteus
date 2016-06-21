#!/usr/bin/env python
"""
Test module for 2D Advection builder
"""
import os
import os.path
import numpy
from proteus.iproteus import *
from proteus import Comm
from proteus import LinearAlgebraTools
import advection_setup_template_2D as tp_2d
comm = Comm.get()
Profiling.logLevel = 7
Profiling.verbose = True

class TestAdvectionConstruction2D():
    """ Run tests to verify the construction of the 2D Laplace operator """
    def __init__(self):
        pass

    def setUp(self):
        """ Initialize the test problem """
        self.advection_object = tp_2d.ns

    def tearDown(self):
        """ Tear down function """
        FileList = ["2D_advection_matrix_test.h5",
                    "2D_advection__matrix_test.xmf",
                    "reference_triangle_2d.ele",
                    "reference_triangle_2d.node",
                    "reference_triangle_2d.poly",
                    "proteus.log",
                    "proteus_default.log"]
        for file in FileList:
            if os.path.isfile(file):
                os.remove(file)

    def build_advection_field(self):
        """ This function builds a dummy advection field for testing (on-reference-element) """
        nElements = self.advection_object.modelList[0].levelModelList[0].mesh.nElements_global
        nQuadPts =  self.advection_object.modelList[0].levelModelList[0].q['x'].shape[1]
        nd =  self.advection_object.modelList[0].levelModelList[0].coefficients.nd
        self.advection_field = numpy.zeros(shape=(nElements,nQuadPts,nd))
        for ele in range(nElements):
            for i,quad_pt in enumerate(self.advection_object.modelList[0].levelModelList[0].q['x'][0]):
                self.advection_field[ele][i][0] = 0.5 * (quad_pt[0] + quad_pt[1])
                self.advection_field[ele][i][1] = 0.25 * (quad_pt[1])

    def test_1(self):
        """ Initial test to check whether this is working """
        self.build_advection_field()
        self.advection_object.modelList[0].levelModelList[0].coefficients.attach_advection_field(self.advection_field)
        import pdb
        pdb.set_trace()
        self.advection_object.modelList[0].levelModelList[0].calculateCoefficients()
        rowptr, colind, nzval = self.advection_object.modelList[0].levelModelList[0].jacobian.getCSRrepresentation()
        self.advection_object.modelList[0].levelModelList[0].scale_dt = False
        self.Asys_rowptr = rowptr.copy()
        self.Asys_colptr = colind.copy()
        self.Asys_nzval = nzval.copy()
        nn = len(self.Asys_rowptr)-1
        self.Asys = LinearAlgebraTools.SparseMatrix(nn,nn,
                                                    self.Asys_nzval.shape[0],
                                                    self.Asys_nzval,
                                                    self.Asys_colptr,
                                                    self.Asys_rowptr)
        self.petsc4py_A = self.advection_object.modelList[0].levelModelList[0].getSpatialJacobian(self.Asys)
        advection_mat = LinearAlgebraTools.superlu_sparse_2_dense(self.petsc4py_A)
        pdb.set_trace()
        comparison_mat = numpy.loadtxt('advection_reference_c0p1_2D.txt')
        assert numpy.allclose(advection_mat,comparison_mat)

if __name__ == "__main__":
    tlc2d = TestAdvectionConstruction2D()
    tlc2d.setUp()
    tlc2d.test_1()
    tlc2d.tearDown()

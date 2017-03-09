#!/usr/bin/env python
"""
Test module for Laplace Matrix builder
"""

from proteus.iproteus import *
from proteus import LinearAlgebraTools
import proteus.test_utils.TestTools

import os.path,sys,inspect
import numpy

proteus.test_utils.TestTools.addSubFolders(inspect.currentframe())
import TPInv_laplace_template_TH_2D_8 as tp_2d

class TestLaplaceConstruction2D(proteus.test_utils.TestTools.SimulationTest):
    """ Run tests to verify the construction of the 2D Laplace operator """

    def setup_method(self,method):
        """ Intialize the test problem """
        reload(tp_2d)
        self.laplace_object = tp_2d.ns
        self._setRelativePath()

    def teardown_method(self,method):
        """ Tear down function """
        FileList = ["Laplace_matrix_test.h5",
                    "Laplace_matrix_test.xmf",
                    "reference_triangle_2d.ele",
                    "reference_triangle_2d.node",
                    "reference_triangle_2d.poly",
                    "proteus.log",
                    "proteus_default.log"]
        self.remove_files(FileList)

    def _setRelativePath(self):
        self._scriptdir = os.path.dirname(__file__)

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
        expected_output = os.path.dirname(os.path.abspath(__file__)) + '/comparison_files/two_phase_invscaled_laplace_TH_8_expected.data'
        comparison_mat = numpy.load(expected_output)
        assert numpy.allclose(laplace_mat,comparison_mat)

if __name__ == "__main__":
    pass

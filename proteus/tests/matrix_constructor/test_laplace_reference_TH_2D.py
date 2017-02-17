#!/usr/bin/env python
"""
Test module for Laplace Matrix
"""

from proteus.iproteus import *
import proteus.test_utils.TestTools

import os,sys,inspect
import numpy

proteus.test_utils.TestTools.addSubFolders(inspect.currentframe())
import laplace_template_TH_2D as L_2d

class TestMassConstruction2D(proteus.test_utils.TestTools.SimulationTest):
    """ Verify construction of 2D Laplce using transport coefficients """

    def setup_method(self,method):
        """ Initialize the test problem """
        reload(L_2d)
        self.Laplace_object = L_2d.ns

    def teardown_method(self,method):
        """ Tear down function """
        FileList = ['Laplace_matrix_test.xmf',
                    'Laplace_matrix_test.h5',
                    'proteus.log',
                    "reference_triangle_2d.ele",
                    "reference_triangle_2d.node",
                    "reference_triangle_2d.face",
                    "reference_triangle_2d.poly"]
        self.remove_files(FileList)
        
    def _setRelativePath(self):
        self._scriptdir = os.path.dirname(__file__)

    def test_1(self):
        """ An initial test of the coefficient class. """
        self.Laplace_object.modelList[0].levelModelList[0].calculateCoefficients()
        rowptr, colind, nzval = self.Laplace_object.modelList[0].levelModelList[0].jacobian.getCSRrepresentation()
        self.Laplace_object.modelList[0].levelModelList[0].scale_dt = False
        self.Asys_rowptr = rowptr.copy()
        self.Asys_colptr = colind.copy()
        self.Asys_nzval = nzval.copy()
        nn = len(self.Asys_rowptr)-1
        self.Asys = LinearAlgebraTools.SparseMatrix(nn,nn,
                                                    self.Asys_nzval.shape[0],
                                                    self.Asys_nzval,
                                                    self.Asys_colptr,
                                                    self.Asys_rowptr)
        self.petsc4py_A = self.Laplace_object.modelList[0].levelModelList[0].getJacobian(self.Asys)
        Laplace_mat = LinearAlgebraTools.superlu_sparse_2_dense(self.petsc4py_A)
        expected_output = os.path.dirname(os.path.abspath(__file__)) + '/comparison_files/Laplace_mat_reference_element_1.npy'
        comparison_mat = numpy.load(expected_output)
        assert numpy.allclose(Laplace_mat,comparison_mat)

    def test_2(self):
        """ Tests the attachMassOperator function in one-level-transport """
        mm = self.Laplace_object.modelList[0].levelModelList[0]
        mm.attachLaplaceOperator()
        Laplace_mat = LinearAlgebraTools.superlu_sparse_2_dense(mm.LaplaceOperator)
        expected_output = os.path.dirname(os.path.abspath(__file__)) + '/comparison_files/Laplace_mat_reference_element_1.npy'
        comparison_mat = numpy.load(expected_output)
        assert numpy.allclose(Laplace_mat,comparison_mat)

if __name__ == '__main__':
    pass

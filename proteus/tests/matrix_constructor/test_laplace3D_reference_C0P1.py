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
import laplace_template_C0P1_3D as tp_3d
from proteus.iproteus import *
from proteus import Comm
from proteus import LinearAlgebraTools
comm = Comm.get()
Profiling.logLevel = 7
Profiling.verbose = True

class TestLaplaceConstruction3D():
    """ Run tests to verify the construction of the 3D Laplace operator """
    @classmethod
    def setup_class(self):
        """ Initialize the test problem """
        self.laplace_object = tp_3d.ns
        self._setRelativePath()
    @classmethod
    def teardown_class(self):
        """ Tear down function """
        FileList = ["Laplace_matrix_test.h5",
                    "Laplace_matrix_test.xmf",
                    "reference_triangle_3d.ele",
                    "reference_triangle_3d.node",
                    "reference_triangle_3d.face",
                    "reference_triangle_3d.poly",
                    "proteus.log",
                    "proteus_default.log"]
        for file in FileList:
            if os.path.isfile(file):
                os.remove(file)
    @classmethod
    def _setRelativePath(self):
        self.scriptdir = os.path.dirname(__file__)

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
        rel_path = "comparison_files/laplace_reference_c0p1_3D.txt"
        comparison_mat = numpy.loadtxt(os.path.join(self.scriptdir,rel_path))
        assert numpy.allclose(laplace_mat,comparison_mat)



if __name__ == "__main__":
    tlc3d = TestLaplaceConstruction3D()
    tlc3d.setUp()
    tlc3d.test_1()
    tlc3d.tearDown()
    

#!/usr/bin/env python
"""
Test module for 2D Quadrilateral Meshes

ARB -- fix this up...
"""
from proteus.iproteus import *
import proteus.test_utils.TestTools # ARB - need to make a descision about this before merge

import os,sys,inspect
import numpy

proteus.test_utils.TestTools.addSubFolders(inspect.currentframe())
import laplace_template_TH_2D_8 as L_2d

class TestLaplaceConstruction2D(proteus.test_utils.TestTools.SimulationTest):
    """ 
    Verify construction of 2D laplace operator with Taylor Hood. Test mesh 
    include 8 elements. 
    """

    def setup_method(self,method):
        """ Initialize the test problem """
        reload(L_2d)
        self.Laplace_object = L_2d.ns
        self._setRelativePath(__file__)

    def test_1(self):
        """ Tests the attachMassOperator function in one-level-transport """
        self.Laplace_object.calculateSolution('tt')
        mm = self.Laplace_object.modelList[0].levelModelList[0]
        op_constructor = LinearSolvers.OperatorConstructor(mm)
        op_constructor.attachLaplaceOperator()
        laplace = LinearAlgebraTools.superlu_sparse_2_dense(op_constructor.LaplaceOperator)
        rel_path = "comparison_files/laplace_TH_mesh.data"
        comparison_mat_p = numpy.load(os.path.join(self.scriptdir,rel_path))
        assert numpy.allclose(laplace,comparison_mat_p,rtol=1e-04,atol=1e-04)

if __name__ == '__main__':
    pass

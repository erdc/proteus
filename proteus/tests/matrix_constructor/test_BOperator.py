#!/usr/bin/env python
"""
"""
from proteus.iproteus import *
import proteus.test_utils.TestTools

import os,sys,inspect
import numpy

proteus.test_utils.TestTools.addSubFolders(inspect.currentframe())
import BOperator_template_TH_2D_BFS as BFS_B_Op

class TestMassConstruction2D(proteus.test_utils.TestTools.SimulationTest):
    """ Verify construction of 2D Mass Matrix using transport coefficients """
 
    def setup_class(self):
        """ Initialize the test problem """
        reload(BFS_B_Op)
        self.B_object = BFS_B_Op.ns
        self._setRelativePath(__file__)

    def test_1(self):
        self.B_object.calculateSolution('tt')
        mm = self.B_object.modelList[0].levelModelList[0]
        op_constructor = LinearSolvers.OperatorConstructor(mm)
        op_constructor.attachBOperator()
        B_full = LinearAlgebraTools.superlu_sparse_2_dense(op_constructor.BOperator)
        B = B_full[0:10,10:64]
        B_transpose = B_full[10:64,0:10]
        B_actual_path = "comparison_files/B_tri_BFS.data"
        B_transpose_actual_path = "comparison_files/B_transpose_tri_BFS.data"
        B_actual = numpy.load(os.path.join(self.scriptdir,B_actual_path))
        B_transpose_actual = numpy.load(os.path.join(self.scriptdir,B_transpose_actual_path))
        assert numpy.allclose(B,B_actual,rtol=1e-04,atol=1e-04)
        assert numpy.allclose(B_transpose,B_transpose_actual,rtol=1e-04,atol=1e-04)
        

if __name__ == '__main__':
    pass

#!/usr/bin/env python
"""
Test module for 2D Quadrilateral Meshes
"""
import proteus.test_utils.TestTools
from proteus.iproteus import *
from proteus import Comm

import os
import sys
import inspect
import numpy

proteus.test_utils.TestTools.addSubFolders( inspect.currentframe() )
import advection_template_pressure_8 as AP_2d

class TestAdvectionConstruction2D_mesh8(proteus.test_utils.TestTools.SimulationTest):
    """ Verify construction of 2D Mass Matrix using transport coefficients """

    def setup_method(self):
        reload(AP_2d)
        self.Advection_object = AP_2d.ns
        self._setRelativePath()

    def teardown_method(self):
        """ Tear down function """
        FileList = ['Advection_matrix_test.xmf',
                    'Advection_matrix_test.h5',
                    'proteus.log',
                    "reference_triangle_2d.ele",
                    "reference_triangle_2d.node",
                    "reference_triangle_2d.face",
                    "reference_triangle_2d.poly",
                    "rdomain.poly"]
        self.remove_files(FileList)

    def _setRelativePath(self):
        self.scriptdir = os.path.dirname(__file__)

    def test_1(self):
        """ Tests the attachAdvectionOperator function """
        mm = self.Advection_object.modelList[0].levelModelList[0]
        op_constructor = LinearSolvers.OperatorConstructor(mm)
        # need to create the advection field
        func_x = lambda x,y : 2*(x+y)
        func_y = lambda x,y : -y
        u_field = numpy.ones((mm.mesh.nElements_global,6))
        v_field = numpy.ones((mm.mesh.nElements_global,6))
        for k in range(mm.mesh.nElements_global):
            for j,pt in enumerate(mm.q['x'][k]):
                u_field[k][j] = func_x(pt[0],pt[1])
                v_field[k][j] = func_y(pt[0],pt[1])       
        advection_field = [u_field,v_field]
        op_constructor.attachAdvectionOperator(advection_field)
        A = LinearAlgebraTools.superlu_sparse_2_dense(op_constructor.AdvectionOperator)
#         numpy.save('advection_pressure_mesh_8',A[0:9,0:9])
        rel_path = "comparison_files/advection_pressure_mesh_8.npy"
        comparison_A = numpy.load(os.path.join(self.scriptdir,rel_path))
        assert numpy.allclose(A[0:9,0:9],comparison_A)


if __name__ == '__main__':
    pass

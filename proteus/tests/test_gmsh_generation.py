from __future__ import division
from past.utils import old_div
import unittest
import numpy.testing as npt
import numpy as np
from subprocess import check_call
from proteus.Profiling import logEvent
from proteus import Comm, Profiling
from proteus import Domain
from proteus import MeshTools
from proteus import SpatialTools as st

comm = Comm.init()
Profiling.procID = comm.rank()
logEvent("Testing Gmsh Mesh Conversion")

class TestGMSH(unittest.TestCase):
    def test_gmsh_generation_2D(self):
        domain = Domain.PlanarStraightLineGraphDomain()
        domain.vertices = [[0., 0., 0.],
                           [5., 0., 0.],
                           [5., 5., 0.],
                           [0., 5., 0.]]
        domain.segments = [[0, 1], [1, 2], [2, 3], [3, 0]]
        domain.facets = [[[0, 1, 2, 3]]]
        domain.writeGeo('gmsh_mesh_test_2D', he_max=0.1)
        gmsh_cmd = "gmsh {0:s} -v 10 -2 -o {1:s} -format msh2".format(domain.geofile+".geo", domain.geofile+".msh")
        check_call(gmsh_cmd, shell=True)
        MeshTools.msh2simplex(domain.geofile, nd=2)
        with open('gmsh_mesh_test_2D.node', 'r') as nodefile:
            npt.assert_equal(nodefile.readline(), '3433 2 0 1\n')
        with open('gmsh_mesh_test_2D.edge', 'r') as edgefile:
            npt.assert_equal(edgefile.readline(), '10096 1\n')
        with open('gmsh_mesh_test_2D.ele', 'r') as elefile:
            npt.assert_equal(elefile.readline(), '6664 3 1\n')

    def test_gmsh_generation_3D(self):
        domain = Domain.PiecewiseLinearComplexDomain()
        cube = st.Cuboid(domain, dim=[2.,2.,2.])
        st.assembleDomain(domain)
        domain.writeGeo('gmsh_mesh_test_3D', he_max=0.1)
        gmsh_cmd = "gmsh {0:s} -v 10 -3 -o {1:s} -format msh2".format(domain.geofile+".geo", domain.geofile+".msh")
        check_call(gmsh_cmd, shell=True)
        MeshTools.msh2simplex(domain.geofile, nd=3)
        # cek disabling exact tests due to cross-platform non-reproducibility
        # with open('gmsh_mesh_test_3D.node', 'r') as nodefile:
        #     npt.assert_equal(nodefile.readline(), '7674 3 0 1\n')
        # with open('gmsh_mesh_test_3D.edge', 'r') as edgefile:
        #     npt.assert_equal(edgefile.readline(), '9384 1\n')
        # with open('gmsh_mesh_test_3D.face', 'r') as facefile:
        #     npt.assert_equal(facefile.readline(), '6256 3 1\n')
        # with open('gmsh_mesh_test_3D.ele', 'r') as elefile:
        #     npt.assert_equal(elefile.readline(), '37473 4 1\n')
        with open('gmsh_mesh_test_3D.node', 'r') as nodefile:
            assert(old_div(abs(int(nodefile.readline().split()[0]) - 7674),7674.0) < .02)
        with open('gmsh_mesh_test_3D.edge', 'r') as edgefile:
            assert(old_div(abs(int(edgefile.readline().split()[0]) - 9384),9384.0) < .02)
        with open('gmsh_mesh_test_3D.face', 'r') as facefile:
            assert(old_div(abs(int(facefile.readline().split()[0]) - 6256),6256.0) < .02)
        with open('gmsh_mesh_test_3D.ele', 'r') as elefile:
            assert(old_div(abs(int(elefile.readline().split()[0]) - 37473),37473.0) < .02)

if __name__ == '__main__':
    unittest.main(verbosity=2)

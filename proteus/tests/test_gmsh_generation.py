import unittest
import numpy.testing as npt
import numpy as np
from subprocess import check_call
from proteus.Profiling import logEvent
from proteus import Comm, Profiling
from proteus import Domain
from proteus import MeshTools

comm = Comm.init()
Profiling.procID = comm.rank()
logEvent("Testing Gmsh Mesh Conversion")

class TestBC(unittest.TestCase):
    def test_gmsh_generation_2D(self):
        domain = Domain.PlanarStraightLineGraphDomain()
        domain.vertices = [[0., 0., 0.],
                           [5., 0., 0.],
                           [5., 5., 0.],
                           [0., 5., 0.]]
        domain.segments = [[0, 1], [1, 2], [2, 3], [3, 0]]
        domain.facets = [[[0, 1, 2, 3]]]
        domain.writeGeo('gmsh_mesh_test', he_max=0.1)
        gmsh_cmd = "time gmsh {0:s} -v 10 -2 -o {1:s} -format msh".format(domain.geofile+".geo", domain.geofile+".msh")
        check_call(gmsh_cmd, shell=True)
        MeshTools.msh2triangle(domain.geofile)
        with open('gmsh_mesh_test.node', 'r') as nodefile:
            npt.assert_equal(nodefile.readline(), '3425 2 0 1\n')
        with open('gmsh_mesh_test.edge', 'r') as edgefile:
            npt.assert_equal(edgefile.readline(), '10072 1\n')
        with open('gmsh_mesh_test.ele', 'r') as elefile:
            npt.assert_equal(elefile.readline(), '6648 3 1\n')


if __name__ == '__main__':
    unittest.main(verbosity=2)

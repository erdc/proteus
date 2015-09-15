from collections import namedtuple

import numpy.testing as npt
import numpy as np
from nose.tools import eq_, ok_
from proteus.MeshTools import (Node,
                               Element,
                               triangleVerticesToNormals,
                               tetrahedronVerticesToNormals,
                               intersectPoints,
                               intersectEdges,
                               intersectPolyhedron,
                               getMeshIntersections)

def test_Node():
    origin_default = Node()
    origin = Node(nodeNumber=0,x=0.0,y=0.0,z=0.0)
    ppp = Node(nodeNumber=1, x= 0.5, y= 0.5, z= 0.5)
    ppm = Node(nodeNumber=2, x= 0.5, y= 0.5, z=-0.5)
    pmp = Node(nodeNumber=3, x= 0.5, y=-0.5, z= 0.5)
    pmm = Node(nodeNumber=4, x= 0.5, y=-0.5, z=-0.5)
    mpp = Node(nodeNumber=5, x=-0.5, y= 0.5, z= 0.5)
    mpm = Node(nodeNumber=6, x=-0.5, y= 0.5, z=-0.5)
    mmp = Node(nodeNumber=7, x=-0.5, y=-0.5, z= 0.5)
    mmm = Node(nodeNumber=8, x=-0.5, y=-0.5, z=-0.5)
    lexicographic_ordering = [ppp, ppm, pmp, pmm, mpp, mpm, mmp, mmm]
    ok_(origin.N == 0)
    ok_(origin.p[0] == 0.0)
    ok_(origin.p[1] == 0.0)
    ok_(origin.p[2] == 0.0)
    ok_(origin.length == 1.0)
    ok_(origin.diameter == 1.0)
    ok_((origin.unitNormal == Node.xUnitVector).all())
    ok_(origin == origin_default)
    ok_(str(origin) == str(origin_default))
    ok_(hash(origin) == hash(origin_default))
    ok_(ppp > origin)
    for  i in range(8):
        p = lexicographic_ordering[i]
        for pg in lexicographic_ordering[:i]:
            ok_(pg > p)
            ok_(pg >= p)
            ok_(pg != p)
        for pl in lexicographic_ordering[i+1:]:
            ok_(pl < p)
            ok_(pl <= p)
            ok_(pl != p)

def test_Element():
    nodesCube= [Node(0,0.0,0.0,0.0),
                Node(1,0.0,1.0,0.0),
                Node(2,1.0,1.0,0.0),
                Node(3,1.0,0.0,0.0),
                Node(4,0.0,0.0,1.0),
                Node(5,0.0,1.0,1.0),
                Node(6,1.0,1.0,1.0),
                Node(7,1.0,0.0,1.0)]
    e0 = Element(elementNumber=1, nodes = nodesCube)
    nodesCube.sort()
    ok_(e0.N == 1)
    for N, n   in enumerate(nodesCube):
        ok_(e0.nodes[N] == n)

def test_intersect_points():
    npt.assert_equal(intersectPoints(([0, 1], [0, 2]), [[0, 1], [0., 1.5]]),
                                          [([0, 1]), ([0., 1.5])])

    npt.assert_equal(intersectPoints(([0, 1], [0, 1]), [[0, 1], [1., 2.]]),
                                          [([0, 1]), None])

def test_intersect_edges():
    # check that overlaps grab furthest point on edge
    npt.assert_equal(intersectEdges([(0, 1), (0, 2)], [[(0, 1), (0, 2)],]),
                     [[0, 2]])

    # check 3D case
    npt.assert_equal(intersectEdges([(3, 1, 1), (3, 1, 2)], [[(3, 1, 1), (3, 1, 4)],],),
                     [[3, 1, 2]])

    # check a proper 3D intersection
    npt.assert_equal(intersectEdges([(1, 1, 0), (1, 1, 2)], [[(5, 5, 5), (-1, -1, -1)],]),
                     [[1, 1, 1]])

def test_intersect_polyhedron():
    # check 2-D triangle intersection
    l = ((-1, 0.5), (1, 0.5))
    p = (((-1, 0), (0, 0)), ((0, -1), (0, 0)), ((1, 1), (0.5, 0.5)))
    npt.assert_equal(intersectPolyhedron(l, p),
                     ((0., 0.5), ((0.5, 0.5))))

    # check 2-D square intersection in 3 space
    l = ((0, 0, 0), (2, 2, 0))
    p = (((-1, 0, 0), (0, 0, 0)),
         ((0, -1, 0), (0, 0, 0)),
         ((1, 0, 0), (1, 1, 0)),
         ((0, 1, 0), (1, 1, 0)))
    npt.assert_equal(intersectPolyhedron(l, p),
                     ((0., 0., 0.), ((1., 1., 0.))))


def test_triangle_vertices_to_normals():
    npt.assert_equal(triangleVerticesToNormals(([0, 0, 0], [0, 1, 0], [1, 0, 0])),
                     [(([-1.,  0.,  0.]), ([0, 0, 0])),
                      (([-0., -1., -0.]), ([0, 0, 0])),
                      (([ 1.,  1.,  0.]), ([0, 1, 0]))])

def test_tetrahedron_vertices_to_normals():
    npt.assert_equal(tetrahedronVerticesToNormals(([0, 0, 0], [0, 1, 0], [1, 0, 0], [0, 0, 1])),
                                                       [(([ 0,  0, -1]), ([0, 0, 0])),
                                                        (([-1,  0,  0]), ([0, 0, 0])),
                                                        (([ 0, -1,  0]), ([0, 0, 0])),
                                                        (([1, 1, 1]), ([0, 1, 0]))])

def test_mesh_intersections():
    """
    Test correctness of intersecting a line and a tetrahedral mesh.
    """

    endpoints = np.asarray(((0, 0, 0), (1, 1, 1)))
    mesh_type = namedtuple('mesh', 'nodeArray elementNodesArray')
    nodeArray = np.asarray(([[ 0. ,  0. ,  0. ],
       [ 0.5,  0. ,  0. ],
       [ 1. ,  0. ,  0. ],
       [ 0. ,  0.5,  0. ],
       [ 0.5,  0.5,  0. ],
       [ 1. ,  0.5,  0. ],
       [ 0. ,  1. ,  0. ],
       [ 0.5,  1. ,  0. ],
       [ 1. ,  1. ,  0. ],
       [ 0. ,  0. ,  0.5],
       [ 0.5,  0. ,  0.5],
       [ 1. ,  0. ,  0.5],
       [ 0. ,  0.5,  0.5],
       [ 0.5,  0.5,  0.5],
       [ 1. ,  0.5,  0.5],
       [ 0. ,  1. ,  0.5],
       [ 0.5,  1. ,  0.5],
       [ 1. ,  1. ,  0.5],
       [ 0. ,  0. ,  1. ],
       [ 0.5,  0. ,  1. ],
       [ 1. ,  0. ,  1. ],
       [ 0. ,  0.5,  1. ],
       [ 0.5,  0.5,  1. ],
       [ 1. ,  0.5,  1. ],
       [ 0. ,  1. ,  1. ],
       [ 0.5,  1. ,  1. ],
       [ 1. ,  1. ,  1. ]]))

    elementNodesArray = np.asarray([[ 0,  1,  4, 13],
       [ 0,  1, 10, 13],
       [ 0,  3,  4, 13],
       [ 0,  3, 12, 13],
       [ 0,  9, 10, 13],
       [ 0,  9, 12, 13],
       [ 1,  2,  5, 14],
       [ 1,  2, 11, 14],
       [ 1,  4,  5, 14],
       [ 1,  4, 13, 14],
       [ 1, 10, 11, 14],
       [ 1, 10, 13, 14],
       [ 3,  4,  7, 16],
       [ 3,  4, 13, 16],
       [ 3,  6,  7, 16],
       [ 3,  6, 15, 16],
       [ 3, 12, 13, 16],
       [ 3, 12, 15, 16],
       [ 4,  5,  8, 17],
       [ 4,  5, 14, 17],
       [ 4,  7,  8, 17],
       [ 4,  7, 16, 17],
       [ 4, 13, 14, 17],
       [ 4, 13, 16, 17],
       [ 9, 10, 13, 22],
       [ 9, 10, 19, 22],
       [ 9, 12, 13, 22],
       [ 9, 12, 21, 22],
       [ 9, 18, 19, 22],
       [ 9, 18, 21, 22],
       [10, 11, 14, 23],
       [10, 11, 20, 23],
       [10, 13, 14, 23],
       [10, 13, 22, 23],
       [10, 19, 20, 23],
       [10, 19, 22, 23],
       [12, 13, 16, 25],
       [12, 13, 22, 25],
       [12, 15, 16, 25],
       [12, 15, 24, 25],
       [12, 21, 22, 25],
       [12, 21, 24, 25],
       [13, 14, 17, 26],
       [13, 14, 23, 26],
       [13, 16, 17, 26],
       [13, 16, 25, 26],
       [13, 22, 23, 26],
       [13, 22, 25, 26]])

    mesh = mesh_type(nodeArray=nodeArray, elementNodesArray=elementNodesArray)
    toPolyhedron = tetrahedronVerticesToNormals
    eq_(getMeshIntersections(mesh, toPolyhedron, endpoints),
        set([((0.5, 0.5, 0.5), (1, 1, 1)), ((0, 0, 0), (0.5, 0.5, 0.5))]))

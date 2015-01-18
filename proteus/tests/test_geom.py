import numpy.testing as npt

from proteus import Geom

def test_intersect_points():
    npt.assert_equal(Geom.intersectPoints(([0, 1], [0, 2]), [[0, 1], [0., 1.5]]),
                                          [([0, 1]), ([0., 1.5])])

    npt.assert_equal(Geom.intersectPoints(([0, 1], [0, 1]), [[0, 1], [1., 2.]]),
                                          [([0, 1]), None])

def test_intersect_edges():
    # check that overlaps grab furthest point on edge
    npt.assert_equal(Geom.intersectEdges([(0, 1), (0, 2)], [[(0, 1), (0, 2)],]),
                     [[0, 2]])

    # check 3D case
    npt.assert_equal(Geom.intersectEdges([(3, 1, 1), (3, 1, 2)], [[(3, 1, 1), (3, 1, 4)],],),
                     [[3, 1, 2]])

    # check a proper 3D intersection
    npt.assert_equal(Geom.intersectEdges([(1, 1, 0), (1, 1, 2)], [[(5, 5, 5), (-1, -1, -1)],]),
                     [[1, 1, 1]])

def test_intersect_polyhedron():
    # check 2-D triangle intersection
    l = ((-1, 0.5), (1, 0.5))
    p = (((-1, 0), (0, 0)), ((0, -1), (0, 0)), ((1, 1), (0.5, 0.5)))
    npt.assert_equal(Geom.intersectPolyhedron(l, p),
                     ((0., 0.5), ((0.5, 0.5))))

    # check 2-D square intersection in 3 space
    l = ((0, 0, 0), (2, 2, 0))
    p = (((-1, 0, 0), (0, 0, 0)),
         ((0, -1, 0), (0, 0, 0)),
         ((1, 0, 0), (1, 1, 0)),
         ((0, 1, 0), (1, 1, 0)))
    npt.assert_equal(Geom.intersectPolyhedron(l, p),
                     ((0., 0., 0.), ((1., 1., 0.))))


def test_triangle_vertices_to_normals():
    npt.assert_equal(Geom.triangleVerticesToNormals(([0, 0, 0], [0, 1, 0], [1, 0, 0])),
                     [(([-1.,  0.,  0.]), ([0, 0, 0])),
                      (([-0., -1., -0.]), ([0, 0, 0])),
                      (([ 1.,  1.,  0.]), ([0, 1, 0]))])

def test_tetrahedron_vertices_to_normals():
    npt.assert_equal(Geom.tetrahedronVerticesToNormals(([0, 0, 0], [0, 1, 0], [1, 0, 0], [0, 0, 1])),
                                                       [(([ 0,  0, -1]), ([0, 0, 0])),
                                                        (([-1,  0,  0]), ([0, 0, 0])),
                                                        (([ 0, -1,  0]), ([0, 0, 0])),
                                                        (([1, 1, 1]), ([0, 1, 0]))])



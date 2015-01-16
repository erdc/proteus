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
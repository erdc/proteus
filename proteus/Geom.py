"""
A set of functions for working with geometric data
"""

import numpy as np

norm = np.linalg.norm

def distance(a, b):
    return norm(b - a)

def intersectPoints(line, points):
    """
    Given a line (defined as two points in three-space), identify all points (defined as points in three space) that
    the line intersects.

    This hasn't been vectorized.
    """

    a, b = line
    a = np.asarray(a)
    b = np.asarray(b)
    distanceAB = distance(a, b)

    def onAB(p):
        p = np.asarray(p)
        eps = 2*np.max((np.max(np.spacing(a)), np.max(np.spacing(b)), np.max(np.spacing(p))))
        distancePA = distance(a, p)
        distancePB = distance(p, b)
        return p if abs(distancePA + distancePB - distanceAB) < eps else None

    return [onAB(p) for p in points]

def intersectEdges(line, edges):
    """
    Given a line (defined as two points in three-space), identify the locations of its intersections with all
    given edges (defined as line segments in three space).  If the line and an edge overlap, the *furthest* point
    along the line (closest to the second point) that is still on the edge is returned.

    This hasn't been vectorized.
    """

    def intersectEdge(line, edge):

        line = np.asarray(line)
        edge = np.asarray(edge)
        a, b = line
        c, d = edge
        v_l = b - a
        v_e = d - c

        vl_cross_ve = np.cross(v_l, v_e)
        mag_vl_cross_ve = norm(vl_cross_ve)

        if mag_vl_cross_ve == 0:
            # lines are parallel, check for overlap
            intersects = intersectPoints(line, edge) + intersectPoints(edge, line)
            # test for an intersect in intersectPoints
            intersect = next((i for i in intersects if i is not None), None)
            if intersect is not None:
                # farthest endpoint is a, so start from there
                closest_endpoint = a
                closest_distance = distance(closest_endpoint, b)
                # could reuse iterator from above, but it's confusing enough as it is :)
                for intersect in intersects:
                    if intersect is None:
                        continue
                    intersect_distance = distance(intersect, b)
                    if intersect_distance < closest_distance:
                        closest_endpoint = intersect
                        closest_distance = intersect_distance
                return closest_endpoint
            else:
                return None

        # lines are not parallel, check for intersection
        vl_cross_ve = np.cross(v_l, v_e)

        # if v_l and v_e intersect, then there is an x that satisfies
        x_vl_cross_ve = np.cross((c - a), v_e)

        # but the two above vectors must be parallel
        if norm(np.cross(vl_cross_ve, x_vl_cross_ve)) > 1e-8:
            return None

        # two lines are parallel, solve for x
        x = norm(x_vl_cross_ve)/norm(vl_cross_ve)

        intersect = a + x*(b-a)

        # and verify intersection is on the line
        points = intersectPoints(line, [intersect])
        assert(len(points) == 1)
        return points[0]

    return [intersectEdge(line, edge) for edge in edges]

def intersectFaces(line, faces):
    """
    Given a line (defined as two points in three-space), identify the locations of all faces (defined as simple
    segments in three space) that the line intersects.
    """
    pass
"""
A set of functions for working with geometric data
"""

from __future__ import division

import numpy as np

norm = np.linalg.norm

def triangleVerticesToNormals(elementVertices):
    """
    Given a set of vertices to a triangle, return normals and a point corresponding to each normal
    """
    elementVertices = np.asarray(elementVertices)

    if norm(elementVertices[:,2]) > 0:
        raise ValueError("Expected triangles in 2D plane, got something else")
    sets = ((0, 1), (0, 2), (1, 2))
    outs = (2, 1, 0)
    faces = []
    rotate = np.asarray(((0., -1., 0.),
                         (1., 0., 0.),
                         (0., 0., 0.)))

    for set, out in zip(sets, outs):
        vertices = elementVertices[[set]]
        ab = vertices[1] - vertices[0]
        v_out = vertices[0] - elementVertices[out]
        normal = rotate.dot(ab)
        # normal should point *away* from remaining point
        if normal.dot(v_out) < 0:
            normal = -1*normal
        faces.append((normal, vertices[0]))
    return faces


def tetrahedronVerticesToNormals(elementVertices):
    """
    Given a set of vertices to a tetrahedron, return normals and a point corresponding to each normal
    """
    elementVertices = np.asarray(elementVertices)

    sets = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))
    outs = (3, 2, 1, 0)

    faces = []

    for set, out in zip(sets, outs):
        vertices = elementVertices[[set]]
        ab = vertices[1] - vertices[0]
        ac = vertices[2] - vertices[0]
        normal = np.cross(ab, ac)
        v_out = vertices[0] - elementVertices[out]
        # normal should point *away* from remaining point
        if normal.dot(v_out) < 0:
            normal = -1*normal
        faces.append((normal, vertices[0]))
    return faces


def distance(a, b):
    return norm(b - a)

def intersectPoints(line, points):
    """
    Given a line segment (defined as two points), identify all points that the line segment intersects.

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
    Given a line segment (defined as two points), identify the locations of its intersections with all
    given edges (defined as line segments).  If the line and an edge overlap, the *furthest* point
    along the line (closest to the second point) that is still on each edge is returned.

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

def intersectPolyhedron(line, polyhedron):
    """
    Given a line (defined as two points), identify the locations that it enters and exits the
    polyhedron (defined as a collection of half-planes in three-space in normal, vertex form)

    If the facets of the polyhedron are in edge form, the normal can be computed by taking the cross product of any
    two non-parallel edges of the facet (in three-space).  Any vertex of the facet will work.

    Implementation of algorithm described here: http://geomalgorithms.com/a13-_intersect-4.html

    This hasn't been vectorized.
    """

    a, b = line
    a, b = np.asarray(a), np.asarray(b)

    if distance(a, b) == 0:
        raise ValueError("Line segment must not have length 0")

    v_l = b - a
    t_e = 0  # location along line entering polyhedron (initial value 0)
    t_l = 1  # location along line leaving polyhedron (initial value 1)

    for plane in polyhedron:
        n, v = plane
        n, v = np.asarray(n), np.asarray(v)
        ndotba = -n.dot(a - v)
        d = n.dot(v_l)
        if d == 0:
            # the line segment is parallel to this face
            if ndotba < 0:
                # the line is outside the face
                return None
            else:
                # the line is in or on the face, ignore this face
                continue
        t = ndotba / d
        if d < 0:
            # segment is entering polyhedron across this facet
            t_e = max(t_e, t)
            if t_e > t_l:
                # segment enters polyhedron after leaving, no intersection
                return None
        else:
            # segment is exiting polyhedron across this facet
            t_l = min(t_l, t)
            if t_l < t_e:
                # segment exits polyhedron before entering, no intersection
                return None

    assert(t_e <= t_l)

    return [a + t_e*v_l, a + t_l*v_l]
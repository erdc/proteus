"""
Module creating predifined or custom shapes. Each shape needs a Domain as
argument (from proteus.Domain). A Domain can contain any number of shapes.
Boundary conditions objects are automatically created for each facet (3D) or
segment (2D) defining the shape.
Classes:
 - Shape: super class, regroups functions common to all shapes
 - Cuboid: creates a 3D cuboid
 - Rectangle: creates a 2D rectangle
 - Custom: creates a custom shape from a given set vertices, facets, etc.


Example
--------
from proteus import Domain
from proteus import SpatialTools as st
import numpy as np

domain = Domain.PlanarStraightLineGraphDomain()
shape1 = st.Rectangle(domain, dim=[0.5, 0.5], coords=[1., 1.])
shape2 = st.Rectangle(domain. dim=[0.3, 0.2], coords=[3., 3.])
shape2.rotate(np.pi/3.)
shape2.BC.left.setNoSlip()

st.buildDomain(domain)
"""

from math import cos, sin, sqrt
import sys
import numpy as np
from proteus import BC as bc
from proteus.Profiling import logEvent as log


class Shape(object):
    """
    Base/super class of all shapes.

    :param domain: domain in which the shape is defined
    """

    def __init__(self, domain, nd=None):
        if nd != domain.nd:
            log('Shape ('+`nd`+'D) and Domain ('+`domain.nd`+'D)' \
                ' have different dimensions!')
            sys.exit()
        self.domain = domain
        domain.shape_list.append(self)
        self.nd = nd
        self.vertices = None
        self.vertexFlags = None
        self.segments = None
        self.segmentFlags = None
        self.facets = None
        self.facetFlags = None
        self.regions = None
        self.regionFlags = None
        self.holes = None
        self.barycenter = np.zeros(3)
        self.coords = None  # Only used for predefined shapes
                            # (can be different from barycenter)
        self.coords_system = np.eye(nd)
        self.b_or = None  # boundary orientation
        self.volume = None
        self.BC_list = []

    def _checkFlags(self, flagSet):
        """
        Checks if flags are set correctly

        :param flagSet: flags to be check (list/array/set)
        """
        flagSet = set(flagSet)
        checkFlag = min(flagSet)
        assert checkFlag == 1, 'Minimum boundary/region tag/flag must be 1'
        for flag in flagSet:
            assert flag == checkFlag, 'Boundary/region tags/flags must be'
            'defined as a suite of numbers with no gap!'
            checkFlag += 1

    def _checkListOfLists(self, list_of_lists):
        """
        Checks if the list of lists has the right dimension

        :param list_of_lists: a list of lists
        """
        assert len(list_of_lists[0]) == self.nd, 'must have be a list of: ' \
            'lists of length ' + self.nd

    def setPosition(self, coords):
        """
        Set position with coords of the barycenter

        :param coords: new set of coordinates for barycenter (list/array)
        """
        old_coords = np.array(self.barycenter)
        if self.domain.nd == 2 and len(old_coords) == 3:
            trans = coords - old_coords[:2]
        else:
            trans = coords - old_coords
        self.translate(trans)

    def setBarycenter(self, barycenter):
        """
        Set barycenter (center of mass) of the shape
        (!) this function does not move the shape

        :param barycenter: global coordinates of barycenter (list/array)
        """
        if self.domain.nd == 2 and len(barycenter) == 2:
            self.barycenter[:2] = barycenter
        else:
            self.barycenter[:] = barycenter

    def setRegions(self, regions, regionFlags=None):
        """
        Sets new regions for the Shape

        :param regions: coordinate of the new region(s) (list/array)
        """
        self._checkListOfLists(regions)
        if regionFlags is not None:
            self._checkFlags(regionFlags)
        else:
            regionFLags = self.regionFlags
        assert len(regions) == len(regionFlags), 'regions and regionFLags'\
            'must have the same length'
        self.regions = np.array(regions)
        self.regionFlags = np.array(regionFlags)

    def setHoles(self, holes):
        """
        Sets a 'hole' in the mesh. The region where the hole is defined will
        not be meshed.

        :param holes: set of coordinates of holes (list/array)
        """
        self._checkListOfLists(holes)
        self.holes = np.array(holes)

    def rotate(self, rot, axis=(0, 0, 1), pivot=None):
        """
        Function to rotate Shape

        :param rot: angle of rotation ()in radians)
        :param axis: axis of rotation (list/array)
        :param pivot: point around which the Shape rotates (list/array)
        -----------------
        Rotated parameters:
        - vertices
        - holes
        - regions
        - local coordinate system
        - boundary orientations
        - coords (if not None)
        - barycenters
        """
        # This function and rotate2D/rotate3D could be optimized
        rot = float(rot)
        nd = self.nd
        if pivot is None:
            pivot = self.barycenter
        if nd == 2:
            pivot = pivot[:2]
            self.vertices[:] = rotation2D(self.vertices, rot, pivot)
            if self.holes is not None:
                self.holes[:] = rotation2D(self.holes, rot, pivot)
            if self.regions is not None:
                self.regions[:] = rotation2D(self.regions, rot, pivot)
            self.barycenter[:2] = rotation2D(self.barycenter[:nd], rot, pivot)
            self.coords_system[:] = rotation2D(self.coords_system, rot,
                                               (0., 0.))
            if self.b_or is not None:
                self.b_or[:] = rotation2D(self.b_or, rot, (0., 0.))
            if self.coords is not None:
                self.coords[:] = rotation2D(self.coords, rot, pivot)
        elif nd == 3:
            self.vertices[:] = rotation3D(self.vertices, rot, axis, pivot)
            if self.holes is not None:
                self.holes[:] = rotation3D(self.holes, rot, axis, pivot)
            if self.regions is not None:
                self.regions[:] = rotation3D(self.regions, rot, axis, pivot)
            self.barycenter[:] = rotation3D(self.barycenter, rot, axis, pivot)
            self.coords_system[:] = rotation3D(self.coords_system, rot, axis,
                                               (0., 0., 0.))
            if self.b_or is not None:
                self.b_or[:] = rotation3D(self.b_or, rot, axis, (0., 0., 0.))
            if self.coords is not None:
                self.coords[:] = rotation3D(self.coords, rot, axis, pivot)

    def translate(self, trans):
        """
        Function to translate Shape

        :param trans: translation values
        -----------------
        Translated parameters:
        - vertices
        - regions
        - coords (if not None)
        - barycenters
        - holes
        """
        self.vertices += trans
        if self.regions is not None:
            self.regions += trans
        if self.coords is not None:
            self.coords += trans
        if self.domain.nd == 2:
            trans2 = (trans[0], trans[1], 0.)
            self.barycenter += trans2
        else:
            self.barycenter += trans
        if self.holes is not None:
            self.holes += trans

    def getPosition(self):
        """
        Returns current position of barycenter
        """
        return self.barycenter

    def getRotation(self):
        """
        Returns local coordinate system relative to global coordinate system
        """
        return self.coords_system


class Cuboid(Shape):
    """
    Class to create a 3D cuboid

    :param domain: domain of the cuboid
    :param dim: dimensions of the cuboid (list/array)
    :param coords: coordinates of the cuboid (list/array)
    """
    count = 0

    def __init__(self, domain, dim=(0., 0., 0.), coords=(0., 0., 0.),
                 barycenter=None):
        super(Cuboid, self).__init__(domain, nd=3)
        self.__class__.count += 1
        self.name = "cuboid" + str(self.__class__.count)
        self.dim = L, W, H = dim  # length, width height
        self.volume = L*W*H
        self.coords = x, y, z = np.array(coords)
        self.vertices = np.array([[x-0.5*L, y-0.5*W, z-0.5*H],
                                  [x-0.5*L, y+0.5*W, z-0.5*H],
                                  [x+0.5*L, y+0.5*W, z-0.5*H],
                                  [x+0.5*L, y-0.5*W, z-0.5*H],
                                  [x-0.5*L, y-0.5*W, z+0.5*H],
                                  [x-0.5*L, y+0.5*W, z+0.5*H],
                                  [x+0.5*L, y+0.5*W, z+0.5*H],
                                  [x+0.5*L, y-0.5*W, z+0.5*H]])
        self.segments = np.array([[0, 1], [1, 2], [2, 3], [3, 0], [4, 5],
                                  [5, 6], [6, 7], [7, 4], [0, 4], [1, 5],
                                  [2, 6], [3, 7]])
        if self.domain.nd == 2:
            self.vertices = np.array([[x-0.5*L, y-0.5*H],
                                      [x+0.5*L, y-0.5*H],
                                      [x+0.5*L, y+0.5*H],
                                      [x-0.5*L, y+0.5*H]])
        self.facets = np.array([[[0, 1, 2, 3]],  # bottom
                                [[0, 1, 5, 4]],  # front
                                [[1, 2, 6, 5]],  # right
                                [[2, 3, 7, 6]],  # back
                                [[3, 0, 4, 7]],  # left
                                [[4, 5, 6, 7]]])  # top
        self.b_or = np.array([[0.,  0., -1.],
                              [-1., 0.,  0.],
                              [0.,  1.,  0.],
                              [1.,  0.,  0.],
                              [0., -1.,  0.],
                              [0.,  0.,  1.]])
        self.regions = np.array([[x, y, z]])
        # defining flags for boundary conditions
        self.facetFlags = np.array([1, 2, 3, 4, 5, 6])
        self.vertexFlags = np.array([1, 1, 1, 1, 6, 6, 6, 6])
        self.segmentFlags = np.array([1, 1, 1, 1, 6, 6, 6, 6, 2, 2, 4, 4])
        self.regionFlags = np.array([1])
        # Initialize (empty) boundary conditions
        self.BC_dict = {'bottom': bc.BoundaryConditions(b_or=self.b_or, b_i=0),
                        'front': bc.BoundaryConditions(b_or=self.b_or, b_i=1),
                        'right': bc.BoundaryConditions(b_or=self.b_or, b_i=2),
                        'back': bc.BoundaryConditions(b_or=self.b_or, b_i=3),
                        'left': bc.BoundaryConditions(b_or=self.b_or, b_i=4),
                        'top': bc.BoundaryConditions(b_or=self.b_or, b_i=5)}
        self.BC_list = [self.BC_dict['bottom'],
                        self.BC_dict['front'],
                        self.BC_dict['right'],
                        self.BC_dict['back'],
                        self.BC_dict['left'],
                        self.BC_dict['top']]
        self.BC = BCContainer(self.BC_dict)
        self.barycenter = np.array(barycenter) or self.coords
        self.It = np.array([[(W**2.+H**2.)/12., 0, 0],
                            [0, (L**2.+H**2.)/12., 0],
                            [0, 0, (W**2.+L**2.)/12.]])

    def setDimensions(self, dim):
        """
        Set dimensions of the shape

        :param dim: new dimensions (list/array)
        """
        self.dim = dim
        L, W, H = dim
        x, y, z = self.coords
        self.vertices[:] = [[x-0.5*L, y-0.5*W, z-0.5*H],
                            [x-0.5*L, y+0.5*W, z-0.5*H],
                            [x+0.5*L, y+0.5*W, z-0.5*H],
                            [x+0.5*L, y-0.5*W, z-0.5*H],
                            [x-0.5*L, y-0.5*W, z+0.5*H],
                            [x-0.5*L, y+0.5*W, z+0.5*H],
                            [x+0.5*L, y+0.5*W, z+0.5*H],
                            [x+0.5*L, y-0.5*W, z+0.5*H]]
        self.volume = L*W*H


class Rectangle(Shape):
    """
    Class to create a rectangle

    :param domain: domain of the rectangle
    :param dim: dimensions of the rectangle (list/array)
    :param coords: coordinates of the rectangle (list/array)
    """
    count = 0

    def __init__(self, domain, dim=(0., 0.), coords=(0., 0.), barycenter=None):
        super(Rectangle, self).__init__(domain, nd=2)
        self.__class__.count += 1
        self.name = "rectangle" + str(self.__class__.count)
        self.dim = L, H = dim  # length, height
        self.coords = x, y = np.array(coords)
        self.vertices = np.array([[x-0.5*L, y-0.5*H],
                                  [x+0.5*L, y-0.5*H],
                                  [x+0.5*L, y+0.5*H],
                                  [x-0.5*L, y+0.5*H]])
        self.segments = np.array([[0, 1], [1, 2], [2, 3], [3, 0]])
        self.barycenter = np.zeros(3)
        if barycenter is not None:
            self.barycenter[0:2] = barycenter[0:2]
        else:
            self.barycenter[0:2] = coords[0:2]
        self.b_or = np.array([[0., -1.],
                              [1., 0.],
                              [0., 1.],
                              [-1., 0.]])
        self.regions = np.array([[x, y]])
        self.segmentFlags = np.array([1, 2, 3, 4])  # bottom, right, top, left
        self.vertexFlags = np.array([1, 1, 3, 3])  # bottom, bottom, top, top
        self.regionFlags = np.array([1])
        self.BC_dict = {'bottom': bc.BoundaryConditions(b_or=self.b_or, b_i=0),
                        'right': bc.BoundaryConditions(b_or=self.b_or, b_i=1),
                        'top': bc.BoundaryConditions(b_or=self.b_or, b_i=2),
                        'left': bc.BoundaryConditions(b_or=self.b_or, b_i=3)}
        self.BC_list = [self.BC_dict['bottom'],
                        self.BC_dict['right'],
                        self.BC_dict['top'],
                        self.BC_dict['left']]
        self.BC = BCContainer(self.BC_dict)
        self.It = L**2+H**2/12

    def setDimensions(self, dim):
        """
        Set dimensions of the shape

        :param dim: new dimensions (list/array)
        """
        self.dim = dim
        L, H = dim
        x, y = self.coords
        self.vertices[:] = [[x-0.5*L, y-0.5*H],
                            [x+0.5*L, y-0.5*H],
                            [x+0.5*L, y+0.5*H],
                            [x-0.5*L, y+0.5*H]]
        self.volume = L*H


class CustomShape(Shape):
    """
    Class to create a custom 2D or 3D shape

    :param domain: domain of the shape
    :param barycenter: barycenter (list/array)
    :param vertices: set of vertices (list/array)
    :param vertexFlags: set of vertex flags (list/array)
    :param segments: set of segments for 2D shape (list/array)
    :param segmentFlags: set of segment flags for 2D shape (list/array)
    :param facets: set of facets for 3D shape (list/array)
    :param facetFlags: set of facet flags for 3D shape (list/array)
    :param holes: set of hole coordinates (list/array)
    :param regions: set of regions of the shape (list/array)
    :param boundaryTags: set of boundary tags to flag shape elements (dict)
    :param boundaryOrientations: set of boundary orientations (list/array)
    """
    count = 0

    def __init__(self, domain, barycenter=None, vertices=None,
                 vertexFlags=None, segments=None, segmentFlags=None,
                 facets=None, facetFlags=None, holes=None, regions=None,
                 regionFlags=None, boundaryTags=None,
                 boundaryOrientations=None):
        print(len(vertices[0]))
        super(CustomShape, self).__init__(domain, nd=len(vertices[0]))
        self.__class__.count += 1
        self.name = "custom" + str(self.__class__.count)
        self._checkFlags(boundaryTags.values())
        self.vertices = np.array(vertices)
        self.vertexFlags = np.array(vertexFlags)
        if segments:
            self.segments = np.array(segments)
            self.segmentFlags = np.array(segmentFlags)
        if facets:
            self.facets = np.array(facets)
            self.facetFlags = np.array(facetFlags)
        if holes is not None:
            self.holes = np.array(holes)
        if regions is not None:
            self._checkFlags(regionFlags)
            self.regions = np.array(regions)
            self.regionFlags = np.array(regionFlags)
        self.BC_dict = {}
        self.BC_list = [None]*len(boundaryTags)
        if boundaryOrientations is not None:
            b_or = []
        else:
            b_or = None
            b_i = None
        for tag, flag in boundaryTags.iteritems():
            b_i = flag-1  # start at index 0
            if boundaryOrientations is not None:
                b_or += [boundaryOrientations[tag]]
            self.BC_dict[tag] = bc.BoundaryConditions(b_or=b_or, b_i=b_i)
            self.BC_list[b_i] = self.BC_dict[tag]
        self.BC = BCContainer(self.BC_dict)
        if barycenter is not None:
            self.barycenter = np.array(barycenter)


class BCContainer(object):
    """
    Creates a class from a dictionary (keys become class variable names)
    """
    def __init__(self, BC_dict):
        self.__dict__ = BC_dict


# --------------------------------------------------------------------------- #
# -------------------------SPATIAL TOOLS FOR SHAPES-------------------------- #
# --------------------------------------------------------------------------- #

def rotation2D(points, rot, pivot=(0., 0.)):
    """
    function to make a set of points/vertices/vectors (arg: points) to rotate
    around a pivot point (arg: pivot)

    :param points: set of 3D points (list or array)
    :param rot: angle of rotation (in radians)
    :param pivot: point around which the set of points rotates (list or array)
    :return points_rot: the rotated set of points (numpy array)
    """
    # function could be optimized
    points = np.array(points)
    rot = float(rot)
    # get coordinates for translation
    x, y = pivot
    # translation matrix
    T = np.array([[1,   0,    0],
                  [0,   1,    0],
                  [-x,  -y,   1]])
    # rotation matrices
    R = np.array([[cos(rot),  sin(rot),  0],
                  [-sin(rot), cos(rot),  0],
                  [0,         0,         1]])
    # full transformation matrix
    M = reduce(np.dot, [T, R, np.linalg.inv(T)])
    # transform points (check also if it is only a 1D array or 2D)
    if points.ndim > 1:
        points_rot = np.ones((len(points), 3))
        points_rot[:, :-1] = points
        points_rot = np.dot(points_rot, M)  # matrix dot product on each vector
        points_rot = points_rot[:, :-1]
    else:
        points_rot = np.ones(3)
        points_rot[:-1] = points
        points_rot = np.dot(points_rot, M)  # matrix dot product on each vector
        points_rot = points_rot[:-1]
    return points_rot


def rotation3D(points, rot, axis=(0., 0., 1.), pivot=(0., 0., 0.)):
    """
    function to make a set of points/vertices/vectors to rotate around an
    arbitrary axis/vector (arg: axis) going through a pivot point.

    :param points: set of 3D points (array)
    :param rot: angle of rotation (in radians)
    :param axis: axis of rotation (list or array)
    :param pivot: point around which the set of points rotates (list or array)
    :return points_rot: the rotated set of points (numpy array)
    """
    # function could be optimized
    points = np.array(points)
    rot = float(rot)
    # get coordinates for translation
    x, y, z = pivot
    # make axis a unity vector
    axis = np.array(axis)
    r = np.linalg.norm(axis)
    axis = axis/r
    # get values for rotation matrix
    cx, cy, cz = axis
    d = sqrt(cy**2+cz**2)
    # rotation matrices
    if d != 0:
        Rx = np.array([[1,         0,        0,    0],
                       [0,         cz/d,     cy/d, 0],
                       [0,         -cy/d,    cz/d, 0],
                       [0,         0,        0,    1]])
    else:  # special case: rotation axis aligned with x axis
        Rx = np.array([[1,         0,        0,    0],
                       [0,         1,        0,    0],
                       [0,         0,        1,    0],
                       [0,         0,        0,    1]])
    Ry = np.array([[d,         0,        cx, 0],
                   [0,         1,        0,   0],
                   [-cx,       0,        d,   0],
                   [0,         0,        0,   1]])
    Rz = np.array([[cos(rot),  sin(rot), 0,   0],
                   [-sin(rot), cos(rot), 0,   0],
                   [0,         0,        1,   0],
                   [0,         0,        0,   1]])
    # translation matrix
    T = np.array([[1,  0,  0,  0],
                  [0,  1,  0,  0],
                  [0,  0,  1,  0],
                  [-x, -y, -z, 1]])
    # full transformation matrix
    inv = np.linalg.inv
    M = reduce(np.dot, [T, Rx, Ry, Rz, inv(Ry), inv(Rx), inv(T)])
    if points.ndim > 1:
        points_rot = np.ones((len(points), 4))
        points_rot[:, :-1] = points
        points_rot = np.dot(points_rot, M)  # matrix dot product on each vector
        points_rot = points_rot[:, :-1]
    else:
        points_rot = np.ones(4)
        points_rot[:-1] = points
        points_rot = np.dot(points_rot, M)  # matrix dot product on each vector
        points_rot = points_rot[:-1]
    return points_rot


def assembleDomain(domain):
    """
    This function sets up everything needed for the domain, meshing, and
    AuxiliaryVariables calculations (if any).
    It should always be called after defining and manipulating all the shapes
    to be attached to the domain.

    :param domain: domain to asssemble
    """
    # reinitialize geometry of domain
    domain.vertices = []
    domain.vertexFlags = []
    domain.segments = []
    domain.segmentFlags = []
    domain.facets = []
    domain.facetFlags = []
    domain.holes = []
    domain.regions = []
    domain.regionFlags = []
    # BC at flag 0
    domain.bc = [bc.BoundaryConditions()]
    domain.bc[0].setParallelFlag0()
    # barycenter at flag 0
    domain.barycenters = np.array([[0., 0., 0.]])
    domain.auxiliaryVariables = []
    start_flag = 0
    start_vertex = 0
    for shape in domain.shape_list:
        # --------------------------- #
        # ----- DOMAIN GEOMETRY ----- #
        # --------------------------- #
        start_flag = len(domain.bc)-1
        start_vertex = len(domain.vertices)
        start_region = len(domain.regions)  # indice 0 ignored
        if domain.regionFlags:
            start_rflag = max(domain.regionFlags)
        else:
            start_rflag = 0
        domain.bc += shape.BC_list
        domain.vertices += list(shape.vertices)
        domain.vertexFlags += list(shape.vertexFlags+start_flag)
        barycenters = np.array([shape.barycenter for bco in shape.BC_list])
        domain.barycenters = np.append(domain.barycenters, barycenters, axis=0)
        if shape.segments is not None:
            domain.segments += list(shape.segments+start_vertex)
            domain.segmentFlags += list(shape.segmentFlags+start_flag)
        if shape.facets is not None:
            domain.facets += list(shape.facets+start_vertex).tolist()
            domain.facetFlags += list(shape.facetFlags+start_flag).tolist()
        if shape.regions is not None:
            domain.regions += list(shape.regions)
            domain.regionFlags += list(shape.regionFlags+start_rflag)
        if shape.holes is not None:
            domain.holes += list(shape.holes)
        domain.getBoundingBox()
    # --------------------------- #
    # ----- MESH GENERATION ----- #
    # --------------------------- #
    mesh = domain.MeshOptions
    if mesh.outputFiles['poly'] is True:
        domain.writePoly(mesh.outputFiles['name'])
    if mesh.outputFiles['ply'] is True:
        domain.writePLY(mesh.outputFiles['name'])
    if mesh.outputFiles['asymptote'] is True:
        domain.writeAsymptote(mesh.outputFiles['name'])
    mesh.setTriangleOptions()
    log("""Mesh generated using: tetgen -%s %s"""  %
        (mesh.triangleOptions, domain.polyfile+".poly"))

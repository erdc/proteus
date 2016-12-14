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
shape2.BC_dict["left"].uOfXT = lambda x, t: 0.

st.assembleDomain(domain)

.. inheritance-diagram:: proteus.SpatialTools
   :parts: 1
"""

from math import cos, sin, sqrt
import sys
import numpy as np
from proteus import BoundaryConditions as bc
from .Profiling import logEvent


class Shape(object):
    """
    Base/super class of all shapes.

    Parameters
    ----------
    domain: proteus.Domain.D_base
        Domain class instance that hold all the geometrical informations and
        boundary conditions of the shape.
    nd: Optional[int]
        Number of dimensions of the shape. If not set, will take the number of
        dimensions of the domain.
    BC_class: Optional[proteus.BoundaryConditions.BC_Base]
        Class to use for boundary conditions (e.g.
        proteus.BoundaryConditions.BC_Base or
        proteus.BoundaryConditions.mprans.BC_RANS).
    """

    def __init__(self, domain, nd=None, BC_class=None):
        if nd != domain.nd:
            logEvent('Shape ('+`nd`+'D) and Domain ('+`domain.nd`+'D)' \
                ' have different dimensions!')
            sys.exit()
        self.Domain = domain
        domain.shape_list.append(self)
        self.nd = nd
        self.BC_class = BC_class or bc.BC_Base
        self.vertices = None
        self.vertexFlags = None
        self.segments = None
        self.segmentFlags = None
        self.facets = None
        self.facetFlags = None
        self.regions = None
        self.volumes = None
        self.regionFlags = None
        self.holes = None
        self.holes_ind = None
        self.barycenter = np.zeros(3)
        self.coords = None  # Only used for predefined shapes
                            # (can be different from barycenter)
        self.coords_system = np.eye(nd)
        self.boundaryTags = None
        self.b_or = None  # boundary orientation
        self.volume = None
        self.children = {}
        self.BC_list = []

    def _checkFlags(self, flagSet):
        """
        Checks if flags are set correctly

        Parameters
        ----------
        flagSet: list
            List of flags.
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

        Parameters
        ----------
        list_of_lists: list
        """
        assert len(list_of_lists[0]) == self.nd, 'must have be a list of: ' \
            'lists of length ' + self.nd

    def _checkNd(self, array):
        """
        Checks if an array is of the same dimension of the Shape instance or
        is of dimension 3

        Parameters
        ----------
        array: array_like
        """
        assert len(array) == self.nd or len(array) == 3, 'wrong dimension'

    def setPosition(self, coords):
        """
        Set position with coords of the barycenter

        Parameters
        ----------
        coords: array_like
            New set of coordinates for barycenter (list/array).
        """
        old_coords = np.array(self.barycenter)
        if self.Domain.nd == 2 and len(old_coords) == 3:
            trans = coords - old_coords[:2]
        else:
            trans = coords - old_coords
        self.translate(trans)

    def setBarycenter(self, barycenter):
        """
        Set barycenter (center of mass) of the shape
        (!) this function does not move the shape

        Parameters
        ----------
        barycenter: array_like
            Global coordinates of barycenter (list/array).
        """
        if self.Domain.nd == 2 and len(barycenter) == 2:
            self.barycenter[:2] = barycenter
        else:
            self.barycenter[:] = barycenter

    def setRegions(self, regions, regionFlags=None):
        """
        Sets new regions for the Shape

        Parameters
        ----------
        regions: array_like
            Array of coordinates of regions.
        regionFlags: Optional[array_like]
            Array of flags.
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

    def setChildShape(self, shape, ind=0):
        """
        Indicates that a shape is contained in this shape
        (this function or setParentShape is necessary for gmsh)

        Parameters
        ----------
        shape: proteus.SpatialTools.Shape
            child of this shape
        ind: int
            parent (current) shape's local index of volume (3D) or facet (2D)
            containing the child shape
        """
        if self.children.has_key(ind):
            self.children[ind] += [shape]
        else:
            self.children[ind] = [shape]

    def setParentShape(self, shape, ind=0):
        """
        Indicates that this shape is within another shape
        (this function or setChildShape is necessary for gmsh)
        Sets the shape within 

        Parameters
        ----------
        shape: proteus.SpatialTools.Shape
            the parent of this shape
        ind: int
            parent shape's local index of volume (3D) or facet (2D) containing
            the child (current) shape
        """
        if shape.children.has_key(ind):
            shape.children[ind] += [self]
        else:
            shape.children[ind] = [self]

    def setHoles(self, holes):
        """
        Sets a 'hole' in the mesh. The region where the hole is defined will
        not be meshed.

        Parameters
        ----------
        holes: array_like
            Array of coordinates of holes (list/array).
        """
        self._checkListOfLists(holes)
        self.holes = np.array(holes)

    def rotate(self, rot, axis=(0, 0, 1), pivot=None):
        """
        Function to rotate Shape

        Parameters
        ----------
        rot: float
            Angle of rotation (in radians).
        axis: Optional[array_like]
            Vector used for rotation. Not necessary for rotation in 2D.
        pivot: Optional[array_list]
            Point around which the shape will rotate. If not set, the
            barycenter will be the center of rotation.

        Notes
        -----
        Rotated attributes:
        - vertices
        - holes
        - regions
        - local coordinate system
        - boundary orientations
        - coords (if not None)
        - barycenter
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

        Parameters
        ----------
        trans: array_like
            Translation values.

        Notes
        -----
        Translated attributes:
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
        if self.Domain.nd == 2:
            trans2 = (trans[0], trans[1], 0.)
            self.barycenter += trans2
        else:
            self.barycenter += trans
        if self.holes is not None:
            self.holes += trans

    def getPosition(self):
        """
        Returns
        -------
        barycenter: array_like
            Current position of barycenter.
        """
        return self.barycenter

    def getRotation(self):
        """
        Returns
        -------
        coords_system: array_like
            Local coordinate system relative to global coordinate system.
        """
        return self.coords_system


class Cuboid(Shape):
    """
    Class to create a 3D cuboid

    Parameters
    ----------
    domain: proteus.Domain.D_base
        Domain class instance that hold all the geometrical informations and
        boundary conditions of the shape.
    dim: Optional[array_like]
        Dimensions of the cuboid.
    coords: Optional[array_like]
        Coordinates of the centroid of the shape.
    barycenter: Optional[array_like]
        Coordinates of the barycenter.
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
        if self.Domain.nd == 2:
            self.vertices = np.array([[x-0.5*L, y-0.5*H],
                                      [x+0.5*L, y-0.5*H],
                                      [x+0.5*L, y+0.5*H],
                                      [x-0.5*L, y+0.5*H]])
        self.facets = np.array([[[0, 1, 2, 3]],  # z-
                                [[1, 2, 6, 5]],  # y+
                                [[2, 3, 7, 6]],  # x+
                                [[3, 0, 4, 7]],  # y-
                                [[0, 1, 5, 4]],  # x-
                                [[4, 5, 6, 7]]])  # z+
        self.b_or = np.array([[0.,  0., -1.],
                              [0.,  1.,  0.],
                              [1.,  0.,  0.],
                              [0., -1.,  0.],
                              [-1., 0.,  0.],
                              [0.,  0.,  1.]])
        self.regions = np.array([[x, y, z]])
        self.volumes = np.array([[[0, 1, 2, 3, 4, 5]]])
        # defining flags for boundary conditions
        self.boundaryTags = bt = {'z-': 1,
                                  'y+': 2,
                                  'x+': 3,
                                  'y-': 4,
                                  'x-': 5,
                                  'z+': 6}
        self.facetFlags = np.array([bt['z-'], bt['y+'], bt['x+'],
                                    bt['y-'], bt['x-'], bt['z+']])
        self.vertexFlags = np.array([bt['z-'], bt['z-'], bt['z-'],
                                     bt['z-'], bt['z+'], bt['z+'],
                                     bt['z+'], bt['z+']])
        self.regionFlags = np.array([1])
        # Initialize (empty) boundary conditions
        self.BC = {'z-': self.BC_class(shape=self, name='z-',
                                       b_or=self.b_or, b_i=0),
                   'y+': self.BC_class(shape=self, name='y+',
                                       b_or=self.b_or, b_i=1),
                   'x+': self.BC_class(shape=self, name='x+',
                                       b_or=self.b_or, b_i=2),
                   'y-': self.BC_class(shape=self, name='y-',
                                       b_or=self.b_or, b_i=3),
                   'x-': self.BC_class(shape=self, name='x-',
                                       b_or=self.b_or, b_i=4),
                   'z+': self.BC_class(shape=self, name='z+',
                                       b_or=self.b_or, b_i=5)}
        self.BC_list = [self.BC['z-'],
                        self.BC['y+'],
                        self.BC['x+'],
                        self.BC['y-'],
                        self.BC['x-'],
                        self.BC['z+']]
        # self.BC = BCContainer(self.BC_dict)
        self.barycenter = np.array(barycenter) or np.array(coords)
        self.It = np.array([[(W**2.+H**2.)/12., 0, 0],
                            [0, (L**2.+H**2.)/12., 0],
                            [0, 0, (W**2.+L**2.)/12.]])

    def setDimensions(self, dim):
        """
        Sets dimensions of the shape.

        Parameters
        ----------

        dim: array_like
            New dimensions of the shape.
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

    Parameters
    ----------
    domain: proteus.Domain.D_base
        Domain class instance that hold all the geometrical informations and
        boundary conditions of the shape.
    dim: Optional[array_like]
        Dimensions of the cuboid.
    coords: Optional[array_like]
        Coordinates of the centroid of the shape.
    barycenter: Optional[array_like]
        Coordinates of the barycenter.
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
        self.facets = np.array([[[0, 1, 2, 3]]])
        self.boundaryTags = bt = {'y-': 1,
                                  'x+': 2,
                                  'y+': 3,
                                  'x-': 4}
        self.segmentFlags = np.array([bt['y-'], bt['x+'], bt['y+'],
                                      bt['x-']])  # y-, x+, y+, x-
        self.vertexFlags = np.array([bt['y-'], bt['y-'], bt['y+'],
                                     bt['y+']])  # y-, y-, y+, y+
        self.regionFlags = np.array([1])
        self.facetFlags = self.regionFlags
        self.BC = {'y-': self.BC_class(shape=self, name='y-',
                                       b_or=self.b_or, b_i=0),
                   'x+': self.BC_class(shape=self, name='x+',
                                       b_or=self.b_or, b_i=1),
                   'y+': self.BC_class(shape=self, name='y+',
                                       b_or=self.b_or, b_i=2),
                   'x-': self.BC_class(shape=self, name='x-',
                                       b_or=self.b_or, b_i=3)}
        self.BC_list = [self.BC['y-'],
                        self.BC['x+'],
                        self.BC['y+'],
                        self.BC['x-']]
        # self.BC = BCContainer(self.BC_dict)
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

    Parameters
    ----------
    domain: proteus.Domain.D_base
        Domain class instance that hold all the geometrical informations and
        boundary conditions of the shape.
    barycenter: Optional[array_like]
        Coordinates of the barycenter.
    vertices: array_like
        Array of vertex coordinates.
    vertexFlags: array_like
        Array of vertex flags (used for boundary conditions)
    segments: array_like
        Array of segments (each defined by indice of 2 vertex).
    segmentFlags: array_like
        Array of segment flags (used for boundary conditions)
    facetss: array_like
        Array of facets (defined by clockwise or counterclockwise loop of
        vertices).
    facetFlags: array_like
        Array of facet flags (used for boundary conditions)
    vertices: array_like
        Array of region coordinates.
    regionFlags: array_like
        Array of region flags (used for boundary conditions)
    volumes: array_like
        Surface loops describing volumes (used for gmsh)
    holes: array_like
        Array of holes coordinates (unmeshed regions)
    boundaryTags: dict
        Dictionary of flags (int) as keys, and tags (e.g. string) for BC.
    boundaryOrientations: Optional[array_like]
        Array of orientation of boundaries. Can be used for BC.
    """
    count = 0

    def __init__(self, domain, barycenter=None, vertices=None,
                 vertexFlags=None, segments=None, segmentFlags=None,
                 facets=None, facetFlags=None, holes=None, regions=None,
                 regionFlags=None, volumes=None, boundaryTags=None,
                 boundaryOrientations=None):
        super(CustomShape, self).__init__(domain, nd=len(vertices[0]))
        self.__class__.count += 1
        self.name = "custom" + str(self.__class__.count)
        self._checkFlags(boundaryTags.values())
        self.boundaryTgs = boundaryTags
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
        if volumes is not None:
            self.volumes = np.array(volumes)
        if regions is not None:
            self._checkFlags(regionFlags)
            self.regions = np.array(regions)
            self.regionFlags = np.array(regionFlags)
        self.BC = {}
        self.BC_list = [None]*len(boundaryTags)
        b_or = [None]*len(boundaryTags)
        self.b_or = b_or
        for tag, flag in boundaryTags.iteritems():
            b_i = flag-1  # start at index 0
            if boundaryOrientations is not None:
                b_or[b_i] = boundaryOrientations[tag]
            self.BC[tag] = self.BC_class(shape=self, name=tag, b_or=b_or, b_i=b_i)
            self.BC_list[b_i] = self.BC[tag]
        # self.BC = BCContainer(self.BC_dict)
        if barycenter is not None:
            self.barycenter = np.array(barycenter)


class ShapeSTL(Shape):
    """
    Class to extract geometrical information from STL file

    Parameters
    ----------
    domain: proteus.Domain.D_base
        Domain class instance that hold all the geometrical informations and
        boundary conditions of the shape.
    filename: string
        Name of the stl file.
    """

    def __init__(self, domain, filename):
        super(ShapeSTL, self).__init__(domain, nd=3)
        self.filename = filename
        self.vertices, self.facets, self.facetnormals = getInfoFromSTL(self.filename)
        self.facetFlags = np.ones(len(self.facets))
        self.vertexFlags = np.ones(len(self.vertices))
        self.volumes = np.array([[[i for i in range(len(self.facets))]]])
        self.boundaryTags = {'stl': 1}
        self.BC = {'stl': self.BC_class(shape=self, name='stl')}
        self.BC_list = [self.BC['stl']]
        # self.BC = BCContainer(self.BC_dict)


def getInfoFromSTL(filename):
    """
    Extracts information from STL file and converts it to a Proteus friendly
    format. Duplicate vertices and segments are removed during the process,
    so the shape is ready for meshing.

    Parameters
    ----------
    filename: name of STL file

    Returns
    -------
    vertices: array_like
        Array of vertices that define STL shape (duplicates removed)
    facets: array_like
        Array of facets (loops of 3 vertices)
    facetnormals: array_like
        normal vertors of each facet
    """
    file = open(filename, 'r')
    facetnormals = []
    facet = []
    facets = []
    vertices = []
    vFlag = 0
    for line in file:
        if "vertex" in line:
            word_list = line.split()
            vertex = (word_list[1], word_list[2], word_list[3])
            vertices += [vertex]
            facet += [vFlag]
            vFlag += 1
        if "facet normal" in line:
            word_list = line.split()
            facetnormals += [[word_list[2], word_list[3], word_list[4]]]
        elif "endfacet" in line:
            facets += [[facet]]
            facet = []
        elif "endsolid" in line:
            pass
        elif "solid" in line:
            word_list = line.split()
            name = word_list[1]
    file.close()
    # vertices_u, inverse = np.unique(vertices, return_inverse=True)
    vertices = np.array(vertices).astype(float)
    facets = np.array(facets).astype(int)
    vertices, inverse = unique_rows(vertices)
    facets = inverse[facets]
    facetnormals = np.array(facetnormals).astype(float)
    return vertices, facets, facetnormals

def unique_rows(arr):
    arr = np.array(arr)
    ca = np.ascontiguousarray(arr).view([('', arr.dtype)] * arr.shape[1])
    unique, indices, inverse = np.unique(ca, return_index=True, return_inverse=True)
    # counts = np.bincount(inverse)
    # sort_indices = np.argsort(counts)[::-1]
    # sorted_arr = arr[indices[sort_indices]]
    # sorted_count = counts[sort_indices]
    return (arr[indices], inverse)

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
    Rotates a set of points/vertices/vectors around a pivotal point in 2D.

    Parameters
    ----------
    points: array_like
        Array of point coordinates to rotate.
    rot: float
        Angle of rotation.
    pivot: array_like
        Pivotal point around which the set of points will be rotated.

    Returns
    -------
    points_rot: array_like
        Rotated set of points.
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
    Rotates a set of points/vertices/vectors around a pivotal point in 3D.

    Parameters
    ----------
    points: array_like
        Array of point coordinates to rotate.
    rot: float
        Angle of rotation.
    axis: array_like
        Axis of rotation.
    pivot: array_like
        Pivotal point around which the set of points will be rotated.

    Returns
    -------
    points_rot: array_like
        Rotated set of points.
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

    Parameters
    ----------
    domain: proteus.Domain.D_base
        Domain class instance that hold all the geometrical informations and
        boundary conditions of the shape.
    """
    # reinitialize geometry of domain
    _assembleGeometry(domain, BC_class=bc.BC_Base)
    _generateMesh(domain)

def _assembleGeometry(domain, BC_class):
    """
    Assembles all the geometrical informations of the shapes attached to a
    domain.

    Parameters
    ----------
    domain: proteus.Domain.D_base
        Domain class instance that hold all the geometrical informations and
        boundary conditions of the shape.
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
    domain.volumes = []
    domain.boundaryTags = {}
    domain.reversed_boundaryTags = {}
    # BC at flag 0
    domain.bc = [BC_class(nd=domain.nd)]
    # domain.bc[0].setNonMaterial()
    # barycenter at flag 0
    domain.barycenters = np.array([[0., 0., 0.]])
    start_flag = 0
    start_vertex = 0
    for shape in domain.shape_list:
        # --------------------------- #
        # ----- DOMAIN GEOMETRY ----- #
        # --------------------------- #
        start_flag = shape.start_flag = len(domain.bc)-1
        start_vertex = shape.start_vertex = len(domain.vertices)
        start_segment = shape.start_segment = len(domain.segments)
        start_facet = shape.start_facet = len(domain.facets)
        if shape.boundaryTags:
            for tag, value in shape.boundaryTags.iteritems():
                domain.boundaryTags[shape.name+'_'+str(tag)] = value+start_flag
                domain.reversed_boundaryTags[value+start_flag] = shape.name+'_'+str(tag)            
        if domain.regionFlags:
            start_rflag = max(domain.regionFlags)
        else:
            start_rflag = 0
        domain.bc += shape.BC_list
        # making copies of shape properties before operations/modifications
        vertices = shape.vertices.copy()
        vertexFlags = shape.vertexFlags.copy()
        if shape.segments is not None:
            segments = shape.segments.copy()
        if shape.facets is not None:
            facets = shape.facets.copy()
        # deleting duplicate vertices and updating segment/facets accordingly
        del_v = 0
        for i_s, vertex in enumerate(shape.vertices):
            if vertex.tolist() in domain.vertices:
                vertices = np.delete(vertices, i_s-del_v, axis=0)
                verticesFlags = np.delete(vertexFlags, i_s-del_v)
                i_s -= del_v
                del_v += 1
                i_d = domain.vertices.index(vertex.tolist())
                if shape.segments is not None:
                    for i in np.nditer(segments, op_flags=['readwrite']):
                        if i > i_s:
                            i[...] -= 1
                        elif i == i_s:
                            i[...] = i_d-start_vertex
                if shape.facets is not None:
                    for i in np.nditer(facets, op_flags=['readwrite']):
                        if i > i_s:
                            i[...] -= 1
                        elif i == i_s:
                            i[...] = i_d-start_vertex
        # adding shape geometry to domain
        domain.vertices += vertices.tolist()
        domain.vertexFlags += (vertexFlags+start_flag).tolist()
        barycenters = np.array([shape.barycenter for bco in shape.BC_list])
        domain.barycenters = np.append(domain.barycenters, barycenters, axis=0)
        if shape.segments is not None:
            domain.segments += (segments+start_vertex).tolist()
            domain.segmentFlags += (shape.segmentFlags+start_flag).tolist()
        if shape.regions is not None:
            domain.regions += (shape.regions).tolist()
            domain.regionFlags += (shape.regionFlags+start_rflag).tolist()
        if shape.facets is not None:
            domain.facets += (facets+start_vertex).tolist()
            if shape.nd == 2:  # facet flags are actually region flags if 2D
                domain.facetFlags += (shape.regionFlags+start_rflag).tolist()
            elif shape.nd == 3:
              domain.facetFlags += (shape.facetFlags+start_flag).tolist()
        if shape.holes is not None:
            domain.holes += (shape.holes).tolist()
    
    domain.holes_ind = []
    if domain.nd == 2 and shape.facets is not None:
        domain.facets = []
        for shape in domain.shape_list:
            if shape.holes_ind is not None:
                    domain.holes_ind += (np.array(shape.holes_ind)+shape.start_facet).tolist()
            if shape.facets is not None:
                facets = (shape.facets+shape.start_vertex).tolist()
                f_to_remove = []
                f_to_add = []
                for i, facet in enumerate(shape.facets):
                    if i in shape.children.keys():
                        for child in shape.children[i]:
                            child_seg = {}
                            child_facs = []
                            for child_fac in child.facets:
                                fac = child_fac[0].tolist()
                                for j, v in enumerate(fac):
                                    if (fac[j-1], fac[j]) in child_seg:
                                        fn = child_seg[(fac[j-1], fac[j])]
                                        f_to_merge = child.facets[fn][0].tolist()
                                        # reverse list
                                        f_to_merge = f_to_merge[::-1] 
                                        # shift lists
                                        f_to_merge = list(deque(f_to_merge).rotate(-f_to_merge.index(fac[j-1])))
                                        fac = list(deque(fac).rotate(-fac.index(fac[j])))
                                        new_f = fac[:-1]+f_to_merge[:-1]
                                        new_f = fac[:fac.index(fac[j-1])]+f_to_merge
                                        f_to_remove += [fn]
                                        f_to_add += [new_f]
                                    elif (fac[j], fac[j-1]) in child_seg:
                                        fn = child_seg[(fac[j], fac[i-1])]
                                        f_to_merge = child.facets[fn][0].tolist()
                                        # shift lists
                                        f_to_merge = list(deque(f_to_merge).rotate(-f_to_merge.index(fac[j-1])))
                                        fac = list(deque(fac).rotate(-fac.index(fac[j])))
                                        new_f = fac[:-1]+f_to_merge[:-1]
                                        new_f = fac[:fac.index(fac[j-1])]+f_to_merge
                                        f_to_remove += [fn, j]
                                        f_to_add += [new_f]
                                    else:
                                        child_seg[(fac[j-1], fac[j])] = j
                            facets2 = [f.tolist()[0] for k, f in enumerate(child.facets) if k not in f_to_remove]
                            if f_to_add:
                                facets2 += [f_to_add]
                        new_facets = np.array(facets2)+child.start_vertex
                        facets[i] += (new_facets).tolist()
                domain.facets += facets


    elif domain.nd == 3 and shape.volumes is not None:
        for shape in domain.shape_list:
            shape.start_volume = len(domain.volumes)
            if shape.holes_ind is not None:
                domain.holes_ind += (np.array(shape.holes_ind)+shape.start_volume).tolist()
            if shape.volumes is not None:
                volumes = (shape.volumes+shape.start_facet).tolist()
                volumes_to_remove = []
                for i, volume in enumerate(volumes):
                    # add holes to volumes
                    if shape.children.has_key(i):
                        for child in shape.children[i]:
                            child_vols = []
                            child_facets = []  # facets in volumes
                            for child_vol in child.volumes:
                                # don't need holes of child volumes, only outer shells:
                                child_vols += [(child_vol[0]+child.start_facet).tolist()]
                            if len(child_vols) > 1:
                                # merge volumes that share a facet
                                inter = np.intersect1d(*child_vols)
                                new_child_vols = list(child_vols)
                                for f in inter:
                                    ind_to_remove = set()
                                    merged = []
                                    for ind, vol in enumerate(child_vols):
                                        if f in vol:
                                            ind_to_remove.add(ind)
                                            merged += vol
                                    for ind in sorted(ind_to_remove, reverse=True):
                                        child_vols.remove(ind)
                                    merged.remove(f)
                                    child_vols += [merged]
                        volume += child_vols
                    if shape.holes_ind is not None and i in shape.holes_ind:
                        volumes_to_remove += [i]
                volumes2 = [vol for i, vol in enumerate(volumes) if i not in volumes_to_remove]
                domain.volumes += volumes2



        domain.getBoundingBox()


def _generateMesh(domain):
    """
    Generates tetgen mesh of domain

    Parameters
    ----------
    domain: proteus.Domain.D_base
        Domain class instance that hold all the geometrical informations and
        boundary conditions of the shape.
    """
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
    logEvent("""Mesh generated using: tetgen -%s %s"""  %
        (mesh.triangleOptions, domain.polyfile+".poly"))


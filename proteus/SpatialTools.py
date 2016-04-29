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
shape2.BC.newBC

st.assembleDomain(domain)
"""

from math import cos, sin, sqrt
import sys
import numpy as np
from proteus import BoundaryConditions as bc
from proteus.Profiling import logEvent as log




class Shape(object):
    """
    Base/super class of all shapes.

    :param domain: domain in which the shape is defined
    """

    def __init__(self, domain, nd=None, BC_class=None):
        if nd != domain.nd:
            log('Shape ({0}D) and Domain ({1}D)'
                ' have different dimensions!'
                .format(nd, domain.nd))
            sys.exit()
        self.Domain = domain
        domain.shape_list.append(self)
        self.MeshOptions = MeshOptions(self)
        self.nd = nd
        self.BC_class = BC_class or bc.BC_Base
        self.name = ''
        self.vertices = None
        self.vertexFlags = None
        self.segments = None
        self.segmentFlags = None
        self.facets = None
        self.facetFlags = None
        self.regions = None
        self.regionFlags = None
        self.volumes = None
        self.holes = None
        self.holes_ind = None
        self.barycenter = np.zeros(3)
        self.coords = None  # Only used for predefined shapes
                            # (can be different from barycenter)
        self.children = {}
        self.coords_system = np.eye(nd)
        self.boundaryTags = None
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

    def _checkNd(self, array):
        """
        Checks if an array is of the same dimension of the Shape instance or
        is of dimension 3
        """
        assert len(array) == self.nd or len(array) == 3, 'wrong dimension'

    def setPosition(self, coords):
        """
        Set position with coords of the barycenter

        :param coords: new set of coordinates for barycenter (list/array)
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

        :param barycenter: global coordinates of barycenter (list/array)
        """
        if self.Domain.nd == 2 and len(barycenter) == 2:
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


    def setParentShape(self, shape, ind=0):
        """
        :param ind: list
        """
        if shape.children.has_key(ind):
            shape.children[ind] += [self]
        else:
            shape.children[ind] = [self]

    def setChildShape(self, shape, ind=0):
        """
        :param shape: shape class instance containing this shape
        :param ind: index of parent shape facet/volume containing this shape
        """
        if self.children.has_key(ind):
            self.children[ind] += [shape]
        else:
            self.children[ind] = [shape]

    def setVolumeHoles(self, ind):
        if self.holes_ind is None:
            self.holes_ind = [ind]
        else:
            self.holes_ind += [ind]

    def setHoles(self, coords):
        """
        Sets a 'hole' in the mesh. The region where the hole is defined will
        not be meshed.

        :param holes: set of coordinates of holes (list/array)
        """
        self._checkListOfLists(holes)
        self.holes = np.array(holes)
        self.holes_ind += [ind]

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
        if self.Domain.nd == 2:
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
        if self.Domain.nd == 2:
            self.vertices = np.array([[x-0.5*L, y-0.5*H],
                                      [x+0.5*L, y-0.5*H],
                                      [x+0.5*L, y+0.5*H],
                                      [x-0.5*L, y+0.5*H]])
        self.facets = np.array([[[0, 1, 2, 3]],  # bottom
                                [[1, 2, 6, 5]],  # front
                                [[2, 3, 7, 6]],  # right
                                [[3, 0, 4, 7]],  # back
                                [[0, 1, 5, 4]],  # left
                                [[4, 5, 6, 7]]])  # top
        self.b_or = np.array([[0.,  0., -1.],
                              [0.,  1.,  0.],
                              [1.,  0.,  0.],
                              [0., -1.,  0.],
                              [-1., 0.,  0.],
                              [0.,  0.,  1.]])
        self.regions = np.array([[x, y, z]])
        self.volumes = np.array([[[0, 1, 2, 3, 4, 5]]])
        # defining flags for boundary conditions
        self.boundaryTags = bt = {'bottom': 1,
                                  'front': 2,
                                  'right': 3,
                                  'back': 4,
                                  'left': 5,
                                  'top': 6}
        self.facetFlags = np.array([bt['bottom'], bt['front'], bt['right'],
                                    bt['back'], bt['left'], bt['top']])
        self.vertexFlags = np.array([bt['bottom'], bt['bottom'], bt['bottom'],
                                     bt['bottom'], bt['top'], bt['top'],
                                     bt['top'], bt['top']])
        self.regionFlags = np.array([1])
        # Initialize (empty) boundary conditions
        self.BC_dict = {'bottom': self.BC_class(shape=self, name='bottom',
                                                 b_or=self.b_or, b_i=0),
                        'front': self.BC_class(shape=self, name='front',
                                                b_or=self.b_or, b_i=1),
                        'right': self.BC_class(shape=self, name='right',
                                                b_or=self.b_or, b_i=2),
                        'back': self.BC_class(shape=self, name='back',
                                               b_or=self.b_or, b_i=3),
                        'left': self.BC_class(shape=self, name='left',
                                               b_or=self.b_or, b_i=4),
                        'top': self.BC_class(shape=self, name='top',
                                              b_or=self.b_or, b_i=5)}
        self.BC_list = [self.BC_dict['bottom'],
                        self.BC_dict['front'],
                        self.BC_dict['right'],
                        self.BC_dict['back'],
                        self.BC_dict['left'],
                        self.BC_dict['top']]
        self.BC = BCContainer(self.BC_dict)
        self.barycenter = np.array(barycenter) or np.array(coords)
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
        self.boundaryTags = bt = {'bottom': 1,
                                  'right': 2,
                                  'top': 3,
                                  'left': 4}
        self.segmentFlags = np.array([bt['bottom'], bt['right'], bt['top'],
                                      bt['left']])  # bottom, right, top, left
        self.vertexFlags = np.array([bt['bottom'], bt['bottom'], bt['top'],
                                     bt['top']])  # bottom, bottom, top, top
        self.regionFlags = np.array([1])
        self.BC_dict = {'bottom': self.BC_class(shape=self, name='bottom',
                                                 b_or=self.b_or, b_i=0),
                        'right': self.BC_class(shape=self, name='right',
                                                b_or=self.b_or, b_i=1),
                        'top': self.BC_class(shape=self, name='top',
                                              b_or=self.b_or, b_i=2),
                        'left': self.BC_class(shape=self, name='left',
                                               b_or=self.b_or, b_i=3)}
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
    :param boundaryOrientations: set of boundary orientations (dict)
    """
    count = 0

    def __init__(self, domain, barycenter=None, vertices=None,
                 vertexFlags=None, segments=None, segmentFlags=None,
                 facets=None, facetFlags=None, volumes=None, holes=None,
                 regions=None, regionFlags=None, boundaryTags=None,
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
        self.BC_dict = {}
        self.BC_list = [None]*len(boundaryTags)
        b_or = [None]*len(boundaryTags)
        self.b_or = b_or
        for tag, flag in boundaryTags.iteritems():
            b_i = flag-1  # start at index 0
            if boundaryOrientations is not None:
                b_or[b_i] = boundaryOrientations[tag]
            self.BC_dict[tag] = self.BC_class(shape=self, name=tag, b_or=b_or, b_i=b_i)
            self.BC_list[b_i] = self.BC_dict[tag]
        self.BC = BCContainer(self.BC_dict)
        if barycenter is not None:
            self.barycenter = np.array(barycenter)


class ShapeSTL(Shape):
    def __init__(self, domain, filename):
        super(ShapeSTL, self).__init__(domain, nd=3)
        self.filename = filename
        self.vertices, self.facets, self.facetnormals = getInfoFromSTL(self.filename)
        self.facetFlags = np.ones(len(self.facets))
        self.vertexFlags = np.ones(len(self.vertices))
        self.volumes = np.array([[[i for i in range(len(self.facets))]]])
        self.boundaryTags = {'stl': 1}
        self.BC_dict = {'stl': self.BC_class(shape=self, name='stl')}
        self.BC_list = [self.BC_dict['stl']]
        self.BC = BCContainer(self.BC_dict)

def getInfoFromSTL(filename):
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

class MeshOptions:
    """
    Mesh options for the domain
    """
    def __init__(self, shape):
        self.Shape = shape
        self.constraints = []

    def _addConstraint(self, entity, cons_type, index, variables):
        entities = ['vertex', 'segment', 'facet', 'global']
        assert entity in entities, 'wrong entity'
        cons_types = ['fixed', 'function', 'TFI', 'box']
        assert cons_type in cons_types, 'wrong constraint type'
        assert isinstance(index, (list, tuple)) or index==None,'must pass a list of index'
        assert isinstance(variables, (dict)), 'variables must be a dictionary'
        self.constraints += [{'entity': entity, 'type': cons_type,
                              'index': index, 'variables':variables}]

    def refineAroundVertex(self, ind, lc_min, lc_max=None, dist_min=None,
                            dist_max=None):
        """
        Refinement (or coarsening) of mesh around vertices.
        tetgen: only lc_min is taken into account
        gmsh: lc_min on segment, lc_max away from vertex, with a transition
              zone between dist_min and dist_max
        :param ind: list of local index of vertices
        :param lc_min: size of element at vertex
        :param lc_max: size of element away from vertex
        :param dist_min: distance away from vertex with lc_min element size
        :param dist_max: distance away from vertex with lc_max element size
        """
        var_dict = {'LcMin': lc_min, 'LcMax': lc_max,
                    'DistMin': dist_min, 'DistMax': dist_max}
        self._addConstraint(entity='vertex', cons_type='fixed',
                            index=ind, variables=var_dict)

    def refineAroundSegment(self, ind, lc_min, lc_max=None, dist_min=None,
                             dist_max=None):
        """
        Refinement (or coarsening) of mesh around segments.
        tetgen: only lc_min is taken into account
        gmsh: lc_min on segment, lc_max away from segment, with a transition
              zone between dist_min and dist_max
        :param ind: list of local index of segments
        :param lc_min: size of element at segment
        :param lc_max: size of element away from segment
        :param dist_min: distance away from segment with lc_min element size
        :param dist_max: distance away from segment with lc_max element size
        """
        var_dict = {'LcMin': lc_min, 'LcMax': lc_max,
                    'DistMin': dist_min, 'DistMax': dist_max}
        self._addConstraint(entity='segment', cons_type='fixed',
                            index=ind, variables=var_dict)

    def refineAroundFacet(self, ind, lc_min, lc_max=None, dist_min=None,
                            dist_max=None):
        """
        Refinement (or coarsening) of mesh around facets.
        tetgen: only lc_min is taken into account
        gmsh: lc_min on segment, lc_max away from facet, with a transition
              zone between dist_min and dist_max
        (!) behaviour can be buggy in gmsh
        :param ind: list of local index of facets
        :param lc_min: size of element at facet
        :param lc_max: size of element away from facet
        :param dist_min: distance away from facet with lc_min element size
        :param dist_max: distance away from facet with lc_max element size
        """
        var_dict = {'LcMin': lc_min, 'LcMax': lc_max,
                    'DistMin': dist_min, 'DistMax': dist_max}
        self._addConstraint(entity='facet', cons_type='fixed',
                            index=ind, variables=var_dict)

    def refineBox(self, lc_in, lc_out, x_min, x_max, y_min, y_max, z_min=None,
                  z_max=None):
        """
        Refinement (or coarsening) of mesh inside a box.
        (!) for gmsh only
        :param lc_in: size of element inside box
        :param lc_out: size of element outside box
        :param x_min: lower limit of x coordinates of box
        :param x_max: upper limit of x coordinates of box
        :param y_min: lower limit of y coordinates of box
        :param y_max: upper limit of y coordinates of box
        :param z_min: lower limit of z coordinates of box
        :param z_max: upper limit of z coordinates of box
        """
        var_dict = {'VIn': lc_in, 'VOut': lc_out, 'XMin': x_min, 'XMax': x_max,
                    'YMin': y_min, 'YMax': y_max, 'ZMin': z_min, 'ZMax': z_max}
        self._addConstraint(entity='global', cons_type='box',
                            index=None, variables=var_dict)

    def setRefinementFunction(self, function):
        """
        Set a function to make the mesh element size vary.
        (!) for gmsh only. Must use MathEval syntax. Can use x, y and z in
            function
        :param function: function that makes mesh vary (string)
        """
        var_dict = {'function': function}
        self._addConstraint(entity='global', cons_type='function',
                            index=None, variables=var_dict)

    def setTransfiniteSegment(self, ind, nb_nodes, prog=1.):
        """
        Sets segment transfinite interpolation. Goes from the segment first
        vertex to its second vertex using the defined progression.
        If TFI should go from segment's second vertex to first vertex,
        use a negative index.
        :param ind: list of local index of segments
        :param nb_nodes: number of nodes on segments
        :param prog: progression parameter for nodes on segments
        """
        var_dict_pos = {'nodes':nb_nodes, 'prog': prog}
        var_dict_neg = {'nodes':nb_nodes, 'prog': 1./prog}
        ind_prog_pos = []
        ind_prog_neg = []
        for ind in indice:
            if ind < 0:
                ind_prog_neg += [abs(ind)]
            else:
                ind_prog_pos += [ind]
        if ind_prog_neg:
            self._addConstraint(entity='segment', cons_type='TFI',
                                index=ind_prog_neg, variables=var_dict_neg)
        if ind_prog_pos:
            self._addConstraint(entity='segment', cons_type='TFI',
                                index=ind_prog_pos, variables=var_dict_pos)

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
    _assembleGeometry(domain, BC_class=bc.BC_Base)
    _generateMesh(domain)

def _assembleGeometry(domain, BC_class):
    mesh = domain.MeshOptions
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
    domain.bc = [BC_class()]
    # domain.bc[0].setNonMaterial()
    # barycenter at flag 0
    domain.barycenters = np.array([[0., 0., 0.]])
    start_flag = 0
    start_vertex = 0
    start_segment = 0
    for shape in domain.shape_list:
        # --------------------------- #
        # ----- DOMAIN GEOMETRY ----- #
        # --------------------------- #
        start_flag = shape.start_flag = len(domain.bc)-1
        start_vertex = shape.start_vertex = len(domain.vertices)
        start_segment = shape.start_segment = len(domain.segments)
        start_facet = shape.start_facet = len(domain.facets)
        start_region = shape.start_region = len(domain.regions)
        if shape.boundaryTags:
            for tag, value in shape.boundaryTags.iteritems():
                domain.boundaryTags[shape.name+'_'+str(tag)] = value+start_flag
                domain.reversed_boundaryTags[value+start_flag] = shape.name+'_'+str(tag)            
            # import operator
            # sorted_tags = sorted(x.items(), key=operator.itemgetter(1))
            # for tag in sorted_tags:
            #     domain.sortedTags = 
            # domain.boundaryTags
        # start_volume = shape.start_volume = len(domain.volumes)
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
        if shape.facets is not None:
            domain.facets += (facets+start_vertex).tolist()
            domain.facetFlags += (shape.facetFlags+start_flag).tolist()
        if shape.regions is not None:
            domain.regions += (shape.regions).tolist()
            domain.regionFlags += (shape.regionFlags+start_rflag).tolist()
        if shape.holes is not None:
            domain.holes += (shape.holes).tolist()
        # adding local shape mesh options to global domain mesh options
        cons = shape.MeshOptions.constraints
        for con in cons:
            domain.MeshOptions.constraints += [con]
            dcon = domain.MeshOptions.constraints[-1]
            if dcon['entity'] == 'vertex':
                dcon['index'] = (np.array(dcon['index'])+start_vertex).tolist()
            if dcon['entity'] == 'segment':
                dcon['index'] = (np.array(dcon['index'])+start_segment).tolist()
            if dcon['entity'] == 'facet':
                print dcon, start_vertex
                dcon['index'] = (np.array(dcon['index'])+start_facet).tolist()
                print dcon
            if dcon['entity'] == 'region':
                dcon['index'] = (np.array(dcon['index'])+start_region).tolist()
        # TFI = np.array(shape.MeshOptions.TFI['segments'])
        # if TFI:
        #     TFI = (TFI+[start_segment, 0, 0]).tolist()
        #     mesh.TFI['segments'] += TFI
        # for line in shape.MeshOptions.TFI['line']:
        #     domain.MeshOptions.TFI['line'] +=
    # for gmsh
    domain.holes_ind = []
    for shape in domain.shape_list:
        shape.start_volume = len(domain.volumes)
        if shape.holes_ind is not None:
            domain.holes_ind += (np.array(shape.holes_ind)+shape.start_volume).tolist()
        if shape.volumes is not None:
            volumes = (shape.volumes+shape.start_facet).tolist()
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
                domain.volumes += volumes


    domain.getBoundingBox()


def _generateMesh(domain):
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


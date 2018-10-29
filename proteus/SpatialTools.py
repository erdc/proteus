"""Module creating predifined or custom shapes. Each shape needs a
Domain as argument (from proteus.Domain). A Domain can contain any
number of shapes.  Boundary conditions objects are automatically
created for each facet (3D) or segment (2D) defining the shape.

Classes:

  * Shape: super class, regroups functions common to all shapes
  * Cuboid: creates a 3D cuboid
  * Rectangle: creates a 2D rectangle
  * Custom: creates a custom shape from a given set vertices, facets, etc.


Example::

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
from __future__ import print_function
from __future__ import division

#from builtins import str
from builtins import range
from past.utils import old_div
from builtins import object
from math import cos, sin, sqrt
import math
import sys
import copy
import numpy as np
from proteus import BoundaryConditions as bc
from .Profiling import logEvent
from subprocess import check_call
from proteus import Comm
from copy import deepcopy
from functools import reduce


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
    count_all = 0

    def __init__(self, domain, nd=None, BC_class=None):
        if nd != domain.nd:
            logEvent('Shape ('+repr(nd)+'D) and Domain ('+repr(domain.nd)+'D)' \
                ' have different dimensions!')
            sys.exit()
        self.Domain = domain
        domain.shape_list.append(self)
        self.__class__.count_all += 1
        self.name = 'shape'+str(self.__class__.count_all)
        self.nd = nd
        self.BC_class = BC_class or bc.BC_Base
        self.vertices = None
        self.vertexFlags = None
        self.vertexFlags_global = None
        self.segments = None
        self.segments_global = None
        self.segmentFlags = None
        self.segmentFlags_global = None
        self.facets = None
        self.facets_global = None
        self.facetFlags = None
        self.facetFlags_global = None
        self.regions = None
        self.volumes = None
        self.regionFlags = None
        self.regionFlags_global = None
        self.holes = None
        self.holes_ind = None
        self.barycenter = np.zeros(3)
        self.coords = None  # Only used for predefined shapes
                            # (can be different from barycenter)
        self.coords_system = np.eye(nd)
        self.boundaryTags = {}
        self.boundaryTags_global = {}
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
        if ind in self.children:
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
        if ind in shape.children:
            shape.children[ind] += [self]
        else:
            shape.children[ind] = [self]

    def setHoles(self, holes, indice=None):
        """
        Sets a 'hole' in the mesh. The region where the hole is defined will
        not be meshed.

        Parameters
        ----------
        holes: array_like
            Array of coordinates of holes (list/array).
        indice: array_like
            Array of index of region where hole is (list/array). Only for gmsh
        """
        self._checkListOfLists(holes)
        self.holes = np.array(holes)
        if indice is not None:
            self.holes_ind = indice

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
        self.volumes = [[[0, 1, 2, 3, 4, 5]]]
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
        if barycenter is None:
            self.barycenter = np.array(coords)
        else:
            self.barycenter = np.array(barycenter)
        self.It = np.array([[old_div((W**2.+H**2.),12.), 0, 0],
                            [0, old_div((L**2.+H**2.),12.), 0],
                            [0, 0, old_div((W**2.+L**2.),12.)]])

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




class Sphere(Shape):
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

    def __init__(self, domain, radius, coords=(0.,0.,0.), barycenter=None,
                 nSectors=10):
        super(Sphere, self).__init__(domain, nd=3)
        self.__class__.count += 1
        self.name = "sphere" + str(self.__class__.count)
        self.radius = radius
        self.coords = np.array(coords)
        self.nSectors = nSectors
        self.constructShape()
        self.volumes = [[[i for i in range(len(self.facets))]]]
        # defining flags for boundary conditions
        self.boundaryTags = bt = {'sphere': 1}
        self.BC = {'sphere': self.BC_class(shape=self, name='sphere')}
        self.BC_list = [self.BC['sphere']]
        # self.BC = BCContainer(self.BC_dict)
        if barycenter is None:
            self.barycenter = np.array(coords)
        else:
            self.barycenter = np.array(barycenter)

    def constructShape(self):
        """
        regualar nx,ny,nz cube packing with padding to boundary
        returns domain size and boundary flags
        """
        coords=(0.,0.,0.)
        n_domain_vertices = 8
        radius = self.radius
        nSectors = self.nSectors
        diameter=2.0*radius
        grain_centers  = {}
        north_poles = {}
        right_poles = {}
        back_poles = {}
        south_poles = {}
        left_poles = {}
        front_poles = {}

        pad= -radius
        noff= 0.
        #cube boundary vertices
        nN=0
        nF=0
        vertices=[]
        facets=[]


        nN = 0
        #cube boundary facets
        # sphereTag = 10
        # boundaries = ['bottom','front','left','right','back','top']
        # boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
        # boundaryTags['sphere']=sphereTag

        nF = 0
        #adding poles of grains to make sure those nodes are unique

        #build sphere
        grain_center = (noff*diameter, noff*diameter, noff*diameter)
        north_pole = nN
        vertices.append([coords[0]+noff*diameter,
                        coords[1]+noff*diameter,
                        coords[2]+noff*diameter+radius])
        nN+=1
        back_pole = nN
        vertices.append([coords[0]+noff*diameter,
                        coords[1]+noff*diameter+radius,
                        coords[2]+noff*diameter])
        nN+=1
        right_pole = nN
        vertices.append([coords[0]+noff*diameter+radius,
                        coords[1]+noff*diameter,
                        coords[2]+noff*diameter])
        nN+=1
        south_pole = nN
        vertices.append([coords[0]+noff*diameter,
                        coords[1]+noff*diameter,
                        coords[2]+noff*diameter-radius])
        nN+=1
        front_pole = nN
        vertices.append([coords[0]+noff*diameter,
                        coords[1]+noff*diameter-radius,
                        coords[2]+noff*diameter])
        nN+=1

        left_pole = nN
        vertices.append([coords[0]+noff*diameter-radius,
                        coords[1]+noff*diameter,
                        coords[2]+noff*diameter])
        nN+=1


        hxi = old_div(radius,(math.sqrt(2.0)*float(nSectors)));
        heta = old_div(radius,(math.sqrt(2.0)*float(nSectors)));
    #now loop over grains
        #top half  sphere nodes
        top_nodes = {}
        for ii in range(2*nSectors+1):
            for jj in range(2*nSectors+1):
                if (ii,jj) == (0,0):
                    top_nodes[(ii,jj)] = left_pole
                elif (ii,jj) == (2*nSectors,0):
                    top_nodes[(ii,jj)] = back_pole
                elif (ii,jj) == (2*nSectors,2*nSectors):
                    top_nodes[(ii,jj)] = right_pole
                elif (ii,jj) == (0,2*nSectors):
                    top_nodes[(ii,jj)] = front_pole
                elif (ii,jj) == (nSectors,nSectors):
                    top_nodes[(ii,jj)] = north_pole
                else:
                    top_nodes[(ii,jj)] = nN
                    #rotate  corners of ref square to line up with poles
                    x0s = (jj-nSectors)*hxi
                    y0s = (ii-nSectors)*heta
                    r0s = math.sqrt(x0s**2 + y0s**2)
                    theta0s = math.atan2(y0s,x0s)
                    theta1s = theta0s - old_div(math.pi,4.0)
                    r1s = r0s
                    x1s = r1s*math.cos(theta1s)
                    y1s = r1s*math.sin(theta1s)
                    #do each quadrant
                    if x1s >= 0.0  and y1s >=0.0:
                        rc = x1s + y1s
                        thetac = theta1s
                    elif x1s < 0.0 and y1s >=0.0:
                        rc = y1s - x1s
                        thetac = theta1s
                    elif x1s <= 0.0 and y1s < 0.0:
                        rc = -(x1s+y1s)
                        thetac = theta1s
                    else:
                        rc = x1s - y1s
                        thetac = theta1s
                    eta = 0.5*math.pi*(radius-rc)/radius
    #                         xc = rc*math.cos(thetac)
    #                         yc = rc*math.sin(thetac)
    #                         zc = math.sqrt(math.fabs(radius**2 - xc**2 - yc**2))
                    xc = radius*math.cos(thetac)*math.cos(eta)
                    yc = radius*math.sin(thetac)*math.cos(eta)
                    zc = radius*math.sin(eta)
                    #zc = math.sqrt(math.fabs(radius**2 - xc**2 - yc**2))
                    #print xc,yc,zc,rc,radius**2 - xc**2 - yc**2
                    #physical coordinates
                    vertices.append([xc+grain_center[0],
                                    yc+grain_center[1],
                                    zc+grain_center[2]])
                    nN+=1
        #bottom half sphere nodes
        bottom_nodes = {}
        for ii in range(2*nSectors+1):
            for jj in range(2*nSectors+1):
                if (ii,jj) == (0,0):
                    bottom_nodes[(ii,jj)] = left_pole
                elif (ii,jj) == (2*nSectors,0):
                    bottom_nodes[(ii,jj)] = back_pole
                elif (ii,jj) == (2*nSectors,2*nSectors):
                    bottom_nodes[(ii,jj)] = right_pole
                elif (ii,jj) == (0,2*nSectors):
                    bottom_nodes[(ii,jj)] = front_pole
                elif (ii,jj) == (nSectors,nSectors):
                    bottom_nodes[(ii,jj)] = south_pole
                elif (ii in [0,2*nSectors] or
                    jj in [0,2*nSectors]):#equator
                    bottom_nodes[(ii,jj)] = top_nodes[(ii,jj)]
                else:
                    bottom_nodes[(ii,jj)] = nN
                    #rotate  corners of ref square to line up with poles
                    x0s = (jj-nSectors)*hxi
                    y0s = (ii-nSectors)*heta
                    r0s = math.sqrt(x0s**2 + y0s**2)
                    theta0s = math.atan2(y0s,x0s)
                    theta1s = theta0s - old_div(math.pi,4.0)
                    r1s = r0s
                    x1s = r1s*math.cos(theta1s)
                    y1s = r1s*math.sin(theta1s)
                    #do each quadrant
                    if x1s >= 0.0  and y1s >=0.0:
                        rc = x1s + y1s
                        thetac = theta1s
                    elif x1s < 0.0 and y1s >=0.0:
                        rc = y1s - x1s
                        thetac = theta1s
                    elif x1s <= 0.0 and y1s < 0.0:
                        rc = -(x1s+y1s)
                        thetac = theta1s
                    else:
                        rc = x1s - y1s
                        thetac = theta1s
                    eta = 0.5*math.pi*(radius-rc)/radius
    #                         xc = rc*math.cos(thetac)
    #                         yc = rc*math.sin(thetac)
    #                         zc = -math.sqrt(math.fabs(radius**2 - xc**2 - yc**2))
                    xc = radius*math.cos(thetac)*math.cos(eta)
                    yc = radius*math.sin(thetac)*math.cos(eta)
                    zc = -radius*math.sin(eta)
                    #print xc,yc,zc
                    #physical coordinates
                    vertices.append([xc+grain_center[0],
                                    yc+grain_center[1],
                                    zc+grain_center[2]])
                    nN+=1
        for ii in range(2*nSectors):
            for jj in range(2*nSectors):
                if ((ii < nSectors and jj < nSectors) or
                    (ii>=nSectors and  jj>=nSectors)):
                    #top
                    facets.append([[ top_nodes[(ii,jj)],top_nodes[(ii+1,jj)],top_nodes[(ii+1,jj+1)]] ])
                    facets.append([[top_nodes[(ii,jj)],top_nodes[(ii+1,jj+1)],top_nodes[(ii,jj+1)]]])
                    nF+=2
                    #bottom
                    facets.append([[ bottom_nodes[(ii,jj)],bottom_nodes[(ii+1,jj)],bottom_nodes[(ii+1,jj+1)]] ])
                    facets.append([[ bottom_nodes[(ii,jj)],bottom_nodes[(ii+1,jj+1)],bottom_nodes[(ii,jj+1)] ]])
                    nF+=2
                else:
                    #top
                    facets.append([[ top_nodes[(ii,jj)],top_nodes[(ii+1,jj)],top_nodes[(ii,jj+1)] ]])
                    facets.append([[ top_nodes[(ii+1,jj)],top_nodes[(ii+1,jj+1)],top_nodes[(ii,jj+1)] ]])
                    nF+=2
                    #bottom
                    facets.append([[ bottom_nodes[(ii,jj)],bottom_nodes[(ii+1,jj)],bottom_nodes[(ii,jj+1)] ]])
                    facets.append([[ bottom_nodes[(ii+1,jj)],bottom_nodes[(ii+1,jj+1)],bottom_nodes[(ii,jj+1)] ]])
                    nF+=2
        self.vertices = np.array(vertices)+self.coords

        self.vertexFlags = np.array([1]*len(vertices))
        self.facets = np.array(facets)
        self.facetFlags = np.array([1]*len(facets))
        self.regions = np.array([self.coords])
        self.regionFlags = np.array([1])


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
        self.It = old_div((L**2+H**2),12)

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

class Cylinder(Shape):
    """
    Class to create a cylinder.

    Parameters
    ----------
    domain: proteus.Domain.D_base
        Domain class instance that hold all the geometrical informations and
        boundary conditions of the shape.
    radius: float
        radius of cylinder.
    height: float
        height of cylinder.
    nPoints: int
        number of points to discretize circles of cylinder.
    coords: Optional[array_like]
        Coordinates of the centroid of the shape.
    barycenter: Optional[array_like]
        Coordinates of the barycenter.
    """
    count = 0
    def __init__(self, domain, radius, height, nPoints, coords=(0.,0.,0.), barycenter=None):
        super(Cylinder, self).__init__(domain, nd=3)
        self.__class__.count += 1
        self.name = "Cylinder" + str(self.__class__.count)
        self.radius = radius
        self.height = height
        self.coords = np.array(coords)
        self.barycenter = np.array(barycenter)
        self.nPoints = nPoints
        self.constructShape()

    def constructShape(self):
        h_offset = np.array([0., 0., self.height])
        arc = 2.*np.pi*self.radius/self.nPoints
        ang = old_div(arc,self.radius)
        vert = []
        facets = []
        segs = []
        for i in range(0, self.nPoints):
            vert += [[self.radius*cos(ang*i),
                      self.radius*sin(ang*i),
                      0]]
            if i > 0:
                segs += [[i-1, i]]
        segs += [[i, 0]]
        segs_bottom = np.array(segs)
        vert_bottom = np.array(vert)
        facets += [[[i for i in range(0, len(vert_bottom))]]]
        vert_top = np.array(vert)+h_offset
        segs_top = np.array(segs)+len(vert) 
        nvb = len(vert_bottom)
        facets += [[[i+nvb for i in range(0, len(vert_top))]]]
        for i in range(len(vert_bottom)-1):
            facets += [[[i, i+1, i+nvb+1, i+nvb]]]
        facets += [[[i+1, 0, nvb, i+1+nvb]]]  # last facet
        self.vertices = np.vstack((vert_bottom, vert_top))-old_div(h_offset,2.)+np.array(self.coords)
        self.segments = np.vstack((segs_bottom, segs_top))
        self.segmentFlags = np.array([1 for i in range(len(segs_bottom))]+[2 for i in range(len(segs_top))])
        self.facets = facets
        self.vertexFlags = np.array([1 for i in range(len(vert_bottom))]+[2 for i in range(len(vert_top))])
        self.facetFlags = np.array([1, 2]+[3 for i in range(len(self.facets)-2)])
        self.boundaryTags = {'z-': 1, 'z+': 2, 'cylinder': 3}
        self.regions = np.array([self.coords])
        self.regionFlags = np.array([1])
        self.volumes = [[[i for i in range(len(self.facets))]]]
        self.BC = {'z-': self.BC_class(shape=self, name='z-'),
                   'z+': self.BC_class(shape=self, name='z+'),
                   'cylinder': self.BC_class(shape=self, name='cylinder')}
        self.BC_list = [self.BC['z-'], self.BC['z+'], self.BC['cylinder']]


class Circle(Shape):
    """
    Class to create a circular shape.

    Parameters
    ----------
    domain: proteus.Domain.D_base
        Domain class instance that hold all the geometrical informations and
        boundary conditions of the shape.
    radius:
        Radius of the circular shape.
    coords: Optional[array_like]
        Coordinates of the centroid of the shape.
    """

    count = 0

    def __init__(self, domain, radius, coords, barycenter, nPoints):
        super(Circle, self).__init__(domain, nd=2)
        self.__class__.count += 1
        self.name = "Circle" + str(self.__class__.count)

        pi=np.pi
        self.radius = radius
        self.coords = xc, yc = np.array(coords)
        self.nPoints = nPoints
        self.arc=2.0*pi*self.radius/self.nPoints
        self.ang=old_div(self.arc,self.radius)
        verti=[]
        segm=[]

        # points
        for iii in range(0,nPoints):
            vert=np.array([xc+self.radius*cos(self.ang*iii),
                           yc+self.radius*sin(self.ang*iii),
                           ])
            verti.append(vert,)
        self.vertices=np.array(verti)
        nVertices=len(self.vertices)

        # segments
        for jjj in range(1,nVertices):
            seg=[jjj-1, jjj]
            segm.append(seg,)
        seg_last=[nVertices-1,0]
        segm.append(seg_last,)
        self.segments=np.array(segm)
        
        # facets
        facets=[]
        for kkk in range(nVertices):
            facets.append(kkk)
        self.facets=np.array([[facets]])

        # barycenter, orientations, region
        if barycenter is not None:
            self.barycenter[0:2] = barycenter[0:2]
        else:
            self.barycenter[0:2] = coords[0:2]

        self.b_or = np.array([[1., 0.],
                              ])
        self.regions = np.array([[xc, yc]])

        # boundary tags
        self.boundaryTags = bt = {'circle':1}
        vertFlag=[]
        segmFlag=[]
        for kkk in range(nVertices):
                vertFlag.append(bt['circle'],)
                segmFlag.append(bt['circle'],)
        self.vertexFlags = np.array(vertFlag)
        self.segmentFlags = np.array(segmFlag)
        self.facetFlags=np.array([1])
        self.regionFlags = np.array([1])

        # boundary list
        self.BC = {'circle': self.BC_class(shape=self, name='circle',
                                       b_or=self.b_or, b_i=0),
                  }
        self.BC_list = [self.BC['circle'],
                       ]

        # set the inertia
        self.It = pi*(radius**4)/2.



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
        Array of holes coordinates (unmeshed regions; for triangle/tetgen only)
    holes: array_like
        Array of holes index (volume or facet; unmeshed regions; for gmsh only)
    boundaryTags: dict
        Dictionary of flags (int) as keys, and tags (e.g. string) for BC.
    boundaryOrientations: Optional[array_like]
        Array of orientation of boundaries. Can be used for BC.
    """
    count = 0

    def __init__(self, domain, barycenter=None, vertices=None,
                 vertexFlags=None, segments=None, segmentFlags=None,
                 facets=None, facetFlags=None, holes=None, holes_ind=None,
                 regions=None, regionFlags=None, volumes=None,
                 boundaryTags=None, boundaryOrientations=None):
        super(CustomShape, self).__init__(domain, nd=len(vertices[0]))
        self.__class__.count += 1
        self.name = "custom" + str(self.__class__.count)
        self._checkFlags(list(boundaryTags.values()))
        self.boundaryTags = boundaryTags
        self.vertices = np.array(vertices)
        self.vertexFlags = np.array(vertexFlags)
        if segments is not None:
            self.segments = np.array(segments)
            self.segmentFlags = np.array(segmentFlags)
        if facets is not None:
            self.facets = facets
            self.facetFlags = np.array(facetFlags)
        if holes is not None:
            self.holes = np.array(holes)
        if holes_ind is not None:
            self.holes_ind = np.array(holes_ind)
        if volumes is not None:
            self.volumes = volumes
        if regions is not None:
            self._checkFlags(regionFlags)
            self.regions = np.array(regions)
            self.regionFlags = np.array(regionFlags)
        self.BC = {}
        self.BC_list = [None]*len(boundaryTags)
        b_or = [None]*len(boundaryTags)
        self.b_or = b_or
        for tag, flag in boundaryTags.items():
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
    count = 0

    def __init__(self, domain, filename):
        super(ShapeSTL, self).__init__(domain, nd=3)
        self.__class__.count += 1
        self.name = 'STL'+str(self.__class__.count)
        self.filename = filename
        self.vertices, self.facets, self.facetnormals = getInfoFromSTL(self.filename)
        self.facetFlags = np.ones(len(self.facets))
        self.vertexFlags = np.ones(len(self.vertices))
        self.volumes = [[[i for i in range(len(self.facets))]]]
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
    axis = old_div(axis,r)
    # get values for rotation matrix
    cx, cy, cz = axis
    d = sqrt(cy**2+cz**2)
    # rotation matrices
    if d != 0:
        Rx = np.array([[1,         0,        0,    0],
                       [0,         old_div(cz,d),     old_div(cy,d), 0],
                       [0,         old_div(-cy,d),    old_div(cz,d), 0],
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
    domain.useSpatialTools = True
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
            shape.boundaryTags_global = {}
            for tag, value in shape.boundaryTags.items():
                new_tag = shape.name+'_'+str(tag)
                new_flag = value+start_flag
                domain.boundaryTags[new_tag] = new_flag
                shape.boundaryTags_global[tag] = new_flag
                domain.reversed_boundaryTags[new_flag] = new_tag
                domain.BCbyFlag[new_flag] = shape.BC[tag]
                domain.BC[new_tag] = shape.BC[tag]
        if domain.regionFlags:
            start_rflag = max(domain.regionFlags)
        else:
            start_rflag = 0
        domain.bc += shape.BC_list
        # making copies of shape properties before operations/modifications
        vertices = copy.deepcopy(shape.vertices)
        vertexFlags = copy.deepcopy(shape.vertexFlags)
        if shape.segments is not None:
            segments = copy.deepcopy(shape.segments)
        if shape.facets is not None:
            facets = copy.deepcopy(shape.facets)
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
                    for facet in shape.facets:
                        for subfacet in facet:
                            for ind, ver in enumerate(subfacet):
                                if ver > i_s:
                                    subfacet[ind] -= 1
                                elif ver == i_s:
                                    subfacet[ind] = i_d-start_vertex
        # adding shape geometry to domain
        domain.vertices += vertices.tolist()
        domain.vertexFlags += (vertexFlags+start_flag).tolist()
        shape.vertexFlags_global = vertexFlags+start_flag
        barycenters = np.array([shape.barycenter for bco in shape.BC_list])
        domain.barycenters = np.append(domain.barycenters, barycenters, axis=0)
        if shape.segments is not None:
            domain.segments += (segments+start_vertex).tolist()
            shape.segments_global = segments+start_vertex
            domain.segmentFlags += (shape.segmentFlags+start_flag).tolist()
            shape.segmentFlags_global = segments+start_flag
        if shape.regions is not None:
            domain.regions += (shape.regions).tolist()
            domain.regionFlags += (shape.regionFlags+start_rflag).tolist()
            shape.regionFlags_global = shape.regionFlags+start_flag
        if shape.facets is not None:
            if type(facets) is np.ndarray:
                facets = facets.tolist()
            else:
                facets = copy.deepcopy(shape.facets)
            for i, facet in enumerate(facets):
                for j, subf in enumerate(facet):
                    for k, v_nb in enumerate(subf):
                        facets[i][j][k] = v_nb+shape.start_vertex
            domain.facets += facets
            shape.facets_global = facets
            if shape.nd == 2:  # facet flags are actually region flags if 2D
                domain.facetFlags += (shape.regionFlags+start_rflag).tolist()
                shape.facetFlags_global = shape.regionFlags+start_rflag
            elif shape.nd == 3:
                domain.facetFlags += (shape.facetFlags+start_flag).tolist()
                shape.facetFlags_global = shape.facetFlags+start_flag
        if shape.holes is not None:
            domain.holes += (shape.holes).tolist()
    
    # 2D holes (only for gmsh)
    domain.holes_ind = []
    if domain.nd == 2 and shape.facets is not None:
        domain.facets = []
        for shape in domain.shape_list:
            if shape.holes_ind is not None:
                    domain.holes_ind += (np.array(shape.holes_ind)+shape.start_facet).tolist()
            if shape.facets is not None:
                facets = deepcopy(shape.facets)
                if type(facets) is np.ndarray:
                    facets = facets.tolist()
                for i, facet in enumerate(facets):
                    for j, subf in enumerate(facet):
                        for k, v_nb in enumerate(subf):
                            facets[i][j][k] = v_nb+shape.start_vertex
                # facets = (shape.facets+shape.start_vertex).tolist()
                f_to_remove = []
                f_to_add = []
                facets_again = deepcopy(shape.facets)
                if type(facets_again) is np.ndarray:
                    facets_again = facets_again.tolist()
                for i, facet in enumerate(facets_again):
                    if i in list(shape.children.keys()):
                        for child in shape.children[i]:
                            child_seg = {}
                            child_facs = deepcopy(child.facets)
                            if type(child_facs) is np.ndarray:
                                child_facs = child_facs.tolist()
                            for child_fac in child_facs:
                                fac = child_fac[0]
                                for j, v in enumerate(fac):
                                    if (fac[j-1], fac[j]) in child_seg:
                                        fn = child_seg[(fac[j-1], fac[j])]
                                        f_to_merge = child.facets[fn][0]
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
                                        f_to_merge = child.facets[fn][0]
                                        # shift lists
                                        f_to_merge = list(deque(f_to_merge).rotate(-f_to_merge.index(fac[j-1])))
                                        fac = list(deque(fac).rotate(-fac.index(fac[j])))
                                        new_f = fac[:-1]+f_to_merge[:-1]
                                        new_f = fac[:fac.index(fac[j-1])]+f_to_merge
                                        f_to_remove += [fn, j]
                                        f_to_add += [new_f]
                                    else:
                                        child_seg[(fac[j-1], fac[j])] = j
                            child_facets = deepcopy(child.facets)
                            if type(child_facets) is np.ndarray:
                                child_facets = child_facets.tolist()
                            facets2 = [f[0] for k, f in enumerate(child_facets) if k not in f_to_remove]
                            if f_to_add:
                                facets2 += [f_to_add]
                        new_facets = facets2
                        for facet_i in range(len(new_facets)):
                            facet = new_facets[facet_i]
                            for v_i in range(len(facet)):
                                print(new_facets[facet_i][v_i])
                                new_facets[facet_i][v_i] += child.start_vertex
                        facets[i] += new_facets
                domain.facets += facets

    # 3D holes (only for gmsh)
    elif domain.nd == 3 and shape.volumes is not None:
        for shape in domain.shape_list:
            shape.start_volume = len(domain.volumes)
            if shape.holes_ind is not None:
                domain.holes_ind += (np.array(shape.holes_ind)+shape.start_volume).tolist()
            if shape.volumes is not None:
                volumes = copy.deepcopy(shape.volumes)
                for volume in volumes:
                    for subvolume in volume:
                        for i in range(len(subvolume)):
                            subvolume[i] += shape.start_facet
                volumes_to_remove = []
                for i, volume in enumerate(volumes):
                    # add holes to volumes
                    if i in shape.children:
                        for child in shape.children[i]:
                            child_vols = []
                            child_facets = []  # facets in volumes
                            for child_vol in child.volumes:
                                # don't need holes of child volumes, only outer shells:
                                child_vol_copy = copy.deepcopy(child_vol[0])
                                for i in range(len(child_vol_copy)):
                                    child_vol_copy[i] += child.start_facet
                                child_vols += [child_vol_copy]
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
    from proteus import Comm
    comm = Comm.get()
    mesh = domain.MeshOptions
    if comm.isMaster():
        if mesh.outputFiles['poly'] is True:
            domain.writePoly(mesh.outputFiles_name)
        if mesh.outputFiles['ply'] is True:
            domain.writePLY(mesh.outputFiles_name)
        if mesh.outputFiles['asymptote'] is True:
            domain.writeAsymptote(mesh.outputFiles_name)
        if mesh.outputFiles['geo'] is True or mesh.use_gmsh is True:
            domain.writeGeo(mesh.outputFiles_name)
    else:
        domain.polyfile=mesh.outputFiles_name
    mesh.setTriangleOptions()


def getGmshPhysicalGroups(geofile):
    boundaryTags = {}
    import re
    # gmsh_cmd = "gmsh {0:s} -0".format(geofile)
    # check_call(gmsh_cmd, shell=True)
    with open(geofile, 'r') as f:
        for line in f:
            if line.startswith("Physical "):
                words = line.lstrip().split('=', 1)
                tagflag = re.search(r'\((.*?)\)', words[0]).group(1).split(',')
                if len(tagflag) == 2:
                    tag = str(tagflag[0][2:-2])
                    boundaryTags[tag] = None  # add empty BC holder
                    flag = int(tagflag[1])
                    print(tagflag)
                else:
                    try:
                        flag = tag = int(tagflag[0])
                    except:
                        flag = tag = tagflag[0]
                boundaryTags[tag] = flag  # add empty BC holder
    return boundaryTags

def extrude2Dto3D(extrusion, vertices, segments, facets, regions=None):
    """
    Extrude 2D geometry attributes in 3D

    Parameters
    ----------
    extrusion: array_like
        value of the extrusion (e.g. [1.4, 0., 0.])
    vertices: array_like
        list of vertices of the 2D shape
    segments: array_like
        list of segments of the 2D shape
    facets: array_like
        list of facets of the 2D shape
    regions: array_like
        list of regions of the 2D shape
    """
    nv = len(vertices)
    ns = len(segments)
    vertices = np.array(vertices)
    vertices = np.vstack((vertices, vertices+extrusion))
    segments_side = [[i, i+nv] for i in range(nv)]
    segments_extruded = copy.deepcopy(segments)
    for segment in segments_extruded:
        segment[0] += ns
        segment[1] += ns
    segments = np.vstack((segments, segments_side, segments_extruded))
    facets_side = [[[i, i+1, i+nv+1, i+nv]] for i in range(nv)]
    facets_side[-1][0][1] = 0
    facets_side[-1][0][2] = nv
    facets_extruded = copy.deepcopy(facets)
    for facet in facets_extruded:
        for subfacet in facet:
            for i in range(len(subfacet)):
                subfacet[i] += nv
    facets = facets+facets_side+facets_extruded
    regions_extruded = copy.deepcopy(regions)
    if regions is not None:
        for region in regions_extruded:
            region += old_div(extrusion,2.)
    return vertices, segments, facets, regions_extruded

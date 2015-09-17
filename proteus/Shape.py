"""
Module creating predifined shapes associated. Each shape needs a domain as argument and the wanted dimensions.
Boundary conditions objects are automatically created for each facet (3D) or segment (2D) defining the shape.
classes:
- Shape: super class, regroups functions common to all shapes
- Cuboid: creates a cuboid 
- Rectangle: creates a rectangle

TO DO:

Add more predifined shape types
"""

import BC as bc
import numpy as np
from math import *

class Shape:
    """
    Class defining a shape
    :param domain: domain in which the shape is defined
    """
    def __init__(self, domain):
        self.domain = domain
        self.vertices = None
        self.vertexFlags = None
        self.segments = None
        self.segmentFlags = None
        self.facets = None
        self.facetFlags = None
        self.regions = None
        self.regionFlags = None
        self.holes = None
        self.volume = None
        self.mass = None
        self.density = None
        self.barycenter = None
        self.cm = None  # centre of mass
        self.bc = None
        self.snv = len(self.domain.vertices)  # total number of vertices in domain when shape.__init__
        self.snf = len(self.domain.facets)
        self.sns = len(self.domain.segments)
        self.snr = len(self.domain.regions)
        self.snbc = len(self.domain.bc)

    def addShape(self):
        """
        Adds Shape information to the domain
        """
        # add new information to the domain
        self.domain.vertices += self.vertices.tolist()
        self.domain.vertexFlags += (self.vertexFlags+self.snv).tolist()
        self.domain.segments += (self.segments+self.snv).tolist()
        self.domain.segmentFlags += (self.segmentFlags+self.sns).tolist()
        if self.domain.nd == 3:
            self.domain.facets += (self.facets+self.snv).tolist()
            self.domain.facetFlags += (self.facetFlags+self.snf).tolist()
        self.domain.holes += self.holes.tolist()
        self.domain.regions += self.regions.tolist()
        self.domain.regionFlags += (self.regionFlags+self.snr).tolist()
        self.domain.bc += self.bc
        self.domain.update()

    def updateDomain(self):
        """
        Updates domain when vertex and region coords have been changed
        """
        nv = len(self.vertices)
        nr = len(self.regions)
        self.domain.vertices[self.snv:self.snv+nv] = self.vertices.tolist()
        self.domain.regions[self.snr:self.snr+nr] = self.regions.tolist()
        self.domain.update()

    def setPosition(self, coords):
        """
        Set position of the Shape
        :arg coords: new set of coordinates for the Shape
        """
        old_coords = self.coords
        self.vertices += coords - old_coords
        self.regions += coords
        self.updateDomain()

    def setDimensions(self, dim):
        """
        Set dimensions of the shape
        :arg dim: new dimensions of the Shape
        """
        self.vertices = self.dimfactor*dim + self.coords 
        
    def setConstraints(self, free_x, free_r):
        """
        Sets constraints on the Shape
        :arg free_x: translational constraints
        :arg free_r: rotational constraints
        """
        self.free_x = free_x
        self.free_r = free_r

    def rotate(self, rot, axis=(0,0,1), pivot=None):
        """
        Function to rotate Shape
        :arg rot: angle of rotation in radians (float)
        :arg axis: axis of rotation (list or array)
        :arg pivot: point around which the Shape rotates
        """
        pivot = pivot or self.coords
        if self.domain.nd == 2:
            self.vertices = rotation2D(self.vertices, rot, pivot)
            self.regions = rotation2D(self.regions, rot, pivot)
        elif self.domain.nd == 3:
            self.vertices = rotation3D(self.vertices, rot, axis, pivot)
            self.holes = rotation3D(self.holes, rot, axis, pivot)
            self.regions = rotation3D(self.regions, rot, axis, pivot)
        self.updateDomain()

    def translate(self, trans):
        """
        Function to translate Shape
        :arg trans: translation values
        """
        self.vertices += trans
        self.updateDomain()

    def setMass(self, mass):
        """
        Set mass of the shape and calculate density
        :arg mass: mass of the Shape
        """
        self.mass = mass
        if self.volume:
            self.density = self.mass/self.volume

    def setDensity(self, density):
        """
        Set density of the Shape and calculate mass
        :arg density: density of the Shape
        """
        self.density = density
        if self.volume:
            self.mass = self.density*self.volume


class Cuboid(Shape):
    """
    Class to create a cuboid
    :arg domain: domain of the cuboid
    :arg dim: dimensions of the cuboid (list or array)
    :arg coords: coordinates of the cuboid (list or array)
    """
    def __init__(self, domain, dim=(0.,0.,0.), coords=(0.,0.,0.)):
        Shape.__init__(self, domain)
        self.dim = list(dim)  # length, width height
        self.volume = dim[0]*dim[1]*dim[2]
        self.coords = list(coords)
        self.coords = x, y, z = coords
        self.dimfactor = np.array([[-0.5, -0.5, -0.5],
                                   [-0.5, +0.5, -0.5],
                                   [+0.5, +0.5, -0.5],
                                   [+0.5, -0.5, -0.5],
                                   [-0.5, -0.5, +0.5],
                                   [-0.5, +0.5, +0.5],
                                   [+0.5, +0.5, +0.5],
                                   [+0.5, -0.5, +0.5]])
        self.vertices = np.array([[x-0.5*dim[0], y-0.5*dim[1], z-0.5*dim[2]],
                                  [x-0.5*dim[0], y+0.5*dim[1], z-0.5*dim[2]],
                                  [x+0.5*dim[0], y+0.5*dim[1], z-0.5*dim[2]],
                                  [x+0.5*dim[0], y-0.5*dim[1], z-0.5*dim[2]],
                                  [x-0.5*dim[0], y-0.5*dim[1], z+0.5*dim[2]],
                                  [x-0.5*dim[0], y+0.5*dim[1], z+0.5*dim[2]],
                                  [x+0.5*dim[0], y+0.5*dim[1], z+0.5*dim[2]],
                                  [x+0.5*dim[0], y-0.5*dim[1], z+0.5*dim[2]]])
        self.segments = np.array([[0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6],
                         [6, 7], [7, 4], [0, 4], [1, 5], [2, 6], [3, 7]])
        self.facets = np.array([[0, 1, 2, 3],  # bottom
                       [0, 1, 5, 4],  # front
                       [1, 2, 6, 5],  # right
                       [2, 3, 7, 6],  # back
                       [3, 0, 4, 7],  # left
                       [4, 5, 6, 7]])  # top
        self.regions = np.array([[x, y, z]])
        self.holes = np.array([coords])
        # defining flags for boundary conditions
        self.facetFlags = np.array([0, 1, 2, 3, 4, 5])  # bottom, front, right, back, left, top
        self.vertexFlags = np.array([0, 0, 0, 0, 5, 5, 5, 5])  # only top and bottom for vertices
        self.segmentFlags = np.array([])
        self.regionFlags = np.array([0])
        # Initialize (empty) boundary conditions
        self.bc_bottom = bc.BoundaryConditions()
        self.bc_front = bc.BoundaryConditions()
        self.bc_right = bc.BoundaryConditions()
        self.bc_back = bc.BoundaryConditions()
        self.bc_left = bc.BoundaryConditions()
        self.bc_top = bc.BoundaryConditions()
        self.bc = [self.bc_bottom, self.bc_front, self.bc_right,
                    self.bc_back, self.bc_left, self.bc_top]
        self.addShape()  # adding shape to domain


class Rectangle(Shape):
    """
    Class to create a rectangle
    :arg domain: domain of the rectangle
    :arg dim: dimensions of the rectangle (list or array)
    :arg coords: coordinates of the rectangle (list or array)
    """
    def __init__(self, domain, dim=(0.,0.), coords=(0.,0.)):
        Shape.__init__(self, domain)
        self.dim = list(dim)
        self.coords = x, y = list(coords)
        self.dimfactor = np.array([[-0.5, -0.5],
                                   [+0.5, -0.5],
                                   [+0.5, +0.5],
                                   [-0.5, +0.5]])
        self.vertices = np.array([[x-0.5*dim[0], y-0.5*dim[1]],
                                  [x+0.5*dim[0], y-0.5*dim[1]],
                                  [x+0.5*dim[0], y+0.5*dim[1]],
                                  [x-0.5*dim[0], y+0.5*dim[1]]])
        self.segments= np.array([[0,1],[1,2],[2,3],[3,0]])
        self.regions = np.array([[x, y]])
        self.segmentFlags = np.array([0, 1, 2, 3]) # bottom, right, top, left
        self.vertexFlags = np.array([0, 0, 2, 2])  # bottom, bottom, top, top
        self.regionFlags = np.array([0])
        # Initialize (empty) boundary conditions
        self.bc_bottom = bc.BoundaryConditions()
        self.bc_right = bc.BoundaryConditions()
        self.bc_left = bc.BoundaryConditions()
        self.bc_top = bc.BoundaryConditions()
        self.bc = [self.bc_bottom, self.bc_right, self.bc_left, self.bc_top]
        self.addShapes()  # adding shape to domain



# ------------------------------------------------------------------------------ #
# --------------------------SPATIAL TOOLS FOR SHAPES---------------------------- #
# ------------------------------------------------------------------------------ #

def rotation2D(points, rot, pivot=(0.,0.)):
    """
    function to make a set of points/vertices/vectors (arg: points) to rotate 
    around a pivot point (arg: pivot) 
    :arg points: set of 3D points (array)
    :arg rot: angle of rotation (in radians) 
    :arg pivot: point around which the set of points rotates (list or array)
    :return points_rot: the rotated set of points
    """
    rot = float(rot)
    # get coordinates for translation
    x, y = pivot
    # translation matrix
    T = np.array([[1,  0,   0],
                  [0,  1,   0],
                  [-x, -y,  1]])
    # rotation matrices
    R = np.array([[cos(rot), -sin(rot), 0],
                  [sin(rot), cos(rot),  0],
                  [0,        0,         1]])
    # full transformation matrix
    M = reduce(np.dot, [T, R, np.linalg.inv(T)])
    # transform points (check also if it is only a 1D array or 2D)
    points_rot = np.ones((len(points),3))
    points_rot[:,:-1] = points
    points_rot = points_rot.T  # set of vectors (transposing array)
    points_rot = np.dot(M, points_rot)  # matrix dot product on each vector
    points_rot = points_rot.T  # transpose back to original array shape
    points_rot = points_rot[:,:2]
    return points_rot

def rotation3D(points, rot, axis=(0.,0.,1.), pivot=(0.,0.,0.)):
    """
    function to make a set of points/vertices/vectors (arg: points) to rotate 
    around an arbitrary axis/vector (arg: axis) going through a pivot point (arg: pivot) 
    :arg points: set of 3D points (array)
    :arg rot: angle of rotation (in radians) 
    :arg axis: axis of rotation (list or array)
    :arg pivot: point around which the set of points rotates (list or array)
    :return points_rot: the rotated set of points
    """
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
                       [0,         cz/d,    -cy/d, 0],
                       [0,         cy/d,    cz/d,  0],
                       [0,         0,        0,    1]])
    else:  # special case: rotation axis aligned with x axis    
        Rx = np.array([[1,         0,        0,    0],
                       [0,         1,        0,    0],
                       [0,         0,        1,    0],
                       [0,         0,        0,    1]])
    Ry = np.array([[d,        0,         -cx, 0],
                   [0,        1,         0,   0],
                   [cx,       0,         d,   0],
                   [0,        0,         0,   1]])
    Rz = np.array([[cos(rot), -sin(rot), 0,   0],
                   [sin(rot), cos(rot),  0,   0],
                   [0,        0,         1,   0],
                   [0,        0,         0,   1]])
    # translation matrix
    T = np.array([[1,  0,  0,  -x],
                  [0,  1,  0,  -y],
                  [0,  0,  1,  -z],
                  [0,  0,  0,  1]])
    # full transformation matrix
    M = reduce(np.dot, [T, Rx, Ry, Rz, np.linalg.inv(Ry), np.linalg.inv(Rx), np.linalg.inv(T)])
    points_rot = np.ones((len(points),4))
    points_rot[:,:-1] = points
    points_rot = points_rot.T  # set of vectors (transposing array)
    points_rot = np.dot(M, points_rot)  # matrix dot product on each vector
    points_rot = points_rot.T  # transpose back to original array shape
    points_rot = points_rot[:,:3]
    return points_rot

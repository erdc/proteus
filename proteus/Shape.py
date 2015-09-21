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
    def __init__(self, domain, dim=None, coords=None):
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
        self.barycentre = None
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
        old_coords = np.array(self.coords)
        trans = coords - old_coords
        self.translate(trans)

    def setDimensions(self, dim):
        """
        Set dimensions of the shape
        :arg dim: new dimensions of the Shape
        """
        self.vertices = self.dimfactor*dim + self.coords
        self.getInertiaTensor()
        self.updateDomain()
        
    def setBarycentre(self, barycentre):
        """
        Set barycentre of the shape, from the reference frame of the shape
        :arg barycentre: coordinates of barycentre from centre of shape
        """
        self.barycentre = np.einsum('ij,i->j', self.coords_system, barycentre)

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
        if pivot is None:
            pivot = self.barycentre
            rotate_barycentre = False
        else:
            rotate_barycentre = True
        if self.domain.nd == 2:
            self.vertices = rotation2D(points=self.vertices, rot=rot, pivot=pivot)
            self.regions = rotation2D(points=self.regions, rot=rot, pivot=pivot)
            self.coords_system = rotation2D(points=self.coords_system, rot=rot, pivot=(0.,0.,0.))
            if rotate_barycentre:
                self.barycentre = rotation2D(points=self.barycentre, rot=rot, pivot=pivot)
        elif self.domain.nd == 3:
            self.vertices = rotation3D(points=self.vertices, rot=rot, axis=axis, pivot=pivot)
            self.holes = rotation3D(points=self.holes, rot=rot, axis=axis, pivot=pivot)
            self.regions = rotation3D(points=self.regions, rot=rot, axis=axis, pivot=pivot)
            self.coords_system = rotation3D(points=self.coords_system, rot=rot, axis=axis, pivot=(0.,0.,0.))
            if rotate_barycentre:
                self.barycentre = rotation3D(points=self.barycentre, rot=rot, axis=axis, pivot=pivot)
        self.updateDomain()

    def translate(self, trans):
        """
        Function to translate Shape
        :arg trans: translation values
        """
        self.vertices += trans
        self.regions += trans
        self.coords += trans
        self.barycentre += trans
        if self.domain.nd == 3:
            self.holes += trans
        self.updateDomain()

    def setMass(self, mass):
        """
        Set mass of the shape and calculate density
        :arg mass: mass of the Shape
        """
        self.mass = float(mass)
        if self.volume:
            self.density = self.mass/self.volume

    def setDensity(self, density):
        """
        Set density of the Shape and calculate mass
        :arg density: density of the Shape
        """
        self.density = float(density)
        if self.volume:
            self.mass = self.density*self.volume

    def getInertia(self, vec=(0.,0.,1.), pivot=None):
        if self.domain.nd == 2:
            I = self.It*self.mass
        elif self.domain.nd == 3:
            # Pivot coords relative to shape centre of mass
            pivot = pivot or self.barycentre
            pivot = pivot-np.array(self.barycentre)
            # making unity vector/axis of rotation
            vec = vx, vy, vz = np.array(vec)
            length_vec = sqrt(vx**2+vy**2+vz**2)
            vec = vec/length_vec
            # vector relative to original position of shape:
            vec = relative_vec(vec, self.coords_system[2])
            cx, cy, cz = vec
            # getting the tensor for calculaing moment of inertia from arbitrary axis 
            vt = np.array([[cx**2, cx*cy, cx*cz],
                           [cx*cy, cy**2, cy*cz],
                           [cx*cz, cy*cz, cz**2]])
            # total moment of inertia
            I = np.einsum('ij,ij->', self.mass*self.It, vt)
        return I


class Cuboid(Shape):
    """
    Class to create a cuboid
    :arg domain: domain of the cuboid
    :arg dim: dimensions of the cuboid (list or array)
    :arg coords: coordinates of the cuboid (list or array)
    """
    def __init__(self, domain, dim=(0.,0.,0.), coords=(0.,0.,0.), barycentre=None):
        Shape.__init__(self, domain)
        self.dim = L, W, H = dim  # length, width height
        self.It = np.array([[(W**2.+H**2.)/12., 0, 0],
                            [0, (L**2.+H**2.)/12., 0],
                            [0, 0, (W**2.+L**2.)/12.]])
        self.volume = dim[0]*dim[1]*dim[2]
        self.coords = x, y, z = coords
        self.coords_system = np.eye(3)
        self.barycentre = barycentre or coords
        self.dimfactor = np.array([[-0.5, -0.5, -0.5],
                                   [-0.5, +0.5, -0.5],
                                   [+0.5, +0.5, -0.5],
                                   [+0.5, -0.5, -0.5],
                                   [-0.5, -0.5, +0.5],
                                   [-0.5, +0.5, +0.5],
                                   [+0.5, +0.5, +0.5],
                                   [+0.5, -0.5, +0.5]])
        self.vertices = np.array([[x-0.5*L, y-0.5*W, z-0.5*H],
                                  [x-0.5*L, y+0.5*W, z-0.5*H],
                                  [x+0.5*L, y+0.5*W, z-0.5*H],
                                  [x+0.5*L, y-0.5*W, z-0.5*H],
                                  [x-0.5*L, y-0.5*W, z+0.5*H],
                                  [x-0.5*L, y+0.5*W, z+0.5*H],
                                  [x+0.5*L, y+0.5*W, z+0.5*H],
                                  [x+0.5*L, y-0.5*W, z+0.5*H]])
        self.segments = np.array([[0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6],
                                  [6, 7], [7, 4], [0, 4], [1, 5], [2, 6], [3, 7]])
        self.facets = np.array([[[0, 1, 2, 3]],  # bottom
                                [[0, 1, 5, 4]],  # front
                                [[1, 2, 6, 5]],  # right
                                [[2, 3, 7, 6]],  # back
                                [[3, 0, 4, 7]],  # left
                                [[4, 5, 6, 7]]])  # top
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

    def getInertiaTensor(self):
        L, W, H = self.dim
        self.It = np.array([[(W**2.+H**2.)/12., 0, 0],
                            [0, (L**2.+H**2.)/12., 0],
                            [0, 0, (W**2.+L**2.)/12.]])
            

class Rectangle(Shape):
    """
    Class to create a rectangle
    :arg domain: domain of the rectangle
    :arg dim: dimensions of the rectangle (list or array)
    :arg coords: coordinates of the rectangle (list or array)
    """
    def __init__(self, domain, dim=(0.,0.), coords=(0.,0.), barycentre=None):
        Shape.__init__(self, domain)
        self.dim = L, H = dim  # length, height
        self.It = (L**2+H**2)/12
        self.coords = x, y = coords
        self.coords_system = np.eye(2)
        self.barycentre = barycentre or coords
        self.dimfactor = np.array([[-0.5, -0.5],
                                   [+0.5, -0.5],
                                   [+0.5, +0.5],
                                   [-0.5, +0.5]])
        self.vertices = np.array([[x-0.5*L, y-0.5*H],
                                  [x+0.5*L, y-0.5*H],
                                  [x+0.5*L, y+0.5*H],
                                  [x-0.5*L, y+0.5*H]])
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

    def getInertiaTensor(self):
        L, H = self.dim
        self.It = np.array([[(L**2.)/3., 0],
                            [0, (H**2.)/3.]])


# ------------------------------------------------------------------------------ #
# --------------------------SPATIAL TOOLS FOR SHAPES---------------------------- #
# ------------------------------------------------------------------------------ #

def rotation2D(points, rot, pivot=(0.,0.)):
    """
    function to make a set of points/vertices/vectors (arg: points) to rotate 
    around a pivot point (arg: pivot) 
    :arg points: set of 3D points (list or array)
    :arg rot: angle of rotation (in radians) 
    :arg pivot: point around which the set of points rotates (list or array)
    :return points_rot: the rotated set of points
    """
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
        points_rot = np.ones((len(points),3))
        points_rot[:,:-1] = points
        points_rot = np.dot(points_rot, M)  # matrix dot product on each vector
        points_rot = points_rot[:,:-1]
    else:
        points_rot = np.ones(3)
        points_rot[:-1] = points
        points_rot = np.dot(points_rot, M)  # matrix dot product on each vector
        points_rot = points_rot[:-1]
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
    print axis
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
    print T
    # full transformation matrix
    M = reduce(np.dot, [T, Rx, Ry, Rz, np.linalg.inv(Ry), np.linalg.inv(Rx), np.linalg.inv(T)])
    print M
    if points.ndim > 1:
        points_rot = np.ones((len(points),4))
        points_rot[:,:-1] = points
        points_rot = np.dot(points_rot, M)  # matrix dot product on each vector
        points_rot = points_rot[:,:-1]
    else:
        points_rot = np.ones(4)
        points_rot[:-1] = points
        points_rot = np.dot(points_rot, M)  # matrix dot product on each vector
        points_rot = points_rot[:-1]
    return points_rot


def relative_vec(vec1, vec0):
    """
    function giving coordinates of a vector relative to another vector
    (projecting vec0 as the z-axis for vec1)
    :arg vec1: vector to get new coordinates
    :arg vec0: vector of reference
    :return: new coordinates of vec1 
    """
    #spherical coords vec0
    x0, y0, z0 = vec0
    r0 = sqrt(x0**2+y0**2+z0**2) # radius from origin
    t0 = atan2(y0,x0) # angle on x-y plane
    p0 = acos(z0/r0) # angle from z-axis
    # spherical coords vec1
    x1, y1, z1 = vec1
    r1 = sqrt(x1**2+y1**2+z1**2)
    t1 = atan2(y1,x1)
    p1 = acos(z1/r1)
    # get new coords for vec1:
    t1_new = t0-t1
    p1_new = p0-p1
    x1_new = r1*sin(p1_new)*cos(t1_new)
    y1_new = r1*sin(p1_new)*sin(t1_new)
    z1_new = r1*cos(p1_new)
    return (x1_new, y1_new, z1_new)

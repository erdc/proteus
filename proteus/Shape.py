"""
Module creating predifined shapes. Each shape needs a proteus.Domain as argument.
Boundary conditions objects are automatically created for each facet (3D) or segment (2D) defining the shape.
classes:
- Shape: super class, regroups functions common to all shapes
- Cuboid: creates a cuboid 
- Rectangle: creates a rectangle
- Custom: creates a custom shape from a given set vertices, facets, etc.
"""

import BC as bc
import numpy as np
from math import cos, sin, sqrt, atan2, acos
from proteus import AuxiliaryVariables, Archiver
from proteus.Profiling import logEvent


class BCContainer(object):
    def __init__(self, BC_dict):
        self.__dict__ = BC_dict


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
        self.barycenter = None
        self.BC_list = None
        self.free_x = (1, 1, 1)
        self.free_r = (1, 1, 1)
        self._snv = len(self.domain.vertices)  # total number of vertices in domain when shape.__init__
        if self.domain.nd == 3:
            self._snf = len(self.domain.facets)
        self._sns = len(self.domain.segments)
        self._snr = len(self.domain.regions)
        self._snbc = len(self.domain.bc)
        self._snh = len(self.domain.holes)

    def _addShape(self):
        """
        Adds Shape information to the domain
        """
        # add new information to the domain
        # get maximum flag defined in domain so far
        # need to add +1 for flags as 0 cannot be used
        flag = 1
        if self.domain.nd == 2:
            if self.domain.vertexFlags or self.domain.segmentFlags:
                flag += max([max(self.domain.vertexFlags),
                             max(self.domain.segmentFlags)])
        elif self.domain.nd == 3:
            if self.domain.vertexFlags or self.domain.segmentFlags or self.domain.facetFlags:
                flag += max([max(self.domain.vertexFlags),
                             max(self.domain.segmentFlags),
                             max(self.domain.facetFlags)])
        self.domain.vertices += self.vertices.tolist()
        self.domain.vertexFlags += (self.vertexFlags+flag).tolist()
        if self.segments is not None:
            self.domain.segments += (self.segments+self._snv).tolist()
            self.domain.segmentFlags += (self.segmentFlags+flag).tolist()
        if self.facets is not None:
            self.domain.facets += (self.facets+self._snv).tolist()
            self.domain.facetFlags += (self.facetFlags+flag).tolist()
        if self.holes is not None:
            self.domain.holes += self.holes.tolist()
        if self.regions is not None:
            self.domain.regions += self.regions.tolist()
            self.domain.regionFlags += (self.regionFlags+flag).tolist()
        if len(self.domain.bc) == 0: # need to add None boundary condition at 0 indice
            self.domain.bc += [bc.BoundaryConditions()]
        self.domain.bc += self.BC_list
        if self.domain.barycenters is not None:
            self.domain.barycenters = np.append(self.domain.barycenters, self.barycenters, axis=0)
        else:
            self.domain.barycenters = np.array([[0., 0., 0.]])
            self.domain.barycenters = np.append(self.domain.barycenters, self.barycenters, axis=0)
        self.domain.update()

    def _updateDomain(self):
        """
        Updates domain when shape has been changed.
        -----------------
        Updated parameters:
        - vertices
        - regions
        - holes
        - barycenters
        """
        nv = len(self.vertices)
        self.domain.vertices[self._snv:self._snv+nv] = self.vertices.tolist()
        if self.regions is not None:
            nr = len(self.regions)
            self.domain.regions[self._snr:self._snr+nr] = self.regions.tolist()
        if self.holes is not None:
            nh = len(self.holes)
            self.domain.holes[self._snh:self._snh+nh]
        if self.domain.nd == 2:
            ns = len(self.segments)
            self.domain.barycenters[self._sns+1:self._sns+1+ns] = self.barycenter
        if self.domain.nd == 3:
            nf = len(self.facets)
            self.domain.barycenters[self._snf+1:self._snf+1+nf] = self.barycenter
        self.domain.update()

    def setPosition(self, coords):
        """
        Set position of the Shape (coords from the barycenter)
        :arg coords: new set of coordinates for the Shape
        """
        old_coords = np.array(self.barycenter)
        if self.domain.nd == 2 and len(old_coords) == 3:
            trans = coords - old_coords[:2]
        self.translate(trans)

    def setBarycenter(self, barycenter):
        """
        Set barycenter (center of mass) of the shape
        :arg barycenter: coordinates of barycenter
        """
    #     new_barycenter = np.einsum('ij,i->j', self.coords_system, barycenter)
    #     # above can be used to rotate barycenter..
        if self.domain.nd == 2:
            if len(barycenter) == 2:
                self.barycenter = np.array([barycenter[0], barycenter[1], 0.])
            elif len(barycenter) == 3:
                self.barycenter = np.array(barycenter)
            self.barycenters[:] = self.barycenter
            self.domain.barycenters[self._sns+1:self._sns+1+len(self.segments)] = self.barycenters
        if self.domain.nd == 3:
            self.barycenter = np.array(barycenter)
            self.barycenters[:] = self.barycenter
            self.domain.barycenters[self._snf+1:self._snf+1+len(self.facets)] = self.barycenters
        self._updateDomain()

    def setConstraints(self, free_x, free_r):
        """
        Sets constraints on the Shape
        :arg free_x: translational constraints
        :arg free_r: rotational constraints
        """
        self.free_x = np.array(free_x)
        self.free_r = np.array(free_r)

    def setRegions(self, regions):
        self.regions = np.array([regions])
        self._updateDomain()
    
    def rotate(self, rot, axis=(0,0,1), pivot=None):
        """
        Function to rotate Shape
        :arg rot: angle of rotation in radians (float)
        :arg axis: axis of rotation (list or array)
        :arg pivot: point around which the Shape rotates
        -----------------
        Rotated parameters:
        - vertices
        - holes
        - regions
        - local coordinate system
        - boundary orientations
        - coordinates
        - barycenters
        """
        nd = self.domain.nd
        if pivot is None:
            pivot = self.barycenter
        if self.domain.nd == 2:
            pivot = pivot[:2]
            self.vertices[:] = rotation2D(points=self.vertices, rot=rot, pivot=pivot)
            if self.holes is not None:
                self.holes[:] = rotation2D(points=self.holes, rot=rot, pivot=pivot)
            if self.regions is not None:
                self.regions[:] = rotation2D(points=self.regions, rot=rot, pivot=pivot)
            self.coords_system[:] = rotation2D(points=self.coords_system, rot=rot, pivot=(0.,0.))
            self.b_or[:] = rotation2D(points=self.b_or, rot=rot, pivot=(0., 0.))
            self.coords[:] = rotation2D(points=self.coords, rot=rot, pivot=pivot)
            self.barycenter[:2] = rotation2D(points=self.barycenter[:nd], rot=rot, pivot=pivot)
        elif self.domain.nd == 3:
            self.vertices[:] = rotation3D(points=self.vertices, rot=rot, axis=axis, pivot=pivot)
            if self.holes is not None:
                self.holes[:] = rotation3D(points=self.holes, rot=rot, axis=axis, pivot=pivot)
            if self.regions is not None:
                self.regions[:] = rotation3D(points=self.regions, rot=rot, axis=axis, pivot=pivot)
            self.coords_system[:] = rotation3D(points=self.coords_system, rot=rot, axis=axis, pivot=(0.,0.,0.))
            self.b_or[:] = rotation3D(points=self.b_or, rot=rot, axis=axis, pivot=(0., 0., 0.))
            self.barycenter[:] = rotation3D(points=self.barycenter, rot=rot, axis=axis, pivot=pivot)
            self.coords[:] = rotation3D(points=self.coords, rot=rot, axis=axis, pivot=pivot)
        self._updateDomain()

    def translate(self, trans):
        """
        Function to translate Shape
        :arg trans: translation values
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
        self._updateDomain()

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

    def getPosition(self):
        return self.barycenter

    def getRotation(self):
        return self.coords_system

    def getInertia(self, vec=(0.,0.,1.), pivot=None):
        if pivot is None:
            pivot = self.barycenter
        # Pivot coords relative to shape centre of mass
        pivot = pivot-np.array(self.barycenter)
        # making unity vector/axis of rotation
        vec = vx, vy, vz = np.array(vec)
        length_vec = sqrt(vx**2+vy**2+vz**2)
        vec = vec/length_vec
        # vector relative to original position of shape:
        if self.domain.nd == 3:
            vec = relative_vec(vec, self.coords_system[2])
        if self.domain.nd == 2:
            vec = (0., 1., 0.)
        cx, cy, cz = vec
        # getting the tensor for calculaing moment of inertia from arbitrary axis 
        vt = np.array([[cx**2, cx*cy, cx*cz],
                        [cx*cy, cy**2, cy*cz],
                        [cx*cz, cy*cz, cz**2]])
        # total moment of inertia
        I = np.einsum('ij,ij->', self.mass*self.It, vt)
        return I

    def setRigidBody(self):
        self.RigidBodyLink = RigidBody(shape=self)

    def setTank(self):
        for bc in self.BC_list:
            bc.setTank()


class Cuboid(Shape):
    """
    Class to create a cuboid
    :arg domain: domain of the cuboid
    :arg dim: dimensions of the cuboid (list or array)
    :arg coords: coordinates of the cuboid (list or array)
    """
    count = 0
    def __init__(self, domain, dim=(0.,0.,0.), coords=(0.,0.,0.), barycenter=None, tank=False, add=True):
        Shape.__init__(self, domain)
        if add is True:
            self.__class__.count += 1
            self.name = "cuboid" + str(self.__class__.count)
        self.dim = L, W, H = dim  # length, width height
        self.volume = L*W*H
        self.coords = x, y, z = np.array(coords)
        self.coords_system = np.eye(3)
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
        self.facetFlags = np.array([0, 1, 2, 3, 4, 5])  # bottom, front, right, back, left, top
        self.vertexFlags = np.array([0, 0, 0, 0, 5, 5, 5, 5])  # only top and bottom for vertices
        self.segmentFlags = np.array([0, 0, 0, 0, 5, 5, 5, 5, 1, 1, 3, 3])
        self.regionFlags = np.array([0])
        # Initialize (empty) boundary conditions
        b_or = self.b_or
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
        self.barycenters = np.array([self.barycenter for facet in self.facets])
        self.It = np.array([[(W**2.+H**2.)/12., 0, 0],
                            [0, (L**2.+H**2.)/12., 0],
                            [0, 0, (W**2.+L**2.)/12.]])
        if add is True:
            self._addShape()  # adding shape to domain

    def _setInertiaTensor(self):
        L, W, H = self.dim
        self.It[:] = [[(W**2.+H**2.)/12., 0, 0],
                      [0, (L**2.+H**2.)/12., 0],
                      [0, 0, (W**2.+L**2.)/12.]]

    def setDimensions(self, dim):
        """
        Set dimensions of the shape
        :arg dim: new dimensions of the Shape
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
        self._setInertiaTensor()
        self._updateDomain()


class Rectangle(Shape):
    """
    Class to create a rectangle
    :arg domain: domain of the rectangle
    :arg dim: dimensions of the rectangle (list or array)
    :arg coords: coordinates of the rectangle (list or array)
    """
    count = 0
    def __init__(self, domain, dim=(0.,0.), coords=(0.,0.), barycenter=None, tank=False, add=True):
        Shape.__init__(self, domain)
        if add is True:
            self.__class__.count += 1
            self.name = "rectangle" + str(self.__class__.count)
        self.dim = L, H = dim  # length, height
        self.coords = x, y = np.array(coords)
        self.coords_system = np.eye(2)
        self.vertices = np.array([[x-0.5*L, y-0.5*H],
                                  [x+0.5*L, y-0.5*H],
                                  [x+0.5*L, y+0.5*H],
                                  [x-0.5*L, y+0.5*H]])
        self.segments = np.array([[0, 1], [1, 2], [2, 3], [3, 0]])  # bottom, right, top, left
        if barycenter is not None:
            if len(barycenter) == 2:
                self.barycenter = np.array([barycenter[0], barycenter[1], 0.])
            if len(barycenter) == 3:
                self.barycenter = np.array([barycenter[0], barycenter[1], barycenter[2]])
        else:
            self.barycenter = np.array([coords[0], coords[1], 0.])
        self.barycenters = np.array([self.barycenter for segment in self.segments])
        self.b_or = np.array([[0., -1.],
                              [1., 0.],
                              [0., 1.],
                              [-1., 0.]])
        self.regions = np.array([[x, y]])
        self.segmentFlags = np.array([0, 1, 2, 3]) # bottom, right, top, left
        self.vertexFlags = np.array([0, 0, 2, 2])  # bottom, bottom, top, top
        self.regionFlags = np.array([0])
        self.BC_dict = {'bottom': bc.BoundaryConditions(b_or=self.b_or, b_i=0),
                        'right': bc.BoundaryConditions(b_or=self.b_or, b_i=1),
                        'top': bc.BoundaryConditions(b_or=self.b_or, b_i=2),
                        'left': bc.BoundaryConditions(b_or=self.b_or, b_i=3)}
        self.BC_list = [self.BC_dict['bottom'],
                        self.BC_dict['right'],
                        self.BC_dict['top'],
                        self.BC_dict['left']]
        self.BC = BCContainer(self.BC_dict)
        self.It = np.array([[(L**2.)/3., 0],
                            [0, (H**2.)/3.]])
        if add is True:
            self._addShape()  # adding shape to domain

    def _setInertiaTensor(self):
        """
        Set the (new) inertia tensor of the shape
        """
        L, H = self.dim
        self.It[:] = [[0., 0, 0.],
                      [0, 0., 0.],
                      [0., 0., (L**2+H**2)/12]]

    def setDimensions(self, dim):
        """
        Set dimensions of the shape
        :arg dim: new dimensions of the Shape
        """
        self.dim = dim
        L, H = dim
        x, y = self.coords
        self.vertices[:] = [[x-0.5*L, y-0.5*H],
                            [x+0.5*L, y-0.5*H],
                            [x+0.5*L, y+0.5*H],
                            [x-0.5*L, y+0.5*H]]
        self.volume = L*H
        self._setInertiaTensor()
        self._updateDomain()


class BodyCuboid(Cuboid):
    """
    Class to create a cuboid rigid body
    :arg domain: domain of the cuboid
    :arg dim: dimensions of the cuboid (list or array)
    :arg coords: coordinates of the cuboid (list or array)
    :arg barycenter: barycenter of the cuboid (list or array)
    """
    count = 0
    def __init__(self, domain, dim=(0.,0.,0.), coords=(0.,0.,0.), barycenter=None):
        Cuboid.__init__(self, domain, dim=dim, coords=coords, barycenter=barycenter, add=False)
        self.__class__.count += 1
        self.name = "body_cuboid" + str(self.__class__.count)
        self.setRigidBody()
        self.holes = np.array([coords])
        for bc in self.BC_list:
            bc.setMoveMesh(self.RigidBodyLink)
        self._addShape()


class BodyRectangle(Rectangle):
    """
    Class to create a rectangle rigid body
    :arg domain: domain of the rectangle
    :arg dim: dimensions of the rectangle (list or array)
    :arg coords: coordinates of the rectangle (list or array)
    :arg barycenter: barycenter of the cuboid (list or array)
    """
    count = 0
    def __init__(self, domain, dim=(0.,0.), coords=(0.,0.), barycenter=None):
        Rectangle.__init__(self, domain, dim=dim, coords=coords, barycenter=barycenter, add=False)
        self.__class__.count += 1
        self.name = "body_rectangle" + str(self.__class__.count)
        self.setRigidBody()
        self.holes = np.array([coords])
        for bc in self.BC_list:
            bc.setMoveMesh(self.RigidBodyLink)
        self.regions = None
        self._addShape()

            
class Tank3D(Cuboid):
    """
    Class to create a 3D tank (cuboid)
    :arg domain: domain of the tank
    :arg dim: dimensions of the tank (list or array)
    """
    count = 0
    def __init__(self, domain, dim=(0.,0.,0.), from_0=True):
        L, W, H = dim
        if from_0 is True:
            x, y, z = L/2., W/2., H/2.
        else:
            x, y, z = 0., 0., 0.
        Cuboid.__init__(self, domain, dim=dim, coords=(x, y, z), add=False)
        self.__class__.count += 1
        self.name = "tank3d" + str(self.__class__.count)
        self.from_0 = from_0
        self.regions = np.array([[x+L/2.-L/100., y+W/2.-W/100., z+H/2.-H/100.]])
        self.barycenter = np.array([0., 0., 0.])
        self.barycenters = np.array([self.barycenter for facet in self.facets])
        for bc in self.BC_list:
            bc.setTank()
        self._addShape()  # adding shape to domain

    def setDimensions(self, dim):
        """
        Set dimensions of the shape
        :arg dim: new dimensions of the Shape
        """
        L, W, H = dim
        if self.from_0 is True:
            x, y, z = L/2., W/2., H/2.
        else:
            x, y, z = 0., 0., 0.
        self.coords = np.array(x, y, z)
        self.vertices = np.array([[x-0.5*L, y-0.5*W, z-0.5*H],
                                  [x-0.5*L, y+0.5*W, z-0.5*H],
                                  [x+0.5*L, y+0.5*W, z-0.5*H],
                                  [x+0.5*L, y-0.5*W, z-0.5*H],
                                  [x-0.5*L, y-0.5*W, z+0.5*H],
                                  [x-0.5*L, y+0.5*W, z+0.5*H],
                                  [x+0.5*L, y+0.5*W, z+0.5*H],
                                  [x+0.5*L, y-0.5*W, z+0.5*H]])
        self.regions = np.array([[x+L/2.-L/100., y+W/2.-W/100., z+H/2.-H/100.]])
        self.volume = L*W*H
        self._updateDomain()


class Tank2D(Rectangle):
    """
    Class to create a 2D tank (rectangle)
    :arg domain: domain of the tank
    :arg dim: dimensions of the tank (list or array)
    """
    count = 0
    def __init__(self, domain, dim=(0.,0.), from_0=True):
        L, H = dim
        if from_0 is True:
            x, y = L/2., H/2.
        else:
            x, y = 0., 0.
        Rectangle.__init__(self, domain, dim=dim, coords=(x, y), add=False)
        self.__class__.count += 1
        self.name = "tank2d" + str(self.__class__.count)
        self.from_0 = from_0
        self.regions = np.array([[L-L/100., H-H/100.]])
        self.barycenter = np.array([0., 0., 0.])
        self.barycenters = np.array([self.barycenter for segment in self.segments])
        for bc in self.BC_list:
            bc.setTank()
        self._addShape()  # adding shape to domain

    def setDimensions(self, dim):
        """
        Set dimensions of the shape
        :arg dim: new dimensions of the Shape
        """
        self.dim = dim
        L, H = dim
        if self.from_0 is True:
            x, y = L/2., H/2.
        else:
            x, y = 0., 0.
        self.vertices = np.array([[x-0.5*L, y-0.5*H],
                                  [x+0.5*L, y-0.5*H],
                                  [x+0.5*L, y+0.5*H],
                                  [x-0.5*L, y+0.5*H]])
        self.coords = np.array(self.dim)/2.
        self.regions = np.array([[x+L/2.-L/100., y+H/2.-H/100.]])
        self._updateDomain()




class CustomShape(Shape):
    """
    Class to create a custom 2D or 3D shape
    :arg domain: domain of the shape
    :arg barycenter: barycenter of the shape (list or array)
    :arg vertices: set of vertices of the shape (list or array)
    :arg facets: set of facets of the shape (list or array)
    :arg segments: set of segments of the shape (list or array)
    :arg regions: set of regions of the shape (list or array)
    """
    count = 0
    def __init__(self, domain, barycenter=None, vertices=None, vertexFlags=None, segments=None, segmentFlags=None, facets=None, facetFlags=None, holes=None, regions=None, regionFlags=None, boundaryTags=None):
        Shape.__init__(self, domain)
        self.__class__.count += 1
        self.name = "custom" + str(self.__class__.count)
        minFlag = min(vertexFlags)
        if segmentFlags:
            minFlag = min(minFlag, min(segmentFlags))
        if facetFlags:
            minFlag = min(minFlag, min(facetFlags))
        flagSet = set()
        self.vertices = np.array(vertices)
        self.vertexFlags = np.array(vertexFlags)-minFlag
        for flag in self.vertexFlags:
            flagSet.add(flag)
        if segments:
            self.segments = np.array(segments)
            self.segmentFlags = np.array(segmentFlags)-minFlag
            for flag in self.segmentFlags:
                    flagSet.add(flag)
        if facets:
            self.facets = np.array(facets)
            self.facetFlags = np.array(facetFlags)-minFlag
            for flag in self.facetFlags:
                flagSet.add(flag)
        previous_flag = -1
        for flag in flagSet:
            assert flag == previous_flag+1, "Flags must be defined as a suite of numbers (e.g. 0, 1, 2, 3, 4 with no gap)!"
            previous_flag = flag
        if holes is not None:
            self.holes = np.array(holes)
        if regions is not None:
            self.regions = np.array(regions)
            self.regionFlags = np.array(regionFlags)
        if barycenter is not None:
            self.barycenter = np.array(barycenter)
        else:
            self.barycenter = np.zeros(domain.nd)
        if self.domain.nd == 2:
            self.barycenters = np.array([self.barycenter for segment in self.segments])
        elif self.domain.nd == 3:
            self.barycenters = np.array([self.barycenter for facet in self.facets])
        self.BC_dict = {}
        self.BC_list = [None]*len(flagSet)
        for tag, index in boundaryTags.iteritems():
            if index-minFlag in flagSet:
                self.BC_dict[tag] = bc.BoundaryConditions()
                self.BC_list[index-minFlag] = self.BC_dict[tag]
        self.BC = BCContainer(self.BC_dict)
        self._addShape()
        
    def _setInertiaTensor(self, It):
        self.It = np.array(It)


class RigidBody(AuxiliaryVariables.AV_base):

    def __init__(self, shape, he=1., cfl_target=0.9, dt_init=0.001):
        self.shape = shape
        self.dt_init = dt_init
        self.he = he
        self.cfl_target = 0.9
        self.last_position = np.array([0., 0., 0.])
        self.rotation_matrix = np.eye(3)
        self.h = np.array([0., 0., 0.])
        self.initialize_values = True

    def step(self, dt):
        nd = self.shape.domain.nd
        # displacement from force
        self.acceleration = self.F/self.shape.mass
        self.velocity = self.last_velocity + self.acceleration*dt
        self.h[:] = self.velocity*dt
        # update barycenters
        self.shape.translate(self.h[:nd])
        i0, i1 = self.nb_start, self.nb_end
        self.barycenters[i0:i1, :] = self.shape.barycenter
        self.position[:] = self.shape.barycenter
        # rotation due to moment
        if sum(self.M) != 0:
            I = self.shape.getInertia(vec=self.M, pivot=self.shape.barycenter)
            assert I != 0, "Zero inertia: inertia tensor (It) was not set correctly! It has to be set as a 3x3 numpy array"
            ang_acc = self.M[:]/I
        else:
            ang_acc = np.array([0., 0., 0.])
        self.angvel[:] = self.last_angvel+ang_acc*dt
        ang_disp = self.angvel*dt
        self.ang = np.linalg.norm(ang_disp)
        if nd == 2 and self.angvel[2] < 0:
            self.ang = -self.ang
        if self.ang != 0.:
            self.shape.rotate(rot=self.ang, axis=self.angvel, pivot=self.shape.barycenter)
            self.rotation[:nd,:nd] = self.shape.coords_system  # this rotation matrix will be used for moveMesh
            self.rotation_matrix[:] = np.dot(np.linalg.inv(self.last_rotation), self.rotation)
        else:
            self.rotation_matrix[:] = np.eye(3)

    def attachModel(self, model, ar):
        self.model = model
        self.ar = ar
        self.writer = Archiver.XdmfWriter()
        self.nd = model.levelModelList[-1].nSpace_global
        m = self.model.levelModelList[-1]
        flagMax = max(m.mesh.elementBoundaryMaterialTypes)
        flagMin = min(m.mesh.elementBoundaryMaterialTypes)
        self.nForces=flagMax+1
        return self

    def calculate_init(self):
        nd = self.shape.domain.nd
        shape = self.shape
        self.position = np.zeros(3)
        self.position[:] = self.shape.barycenter.copy()
        self.last_position[:] = self.position
        self.velocity = np.zeros(3, 'd')
        self.last_velocity = np.zeros(3, 'd')
        self.rotation = np.eye(3)
        self.rotation[:nd,:nd] = shape.coords_system
        self.last_rotation = np.eye(3)
        self.last_rotation[:nd,:nd] = shape.coords_system
        self.F = np.zeros(3, 'd')
        self.M = np.zeros(3, 'd')
        self.last_F = np.zeros(3, 'd')
        self.last_M = np.zeros(3, 'd')
        self.ang = 0.
        self.barycenters = shape.domain.barycenters
        self.angvel = np.zeros(3, 'd')
        self.last_angvel = np.zeros(3, 'd')
        if nd == 2:
            self.nb_start = self.shape._sns+1
            self.nb_end = self.nb_start+len(self.shape.segments)
        if nd == 3:
            self.nb_start = self.shape._snf+1  # must skip indice 0 for forces / barycenters
            self.nb_end = self.nb_start + len(self.shape.facets)
        if nd == 2:
            self.Fg = self.shape.mass*np.array([0., -9.81, 0.])
        if nd == 3:
            self.Fg = self.shape.mass*np.array([0., 0., -9.81])

    def calculate(self):
        if self.initialize_values is True:
            self.calculate_init()
            self.initialize_values = False
        import copy
        self.last_position[:] = self.position
        self.last_velocity[:] = self.velocity
        self.last_rotation[:] = self.rotation
        self.last_angvel[:] = self.angvel
        # store forces
        self.last_F[:] = self.F
        self.last_M[:] = self.M
        # self.last_rotation_inv = np.linalg.inv(self.last_rotation)
        try:
            dt = self.model.levelModelList[-1].dt_last
        except:
            dt = self.dt_init
        i0, i1 = self.nb_start, self.nb_end
        F = np.sum(self.model.levelModelList[-1].coefficients.netForces_p[i0:i1,:] + self.model.levelModelList[-1].coefficients.netForces_v[i0:i1,:], axis=0) + self.Fg
        M = np.sum(self.model.levelModelList[-1].coefficients.netMoments[i0:i1,:], axis=0)
        logEvent(`self.model.levelModelList[-1].coefficients.netForces_p+ self.model.levelModelList[-1].coefficients.netForces_v`)
        self.F[:] = F*self.shape.free_x
        self.M[:] = M*self.shape.free_r
        logEvent("==============================================================")
        logEvent("================= Floating Body Calculation ==================")
        logEvent("==============================================================")
        logEvent("Name: " + `self.shape.name`)
        logEvent("==============================================================")
        logEvent("[proteus]     t=%1.5fsec" % (self.model.stepController.t_model_last))
        logEvent("[proteus]    dt=%1.5fsec" % (dt))
        logEvent("[proteus]     F=(%21.16e, %21.16e, %21.16e)" % (F[0], F[1], F[2]))
        logEvent("[proteus] F*DOF=(%21.16e, %21.16e, %21.16e)" % (self.F[0], self.F[1], self.F[2]))
        logEvent("[proteus]     M=(%21.16e, %21.16e, %21.16e)" % (M[0], M[1], M[2]))
        logEvent("[proteus] M*DOF=(%21.16e, %21.16e, %21.16e)" % (self.M[0], self.M[1], self.M[2]))
        logEvent("[body] ============= Pre-calculation attributes  =============")
        logEvent("[body]      t=%1.5fsec" % (self.model.stepController.t_model_last))
        logEvent("[body]    pos=(%21.16e, %21.16e, %21.16e)" % (self.position[0], self.position[1], self.position[2]))
        logEvent("[body]    vel=(%21.16e, %21.16e, %21.16e)" % (self.velocity[0], self.velocity[1], self.velocity[2]))
        self.step(dt)
        logEvent("[body] ============= Post-calculation attributes =============")
        logEvent("[body]      t=%1.5fsec" % (self.model.stepController.t_model))
        logEvent("[body]    pos=(%21.16e, %21.16e, %21.16e)" % (self.position[0], self.position[1], self.position[2]))
        logEvent("[body]    vel=(%21.16e, %21.16e, %21.16e)" % (self.velocity[0], self.velocity[1], self.velocity[2]))
        logEvent("[body]      h=(%21.16e, %21.16e, %21.16e)" % (self.h[0], self.h[1], self.h[2]))
        logEvent("[body]  r vel=(%21.16e, %21.16e, %21.16e)" % (self.angvel[0], self.angvel[1], self.angvel[2]))
        if sum(self.angvel) != 0:
            rot_axis=self.angvel/np.linalg.norm(self.angvel)
        else:
            rot_axis=(0., 0., 0.)
        logEvent("[body] r axis=(%21.16e, %21.16e, %21.16e)" % (rot_axis[0], rot_axis[1], rot_axis[2]))
        logEvent("[body]  r ang=(%21.16e)" % (self.ang))
        logEvent("[body] r matrix:")
        logEvent("[body] " + `self.rotation_matrix[0].tolist()`)
        logEvent("[body] " + `self.rotation_matrix[1].tolist()`)
        logEvent("[body] " + `self.rotation_matrix[2].tolist()`)
        logEvent("==============================================================")







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
    :return points_rot: the rotated set of points (numpy array)
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
    :return points_rot: the rotated set of points (numpy array)
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
    M = reduce(np.dot, [T, Rx, Ry, Rz, np.linalg.inv(Ry), np.linalg.inv(Rx), np.linalg.inv(T)])
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

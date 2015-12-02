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
from proteus import AuxiliaryVariables, Archiver, Comm, Profiling
from proteus.Profiling import logEvent
from itertools import compress
import csv
import os


class BCContainer(object):
    def __init__(self, BC_dict):
        self.__dict__ = BC_dict


class Shape:
    """
    Class defining a shape
    :param domain: domain in which the shape is defined
    """

    def __init__(self, domain, dim=None, coords=None):
        self.instance_i = len(domain.shape_list)
        domain.shape_list.append(self)
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
        self.record_values = False
        self.BC_list = None
        self.free_x = (1, 1, 1)
        self.free_r = (1, 1, 1)
        self._snv = None  # total number of vertices in domain when shape.__init__
        self._snf = None
        self._sns = None
        self._snr = None
        self._snh = None
        self._snbc = None

    def _addShape(self):
        """
        Adds Shape information to the domain
        """
        # add new information to the domain
        # get maximum flag defined in domain so far
        # need to add +1 for flags as 0 cannot be used
        self._snbc = len(self.domain.bc)
        self._enbc = self._snbc + len(self.BC_list)
        self.domain.bc += self.BC_list
        self.domain.barycenters = np.append(self.domain.barycenters, self.barycenters, axis=0)
        self._snv = len(self.domain.vertices)  # total number of vertices in domain when shape.__init__
        self._env = self._snv + len(self.vertices)  # total number of vertices in domain when shape.__init__
        flag = self._snbc-1
        self.domain.vertices += self.vertices.tolist()
        self.domain.vertexFlags += (self.vertexFlags+flag).tolist()
        if self.segments is not None:
            self._sns = len(self.domain.segments)
            self._ens = self._sns + len(self.segments)
            self.domain.segments += (self.segments+self._snv).tolist()
            self.domain.segmentFlags += (self.segmentFlags+flag).tolist()
        if self.facets is not None:
            self._snf = len(self.domain.facets)
            self._enf = self._snf + len(self.facets)
            self.domain.facets += (self.facets+self._snv).tolist()
            self.domain.facetFlags += (self.facetFlags+flag).tolist()
        if self.holes is not None:
            self._snh = len(self.domain.holes)
            self._enh = self._snh + len(self.holes)
            self.domain.holes += self.holes.tolist()
        if self.regions is not None:
            self._snr = len(self.domain.regions)
            self._enr = self._snr + len(self.regions)
            self.domain.regions += self.regions.tolist()
            self.domain.regionFlags += (self.regionFlags+flag).tolist()
        self.domain.update()

    def _resetShape(self):
        """
        Updates domain when shape has been changed.
        This function deletes and resets values of the shape in the domain using their attribute indice.
        The shapes that were defined after the current shape have their attribute indice updated.
        -----------------
        Updated parameters:
        - vertices
        - segments (if any)
        - facets (if any)
        - regions (if any)
        - holes (if any)
        - boundary conditions
        - barycenters
        """
        domain = self.domain
        flag = self._snbc-1
        # vertices
        i0, i1 = self._snv, self._env
        del domain.vertices[i0:i1]
        del domain.vertexFlags[i0:i1]
        domain.vertices[i0:i0] = self.vertices.tolist()
        domain.vertexFlags[i0:i0] = (self.vertexFlags+flag).tolist()
        self._env = i0 + len(self.vertices)
        vdiff = self._env - i1
        # segments
        if self._sns is not None:
            i0, i1 = self._sns, self._ens
            del domain.segments[i0:i1]
            del domain.segmentFlags[i0:i1]
            domain.segments[i0:i0] = (self.segments+self._snv).tolist()
            domain.segmentFlags[i0:i0] = (self.segmentFlags+flag).tolist()
            self._ens = i0 + len(self.segments)
            sdiff = self._ens - i1
        # facets
        if self._snf is not None:
            i0, i1 = self._snf, self._enf
            del domain.facets[i0:i1]
            del domain.facetFlags[i0:i1]
            domain.facets[i0:i0] = (self.facets+self._snv).tolist()
            domain.facetFlags[i0:i0] = (self.facetFlags+flag).tolist()
            self._enf = i0 + len(self.facets)
            fdiff = self._enf - i1
        # regions
        if self._snr is not None:
            i0, i1 = self._snr, self._enr
            del domain.regions[i0:i1]
            del domain.regionFlags[i0:i1]
            domain.regions[i0:i0] = self.regions.tolist()
            domain.regionFlags[i0:i0] = (self.regionFlags+self._snr).tolist()
            self._enr = i0 + len(self.regions)
            rdiff = self._enr - i1
        # holes
        if self._snh is not None:
            i0, i1 = self._snh, self._enh
            del domain.holes[i0:i1]
            domain.holes[i0:i0] = self.holes.tolist()
            self._enh = i0 + len(self.holes)
            hdiff = self._enh - i1
        # boundary conditions and barycenters
        i0, i1 = self._snbc, self._enbc
        del domain.bc[i0:i1]
        # this is a numpy array!!
        ilist = [i0 + i for i in range(i1-i0)]
        domain.barycenters = np.delete(domain.barycenters, ilist, axis=0)
        domain.bc[i0:i0] = self.BC_list
        domain.barycenters = np.insert(domain.barycenters, i0, self.barycenters, axis=0)
        self._enbc = i0 + len(self.BC_list)
        # updating the attributes of the shape defined after current shape
        bcdiff = self._enbc - i1
        for ind, inst in enumerate(self.domain.shape_list):
            if ind > self.instance_i:
                inst._snv += vdiff
                inst._env += vdiff
                inst._snbc += bcdiff
                inst._enbc += bcdiff
                if inst.segments is not None:
                    inst._sns += sdiff
                    inst._ens += sdiff
                    inst.segments += vdiff
                    inst.segmentFlags += bcdiff
                if inst.facets is not None:
                    inst._snf += fdiff
                    inst._enf += fdiff
                    inst.facets += vdiff
                    inst.facetFlags += bcdiff
                if inst.regions is not None:
                    inst._snr += rdiff
                    inst._enr += rdiff
                    inst.regionFlags += bcdiff
                if inst.holes is not None:
                    inst._snh += fdiff
                    inst._enh += fdiff

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
        self.domain.vertices[self._snv:self._env] = self.vertices.tolist()
        if self.regions is not None:
            self.domain.regions[self._snr:self._enr] = self.regions.tolist()
        if self.holes is not None:
            self.domain.holes[self._snh:self._enh]
        self.domain.barycenters[self._snbc:self._enbc] = self.barycenter
        self.domain.update()

    def setPosition(self, coords):
        """
        Set position of the Shape (coords from the barycenter)
        :arg coords: new set of coordinates for the Shape
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
        :arg barycenter: coordinates of barycenter
        """
    #     new_barycenter = np.einsum('ij,i->j', self.coords_system, barycenter)
    #     # above can be used to rotate barycenter..
        if self.domain.nd == 2:
            if len(barycenter) == 2:
                self.barycenter = np.array([barycenter[0], barycenter[1], 0.])
            elif len(barycenter) == 3:
                self.barycenter = np.array(barycenter)
        if self.domain.nd == 3:
            self.barycenter = np.array(barycenter)
        self.barycenters[:] = self.barycenter
        self.domain.barycenters[self._snbc:self._enbc] = self.barycenters
        self._resetShape()

    def setConstraints(self, free_x, free_r):
        """
        Sets constraints on the Shape
        :arg free_x: translational constraints
        :arg free_r: rotational constraints
        """
        self.free_x = np.array(free_x)
        self.free_r = np.array(free_r)

    def setRegions(self, regions):
        """
        Sets new regions for the Shape
        :arg regions: coordinate of the new region(s) (list or array)
        """
        self.regions = np.array([regions])
        if not hasattr(self, '_snr'):
            self._snr = len(self.domain.regions)
        self._enr = self._snr + len(self.regions)
        self._resetShape()

    def rotate(self, rot, axis=(0, 0, 1), pivot=None):
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

    def getInertia(self, vec=(0., 0., 1.), pivot=None):
        try:
            self.It
        except NameError:
            print "This shape does not have its inertia tensor defined! (", self.name, ")"
        if pivot is None:
            pivot = self.barycenter
        # Pivot coords relative to shape centre of mass
        pivot = pivot-np.array(self.barycenter)
        # making unity vector/axis of rotation
        vec = vx, vy, vz = np.array(vec)
        length_vec = sqrt(vx**2+vy**2+vz**2)
        vec = vec/length_vec
        # vector relative to original position of shape:
        if self.domain.nd == 2:
            I = self.It*self.mass
        elif self.domain.nd == 3:
            vec = relative_vec(vec, self.coords_system[2])
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
        self.holes = np.array([self.coords])
        for boundcond in self.BC_list:
            boundcond.setMoveMesh(self.RigidBodyLink)

    def setTank(self):
        for boundcond in self.BC_list:
            boundcond.setTank()

    def setRecordValues(self, all_values=False, time=True, pos=False,
                        pos_x=False, pos_y=False, pos_z=False, rot=False,
                        rot_x=False, rot_y=False, rot_z=False, F=False,
                        Fx=False, Fy=False, Fz=False, M=False, Mx=False,
                        My=False, Mz=False, inertia=False, vel=False,
                        vel_x=False, vel_y=False, vel_z=False, acc=False,
                        acc_x=False, acc_y=False, acc_z=False, filename=None):
        """
        values to be recorded in a csv file (for rigid bodies)
        """
        self.record_values = True
        if pos is True:
            pos_x = pos_y = pos_z = True
        if rot is True:
            rot_x = rot_y = rot_z = True
        if F is True:
            Fx = Fy = Fz = True
        if M is True:
            Mx = My = Mz = True
        if vel is True:
            vel_x = vel_y = vel_z = True
        if acc is True:
            acc_x = acc_y = acc_z = True
        self.record_bool = [time, pos, pos_x, pos_y, pos_z, rot, rot_x, rot_y,
                            rot_z, F, Fx, Fy, Fz, M, Mx, My, Mz, inertia, vel_x, vel_y,
                            vel_z, acc_x, acc_y, acc_z]
        if all_values is True:
            self.record_bool = [True for value in self.record_bool]
        self.record_names = ['time', 'pos_x', 'pos_y', 'pos_z',
                             'rot_x', 'rot_y', 'rot_z', 'Fx', 'Fy', 'Fz',
                             'Mx', 'My', 'Mz', 'inertia', 'vel_x', 'vel_y', 'vel_z',
                             'acc_x', 'acc_y', 'acc_z']
        self.record_names = list(compress(self.record_names, self.record_bool))
        if filename is None:
            self.record_filename = 'record_' + self.name + '.csv'
        else:
            self.record_filename = filename + '.csv'


class Cuboid(Shape):
    """
    Class to create a cuboid
    :arg domain: domain of the cuboid
    :arg dim: dimensions of the cuboid (list or array)
    :arg coords: coordinates of the cuboid (list or array)
    """
    count = 0
    def __init__(self, domain, dim=(0., 0., 0.), coords=(0., 0., 0.),
                 barycenter=None, tank=False, add=True):
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
        self.facetFlags = np.array([1, 2, 3, 4, 5, 6])  # bottom, front, right, back, left, top
        self.vertexFlags = np.array([1, 1, 1, 1, 6, 6, 6, 6])  # only top and bottom for vertices
        self.segmentFlags = np.array([1, 1, 1, 1, 6, 6, 6, 6, 2, 2, 4, 4])
        self.regionFlags = np.array([1])
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
    def __init__(self, domain, dim=(0., 0.), coords=(0., 0.), barycenter=None,
                 tank=False, add=True):
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
        if add is True:
            self._addShape()  # adding shape to domain

    def _setInertiaTensor(self):
        """
        Set the (new) inertia tensor of the shape
        """
        L, H = self.dim
        self.It[:] = (L**2+H**2)/12

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

    def __init__(self, domain, dim=(0., 0., 0.), coords=(0., 0., 0.),
                 barycenter=None):
        Cuboid.__init__(self, domain, dim=dim, coords=coords,
                        barycenter=barycenter, add=False)
        self.__class__.count += 1
        self.name = "body_cuboid" + str(self.__class__.count)
        self.setRigidBody()
        self.regions = None
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

    def __init__(self, domain, dim=(0., 0.), coords=(0., 0.), barycenter=None):
        Rectangle.__init__(self, domain, dim=dim, coords=coords,
                           barycenter=barycenter, add=False)
        self.__class__.count += 1
        self.name = "body_rectangle" + str(self.__class__.count)
        self.setRigidBody()
        self.regions = None
        self._addShape()


class Tank3D(Shape):
    """
    Class to create a 3D tank (cuboid)
    :arg domain: domain of the tank
    :arg dim: dimensions of the tank (list or array)
    """
    count = 0

    def __init__(self, domain, dim=(0., 0., 0.), from_0=True, leftSponge=None,
                 rightSponge=None, frontSponge=None, backSponge=None):
        Shape.__init__(self, domain)
        self.__class__.count += 1
        self.name = "tank3d" + str(self.__class__.count)
        self.from_0 = from_0
        self.holes = None
        self.coords_system = np.eye(3)
        self.leftSponge = leftSponge
        self.rightSponge = rightSponge
        self.backSponge = backSponge
        self.frontSponge = frontSponge
        self.boundaryTags = {'bottom': 1,
                             'front': 2,
                             'right': 3,
                             'back': 4,
                             'left': 5,
                             'top': 6,
                             'sponge': 7}
        self.b_or = np.array([[0.,  0., -1.],
                              [-1., 0.,  0.],
                              [0.,  1.,  0.],
                              [1.,  0.,  0.],
                              [0., -1.,  0.],
                              [0.,  0.,  1.]])
        self.BC_dict = {'bottom': bc.BoundaryConditions(b_or=self.b_or, b_i=0),
                        'front': bc.BoundaryConditions(b_or=self.b_or, b_i=1),
                        'right': bc.BoundaryConditions(b_or=self.b_or, b_i=2),
                        'back': bc.BoundaryConditions(b_or=self.b_or, b_i=3),
                        'left': bc.BoundaryConditions(b_or=self.b_or, b_i=4),
                        'top': bc.BoundaryConditions(b_or=self.b_or, b_i=5),
                        'sponge': bc.BoundaryConditions()}
        self.BC_list = [self.BC_dict['bottom'],
                        self.BC_dict['front'],
                        self.BC_dict['right'],
                        self.BC_dict['back'],
                        self.BC_dict['left'],
                        self.BC_dict['top'],
                        self.BC_dict['sponge']]
        self.BC = BCContainer(self.BC_dict)
        for i in range(0, 5):
            self.BC_list[i].setTank()
        self.barycenter = np.array([0., 0., 0.])
        self.barycenters = np.array([self.barycenter for bt in self.boundaryTags])
        self.setDimensions(dim)
        self.porosityTypes = np.ones(len(self.regionFlags))
        self.dragAlphaTypes = np.zeros(len(self.regionFlags))
        self.dragBetaTypes = np.zeros(len(self.regionFlags))
        self.epsFact_solid = np.zeros(len(self.regionFlags))
        self.zones = {}
        self.RelaxationZones = RelaxationZoneWaveGenerator(self.zones, shape=self)
        self._addShape()  # adding shape to domain

    def setSponge(self, left=None, right=None, back=None, front=None):
        ''' Not working yet! '''
        self.leftSponge = left
        self.rightSponge = right
        self.backSponge = back
        self.frontSponge = front
        self.setDimensions(self.dim)
        self._resetShape()

    def setDimensions(self, dim):
        L, W, H = dim
        self.dim = dim
        if self.from_0 is True:
            x, y, z = L/2., W/2., H/2.
        else:
            x, y, z = 0., 0., 0.
        self.coords = [x, y, z]
        x0, x1 = x-0.5*L, x+0.5*L
        y0, y1 = y-0.5*H, y+0.5*H
        z0, z1 = z-0.5*L, z+0.5*L
        # ---------------------------------------------
        # first add all vecors, facets, regions at the bottom
        # ---------------------------------------------
        bt = self.boundaryTags
        leftSponge = self.leftSponge or 0.
        backSponge = self.backSponge or 0.
        rightSponge =self.rightSponge or 0.
        frontSponge = self.frontSponge or 0.
        vertices = [[x0+frontSponge, y0+leftSponge, z0],
                    [x1-backSponge, y0+leftSponge, z0],
                    [x1-backSponge, y1-rightSponge, z0],
                    [x0+frontSponge, y1-rightSponge, z0]]
        vertexFlags = [bt['bottom'], bt['bottom'], bt['bottom'], bt['bottom']]
        segments = [[0, 1], [1, 2], [2, 3], [3, 0]]
        segmentFlags = [bt['bottom'], bt['bottom'], bt['bottom'], bt['bottom']]
        facets = [[[0, 1, 2, 3]]]
        facetFlags = [bt['bottom']]
        regions = [[((x0+frontSponge)+(x1-backSponge))/2., ((y0+leftSponge)+(y1-rightSponge))/2., (z0+z1)/2.]]
        regionFlags = [1]
        self.regionIndice = {0: 'tank'}
        v_i = 4  # index of next vector to add
        r_i = 1  # index of next region to add
        nb_sponge = 0  # number of sponge layers defined

        if leftSponge:
            vertices += [[x0+frontSponge, y0, z0],
                         [x1-backSponge, y0, z0]]
            segments += [[0, v_i], [v_i, v_i+1], [v_i+1, 1]]
            facets += [[[0, 1, v_i+1, v_i]]]
            regions += [[((x0+frontSponge)+(x1-backSponge))/2., (y0+(y0+leftSponge))/2., (z0+z1)/2.]]
            self.regionIndice[r_i] = 'leftSponge'
            regionFlags += [r_i+1]
            v_i += 2  # 2 vertices were added
            r_i += 1  # 1 region was added
            nb_sponge += 1
        if backSponge:
            vertices += [[x1, y0+leftSponge, z0],
                         [x1, y1-rightSponge, z0]]
            segments += [[1, v_i], [v_i, v_i+1], [v_i+1, 2]]
            facets += [[[1, 2, v_i+1, v_i]]]
            regions += [[((x1-backSponge)+x1)/2., ((y0+leftSponge)+(y1-rightSponge))/2., (z0+z1)/2.]]
            self.regionIndice[r_i] = 'backSponge'
            regionFlags += [r_i+1]
            v_i += 2
            r_i += 1
            nb_sponge += 1
        if rightSponge:
            vertices += [[x1-backSponge, y1, z0],
                         [x0+frontSponge, y1, z0]]
            segments += [[2, v_i], [v_i, v_i+1], [v_i+1, 3]]
            facets += [[[2, 3, v_i+1, v_i]]]
            regions += [[((x0+frontSponge)+(x1-backSponge))/2., (y1+(y1-rightSponge))/2., (z0+z1)/2.]]
            self.regionIndice[r_i] = 'rightSponge'
            regionFlags += [r_i+1]
            v_i += 2
            r_i += 1
            nb_sponge += 1
        if frontSponge:
            vertices += [[x0, y1-rightSponge, z0],
                         [x0, y0+leftSponge, z0]]
            segments += [[3, v_i], [v_i, v_i+1], [v_i+1, 0]]
            facets += [[[3, 0, v_i+1, v_i]]]
            regions += [[((x0+frontSponge)+x0)/2., ((y0+leftSponge)+(y1-rightSponge))/2., (z0+z1)/2.]]
            self.regionIndice[r_i] = 'frontSponge'
            regionFlags += [r_i+1]
            v_i += 2
            r_i += 1
            nb_sponge += 1
        # all flags as bottom flags
        for i in range(nb_sponge):
            vertexFlags += [bt['bottom'], bt['bottom']]
            segmentFlags += [bt['bottom'], bt['bottom'], bt['bottom']]
            facetFlags += [bt['bottom']]
        nb_corner = 0  # number of new corners
        if leftSponge and backSponge:
            vertices += [[x1, y0, z0]]
            segments += [[5+nb_corner, v_i], [v_i, 6+nb_corner]]
            facets += [[[1, 5+nb_corner, v_i, 6+nb_corner]]]
            regions += [[(x1+(x1-backSponge))/2., (y0+(y0+leftSponge))/2., (z0+z1)/2.]]
            self.regionIndice[r_i] = 'back_left_Sponge'
            regionFlags += [r_i+1]
            nb_corner += 2
            v_i += 1
            r_i += 1
        if backSponge and rightSponge:
            vertices += [[x1, y1, z0]]
            segments += [[5+nb_corner, v_i], [v_i, 6+nb_corner]]
            facets += [[[2, 5+nb_corner, v_i, 6+nb_corner]]]
            regions += [[(x1+(x1-backSponge))/2., (y1+(y1-rightSponge))/2., (z0+z1)/2.]]
            self.regionIndice[r_i] = 'back_right_Sponge'
            regionFlags += [r_i+1]
            nb_corner += 2
            v_i += 1
            r_i += 1
        if rightSponge and frontSponge:
            vertices += [[x0, y1, z0]]
            segments += [[5+nb_corner, v_i], [v_i, 6+nb_corner]]
            facets += [[[3, 5+nb_corner, v_i, 6+nb_corner]]]
            regions += [[(x0+(x0+frontSponge))/2., (y1+(y1-rightSponge))/2., (z0+z1)/2.]]
            self.regionIndice[r_i] = 'front_right_Sponge'
            regionFlags += [r_i+1]
            nb_corner += 2
            v_i += 1
            r_i += 1
        if frontSponge and leftSponge:
            vertices += [[x0, y0, z0]]
            segments += [[5+nb_corner, v_i], [v_i, 4]]
            facets += [[[0, 5+nb_corner, v_i, 4]]]
            regions += [[(x0+(x0+frontSponge))/2., (y0+(y0+leftSponge))/2., (z0+z1)/2.]]
            self.regionIndice[r_i] = 'front_left_Sponge'
            regionFlags += [r_i+1]
            nb_corner += 2
            v_i += 1
            r_i += 1
        for i in range(nb_corner/2):
            vertexFlags += [bt['bottom']]
            segmentFlags += [bt['bottom'], bt['bottom'], bt['bottom']]
            facetFlags += [bt['bottom']]
        # ---------------------------------------------
        # Then add the rest of the vectors (top) by symmetry
        # ---------------------------------------------
        # copying list of bottom segments to get top and side segments
        segments_bottom = segments[:]
        facets_bottom = facets[:]
        # getting top
        vertexFlags += [bt['top'] for i in range(len(vertices))]
        segmentFlags += [bt['top'] for i in range(len(segments))]
        facetFlags += [bt['top'] for i in range(len(facets))]
        vertices_top = np.array(vertices)
        vertices_top[:,2] = z1
        vertices += vertices_top.tolist()
        segments_top = np.array(segments)
        segments_top += v_i
        segments += segments_top.tolist()
        facets_top = np.array(facets)
        facets_top += v_i
        facets += facets_top.tolist()
        facets_side = []
        # getting sides
        for s in segments_bottom:  # for vertical facets
            facets += [[[s[0], s[1], s[1]+v_i, s[0]+v_i]]]
            if vertices[s[0]][0] == vertices[s[1]][0] == x0:
                facetFlags += [bt['front']]
            elif vertices[s[0]][0] == vertices[s[1]][0] == x1:
                facetFlags += [bt['back']]
            elif vertices[s[0]][1] ==  vertices[s[1]][1] == y0:
                facetFlags += [bt['left']]
            elif vertices[s[0]][1] == vertices[s[1]][1] == y1:
                facetFlags += [bt['right']]
            else:
                facetFlags += [bt['sponge']]
        for i in range(v_i):  # for vertical segments
            segments += [[i, i+v_i]]
            if vertices[i][0] == vertices[i+v_i][0] == x0:
                segmentFlags += [bt['front']]
            elif vertices[i][0] == vertices[i+v_i][0] == x1:
                segmentFlags += [bt['back']]
            elif vertices[i][1] == vertices[i+v_i][1] == y0:
                segmentFlags += [bt['left']]
            elif vertices[i][1] == vertices[i+v_i][1] == y1:
                segmentFlags += [bt['right']]
            else:
                segmentFlags += [bt['sponge']]
        self.vertices = np.array(vertices)
        self.vertices = np.dot(self.vertices, self.coords_system)
        self.vertexFlags = np.array(vertexFlags)
        self.segments = np.array(segments)
        self.segmentFlags = np.array(segmentFlags)
        self.facets = np.array(facets)
        self.facetFlags = np.array(facetFlags)
        self.regions = np.array(regions)
        self.regionFlags = np.array(regionFlags)
        if self._snbc is not None:
            self._resetShape()





    def setAbsorptionZones(self, allSponge=False, left=False, right=False, front=False,
                           back=False, front_left=False, front_right=False,
                           back_left=False, back_right=False):
        self.abs_zones = {'leftSponge': left,
                          'rightSponge': right,
                          'frontSponge': front,
                          'backSponge': back,
                          'front_left_Sponge': front_left,
                          'front_right_Sponge': front_right,
                          'back_left_Sponge': back_left,
                          'back_right_Sponge': back_right}
        if allSponge is True:
            for key in self.abs_zones:
                self.abs_zones[key] = True
        for key, val in self.abs_zones.iteritems():
            if val is True:
                ind = self.regionIndice[key]
                key = ind + (self._snr + 1)
                self.zones[key] = RelaxationZone(self.regions[ind, 0],
                                                 1.,
                                                 lambda x, t: 0.,
                                                 lambda x, t: 0.,
                                                 lambda x, t: 0.,)
                self.dragAlphaTypes[ind] = 0.5/1.005e-6
                self.epsFact_solid[ind] = self.leftSponge/2.


class Tank2D(Rectangle):
    """
    Class to create a 2D tank (rectangle)
    :arg domain: domain of the tank
    :arg dim: dimensions of the tank (list or array)
    :leftSponge: width of left sponge (float)
    :rightSponge: width of right sponge (float)
    """
    count = 0
    def __init__(self, domain, dim=(0.,0.), leftSponge=None, rightSponge=None,
                 from_0=True):
        L, H = dim
        if from_0 is True:
            x, y = L/2., H/2.
        else:
            x, y = 0., 0.
        Rectangle.__init__(self, domain, dim=dim, coords=(x, y), add=False)
        self.__class__.count += 1
        self.name = "tank2d" + str(self.__class__.count)
        self.from_0 = from_0
        self.regions = np.array([[L/2., H/2.]])
        for boundcond in self.BC_list:
            boundcond.setTank()
        extra_vertices = []
        extra_vertexFlags = []
        extra_regions = []
        extra_regionFlags = []
        self.regionIndice = {'tank': 0}
        ind_region = 1
        x0 = x-0.5*L
        x1 = x+0.5*L
        y0 = y-0.5*H
        y1 = y+0.5*H
        self.leftSponge = leftSponge
        self.rightSponge = rightSponge
        if leftSponge is not None:
            extra_vertices += [[x0+leftSponge, y0], [x0+leftSponge, y1]]
            extra_vertexFlags += [1, 3]
            extra_regions += [[(x0+leftSponge)/2., (y0+y1)/2.]]
            self.regionIndice['leftSponge'] = ind_region
            ind_region += 1
            extra_regionFlags += [ind_region]
            self.BC_list += [bc.BoundaryConditions()]
        if rightSponge is not None:
            extra_vertices += [[x1-rightSponge, y0], [x1-rightSponge, y1]]
            extra_vertexFlags += [1, 3]
            extra_regions += [[((x1-rightSponge)+x1)/2., (y0+y1)/2.]]
            self.regionIndice['rightSponge'] = ind_region
            ind_region += 1
            extra_regionFlags += [ind_region]
            self.BC_list += [bc.BoundaryConditions()]
        # getting the right segments if sponge layers are defined
        if leftSponge is not None and rightSponge is not None:
            self.segments = np.array([[0, 4], [4, 6], [6, 1], [1, 2], [2, 7],
                                      [7, 5], [5, 3], [3, 0], [4, 5], [6, 7]])
            self.segmentFlags = np.array([1, 1, 1, 2, 3, 3, 3, 4, 5, 6])
            self.regions[0] = [(x0+leftSponge)+(x1-rightSponge)/2., (y0+y1)/2.]
        elif leftSponge is not None or rightSponge is not None:
            self.segments = np.array([[0, 4], [4, 1], [1, 2], [2, 5], [5, 3],
                                      [3, 0], [4, 5]])
            self.segmentFlags = np.array([1, 1, 2, 3, 3, 4, 5])
            if leftSponge is not None:
                self.regions[0] = [((x0+leftSponge)+x1)/2., (y0+y1)/2.]
            elif rightSponge is not None:
                self.regions[0] = [(x0+(x1-rightSponge))/2., (y0+y1)/2.]
        # need to check that original region is not in new sponge regions!
        if len(extra_vertices):
            self.vertices = np.append(self.vertices, extra_vertices, axis=0)
            self.vertexFlags = np.append(self.vertexFlags, extra_vertexFlags, axis=0)
        if len(extra_regions):
            self.regions = np.append(self.regions, extra_regions, axis=0)
            self.regionFlags = np.append(self.regionFlags, extra_regionFlags, axis=0)
        self.porosityTypes = np.ones(len(self.regionFlags))
        self.dragAlphaTypes = np.zeros(len(self.regionFlags))
        self.dragBetaTypes = np.zeros(len(self.regionFlags))
        self.epsFact_solid = np.zeros(len(self.regionFlags))
        self.zones = {}
        self.RelaxationZones = RelaxationZoneWaveGenerator(self.zones, shape=self)
        self.barycenter = np.array([0., 0., 0.])
        self.barycenters = np.array([self.barycenter for i in range(max(self.segmentFlags))])
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
            self.coords[:] = [x, y]
        else:
            x, y = self.coords
        x0, x1 = x-0.5*L, x+0.5*L
        y0, y1 = y-0.5*H, y+0.5*H
        vertices = [[x0, y0],
                    [x1, y0],
                    [x1, y1],
                    [x0, y1]]
        leftSponge = self.leftSponge
        rightSponge = self.rightSponge
        if leftSponge is not None:
            vertices += [[x0+leftSponge, y0],
                         [x0+leftSponge, y1]]
            regions = [[((x1-leftSponge)+x1)/2., (y0+y1)/2.],
                       [(x0+leftSponge)/2., (y0+y1)/2.]]
        if rightSponge is not None:
            vertices += [[x1-rightSponge, y0],
                         [x1-rightSponge, y1]]
            regions = [[(x0+rightSponge)/2., (y0+y1)/2.],
                       [((x1-rightSponge)+x1)/2., (y0+y1)/2.]]
        if rightSponge is not None and leftSponge is not None:
            regions = [[((x0+leftSponge)+(x1-rightSponge))/2., (y0+y1)/2.],
                       [(x0+leftSponge)/2., (y0+y1)/2.],
                       [((x1-rightSponge)+x1)/2., (y0+y1)/2.]]
        self.vertices[:] = vertices
        self.regions[:] = regions
        self._updateDomain()

    def setAbsorptionZones(self, left=False, right=False):
        self.leftSpongeAbs = left
        self.rightSpongeAbs = right
        if self.leftSpongeAbs is True:
            ind = self.regionIndice['leftSponge']
            key = ind + self._snr + 1
            self.zones[key] = RelaxationZone(self.regions[ind, 0],
                                             1.,
                                             lambda x, t: 0.,
                                             lambda x, t: 0.,
                                             lambda x, t: 0.,)
            self.dragAlphaTypes[ind] = 0.5/1.005e-6
            self.epsFact_solid[ind] = self.leftSponge/2.
        if self.rightSpongeAbs is True:
            ind = self.regionIndice['rightSponge']
            key = ind + self._snr + 1
            self.zones[key] = RelaxationZone(self.regions[ind, 0],
                                             1.,
                                             lambda x, t: 0.,
                                             lambda x, t: 0.,
                                             lambda x, t: 0.,)
            self.dragAlphaTypes[ind] = 0.5/1.005e-6
            self.epsFact_solid[ind] = self.rightSponge/2.


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

    def __init__(self, domain, barycenter=None, vertices=None,
                 vertexFlags=None, segments=None, segmentFlags=None,
                 facets=None, facetFlags=None, holes=None, regions=None,
                 regionFlags=None, boundaryTags=None, boundaryOrientations=None):
        Shape.__init__(self, domain)
        self.__class__.count += 1
        self.name = "custom" + str(self.__class__.count)
        flagSet = set()
        for tag, value in boundaryTags.iteritems():
            flagSet.add(value)
        print flagSet
        minFlag = min(flagSet)
        previous_flag = minFlag-1
        for flag in flagSet:
            print flag
            assert flag == previous_flag+1, "Flags must be defined as a suite of numbers (e.g. 0, 1, 2, 3, 4 with no gap)!"
            previous_flag = flag
        self.vertices = np.array(vertices)
        self.vertexFlags = np.array(vertexFlags)-(minFlag+1)
        if segments:
            self.segments = np.array(segments)
            self.segmentFlags = np.array(segmentFlags)-(minFlag+1)
        if facets:
            self.facets = np.array(facets)
            self.facetFlags = np.array(facetFlags)-(minFlag+1)
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
        if boundaryOrientations is not None:
            b_or = []
        else:
            b_or = None
            b_i = None
        for tag, index in boundaryTags.iteritems():
            if boundaryOrientations is not None:
                b_or += [boundaryOrientations[tag]]
                b_i = index-minFlag
            self.BC_dict[tag] = bc.BoundaryConditions(b_or=b_or, b_i=b_i)
            self.BC_list[index-minFlag] = self.BC_dict[tag]
        self.BC = BCContainer(self.BC_dict)
        self._addShape()

    def _setInertiaTensor(self, It):
        self.It = np.array(It)


class RigidBody(AuxiliaryVariables.AV_base):

    def __init__(self, shape, he=1., cfl_target=0.9, dt_init=0.001):
        self.shape = shape
        shape.domain.auxiliaryVariables += [self]
        self.dt_init = dt_init
        self.he = he
        self.cfl_target = 0.9
        self.last_position = np.array([0., 0., 0.])
        self.rotation_matrix = np.eye(3)
        self.h = np.array([0., 0., 0.])

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
            self.inertia = self.shape.getInertia(vec=self.M, pivot=self.shape.barycenter)
            assert self.inertia != 0, "Zero inertia: inertia tensor (It) was not set correctly!"
            ang_acc = self.M[:]/self.inertia
        else:
            self.inertia = None
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
        if self.shape.record_values is True:
            self.recordValues()

    def recordValues(self):
        time = self.model.stepController.t_model_last-self.model.levelModelList[-1].dt_last
        self.record_time = time
        pos_x, pos_y, pos_z = self.last_position
        rot = self.last_rotation
        rot_x = atan2(rot[2, 1], rot[1, 2])
        rot_y = atan2(-rot[0, 2], sqrt(rot[2, 1]**2+rot[2, 2]**2))
        rot_z = atan2(rot[1, 0], rot[0, 0])
        Fx, Fy, Fz = self.F
        Mx, My, Mz = self.M
        inertia = self.inertia
        vel_x, vel_y, vel_z = self.velocity
        acc_x, acc_y, acc_z = self.acceleration
        values = [time, pos_x, pos_y, pos_z, rot_x, rot_y,
                  rot_z, Fx, Fy, Fz, Mx, My, Mz, inertia,
                  vel_x, vel_y, vel_z, acc_x, acc_y, acc_z]
        values_towrite = list(compress(values, self.shape.record_bool))
        comm = Comm.get()
        if  comm.isMaster():
            if self.shape.record_values is True:
                with open(self.record_file, 'a') as csvfile:
                    writer = csv.writer(csvfile, delimiter=',')
                    writer.writerow(values_towrite)

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
        self.acceleration = np.zeros(3, 'd')
        self.last_acceleration = np.zeros(3, 'd')
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
        self.nb_start = self.shape._snbc
        self.nb_end = self.shape._enbc
        if nd == 2:
            self.Fg = self.shape.mass*np.array([0., -9.81, 0.])
        if nd == 3:
            self.Fg = self.shape.mass*np.array([0., 0., -9.81])
        comm = Comm.get()
        if  comm.isMaster():
            if self.shape.record_values is True:
                self.record_file = os.path.join(Profiling.logDir, self.shape.record_filename)
                with open(self.record_file, 'w') as csvfile:
                    writer = csv.writer(csvfile, delimiter=',')
                    writer.writerow(self.shape.record_names)

    def calculate(self):
        self.last_position[:] = self.position
        self.last_velocity[:] = self.velocity
        self.last_acceleration[:] = self.acceleration
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
        self.F[:] = F*self.shape.free_x
        self.M[:] = M*self.shape.free_r
        logEvent("========================================================================")
        logEvent("======================= Rigid Body Calculation =========================")
        logEvent("========================================================================")
        logEvent("Name: " + `self.shape.name`)
        logEvent("========================================================================")
        logEvent("[proteus]     t=%1.5fsec to t=%1.5fsec" % (self.model.stepController.t_model_last-dt, self.model.stepController.t_model_last))
        logEvent("[proteus]    dt=%1.5fsec" % (dt))
        logEvent("[body] ================== Pre-calculation attributes  ==================")
        logEvent("[proteus]     t=%1.5fsec" % (self.model.stepController.t_model_last-dt))
        logEvent("[proteus]     F=(% 21.16e, % 21.16e, % 21.16e)" % (F[0], F[1], F[2]))
        logEvent("[proteus] F*DOF=(% 21.16e, % 21.16e, % 21.16e)" % (self.F[0], self.F[1], self.F[2]))
        logEvent("[proteus]     M=(% 21.16e, % 21.16e, % 21.16e)" % (M[0], M[1], M[2]))
        logEvent("[proteus] M*DOF=(% 21.16e, % 21.16e, % 21.16e)" % (self.M[0], self.M[1], self.M[2]))
        logEvent("[body]        h=(% 21.16e, % 21.16e, % 21.16e)" % (self.h[0], self.h[1], self.h[2]))
        logEvent("[body]      pos=(% 21.16e, % 21.16e, % 21.16e)" % (self.last_position[0], self.last_position[1], self.last_position[2]))
        logEvent("[body]      vel=(% 21.16e, % 21.16e, % 21.16e)" % (self.last_velocity[0], self.last_velocity[1], self.last_velocity[2]))
        self.step(dt)
        logEvent("[body] ================== Post-calculation attributes ==================")
        logEvent("[body]        t=%1.5fsec" % (self.model.stepController.t_model_last))
        logEvent("[body]      pos=(% 21.16e, % 21.16e, % 21.16e)" % (self.position[0], self.position[1], self.position[2]))
        logEvent("[body]      vel=(% 21.16e, % 21.16e, % 21.16e)" % (self.velocity[0], self.velocity[1], self.velocity[2]))
        logEvent("[body]    r vel=(% 21.16e, % 21.16e, % 21.16e)" % (self.angvel[0], self.angvel[1], self.angvel[2]))
        if sum(self.angvel) != 0:
            rot_axis=self.angvel/np.linalg.norm(self.angvel)
        else:
            rot_axis=(0., 0., 0.)
        rot = self.rotation
        rot_x = -atan2(rot[2, 1], rot[1, 2])
        rot_y = -atan2(-rot[0, 2], sqrt(rot[2, 1]**2+rot[2, 2]**2))
        rot_z = -atan2(rot[1, 0], rot[0, 0])
        logEvent("[body]   r axis=(% 21.16e, % 21.16e, % 21.16e)" % (rot_axis[0], rot_axis[1], rot_axis[2]))
        logEvent("[body]    r ang=(% 21.16e)" % (self.ang))
        logEvent("[body]      rot=(% 21.16e, % 21.16e, % 21.16e)" % (rot_x, rot_y, rot_z))
        logEvent("========================================================================")







# ------------------------------------------------------------------------------ #
# --------------------------SPATIAL TOOLS FOR SHAPES---------------------------- #
# ------------------------------------------------------------------------------ #


def rotation2D(points, rot, pivot=(0., 0.)):
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





class RelaxationZone:
    def __init__(self, center_x, sign, u, v, w):
        self.center_x = center_x
        self.sign = sign
        self.u = u
        self.v = v
        self.w = w


class RelaxationZoneWaveGenerator(AuxiliaryVariables.AV_base):
    """ Prescribe a velocity penalty scaling in a material zone via a Darcy-Forchheimer penalty

    :param zones: A dictionary mapping integer material types to Zones, where a Zone is a named tuple
    specifying the x coordinate of the zone center and the velocity components
    """
    def __init__(self, zones, shape):
        assert isinstance(zones,dict)
        self.zones = zones
        self.shape = shape
        shape.domain.auxiliaryVariables += [self]
    def calculate(self):
        for l, m in enumerate(self.model.levelModelList):
            for eN in range(m.coefficients.q_phi.shape[0]):
                mType = m.mesh.elementMaterialTypes[eN]
                if self.zones.has_key(mType):
                    for k in range(m.coefficients.q_phi.shape[1]):
                        t = m.timeIntegration.t
                        x = m.q['x'][eN,k]
                        m.coefficients.q_phi_solid[eN,k] = self.zones[mType].sign*(self.zones[mType].center_x - x[0])
                        m.coefficients.q_velocity_solid[eN,k,0] = self.zones[mType].u(x,t)
                        m.coefficients.q_velocity_solid[eN,k,1] = self.zones[mType].v(x,t)
                        if self.shape.domain.nd > 2:
                            m.coefficients.q_velocity_solid[eN,k,2] = self.zones[mType].w(x,t)
        m.q['phi_solid'] = m.coefficients.q_phi_solid
        m.q['velocity_solid'] = m.coefficients.q_velocity_solid

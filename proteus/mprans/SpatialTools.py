"""
This module adds functionality to proteus.SpatialTools module by enabling
two-phase flow functionality such as converting shapes to moving rigid bodies,
or adding wave absorption and generation zones.


Example
-------
from proteus import Domain
from proteus.mprans import SpatialTools as st
import numpy as np

domain = Domain.PlanarStraightLineGraphDomain()
tank = st.Tank2D(domain. dim=[4., 4.])
tank.setSponge(left=0.4)
tank.setAbsorptionZones(left=true)
shape = st.Rectangle(domain, dim=[0.5, 0.5], coords=[1., 1.])
shape.setRigidBody()
shape2.rotate(np.pi/3.)

st.buildDomain(domain)
"""

from math import cos, sin, sqrt, atan2, acos
from itertools import compress
import csv
import os
import numpy as np
from proteus import AuxiliaryVariables, Archiver, Comm, Profiling
from proteus.Profiling import logEvent as log
from proteus import BC as bc
from proteus.SpatialTools import (Shape, BCContainer, Cuboid, Rectangle,
                                  CustomShape)


class ShapeRANS(Shape):
    """
    Super class of shapes defined below. All shapes will have the arguments and
    functions defined here + the ones from Shape of proteus.SpatialTools
    """

    def __init__(self, domain, nd):
        super(ShapeRANS, self).__init__(domain, nd)
        self.mass = None
        self.density = None
        self.free_x = (1, 1, 1)
        self.free_r = (1, 1, 1)
        self.record_values = False
        self.zones = {}  # for absorption/generation zones
        self.auxiliaryVariables = {}  # list of auxvar attached to shape
        self.It = None  # inertia tensor

    def _attachAuxiliaryVariable(self, key):
        """
        Attaches an auxiliary variable to the auxiliaryVariables dictionary of
        the shape (used in buildDomain function)
        (!) should not be used manually
        """
        if key not in self.auxiliaryVariables:
            if key == 'RigidBody':
                self.auxiliaryVariables[key] = True
            if key == 'RelaxZones':
                self.auxiliaryVariables[key] = self.zones

    def setRigidBody(self, holes=None):
        """
        Makes the shape a rigid body

        :param holes: set of hole coordinates (the interior of the rigid body
                      will not be meshed) (list/array)
        """
        self._attachAuxiliaryVariable('RigidBody')
        if holes is None:
            self.holes = np.array([self.barycenter[:self.nd]])
        else:
            self._checkListOfLists(holes)
            self.holes = np.array(holes)

    def setTank(self):
        """
        Sets tank boundary conditions (for moving domain).
        """
        for boundcond in self.BC_list:
            boundcond.setTank()

    def setConstraints(self, free_x, free_r):
        """
        Sets constraints on the Shape (for moving bodies)

        :param free_x: translational constraints
        :param free_r: rotational constraints
        """
        self.free_x = np.array(free_x)
        self.free_r = np.array(free_r)

    def setMass(self, mass):
        """
        Set mass of the shape and calculate density

        :param mass: mass (int/float)
        """
        self.mass = float(mass)
        if self.volume:
            self.density = self.mass/self.volume

    def setDensity(self, density):
        """
        Set density and calculate mass

        :param density: density (int/float)
        """
        self.density = float(density)
        if self.volume:
            self.mass = self.density*self.volume

    def _setInertiaTensor(self, It):
        """
        Set the inertia tensor of the shape

        :param It: inertia tensor (list/array)
        """
        It = np.array(It)
        if self.nd == 2:
            assert isinstance(It, float), 'the inertia tensor of a 2D shape ' \
                'must be a float'
        if self.nd == 3:
            assert It.shape == (3, 3), 'the inertia tensor of a 3D shape ' \
                'must have a (3, 3) shape'
        self.It = It

    def getInertia(self, vec=(0., 0., 1.), pivot=None):
        """
        Gives the inertia of the shape from an axis and a pivot
        (!) The inertia tensor of the shape must be set
        """
        assert self.It is not None, 'No inertia tensor! (' + self.name + ')'
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
            # getting the tensor for calculaing moment of inertia
            # from arbitrary axis
            vt = np.array([[cx**2, cx*cy, cx*cz],
                           [cx*cy, cy**2, cy*cz],
                           [cx*cz, cy*cz, cz**2]])
            # total moment of inertia
            I = np.einsum('ij,ij->', self.mass*self.It, vt)
        return I

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
                            rot_z, F, Fx, Fy, Fz, M, Mx, My, Mz, inertia,
                            vel_x, vel_y, vel_z, acc_x, acc_y, acc_z]
        if all_values is True:
            self.record_bool = [True for value in self.record_bool]
        self.record_names = ['time', 'pos_x', 'pos_y', 'pos_z', 'rot_x',
                             'rot_y', 'rot_z', 'Fx', 'Fy', 'Fz', 'Mx', 'My',
                             'Mz', 'inertia', 'vel_x', 'vel_y', 'vel_z',
                             'acc_x', 'acc_y', 'acc_z']
        self.record_names = list(compress(self.record_names, self.record_bool))
        if filename is None:
            self.record_filename = 'record_' + self.name + '.csv'
        else:
            self.record_filename = filename + '.csv'

    def setAbsorptionZones(self, flags, epsFact_solid, sign=None,
                           center_x=None, dragAlphaTypes=0.5/1.005e-6,
                           dragBetaTypes=0., porosityTypes=1.):
        """
        Sets a region (given the local index) to an absorption zone

        :param flags: local flags of the region. Can be an integer or a list.
        :param epsFact_solid: half of absorption zone length
        :param sign: direction vector from the boundary to the sponge (x-axis)
        :param center_x: center of the absorption zone
        :param dragAlpha
        """
        self._attachAuxiliaryVariable('RelaxZones')
        waves = None
        windSpeed = 0.
        if isinstance(flags, int):
            flags = [flags]
            sign = [sign]
            epsFact_solid = [epsFact_solid]
            center_x = [center_x]
            dragAlphaTypes = [dragAlphaTypes]
            dragBetaTypes = [dragBetaTypes]
            porosityTypes = [porosityTypes]
        for i, flag in enumerate(flags):
            self.zones[flag] = bc.RelaxationZone(domain=self.domain,
                                                 zone_type='absorption',
                                                 sign=sign[i],
                                                 center_x=center_x[i],
                                                 waves=waves,
                                                 windSpeed=windSpeed,
                                                 epsFact_solid=epsFact_solid[i],
                                                 dragAlphaTypes=dragAlphaTypes[i],
                                                 dragBetaTypes=dragBetaTypes[i],
                                                 porosityTypes=porosityTypes[i])

    def setGenerationZones(self, flags, epsFact_solid, sign, center_x, waves,
                           windSpeed=(0., 0., 0.), dragAlphaTypes=0.5/1.005e-6,
                           dragBetaTypes=0., porosityTypes=1.):
        """
        Sets a region (given the local index) to a generation zone

        :param flags: local flags of the region. Can be an integer or a list.
        :param epsFact_solid: absorption zone length (length of region/2)
        """
        self._attachAuxiliaryVariable('RelaxZones')
        if isinstance(flags, int):
            flags = [flags]
            sign = [sign]
            waves = [waves]
            windSpeed = [windSpeed]
            epsFact_solid = [epsFact_solid]
            center_x = [center_x]
            dragAlphaTypes = [dragAlphaTypes]
            dragBetaTypes = [dragBetaTypes]
            porosityTypes = [porosityTypes]
        for i, flag in enumerate(flags):
            self.zones[flag] = bc.RelaxationZone(domain=self.domain,
                                                 zone_type='generation',
                                                 sign=sign[i],
                                                 center_x=center_x[i],
                                                 waves=waves[i],
                                                 windSpeed=windSpeed[i],
                                                 epsFact_solid=epsFact_solid[i],
                                                 dragAlphaTypes=dragAlphaTypes[i],
                                                 dragBetaTypes=dragBetaTypes[i],
                                                 porosityTypes=porosityTypes[i])


# -----------------------------------------------------------------------------
# ADDING FUNCTIONALITY TO SHAPE FROM proteus.SpatialTools
# -----------------------------------------------------------------------------

# reassigning base/super class to access all functions from ShapeRANS and Shape
Rectangle.__bases__ = (ShapeRANS,)
Cuboid.__bases__ = (ShapeRANS,)
CustomShape.__bases__ = (ShapeRANS,)

# adding extra functionality to predefined shapes

def _CuboidsetInertiaTensor(self):
    """
    Sets the inertia tensor of the cuboid
    (!) should not be used manually
    """
    L, W, H = self.dim
    self.It = [[(W**2.+H**2.)/12., 0, 0],
               [0, (L**2.+H**2.)/12., 0],
               [0, 0, (W**2.+L**2.)/12.]]

Cuboid._setInertiaTensor = _CuboidsetInertiaTensor

def _RectanglesetInertiaTensor(self):
    """
    Sets the inertia tensor of the rectangle
    (!) should not be used manually
    """
    L, H = self.dim
    self.It = (L**2+H**2)/12

Rectangle._setInertiaTensor = _RectanglesetInertiaTensor


# -----------------------------------------------------------------------------
# DEFINING NEW SHAPES TYPES
# -----------------------------------------------------------------------------

class Tank3D(ShapeRANS):
    """
    Class to create a 3D tank (cuboid)

    :param domain: domain of the tank
    :param dim: dimensions of the tank (list or array)
    """
    count = 0

    def __init__(self, domain, dim=(0., 0., 0.), from_0=True, leftSponge=None,
                 rightSponge=None, frontSponge=None, backSponge=None):
        super(Tank3D, self).__init__(domain, nd=3)
        self.__class__.count += 1
        self.name = "tank3d" + str(self.__class__.count)
        self.from_0 = from_0
        self.holes = None
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
        for i in range(6):
            self.BC_list[i].setTank()
        self.barycenter = np.array([0., 0., 0.])
        self.setDimensions(dim)
        self.porosityTypes = np.ones(len(self.regionFlags))
        self.dragAlphaTypes = np.zeros(len(self.regionFlags))
        self.dragBetaTypes = np.zeros(len(self.regionFlags))
        self.epsFact_solid = np.zeros(len(self.regionFlags))

    def setSponge(self, left=None, right=None, back=None, front=None):
        """
        Set length of sponge layers of the tank.
        (!) Sponge layers expand inwards.

        :param left: length of the left sponge (int/float)
        :param right: length of the right sponge (int/float)
        :param back: length of the back sponge (int/float)
        :param front: length of the front sponge (int/float)
        """
        self.leftSponge = left
        self.rightSponge = right
        self.backSponge = back
        self.frontSponge = front
        self.setDimensions(self.dim)

    def setDimensions(self, dim):
        """
        Set dimension of the tank

        :param dim: new dimensions (list/array)
        """
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
        rightSponge = self.rightSponge or 0.
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
        regions = [[((x0+frontSponge)+(x1-backSponge))/2.,
                    ((y0+leftSponge)+(y1-rightSponge))/2.,
                    (z0+z1)/2.]]
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
            regions += [[((x0+frontSponge)+(x1-backSponge))/2.,
                         (y0+(y0+leftSponge))/2.,
                         (z0+z1)/2.]]
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
            regions += [[((x1-backSponge)+x1)/2.,
                         ((y0+leftSponge)+(y1-rightSponge))/2.,
                         (z0+z1)/2.]]
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
            regions += [[((x0+frontSponge)+(x1-backSponge))/2.,
                         (y1+(y1-rightSponge))/2.,
                         (z0+z1)/2.]]
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
            regions += [[((x0+frontSponge)+x0)/2.,
                         ((y0+leftSponge)+(y1-rightSponge))/2.,
                         (z0+z1)/2.]]
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
            regions += [[(x1+(x1-backSponge))/2.,
                         (y0+(y0+leftSponge))/2.,
                         (z0+z1)/2.]]
            self.regionIndice[r_i] = 'back_left_Sponge'
            regionFlags += [r_i+1]
            nb_corner += 2
            v_i += 1
            r_i += 1
        if backSponge and rightSponge:
            vertices += [[x1, y1, z0]]
            segments += [[5+nb_corner, v_i], [v_i, 6+nb_corner]]
            facets += [[[2, 5+nb_corner, v_i, 6+nb_corner]]]
            regions += [[(x1+(x1-backSponge))/2.,
                         (y1+(y1-rightSponge))/2.,
                         (z0+z1)/2.]]
            self.regionIndice[r_i] = 'back_right_Sponge'
            regionFlags += [r_i+1]
            nb_corner += 2
            v_i += 1
            r_i += 1
        if rightSponge and frontSponge:
            vertices += [[x0, y1, z0]]
            segments += [[5+nb_corner, v_i], [v_i, 6+nb_corner]]
            facets += [[[3, 5+nb_corner, v_i, 6+nb_corner]]]
            regions += [[(x0+(x0+frontSponge))/2.,
                         (y1+(y1-rightSponge))/2.,
                         (z0+z1)/2.]]
            self.regionIndice[r_i] = 'front_right_Sponge'
            regionFlags += [r_i+1]
            nb_corner += 2
            v_i += 1
            r_i += 1
        if frontSponge and leftSponge:
            vertices += [[x0, y0, z0]]
            segments += [[5+nb_corner, v_i], [v_i, 4]]
            facets += [[[0, 5+nb_corner, v_i, 4]]]
            regions += [[(x0+(x0+frontSponge))/2.,
                         (y0+(y0+leftSponge))/2.,
                         (z0+z1)/2.]]
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
        # getting top
        vertexFlags += [bt['top'] for i in range(len(vertices))]
        segmentFlags += [bt['top'] for i in range(len(segments))]
        facetFlags += [bt['top'] for i in range(len(facets))]
        vertices_top = np.array(vertices)
        vertices_top[:, 2] = z1
        vertices += vertices_top.tolist()
        segments_top = np.array(segments)
        segments_top += v_i
        segments += segments_top.tolist()
        facets_top = np.array(facets)
        facets_top += v_i
        facets += facets_top.tolist()
        # getting sides
        for s in segments_bottom:  # for vertical facets
            facets += [[[s[0], s[1], s[1]+v_i, s[0]+v_i]]]
            if vertices[s[0]][0] == vertices[s[1]][0] == x0:
                facetFlags += [bt['front']]
            elif vertices[s[0]][0] == vertices[s[1]][0] == x1:
                facetFlags += [bt['back']]
            elif vertices[s[0]][1] == vertices[s[1]][1] == y0:
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

    # def setAbsorptionZones(self, epsFact_solid, allSponge=False, left=False,
    #                        right=False, front=False, back=False,
    #                        front_left=False, front_right=False,
    #                        back_left=False, back_right=False,
    #                        dragAlphaTypes=0.5/1.005e-6, dragBetaTypes=0.,
    #                        porosityTypes=1.):
    #     self.abs_zones = {'leftSponge': left,
    #                       'rightSponge': right,
    #                       'frontSponge': front,
    #                       'backSponge': back,
    #                       'front_left_Sponge': front_left,
    #                       'front_right_Sponge': front_right,
    #                       'back_left_Sponge': back_left,
    #                       'back_right_Sponge': back_right}
    #     if allSponge is True:
    #         for key in self.abs_zones:
    #             self.abs_zones[key] = True
    #     indice = []
    #     for key, value in self.abs_zones:
    #         if value is True:
    #             indice += [self.regionIndice[key]]
    #     print('Tank3D absorption zones not implemented yet!')
    #     sys.exit()


class Tank2D(ShapeRANS):
    """
    Class to create a 2D tank (rectangle)

    :param domain: domain of the tank
    :param dim: dimensions of the tank (list or array)
    :param leftSponge: width of left sponge (float)
    :param rightSponge: width of right sponge (float)
    """
    count = 0

    def __init__(self, domain, dim=(0., 0.), leftSponge=None, rightSponge=None,
                 from_0=True):
        super(Tank2D, self).__init__(domain, nd=2)
        self.__class__.count += 1
        self.name = "tank2d" + str(self.__class__.count)
        self.from_0 = from_0
        self.leftSponge = leftSponge
        self.rightSponge = rightSponge
        self.leftAbs = False
        self.rightAbs = False
        self.boundaryTags = {'bottom': 1,
                             'right': 2,
                             'top': 3,
                             'left': 4,
                             'sponge': 5}
        self.b_or = np.array([[0., -1.],
                              [1., 0.],
                              [0., 1.],
                              [-1., 0.]])
        self.BC_dict = {'bottom': bc.BoundaryConditions(b_or=self.b_or, b_i=0),
                        'right': bc.BoundaryConditions(b_or=self.b_or, b_i=1),
                        'top': bc.BoundaryConditions(b_or=self.b_or, b_i=2),
                        'left': bc.BoundaryConditions(b_or=self.b_or, b_i=3),
                        'sponge': bc.BoundaryConditions()}
        self.BC_list = [self.BC_dict['bottom'],
                        self.BC_dict['right'],
                        self.BC_dict['top'],
                        self.BC_dict['left'],
                        self.BC_dict['sponge']]
        self.BC = BCContainer(self.BC_dict)
        for i in range(4):
            self.BC_list[i].setTank()
        self.setDimensions(dim)
        self.porosityTypes = np.ones(len(self.regionFlags))
        self.dragAlphaTypes = np.zeros(len(self.regionFlags))
        self.dragBetaTypes = np.zeros(len(self.regionFlags))
        self.epsFact_solid = np.zeros(len(self.regionFlags))

    def setDimensions(self, dim):
        """
        Set dimension of the tank

        :param dim: new dimensions (list/array)
        """
        self.dim = dim
        L, H = dim
        if self.from_0 is True:
            x, y = L/2., H/2.
        else:
            x, y = 0., 0.
        self.coords = [x, y]
        x0, x1 = x-0.5*L, x+0.5*L
        y0, y1 = y-0.5*H, y+0.5*H
        bt = self.boundaryTags
        # add attributes
        leftSponge = self.leftSponge or 0.
        rightSponge = self.rightSponge or 0.
        regions_y = y1-y1/100.  # y coord of regions
        vertices = [[x-0.5*L, y-0.5*H],
                    [x+0.5*L, y-0.5*H],
                    [x+0.5*L, y+0.5*H],
                    [x-0.5*L, y+0.5*H]]
        vertexFlags = [bt['bottom'], bt['bottom'], bt['top'], bt['top']]
        segments = [[0, 1], [1, 2], [2, 3], [3, 0]]
        segmentFlags = [1, 2, 3, 4]  # bottom, right, top, left
        regions = [[(x0+leftSponge+x1-rightSponge)/2., regions_y]]
        regionFlags = [1]
        self.regionIndice = {'tank': 0}
        ind_region = 1
        if leftSponge:
            vertices += [[x0+leftSponge, y0], [x0+leftSponge, y1]]
            vertexFlags += [bt['bottom'], bt['top']]
            regions += [[(x0+leftSponge)/2., regions_y]]
            self.regionIndice['leftSponge'] = ind_region
            ind_region += 1
            regionFlags += [ind_region]
        if rightSponge:
            vertices += [[x1-rightSponge, y0], [x1-rightSponge, y1]]
            vertexFlags += [bt['bottom'], bt['top']]
            regions += [[((x1-rightSponge)+x1)/2., regions_y]]
            self.regionIndice['rightSponge'] = ind_region
            ind_region += 1
            regionFlags += [ind_region]
        # getting the right segments if sponge layers are defined
        if leftSponge and rightSponge:
            segments = [[0, 4], [4, 6], [6, 1], [1, 2], [2, 7],
                        [7, 5], [5, 3], [3, 0], [4, 5], [6, 7]]
            segmentFlags = [1, 1, 1, 2, 3, 3, 3, 4, 5, 5]
        elif leftSponge or rightSponge:
            segments = [[0, 4], [4, 1], [1, 2], [2, 5], [5, 3], [3, 0], [4, 5]]
            segmentFlags = [1, 1, 2, 3, 3, 4, 5]
        else:
            segments = [[0, 1], [1, 2], [2, 3], [3, 0]]
            segmentFlags = [1, 2, 3, 4]  # bottom, right, top, left
        # need to check that original region is not in new sponge regions!
        self.vertices = np.array(vertices)
        self.vertexFlags = np.array(vertexFlags)
        self.segments = np.array(segments)
        self.segmentFlags = np.array(segmentFlags)
        self.regions = np.array(regions)
        self.regionFlags = np.array(regionFlags)

    def setSponge(self, left=None, right=None):
        """
        Set length of sponge layers of the tank.
        (!) Sponge layers expand inwards.

        :param leftSponge: length of left sponge (float)
        :param rightSponge: length of right sponge (float)
        """
        self.leftSponge = left
        self.rightSponge = right
        self.setDimensions(self.dim)

    def setAbsorptionZones(self, left=False, right=False,
                           dragAlphaTypes=0.5/1.005e-6, dragBetaTypes=0.,
                           porosityTypes=1.):
        self.leftSpongeAbs = left
        self.rightSpongeAbs = right
        waves = None
        windSpeed = 0.
        if left or right:
            self._attachAuxiliaryVariable('RelaxZones')
        if self.leftSpongeAbs is True:
            ind = self.regionIndice['leftSponge']
            flag = self.regionFlags[ind]
            center_x = self.coords[0]-0.5*self.dim[0]+self.leftSponge/2.
            epsFact_solid = self.leftSponge/2.
            sign = 1.
            self.zones[flag] = bc.RelaxationZone(domain=self.domain,
                                                 zone_type='absorption',
                                                 sign=sign, center_x=center_x,
                                                 waves=waves,
                                                 windSpeed=windSpeed,
                                                 epsFact_solid=epsFact_solid,
                                                 dragAlphaTypes=dragAlphaTypes,
                                                 dragBetaTypes=dragBetaTypes,
                                                 porosityTypes=porosityTypes)
        if self.rightSpongeAbs is True:
            ind = self.regionIndice['rightSponge']
            flag = self.regionFlags[ind]
            center_x = self.coords[0]+0.5*self.dim[0]-self.rightSponge/2.
            epsFact_solid = self.rightSponge/2.
            sign = -1.
            self.zones[flag] = bc.RelaxationZone(domain=self.domain,
                                                 zone_type='absorption',
                                                 sign=sign, center_x=center_x,
                                                 waves=waves,
                                                 windSpeed=windSpeed,
                                                 epsFact_solid=epsFact_solid,
                                                 dragAlphaTypes=dragAlphaTypes,
                                                 dragBetaTypes=dragBetaTypes,
                                                 porosityTypes=porosityTypes)


class RigidBody(AuxiliaryVariables.AV_base):

    def __init__(self, shape, he=1., cfl_target=0.9, dt_init=0.001):
        self.shape = shape
        # if isinstance(shape, (Rectangle, Cuboid)):
        #     shape._setInertiaTensor()
        self.dt_init = dt_init
        self.he = he
        self.cfl_target = 0.9
        self.last_position = np.array([0., 0., 0.])
        self.rotation_matrix = np.eye(3)
        self.h = np.array([0., 0., 0.])
        self.barycenter = np.zeros(3)
        self.i_start = None  # will be retrieved from setValues() of Domain
        self.i_end = None  # will be retrieved from setValues() of Domain

    def step(self, dt):
        nd = self.shape.domain.nd
        # acceleration from force
        self.acceleration = self.F/self.shape.mass
        # angular acceleration from moment
        if sum(self.M) != 0:
            self.inertia = self.shape.getInertia(self.M, self.shape.barycenter)
            assert self.inertia != 0, 'Zero inertia: inertia tensor (It)' \
                                      'was not set correctly!'
            ang_acc = self.M[:]/self.inertia
        else:
            self.inertia = None
            ang_acc = np.array([0., 0., 0.])
        # substeps for smoother motion between timesteps
        ang_disp = 0
        substeps = 20
        dt_sub = dt/float(substeps)
        self.h = np.zeros(3)
        for i in range(substeps):
            # displacement
            self.velocity += self.acceleration*dt_sub
            self.h[:] += self.velocity*dt_sub
            # rotation
            self.angvel += ang_acc*dt_sub
            ang_disp += self.angvel*dt_sub
        # translate
        self.shape.translate(self.h[:nd])
        # rotate
        self.ang = np.linalg.norm(ang_disp)
        if nd == 2 and self.angvel[2] < 0:
            self.ang = -self.ang
        if self.ang != 0.:
            self.shape.rotate(self.ang, self.angvel, self.shape.barycenter)
            self.rotation[:nd, :nd] = self.shape.coords_system
            self.rotation_matrix[:] = np.dot(np.linalg.inv(self.last_rotation),
                                             self.rotation)
        else:
            self.rotation_matrix[:] = np.eye(3)
        self.barycenter[:] = self.shape.barycenter
        self.position[:] = self.shape.barycenter
        if self.shape.record_values is True:
            self.recordValues()

    def recordValues(self):
        t_last = self.model.stepController.t_model_last
        dt_last = self.model.levelModelList[-1].dt_last
        time = t_last-dt_last
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
        if comm.isMaster():
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
        # flagMin = min(m.mesh.elementBoundaryMaterialTypes)
        self.nForces = flagMax+1
        return self

    def calculate_init(self):
        """
        Function called at the very beginning of the simulation by proteus.
        (!) name of the function has to be calculate_init()
        """
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
        self.rotation[:nd, :nd] = shape.coords_system
        self.last_rotation = np.eye(3)
        self.last_rotation[:nd, :nd] = shape.coords_system
        self.F = np.zeros(3, 'd')
        self.M = np.zeros(3, 'd')
        self.last_F = np.zeros(3, 'd')
        self.last_M = np.zeros(3, 'd')
        self.ang = 0.
        self.barycenter = self.shape.barycenter
        self.angvel = np.zeros(3, 'd')
        self.last_angvel = np.zeros(3, 'd')
        if nd == 2:
            self.Fg = self.shape.mass*np.array([0., -9.81, 0.])
        if nd == 3:
            self.Fg = self.shape.mass*np.array([0., 0., -9.81])
        comm = Comm.get()
        if comm.isMaster():
            if self.shape.record_values is True:
                self.record_file = os.path.join(Profiling.logDir,
                                                self.shape.record_filename)
                with open(self.record_file, 'w') as csvfile:
                    writer = csv.writer(csvfile, delimiter=',')
                    writer.writerow(self.shape.record_names)

    def calculate(self):
        """
        Function called at each time step by proteus.
        (!) name of the function has to be calculate()
        """
        # store previous values
        self.last_position[:] = self.position
        self.last_velocity[:] = self.velocity
        self.last_acceleration[:] = self.acceleration
        self.last_rotation[:] = self.rotation
        self.last_angvel[:] = self.angvel
        self.last_F[:] = self.F
        self.last_M[:] = self.M
        # for first time step
        try:
            dt = self.model.levelModelList[-1].dt_last
        except:
            dt = self.dt_init
        # update forces and moments for current body/shape
        i0, i1 = self.i_start, self.i_end
        # get forces
        F_p = self.model.levelModelList[-1].coefficients.netForces_p[i0:i1, :]
        F_v = self.model.levelModelList[-1].coefficients.netForces_v[i0:i1, :]
        F_g = self.Fg
        F = np.sum(F_p + F_v, axis=0) + F_g
        # get moments
        M_t = self.model.levelModelList[-1].coefficients.netMoments[i0:i1, :]
        M = np.sum(M_t, axis=0)
        # store F and M with DOF constraints to body
        self.F[:] = F*self.shape.free_x
        self.M[:] = M*self.shape.free_r
        # calculate new properties
        self.step(dt)
        # log values
        t_previous = self.model.stepController.t_model_last-dt
        t_current = self.model.stepController.t_model_last
        h = self.h
        last_pos, pos = self.last_position, self.position
        last_vel, vel = self.last_velocity, self.velocity
        rot = self.rotation
        rot_x = -atan2(rot[2, 1], rot[1, 2])
        rot_y = -atan2(-rot[0, 2], sqrt(rot[2, 1]**2+rot[2, 2]**2))
        rot_z = np.degrees(-atan2(rot[1, 0], rot[0, 0]))
        log("================================================================")
        log("=================== Rigid Body Calculation =====================")
        log("================================================================")
        log("Name: " + `self.shape.name`)
        log("================================================================")
        log("[proteus]     t=%1.5fsec to t=%1.5fsec" % \
            (t_previous, t_current))
        log("[proteus]    dt=%1.5fsec" % (dt))
        log("[body] ============== Pre-calculation attributes  ==============")
        log("[proteus]     t=%1.5fsec" % (t_previous))
        log("[proteus]     F=(% 12.7e, % 12.7e, % 12.7e)" % (F[0], F[1], F[2]))
        log("[proteus] F*DOF=(% 12.7e, % 12.7e, % 12.7e)" % (F[0], F[1], F[2]))
        log("[proteus]     M=(% 12.7e, % 12.7e, % 12.7e)" % (M[0], M[1], M[2]))
        log("[proteus] M*DOF=(% 12.7e, % 12.7e, % 12.7e)" % (M[0], M[1], M[2]))
        log("[body]      pos=(% 12.7e, % 12.7e, % 12.7e)" % \
            (last_pos[0], last_pos[1], last_pos[2]))
        log("[body]      vel=(% 12.7e, % 12.7e, % 12.7e)" % \
            (last_vel[0], last_vel[1], last_vel[2]))
        log("[body] ===============Post-calculation attributes ==============")
        log("[body]        t=%1.5fsec" % (t_current))
        log("[body]        h=(% 12.7e, % 12.7e, % 12.7e)" % (h[0], h[1], h[2]))
        log("[body]      pos=(% 12.7e, % 12.7e, % 12.7e)" % \
            (pos[0], pos[1], pos[2]))
        log("[body]      vel=(% 12.7e, % 12.7e, % 12.7e)" % \
            (vel[0], vel[1], vel[2]))
        log("[body]      rot=(% 12.7e, % 12.7e, % 12.7e)" % \
            (rot_x, rot_y, rot_z))
        log("================================================================")


def relative_vec(vec1, vec0):
    """
    function giving coordinates of a vector relative to another vector
    (projecting vec0 as the z-axis for vec1)

    :param vec1: vector to get new coordinates
    :param vec0: vector of reference
    :return: new coordinates of vec1
    """
    # spherical coords vec0
    x0, y0, z0 = vec0
    r0 = sqrt(x0**2+y0**2+z0**2)  # radius from origin
    t0 = atan2(y0, x0)  # angle on x-y plane
    p0 = acos(z0/r0)  # angle from z-axis
    # spherical coords vec1
    x1, y1, z1 = vec1
    r1 = sqrt(x1**2+y1**2+z1**2)
    t1 = atan2(y1, x1)
    p1 = acos(z1/r1)
    # get new coords for vec1:
    t1_new = t0-t1
    p1_new = p0-p1
    x1_new = r1*sin(p1_new)*cos(t1_new)
    y1_new = r1*sin(p1_new)*sin(t1_new)
    z1_new = r1*cos(p1_new)
    return (x1_new, y1_new, z1_new)


def assembleDomain(domain):
    """
    This function sets up everything needed for the domain, meshing, and
    AuxiliaryVariables calculations (if any).
    It should always be called after defining and manipulating all the shapes
    to be attached to the domain.

    :param domain: domain to assemble
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
    zones_global = {}
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
            domain.facets += list(shape.facets+start_vertex)
            domain.facetFlags += list(shape.facetFlags+start_flag)
        if shape.regions is not None:
            domain.regions += list(shape.regions)
            domain.regionFlags += list(shape.regionFlags+start_rflag)
        if shape.holes is not None:
            domain.holes += list(shape.holes)
        domain.getBoundingBox()
        # --------------------------- #
        # --- AUXILIARY VARIABLES --- #
        # --------------------------- #
        aux = domain.auxiliaryVariables
        # ----------------------------
        # RIGID BODIES
        if 'RigidBody' in shape.auxiliaryVariables.keys():
            aux += [RigidBody(shape)]
            # fixing mesh on rigid body
            for boundcond in shape.BC_list:
                boundcond.setMoveMesh(aux[-1])
            # update the indice for force/moment calculations
            aux[-1].i_start = start_flag+1
            aux[-1].i_end = start_flag+1+len(shape.BC_list)
        # ----------------------------
        # ABSORPTION/GENERATION ZONES
        if 'RelaxZones' in shape.auxiliaryVariables.keys():
            if not zones_global:
                aux += [bc.RelaxationZoneWaveGenerator(zones_global,
                                                       domain.nd)]
            # create arrays of default values
            domain.porosityTypes = np.ones(len(domain.regionFlags)+1)
            domain.dragAlphaTypes = np.zeros(len(domain.regionFlags)+1)
            domain.dragBetaTypes = np.zeros(len(domain.regionFlags)+1)
            domain.epsFact_solid = np.zeros(len(domain.regionFlags)+1)
            i0 = start_region+1
            for flag, zone in shape.zones.iteritems():
                ind = [i for i, f in enumerate(shape.regionFlags) if f == flag]
                for i1 in ind:
                    domain.porosityTypes[i0+i1] = zone.porosityTypes
                    domain.dragAlphaTypes[i0+i1] = zone.dragAlphaTypes
                    domain.dragBetaTypes[i0+i1] = zone.dragBetaTypes
                    domain.epsFact_solid[i0+i1] = zone.epsFact_solid
                # update dict with global key instead of local key
                key = flag+start_rflag
                zones_global[key] = zone
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

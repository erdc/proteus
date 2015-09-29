"""
Module creating predifined shapes. Each shape needs a proteus.Domain as argument.
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
from math import cos, sin, sqrt
from proteus import AuxiliaryVariables, Archiver

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
        self.free_x = (1, 1, 1)
        self.free_r = (1, 1, 1)
        self.snv = len(self.domain.vertices)  # total number of vertices in domain when shape.__init__
        self.snf = len(self.domain.facets)
        self.sns = len(self.domain.segments)
        self.snr = len(self.domain.regions)
        self.snbc = len(self.domain.bc)

    def _addShape(self):
        """
        Adds Shape information to the domain
        """
        # add new information to the domain
        # get maximum flag defined in domain so far
        # need to add +1 for flags as 0 cannot be used
        if self.domain.vertexFlags or self.domain.segmentFlags or self.domain.facetFlags:
            flag = max([max(self.domain.vertexFlags), max(self.domain.segmentFlags), max(self.domain.facetFlags)])
            flag += 1
        else:
            flag = 1
        self.domain.vertices += self.vertices.tolist()
        self.domain.vertexFlags += (self.vertexFlags+flag).tolist()
        self.domain.segments += (self.segments+self.snv).tolist()
        self.domain.segmentFlags += (self.segmentFlags+flag).tolist()
        if self.domain.nd == 3:
            self.domain.facets += (self.facets+self.snv).tolist()
            self.domain.facetFlags += (self.facetFlags+flag).tolist()
            if len(self.holes) != 0:
                self.domain.holes += self.holes.tolist()
        self.domain.regions += self.regions.tolist()
        self.domain.regionFlags += (self.regionFlags+flag).tolist()
        if len(self.domain.bc) == 0: # need to add None boundary condition at 0 indice
            self.domain.bc += [bc.BoundaryConditions()]
        self.domain.bc += self.bc
        if self.barycentre:
            if self.domain.barycenters is not None:
                self.domain.barycenters = np.append(self.domain.barycenters, self.barycenters, axis=0)  # need to change to array
            else:
                self.domain.barycenters = self.barycenters
        self.domain.update()

    def _updateDomain(self):
        """
        Updates domain when vertex and region coords have been changed
        """
        nv = len(self.vertices)
        nr = len(self.regions)
        nf = len(self.facets)
        self.domain.vertices[self.snv:self.snv+nv] = self.vertices.tolist()
        self.domain.regions[self.snr:self.snr+nr] = self.regions.tolist()
        self.domain.barycenters[self.snf:self.snf+nf] = self.barycentre
        self.domain.update()

    def setPosition(self, coords):
        """
        Set position of the Shape
        :arg coords: new set of coordinates for the Shape
        """
        old_coords = np.array(self.coords)
        trans = coords - old_coords
        self.translate(trans)

    def setBarycentre(self, barycentre):
        """
        Set barycentre of the shape, from the reference frame of the shape
        :arg barycentre: coordinates of barycentre from centre of shape
        """
        self.barycentre = np.einsum('ij,i->j', self.coords_system, barycentre)
        self.domain.barycenters[:] = self.barycentre

    def setConstraints(self, free_x, free_r):
        """
        Sets constraints on the Shape
        :arg free_x: translational constraints
        :arg free_r: rotational constraints
        """
        self.free_x = free_x
        self.free_r = free_r

    def setRegions(self, regions):
        self.regions = regions
        self._updateDomain
    
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
            self.vertices[:] = rotation2D(points=self.vertices, rot=rot, pivot=pivot)
            self.regions[:] = rotation2D(points=self.regions, rot=rot, pivot=pivot)
            self.coords_system[:] = rotation2D(points=self.coords_system, rot=rot, pivot=(0.,0.,0.))
            self.segments_orientation[:] = rotation2D(points=self.facets_orientation, rot=rot, pivot=(0., 0., 0.))
            if rotate_barycentre:
                self.barycentre[:] = rotation2D(points=self.barycentre, rot=rot, pivot=pivot)
        elif self.domain.nd == 3:
            self.vertices[:] = rotation3D(points=self.vertices, rot=rot, axis=axis, pivot=pivot)
            if len(self.holes) != 0:
                self.holes[:] = rotation3D(points=self.holes, rot=rot, axis=axis, pivot=pivot)
            self.regions[:] = rotation3D(points=self.regions, rot=rot, axis=axis, pivot=pivot)
            self.coords_system[:] = rotation3D(points=self.coords_system, rot=rot, axis=axis, pivot=(0.,0.,0.))
            self.facets_orientation[:] = rotation3D(points=self.facets_orientation, rot=rot, axis=axis, pivot=(0., 0., 0.))
            if rotate_barycentre:
                self.barycentre[:] = rotation3D(points=self.barycentre, rot=rot, axis=axis, pivot=pivot)
        self._updateDomain()

    def translate(self, trans):
        """
        Function to translate Shape
        :arg trans: translation values
        """
        self.vertices += trans
        self.regions += trans
        self.coords += trans
        self.barycentre += trans
        if self.domain.nd == 3 and len(self.holes) != 0:
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
        return self.coords

    def getRotation(self):
        return self.coords_system

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

    def setRigidBody(self):
        self.RigidBody = RigidBody(shape=self)

    def setTank(self):
        for bc in self.bc:
            bc.setTank()
        self.RigidBody = RigidBody(shape=self)


class Cuboid(Shape):
    """
    Class to create a cuboid
    :arg domain: domain of the cuboid
    :arg dim: dimensions of the cuboid (list or array)
    :arg coords: coordinates of the cuboid (list or array)
    """
    def __init__(self, domain, dim=(0.,0.,0.), coords=(0.,0.,0.), barycentre=None, tank=False):
        Shape.__init__(self, domain)
        self.dim = L, W, H = dim  # length, width height
        self.volume = L*W*H
        self.coords = x, y, z = coords
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
        self.facets = np.array([[[0, 1, 2, 3]],  # bottom
                                [[0, 1, 5, 4]],  # front
                                [[1, 2, 6, 5]],  # right
                                [[2, 3, 7, 6]],  # back
                                [[3, 0, 4, 7]],  # left
                                [[4, 5, 6, 7]]])  # top
        self.facets_orientation = np.array([[0,  0, -1],
                                            [-1, 0,  0],
                                            [0,  1,  0],
                                            [1,  0,  0],
                                            [0, -1,  0],
                                            [0,  0,  1]])
        self.regions = np.array([[x, y, z]])
        if not tank:
            self.holes = np.array([coords])
        # defining flags for boundary conditions
        self.facetFlags = np.array([0, 1, 2, 3, 4, 5])  # bottom, front, right, back, left, top
        self.vertexFlags = np.array([0, 0, 0, 0, 5, 5, 5, 5])  # only top and bottom for vertices
        self.segmentFlags = np.array([0, 0, 0, 0, 5, 5, 5, 5, 1, 1, 3, 3])
        self.regionFlags = np.array([0])
        # Initialize (empty) boundary conditions
        b_or = self.facets_orientation
        self.bc_bottom = bc.BoundaryConditions(b_or=b_or, b_i=0)
        self.bc_front = bc.BoundaryConditions(b_or=b_or, b_i=1)
        self.bc_right = bc.BoundaryConditions(b_or=b_or, b_i=2)
        self.bc_back = bc.BoundaryConditions(b_or=b_or, b_i=3)
        self.bc_left = bc.BoundaryConditions(b_or=b_or, b_i=4)
        self.bc_top = bc.BoundaryConditions(b_or=b_or, b_i=5)
        self.bc = [self.bc_bottom, self.bc_front, self.bc_right,
                    self.bc_back, self.bc_left, self.bc_top]
        self.barycentre = barycentre or self.coords
        self.barycenters = np.array([self.barycentre for facet in self.facets])
        self.It = np.array([[(W**2.+H**2.)/12., 0, 0],
                            [0, (L**2.+H**2.)/12., 0],
                            [0, 0, (W**2.+L**2.)/12.]])
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
    def __init__(self, domain, dim=(0.,0.), coords=(0.,0.), barycentre=None):
        Shape.__init__(self, domain)
        self.dim = L, H = dim  # length, height
        self.coords = x, y = coords
        self.coords_system = np.eye(2)
        self.barycentre = barycentre or coords
        self.vertices = np.array([[x-0.5*L, y-0.5*H],
                                  [x+0.5*L, y-0.5*H],
                                  [x+0.5*L, y+0.5*H],
                                  [x-0.5*L, y+0.5*H]])
        self.segments= np.array([[0,1],[1,2],[2,3],[3,0]])  # bottom, right, top, left
        self.segments_orientation = [[0, -1],
                                     [1, 0],
                                     [0, 1],
                                     [0, -1]]
        self.barycenters = np.array([self.barycentre for segment in self.sgments])
        self.regions = np.array([[x, y]])
        self.segmentFlags = np.array([0, 1, 2, 3]) # bottom, right, top, left
        self.vertexFlags = np.array([0, 0, 2, 2])  # bottom, bottom, top, top
        self.regionFlags = np.array([0])
        # Initialize (empty) boundary conditions
        self.bc_bottom = bc.BoundaryConditions(b_or=self.segments_orientation, b_i=0)
        self.bc_right = bc.BoundaryConditions(b_or=self.segments_orientation, b_i=1)
        self.bc_top = bc.BoundaryConditions(b_or=self.segments_orientation, b_i=2)
        self.bc_left = bc.BoundaryConditions(b_or=self.segments_orientation, b_i=3)
        self.bc = [self.bc_bottom, self.bc_right, self.bc_top, self.bc_left]
        self.It = np.array([[(L**2.)/3., 0],
                            [0, (H**2.)/3.]])
        self._addShapes()  # adding shape to domain

    def _setInertiaTensor(self):
        """
        Set the (new) inertia tensor of the shape
        """
        L, H = self.dim
        self.It[:] = [[(L**2.)/3., 0],
                      [0, (H**2.)/3.]]

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






class RigidBody(AuxiliaryVariables.AV_base):

    def __init__(self, shape, he=1., cfl_target=0.9, dt_init=0.001):
        self.shape = shape
        self.dt_init = dt_init
        self.he = he
        self.cfl_target = 0.9
        self.position = self.shape.barycentre
        nd = self.shape.domain.nd
        self.velocity = np.zeros(nd)
        self.last_velocity = np.zeros(nd)
        self.h = np.zeros(nd)
        self.rotation = np.eye(nd)
        self.last_rotation = np.eye(nd)
        self.F = np.zeros(nd, 'd')
        self.M = np.zeros(nd, 'd')
        self.last_F = np.zeros(nd, 'd')
        self.last_M = np.zeros(nd, 'd')
        self.barycenters = np.zeros((len(self.shape.domain.facets), nd), 'd')
        startnf = self.shape.snf
        endnf = startnf + len(self.shape.facets) + 1
        self.barycenters[startnf:endnf,:] = self.shape.barycentre

    def step(self, dt):
        # displacement from force
        self.acceleration = self.F/self.shape.mass
        self.velocity += self.acceleration*dt
        displacement = self.velocity*dt
        self.shape.translate(displacement)        # rotation due to moment
        I = self.shape.getInertia(vec=self.M, pivot=self.shape.barycentre)
        acc_ang = self.M/I
        velocity_ang += acc_ang*dt
        ang_disp = velocity_ang*dt
        self.ang = np.linalg.norm(ang_disp)
        self.shape.rotate(rot=self.ang, axis=ang_vel, pivot=self.shape.barycentre)

    def attachModel(self, model, ar):
        self.model = model
        self.ar = ar
        self.writer = Archiver.XdmfWriter()
        self.nd = model.levelModelList[-1].nSpace_global
        m = self.model.levelModelList[-1]
        flagMax = max(m.mesh.elementBoundaryMaterialTypes)
        flagMin = min(m.mesh.elementBoundaryMaterialTypes)
        # assert(flagMin >= 0)
        # assert(flagMax <= 7)
        # self.nForces=flagMax+1
        # assert(self.nForces <= 8)
        return self

    def get_u(self):
        return self.last_velocity[0]

    def get_v(self):
        return self.last_velocity[1]

    def get_w(self):
        return self.last_velocity[2]

    def getVelocity(self):
        return self.last_velocity

    def calculate_init(self):
        self.last_F = None
        self.calculate()

    def calculate(self):
        self.last_velocity = copy.deepcopy(self.velocity)
        self.last_position = copy.deepcopy(self.position)
        self.last_rotation = self.rotation.copy()
        # self.last_rotation_inv = np.linalg.inv(self.last_rotation)
        import copy
        try:
            dt = self.model.levelModelList[-1].dt_last
        except:
            dt = self.dt_init
        startnf = self.shape.snf
        endnf = startnf + len(self.shape.facets)
        F = np.sum(self.model.levelModelList[-1].coefficients.netForces_p[startnf:endnf,:], axis=0) + self.model.levelModelList[-1].coefficients.netForces_v[startnf:endnf,:];
        M = np.sum(self.model.levelModelList[-1].coefficients.netMoments[startnf:endnf,:], axis=0)
        logEvent("x Force " +`self.model.stepController.t_model_last`+" "+`F[0]`)
        logEvent("y Force " +`self.model.stepController.t_model_last`+" "+`F[1]`)
        logEvent("z Force " +`self.model.stepController.t_model_last`+" "+`F[2]`)
        logEvent("x Moment " +`self.model.stepController.t_model_last`+" "+`M[0]`)
        logEvent("y Moment " +`self.model.stepController.t_model_last`+" "+`M[1]`)
        logEvent("z Moment " +`self.model.stepController.t_model_last`+" "+`M[2]`)
        logEvent("dt " +`dt`)
        nSteps = 10
        Fstar = F
        Mstar = M
        for i in range(nSteps):
            self.F = Fstar*self.shape.free_x
            self.M = Mstar*self.shape.free_r
            self.step(dt/float(nSteps))
        # store forces
        self.last_F[:] = F
        self.last_M[:] = M
        x,y,z = self.shape.barycentre
        u,v,w = self.getVelocity()
        # update barycenters
        self.barycenters[startnf:endnf, :] = [x, y, z]
        self.position = (x,y,z)
        self.velocity = (u,v,w)
        self.rotation = self.shape.getRotation()
        self.h = np.array(self.position)-self.last_position
        # set new boundary conditions for moveMesh
        for bc in self.shape.bc:
            bc.setMoveMesh(self.h, self.last_rotation, self.rotation, self.last_position)
        logEvent("%1.2fsec: pos=(%21.16e, %21.16e, %21.16e) vel=(%21.16e, %21.16e, %21.16e) h=(%21.16e, %21.16e, %21.16e)" % (self.model.stepController.t_model_last,
                                                                                                                            self.position[0],
                                                                                                                            self.position[1],
                                                                                                                            self.position[2],
                                                                                                                            self.velocity[0],
                                                                                                                            self.velocity[1],
                                                                                                                            self.velocity[2],
                                                                                                                            self.h[0],
                                                                                                                            self.h[1],
                                                                                                                            self.h[2]))







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

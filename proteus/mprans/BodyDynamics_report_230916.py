from math import cos, sin, sqrt, atan2, acos, asin
import csv
import os
import numpy as np
from proteus import AuxiliaryVariables, Archiver, Comm, Profiling
from proteus.Profiling import logEvent as logEvent



class RigidBody(AuxiliaryVariables.AV_base, object):
    """
    Auxiliary variable used to calculate attributes of an associated shape
    class instance acting as a rigid body. To set a shape as a rigid body, use
    shape.setRigidBody(). The class instance is created automatically when
    shape.setRigidBody() has been called and after calling assembleDomain().

    Parameters
    ----------
    shape: proteus.mprans.SpatialTools.Shape_RANS
        Class instance of the shape associated to the rigid body calculations.
    cfl_target: Optional[float]
        UNUSED (to implement), sets the maximum displacement of the body
        allowed per time step.
    dt_init: float
        first time step of the simulation.
    """

    def __init__(self, shape, cfl_target=0.9, dt_init=0.001):
        self.Shape = shape
        # if isinstance(shape, (Rectangle, Cuboid)):
        #     shape._setInertiaTensor()
        self.dt_init = dt_init
        self.cfl_target = 0.9
        self.last_position = np.array([0., 0., 0.])
        self.rotation_matrix = np.eye(3)
        self.h = np.array([0., 0., 0.])
        self.barycenter = np.zeros(3)
        self.i_start = None  # will be retrieved from setValues() of Domain
        self.i_end = None  # will be retrieved from setValues() of Domain
        self.It = self.Shape.It
        if 'RigidBody' not in shape.auxiliaryVariables:
            shape.auxiliaryVariables['RigidBody'] = self

    def attachModel(self, model, ar):
        """
        Attaches model to auxiliary variable
        """
        self.model = model
        self.ar = ar
        self.writer = Archiver.XdmfWriter()
        self.nd = model.levelModelList[-1].nSpace_global
        m = self.model.levelModelList[-1]
        flagMax = max(m.mesh.elementBoundaryMaterialTypes)
        # flagMin = min(m.mesh.elementBoundaryMaterialTypes)
        self.nForces = flagMax+1
        return self


    def getPressureForces(self):
        i0, i1 = self.i_start, self.i_end
        return self.model.levelModelList[-1].coefficients.netForces_p[i0:i1, :]

    def getShearForces(self):
        i0, i1 = self.i_start, self.i_end
        return self.model.levelModelList[-1].coefficients.netForces_v[i0:i1, :]

    def getGravityForce(self):
        return self.Fg

    def getMoments(self):
        i0, i1 = self.i_start, self.i_end
        return self.model.levelModelList[-1].coefficients.netMoments[i0:i1, :]

    def getTotalMoments(self):
        M = self.getMoments()
        M_t = np.sum(M, axis=0)
        return M_t

    def getTotalForce(self):
        F_p = self.getPressureForces()
        F_v = self.getShearForces()
        F_g = self.Fg
        F_t = np.sum(F_p + F_v, axis=0) + F_g
        return F_t

    def getAcceleration(self):
        a = self.F/self.mass
        return a

    def getAngularAcceleration(self):
        if sum(self.M) != 0:
            self.inertia = self.getInertia(self.M, self.Shape.barycenter)
            assert self.inertia != 0, 'Zero inertia: inertia tensor (It)' \
                                      'was not set correctly!'
            self.ang_acc = self.M[:]/self.inertia
        else:
            self.inertia = None
            self.ang_acc = np.array([0., 0., 0.])

    def getDisplacement(self,dt):
        # acceleration from force
        self.acceleration = self.getAcceleration()
        # substeps for smoother motion between timesteps
        dt_sub = dt/float(substeps)
        for i in range(substeps):
            # displacement
            self.h[:], self.velocity[:] = forward_euler(p0=self.h, v0=self.velocity,
                                                             a=self.acceleration, dt=dt_sub)

    def getAngularDisplacement(self):
        # angular acceleration from moment
        self.ang_acc = self.getAngularAcceleration()
        dt_sub = dt/float(substeps)
        self.ang_disp = np.array([0.,0.,0.])
        for i in range(substeps):
            # rotation
            self.ang_disp[:], self.ang_vel[:] = forward_euler(p0=self.ang_disp, v0=self.ang_vel,
                                                         a=self.ang_acc, dt=dt_sub)

    def step(self, dt, substeps=20):
        """
        Step for rigid body calculations in Python

        Parameters
        ----------
        dt: float
            time step
        """
        nd = self.Shape.Domain.nd
        # reinitialise displacement values
        ang_disp = 0
        self.ang_disp = np.zeros(3)
        self.h[:] = np.zeros(3)

        # Translational motion calculation
        self.getDisplacement(dt)
        # Rotational motion calculation
        self.getAngularDisplacement(dt)

        # translate
        self.Shape.translate(self.h[:nd])
        # rotate
        self.ang = np.linalg.norm(self.ang_disp)
        if nd == 2 and self.ang_vel[2] < 0:
            self.ang = -self.ang
        if self.ang != 0.:
            if self.overturning != True:
                self.Shape.rotate(self.ang, self.ang_vel, self.Shape.barycenter)
            else:
                self.Shape.rotate(self.ang, self.ang_vel, self.pivot_friction)
            self.rotation[:nd, :nd] = self.Shape.coords_system
            self.rotation_matrix[:] = np.dot(np.linalg.inv(self.last_rotation),
                                             self.rotation)
        else:
            self.rotation_matrix[:] = np.eye(3)
        self.barycenter[:] = self.Shape.barycenter
        self.position[:] = self.Shape.barycenter
        if self.record_values is True:
            self.recordValues()

    def recordValues(self):
        """
        Records values of rigid body attributes at each time step in a csv file.
        """
        comm = Comm.get()
        if comm.isMaster():
            t_last = self.model.stepController.t_model_last
            dt_last = self.model.levelModelList[-1].dt_last
            values_towrite = []
            t = t_last-dt_last
            if t == 0:
                headers = []
                if self.record_dict['time'] is True:
                    headers += ['t']
                if self.record_dict['pos'] is True:
                    headers += ['x', 'y', 'z']
                if self.record_dict['rot'] is True:
                    headers += ['rx', 'ry', 'rz']
                if self.record_dict['ang_disp'] is True:
                    headers += ['ang_dispx', 'ang_dispy', 'ang_dispz']
                if self.record_dict['F'] is True:
                    headers += ['Fx', 'Fy', 'Fz']
                if self.record_dict['M'] is True:
                    headers += ['Mx', 'My', 'Mz']
                if self.record_dict['inertia'] is True:
                    headers += ['inertia']
                if self.record_dict['vel'] is True:
                    headers += ['vel_x', 'vel_y', 'vel_z']
                if self.record_dict['acc'] is True:
                    headers += ['acc_x', 'acc_y', 'acc_z']
                if self.record_dict['ang_vel'] is True:
                    headers += ['ang_velx', 'ang_vely', 'ang_velz']
                if self.record_dict['ang_acc'] is True:
                    headers += ['ang_accx', 'ang_accy', 'ang_accz']
                with open(self.record_file, 'w') as csvfile:
                    writer = csv.writer(csvfile, delimiter=',')
                    writer.writerow(headers)
            if self.record_dict['time'] is True:
                t = t_last-dt_last
                values_towrite += [t]
            if self.record_dict['pos'] is True:
                x, y, z = self.last_position
                values_towrite += [x, y, z]
            if self.record_dict['rot'] is True:
                rot = self.last_rotation
                rx = atan2(rot[1, 2], rot[2, 2])
                ry = -asin(rot[0, 2])
                rz = atan2(rot[0, 1], rot[0, 0])
                values_towrite += [rx, ry, rz]
            if self.record_dict['ang_disp'] is True:
                ang_dispx, ang_dispy, ang_dispz = self.ang_disp
                values_towrite += [ang_dispx, ang_dispy, ang_dispz]
            if self.record_dict['F'] is True:
                Fx, Fy, Fz = self.F
                values_towrite += [Fx, Fy, Fz]
            if self.record_dict['M'] is True:
                Mx, My, Mz = self.M
                values_towrite += [Mx, My, Mz]
            if self.record_dict['inertia'] is True:
                values_towrite += [self.inertia]
            if self.record_dict['vel'] is True:
                vel_x, vel_y, vel_z = self.velocity
                values_towrite += [vel_x, vel_y, vel_z]
            if self.record_dict['acc'] is True:
                acc_x, acc_y, acc_z = self.acceleration
                values_towrite += [acc_x, acc_y, acc_z]
            if self.record_dict['ang_vel'] is True:
                ang_velx, ang_vely, ang_velz = self.ang_vel
                values_towrite += [ang_velx, ang_vely, ang_velz]
            if self.record_dict['ang_acc'] is True:
                ang_accx, ang_accy, ang_accz = self.ang_acc
                values_towrite += [ang_accx, ang_accy, ang_accz]
            with open(self.record_file, 'a') as csvfile:
                writer = csv.writer(csvfile, delimiter=',')
                writer.writerow(values_towrite)


    def setConstraints(self, free_x, free_r):
        """
        Sets constraints on the Shape (for moving bodies)

        Parameters
        ----------
        free_x: array_like
            Translational constraints.
        free_r: array_like
            Rotational constraints.
        """
        self.free_x = np.array(free_x)
        self.free_r = np.array(free_r)

    def setMass(self, mass):
        """
        Set mass of the shape.

        Parameters
        ----------
        mass: float
            mass of the body
        """
        self.mass = float(mass)

    def setDensity(self, density):
        """
        Set density.

        Parameters
        ----------
        density: float
            Density of the shape
        """
        self.density = float(density)

    def setInertiaTensor(self, It):
        """
        Set the inertia tensor of the shape

        Parameters
        ----------
        It: array_like, float
            Inertia tensor of the body (3x3 array in 3D, float in 2D)

        Notes
        -----
        The inertia tensor should not be already scaled with the mass of the
        shape.
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

        Parameters
        ----------
        vec: array_like
            Vector around which the body rotates.
        pivot: Optional[array_like]
            Pivotal point around which the body rotates. If not set, it will
            be the barycenter coordinates

        Returns
        -------
        I: float
            inertia of the mass

        Notes
        -----
        The inertia is calculated relative to the coordinate system of the
        shape (self.coords_system). If the shape was not initialised with a
        position corresponding to its inertia tensor (e.g. shape was already
        rotated when initialised), set the coordinate system accordingly
        before calling this function
        """
        assert self.It is not None, 'No inertia tensor! (' + self.name + ')'
        if pivot is None:
            pivot = self.barycenter
        # Pivot coords relative to shape centre of mass
        pivot = pivot-np.array(self.barycenter)
        # making unity vector/axis of rotation
        vec = vx, vy, vz = np.array(vec)
        length_vec = np.sqrt(vx**2+vy**2+vz**2)
        vec = vec/length_vec
        if self.Shape.Domain.nd == 2:
            I = self.It*self.mass
        elif self.Shape.Domain.nd == 3:
            # vector relative to original position of shape:
            vec = np.dot(vec, np.linalg.inv(self.coords_system))
            cx, cy, cz = vec
            # getting the tensor for calculaing moment of inertia
            # from arbitrary axis
            vt = np.array([[cx**2, cx*cy, cx*cz],
                           [cx*cy, cy**2, cy*cz],
                           [cx*cz, cy*cz, cz**2]])
            # total moment of inertia
            I = np.einsum('ij,ij->', self.mass*self.It, vt)
        return I

    def setRecordValues(self, filename=None, all_values=False, time=True,
                        pos=False, rot=False, ang_disp=False, F=False, M=False, inertia=False,
                        vel=False, acc=False, ang_vel=False, ang_acc=False):
        """
        Sets the rigid body attributes that are to be recorded in a csv file
        during the simulation.

        Parameters
        ----------
        filename: Optional[string]
            Name of file, if not set, the file will be named as follows:
            'record_[shape.name].csv'
        all_values: bool
            Set to True to record all values listed below.
        time: bool
            Time of recorded row (default: True).
        pos: bool
            Position of body (default: False. Set to True to record).
        rot: bool
            Rotation of body (default: False. Set to True to record).
        ang_disp: array
            Angular displecement calculated during rigid body calculation step.
            Applied on the body in order to make it rotating.
        F: bool
            Forces applied on body (default: False. Set to True to record).
        M: bool
            Moments applied on body (default: False. Set to True to record).
        inertia: bool
            Inertia of body (default: False. Set to True to record).
        vel: bool
            Velocity of body (default: False. Set to True to record).
        acc: bool
            Acceleration of body (default: False. Set to True to record).
        ang_vel: array
            Angular velocity of body (default: False. Set to True to record).
        ang_acc: bool
            Angular acceleration of body (default: False. Set to True to record).
        """
        self.record_values = True
        if pos is True:
            x = y = z = True
        if rot is True:
            rot_x = rot_y = rot_z = True
        if ang_disp is True:
            ang_dispx =  ang_dispy = ang_dispz = True
        if F is True:
            Fx = Fy = Fz = True
        if M is True:
            Mx = My = Mz = True
        if vel is True:
            vel_x = vel_y = vel_z = True
        if acc is True:
            acc_x = acc_y = acc_z = True
        if ang_vel is True:
            ang_velx = ang_vely = ang_velz = True
        if ang_acc is True:
            ang_accx = ang_accy = ang_accz = True

        self.record_dict = {'time':time, 'pos': pos, 'rot':rot, 'ang_disp':ang_disp, 'F':F, 'M':M,
                            'inertia': inertia, 'vel': vel, 'acc': acc, 'ang_vel': ang_vel, 'ang_acc': ang_acc}

        if all_values is True:
            for key in self.record_dict:
                self.record_dict[key] = True
        if filename is None:
            self.record_filename = 'record_' + self.name + '.csv'
        else:
            self.record_filename = filename + '.csv'



##################################################
# Calculation step
##################################################

    def calculate_init(self):
        """
        Function called at the very beginning of the simulation by proteus.
        """
        nd = self.Shape.Domain.nd
        shape = self.Shape

        self.position = np.zeros(3)
        self.position[:] = self.Shape.barycenter.copy()
        self.last_position[:] = self.position
        self.velocity = np.zeros(3, 'd')
        self.last_velocity = np.zeros(3, 'd')
        self.acceleration = np.zeros(3, 'd')
        self.last_acceleration = np.zeros(3, 'd')

        self.ang_disp = np.zeros(3, 'd')
        self.rotation = np.eye(3)
        self.rotation[:nd, :nd] = shape.coords_system
        self.last_rotation = np.eye(3)
        self.last_rotation[:nd, :nd] = shape.coords_system
        self.ang_vel = np.zeros(3, 'd')
        self.last_ang_vel = np.zeros(3, 'd')
        self.ang_acc = np.zeros(3, 'd')
        self.last_ang_acc = np.zeros(3, 'd')
        self.F = np.zeros(3, 'd')
        self.M = np.zeros(3, 'd')
        self.last_F = np.zeros(3, 'd')
        self.last_M = np.zeros(3, 'd')
        self.ang = 0.
        self.barycenter = self.Shape.barycenter

        if nd == 2:
            self.Fg = self.mass*np.array([0., -9.81, 0.])
        if nd == 3:
            self.Fg = self.mass*np.array([0., 0., -9.81])
        if self.record_values is True:
            self.record_file = os.path.join(Profiling.logDir,
                                            self.record_filename)


    def calculate(self):
        """
        Function called at each time step by proteus.
        """
        # store previous values
        self.last_position[:] = self.position
        self.last_velocity[:] = self.velocity
        self.last_acceleration[:] = self.acceleration
        self.last_rotation[:] = self.rotation
        self.last_ang_acc[:] = self.ang_acc
        self.last_ang_vel[:] = self.ang_vel
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
        F = self.getTotalForce()
        M = self.getTotalMoments()
        self.F[:] = F2 = F*self.free_x
        # get moments
        self.M[:] = M2 = M*self.free_r
        # store F and M with DOF constraints to body
        # calculate new properties
        self.step(dt)
        # log values
        t_previous = self.model.stepController.t_model_last-dt
        t_current = self.model.stepController.t_model_last
        h = self.h
        last_pos, pos = self.last_position, self.position
        last_vel, vel = self.last_velocity, self.velocity
        rot = self.rotation
        rot_x = atan2(rot[1, 2], rot[2, 2])
        rot_y = -asin(rot[0, 2])
        rot_z = atan2(rot[0, 1], rot[0, 0])
        logEvent("================================================================")
        logEvent("=================== Rigid Body Calculation =====================")
        logEvent("================================================================")
        logEvent("Name: " + `self.Shape.name`)
        logEvent("================================================================")
        logEvent("[proteus]     t=%1.5fsec to t=%1.5fsec" % \
            (t_previous, t_current))
        logEvent("[proteus]    dt=%1.5fsec" % (dt))
        logEvent("[body] ============== Pre-calculation attributes  ==============")
        logEvent("[proteus]     t=%1.5fsec" % (t_previous))
        logEvent("[proteus]     F=(% 12.7e, % 12.7e, % 12.7e)" % (F[0], F[1], F[2]))
        logEvent("[proteus] F*DOF=(% 12.7e, % 12.7e, % 12.7e)" % (F2[0], F2[1], F2[2]))
        logEvent("[proteus]     M=(% 12.7e, % 12.7e, % 12.7e)" % (M[0], M[1], M[2]))
        logEvent("[proteus] M*DOF=(% 12.7e, % 12.7e, % 12.7e)" % (M2[0], M2[1], M2[2]))
        logEvent("[body]      pos=(% 12.7e, % 12.7e, % 12.7e)" % \
            (last_pos[0], last_pos[1], last_pos[2]))
        logEvent("[body]      vel=(% 12.7e, % 12.7e, % 12.7e)" % \
            (last_vel[0], last_vel[1], last_vel[2]))
        logEvent("[body] ===============Post-calculation attributes ==============")
        logEvent("[body]        t=%1.5fsec" % (t_current))
        logEvent("[body]        h=(% 12.7e, % 12.7e, % 12.7e)" % (h[0], h[1], h[2]))
        logEvent("[body]      pos=(% 12.7e, % 12.7e, % 12.7e)" % \
            (pos[0], pos[1], pos[2]))
        logEvent("[body]      vel=(% 12.7e, % 12.7e, % 12.7e)" % \
            (vel[0], vel[1], vel[2]))
        logEvent("[body]      rot=(% 12.7e, % 12.7e, % 12.7e)" % \
            (rot_x, rot_y, rot_z))
        logEvent("================================================================")


##################################################


def forward_euler(p0, v0, a, dt):
    v1 = v0+a*dt
    p1 = p0+v1*dt
    return p1, v1

def leapfrog(p0, v0, a, dt):
    v1 = v0+a*dt
    p1 = p0+v0*dt+0.5*a*dt**2
    return p1, v1


########################################################################################################################################################################################################
########################################################################################################################################################################################################

class CaissonBody(RigidBody):
    """
    Sub-class to create a caisson rigid body.
    """

    def __init__(self, shape):
        super(CaissonBody, self).__init__(shape)
        self.Shape = shape
        self.barycenter = shape.barycenter


    def getDisplacement(self, dt):
        # acceleration from force
        self.acceleration = self.getAcceleration()
        # substeps for smoother motion between timesteps
        dt_sub = dt/float(substeps)
        for i in range(substeps):
            # displacement
            #self.velocity += self.acceleration*dt_sub
            #self.h += self.velocity*dt_sub
            self.h[:], self.velocity[:] = forward_euler(p0=self.h, v0=self.velocity,
                                                             a=self.acceleration, dt=dt_sub)



    def getAngularDisplacement(self, dt):
        # angular acceleration from moment
        self.ang_acc = self.getAngularAcceleration()
        dt_sub = dt/float(substeps)
        self.ang_disp = np.array([0.,0.,0.])
        for i in range(substeps):
            # rotation
            #self.ang_vel += self.ang_acc*dt_sub
            #self.ang_disp += self.ang_vel*dt_sub
            self.ang_disp[:], self.ang_vel[:] = forward_euler(p0=self.ang_disp, v0=self.ang_vel,
                                                         a=self.ang_acc, dt=dt_sub)


    def step(self, dt, substeps=20):
        """
        Step for rigid body calculations in Python

        Parameters
        ----------
        dt: float
            time step
        """
        nd = self.Shape.Domain.nd
        # reinitialise displacement values
        self.ang_disp[:] = np.zeros(3)
        self.h[:] = np.zeros(3)

        # Translational motion calculation
        if self.friction != True:
            self.getDisplacement(dt)
        else:
            self.friction_module(dt)

        # Rotational motion calculation
        if self.overturning != True:
            self.getAngularDisplacement(dt)
        else:
            self.overturning_module(dt)

        # translate
        self.Shape.translate(self.h[:nd])
        # rotate
        if nd ==2:
            self.ang = self.ang_disp[2]
        else:
            self.ang = np.linalg.norm(self.ang_disp)
        if self.ang != 0.:
            if self.overturning != True:
                self.Shape.rotate(self.ang, self.ang_vel, self.Shape.barycenter)
            else:
                self.Shape.rotate(self.ang, self.ang_vel, self.pivot_friction)
            self.rotation[:nd, :nd] = self.Shape.coords_system
            self.rotation_matrix[:] = np.dot(np.linalg.inv(self.last_rotation),
                                             self.rotation)
        else:
            self.rotation_matrix[:] = np.eye(3)
        self.barycenter[:] = self.Shape.barycenter
        self.position[:] = self.Shape.barycenter

        # update vertices for friction and overturning modules
        self.cV[0] = self.Shape.vertices[0]
        self.cV[1] = self.Shape.vertices[1]
        self.cV[2] = self.Shape.vertices[2]
        self.cV[3] = self.Shape.vertices[3]

        #import pdb; pdb.set_trace()

        # record values
        if self.record_values is True:
            self.recordValues()

    def recordValues(self):
        """
        Records values of rigid body attributes at each time step in a csv file.
        """
        comm = Comm.get()
        if comm.isMaster():
            t_last = self.model.stepController.t_model_last
            dt_last = self.model.levelModelList[-1].dt_last
            values_towrite = []
            t = t_last-dt_last
            if t == 0:
                headers = []
                if self.record_dict['time'] is True:
                    headers += ['t']
                if self.record_dict['pos'] is True:
                    headers += ['x', 'y', 'z']
                if self.record_dict['caissonVertices'] is True:
                    headers += ['1', '2', '3', '4']
                if self.record_dict['rot'] is True:
                    headers += ['rx', 'ry', 'rz']
                if self.record_dict['ang_disp'] is True:
                    headers += ['ang_dispx', 'ang_dispy', 'ang_dispz']
                if self.record_dict['F'] is True:
                    headers += ['Fx', 'Fy', 'Fz']
                if self.record_dict['M'] is True:
                    headers += ['Mx', 'My', 'Mz']
                if self.record_dict['inertia'] is True:
                    headers += ['inertia']
                if self.record_dict['vel'] is True:
                    headers += ['vel_x', 'vel_y', 'vel_z']
                if self.record_dict['acc'] is True:
                    headers += ['acc_x', 'acc_y', 'acc_z']
                if self.record_dict['ang_vel'] is True:
                    headers += ['ang_velx', 'ang_vely', 'ang_velz']
                if self.record_dict['ang_acc'] is True:
                    headers += ['ang_accx', 'ang_accy', 'ang_accz']
                if self.record_dict['rp'] is True:
                    headers += ['rpx', 'rpy', 'rpz']
                if self.record_dict['Msoil'] is True:
                    headers += ['Msoil']
                with open(self.record_file, 'w') as csvfile:
                    writer = csv.writer(csvfile, delimiter=',')
                    writer.writerow(headers)
            if self.record_dict['time'] is True:
                t = t_last-dt_last
                values_towrite += [t]
            if self.record_dict['pos'] is True:
                x, y, z = self.last_position
                values_towrite += [x, y, z]
            if self.record_dict['caissonVertices'] is True:
                cV1, cV2, cV3, cV4 = self.cV_last
                values_towrite += [cV1, cV2, cV3, cV4]
            if self.record_dict['rot'] is True:
                rot = self.last_rotation
                rx = atan2(rot[1, 2], rot[2, 2])
                ry = -asin(rot[0, 2])
                rz = atan2(rot[0, 1], rot[0, 0])
                values_towrite += [rx, ry, rz]
            if self.record_dict['ang_disp'] is True:
                ang_dispx, ang_dispy, ang_dispz = self.ang_disp
                values_towrite += [ang_dispx, ang_dispy, ang_dispz]
            if self.record_dict['F'] is True:
                Fx, Fy, Fz = self.F
                values_towrite += [Fx, Fy, Fz]
            if self.record_dict['M'] is True:
                Mx, My, Mz = self.M
                values_towrite += [Mx, My, Mz]
            if self.record_dict['inertia'] is True:
                values_towrite += [self.inertia]
            if self.record_dict['vel'] is True:
                vel_x, vel_y, vel_z = self.velocity
                values_towrite += [vel_x, vel_y, vel_z]
            if self.record_dict['acc'] is True:
                acc_x, acc_y, acc_z = self.acceleration
                values_towrite += [acc_x, acc_y, acc_z]
            if self.record_dict['ang_vel'] is True:
                ang_velx, ang_vely, ang_velz = self.ang_vel
                values_towrite += [ang_velx, ang_vely, ang_velz]
            if self.record_dict['ang_acc'] is True:
                ang_accx, ang_accy, ang_accz = self.ang_acc
                values_towrite += [ang_accx, ang_accy, ang_accz]
            if self.record_dict['rp'] is True:
                rpx, rpy, rpz = self.rp
                values_towrite += [rpx, rpy, rpz]
            if self.record_dict['Msoil'] is True:
                Msoil = self.Msoil
                values_towrite += [Msoil]
            with open(self.record_file, 'a') as csvfile:
                writer = csv.writer(csvfile, delimiter=',')
                writer.writerow(values_towrite)

    def getInertia(self, vec=(0., 0., 1.), pivot=None):
        """
        Gives the inertia of the shape from an axis and a pivot

        Parameters
        ----------
        vec: array_like
            Vector around which the body rotates.
        pivot: Optional[array_like]
            Pivotal point around which the body rotates. If not set, it will
            be the barycenter coordinates

        Returns
        -------
        I: float
            inertia of the mass

        Notes
        -----
        The inertia is calculated relative to the coordinate system of the
        shape (self.coords_system). If the shape was not initialised with a
        position corresponding to its inertia tensor (e.g. shape was already
        rotated when initialised), set the coordinate system accordingly
        before calling this function
        """
        assert self.It is not None, 'No inertia tensor! (' + self.name + ')'
        if pivot is None:
            pivot = self.barycenter
        # Pivot coords relative to shape centre of mass
        distance = pivot-np.array(self.barycenter)
        # making unity vector/axis of rotation
        vec = vx, vy, vz = np.array(vec)
        length_vec = np.sqrt(vx**2+vy**2+vz**2)
        vec = vec/length_vec
        if self.Shape.Domain.nd == 2:
            L, H = self.Shape.dim
            dx, dy, dz = distance                             # To calculate It when pivot != barycenter
            Ix = float((H**2)/12.) + float(dy**2)             # To calculate x-component of It when pivot != barycenter
            Iy = float((L**2)/12.) + float(dx**2)             # To calculate y-component of It when pivot != barycenter
            self.It = Ix + Iy
            I = self.It*self.mass
        elif self.Shape.Domain.nd == 3:
            # vector relative to original position of shape:
            vec = np.dot(vec, np.linalg.inv(self.coords_system))
            cx, cy, cz = vec
            # getting the tensor for calculaing moment of inertia
            # from arbitrary axis
            vt = np.array([[cx**2, cx*cy, cx*cz],
                           [cx*cy, cy**2, cy*cz],
                           [cx*cz, cy*cz, cz**2]])
            # total moment of inertia
            I = np.einsum('ij,ij->', self.mass*self.It, vt)
        return I

    def setRecordValues(self, filename=None, all_values=False, time=True,
                        pos=False, caissonVertices=None, rot=False, ang_disp=False, F=False, M=False, inertia=False,
                        vel=False, acc=False, ang_vel=False, ang_acc=False, rp=False, Msoil=False):
        """
        Sets the rigid body attributes that are to be recorded in a csv file
        during the simulation.

        Parameters
        ----------
        filename: Optional[string]
            Name of file, if not set, the file will be named as follows:
            'record_[shape.name].csv'
        all_values: bool
            Set to True to record all values listed below.
        time: bool
            Time of recorded row (default: True).
        pos: bool
            Position of body (default: False. Set to True to record).
        caissonVertices: list of points
            Caisson vertices position (default: False. Set to True to record).
        rot: bool
            Rotation of body (default: False. Set to True to record).
        ang_disp: array
            Angular displecement calculated during rigid body calculation step.
            Applied on the body in order to make it rotating.
        F: bool
            Forces applied on body (default: False. Set to True to record).
        M: bool
            Moments applied on body (default: False. Set to True to record).
        inertia: bool
            Inertia of body (default: False. Set to True to record).
        vel: bool
            Velocity of body (default: False. Set to True to record).
        acc: bool
            Acceleration of body (default: False. Set to True to record).
        ang_vel: array
            Angular velocity of body (default: False. Set to True to record).
        ang_acc: array
            Angular acceleration of body (default: False. Set to True to record).
        rp: array
            Distance between barycenter and pivot (rotation centre).
            It's for moment transformation.
        Msoil: array
            Moment offered by soil.
        """
        self.record_values = True
        if pos is True:
            x = y = z = True
        if caissonVertices is True:
            cV1 = cV2 = cV3 = cV4 = True
        if rot is True:
            rot_x = rot_y = rot_z = True
        if ang_disp is True:
            ang_dispx =  ang_dispy = ang_dispz = True
        if F is True:
            Fx = Fy = Fz = True
        if M is True:
            Mx = My = Mz = True
        if vel is True:
            vel_x = vel_y = vel_z = True
        if acc is True:
            acc_x = acc_y = acc_z = True
        if ang_vel is True:
            ang_velx = ang_vely = ang_velz = True
        if ang_acc is True:
            ang_accx = ang_accy = ang_accz = True
        if rp is True:
            rpx, rpy, rpz = True
        if Msoil is True:
            Msoil = True

        self.record_dict = {'time':time, 'pos': pos, 'caissonVertices':caissonVertices, 'rot':rot, 'ang_disp':ang_disp, 'F':F, 'M':M,
                            'inertia': inertia, 'vel': vel, 'acc': acc, 'ang_vel': ang_vel, 'ang_acc': ang_acc, 'rp': rp, 'Msoil':Msoil}

        if all_values is True:
            for key in self.record_dict:
                self.record_dict[key] = True
        if filename is None:
            self.record_filename = 'record_' + self.name + '.csv'
        else:
            self.record_filename = filename + '.csv'


#  Friction, overturning and soil modules

    def setFriction(self, friction, m_static, m_dynamic, tolerance, grainSize):
        """
        Sets material properties for sliding and overturning modules

        Parameters
        ----------
        friction: string
            If True, friction module is switched on.
        m_static: float
            Static friction parameter (see Coulomb equation).
        m_dynamic: float
            Dynamic friction parameter.
        tolerance: float
            It's used to check if the body is in rotated state or not.
            It imposes a tolerance limit for the difference between the vertical coordinates of the 2 bottom vertices of the rigid body.
        grainSize: float
            Typical grain size of the rubble mound (if exists!) under the caisson.
            It offers an extra check of the body position.
            If the rigid body lower point position is higher than this value, it is a floating body.
        """
        self.friction = friction
        self.m_static = m_static
        self.m_dynamic = m_dynamic
        self.tolerance = tolerance
        self.grainSize = grainSize



    def setOverturning(self, overturning):
        """
        Sets material properties for sliding and overturning modules

        Parameters
        ----------
        overturning: string
            If True, overturning module is switched on.
        """
        self.overturning = overturning



    def setSprings(self, springs, Kx, Ky, Krot, C, Crot):
        """
        Sets a system of uniform springs to model soil's reactions (for moving bodies)

        Parameters
        ----------
        spring: string
            If True, spring module is switched on.
        Kx: float
            horizontal stiffness
        Ky: float
            vertical stiffness
        Krot: float
            rotational stiffness
        C: float
            damping parameter
        Crot: float
            rotational damping parameter
        """
        self.springs = springs
        self.Kx = Kx
        self.Ky = Ky
        self.Krot = Krot
        self.C = C
        self.Crot = Crot



    def soil_force(self, reaction):
        """
        Calculate force reactions offered  by soil foundation (springs model).

        Parameters
        ----------
        reaction : Reaction flag.
            fx = Horizontal force is returned.
            fy = Vertical force is returned.
        """
        # spring
        Lcais, Hcais = self.Shape.dim
        Kx = self.Ky        # Pa
        Ky = self.Ky        # Pa
        Krot = self.Krot    # N
        C = self.C          # Pa s
        Crot = self.Crot
        l_caisson = Lcais * 0.5
        pos_x, pos_y, pos_z = self.last_position

      ### Horizontal reaction
        dx = (pos_x-self.init_barycenter[0])/self.init_barycenter[0]
        fsx = -Kx*dx*Lcais                 # spring
        fdx = -C*self.velocity[0]          # damper
        fx = fsx + fdx

      ### Vertical reaction
        dy = (pos_y-self.init_barycenter[1])/self.init_barycenter[1]
        fsy = -Ky*dy*Lcais                 # spring
        fdy = -10*C*self.velocity[1]       # damper
        fy = fsy + fdy
        if reaction == 'fx':
            return fx
        elif reaction == 'fy':
            return fy



    def soil_moment(self, h):
        """
        Calculate moment reactions offered  by soil foundation (springs model).

        Parameters
        ----------
        h : Caisson's vertex flag.
            The choice is between vertices at the bottom layer in contact with the soil foundation.
            It's used to check position of the caisson and spring deformation.
            h = 0 represents body in unrotated position.
            h = 1 represents body in rotated position, vertex_1 inside porous medium, positive rotation.
            h = 2 represents body in rotated position, vertex_2 inside porous medium, negative rotation.
        """
        # spring
        Lcais, Hcais = self.Shape.dim
        Kx = self.Kx        # Pa
        Ky = self.Ky        # Pa
        Krot = self.Krot    # N
        C = self.C          # Pa s
        Crot = self.Crot
        l_caisson = Lcais * 0.5
        pos_x, pos_y, pos_z = self.last_position

      ### Rotational reaction
        if Krot != 0.0:
            rot = self.rotation
            rot_z = atan2(rot[0, 1], rot[0, 0])    # rotation from position
            msz = -Krot*rot_z                      # spring
        else:
        # ------ left side
            Def1 = (self.cV_last[0][1]-self.cV_init[0][1])/self.cV_init[0][1] # non-dimensional
            # extra check
            if Def1 > 0.0:
                Def1=0.0
            Fspr1 = -Ky*Def1*l_caisson/2. # semplified version, using costant Ky, formula is for the triangle area
            lever_arm1 = -l_caisson*2./3.
            Mspr1 = Fspr1*lever_arm1
        # ------ right side
            Def2 = (self.cV_last[1][1]-self.cV_init[1][1])/self.cV_init[1][1] # non-dimensional
            # extra check
            if Def2 > 0.0:
                Def2=0.0
            Fspr2 = -Ky*Def2*l_caisson/2. # semplified version, using costant Ky, formula is for the triangle area
            lever_arm2 = l_caisson*2./3.
            Mspr2 = Fspr2*lever_arm2
            # ----- springs total moment
            if h == 0:
                Mspr1 = 0.0
                Mspr2 = 0.0
            elif h == 1:
                Mspr2 = 0.0
            elif h == 2:
                Mspr1 = 0.0
            msz = Mspr1 + Mspr2             # spring
        if h==0:
            Crot=0.0

        mdz = -Crot*self.last_ang_vel[2]          # damper
        return msz, mdz



    def friction_module(self,dt):
        """
        Calculate sliding motion modelling frictional force.

        Parameters
        ----------
        dt : Time step.
        """
        nd = self.Shape.Domain.nd
        substeps = 20
        dt_sub = dt/float(substeps)
        # movement_functions for friction test cases
        Fx, Fy, Fz = self.F
        eps = 0.000000000000000001 # to avoid 0/0
        mass = self.mass
        sign_static = Fx/(abs(Fx)+eps)
        sign_dynamic = self.last_velocity[0]/(abs(self.last_velocity[0])+eps)
        if nd==2:
            g=np.array([0.,-9.81,0.])
            Fv = Fy
            gv = g[1]
        if nd==3:
            g=np.array([0.,0.,-9.81])
            Fv = Fz
            gv = g[2]
        acceleration = np.zeros(3)

        #---------------------------------------------------------------
        def static_case(self, sign, Fx, Fv, mass, m):
            """
            Set a static friction.

            Parameters
            ----------
            sign : It's function of horizontal force.
                It's used to calculate frictional force.
            Fx : Total horizontal force from rigid body calculation (wave loading).
            Fy : Total vertical force from rigid body calculation (wave loading + weight of the body).
            mass : Mass of the rigid body.
            m : static friction factor.
            """
            Ftan = -sign*m*abs(Fv)
            if abs(Fx)<abs(Ftan):
                self.acceleration = np.zeros(3)
                self.velocity = np.zeros(3)
                self.h[:] = np.zeros(3)
            else:
                fs = self.soil_force(reaction='fx')
                Fh = Fx+Ftan#+fs
                self.acceleration[0] = Fh/mass
                self.acceleration[1] = 0.0
                self.acceleration[2] = 0.0
                self.velocity += self.acceleration*dt
                self.h += self.velocity*dt
                    #self.h[:], self.velocity[:] = forward_euler(p0=self.h, v0=self.velocity,
                    #                                         a=self.acceleration, dt=dt_sub)
            self.fromDynamic_toStatic = False

        #---------------------------------------------------------------
        def dynamic_case(self, sign, Fx, Fv, mass, m):
            """
            Set a dynamic friction.

            Parameters
            ----------
            sign : It's function of horizontal force.
                It's used to calculate frictional force.
            Fx : Total horizontal force from rigid body calculation (wave loading).
            Fy : Total vertical force from rigid body calculation (wave loading + weight of the body).
            mass : Mass of the rigid body.
            m : dynamic friction factor.
            """
            Ftan = -sign*m*abs(Fv)
            fs = self.soil_force(reaction='fx')
            Fh = Fx+Ftan#+fs
            self.acceleration[0] = Fh/mass
            self.acceleration[1] = 0.0
            self.acceleration[2] = 0.0
            for i in range(substeps):
                self.velocity += self.acceleration*dt_sub
                self.h += self.velocity*dt_sub
                #self.h[:], self.velocity[:] = forward_euler(p0=self.h, v0=self.velocity,
                #                                             a=self.acceleration, dt=dt_sub)
                # When velocity changes sign, it means that 0-condition is passed
                # Loop must start from static case again
                self.fromDynamic_toStatic = False
                velCheck=self.velocity[0]
                last_velCheck=self.last_velocity[0]
                if (last_velCheck*velCheck) < 0.0:
                    self.fromDynamic_toStatic = True
                    break



        #---------------------------------------------------------------

        if (Fv*gv)>0:
        #--- Friction module, static case
            if self.last_velocity[0] == 0.0 or self.last_fromDynamic_toStatic ==True:
                static_case(self, sign_static, Fx, Fv, mass, m=self.m_static)
        #--- Friction module, dynamic case
            else :
                dynamic_case(self, sign_dynamic, Fx, Fv, mass, m=self.m_dynamic)

        if (Fv*gv)<0:
        #--- Floating module, static case
            if self.last_velocity[0] == 0.0 or self.last_fromDynamic_toStatic ==True:
                static_case(self, sign_static, Fx, Fv, mass, m=0.0)
        #--- Floating module, dynamic case
            else :
                dynamic_case(self, sign_dynamic, Fx, Fv, mass, m=0.0)



    def overturning_module(self,dt):
        """
        Calculate overturning motion modelling soil foundation reactions.

        Parameters
        ----------
        dt : Time step.
        """
        nd = self.Shape.Domain.nd
        substeps = 20
        dt_sub = dt/float(substeps)


        #---------------------------------------------------------------
        def calculate_rotation(self, floating, h):
            """
            Calculate rotation.

            Parameters
            ----------
            floating : enable floating body calculation with NO restrictions offered by friction or springs.
            h : Caisson's vertex flag.
                The choice is between vertices at the bottom layer in contact with the soil foundation.
                It's used to check position of the caisson and spring deformation.
                h = 0 represents body in unrotated position.
                h = 1 represents body in rotated position, vertex_1 inside porous medium, positive rotation.
                h = 2 represents body in rotated position, vertex_2 inside porous medium, negative rotation.
            """

            if floating == True:
                self.pivot_friction = self.Shape.barycenter
            else: # medium point in the bottom of the caisson between vertex1 and vertex2
                self.pivot_friction = np.array((0.5*(self.cV[0][0]+self.cV[1][0]), 0.5*(self.cV[0][1]+self.cV[1][1]), 0.0), dtype=float)

            # angular acceleration from moment
            self.rp = (self.pivot_friction-self.Shape.barycenter)
            rpx, rpy, rpz = self.rp
            Fx, Fy, Fz = self.F
            Mpivot = np.array([(rpy*Fz-rpz*Fy), -(rpx*Fz-rpz*Fx), (rpx*Fy-rpy*Fx)]) # moment transformation calculated in pivot
            Mp = self.M - Mpivot                                                    # moment transformation

            # moment equilibrium -- M = I*angAcc + C*angVel + K*angDispl
            # all these terms can be written in function of the rotation (angDispl)

            if floating == True or self.springs==False:
                Kx = 0.0
                Ky = 0.0
                Krot = 0.0
                C = 0.0
                Crot = 0.0
            else:
                Kx = self.Kx
                Ky = self.Ky
                Krot = self.Krot
                C = self.C
                Crot = self.Crot

            self.inertia = self.getInertia(self.Mp, self.pivot_friction)
            assert self.inertia != 0, 'Zero inertia: inertia tensor (It)' \
                                      'was not set correctly!'
            rot = self.rotation
            phi0 =  atan2(rot[0, 1], rot[0, 0])    # rotation from position, last timestep
            # rotation (phi) calculation. Inversion of the moment equilibrium
            num = Mp[2]*(dt**2) + Crot*dt*phi0 + self.inertia*(phi0+self.last_ang_vel[2]*dt)
            den = Crot*dt + Krot*(dt**2) + self.inertia
            phi = float(num/den)

            # velocity and acceleration updated from rotation calculation
            self.ang_disp[2] = phi-phi0
            self.ang_vel[2]  = self.ang_disp[2] / dt
            self.ang_acc[2]  = (self.ang_vel[2] - self.last_ang_vel[2]) / dt

        #---------------------------------------------------------------

        # Check position and then call calculate rotation
        # h1 and h2 in the same horizontal line, flat
        if abs(self.cV_last[0][1]-self.cV_last[1][1]) < self.tolerance:

            # caisson inside rubble mound
            if self.cV_last[0][1] < self.cV_init[0][1]:
                transl_x = 0.0
                transl_y = self.cV_init[0][1] - self.cV_last[0][1]
                transl_z = 0.0
                transl = (transl_x, transl_y)
                self.Shape.translate(transl)
                calculate_rotation(self, floating=False, h=0)     # ----> caisson under the rubble mound: I move it again on the mound!!!

            # caisson on the rubble mound (still in contact)
            elif self.cV_last[0][1] == self.cV_init[0][1] or (self.cV_last[0][1] - self.cV_init[0][1]) < self.grainSize:
                calculate_rotation(self, floating=False, h=0)     # ----> caisson on the rubble mound: I calculate normally rotation

            # floating caisson (no contact)
            else:
                calculate_rotation(self, floating=True, h=0)     # ----> caisson up the rubble mound: FLOATING CASE I calculate normally rotation on barycenter!!!


        # h1 and h2 NOT in a horizontal line, NOT flat
        else:
            if self.cV_last[1][1] < self.cV_last[0][1]:
                if self.cV_last[1][1] == self.cV_init[1][1] or (self.cV_last[1][1] - self.cV_init[1][1]) < self.grainSize:
                    calculate_rotation(self, floating=False, h=2)        # ----> caisson on the rubble mound: I calculate normally rotation
                else:
                    calculate_rotation(self, floating=True, h=0)     # ----> caisson up the rubble mound: FLOATING CASE I calculate normally rotation on barycenter!!!

            elif self.cV_last[0][1] < self.cV_last[1][1]:
                if self.cV_last[0][1] == self.cV_init[0][1] or (self.cV_last[0][1] - self.cV_init[0][1]) < self.grainSize:
                    calculate_rotation(self, floating=False, h=1)        # ----> caisson on the rubble mound: I calculate normally rotation
                else:
                    calculate_rotation(self, floating=True, h=0)     # ----> caisson up the rubble mound: FLOATING CASE I calculate normally rotation on barycenter!!!

        # check variables
        #import pdb; pdb.set_trace()



# Calculation step

    def calculate_init(self):
        """
        Function called at the very beginning of the simulation by proteus.
        """
        nd = self.Shape.Domain.nd
        shape = self.Shape

        self.position = np.zeros(3)
        self.position[:] = self.Shape.barycenter.copy()
        self.last_position[:] = self.position
        self.velocity = np.zeros(3, 'd')
        self.last_velocity = np.zeros(3, 'd')
        self.acceleration = np.zeros(3, 'd')
        self.last_acceleration = np.zeros(3, 'd')

        self.ang_disp = np.zeros(3, 'd')
        self.last_ang_disp = np.zeros(3, 'd')
        self.rotation = np.eye(3)
        self.rotation[:nd, :nd] = shape.coords_system
        self.last_rotation = np.eye(3)
        self.last_rotation[:nd, :nd] = shape.coords_system
        self.ang_vel = np.zeros(3, 'd')
        self.last_ang_vel = np.zeros(3, 'd')
        self.ang_acc = np.zeros(3, 'd')
        self.last_ang_acc = np.zeros(3, 'd')
        self.F = np.zeros(3, 'd')
        self.M = np.zeros(3, 'd')
        self.last_F = np.zeros(3, 'd')
        self.last_M = np.zeros(3, 'd')
        self.ang = 0.
        self.barycenter = self.Shape.barycenter

        if nd == 2:
            self.Fg = self.mass*np.array([0., -9.81, 0.])
        if nd == 3:
            self.Fg = self.mass*np.array([0., 0., -9.81])
        if self.record_values is True:
            self.record_file = os.path.join(Profiling.logDir,
                                            self.record_filename)
        # Initial position of the 2D caisson vertices
        self.cV_init = np.array([(self.Shape.vertices[0][0], self.Shape.vertices[0][1]),
                                 (self.Shape.vertices[1][0], self.Shape.vertices[1][1]),
                                 (self.Shape.vertices[2][0], self.Shape.vertices[2][1]),
                                 (self.Shape.vertices[3][0], self.Shape.vertices[3][1])])
        # Position of the 2D caisson vertices
        self.cV =  np.array([(self.Shape.vertices[0][0], self.Shape.vertices[0][1]),
                             (self.Shape.vertices[1][0], self.Shape.vertices[1][1]),
                             (self.Shape.vertices[2][0], self.Shape.vertices[2][1]),
                             (self.Shape.vertices[3][0], self.Shape.vertices[3][1])])
        # Last position of the 2D caisson vertices
        self.cV_last =  np.array([(self.Shape.vertices[0][0], self.Shape.vertices[0][1]),
                                 (self.Shape.vertices[1][0], self.Shape.vertices[1][1]),
                                 (self.Shape.vertices[2][0], self.Shape.vertices[2][1]),
                                 (self.Shape.vertices[3][0], self.Shape.vertices[3][1])])
        # friciton module parameter used for switching to dynamic motion cases
        self.fromDynamic_toStatic = False
        self.last_fromDynamic_toStatic = False
        # friction and overturning parameters to be initialised
        self.pivot_friction = np.zeros(3)
        self.last_pivot_friction = np.zeros(3)
        self.Ftan = 0.0
        self.last_Ftan = 0.0
        self.Mp = np.zeros(3, 'd')
        self.last_Mp = np.zeros(3, 'd')
        self.rp = np.zeros(3, 'd')
        self.last_rp = np.zeros(3, 'd')
        self.init_barycenter = self.Shape.Domain.barycenters[1]
        self.Msoil = np.zeros(3, 'd')
        self.last_Msoil = np.zeros(3, 'd')

    def calculate(self):
        """
        Function called at each time step by proteus.
        """
        # store previous values
        self.last_position[:] = self.position
        self.last_velocity[:] = self.velocity
        self.last_acceleration[:] = self.acceleration
        self.last_rotation[:] = self.rotation
        self.last_ang_acc[:] = self.ang_acc
        self.last_ang_vel[:] = self.ang_vel
        self.last_ang_disp[:] = self.ang_disp
        self.last_F[:] = self.F
        self.last_M[:] = self.M
        # friciton and overturning
        self.cV_last[:] = self.cV
        self.last_fromDynamic_toStatic = self.fromDynamic_toStatic
        self.last_pivot_friction = self.pivot_friction
        self.last_Mp[:] = self.Mp
        self.last_rp = self.rp
        self.last_Msoil[:] = self.Msoil
        # for first time step
        try:
            dt = self.model.levelModelList[-1].dt_last
        except:
            dt = self.dt_init
        # update forces and moments for current body/shape
        i0, i1 = self.i_start, self.i_end
        # get forces
        F = self.getTotalForce()
        M = self.getTotalMoments()
        self.F[:] = F2 = F*self.free_x
        # get moments
        self.M[:] = M2 = M*self.free_r
        # store F and M with DOF constraints to body
        # calculate new properties
        self.step(dt)
        # log values
        t_previous = self.model.stepController.t_model_last-dt
        t_current = self.model.stepController.t_model_last
        h = self.h
        last_pos, pos = self.last_position, self.position
        last_vel, vel = self.last_velocity, self.velocity
        rot = self.rotation
        rot_x = atan2(rot[1, 2], rot[2, 2])
        rot_y = -asin(rot[0, 2])
        rot_z = atan2(rot[0, 1], rot[0, 0])
        logEvent("================================================================")
        logEvent("=================== Rigid Body Calculation =====================")
        logEvent("================================================================")
        logEvent("Name: " + `self.Shape.name`)
        logEvent("================================================================")
        logEvent("[proteus]     t=%1.5fsec to t=%1.5fsec" % \
            (t_previous, t_current))
        logEvent("[proteus]    dt=%1.5fsec" % (dt))
        logEvent("[body] ============== Pre-calculation attributes  ==============")
        logEvent("[proteus]     t=%1.5fsec" % (t_previous))
        logEvent("[proteus]     F=(% 12.7e, % 12.7e, % 12.7e)" % (F[0], F[1], F[2]))
        logEvent("[proteus] F*DOF=(% 12.7e, % 12.7e, % 12.7e)" % (F2[0], F2[1], F2[2]))
        logEvent("[proteus]     M=(% 12.7e, % 12.7e, % 12.7e)" % (M[0], M[1], M[2]))
        logEvent("[proteus] M*DOF=(% 12.7e, % 12.7e, % 12.7e)" % (M2[0], M2[1], M2[2]))
        logEvent("[body]      pos=(% 12.7e, % 12.7e, % 12.7e)" % \
            (last_pos[0], last_pos[1], last_pos[2]))
        logEvent("[body]      vel=(% 12.7e, % 12.7e, % 12.7e)" % \
            (last_vel[0], last_vel[1], last_vel[2]))
        logEvent("[body] ===============Post-calculation attributes ==============")
        logEvent("[body]        t=%1.5fsec" % (t_current))
        logEvent("[body]        h=(% 12.7e, % 12.7e, % 12.7e)" % (h[0], h[1], h[2]))
        logEvent("[body]      pos=(% 12.7e, % 12.7e, % 12.7e)" % \
            (pos[0], pos[1], pos[2]))
        logEvent("[body]      vel=(% 12.7e, % 12.7e, % 12.7e)" % \
            (vel[0], vel[1], vel[2]))
        logEvent("[body]      rot=(% 12.7e, % 12.7e, % 12.7e)" % \
            (rot_x, rot_y, rot_z))
        logEvent("================================================================")




import numpy as np
from proteus import AuxiliaryVariables


class RigidBody(AuxiliaryVariables.AV_base):
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
        self.rotation = np.eye(3)
        self.rotation[:nd, :nd] = shape.coords_system
        self.last_rotation = np.eye(3)
        self.last_rotation[:nd, :nd] = shape.coords_system
        self.F = np.zeros(3, 'd')
        self.M = np.zeros(3, 'd')
        self.last_F = np.zeros(3, 'd')
        self.last_M = np.zeros(3, 'd')
        self.ang = 0.
        self.barycenter = self.Shape.barycenter
        self.angvel = np.zeros(3, 'd')
        self.last_angvel = np.zeros(3, 'd')
        if nd == 2:
            self.Fg = self.mass*np.array([0., -9.81, 0.])
        if nd == 3:
            self.Fg = self.mass*np.array([0., 0., -9.81])
        if self.record_values is True:
            self.record_file = os.path.join(Profiling.logDir,
                                            self.record_filename)

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
        M = getMoments()
        M_t = np.sum(M)
        return M_t

    def getTotalForce(self):
        F_p = getPressureForces()
        F_v = getShearForces()
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
            ang_acc = self.M[:]/self.inertia
        else:
            self.inertia = None
            ang_acc = np.array([0., 0., 0.])
        return ang_acc

    def calculate(self):
        """
        Function called at each time step by proteus.
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
        F = getTotalForce()
        M = getTotalMoments()
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
        self.h[:] = np.zeros(3)
        # acceleration from force
        self.acceleration = getAcceleration()
        # angular acceleration from moment
        ang_acc = getAngularAcceleration()
        # substeps for smoother motion between timesteps
        dt_sub = dt/float(substeps)
        for i in range(substeps):
            # displacement
            self.h[:], self.velocity[:] = forward_euler(p0=self.h, v0=self.velocity,
                                                        a=self.acceleration, dt=dt_sub)
            # rotation
            ang_disp, self.angvel[:] = forward_euler(p0=ang_disp, v0=self.angvel,
                                                     a=ang_acc, dt=dt_sub)
        # translate
        self.Shape.translate(self.h[:nd])
        # rotate
        self.ang = np.linalg.norm(ang_disp)
        if nd == 2 and self.angvel[2] < 0:
            self.ang = -self.ang
        if self.ang != 0.:
            self.Shape.rotate(self.ang, self.angvel, self.Shape.barycenter)
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
        Set mass of the shape and calculate density if volume is defined.

        Parameters
        ----------
        mass: float
            mass of the body
        """
        self.mass = float(mass)
        if self.volume:
            self.density = self.mass/self.volume

    def setDensity(self, density):
        """
        Set density and calculate mass is volume is defined.

        Parameters
        ----------
        density: float
            Density of the shape
        """
        self.density = float(density)
        if self.volume:
            self.mass = self.density*self.volume

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
        length_vec = sqrt(vx**2+vy**2+vz**2)
        vec = vec/length_vec
        if self.Domain.nd == 2:
            I = self.It*self.mass
        elif self.Domain.nd == 3:
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
                        pos=False, rot=False, F=False, M=False, inertia=False,
                        vel=False, acc=False):
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

        """
        self.record_values = True
        if pos is True:
            x = y = z = True
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
        self.record_dict = {'time':time, 'pos': pos, 'rot':rot, 'F':F, 'M':M,
                            'inertia': inertia, 'vel': vel, 'acc': acc}
        if all_values is True:
            for key in self.record_dict:
                self.record_dict[key] = True
        if filename is None:
            self.record_filename = 'record_' + self.name + '.csv'
        else:
            self.record_filename = filename + '.csv'



def forward_euler(p0, v0, a, dt):
    v1 = v0+a*dt
    p1 = p0+v1*dt
    return p1, v1

def leapfrog(p0, v0, a, dt):
    p1 = p0+v0*dt+0.5*a*dt**2
    v1 = v0+a*dt
    return p1, v1

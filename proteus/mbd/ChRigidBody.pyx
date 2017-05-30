# distutils: language = c++

import os
import sys
import csv
cimport numpy as np
import numpy as np
import mpi4py as MPI
from scipy import spatial
from proteus import AuxiliaryVariables, Archiver, Comm, Profiling
from proteus import SpatialTools as st
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool
from libcpp.memory cimport (shared_ptr,
                            make_shared)
from collections import OrderedDict
from cython.operator cimport dereference as deref
# chrono C++ headers
cimport ChronoHeaders as ch
# chrono Python headers
from proteus.mbd cimport pyChronoCore as pych


cdef extern from "ChMoorings.h":
    cdef cppclass cppMesh:
        shared_ptr[ch.ChMesh] mesh
        void SetAutomaticGravity(bool val)
    cppMesh * newMesh(ch.ChSystemSMC&, shared_ptr[ch.ChMesh])
    cdef cppclass cppCable:
        ch.ChSystemSMC& system
        ch.ChMesh& mesh
        vector[shared_ptr[ch.ChNodeFEAxyzDD]] nodes
        vector[shared_ptr[ch.ChElementBeamANCF]] elems
        double L0
        double length
        int nb_elems
        vector[ch.ChVector] mvecs
        void buildNodes()
        void buildMaterials()
        void buildElements()
        void buildMesh()
        void updateDragForces()
        void updateBuoyancyForces()
        void setDragCoefficients(double Cd_axial, double Cd_normal)
        void setAddedMassCoefficients(double Cd_axial, double Cd_normal)
    cdef cppclass cppMultiSegmentedCable:
        ch.ChSystemSMC& system
        ch.ChMesh& mesh
        vector[shared_ptr[cppCable]] cables
        vector[shared_ptr[ch.ChNodeFEAxyzDD]] nodes
        vector[shared_ptr[ch.ChElementBeamANCF]] elems
        shared_ptr[ch.ChLinkPointFrame] constraint_back
        shared_ptr[ch.ChLinkPointFrame] constraint_front
        void buildNodes()
        void buildCable()
        # void setVelocityAtNodes(double* fluid_velocity)
        void attachFrontNodeToBody(shared_ptr[ch.ChBody])
        void attachBackNodeToBody(shared_ptr[ch.ChBody])
        void updateDragForces()
        void updateAddedMassForces()
        void applyForces()
        void updateBuoyancyForces()
        void setFluidVelocityAtNodes(vector[ch.ChVector] fluid_velocity)
        void setFluidAccelerationAtNodes(vector[ch.ChVector] fluid_acceleration)
        void setFluidDensityAtNodes(vector[double] fluid_density)
        void setContactMaterial(shared_ptr[ch.ChMaterialSurfaceSMC] material)
        ch.ChVector getTensionElement(int i)
    cppMultiSegmentedCable * newMoorings(ch.ChSystemSMC& system,
                                         shared_ptr[ch.ChMesh] mesh,
                                         vector[double] length,
                                         vector[int] nb_elems,
                                         vector[double] d,
                                         vector[double] rho,
                                         vector[double] E,
                                         string beam_type
        )
    cdef cppclass cppSurfaceBoxNodesCloud:
        ch.ChSystemSMC& system
        ch.ChVector position
        ch.ChVector dimensions
        shared_ptr[ch.ChBodyEasyBox] body;
        void setNodesSize(double size)
    cppSurfaceBoxNodesCloud * newSurfaceBoxNodesCloud(ch.ChSystemSMC& system,
                                                shared_ptr[ch.ChMesh] mesh,
                                                ch.ChVector position,
                                                ch.ChVector dimensions)


cdef extern from "ChRigidBody.h":
    cdef cppclass cppSystem:
        ch.ChSystemSMC system
        void DoStepDynamics(dt)
        void step(double proteus_dt, int n_substeps)
        void recordBodyList()
        void setChTimeStep(double dt)
        void setGravity(double* gravity)
        void setDirectory(string directory)
        void setTimestepperType(string tstype, bool verbose)
    cppSystem * newSystem(double* gravity)
    cdef cppclass cppRigidBody:
        shared_ptr[ch.ChBody] body
        double mass
        ch.ChVector pos
        ch.ChVector pos_last
        ch.ChVector vel
        ch.ChVector vel_last
        ch.ChVector acc
        ch.ChVector acc_last
        ch.ChVector angvel
        ch.ChVector angvel_last
        ch.ChVector angacc
        ch.ChVector angacc_last
        # ChVector inertia
        ch.ChMatrix33 rotm
        ch.ChMatrix33 rotm_last
        ch.ChQuaternion rotq
        ch.ChQuaternion rotq_last
        # double* free_x
        # double* free_r
        ch.ChVector F
        ch.ChVector F_last
        ch.ChVector M
        ch.ChVector M_last
        cppRigidBody(cppSystem* system)
        void prestep(double* force, double* torque)
        void poststep()
        ch.ChVector hxyz(double* x, double dt)
        double hx(double* x, double dt)
        double hy(double* x, double dt)
        double hz(double* x, double dt)
        void addSpring(double stiffness,
                       double damping,
                       double* fairlead,
                       double* anchor,
                       double rest_length)
        void addPrismaticLinksWithSpring(double* pris1,
                                         double* pris2,
                                         double stiffness,
                                         double damping,
                                         double rest_length);
        void setRotation(double* quat)
        void setPosition(double* pos)
        void setConstraints(double* free_x, double* free_r)
        void setInertiaXX(double* inertia)
        void setName(string name)
        void setPrescribedMotionPoly(double coeff1)
        void setPrescribedMotionSine(double a, double f)


    cppRigidBody * newRigidBody(cppSystem* system)

cdef class ProtChBody:
    cdef cppRigidBody * thisptr
    cdef ch.ChQuaternion rotation
    cdef ch.ChQuaternion rotation_last
    cdef public:
      str record_file
      object model
      ProtChSystem ProtChSystem
      object Shape
      int nd, i_start, i_end
      double dt
      double width_2D
      object record_dict
      object prescribed_motion_function
      pych.ChBody ChBody
      np.ndarray position
      np.ndarray position_last
      np.ndarray F
      np.ndarray M
      np.ndarray F_last
      np.ndarray M_last
      np.ndarray F_prot
      np.ndarray M_prot
      np.ndarray F_prot_last
      np.ndarray M_prot_last
      np.ndarray F_applied
      np.ndarray M_applied
      np.ndarray F_applied_last
      np.ndarray M_applied_last
      np.ndarray acceleration
      np.ndarray acceleration_last
      np.ndarray velocity
      np.ndarray velocity_last
      np.ndarray ang_acceleration_last
      np.ndarray ang_acceleration
      np.ndarray ang_velocity_last
      np.ndarray ang_velocity
      double ang_vel_norm
      double ang_vel_norm_last
      np.ndarray barycenter0
      np.ndarray rotation_init
      np.ndarray rotm
      np.ndarray rotm_last
      np.ndarray rotq
      np.ndarray rotq_last
      string name
      # np.ndarray free_r
      # np.ndarray free_x
    def __cinit__(self,
                  ProtChSystem system,
                  shape=None,
                  nd=None):
        self.ProtChSystem = system
        self.thisptr = newRigidBody(system.thisptr)
        self.ChBody = pych.ChBody()
        self.thisptr.body = self.ChBody.sharedptr_chbody  # give pointer to cpp class
        self.ProtChSystem.thisptr.system.AddBody(self.thisptr.body)
        self.Shape = None
        self.setRotation(np.array([1.,0.,0.,0.]))  # initialise rotation (nan otherwise)
        self.attachShape(shape)  # attach shape (if any)
        self.ProtChSystem.addSubcomponent(self)  # add body to system (for pre and post steps)
        self.record_dict = OrderedDict()
        self.F_prot = np.zeros(3)  # initialise empty Proteus force
        self.M_prot = np.zeros(3)  # initialise empty Proteus moment
        self.F_applied = np.zeros(3)  # initialise empty Applied force
        self.M_applied = np.zeros(3)  # initialise empty Applied moment
        self.prescribed_motion_function = None

        self.velocity = np.zeros(3)
        self.velocity_last = np.zeros(3)
        self.ang_velocity = np.zeros(3)
        self.position_last = self.position
        self.ang_vel_norm = 0.  # used for mesh disp prediction
        self.ang_vel_norm_last = 0.  # used for mesh disp prediction
        # if self.nd is None:
        #     assert nd is not None, "must set nd if SpatialTools.Shape is not passed"
        #     self.nd = nd

    def attachShape(self, shape):
        """Attach proteus.SpatialTools shape to body.
        Used for automatic calculation of external forces from Proteus.
        Called automatically when creating a body and passing a shape
        instance.

        Parameters
        ----------
        shape: SpatialTools.Shape
            instance of shape from proteus.SpatialTools or
            proteus.mprans.SpatialTools
        """
        if shape is not None:
            assert self.Shape is None, 'Shape '+self.Shape.name+' was already attached'
            self.Shape = shape
            if 'ChRigidBody' not in shape.auxiliaryVariables:
                shape.auxiliaryVariables['ChRigidBody'] = self
                self.setName(shape.name)
            self.nd = shape.Domain.nd
            self.SetPosition(shape.barycenter)

    def SetBodyFixed(self, bool state):
        """Fix body in space

        Parameters
        ----------
        state: bool
            body fixed if True, body free if False
        """
        deref(self.thisptr.body).SetBodyFixed(state)

    def setWidth2D(self, width):
        """Sets width of 2D body (for forces and moments calculation)

        Parameters
        ----------
        width: float
            width of the body
        """
        self.width_2D = width

    def set_indices(self, i_start, i_end):
        """Sets the flags of the boundaries of the body
        numbers must be gloabal (from domain.segmentFlags or
        domain.facetFlags) and a range from i_start to i_end.

        Parameters
        ----------
        i_start: int
            first global flag of body boundaries
        i_end: int
            last global flag (+1) of body boundaries
        """
        self.i_start = i_start
        self.i_end = i_end

    def attachAuxiliaryVariables(self,avDict):
        pass

    def setInertiaXX(self, np.ndarray inertia):
        self.thisptr.setInertiaXX(<double*> inertia.data)

    def setInitialRot(self, rot):
        cdef np.ndarray zeros = np.zeros(3)
        self.rotation_init = rot
        self.thisptr.prestep(<double*> zeros.data,
                             <double*> zeros.data)
        if self.rotation_init is not None:
            Profiling.logEvent("$$$$$ SETTING ROT")
            self.thisptr.setRotation(<double*> self.rotation_init.data)
        self.thisptr.poststep()

    def hxyz(self, np.ndarray x, double t):
        cdef np.ndarray h
        cdef np.ndarray xx
        cdef double ang, and_last
        cdef double dt_half, dt_half_last
        cdef np.ndarray d_tra, d_tra_last # translational displacements
        cdef np.ndarray d_rot, d_rot_last # rotational displacements
        cdef np.ndarray h_body  # displacement from body
        cdef ch.ChVector h_body_vec
        h = np.zeros(3)
        if self.ProtChSystem.scheme == "CSS":
            h_body_vec = self.thisptr.hxyz(<double*> x.data, t)
            h_body = pych.ChVector_to_npArray(h_body_vec)
            h += h_body
        elif self.ProtChSystem.scheme == "ISS":
            dt_half = self.ProtChSystem.dt_fluid/2.
            dt_half_last = self.ProtChSystem.dt_fluid_last/2.
            # remove previous prediction
            # translate back first
            d_tra_last = -self.velocity_last*dt_half_last
            h += d_tra_last
            # rotate back
            ang_last = -self.ang_vel_norm*dt_half_last
            if ang > 0:
                d_rot_last = (st.rotation3D(points=x+h,  # (translated back)
                                            rot=ang_last,
                                            axis=self.ang_velocity_last,
                                            pivot=self.position_last)
                              -self.position_last)
                h += d_rot_last
            # add rigid body displacement
            xx = x+h # previous position of body
            h_body_vec = self.thisptr.hxyz(<double*> xx.data, t)
            h_body = pych.ChVector_to_npArray(h_body_vec)
            h += h_body
            # add current prediction
            # rotate first
            ang = self.ang_vel_norm*dt_half
            if ang > 0:
                d_rot = (st.rotation3D(points=x+h,
                                       rot=ang,
                                       axis=self.ang_velocity,
                                       pivot=self.position)
                         -self.position)
                h += d_rot
            # translate
            d_tra = self.velocity*dt_half
            h += d_tra
        return h

    def hx(self, np.ndarray x, double t):
        """BC function for mesh nodes displacement (x component)

        Parameters
        ----------
        x: array_like
            coordinates of the node before displacement
        t: double
            simulation time
        """
        return self.hxyz(x, t)[0]

    def hy(self, np.ndarray x, double t):
        """BC function for mesh nodes displacement (y component)

        Parameters
        ----------
        x: array_like
            coordinates of the node before displacement
        t: double
            simulation time
        """
        return self.hxyz(x, t)[1]

    def hz(self, np.ndarray x, double t):
        """BC function for mesh nodes displacement (z component)

        Parameters
        ----------
        x: array_like
            coordinates of the node before displacement
        t: double
            simulation time
        """
        return self.hxyz(x, t)[2]

    def addSpring(self, double stiffness, double damping, np.ndarray fairlead,
                  np.ndarray anchor, double rest_length):
        self.thisptr.addSpring(stiffness, damping, <double*> fairlead.data,
                               <double*> anchor.data, rest_length)

    def SetPosition(self, np.ndarray position):
        """Sets position of body manually

        Parameters
        ----------
        position: array_like
            new position of body (must be array of length 3)
        """
        self.position = position
        self.thisptr.setPosition(<double*> position.data)

    def SetPosition(self, np.ndarray position):
        """Sets position of body manually

        Parameters
        ----------
        position: array_like
            new position of body (must be array of length 3)
        """
        self.thisptr.setPosition(<double*> position.data)

    def setRotation(self, np.ndarray quaternion):
        """Sets rotation of body manually

        Parameters
        ----------
        rotation: array_like
            new rotation of body (quaternion: must be array of length 4)
        """
        self.thisptr.setRotation(<double*> quaternion.data)

    def setConstraints(self, np.ndarray free_x, np.ndarray free_r):
        """Sets constraints on the body
        (!) Only acts on Proteus and gravity forces

        Parameters
        ----------
        free_x: array_like
            Translational constraints.
        free_r: array_like
            Rotational constraints.
        """
        self.thisptr.setConstraints(<double*> free_x.data, <double*> free_r.data)

    def SetMass(self, double mass):
        """Sets mass of body.

        Parameters
        ----------
        mass: double
            mass of the body
        """
        deref(self.thisptr.body).SetMass(mass)

    def getPressureForces(self):
        """Gives pressure forces from fluid (Proteus) acting on body.
        (!) Only works during proteus simulation

        Returns
        -------
        F_p: array_like
            pressure forces (x, y, z) as provided by Proteus
        """
        i0, i1 = self.i_start, self.i_end
        F_p = self.ProtChSystem.model.levelModelList[-1].coefficients.netForces_p[i0:i1, :]
        F_t = np.sum(F_p, axis=0)
        return F_t

    def getShearForces(self):
        """Gives shear forces from fluid (Proteus) acting on body
        (!) Only works during proteus simulation

        Returns
        -------
        F_v: array_like
            shear forces (x, y, z) as provided by Proteus
        """
        i0, i1 = self.i_start, self.i_end
        F_v = self.ProtChSystem.model.levelModelList[-1].coefficients.netForces_v[i0:i1, :]
        F_t = np.sum(F_v, axis=0)
        return F_t

    def getMoments(self):
        """Gives moments from fluid (Proteus) acting on body
        (!) Only works during proteus simulation

        Returns
        -------
        M: array_like
            moments (x, y, z) as provided by Proteus
        """
        i0, i1 = self.i_start, self.i_end
        M = self.ProtChSystem.model.levelModelList[-1].coefficients.netMoments[i0:i1, :]
        M_t = np.sum(M, axis=0)
        # !!!!!!!!!!!! UPDATE BARYCENTER !!!!!!!!!!!!
        Fx, Fy, Fz = self.F_prot
        rx, ry, rz = self.barycenter0-self.ChBody.GetPos()
        Mp = np.array([ry*Fz-rz*Fy, -(rx*Fz-rz*Fx), (rx*Fy-ry*Fx)])
        M_t += Mp
        return M_t

    def SetPos(self, double[:] pos):
        """Sets current position of body

        Parameters
        ----------
        pos: array_like
            position of body (array of length 3)
        """
        assert len(pos) == 3, 'Position aray must be of length 3'
        cdef ch.ChVector pos_vec
        pos_vec = ch.ChVector[double](pos[0], pos[1], pos[2])
        deref(self.thisptr.body).SetPos(pos_vec)

    def SetRot(self, double[:] rot):
        """Sets current rotation (quaternion) of body

        Parameters
        ----------
        rot: array_like
           rotation of body (array of length 4)
        """
        assert len(rot) == 4, 'Position aray must be of length 4'
        cdef ch.ChQuaternion rot_vec
        rot_vec = ch.ChQuaternion[double](rot[0], rot[1], rot[2], rot[3])
        deref(self.thisptr.body).SetRot(rot_vec)

    def getRotationMatrix(self):
        """Gives current rotation (matrix) of body

        Returns
        -------
        rot: array_like
            current rotation (matrix) of body
        """
        x0 = self.thisptr.rotm.Get_A_Xaxis().x()
        x1 = self.thisptr.rotm.Get_A_Xaxis().y()
        x2 = self.thisptr.rotm.Get_A_Xaxis().z()
        y0 = self.thisptr.rotm.Get_A_Yaxis().x()
        y1 = self.thisptr.rotm.Get_A_Yaxis().y()
        y2 = self.thisptr.rotm.Get_A_Yaxis().z()
        z0 = self.thisptr.rotm.Get_A_Zaxis().x()
        z1 = self.thisptr.rotm.Get_A_Zaxis().y()
        z2 = self.thisptr.rotm.Get_A_Zaxis().z()
        rot = np.array([x0, x1, x2],
                       [y0, y1, y2],
                       [z0, z1, z2])
        return rot

    def prestep(self):
        """Called before Chrono system step.
        Sets external forces automatically from Proteus solution.
        """


        self.storeValues()
        if self.ProtChSystem.model is not None:
            self.F_prot = self.getPressureForces()+self.getShearForces()
            self.M_prot = self.getMoments()
            if self.ProtChSystem.first_step is True:
                # just apply initial conditions for 1st time step
                self.F_applied = self.F_prot
                self.M_applied = self.M_prot
            else:
                # actual force applied to body
                if self.ProtChSystem.scheme == "CSS":
                    self.F_applied = self.F_prot
                    self.M_applied = self.M_prot
                elif self.ProtChSystem.scheme == "ISS":
                    F_bar = (self.F_prot+self.F_prot_last)/2.
                    M_bar = (self.M_prot+self.M_prot_last)/2.
                    self.F_applied = F_bar
                    self.M_applied = M_bar
                    # self.F_applied = self.F_prot
                    # self.M_applied = self.M_prot
                    # self.F_applied = 2*F_bar - self.F_applied_last
                    # self.M_applied = 2*M_bar - self.M_applied_last
        self.setExternalForces(self.F_applied, self.M_applied)

    def setExternalForces(self, np.ndarray forces, np.ndarray moments):
        """Sets external forces to body.
        Called during prestep or can be called manually. If called manually,
        must be a Chrono only simulation.

        Parameters
        ----------
        forces: array_like
            forces array (length 3)
        moments: array_like
            moments array (length 3)
        """
        self.thisptr.prestep(<double*> forces.data,
                             <double*> moments.data)

    def poststep(self):
        """Called after Chrono system step.
        Records values to csv, broadcast new position and rotation from
        calculating processor to all processors for moving mesh BC.
        """
        if self.prescribed_motion_function is not None:
            new_x = self.callPrescribedMotion(self.ProtChSystem.model.stepController.t_model_last)
            self.thisptr.setPosition(<double*> new_x.data)
        self.thisptr.poststep()
        self.getValues()
        comm = Comm.get().comm.tompi4py()
        cdef ch.ChQuaternion rotq
        cdef ch.ChQuaternion rotq_last
        cdef ch.ChVector pos
        cdef ch.ChVector pos_last
        cdef double e0, e1, e2, e3, e0_last, e1_last, e2_last, e3_last
        cdef double posx, posy, posz, posx_last, posy_last, posz_last
        if self.ProtChSystem.parallel_mode is True:
            # need to broadcast values to all processors on the C++ side
            # comm.Barrier()
            self.rotq = comm.bcast(self.rotq,
                                  self.ProtChSystem.chrono_processor)
            self.rotq_last = comm.bcast(self.rotq_last,
                                       self.ProtChSystem.chrono_processor)
            self.position = comm.bcast(self.position,
                                  self.ProtChSystem.chrono_processor)
            self.position_last = comm.bcast(self.position_last,
                                       self.ProtChSystem.chrono_processor)
            if comm.rank == 1:
                self._recordValues()
        else:
            if comm.rank == 0:
                self._recordValues()
        # need to pass position and rotation values to C++ side
        # needed for transformations when calling hx, hy, hz, hxyz
        e0, e1, e2, e3 = self.rotq
        e0_last, e1_last, e2_last, e3_last = self.rotq_last
        posx, posy, posz = self.position
        posx_last, posy_last, posz_last = self.position_last
        pos = ch.ChVector[double](posx, posy, posz)
        pos_last = ch.ChVector[double](posx_last, posy_last, posz_last)
        rotq = ch.ChQuaternion[double](e0, e1, e2, e3)
        rotq_last = ch.ChQuaternion[double](e0_last, e1_last, e2_last, e3_last)
        self.thisptr.rotq = rotq
        self.thisptr.rotq_last = rotq_last
        self.thisptr.pos = pos
        self.thisptr.pos_last = pos_last

    def calculate_init(self):
        """Called from self.ProtChSystem.calculate_init()
        before simulation starts
        """
        # barycenter0 used for moment calculations
        if self.Shape is not None:
            self.barycenter0 = self.Shape.barycenter.copy()
        # get the initial values for F and M
        cdef np.ndarray zeros = np.zeros(3)
        self.setExternalForces(zeros, zeros)
        self.thisptr.poststep()
        # get first, store then on initial time step
        self.getValues()
        self.storeValues()
        # self.thisptr.setRotation(<double*> self.rotation_init.data)
        #

    def calculate(self):
        pass

    def setPrescribedMotionSine(self, double a, double f):
        """Sets sinusoidal prescribed motion for body

        Parameters
        ----------
        a: double
            amplitude of sinusoid
        f: double
            frequency of sinusoid
        """
        self.thisptr.setPrescribedMotionSine(a, f)

    def setPrescribedMotionPoly(self, double coeff1):
        """Sets polynomial prescribed motion for body

        Parameters
        ----------
        coeff1: double
            coeff1 of polynomial
        """
        self.thisptr.setPrescribedMotionPoly(coeff1)

    def setPrescribedMotion(self, function):
        """Sets custom prescribed motion function
        (!) should be preferably set only if body is free and not
        linked to other bodies or other elements (such as moorings)
        as this is used only for setting moving mesh BC by enforcing
        position of the body at each time step.
        Use setPrescribedMotionPoly or setPrescribedMotionSine for
        functions that are safe to use with a body linked with other
        elements.

        Parameters
        ----------
        function:
            must be a function of time returning an array of the
            absolute position of the body (numpy array of length 3:
            x, y, z)
        """
        self.prescribed_motion_function = function

    cdef np.ndarray callPrescribedMotion(self, double t):
        return self.prescribed_motion_function(t)

    def storeValues(self):
        self.velocity_last = self.velocity
        self.position_last = self.position
        self.acceleration_last = self.acceleration
        self.rotq_last = self.rotq
        self.rotm_last = self.rotm
        self.ang_acceleration_last = self.ang_acceleration
        self.ang_velocity_last = self.ang_velocity
        self.F_last = self.F
        self.M_last = self.M
        self.ang_vel_norm_last = self.ang_vel_norm
        # force from fluid at current time step
        self.F_prot_last = np.array(self.F_prot)
        self.M_prot_last = np.array(self.M_prot)
        self.F_applied_last = np.array(self.F_applied)
        self.M_applied_last = np.array(self.M_applied)
        if self.ProtChSystem.parallel_mode is True:
            comm = Comm.get().comm.tompi4py()
            self.position_last = comm.bcast(self.position_last,
                                            self.ProtChSystem.chrono_processor)
            self.velocity_last = comm.bcast(self.velocity_last,
                                            self.ProtChSystem.chrono_processor)
            self.acceleration_last = comm.bcast(self.acceleration_last,
                                                self.ProtChSystem.chrono_processor)
            self.rotq_last = comm.bcast(self.rotq_last,
                                        self.ProtChSystem.chrono_processor)
            self.rotm_last = comm.bcast(self.rotm_last,
                                        self.ProtChSystem.chrono_processor)
            self.ang_velocity_last = comm.bcast(self.ang_velocity_last,
                                                self.ProtChSystem.chrono_processor)
            self.ang_acceleration_last = comm.bcast(self.ang_acceleration_last,
                                                    self.ProtChSystem.chrono_processor)
            self.ang_vel_norm_last = comm.bcast(self.ang_vel_norm_last,
                                                self.ProtChSystem.chrono_processor)
            self.F_prot_last = comm.bcast(self.F_prot_last,
                                          self.ProtChSystem.chrono_processor)
            self.M_prot_last = comm.bcast(self.M_prot_last,
                                          self.ProtChSystem.chrono_processor)

    def getValues(self):
        """Get values (pos, vel, acc, etc.) from C++ to python
        """
        # position
        self.position = self.ChBody.GetPos()
        # rotation
        self.rotq = self.ChBody.GetRot()
        self.rotm = self.ChBody.GetA()
        # acceleration
        self.acceleration = self.ChBody.GetPos_dtdt()
        #velocity
        self.velocity = self.ChBody.GetPos_dt()
        # angular acceleration
        self.ang_acceleration = self.ChBody.GetWacc_loc()
        # angular velocity
        self.ang_velocity = self.ChBody.GetWvel_loc()
        # norm of angular velocity
        self.ang_vel_norm = np.sqrt(self.ang_velocity[0]**2
                                    +self.ang_velocity[1]**2
                                    +self.ang_velocity[2]**2)
        # force
        self.F = pych.ChVector_to_npArray(self.thisptr.F)
        # moment
        self.M = pych.ChVector_to_npArray(self.thisptr.M)
        if self.ProtChSystem.parallel_mode is True:
            comm = Comm.get().comm.tompi4py()
            self.position = comm.bcast(self.position, self.ProtChSystem.chrono_processor)
            self.velocity = comm.bcast(self.velocity, self.ProtChSystem.chrono_processor)
            self.acceleration = comm.bcast(self.acceleration, self.ProtChSystem.chrono_processor)
            self.rotq = comm.bcast(self.rotq, self.ProtChSystem.chrono_processor)
            self.rotm = comm.bcast(self.rotm, self.ProtChSystem.chrono_processor)
            self.ang_velocity = comm.bcast(self.ang_velocity, self.ProtChSystem.chrono_processor)
            self.ang_acceleration = comm.bcast(self.ang_acceleration, self.ProtChSystem.chrono_processor)
            self.ang_vel_norm = comm.bcast(self.ang_vel_norm, self.ProtChSystem.chrono_processor)


    def setRecordValues(self, filename=None, all_values=False, pos=False,
                        rot=False, ang_disp=False, F=False, M=False,
                        inertia=False, vel=False, acc=False, ang_vel=False, ang_acc=False):
        """
        Sets the body attributes that are to be recorded in a csv file
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
        Notes
        -----
        To add another value manually, add to dictionary self.record_dict:
        key: header of the column in .csv
        value: list of length 2: [variable name, index within variable]
                                                 (if no index, use None)
        e.g. self.record_dict['m']['mass', None]
        """
        if all_values is True:
            pos = rot = F = M = acc = vel = ang_acc = ang_vel = True
        if pos is True:
            self.record_dict['x'] = ['position_last', 0]
            self.record_dict['y'] = ['position_last', 1]
            self.record_dict['z'] = ['position_last', 2]
        if rot is True:
            self.record_dict['rotq_e0'] = ['rotq_last', 0]
            self.record_dict['rotq_e1'] = ['rotq_last', 1]
            self.record_dict['rotq_e2'] = ['rotq_last', 2]
            self.record_dict['rotq_e3'] = ['rotq_last', 3]
        if F is True:
            self.record_dict['Fx'] = ['F', 0]
            self.record_dict['Fy'] = ['F', 1]
            self.record_dict['Fz'] = ['F', 2]
            self.record_dict['Fx_prot'] = ['F_prot', 0]
            self.record_dict['Fy_prot'] = ['F_prot', 1]
            self.record_dict['Fz_prot'] = ['F_prot', 2]
            Fx = Fy = Fz = True
        if M is True:
            self.record_dict['Mx'] = ['M', 0]
            self.record_dict['My'] = ['M', 1]
            self.record_dict['Mz'] = ['M', 2]
            self.record_dict['Mx_prot'] = ['M_prot', 0]
            self.record_dict['My_prot'] = ['M_prot', 1]
            self.record_dict['Mz_prot'] = ['M_prot', 2]
        if acc is True:
            self.record_dict['ax'] = ['acceleration_last', 0]
            self.record_dict['ay'] = ['acceleration_last', 1]
            self.record_dict['az'] = ['acceleration_last', 2]
        if vel is True:
            self.record_dict['ux'] = ['velocity_last', 0]
            self.record_dict['uy'] = ['velocity_last', 1]
            self.record_dict['uz'] = ['velocity_last', 2]
        if ang_acc is True:
            self.record_dict['ang_ax'] = ['ang_acceleration_last', 0]
            self.record_dict['ang_ay'] = ['ang_acceleration_last', 1]
            self.record_dict['ang_az'] = ['ang_acceleration_last', 2]
        if ang_vel is True:
            self.record_dict['ang_ux'] = ['ang_velocity_last', 0]
            self.record_dict['ang_uy'] = ['ang_velocity_last', 1]
            self.record_dict['ang_uz'] = ['ang_velocity_last', 2]
        if inertia is True:
            self.record_dict['intertia'] = ['inertia', None]

    def _recordValues(self):
        """Records values of body attributes in a csv file.
        """
        if self.Shape is not None:
            self.record_file = os.path.join(Profiling.logDir, 'record_' + self.Shape.name + '.csv')
        else:
            self.record_file = os.path.join(Profiling.logDir, 'record_' + 'body' + '.csv')
        if self.ProtChSystem.model is not None:
            t_last = self.ProtChSystem.model.stepController.t_model_last
            try:
                dt_last = self.ProtChSystem.model.levelModelList[-1].dt_last
            except:
                dt_last = 0
            t = t_last-dt_last
        else:
            t = self.ProtChSystem.thisptr.system.GetChTime()
        t_prot = Profiling.time()-Profiling.startTime
        values_towrite = [t, t_prot]
        if t == 0:
            headers = ['t', 't_prot']
            for key in self.record_dict:
                headers += [key]
            with open(self.record_file, 'w') as csvfile:
                writer = csv.writer(csvfile, delimiter=',')
                writer.writerow(headers)
        for key, val in self.record_dict.iteritems():
            if val[1] is not None:
                values_towrite += [getattr(self, val[0])[val[1]]]
            else:
                values_towrite += [getattr(self, val[0])]
        with open(self.record_file, 'a') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerow(values_towrite)

    def addPrismaticLinksWithSpring(self, np.ndarray pris1,
                                    np.ndarray pris2, double stiffness, double damping,
                                    double rest_length):
        """
        fairlead: barycenter coords
        pris: absolute coords
        pris1-------fairlead(barycenter)
        |
        |
        |
        |
        pris2
        """
        self.thisptr.addPrismaticLinksWithSpring(<double*> pris1.data, 
                                                 <double*> pris2.data,
                                                 stiffness,
                                                 damping,
                                                 rest_length)

    def setName(self, string name):
        """Sets name of body (used for csv file)

        Parameters
        ----------
        name: str
            name of the body
        """
        self.name = name
        self.thisptr.setName(name)


cdef class ProtChSystem:
    cdef cppSystem * thisptr
    cdef public object model
    cdef object subcomponents
    cdef public double dt_init
    cdef double proteus_dt
    cdef double proteus_dt_last
    cdef double chrono_dt
    cdef string directory
    cdef object u
    cdef int nd
    cdef object femSpace_velocity
    cdef object femSpace_pressure
    cdef object nodes_kdtree
    cdef int min_nb_steps
    cdef double dt_fluid
    cdef double dt_fluid_last
    cdef public:
        bool build_kdtree
        bool parallel_mode
        int chrono_processor
        bool first_step
        string scheme

    def __cinit__(self, np.ndarray gravity, int nd=3, dt_init=0.001):
        self.thisptr = newSystem(<double*> gravity.data)
        self.subcomponents = []
        self.dt_init = 0.001
        self.model = None
        self.nd = nd
        self.build_kdtree = False
        comm = Comm.get().comm.tompi4py()
        if comm.Get_size() > 1:
            parallel_mode = True
            chrono_processor = 1
        else:
            parallel_mode = False
            chrono_processor = 0
        self.parallel_mode = parallel_mode
        self.chrono_processor = chrono_processor
        self.min_nb_steps = 1  # minimum number of chrono substeps
        self.proteus_dt = 0.
        self.first_step = True  # just to know if first step
        self.scheme = "CSS"
        self.dt_fluid = 0

    def GetChTime(self):
        """Gives time of Chrono system simulation

        Returns
        -------
        time: double
            time of chrono system
        """
        time = self.thisptr.system.GetChTime()
        return time

    def setCouplingScheme(self, string scheme):
        assert scheme == "CSS" or scheme == "ISS", "Coupling scheme requested unknown"
        self.scheme = scheme

    def attachModel(self, model, ar):
        """Attaches Proteus model to auxiliary variable
        """
        self.model = model
        return self

    def attachAuxiliaryVariables(self,avDict):
        pass

    def setMinimumSubsteps(self, int nb):
        """Sets the minimum nb of chrono substeps per proteus step
        if prot_dt=0.001 and ch_dt=0.002, there will be <nb>
        substeps of chrono instead of just 1.

        Parameters
        ----------
        nb: int
            Minimum number of chrono substeps.
        """
        self.min_nb_steps = nb

    def step(self, dt):
        self.dt_fluid_last = self.dt_fluid
        self.dt_fluid = dt
        if self.scheme == "ISS":
            dt = (self.dt_fluid+self.dt_fluid_last)/2.
        # calculate number of time steps
        nb_steps = max(int(dt/self.chrono_dt), 1)
        if nb_steps < self.min_nb_steps:
            nb_steps = self.min_nb_steps
        # solve Chrono system
        Profiling.logEvent('Solving Chrono system from t='
                           +str(self.thisptr.system.GetChTime())
                           +' with dt='+str(dt)
                           +' and '+str(nb_steps)+' substeps')


        comm = Comm.get().comm.tompi4py()
        if comm.rank == self.chrono_processor:
            self.thisptr.step(<double> dt, n_substeps=nb_steps)
        Profiling.logEvent('Solved Chrono system to t='
                           +str(self.thisptr.system.GetChTime()))

    def calculate(self, proteus_dt=None):
        """Does chrono system calculation for a Proteus time step
        Calls prestep and poststep on all subcomponents (bodies,
        moorings, etc) attached to the system.

        Parameters
        ----------
        proteus_dt: Optional[proteus_dt]
            Manually sets a time step. The time step is set
            automatically when coupled with a Proteus simulation
        """
        self.proteus_dt_last = self.proteus_dt
        if self.model is not None:
            try:
                self.proteus_dt = self.model.levelModelList[-1].dt_last
            except:
                self.proteus_dt = self.dt_init
        elif proteus_dt is not None:
            self.proteus_dt = proteus_dt
        else:
            sys.exit('no time step set')
        if self.model is not None and self.build_kdtree is True:
            self.nodes_kdtree = spatial.cKDTree(self.model.levelModelList[-1].mesh.nodeArray)
        for s in self.subcomponents:
            s.prestep()
        self.step(self.proteus_dt)
        for s in self.subcomponents:
            s.poststep()
        if self.first_step is True:  # check if first step passed
            self.first_step = False

    def calculate_init(self):
        """Does chrono system initialisation
        (!) Must be called before the first calculate() call.
        Calls calculate_init and poststep on all subcomponents
        (bodies, moorings, etc) attached to the system.
        """
        self.directory = str(Profiling.logDir)+'/'
        self.thisptr.setDirectory(self.directory)
        if self.model is not None and self.build_kdtree is True:
            self.u = self.model.levelModelList[-1].u
            # finite element space (! linear for p, quadratic for velocity)
            self.femSpace_velocity = self.u[1].femSpace
            self.femSpace_pressure = self.u[0].femSpace
            self.nodes_kdtree = spatial.cKDTree(self.model.levelModelList[-1].mesh.nodeArray)
        for s in self.subcomponents:
            s.calculate_init()
        self.thisptr.system.SetupInitial()
        for s in self.subcomponents:
            s.poststep()

    def setTimestepperType(self, string tstype, bool verbose=False):
        """Change timestepper (default: Euler)

        Parameters
        ----------
        tstype: str
            type of timestepper ('Euler' or 'HHT')
        """
        tstypes = ["Euler", "HHT"]
        assert str(tstype) in tstypes, str(tstype)+" not a valid choice."
        self.thisptr.setTimestepperType(tstype, verbose)

    def setTimeStep(self, double dt):
        """Sets time step for Chrono solver.
        Calculations in Chrono will use this time step within the
        Proteus time step (if bigger)

        Parameters
        ----------
        dt: float
            Chrono time step size
        """
        self.chrono_dt = dt
        self.thisptr.setChTimeStep(dt)

    def setGravity(self, np.ndarray gravity):
        """Sets gravity acceleration of the Chrono system

        Parameters
        ----------
        gravity: array_like
            Gravity acceleration (array of length 3)
        """
        self.thisptr.setGravity(<double*> gravity.data)

    def addSubcomponent(self, subcomponent):
        """Adds subcomponent to system
        calculate_init() of subcomponent called before initial timestep
        prestep() and poststep of subcomponent() called at all timestep

        Parameters
        ----------
        subcomponent: class instance
            class instance of subcomponent
        """
        self.subcomponents += [subcomponent]

    def recordBodyList(self):
        comm = Comm.get().comm.tompi4py()
        if self.parallel_mode is True:
            if comm.rank == self.chrono_processor:
                self.thisptr.recordBodyList()
        else:
            if comm.rank == 0:
                self.thisptr.recordBodyList()

    def findELementContainingCoords(self, coords):
        """
        Parameters
        ----------
        coords: array_like
            global coordinates to look for

        Returns
        -------
        xi:
            local coordinates 
        eN: int
            (local) element number
        rank: int
            processor rank containing element
        """
        comm = Comm.get().comm.tompi4py()
        # get nearest node on each processor
        nearest_node, nearest_node_distance = getLocalNearestNode(coords, self.nodes_kdtree)
        # look for element containing coords on each processor (if it exists)
        local_element = getLocalElement(self.femSpace_velocity, coords, nearest_node)
        # check which processor has element (if any)
        haveElement = int(local_element is not None)
        if haveElement:
            owning_proc = comm.rank
        if local_element:
            # get local coords
            xi0 = self.femSpace_velocity.elementMaps.getInverseValue(local_element, coords)
            xi = comm.bcast(xi0, owning_proc)
            element = comm.bcast(local_element, owning_proc)
            rank = comm.bcast(owning_proc, owning_proc)
        return xi, element, rank

    def getFluidVelocityLocalCoords(self, xi, element, rank):
        """
        Parameters
        ----------
        xi: 
            local coords in element
        element: int
            element number (local to processor 'rank')
        rank: int
            rank of processor owning the element
        """
        comm = Comm.get().comm.tompi4py()
        if comm.rank == rank:
            u = self.u[1].getValue(element, xi)
            u = comm.bcast(u, rank)
            v = self.u[2].getValue(element, xi)
            v = comm.bcast(v, rank)
            if self.nd > 2:
                w = self.u[3].getValue(element, xi)
            if self.nd <= 2:
                w = 0
            # broadcast to all processors
            w = comm.bcast(w, rank)
        return u, v, w

    def getFluidVelocityGradientLocalCoords(self, xi, element, rank):
        comm = Comm.get().comm.tompi4py()
        if comm.rank == rank:
            u_grad = self.u_grad[1].getGradientValue(element, xi)
            u_grad = comm.bcast(u_grad, rank)
            v_grad = self.u[2].getGradientValue(element, xi)
            v_grad = comm.bcast(v_grad, rank)
            if self.nd > 2:
                w_grad = self.u[3].getGradientValue(element, xi)
            if self.nd <= 2:
                w_grad = 0
            # broadcast to all processors
            w_grad = comm.bcast(w_grad, rank)
        return u_grad, v_grad, w_grad

    def findFluidVelocityAtCoords(self, coords):
        """Finds solution from NS for velocity of fluid at given coordinates

        Parameters
        ----------
        coords: array_like
            coordinates at which velocity solution is sought

        Returns
        -------
        u: float
            velocity in the x direction
        v: float
            velocity in the y direction
        w: float
            velocity in the z direction (0 if 2D)
        """
        comm = Comm.get().comm.tompi4py()
        # get nearest node on each processor
        nearest_node, nearest_node_distance = getLocalNearestNode(coords, self.nodes_kdtree)
        # look for element containing coords on each processor (if it exists)
        local_element = getLocalElement(self.femSpace_velocity, coords, nearest_node)
        # check which processor has element (if any)
        haveElement = int(local_element is not None)
        if haveElement:
            owning_proc = comm.rank
        if local_element:
            # NEXT LINE TO CHANGE
            nd = self.nd
            # get local coords
            xi = self.femSpace_velocity.elementMaps.getInverseValue(local_element, coords)
            # get solution
            u = self.u[1].getValue(local_element, xi)
            v = self.u[2].getValue(local_element, xi)
            if nd > 2:
                w = self.u[3].getValue(local_element, xi)
            # broadcast to all processors
            u = comm.bcast(u, owning_proc)
            v = comm.bcast(v, owning_proc)
            if nd > 2:
                w = comm.bcast(w, owning_proc)
            if nd <= 2:
                w = 0
        else:
            sys.exit('{coords} outside of domain'.format(coords=str(coords)))
        return u, v, w

# ctypedef np.ndarray vecarray(ChVector)

# ctypedef np.ndarray (*ChVector_to_npArray) (ChVector)

#def testx():
#    cdef ChSystem system = ChSystem()
#    cdef ChBody bod = ChBody()
#    cdef ChVector oo = ChVector[double](2.,3.,4.)
#    bod.SetPos_dt(oo)
#    cdef ChVector& gg = bod.GetPos_dt()
#    print(gg.x, gg.y, gg.z)


cdef class Mesh:
    cdef cppMesh * thisptr
    def __cinit__(self, ProtChSystem system):
        cdef shared_ptr[ch.ChMesh] mesh = make_shared[ch.ChMesh]()
        self.thisptr = newMesh(system.thisptr.system, mesh)
    def setAutomaticGravity(self, bool val):
        self.thisptr.SetAutomaticGravity(val)

cdef class SurfaceBoxNodesCloud:
    cdef cppSurfaceBoxNodesCloud * thisptr
    def __cinit__(self, ProtChSystem system, Mesh mesh, np.ndarray position, np.ndarray dimensions):
        cdef ch.ChVector[double] pos = ch.ChVector[double](position[0], position[1], position[2])
        cdef ch.ChVector[double] dim = ch.ChVector[double](dimensions[0], dimensions[1], dimensions[2])
        self.thisptr = newSurfaceBoxNodesCloud(system.thisptr.system,
                                               mesh.thisptr.mesh,
                                               pos,
                                               dim)
        # self.System.addBody(self)
    def setNodesSize(self, double size):
        self.thisptr.setNodesSize(size)

# def linkMoorings(System system, Moorings mooringA, int nodeA_ind, Moorings mooringB, int nodeB_ind):
#     """
#     Parameters
#     ----------
#     system: System
#         class instance of the System where to add link
#     mooringA: Moorings
#         class instance of first mooring cable (A)
#     nodeA_ind: int
#         index of node to be used as link on mooring cable A
#     mooringB: Moorings
#         class instance of second mooring cable (B)
#     nodeB_ind: int
#         index of node to be used as link on mooring cable B
#     """
#     cdef shared_ptr[ch.ChLinkPointPoint] link = make_shared[ch.ChLinkPointPoint]()
#     deref(link).Initialize(<shared_ptr[ch.ChNodeFEAxyz]> mooringA.thisptr.nodes[nodeA_ind], <shared_ptr[ch.ChNodeFEAxyz]> mooringB.thisptr.nodes[nodeB_ind])
#     system.thisptr.system.Add(<shared_ptr[ch.ChPhysicsItem]> link)


cdef class ProtChMoorings:
    """Class for building mooring cables

    Parameters
    ----------
    system: System
        Class instance of the system.
    mesh: Mesh
        Class instance of the mesh.
    length: np.ndarray
        Length of cable segments. Must be an array, if the cable only
        has one type of segment (e.g. only one type of chain), an
        array of length 1 can be passed.
    nb_elems: np.ndarray
        Number of elements per segments.
    d: np.ndarray
        Diameter of segments.
    rho: np.ndarray
        Density of segments.
    E: np.ndarray
        Young's modulus of segments.
    beam_type: str
        Type of elements (default: "CableANCF").
    """
    cdef cppMultiSegmentedCable * thisptr
    cdef public:
      # vector[pyfea.ChNodeFEAxyzD] nodes
      # vector[pyfea.ChElementCableANCF] elems
      # vector[ProtChCable] cables
      str record_file
      object model
      ProtChSystem ProtChSystem
      object Mesh
      int nd
      object nodes_function
      ProtChBody body_front
      ProtChBody body_back
      bool front_body
      bool back_body
      bool nodes_built
      bool external_forces_from_ns
      bool external_forces_manual
      np.ndarray fluid_density_array
      np.ndarray fluid_velocity_array
      np.ndarray fluid_velocity_array_previous
      np.ndarray fluid_acceleration_array
      string name
      string beam_type
    def __cinit__(self,
                  ProtChSystem system,
                  Mesh mesh,
                  double[:] length,
                  np.ndarray nb_elems,
                  double[:] d,
                  double[:] rho,
                  double[:] E,
                  string beam_type="CableANCF"):
        check_arrays = [len(length), len(nb_elems), len(d), len(rho), len(E)]
        assert all(v == len(length) for v in check_arrays), 'arrays are not of same length'
        self.ProtChSystem = system
        self.ProtChSystem.addSubcomponent(self)
        self.nd = self.ProtChSystem.nd
        self.Mesh = mesh
        self.beam_type = beam_type
        cdef vector[double] vec_length
        cdef vector[int] vec_nb_elems
        cdef vector[double] vec_d
        cdef vector[double] vec_rho
        cdef vector[double] vec_E
        for i in range(len(length)):
            vec_length.push_back(length[i])
            vec_nb_elems.push_back(nb_elems[i])
            vec_d.push_back(d[i])
            vec_rho.push_back(rho[i])
            vec_E.push_back(E[i])
        self.thisptr = newMoorings(system.thisptr.system,
                                   mesh.thisptr.mesh,
                                   vec_length,
                                   vec_nb_elems,
                                   vec_d,
                                   vec_rho,
                                   vec_E,
                                   beam_type
                                   )
        self.nodes_function = lambda s: (s, s, s)
        self.nodes_built = False
        self.name = 'record_moorings'
        self.external_forces_from_ns = False
        self.external_forces_manual = False

    def setName(self, string name):
        """Sets name of cable, used for csv file

        Parameters
        ----------
        name: str
            Name of cable.
        """
        self.name = name

    def _recordValues(self):
        """Records values in csv files
        """
        if self.ProtChSystem.model is not None:
            t_last = self.ProtChSystem.model.stepController.t_model_last
            try:
                dt_last = self.ProtChSystem.model.levelModelList[-1].dt_last
            except:
                dt_last = 0
            t = t_last-dt_last
        else:
            t = self.ProtChSystem.thisptr.system.GetChTime()
        t_prot = Profiling.time()-Profiling.startTime
        self.record_file = os.path.join(Profiling.logDir, self.name+'_pos.csv')
        if t == 0:
            headers = ['t', 't_prot']
            for i in range(self.thisptr.nodes.size()):
                headers += ['x'+str(i), 'y'+str(i), 'z'+str(i)]
            with open(self.record_file, 'w') as csvfile:
                writer = csv.writer(csvfile, delimiter=',')
                writer.writerow(headers)
        row = [t, t_prot]
        positions = self.getNodesPosition()
        for pos in positions:
            row += [pos[0], pos[1], pos[2]]
        with open(self.record_file, 'a') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerow(row)
        self.record_file = os.path.join(Profiling.logDir, self.name+'_T.csv')
        if t == 0:
            headers = ['t', 't_prot']
            # for i in range(self.thisptr.nodes.size()):
            #     headers += ['Tbx'+str(i), 'Tby'+str(i), 'Tbz'+str(i), 'Tfx'+str(i), 'Tfy'+str(i), 'Tfz'+str(i)]
            headers += ['Tb0', 'Tb1', 'Tb2', 'Tf0', 'Tf1', 'Tf2']
            with open(self.record_file, 'w') as csvfile:
                writer = csv.writer(csvfile, delimiter=',')
                writer.writerow(headers)
        row = [t, t_prot]
        Tb = self.getTensionBack()
        Tf = self.getTensionFront()
        row += [Tb[0], Tb[1], Tb[2], Tf[0], Tf[1], Tf[2]]
        with open(self.record_file, 'a') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerow(row)

    def getTensionBack(self):
        """
        """
        cdef ch.ChVector T
        if self.thisptr.constraint_back:
            T = deref(self.thisptr.constraint_back).GetReactionOnNode()
            return pych.ChVector_to_npArray(T)
        else:
            return np.zeros(3)

    def getTensionFront(self):
        """
        """
        cdef ch.ChVector T
        if self.thisptr.constraint_front:
            T = deref(self.thisptr.constraint_front).GetReactionOnNode()
            return pych.ChVector_to_npArray(T)
        else:
            return np.zeros(3)

    def calculate_init(self):
        # build position vector of nodes (for each segment)
        # self.setNodesPosition()
        # build cable (nodes, elements, etc)
        self.thisptr.buildCable()
        nb_nodes = self.thisptr.nodes.size()
        if self.fluid_velocity_array is None:
            self.fluid_velocity_array = np.zeros((nb_nodes, 3))
            self.fluid_velocity_array_previous = np.zeros((nb_nodes, 3))
        if self.fluid_acceleration_array is None:
            self.fluid_acceleration_array = np.zeros((nb_nodes, 3))
        if self.fluid_density_array is None:
            self.fluid_density_array = np.zeros(nb_nodes)


    def prestep(self):
        """Sets external forces on the cable (if any)
        """
        if self.ProtChSystem.model is not None and self.external_forces_from_ns is True:
            Profiling.logEvent('moorings extern forces')
            self.setExternalForces()
        elif self.external_forces_manual is not None:
            self.setExternalForces()

    def poststep(self):
        """Records values
        """
        comm = Comm.get().comm.tompi4py()
        if comm.rank == self.ProtChSystem.chrono_processor:
            self._recordValues()

    def setNodesPositionFunction(self, function):
        """Function to build nodes

        Parameters
        ----------
        function:
            Must be a function taking one argument (e.g. distance
            along cable) and returning 3 arguments (x, y, z) coords.
        """
        self.nodes_function = function

    def fixFrontNode(self, bool fixed):
        """Fix front node of cable

        Parameters
        ----------
        fixed: bool
            Fixes node if True
        """
        assert self.nodes_built is True, 'call buildNodes() before calling this function'
        deref(self.thisptr.nodes.front()).SetFixed(fixed)

    def fixBackNode(self, bool fixed):
        """Fix back node of cable

        Parameters
        ----------
        fixed: bool
            Fixes node if True
        """
        assert self.nodes_built is True, 'call buildNodes() before calling this function'
        deref(self.thisptr.nodes.back()).SetFixed(fixed)

    def attachBackNodeToBody(self, ProtChBody body):
        """Attaches back node to a body with ChLinkLockLock

        Parameters
        ----------
        body: ProtChBody
            body to which the node will be attached
        """
        assert self.nodes_built is True, 'call buildNodes() before calling this function'
        self.thisptr.attachBackNodeToBody(body.thisptr.body)

    def attachFrontNodeToBody(self, ProtChBody body):
        """Attaches front node to a body with ChLinkLockLock

        Parameters
        ----------
        body: ProtChBody
            body to which the node will be attached
        """
        assert self.nodes_built is True, 'call buildNodes() before calling this function'
        self.thisptr.attachFrontNodeToBody(body.thisptr.body)

    def getElementMass(self, int i=0):
        mass = deref(self.thisptr.elems[i]).GetMass()
        return mass

    def getTensionElement(self, int i=0):
        cdef ch.ChVector[double] F
        F = self.thisptr.getTensionElement(i)
        return np.array([F.x(), F.y(), F.z()])

    def getLengthElems(self):
        lengths = np.zeros(self.thisptr.elems.size())
        for i in range(self.thisptr.elems.size()):
            lengths[i] = deref(self.thisptr.elems[i]).GetLengthX()
        return lengths

    def setDragCoefficients(self, double tangential, double normal, int segment_nb):
        """Sets drag coefficients of cable

        Parameters
        ----------
        tangential: double
            Tangential drag coefficient.
        normal: double
            Normal drag coefficient.
        segment_nb: int
            Segment number to which these coefficients apply.
        """
        deref(self.thisptr.cables[segment_nb]).setDragCoefficients(tangential, normal)

    def setAddedMassCoefficients(self, double tangential, double normal, int segment_nb):
        """Sets added mass coefficients of cable

        Parameters
        ----------
        tangential: double
            Tangential added mass coefficient.
        normal: double
            Normal added mass coefficient.
        segment_nb: int
            Segment number to which these coefficients apply.
        """
        deref(self.thisptr.cables[segment_nb]).setAddedMassCoefficients(tangential, normal)

    def setNodesPosition(self):
        """Builds the nodes of the cable.

        (!) Must be called after setNodesPositionFunction()
        """
        cdef ch.ChVector[double] vec
        for i in range(self.thisptr.cables.size()):
            deref(self.thisptr.cables[i]).mvecs.clear()
            L0 = deref(self.thisptr.cables[i]).L0
            L = deref(self.thisptr.cables[i]).length
            nb_elems = deref(self.thisptr.cables[i]).nb_elems
            if self.beam_type == "BeamANCF":
                nb_nodes = nb_elems*2+1
            elif self.beam_type == "CableANCF":
                nb_nodes = nb_elems+1
            else:
                print("set element type")
                sys.exit()
            ds = L/(nb_nodes-1)
            for j in range(nb_nodes):
                x, y, z = self.nodes_function(L0+ds*j)
                vec = ch.ChVector[double](x, y, z)
                deref(self.thisptr.cables[i]).mvecs.push_back(vec)
        self.buildNodes()

    def buildNodes(self):
        self.thisptr.buildNodes()
        self.nodes_built = True

    def getNodesPosition(self):
        """Gives array of nodes position

        Returns
        -------
        pos: np.ndarray
            Array of nodes position.
        """
        pos = np.zeros(( self.thisptr.nodes.size(),3 ))
        for i in range(self.thisptr.nodes.size()):
            vec = deref(self.thisptr.nodes[i]).GetPos()
            pos[i] = [vec.x(), vec.y(), vec.z()]
        return pos

    def getNodesD(self):
        """Gives direction of nodes

        Returns
        -------
        dire: np.ndarray
            Array of nodes direction.
        """
        dire = np.zeros(( self.thisptr.nodes.size(),3 ))
        for i in range(self.thisptr.nodes.size()):
            vec = deref(self.thisptr.nodes[i]).GetD()
            dire[i] = [vec.x(), vec.y(), vec.z()]
        return dire

    def getNodesDD(self):
        """(!) Only for BeamANCF
        """
        pos = np.zeros(( self.thisptr.nodes.size(),3 ))
        for i in range(self.thisptr.nodes.size()):
            vec = deref(self.thisptr.nodes[i]).GetDD()
            pos[i] = [vec.x(), vec.y(), vec.z()]
        return pos

    def setContactMaterial(self, pych.ChMaterialSurfaceSMC mat):
        """Sets contact material of the cable

        Parameters
        ----------
        mat: ChMaterialSurfaceSMC
            Material of cable.
        """
        self.thisptr.setContactMaterial(mat.sharedptr)

    def setExternalForces(self, fluid_velocity_array=None, fluid_density_array=None,
                          fluid_acceleration_array=None):
        """
        Sets external forces acting on cables
        Pass fluid velocity_array as argument only for debugging (must be an array as long as the number of nodes)
        """
        # get velocity at nodes
        # cdef np.ndarray fluid_velocity = np.zeros((len(self.thisptr.nodes.size()), 3))
        self.fluid_velocity_array_previous[:] = self.fluid_velocity_array
        if fluid_velocity_array is not None:
            self.fluid_velocity_array = fluid_velocity_array
        if fluid_density_array is not None:
            self.fluid_density_array = fluid_density_array
        if fluid_acceleration_array is not None:
            self.fluid_acceleration_array = fluid_acceleration_array
        cdef vector[ch.ChVector[double]] fluid_velocity
        cdef vector[ch.ChVector[double]] fluid_acceleration
        cdef ch.ChVector[double] vel
        cdef ch.ChVector[double] acc
        cdef vector[double] fluid_density
        cdef double dens
        for i in range(self.thisptr.nodes.size()):
            if self.ProtChSystem.model is not None and self.external_forces_from_ns is True:
                vec = deref(self.thisptr.nodes[i]).GetPos()
                x = vec.x()
                y = vec.y()
                z = vec.z()
                coords = np.array([x, y, z])
                vel_arr = np.zeros(3)
                vel_grad_arr = np.zeros(3)
                xi, el, rank = self.ProtChSystem.findElementContainingCoords(coords[:self.nd])
                vel_arr[:] = self.ProtChSystem.getFluidVelocityLocalCoords(xi, el, rank)
                #vel_grad_arr[:] = self.ProtChSystem.getFluidVelocityGradientLocalCoords(xi, el, rank)
                # acc = du/dt+u.grad(u)
                #acc_arr = (vel_arr-fluid_velocity_array_previous[i])/dt+vel_arr*vel_grad_arr
                #arr[:self.nd] = self.ProtChSystem.findFluidVelocityAtCoords(coords[:self.nd])
                self.fluid_velocity_array[i] = vel_arr
                vel = ch.ChVector[double](vel_arr[0], vel_arr[1], vel_arr[2])
                fluid_velocity.push_back(vel)
            else:
                vel = ch.ChVector[double](self.fluid_velocity_array[i][0], self.fluid_velocity_array[i][1], self.fluid_velocity_array[i][2])
                fluid_velocity.push_back(vel)
                acc = ch.ChVector[double](self.fluid_acceleration_array[i][0], self.fluid_acceleration_array[i][1], self.fluid_acceleration_array[i][2])
                fluid_acceleration.push_back(acc)
                dens = self.fluid_density_array[i]
                fluid_density.push_back(dens)
        self.thisptr.setFluidAccelerationAtNodes(fluid_acceleration)
        self.thisptr.setFluidVelocityAtNodes(fluid_velocity)
        self.thisptr.setFluidDensityAtNodes(fluid_density)
            # update drag forces
        self.thisptr.updateDragForces()
        self.thisptr.updateAddedMassForces()
        self.thisptr.applyForces()
        # update buoyancy forces
        # self.thisptr.updateBuoyancyForces()
        # update added mass forces
        # self.thisptr.updateAddedMassForces()

    def setFluidDensityAtNodes(self, np.ndarray density_array):
        cdef vector[double] fluid_density
        self.fluid_density_array = density_array
        cdef double dens
        for d in density_array:
            fluid_density.push_back(d)
        self.thisptr.setFluidDensityAtNodes(fluid_density)

    def setFluidVelocityAtNodes(self, np.ndarray velocity_array):
        cdef vector[ch.ChVector[double]] fluid_velocity
        cdef ch.ChVector[double] vel
        self.fluid_velocity_array = velocity_array
        for v in velocity_array:
            vel = ch.ChVector[double](v[0], v[1], v[2])
            fluid_velocity.push_back(vel)
        self.thisptr.setFluidVelocityAtNodes(fluid_velocity)

    def setFluidAccelerationAtNodes(self, np.ndarray acceleration_array):
        cdef vector[ch.ChVector[double]] fluid_acceleration
        cdef ch.ChVector[double] acc
        self.fluid_acceleration_array = acceleration_array
        for a in acceleration_array:
            acc = ch.ChVector[double](a[0], a[1], a[2])
            fluid_acceleration.push_back(acc)
        self.thisptr.setFluidAccelerationAtNodes(fluid_acceleration)


def getLocalNearestNode(coords, kdtree):
    """Finds nearest node to coordinates (local)
    Parameters
    ----------
    coords: array_like
        coordinates from which to find nearest node
    kdtree: scipy.spatial.cKDTree
        instance of scipy kdtree

    Returns
    -------
    node: int
        nearest node index
    distance: float
        distance to nearest node
    """
    # determine local nearest node distance
    distance, node = kdtree.query(coords)
    return node, distance

def getLocalElement(femSpace, coords, node):
    """Given coordinates and its nearest node, determine if it is on a
    local element.

    Parameters
    ----------
    femSpace: object
        finite element space
    coords: array_like
        coordinates from which to element
    node: int
        nearest node index

    Returns
    -------
    eN: int or None
        local index of element (None if not found)
    """
    patchBoundaryNodes=set()
    checkedElements=[]
    # nodeElementOffsets give the indices to get the elements sharing the node
    for eOffset in range(femSpace.mesh.nodeElementOffsets[node], femSpace.mesh.nodeElementOffsets[node + 1]):
        eN = femSpace.mesh.nodeElementsArray[eOffset]
        checkedElements.append(eN)
        # union of set
        patchBoundaryNodes|=set(femSpace.mesh.elementNodesArray[eN])
        # evaluate the inverse map for element eN (global to local)
        xi = femSpace.elementMaps.getInverseValue(eN, coords)
        #J = femSpace.elementMaps.getJacobianValues(eN, )
        # query whether xi lies within the reference element
        if femSpace.elementMaps.referenceElement.onElement(xi):
            return eN
    # extra loop if case coords is in neighbour element
    for node in patchBoundaryNodes:
        for eOffset in range(femSpace.mesh.nodeElementOffsets[node], femSpace.mesh.nodeElementOffsets[node + 1]):
            eN = femSpace.mesh.nodeElementsArray[eOffset]
            if eN not in checkedElements:
                checkedElements.append(eN)
                # evaluate the inverse map for element eN
                xi = femSpace.elementMaps.getInverseValue(eN, coords)
                # query whether xi lies within the reference element
                if femSpace.elementMaps.referenceElement.onElement(xi):
                    return eN
    # no elements found
    return None


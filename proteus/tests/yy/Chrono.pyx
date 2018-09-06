import numpy
cimport numpy
from proteus import AuxiliaryVariables, Archiver
from proteus.Profiling import  logEvent
from cython.view cimport array as cvarray
# import ode
import numpy as np
# Import the C-level symbols of numpy
cimport numpy as np
from libcpp cimport bool

cdef extern from "Chrono.h":
    cdef cppclass cppMBDModel:
        void step(double* forces, double* torques, double dt)
        const int num_particles
        double *diam_
        double *pos_
        double *vel_
        double *angular_vel_
        void writeThisFrame()

    cppMBDModel* newMBDModel(double m_timeStep,
                             double* m_container_center,
                             double* m_container_dims,
                             double m_particles_density,
                             double m_particles_diameter,
                             double* m_gravity,
                             int nParticle,
                             double* ball_center,
                             double* ball_radius)

cdef class MBDModel:
    cdef cppMBDModel* thisptr
    cdef object model
    cdef double x,y,z

#     cdef object np_pos
#     cdef object np_normal
#     cdef object n_coor
#     cdef object numnodes

    def __cinit__(self,
                  double m_timeStep,
                  numpy.ndarray m_container_center,
                  numpy.ndarray m_container_dims,
                  double m_particles_density,
                  double m_particles_diameter,
                  numpy.ndarray m_gravity,
                  int nParticle,
                  numpy.ndarray ball_center,
                  numpy.ndarray ball_radius,
                  ):
       
        self.thisptr =  newMBDModel(m_timeStep,
                                    <double*> m_container_center.data,
                                    <double*> m_container_dims.data,
                                    m_particles_density,
                                    m_particles_diameter,
                                    <double*> m_gravity.data,
                                    nParticle,
                                    <double*> ball_center.data,
                                    <double*> ball_radius.data
                                     )

#         self.n_coor=np.zeros((self.thisptr.num_particles,3))
        #self.calcNodalInfo()

    def getNumParticles(self):
        return self.thisptr.num_particles

    def attachModel(self,model,ar):
        self.model=model
        return self

    def get_u(self):
        return 0
    def get_v(self):
        return 0
    def get_w(self):
        return 0
        
    def calculate_init(self):
        self.calculate()

    def SyncData(self,numpy.ndarray new_coors):
        #self.thisptr.prepareData(<double*> new_coors.data)
        return

#     def d_N_IBM(self, numpy.ndarray x):
#         import numpy as np
#         cdef double dN[4]
#         self.thisptr.calc_d_N_IBM(<double*> x.data,dN)
#         dist=dN[0]
#         normal=np.array([dN[1],dN[2],dN[3]])
#         return dist, normal

    def getParticlesDiameter(self, i):
        return self.thisptr.diam_[i]

    def get_Angular_Vel(self, i):
        vel=np.array([ self.thisptr.angular_vel_[3*i+0],self.thisptr.angular_vel_[3*i+1], self.thisptr.angular_vel_[3*i+2] ])
        return vel

    def get_Vel(self, i):
        vel=np.array([ self.thisptr.vel_[3*i+0],self.thisptr.vel_[3*i+1], self.thisptr.vel_[3*i+2] ])
        return vel

    def get_Pos(self, i):
        pos=np.array([ self.thisptr.pos_[3*i+0],self.thisptr.pos_[3*i+1], self.thisptr.pos_[3*i+2] ])
        return pos

    def step(self,
             numpy.ndarray force,
             numpy.ndarray torque,
             double dt):
        self.thisptr.step(<double*> force.data,
                          <double*> torque.data,
                          dt)

# 
#     def get_par_pos(self):
#         import numpy as np
#         np_pos=np.zeros((self.thisptr.num_particles,3))
#         for i in range(self.thisptr.num_particles):
#             self.thisptr.get_pos(i,np_pos[i,0],np_pos[i,1],np_pos[i,1])
#         return np_pos
# 
# 
#     def get_par_vel(self):
#         import numpy as np
#         np_vel=np.zeros((self.thisptr.num_particles,3))
#         for i in range(self.thisptr.num_particles):
#             self.thisptr.get_vel(i,np_vel[i,0],np_vel[i,1],np_vel[i,1])
# 
#         return np_vel

    def calculate(self, 
                  numpy.ndarray Forces,
                  numpy.ndarray torques,
                  double dt):
        logEvent("Calling chrono with dt " +`dt`)
        self.step(Forces,torques,dt)

    def writeFrame(self):
        self.thisptr.writeThisFrame()

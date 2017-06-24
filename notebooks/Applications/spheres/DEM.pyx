import numpy
cimport numpy
from proteus import AuxiliaryVariables, Archiver
from proteus.Profiling import  logEvent

cdef extern from "DEM.h":
    cdef cppclass cppDEM:
        int get_nSolids()
        void step(double* force, double* torque, double dt)
        double hx(double* x, double dt)
        double hy(double* x, double dt)
        double hz(double* x, double dt)
        void set_center_array(double* center)
        void set_radius_array(double* radius)
    cppDEM* newDEM(const char* config_in)

cdef class DEM:
    cdef cppDEM* thisptr
    cdef object model
    cdef numpy.ndarray center_array
    cdef numpy.ndarray radius_array
    def __cinit__(self,
                  config_in):
        self.thisptr = newDEM(config_in)
        self.center_array=numpy.zeros((self.thisptr.get_nSolids(),3),'d')
        self.radius_array=numpy.zeros((self.thisptr.get_nSolids()),'d')
    def attachModel(self,model,ar):
        self.model=model
        return self
    def get_radius_array(self):
        return self.radius_array
    def get_center_array(self):
        return self.center_array
    def get_u(self):
        return self.last_velocity[0]
    def get_v(self):
        return self.last_velocity[1]
    def get_w(self):
        return self.last_velocity[2]
    def calculate_init(self):
        self.calculate()
    def hx(self, numpy.ndarray x, double t):
        return self.thisptr.hx(<double*> x.data,t)
    def hy(self, numpy.ndarray x, double t):
        return self.thisptr.hy(<double*> x.data,t)
    def hz(self, numpy.ndarray x, double t):
        return self.thisptr.hz(<double*> x.data,t)
    def step(self,
             numpy.ndarray force,
             numpy.ndarray torque,
             double dt):
        self.thisptr.step(<double*> force.data,
                          <double*> torque.data,
                          dt)
        self.thisptr.set_center_array(<double*> self.center_array.data)
        self.thisptr.set_radius_array(<double*> self.radius_array.data)
    def calculate(self):
        from numpy.linalg import inv
        import copy
        try:
            dt = self.model.levelModelList[-1].dt_last
        except:
            dt = 1.0e-8
        #cdef numpy.ndarray F;
        F = self.model.levelModelList[-1].coefficients.netForces_p[7,:] + self.model.levelModelList[-1].coefficients.netForces_v[7,:];
        #cdef numpy.ndarray M;
        M = self.model.levelModelList[-1].coefficients.netMoments[7,:]
        logEvent("x Force " +`self.model.stepController.t_model_last`+" "+`F[0]`)
        logEvent("y Force " +`self.model.stepController.t_model_last`+" "+`F[1]`)
        logEvent("z Force " +`self.model.stepController.t_model_last`+" "+`F[2]`)
        logEvent("x Moment " +`self.model.stepController.t_model_last`+" "+`M[0]`)
        logEvent("y Moment " +`self.model.stepController.t_model_last`+" "+`M[1]`)
        logEvent("z Moment " +`self.model.stepController.t_model_last`+" "+`M[2]`)
        logEvent("dt " +`dt`)
        scriptMotion=False
        linearOnly=False
        self.step( F,
                   M,
                   dt)

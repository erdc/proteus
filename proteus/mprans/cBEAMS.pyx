import numpy
cimport numpy
from proteus import *
from proteus.Transport import *
from proteus.Transport import OneLevelTransport

cdef extern from "BEAMS.h" namespace "proteus":
    cdef cppclass BEAMS_base:
        void calculateBeams(int nElements_global,
				   double rho_0,
				   double rho_1,
				   double* phi,
				   double* q_x,
				   double* q_velocity,
				   double* q_dV,
				   double* q_dragBeam1,
				   double* q_dragBeam2,
				   double* q_dragBeam3,
				   int nBeams,
				int nBeamElements,
				int beam_quadOrder,
				double beam_Cd,
				double* beamRadius,
				double* xq,
				double* yq,
				double* zq,
				double* Beam_h,
				double* dV_beam,
				double* q1,
				double* q2,
				double* q3,
				double* vel_avg,
				double* netBeamDrag,
				int* beamIsLocal)       
     
    BEAMS_base* newBEAMS(int nSpaceIn,
                           int nQuadraturePoints_elementIn,
                           int nDOF_mesh_trial_elementIn,
                           int nDOF_trial_elementIn,
                           int nDOF_test_elementIn,
                           int nQuadraturePoints_elementBoundaryIn,
                           int CompKernelFlag)

cdef class cBEAMS_base:
   cdef BEAMS_base* thisptr
   def __cinit__(self,
                 int nSpaceIn,
                 int nQuadraturePoints_elementIn,
                 int nDOF_mesh_trial_elementIn,
                 int nDOF_trial_elementIn,
                 int nDOF_test_elementIn,
                 int nQuadraturePoints_elementBoundaryIn,
                 int CompKernelFlag):
       self.thisptr = newBEAMS(nSpaceIn,
                              nQuadraturePoints_elementIn,
                              nDOF_mesh_trial_elementIn,
                              nDOF_trial_elementIn,
                              nDOF_test_elementIn,
                              nQuadraturePoints_elementBoundaryIn,
                              CompKernelFlag)
   def __dealloc__(self):
       del self.thisptr
   		 


   def calculateBeams(self,
                         int nElements_global,
                         double rho_0,
                         double rho_1,
                         numpy.ndarray phi,
                         numpy.ndarray q_x,
                         numpy.ndarray q_velocity,
			 numpy.ndarray q_dV,
			 numpy.ndarray q_dragBeam1,
			numpy.ndarray q_dragBeam2,
				   numpy.ndarray q_dragBeam3,
				   int nBeams,
			 int nBeamElements,
			 int beam_quadOrder,
			 double beam_Cd,
			 numpy.ndarray beamRadius,
			 numpy.ndarray xq,
			 numpy.ndarray yq,
			 numpy.ndarray zq,
			 numpy.ndarray Beam_h,
			 numpy.ndarray dV_beam,
			 numpy.ndarray q1,
			 numpy.ndarray q2,
			 numpy.ndarray q3,
			 numpy.ndarray vel_avg,
			 numpy.ndarray netBeamDrag,
			 numpy.ndarray beamIsLocal):			 
       self.thisptr.calculateBeams(nElements_global,
                                       rho_0,
                                       rho_1,
                                       <double*> phi.data,
                                       <double*> q_x.data,
                                       <double*> q_velocity.data,
				       <double*> q_dV.data,
				       <double*> q_dragBeam1.data,
				   <double*> q_dragBeam2.data,
				   <double*> q_dragBeam3.data,
				   nBeams,
				       nBeamElements,
				       beam_quadOrder,
				       beam_Cd,
				       <double*> beamRadius.data,
				       <double*> xq.data,
				       <double*> yq.data,
				       <double*> zq.data,
				       <double*> Beam_h.data,
				       <double*> dV_beam.data,
				       <double*> q1.data,
				       <double*> q2.data,
				       <double*> q3.data,
				       <double*> vel_avg.data,
				       <double*> netBeamDrag.data,
				       <int*> beamIsLocal.data)
			       
   
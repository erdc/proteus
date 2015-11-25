import numpy
cimport numpy
from cpython cimport array
#pull in C++ definitions and declare interface
cdef extern from "mprans/SedClosure.h" namespace "proteus":
    cdef cppclass cppHsuSedStress:
        cppHsuSedStress(
            double parameterIn,  # Dummy parameter for first function
		 double aDarcy, # darcy parameter for drag term. Default value from Ergun (1952) is 150
		 double betaForch, # forchheimer parameter for drag term. Default value from Ergun (1952) is 1.75
		 double grain, # Grain size, default assumed as d50
		 double packFraction #Critical volume fraction for switching the drag relation 0.2 by default, see Chen and Hsu 2014
		 )
        double dummy(double porosity,double pf)
        double granularDrag(
                            double sedF, # Sediment fraction
                            double* uFluid, #Fluid velocity
                            double* uSolid, #Sediment velocity
                            double nu #Kinematic viscosity
                           )
#define the way we want to present to Python
cdef class HsuSedStress:
    cdef  cppHsuSedStress* thisptr
    def __cinit__(self, parameterIn, aDarcy, betaForch, grain, packFactor ):
        self.thisptr = new cppHsuSedStress( parameterIn,aDarcy, betaForch, grain, packFactor)
    def __dealloc__(self):
        del self.thisptr
    def dummy(self, porosity, pf):
        return self.thisptr.dummy(porosity,pf)
    def granularDrag(self, 
                     sedF, 
                     numpy.ndarray uFluid, 
                     numpy.ndarray uSolid, 
                     nu):
        return self.thisptr.granularDrag(sedF, 
                                  < double * > uFluid.data,
                                  < double * > uSolid.data, 
                                  nu)

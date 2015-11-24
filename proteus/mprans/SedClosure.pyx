import numpy
cimport numpy

#pull in C++ definitions and declare interface
cdef extern from "mprans/SedClosure.h" namespace "proteus":
    cdef cppclass cppHsuSedStress:
        cppHsuSedStress(double parameterIn)
        double M_sf_x(double porosity,double pf)

#define the way we want to present to Python
cdef class HsuSedStress:
    cdef  cppHsuSedStress* thisptr
    def __cinit__(self, parameterIn=3.0):
        self.thisptr = new cppHsuSedStress(parameterIn)
    def __dealloc__(self):
        del self.thisptr
    def M_sf_x(self, porosity, pf):
        return self.thisptr.M_sf_x(porosity,pf)

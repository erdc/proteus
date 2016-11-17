#!python
#cython: embedsignature=True
#cython: profile=True

cimport cython
from libc.math cimport tanh, sqrt, exp, log, sin, cos, cosh, sinh, M_PI #from math import pi, tanh, sqrt, exp, log, sin, cos, cosh, sinh
cimport numpy as np


cdef extern from "WaveTools.h" namespace "proteus":
    cdef const int nDim
    cdef double __cpp_eta_mode(double* x, double t, double* kDir, double omega, double phi, double amplitude)
    cdef double* __cpp_vel_mode(double* x, double t, double* kDir, double kAbs, double omega, double phi, double amplitude, double mwl, double depth, double* waveDir, double* vDir, double sinhL)
    cdef double __cpp_etaFenton(double* x, double t, double* kDir, double kAbs, double omega, double phi0, double amplitude, int Nf, double* Ycoeff)
    cdef double* __cpp_uFenton(double* x, double t, double* kDir, double kAbs, double omega, double phi0, double amplitude, double mwl, double depth,double gAbs, int Nf, double* Bcoeff, double* mV, double* waveDir, double* vDir, double* sinhF, double* tanhF)


# pointer to eta function
ctypedef double (*cfeta) (MonochromaticWaves, double* , double )  

# pointer to velocity function
ctypedef double* (*cfvel) (MonochromaticWaves, double* , double )

ctypedef double[:] double_memview1



cdef class  MonochromaticWaves:
    cdef np.ndarray g
    cdef np.ndarray waveDir
    cdef np.ndarray vDir
    cdef double gAbs
    cdef double phi
    cdef double depth
    cdef double omega
    cdef double k
    cdef double phi0
    cdef double sinhL
    cdef int Nf
    cdef np.ndarray Ycoeff
    cdef np.ndarray Bcoeff
    cdef np.ndarray kDir
    cdef np.ndarray sinhF
    cdef np.ndarray tanhF
    cdef double amplitude
    cdef np.ndarray mV    
    cdef double* kDir_
    cdef double* waveDir_
    cdef double* vDir_
    cdef double* mV_    
    cdef double* Ycoeff_
    cdef double* Bcoeff_
    cdef double* sinhF_
    cdef double* tanhF_
    cdef double[3] kDir_c
    cdef double[3] waveDir_c
    cdef double[3] vDir_c
    cdef double[3] mV_c    
    cdef double[1000] Ycoeff_c
    cdef double[1000] Bcoeff_c
    cdef double[1000] sinh_c
    cdef double[1000] tanh_c
    cdef public:
        double wavelength
        double mwl
    cdef  cfeta _cpp_eta
    cdef  cfvel _cpp_u
    cdef object waveType
    cdef double etaLinear(self, double* x, double t)
    cdef double etaFenton(self, double* x, double t)
    cdef double* uLinear(self, double* x, double t)
    cdef double* uFenton(self, double* x, double t)
'''
cdef class RandomWaves:
    cdef np.ndarray g
    cdef np.ndarray waveDir
    cdef np.ndarray vDir
    cdef double gAbs
    cdef double phi
    cdef double depth
    cdef double omega
    cdef double k
    cdef double phi0
    cdef double sinhL
    cdef int Nf
    cdef np.ndarray Ycoeff
    cdef np.ndarray Bcoeff
    cdef np.ndarray kDir
    cdef np.ndarray sinhF
    cdef np.ndarray tanhF
    cdef double amplitude
    cdef np.ndarray mV    
    cdef double* kDir_
    cdef double* waveDir_
    cdef double* vDir_
    cdef double* mV_    
    cdef double* Ycoeff_
    cdef double* Bcoeff_
    cdef double* sinhF_
    cdef double* tanhF_
    cdef double[3] kDir_c
    cdef double[3] waveDir_c
    cdef double[3] vDir_c
    cdef double[3] mV_c    
    cdef double[1000] Ycoeff_c
    cdef double[1000] Bcoeff_c
    cdef double[1000] sinh_c
    cdef double[1000] tanh_c
    cdef public:
        double wavelength
        double mwl
    cdef  cfeta _cpp_eta
    cdef  cfvel _cpp_u
    cdef object waveType
    cdef double etaLinear(self, double* x, double t)
    cdef double etaFenton(self, double* x, double t)
    cdef double* uLinear(self, double* x, double t)
    cdef double* uFenton(self, double* x, double t)
'''

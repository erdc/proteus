#!python
#cython: embedsignature=True
#cython: profile=True

cimport cython
from libc.math cimport tanh, sqrt, exp, log, sin, cos, cosh, sinh, M_PI #from math import pi, tanh, sqrt, exp, log, sin, cos, cosh, sinh
cimport numpy as np
from libcpp cimport bool

cdef extern from "WaveTools.h" namespace "proteus":
    cdef const int nDim
    cdef double fastcosh(double k, double Z,  bool cosh)
    cdef double fastcos(double phase)
    cdef double __cpp_eta_mode(double* x, double t, double* kDir, double omega, double phi, double amplitude)
    cdef double* __cpp_vel_mode(double* x, double t, double* kDir, double kAbs, double omega, double phi, double amplitude, double mwl, double depth, double* waveDir, double* vDir, double tanhkd)
    cdef double __cpp_etaFenton(double* x, double t, double* kDir, double kAbs, double omega, double phi0, double amplitude, int Nf, double* Ycoeff)
    cdef double* __cpp_uFenton(double* x, double t, double* kDir, double kAbs, double omega, double phi0, double amplitude, double mwl, double depth,double gAbs, int Nf, double* Bcoeff, double* mV, double* waveDir, double* vDir, double* tanhF)
    cdef double __cpp_etaRandom(double* x, double t, double* kDir, double* omega, double* phi, double* amplitude, int N)
    cdef double* __cpp_uRandom(double* x,double t,double* kDir,double* kAbs, double* omega, double* phi, double* amplitude, double mwl, double depth, int N, double* waveDir, double* vDir, double* tanhKd )
    cdef double* __cpp_uDir(double* x,double t,double* kDir,double* kAbs, double* omega, double* phi, double* amplitude, double mwl, double depth, int N, double* waveDir, double* vDir, double* tanhKd )
    cdef int __cpp_findWindow(double t, double handover, double t0, double Twindow, int Nwindows, double* windows_handover)
    cdef double __cpp_etaDirect(double* x, double* x0, double t, double* kDir, double* omega, double* phi, double* amplitude, int N)
    cdef double* __cpp_uDirect(double* x,double* x0,double t,double* kDir,double* kAbs, double* omega, double* phi, double* amplitude, double mwl, double depth, int N, double* waveDir, double* vDir, double* tanhKd )
    cdef double __cpp_etaWindow(double* x, double* x0, double t, double* t0, double* kDir, double* omega, double* phi, double* amplitude, int N, int Nw)
    cdef  double* __cpp_uWindow(double* x, double* x0, double t, double* T0, double* kDir, double* kAbs, double* omega, double* phi, double* amplitude, double mwl, double depth, int N,int Nw, double* waveDir, double* vDir, double* tanhKd )

# pointer to eta function
ctypedef double (*cfeta) (MonochromaticWaves, double* , double )  

# pointer to velocity function
ctypedef double* (*cfvel) (MonochromaticWaves, double* , double )



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
    cdef double tanhL
    cdef int Nf
    cdef np.ndarray Ycoeff
    cdef np.ndarray Bcoeff
    cdef np.ndarray kDir
    cdef np.ndarray tanhF
    cdef double amplitude
    cdef np.ndarray mV    
    cdef double* kDir_
    cdef double* waveDir_
    cdef double* vDir_
    cdef double* mV_    
    cdef double* Ycoeff_
    cdef double* Bcoeff_
    cdef double* tanhF_
    cdef double[3] kDir_c
    cdef double[3] waveDir_c
    cdef double[3] vDir_c
    cdef double[3] mV_c    
    cdef double[1000] Ycoeff_c
    cdef double[1000] Bcoeff_c
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

cdef class RandomWaves:
    cdef np.ndarray g
    cdef np.ndarray waveDir
    cdef np.ndarray vDir
    cdef double gAbs
    cdef double Hs
    cdef double Tp
    cdef double fp
    cdef double bandFactor
    cdef object phi
    cdef double depth
    cdef double df
    cdef int N
    cdef np.ndarray kDir
    cdef np.ndarray ai
    cdef double* waveDir_
    cdef double* vDir_
    cdef double* tanh_
    cdef double* ai_
    cdef double* phi_
    cdef double* omega_
    cdef double* kDir_
    cdef double* ki_
    cdef double[3] waveDir_c
    cdef double[3] vDir_c
    cdef double[30000] kDir_c
    cdef double[10000] omega_c
    cdef double[10000] ki_c
    cdef double[10000] ai_c
    cdef double[10000] tanh_c
    cdef double[10000] phi_c
    cdef public:
        double mwl
        np.ndarray fi
        np.ndarray fim
        np.ndarray Si_Jm
        np.ndarray ki
        np.ndarray omega
        np.ndarray tanhF
    cdef double _cpp_eta(self, double* x, double t)
    cdef double* _cpp_u(self, double* x, double t)



cdef class MultiSpectraRandomWaves:
    cdef int Nall
    cdef int N
    cdef np.ndarray g
    cdef np.ndarray vDir
    cdef np.ndarray waveDir
    cdef np.ndarray omegaM
    cdef np.ndarray phiM
    cdef np.ndarray kiM
    cdef np.ndarray kDirM
    cdef np.ndarray tanhFM
    cdef np.ndarray aiM
    cdef double * vDir_
    cdef double* omegaM_
    cdef double* phiM_
    cdef double* kiM_
    cdef double* kDirM_
    cdef double* tanhM_
    cdef double* waveDirM_
    cdef double* aiM_
    cdef double[3] vDir_c
    cdef double[30000] kDir_cM
    cdef double[30000] waveDir_cM
    cdef double[10000] omega_cM
    cdef double[10000] ki_cM
    cdef double[10000] ai_cM
    cdef double[10000] tanh_cM
    cdef double[10000] phi_cM
    cdef public:
        double mwl
        double depth
    cdef double _cpp_eta(self, double* x, double t)
    cdef double* _cpp_u(self, double* x, double t)

cdef class DirectionalWaves:
    cdef int Nall
    cdef int Mtot
    cdef int N
    cdef np.ndarray vDir
    cdef np.ndarray omega
    cdef np.ndarray tanh
    cdef np.ndarray waveDir0
    cdef np.ndarray waveDirs
    cdef np.ndarray phiDirs
    cdef np.ndarray aiDirs
    cdef np.ndarray ki
    cdef np.ndarray kDirs
    cdef np.ndarray tanhF
    cdef double * vDir_
    cdef double* omega_
    cdef double* phi_
    cdef double* ki_
    cdef double* kDir_
    cdef double* tanh_
    cdef double* waveDir_
    cdef double* ai_
    cdef double[3] vDir_c
    cdef double[300000] kDir_c
    cdef double[300000] waveDir_c
    cdef double[100000] omega_c
    cdef double[100000] ki_c
    cdef double[100000] ai_c
    cdef double[100000] tanh_c
    cdef double[100000] phi_c
    cdef public:
        double mwl
        double depth
    cdef double _cpp_eta(self, double* x, double t)
    cdef double* _cpp_u(self, double* x, double t)


# pointer to eta function
ctypedef double (*cfeta2) (TimeSeries, double* , double )  

# pointer to velocity function
ctypedef double* (*cfvel2) (TimeSeries, double* , double )

cdef class  TimeSeries:
    cdef np.ndarray g
    cdef np.ndarray waveDir
    cdef np.ndarray vDir
    cdef double gAbs
    cdef double depth
    cdef int N
    cdef int Nall
    cdef np.ndarray x0
    cdef np.ndarray kDir
    cdef np.ndarray tanhF
    cdef np.ndarray time
    cdef np.ndarray etaS
    cdef np.ndarray ai
    cdef np.ndarray omega
    cdef np.ndarray phi
    cdef np.ndarray ki
    cdef int Nf
    cdef int Nwaves
    cdef double Tm
    cdef double overlap
    cdef double cutoff
    cdef double setup
    cdef double handover
    cdef double Twindow
    cdef double Tlag
    cdef double Toverlap
    cdef int Nwindows
    cdef list windows_handover
    cdef list windows_rec
    cdef list decompose_window
    cdef double dt
    cdef double t0
    cdef double tlength
    cdef np.ndarray mV    
    cdef double* kDir_
    cdef double* waveDir_
    cdef double* vDir_
    cdef double* tanh_
    cdef double* whand_
    cdef double* ai_
    cdef double* omega_
    cdef double* phi_
    cdef double* ki_
    cdef double* T0_
    cdef double[3] x0_
    cdef double[3000000] kDir_c
    cdef double[1000000] ki_c
    cdef double[1000000] ai_c
    cdef double[1000000] omega_c
    cdef double[1000000] phi_c
    cdef double[1000000] tanh_c
    cdef double[3] waveDir_c
    cdef double[3] vDir_c
    cdef double[3] x0_c
    cdef double[1000000] whand_c
    cdef double[1000000] T0
    cdef public:
        double wavelength
        double mwl
        object eta
        object u
    cdef  cfeta2 _cpp_eta
    cdef  cfvel2 _cpp_u
    cdef double _cpp_etaDirect(self, double* x, double t) 
    cdef double _cpp_etaWindow(self, double* x, double t) 
    cdef double* _cpp_uDirect(self, double* x, double t) 
    cdef double* _cpp_uWindow(self, double* x, double t) 

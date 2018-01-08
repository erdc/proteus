#!python
#cython: embedsignature=True
#cython: profile=True

cimport cython
from libc.math cimport tanh, sqrt, exp, log, sin, cos, cosh, sinh, M_PI #from math import pi, tanh, sqrt, exp, log, sin, cos, cosh, sinh
cimport numpy as np
from libcpp cimport bool

cdef extern from "WaveTools.h" namespace "proteus":
    cdef const int nDim
    cdef double fastcosh(double * hype,double k, double Z, bool fast)
    cdef double fastcos(double phase,bool fast)
    cdef double __cpp_eta_mode(double* x, double t, double* kDir, double omega, double phi, double amplitude,bool fast)
    cdef double __cpp_etaFenton(double* x, double t, double* kDir, double kAbs, double omega, double phi0, double amplitude, int Nf, double* Ycoeff,bool fast)
    cdef void __cpp_uFenton(double * U,double* x, double t, double* kDir, double kAbs, double omega, double phi0, double amplitude, double mwl, double depth,double gAbs, int Nf, double* Bcoeff, double* mV, double* waveDir, double* vDir, double* tanhF,bool fast)
    cdef double __cpp_etaRandom(double* x, double t, double* kDir, double* omega, double* phi, double* amplitude, int N,bool fast)
    cdef void __cpp_uRandom(double * U, double* x,double t,double* kDir,double* kAbs, double* omega, double* phi, double* amplitude, double mwl, double depth, int N, double* waveDir, double* vDir, double* tanhKd, double gAbs ,bool fast)
    cdef void __cpp_uDir(double* U, double* x,double t,double* kDir,double* kAbs, double* omega, double* phi, double* amplitude, double mwl, double depth, int N, double* waveDir, double* vDir, double* tanhKd , double gAbs,bool fast)
    cdef int __cpp_findWindow(double t, double handover, double t0, double Twindow, int Nwindows, double* windows_handover)
    cdef double __cpp_etaDirect(double* x, double* x0, double t, double* kDir, double* omega, double* phi, double* amplitude, int N,bool fast)
    cdef void __cpp_uDirect(double* U, double* x,double* x0,double t,double* kDir,double* kAbs, double* omega, double* phi, double* amplitude, double mwl, double depth, int N, double* waveDir, double* vDir, double* tanhKd, double gAbs ,bool fast)
    cdef double __cpp_etaWindow(double* x, double* x0, double t, double* t0, double* kDir, double* omega, double* phi, double* amplitude, int N, int Nw,bool fast)
    cdef  void __cpp_uWindow(double* U,double* x, double* x0, double t, double* T0, double* kDir, double* kAbs, double* omega, double* phi, double* amplitude, double mwl, double depth, int N,int Nw, double* waveDir, double* vDir, double* tanhKd , double gAbs,bool fast)
    cdef double __cpp_eta2nd(double* x, double t, double* kDir, double* ki, double* omega, double* phi, double* amplitude, int N, double* sinhKd, double* tanhKd,bool fast)
    cdef double __cpp_eta_short(double* x, double t, double* kDir, double* ki, double* omega, double* phi, double* amplitude, int N, double* sinhKd, double* tanhKd, double gAbs,bool fast)
    cdef double __cpp_eta_long(double* x, double t, double* kDir, double* ki, double* omega, double* phi, double* amplitude, int N, double* sinhKd, double* tanhKd, double gAbs,bool fast)
    cdef void __cpp_vel_mode_p(double* U, double * x, double t, double *kDir,double kAbs, double omega, double phi, double amplitude,double mwl, double depth, double *waveDir, double *vDir, double tanhkd, double gAbs,bool fast)

# pointer to eta function
ctypedef double (*cfeta) (MonochromaticWaves, double* , double )  

# pointer to velocity function
ctypedef void (*cfvel) (MonochromaticWaves, double*, double* , double )


cdef class  SolitaryWave:
    cdef double H,gAbs,c,depth,K,d2,d3
    cdef public:
        double mwl
    cdef np.ndarray g,waveDir,vDir,trans
    cdef bool fast

cdef class  MonochromaticWaves:
    cdef bool fast
    cdef np.ndarray g,waveDir,vDir,Ycoeff,Bcoeff,kDir,tanhF,mV
    cdef double gAbs,phi,depth,omega,k,phi0,tanhL,amplitude
    cdef int Nf					
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
    cdef void uLinear(self, double* U, double* x, double t)
    cdef void uFenton(self, double* U, double* x, double t)

cdef class RandomWaves:
    cdef bool fast
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
        double mwl,depth,gAbs,Tlag,Hs,Tp,fp,bandFactor,df
        int N
        np.ndarray fi,fim,Si_Jm,ki,omega,tanhF,g,waveDir,vDir,kDir,ai
        cdef object phi
    cdef double _cpp_eta(self , double* x, double t)
    cdef void _cpp_u(self, double *U, double* x, double t)



cdef class MultiSpectraRandomWaves:
    cdef bool fast
    cdef double gAbs
    cdef int Nall,N
    cdef np.ndarray g,vDir,waveDir,omegaM,phiM,kiM,kDirM,tanhFM,aiM
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
        double mwl,depth
    cdef double _cpp_eta(self, double* x, double t)
    cdef void _cpp_u(self, double* U, double* x, double t)

cdef class DirectionalWaves:
    cdef bool fast
    cdef double gAbs
    cdef int Nall,Mtot,N
    cdef np.ndarray vDir,omega,tanh,waveDir0,waveDirs,phiDirs,aiDirs,ki,kDirs,tanhF
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
        double mwl,depth
    cdef double _cpp_eta(self, double* x, double t)
    cdef void _cpp_u(self, double* U, double* x, double t)


# pointer to eta function
ctypedef double (*cfeta2) (TimeSeries, double* , double )  

# pointer to velocity function
ctypedef void (*cfvel2) (TimeSeries, double*, double* , double )

cdef class  TimeSeries:
    cdef bool fast,rec_direct
    cdef np.ndarray g,waveDir,vDir,x0,kDir,tanhF,time,etaS,ai,omega,phi,ki
    cdef double gAbs,depth,Tm,overlap,cutoff,setup,handover,Twindow,Tlag,Toverlap,dt,t0,tlength
    cdef int N,Nall,Nf,Nwaves,Nwindows
    cdef list windows_handover,windows_rec,decompose_window
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
        double wavelength,mwl
        object eta,u
    cdef cfeta2 _cpp_eta
    cdef cfvel2 _cpp_u
    cdef double _cpp_etaDirect(self, double* x, double t) 
    cdef double _cpp_etaWindow(self, double* x, double t) 
    cdef void _cpp_uDirect(self, double* U,double* x, double t) 
    cdef void _cpp_uWindow(self, double* U, double* x, double t) 

cdef class RandomNLWaves:
    cdef bool fast
    cdef np.ndarray omega,ki,kDir,phi,tanhKd,sinhKd,waveDir,ai
    cdef int N
    cdef double depth,gAbs
    cdef double* tanhKd_
    cdef double* sinhKd_
    cdef double* ai_
    cdef double* phi_
    cdef double* omega_
    cdef double* kDir_
    cdef double* ki_
    cdef double[30000] kDir_c
    cdef double[10000] omega_c
    cdef double[10000] ki_c
    cdef double[10000] ai_c
    cdef double[10000] tanh_c
    cdef double[10000] sinh_c
    cdef double[10000] phi_c
    cdef double _cpp_eta_2ndOrder(self,double* x, double t)
    cdef double _cpp_eta_short(self,double* x, double t)
    cdef double _cpp_eta_long(self,double* x, double t)
    cdef public:
        object eta
        object u
        object eta_linear
        

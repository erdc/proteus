cimport numpy
cdef extern from "math.h":
   double fabs(double x)
cdef extern from "transportCoefficients.h":
   double smoothedDirac(double eps,double phi)
   double smoothedHeaviside(double eps,double phi)

cdef inline double double_max(double a, double b): return a if a >= b else b
cdef inline double double_min(double a, double b): return a if a <= b else b

ctypedef numpy.double_t DTYPE_t
from math import sin,asin,cosh,sinh,cos,sqrt,pi
def monochromaticWave(int nElements, int nQuadraturePoints_element, double waveHeight, double waveCelerity, double waveFrequency,
                      double source_x0, double source_x1, double source_y0, double source_y1,
                      double eps,
                      numpy.ndarray[DTYPE_t,ndim=3] x,numpy.ndarray[DTYPE_t,ndim=2] r,double t):
    cdef int ie, k
    cdef double factor,source_area
    cdef double distance_x,distance_y,delta,dx_source,dy_source
    cdef double tmp_delta
    dx_source   = source_x1-source_x0
    dy_source   = source_y1-source_y0
    source_area = dx_source*dy_source
    factor = waveHeight/source_area*waveCelerity*sin(waveFrequency*t)#/sourceVolume
    #mwf debug
    #print "monochromaticWave t=%s x_0=%s x_1=%s y_0=%s y_1=%s factor=%s " % (t,source_x0,source_x1,source_y0,source_y1,factor)
    for ie in range(nElements):
        for k in range(nQuadraturePoints_element):
            #if (source_y0 <= x[ie,k,1] and x[ie,k,1] <= source_y1):
                #delta = smoothedDirac(eps,distance_x)
            distance_x = fabs(x[ie,k,0]-0.5*(source_x0+source_x1)) - 0.5*dx_source
            distance_y = fabs(x[ie,k,1]-0.5*(source_y0+source_y1)) - 0.5*dy_source
            tmp_delta =  (1.0-smoothedHeaviside(eps,distance_x))*(1.0-smoothedHeaviside(eps,distance_y))
            delta = tmp_delta #try quadratic instead tmp_delta 
            r[ie,k] -= delta*factor
            #mwf debug
            #tmp_delta = fabs(x[ie,k,0]-0.5*(source_x0+source_x1)) - 0.5*dy_source
            #if (1.0-smoothedHeaviside(eps,distance_x))*(1.0-smoothedHeaviside(eps,tmp_delta)) > 0.0 and delta < 1.0e-4:
            #    print "monochromaticWave t=%s x=%s y=%s factor=%s delta=%s tmp_delta= %s orig_delta= %s r[%s,%s]= %s" % (t,x[ie,k,0],x[ie,k,1],factor,delta,
            #                                                                                                              tmp_delta,(1.0-smoothedHeaviside(eps,distance_x))*(1.0-smoothedHeaviside(eps,tmp_delta)),
            #                                                                                                              ie,k,r[ie,k])

def secondOrderStokesWave(int nElements, int nQuadraturePoints_element, 
                          double waveHeight, double waveCelerity, double waveFrequency, double waveNumber,
                          double waterDepth,
                          double source_x0, double source_x1, double source_y0, double source_y1,
                          double eps,
                          numpy.ndarray[DTYPE_t,ndim=3] x,numpy.ndarray[DTYPE_t,ndim=2] r,double t):
    cdef int ie, k
    cdef double p_s,a_s,b_s,kd,term1,term2,term3
    cdef double distance_x,distance_y,delta,dx_source,dy_source,source_area,tmp_delta
    dx_source   = source_x1-source_x0
    dy_source   = source_y1-source_y0
    source_area = dx_source*dy_source

    kd  = waveNumber*waterDepth
    a_s = waveHeight*0.5 
    b_s = waveHeight*waveHeight*waveNumber*cosh(kd)*(2.0 + cosh(2.*kd)/(16.0+sinh(kd)**3))
    term1 = -a_s + sqrt(a_s*a_s + 8.0*b_s*b_s)/(4.0*b_s)
    p_s = asin(term1)
    term1 = waveCelerity*waveHeight*cos(pi*0.5 - waveFrequency*t - p_s)
    term2 = waveCelerity*waveHeight*waveHeight*cosh(kd)/(8.0*sinh(kd)**3)
    term3 = 2.0 + cosh(2.0*kd)*cos(2.0*(pi*0.5 - waveFrequency*t - p_s))
    factor = (term1 + term2*term3)/source_area
    for ie in range(nElements):
        for k in range(nQuadraturePoints_element):
            #if (source_y0 <= x[ie,k,1] and x[ie,k,1] <= source_y1):
            #    distance_x = x[ie,k,0]-0.5*(source_x0+source_x1)
            #    delta = smoothedDirac(eps,distance_x)
            distance_x = fabs(x[ie,k,0]-0.5*(source_x0+source_x1)) - 0.5*dx_source
            distance_y = fabs(x[ie,k,1]-0.5*(source_y0+source_y1)) - 0.5*dy_source
            tmp_delta = (1.0-smoothedHeaviside(eps,distance_x))*(1.0-smoothedHeaviside(eps,distance_y))
            delta     = tmp_delta #try quadratic instead tmp_delta?
            r[ie,k] -= factor*delta

def solitaryWave(int nElements, int nQuadraturePoints_element, 
                 double waveHeight, double waveCelerity, double waveFrequency,
                 double waterDepth,
                 double source_x0, double source_x1, double source_y0, double source_y1,
                 double eps,
                 numpy.ndarray[DTYPE_t,ndim=3] x,numpy.ndarray[DTYPE_t,ndim=2] r,double t):
    cdef int ie, k
    cdef double factor,term1,x_s,term2
    cdef double distance_x,distance_y,delta,dx_source,dy_source,source_area,tmp_delta
    dx_source   = source_x1-source_x0
    dy_source   = source_y1-source_y0
    source_area = dx_source*dy_source

    x_s  = 4.0*waveHeight/sqrt(waveHeight/waterDepth)
    term1= sqrt(3.0*waveHeight/(4.0*waterDepth*waterDepth*waterDepth))*(x_s - waveCelerity*t)
    term1= double_max(term1,-80.0)
    term1= double_min(term1, 80.0)
    #mwf debug
    #print "solitaryWave; waveHeight=%s waterDepth=%s x_s=%s waveCelerity= %s t=%s term1= %s" % (waveHeight,
    #                                                                                            waterDepth,
    #                                                                                            x_s,
    #                                                                                            waveCelerity,
    #                                                                                            t,
    #                                                                                            term1)
    term2= 1.0/(cosh(term1)+1.0e-12)
    factor = waveHeight*waveCelerity*term2*term2/source_area
    for ie in range(nElements):
        for k in range(nQuadraturePoints_element):
            #if (source_y0 <= x[ie,k,1] and x[ie,k,1] <= source_y1):
            #    distance_x = x[ie,k,0]-0.5*(source_x0+source_x1)
            #    delta = smoothedDirac(eps,distance_x)
            distance_x = fabs(x[ie,k,0]-0.5*(source_x0+source_x1)) - 0.5*dx_source
            distance_y = fabs(x[ie,k,1]-0.5*(source_y0+source_y1)) - 0.5*dy_source
            tmp_delta = (1.0-smoothedHeaviside(eps,distance_x))*(1.0-smoothedHeaviside(eps,distance_y))
            delta     = tmp_delta #try quadratic instead tmp_delta
            r[ie,k] -= delta*factor

def monochromaticWave3d(int nElements, int nQuadraturePoints_element, double waveHeight, double waveCelerity, double waveFrequency,
                        double source_x0, double source_x1, double source_y0, double source_y1, double source_z0, double source_z1,
                        double eps,
                        numpy.ndarray[DTYPE_t,ndim=3] x,numpy.ndarray[DTYPE_t,ndim=2] r,double t):
    cdef int ie, k
    cdef double factor,source_volume
    cdef double distance_x,distance_y,distance_z,delta,dx_source,dy_source,dz_source
    dx_source   = source_x1-source_x0
    dy_source   = source_y1-source_y0
    dz_source   = source_z1-source_z0
    source_volume = dx_source*dy_source*dz_source
    factor = waveHeight/source_volume*waveCelerity*sin(waveFrequency*t)#
    #mwf debug
    #print "monochromaticWave t=%s x_0=%s x_1=%s y_0=%s y_1=%s factor=%s " % (t,source_x0,source_x1,source_y0,source_y1,factor)
    for ie in range(nElements):
        for k in range(nQuadraturePoints_element):
            #if (source_y0 <= x[ie,k,1] and x[ie,k,1] <= source_y1):
                #delta = smoothedDirac(eps,distance_x)
            distance_x = fabs(x[ie,k,0]-0.5*(source_x0+source_x1)) - 0.5*dx_source
            distance_y = fabs(x[ie,k,1]-0.5*(source_y0+source_y1)) - 0.5*dy_source
            distance_z = fabs(x[ie,k,2]-0.5*(source_z0+source_z1)) - 0.5*dz_source
            delta = (1.0-smoothedHeaviside(eps,distance_x))*(1.0-smoothedHeaviside(eps,distance_y))*(1.0-smoothedHeaviside(eps,distance_z))
            r[ie,k] -= delta*factor
                #mwf debug
                #print "monochromaticWave t=%s x=%s y=%s z=%s factor=%s delta=%s r[%s,%s]= %s" % (t,x[ie,k,0],x[ie,k,1],x[ie,k,2],factor,delta,ie,k,r[ie,k])

def secondOrderStokesWave3d(int nElements, int nQuadraturePoints_element, 
                            double waveHeight, double waveCelerity, double waveFrequency, double waveNumber,
                            double waterDepth,
                            double source_x0, double source_x1, double source_y0, double source_y1, double source_z0, double source_z1,
                            double eps,
                            numpy.ndarray[DTYPE_t,ndim=3] x,numpy.ndarray[DTYPE_t,ndim=2] r,double t):
    cdef int ie, k
    cdef double p_s,a_s,b_s,kd,term1,term2,term3
    cdef double distance_x,distance_y,distance_z,delta,dx_source,dy_source,dz_source,source_volume
    dx_source   = source_x1-source_x0
    dy_source   = source_y1-source_y0
    dz_source   = source_z1-source_z0
    source_volume = dx_source*dy_source*dz_source

    kd  = waveNumber*waterDepth
    a_s = waveHeight*0.5 
    b_s = waveHeight*waveHeight*waveNumber*cosh(kd)*(2.0 + cosh(2.*kd)/(16.0+sinh(kd)**3))
    term1 = -a_s + sqrt(a_s*a_s + 8.0*b_s*b_s)/(4.0*b_s)
    p_s = asin(term1)
    term1 = waveCelerity*waveHeight*cos(pi*0.5 - waveFrequency*t - p_s)
    term2 = waveCelerity*waveHeight*waveHeight*cosh(kd)/(8.0*sinh(kd)**3)
    term3 = 2.0 + cosh(2.0*kd)*cos(2.0*(pi*0.5 - waveFrequency*t - p_s))
    factor = (term1 + term2*term3)/source_volume
    for ie in range(nElements):
        for k in range(nQuadraturePoints_element):
            #if (source_y0 <= x[ie,k,1] and x[ie,k,1] <= source_y1):
            #    distance_x = x[ie,k,0]-0.5*(source_x0+source_x1)
            #    delta = smoothedDirac(eps,distance_x)
            distance_x = fabs(x[ie,k,0]-0.5*(source_x0+source_x1)) - 0.5*dx_source
            distance_y = fabs(x[ie,k,1]-0.5*(source_y0+source_y1)) - 0.5*dy_source
            distance_z = fabs(x[ie,k,2]-0.5*(source_z0+source_z1)) - 0.5*dz_source
            delta = (1.0-smoothedHeaviside(eps,distance_x))*(1.0-smoothedHeaviside(eps,distance_y))*(1.0-smoothedHeaviside(eps,distance_z))
            r[ie,k] -= factor*delta

def solitaryWave3d(int nElements, int nQuadraturePoints_element, 
                   double waveHeight, double waveCelerity, double waveFrequency,
                   double waterDepth,
                   double source_x0, double source_x1, double source_y0, double source_y1, double source_z0, double source_z1,
                   double eps,
                   numpy.ndarray[DTYPE_t,ndim=3] x,numpy.ndarray[DTYPE_t,ndim=2] r,double t):
    cdef int ie, k
    cdef double factor,term1,x_s,term2
    cdef double distance_x,distance_y,distance_z,delta,dx_source,dy_source,dz_source,source_volume
    dx_source   = source_x1-source_x0
    dy_source   = source_y1-source_y0
    dz_source   = source_z1-source_z0
    source_volume = dx_source*dy_source*dz_source

    x_s  = 4.0*waveHeight/sqrt(waveHeight/waterDepth)
    term1= sqrt(3.0*waveHeight/(4.0*waterDepth*waterDepth*waterDepth))*(x_s - waveCelerity*t)
    term1= double_max(term1,-80.0)
    term1= double_min(term1, 80.0)
    #mwf debug
    #print "solitaryWave; waveHeight=%s waterDepth=%s x_s=%s waveCelerity= %s t=%s term1= %s" % (waveHeight,
    #                                                                                            waterDepth,
    #                                                                                            x_s,
    #                                                                                            waveCelerity,
    #                                                                                            t,
    #                                                                                            term1)
    term2= 1.0/(cosh(term1)+1.0e-12)
    factor = waveHeight*waveCelerity*term2*term2/source_volume
    for ie in range(nElements):
        for k in range(nQuadraturePoints_element):
            #if (source_y0 <= x[ie,k,1] and x[ie,k,1] <= source_y1):
            #    distance_x = x[ie,k,0]-0.5*(source_x0+source_x1)
            #    delta = smoothedDirac(eps,distance_x)
            distance_x = fabs(x[ie,k,0]-0.5*(source_x0+source_x1)) - 0.5*dx_source
            distance_y = fabs(x[ie,k,1]-0.5*(source_y0+source_y1)) - 0.5*dy_source
            distance_z = fabs(x[ie,k,2]-0.5*(source_z0+source_z1)) - 0.5*dz_source
            delta = (1.0-smoothedHeaviside(eps,distance_x))*(1.0-smoothedHeaviside(eps,distance_y))*(1.0-smoothedHeaviside(eps,distance_z))
            r[ie,k] -= delta*factor


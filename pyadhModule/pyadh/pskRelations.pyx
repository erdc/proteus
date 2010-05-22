import numpy
cimport numpy

cdef extern from "math.h":
   double fabs(double x)
   double pow(double x, double y)
   double sqrt(double x)
   double exp(double x)
cdef inline double double_max(double a, double b): return a if a >= b else b
cdef inline double double_min(double a, double b): return a if a <= b else b

ctypedef numpy.double_t DTYPE_t
ctypedef int ITYPE_t

########################################################################
#play with cython implementations of psk relations rather than wrapping?
########################################################################



cdef class cPskRelation:
#Try to bring psk relation hierarchy into cython and allow more direct access?

  #only accessible from cython unless add public or readonly tag
  cdef readonly double Se, dSe_dSw, Sw_min, Sw_max, krw, dkrw, krn, dkrn, psic, dpsic

  def __init__(self, double Sw_min=0.0, double Sw_max=1.0):
      self.Sw_min=Sw_min; self.Sw_max=Sw_max

  cdef void calc(self,double Sw):
      self.calc_Se(Sw)
      self.krw = self.Se
      self.dkrw = self.dSe_dSw
    
      self.krn  = (1.0-self.Se)
      self.dkrn = -self.dSe_dSw
      
      self.psic  = self.Se
      self.dpsic = self.dSe_dSw   
      
  cdef void calc_Se(self,double Sw):
      if Sw < self.Sw_min: 
          self.Se = 0.0; self.dSe_dSw = 0.0
      elif Sw > self.Sw_max:
          self.Se = 1.0; self.dSe_dSw = 0.0
      else:
          self.Se = (Sw - self.Sw_min)/(self.Sw_max-self.Sw_min)
          self.dSe_dSw = 1.0/(self.Sw_max-self.Sw_min)

cdef class cSimplePsk(cPskRelation):
#quadratic relative perms
  def __init__(self, double Sw_min=0.0, double Sw_max=1.0):
      cPskRelation.__init__(self,Sw_min,Sw_max)

  cdef void calc(self,double Sw):
      self.calc_Se(Sw)
      self.krw = self.Se*self.Se
      self.dkrw = 2.0*self.Se*self.dSe_dSw;
    
      self.krn  = (1.0-self.Se)*(1.0-self.Se);
      self.dkrn = 2.0*(self.Se-1.0)*self.dSe_dSw;
      
      self.psic  = -self.Se;
      self.dpsic = -self.dSe_dSw;   


cdef class cVGM(cPskRelation):
#van Genuchten-mualem
  cdef double alpha,n,m,eps_se
  def __init__(self, double Sw_min=0.0, double Sw_max=1.0,alpha=1.0,m=0.5,eps_se=1.0e-2):
      cPskRelation.__init__(self,Sw_min,Sw_max)
      self.alpha=alpha; self.m=m; self.n=1.0/(1.0-self.m)
      self.eps_se=eps_se
  cdef void calc(self,double Sw):
      cdef double Seovmmo,Seovm,Seovmmh,S,Smmo,S2mmo,Sm,S2m,pcesub1,pcesub2,Se_eps
      cdef double m = self.m, n=self.n, alpha=self.alpha,Se,dSe_dSw
      self.calc_Se(Sw)
      Se_eps = double_max(self.Se,self.eps_se)
      Se_eps = double_min(Se_eps,1.0-self.eps_se)
      Se = self.Se; dSe_dSw = self.dSe_dSw

      Seovmmo  = pow(Se_eps,((1.0/m)-1.0));
      Seovm    = Seovmmo*Se_eps;
      Seovmmh  = pow(Se_eps,((1.0/m)-0.5)); 
      S       = 1.0 - Seovm;
      Smmo    = pow(S,m-1.0);
      S2mmo   = pow(S,2.0*m-1.0);
      Sm      = Smmo*S;
      S2m     = S2mmo*S; 
      
      pcesub1 = pow(((1.0/Seovm)-1.0),((1.0/n)-1.0));
      pcesub2 = pcesub1*((1.0/Seovm)-1.0);
      #
      self.krw  = sqrt(Se)*(1.0-Sm)*(1.0-Sm);
      self.dkrw = (0.5*(1.0/sqrt(Se_eps))*(1.0-Sm)*(1.0-Sm) + 2.0*(1-Sm)*Smmo*Seovmmh)*dSe_dSw;
  
      self.krn  = sqrt(1.0-Se)*S2m;
      self.dkrn = (-0.5*(1.0/sqrt(1.0-Se_eps))*S2m - 2.0*sqrt(1.0-Se)*(S2mmo)*Seovmmo)*dSe_dSw; 
      
      self.psic  = pow((pow(Se_eps,(-1.0/m)) - 1.0),(1.0/n))/alpha;
      self.dpsic = ((-1.0/(alpha*n*m))*(pow((pow(Se_eps,-1.0/m)-1.0),-m))*(pow(Se_eps,((-1.0/m)-1.0))))*dSe_dSw;  

#########################
cdef class cDensityRelation:
  cdef readonly double rho,drho
  
  def __init__(self,rho=1.0):
      self.rho=rho; self.drho=0.0

  cdef void calc(self, double psi):
      pass

cdef class cExponentialDensity(cDensityRelation):
   cdef readonly rho_0,psi_0,beta
   
   def __init__(self,rho_0=1.0,psi_0=0.0,beta=0.0):
       cDensityRelation.__init__(self,rho_0)
       self.rho_0=rho_0; self.psi_0=psi_0; self.beta=beta
       
   cdef void calc(self,double psi):
       self.rho=self.rho_0*exp(self.beta*(psi-self.psi_0))
       self.drho=self.beta*self.rho

cdef class cIdealGasDensity(cDensityRelation):
  cdef readonly double T,W,R,convFactor,WoverRT,rho_0,psi_0 

  def __init__(self,T,W,R,convFactor,rho_0,psi_0):
      self.T=T; self.W=W; self.R=R; self.convFactor=convFactor; self.rho_0=rho_0; self.psi_0=psi_0
      self.WoverRT = self.convFactor*self.W/(self.R*self.T)
  cdef void calc(self, double psi):
      self.rho = self.rho_0 + (psi-self.psi_0)*self.WoverRT
      self.drho= self.WoverRT


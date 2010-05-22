from math import *
import numpy
#John Chrispell, Summer 07
import pyadh
from pyadh.TransportCoefficients import *

class Nonlinear_VA_linearDR_Coefficients(TC_base):
    """
    This class implements constant coefficients with no cross-component diffusion:
    
    (Mu)_t + \deld (Bu - A \grad u) + C u = 0 
    """
#    from ctransportCoefficients import linearADR_ConstantCoefficientsEvaluate
    def __init__(self,nc=1,q=[1.0],a=1.0e-3,mu_r=0.5):
    	self.nc=nc
	self.variableNames=['S']
        mass={}
        advection={}
        diffusion={}
        potential={}
        reaction={}
        hamiltonian={}
        for i in range(nc):
            mass[i]      = {i:'linear'}
            advection[i] = {i:'nonlinear'}
            diffusion[i] = {i: {i:'constant'}}
            potential[i] = {i: 'u'}
            reaction[i]  = {i: 'linear'}
        TC_base.__init__(self,
                         nc,
                         mass,
                         advection,
                         diffusion,
                         potential,
                         reaction,
                         hamiltonian)
        self.q=q
	self.a=a
	self.mu_r = mu_r #mobility ratio more or less
    # Simple Model Definitions (Buckley Leverett)	
    def krw_simp(self,u):
    	return u**2
    def krn_simp(self,u):
        return (1.0-u)**2
    def dkrw_simp(self,u):
    	return 2.0*u
    def dkrn_simp(self,u):
    	return u-1.0
    	
    def f(self,krw,krn):
	return krw/(krw+self.mu_r*krn)
    def df(self,krw,dkrw,krn,dkrn):
	return (dkrw*(krw+self.mu_r*krn) - krw*(dkrw+self.mu_r*dkrn)) / (krw+self.mu_r*krn)**2
	
    # Set the Values of krw and krn 
    #    Choices: *_simp
    
    krw = krw_simp
    dkrw = dkrw_simp
    krn = krn_simp
    dkrn = dkrn_simp
    
    def evaluate(self,t,c):
    	#loop over components of the equation system
        for  i  in range(self.nc):
		space_dim = c[('f',i)].shape[-1]
		#loop over the points at which the coefficients are evaluated
		for pi in range(len(c[('u',0)].flat)):
			#compute the p-s-k relations
			u_pi=c[('u',i)].flat[pi]
			krw_pi = self.krw(u_pi)
			krn_pi = self.krn(u_pi)
			dkrw_pi = self.dkrw(u_pi)
			dkrn_pi = self.dkrn(u_pi)
			#now set the coefficients of the ADR equation
			c[('m',i)].flat[pi] = u_pi
			c[('dm',i,i)].flat[pi] = 1.0
			#loop over the space dimension of the problem
			for I in range(space_dim):
				c[('f',i)].flat[pi*space_dim+I] = self.q[I]*self.f(krw_pi,krn_pi)	
				c[('df',i,i)].flat[pi*space_dim+I] = self.q[I]*self.df(krw_pi,dkrw_pi,krn_pi,dkrn_pi)	
				c[('a',i,i)].flat[pi*space_dim**2+I*space_dim+I] = self.a
			#since a is constant, no need to eval da
			#since phi=u, no need to eval phi or dphi
			#r is zero in this problem for now
			
			
			
			
			
			
#encapsulate analytical solver for BuckleyLeverett Riemann problems
from pyadh.ObjectiveFunctions import OsherFuncCoef
from pyadh.Optimizers import fminbound
class Buckley_Leverett_RiemannSoln:
    def __init__(self,coefficients,uLeft=1.0,uRight=0.0,t0=0.0,x0=0.0,T=0.5,ftol=1.0e-8):
        #ought to start with range and spline for eval later I guess
        self.coefficients = coefficients
        self.uL = uLeft; self.uR = uRight;
        self.t0 = t0; self.x0 = x0;
        self.T  = T;  self.ftol = ftol
        self.riemF = OsherFuncCoef(self.uL,self.uR,self.coefficients,self.T-self.t0,self.x0)
        self.solver= fminbound(self.riemF,self.ftol)

    def uOfXT(self,x,t):
        if abs(t-self.t0) < 1.0e-7:
            if x[0]-self.x0 <= 0.0:
                return self.uL
            return self.uR
        self.riemF.xi = (x[0]-self.x0)/(t-self.t0)
        u,f = self.solver.solve(Guess_x = 0.5*(self.uL+self.uR))
        #mwf debug
        #import pdb
        #pdb.set_trace()
        return u
    def uOfX(self,x):
        return self.uOfXT(x,self.T)
    

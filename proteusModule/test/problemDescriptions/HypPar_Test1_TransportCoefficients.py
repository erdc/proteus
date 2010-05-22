from math import *
import numpy
#John Chrispell, Spring 08
import pyadh
from pyadh.TransportCoefficients import *

class HypPar_Test1_Coefficients(TC_base):
    """
    This class implements constant coefficients with no cross-component diffusion:
    
    (Mu)_t + \deld (Bu - A \grad u) + C u = 0 
    """
#   from ctransportCoefficients import linearADR_ConstantCoefficientsEvaluate
    def __init__(self,nc=1,q=[1.0],alpha = 0.05,
                                   rFuncClass = None, 
				   setParamsFunc = 1.0):
    	self.nc=nc
	self.variableNames=['u']
        mass={}
        advection={}
        diffusion={}
        potential={}
        reaction={}
        hamiltonian={}
        for i in range(nc):
            mass[i]      = {i:'linear'}
            advection[i] = {i:'linear'} # May set to nonlinear in future. 
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
	self.alpha = alpha	
	
	if rFuncClass != None:
	    self.rFuncClass = rFuncClass()
	if rFuncClass == None:
	    self.rFuncClass = None
	self.setParamsFunc = setParamsFunc
	
    def initializeElementQuadrature(self,t,cq):
	self.alpha_q= Numeric.zeros(cq[('u',0)].shape,Numeric.Float)
	self.setParamsFunc(cq['x'],
                           self.alpha_q)
	
    def initializeElementBoundaryQuadrature(self,t,cebq,cebq_global):
        self.alpha_ebq= Numeric.zeros(cebq[('u',0)].shape,Numeric.Float)
        self.setParamsFunc(cebq['x'],
                           self.alpha_ebq)
			   
    def evaluate(self,t,c):
    
        if c[('u',0)].shape == self.alpha_q.shape:
	    alpha = self.alpha_q
	else:
	    alpha = self.alpha_ebq
    
    	#loop over components of the equation system
        for  i  in range(self.nc):
	    space_dim = c[('f',i)].shape[-1]
	    #print "space_dim = c[('f',i)].shape[-1]= ", c[('f',i)].shape[-1] 
	    #loop over the points at which the coefficients are evaluated
	    
	    for pi in range(len(c[('u',0)].flat)):
		u_pi=c[('u',i)].flat[pi]
		alpha_pi = alpha.flat[pi]
		#now set the coefficients of the ADR equation
		c[('m',i)].flat[pi] = u_pi
		c[('dm',i,i)].flat[pi] = 1.0
		if self.rFuncClass != None:
		    RHS_value  = self.rFuncClass.rOfXT(c['x'].flat[pi*3:(pi+1)*3],t)    
		    c[('r',i)].flat[pi]    = 0.0 - RHS_value
		    #c[('dr',i,i)].flat[pi] = 0.0  
		if self.rFuncClass == None: 
		    print "+++++++++++++++++++++++++++++++++++++++++++++++ Boing +++++++" 
		    c[('r',i)].flat[pi] = u_pi  
		    c[('dr',i,i)].flat[pi] = 1.0 
		
		#loop over the space dimension of the problem
		for I in range(space_dim):
		    c[('f',i)].flat[pi*space_dim+I] = self.q[I]*u_pi #self.f(krw_pi,krn_pi)	
		    c[('df',i,i)].flat[pi*space_dim+I] = self.q[I]	
		    c[('a',i,i)].flat[pi*space_dim**2+I*space_dim+I] = alpha_pi
		    #print "-----------------------------------------------------------------> alpha_pi", alpha_pi

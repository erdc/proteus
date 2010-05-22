from math import *
import numpy
#John Chrispell, Spring 08
import pyadh
from pyadh.TransportCoefficients import *

class HypPar_Coefficients(TC_base):
    """
    This class implements constant coefficients with no cross-component diffusion:
    
    (Mu)_t + \deld (Bu - A \grad u) + C u = 0 
    """
#   from ctransportCoefficients import linearADR_ConstantCoefficientsEvaluate
    def __init__(self,nc=1,q=[1.0],alpha = 0.05,
                                   delta = 1.0, 
                                   phase = "ADR_start",
				   rFuncClass = None, 
				   setParamsFunc = 1.0):
    	self.nc=nc
	if phase == 'ADR_start':
	    variableNames=['u']
	    
	if phase == 'ADR_finish':
	    variableNames=['ufinal']
	    
	if phase == 'correctionOne':
	    variableNames=['ubar']
	    	
	if phase == 'correctionTwo':
	    variableNames=['udoublebar']
	#else: 
	#    variableNames=['Not Set']
	   
        mass={}
        advection={}
        diffusion={}
        potential={}
        reaction={}
        hamiltonian={}
        for i in range(nc):
            if phase in ['ADR_start','ADR_finish']:
                mass[i]      = {i:'linear'}
                advection[i] = {i:'linear'} # May set to nonlinear in future. 
                diffusion[i] = {i: {i:'constant'}}
                potential[i] = {i: 'u'}
                reaction[i]  = {i: 'linear'}
            else:
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
                         hamiltonian,
			 variableNames)		 
        self.q=q
	self.alpha = alpha
	self.phase = phase
	self.delta = delta
	
	if rFuncClass != None:
	    self.rFuncClass = rFuncClass()
	if rFuncClass == None:
	    self.rFuncClass = None
	    
	self.setParamsFunc = setParamsFunc
	self.ubarValue = [0.0]
	self.uValue = [0.0]
	
	# adding stuff for a split operator
	if phase == "ADR_start":
	    self.nModel = 2
	    self.q_ubar = None #Numeric.array(self.ubarValue)    # ubar at quadrature points
	    self.ebq_ubar = None #Numeric.array(self.ubarValue)  # ubar at element boundary quadrature points
	
	if phase == "ADR_finish":
	    self.nModel = 2
	    self.q_ubar = None #Numeric.array(self.ubarValue)    # ubar at quadrature points
	    self.ebq_ubar = None #Numeric.array(self.ubarValue)  # ubar at element boundary quadrature points

	if phase == "correctionOne":
	    self.nModel = 0
	    self.q_u   = Numeric.array(self.uValue) # solution u at quadrature points
	    self.ebq_u = Numeric.array(self.uValue) # solution u at element boundary quadrature points	
	
	if phase == "correctionTwo":
	    self.nModel = 1
	    self.q_u   = Numeric.array(self.uValue) # solution u at quadrature points
	    self.ebq_u = Numeric.array(self.uValue) # solution u at element boundary quadrature points
	    
    def attachModels(self,modelList):
        self.startModel = modelList[0]
        self.finishModel=modelList[-1]
        if self.phase == "ADR_start":
	    #print "----------------------------------------------------------> ADR_Start" 
	    self.q_ubar   = modelList[self.nModel].q[('u',0)]
	    self.ebq_ubar = modelList[self.nModel].ebq[('u',0)]
	if self.phase == "correctionOne":
	    #print "----------------------------------------------------------> correction one" 
	    self.q_u = modelList[self.nModel].q[('u',0)]
	    self.ebq_u = modelList[self.nModel].ebq[('u',0)]
	if self.phase == "correctionTwo":
	    #print "----------------------------------------------------------> correction two" 
	    self.q_u = modelList[self.nModel].q[('u',0)]
	    self.ebq_u = modelList[self.nModel].ebq[('u',0)]
	    #print self.q_u
	    #print self.nModel
	if self.phase == "ADR_finish":
	    #print "----------------------------------------------------------> ADR_finish"
	    self.q_ubar   = modelList[self.nModel].q[('u',0)]
	    self.ebq_ubar = modelList[self.nModel].ebq[('u',0)]   
	     
    def initializeMesh(self,mesh):
        pass
	
    def initializeElementQuadrature(self,t,cq):
	self.alpha_q= Numeric.zeros(cq[('u',0)].shape,Numeric.Float)
	self.setParamsFunc(cq['x'],
                           self.alpha_q)
	# set up dummy values incase we arn't running the other model. 
	if self.phase == "ADR_start":
	    self.q_ubar = Numeric.ones(cq[('u',0)].shape,Numeric.Float)
	    self.q_ubar.flat[:] = self.ubarValue	
	if self.phase == "ADR_finish":
	    self.q_ubar = Numeric.ones(cq[('u',0)].shape,Numeric.Float)
	    self.q_ubar.flat[:] = self.ubarValue
	if self.phase == "correctionOne":
	    self.q_u = Numeric.ones(cq[('u',0)].shape,Numeric.Float)
	    self.q_u.flat[:] = self.uValue	
	if self.phase == "correctionTwo":
	    self.q_u = Numeric.ones(cq[('u',0)].shape,Numeric.Float)
	    self.q_u.flat[:] = self.uValue
	
    def initializeElementBoundaryQuadrature(self,t,cebq,cebq_global):
        self.alpha_ebq= Numeric.zeros(cebq[('u',0)].shape,Numeric.Float)
        self.setParamsFunc(cebq['x'],
                           self.alpha_ebq)
	# set up dummy values incase we arn't running the other model. 
	if (self.phase == "ADR_start" or self.phase == "ADR_finish"):
	    self.ebq_ubar = Numeric.ones(cebq[('u',0)].shape,Numeric.Float)
	    self.ebq_ubar.flat[:] = self.ubarValue 
	if (self.phase == "correctionOne" or self.phase == "correctionTwo"):
	    self.ebq_u = Numeric.ones(cebq[('u',0)].shape,Numeric.Float)
	    self.ebq_u.flat == self.uValue
	
	   
    def preStep(self,t,firstStep=False):
        if self.phase == 'ADR_start':
            self.startModel.u[0].dof.flat[:]  = self.finishModel.u[0].dof.flat[:]
            self.startModel.calculateCoefficients()
            self.startModel.calculateElementResidual()
            self.startModel.updateTimeHistory(t,resetFromDOF=True)
        copyInstructions = {}#{'copy_uList':True,'uList_model':self.nModelId}
        return copyInstructions
    def evaluate(self,t,c):
    	#c['x']
        if((self.phase == 'ADR_start') or (self.phase == 'ADR_finish')):
	    #print "------------------------------------------------------------- Phase = ", self.phase
            if c[('u',0)].shape == self.alpha_q.shape:
	        alpha = self.alpha_q
		ubar = self.q_ubar
		#print "ubar = ", ubar
	    else:
	        alpha = self.alpha_ebq
	        ubar  = self.ebq_ubar
	    #print ubar    
	    #raw_input('ADR wait')
    	    #loop over components of the equation system
            for  i  in range(self.nc):
	        space_dim = c[('f',i)].shape[-1]
	        #print "space_dim = c[('f',i)].shape[-1]= ", c[('f',i)].shape[-1] 
	        #loop over the points at which the coefficients are evaluated
	        nSpace=c['x'].shape[-1]
		for pi in range(len(c[('u',0)].flat)):
		    #print "pi = ", pi
		    u_pi=c[('u',i)].flat[pi]
		    alpha_pi = alpha.flat[pi]
		    ubar_pi  = ubar.flat[pi]
		    #print "Phase =", self.phase, " ubar_pi = ", ubar_pi
		    #print "u_pi =", u_pi
		    #now set the coefficients of the ADR equation
		    c[('m',i)].flat[pi] = u_pi
		    c[('dm',i,i)].flat[pi] = 1.0		    
		    if self.rFuncClass != None:
			RHS_value  = self.rFuncClass.rOfXT(c['x'].flat[pi*3:(pi+1)*3],t) 
		        c[('r',i)].flat[pi] = u_pi - ubar_pi - RHS_value
		        c[('dr',i,i)].flat[pi] = 1.0      
		    if self.rFuncClass == None: 
		        c[('r',i)].flat[pi] = u_pi - ubar_pi 
		        c[('dr',i,i)].flat[pi] = 1.0                     
# 		        c[('r',i)].flat[pi] =0.0
# 		        c[('dr',i,i)].flat[pi] = 0.0
		    #loop over the space dimension of the problem
		    for I in range(space_dim):
		        c[('f',i)].flat[pi*space_dim+I] = self.q[I]*1.0*u_pi  
		        c[('df',i,i)].flat[pi*space_dim+I] = self.q[I]*1.0
		        c[('a',i,i)].flat[pi*space_dim**2+I*space_dim+I] = alpha_pi 
        if((self.phase == 'correctionOne') or (self.phase == 'correctionTwo')):
	    #print "------------------------------------------------------------- Phase = ", self.phase
            if c[('u',0)].shape == self.q_u.shape:
	        sol_u = self.q_u   
	    else:
	        sol_u = self.ebq_u
    	    #print sol_u
	    #raw_input('wait in correction')
    	    #loop over components of the equation system
            for  i  in range(self.nc):
	        space_dim = c[('a',i,i)].shape[-1]
	        #loop over the points at which the coefficients are evaluated
	        for pi in range(len(c[('u',0)].flat)):
		    ubar_pi  = c[('u',i)].flat[pi]
		    sol_u_pi = sol_u.flat[pi]
		    #print "ubar_pi = ", ubar_pi
		    #now set the coefficients of the correction equation
		    #c[('m',i)].flat[pi] = ubar_pi
		    #c[('dm',i,i)].flat[pi] = 1.0
		    c[('r',i)].flat[pi] = ubar_pi - sol_u_pi
		    c[('dr',i,i)].flat[pi] = 1.0
		    #loop over the space dimension of the problem
		    for I in range(space_dim):
		        #c[('f',i)].flat[pi*space_dim+I] = 0.0 
		        c[('a',i,i)].flat[pi*space_dim**2+I*space_dim+I] = 0#self.delta
		        

from proteus.TransportCoefficients import *
import numpy as np


class Transport_Dilute_Dispersion(LinearVADR_ConstantCoefficients):
	from proteus.dilutetransport.cdilutetransportCoefficients import DiluteDispersionEvaluate
	def __init__(self, nc=1, poro = 0, diff=0, alpha_L=0, diff_ModelId = None, adv_ModelId = None, mom_ModelId = None):
		mass={}
		advection={}
		diffusion={}
		potential={}
		reaction={}
		hamiltonian={}
		for i in range(nc):
			mass[i]      = {i:'nonlinear'}
			advection[i] = {i:'constant'}
			diffusion[i] = {i: {i:'linear'}}
			potential[i] = {i: 'u'}
			reaction[i]  = {i: 'constant'}
		TC_base.__init__(self,nc,mass,advection,diffusion,potential,reaction,hamiltonian)
		self.q={}; self.ebqe={}; self.ebq ={}; self.ebq_global = {}
		self.diff = diff
		self.poro = poro
		self.alpha_L = alpha_L
		self.adv_ModelId  = adv_ModelId
		self.diff_ModelId = diff_ModelId
		self.mom_ModelId  = mom_ModelId
		self.variableNames=['u_diff']
		self.nd = 1

	def attachModels(self,modelList):
		self.mom_Model = modelList[self.mom_ModelId]
		self.adv_Model = modelList[self.adv_ModelId]
		self.diff_Model = modelList[self.diff_ModelId]
		self.q_vel  = self.mom_Model.q[('velocity',0)]/self.poro
		self.ebqe_vel = self.mom_Model.ebqe[('velocity',0)]/self.poro
		if self.mom_Model.ebq.has_key(('velocity',0)):
			self.ebq_vel = self.mom_Model.ebq[('velocity',0)]/self.poro
		if self.mom_Model.ebq_global.has_key(('velocity',0)):
			self.ebq_global_vel = self.mom_Model.ebq_global[('velocity',0)] /self.poro


	def initializeElementQuadrature(self,t,cq):
		self.q_x_shape = cq['x'].shape
		self.q_vel = numpy.zeros((cq['x'].shape[0],cq['x'].shape[1],self.nd),'d')

	def initializeElementBoundaryQuadrature(self,t,cebq,cebq_global):
		self.ebq_global_x_shape = cebq_global['x'].shape
		self.ebq_x_shape = cebq['x'].shape
		self.ebq_vel = numpy.zeros((cebq['x'].shape[0],cebq['x'].shape[1],cebq['x'].shape[2],self.nd),'d')
		self.ebq_global_vel = numpy.zeros((cebq_global['x'].shape[0],cebq_global['x'].shape[1],self.nd),'d')

	def initializeGlobalExteriorElementBoundaryQuadrature(self,t,cebqe):
		self.ebqe_x_shape = cebqe['x'].shape
		self.ebqe_vel = numpy.zeros((cebqe['x'].shape[0],cebqe['x'].shape[1],self.nd),'d')

 	def preStep(self,t,firstStep=True):
# 		#import pdb; pdb.set_trace()
 		self.diff_Model.timeIntegration.m_tmp = self.adv_Model.timeIntegration.m_tmp
 		self.diff_Model.timeIntegration.m_last = self.adv_Model.timeIntegration.m_tmp
		self.diff_Model.u[0].dof = self.adv_Model.u[0].dof
 		self.diff_Model.calculateCoefficients()
 		self.diff_Model.calculateElementResidual()
 		self.diff_Model.timeIntegration.updateTimeHistory(resetFromDOF=True)
 		self.diff_Model.timeIntegration.resetTimeHistory(resetFromDOF=True)
 		self.diff_Model.updateTimeHistory(t,resetFromDOF=True)
 		copyInstructions = {'copy_uList':True,'uList_model':self.adv_ModelId}
 		copyInstructions = {'reset_uList':True}
 		return copyInstructions

	def evaluate(self,t,c): 
		#import pdb; pdb.set_trace()
		if self.q_vel.shape == c[('df',0,0)].shape:
			velocity = self.q_vel
		elif self.ebqe_vel.shape == c[('df',0,0)].shape:
			velocity = self.ebqe_vel
		elif self.ebq_vel.shape == c[('df',0,0)].shape:
			velocity = self.ebq_vel
		else:
			print "no v---------------------"
			raise RuntimeError
		self.DiluteDispersionEvaluate(self.poro,
                                    self.diff,
                                    self.alpha_L,
                                    c[('u',0)],
                                    c[('m',0)],
                                    c[('dm',0,0)],
                                    c[('f',0)],
                                    c[('df',0,0)],
                                    c[('a',0,0)],
                                    c[('da',0,0,0)],
                                    velocity)


class Porous_Media_Flow(TC_base):
	from proteus.dilutetransport.cdilutetransportCoefficients import ConMassFluidEvaluate
	def __init__(self,nc=1,nd=1,diff_ModelId = None, adv_ModelId = None, mom_ModelId = None, K = 0.0 , grav = -980., L = 0):
		mass={}
		advection={}
		diffusion={}
		potential={}
		reaction={}
		hamiltonian={}
		self.nd=nd
		mass = {0:{0:'linear'}}
		advection = {0:{0:'linear'}}
		diffusion = {0:{0:{0:'nonlinear'}}}
		potential = {0:{0:'nonlinear'}}
		reaction = {0:{0:'constant'}}
		TC_base.__init__(self,
                         nc,
                         mass,
                         advection,
                         diffusion,
                         potential,
                         reaction,
                         hamiltonian)
		self.adv_ModelId  = adv_ModelId
		self.diff_ModelId = diff_ModelId
		self.mom_ModelId  = mom_ModelId
		self.K = K
		self.grav = grav
		self.L = L
		self.variableNames=['pressure']
        
	def attachModels(self,modelList):
		self.mom_Model = modelList[self.mom_ModelId]
		if (self.adv_ModelId is not None):
			self.adv_Model = modelList[self.adv_ModelId]
		self.diff_Model = modelList[self.diff_ModelId]
		self.q_c  = self.diff_Model.q[('u',0)]
		self.ebqe_c = self.diff_Model.ebqe[('u',0)]
		if self.diff_Model.phi_ip.has_key(('u',0)):
			self.cip_c = self.diff_Model.phi_ip[('u',0)]
		if self.diff_Model.ebq.has_key(('u',0)):
			self.ebq_c = self.diff_Model.ebq[('u',0)]
		if self.diff_Model.ebq_global.has_key(('u',0)):
			self.ebq_global_c = self.diff_Model.ebq_global[('u',0)] 

	def initializeElementQuadrature(self,t,cq):
		self.q_c = np.zeros(cq[('u',0)].shape,'d')

	def initializeGeneralizedInterpolationPointQuadrature(self,t,cip):
		self.cip_c = numpy.zeros((cip['x'].shape[0],cip['x'].shape[1]),'d') 

	def initializeElementBoundaryQuadrature(self,t,cebq,cebq_global):
		ebq_shape = cebq['x'].shape[:-1];         
		ebq_global_shape = cebq_global['x'].shape[:-1]
		self.ebq_c = np.zeros(ebq_shape,'d');  
		self.ebq_global_c = np.zeros(ebq_global_shape,'d')


	def initializeGlobalExteriorElementBoundaryQuadrature(self,t,cebqe):
		self.ebqe_c = np.zeros(cebqe[('u',0)].shape,'d')

	def evaluate(self,t,c):
		if c[('u',0)].shape == self.diff_Model.q[('u',0)].shape:
			mass_frac = self.q_c;
		elif c[('u',0)].shape == self.diff_Model.ebqe[('u',0)].shape:
			mass_frac  = self.ebqe_c;
		elif c[('u',0)].shape == self.diff_Model.ebq[('u',0)].shape:
			mass_frac  = self.ebq_c;
		elif c[('u',0)].shape == self.diff_Model.phi_ip[('u',0)].shape:
			mass_frac  = self.cip_c;
		else:
			raise NotImplementedError

		self.ConMassFluidEvaluate(self.K,
                                  self.grav,
                                  self.L,
                                  c[('u',0)],
                                  c[('m',0)],
                                  c[('dm',0,0)],
                                  c[('f',0)],
                                  c[('df',0,0)],
                                  c[('phi',0)],
                                  c[('dphi',0,0)],
                                  c[('a',0,0)],
                                  c[('r',0)],
                                  c[('x')],
                                  mass_frac)


class Transport_Dilute(LinearVADR_ConstantCoefficients):
	from proteus.dilutetransport.cdilutetransportCoefficients import DiluteEvaluate
	def __init__(self, nc=1, poro = 0, mass=0, diff=0, alpha_L=0, mom_ModelId = None, diff_ModelId = None):
		mass={}
		advection={}
		diffusion={}
		potential={}
		reaction={}
		hamiltonian={}
		for i in range(nc):
			mass[i]      = {i:'linear'}
			advection[i] = {i:'constant'}
			diffusion[i] = {i: {i:'constant'}}
			potential[i] = {i: 'u'}
			reaction[i]  = {i: 'constant'}
		TC_base.__init__(self,nc,mass,advection,diffusion,potential,reaction,hamiltonian)
		self.q={}; self.ebqe={}; self.ebq ={}; self.ebq_global = {}
		self.mass = mass
		self.diff = diff
		self.poro = poro
		self.alpha_L = alpha_L
		self.mom_ModelId  = mom_ModelId
		self.diff_ModelId = diff_ModelId
		self.variableNames=['u']
		self.nd = 1

	def attachModels(self,modelList):
		self.mom_Model = modelList[self.mom_ModelId]
		self.diff_Model = modelList[self.diff_ModelId]
		self.q_vel  = self.mom_Model.q[('velocity',0)]/self.poro
		self.ebqe_vel = self.mom_Model.ebqe[('velocity',0)]/self.poro
		if self.mom_Model.ebq.has_key(('velocity',0)):
			self.ebq_vel = self.mom_Model.ebq[('velocity',0)]/self.poro
		if self.mom_Model.ebq_global.has_key(('velocity',0)):
			self.ebq_global_vel = self.mom_Model.ebq_global[('velocity',0)] /self.poro


	def initializeElementQuadrature(self,t,cq):
		self.q_x_shape = cq['x'].shape
		self.q_vel = numpy.zeros((cq['x'].shape[0],cq['x'].shape[1],self.nd),'d')

	def initializeElementBoundaryQuadrature(self,t,cebq,cebq_global):
		self.ebq_global_x_shape = cebq_global['x'].shape
		self.ebq_x_shape = cebq['x'].shape
		self.ebq_vel = numpy.zeros((cebq['x'].shape[0],cebq['x'].shape[1],cebq['x'].shape[2],self.nd),'d')
		self.ebq_global_vel = numpy.zeros((cebq_global['x'].shape[0],cebq_global['x'].shape[1],self.nd),'d')

	def initializeGlobalExteriorElementBoundaryQuadrature(self,t,cebqe):
		self.ebqe_x_shape = cebqe['x'].shape
		self.ebqe_vel = numpy.zeros((cebqe['x'].shape[0],cebqe['x'].shape[1],self.nd),'d')


	def evaluate(self,t,c): 
		#import pdb; pdb.set_trace()
		if self.q_vel.shape == c[('df',0,0)].shape:
			velocity = self.q_vel
		elif self.ebqe_vel.shape == c[('df',0,0)].shape:
			velocity = self.ebqe_vel
		elif self.ebq_vel.shape == c[('df',0,0)].shape:
			velocity = self.ebq_vel
		else:
			print "no v---------------------"
			raise RuntimeError
		self.DiluteEvaluate(self.poro,
                                    self.diff,
                                    self.alpha_L,
                                    c[('u',0)],
                                    c[('m',0)],
                                    c[('dm',0,0)],
                                    c[('f',0)],
                                    c[('df',0,0)],
                                    c[('a',0,0)],
                                    c[('da',0,0,0)],
                                    velocity)



class Transport_Advection(LinearVADR_ConstantCoefficients):
	from proteus.dilutetransport.cdilutetransportCoefficients import AdvectionEvaluate
	def __init__(self, nc=1, poro = 0, mass=0, diff=0, alpha_L=0, mom_ModelId = None, adv_ModelId = None, diff_ModelId = None):
		mass={}
		advection={}
		diffusion={}
		potential={}
		reaction={}
		hamiltonian={}
		for i in range(nc):
			mass[i]      = {i:'linear'}
			advection[i] = {i:'constant'}
			diffusion[i] = {i: {i:'constant'}}
			potential[i] = {i: 'u'}
			reaction[i]  = {i: 'constant'}
		TC_base.__init__(self,nc,mass,advection,diffusion,potential,reaction,hamiltonian)
		self.q={}; self.ebqe={}; self.ebq ={}; self.ebq_global = {}
		self.mass = mass
		self.diff = diff
		self.poro = poro
		self.alpha_L = alpha_L
		self.mom_ModelId  = mom_ModelId
		self.adv_ModelId  = adv_ModelId
		self.diff_ModelId = diff_ModelId
		self.variableNames=['u_adv']
		self.nd = 1

	def attachModels(self,modelList):
		self.mom_Model = modelList[self.mom_ModelId]
		self.adv_Model = modelList[self.adv_ModelId]
		self.diff_Model = modelList[self.diff_ModelId]
		self.q_vel  = self.mom_Model.q[('velocity',0)]/self.poro
		self.ebqe_vel = self.mom_Model.ebqe[('velocity',0)]/self.poro
		if self.mom_Model.ebq.has_key(('velocity',0)):
			self.ebq_vel = self.mom_Model.ebq[('velocity',0)]/self.poro
		if self.mom_Model.ebq_global.has_key(('velocity',0)):
			self.ebq_global_vel = self.mom_Model.ebq_global[('velocity',0)] /self.poro


	def preStep(self,t,firstStep=False):
		#import pdb; pdb.set_trace()
		self.adv_Model.u[0].dof = self.diff_Model.u[0].dof
		self.adv_Model.timeIntegration.m_last = self.diff_Model.timeIntegration.m_last
		#self.adv_Model.calculateCoefficients()
		#self.adv_Model.calculateElementResidual()
		#self.adv_Model.timeIntegration.updateTimeHistory(resetFromDOF=False)
		#self.adv_Model.timeIntegration.resetTimeHistory(resetFromDOF=True)
		#self.adv_Model.updateTimeHistory(t,resetFromDOF=True)
		copyInstructions = {'copy_uList':True,'uList_model':self.adv_ModelId}
		#copyInstructions = {'reset_uList':True}
		return copyInstructions

	def initializeElementQuadrature(self,t,cq):
		self.q_x_shape = cq['x'].shape
		self.q_vel = numpy.zeros((cq['x'].shape[0],cq['x'].shape[1],self.nd),'d')

	def initializeElementBoundaryQuadrature(self,t,cebq,cebq_global):
		self.ebq_global_x_shape = cebq_global['x'].shape
		self.ebq_x_shape = cebq['x'].shape
		self.ebq_vel = numpy.zeros((cebq['x'].shape[0],cebq['x'].shape[1],cebq['x'].shape[2],self.nd),'d')
		self.ebq_global_vel = numpy.zeros((cebq_global['x'].shape[0],cebq_global['x'].shape[1],self.nd),'d')

	def initializeGlobalExteriorElementBoundaryQuadrature(self,t,cebqe):
		self.ebqe_x_shape = cebqe['x'].shape
		self.ebqe_vel = numpy.zeros((cebqe['x'].shape[0],cebqe['x'].shape[1],self.nd),'d')

	#def postStep(self,t,firstStep=False):
	#	#import pdb; pdb.set_trace()
	#	print self.adv_Model.u[0].dof#pass#self.adv_Model.q['dV_last'][:] = self.adv_Model.q['dV']

	def evaluate(self,t,c): 
		#import pdb; pdb.set_trace()
		if self.q_vel.shape == c[('df',0,0)].shape:
			velocity = self.q_vel
		elif self.ebqe_vel.shape == c[('df',0,0)].shape:
			velocity = self.ebqe_vel
		elif self.ebq_vel.shape == c[('df',0,0)].shape:
			velocity = self.ebq_vel
		else:
			print "no v---------------------"
			raise RuntimeError
		self.AdvectionEvaluate(self.poro,
                                    self.diff,
                                    self.alpha_L,
                                    c[('u',0)],
                                    c[('m',0)],
                                    c[('dm',0,0)],
                                    c[('f',0)],
                                    c[('df',0,0)],
                                    c[('a',0,0)],
                                    c[('da',0,0,0)],
                                    velocity)










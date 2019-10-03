from proteus.TransportCoefficients import *
import numpy as np

class Transport_NonDilute_Dispersion(LinearVADR_ConstantCoefficients):
	from proteus.nondilutetransport.cnondilutetransportCoefficients import NonDiluteDispersionEvaluate
	from proteus.nondilutetransport.cnondilutetransportCoefficients import NonDilutePhiDispersionEvaluate
	def __init__(self, nc = 2, poro = 0, diff=0, alpha_L=0, mom_ModelId = None, diff_ModelId = None, adv_ModelId = None, mom_ModelId2 = None, adv_ModelId2 = None, ModelId = None, Alt_Split = False, Mom_Split = False, R = 0, theta = 0, MW_a = 0, MW_b = 0, beta1 = 0, beta2 = 0):
		mass={}
		advection={}
		diffusion={}
		potential={}
		reaction={}
		hamiltonian={}
		mass = {0:{0:'nonlinear'},
                1:{1:'constant'}}
		advection = {0:{0:'constant'},
                     1:{1:'constant'}}
		diffusion = {0:{0:{0:'nonlinear',1:'nonlinear'},
                        1:{0:'nonlinear',1:'nonlinear'}},
                     1:{0:{0:'nonlinear',1:'nonlinear'},
                        1:{0:'nonlinear',1:'nonlinear'}}}
		potential = {0:{0:'nonlinear'},
                     1:{0:'nonlinear',
                        1:'nonlinear'}}
		reaction = {0:{0:'constant'},
                    1:{0:'nonlinear',1:'nonlinear'}}
		TC_base.__init__(self,nc,mass,advection,diffusion,potential,reaction,hamiltonian,useSparseDiffusion = True)
		self.q={}; self.ebqe={}; self.ebq ={}; self.ebq_global = {}; self.phi_ip = {}
		self.q_mf_old = {}; self.ebqe_mf_old = {}; self.cip_mf_old = {}; self.ebq_mf_old = {}; self.ebq_global_mf_old = {}
		self.dof_mf_old = {}
		self.diff = diff
		self.poro = poro
		self.alpha_L = alpha_L
		self.R = R
		self.theta = theta
		self.MW_a = MW_a
		self.MW_b = MW_b
		self.beta1 = beta1
		self.beta2 = beta2
		self.mom_ModelId  = mom_ModelId
		self.diff_ModelId = diff_ModelId
		self.adv_ModelId = adv_ModelId
		self.Alt_Split = Alt_Split
		self.Mom_Split = Mom_Split
		if self.Alt_Split:
			if self.Mom_Split:
				self.mom_ModelId2 = mom_ModelId2
			self.adv_ModelId2 = adv_ModelId2
		self.variableNames=['u','act']
		self.nd = 1
		self.nc = nc

	def attachModels(self,modelList):
		self.mom_Model = modelList[self.mom_ModelId]
		self.diff_Model = modelList[self.diff_ModelId]
		self.adv_Model = modelList[self.adv_ModelId]
		if self.Alt_Split:
			self.adv_Model2 = modelList[self.adv_ModelId2]
			if self.Mom_Split:
				self.mom_Model2 = modelList[self.mom_ModelId2]
		self.q_p  = self.mom_Model.q[('u',0)]
		self.ebqe_p = self.mom_Model.ebqe[('u',0)]
		if self.mom_Model.ebq.has_key(('u',0)):
			self.ebq_p = self.mom_Model.ebq[('u',0)]
		if self.mom_Model.ebq_global.has_key(('u',0)):
			self.ebq_global_p = self.mom_Model.ebq_global[('u',0)]
		self.diff_Model.q_mf_old = np.copy(self.diff_Model.q[('u',0)])
		self.diff_Model.q_mf = np.copy(self.diff_Model.q[('u',0)])
		self.diff_Model.cip_mf_old = np.copy(self.diff_Model.phi_ip[('u',0)])
		self.diff_Model.cip_mf = np.copy(self.diff_Model.phi_ip[('u',0)])
		self.diff_Model.ebqe_mf_old = np.copy(self.diff_Model.ebqe[('u',0)])
		self.diff_Model.ebqe_mf = np.copy(self.diff_Model.ebqe[('u',0)])
		if self.diff_Model.ebq.has_key(('u',0)):
			self.diff_Model.ebq_mf_old = np.copy(self.diff_Model.ebq[('u',0)])
			self.diff_Model.ebq_mf = np.copy(self.diff_Model.ebq[('u',0)])
		if self.diff_Model.ebq_global.has_key(('u',0)):
			self.diff_Model.ebq_global_mf_old = np.copy(self.diff_Model.ebq_global[('u',0)])
			self.diff_Model.ebq_global_mf = np.copy(self.diff_Model.ebq_global[('u',0)])


	def preStep(self,t,firstStep=True):
		if self.Alt_Split:
			#import pdb; pdb.set_trace()
			self.diff_Model.q_mf_old = np.copy(self.diff_Model.q[('u',0)])
			self.diff_Model.cip_mf_old = np.copy(self.diff_Model.phi_ip[('u',0)])
			self.diff_Model.ebqe_mf_old = np.copy(self.diff_Model.ebqe[('u',0)])
			if self.diff_Model.ebq.has_key(('u',0)):
				self.diff_Model.ebq_mf_old = np.copy(self.diff_Model.ebq[('u',0)])
			if self.diff_Model.ebq_global.has_key(('u',0)):
				self.diff_Model.ebq_global_mf_old = np.copy(self.diff_Model.ebq_globa[('u',0)])
		else:
			self.diff_Model.q_mf_old = np.copy(self.diff_Model.q[('u',0)])
			self.diff_Model.cip_mf_old = np.copy(self.diff_Model.phi_ip[('u',0)])
			self.diff_Model.ebqe_mf_old = np.copy(self.diff_Model.ebqe[('u',0)])
			if self.diff_Model.ebq.has_key(('u',0)):
				self.diff_Model.ebq_mf_old = np.copy(self.diff_Model.ebq[('u',0)])
			if self.diff_Model.ebq_global.has_key(('u',0)):
				self.diff_Model.ebq_global_mf_old = np.copy(self.diff_Model.ebq_globa[('u',0)])

		self.diff_Model.u[0].dof = np.copy(self.adv_Model.u[0].dof)
 		self.diff_Model.calculateCoefficients()
 		self.diff_Model.calculateElementResidual()
 		self.diff_Model.timeIntegration.updateTimeHistory(resetFromDOF=True)
 		self.diff_Model.timeIntegration.resetTimeHistory(resetFromDOF=True)
 		self.diff_Model.updateTimeHistory(t,resetFromDOF=True)
 		copyInstructions = {'copy_uList':True,'uList_model':self.diff_ModelId}
 		copyInstructions = {'reset_uList':True}
 		return copyInstructions

	def postStep(self,t,firstStep=True):
		self.diff_Model.q_mf = np.copy(self.diff_Model.q[('u',0)])
		self.diff_Model.cip_mf = np.copy(self.diff_Model.phi_ip[('u',0)])
		self.diff_Model.ebqe_mf = np.copy(self.diff_Model.ebqe[('u',0)])
		if self.diff_Model.ebq.has_key(('u',0)):
			self.diff_Model.ebq_mf = np.copy(self.diff_Model.ebq[('u',0)])
		if self.diff_Model.ebq_global.has_key(('u',0)):
			self.diff_Model.ebq_global_mf = np.copy(self.diff_Model.ebq_globa[('u',0)])



	def initializeElementQuadrature(self,t,cq):
		self.q_x_shape = cq['x'].shape
		self.q_vel = numpy.zeros((cq['x'].shape[0],cq['x'].shape[1],self.nd),'d')

	def initializeElementBoundaryQuadrature(self,t,cebq,cebq_global):
		self.ebq_global_x_shape = cebq_global['x'].shape
		self.ebq_x_shape = cebq['x'].shape
		self.ebq_vel = numpy.zeros((cebq['x'].shape[0],cebq['x'].shape[1],cebq['x'].shape[2],self.nd),'d')
		self.ebq_global_vel = numpy.zeros((cebq_global['x'].shape[0],cebq_global['x'].shape[1],self.nd),'d')

	def initializeGeneralizedInterpolationPointQuadrature(self,t,cip):
		self.cip_vel = numpy.zeros((cip['x'].shape[0],cip['x'].shape[1],self.nd),'d')

	def initializeGlobalExteriorElementBoundaryQuadrature(self,t,cebqe):
		self.ebqe_x_shape = cebqe['x'].shape
		self.ebqe_vel = numpy.zeros((cebqe['x'].shape[0],cebqe['x'].shape[1],self.nd),'d')

	def evaluate(self,t,c):
		phi_eval = False
		if self.q_vel.shape == c[('df',0,0)].shape:
			velocity = self.mom_Model.q[('velocity',0)]/self.poro
			pressure = self.q_p;
			w_old = np.zeros_like(c[('u',0)])
			if len(self.diff_Model.q_mf_old) < 1:
				w_old = c[('u',0)]
			else:
				w_old = self.diff_Model.q_mf_old
		elif self.ebqe_vel.shape == c[('df',0,0)].shape:
			velocity = self.mom_Model.ebqe[('velocity',0)]/self.poro
			pressure = self.ebqe_p
			w_old = np.zeros_like(c[('u',0)])
			if len(self.diff_Model.ebqe_mf_old) < 1:
				w_old = c[('u',0)]
			else:
				w_old = self.diff_Model.ebqe_mf_old
		elif self.ebq_vel.shape == c[('df',0,0)].shape:
			velocity = self.mom_Model.ebq[('velocity',0)]/self.poro
			pressure = self.ebq_p
			w_old = np.zeros_like(c[('u',0)])
			if len(self.diff_Model.ebq_mf_old) < 1:
				w_old = c[('u',0)]
			else:
				w_old = self.diff_Model.ebq_mf_old
		elif self.cip_vel.shape == c[('df',0,0)].shape:
			phi_eval = True
		else:
			print "no v---------------------"
			raise RuntimeError
		if phi_eval:
			self.NonDilutePhiDispersionEvaluate(self.poro,
                               self.diff,
                               self.alpha_L,
                               self.R,
                               self.theta,
                               self.MW_a,
                               self.MW_b,
                               self.beta1,
                               self.beta2,
                               c[('u',0)],
                               c[('u',1)],
                               c[('phi',0)],
                               c[('dphi',0,0)],
                               c[('dphi',0,1)],
                               c[('phi',1)],
                               c[('dphi',1,0)],
                               c[('dphi',1,1)],
                               c[('f',0)])
		else:
			self.NonDiluteDispersionEvaluate(self.poro,
                               self.diff,
                               self.alpha_L,
                               self.R,
                               self.theta,
                               self.MW_a,
                               self.MW_b,
                               self.beta1,
                               self.beta2,
                               c[('u',0)],
                               c[('u',1)],
                               c[('m',0)],
                               c[('dm',0,0)],
                               c[('f',0)],
                               c[('df',0,0)],
                               c[('r',1)],
                               c[('dr',1,0)],
                               c[('dr',1,1)],
                               c[('a',0,0)],
                               c[('da',0,0,0)],
                               c[('a',0,1)],
                               c[('da',0,1,0)],
                               c[('da',0,1,1)],
                               velocity,
                               pressure,
                               w_old,
                               c[('grad(u)',0)],
                               c[('grad(u)',1)],
                               c[('phi',0)],
                               c[('dphi',0,0)],
                               c[('dphi',0,1)],
                               c[('phi',1)],
                               c[('dphi',1,0)],
                               c[('dphi',1,1)],
                               c[('x')])


class Porous_Media_Flow(TC_base):
	from proteus.nondilutetransport.cnondilutetransportCoefficients import ConMassFluidEvaluate
	def __init__(self,nc=1,nd=1,diff_ModelId = None, adv_ModelId = None, mom_ModelId = None, mom_ModelId2 = None, adv_ModelId2 = None, ModelId = None, Alt_Split = False, Mom_Split = False, K = 0.0 , grav = -980., poro = 0.0, L = 0.0):
		mass={}
		advection={}
		diffusion={}
		potential={}
		reaction={}
		hamiltonian={}
		self.nd=nd
		mass = {0:{0:'constant'}}
		advection = {0:{0:'constant'}}
		diffusion = {0:{0:{0:'constant'}}}
		potential = {0:{0:'u'}}
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
		self.Alt_Split = Alt_Split
		self.Mom_Split = Mom_Split
		if self.Alt_Split:
			if self.Mom_Split:
				self.mom_ModelId2 = mom_ModelId2
			self.adv_ModelId2 = adv_ModelId2
		self.ModelId = ModelId
		self.K = K
		self.grav = grav
		self.poro = poro
		self.L = L
		self.variableNames=['pressure']
		self.dt = 0.0
		self.nc = nc

	def attachModels(self,modelList):
		self.mom_Model = modelList[self.mom_ModelId]
		if (self.diff_ModelId is not None):
			self.diff_Model = modelList[self.diff_ModelId]
		if (self.adv_ModelId is not None):
				self.adv_Model = modelList[self.adv_ModelId]
		if self.Alt_Split:
			self.adv_Model2 = modelList[self.adv_ModelId2]
			if self.Mom_Split:
				self.mom_Model2 = modelList[self.mom_ModelId2]

		if self.Alt_Split:
			if self.ModelId == self.mom_ModelId:
				self.oldTime_Model = self.adv_Model2
				self.currentTime_Model = self.adv_Model2
			elif self.ModelId == self.mom_ModelId2:
				self.oldTime_Model = self.diff_Model
				self.currentTime_Model = self.diff_Model
		else:
			self.oldTime_Model = self.diff_Model
			self.currentTime_Model = self.diff_Model

		self.q_c  = self.currentTime_Model.q[('u',0)]
		self.ebqe_c = self.currentTime_Model.ebqe[('u',0)]
		if self.currentTime_Model.phi_ip.has_key(('u',0)):
			self.cip_c = self.currentTime_Model.phi_ip[('u',0)]
		if self.currentTime_Model.ebq.has_key(('u',0)):
			self.ebq_c = self.currentTime_Model.ebq[('u',0)]
		if self.currentTime_Model.ebq_global.has_key(('u',0)):
			self.ebq_global_c = self.currentTime_Model.ebq_global[('u',0)]

	def preStep(self,t,firstStep=True):
		self.dt = self.currentTime_Model.timeIntegration.dt

	# def postStep(self,t,firstStep=True):
	# 	print self.ModelId
	# 	print np.max(self.mom_Model.q['velocity',0]-self.mom_Model2.q['velocity',0])

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
		#import pdb; pdb.set_trace()
		if c[('u',0)].shape == self.currentTime_Model.q[('u',0)].shape:
			mass_frac = self.currentTime_Model.q_mf;
			if ( not hasattr(self.oldTime_Model,'q_mf_old')  or len(self.oldTime_Model.q_mf_old) < 1):
				mass_frac_old = mass_frac
 			else:
				mass_frac_old = self.oldTime_Model.q_mf_old
				#import pdb; pdb.set_trace()
		elif c[('u',0)].shape == self.currentTime_Model.ebqe[('u',0)].shape:
			mass_frac  = self.currentTime_Model.ebqe_mf;
			if (not hasattr(self.oldTime_Model,'ebqe_mf_old') or len(self.oldTime_Model.ebqe_mf_old) < 1):
				mass_frac_old = mass_frac
 			else:
				mass_frac_old  = self.oldTime_Model.ebqe_mf_old
		elif c[('u',0)].shape == self.currentTime_Model.ebq[('u',0)].shape:
			mass_frac  = self.currentTime_Model.ebq_mf;
			if (not hasattr(self.oldTime_Model,'ebq_mf_old') or len(self.oldTime_Model.ebq_mf_old) < 1):
				mass_frac_old = mass_frac
 			else:
				mass_frac_old  = self.oldTime_Model.ebq_mf_old
		elif c[('u',0)].shape == self.currentTime_Model.phi_ip[('u',0)].shape:
			mass_frac  = self.currentTime_Model.cip_mf;
			if (not hasattr(self.oldTime_Model,'cip_mf_old') or len(self.oldTime_Model.cip_mf_old) < 1):
				mass_frac_old = mass_frac
 			else:
				mass_frac_old  = self.oldTime_Model.cip_mf_old
		else:
			raise NotImplementedError
		elements = mass_frac.shape[0]
		points = mass_frac.shape[1]
		dof_length = self.currentTime_Model.u[0].dof.shape[0]
		self.ConMassFluidEvaluate(elements,points,dof_length,
		                          self.currentTime_Model.u[0].dof,
		                          self.K,
                                  self.grav,
                                  self.dt,
                                  self.poro,
                                  self.L,
                                  c[('u',0)],
                                  c[('m',0)],
                                  c[('dm',0,0)],
                                  c[('f',0)],
                                  c[('df',0,0)],
                                  c[('phi',0)],
                                  c[('dphi',0,0)],
                                  c[('a',0,0)],
                                  c[('da',0,0,0)],
                                  c[('r',0)],
                                  c[('x')],
                                  mass_frac,
                                  mass_frac_old)

class Transport_NonDilute(LinearVADR_ConstantCoefficients):
	from proteus.nondilutetransport.cnondilutetransportCoefficients import NonDiluteEvaluate
	from proteus.nondilutetransport.cnondilutetransportCoefficients import NonDilutePhiEvaluate
	def __init__(self, nc = 2, poro = 0, diff=0, alpha_L=0, mom_ModelId = None, diff_ModelId = None, R = 0, theta = 0, MW_a = 0, MW_b = 0, beta1 = 0, beta2 = 0):
		mass={}
		advection={}
		diffusion={}
		potential={}
		reaction={}
		hamiltonian={}
		mass = {0:{0:'nonlinear'},
                1:{1:'linear'}}
		advection = {0:{0:'nonlinear'},
                     1:{1:'constant'}}
		diffusion = {0:{0:{0:'nonlinear',1:'constant'},
                        1:{0:'nonlinear',1:'nonlinear'}},
                     1:{0:{0:'constant',1:'constant'},
                        1:{0:'constant',1:'constant'}}}
		potential = {0:{0:'u'},
                     1:{0:'nonlinear',
                        1:'nonlinear'}}
		reaction = {0:{0:'constant'},
                    1:{0:'nonlinear',1:'nonlinear'}}
		TC_base.__init__(self,nc,mass,advection,diffusion,potential,reaction,hamiltonian,useSparseDiffusion = True)
		self.q={}; self.ebqe={}; self.ebq ={}; self.ebq_global = {}; self.phi_ip = {}
		self.q_mf_old = {}; self.ebqe_mf_old = {}; self.cip_mf_old = {}; self.ebq_mf_old = {}; self.ebq_global_mf_old = {}
		self.diff = diff
		self.poro = poro
		self.alpha_L = alpha_L
		self.R = R
		self.theta = theta
		self.MW_a = MW_a
		self.MW_b = MW_b
		self.beta1 = beta1
		self.beta2 = beta2
		self.mom_ModelId  = mom_ModelId
		self.diff_ModelId = diff_ModelId
		self.variableNames=['u','act']
		self.nd = 1
		self.nc = nc

	def attachModels(self,modelList):
		self.mom_Model = modelList[self.mom_ModelId]
		self.diff_Model = modelList[self.diff_ModelId]
		self.q_p  = self.mom_Model.q[('grad(u)',0)]
		self.ebqe_p = self.mom_Model.ebqe[('grad(u)',0)]
		if self.mom_Model.ebq.has_key(('grad(u)',0)):
			self.ebq_p = self.mom_Model.ebq[('grad(u)',0)]
		if self.mom_Model.ebq_global.has_key(('grad(u)',0)):
			self.ebq_global_p = self.mom_Model.ebq_global[('grad(u)',0)]
		self.diff_Model.q_mf_old = np.copy(self.diff_Model.q[('u',0)])
		self.diff_Model.q_mf = np.copy(self.diff_Model.q[('u',0)])
		self.diff_Model.cip_mf_old = np.copy(self.diff_Model.phi_ip[('u',0)])
		self.diff_Model.cip_mf = np.copy(self.diff_Model.phi_ip[('u',0)])
		self.diff_Model.ebqe_mf_old = np.copy(self.diff_Model.ebqe[('u',0)])
		self.diff_Model.ebqe_mf = np.copy(self.diff_Model.ebqe[('u',0)])
		if self.diff_Model.ebq.has_key(('u',0)):
			self.diff_Model.ebq_mf_old = np.copy(self.diff_Model.ebq[('u',0)])
			self.diff_Model.ebq_mf = np.copy(self.diff_Model.ebq[('u',0)])
		if self.diff_Model.ebq_global.has_key(('u',0)):
			self.diff_Model.ebq_global_mf_old = np.copy(self.diff_Model.ebq_global[('u',0)])
			self.diff_Model.ebq_global_mf = np.copy(self.diff_Model.ebq_global[('u',0)])

	def preStep(self,t,firstStep=True):
		self.diff_Model.q_mf_old  = np.copy(self.diff_Model.q[('u',0)])
		self.diff_Model.ebqe_mf_old = np.copy(self.diff_Model.ebqe[('u',0)])
		if self.diff_Model.phi_ip.has_key(('u',0)):
			self.diff_Model.cip_mf_old = np.copy(self.diff_Model.phi_ip[('u',0)])
		if self.diff_Model.ebq.has_key(('u',0)):
			self.diff_Model.ebq_mf_old = np.copy(self.diff_Model.ebq[('u',0)])
		if self.diff_Model.ebq_global.has_key(('u',0)):
			self.diff_Model.ebq_global_mf_old = np.copy(self.diff_Model.ebq_global[('u',0)] )

	def postStep(self,t,firstStep=True):
		self.diff_Model.q_mf = np.copy(self.diff_Model.q[('u',0)])
		self.diff_Model.cip_mf = np.copy(self.diff_Model.phi_ip[('u',0)])
		self.diff_Model.ebqe_mf = np.copy(self.diff_Model.ebqe[('u',0)])
		if self.diff_Model.ebq.has_key(('u',0)):
			self.diff_Model.ebq_mf = np.copy(self.diff_Model.ebq[('u',0)])
		if self.diff_Model.ebq_global.has_key(('u',0)):
			self.diff_Model.ebq_global_mf = np.copy(self.diff_Model.ebq_globa[('u',0)])

	def initializeElementQuadrature(self,t,cq):
		self.q_x_shape = cq['x'].shape
		self.q_vel = numpy.zeros((cq['x'].shape[0],cq['x'].shape[1],self.nd),'d')

	def initializeElementBoundaryQuadrature(self,t,cebq,cebq_global):
		self.ebq_global_x_shape = cebq_global['x'].shape
		self.ebq_x_shape = cebq['x'].shape
		self.ebq_vel = numpy.zeros((cebq['x'].shape[0],cebq['x'].shape[1],cebq['x'].shape[2],self.nd),'d')
		self.ebq_global_vel = numpy.zeros((cebq_global['x'].shape[0],cebq_global['x'].shape[1],self.nd),'d')

	def initializeGeneralizedInterpolationPointQuadrature(self,t,cip):
		self.cip_vel = numpy.zeros((cip['x'].shape[0],cip['x'].shape[1],self.nd),'d')

	def initializeGlobalExteriorElementBoundaryQuadrature(self,t,cebqe):
		self.ebqe_x_shape = cebqe['x'].shape
		self.ebqe_vel = numpy.zeros((cebqe['x'].shape[0],cebqe['x'].shape[1],self.nd),'d')

	def evaluate(self,t,c):
		phi_eval = False;
		if self.q_vel.shape == c[('df',0,0)].shape:
			velocity = self.mom_Model.q[('velocity',0)]/self.poro; #import pdb; pdb.set_trace()
			pressure = self.q_p;
		elif self.ebqe_vel.shape == c[('df',0,0)].shape:
			velocity = self.mom_Model.ebqe[('velocity',0)]/self.poro
			pressure = self.ebqe_p
		elif self.ebq_vel.shape == c[('df',0,0)].shape:
			velocity = self.mom_Model.ebq[('velocity',0)]/self.poro
			pressure = self.ebq_p
		elif self.cip_vel.shape == c[('df',0,0)].shape:
			velocity = self.cip_vel
			phi_eval = True
		else:
			print "no v---------------------"
			raise RuntimeError
		if phi_eval:
			self.NonDilutePhiEvaluate(self.poro,
                               self.diff,
                               self.alpha_L,
                               self.R,
                               self.theta,
                               self.MW_a,
                               self.MW_b,
                               self.beta1,
                               self.beta2,
                               c[('u',0)],
                               c[('u',1)],
                               c[('phi',0)],
                               c[('dphi',0,0)],
                               c[('dphi',0,1)],
                               c[('phi',1)],
                               c[('dphi',1,0)],
                               c[('dphi',1,1)],
                               c[('f',0)])
		else:
			self.NonDiluteEvaluate(self.poro,
                               self.diff,
                               self.alpha_L,
                               self.R,
                               self.theta,
                               self.MW_a,
                               self.MW_b,
                               self.beta1,
                               self.beta2,
                               c[('u',0)],
                               c[('u',1)],
                               c[('m',0)],
                               c[('dm',0,0)],
                               c[('f',0)],
                               c[('df',0,0)],
                               c[('r',1)],
                               c[('dr',1,0)],
                               c[('dr',1,1)],
                               c[('a',0,0)],
                               c[('da',0,0,0)],
                               c[('a',0,1)],
                               c[('da',0,1,0)],
                               c[('da',0,1,1)],
                               velocity,
                               pressure,
                               c[('grad(u)',0)],
                               c[('grad(u)',1)],
                               c[('phi',0)],
                               c[('dphi',0,0)],
                               c[('dphi',0,1)],
                               c[('phi',1)],
                               c[('dphi',1,0)],
                               c[('dphi',1,1)],
                               c[('x')])



class Transport_Advection(LinearVADR_ConstantCoefficients):
	from proteus.nondilutetransport.cnondilutetransportCoefficients import AdvectionEvaluate
	def __init__(self, nc=1, poro = 0, diff=0, alpha_L=0, mom_ModelId = None, adv_ModelId = None, diff_ModelId = None):
		mass={}
		advection={}
		diffusion={}
		potential={}
		reaction={}
		hamiltonian={}
		for i in range(nc):
			mass[i]      = {i:'nonlinear'}
			advection[i] = {i:'nonlinear'}
			diffusion[i] = {i: {i:'constant'}}
			potential[i] = {i: 'u'}
			reaction[i]  = {i: 'constant'}
		TC_base.__init__(self,nc,mass,advection,diffusion,potential,reaction,hamiltonian)
		self.q={}; self.ebqe={}; self.ebq ={}; self.ebq_global = {}
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
		self.adv_Model.u[0].dof = self.diff_Model.u[0].dof
		#self.adv_Model.timeIntegration.m_last[0] = self.diff_Model.timeIntegration.m_last[0]
		copyInstructions = {'copy_uList':True,'uList_model':self.adv_ModelId}
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


	def evaluate(self,t,c):
		#import pdb; pdb.set_trace()
		if self.q_vel.shape == c[('df',0,0)].shape:
			velocity = self.mom_Model.q[('velocity',0)]/self.poro
		elif self.ebqe_vel.shape == c[('df',0,0)].shape:
			velocity = self.mom_Model.ebqe[('velocity',0)]/self.poro
		elif self.ebq_vel.shape == c[('df',0,0)].shape:
			velocity = self.mom_Model.ebq[('velocity',0)]/self.poro
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


### Fully Coupled Model - Error with boundary conditions
class Transport_NonDilute_TEST(LinearVADR_ConstantCoefficients):
	from proteus.nondilutetransport.cnondilutetransportCoefficients import NonDiluteTESTEvaluate
	from proteus.nondilutetransport.cnondilutetransportCoefficients import NonDilutePhiTESTEvaluate
	def __init__(self, nc = 2, poro = 0, grav= 0, K = 0, L = 0, diff=0, alpha_L=0,  R = 0, theta = 0, MW_a = 0, MW_b = 0, beta1 = 0, beta2 = 0):
		mass={}
		advection={}
		diffusion={}
		potential={}
		reaction={}
		hamiltonian={}
		mass = {0:{0:'constant',
                   1:'nonlinear'},
                1:{1:'nonlinear'}}
		advection = {0:{0:'constant'},
                     1:{0:'nonlinear',1:'nonlinear'}}
		diffusion = {0:{0:{0:'constant',
                           1:'nonlinear'},
                        1:{0:'constant',
                           1:'constant'}},
                     1:{0:{0:'constant',
                           1:'nonlinear'},
                        1:{0:'constant',
                           1:'nonlinear'}}}
		potential = {0:{0:'nonlinear',1:'nonlinear'},
                     1:{1:'nonlinear'}}
		reaction = {0:{0:'constant'},
                    1:{1:'constant'}}
		TC_base.__init__(self,nc,mass,advection,diffusion,potential,reaction,hamiltonian,useSparseDiffusion = True)
		self.poro = poro
		self.grav = grav
		self.K = K
		self.L = L
		self.diff = diff
		self.alpha_L = alpha_L
		self.R = R
		self.theta = theta
		self.MW_a = MW_a
		self.MW_b = MW_b
		self.beta1 = beta1
		self.beta2 = beta2
		self.variableNames=['press','w']
		self.nd = 1
		self.nc = 2
		self.upwind_potential = {0:[0],1:[0]}

	def initializeElementQuadrature(self,t,cq):
		self.q_x_shape = cq['x'].shape
		self.q_vel = numpy.zeros((cq['x'].shape[0],cq['x'].shape[1],self.nd),'d')

	def initializeElementBoundaryQuadrature(self,t,cebq,cebq_global):
		self.ebq_global_x_shape = cebq_global['x'].shape
		self.ebq_x_shape = cebq['x'].shape
		self.ebq_vel = numpy.zeros((cebq['x'].shape[0],cebq['x'].shape[1],cebq['x'].shape[2],self.nd),'d')
		self.ebq_global_vel = numpy.zeros((cebq_global['x'].shape[0],cebq_global['x'].shape[1],self.nd),'d')

	def initializeGeneralizedInterpolationPointQuadrature(self,t,cip):
		self.cip_vel = numpy.zeros((cip['x'].shape[0],cip['x'].shape[1],self.nd),'d')

	def initializeGlobalExteriorElementBoundaryQuadrature(self,t,cebqe):
		self.ebqe_x_shape = cebqe['x'].shape
		self.ebqe_vel = numpy.zeros((cebqe['x'].shape[0],cebqe['x'].shape[1],self.nd),'d')



	def evaluate(self,t,c):
		#import pdb; pdb.set_trace()
		phi_eval = False
		if self.q_vel.shape == c[('df',0,0)].shape:
			pass
		elif self.ebqe_vel.shape == c[('df',0,0)].shape:
			pass
		elif self.ebq_vel.shape == c[('df',0,0)].shape:
			pass
		elif self.cip_vel.shape == c[('df',0,0)].shape:
			phi_eval = True;
		else:
			print "no v---------------------"
			raise RuntimeError
		if not phi_eval:
			c[('upwind_velocity',0)] = numpy.zeros(c[('grad(u)',0)].shape,'d')
			self.NonDiluteTESTEvaluate(self.poro,
                               self.grav,
                               self.K,
                               self.L,
                               self.diff,
                               self.alpha_L,
                               self.R,
                               self.theta,
                               self.MW_a,
                               self.MW_b,
                               self.beta1,
                               self.beta2,
                               c[('u',0)],
                               c[('u',1)],
                               c[('m',0)],
                               c[('m',1)],
                               c[('dm',0,1)],
                               c[('dm',1,1)],
                               c[('phi',0)],
                               c[('phi',1)],
                               c[('dphi',0,0)],
                               c[('dphi',0,1)],
                               c[('dphi',1,0)],
                               c[('dphi',1,1)],
                               c[('f',1)],
                               c[('df',1,0)],
                               c[('df',1,1)],
                               c[('a',0,0)],
                               c[('a',1,0)],
                               c[('a',1,1)],
                               c[('da',0,0,1)],
                               c[('da',1,0,1)],
                               c[('da',1,1,1)],
                               c[('x')],
                               c[('velocity',0)])
			space_dim = c[('f',0)].shape[-1]
			for k in range(len(c[('u',0)].flat)):
				c[('upwind_velocity',0)].flat[k*space_dim] = -c[('dphi',0,0)].flat[k]*c[('grad(u)',0)].flat[k*space_dim]-c[('dphi',0,1)].flat[k]*c[('grad(u)',1)].flat[k*space_dim]
			c[('upwind_velocity',1)] = c[('upwind_velocity',0)]



		else:
			self.NonDilutePhiTESTEvaluate(self.poro,
                               self.grav,
                               self.L,
                               c[('u',0)],
                               c[('u',1)],
                               c[('phi',0)],
                               c[('phi',1)],
                               c[('dphi',0,0)],
                               c[('dphi',0,1)],
                               c[('dphi',1,0)],
                               c[('dphi',1,1)],
                               c[('f',0)],
                               c[('x')])

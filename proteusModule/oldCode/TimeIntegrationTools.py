from LinearAlgebra import *
import sys
## \defgroup TimeIntegrationTools TimeIntegrationTools
#
# @{

"""
Tools for approximating time-dependent PDE's
"""

class TimeIntegration:
    """
    The base class for time integration methods of
    the scalar transport equation. TimeIntegration
    objects generate the following multistep approximations

    m_t = a_{m,0} m(u_n) + \sum_{k=1} a_{m,k} m_{n-k}
    u_* = a_{u,0) u_n    + \sum_{k=1} a_{u,k} u_{n-k}
    f_* = a_{f,0} f(u_*) + \sum_{k=1} a_{f,k} f_{n-k}
    ...
    
    The base class implements no time integration
    method and results in the steady state scalar
    transport equation:

    \deld ( f - a \grad phi) + r = 0
    """
    def __init__(self,scalarTransport):
        """
        Set flags that indicate that all terms
        are implicit.
        """
        self.u = None
        self.duStar_du = None
        self.DT = 1.0
        self.massIsImplicit = True
        self.advectionIsImplicit = True
        self.diffusionIsImplicit = True
        self.reactionIsImplicit = True
        self.stabilizationIsImplicit = True
        self.shockCapturingIsImplicit = True
    def calculateU(self,u):
        """
        Generate u_*
        """
        self.u = u
    def updateMass(self,m,mt,dm,dmt):
        """
        Calculate m_t, and dm_t/du given m,dm.
        """
        pass
    def updateAdvection(self,f,df):
        """
        Given f,df, update with f^*,df^*/du,...
        """
        pass
    def updateDivF(self,div_f,ddiv_f):
        """
        Given div_f,ddiv_f, update with div_f^*,ddiv_f^*/du,...
        """
        pass
    def updateDivA(self,div_a,ddiv_a):
        """
        Given div_a,ddiv_a, update with div_a^*,ddiv_a^*/du,...
        """
        pass
    def updateAdjoint(self,Lstar_w):
        pass
    def updateDiffusion(self,a,da):
        """
        Given a,da, update with a^*,da^*/du,...
        """
        pass
    def updateGradients(self,dphi,grad_phi_x_grad_w,grad_u_x_grad_w):
        """
        Given grad_u_x_grad_w,grad_phi_x_grad_w update  with grad_u_x_grad_w*, grad_phi_x_grad_w*
        """
        pass
    def updateReaction(self,r,dr):
        """
        Given r,dr, update with r^*,dr^*/du,...
        """
        pass
    def updateStabilization(self,tau):
        """
        Given \tau, update with \tau_*
        """
        pass
    def updateShockCapturing(self,numDiff):
        """
        Given numDiff, update with numDiff_*
        """
        pass
    def updatePDEResidual(self,pdeResidual,dpdeResidual):
        """
        Given R,DR, update with R^*,dR^*/du
        """
        pass
    def chooseDT(self):
        """
        Modify self.DT
        """
        pass
    def updateTimeHistory(self):
        """
        Push necessary information into time history arrays
        """
        pass
    #mwf add
    def updateStage(self):
        """
        increment stage counters and push necessary information
        into stage arrays
        """
        pass
    #mwf add
    def setInitialStageValues(self):
        """
        set the all stage values assuming this is the first step
        after a problem reset
        """
        pass
    #mwf add
    def updateNumericalFlux(self,f,d):
        """
        Given advective and diffusive numerical fluxes and derivative update
        """
        pass
    def updateNumericalFluxJacobian(self,df):
        """
        Given advective and diffusive numerical fluxes and derivative update
        """
        pass
    
NoIntegration = TimeIntegration

class BackwardEuler:
    def __init__(self,scalarTransport):
        self.runCFL=0.1
        self.u = None
        self.duStar_du = None
        self.scalarTransport = scalarTransport
        self.cfl = scalarTransport.q['cfl']
        self.DT = 1.0
        self.massIsImplicit = True
        self.advectionIsImplicit = True
        self.diffusionIsImplicit = True
        self.reactionIsImplicit = True
        self.stabilizationIsImplicit = True
        self.shockCapturingIsImplicit = True
        self.m_last = Numeric.array(scalarTransport.q['m'])
        if scalarTransport.q.has_key('tau'):
            self.tau_last = Numeric.array(scalarTransport.q['tau'])
            self.tau_tmp = Numeric.array(scalarTransport.q['tau'])
        else:
            self.tau_last = None
            self.tau_tmp = None
        if scalarTransport.q.has_key('numDiff'):
            self.numDiff_last = Numeric.array(scalarTransport.q['numDiff'])
            self.numDiff_tmp = Numeric.array(scalarTransport.q['numDiff'])
        else:
            self.numDiff_last = None
            self.numDiff_tmp = None
    def calculateU(self,u):
        self.u = u
    def updateMass(self,m,mt,dm,dmt):
        self.m_tmp = m
        mt.flat[:]=m.flat
        mt-=self.m_last
        mt/=self.DT
        dmt.flat[:]=dm.flat
        dmt/=self.DT
    def updateAdvection(self,f,df):
        pass
    def updateDivF(self,div_f,ddiv_f):
        pass
    def updateDivA(self,div_a,ddiv_a):
        pass
    def updateAdjoint(self,Lstar_w):
        pass
    def updateDiffusion(self,a,da):
        pass
    def updateGradients(self,dphi,grad_phi_x_grad_w,grad_u_x_grad_w):
        pass
    def updateReaction(self,r,dr):
        pass
    def updateStabilization(self,tau):
        self.tau_tmp.flat[:] = tau.flat
        tau.flat[:] = self.tau_last.flat
    def updateShockCapturing(self,numDiff):
        self.numDiff_tmp.flat[:] = numDiff.flat
        numDiff.flat[:] = self.numDiff_last.flat
    def updatePDEResidual(self,pdeResidual,dpdeResidual):
        pass
    def chooseDT(self):
        self.DT = self.runCFL/max(max(self.cfl.flat),1.0e-8)
    def updateTimeHistory(self):
        self.m_last.flat[:] = self.m_tmp.flat
        if self.tau_last != None:
            self.tau_last.flat[:] = self.tau_tmp.flat
        if self.numDiff_last != None:
            self.numDiff_last.flat[:] = self.numDiff_tmp.flat
     #mwf add
    def updateStage(self):
        """
        increment stage counters and push necessary information
        into stage arrays
        """
        pass
    #mwf add
    def setInitialStageValues(self):
        """
        set the all stage values assuming this is the first step
        after a problem reset
        """
        pass
    #mwf add
    def updateNumericalFlux(self,f,d):
        """
        Given advective and diffusive numerical fluxes and derivative update
        """
        pass
    def updateNumericalFluxJacobian(self,df):
        """
        Given advective and diffusive numerical fluxes and derivative update
        """
        pass
    
class InnerTheta:
    """
    Broken,need to fix u/mass mismatch
    """
    def __init__(self,scalarTransport,theta=0.5):
        self.runCFL=0.1
        self.theta=theta
        self.duStar_du = self.theta
        self.u  = Vec(scalarTransport.nFreeDOF_global)
        self.oneMtheta_u_last = Vec(scalarTransport.nFreeDOF_global)
        self.u_tmp = Vec(scalarTransport.nFreeDOF_global)
        self.scalarTransport = scalarTransport
        self.cfl = scalarTransport.q['cfl']
        self.DT = 1.0
        self.massIsImplicit = True
        self.advectionIsImplicit = True
        self.diffusionIsImplicit = True
        self.reactionIsImplicit = True
        self.stabilizationIsImplicit = True
        self.shockCapturingIsImplicit = True
        self.m_last = Numeric.array(scalarTransport.q['m'])
        if scalarTransport.q.has_key('tau'):
            self.tau_last = Numeric.array(scalarTransport.q['tau'])
            self.tau_tmp = Numeric.array(scalarTransport.q['tau'])
        else:
            self.tau_last = None
            self.tau_tmp = None
        if scalarTransport.q.has_key('numDiff'):
            self.numDiff_last = Numeric.array(scalarTransport.q['numDiff'])
            self.numDiff_tmp = Numeric.array(scalarTransport.q['numDiff'])
        else:
            self.numDiff_last = None
            self.numDiff_tmp = None
    def calculateU(self,u):
        self.u_tmp[:]=u
        u*=self.theta
        u+=self.oneMtheta_u_last
    def updateMass(self,m,mt,dm,dmt):
        self.m_tmp = m
        mt.flat[:]=m.flat
        mt-=self.m_last
        mt/=self.DT
        dmt.flat[:]=dm.flat
        dmt/=self.DT
    def updateAdvection(self,f,df):
        pass
    def updateDivF(self,div_f,ddiv_f):
        pass
    def updateDivA(self,div_a,ddiv_a):
        pass
    def updateDiffusion(self,a,da):
        pass
    def updateGradients(self,dphi,grad_phi_x_grad_w,grad_u_x_grad_w):
        pass
    def updateReaction(self,r,dr):
        pass
    def updateStabilization(self,tau):
        self.tau_tmp.flat[:] = tau.flat
        tau.flat[:] = self.tau_last.flat
    def updateShockCapturing(self,numDiff):
        self.numDiff_tmp.flat[:] = numDiff.flat
        numDiff.flat[:] = self.numDiff_last.flat
    def updatePDEResidual(self,pdeResidual,dpdeResidual):
        pass
    def chooseDT(self):
        self.DT = self.runCFL/max(max(self.cfl.flat),1.0e-8)
    def updateTimeHistory(self):
        self.m_last.flat[:] = self.m_tmp.flat
        if self.tau_last != None:
            self.tau_last.flat[:] = self.tau_tmp.flat
        if self.numDiff_last != None:
            self.numDiff_last.flat[:] = self.numDiff_tmp.flat
        self.oneMtheta_u_last[:] = self.u_tmp.flat
        self.oneMtheta_u_last *= (1.0 - self.theta)
    #mwf add
    def updateStage(self):
        """
        increment stage counters and push necessary information
        into stage arrays
        """
        pass
    #mwf add
    def setInitialStageValues(self):
        """
        set the all stage values assuming this is the first step
        after a problem reset
        """
        pass
    #mwf add
    def updateNumericalFlux(self,f,d):
        """
        Given advective and diffusive numerical fluxes and derivative update
        """
        pass
    def updateNumericalFluxJacobian(self,df):
        """
        Given advective and diffusive numerical fluxes and derivative update
        """
        pass
    
class ForwardEuler:
    def __init__(self,scalarTransport):
        self.runCFL=0.1
        self.u = None
        self.duStar_du = None
        self.scalarTransport = scalarTransport
        self.cfl = scalarTransport.q['cfl']
        self.DT = 1.0
        self.massIsImplicit = True
        self.advectionIsImplicit = False
        self.diffusionIsImplicit = False
        self.reactionIsImplicit = False
        self.stabilizationIsImplicit = True
        self.shockCapturingIsImplicit = False
        self.m_last = Numeric.array(scalarTransport.q['m'])
        self.f_last = Numeric.array(scalarTransport.q['f'])
        self.div_f_last = Numeric.array(scalarTransport.q['div(f)'])
        if scalarTransport.q.has_key('Lstar*w'):
            self.Lstar_w_last = Numeric.array(scalarTransport.q['Lstar*w'])
        else:
            self.Lstar_w_last = None
        self.a_last = Numeric.array(scalarTransport.q['a'])
        self.r_last = Numeric.array(scalarTransport.q['r'])
        if scalarTransport.q.has_key('tau'):
            self.tau_last = Numeric.array(scalarTransport.q['tau'])
            self.tau_tmp = Numeric.array(scalarTransport.q['tau'])
        else:
            self.tau_last = None
            self.tau_tmp = None
        if scalarTransport.q.has_key('numDiff'):
            self.numDiff_last = Numeric.array(scalarTransport.q['numDiff'])
            self.numDiff_tmp = Numeric.array(scalarTransport.q['numDiff'])
        else:
            self.numDiff_last = None
            self.numDiff_tmp = None
        self.m_tmp = Numeric.array(scalarTransport.q['m'])
        self.f_tmp = Numeric.array(scalarTransport.q['f'])
        self.div_f_tmp = Numeric.array(scalarTransport.q['div(f)'])
        self.div_a_tmp = Numeric.array(scalarTransport.q['div(a)'])
        if scalarTransport.q.has_key('Lstar*w'):
            self.Lstar_w_tmp = Numeric.array(scalarTransport.q['Lstar*w'])
        else:
            self.Lstar_w_tmp = None
	if scalarTransport.ebq_global.has_key('advectiveFlux*dx_f'):
	    self.advectiveFlux_tmp = Numeric.array(scalarTransport.ebq_global['advectiveFlux*dx_f'])
	    self.advectiveFlux_last = Numeric.array(scalarTransport.ebq_global['advectiveFlux*dx_f'])
	else:
	    self.advectiveFlux_tmp = None
	    self.advectiveFlux_last = None
	if scalarTransport.ebq_global.has_key('diffusiveFlux*dx_a'):
	    self.diffusiveFlux_tmp = Numeric.array(scalarTransport.ebq_global['diffusiveFlux*dx_a'])
	    self.diffusiveFlux_last = Numeric.array(scalarTransport.ebq_global['diffusiveFlux*dx_a'])
	else:
	    self.diffusiveFlux_tmp = None
	    self.diffusiveFlux_last = None
        self.grad_u_x_grad_w_last = Numeric.array(scalarTransport.q['grad(u)*grad(w)'])
        self.grad_phi_x_grad_w_last = Numeric.array(scalarTransport.q['grad(phi)*grad(w)'])
        self.grad_u_x_grad_w_tmp = Numeric.array(scalarTransport.q['grad(u)*grad(w)'])
        self.grad_phi_x_grad_w_tmp = Numeric.array(scalarTransport.q['grad(phi)*grad(w)'])
        self.a_tmp = Numeric.array(scalarTransport.q['a'])
        self.r_tmp = Numeric.array(scalarTransport.q['r'])
    def calculateU(self,u):
        self.u = u
    def updateMass(self,m,mt,dm,dmt):
        self.m_tmp[:] = m
        mt.flat[:]=m.flat
        mt-=self.m_last
        mt/=self.DT
        dmt.flat[:]=dm.flat
        dmt/=self.DT
    def updateAdvection(self,f,df):
        self.f_tmp.flat[:] = f.flat
        f.flat[:] = self.f_last.flat
        df.flat[:]=0.0
    def updateDivF(self,div_f,ddiv_f):
        self.div_f_tmp.flat[:] = div_f.flat
        div_f.flat[:] = self.div_f_last.flat
        ddiv_f[:]=0.0
    def updateDivA(self,div_a,ddiv_a):
        self.div_a_tmp.flat[:] = div_a.flat
        div_a.flat[:] = self.div_a_last.flat
        ddiv_a[:]=0.0
    def updateAdjoint(self,Lstar_w):
        self.Lstar_w_tmp.flat[:] = Lstar_w.flat
        Lstar_w.flat[:] = self.Lstar_w_last.flat
    def updateDiffusion(self,a,da):
        self.a_tmp.flat[:] = a.flat
        a.flat[:] = self.a_last.flat
        da.flat[:] = 0.0
    def updateGradients(self,dphi,grad_phi_x_grad_w,grad_u_x_grad_w):
        dphi.flat[:]=0.0
        self.grad_phi_x_grad_w_tmp.flat[:] = grad_phi_x_grad_w.flat
        grad_phi_x_grad_w.flat[:] = self.grad_phi_x_grad_w_last.flat
        self.grad_u_x_grad_w_tmp.flat[:] = grad_u_x_grad_w.flat
        grad_u_x_grad_w.flat[:] = self.grad_u_x_grad_w_last.flat
    def updateReaction(self,r,dr):
        self.r_tmp.flat[:] = r.flat
        r.flat[:] = self.r_last.flat
        dr.flat[:] = 0.0
    def updateStabilization(self,tau):
        self.tau_tmp.flat[:] = tau.flat
        tau.flat[:] = self.tau_last.flat
    def updateShockCapturing(self,numDiff):
        self.numDiff_tmp.flat[:] = numDiff.flat
        numDiff.flat[:] = self.numDiff_last.flat
    def chooseDT(self):
        self.DT = self.runCFL/max(max(self.cfl.flat),1.0e-8)
    def updateTimeHistory(self):
        self.m_last.flat[:] = self.m_tmp.flat
        self.f_last.flat[:] = self.f_tmp.flat
        self.div_f_last.flat[:] = self.div_f_tmp.flat
        self.div_a_last.flat[:] = self.div_a_tmp.flat
        if self.Lstar_w_last != None:
            self.Lstar_w_last.flat[:] = self.Lstar_w_tmp.flat
        self.a_last.flat[:] = self.a_tmp.flat
        self.grad_phi_x_grad_w_last.flat[:] = self.grad_phi_x_grad_w_tmp.flat
        self.grad_u_x_grad_w_last.flat[:] = self.grad_u_x_grad_w_tmp.flat
        self.r_last.flat[:] = self.r_tmp.flat
        if self.tau_last != None:
            self.tau_last.flat[:] = self.tau_tmp.flat
        if self.numDiff_last != None:
            self.numDiff_last.flat[:] = self.numDiff_tmp.flat
	if self.advectiveFlux_last != None:
	    self.advectiveFlux_last.flat[:] = self.advectiveFlux_tmp.flat
	if self.diffusiveFlux_last != None:
	    self.diffusiveFlux_last.flat[:] = self.diffusiveFlux_tmp.flat
    def updateStage(self):
        """
        increment stage counters and push necessary information
        into stage arrays
        """
        pass
    def setInitialStageValues(self):
        """
        set the all stage values assuming this is the first step
        after a problem reset
        """
        pass
    def updateNumericalFlux(self,advectiveFlux,diffusiveFlux):
        """
        Given advective and diffusive numerical fluxes update
        """
	self.advectiveFlux_tmp.flat[:]=advectiveFlux.flat
	advectiveFlux.flat[:]=self.advectiveFlux_last.flat
	self.diffusiveFlux_tmp.flat[:]=diffusiveFlux.flat
	diffusiveFlux.flat[:]=self.diffusiveFlux_last.flat
    def updateNumericalFluxJacobian(self,numericalFluxJacobian):
        """
        Given advective and diffusive numerical fluxes and derivative update
        """
	numericalFluxJacobian.flat[:]=0.0
    
class OuterTheta:
    def __init__(self,scalarTransport,
                 thetaAdvection=0.5,
                 thetaDiffusion=0.5,
                 thetaReaction=0.5):
        self.runCFL = 0.1
        self.u = None
        self.duStar_du = None
        self.thetaAdvection=thetaAdvection
        self.thetaDiffusion=thetaDiffusion
        self.thetaReaction=thetaReaction
        self.scalarTransport = scalarTransport
        self.cfl = scalarTransport.q['cfl']
        self.DT = 1.0
        self.massIsImplicit = True
        self.advectionIsImplicit = True
        self.diffusionIsImplicit = True
        self.reactionIsImplicit = True
        self.stabilizationIsImplicit = True
        self.shockCapturingIsImplicit = True
        self.m_last = Numeric.array(scalarTransport.q['m'])
        self.f_last = Numeric.array(scalarTransport.q['f'])
        self.div_f_last = Numeric.array(scalarTransport.q['div(f)'])
        self.div_a_last = Numeric.array(scalarTransport.q['div(a)'])
        if scalarTransport.q.has_key('Lstar*w'):
            self.Lstar_w_last = Numeric.array(scalarTransport.q['Lstar*w'])
        else:
            self.Lstar_w_last = None
        self.a_last = Numeric.array(scalarTransport.q['a'])
        self.r_last = Numeric.array(scalarTransport.q['r'])
        if scalarTransport.q.has_key('tau'):
            self.tau_last = Numeric.array(scalarTransport.q['tau'])
            self.tau_tmp = Numeric.array(scalarTransport.q['tau'])
        else:
            self.tau_last = None
            self.tau_tmp = None
        if scalarTransport.q.has_key('numDiff'):
            self.numDiff_last = Numeric.array(scalarTransport.q['numDiff'])
            self.numDiff_tmp = Numeric.array(scalarTransport.q['numDiff'])
        else:
            self.numDiff_last = None
            self.numDiff_tmp = None
        self.m_tmp = Numeric.array(scalarTransport.q['m'])
        self.f_tmp = Numeric.array(scalarTransport.q['f'])
        self.div_f_tmp = Numeric.array(scalarTransport.q['div(f)'])
        self.div_a_tmp = Numeric.array(scalarTransport.q['div(a)'])
        if scalarTransport.q.has_key('Lstar*w'):
            self.Lstar_w_tmp = Numeric.array(scalarTransport.q['Lstar*w'])
        else:
            self.Lstar_w_tmp = None
        if scalarTransport.ebq_global.has_key('advectiveFlux*dx_f'):
            self.advectiveFlux_tmp = Numeric.array(scalarTransport.ebq_global['advectiveFlux*dx_f'])
            self.advectiveFlux_last = Numeric.array(scalarTransport.ebq_global['advectiveFlux*dx_f'])
        else:
            self.advectiveFlux_tmp = None
            self.advectiveFlux_last = None
        if scalarTransport.ebq_global.has_key('diffusiveFlux*dx_a'):
            self.diffusiveFlux_tmp = Numeric.array(scalarTransport.ebq_global['diffusiveFlux*dx_a'])
            self.diffusiveFlux_last = Numeric.array(scalarTransport.ebq_global['diffusiveFlux*dx_a'])
        else:
            self.diffusiveFlux_tmp = None
            self.diffusiveFlux_last = None
        self.grad_u_x_grad_w_last = Numeric.array(scalarTransport.q['grad(u)*grad(w)'])
        self.grad_phi_x_grad_w_last = Numeric.array(scalarTransport.q['grad(phi)*grad(w)'])
        self.grad_u_x_grad_w_tmp = Numeric.array(scalarTransport.q['grad(u)*grad(w)'])
        self.grad_phi_x_grad_w_tmp = Numeric.array(scalarTransport.q['grad(phi)*grad(w)'])
        self.a_tmp = Numeric.array(scalarTransport.q['a'])
        self.r_tmp = Numeric.array(scalarTransport.q['r'])
    def calculateU(self,u):
        self.u = u
    def updateMass(self,m,mt,dm,dmt):
        self.m_tmp[:] = m
        mt.flat[:]=m.flat
        mt-=self.m_last
        mt/=self.DT
        dmt.flat[:]=dm.flat
        dmt/=self.DT
    def updateAdvection(self,f,df):
        self.f_tmp.flat[:] = f.flat
        f *= self.thetaAdvection
        f.flat[:] += (1.0-self.thetaAdvection)*self.f_last.flat
        df *= self.thetaAdvection
    def updateDivF(self,div_f,ddiv_f):
        self.div_f_tmp.flat[:] = div_f.flat
        div_f *=self.thetaAdvection
        div_f.flat[:] += (1.0-self.thetaAdvection)*self.div_f_last.flat
        ddiv_f[:]*=self.thetaAdvection
    def updateDivA(self,div_a,ddiv_a):
        self.div_a_tmp.flat[:] = div_a.flat
        div_a *=self.thetaAdvection
        div_a.flat[:] += (1.0-self.thetaAdvection)*self.div_a_last.flat
        ddiv_a[:]*=self.thetaDiffusion
    def updateAdjoint(self,Lstar_w):
        self.Lstar_w_tmp.flat[:] = Lstar_w.flat
        Lstar_w *=self.thetaAdvection
        Lstar_w.flat[:] += (1.0-self.thetaAdvection)*self.Lstar_w_last.flat
    def updateDiffusion(self,a,da):
        self.a_tmp.flat[:] = a.flat
        a.flat[:] = self.a_last.flat
    def updateGradients(self,dphi,grad_phi_x_grad_w,grad_u_x_grad_w):
        dphi.flat[:] *= self.thetaDiffusion
        self.grad_phi_x_grad_w_tmp.flat[:] = grad_phi_x_grad_w.flat
        grad_phi_x_grad_w.flat[:] *= self.thetaDiffusion
        grad_phi_x_grad_w.flat[:] += (1.0-self.thetaDiffusion)*self.grad_phi_x_grad_w_last.flat
        self.grad_u_x_grad_w_tmp.flat[:] = grad_u_x_grad_w.flat
        grad_u_x_grad_w.flat[:] *= self.thetaDiffusion
        grad_u_x_grad_w.flat[:] += (1.0-self.thetaDiffusion)*self.grad_u_x_grad_w_last.flat
    def updateReaction(self,r,dr):
        self.r_tmp.flat[:] = r.flat
        r *= self.thetaReaction
        r.flat[:] += (1.0 - self.thetaReaction)*self.r_last.flat
        dr.flat[:] *=self.thetaReaction
    def updateStabilization(self,tau):
        self.tau_tmp.flat[:] = tau.flat
        tau.flat[:] = self.tau_last.flat
    def updateShockCapturing(self,numDiff):
        self.numDiff_tmp.flat[:] = numDiff.flat
        numDiff.flat[:] = self.numDiff_last.flat
    def chooseDT(self):
        self.DT = self.runCFL/max(max(self.cfl.flat),1.0e-8)
    def updateTimeHistory(self):
        self.m_last.flat[:] = self.m_tmp.flat
        self.f_last.flat[:] = self.f_tmp.flat
        self.div_f_last.flat[:] = self.div_f_tmp.flat
        self.div_a_last.flat[:] = self.div_a_tmp.flat
        if self.Lstar_w_last != None:
            self.Lstar_w_last.flat[:] = self.Lstar_w_tmp.flat
        self.a_last.flat[:] = self.a_tmp.flat
        self.grad_phi_x_grad_w_last.flat[:] = self.grad_phi_x_grad_w_tmp.flat
        self.grad_u_x_grad_w_last.flat[:] = self.grad_u_x_grad_w_tmp.flat
        self.r_last.flat[:] = self.r_tmp.flat
        if self.tau_last != None:
            self.tau_last.flat[:] = self.tau_tmp.flat
        if self.numDiff_last != None:
            self.numDiff_last.flat[:] = self.numDiff_tmp.flat
	if self.advectiveFlux_last != None:
	    self.advectiveFlux_last.flat[:] = self.advectiveFlux_tmp.flat
	if self.diffusiveFlux_last != None:
	    self.diffusiveFlux_last.flat[:] = self.diffusiveFlux_tmp.flat
    #mwf add
    def updateStage(self):
        """
        increment stage counters and push necessary information
        into stage arrays
        """
        pass
    #mwf add
    def setInitialStageValues(self):
        """
        set the all stage values assuming this is the first step
        after a problem reset
        """
        pass
    #mwf add
    def updateNumericalFlux(self,advectiveFlux,diffusiveFlux):
        """
        Given advective and diffusive numerical fluxes and derivative update
        """
 	#TODO fix numerical flux jacobian or place where these are lagged/avg'd
	self.advectiveFlux_tmp.flat[:] = advectiveFlux.flat
        advectiveFlux.flat[:] *= self.thetaAdvection
	advectiveFlux.flat[:] += (1.0 - self.thetaAdvection)*self.advectiveFlux_last.flat
	self.diffusiveFlux_tmp.flat[:] = diffusiveFlux.flat
        diffusiveFlux.flat[:] *= self.thetaDiffusion
	diffusiveFlux.flat[:] += (1.0 - self.thetaDiffusion)*self.diffusiveFlux_last.flat
    def updateNumericalFluxJacobian(self,df):
        """
        Given advective and diffusive numerical fluxes and derivative update
        """
        df.flat[:] *= self.thetaAdvection
class IMEX(OuterTheta):
    def __init__(self,scalarTransport):
        OuterTheta.__init__(self,scalarTransport,thetaAdvection = 0.0, thetaDiffusion = 1.0)
        self.advectionIsImplicit = False

"""
define set of explicit strong stability preserving (SSP) Runge Kutta
schemes as described in Gottlieb and Shu notes. These are SSP for
linear spatial discretizations that also have a few caveats similar to
coercivity I believe.

TODO:
  get working for new version of code

  mwf
"""

class SSPRKintegration(TimeIntegration):
    """
    implementation of Explicit RK schemes for ODE that are SSP for
    linear operators (aka) spatial discretizations in 

    \od{\vec y}{t} = \mat{L}\vec y

    See Gottlieb, Shu, Tadmor siam review article and my notes

    This is implemented following Gottlieb and Shu etal formulation rather
    than more traditional one used in p2mesh code.

    It also is inteded for use in a Rothe framework rather than MOL.

    """
    def __init__(self,scalarTransport,order=2,runCFL=0.1):
        """
        just setup stage coefficient arrays based on order

        could use recursive formula if want to fill in entire
        alpha matrix (lower tridiagonal)
        """
        TimeIntegration.__init__(self,scalarTransport)
        #change base class info for explicit RK methods
        self.massIsImplicit = True #assume linearly implicit disc.
        self.advectionIsImplicit = False
        self.diffusionIsImplicit = False
        self.reactionIsImplicit  = False
        self.stabilizationIsImplicit  = False
        self.shockCapturingIsImplicit = False
        #
        self.order  = order #order of approximation
        self.nstages= order #number of stages total
        self.lstage = 0     #last stage completed
        self.runCFL = runCFL
        self.DT     = 1.0
        self.alpha  = self.getAlphaCoefs(order)
        self.u      = None  #"current" value of solution
        self.cfl    = scalarTransport.q['cfl']
        #types of terms have to account for in system solving
        eqTermsAll = ['m','f','div(f)','div(a)','Lstar*w','a','r','tau',
                      'numDiff','grad(u)*grad(w)','grad(phi)*grad(w)',
                      'advectiveFlux*dx_f','diffusiveFlux*dx_a']
        self.eqTerms = []
        for term in eqTermsAll:
            if scalarTransport.q.has_key(term):
                self.eqTerms.append(term)
            elif scalarTransport.ebq_global.has_key(term):
                self.eqTerms.append(term)
            #end term in scalarTransport
        #end for all possible terms allowed
        #mwf debug
        print """SSPRK ctor order= %d nstages= %d\n
        """ % (self.order,self.nstages)
        #mwf debug
        #print """SSPRK ctor order= %d nstages= %d\n eqTerms= %s
        #      """ % (self.order,self.nstages,self.eqTerms)
        #variables needed for storing solution at time levels
        self.stageValues = [] #stage values for u, u[0]= u^n not U^1
        for i in range(self.nstages+1):
            self.stageValues.append({})
            #mwf debug
            #print """SSPRK ctor stageValues for stages %d """ % i
            #terms in quadrature array
            for term in self.eqTerms:
                if scalarTransport.q.has_key(term):
                    self.stageValues[i][term] = Numeric.array(
                        scalarTransport.q[term])
                elif scalarTransport.ebq_global.has_key(term):
                    self.stageValues[i][term] = Numeric.array(
                        scalarTransport.ebq_global[term])
                #end if on 
            #end for term
        #end for i

    #end init

    def getAlphaCoefs(self,order):
        """
        alpha matrix in Gottlieb, Shu, Tadmor review paper
        """
        alpha = Numeric.zeros(self.order,Numeric.Float)
        if order == 1:
            alpha[0] = 1.0
        elif order == 2:
            alpha[0] = 0.5;     alpha[1] = 0.5
        elif order == 3:
            alpha[0] = 1./3.;   alpha[1] = 1./2.; alpha[2] = 1./6.
        elif order == 4:
            alpha[0] = 3./8.;   alpha[1] = 1./3.; alpha[2] = 1./4.;
            alpha[3] = 1./24.;
        elif order == 5:
            alpha[0] = 11./30.; alpha[1] = 3./8.; alpha[2] = 1./6.;
            alpha[3] = 1./12.;  alpha[4] = 1./120.
        else:
            print 'order= ',order,' not implemented. Quitting !!!'
            sys.exit(1)
        #end if on order
        return alpha

    def updateStage(self):
        """
        increment stage counter by 1.
        lstage here is last stage completed
        """
        self.lstage +=1
        if self.lstage < 0 or self.lstage > self.nstages:
            print """stage= %d out of allowed bounds [0,%d]""" % (self.lstage,self.nstages)
        #end if
        #mwf debug
        #print """SSPRK updateStage finished lstage= %d m= %s""" % (self.lstage,
        #                                                          self.stageValues[self.lstage]['m'])
        #if self.lstage < self.nstages:
        #    print """SSPRK updateStage lstage+1= %d m= %s """ % (self.lstage+1,
        #                                                         self.stageValues[self.lstage+1]['m'])
    def updateMass(self,m,mt,dm,dmt):
        """
        Calculate m_t, and dm_t/du given m,dm.
        """
        #except for last stage all stage values can be had
        #by solving problem with dm/dt = (m^k-m^{k-1})/dt
        #mwf debug
        #print """SSPRK lstage= %d DT = %g updateMass in m= \n%s""" % (self.lstage,self.DT,m)  
        #print """SSPRK updateMass in stage M[%d]=%s""" % (self.lstage,
        #                                                  self.stageValues[self.lstage]['m'])
        #print """SSPRK updateMass in stage M[%d]=%s""" % (self.lstage+1,
        #                                                  self.stageValues[self.lstage+1]['m'])

        
        if self.DT <= 1.0e-24:
            print 'WARNING DT= ',self.DT,' too small in updateMass quitting '
            sys.exit(1)
        #end if
        self.stageValues[self.lstage+1]['m'].flat[:] = m.flat[:]
        mt.flat[:]=m.flat[:]
        mtmp = Numeric.zeros(m.shape,Numeric.Float)
        if self.lstage == self.nstages-1:
            for i in range(self.nstages): #go from 0 to s-1
                mtmp.flat[:]+=self.alpha[i]*self.stageValues[i]['m'].flat[:]
            #end i loop through stages
        else:
            mtmp.flat[:]=self.stageValues[self.lstage]['m'].flat[:] 
        #end if last stage
        mt-=mtmp
        mt/=self.DT
        dmt.flat[:]=dm.flat[:]
        dmt/=self.DT
        #mwf debug
        #print """SSPRK updateMass out stage M[%d]=%s""" % (self.lstage,
        #                                                   self.stageValues[self.lstage]['m'])
        #print """SSPRK updateMass out stage M[%d]=%s""" % (self.lstage+1,
        #                                                   self.stageValues[self.lstage+1]['m'])
        #print """SSPRK updateMass out mt= %s """ % mt
        
    def updateAdvection(self,f,df):
        """
        Given f,df, update with f^*,df^*/du,...
        """
        term = 'f'
        if term in self.eqTerms:
            self.stageValues[self.lstage+1][term].flat[:]=f.flat[:]
            if self.lstage == self.nstages-1:
                f.flat[:] = self.alpha[self.nstages-1]*self.stageValues[self.lstage][term].flat[:]
            else:
                f.flat[:] = self.stageValues[self.lstage][term].flat[:]
            #end if last stage
            df.flat[:]= 0.0
        #end if in eqTerms
      
    def updateDivF(self,div_f,ddiv_f):
        """
        Given dif_f,ddif_f, update with dif_f^*,ddif_f^*/du,...
        """
        term = 'div(f)'
        if term in self.eqTerms:
            self.stageValues[self.lstage+1][term].flat[:]=div_f.flat[:]
            if self.lstage == self.nstages-1:
                div_f.flat[:] = self.alpha[self.nstages-1]*self.stageValues[self.lstage][term].flat[:]
            else:
                div_f.flat[:] = self.stageValues[self.lstage][term].flat[:]
            #end if last stage
            ddiv_f.flat[:]= 0.0
        #end term found
    def updateDivA(self,div_a,ddiv_a):
        """
        Given div_a,ddiv_a, update with div_a^*,ddiv_a^*/du,...
        """
        term = 'div(a)'
        if term in self.eqTerms:
            self.stageValues[self.lstage+1][term].flat[:]=div_a.flat[:]
            if self.lstage == self.nstages-1:
                div_a.flat[:] = self.alpha[self.nstages-1]*self.stageValues[self.lstage][term].flat[:]
            else:
                div_a.flat[:] = self.stageValues[self.lstage][term].flat[:]
            #end if last stage
            ddiv_a.flat[:]= 0.0
        #end term found
    def updateAdjoint(self,Lstar_w):
        term = 'Lstar*w'
        if term in self.eqTerms:
            self.stageValues[self.lstage+1][term].flat[:]=Lstar_w.flat[:]
            if self.lstage == self.nstages-1:
                Lstar_w.flat[:] = self.alpha[self.nstages-1]*self.stageValues[self.lstage][term].flat[:]
            else:
                Lstar_w.flat[:] = self.stageValues[self.lstage][term].flat[:]
            #end if last stage
        #end if term found
    def updateDiffusion(self,a,da):
        """
        Given a,da, update with a^*,da^*/du,...
        """
        term = 'a'
        if term in self.eqTerms:
            self.stageValues[self.lstage+1][term].flat[:]=a.flat[:]
            if self.lstage == self.nstages-1:
                a.flat[:] = self.alpha[self.nstages-1]*self.stageValues[self.lstage][term].flat[:]
            else:
                a.flat[:] = self.stageValues[self.lstage][term].flat[:]
            #end if last stage
            da.flat[:]= 0.0
        #end if term found
    def updateGradients(self,dphi,grad_phi_x_grad_w,grad_u_x_grad_w):
        term = 'grad(phi)*grad(w)'
        if term in self.eqTerms:
            self.stageValues[self.lstage+1][term].flat[:]=grad_phi_x_grad_w.flat[:]
            #don't do anything different to phi since a gets multiplied by dt coefs
            grad_phi_x_grad_w.flat[:] = self.stageValues[self.lstage][term].flat[:]
            dphi.flat[:] = 0.0
        term = 'grad(u)*grad(w)'
        if term in self.eqTerms:
            self.stageValues[self.lstage+1][term].flat[:]=grad_u_x_grad_w.flat[:]
            #don't do anything different to phi since a gets multiplied by dt coefs
            grad_u_x_grad_w.flat[:] = self.stageValues[self.lstage][term].flat[:]
        #end if term found
    def updateReaction(self,r,dr):
        """
        Given r,dr, update with r^*,dr^*/du,...
        """
        term = 'r'
        if term in self.eqTerms:
            #mwf debug
            #print """SSPRK lstage= %d updateRxn in r= %s""" % (self.lstage,r)  
            #print """SSPRK updateRxn in stage R[%d]=%s""" % (self.lstage,
            #                                                 self.stageValues[self.lstage][term])
            #print """SSPRK updateMass in stage R[%d]=%s""" % (self.lstage+1,
            #                                                  self.stageValues[self.lstage+1][term])

            self.stageValues[self.lstage+1][term].flat[:]=r.flat[:]
            if self.lstage == self.nstages-1:
                r.flat[:] = self.alpha[self.nstages-1]*self.stageValues[self.lstage][term].flat[:]
            else:
                r.flat[:] = self.stageValues[self.lstage][term].flat[:]
            #end if on last stage
            #mwf debug
            #print """SSPRK updateRxn out stage R[%d]=%s""" % (self.lstage,
            #                                                 self.stageValues[self.lstage][term])
            #print """SSPRK updateMass out stage R[%d]=%s""" % (self.lstage+1,
            #                                               self.stageValues[self.lstage+1][term])

            dr.flat[:]= 0.0
        #end if ferm found

    def updateStabilization(self,tau):
        """
        Given \tau, update with \tau_*
        """
        term = 'tau'
        if term in self.eqTerms:
            self.stageValues[self.lstage+1][term].flat[:]=tau.flat[:]
            if self.lstage == self.nstages-1:
                tau.flat[:] = self.alpha[self.nstages-1]*self.stageValues[self.lstage][term].flat[:]
            else:
                tau.flat[:] = self.stageValues[self.lstage][term].flat[:]
            #end if last stage
        #end if

    def updateShockCapturing(self,numDiff):
        """
        Given numDiff, update with numDiff_*
        """
        term = 'numDiff'
        if term in self.eqTerms:
            self.stageValues[self.lstage+1][term].flat[:]=numDiff.flat[:]
            if self.lstage == self.nstages-1:
                numDiff.flat[:] = self.alpha[self.nstages-1]*self.stageValues[self.lstage][term].flat[:]
            else:
                numDiff.flat[:] = self.stageValues[self.lstage][term].flat[:]
            #end if last stage
        #end if term found

    def chooseDT(self):
        """
        Modify self.DT
        """
        #mwf fix: cfl is always zero right now for scalar transport problem
        
        self.DT = self.runCFL/max(max(self.cfl.flat),1.0e-6)

    def setInitialStageValues(self):
        """
        setup all the stage values if we assume that the initial condition
        has been loaded into lstage+1
        """
        for i in range(self.nstages+1):
            if i != self.lstage+1:
                for term in self.eqTerms:
                    self.stageValues[i][term].flat[:] = self.stageValues[self.lstage+1][term].flat[:]
                #end for terms
            #end i not equal to lstage+1
        #end for i
        self.lstage = 0

    #end updateTimeHistory

    def updateTimeHistory(self):
        """
        setup for next time step, cycle U^s --> U^0
        """
        for term in self.eqTerms:
            self.stageValues[0][term].flat[:] = self.stageValues[self.nstages][term].flat[:]
        #end for
        self.lstage = 0

    #end updateTimeHistory
    #cek modified, mwf remove df too
    def updateNumericalFlux(self,f,d):
        """
        Given advective and diffusive numerical fluxes and derivative update with correct stage value
        """
        term = 'advectiveFlux*dx_f'
        if term in self.eqTerms:
            self.stageValues[self.lstage+1][term].flat[:]=f.flat[:]
            if self.lstage == self.nstages-1:
                f.flat[:] = self.alpha[self.nstages-1]*self.stageValues[self.lstage][term].flat[:]
            else:
                f.flat[:] = self.stageValues[self.lstage][term].flat[:]
            #end if last stage
        #end if in eqTerms
        term = 'diffusiveFlux*dx_a'
        if term in self.eqTerms:
            self.stageValues[self.lstage+1][term].flat[:]=d.flat[:]
            if self.lstage == self.nstages-1:
                d.flat[:] = self.alpha[self.nstages-1]*self.stageValues[self.lstage][term].flat[:]
            else:
                d.flat[:] = self.stageValues[self.lstage][term].flat[:]
            #end if last stage
        #end if in eqTerms
    def updateNumericalFluxJacobian(self,df):
        """
        Given flux derivative update with correct stage value
        """
        df.flat[:]= 0.0
      
## @}    

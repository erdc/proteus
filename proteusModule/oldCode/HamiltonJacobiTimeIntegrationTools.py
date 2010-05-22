from LinearAlgebra import *
import sys
"""
Tools for discretizing time-dependent Hamilton-Jacobi equations.
"""

## \defgroup TimeIntegration
#
# Tools for discretizing time-dependent Hamilton-Jacobi equations.                                       #
# @{   

class TimeIntegration:
    def __init__(self,hamiltonJacobi):
        self.u = None
        self.duStar_du = None
        self.DT = 1.0
        self.massIsImplicit = True
        self.hamiltonianIsImplicit = True
        self.stabilizationIsImplicit = True
        self.shockCapturingIsImplicit = True
    def calculateU(self,u):
        self.u = u
    def updateMass(self,m,mt,dm,dmt):
        pass
    def updateHamiltonian(self,h,dh,dhdu,grad_u):
        pass
    def updateAdjoint(self,Lstar_w):
        pass
    def updateStabilization(self,tau):
        pass
    def updateShockCapturing(self,numDiff):
        pass
    def updateGradients(self,grad_u_x_grad_w):
        """
        Given grad_u_x_grad_w,grad_phi_x_grad_w update  with grad_u_x_grad_w*, grad_phi_x_grad_w*
        """
        pass
    def chooseDT(self):
        pass
    def updateTimeHistory(self):
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
    def __init__(self,hamiltonJacobi):
        self.runCFL=0.1
        self.u = None
        self.duStar_du = None
        self.hamiltonJacobi = hamiltonJacobi
        self.cfl = hamiltonJacobi.q['cfl']
        self.DT = 1.0
        self.massIsImplicit = True
        self.hamiltonianIsImplicit = True
        self.stabilizationIsImplicit = True
        self.shockCapturingIsImplicit = True
        self.m_last = Numeric.array(hamiltonJacobi.q['m'])
        self.tau_last = Numeric.array(hamiltonJacobi.q['tau'])
        self.numDiff_last = Numeric.array(hamiltonJacobi.q['numDiff'])
        self.tau_tmp = Numeric.array(hamiltonJacobi.q['tau'])
        self.numDiff_tmp = Numeric.array(hamiltonJacobi.q['numDiff'])
    def calculateU(self,u):
        self.u = u
    def updateMass(self,m,mt,dm,dmt):
        self.m_tmp = m
        mt.flat[:]=m.flat
        mt-=self.m_last
        mt/=self.DT
        dmt.flat[:]=dm.flat
        dmt/=self.DT
    def updateHamiltonian(self,h,dh,dhdu,grad_u):
        dhdu.flat[:]=dh.flat
    def updateAdjoint(self,Lstar_w):
        pass
    def updateStabilization(self,tau):
        self.tau_tmp.flat[:] = tau.flat
        tau.flat[:] = self.tau_last.flat
    def updateShockCapturing(self,numDiff):
        self.numDiff_tmp.flat[:] = numDiff.flat
        numDiff.flat[:] = self.numDiff_last.flat
    def updateGradients(self,grad_u_x_grad_w):
        pass
    def chooseDT(self):
        self.DT = self.runCFL/max(max(self.cfl.flat),1.0e-8)
    def updateTimeHistory(self):
        self.m_last.flat[:] = self.m_tmp.flat
        self.tau_last.flat[:]=self.tau_tmp.flat
        self.numDiff_last.flat[:]=self.numDiff_tmp.flat
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
    def __init__(self,hamiltonJacobi):
        self.runCFL=0.1
        self.u = None
        self.duStar_du = None
        self.hamiltonJacobi = hamiltonJacobi
        self.cfl = hamiltonJacobi.q['cfl']
        self.DT = 1.0
        self.massIsImplicit = True
        self.hamiltonianIsImplicit = False
        self.stabilizationIsImplicit = True
        self.shockCapturingIsImplicit = False
        self.m_last = Numeric.array(hamiltonJacobi.q['m'])
        self.h_last = Numeric.array(hamiltonJacobi.q['h'])
        self.dh_last = Numeric.array(hamiltonJacobi.q['dh'])
        self.grad_u_last = Numeric.array(hamiltonJacobi.q['grad(u)'])
        self.tau_last = Numeric.array(hamiltonJacobi.q['tau'])
        self.numDiff_last = Numeric.array(hamiltonJacobi.q['numDiff'])
        self.m_tmp = Numeric.array(hamiltonJacobi.q['m'])
        self.h_tmp = Numeric.array(hamiltonJacobi.q['h'])
        self.dh_tmp = Numeric.array(hamiltonJacobi.q['dh'])
        self.grad_u_tmp = Numeric.array(hamiltonJacobi.q['grad(u)'])
        self.tau_tmp = Numeric.array(hamiltonJacobi.q['tau'])
        self.numDiff_tmp = Numeric.array(hamiltonJacobi.q['numDiff'])
        if hamiltonJacobi.q.has_key('Lstar*w'):
            self.Lstar_w_last = Numeric.array(hamiltonJacobi.q['Lstar*w'])
            self.Lstar_w_tmp = Numeric.array(hamiltonJacobi.q['Lstar*w'])
        else:
            self.Lstar_w_last = None
            self.Lstar_w_tmp  = None
        self.grad_u_x_grad_w_last = Numeric.array(hamiltonJacobi.q['grad(u)*grad(w)'])
        self.grad_u_x_grad_w_tmp = Numeric.array(hamiltonJacobi.q['grad(u)*grad(w)'])
    def calculateU(self,u):
        self.u = u
    def updateMass(self,m,mt,dm,dmt):
        self.m_tmp = m
        mt.flat[:]=m.flat
        mt-=self.m_last
        mt/=self.DT
        dmt.flat[:]=dm.flat
        dmt/=self.DT
    def updateHamiltonian(self,h,dh,dhdu,grad_u):
        self.h_tmp.flat[:] = h.flat
        h.flat[:] = self.h_last.flat
        self.dh_tmp.flat[:] = dh.flat
        dh.flat[:]=self.dh_last.flat
        dhdu.flat[:] = 0.0
        self.grad_u_tmp.flat[:] = grad_u.flat
        grad_u.flat[:] = self.grad_u_last.flat
    def updateAdjoint(self,Lstar_w):
        self.Lstar_w_tmp.flat[:] = Lstar_w.flat
        Lstar_w.flat[:] = self.Lstar_w_last.flat
    def updateStabilization(self,tau):
        self.tau_tmp.flat[:] = tau.flat
        tau.flat[:] = self.tau_last.flat
    def updateShockCapturing(self,numDiff):
        self.numDiff_tmp.flat[:] = numDiff.flat
        numDiff.flat[:] = self.numDiff_last.flat
    def updateGradients(self,grad_u_x_grad_w):
        grad_u_x_grad_w.flat[:] = self.grad_u_x_grad_w_last.flat
    def chooseDT(self):
        self.DT = self.runCFL/max(max(self.cfl.flat),1.0e-8)
    def updateTimeHistory(self):
        self.m_last.flat[:] = self.m_tmp.flat
        self.h_last.flat[:] = self.h_tmp.flat
        self.dh_last.flat[:] = self.dh_tmp.flat
        self.Lstar_w_last.flat[:] = self.Lstar_w_tmp.flat
        self.grad_u_last.flat[:] = self.grad_u_tmp.flat
        self.tau_last.flat[:] = self.tau_tmp.flat
        self.numDiff_last.flat[:] = self.numDiff_tmp.flat
        self.grad_u_x_grad_w_last.flat[:] = self.grad_u_x_grad_w_tmp.flat
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

class OuterTheta:
    def __init__(self,hamiltonJacobi,
                 theta=0.5):
        self.runCFL = 0.1
        self.u = None
        self.duStar_du = None
        self.theta=theta
        self.hamiltonJacobi = hamiltonJacobi
        self.cfl = hamiltonJacobi.q['cfl']
        self.DT = 1.0
        self.massIsImplicit = True
        self.hamiltonianIsImplicit = True
        self.stabilizationIsImplicit = True
        self.shockCapturingIsImplicit = True
        self.m_last = Numeric.array(hamiltonJacobi.q['m'])
        self.h_last = Numeric.array(hamiltonJacobi.q['h'])
        self.dh_last = Numeric.array(hamiltonJacobi.q['dh'])
        self.grad_u_last = Numeric.array(hamiltonJacobi.q['grad(u)'])
        self.tau_last = Numeric.array(hamiltonJacobi.q['tau'])
        self.numDiff_last = Numeric.array(hamiltonJacobi.q['numDiff'])
        self.m_tmp = Numeric.array(hamiltonJacobi.q['m'])
        self.h_tmp = Numeric.array(hamiltonJacobi.q['h'])
        self.dh_tmp = Numeric.array(hamiltonJacobi.q['dh'])
        self.grad_u_tmp = Numeric.array(hamiltonJacobi.q['grad(u)'])
        self.tau_tmp = Numeric.array(hamiltonJacobi.q['tau'])
        self.numDiff_tmp = Numeric.array(hamiltonJacobi.q['numDiff'])
    def calculateU(self,u):
        self.u = u
    def updateMass(self,m,mt,dm,dmt):
        self.m_tmp[:] = m
        mt.flat[:]=m.flat
        mt-=self.m_last
        mt/=self.DT
        dmt.flat[:]=dm.flat
        dmt/=self.DT
    def updateHamiltonian(self,h,dh,dhdu,grad_u):
        self.h_tmp.flat[:] = h.flat
        h *= self.theta
        h.flat[:] += (1.0-self.theta)*self.h_last.flat
        self.dh_tmp.flat[:] = dh.flat
        dh *=self.theta
        dhdu.flat[:] = dh.flat
        dh.flat[:] += (1.0-self.theta)*self.dh_last.flat
        self.grad_u_tmp.flat[:] = grad_u.flat
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
        self.h_last.flat[:] = self.h_tmp.flat
        self.dh_last.flat[:] = self.dh_tmp.flat
        self.grad_u_last.flat[:] = self.grad_u_tmp.flat
        self.tau_last.flat[:] = self.tau_tmp.flat
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

## @}

# A type of -*- python -*- file
"""
An optimized Advection-Diffusion-Reaction module
"""

import numpy
cimport numpy
from math import fabs
import proteus
from proteus import  cfemIntegrals, Quadrature
from proteus.NonlinearSolvers import NonlinearEquation
from proteus.FemTools import (DOFBoundaryConditions,
                              FluxBoundaryConditions,
                              C0_AffineLinearOnSimplexWithNodalBasis)
from proteus.Comm import globalMax, globalSum
from proteus.Profiling import  memory
from proteus.Profiling import  logEvent as log
from proteus.Transport import OneLevelTransport
from proteus.TransportCoefficients import TC_base
from proteus.SubgridError import SGE_base
from proteus.ShockCapturing import  ShockCapturing_base

cdef extern from "ADR.h" namespace "proteus":
    cdef cppclass cppADR_base:
        void calculateResidual(double* mesh_trial_ref,
                               double* mesh_grad_trial_ref,
                               double* mesh_dof,
                               int* mesh_l2g,
                               double* dV_ref,
                               double* u_trial_ref,
                               double* u_grad_trial_ref,
                               double* u_test_ref,
                               double* u_grad_test_ref,
                               double* elementDiameter,
                               double* cfl,
                               double CT_sge,
                               double sc_uref,
                               double sc_alpha,
                               double useMetrics,
                               double* mesh_trial_trace_ref,
                               double* mesh_grad_trial_trace_ref,
                               double* dS_ref,
                               double* u_trial_trace_ref,
                               double* u_grad_trial_trace_ref,
                               double* u_test_trace_ref,
                               double* u_grad_test_trace_ref,
                               double* normal_ref,
                               double* boundaryJac_ref,
                               int nElements_global,
                               int* u_l2g,
                               double* u_dof,
                               int* sdInfo_u_u_rowptr,
                               int* sdInfo_u_u_colind,
                               double* q_a,
			       double* q_v,
                               double* q_r,
			       int lag_shockCapturing,
			       double shockCapturingDiffusion,
			       double* q_numDiff_u,
			       double* q_numDiff_u_last,
                               int offset_u, int stride_u,
                               double* globalResidual,
                               int nExteriorElementBoundaries_global,
                               int* exteriorElementBoundariesArray,
                               int* elementBoundaryElementsArray,
                               int* elementBoundaryLocalElementBoundariesArray,
                               double* ebqe_a,
			       double* ebqe_v,
                               int* isDOFBoundary_u,
                               double* ebqe_bc_u_ext,
                               int* isFluxBoundary_u,
			       int* isAdvectiveFluxBoundary_u,
                               double* ebqe_bc_flux_u_ext,
			       double* ebqe_bc_advectiveFlux_u_ext,
                               double* ebqe_penalty,
                               double adjoint_sigma)
        void calculateJacobian(double* mesh_trial_ref,
                               double* mesh_grad_trial_ref,
                               double* mesh_dof,
                               int* mesh_l2g,
                               double* dV_ref,
                               double* u_trial_ref,
                               double* u_grad_trial_ref,
                               double* u_test_ref,
                               double* u_grad_test_ref,
                               double* elementDiameter,
                               double* cfl,
                               double CT_sge,
                               double sc_uref,
                               double sc_alpha,
                               double useMetrics,
                               double* mesh_trial_trace_ref,
                               double* mesh_grad_trial_trace_ref,
                               double* dS_ref,
                               double* u_trial_trace_ref,
                               double* u_grad_trial_trace_ref,
                               double* u_test_trace_ref,
                               double* u_grad_test_trace_ref,
                               double* normal_ref,
                               double* boundaryJac_ref,
                               int nElements_global,
                               int* u_l2g,
                               double* u_dof,
                               int* sdInfo_u_u_rowptr,
                               int* sdInfo_u_u_colind,
                               double* q_a,
			       double* q_v,
                               double* q_r,
			       int lag_shockCapturing,
			       double shockCapturingDiffusion,
			       double* q_numDiff_u,
			       double* q_numDiff_u_last,
                               int* csrRowIndeces_u_u,
                               int* csrColumnOffsets_u_u,
                               double* globalJacobian,
                               int nExteriorElementBoundaries_global,
                               int* exteriorElementBoundariesArray,
                               int* elementBoundaryElementsArray,
                               int* elementBoundaryLocalElementBoundariesArray,
                               double* ebqe_a,
			       double* ebqe_v,
                               int* isDOFBoundary_u,
                               double* ebqe_bc_u_ext,
                               int* isFluxBoundary_u,
			       int* isAdvectiveFluxBoundary_u,
                               double* ebqe_bc_flux_u_ext,
			       double* ebqe_bc_advectiveFlux_u_ext,
                               int* csrColumnOffsets_eb_u_u,
                               double* ebqe_penalty,
                               double adjoint_sigma)
    cppADR_base* newADR(int nSpaceIn,
                        int nQuadraturePoints_elementIn,
                        int nDOF_mesh_trial_elementIn,
                        int nDOF_trial_elementIn,
                        int nDOF_test_elementIn,
                        int nQuadraturePoints_elementBoundaryIn,
                        int CompKernelFlag)


cdef class ADR:
    """
    Optimized ADR member functions
    """

    cdef cppADR_base* thisptr
    def __cinit__(self,
                  int nSpaceIn,
                  int nQuadraturePoints_elementIn,
                  int nDOF_mesh_trial_elementIn,
                  int nDOF_trial_elementIn,
                  int nDOF_test_elementIn,
                  int nQuadraturePoints_elementBoundaryIn,
                  int CompKernelFlag):
        self.thisptr = newADR(nSpaceIn,
                              nQuadraturePoints_elementIn,
                              nDOF_mesh_trial_elementIn,
                              nDOF_trial_elementIn,
                              nDOF_test_elementIn,
                              nQuadraturePoints_elementBoundaryIn,
                              CompKernelFlag)
    def __dealloc__(self):
        del self.thisptr
    def calculateResidual(self,
                          numpy.ndarray mesh_trial_ref,
                          numpy.ndarray mesh_grad_trial_ref,
                          numpy.ndarray mesh_dof,
                          numpy.ndarray mesh_l2g,
                          numpy.ndarray dV_ref,
                          numpy.ndarray u_trial_ref,
                          numpy.ndarray u_grad_trial_ref,
                          numpy.ndarray u_test_ref,
                          numpy.ndarray u_grad_test_ref,
                          numpy.ndarray elementDiameter,
                          numpy.ndarray cfl,
                          double CT_sge,
                          double sc_uref,
                          double sc_alpha,
                          double useMetrics,
                          numpy.ndarray mesh_trial_trace_ref,
                          numpy.ndarray mesh_grad_trial_trace_ref,
                          numpy.ndarray dS_ref,
                          numpy.ndarray u_trial_trace_ref,
                          numpy.ndarray u_grad_trial_trace_ref,
                          numpy.ndarray u_test_trace_ref,
                          numpy.ndarray u_grad_test_trace_ref,
                          numpy.ndarray normal_ref,
                          numpy.ndarray boundaryJac_ref,
                          int nElements_global,
                          numpy.ndarray u_l2g,
                          numpy.ndarray u_dof,
                          numpy.ndarray sdInfo_u_u_rowptr,
                          numpy.ndarray sdInfo_u_u_colind,
                          numpy.ndarray q_a,
                          numpy.ndarray q_v,
                          numpy.ndarray q_r,
                          int lag_shockCapturing,
                          double shockCapturingDiffusion,
                          numpy.ndarray q_numDiff_u,
                          numpy.ndarray q_numDiff_u_last,
                          int offset_u,
                          int stride_u,
                          numpy.ndarray globalResidual,
                          int nExteriorElementBoundaries_global,
                          numpy.ndarray exteriorElementBoundariesArray,
                          numpy.ndarray elementBoundaryElementsArray,
                          numpy.ndarray elementBoundaryLocalElementBoundariesArray,
                          numpy.ndarray ebqe_a,
                          numpy.ndarray ebqe_v,
                          numpy.ndarray isDOFBoundary_u,
                          numpy.ndarray ebqe_bc_u_ext,
                          numpy.ndarray isFluxBoundary_u,
                          numpy.ndarray isAdvectiveFluxBoundary_u,
                          numpy.ndarray ebqe_bc_flux_u_ext,
                          numpy.ndarray ebqe_bc_advectiveFlux_u_ext,
                          numpy.ndarray ebqe_penalty,
                          double adjoint_sigma):
        self.thisptr.calculateResidual(<double*> mesh_trial_ref.data,
                                       <double*> mesh_grad_trial_ref.data,
                                       <double*> mesh_dof.data,
                                       <int*> mesh_l2g.data,
                                       <double*> dV_ref.data,
                                       <double*> u_trial_ref.data,
                                       <double*> u_grad_trial_ref.data,
                                       <double*> u_test_ref.data,
                                       <double*> u_grad_test_ref.data,
                                       <double*> elementDiameter.data,
                                       <double*> cfl.data,
                                       CT_sge,
                                       sc_uref,
                                       sc_alpha,
                                       useMetrics,
                                       <double*> mesh_trial_trace_ref.data,
                                       <double*> mesh_grad_trial_trace_ref.data,
                                       <double*> dS_ref.data,
                                       <double*> u_trial_trace_ref.data,
                                       <double*> u_grad_trial_trace_ref.data,
                                       <double*> u_test_trace_ref.data,
                                       <double*> u_grad_test_trace_ref.data,
                                       <double*> normal_ref.data,
                                       <double*> boundaryJac_ref.data,
                                       nElements_global,
                                       <int*> u_l2g.data,
                                       <double*> u_dof.data,
                                       <int*> sdInfo_u_u_rowptr.data,
                                       <int*> sdInfo_u_u_colind.data,
                                       <double*> q_a.data,
                                       <double*> q_v.data,
				       <double*> q_r.data,
				       lag_shockCapturing,
				       shockCapturingDiffusion,
				       <double*> q_numDiff_u.data,
				       <double*> q_numDiff_u_last.data,
                                       offset_u,
                                       stride_u,
                                       <double*> globalResidual.data,
                                       nExteriorElementBoundaries_global,
                                       <int*> exteriorElementBoundariesArray.data,
                                       <int*> elementBoundaryElementsArray.data,
                                       <int*> elementBoundaryLocalElementBoundariesArray.data,
                                       <double*> ebqe_a.data,
				       <double*> ebqe_v.data,
                                       <int*> isDOFBoundary_u.data,
                                       <double*> ebqe_bc_u_ext.data,
                                       <int*> isFluxBoundary_u.data,
				       <int*> isAdvectiveFluxBoundary_u.data,
                                       <double*> ebqe_bc_flux_u_ext.data,
				       <double*> ebqe_bc_advectiveFlux_u_ext.data,
                                       <double*> ebqe_penalty.data,
                                       adjoint_sigma)
    def calculateJacobian(self,
                          numpy.ndarray mesh_trial_ref,
                          numpy.ndarray mesh_grad_trial_ref,
                          numpy.ndarray mesh_dof,
                          numpy.ndarray mesh_l2g,
                          numpy.ndarray dV_ref,
                          numpy.ndarray u_trial_ref,
                          numpy.ndarray u_grad_trial_ref,
                          numpy.ndarray u_test_ref,
                          numpy.ndarray u_grad_test_ref,
                          numpy.ndarray elementDiameter,
                          numpy.ndarray cfl,
                          double CT_sge,
                          double sc_uref,
                          double sc_alpha,
                          double useMetrics,
                          numpy.ndarray mesh_trial_trace_ref,
                          numpy.ndarray mesh_grad_trial_trace_ref,
                          numpy.ndarray dS_ref,
                          numpy.ndarray u_trial_trace_ref,
                          numpy.ndarray u_grad_trial_trace_ref,
                          numpy.ndarray u_test_trace_ref,
                          numpy.ndarray u_grad_test_trace_ref,
                          numpy.ndarray normal_ref,
                          numpy.ndarray boundaryJac_ref,
                          int nElements_global,
                          numpy.ndarray u_l2g,
                          numpy.ndarray u_dof,
                          numpy.ndarray sdInfo_u_u_rowptr,
                          numpy.ndarray sdInfo_u_u_colind,
                          numpy.ndarray q_a,
                          numpy.ndarray q_v,
                          numpy.ndarray q_r,
                          int lag_shockCapturing,
                          double shockCapturingDiffusion,
                          numpy.ndarray q_numDiff_u,
                          numpy.ndarray q_numDiff_u_last,
                          numpy.ndarray csrRowIndeces_u_u,
                          numpy.ndarray csrColumnOffsets_u_u,
                          globalJacobian,
                          int nExteriorElementBoundaries_global,
                          numpy.ndarray exteriorElementBoundariesArray,
                          numpy.ndarray elementBoundaryElementsArray,
                          numpy.ndarray elementBoundaryLocalElementBoundariesArray,
                          numpy.ndarray ebqe_a,
                          numpy.ndarray ebqe_v,
                          numpy.ndarray isDOFBoundary_u,
                          numpy.ndarray ebqe_bc_u_ext,
                          numpy.ndarray isFluxBoundary_u,
                          numpy.ndarray isAdvectiveFluxBoundary_u,
                          numpy.ndarray ebqe_bc_flux_u_ext,
                          numpy.ndarray ebqe_bc_advectiveFlux_u_ext,
                          numpy.ndarray csrColumnOffsets_eb_u_u,
                          numpy.ndarray ebqe_penalty,
                          double adjoint_sigma):
        """
        Optimized jacobian calculation
        """
        cdef numpy.ndarray rowptr,colind,globalJacobian_a
        (rowptr,colind,globalJacobian_a) = globalJacobian.getCSRrepresentation()
        self.thisptr.calculateJacobian(<double*> mesh_trial_ref.data,
                                       <double*> mesh_grad_trial_ref.data,
                                       <double*> mesh_dof.data,
                                       <int*> mesh_l2g.data,
                                       <double*> dV_ref.data,
                                       <double*> u_trial_ref.data,
                                       <double*> u_grad_trial_ref.data,
                                       <double*> u_test_ref.data,
                                       <double*> u_grad_test_ref.data,
                                       <double*> elementDiameter.data,
                                       <double*> cfl.data,
                                       CT_sge,
                                       sc_uref,
                                       sc_alpha,
                                       useMetrics,
                                       <double*> mesh_trial_trace_ref.data,
                                       <double*> mesh_grad_trial_trace_ref.data,
                                       <double*> dS_ref.data,
                                       <double*> u_trial_trace_ref.data,
                                       <double*> u_grad_trial_trace_ref.data,
                                       <double*> u_test_trace_ref.data,
                                       <double*> u_grad_test_trace_ref.data,
                                       <double*> normal_ref.data,
                                       <double*> boundaryJac_ref.data,
                                       nElements_global,
                                       <int*> u_l2g.data,
                                       <double*> u_dof.data,
                                       <int*> sdInfo_u_u_rowptr.data,
                                       <int*> sdInfo_u_u_colind.data,
                                       <double*> q_a.data,
                                       <double*> q_v.data,
				       <double*> q_r.data,
                                       lag_shockCapturing,
                                       shockCapturingDiffusion,
                                       <double*> q_numDiff_u.data,
                                       <double*> q_numDiff_u_last.data,
                                       <int*> csrRowIndeces_u_u.data,
                                       <int*> csrColumnOffsets_u_u.data,
                                       <double*> globalJacobian_a.data,
                                       nExteriorElementBoundaries_global,
                                       <int*> exteriorElementBoundariesArray.data,
                                       <int*> elementBoundaryElementsArray.data,
                                       <int*> elementBoundaryLocalElementBoundariesArray.data,
                                       <double*> ebqe_a.data,
                                       <double*> ebqe_v.data,
                                       <int*> isDOFBoundary_u.data,
                                       <double*> ebqe_bc_u_ext.data,
                                       <int*> isFluxBoundary_u.data,
                                       <int*> isAdvectiveFluxBoundary_u.data,
                                       <double*> ebqe_bc_flux_u_ext.data,
                                       <double*> ebqe_bc_advectiveFlux_u_ext.data,
                                       <int*> csrColumnOffsets_eb_u_u.data,
                                       <double*> ebqe_penalty.data,
                                       adjoint_sigma)


class SubgridError(SGE_base):
    """
    SubgridError approximation for ADR equations

    .. inheritance-diagram:: SubgridError
       :parts: 2
    """

    def __init__(self,coefficients,nd):
        SGE_base.__init__(self,coefficients,nd,lag=False)
    def initializeElementQuadrature(self,mesh,t,cq):
        pass
    def updateSubgridErrorHistory(self,initializationPhase=False):
        pass
    def calculateSubgridError(self,q):
        pass

class ShockCapturing(ShockCapturing_base):
    """
    Residual-based shock capturing for ADR equations

    .. inheritance-diagram:: ShockCapturing
       :parts: 2
    """
    def __init__(self,coefficients,nd,shockCapturingFactor=0.25,lag=True,nStepsToDelay=None):
        ShockCapturing_base.__init__(self,coefficients,nd,shockCapturingFactor,lag)
        self.nStepsToDelay = nStepsToDelay
        self.nSteps=0
        if self.lag:
            log("ADR.ShockCapturing: lagging requested but must lag the first step; switching lagging off and delaying")
            self.nStepsToDelay=1
            self.lag=False
    def initializeElementQuadrature(self,mesh,t,cq):
        self.mesh=mesh
        self.numDiff=[]
        self.numDiff_last=[]
        for ci in range(self.nc):
            self.numDiff.append(cq[('numDiff',ci,ci)])
            self.numDiff_last.append(cq[('numDiff',ci,ci)])
    def updateShockCapturingHistory(self):
        self.nSteps += 1
        if self.lag:
            for ci in range(self.nc):
                self.numDiff_last[ci][:] = self.numDiff[ci]
        if self.lag == False and self.nStepsToDelay is not None and self.nSteps > self.nStepsToDelay:
            log("ADR.ShockCapturing: switched to lagged shock capturing")
            self.lag = True
            self.numDiff_last=[]
            for ci in range(self.nc):
                self.numDiff_last.append(self.numDiff[ci].copy())
        log("VOF: max numDiff %e" % (globalMax(self.numDiff_last[0].max()),))

class NumericalFlux_IIPG(proteus.NumericalFlux.Advection_DiagonalUpwind_Diffusion_IIPG_exterior):
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        proteus.NumericalFlux.Advection_DiagonalUpwind_Diffusion_IIPG_exterior.__init__(self,vt,getPointwiseBoundaryConditions,
                                                                                        getAdvectiveFluxBoundaryConditions,
                                                                                        getDiffusiveFluxBoundaryConditions)

class NumericalFlux_SIPG(proteus.NumericalFlux.Advection_DiagonalUpwind_Diffusion_SIPG_exterior):
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        proteus.NumericalFlux.Advection_DiagonalUpwind_Diffusion_SIPG_exterior.__init__(self,vt,getPointwiseBoundaryConditions,
                                                                                        getAdvectiveFluxBoundaryConditions,
                                                                                        getDiffusiveFluxBoundaryConditions)

class NumericalFlux_NIPG(proteus.NumericalFlux.Advection_DiagonalUpwind_Diffusion_NIPG_exterior):
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        proteus.NumericalFlux.Advection_DiagonalUpwind_Diffusion_NIPG_exterior.__init__(self,vt,getPointwiseBoundaryConditions,
                                                                                        getAdvectiveFluxBoundaryConditions,
                                                                                        getDiffusiveFluxBoundaryConditions)

NumericalFlux = NumericalFlux_SIPG

class Coefficients(TC_base):
    """
    Coefficients of linear ADR equations

    .. inheritance-diagram:: Coefficients
       :parts: 2
    """
    from proteus.ctransportCoefficients import L2projectEvaluate
    def __init__(self,aOfX,fOfX,velocity=None,nc=1,nd=2,l2proj=None,
                 timeVaryingCoefficients=False,
                 forceStrongDirichlet=False,
                 useMetrics=0.0,
                 sc_uref=1.0,
                 sc_beta=1.0):
        self.useMetrics = useMetrics
        self.forceStrongDirichlet=forceStrongDirichlet
        self.aOfX = aOfX
        self.fOfX = fOfX
        self.velocity=velocity
        self.nd = nd
        self.l2proj = l2proj
        self.timeVaryingCoefficients=timeVaryingCoefficients
        self.sc_uref=sc_uref
        self.sc_beta=sc_beta
        mass = {}
        advection = {}
        diffusion = {}
        potential = {}
        reaction  = {}
        hamiltonian = {}
        for i in range(nc):
            diffusion[i] = {i : {i:'constant'}}
            reaction[i]  = {i : 'constant'}
            advection[i] = {i : 'linear'} #now include for gravity type terms
            potential[i] = {i : 'u'}
        #end i
        sdInfo = {(0,0):(numpy.arange(start=0,stop=self.nd**2+1,step=self.nd,dtype='i'),
                         numpy.array([range(self.nd) for row in range(self.nd)],dtype='i'))}
        TC_base.__init__(self,
                         nc,
                         mass,
                         advection,
                         diffusion,
                         potential,
                         reaction,
                         hamiltonian,
                         variableNames=['u'],
                         sparseDiffusionTensors=sdInfo,
                         useSparseDiffusion=True,
                         movingDomain=False)
    def initializeElementQuadrature(self,t,cq):
        nd = self.nd
        for ci in range(self.nc):
            if cq.has_key(('df',ci,ci)):
                if self.velocity is not None:
                    cq[('df',ci,ci)][...,:] = self.velocity
                else:
                    cq[('df',ci,ci)].flat[:] = 0.0
            for i in range(len(cq[('r',ci)].flat)):
                cq[('r',ci)].flat[i] = -self.fOfX[ci](cq['x'].flat[3*i:3*(i+1)])
                cq[('a',ci,ci)].flat[nd*nd*i:nd*nd*(i+1)] = self.aOfX[ci](cq['x'].flat[3*i:3*(i+1)]).flat
    def initializeElementBoundaryQuadrature(self,t,cebq,cebq_global):
        nd = self.nd
        for c in [cebq,cebq_global]:
            for ci in range(self.nc):
                if c.has_key(('df',ci,ci)):
                    if self.velocity is not None:
                        c[('df',ci,ci)][...,:] = self.velocity
                    else:
                        c[('df',ci,ci)].flat[:] = 0.0
                if c.has_key(('r',ci)) and c.has_key(('a',ci,ci)):
                    for i in range(len(c[('u',ci)].flat)):
                        c[('r',ci)].flat[i] = -self.fOfX[ci](c['x'].flat[3*i:3*(i+1)])
                        c[('a',ci,ci)].flat[nd*nd*i:nd*nd*(i+1)] = self.aOfX[ci](c['x'].flat[3*i:3*(i+1)]).flat
    def initializeGlobalExteriorElementBoundaryQuadrature(self,t,cebqe):
        nd = self.nd
        for c in [cebqe]:
            for ci in range(self.nc):
                if c.has_key(('df',ci,ci)):
                    if self.velocity is not None:
                        c[('df',ci,ci)][...,:] = self.velocity
                    else:
                        c[('df',ci,ci)].flat[:] = 0.0
                if c.has_key(('r',ci)) and c.has_key(('a',ci,ci)):
                    for i in range(len(c[('u',ci)].flat)):
                        c[('r',ci)].flat[i] = -self.fOfX[ci](c['x'].flat[3*i:3*(i+1)])
                        c[('a',ci,ci)].flat[nd*nd*i:nd*nd*(i+1)] = self.aOfX[ci](c['x'].flat[3*i:3*(i+1)]).flat
    def evaluate(self,t,c):
        if self.timeVaryingCoefficients:
            nd = self.nd
            for ci in range(self.nc):
                if self.velocity is not None:
                    c[('df',ci,ci)][...,:] = self.velocity
                else:
                    c[('df',ci,ci)].flat[:] = 0.0
                for i in range(len(c[('r',ci)].flat)):
                    c[('r',ci)].flat[i] = -self.fOfX[ci](c['x'].flat[3*i:3*(i+1)])
                    c[('a',ci,ci)].flat[nd*nd*i:nd*nd*(i+1)] = self.aOfX[ci](c['x'].flat[3*i:3*(i+1)]).flat

class LevelModel(proteus.Transport.OneLevelTransport):
    """
    Optimized LevelModel for ADR equations

    .. inheritance-diagram:: LevelModel
       :parts: 2
    """
    nCalls=0
    def __init__(self,
                 uDict,
                 phiDict,
                 testSpaceDict,
                 matType,
                 dofBoundaryConditionsDict,
                 dofBoundaryConditionsSetterDict,
                 coefficients,
                 elementQuadrature,
                 elementBoundaryQuadrature,
                 fluxBoundaryConditionsDict=None,
                 advectiveFluxBoundaryConditionsSetterDict=None,
                 diffusiveFluxBoundaryConditionsSetterDictDict=None,
                 stressTraceBoundaryConditionsSetterDict=None,
                 stabilization=None,
                 shockCapturing=None,
                 conservativeFluxDict=None,
                 numericalFluxType=None,
                 TimeIntegrationClass=None,
                 massLumping=False,
                 reactionLumping=False,
                 options=None,
                 name='defaultName',
                 reuse_trial_and_test_quadrature=True,
                 sd = True,
                 movingDomain=False):#,
        from proteus import Comm
        #
        #set the objects describing the method and boundary conditions
        #
        self.movingDomain=movingDomain
        self.tLast_mesh=None
        #
        self.name=name
        self.sd=sd
        self.Hess=False
        self.lowmem=True
        self.timeTerm=True#allow turning off  the  time derivative
        #self.lowmem=False
        self.testIsTrial=True
        self.phiTrialIsTrial=True
        self.u = uDict
        self.ua = {}#analytical solutions
        self.phi  = phiDict
        self.dphi={}
        self.matType = matType
        #mwf try to reuse test and trial information across components if spaces are the same
        self.reuse_test_trial_quadrature = reuse_trial_and_test_quadrature#True#False
        if self.reuse_test_trial_quadrature:
            for ci in range(1,coefficients.nc):
                assert self.u[ci].femSpace.__class__.__name__ == self.u[0].femSpace.__class__.__name__, "to reuse_test_trial_quad all femSpaces must be the same!"
        ## Simplicial Mesh
        self.mesh = self.u[0].femSpace.mesh #assume the same mesh for  all components for now
        self.testSpace = testSpaceDict
        self.dirichletConditions = dofBoundaryConditionsDict
        self.dirichletNodeSetList=None #explicit Dirichlet  conditions for now, no Dirichlet BC constraints
        self.coefficients = coefficients
        self.coefficients.initializeMesh(self.mesh)
        self.nc = self.coefficients.nc
        self.stabilization = stabilization
        self.shockCapturing = shockCapturing
        self.conservativeFlux = conservativeFluxDict #no velocity post-processing for now
        self.fluxBoundaryConditions=fluxBoundaryConditionsDict
        self.advectiveFluxBoundaryConditionsSetterDict=advectiveFluxBoundaryConditionsSetterDict
        self.diffusiveFluxBoundaryConditionsSetterDictDict = diffusiveFluxBoundaryConditionsSetterDictDict
        #determine whether  the stabilization term is nonlinear
        self.stabilizationIsNonlinear = False
        #cek come back
        if self.stabilization is not None:
            for ci in range(self.nc):
                if coefficients.mass.has_key(ci):
                    for flag in coefficients.mass[ci].values():
                        if flag == 'nonlinear':
                            self.stabilizationIsNonlinear=True
                if  coefficients.advection.has_key(ci):
                    for  flag  in coefficients.advection[ci].values():
                        if flag == 'nonlinear':
                            self.stabilizationIsNonlinear=True
                if  coefficients.diffusion.has_key(ci):
                    for diffusionDict in coefficients.diffusion[ci].values():
                        for  flag  in diffusionDict.values():
                            if flag != 'constant':
                                self.stabilizationIsNonlinear=True
                if  coefficients.potential.has_key(ci):
                     for flag in coefficients.potential[ci].values():
                        if  flag == 'nonlinear':
                            self.stabilizationIsNonlinear=True
                if coefficients.reaction.has_key(ci):
                    for flag in coefficients.reaction[ci].values():
                        if  flag == 'nonlinear':
                            self.stabilizationIsNonlinear=True
                if coefficients.hamiltonian.has_key(ci):
                    for flag in coefficients.hamiltonian[ci].values():
                        if  flag == 'nonlinear':
                            self.stabilizationIsNonlinear=True
        #determine if we need element boundary storage
        self.elementBoundaryIntegrals = {}
        for ci  in range(self.nc):
            self.elementBoundaryIntegrals[ci] = ((self.conservativeFlux is not None) or
                                                 (numericalFluxType is not None) or
                                                 (self.fluxBoundaryConditions[ci] == 'outFlow') or
                                                 (self.fluxBoundaryConditions[ci] == 'mixedFlow') or
                                                 (self.fluxBoundaryConditions[ci] == 'setFlow'))
        #
        #calculate some dimensions
        #
        self.nSpace_global    = self.u[0].femSpace.nSpace_global #assume same space dim for all variables
        self.nDOF_trial_element     = [u_j.femSpace.max_nDOF_element for  u_j in self.u.values()]
        self.nDOF_phi_trial_element     = [phi_k.femSpace.max_nDOF_element for  phi_k in self.phi.values()]
        self.n_phi_ip_element = [phi_k.femSpace.referenceFiniteElement.interpolationConditions.nQuadraturePoints for  phi_k in self.phi.values()]
        self.nDOF_test_element     = [femSpace.max_nDOF_element for femSpace in self.testSpace.values()]
        self.nFreeDOF_global  = [dc.nFreeDOF_global for dc in self.dirichletConditions.values()]
        self.nVDOF_element    = sum(self.nDOF_trial_element)
        self.nFreeVDOF_global = sum(self.nFreeDOF_global)
        #
        NonlinearEquation.__init__(self,self.nFreeVDOF_global)
        #
        #build the quadrature point dictionaries from the input (this
        #is just for convenience so that the input doesn't have to be
        #complete)
        #
        elementQuadratureDict={}
        elemQuadIsDict = isinstance(elementQuadrature,dict)
        if elemQuadIsDict: #set terms manually
            for I in self.coefficients.elementIntegralKeys:
                if elementQuadrature.has_key(I):
                    elementQuadratureDict[I] = elementQuadrature[I]
                else:
                    elementQuadratureDict[I] = elementQuadrature['default']
        else:
            for I in self.coefficients.elementIntegralKeys:
                elementQuadratureDict[I] = elementQuadrature
        if self.stabilization is not None:
            for I in self.coefficients.elementIntegralKeys:
                if elemQuadIsDict:
                    if elementQuadrature.has_key(I):
                        elementQuadratureDict[('stab',)+I[1:]] = elementQuadrature[I]
                    else:
                        elementQuadratureDict[('stab',)+I[1:]] = elementQuadrature['default']
                else:
                    elementQuadratureDict[('stab',)+I[1:]] = elementQuadrature
        if self.shockCapturing is not None:
            for ci in self.shockCapturing.components:
                if elemQuadIsDict:
                    if elementQuadrature.has_key(('numDiff',ci,ci)):
                        elementQuadratureDict[('numDiff',ci,ci)] = elementQuadrature[('numDiff',ci,ci)]
                    else:
                        elementQuadratureDict[('numDiff',ci,ci)] = elementQuadrature['default']
                else:
                    elementQuadratureDict[('numDiff',ci,ci)] = elementQuadrature
        if massLumping:
            for ci in self.coefficients.mass.keys():
                elementQuadratureDict[('m',ci)] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global,1)
            for I in self.coefficients.elementIntegralKeys:
                elementQuadratureDict[('stab',)+I[1:]] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global,1)
        if reactionLumping:
            for ci in self.coefficients.mass.keys():
                elementQuadratureDict[('r',ci)] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global,1)
            for I in self.coefficients.elementIntegralKeys:
                elementQuadratureDict[('stab',)+I[1:]] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global,1)
        elementBoundaryQuadratureDict={}
        if isinstance(elementBoundaryQuadrature,dict): #set terms manually
            for I in self.coefficients.elementBoundaryIntegralKeys:
                if elementBoundaryQuadrature.has_key(I):
                    elementBoundaryQuadratureDict[I] = elementBoundaryQuadrature[I]
                else:
                    elementBoundaryQuadratureDict[I] = elementBoundaryQuadrature['default']
        else:
            for I in self.coefficients.elementBoundaryIntegralKeys:
                elementBoundaryQuadratureDict[I] = elementBoundaryQuadrature
        #
        # find the union of all element quadrature points and
        # build a quadrature rule for each integral that has a
        # weight at each point in the union
        #mwf include tag telling me which indices are which quadrature rule?
        (self.elementQuadraturePoints,self.elementQuadratureWeights,
         self.elementQuadratureRuleIndeces) = Quadrature.buildUnion(elementQuadratureDict)
        self.nQuadraturePoints_element = self.elementQuadraturePoints.shape[0]
        self.nQuadraturePoints_global = self.nQuadraturePoints_element*self.mesh.nElements_global
        #
        #Repeat the same thing for the element boundary quadrature
        #
        (self.elementBoundaryQuadraturePoints,
         self.elementBoundaryQuadratureWeights,
         self.elementBoundaryQuadratureRuleIndeces) = Quadrature.buildUnion(elementBoundaryQuadratureDict)
        self.nElementBoundaryQuadraturePoints_elementBoundary = self.elementBoundaryQuadraturePoints.shape[0]
        self.nElementBoundaryQuadraturePoints_global = (self.mesh.nElements_global*
                                                        self.mesh.nElementBoundaries_element*
                                                        self.nElementBoundaryQuadraturePoints_elementBoundary)
        if type(self.u[0].femSpace) == C0_AffineLinearOnSimplexWithNodalBasis:
            if self.nSpace_global == 3:
                assert(self.nQuadraturePoints_element == 5)
            elif self.nSpace_global == 2:
                assert(self.nQuadraturePoints_element == 6)
            elif self.nSpace_global == 1:
                assert(self.nQuadraturePoints_element == 3)

            if self.nSpace_global == 3:
                assert(self.nElementBoundaryQuadraturePoints_elementBoundary == 4)
            elif self.nSpace_global == 2:
                assert(self.nElementBoundaryQuadraturePoints_elementBoundary == 4)
            elif self.nSpace_global == 1:
                assert(self.nElementBoundaryQuadraturePoints_elementBoundary == 1)

        #pdb.set_trace()
        #
        #simplified allocations for test==trial and also check if space is mixed or not
        #
        self.q={}
        self.ebq={}
        self.ebq_global={}
        self.ebqe={}
        self.phi_ip={}
        #mesh
        self.q['x'] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,3),'d')
        self.ebqe['x'] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,3),'d')
        self.q[('u',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('grad(u)',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nSpace_global),'d')
        self.q[('a',0,0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.coefficients.sdInfo[(0,0)][0][-1]),'d')
        self.q[('df',0,0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nSpace_global),'d')
        self.q[('r',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('cfl',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('numDiff',0,0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')

        self.ebqe['penalty'] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('u',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('diffusiveFlux_bc_flag',0,0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'i')
        self.ebqe[('diffusiveFlux_bc',0,0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('advectiveFlux_bc_flag',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'i')
        self.ebqe[('advectiveFlux_bc',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('grad(u)',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),'d')
        self.ebqe[('a',0,0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.coefficients.sdInfo[(0,0)][0][-1]),'d')
        self.ebqe[('df',0,0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),'d')
        self.ebqe[('r',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')

        self.points_elementBoundaryQuadrature= set()
        self.scalars_elementBoundaryQuadrature= set([('u',ci) for ci in range(self.nc)])
        self.vectors_elementBoundaryQuadrature= set()
        self.tensors_elementBoundaryQuadrature= set()
        log(memory("element and element boundary Jacobians","OneLevelTransport"),level=4)
        self.inflowBoundaryBC = {}
        self.inflowBoundaryBC_values = {}
        self.inflowFlux = {}
        for cj in range(self.nc):
             self.inflowBoundaryBC[cj] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,),'i')
             self.inflowBoundaryBC_values[cj] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nDOF_trial_element[cj]),'d')
             self.inflowFlux[cj] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.internalNodes = set(range(self.mesh.nNodes_global))
        #identify the internal nodes this is ought to be in mesh
        ##\todo move this to mesh
        for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
            ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
            eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
            ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
            for i in range(self.mesh.nNodes_element):
                if i != ebN_element:
                    I = self.mesh.elementNodesArray[eN_global,i]
                    self.internalNodes -= set([I])
        self.nNodes_internal = len(self.internalNodes)
        self.internalNodesArray=numpy.zeros((self.nNodes_internal,),'i')
        for nI,n in enumerate(self.internalNodes):
            self.internalNodesArray[nI]=n
        #
        del self.internalNodes
        self.internalNodes = None
        log("Updating local to global mappings",2)
        self.updateLocal2Global()
        log("Building time integration object",2)
        log(memory("inflowBC, internalNodes,updateLocal2Global","OneLevelTransport"),level=4)
        #mwf for interpolating subgrid error for gradients etc
        if self.stabilization and self.stabilization.usesGradientStabilization:
            self.timeIntegration = TimeIntegrationClass(self,integrateInterpolationPoints=True)
        else:
             self.timeIntegration = TimeIntegrationClass(self)

        if options is not None:
            self.timeIntegration.setFromOptions(options)
        log(memory("TimeIntegration","OneLevelTransport"),level=4)
        log("Calculating numerical quadrature formulas",2)
        self.calculateQuadrature()

        comm = Comm.get()
        self.comm=comm
        if comm.size() > 1:
            assert numericalFluxType is not None and numericalFluxType.useWeakDirichletConditions,"You must use a numerical flux to apply weak boundary conditions for parallel runs"

        self.setupFieldStrides()

        log(memory("stride+offset","OneLevelTransport"),level=4)
        if numericalFluxType is not None:
            if options is None or options.periodicDirichletConditions is None:
                self.numericalFlux = numericalFluxType(self,
                                                       dofBoundaryConditionsSetterDict,
                                                       advectiveFluxBoundaryConditionsSetterDict,
                                                       diffusiveFluxBoundaryConditionsSetterDictDict)
            else:
                self.numericalFlux = numericalFluxType(self,
                                                       dofBoundaryConditionsSetterDict,
                                                       advectiveFluxBoundaryConditionsSetterDict,
                                                       diffusiveFluxBoundaryConditionsSetterDictDict,
                                                       options.periodicDirichletConditions)
        else:
            self.numericalFlux = None
        #set penalty terms
        #cek todo move into numerical flux initialization
        if self.ebq_global.has_key('penalty'):
            for ebN in range(self.mesh.nElementBoundaries_global):
                for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                    self.ebq_global['penalty'][ebN,k] = self.numericalFlux.penalty_constant/(self.mesh.elementBoundaryDiametersArray[ebN]**self.numericalFlux.penalty_power)
        #penalty term
        #cek move  to Numerical flux initialization
        if self.ebqe.has_key('penalty'):
            for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
                ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
                for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                    self.ebqe['penalty'][ebNE,k] = self.numericalFlux.penalty_constant/self.mesh.elementBoundaryDiametersArray[ebN]**self.numericalFlux.penalty_power
        log(memory("numericalFlux","OneLevelTransport"),level=4)
        self.elementEffectiveDiametersArray  = self.mesh.elementInnerDiametersArray
        #use post processing tools to get conservative fluxes, None by default
        from proteus import PostProcessingTools
        self.velocityPostProcessor = PostProcessingTools.VelocityPostProcessingChooser(self)
        log(memory("velocity postprocessor","OneLevelTransport"),level=4)
        #helper for writing out data storage
        from proteus import Archiver
        self.elementQuadratureDictionaryWriter = Archiver.XdmfWriter()
        self.elementBoundaryQuadratureDictionaryWriter = Archiver.XdmfWriter()
        self.exteriorElementBoundaryQuadratureDictionaryWriter = Archiver.XdmfWriter()
        for ci,fbcObject  in self.fluxBoundaryConditionsObjectsDict.iteritems():
            self.ebqe[('advectiveFlux_bc_flag',ci)] = numpy.zeros(self.ebqe[('advectiveFlux_bc',ci)].shape,'i')
            for t,g in fbcObject.advectiveFluxBoundaryConditionsDict.iteritems():
                if self.coefficients.advection.has_key(ci):
                    self.ebqe[('advectiveFlux_bc',ci)][t[0],t[1]] = g(self.ebqe[('x')][t[0],t[1]],self.timeIntegration.t)
                    self.ebqe[('advectiveFlux_bc_flag',ci)][t[0],t[1]] = 1
            for ck,diffusiveFluxBoundaryConditionsDict in fbcObject.diffusiveFluxBoundaryConditionsDictDict.iteritems():
                self.ebqe[('diffusiveFlux_bc_flag',ck,ci)] = numpy.zeros(self.ebqe[('diffusiveFlux_bc',ck,ci)].shape,'i')
                for t,g in diffusiveFluxBoundaryConditionsDict.iteritems():
                    self.ebqe[('diffusiveFlux_bc',ck,ci)][t[0],t[1]] = g(self.ebqe[('x')][t[0],t[1]],self.timeIntegration.t)
                    self.ebqe[('diffusiveFlux_bc_flag',ck,ci)][t[0],t[1]] = 1
        self.numericalFlux.setDirichletValues(self.ebqe)
        if self.movingDomain:
            self.MOVING_DOMAIN=1.0
        else:
            self.MOVING_DOMAIN=0.0
        #cek hack
        self.movingDomain=False
        self.MOVING_DOMAIN=0.0
        if self.mesh.nodeVelocityArray is None:
            self.mesh.nodeVelocityArray = numpy.zeros(self.mesh.nodeArray.shape,'d')
        #cek/ido todo replace python loops in modules with optimized code if possible/necessary
        self.forceStrongConditions=coefficients.forceStrongDirichlet
        self.dirichletConditionsForceDOF = {}
        if self.forceStrongConditions:
            for cj in range(self.nc):
                self.dirichletConditionsForceDOF[cj] = DOFBoundaryConditions(self.u[cj].femSpace,dofBoundaryConditionsSetterDict[cj],weakDirichletConditions=False)
        compKernelFlag = 0
        self.adr = ADR(self.nSpace_global,
                        self.nQuadraturePoints_element,
                        self.u[0].femSpace.elementMaps.localFunctionSpace.dim,
                        self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.dim,
                        self.testSpace[0].referenceFiniteElement.localFunctionSpace.dim,
                        self.nElementBoundaryQuadraturePoints_elementBoundary,
                        compKernelFlag)
    #mwf these are getting called by redistancing classes,
    def calculateCoefficients(self):
        pass
    def getResidual(self,u,r):
        import pdb
        import copy
        """
        Calculate the element residuals and add in to the global residual
        """
        r.fill(0.0)
        #Load the unknowns into the finite element dof
        self.setUnknowns(u)


        #no flux boundary conditions
        self.adr.calculateResidual(#element
            self.u[0].femSpace.elementMaps.psi,
            self.u[0].femSpace.elementMaps.grad_psi,
            self.mesh.nodeArray,
            self.mesh.elementNodesArray,
            self.elementQuadratureWeights[('u',0)],
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            self.mesh.elementDiametersArray,
            self.q[('cfl',0)],
            self.shockCapturing.shockCapturingFactor,
	    self.coefficients.sc_uref,
	    self.coefficients.sc_beta,
	    self.coefficients.useMetrics,
            #element boundary
            self.u[0].femSpace.elementMaps.psi_trace,
            self.u[0].femSpace.elementMaps.grad_psi_trace,
            self.elementBoundaryQuadratureWeights[('u',0)],
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.elementMaps.boundaryNormals,
            self.u[0].femSpace.elementMaps.boundaryJacobians,
            #physics
            self.mesh.nElements_global,
            self.u[0].femSpace.dofMap.l2g,
            self.u[0].dof,
            self.coefficients.sdInfo[(0,0)][0],
            self.coefficients.sdInfo[(0,0)][1],
            self.q[('a',0,0)],
            self.q[('df',0,0)],
            self.q[('r',0)],
            self.shockCapturing.lag,
            self.shockCapturing.shockCapturingFactor,
            self.shockCapturing.numDiff[0],
            self.shockCapturing.numDiff_last[0],
            self.offset[0],
            self.stride[0],
            r,
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray,
            self.ebqe[('a',0,0)],
            self.ebqe[('df',0,0)],
            self.numericalFlux.isDOFBoundary[0],
            self.numericalFlux.ebqe[('u',0)],
            self.ebqe[('diffusiveFlux_bc_flag',0,0)],
            self.ebqe[('advectiveFlux_bc_flag',0)],
            self.ebqe[('diffusiveFlux_bc',0,0)],
            self.ebqe[('advectiveFlux_bc',0)],
            self.ebqe['penalty'],
            self.numericalFlux.boundaryAdjoint_sigma)
        log("Global residual",level=9,data=r)
        self.coefficients.massConservationError = fabs(globalSum(sum(r.flat[:self.mesh.nElements_owned])))
        log("   Mass Conservation Error",level=3,data=self.coefficients.massConservationError)
        self.nonlinear_function_evaluations += 1
    def getJacobian(self,jacobian):
        #import superluWrappers
        #import numpy
        import pdb
        cfemIntegrals.zeroJacobian_CSR(self.nNonzerosInJacobian,jacobian)
        self.adr.calculateJacobian(#element
            self.u[0].femSpace.elementMaps.psi,
            self.u[0].femSpace.elementMaps.grad_psi,
            self.mesh.nodeArray,
            self.mesh.elementNodesArray,
            self.elementQuadratureWeights[('u',0)],
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            self.mesh.elementDiametersArray,
            self.q[('cfl',0)],
            self.shockCapturing.shockCapturingFactor,
	    self.coefficients.sc_uref,
	    self.coefficients.sc_beta,
	    self.coefficients.useMetrics,
            #element boundary
            self.u[0].femSpace.elementMaps.psi_trace,
            self.u[0].femSpace.elementMaps.grad_psi_trace,
            self.elementBoundaryQuadratureWeights[('u',0)],
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.elementMaps.boundaryNormals,
            self.u[0].femSpace.elementMaps.boundaryJacobians,
            self.mesh.nElements_global,
            self.u[0].femSpace.dofMap.l2g,
            self.u[0].dof,
            self.coefficients.sdInfo[(0,0)][0],
            self.coefficients.sdInfo[(0,0)][1],
            self.q[('a',0,0)],
            self.q[('df',0,0)],
            self.q[('r',0)],
            self.shockCapturing.lag,
            self.shockCapturing.shockCapturingFactor,
            self.shockCapturing.numDiff[0],
            self.shockCapturing.numDiff_last[0],
            self.csrRowIndeces[(0,0)],
            self.csrColumnOffsets[(0,0)],
            jacobian,
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray,
            self.ebqe[('a',0,0)],
            self.ebqe[('df',0,0)],
            self.numericalFlux.isDOFBoundary[0],
            self.numericalFlux.ebqe[('u',0)],
            self.ebqe[('diffusiveFlux_bc_flag',0,0)],
            self.ebqe[('advectiveFlux_bc_flag',0)],
            self.ebqe[('diffusiveFlux_bc',0,0)],
            self.ebqe[('advectiveFlux_bc',0)],
            self.csrColumnOffsets_eb[(0,0)],
            self.ebqe['penalty'],
            self.numericalFlux.boundaryAdjoint_sigma)
        log("Jacobian ",level=10,data=jacobian)
        #mwf decide if this is reasonable for solver statistics
        self.nonlinear_function_jacobian_evaluations += 1
        return jacobian
    def calculateElementQuadrature(self):
        """
        Calculate the physical location and weights of the quadrature rules
        and the shape information at the quadrature points.

        This function should be called only when the mesh changes.
        """
        self.u[0].femSpace.elementMaps.getValues(self.elementQuadraturePoints,
                                                 self.q['x'])
        self.u[0].femSpace.elementMaps.getBasisValuesRef(self.elementQuadraturePoints)
        self.u[0].femSpace.elementMaps.getBasisGradientValuesRef(self.elementQuadraturePoints)
        self.u[0].femSpace.getBasisValuesRef(self.elementQuadraturePoints)
        self.u[0].femSpace.getBasisGradientValuesRef(self.elementQuadraturePoints)
        self.coefficients.initializeElementQuadrature(self.timeIntegration.t,self.q)
        if self.stabilization is not None:
            self.stabilization.initializeElementQuadrature(self.mesh,self.timeIntegration.t,self.q)
            self.stabilization.initializeTimeIntegration(self.timeIntegration)
        if self.shockCapturing is not None:
            self.shockCapturing.initializeElementQuadrature(self.mesh,self.timeIntegration.t,self.q)
    def calculateElementBoundaryQuadrature(self):
        pass
    def calculateExteriorElementBoundaryQuadrature(self):
        self.u[0].femSpace.elementMaps.getBasisValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.elementMaps.getBasisGradientValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.getBasisValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.getBasisGradientValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.elementMaps.getValuesGlobalExteriorTrace(self.elementBoundaryQuadraturePoints,
                                                                    self.ebqe['x'])
        self.fluxBoundaryConditionsObjectsDict = dict([(cj,FluxBoundaryConditions(self.mesh,
                                                                                  self.nElementBoundaryQuadraturePoints_elementBoundary,
                                                                                  self.ebqe[('x')],
                                                                                  self.advectiveFluxBoundaryConditionsSetterDict[cj],
                                                                                  self.diffusiveFluxBoundaryConditionsSetterDictDict[cj]))
                                                       for cj in self.advectiveFluxBoundaryConditionsSetterDict.keys()])
        self.coefficients.initializeGlobalExteriorElementBoundaryQuadrature(self.timeIntegration.t,self.ebqe)
    def estimate_mt(self):
        pass
    def calculateAuxiliaryQuantitiesAfterStep(self):
        pass
    def calculateSolutionAtQuadrature(self):
        pass

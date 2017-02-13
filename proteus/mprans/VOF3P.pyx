# A type of -*- python -*- file
"""
An optimized volume-of-fluid  transport module
"""
import numpy
cimport numpy
from math import fabs
import proteus
from proteus import LinearAlgebraTools 
from proteus import cfemIntegrals, Quadrature, Norms, Comm
from proteus.NonlinearSolvers import NonlinearEquation
from proteus.FemTools import (DOFBoundaryConditions,
                              FluxBoundaryConditions,
                              C0_AffineLinearOnSimplexWithNodalBasis)
from proteus.flcbdfWrappers import globalMax
from proteus.Profiling import memory
from proteus.Profiling import logEvent as log
from proteus.Transport import OneLevelTransport
from proteus.TransportCoefficients import TC_base
from proteus.SubgridError import SGE_base
from proteus.ShockCapturing import ShockCapturing_base

cdef extern from "mprans/VOF3P.h" namespace "proteus":
    cdef cppclass cppVOF3P_base:
        void FCTStep(double dt, 
	             int NNZ,
		     int numDOFs,
		     double* lumped_mass_matrix, 
		     double* soln, 
		     double* solH, 
		     double* flux_plus_dLij_times_soln, 
		     int* csrRowIndeces_DofLoops, 
		     int* csrColumnOffsets_DofLoops, 
		     double* MassMatrix, 
		     double* dL_minus_dC,
                     double* min_u_bc,
                     double* max_u_bc) 
        void getInflowDOFs(double* mesh_dof,            
                           int* mesh_l2g,
                           double* mesh_trial_trace_ref,
                           double* mesh_grad_trial_trace_ref,
                           double* normal_ref,
                           double* boundaryJac_ref,
                           int* u_l2g, 
                           int offset_u, 
                           int stride_u, 
                           int nExteriorElementBoundaries_global,
                           int* exteriorElementBoundariesArray,
                           int* elementBoundaryElementsArray,
                           int* elementBoundaryLocalElementBoundariesArray,
                           double* ebqe_velocity_ext,
                           double* inflow_DOFs)
        void calculateResidual(double * mesh_trial_ref,
                               double * mesh_grad_trial_ref,
                               double * mesh_dof,
                               double * meshVelocity_dof,
                               double MOVING_DOMAIN,
                               int * mesh_l2g,
                               double * dV_ref,
                               double * u_trial_ref,
                               double * u_grad_trial_ref,
                               double * u_test_ref,
                               double * u_grad_test_ref,
                               int nDOF_vel_trial_element, 
                               int * vel_l2g,
                               double * vel_grad_trial_ref,
                               double * mesh_trial_trace_ref,
                               double * mesh_grad_trial_trace_ref,
                               double * dS_ref,
                               double * u_trial_trace_ref,
                               double * u_grad_trial_trace_ref,
                               double * u_test_trace_ref,
                               double * u_grad_test_trace_ref,
                               double * normal_ref,
                               double * boundaryJac_ref,
                               int nElements_global,
                               double useMetrics,
                               double alphaBDF,
                               int lag_shockCapturing,
                               double shockCapturingDiffusion,
                               double sc_uref, double sc_alpha,
                               const double * q_vos,
                               int * u_l2g,
                               double * elementDiameter,
                               double * u_dof, 
                               double * u_dof_old,
                               double * u_dof_old_old,
                               double * velx_tn_dof,
                               double * vely_tn_dof,
                               double * velz_tn_dof,
                               double * velocity,
                               double * div_velocity,
                               double * q_m,
                               double * q_u,
                               double * q_m_betaBDF,
                               double * q_dV,
                               double * q_dV_last,
                               double * cfl,
                               double * q_numDiff_u,
                               double * q_numDiff_u_last,
                               int offset_u, int stride_u,
                               double * globalResidual,
                               int nExteriorElementBoundaries_global,
                               int * exteriorElementBoundariesArray,
                               int * elementBoundaryElementsArray,
                               int * elementBoundaryLocalElementBoundariesArray,
                               double * ebqe_velocity_ext,
                               double * ebqe_div_velocity_ext,
                               const double * ebqe_vos_ext,
                               int * isDOFBoundary_u,
                               double * ebqe_bc_u_ext,
                               int * isFluxBoundary_u,
                               double * ebqe_bc_flux_u_ext,
                               double * ebqe_phi, double epsFact,
                               double * ebqe_u,
                               double * ebqe_flux,
                               int EDGE_VISCOSITY,
                               int ENTROPY_VISCOSITY, 
                               double cE,
                               double cMax, 
                               double cK,
                               double uL, 
                               double uR, 
                               int numDOFs,                              
                               int NNZ,
                               int* csrRowIndeces_DofLoops,
                               int* csrColumnOffsets_DofLoops,
                               int* csrRowIndeces_CellLoops,
                               int* csrColumnOffsets_CellLoops,
                               int * csrColumnOffsets_eb_CellLoops,
                               double * Cx,
                               double * Cy,
                               double * Cz,
                               double * CTx,
                               double * CTy,
                               double * CTz,
                               int POWER_SMOOTHNESS_INDICATOR, 
                               int LUMPED_MASS_MATRIX, 
                               double * flux_plus_dLij_times_soln,
                               double * dL_minus_dC,
                               double * min_u_bc,
                               double * max_u_bc,
                               double * quantDOFs)
        void calculateJacobian(double * mesh_trial_ref,
                               double * mesh_grad_trial_ref,
                               double * mesh_dof,
                               double * mesh_velocity_dof,
                               double MOVING_DOMAIN,
                               int * mesh_l2g,
                               double * dV_ref,
                               double * u_trial_ref,
                               double * u_grad_trial_ref,
                               double * u_test_ref,
                               double * u_grad_test_ref,
                               double * mesh_trial_trace_ref,
                               double * mesh_grad_trial_trace_ref,
                               double * dS_ref,
                               double * u_trial_trace_ref,
                               double * u_grad_trial_trace_ref,
                               double * u_test_trace_ref,
                               double * u_grad_test_trace_ref,
                               double * normal_ref,
                               double * boundaryJac_ref,
                               int nElements_global,
                               double useMetrics,
                               double alphaBDF,
                               int lag_shockCapturing,
                               double shockCapturingDiffusion,
                               const double * q_vos,
                               int * u_l2g,
                               double * elementDiameter,
                               double * u_dof,
                               double * velocity,
                               double * q_m_betaBDF,
                               double * cfl,
                               double * q_numDiff_u_last,
                               int * csrRowIndeces_u_u,
                               int * csrColumnOffsets_u_u,
                               double * globalJacobian,
                               int nExteriorElementBoundaries_global,
                               int * exteriorElementBoundariesArray,
                               int * elementBoundaryElementsArray,
                               int * elementBoundaryLocalElementBoundariesArray,
                               double * ebqe_velocity_ext,
                               const double * ebqe_vos_ext,
                               int * isDOFBoundary_u,
                               double * ebqe_bc_u_ext,
                               int * isFluxBoundary_u,
                               double * ebqe_bc_flux_u_ext,
                               int * csrColumnOffsets_eb_u_u,
                               int EDGE_VISCOSITY,
                               int ENTROPY_VISCOSITY,
                               int LUMPED_MASS_MATRIX)
    cppVOF3P_base* newVOF3P(int nSpaceIn,
                            int nQuadraturePoints_elementIn,
                            int nDOF_mesh_trial_elementIn,
                            int nDOF_trial_elementIn,
                            int nDOF_test_elementIn,
                            int nQuadraturePoints_elementBoundaryIn,
                            int CompKernelFlag)

cdef class VOF3P:
    """
    Optimized VOF3P member functions
    """

    cdef cppVOF3P_base* thisptr
    def __cinit__(self,
                  int nSpaceIn,
                  int nQuadraturePoints_elementIn,
                  int nDOF_mesh_trial_elementIn,
                  int nDOF_trial_elementIn,
                  int nDOF_test_elementIn,
                  int nQuadraturePoints_elementBoundaryIn,
                  int CompKernelFlag):
        self.thisptr = newVOF3P(nSpaceIn,
                                nQuadraturePoints_elementIn,
                                nDOF_mesh_trial_elementIn,
                                nDOF_trial_elementIn,
                                nDOF_test_elementIn,
                                nQuadraturePoints_elementBoundaryIn,
                                CompKernelFlag)
    def __dealloc__(self):
        del self.thisptr
    def FCTStep(self, 
                double dt, 
                int NNZ,
                int numDOFs,
                numpy.ndarray lumped_mass_matrix, 
                numpy.ndarray soln, 
                numpy.ndarray solH, 
                numpy.ndarray flux_plus_dLij_times_soln, 
                numpy.ndarray csrRowIndeces_DofLoops, 
                numpy.ndarray csrColumnOffsets_DofLoops, 
                numpy.ndarray MassMatrix, 
                numpy.ndarray dL_minus_dC,
                numpy.ndarray min_u_bc,
                numpy.ndarray max_u_bc):
        self.thisptr.FCTStep(dt, 
                             NNZ,
                             numDOFs,
                             <double*> lumped_mass_matrix.data, 
                             <double*> soln.data, 
                             <double*> solH.data,
                             <double*> flux_plus_dLij_times_soln.data,
                             <int*> csrRowIndeces_DofLoops.data,
                             <int*> csrColumnOffsets_DofLoops.data,
                             <double*> MassMatrix.data,
                             <double*> dL_minus_dC.data,
                             <double*> min_u_bc.data,
                             <double*> max_u_bc.data)
    def getInflowDOFs(self, 
                      numpy.ndarray mesh_dof,            
                      numpy.ndarray mesh_l2g,
                      numpy.ndarray mesh_trial_trace_ref,
                      numpy.ndarray mesh_grad_trial_trace_ref,
                      numpy.ndarray normal_ref,
                      numpy.ndarray boundaryJac_ref,
                      numpy.ndarray u_l2g, 
                      int offset_u, 
                      int stride_u, 
                      int nExteriorElementBoundaries_global,
                      numpy.ndarray exteriorElementBoundariesArray,
                      numpy.ndarray elementBoundaryElementsArray,
                      numpy.ndarray elementBoundaryLocalElementBoundariesArray,
                      numpy.ndarray ebqe_velocity_ext,
                      numpy.ndarray inflow_DOFs):
        self.thisptr.getInflowDOFs(<double*> mesh_dof.data,
                                    <int*> mesh_l2g.data,
                                    <double*> mesh_trial_trace_ref.data,
                                    <double*> mesh_grad_trial_trace_ref.data,
                                    <double*> normal_ref.data,
                                    <double*> boundaryJac_ref.data,
                                    <int*> u_l2g.data,
                                    offset_u, 
                                    stride_u, 
                                    nExteriorElementBoundaries_global,
                                    <int*> exteriorElementBoundariesArray.data,
                                    <int*> elementBoundaryElementsArray.data,
                                    <int*> elementBoundaryLocalElementBoundariesArray.data,
                                    <double*> ebqe_velocity_ext.data,
                                    <double*> inflow_DOFs.data)
    def calculateResidual(self,
                          numpy.ndarray mesh_trial_ref,
                          numpy.ndarray mesh_grad_trial_ref,
                          numpy.ndarray mesh_dof,
                          numpy.ndarray meshVelocity_dof,
                          double MOVING_DOMAIN,
                          numpy.ndarray mesh_l2g,
                          numpy.ndarray dV_ref,
                          numpy.ndarray u_trial_ref,
                          numpy.ndarray u_grad_trial_ref,
                          numpy.ndarray u_test_ref,
                          numpy.ndarray u_grad_test_ref,
                          int nDOF_vel_trial_element, 
                          numpy.ndarray vel_l2g,
                          numpy.ndarray vel_grad_trial_ref,
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
                          double useMetrics,
                          double alphaBDF,
                          int lag_shockCapturing,
                          double shockCapturingDiffusion,
                          double sc_uref, double sc_alpha,
                          numpy.ndarray q_vos,
                          numpy.ndarray u_l2g,
                          numpy.ndarray elementDiameter,
                          numpy.ndarray u_dof, 
                          numpy.ndarray u_dof_old,
                          numpy.ndarray u_dof_old_old,
                          numpy.ndarray velx_tn_dof,
                          numpy.ndarray vely_tn_dof,
                          numpy.ndarray velz_tn_dof,
                          numpy.ndarray velocity,
                          numpy.ndarray div_velocity,
                          numpy.ndarray q_m,
                          numpy.ndarray q_u,
                          numpy.ndarray q_m_betaBDF,
                          numpy.ndarray q_dV,
                          numpy.ndarray q_dV_last,
                          numpy.ndarray cfl,
                          numpy.ndarray q_numDiff_u,
                          numpy.ndarray q_numDiff_u_last,
                          int offset_u, int stride_u,
                          numpy.ndarray globalResidual,
                          int nExteriorElementBoundaries_global,
                          numpy.ndarray exteriorElementBoundariesArray,
                          numpy.ndarray elementBoundaryElementsArray,
                          numpy.ndarray elementBoundaryLocalElementBoundariesArray,
                          numpy.ndarray ebqe_velocity_ext,
                          numpy.ndarray ebqe_div_velocity_ext,
                          numpy.ndarray ebqe_vos_ext,
                          numpy.ndarray isDOFBoundary_u,
                          numpy.ndarray ebqe_bc_u_ext,
                          numpy.ndarray isFluxBoundary_u,
                          numpy.ndarray ebqe_bc_flux_u_ext,
                          numpy.ndarray ebqe_phi, double epsFact,
                          numpy.ndarray ebqe_u,
                          numpy.ndarray ebqe_flux, 
                          int EDGE_VISCOSITY, 
                          int ENTROPY_VISCOSITY, 
                          double cE, 
                          double cMax, 
                          double cK,
                          double uL, 
                          double uR,
                          int numDOFs,
                          int NNZ,
                          numpy.ndarray csrRowIndeces_DofLoops,
                          numpy.ndarray csrColumnOffsets_DofLoops,
                          numpy.ndarray csrRowIndeces_CellLoops,
                          numpy.ndarray csrColumnOffsets_CellLoops,
                          numpy.ndarray csrColumnOffsets_eb_CellLoops,
                          numpy.ndarray Cx,
                          numpy.ndarray Cy,
                          numpy.ndarray Cz,
                          numpy.ndarray CTx,
                          numpy.ndarray CTy,
                          numpy.ndarray CTz,
                          int POWER_SMOOTHNESS_INDICATOR, 
                          int LUMPED_MASS_MATRIX, 
                          numpy.ndarray flux_plus_dLij_times_soln,
                          numpy.ndarray dL_minus_dC, 
                          numpy.ndarray min_u_bc,
                          numpy.ndarray max_u_bc,
                          numpy.ndarray quantDOFs):
        self.thisptr.calculateResidual(<double*> mesh_trial_ref.data,
                                       <double*> mesh_grad_trial_ref.data,
                                       <double*> mesh_dof.data,
                                       <double*> meshVelocity_dof.data,
                                       MOVING_DOMAIN,
                                       <int*> mesh_l2g.data,
                                       <double*> dV_ref.data,
                                       <double*> u_trial_ref.data,
                                       <double*> u_grad_trial_ref.data,
                                       <double*> u_test_ref.data,
                                       <double*> u_grad_test_ref.data,
                                       nDOF_vel_trial_element,
                                       <int*> vel_l2g.data,
                                       <double*> vel_grad_trial_ref.data,
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
                                       useMetrics,
                                       alphaBDF,
                                       lag_shockCapturing,
                                       shockCapturingDiffusion,
                                       sc_uref,
                                       sc_alpha,
                                       <double*> q_vos.data,
                                       <int*> u_l2g.data,
                                       <double*> elementDiameter.data,
                                       <double*> u_dof.data,
                                       <double*> u_dof_old.data,
                                       <double*> u_dof_old_old.data,
                                       <double*> velx_tn_dof.data,
                                       <double*> vely_tn_dof.data,
                                       <double*> velz_tn_dof.data,
                                       <double*> velocity.data,
                                       <double*> div_velocity.data,
                                       <double*> q_m.data,
                                       <double*> q_u.data,
                                       <double*> q_m_betaBDF.data,
                                       <double*> q_dV.data,
                                       <double*> q_dV_last.data,
                                       <double*> cfl.data,
                                       <double*> q_numDiff_u.data,
                                       <double*> q_numDiff_u_last.data,
                                       offset_u,
                                       stride_u,
                                       <double*> globalResidual.data,
                                       nExteriorElementBoundaries_global,
                                       <int*> exteriorElementBoundariesArray.data,
                                       <int*> elementBoundaryElementsArray.data,
                                       <int*> elementBoundaryLocalElementBoundariesArray.data,
                                       <double*> ebqe_velocity_ext.data,
                                       <double*> ebqe_div_velocity_ext.data,
                                       <double*> ebqe_vos_ext.data,
                                       <int*> isDOFBoundary_u.data,
                                       <double*> ebqe_bc_u_ext.data,
                                       <int*> isFluxBoundary_u.data,
                                       <double*> ebqe_bc_flux_u_ext.data,
                                       <double*> ebqe_phi.data,
                                       epsFact,
                                       <double*> ebqe_u.data,
                                       <double*> ebqe_flux.data,
                                       EDGE_VISCOSITY,
                                       ENTROPY_VISCOSITY, 
                                       cE, 
                                       cMax, 
                                       cK,
                                       uL, 
                                       uR, 
                                       numDOFs,
                                       NNZ,
				       <int*> csrRowIndeces_DofLoops.data,
				       <int*> csrColumnOffsets_DofLoops.data,
				       <int*> csrRowIndeces_CellLoops.data,
				       <int*> csrColumnOffsets_CellLoops.data,
                                       <int*> csrColumnOffsets_eb_CellLoops.data,
                                       <double*> Cx.data,
                                       <double*> Cy.data,
                                       <double*> Cz.data,
                                       <double*> CTx.data,
                                       <double*> CTy.data,
                                       <double*> CTz.data,
                                       POWER_SMOOTHNESS_INDICATOR, 
                                       LUMPED_MASS_MATRIX, 
                                       <double*> flux_plus_dLij_times_soln.data,
                                       <double*> dL_minus_dC.data, 
                                       <double*> min_u_bc.data,
                                       <double*> max_u_bc.data,
                                       <double*> quantDOFs.data)
    def calculateJacobian(self,
                          numpy.ndarray mesh_trial_ref,
                          numpy.ndarray mesh_grad_trial_ref,
                          numpy.ndarray mesh_dof,
                          numpy.ndarray mesh_velocity_dof,
                          double MOVING_DOMAIN,
                          numpy.ndarray mesh_l2g,
                          numpy.ndarray dV_ref,
                          numpy.ndarray u_trial_ref,
                          numpy.ndarray u_grad_trial_ref,
                          numpy.ndarray u_test_ref,
                          numpy.ndarray u_grad_test_ref,
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
                          double useMetrics,
                          double alphaBDF,
                          int lag_shockCapturing,
                          double shockCapturingDiffusion,
                          numpy.ndarray q_vos,
                          numpy.ndarray u_l2g,
                          numpy.ndarray elementDiameter,
                          numpy.ndarray u_dof,
                          numpy.ndarray velocity,
                          numpy.ndarray q_m_betaBDF,
                          numpy.ndarray cfl,
                          numpy.ndarray q_numDiff_u_last,
                          numpy.ndarray csrRowIndeces_u_u,
                          numpy.ndarray csrColumnOffsets_u_u,
                          globalJacobian,
                          int nExteriorElementBoundaries_global,
                          numpy.ndarray exteriorElementBoundariesArray,
                          numpy.ndarray elementBoundaryElementsArray,
                          numpy.ndarray elementBoundaryLocalElementBoundariesArray,
                          numpy.ndarray ebqe_velocity_ext,
                          numpy.ndarray ebqe_vos_ext,
                          numpy.ndarray isDOFBoundary_u,
                          numpy.ndarray ebqe_bc_u_ext,
                          numpy.ndarray isFluxBoundary_u,
                          numpy.ndarray ebqe_bc_flux_u_ext,
                          numpy.ndarray csrColumnOffsets_eb_u_u,
                          int EDGE_VISCOSITY, 
                          int ENTROPY_VISCOSITY,
                          int LUMPED_MASS_MATRIX):
        """
        Optimized jacobian calculation
        """
        cdef numpy.ndarray rowptr,colind,globalJacobian_a
        (rowptr,colind,globalJacobian_a) = globalJacobian.getCSRrepresentation()
        self.thisptr.calculateJacobian(<double*> mesh_trial_ref.data,
                                       <double*> mesh_grad_trial_ref.data,
                                       <double*> mesh_dof.data,
                                       <double*> mesh_velocity_dof.data,
                                       MOVING_DOMAIN,
                                       <int*> mesh_l2g.data,
                                       <double*> dV_ref.data,
                                       <double*> u_trial_ref.data,
                                       <double*> u_grad_trial_ref.data,
                                       <double*> u_test_ref.data,
                                       <double*> u_grad_test_ref.data,
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
                                       useMetrics,
                                       alphaBDF,
                                       lag_shockCapturing,
                                       shockCapturingDiffusion,
                                       <double*> q_vos.data,
                                       <int*> u_l2g.data,
                                       <double*> elementDiameter.data,
                                       <double*> u_dof.data,
                                       <double*> velocity.data,
                                       <double*> q_m_betaBDF.data,
                                       <double*> cfl.data,
                                       <double*> q_numDiff_u_last.data,
                                       <int*> csrRowIndeces_u_u.data,
                                       <int*> csrColumnOffsets_u_u.data,
                                       <double*> globalJacobian_a.data,
                                       nExteriorElementBoundaries_global,
                                       <int*> exteriorElementBoundariesArray.data,
                                       <int*> elementBoundaryElementsArray.data,
                                       <int*> elementBoundaryLocalElementBoundariesArray.data,
                                       <double*> ebqe_velocity_ext.data,
                                       <double*> ebqe_vos_ext.data,
                                       <int*> isDOFBoundary_u.data,
                                       <double*> ebqe_bc_u_ext.data,
                                       <int*> isFluxBoundary_u.data,
                                       <double*> ebqe_bc_flux_u_ext.data,
                                       <int*> csrColumnOffsets_eb_u_u.data,
                                       EDGE_VISCOSITY, 
                                       ENTROPY_VISCOSITY,
                                       LUMPED_MASS_MATRIX)

class SubgridError(SGE_base):

    def __init__(self, coefficients, nd):
        proteus.SubgridError.SGE_base.__init__(
            self, coefficients, nd, lag=False)

    def initializeElementQuadrature(self, mesh, t, cq):
        pass

    def updateSubgridErrorHistory(self, initializationPhase=False):
        pass

    def calculateSubgridError(self, q):
        pass


class ShockCapturing(ShockCapturing_base):

    def __init__(
            self,
            coefficients,
            nd,
            shockCapturingFactor=0.25,
            lag=True,
            nStepsToDelay=None):
        proteus.ShockCapturing.ShockCapturing_base.__init__(
            self, coefficients, nd, shockCapturingFactor, lag)
        self.nStepsToDelay = nStepsToDelay
        self.nSteps = 0
        if self.lag:
            log("VOF3P.ShockCapturing: lagging requested but must lag the first step; switching lagging off and delaying")
            self.nStepsToDelay = 1
            self.lag = False

    def initializeElementQuadrature(self, mesh, t, cq):
        self.mesh = mesh
        self.numDiff = []
        self.numDiff_last = []
        for ci in range(self.nc):
            self.numDiff.append(cq[('numDiff', ci, ci)])
            self.numDiff_last.append(cq[('numDiff', ci, ci)])

    def updateShockCapturingHistory(self):
        self.nSteps += 1
        if self.lag:
            for ci in range(self.nc):
                self.numDiff_last[ci][:] = self.numDiff[ci]
        if self.lag == False and self.nStepsToDelay != None and self.nSteps > self.nStepsToDelay:
            log("VOF3P.ShockCapturing: switched to lagged shock capturing")
            self.lag = True
            self.numDiff_last = []
            for ci in range(self.nc):
                self.numDiff_last.append(self.numDiff[ci].copy())
        log("VOF3P: max numDiff %e" % (globalMax(self.numDiff_last[0].max()),))


class NumericalFlux(
        proteus.NumericalFlux.Advection_DiagonalUpwind_Diffusion_IIPG_exterior):

    def __init__(self, vt, getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        proteus.NumericalFlux.Advection_DiagonalUpwind_Diffusion_IIPG_exterior.__init__(
            self,
            vt,
            getPointwiseBoundaryConditions,
            getAdvectiveFluxBoundaryConditions,
            getDiffusiveFluxBoundaryConditions)


class Coefficients(proteus.TransportCoefficients.TC_base):
    from proteus.ctransportCoefficients import VOFCoefficientsEvaluate
    from proteus.UnstructuredFMMandFSWsolvers import FMMEikonalSolver, FSWEikonalSolver
    from proteus.NonlinearSolvers import EikonalSolver
    from proteus.ctransportCoefficients import VolumeAveragedVOFCoefficientsEvaluate
    from proteus.cfemIntegrals import copyExteriorElementBoundaryValuesFromElementBoundaryValues

    def __init__(
            self,
            EDGE_VISCOSITY=0,
            ENTROPY_VISCOSITY=0,
            FCT=0,
            POWER_SMOOTHNESS_INDICATOR=2,
            LUMPED_MASS_MATRIX=0,
            # FOR LOG BASED ENTROPY FUNCTION
            uL=0.0, 
            uR=1.0,
            # FOR ARTIFICIAL COMPRESSION
            cK=0.25,
            # FOR IMPOSING DIRICHLET BCs STRONGLY
            forceStrongConditions=0,
            # FOR ELEMENT BASED ENTROPY VISCOSITY
            cMax=0.1,
            cE=1.0,
            LS_model=None,
            V_model=0,
            RD_model=None,
            ME_model=1,
            VOS_model=None,
            EikonalSolverFlag=0,
            checkMass=True,
            epsFact=0.0,
            useMetrics=0.0,
            sc_uref=1.0,
            sc_beta=1.0,
            setParamsFunc=None,
            movingDomain=False,
            set_vos=None):
        self.set_vos=set_vos
        self.useMetrics = useMetrics
        self.variableNames = ['vof']
        nc = 1
        mass = {0: {0: 'linear'}}
        advection = {0: {0: 'linear'}}
        hamiltonian = {}
        diffusion = {}
        potential = {}
        reaction = {}
        TC_base.__init__(self,
                         nc,
                         mass,
                         advection,
                         diffusion,
                         potential,
                         reaction,
                         hamiltonian,
                         self.variableNames,
                         movingDomain=movingDomain)
        self.epsFact = epsFact
        self.V_model = V_model
        self.modelIndex = ME_model
        self.RD_modelIndex = RD_model
        self.LS_modelIndex = LS_model
        self.VOS_model = VOS_model
        self.sc_uref = sc_uref
        self.sc_beta = sc_beta
        # mwf added
        self.eikonalSolverFlag = EikonalSolverFlag
        if self.eikonalSolverFlag >= 1:  # FMM
            assert self.RD_modelIndex < 0, "no redistance with eikonal solver too"
        self.checkMass = checkMass
        # VRANS
        self.q_vos = None
        self.ebqe_vos = None
        self.setParamsFunc = setParamsFunc
        self.flowCoefficients = None
        self.movingDomain = movingDomain
        # EDGE BASED (AND ENTROPY) VISCOSITY 
        self.EDGE_VISCOSITY=EDGE_VISCOSITY
        self.ENTROPY_VISCOSITY=ENTROPY_VISCOSITY
        self.POWER_SMOOTHNESS_INDICATOR=POWER_SMOOTHNESS_INDICATOR
        self.LUMPED_MASS_MATRIX=LUMPED_MASS_MATRIX
        self.FCT=FCT
        self.uL=uL
        self.uR=uR
        self.cK=cK
        self.forceStrongConditions=forceStrongConditions
        self.cMax=cMax
        self.cE=cE        

    def initializeMesh(self, mesh):
        self.eps = self.epsFact * mesh.h

    def attachModels(self, modelList):
        # self
        self.model = modelList[self.modelIndex]

        self.u_dof_old = numpy.copy(self.model.u[0].dof)
        self.u_dof_old_old = numpy.copy(self.model.u[0].dof)

        # redistanced level set
        if self.RD_modelIndex is not None:
            self.rdModel = modelList[self.RD_modelIndex]
        # level set
        if self.LS_modelIndex is not None:
            self.lsModel = modelList[self.LS_modelIndex]
            self.q_phi = modelList[self.LS_modelIndex].q[('u', 0)]
            self.ebqe_phi = modelList[self.LS_modelIndex].ebqe[('u', 0)]
            if modelList[self.LS_modelIndex].ebq.has_key(('u', 0)):
                self.ebq_phi = modelList[self.LS_modelIndex].ebq[('u', 0)]
        else:
            self.ebqe_phi = numpy.zeros(
                self.model.ebqe[
                    ('u', 0)].shape, 'd')  # cek hack, we don't need this
        # flow model
        # print "flow model
        # index------------",self.V_model,modelList[self.V_model].q.has_key(('velocity',0))
        
        if self.V_model is not None:
            self.velx_tn_dof = modelList[self.V_model].u[0].dof #NOTE: this is not a copy. It is a reference 
            self.vely_tn_dof = modelList[self.V_model].u[1].dof
            self.velz_tn_dof = modelList[self.V_model].u[2].dof

            if modelList[self.V_model].q.has_key(('velocity', 0)):
                self.q_v = modelList[self.V_model].q[('velocity', 0)]
                self.ebqe_v = modelList[self.V_model].ebqe[('velocity', 0)]
            else:
                self.q_v = modelList[self.V_model].q[('f', 0)]
                self.ebqe_v = modelList[self.V_model].ebqe[('f', 0)]
            if modelList[self.V_model].ebq.has_key(('velocity', 0)):
                self.ebq_v = modelList[self.V_model].ebq[('velocity', 0)]
            else:
                if modelList[self.V_model].ebq.has_key(('f', 0)):
                    self.ebq_v = modelList[self.V_model].ebq[('f', 0)]

            # Take divergence from velocity model 
            self.q_div_velocity = modelList[self.V_model].q['div_velocity']
            self.ebqe_div_velocity = modelList[self.V_model].ebqe['div_velocity']
            # Take trial functions and their grads on vel space
            self.nDOF_vel_trial_element = modelList[self.V_model].nDOF_trial_element[0]
            self.vel_l2g = modelList[self.V_model].u[0].femSpace.dofMap.l2g            
            self.velSpace_grad_psi = modelList[self.V_model].u[0].femSpace.grad_psi
        else: #Then it is assumed that the velocity is in the same space as the VOF solution
            self.velx_tn_dof = numpy.zeros(self.model.u[0].dof.shape,'d')
            self.vely_tn_dof = numpy.zeros(self.model.u[0].dof.shape,'d')
            self.velz_tn_dof = numpy.zeros(self.model.u[0].dof.shape,'d')

            # If no velocity model I assume the vel field is div free 
            self.q_div_velocity = numpy.zeros(self.model.q[('u', 0)].shape,'d')
            self.ebqe_div_velocity = numpy.zeros(self.model.ebqe[('u', 0)].shape,'d')

            # Take trial functions and their grads on vel space
            self.nDOF_vel_trial_element = self.model.nDOF_trial_element[0]
            self.vel_l2g = self.model.u[0].femSpace.dofMap.l2g
            self.velSpace_grad_psi = self.model.u[0].femSpace.grad_psi

        if self.eikonalSolverFlag == 2:  # FSW
            self.resDummy = numpy.zeros(self.model.u[0].dof.shape, 'd')
            eikonalSolverType = self.FSWEikonalSolver
            self.eikonalSolver = self.EikonalSolver(eikonalSolverType,
                                                    self.model,
                                                    relativeTolerance=0.0, absoluteTolerance=1.0e-12,
                                                    frontTolerance=1.0e-4,  # default 1.0e-4
                                                    frontInitType='frontIntersection',  # 'frontIntersection',#or 'magnitudeOnly'
                                                    useLocalPWLreconstruction=False)
        elif self.eikonalSolverFlag == 1:  # FMM
            self.resDummy = numpy.zeros(self.model.u[0].dof.shape, 'd')
            eikonalSolverType = self.FMMEikonalSolver
            self.eikonalSolver = self.EikonalSolver(eikonalSolverType,
                                                    self.model,
                                                    frontTolerance=1.0e-4,  # default 1.0e-4
                                                    frontInitType='frontIntersection',  # 'frontIntersection',#or 'magnitudeOnly'
                                                    useLocalPWLreconstruction=False)
        # if self.checkMass:
        #     self.m_pre = Norms.scalarDomainIntegral(self.model.q['dV'],
        #                                              self.model.q[('m',0)],
        #                                              self.model.mesh.nElements_owned)
        #     log("Attach Models VOF: Phase  0 mass after VOF step = %12.5e" % (self.m_pre,),level=2)
        #     self.m_post = Norms.scalarDomainIntegral(self.model.q['dV'],
        #                                              self.model.q[('m',0)],
        #                                              self.model.mesh.nElements_owned)
        #     log("Attach Models VOF: Phase  0 mass after VOF step = %12.5e" % (self.m_post,),level=2)
        #     if self.model.ebqe.has_key(('advectiveFlux',0)):
        #         self.fluxIntegral = Norms.fluxDomainBoundaryIntegral(self.model.ebqe['dS'],
        #                                                              self.model.ebqe[('advectiveFlux',0)],
        #                                                              self.model.mesh)
        #         log("Attach Models VOF: Phase  0 mass conservation after VOF step = %12.5e" % (self.m_post - self.m_pre + self.model.timeIntegration.dt*self.fluxIntegral,),level=2)
        # VRANS
        self.flowCoefficients = modelList[self.V_model].coefficients
        if self.VOS_model is not None:
            self.q_vos = modelList[self.VOS_model].q[('u',0)]
            self.ebqe_vos = modelList[self.VOS_model].ebqe[('u',0)]
        else:
            self.q_porosity = numpy.ones(
                modelList[
                    self.modelIndex].q[
                    ('u', 0)].shape, 'd')
            if self.setParamsFunc is not None:
                self.setParamsFunc(
                    modelList[
                        self.modelIndex].q['x'],
                    self.q_porosity)
            self.q_vos = 1.0 - self.q_porosity
            if hasattr(self.flowCoefficients, 'ebq_porosity'):
                self.ebq_porosity = self.flowCoefficients.ebq_porosity
                self.ebq_vos = 1.0 - self.ebq_porosity
            elif modelList[self.modelIndex].ebq.has_key(('u', 0)):
                self.ebq_porosity = numpy.ones(
                    modelList[
                        self.modelIndex].ebq[
                        ('u', 0)].shape, 'd')
                if self.setParamsFunc is not None:
                    self.setParamsFunc(
                        modelList[
                            self.modelIndex].ebq['x'],
                        self.ebq_porosity)
                self.ebq_vos = 1.0 - self.ebq_porosity
            if hasattr(self.flowCoefficients, 'ebqe_porosity'):
                self.ebqe_porosity = self.flowCoefficients.ebqe_porosity
                self.ebqe_vos = 1.0-self.ebqe_porosity
            else:
                self.ebqe_porosity = numpy.ones(
                    modelList[
                        self.LS_modelIndex].ebqe[
                        ('u', 0)].shape, 'd')
                if self.setParamsFunc is not None:
                    self.setParamsFunc(
                        modelList[
                            self.LS_modelIndex].ebqe['x'],
                        self.ebqe_porosity)
                self.ebqe_vos = 1.0-self.ebqe_porosity

    def initializeElementQuadrature(self, t, cq):
        if self.V_model is None:
            self.q_v = numpy.ones(cq[('f', 0)].shape, 'd')
        # VRANS
        self.q_porosity = numpy.ones(cq[('u', 0)].shape, 'd')

    def initializeElementBoundaryQuadrature(self, t, cebq, cebq_global):
        if self.V_model is None:
            self.ebq_v = numpy.ones(cebq[('f', 0)].shape, 'd')
        # VRANS
        self.ebq_porosity = numpy.ones(cebq[('u', 0)].shape, 'd')

    def initializeGlobalExteriorElementBoundaryQuadrature(self, t, cebqe):
        if self.V_model is None:
            self.ebqe_v = numpy.ones(cebqe[('f', 0)].shape, 'd')
        # VRANS
        self.ebqe_porosity = numpy.ones(cebqe[('u', 0)].shape, 'd')

    def preStep(self, t, firstStep=False):
        self.u_dof_old_old = numpy.copy(self.u_dof_old)
        self.u_dof_old = numpy.copy(self.model.u[0].dof)
        if self.checkMass:
            self.m_pre = Norms.scalarDomainIntegral(
                self.model.q['dV_last'], self.model.q[
                    ('m', 0)], self.model.mesh.nElements_owned)
            log("Phase  0 mass before VOF3P step = %12.5e" %
                (self.m_pre,), level=2)
        #     self.m_last = Norms.scalarDomainIntegral(self.model.q['dV'],
        #                                              self.model.timeIntegration.m_last[0],
        #                                              self.model.mesh.nElements_owned)
        #     log("Phase  0 mass before VOF3P (m_last) step = %12.5e" % (self.m_last,),level=2)
        copyInstructions = {}
        return copyInstructions

    def postStep(self, t, firstStep=False):
        if (self.FCT==1):
            self.model.FCTStep()
        self.model.q['dV_last'][:] = self.model.q['dV']
        if self.checkMass:
            self.m_post = Norms.scalarDomainIntegral(
                self.model.q['dV'], self.model.q[
                    ('m', 0)], self.model.mesh.nElements_owned)
            log("Phase  0 mass after VOF3P step = %12.5e" %
                (self.m_post,), level=2)
            # self.fluxIntegral = Norms.fluxDomainBoundaryIntegral(self.model.ebqe['dS'],
            #                                                     self.model.ebqe[('advectiveFlux',0)],
            #                                                     self.model.mesh)
            #log("Phase  0 mass flux boundary integral after VOF step = %12.5e" % (self.fluxIntegral,),level=2)
            #log("Phase  0 mass conservation after VOF step = %12.5e" % (self.m_post - self.m_last + self.model.timeIntegration.dt*self.fluxIntegral,),level=2)
            # divergence = Norms.fluxDomainBoundaryIntegralFromVector(self.model.ebqe['dS'],
            #                                                        self.ebqe_v,
            #                                                        self.model.ebqe['n'],
            #                                                        self.model.mesh)
            #log("Divergence = %12.5e" % (divergence,),level=2)
        copyInstructions = {}
        return copyInstructions

    def updateToMovingDomain(self, t, c):
        # in a moving domain simulation the velocity coming in is already for
        # the moving domain
        pass

    def evaluate(self, t, c):
        # mwf debug
        # print "VOF3Pcoeficients eval t=%s " % t
        if c[('f', 0)].shape == self.q_v.shape:
            v = self.q_v
            phi = self.q_phi
            porosity = self.q_porosity
        elif c[('f', 0)].shape == self.ebqe_v.shape:
            v = self.ebqe_v
            phi = self.ebqe_phi
            porosity = self.ebq_porosity
        elif ((self.ebq_v is not None and self.ebq_phi is not None) and c[('f', 0)].shape == self.ebq_v.shape):
            v = self.ebq_v
            phi = self.ebq_phi
            porosity = self.ebq_porosity
        else:
            v = None
            phi = None
            porosity = None
        if v is not None:
            # self.VOF3PCoefficientsEvaluate(self.eps,
            #                              v,
            #                              phi,
            #                              c[('u',0)],
            #                              c[('m',0)],
            #                              c[('dm',0,0)],
            #                              c[('f',0)],
            #                              c[('df',0,0)])
            self.VolumeAveragedVOFCoefficientsEvaluate(self.eps,
                                                       v,
                                                       phi,
                                                       porosity,
                                                       c[('u', 0)],
                                                       c[('m', 0)],
                                                       c[('dm', 0, 0)],
                                                       c[('f', 0)],
                                                       c[('df', 0, 0)])
        # if self.checkMass:
        #     log("Phase  0 mass in eavl = %12.5e" % (Norms.scalarDomainIntegral(self.model.q['dV'],
        #                                                                        self.model.q[('m',0)],
        # self.model.mesh.nElements_owned),),level=2)


class LevelModel(proteus.Transport.OneLevelTransport):
    nCalls = 0

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
                 sd=True,
                 movingDomain=False):
        #
        # set the objects describing the method and boundary conditions
        #
        self.movingDomain = movingDomain
        self.tLast_mesh = None
        #
        self.name = name
        self.sd = sd
        self.Hess = False
        self.lowmem = True
        self.timeTerm = True  # allow turning off  the  time derivative
        # self.lowmem=False
        self.testIsTrial = True
        self.phiTrialIsTrial = True
        self.u = uDict
        self.ua = {}  # analytical solutions
        self.phi = phiDict
        self.dphi = {}
        self.matType = matType
        # mwf try to reuse test and trial information across components if
        # spaces are the same
        self.reuse_test_trial_quadrature = reuse_trial_and_test_quadrature  # True#False
        if self.reuse_test_trial_quadrature:
            for ci in range(1, coefficients.nc):
                assert self.u[ci].femSpace.__class__.__name__ == self.u[
                    0].femSpace.__class__.__name__, "to reuse_test_trial_quad all femSpaces must be the same!"
        # Simplicial Mesh
        # assume the same mesh for  all components for now
        self.mesh = self.u[0].femSpace.mesh
        self.testSpace = testSpaceDict
        self.dirichletConditions = dofBoundaryConditionsDict
        # explicit Dirichlet  conditions for now, no Dirichlet BC constraints
        self.dirichletNodeSetList = None
        self.coefficients = coefficients
        self.coefficients.initializeMesh(self.mesh)
        self.nc = self.coefficients.nc
        self.stabilization = stabilization
        self.shockCapturing = shockCapturing
        # no velocity post-processing for now
        self.conservativeFlux = conservativeFluxDict
        self.fluxBoundaryConditions = fluxBoundaryConditionsDict
        self.advectiveFluxBoundaryConditionsSetterDict = advectiveFluxBoundaryConditionsSetterDict
        self.diffusiveFluxBoundaryConditionsSetterDictDict = diffusiveFluxBoundaryConditionsSetterDictDict
        # determine whether  the stabilization term is nonlinear
        self.stabilizationIsNonlinear = False
        # cek come back
        if self.stabilization is not None:
            for ci in range(self.nc):
                if coefficients.mass.has_key(ci):
                    for flag in coefficients.mass[ci].values():
                        if flag == 'nonlinear':
                            self.stabilizationIsNonlinear = True
                if coefficients.advection.has_key(ci):
                    for flag in coefficients.advection[ci].values():
                        if flag == 'nonlinear':
                            self.stabilizationIsNonlinear = True
                if coefficients.diffusion.has_key(ci):
                    for diffusionDict in coefficients.diffusion[ci].values():
                        for flag in diffusionDict.values():
                            if flag != 'constant':
                                self.stabilizationIsNonlinear = True
                if coefficients.potential.has_key(ci):
                    for flag in coefficients.potential[ci].values():
                        if flag == 'nonlinear':
                            self.stabilizationIsNonlinear = True
                if coefficients.reaction.has_key(ci):
                    for flag in coefficients.reaction[ci].values():
                        if flag == 'nonlinear':
                            self.stabilizationIsNonlinear = True
                if coefficients.hamiltonian.has_key(ci):
                    for flag in coefficients.hamiltonian[ci].values():
                        if flag == 'nonlinear':
                            self.stabilizationIsNonlinear = True
        # determine if we need element boundary storage
        self.elementBoundaryIntegrals = {}
        for ci in range(self.nc):
            self.elementBoundaryIntegrals[ci] = (
                (self.conservativeFlux is not None) or (
                    numericalFluxType is not None) or (
                    self.fluxBoundaryConditions[ci] == 'outFlow') or (
                    self.fluxBoundaryConditions[ci] == 'mixedFlow') or (
                    self.fluxBoundaryConditions[ci] == 'setFlow'))
        #
        # calculate some dimensions
        #
        # assume same space dim for all variables
        self.nSpace_global = self.u[0].femSpace.nSpace_global
        self.nDOF_trial_element = [
            u_j.femSpace.max_nDOF_element for u_j in self.u.values()]
        self.nDOF_phi_trial_element = [
            phi_k.femSpace.max_nDOF_element for phi_k in self.phi.values()]
        self.n_phi_ip_element = [
            phi_k.femSpace.referenceFiniteElement.interpolationConditions.nQuadraturePoints for phi_k in self.phi.values()]
        self.nDOF_test_element = [
            femSpace.max_nDOF_element for femSpace in self.testSpace.values()]
        self.nFreeDOF_global = [
            dc.nFreeDOF_global for dc in self.dirichletConditions.values()]
        self.nVDOF_element = sum(self.nDOF_trial_element)
        self.nFreeVDOF_global = sum(self.nFreeDOF_global)
        #
        NonlinearEquation.__init__(self, self.nFreeVDOF_global)
        #
        # build the quadrature point dictionaries from the input (this
        # is just for convenience so that the input doesn't have to be
        # complete)
        #

        elementQuadratureDict = {}
        elemQuadIsDict = isinstance(elementQuadrature, dict)
        if elemQuadIsDict:  # set terms manually
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
                        elementQuadratureDict[
                            ('stab',) + I[1:]] = elementQuadrature[I]
                    else:
                        elementQuadratureDict[
                            ('stab',) + I[1:]] = elementQuadrature['default']
                else:
                    elementQuadratureDict[
                        ('stab',) + I[1:]] = elementQuadrature
        if self.shockCapturing is not None:
            for ci in self.shockCapturing.components:
                if elemQuadIsDict:
                    if elementQuadrature.has_key(('numDiff', ci, ci)):
                        elementQuadratureDict[('numDiff', ci, ci)] = elementQuadrature[
                            ('numDiff', ci, ci)]
                    else:
                        elementQuadratureDict[('numDiff', ci, ci)] = elementQuadrature[
                            'default']
                else:
                    elementQuadratureDict[
                        ('numDiff', ci, ci)] = elementQuadrature
        if massLumping:
            for ci in self.coefficients.mass.keys():
                elementQuadratureDict[('m', ci)] = Quadrature.SimplexLobattoQuadrature(
                    self.nSpace_global, 1)
            for I in self.coefficients.elementIntegralKeys:
                elementQuadratureDict[
                    ('stab',) + I[1:]] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global, 1)
        if reactionLumping:
            for ci in self.coefficients.mass.keys():
                elementQuadratureDict[('r', ci)] = Quadrature.SimplexLobattoQuadrature(
                    self.nSpace_global, 1)
            for I in self.coefficients.elementIntegralKeys:
                elementQuadratureDict[
                    ('stab',) + I[1:]] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global, 1)
        elementBoundaryQuadratureDict = {}
        if isinstance(elementBoundaryQuadrature, dict):  # set terms manually
            for I in self.coefficients.elementBoundaryIntegralKeys:
                if elementBoundaryQuadrature.has_key(I):
                    elementBoundaryQuadratureDict[
                        I] = elementBoundaryQuadrature[I]
                else:
                    elementBoundaryQuadratureDict[
                        I] = elementBoundaryQuadrature['default']
        else:
            for I in self.coefficients.elementBoundaryIntegralKeys:
                elementBoundaryQuadratureDict[I] = elementBoundaryQuadrature
        #
        # find the union of all element quadrature points and
        # build a quadrature rule for each integral that has a
        # weight at each point in the union
        # mwf include tag telling me which indices are which quadrature rule?
        (self.elementQuadraturePoints, self.elementQuadratureWeights,
         self.elementQuadratureRuleIndeces) = Quadrature.buildUnion(elementQuadratureDict)
        self.nQuadraturePoints_element = self.elementQuadraturePoints.shape[0]
        self.nQuadraturePoints_global = self.nQuadraturePoints_element * \
            self.mesh.nElements_global
        #
        # Repeat the same thing for the element boundary quadrature
        #
        (self.elementBoundaryQuadraturePoints, self.elementBoundaryQuadratureWeights,
         self.elementBoundaryQuadratureRuleIndeces) = Quadrature.buildUnion(elementBoundaryQuadratureDict)
        self.nElementBoundaryQuadraturePoints_elementBoundary = self.elementBoundaryQuadraturePoints.shape[
            0]
        self.nElementBoundaryQuadraturePoints_global = (
            self.mesh.nElements_global *
            self.mesh.nElementBoundaries_element *
            self.nElementBoundaryQuadraturePoints_elementBoundary)
#        if isinstance(self.u[0].femSpace,C0_AffineLinearOnSimplexWithNodalBasis):
#            print self.nQuadraturePoints_element
#            if self.nSpace_global == 3:
#                assert(self.nQuadraturePoints_element == 5)
#            elif self.nSpace_global == 2:
#                assert(self.nQuadraturePoints_element == 6)
#            elif self.nSpace_global == 1:
#                assert(self.nQuadraturePoints_element == 3)
#
#            print self.nElementBoundaryQuadraturePoints_elementBoundary
#            if self.nSpace_global == 3:
#                assert(self.nElementBoundaryQuadraturePoints_elementBoundary == 4)
#            elif self.nSpace_global == 2:
#                assert(self.nElementBoundaryQuadraturePoints_elementBoundary == 4)
#            elif self.nSpace_global == 1:
#                assert(self.nElementBoundaryQuadraturePoints_elementBoundary == 1)

        #
        # storage dictionaries
        self.scalars_element = set()
        #
        # simplified allocations for test==trial and also check if space is mixed or not
        #
        self.q = {}
        self.ebq = {}
        self.ebq_global = {}
        self.ebqe = {}
        self.phi_ip = {}
        # mesh
        self.q['x'] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,3),'d')
        self.ebqe['x'] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary,
             3),
            'd')
        self.q[('u', 0)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('dV_u', 0)] = (1.0 / self.mesh.nElements_global) * \
            numpy.ones((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[
            ('grad(u)',
             0)] = numpy.zeros(
            (self.mesh.nElements_global,
             self.nQuadraturePoints_element,
             self.nSpace_global),
            'd')
        self.q[('m', 0)] = self.q[('u', 0)]
        self.q[('m_last', 0)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('mt', 0)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q['dV'] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q['dV_last'] = -1000 * \
            numpy.ones((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('m_tmp', 0)] = self.q[('u', 0)]
        self.q[('cfl', 0)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('numDiff', 0, 0)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.ebqe[
            ('u',
             0)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary),
            'd')
        self.ebqe[
            ('grad(u)',
             0)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary,
             self.nSpace_global),
            'd')
        self.ebqe[
            ('advectiveFlux_bc_flag',
             0)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary),
            'i')
        self.ebqe[
            ('advectiveFlux_bc',
             0)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary),
            'd')
        self.ebqe[
            ('advectiveFlux',
             0)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary),
            'd')
        self.points_elementBoundaryQuadrature = set()
        self.scalars_elementBoundaryQuadrature = set(
            [('u', ci) for ci in range(self.nc)])
        self.vectors_elementBoundaryQuadrature = set()
        self.tensors_elementBoundaryQuadrature = set()
        self.inflowBoundaryBC = {}
        self.inflowBoundaryBC_values = {}
        self.inflowFlux = {}
        for cj in range(self.nc):
            self.inflowBoundaryBC[cj] = numpy.zeros(
                (self.mesh.nExteriorElementBoundaries_global,), 'i')
            self.inflowBoundaryBC_values[cj] = numpy.zeros(
                (self.mesh.nExteriorElementBoundaries_global, self.nDOF_trial_element[cj]), 'd')
            self.inflowFlux[cj] = numpy.zeros(
                (self.mesh.nExteriorElementBoundaries_global,
                 self.nElementBoundaryQuadraturePoints_elementBoundary),
                'd')
        self.internalNodes = set(range(self.mesh.nNodes_global))
        # identify the internal nodes this is ought to be in mesh
        # \todo move this to mesh
        for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
            ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
            eN_global = self.mesh.elementBoundaryElementsArray[ebN, 0]
            ebN_element = self.mesh.elementBoundaryLocalElementBoundariesArray[
                ebN, 0]
            for i in range(self.mesh.nNodes_element):
                if i != ebN_element:
                    I = self.mesh.elementNodesArray[eN_global, i]
                    self.internalNodes -= set([I])
        self.nNodes_internal = len(self.internalNodes)
        self.internalNodesArray = numpy.zeros((self.nNodes_internal,), 'i')
        for nI, n in enumerate(self.internalNodes):
            self.internalNodesArray[nI] = n
        #
        del self.internalNodes
        self.internalNodes = None
        log("Updating local to global mappings", 2)
        self.updateLocal2Global()
        log("Building time integration object", 2)
        log(memory("inflowBC, internalNodes,updateLocal2Global",
                   "OneLevelTransport"), level=4)
        # mwf for interpolating subgrid error for gradients etc
        if self.stabilization and self.stabilization.usesGradientStabilization:
            self.timeIntegration = TimeIntegrationClass(
                self, integrateInterpolationPoints=True)
        else:
            self.timeIntegration = TimeIntegrationClass(self)

        if options is not None:
            self.timeIntegration.setFromOptions(options)
        log(memory("TimeIntegration", "OneLevelTransport"), level=4)
        log("Calculating numerical quadrature formulas", 2)
        self.calculateQuadrature()
        self.setupFieldStrides()

        #cek adding empty data member for low order numerical viscosity structures here for now
        self.ML=None #lumped mass matrix
        self.MC_global=None #consistent mass matrix
        self.cterm_global=None
        self.cterm_transpose_global=None
        # dL_global and dC_global are not the full matrices but just the CSR arrays containing the non zero entries
        self.flux_plus_dLij_times_soln=None
        self.dL_minus_dC=None
        self.min_u_bc=None
        self.max_u_bc=None
        self.inflow_DOFs=None
        # Aux quantity at DOFs to be filled by optimized code (MQL)
        self.quantDOFs=None

        comm = Comm.get()
        self.comm = comm
        if comm.size() > 1:
            assert numericalFluxType is not None and numericalFluxType.useWeakDirichletConditions, "You must use a numerical flux to apply weak boundary conditions for parallel runs"

        log(memory("stride+offset", "OneLevelTransport"), level=4)
        if numericalFluxType is not None:
            if options is None or options.periodicDirichletConditions is None:
                self.numericalFlux = numericalFluxType(
                    self,
                    dofBoundaryConditionsSetterDict,
                    advectiveFluxBoundaryConditionsSetterDict,
                    diffusiveFluxBoundaryConditionsSetterDictDict)
            else:
                self.numericalFlux = numericalFluxType(
                    self,
                    dofBoundaryConditionsSetterDict,
                    advectiveFluxBoundaryConditionsSetterDict,
                    diffusiveFluxBoundaryConditionsSetterDictDict,
                    options.periodicDirichletConditions)
        else:
            self.numericalFlux = None
        # set penalty terms
        # cek todo move into numerical flux initialization
        if self.ebq_global.has_key('penalty'):
            for ebN in range(self.mesh.nElementBoundaries_global):
                for k in range(
                        self.nElementBoundaryQuadraturePoints_elementBoundary):
                    self.ebq_global['penalty'][ebN, k] = self.numericalFlux.penalty_constant / (
                        self.mesh.elementBoundaryDiametersArray[ebN]**self.numericalFlux.penalty_power)
        # penalty term
        # cek move  to Numerical flux initialization
        if self.ebqe.has_key('penalty'):
            for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
                ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
                for k in range(
                        self.nElementBoundaryQuadraturePoints_elementBoundary):
                    self.ebqe['penalty'][ebNE, k] = self.numericalFlux.penalty_constant / \
                        self.mesh.elementBoundaryDiametersArray[ebN]**self.numericalFlux.penalty_power
        log(memory("numericalFlux", "OneLevelTransport"), level=4)
        self.elementEffectiveDiametersArray = self.mesh.elementInnerDiametersArray
        # use post processing tools to get conservative fluxes, None by default
        from proteus import PostProcessingTools
        self.velocityPostProcessor = PostProcessingTools.VelocityPostProcessingChooser(
            self)
        log(memory("velocity postprocessor", "OneLevelTransport"), level=4)
        # helper for writing out data storage
        from proteus import Archiver
        self.elementQuadratureDictionaryWriter = Archiver.XdmfWriter()
        self.elementBoundaryQuadratureDictionaryWriter = Archiver.XdmfWriter()
        self.exteriorElementBoundaryQuadratureDictionaryWriter = Archiver.XdmfWriter()
        # TODO get rid of this
        for ci, fbcObject in self.fluxBoundaryConditionsObjectsDict.iteritems():
            self.ebqe[('advectiveFlux_bc_flag', ci)] = numpy.zeros(
                self.ebqe[('advectiveFlux_bc', ci)].shape, 'i')
            for t, g in fbcObject.advectiveFluxBoundaryConditionsDict.iteritems():
                if self.coefficients.advection.has_key(ci):
                    self.ebqe[
                        ('advectiveFlux_bc', ci)][
                        t[0], t[1]] = g(
                        self.ebqe[
                            ('x')][
                            t[0], t[1]], self.timeIntegration.t)
                    self.ebqe[('advectiveFlux_bc_flag', ci)][t[0], t[1]] = 1

        if hasattr(self.numericalFlux, 'setDirichletValues'):
            self.numericalFlux.setDirichletValues(self.ebqe)
        if not hasattr(self.numericalFlux, 'isDOFBoundary'):
            self.numericalFlux.isDOFBoundary = {
                0: numpy.zeros(self.ebqe[('u', 0)].shape, 'i')}
        if not hasattr(self.numericalFlux, 'ebqe'):
            self.numericalFlux.ebqe = {
                ('u', 0): numpy.zeros(self.ebqe[('u', 0)].shape, 'd')}
        # TODO how to handle redistancing calls for
        # calculateCoefficients,calculateElementResidual etc
        self.globalResidualDummy = None
        compKernelFlag = 0
        self.vof = VOF3P(
            self.nSpace_global,
            self.nQuadraturePoints_element,
            self.u[0].femSpace.elementMaps.localFunctionSpace.dim,
            self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.dim,
            self.testSpace[0].referenceFiniteElement.localFunctionSpace.dim,
            self.nElementBoundaryQuadraturePoints_elementBoundary,
            compKernelFlag)

        self.forceStrongConditions = self.coefficients.forceStrongConditions
        if self.forceStrongConditions:
            self.dirichletConditionsForceDOF = DOFBoundaryConditions(
                self.u[0].femSpace,
                dofBoundaryConditionsSetterDict[0],
                weakDirichletConditions=False)

        if self.movingDomain:
            self.MOVING_DOMAIN = 1.0
        else:
            self.MOVING_DOMAIN = 0.0
        if self.mesh.nodeVelocityArray is None:
            self.mesh.nodeVelocityArray = numpy.zeros(
                self.mesh.nodeArray.shape, 'd')
    # mwf these are getting called by redistancing classes,

    def FCTStep(self):
        rowptr, colind, MassMatrix = self.MC_global.getCSRrepresentation()
        self.vof.FCTStep(self.timeIntegration.dt, 
                         self.nnz, #number of non zero entries 
                         len(rowptr)-1, #number of DOFs
                         self.ML, #Lumped mass matrix
                         self.coefficients.u_dof_old, #soln
                         self.u[0].dof, #solH
                         self.flux_plus_dLij_times_soln, 
                         rowptr, #Row indices for Sparsity Pattern (convenient for DOF loops)
                         colind, #Column indices for Sparsity Pattern (convenient for DOF loops)
                         MassMatrix, 
                         self.dL_minus_dC,
                         self.min_u_bc,
                         self.max_u_bc)
    def calculateCoefficients(self):
        pass

    def calculateElementResidual(self):
        if self.globalResidualDummy is not None:
            self.getResidual(self.u[0].dof, self.globalResidualDummy)

    def getResidual(self, u, r):
        import pdb
        import copy
        """
        Calculate the element residuals and add in to the global residual
        """

        #COMPUTE C MATRIX 
        if self.cterm_global is None:
            #since we only need cterm_global to persist, we can drop the other self.'s
            self.cterm={}
            self.cterm_a={}
            self.cterm_global={}
            self.cterm_transpose={}
            self.cterm_a_transpose={}
            self.cterm_global_transpose={}
            rowptr, colind, nzval = self.jacobian.getCSRrepresentation()
            nnz = nzval.shape[-1] #number of non-zero entries in sparse matrix
            di = self.q[('grad(u)',0)].copy() #direction of derivative
            # JACOBIANS (FOR ELEMENT TRANSFORMATION)
            self.q[('J')] = numpy.zeros((self.mesh.nElements_global,
                                         self.nQuadraturePoints_element,
                                         self.nSpace_global,
                                         self.nSpace_global),
                                        'd')
            self.q[('inverse(J)')] = numpy.zeros((self.mesh.nElements_global,
                                                  self.nQuadraturePoints_element,
                                                  self.nSpace_global,
                                                  self.nSpace_global),
                                                 'd')
            self.q[('det(J)')] = numpy.zeros((self.mesh.nElements_global,
                                              self.nQuadraturePoints_element),
                                             'd')
            self.u[0].femSpace.elementMaps.getJacobianValues(self.elementQuadraturePoints,
                                                             self.q['J'],
                                                             self.q['inverse(J)'],
                                                             self.q['det(J)'])
            self.q['abs(det(J))'] = numpy.abs(self.q['det(J)'])
            # SHAPE FUNCTIONS
            self.q[('w',0)] = numpy.zeros((self.mesh.nElements_global,
                                           self.nQuadraturePoints_element,
                                           self.nDOF_test_element[0]),
                                          'd')
            self.q[('w*dV_m',0)] = self.q[('w',0)].copy()
            self.u[0].femSpace.getBasisValues(self.elementQuadraturePoints, self.q[('w',0)])
            cfemIntegrals.calculateWeightedShape(self.elementQuadratureWeights[('u',0)],
                                                 self.q['abs(det(J))'],
                                                 self.q[('w',0)],
                                                 self.q[('w*dV_m',0)])
            #### GRADIENT OF TEST FUNCTIONS 
            self.q[('grad(w)',0)] = numpy.zeros((self.mesh.nElements_global,
                                                 self.nQuadraturePoints_element,
                                                 self.nDOF_test_element[0],
                                                 self.nSpace_global),
                                                'd')
            self.u[0].femSpace.getBasisGradientValues(self.elementQuadraturePoints,
                                                      self.q['inverse(J)'],
                                                      self.q[('grad(w)',0)])
            self.q[('grad(w)*dV_f',0)] = numpy.zeros((self.mesh.nElements_global,
                                                      self.nQuadraturePoints_element,
                                                      self.nDOF_test_element[0],
                                                      self.nSpace_global),
                                                     'd')
            cfemIntegrals.calculateWeightedShapeGradients(self.elementQuadratureWeights[('u',0)],
                                                          self.q['abs(det(J))'],
                                                          self.q[('grad(w)',0)],
                                                          self.q[('grad(w)*dV_f',0)])
            #
            #lumped mass matrix
            #
            #assume a linear mass term
            dm = numpy.ones(self.q[('u',0)].shape,'d')
            elementMassMatrix = numpy.zeros((self.mesh.nElements_global,
                                             self.nDOF_test_element[0],
                                             self.nDOF_trial_element[0]),'d')
            cfemIntegrals.updateMassJacobian_weak_lowmem(dm,
                                                         self.q[('w',0)],
                                                         self.q[('w*dV_m',0)],
                                                         elementMassMatrix)
            self.MC_a = nzval.copy()
            self.MC_global = LinearAlgebraTools.SparseMat(self.nFreeDOF_global[0],
                                                          self.nFreeDOF_global[0],
                                                          nnz,
                                                          self.MC_a,
                                                          colind,
                                                          rowptr)
            cfemIntegrals.zeroJacobian_CSR(self.nnz, self.MC_global)
            cfemIntegrals.updateGlobalJacobianFromElementJacobian_CSR(self.l2g[0]['nFreeDOF'],
                                                                      self.l2g[0]['freeLocal'],
                                                                      self.l2g[0]['nFreeDOF'],
                                                                      self.l2g[0]['freeLocal'],
                                                                      self.csrRowIndeces[(0,0)],
                                                                      self.csrColumnOffsets[(0,0)],
                                                                      elementMassMatrix,
                                                                      self.MC_global)
            self.ML = numpy.zeros((self.nFreeDOF_global[0],),'d')
            for i in range(self.nFreeDOF_global[0]):
                self.ML[i] = self.MC_a[rowptr[i]:rowptr[i+1]].sum()
            numpy.testing.assert_almost_equal(self.ML.sum(), self.mesh.volume, err_msg="Trace of lumped mass matrix should be the domain volume",verbose=True)
            for d in range(self.nSpace_global): #spatial dimensions
                #C matrices
                self.cterm[d] = numpy.zeros((self.mesh.nElements_global,
                                             self.nDOF_test_element[0],
                                             self.nDOF_trial_element[0]),'d')
                self.cterm_a[d] = nzval.copy()
                self.cterm_global[d] = LinearAlgebraTools.SparseMat(self.nFreeDOF_global[0],
                                                                    self.nFreeDOF_global[0],
                                                                    nnz,
                                                                    self.cterm_a[d],
                                                                    colind,
                                                                    rowptr)
                cfemIntegrals.zeroJacobian_CSR(self.nnz, self.cterm_global[d])
                di[:] = 0.0
                di[...,d] = 1.0
                cfemIntegrals.updateHamiltonianJacobian_weak_lowmem(di,
                                                                    self.q[('grad(w)*dV_f',0)],
                                                                    self.q[('w',0)],
                                                                    self.cterm[d]) # int[(di*grad(wj))*wi*dV]
                cfemIntegrals.updateGlobalJacobianFromElementJacobian_CSR(self.l2g[0]['nFreeDOF'],
                                                                          self.l2g[0]['freeLocal'],
                                                                          self.l2g[0]['nFreeDOF'],
                                                                          self.l2g[0]['freeLocal'],
                                                                          self.csrRowIndeces[(0,0)],
                                                                          self.csrColumnOffsets[(0,0)],
                                                                          self.cterm[d],
                                                                          self.cterm_global[d])
                #C Transpose matrices
                self.cterm_transpose[d] = numpy.zeros((self.mesh.nElements_global,
                                                       self.nDOF_test_element[0],
                                                       self.nDOF_trial_element[0]),'d')
                self.cterm_a_transpose[d] = nzval.copy()
                self.cterm_global_transpose[d] = LinearAlgebraTools.SparseMat(self.nFreeDOF_global[0],
                                                                              self.nFreeDOF_global[0],
                                                                              nnz,
                                                                              self.cterm_a_transpose[d],
                                                                              colind,
                                                                              rowptr)
                cfemIntegrals.zeroJacobian_CSR(self.nnz, self.cterm_global_transpose[d])
                di[:] = 0.0
                di[...,d] = -1.0
                cfemIntegrals.updateAdvectionJacobian_weak_lowmem(di,
                                                                  self.q[('w',0)], 
                                                                  self.q[('grad(w)*dV_f',0)], 
                                                                  self.cterm_transpose[d]) # -int[(-di*grad(wi))*wj*dV]
                cfemIntegrals.updateGlobalJacobianFromElementJacobian_CSR(self.l2g[0]['nFreeDOF'],
                                                                          self.l2g[0]['freeLocal'],
                                                                          self.l2g[0]['nFreeDOF'],
                                                                          self.l2g[0]['freeLocal'],
                                                                          self.csrRowIndeces[(0,0)],
                                                                          self.csrColumnOffsets[(0,0)],
                                                                          self.cterm_transpose[d],
                                                                          self.cterm_global_transpose[d])

        rowptr, colind, Cx = self.cterm_global[0].getCSRrepresentation()
        rowptr, colind, Cy = self.cterm_global[1].getCSRrepresentation()        
        Cz = numpy.zeros(Cx.shape,'d')
        if (self.nSpace_global==3):
            rowptr, colind, Cz = self.cterm_global[2].getCSRrepresentation()        
        rowptr, colind, CTx = self.cterm_global_transpose[0].getCSRrepresentation()
        rowptr, colind, CTy = self.cterm_global_transpose[1].getCSRrepresentation()
        CTz = numpy.zeros(Cx.shape,'d')
        if (self.nSpace_global==3):
            rowptr, colind, CTz = self.cterm_global_transpose[2].getCSRrepresentation()        
        # This is dummy. I just care about the csr structure of the sparse matrix
        self.dL_minus_dC = numpy.zeros(Cx.shape,'d')
        self.min_u_bc = numpy.zeros(self.u[0].dof.shape,'d')
        self.max_u_bc = numpy.zeros(self.u[0].dof.shape,'d')
        self.min_u_bc.fill(1E10);
        self.max_u_bc.fill(-1E10);
        self.flux_plus_dLij_times_soln = numpy.zeros(self.u[0].dof.shape,'d')
        self.inflow_DOFs = numpy.zeros(self.u[0].dof.shape,'d')
        self.quantDOFs = numpy.zeros(self.u[0].dof.shape,'d')
        #
        #cek end computationa of cterm_global
        #
        #cek showing mquezada an example of using cterm_global sparse matrix
        #calculation y = c*x where x==1
        #direction=0
        #rowptr, colind, c = self.cterm_global[direction].getCSRrepresentation()
        #y = numpy.zeros((self.nFreeDOF_global[0],),'d')
        #x = numpy.ones((self.nFreeDOF_global[0],),'d')
        #ij=0
        #for i in range(self.nFreeDOF_global[0]):
        #    for offset in range(rowptr[i],rowptr[i+1]):
        #        j = colind[offset]
        #        y[i] += c[ij]*x[j]
        #        ij+=1

        if self.coefficients.set_vos:
            self.coefficients.set_vos(self.q['x'],self.coefficients.q_vos)
        r.fill(0.0)
        # Load the unknowns into the finite element dof
        self.timeIntegration.calculateCoefs()
        # print "***************max/min(m_last)*********************",max(self.timeIntegration.m_last[0].flat[:]),min(self.timeIntegration.m_last[0].flat[:])
        # print
        # "***************max/min(m_last)*********************",max(-self.timeIntegration.dt*self.timeIntegration.beta_bdf[0].flat[:]),min(-self.timeIntegration.dt*self.timeIntegration.beta_bdf[0].flat[:]),
        self.timeIntegration.calculateU(u)
        self.setUnknowns(self.timeIntegration.u)
        # cek can put in logic to skip of BC's don't depend on t or u
        # Dirichlet boundary conditions
        # if hasattr(self.numericalFlux,'setDirichletValues'):
        self.numericalFlux.setDirichletValues(self.ebqe)
        # flux boundary conditions
        for t, g in self.fluxBoundaryConditionsObjectsDict[
                0].advectiveFluxBoundaryConditionsDict.iteritems():
            self.ebqe[
                ('advectiveFlux_bc', 0)][
                t[0], t[1]] = g(
                self.ebqe[
                    ('x')][
                    t[0], t[1]], self.timeIntegration.t)
            self.ebqe[('advectiveFlux_bc_flag', 0)][t[0], t[1]] = 1

        # MQL: Find out the DOFs that are inflow
        if self.forceStrongConditions:
            self.vof.getInflowDOFs(self.mesh.nodeArray,
                                   self.mesh.elementNodesArray,
                                   self.u[0].femSpace.elementMaps.psi_trace,
                                   self.u[0].femSpace.elementMaps.grad_psi_trace,
                                   self.u[0].femSpace.elementMaps.boundaryNormals,
                                   self.u[0].femSpace.elementMaps.boundaryJacobians,
                                   self.u[0].femSpace.dofMap.l2g,
                                   self.offset[0], 
                                   self.stride[0],
                                   self.mesh.nExteriorElementBoundaries_global,
                                   self.mesh.exteriorElementBoundariesArray,
                                   self.mesh.elementBoundaryElementsArray,
                                   self.mesh.elementBoundaryLocalElementBoundariesArray,
                                   self.coefficients.ebqe_v,
                                   self.inflow_DOFs)
            
        if self.forceStrongConditions:
            for dofN, g in self.dirichletConditionsForceDOF.DOFBoundaryConditionsDict.iteritems():
                #if self.coefficients.vely_tn_dof[dofN] < 0: #HACKED FOR THIS PROBLEM (MQL) TMP
                if self.inflow_DOFs[dofN] == 1.0: 
                    self.u[0].dof[dofN] = g(
                        self.dirichletConditionsForceDOF.DOFBoundaryPointDict[dofN],
                        self.timeIntegration.t)
        assert (self.coefficients.q_porosity == 1).all()

        self.vof.calculateResidual(  # element
            self.u[0].femSpace.elementMaps.psi,
            self.u[0].femSpace.elementMaps.grad_psi,
            self.mesh.nodeArray,
            self.mesh.nodeVelocityArray,
            self.MOVING_DOMAIN,
            self.mesh.elementNodesArray,
            self.elementQuadratureWeights[('u', 0)],
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            # trial grads on vel Space
            self.coefficients.nDOF_vel_trial_element,
            self.coefficients.vel_l2g,
            self.coefficients.velSpace_grad_psi,
            # element boundary
            self.u[0].femSpace.elementMaps.psi_trace,
            self.u[0].femSpace.elementMaps.grad_psi_trace,
            self.elementBoundaryQuadratureWeights[('u', 0)],
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.elementMaps.boundaryNormals,
            self.u[0].femSpace.elementMaps.boundaryJacobians,
            # physics
            self.mesh.nElements_global,
            self.coefficients.useMetrics,
            self.timeIntegration.alpha_bdf,
            self.shockCapturing.lag,
            self.shockCapturing.shockCapturingFactor,
            self.coefficients.sc_uref,
            self.coefficients.sc_beta,
            # VRANS start
            self.coefficients.q_vos,
            # VRANS end
            self.u[0].femSpace.dofMap.l2g,
            self.mesh.elementDiametersArray,
            self.u[0].dof,
            self.coefficients.u_dof_old,
            self.coefficients.u_dof_old_old,
            self.coefficients.velx_tn_dof, 
            self.coefficients.vely_tn_dof, 
            self.coefficients.velz_tn_dof, 
            self.coefficients.q_v,
            self.coefficients.q_div_velocity,
            self.timeIntegration.m_tmp[0],
            self.q[('u', 0)],
            self.timeIntegration.beta_bdf[0],
            self.q['dV'],
            self.q['dV_last'],
            self.q[('cfl', 0)],
            self.shockCapturing.numDiff[0],
            self.shockCapturing.numDiff_last[0],
            self.offset[0], self.stride[0],
            r,
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray,
            self.coefficients.ebqe_v,
            self.coefficients.ebqe_div_velocity,
            # VRANS start
            self.coefficients.ebqe_vos,
            # VRANS end
            self.numericalFlux.isDOFBoundary[0],
            self.numericalFlux.ebqe[('u', 0)],
            self.ebqe[('advectiveFlux_bc_flag', 0)],
            self.ebqe[('advectiveFlux_bc', 0)],
            self.coefficients.ebqe_phi, self.coefficients.epsFact,
            self.ebqe[('u', 0)],
            self.ebqe[('advectiveFlux', 0)],             
            # PARAMETERS FOR EDGE_VISCOSITY
            self.coefficients.EDGE_VISCOSITY,
            self.coefficients.ENTROPY_VISCOSITY, 
            # PARAMETERS FOR ELEMENT BASED ENTROPY VISCOSITY
            self.coefficients.cE, 
            self.coefficients.cMax, 
            self.coefficients.cK,
            # PARAMETERS FOR LOG BASED ENTROPY FUNCTION 
            self.coefficients.uL,
            self.coefficients.uR,
            len(rowptr)-1, #num of DOFs
            self.nnz,
            rowptr, #Row indices for Sparsity Pattern (convenient for DOF loops)
            colind, #Column indices for Sparsity Pattern (convenient for DOF loops)
            self.csrRowIndeces[(0,0)], #row indices (convenient for element loops) 
            self.csrColumnOffsets[(0,0)], #column indices (convenient for element loops)
            self.csrColumnOffsets_eb[(0, 0)], #indices for boundary terms
            Cx, 
            Cy, 
            Cz, 
            CTx,
            CTy,
            CTz,
            # PARAMETERS FOR 1st or 2nd ORDER MPP METHOD 
            self.coefficients.POWER_SMOOTHNESS_INDICATOR,
            self.coefficients.LUMPED_MASS_MATRIX,
            # FLUX CORRECTED TRANSPORT 
            self.flux_plus_dLij_times_soln,
            self.dL_minus_dC, 
            self.min_u_bc,
            self.max_u_bc,
            self.quantDOFs) #TMP

        if self.forceStrongConditions:
            for dofN, g in self.dirichletConditionsForceDOF.DOFBoundaryConditionsDict.iteritems():
                #if self.coefficients.vely_tn_dof[dofN] < 0: #HACKED FOR THIS PROBLEM (MQL)
                if self.inflow_DOFs[dofN] == 1.0: #HACKED FOR THIS PROBLEM (MQL) TMP
                    r[dofN] = 0
        if self.stabilization:
            self.stabilization.accumulateSubgridMassHistory(self.q)
        log("Global residual", level=9, data=r)
        # mwf debug
        # pdb.set_trace()
        # mwf decide if this is reasonable for keeping solver statistics
        self.nonlinear_function_evaluations += 1
        if self.globalResidualDummy is None:
            self.globalResidualDummy = numpy.zeros(r.shape, 'd')

    def getJacobian(self, jacobian):
        cfemIntegrals.zeroJacobian_CSR(self.nNonzerosInJacobian,
                                       jacobian)
        self.vof.calculateJacobian(  # element
            self.u[0].femSpace.elementMaps.psi,
            self.u[0].femSpace.elementMaps.grad_psi,
            self.mesh.nodeArray,
            self.mesh.nodeVelocityArray,
            self.MOVING_DOMAIN,
            self.mesh.elementNodesArray,
            self.elementQuadratureWeights[('u', 0)],
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            # element boundary
            self.u[0].femSpace.elementMaps.psi_trace,
            self.u[0].femSpace.elementMaps.grad_psi_trace,
            self.elementBoundaryQuadratureWeights[('u', 0)],
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.elementMaps.boundaryNormals,
            self.u[0].femSpace.elementMaps.boundaryJacobians,
            self.mesh.nElements_global,
            self.coefficients.useMetrics,
            self.timeIntegration.alpha_bdf,
            self.shockCapturing.lag,
            self.shockCapturing.shockCapturingFactor,
            # VRANS start
            self.coefficients.q_vos,
            # VRANS end
            self.u[0].femSpace.dofMap.l2g,
            self.mesh.elementDiametersArray,
            self.u[0].dof,
            self.coefficients.q_v,
            self.timeIntegration.beta_bdf[0],
            self.q[('cfl', 0)],
            self.shockCapturing.numDiff_last[0],
            self.csrRowIndeces[(0, 0)], self.csrColumnOffsets[(0, 0)],
            jacobian,
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray,
            self.coefficients.ebqe_v,
            # VRANS start
            self.coefficients.ebqe_vos,
            # VRANS end
            self.numericalFlux.isDOFBoundary[0],
            self.numericalFlux.ebqe[('u', 0)],
            self.ebqe[('advectiveFlux_bc_flag', 0)],
            self.ebqe[('advectiveFlux_bc', 0)],
            self.csrColumnOffsets_eb[(0, 0)], 
            self.coefficients.EDGE_VISCOSITY,
            self.coefficients.ENTROPY_VISCOSITY,
            self.coefficients.LUMPED_MASS_MATRIX)
        # Load the Dirichlet conditions directly into residual
        if self.forceStrongConditions:
            scaling = 1.0  # probably want to add some scaling to match non-dirichlet diagonals in linear system
            for dofN in self.dirichletConditionsForceDOF.DOFBoundaryConditionsDict.keys():
                global_dofN = dofN
                #if self.coefficients.vely_tn_dof[dofN] < 0: #HACKED FOR THIS PROBLEM (MQL)
                if self.inflow_DOFs[dofN] == 1.0: #HACKED FOR THIS PROBLEM (MQL) TMP
                    for i in range(
                        self.rowptr[global_dofN],
                        self.rowptr[
                            global_dofN + 1]):
                        if (self.colind[i] == global_dofN):
                            # print "RBLES forcing residual cj = %s dofN= %s
                            # global_dofN= %s was self.nzval[i]= %s now =%s " %
                            # (cj,dofN,global_dofN,self.nzval[i],scaling)
                            self.nzval[i] = scaling
                        else:
                            self.nzval[i] = 0.0
                            # print "RBLES zeroing residual cj = %s dofN= %s
                            # global_dofN= %s " % (cj,dofN,global_dofN)
        log("Jacobian ", level=10, data=jacobian)
        # mwf decide if this is reasonable for solver statistics
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
        self.u[0].femSpace.elementMaps.getBasisValuesRef(
            self.elementQuadraturePoints)
        self.u[0].femSpace.elementMaps.getBasisGradientValuesRef(
            self.elementQuadraturePoints)
        self.u[0].femSpace.getBasisValuesRef(self.elementQuadraturePoints)
        self.u[0].femSpace.getBasisGradientValuesRef(
            self.elementQuadraturePoints)
        self.coefficients.initializeElementQuadrature(
            self.timeIntegration.t, self.q)
        if self.stabilization is not None:
            self.stabilization.initializeElementQuadrature(
                self.mesh, self.timeIntegration.t, self.q)
            self.stabilization.initializeTimeIntegration(self.timeIntegration)
        if self.shockCapturing is not None:
            self.shockCapturing.initializeElementQuadrature(
                self.mesh, self.timeIntegration.t, self.q)

    def calculateElementBoundaryQuadrature(self):
        pass

    def calculateExteriorElementBoundaryQuadrature(self):
        """
        Calculate the physical location and weights of the quadrature rules
        and the shape information at the quadrature points on global element boundaries.

        This function should be called only when the mesh changes.
        """
        #
        # get physical locations of element boundary quadrature points
        #
        # assume all components live on the same mesh
        self.u[0].femSpace.elementMaps.getBasisValuesTraceRef(
            self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.elementMaps.getBasisGradientValuesTraceRef(
            self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.getBasisValuesTraceRef(
            self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.getBasisGradientValuesTraceRef(
            self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.elementMaps.getValuesGlobalExteriorTrace(
            self.elementBoundaryQuadraturePoints, self.ebqe['x'])
        self.fluxBoundaryConditionsObjectsDict = dict([(cj, FluxBoundaryConditions(self.mesh,
                                                                                   self.nElementBoundaryQuadraturePoints_elementBoundary,
                                                                                   self.ebqe[('x')],
                                                                                   self.advectiveFluxBoundaryConditionsSetterDict[cj],
                                                                                   self.diffusiveFluxBoundaryConditionsSetterDictDict[cj]))
                                                       for cj in self.advectiveFluxBoundaryConditionsSetterDict.keys()])
        self.coefficients.initializeGlobalExteriorElementBoundaryQuadrature(
            self.timeIntegration.t, self.ebqe)

    def estimate_mt(self):
        pass

    def calculateSolutionAtQuadrature(self):
        pass

    def calculateAuxiliaryQuantitiesAfterStep(self):
        pass

    def updateAfterMeshMotion(self):
        pass

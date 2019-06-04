#ifndef SW2DCV_H
#define SW2DCV_H
#include <cmath>
#include <iostream>
#include <valarray>
#include "CompKernel.h"
#include "ModelFactory.h"
#include <assert.h>

//cek todo
//2. Get stabilization right
//3. Add Riemann solvers for external flux
//4. Add Riemann solvers for internal flux and DG terms
//5. Try other choices of variables h,hu,hv, Bova-Carey symmetrization?

#define GLOBAL_FCT 0
#define POWER_SMOOTHNESS_INDICATOR 2
#define VEL_FIX_POWER 2.
#define REESTIMATE_MAX_EDGE_BASED_CFL 1

namespace proteus
{
  // FOR CELL BASED ENTROPY VISCOSITY
  inline double ENTROPY(const double& g, const double& h, const double& hu, const double& hv, const double& z, const double& one_over_hReg){
    return 0.5*(g*h*h + one_over_hReg*(hu*hu+hv*hv) + 2.*g*h*z);
  }
  inline double DENTROPY_DH(const double& g, const double& h, const double& hu, const double& hv, const double& z, const double& one_over_hReg){
    return g*h - 0.5*(hu*hu+hv*hv)*std::pow(one_over_hReg,2) + g*z;
  }
  inline double DENTROPY_DHU(const double& g, const double& h, const double& hu, const double& hv, const double& z, const double& one_over_hReg){
    return hu*one_over_hReg;
  }
  inline double DENTROPY_DHV(const double& g, const double& h, const double& hu, const double& hv, const double& z, const double& one_over_hReg){
    return hv*one_over_hReg;
  }
  inline double ENTROPY_FLUX1(const double& g, const double& h, const double& hu, const double& hv, const double& z, const double& one_over_hReg){
    return (ENTROPY(g,h,hu,hv,z,one_over_hReg) + 0.5*g*h*h + g*h*z)*hu*one_over_hReg;
  }
  inline double ENTROPY_FLUX2(const double& g, const double& h, const double& hu, const double& hv, const double& z, const double& one_over_hReg){
    return (ENTROPY(g,h,hu,hv,z,one_over_hReg) + 0.5*g*h*h + g*h*z)*hv*one_over_hReg;
  }
  // FOR ESTIMATING MAX WAVE SPEEDS
  inline double f(const double& g, const double& h, const double& hZ){
    return ( (h <= hZ) ? 2.*(sqrt(g*h)-sqrt(g*hZ)) : (h-hZ)*sqrt(0.5*g*(h+hZ)/h/hZ) );
  }
  inline double phi(const double& g, const double& h, const double& hL, const double& hR, const double& uL, const double& uR){
    return ( f(g,h,hL) + f(g,h,hR) + uR - uL );
  }
  inline double fp(const double& g, const double& h, const double& hZ){
    return ( (h <= hZ) ? sqrt(g/h) : g*(2*h*h+h*hZ+hZ*hZ)/(2*sqrt(2*g)*h*h*hZ*sqrt(1/h+1/hZ)) );
  }
  inline double phip(const double& g, const double& h, const double& hL, const double& hR){
    return ( fp(g,h,hL) + fp(g,h,hR) );
  }
  inline double nu1(const double& g, const double& hStar, const double& hL, const double& uL){
    return ( uL - sqrt(g*hL)*sqrt( (1+fmax((hStar-hL)/2/hL,0.0)) * (1+fmax((hStar-hL)/hL,0.)) ) );
  }
  inline double nu3(const double& g, const double& hStar, const double& hR, const double& uR){
    return ( uR + sqrt(g*hR)*sqrt( (1+fmax((hStar-hR)/2/hR,0.0)) * (1+fmax((hStar-hR)/hR,0.)) ) );
  }
  inline double phiDiff(const double& g, const double& h1k, const double& h2k, const double& hL, const double& hR, const double& uL, const double& uR){
    return ( (phi(g,h2k,hL,hR,uL,uR) - phi(g,h1k,hL,hR,uL,uR))/(h2k-h1k) );
  }
  inline double phiDDiff1(const double& g, const double& h1k, const double& h2k, const double& hL, const double& hR, const double& uL, const double& uR){
  return ( (phiDiff(g,h1k,h2k,hL,hR,uL,uR) - phip(g,h1k,hL,hR))/(h2k-h1k) );
  }
  inline double phiDDiff2(const double& g, const double& h1k, const double& h2k, const double& hL, const double& hR, const double& uL, const double& uR){
  return ( (phip(g,h2k,hL,hR) - phiDiff(g,h1k,h2k,hL,hR,uL,uR))/(h2k-h1k) );
  }
  inline double hStarLFromQuadPhiFromAbove(const double& g, const double& hStarL, const double& hStarR, const double& hL, const double& hR, const double& uL, const double& uR){
    return ( hStarL-2*phi(g,hStarL,hL,hR,uL,uR)/(phip(g,hStarL,hL,hR)+sqrt(std::pow(phip(g,hStarL,hL,hR),2)-4*phi(g,hStarL,hL,hR,uL,uR)*phiDDiff1(g,hStarL,hStarR,hL,hR,uL,uR))) );
  }
  inline double hStarRFromQuadPhiFromBelow(const double& g, const double& hStarL, const double& hStarR, const double& hL, const double& hR, const double& uL, const double& uR){
    return ( hStarR-2*phi(g,hStarR,hL,hR,uL,uR)/(phip(g,hStarR,hL,hR)+sqrt(std::pow(phip(g,hStarR,hL,hR),2)-4*phi(g,hStarR,hL,hR,uL,uR)*phiDDiff2(g,hStarL,hStarR,hL,hR,uL,uR))) );
  }
}

namespace proteus
{
  class SW2DCV_base
  {
  public:
    std::valarray<double> Rneg, Rpos, hLow, huLow, hvLow, Kmax;
    virtual ~SW2DCV_base(){}
    virtual void FCTStep(double dt,
			 int NNZ, //number on non-zero entries on sparsity pattern
			 int numDOFs, //number of DOFs
			 double* lumped_mass_matrix, //lumped mass matrix (as vector)
			 double* h_old, //DOFs of solution at last stage
			 double* hu_old,
			 double* hv_old,
			 double* b_dof,
			 double* high_order_hnp1, //DOFs of high order solution at tnp1
			 double* high_order_hunp1,
			 double* high_order_hvnp1,
			 double* extendedSourceTerm_hu,
			 double* extendedSourceTerm_hv,
			 double* limited_hnp1,
			 double* limited_hunp1,
			 double* limited_hvnp1,
			 int* csrRowIndeces_DofLoops, //csr row indeces
			 int* csrColumnOffsets_DofLoops, //csr column offsets
			 double* MassMatrix, //mass matrix
			 double* dH_minus_dL,
			 double* muH_minus_muL,
			 double hEps,
			 double* hReg,
			 int LUMPED_MASS_MATRIX,
			 // LOCAL LIMITING
			 double* dLow,
			 double* hBT,
			 double* huBT,
			 double* hvBT
			 )=0;
    virtual void convexLimiting(double dt,
				int NNZ, //number on non-zero entries on sparsity pattern
				int numDOFs, //number of DOFs
				double* lumped_mass_matrix, //lumped mass matrix (as vector)
				double* h_old, //DOFs of solution at last stage
				double* hu_old,
				double* hv_old,
				double* b_dof,
				double* high_order_hnp1, //DOFs of high order solution at tnp1
				double* high_order_hunp1,
				double* high_order_hvnp1,
				double* extendedSourceTerm_hu,
				double* extendedSourceTerm_hv,
				double* limited_hnp1,
				double* limited_hunp1,
				double* limited_hvnp1,
				int* csrRowIndeces_DofLoops, //csr row indeces
				int* csrColumnOffsets_DofLoops, //csr column offsets
				double* MassMatrix, //mass matrix
				double* dH_minus_dL,
				double* muH_minus_muL,
				double hEps,
				double* hReg,
				int LUMPED_MASS_MATRIX,
				// LOCAL LIMITING
				double* dLow,
				double* hBT,
				double* huBT,
				double* hvBT
				)=0;
    virtual double calculateEdgeBasedCFL(double g,
					 int numDOFsPerEqn, //number of DOFs
					 double* lumped_mass_matrix, //lumped mass matrix (as vector)
					 double* h_old, //DOFs of solution at last stage
					 double* hu_old,
					 double* hv_old,
					 double* b_dof,
					 int* csrRowIndeces_DofLoops, //csr row indeces
					 int* csrColumnOffsets_DofLoops, //csr column offsets
					 double hEps,
					 double* hReg,
					 double* Cx,
					 double* Cy,
					 double* CTx,
					 double* CTy,
					 double* dLow,
					 double run_cfl,
					 double* edge_based_cfl,
					 int debug
					 )=0;
    virtual void calculateResidual(// last EDGE BASED version
                                   double* mesh_trial_ref,
                                   double* mesh_grad_trial_ref,
                                   double* mesh_dof,
                                   double* mesh_velocity_dof,
                                   double MOVING_DOMAIN,//0 or 1
                                   int* mesh_l2g,
                                   double* dV_ref,
                                   double* h_trial_ref,
                                   double* h_grad_trial_ref,
                                   double* h_test_ref,
                                   double* h_grad_test_ref,
                                   double* vel_trial_ref,
                                   double* vel_grad_trial_ref,
                                   double* vel_test_ref,
                                   double* vel_grad_test_ref,
                                   //element boundary
                                   double* mesh_trial_trace_ref,
                                   double* mesh_grad_trial_trace_ref,
                                   double* dS_ref,
                                   double* h_trial_trace_ref,
                                   double* h_grad_trial_trace_ref,
                                   double* h_test_trace_ref,
                                   double* h_grad_test_trace_ref,
                                   double* vel_trial_trace_ref,
                                   double* vel_grad_trial_trace_ref,
                                   double* vel_test_trace_ref,
                                   double* vel_grad_test_trace_ref,
                                   double* normal_ref,
                                   double* boundaryJac_ref,
                                   //physics
                                   double* elementDiameter,
                                   int nElements_global,
                                   double useRBLES,
                                   double useMetrics,
                                   double alphaBDF,
                                   double nu,
                                   double g,
                                   int* h_l2g,
                                   int* vel_l2g,
                                   double* h_dof_old,
                                   double* hu_dof_old,
                                   double* hv_dof_old,
                                   double* b_dof,
                                   double* h_dof,
                                   double* hu_dof,
                                   double* hv_dof,
                                   double* h_dof_sge,
                                   double* hu_dof_sge,
                                   double* hv_dof_sge,
                                   double* q_mass_acc,
                                   double* q_mom_hu_acc,
                                   double* q_mom_hv_acc,
                                   double* q_mass_adv,
                                   double* q_mass_acc_beta_bdf,
                                   double* q_mom_hu_acc_beta_bdf,
                                   double* q_mom_hv_acc_beta_bdf,
                                   double* q_cfl,
                                   int* sdInfo_hu_hu_rowptr,
                                   int* sdInfo_hu_hu_colind,
                                   int* sdInfo_hu_hv_rowptr,
                                   int* sdInfo_hu_hv_colind,
                                   int* sdInfo_hv_hv_rowptr,
                                   int* sdInfo_hv_hv_colind,
                                   int* sdInfo_hv_hu_rowptr,
                                   int* sdInfo_hv_hu_colind,
                                   int offset_h,
                                   int offset_hu,
                                   int offset_hv,
                                   int stride_h,
                                   int stride_hu,
                                   int stride_hv,
                                   double* globalResidual,
                                   int nExteriorElementBoundaries_global,
                                   int* exteriorElementBoundariesArray,
                                   int* elementBoundaryElementsArray,
                                   int* elementBoundaryLocalElementBoundariesArray,
                                   int* isDOFBoundary_h,
                                   int* isDOFBoundary_hu,
                                   int* isDOFBoundary_hv,
                                   int* isAdvectiveFluxBoundary_h,
                                   int* isAdvectiveFluxBoundary_hu,
                                   int* isAdvectiveFluxBoundary_hv,
                                   int* isDiffusiveFluxBoundary_hu,
                                   int* isDiffusiveFluxBoundary_hv,
                                   double* ebqe_bc_h_ext,
                                   double* ebqe_bc_flux_mass_ext,
                                   double* ebqe_bc_flux_mom_hu_adv_ext,
                                   double* ebqe_bc_flux_mom_hv_adv_ext,
                                   double* ebqe_bc_hu_ext,
                                   double* ebqe_bc_flux_hu_diff_ext,
                                   double* ebqe_penalty_ext,
                                   double* ebqe_bc_hv_ext,
                                   double* ebqe_bc_flux_hv_diff_ext,
                                   double* q_velocity,
                                   double* ebqe_velocity,
                                   double* flux,
                                   double* elementResidual_h,
                                   // C matrices
                                   double* Cx,
                                   double* Cy,
                                   double* CTx,
                                   double* CTy,
                                   // PARAMETERS FOR EDGE BASED STABILIZATION
                                   int numDOFsPerEqn,
                                   int NNZ,
                                   int* csrRowIndeces_DofLoops,
                                   int* csrColumnOffsets_DofLoops,
                                   // LUMPED MASS MATRIX
                                   double* lumped_mass_matrix,
                                   double cfl_run,
                                   double hEps,
                                   double* hReg,
                                   // SAVE SOLUTION (mql)
                                   double* hnp1_at_quad_point,
                                   double* hunp1_at_quad_point,
                                   double* hvnp1_at_quad_point,
                                   // TO COMPUTE LOW ORDER
				   double* extendedSourceTerm_hu,
				   double* extendedSourceTerm_hv,
                                   // FOR FCT
                                   double* dH_minus_dL,
                                   double* muH_minus_muL,
                                   double cE,
                                   int LUMPED_MASS_MATRIX,
                                   double dt,
                                   int LINEAR_FRICTION,
                                   double mannings,
                                   // Quant of interests
                                   double* quantDOFs,
                                   int SECOND_CALL_CALCULATE_RESIDUAL,
                                   // NORMAL COMPONENTS
                                   int COMPUTE_NORMALS,
                                   double* normalx,
                                   double* normaly,
                                   double* dLow,
				   double* hBT,
				   double* huBT,
				   double* hvBT,
                                   int lstage
                                   )=0;
    virtual void calculateMassMatrix(//element
                                     double* mesh_trial_ref,
                                     double* mesh_grad_trial_ref,
                                     double* mesh_dof,
                                     double* mesh_velocity_dof,
                                     double MOVING_DOMAIN,
                                     int* mesh_l2g,
                                     double* dV_ref,
                                     double* h_trial_ref,
                                     double* h_grad_trial_ref,
                                     double* h_test_ref,
                                     double* h_grad_test_ref,
                                     double* vel_trial_ref,
                                     double* vel_grad_trial_ref,
                                     double* vel_test_ref,
                                     double* vel_grad_test_ref,
                                     //element boundary
                                     double* mesh_trial_trace_ref,
                                     double* mesh_grad_trial_trace_ref,
                                     double* dS_ref,
                                     double* h_trial_trace_ref,
                                     double* h_grad_trial_trace_ref,
                                     double* h_test_trace_ref,
                                     double* h_grad_test_trace_ref,
                                     double* vel_trial_trace_ref,
                                     double* vel_grad_trial_trace_ref,
                                     double* vel_test_trace_ref,
                                     double* vel_grad_test_trace_ref,
                                     double* normal_ref,
                                     double* boundaryJac_ref,
                                     //physics
                                     double* elementDiameter,
                                     int nElements_global,
                                     double useRBLES,
                                     double useMetrics,
                                     double alphaBDF,
                                     double nu,
                                     double g,
                                     int* h_l2g,
                                     int* vel_l2g,
                                     double* b_dof,
                                     double* h_dof,
                                     double* hu_dof,
                                     double* hv_dof,
                                     double* h_dof_sge,
                                     double* hu_dof_sge,
                                     double* hv_dof_sge,
                                     double* q_mass_acc_beta_bdf,
                                     double* q_mom_hu_acc_beta_bdf,
                                     double* q_mom_hv_acc_beta_bdf,
                                     double* q_cfl,
                                     int* sdInfo_hu_hu_rowptr,
                                     int* sdInfo_hu_hu_colind,
                                     int* sdInfo_hu_hv_rowptr,
                                     int* sdInfo_hu_hv_colind,
                                     int* sdInfo_hv_hv_rowptr,
                                     int* sdInfo_hv_hv_colind,
                                     int* sdInfo_hv_hu_rowptr,
                                     int* sdInfo_hv_hu_colind,
                                     int* csrRowIndeces_h_h,
                                     int* csrColumnOffsets_h_h,
                                     int* csrRowIndeces_h_hu,
                                     int* csrColumnOffsets_h_hu,
                                     int* csrRowIndeces_h_hv,
                                     int* csrColumnOffsets_h_hv,
                                     int* csrRowIndeces_hu_h,
                                     int* csrColumnOffsets_hu_h,
                                     int* csrRowIndeces_hu_hu,
                                     int* csrColumnOffsets_hu_hu,
                                     int* csrRowIndeces_hu_hv,
                                     int* csrColumnOffsets_hu_hv,
                                     int* csrRowIndeces_hv_h,
                                     int* csrColumnOffsets_hv_h,
                                     int* csrRowIndeces_hv_hu,
                                     int* csrColumnOffsets_hv_hu,
                                     int* csrRowIndeces_hv_hv,
                                     int* csrColumnOffsets_hv_hv,
                                     double* globalJacobian,
                                     int nExteriorElementBoundaries_global,
                                     int* exteriorElementBoundariesArray,
                                     int* elementBoundaryElementsArray,
                                     int* elementBoundaryLocalElementBoundariesArray,
                                     int* isDOFBoundary_h,
                                     int* isDOFBoundary_hu,
                                     int* isDOFBoundary_hv,
                                     int* isAdvectiveFluxBoundary_h,
                                     int* isAdvectiveFluxBoundary_hu,
                                     int* isAdvectiveFluxBoundary_hv,
                                     int* isDiffusiveFluxBoundary_hu,
                                     int* isDiffusiveFluxBoundary_hv,
                                     double* ebqe_bc_h_ext,
                                     double* ebqe_bc_flux_mass_ext,
                                     double* ebqe_bc_flux_mom_hu_adv_ext,
                                     double* ebqe_bc_flux_mom_hv_adv_ext,
                                     double* ebqe_bc_hu_ext,
                                     double* ebqe_bc_flux_hu_diff_ext,
                                     double* ebqe_penalty_ext,
                                     double* ebqe_bc_hv_ext,
                                     double* ebqe_bc_flux_hv_diff_ext,
                                     int* csrColumnOffsets_eb_h_h,
                                     int* csrColumnOffsets_eb_h_hu,
                                     int* csrColumnOffsets_eb_h_hv,
                                     int* csrColumnOffsets_eb_hu_h,
                                     int* csrColumnOffsets_eb_hu_hu,
                                     int* csrColumnOffsets_eb_hu_hv,
                                     int* csrColumnOffsets_eb_hv_h,
                                     int* csrColumnOffsets_eb_hv_hu,
                                     int* csrColumnOffsets_eb_hv_hv,
                                     double dt)=0;
  virtual void calculateLumpedMassMatrix(//element
                                         double* mesh_trial_ref,
                                         double* mesh_grad_trial_ref,
                                         double* mesh_dof,
                                         double* mesh_velocity_dof,
                                         double MOVING_DOMAIN,
                                         int* mesh_l2g,
                                         double* dV_ref,
                                         double* h_trial_ref,
                                         double* h_grad_trial_ref,
                                         double* h_test_ref,
                                         double* h_grad_test_ref,
                                         double* vel_trial_ref,
                                         double* vel_grad_trial_ref,
                                         double* vel_test_ref,
                                         double* vel_grad_test_ref,
                                         //element boundary
                                         double* mesh_trial_trace_ref,
                                         double* mesh_grad_trial_trace_ref,
                                         double* dS_ref,
                                         double* h_trial_trace_ref,
                                         double* h_grad_trial_trace_ref,
                                         double* h_test_trace_ref,
                                         double* h_grad_test_trace_ref,
                                         double* vel_trial_trace_ref,
                                         double* vel_grad_trial_trace_ref,
                                         double* vel_test_trace_ref,
                                         double* vel_grad_test_trace_ref,
                                         double* normal_ref,
                                         double* boundaryJac_ref,
                                         //physics
                                         double* elementDiameter,
                                         int nElements_global,
                                         double useRBLES,
                                         double useMetrics,
                                         double alphaBDF,
                                         double nu,
                                         double g,
                                         int* h_l2g,
                                         int* vel_l2g,
                                         double* b_dof,
                                         double* h_dof,
                                         double* hu_dof,
                                         double* hv_dof,
                                         double* h_dof_sge,
                                         double* hu_dof_sge,
                                         double* hv_dof_sge,
                                         double* q_mass_acc_beta_bdf,
                                         double* q_mom_hu_acc_beta_bdf,
                                         double* q_mom_hv_acc_beta_bdf,
                                         double* q_cfl,
                                         int* sdInfo_hu_hu_rowptr,
                                         int* sdInfo_hu_hu_colind,
                                         int* sdInfo_hu_hv_rowptr,
                                         int* sdInfo_hu_hv_colind,
                                         int* sdInfo_hv_hv_rowptr,
                                         int* sdInfo_hv_hv_colind,
                                         int* sdInfo_hv_hu_rowptr,
                                         int* sdInfo_hv_hu_colind,
                                         int* csrRowIndeces_h_h,
                                         int* csrColumnOffsets_h_h,
                                         int* csrRowIndeces_h_hu,
                                         int* csrColumnOffsets_h_hu,
                                         int* csrRowIndeces_h_hv,
                                         int* csrColumnOffsets_h_hv,
                                         int* csrRowIndeces_hu_h,
                                         int* csrColumnOffsets_hu_h,
                                         int* csrRowIndeces_hu_hu,
                                         int* csrColumnOffsets_hu_hu,
                                         int* csrRowIndeces_hu_hv,
                                         int* csrColumnOffsets_hu_hv,
                                         int* csrRowIndeces_hv_h,
                                         int* csrColumnOffsets_hv_h,
                                         int* csrRowIndeces_hv_hu,
                                         int* csrColumnOffsets_hv_hu,
                                         int* csrRowIndeces_hv_hv,
                                         int* csrColumnOffsets_hv_hv,
                                         double* globalJacobian,
                                         int nExteriorElementBoundaries_global,
                                         int* exteriorElementBoundariesArray,
                                         int* elementBoundaryElementsArray,
                                         int* elementBoundaryLocalElementBoundariesArray,
                                         int* isDOFBoundary_h,
                                         int* isDOFBoundary_hu,
                                         int* isDOFBoundary_hv,
                                         int* isAdvectiveFluxBoundary_h,
                                         int* isAdvectiveFluxBoundary_hu,
                                         int* isAdvectiveFluxBoundary_hv,
                                         int* isDiffusiveFluxBoundary_hu,
                                         int* isDiffusiveFluxBoundary_hv,
                                         double* ebqe_bc_h_ext,
                                         double* ebqe_bc_flux_mass_ext,
                                         double* ebqe_bc_flux_mom_hu_adv_ext,
                                         double* ebqe_bc_flux_mom_hv_adv_ext,
                                         double* ebqe_bc_hu_ext,
                                         double* ebqe_bc_flux_hu_diff_ext,
                                         double* ebqe_penalty_ext,
                                         double* ebqe_bc_hv_ext,
                                         double* ebqe_bc_flux_hv_diff_ext,
                                         int* csrColumnOffsets_eb_h_h,
                                         int* csrColumnOffsets_eb_h_hu,
                                         int* csrColumnOffsets_eb_h_hv,
                                         int* csrColumnOffsets_eb_hu_h,
                                         int* csrColumnOffsets_eb_hu_hu,
                                         int* csrColumnOffsets_eb_hu_hv,
                                         int* csrColumnOffsets_eb_hv_h,
                                         int* csrColumnOffsets_eb_hv_hu,
                                         int* csrColumnOffsets_eb_hv_hv,
                                         double dt)=0;
  };

  template<class CompKernelType,
           int nSpace,
           int nQuadraturePoints_element,
           int nDOF_mesh_trial_element,
           int nDOF_trial_element,
           int nDOF_test_element,
           int nQuadraturePoints_elementBoundary>
  class SW2DCV : public SW2DCV_base
  {
  public:
    const int nDOF_test_X_trial_element;
    CompKernelType ck;
    SW2DCV():
      nDOF_test_X_trial_element(nDOF_test_element*nDOF_trial_element),
      ck()
    {
      std::cout<<"Constructing SW2DCV<CompKernelTemplate<"
               <<nSpace<<","
               <<nQuadraturePoints_element<<","
               <<nDOF_mesh_trial_element<<","
               <<nDOF_trial_element<<","
               <<nDOF_test_element<<","
               <<nQuadraturePoints_elementBoundary<<">());"
               <<std::endl<<std::flush;
    }

    inline
      double maxWaveSpeedSharpInitialGuess(double g, double nx, double ny,
                                           double hL, double huL, double hvL,
                                           double hR, double huR, double hvR,
                                           double hEpsL, double hEpsR,
                                           bool debugging)
    {
      double lambda1, lambda3;
      //1-eigenvalue: uL-sqrt(g*hL)
      //3-eigenvalue: uR+sqrt(g*hR)

      double hVelL = nx*huL + ny*hvL;
      double hVelR = nx*huR + ny*hvR;
      double velL = 2*hL/(hL*hL+std::pow(fmax(hL,hEpsL),2))*hVelL;
      double velR = 2*hR/(hR*hR+std::pow(fmax(hR,hEpsR),2))*hVelR;

      if (debugging)
        std::cout << "hL, hR, hVelL, hVelR, velL, velR: "
                  << hL << "\t"
                  << hR << "\t"
                  << hVelL <<  "\t"
                  << hVelR <<  "\t"
                  << velL <<  "\t"
                  << velR <<  "\t"
                  << std::endl;
      // CHECK IF BOTH STATES ARE DRY:
      if (hL==0 && hR==0)
        {
          lambda1=0.;
          lambda3=0.;
        }
      else if (hL==0) // left dry state
        {
          lambda1 = velR-2*sqrt(g*hR);
          lambda3 = velR+sqrt(g*hR);
          if (debugging)
            {
              std::cout << "hL=0" << std::endl;
              std::cout << lambda1 << "\t" << lambda3 << std::endl;
            }
        }
      else if (hR==0) // right dry state
        {
          lambda1 = velL-sqrt(g*hL);
          lambda3 = velL+2*sqrt(g*hL);
          if (debugging)
            {
              std::cout << "hR=0" << std::endl;
              std::cout << lambda1 << "\t" << lambda3 << std::endl;
            }
        }
      else // both states are wet
        {
          double x0 = std::pow(2.*sqrt(2.)-1.,2.);
          double hMin = fmin(hL,hR);
          double hMax = fmax(hL,hR);

          double hStar;
          double fMin = phi(g,x0*hMin,hL,hR,velL,velR);
          double fMax = phi(g,x0*hMax,hL,hR,velL,velR);

          if (debugging)
            std::cout << "hMin, hMax, fMin, fMax: "
                      << hMin << ", " << hMax << ", "
                      << fMin << ", " << fMax
                      << std::endl;

          if (0 <= fMin)
            {
              hStar = std::pow(fmax(0.,velL-velR+2*sqrt(g)*(sqrt(hL)+sqrt(hR))),2)/16./g;
              if (debugging)
                std::cout << "**********... THIS IS A RAREFACTION"
                          << std::endl;
              if (debugging)
                {
                  std::cout << "h* = " << hStar << std::endl;
                  lambda1 = nu1(g,hStar,hL,velL);
                  lambda3 = nu3(g,hStar,hR,velR);
                  std::cout << "lambda1, lambda3: " << lambda1 << ", " << lambda3 << std::endl;
                }
            }
          else if (0 <= fMax)
            hStar = std::pow(-sqrt(2*hMin)+sqrt(3*hMin+2*sqrt(2*hMin*hMax)+sqrt(2./g)*(velL-velR)*sqrt(hMin)),2);
          else // fMax < 0
            hStar = sqrt(hMin*hMax)*(1+(sqrt(2)*(velL-velR))/(sqrt(g*hMin)+sqrt(g*hMax)));
          // Compute max wave speed based on hStar0
          lambda1 = nu1(g,hStar,hL,velL);
          lambda3 = nu3(g,hStar,hR,velR);
        }
      if (debugging)
	{
	  std::cout << "lambda1, lambda3: " << lambda1 << ", " << lambda3 << std::endl;
	  if (isinf(lambda1)==1 || isinf(lambda3)==1)
	    abort();
	}
      //return fmax(fmax(0.,-lambda1), fmax(0,lambda3));
      return fmax(lambda1, lambda3);
    }

    inline
      double maxWaveSpeedIterativeProcess(double g, double nx, double ny,
                                          double hL, double huL, double hvL,
                                          double hR, double huR, double hvR,
                                          double hEpsL, double hEpsR,
                                          bool verbose)
    {
      double tol = 1E-15;
      //1-eigenvalue: uL-sqrt(g*hL)
      //3-eigenvalue: uR+sqrt(g*hR)

      double hVelL = nx*huL + ny*hvL;
      double hVelR = nx*huR + ny*hvR;
      double velL = 2*hL/(hL*hL+std::pow(fmax(hL,hEpsL),2))*hVelL;
      double velR = 2*hR/(hR*hR+std::pow(fmax(hR,hEpsR),2))*hVelR;

      double lambda1, lambda3;

      // CHECK IF BOTH STATES ARE DRY:
      if (hL==0 && hR==0)
        {
          lambda1=0.;
          lambda3=0.;
          return 0.;
        }
      else if (hL==0) // left dry state
        {
          lambda1 = velR-2*sqrt(g*hR);
          lambda3 = velR+sqrt(g*hR);
          return fmax(fabs(lambda1),fabs(lambda3));
        }
      else if (hR==0) // right dry state
        {
          lambda1 = velL-sqrt(g*hL);
          lambda3 = velL+2*sqrt(g*hL);
          return fmax(fabs(lambda1),fabs(lambda3));
        }
      else
        {
          ////////////////////
          // ESTIMATE hStar //
          ////////////////////
          // Initial estimate of hStar0 from above.
          // This is computed via phiR(h) >= phi(h) ---> hStar0 >= hStar
          double hStar0 = 1;
          double hStar = hStar0;

          /////////////////////////////////
          // ALGORITHM 1: Initialization //
          /////////////////////////////////
          // Requires: tol
          // Ensures: hStarL, hStarR
          double hStarL, hStarR;
          double hMin = fmin(hL,hR);
          double hMax = fmin(hL,hR);
          double phiMin = phi(g,hMin,hL,hR,velL,velR);
          double phiMax = phi(g,hMax,hL,hR,velL,velR);
          if (0 <= phiMin)
            {
              // This is a 1- and 3-rarefactions situation. We know the solution in this case
              lambda1 = velL - sqrt(g*hL);
              lambda3 = velR + sqrt(g*hR);

              std::cout << "lambda Min, lambda Max: "
                        << lambda1 << ", " << lambda3
                        << std::endl;

              return fmax(fabs(lambda1),fabs(lambda3));
            }
          if (phiMax == 0) // if hMax "hits" hStar (very unlikely)
            {
              hStar = hMax;
              lambda1 = nu1(g,hStar,hL,velL);
              lambda3 = nu3(g,hStar,hR,velR);
              return fmax(fabs(lambda1),fabs(lambda3));
            }
          double hStarTwoRarefactions = std::pow(velL-velR+2*sqrt(g)*(sqrt(hL)+sqrt(hR)),2)/16/g;
          if (phiMax < 0) // This is a 1- and 3-shock situation
            {
              hStarL = hMax;
              hStarR = hStarTwoRarefactions;
            }
          else // Here we have one shock and one rarefaction
            {
              hStarL = hMin;
              hStarR = fmin(hMax,hStarTwoRarefactions);
            }

          // improve estimate from below via one newton step (not required)
          hStarL = fmax(hStarL,hStarR-phi(g,hStarR,hL,hR,velL,velR)/phip(g,hStarR,hL,hR));
          // COMPUTE lambdaMin0 and lambdaMax0
          double nu11 = nu1(g,hStarR,hL,velL);
          double nu12 = nu1(g,hStarL,hL,velL);
          double nu31 = nu3(g,hStarL,hR,velR);
          double nu32 = nu3(g,hStarR,hR,velR);

          double lambdaMin = fmax(fmax(nu31,0), fmax(-nu12,0));
          double lambdaMax = fmax(fmax(nu32,0), fmax(-nu11,0));

          if (verbose)
            {
              std::cout << "hStarL, hStarR: " << hStarL << ", " << hStarR << "\t"
                        << "lambda Min, lambda Max: "
                        << lambdaMin << ", " << lambdaMax
                        << std::endl;
            }
          // CHECK IF TOL IS SATISFIED. O.W. GOT TO ALGORITHM 2 //
          if (lambdaMin > 0 && lambdaMax/lambdaMin - 1 <= tol)
            return lambdaMax;
          else // Proceed to algorithm 2
            {
              ///////////////////////////////////////////
              // ALGORITHM 2: ESTIMATION OF LAMBDA MAX //
              ///////////////////////////////////////////
              // Requires: hStarL, hStarR
              // Ensures: lambdaMax
              int aux_counter = 0;
              while (true)
                {
                  aux_counter++;
                  // Start having lambdaMin and lambdaMax
                  // Check if current lambdaMin and lambdaMax satisfy the tolerance
                  if (verbose)
                    {
                      std::cout << lambdaMin << ", " << lambdaMax << std::endl;
                    }
                  if (lambdaMin > 0 && lambdaMax/lambdaMin - 1 <= tol)
                    return lambdaMax;
                  // Check for round off error
                  if (phi(g,hStarL,hL,hR,velL,velR) > 0 || phi(g,hStarR,hL,hR,velL,velR) < 0)
                    return lambdaMax;

                  // Save old estimates of hStar
                  double hStarL_old = hStarL;
                  double hStarR_old = hStarR;
                  // Compute new estimates on hStarL and hStarR
                  // NOTE (MQL): hStarL and hStarR must be computed using the old values
                  hStarL = hStarLFromQuadPhiFromAbove(g,hStarL_old,hStarR_old,hL,hR,velL,velR);
                  hStarR = hStarRFromQuadPhiFromBelow(g,hStarL_old,hStarR_old,hL,hR,velL,velR);

                  // Compute lambdaMax and lambdaMin
                  nu11 = nu1(g,hStarR,hL,velL);
                  nu12 = nu1(g,hStarL,hL,velL);
                  nu31 = nu3(g,hStarL,hR,velR);
                  nu32 = nu3(g,hStarR,hR,velR);

                  lambdaMin = fmax(fmax(nu31,0), fmax(-nu12,0));
                  lambdaMax = fmax(fmax(nu32,0), fmax(-nu11,0));

                  if (aux_counter>1000) //TMP
                    {
                      std::cout << "**** AUX COUNTER > 1000... aborting!" << std::endl;
                      std::cout << "**** Initial guess hStar: " << hStar0 << std::endl;

                      hStar = hStar0;
                      lambda1 = nu1(g,hStar,hL,velL);
                      lambda3 = nu3(g,hStar,hR,velR);
                      std::cout << "**** Initial estimate of max wave speed: "
                                << fmax(fabs(lambda1),fabs(lambda3)) << std::endl;

                      abort();
                    }
                  //else
                  //{
                  //  std::cout << "*****... AUX COUNTER: " << aux_counter << std::endl; //TMP
                  //}
                }
            }
        }
    }

    inline
      void calculateCFL(const double& elementDiameter,
                        const double& g,
                        const double& h,
                        const double& hu,
                        const double& hv,
                        const double hEps,
                        double& cfl)
    {
      double cflx, cfly, c=sqrt(fmax(g*hEps,g*h));
      double u = 2*h/(h*h+std::pow(fmax(h,hEps),2))*hu;
      double v = 2*h/(h*h+std::pow(fmax(h,hEps),2))*hv;

      if (u > 0.0)
        cflx = (u+c)/elementDiameter;
      else
        cflx = fabs(u-c)/elementDiameter;

      if (v > 0.0)
        cfly = (v+c)/elementDiameter;
      else
        cfly = fabs(v-c)/elementDiameter;
      cfl = sqrt(cflx*cflx+cfly*cfly);//hack, conservative estimate
    }

    void FCTStep(double dt,
		 int NNZ, //number on non-zero entries on sparsity pattern
		 int numDOFs, //number of DOFs
		 double* lumped_mass_matrix, //lumped mass matrix (as vector))
		 double* h_old, //DOFs of solution at last stage
		 double* hu_old,
		 double* hv_old,
		 double* b_dof,
		 double* high_order_hnp1, //DOFs of high order solution at tnp1
		 double* high_order_hunp1,
		 double* high_order_hvnp1,
		 double* extendedSourceTerm_hu,
		 double* extendedSourceTerm_hv,
		 double* limited_hnp1,
		 double* limited_hunp1,
		 double* limited_hvnp1,
		 int* csrRowIndeces_DofLoops, //csr row indeces
		 int* csrColumnOffsets_DofLoops, //csr column offsets
		 double* MassMatrix, //mass matrix
		 double* dH_minus_dL,
		 double* muH_minus_muL,
		 double hEps,
		 double* hReg,
		 int LUMPED_MASS_MATRIX,
		 double* dLow,
		 double* hBT,
		 double* huBT,
		 double* hvBT)
    {
      Rneg.resize(numDOFs,0.0);
      Rpos.resize(numDOFs,0.0);
      hLow.resize(numDOFs,0.0);
      huLow.resize(numDOFs,0.0);
      hvLow.resize(numDOFs,0.0);

      ////////////////////////
      // FIRST LOOP in DOFs //
      ////////////////////////
      int ij=0;
      for (int i=0; i<numDOFs; i++)
        {
          //read some vectors
          double high_order_hnp1i  = high_order_hnp1[i];
          double hi = h_old[i];
          double hui = hu_old[i];
          double hvi = hv_old[i];
          double Zi = b_dof[i];
          double mi = lumped_mass_matrix[i];

	  double hiMin = hi;
	  double hiMax = hi;
          double Pnegi=0., Pposi=0.;

	  // LOW ORDER SOLUTION without extended source term. Eqn 6.23
	  hLow[i]  = hi;
	  huLow[i] = hui;
	  hvLow[i] = hvi;

          // LOOP OVER THE SPARSITY PATTERN (j-LOOP)//
          for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
            {
              int j = csrColumnOffsets_DofLoops[offset];
              // read some vectors
              double hj = h_old[j];
              double Zj = b_dof[j];

              // COMPUTE STAR SOLUTION // hStar, huStar and hvStar
              double hStarij  = fmax(0., hi + Zi - fmax(Zi,Zj));
              double hStarji  = fmax(0., hj + Zj - fmax(Zi,Zj));

              // i-th row of flux correction matrix
              double ML_minus_MC = (LUMPED_MASS_MATRIX == 1 ? 0. :
				    (i==j ? 1. : 0.)*mi - MassMatrix[ij]);
              double FluxCorrectionMatrix1 =
                ML_minus_MC*(high_order_hnp1[j]-hj - (high_order_hnp1i-hi))
                + dt*(dH_minus_dL[ij]-muH_minus_muL[ij])*(hStarji-hStarij)
                + dt*muH_minus_muL[ij]*(hj-hi);

              // COMPUTE P VECTORS //
              Pnegi += FluxCorrectionMatrix1*((FluxCorrectionMatrix1 < 0) ? 1. : 0.);
	      Pposi += FluxCorrectionMatrix1*((FluxCorrectionMatrix1 > 0) ? 1. : 0.);

	      // COMPUTE LOCAL BOUNDS //
	      hiMin = std::min(hiMin, hBT[ij]);
	      hiMax = std::max(hiMax, hBT[ij]);

	      // COMPUTE LOW ORDER SOLUTION (WITHOUT EXTENDED SOURCE) //
	      if (i!=j)
		{
		  hLow[i]  += hi*(-dt/mi*2*dLow[ij])  + dt/mi*(2*dLow[ij]*hBT[ij]);
		  huLow[i] += hui*(-dt/mi*2*dLow[ij]) + dt/mi*(2*dLow[ij]*huBT[ij]);
		  hvLow[i] += hvi*(-dt/mi*2*dLow[ij]) + dt/mi*(2*dLow[ij]*hvBT[ij]);
		}
              //update ij
              ij+=1;
            }
	  // clean up hLow from round off error
	  if (hLow[i] < hEps)
	    hLow[i]=0;
          ///////////////////////
          // COMPUTE Q VECTORS //
          ///////////////////////
	  if (GLOBAL_FCT==1)
	    hiMin=0;
          double Qnegi = mi*(hiMin-hLow[i]);
	  double Qposi = mi*(hiMax-hLow[i]);

          ///////////////////////
          // COMPUTE R VECTORS //
          ///////////////////////
          if (high_order_hnp1[i] <= hReg[i]) //hEps
	    {
	      Rneg[i] = 0.;
	      Rpos[i] = 0.;
	    }
          else
	    {
	      Rneg[i] = ((Pnegi==0) ? 1. : std::min(1.0,Qnegi/Pnegi));
	      Rpos[i] = ((Pposi==0) ? 1. : std::min(1.0,Qposi/Pposi));
	    }
        } // i DOFs

      //////////////////////
      // COMPUTE LIMITERS //
      //////////////////////
      ij=0;
      for (int i=0; i<numDOFs; i++)
        {
          //read some vectors
          double high_order_hnp1i  = high_order_hnp1[i];
          double high_order_hunp1i = high_order_hunp1[i];
          double high_order_hvnp1i = high_order_hvnp1[i];
          double hi = h_old[i];
          double huni = hu_old[i];
          double hvni = hv_old[i];
          double Zi = b_dof[i];
          double mi = lumped_mass_matrix[i];
          double one_over_hiReg = 2*hi/(hi*hi+std::pow(fmax(hi,hEps),2)); //hEps

          double ith_Limiter_times_FluxCorrectionMatrix1 = 0.;
          double ith_Limiter_times_FluxCorrectionMatrix2 = 0.;
          double ith_Limiter_times_FluxCorrectionMatrix3 = 0.;

          // LOOP OVER THE SPARSITY PATTERN (j-LOOP)//
          for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
            {
              int j = csrColumnOffsets_DofLoops[offset];
              // read some vectors
              double hj = h_old[j];
              double hunj = hu_old[j];
              double hvnj = hv_old[j];
              double Zj = b_dof[j];
              double one_over_hjReg = 2*hj/(hj*hj+std::pow(fmax(hj,hEps),2)); //hEps

              // COMPUTE STAR SOLUTION // hStar, huStar and hvStar
              double hStarij  = fmax(0., hi + Zi - fmax(Zi,Zj));
              double huStarij = huni*hStarij*one_over_hiReg;
              double hvStarij = hvni*hStarij*one_over_hiReg;

              double hStarji  = fmax(0., hj + Zj - fmax(Zi,Zj));
              double huStarji = hunj*hStarji*one_over_hjReg;
              double hvStarji = hvnj*hStarji*one_over_hjReg;

              // COMPUTE FLUX CORRECTION MATRICES
              double ML_minus_MC = (LUMPED_MASS_MATRIX == 1 ? 0. :
				    (i==j ? 1. : 0.)*mi - MassMatrix[ij]);
              double FluxCorrectionMatrix1 =
                ML_minus_MC*(high_order_hnp1[j]-hj - (high_order_hnp1i-hi))
                + dt*(dH_minus_dL[ij]-muH_minus_muL[ij])*(hStarji-hStarij)
                + dt*muH_minus_muL[ij]*(hj-hi);

              double FluxCorrectionMatrix2 =
                ML_minus_MC*(high_order_hunp1[j]-hunj - (high_order_hunp1i-huni))
                + dt*(dH_minus_dL[ij]-muH_minus_muL[ij])*(huStarji-huStarij)
                + dt*muH_minus_muL[ij]*(hunj-huni);

              double FluxCorrectionMatrix3 =
                ML_minus_MC*(high_order_hvnp1[j]-hvnj - (high_order_hvnp1i-hvni))
                + dt*(dH_minus_dL[ij]-muH_minus_muL[ij])*(hvStarji-hvStarij)
                + dt*muH_minus_muL[ij]*(hvnj-hvni);

              // COMPUTE LIMITER // based on water height
	      double Lij = 0.;
	      if (FluxCorrectionMatrix1 >= 0)
		Lij = std::min(Rpos[i],Rneg[j]);
	      else
		Lij = std::min(Rneg[i],Rpos[j]);
	      if (GLOBAL_FCT==1)
		Lij = (FluxCorrectionMatrix1 >= 0. ? std::min(1.,Rneg[j]) : std::min(Rneg[i],1.));

	      // COMPUTE LIMITED FLUX //
              ith_Limiter_times_FluxCorrectionMatrix1 += Lij*FluxCorrectionMatrix1;
              ith_Limiter_times_FluxCorrectionMatrix2 += Lij*FluxCorrectionMatrix2;
              ith_Limiter_times_FluxCorrectionMatrix3 += Lij*FluxCorrectionMatrix3;

              //update ij
              ij+=1;
            }
          double one_over_mi = 1.0/lumped_mass_matrix[i];
	  limited_hnp1[i]  = hLow[i] + one_over_mi*ith_Limiter_times_FluxCorrectionMatrix1;
          limited_hunp1[i] = ((huLow[i] - dt/mi*extendedSourceTerm_hu[i]) // low_order_hunp1+...
	  		      +one_over_mi*ith_Limiter_times_FluxCorrectionMatrix2);
          limited_hvnp1[i] = ((hvLow[i] - dt/mi*extendedSourceTerm_hv[i]) // low_order_hvnp1+...
			      +one_over_mi*ith_Limiter_times_FluxCorrectionMatrix3);

          if (limited_hnp1[i] < -hEps && dt < 1.0)
            {
              std::cout << "Limited water height is negative: "
			<< "hLow: " << hLow[i] << "\t"
                        << "hHigh: " << limited_hnp1[i] << "\t"
                        << " ... aborting!" << std::endl;
              abort();
            }
          else
            {
	      // clean up uHigh from round off error
	      if (limited_hnp1[i] < hEps)
		limited_hnp1[i] = 0;
              //double aux = fmax(limited_hnp1[i],hEps); // hEps
              double aux = fmax(limited_hnp1[i],hReg[i]); // hReg makes the code more robust
              limited_hunp1[i] *=
		2*std::pow(limited_hnp1[i],VEL_FIX_POWER)
		/(std::pow(limited_hnp1[i],VEL_FIX_POWER)+std::pow(aux,VEL_FIX_POWER));
              limited_hvnp1[i] *=
		2*std::pow(limited_hnp1[i],VEL_FIX_POWER)
		/(std::pow(limited_hnp1[i],VEL_FIX_POWER)+std::pow(aux,VEL_FIX_POWER));
            }
        }
    }

    void convexLimiting(double dt,
			int NNZ, //number on non-zero entries on sparsity pattern
			int numDOFs, //number of DOFs
			double* lumped_mass_matrix, //lumped mass matrix (as vector))
			double* h_old, //DOFs of solution at last stage
			double* hu_old,
			double* hv_old,
			double* b_dof,
			double* high_order_hnp1, //DOFs of high order solution at tnp1
			double* high_order_hunp1,
			double* high_order_hvnp1,
			double* extendedSourceTerm_hu,
			double* extendedSourceTerm_hv,
			double* limited_hnp1,
			double* limited_hunp1,
			double* limited_hvnp1,
			int* csrRowIndeces_DofLoops, //csr row indeces
			int* csrColumnOffsets_DofLoops, //csr column offsets
			double* MassMatrix, //mass matrix
			double* dH_minus_dL,
			double* muH_minus_muL,
			double hEps,
			double* hReg,
			int LUMPED_MASS_MATRIX,
			double* dLow,
			double* hBT,
			double* huBT,
			double* hvBT)
    {
      Rneg.resize(numDOFs,0.0);
      Rpos.resize(numDOFs,0.0);
      hLow.resize(numDOFs,0.0);
      huLow.resize(numDOFs,0.0);
      hvLow.resize(numDOFs,0.0);
      Kmax.resize(numDOFs,0.0);

      ////////////////////////
      // FIRST LOOP in DOFs //
      ////////////////////////
      int ij=0;
      for (int i=0; i<numDOFs; i++)
        {
          //read some vectors
          double high_order_hnp1i  = high_order_hnp1[i];
          double hi = h_old[i];
	  double hui = hu_old[i];
	  double hvi = hv_old[i];
          double Zi = b_dof[i];
          double mi = lumped_mass_matrix[i];

	  double hiMin = hi;
	  double hiMax = hi;
          double Pnegi=0., Pposi=0.;

	  // LOW ORDER SOLUTION without extended source term. Eqn 6.23
	  hLow[i]  = hi;
	  huLow[i] = hui;
	  hvLow[i] = hvi;
	  Kmax[i] = 0;

          // LOOP OVER THE SPARSITY PATTERN (j-LOOP)//
          for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
            {
              int j = csrColumnOffsets_DofLoops[offset];
	      double psi_ij = 0;
	      if (hBT != 0)
		psi_ij = 1/(2*hBT[ij])*(huBT[ij]*huBT[ij] + huBT[ij]*huBT[ij]); // Eqn (6.29)
	      Kmax[i] = fmax(psi_ij,Kmax[i]);

              // read some vectors
              double hj = h_old[j];
              double Zj = b_dof[j];

              // COMPUTE STAR SOLUTION // hStar, huStar and hvStar
              double hStarij  = fmax(0., hi + Zi - fmax(Zi,Zj));
              double hStarji  = fmax(0., hj + Zj - fmax(Zi,Zj));

              // i-th row of flux correction matrix
              double ML_minus_MC = (LUMPED_MASS_MATRIX == 1 ? 0. :
				    (i==j ? 1. : 0.)*mi - MassMatrix[ij]);
              double FluxCorrectionMatrix1 =
                ML_minus_MC*(high_order_hnp1[j]-hj - (high_order_hnp1i-hi))
                + dt*(dH_minus_dL[ij]-muH_minus_muL[ij])*(hStarji-hStarij)
                + dt*muH_minus_muL[ij]*(hj-hi);

              // COMPUTE P VECTORS //
              Pnegi += FluxCorrectionMatrix1*((FluxCorrectionMatrix1 < 0) ? 1. : 0.);
	      Pposi += FluxCorrectionMatrix1*((FluxCorrectionMatrix1 > 0) ? 1. : 0.);

	      // COMPUTE LOCAL BOUNDS //
	      hiMin = std::min(hiMin, hBT[ij]);
	      hiMax = std::max(hiMax, hBT[ij]);

	      // COMPUTE LOW ORDER SOLUTION (WITHOUT EXTENDED SOURCE) //
	      if (i!=j)
		{
		  hLow[i]  += hi*(-dt/mi*2*dLow[ij])  + dt/mi*(2*dLow[ij]*hBT[ij]);
		  huLow[i] += hui*(-dt/mi*2*dLow[ij]) + dt/mi*(2*dLow[ij]*huBT[ij]);
		  hvLow[i] += hvi*(-dt/mi*2*dLow[ij]) + dt/mi*(2*dLow[ij]*hvBT[ij]);
		}
              // UPDATE ij //
              ij+=1;
            }
	  // clean up hLow from round off error
	  if (hLow[i] < hEps)
	    hLow[i]=0;
          ///////////////////////
          // COMPUTE Q VECTORS //
          ///////////////////////
          double Qnegi = mi*(hiMin-hLow[i]);
	  double Qposi = mi*(hiMax-hLow[i]);

          ///////////////////////
          // COMPUTE R VECTORS //
          ///////////////////////
          if (high_order_hnp1[i] <= hReg[i]) //hEps
	    {
	      Rneg[i] = 0.;
	      Rpos[i] = 0.;
	    }
          else
	    {
	      Rneg[i] = ((Pnegi==0) ? 1. : std::min(1.0,Qnegi/Pnegi));
	      Rpos[i] = ((Pposi==0) ? 1. : std::min(1.0,Qposi/Pposi));
	    }
        } // i DOFs

      //////////////////////
      // COMPUTE LIMITERS //
      //////////////////////
      ij=0;
      for (int i=0; i<numDOFs; i++)
        {
          //read some vectors
          double high_order_hnp1i  = high_order_hnp1[i];
          double high_order_hunp1i = high_order_hunp1[i];
          double high_order_hvnp1i = high_order_hvnp1[i];
          double hi = h_old[i];
          double huni = hu_old[i];
          double hvni = hv_old[i];
          double Zi = b_dof[i];
          double mi = lumped_mass_matrix[i];
          double one_over_hiReg = 2*hi/(hi*hi+std::pow(fmax(hi,hEps),2)); //hEps

          double ith_Limiter_times_FluxCorrectionMatrix1 = 0.;
          double ith_Limiter_times_FluxCorrectionMatrix2 = 0.;
          double ith_Limiter_times_FluxCorrectionMatrix3 = 0.;

	  double ci = Kmax[i]*hLow[i]-0.5*(huLow[i]*huLow[i]+hvLow[i]*hvLow[i]); // for conv. lim.

          // LOOP OVER THE SPARSITY PATTERN (j-LOOP)//
          for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
            {
              int j = csrColumnOffsets_DofLoops[offset];
              // read some vectors
              double hj = h_old[j];
              double hunj = hu_old[j];
              double hvnj = hv_old[j];
              double Zj = b_dof[j];
              double one_over_hjReg = 2*hj/(hj*hj+std::pow(fmax(hj,hEps),2)); //hEps

              // COMPUTE STAR SOLUTION // hStar, huStar and hvStar
              double hStarij  = fmax(0., hi + Zi - fmax(Zi,Zj));
              double huStarij = huni*hStarij*one_over_hiReg;
              double hvStarij = hvni*hStarij*one_over_hiReg;

              double hStarji  = fmax(0., hj + Zj - fmax(Zi,Zj));
              double huStarji = hunj*hStarji*one_over_hjReg;
              double hvStarji = hvnj*hStarji*one_over_hjReg;

              // COMPUTE FLUX CORRECTION MATRICES
              double ML_minus_MC = (LUMPED_MASS_MATRIX == 1 ? 0. :
				    (i==j ? 1. : 0.)*mi - MassMatrix[ij]);
              double FluxCorrectionMatrix1 =
                ML_minus_MC*(high_order_hnp1[j]-hj - (high_order_hnp1i-hi))
                + dt*(dH_minus_dL[ij]-muH_minus_muL[ij])*(hStarji-hStarij)
                + dt*muH_minus_muL[ij]*(hj-hi);

              double FluxCorrectionMatrix2 =
                ML_minus_MC*(high_order_hunp1[j]-hunj - (high_order_hunp1i-huni))
                + dt*(dH_minus_dL[ij]-muH_minus_muL[ij])*(huStarji-huStarij)
                + dt*muH_minus_muL[ij]*(hunj-huni);

              double FluxCorrectionMatrix3 =
                ML_minus_MC*(high_order_hvnp1[j]-hvnj - (high_order_hvnp1i-hvni))
                + dt*(dH_minus_dL[ij]-muH_minus_muL[ij])*(hvStarji-hvStarij)
                + dt*muH_minus_muL[ij]*(hvnj-hvni);

	      // compute limiter based on water height
	      double Lij = 0.;
	      if (FluxCorrectionMatrix1 >= 0)
		Lij = std::min(Rpos[i],Rneg[j]);
	      else
		Lij = std::min(Rneg[i],Rpos[j]);

	      // CONVEX LIMITING // for kinetic energy
	      // root of ith-DOF
	      double lambdaj = csrRowIndeces_DofLoops[i+1]-csrRowIndeces_DofLoops[i]-1;
	      double Ph_ij = FluxCorrectionMatrix1/mi/lambdaj;
	      double Phu_ij = FluxCorrectionMatrix2/mi/lambdaj;
	      double Phv_ij = FluxCorrectionMatrix3/mi/lambdaj;

	      double ai = -0.5*(Phu_ij*Phu_ij+Phv_ij*Phv_ij);
	      double bi = Kmax[i]*Ph_ij-(huLow[i]*Phu_ij + hvLow[i]*Phv_ij);

	      double r1 = ai==0 ? (bi==0 ? 1. : -ci/bi) : (-bi + std::sqrt(bi*bi-4*ai*ci))/2./ai;
	      double r2 = ai==0 ? (bi==0 ? 1. : -ci/bi) : (-bi - std::sqrt(bi*bi-4*ai*ci))/2./ai;
	      if (r1<0 && r2<0)
		{
		  r1=1.;
		  r2=1.;
		}
	      double ri = fabs(fmax(r1,r2));

	      // root of jth-DOF (To compute transpose component)
	      double lambdai = csrRowIndeces_DofLoops[j+1]-csrRowIndeces_DofLoops[j]-1;
	      double mj = lumped_mass_matrix[j];
	      double cj = Kmax[j]*hLow[j]-0.5*(huLow[j]*huLow[j]+hvLow[j]*hvLow[j]);
	      double Ph_ji  = -FluxCorrectionMatrix1/mj/lambdai; //Aij=-Aji
	      double Phu_ji = -FluxCorrectionMatrix2/mj/lambdai;
	      double Phv_ji = -FluxCorrectionMatrix3/mj/lambdai;
	      double aj = -0.5*(Phu_ji*Phu_ji+Phv_ji*Phv_ji);
	      double bj = Kmax[j]*Ph_ji-(huLow[j]*Phu_ji + hvLow[j]*Phv_ji);

	      r1 = aj==0 ? (bj==0 ? 1. : -cj/bj) : (-bj + std::sqrt(bj*bj-4*aj*cj))/2./aj;
	      r2 = aj==0 ? (bj==0 ? 1. : -cj/bj) : (-bj - std::sqrt(bj*bj-4*aj*cj))/2./aj;
	      if (r1<0 && r2<0)
		{
		  r1=1.;
		  r2=1.;
		}
	      double rj = fabs(fmax(r1,r2));

	      // COMPUTE LIMITER //
	      Lij = fmin(fmin(ri,Lij),fmin(rj,Lij)); //Lij=Lji

	      // COMPUTE LIMITED FLUX //
              ith_Limiter_times_FluxCorrectionMatrix1 += Lij*FluxCorrectionMatrix1;
              ith_Limiter_times_FluxCorrectionMatrix2 += Lij*FluxCorrectionMatrix2;
              ith_Limiter_times_FluxCorrectionMatrix3 += Lij*FluxCorrectionMatrix3;

              //update ij
              ij+=1;
            }
          double one_over_mi = 1.0/lumped_mass_matrix[i];
          limited_hnp1[i]  = hLow[i] + one_over_mi*ith_Limiter_times_FluxCorrectionMatrix1;
          limited_hunp1[i] = ((huLow[i] - dt/mi*extendedSourceTerm_hu[i]) // low_order_hunp1+...
			      +one_over_mi*ith_Limiter_times_FluxCorrectionMatrix2);
          limited_hvnp1[i] = ((hvLow[i] - dt/mi*extendedSourceTerm_hv[i]) // low_order_hvnp1+...
			      +one_over_mi*ith_Limiter_times_FluxCorrectionMatrix3);

          if (limited_hnp1[i] < -hEps && dt < 1.0)
	    {
              std::cout << "Limited water height is negative: "
			<< "hLow: " << hLow[i] << "\t"
                        << "hHigh: " << limited_hnp1[i] << "\t"
                        << " ... aborting!" << std::endl;
              abort();
	    }
          else
	    {
	      // clean up uHigh from round off error
	      if (limited_hnp1[i] < hEps)
		limited_hnp1[i] = 0;
              //double aux = fmax(limited_hnp1[i],hEps); // hEps
              double aux = fmax(limited_hnp1[i],hReg[i]); // hReg makes the code more robust
              limited_hunp1[i] *=
		2*std::pow(limited_hnp1[i],VEL_FIX_POWER)
		/(std::pow(limited_hnp1[i],VEL_FIX_POWER)+std::pow(aux,VEL_FIX_POWER));
              limited_hvnp1[i] *=
		2*std::pow(limited_hnp1[i],VEL_FIX_POWER)
		/(std::pow(limited_hnp1[i],VEL_FIX_POWER)+std::pow(aux,VEL_FIX_POWER));
	    }
        }
    }

    double calculateEdgeBasedCFL(double g,
                                 int numDOFsPerEqn, //number of DOFs
                                 double* lumped_mass_matrix, //lumped mass matrix (as vector))
                                 double* h_dof_old, //DOFs of solution at last stage
                                 double* hu_dof_old,
                                 double* hv_dof_old,
                                 double* b_dof,
                                 int* csrRowIndeces_DofLoops, //csr row indeces
                                 int* csrColumnOffsets_DofLoops, //csr column offsets
                                 double hEps,
                                 double* hReg,
                                 double* Cx,
                                 double* Cy,
                                 double* CTx,
                                 double* CTy,
                                 double* dLow,
                                 double run_cfl,
                                 double* edge_based_cfl,
				 int debug)
    {
      std::valarray<double> psi(numDOFsPerEqn);
      double max_edge_based_cfl = 0.;
      int ij=0;
      for (int i=0; i<numDOFsPerEqn; i++)
        {
          double hi = h_dof_old[i]; // solution at time tn for the ith DOF
          double hui = hu_dof_old[i];
          double hvi = hv_dof_old[i];
          double dLowii = 0.;

          double alphai;
          double alpha_numerator = 0.;
          double alpha_denominator = 0.;
          for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
            { //loop in j (sparsity pattern)
              int j = csrColumnOffsets_DofLoops[offset];
              double hj = h_dof_old[j]; // solution at time tn for the jth DOF
              double huj = hu_dof_old[j];
              double hvj = hv_dof_old[j];

              if (i != j)
                {
                  ////////////////////////
                  // DISSIPATIVE MATRIX //
                  ////////////////////////
                  double cij_norm = sqrt(Cx[ij]*Cx[ij] + Cy[ij]*Cy[ij]);
                  double cji_norm = sqrt(CTx[ij]*CTx[ij] + CTy[ij]*CTy[ij]);
                  double nxij = Cx[ij]/cij_norm, nyij = Cy[ij]/cij_norm;
                  double nxji = CTx[ij]/cji_norm, nyji = CTy[ij]/cji_norm;
                  dLow[ij] = fmax(maxWaveSpeedSharpInitialGuess(g,nxij,nyij,
                                                                hi,hui,hvi,
                                                                hj,huj,hvj,
                                                                hEps,hEps,debug)*cij_norm, //hEps
                                  maxWaveSpeedSharpInitialGuess(g,nxji,nyji,
                                                                hj,huj,hvj,
                                                                hi,hui,hvi,
                                                                hEps,hEps,debug)*cji_norm); //hEps
                  dLowii -= dLow[ij];

                  // FOR SMOOTHNESS INDICATOR //
                  alpha_numerator += hj - hi;
                  alpha_denominator += fabs(hj - hi);
                }
              else
                dLow[ij] = 0.;
              //update ij
              ij+=1;
            }
          //////////////////////////////
          // CALCULATE EDGE BASED CFL //
          //////////////////////////////
          double mi = lumped_mass_matrix[i];
          edge_based_cfl[i] = 2*fabs(dLowii)/mi;
          max_edge_based_cfl = fmax(max_edge_based_cfl,edge_based_cfl[i]);

          //////////////////////////////////
          // COMPUTE SMOOTHNESS INDICATOR //
          //////////////////////////////////
          if (hi <= hReg[i]) //hEps, hReg makes the method more robust
            alphai = 1.;
          else
            {
              if (fabs(alpha_numerator) <= hEps) //hEps. Force alphai=0 in constant states. This for well balancing wrt friction
                alphai = 0.;
              else
                alphai = fabs(alpha_numerator)/(alpha_denominator+1E-15);
            }
          if (POWER_SMOOTHNESS_INDICATOR==0)
            psi[i] = 1.0;
          else
            psi[i] = std::pow(alphai,POWER_SMOOTHNESS_INDICATOR); //NOTE: alpha^2 in the paper
        }

      if (REESTIMATE_MAX_EDGE_BASED_CFL==1)
        {
          // CALCULATE FIRST GUESS dt //
          double dt = run_cfl/max_edge_based_cfl;
          ij=0;
          for (int i=0; i<numDOFsPerEqn; i++)
            {
              double hi = h_dof_old[i]; // solution at time tn for the ith DOF
              double hui = hu_dof_old[i];
              double hvi = hv_dof_old[i];
              double Zi = b_dof[i];
              double one_over_hiReg = 2*hi/(hi*hi+std::pow(fmax(hi,hEps),2)); // hEps
              double ui = hui*one_over_hiReg;
              double vi = hvi*one_over_hiReg;
              // flux and stabilization variables to compute low order solution
              double ith_flux_term1 = 0.;
              double ith_dLij_minus_muLij_times_hStarStates = 0.;
              double ith_muLij_times_hStates = 0.;

              for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
                {
                  int j = csrColumnOffsets_DofLoops[offset];
                  double hj = h_dof_old[j]; // solution at time tn for the jth DOF
                  double huj = hu_dof_old[j];
                  double hvj = hv_dof_old[j];
                  double Zj = b_dof[j];
                  double one_over_hjReg = 2*hj/(hj*hj+std::pow(fmax(hj,hEps),2)); //hEps
                  double uj = huj*one_over_hjReg;
                  double vj = hvj*one_over_hjReg;

                  // Star states for water height
                  double dLij, muLij, muLowij;
                  double hStarij  = fmax(0., hi + Zi - fmax(Zi,Zj));
                  double hStarji  = fmax(0., hj + Zj - fmax(Zi,Zj));

                  // compute flux term
                  ith_flux_term1 += huj*Cx[ij] + hvj*Cy[ij]; // f1*C
                  if (i != j)
                    {
                      dLij = dLow[ij];//*fmax(psi[i],psi[j]); // enhance the order to 2nd order. No EV

                      muLowij = fmax(fmax(0.,-(ui*Cx[ij] + vi*Cy[ij])),fmax(0,(uj*Cx[ij] + vj*Cy[ij])));
                      muLij = muLowij;//*fmax(psi[i],psi[j]); // enhance the order to 2nd order. No EV
                      // compute dissipative terms
                      ith_dLij_minus_muLij_times_hStarStates  += (dLij - muLij)*(hStarji-hStarij);
                      ith_muLij_times_hStates  += muLij*(hj-hi);
                    }
                  // update ij
                  ij += 1;
                }

              double mi = lumped_mass_matrix[i];
              double low_order_hnp1 = hi - dt/mi*(ith_flux_term1
                                                  - ith_dLij_minus_muLij_times_hStarStates
                                                  - ith_muLij_times_hStates);
              while (low_order_hnp1 < -1E-14 && dt < 1.0)
                { // Water height is negative. Recalculate dt
                  std::cout << "********** ... Reducing dt from original estimate to achieve positivity... **********" << std::endl;
                  dt /= 2.;
                  low_order_hnp1 = hi - dt/mi*(ith_flux_term1
                                               - ith_dLij_minus_muLij_times_hStarStates
                                               - ith_muLij_times_hStates);
                }
            }
          // new max_edge_based_cfl
          max_edge_based_cfl = run_cfl/dt;
        }
      return max_edge_based_cfl;
    }

    void calculateResidual(// last EDGE BASED version
                           double* mesh_trial_ref,
                           double* mesh_grad_trial_ref,
                           double* mesh_dof,
                           double* mesh_velocity_dof,
                           double MOVING_DOMAIN,
                           int* mesh_l2g,
                           double* dV_ref,
                           double* h_trial_ref,
                           double* h_grad_trial_ref,
                           double* h_test_ref,
                           double* h_grad_test_ref,
                           double* vel_trial_ref,
                           double* vel_grad_trial_ref,
                           double* vel_test_ref,
                           double* vel_grad_test_ref,
                           //element boundary
                           double* mesh_trial_trace_ref,
                           double* mesh_grad_trial_trace_ref,
                           double* dS_ref,
                           double* h_trial_trace_ref,
                           double* h_grad_trial_trace_ref,
                           double* h_test_trace_ref,
                           double* h_grad_test_trace_ref,
                           double* vel_trial_trace_ref,
                           double* vel_grad_trial_trace_ref,
                           double* vel_test_trace_ref,
                           double* vel_grad_test_trace_ref,
                           double* normal_ref,
                           double* boundaryJac_ref,
                           //physics
                           double* elementDiameter,
                           int nElements_global,
                           double useRBLES,
                           double useMetrics,
                           double alphaBDF,
                           double nu,
                           double g,
                           int* h_l2g,
                           int* vel_l2g,
                           double* h_dof_old,
                           double* hu_dof_old,
                           double* hv_dof_old,
                           double* b_dof,
                           double* h_dof,
                           double* hu_dof,
                           double* hv_dof,
                           double* h_dof_sge,
                           double* hu_dof_sge,
                           double* hv_dof_sge,
                           double* q_mass_acc,
                           double* q_mom_hu_acc,
                           double* q_mom_hv_acc,
                           double* q_mass_adv,
                           double* q_mass_acc_beta_bdf,
                           double* q_mom_hu_acc_beta_bdf,
                           double* q_mom_hv_acc_beta_bdf,
                           double* q_cfl,
                           int* sdInfo_hu_hu_rowptr,
                           int* sdInfo_hu_hu_colind,
                           int* sdInfo_hu_hv_rowptr,
                           int* sdInfo_hu_hv_colind,
                           int* sdInfo_hv_hv_rowptr,
                           int* sdInfo_hv_hv_colind,
                           int* sdInfo_hv_hu_rowptr,
                           int* sdInfo_hv_hu_colind,
                           int offset_h,
                           int offset_hu,
                           int offset_hv,
                           int stride_h,
                           int stride_hu,
                           int stride_hv,
                           double* globalResidual,
                           int nExteriorElementBoundaries_global,
                           int* exteriorElementBoundariesArray,
                           int* elementBoundaryElementsArray,
                           int* elementBoundaryLocalElementBoundariesArray,
                           int* isDOFBoundary_h,
                           int* isDOFBoundary_hu,
                           int* isDOFBoundary_hv,
                           int* isAdvectiveFluxBoundary_h,
                           int* isAdvectiveFluxBoundary_hu,
                           int* isAdvectiveFluxBoundary_hv,
                           int* isDiffusiveFluxBoundary_hu,
                           int* isDiffusiveFluxBoundary_hv,
                           double* ebqe_bc_h_ext,
                           double* ebqe_bc_flux_mass_ext,
                           double* ebqe_bc_flux_mom_hu_adv_ext,
                           double* ebqe_bc_flux_mom_hv_adv_ext,
                           double* ebqe_bc_hu_ext,
                           double* ebqe_bc_flux_hu_diff_ext,
                           double* ebqe_penalty_ext,
                           double* ebqe_bc_hv_ext,
                           double* ebqe_bc_flux_hv_diff_ext,
                           double* q_velocity,
                           double* ebqe_velocity,
                           double* flux,
                           double* elementResidual_h_save,
                           // C matrices
                           double* Cx,
                           double* Cy,
                           double* CTx,
                           double* CTy,
                           // PARAMETERS FOR EDGE BASED STABILIZATION
                           int numDOFsPerEqn,
                           int NNZ,
                           int* csrRowIndeces_DofLoops,
                           int* csrColumnOffsets_DofLoops,
                           // LUMPED MASS MATRIX
                           double* lumped_mass_matrix,
                           double cfl_run,
                           double hEps,
                           double* hReg,
                           // SAVE SOLUTION (mql)
                           double* hnp1_at_quad_point,
                           double* hunp1_at_quad_point,
                           double* hvnp1_at_quad_point,
                           // TO COMPUTE LOW ORDER
			   double* extendedSourceTerm_hu,
			   double* extendedSourceTerm_hv,
                           // FOR FCT
                           double* dH_minus_dL,
                           double* muH_minus_muL,
                           double cE,
                           int LUMPED_MASS_MATRIX,
                           double dt,
                           int LINEAR_FRICTION,
                           double mannings,
                           // Quant of interests
                           double* quantDOFs,
                           int SECOND_CALL_CALCULATE_RESIDUAL,
                           // NORMAL COMPONENTS
                           int COMPUTE_NORMALS,
                           double* normalx,
                           double* normaly,
                           double* dLow,
			   // LOCAL LIMITING
			   double* hBT,
			   double* huBT,
			   double* hvBT,
                           int lstage)
    {
      //FOR FRICTION//
      double n2 = std::pow(mannings,2.);
      double gamma=4./3;
      // mql. This parameter relaxes the cfl restriction.
      // It is now conservative. I might change it after the local bounds are implemented
      double xi=10.;

      //////////////////////////////////////
      // ********** CELL LOOPS ********** //
      //////////////////////////////////////
      // To compute:
      //      * Time derivative term
      //      * Cell based CFL
      //      * Velocity and soln at quad points (for other models)
      for(int eN=0;eN<nElements_global;eN++)
        {
          //declare local storage for element residual and initialize
          register double
            elementResidual_h[nDOF_test_element],
            elementResidual_hu[nDOF_test_element],
            elementResidual_hv[nDOF_test_element];

          for (int i=0;i<nDOF_test_element;i++)
            {
              elementResidual_h[i]=0.0;
              elementResidual_hu[i]=0.0;
              elementResidual_hv[i]=0.0;
            }
          //
          //loop over quadrature points and compute integrands
          //
          for(int k=0;k<nQuadraturePoints_element;k++)
            {
              //compute indices and declare local storage
              register int eN_k = eN*nQuadraturePoints_element+k,
                eN_k_nSpace = eN_k*nSpace,
                eN_nDOF_trial_element = eN*nDOF_trial_element;
              register double
                h=0.0,hu=0.0,hv=0.0, // solution at current time
                h_old=0.0,hu_old=0.0,hv_old=0.0, // solution at lstage
                jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace],
                h_test_dV[nDOF_trial_element],
                dV,x,y,xt,yt;
              //get jacobian, etc for mapping reference element
              ck.calculateMapping_element(eN,
                                          k,
                                          mesh_dof,
                                          mesh_l2g,
                                          mesh_trial_ref,
                                          mesh_grad_trial_ref,
                                          jac,
                                          jacDet,
                                          jacInv,
                                          x,y);
              //get the physical integration weight
              dV = fabs(jacDet)*dV_ref[k];
              //get the solution at current time
              ck.valFromDOF(h_dof,
			    &h_l2g[eN_nDOF_trial_element],
			    &h_trial_ref[k*nDOF_trial_element],
			    h);
              ck.valFromDOF(hu_dof,
			    &vel_l2g[eN_nDOF_trial_element],
			    &vel_trial_ref[k*nDOF_trial_element],
			    hu);
              ck.valFromDOF(hv_dof,
			    &vel_l2g[eN_nDOF_trial_element],
			    &vel_trial_ref[k*nDOF_trial_element],
			    hv);
              // get the solution at the lstage
              ck.valFromDOF(h_dof_old,
			    &h_l2g[eN_nDOF_trial_element],
			    &h_trial_ref[k*nDOF_trial_element],
			    h_old);
              ck.valFromDOF(hu_dof_old,
			    &vel_l2g[eN_nDOF_trial_element],
			    &vel_trial_ref[k*nDOF_trial_element],
			    hu_old);
              ck.valFromDOF(hv_dof_old,
			    &vel_l2g[eN_nDOF_trial_element],
			    &vel_trial_ref[k*nDOF_trial_element],
			    hv_old);
              // calculate cell based CFL to keep a reference
              calculateCFL(elementDiameter[eN],
                           g,
                           h_old,
                           hu_old,
                           hv_old,
                           hEps,
                           q_cfl[eN_k]);
              //precalculate test function products with integration weights
              for (int j=0;j<nDOF_trial_element;j++)
                h_test_dV[j] = h_test_ref[k*nDOF_trial_element+j]*dV;

              // SAVE VELOCITY // at quadrature points for other models to use
              q_velocity[eN_k_nSpace+0] = 2*h/(h*h+std::pow(fmax(h,hEps),2))*hu;
              q_velocity[eN_k_nSpace+1] = 2*h/(h*h+std::pow(fmax(h,hEps),2))*hv;
              hnp1_at_quad_point[eN_k] = h;
              hunp1_at_quad_point[eN_k] = hu;
              hvnp1_at_quad_point[eN_k] = hv;

              for(int i=0;i<nDOF_test_element;i++)
                {
                  // compute time derivative part of global residual. NOTE: no lumping
                  elementResidual_h[i]  += (h  - h_old)*h_test_dV[i];
                  elementResidual_hu[i] += (hu - hu_old)*h_test_dV[i];
                  elementResidual_hv[i] += (hv - hv_old)*h_test_dV[i];
                }
            }
          // distribute
          for(int i=0;i<nDOF_test_element;i++)
            {
              register int eN_i=eN*nDOF_test_element+i;
              int h_gi = h_l2g[eN_i]; //global i-th index for h
              int vel_gi = vel_l2g[eN_i]; //global i-th index for velocities

              // distribute time derivative to global residual
              globalResidual[offset_h+stride_h*h_gi]  += elementResidual_h[i];
              globalResidual[offset_hu+stride_hu*vel_gi] += elementResidual_hu[i];
              globalResidual[offset_hv+stride_hv*vel_gi] += elementResidual_hv[i];
            }
        }
      // ********** END OF CELL LOOPS ********** //

      if (SECOND_CALL_CALCULATE_RESIDUAL==0) // This is to save some time
        {
          //////////////////////////////////////////////
          // ********** FIRST LOOP ON DOFs ********** //
          //////////////////////////////////////////////
	  // To compute:
	  //     * Entropy at i-th node
	  std::valarray<double> eta(numDOFsPerEqn);
          for (int i=0; i<numDOFsPerEqn; i++)
	    {
	      // COMPUTE ENTROPY. NOTE: WE CONSIDER A FLAT BOTTOM
	      double hi = h_dof_old[i];
	      double one_over_hiReg = 2*hi/(hi*hi+std::pow(fmax(hi,hEps),2)); //hEps
	      eta[i] = ENTROPY(g,hi,hu_dof_old[i],hv_dof_old[i],0.,one_over_hiReg);
	    }
          // ********** END OF FIRST LOOP ON DOFs ********** //

          ///////////////////////////////////////////////
          // ********** SECOND LOOP ON DOFs ********** //
          ///////////////////////////////////////////////
	  // To compute:
	  //     * Hyperbolic part of the flux
	  //     * Extended source term (eqn 6.19)
	  //     * Smoothness indicator
	  //     * global entropy residual
          int ij = 0;
	  std::valarray<double> hyp_flux_h(numDOFsPerEqn), hyp_flux_hu(numDOFsPerEqn), hyp_flux_hv(numDOFsPerEqn),
	    global_entropy_residual(numDOFsPerEqn), psi(numDOFsPerEqn), etaMax(numDOFsPerEqn), etaMin(numDOFsPerEqn);
	  for (int i=0; i<numDOFsPerEqn; i++)
            {
              double hi = h_dof_old[i]; // solution at time tn for the ith DOF
              double hui = hu_dof_old[i];
              double hvi = hv_dof_old[i];
              double one_over_hiReg = 2*hi/(hi*hi+std::pow(fmax(hi,hEps),2)); //hEps
	      double ui = hui*one_over_hiReg;
              double vi = hvi*one_over_hiReg;
	      double mi = lumped_mass_matrix[i];

              // For eta min and max
              etaMax[i] = fabs(eta[i]);
              etaMin[i] = fabs(eta[i]);

	      // COMPUTE EXTENDED SOURCE TERM //
	      // FRICTION //
	      if (LINEAR_FRICTION==1)
		{
		  extendedSourceTerm_hu[i] = mannings*hui*mi;
		  extendedSourceTerm_hv[i] = mannings*hvi*mi;
		}
	      else
		{
		  double veli_norm = std::sqrt(ui*ui+vi*vi);
		  double hi_to_the_gamma = std::pow(hi,gamma);
		  double friction_aux =
		    veli_norm == 0. ? 0.:(2*g*n2*veli_norm*mi/
					  (hi_to_the_gamma +
					   fmax(hi_to_the_gamma,xi*g*n2*dt*veli_norm)));
		  extendedSourceTerm_hu[i] = friction_aux*hui;
		  extendedSourceTerm_hv[i] = friction_aux*hvi;
		}

	      // HYPERBOLIC FLUXES //
	      hyp_flux_h[i]=0;
	      hyp_flux_hu[i]=0;
	      hyp_flux_hv[i]=0;

	      // FOR ENTROPY RESIDUAL //
              double ith_flux_term1=0., ith_flux_term2=0., ith_flux_term3=0.;
              double entropy_flux=0.;

	      // FOR SMOOTHNESS INDICATOR //
              double alphai; // smoothness indicator of solution
              double alpha_numerator = 0;
              double alpha_denominator = 0;

              for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
                { //loop in j (sparsity pattern)
                  int j = csrColumnOffsets_DofLoops[offset];
                  double hj = h_dof_old[j]; // solution at time tn for the jth DOF
                  double huj = hu_dof_old[j];
                  double hvj = hv_dof_old[j];
                  double one_over_hjReg = 2*hj/(hj*hj+std::pow(fmax(hj,hEps),2)); //hEps
                  double uj = huj*one_over_hjReg;
                  double vj = hvj*one_over_hjReg;
		  double Zj = b_dof[j];

		  // auxiliary functions to compute fluxes
		  double aux_h = huj*Cx[ij] + hvj*Cy[ij]; // f1*C = hj*(uj*Cx[ij] + vj*Cy[ij]);
		  double aux_hu = uj*huj*Cx[ij] + uj*hvj*Cy[ij];
		  double aux_hv = vj*huj*Cx[ij] + vj*hvj*Cy[ij];

		  // HYPERBOLIC FLUX //
		  hyp_flux_h[i]  += aux_h;
		  hyp_flux_hu[i] += aux_hu;
		  hyp_flux_hv[i] += aux_hv;

		  // EXTENDED SOURCE //
		  extendedSourceTerm_hu[i] += g*hi*(hj+Zj)*Cx[ij];
		  extendedSourceTerm_hv[i] += g*hi*(hj+Zj)*Cy[ij];

		  // flux for entropy
                  ith_flux_term1 += aux_h;
                  ith_flux_term2 += aux_hu + g*hi*(hj+0.)*Cx[ij]; // NOTE: Zj=0
                  ith_flux_term3 += aux_hv + g*hi*(hj+0.)*Cy[ij]; // NOTE: Zj=0

                  // NOTE: WE CONSIDER FLAT BOTTOM
                  entropy_flux += ( Cx[ij]*ENTROPY_FLUX1(g,hj,huj,hvj,0.,one_over_hjReg) +
                                    Cy[ij]*ENTROPY_FLUX2(g,hj,huj,hvj,0.,one_over_hjReg) );

                  // COMPUTE ETA MIN AND ETA MAX //
                  etaMax[i] = fmax(etaMax[i],fabs(eta[j]));
                  etaMin[i] = fmin(etaMin[i],fabs(eta[j]));

                  // FOR SMOOTHNESS INDICATOR //
                  alpha_numerator += hj - hi;
                  alpha_denominator += fabs(hj - hi);

                  //update ij
                  ij+=1;
                }
	      // COMPUTE ENTROPY RESIDUAL //
              double one_over_entNormFactori = 2./(etaMax[i]-etaMin[i]+1E-15);
              double eta_prime1 = DENTROPY_DH(g,hi,hui,hvi,0.,one_over_hiReg); // NOTE: FLAT BOTTOM
              double eta_prime2 = DENTROPY_DHU(g,hi,hui,hvi,0.,one_over_hiReg);
              double eta_prime3 = DENTROPY_DHV(g,hi,hui,hvi,0.,one_over_hiReg);
              global_entropy_residual[i] = one_over_entNormFactori*
		fabs(entropy_flux -(ith_flux_term1*eta_prime1 +
				    ith_flux_term2*eta_prime2 +
				    ith_flux_term3*eta_prime3));

              // COMPUTE SMOOTHNESS INDICATOR //
              if (hi <= hReg[i]) //hEps, hReg makes the method more robust
                {
                  alphai = 1.;
                  global_entropy_residual[i] = 1E10;
                }
              else
                {
                  if (fabs(alpha_numerator) <= hEps) //hEps. Force alphai=0 in constant states
                    alphai = 0.;
                  else
                    alphai = fabs(alpha_numerator)/(alpha_denominator+1E-15);
                }
              if (POWER_SMOOTHNESS_INDICATOR==0)
                psi[i] = 1.0;
              else
                psi[i] = std::pow(alphai,POWER_SMOOTHNESS_INDICATOR); //NOTE: alpha^2 in the paper
            }
          // ********** END OF 2nd LOOP ON DOFS ********** //

          /////////////////////////////////////////////
          // ********** MAIN LOOP ON DOFs ********** // to compute flux and dissipative terms
          /////////////////////////////////////////////
          ij = 0;
          for (int i=0; i<numDOFsPerEqn; i++)
            {
              double hi = h_dof_old[i];
              double hui = hu_dof_old[i];
              double hvi = hv_dof_old[i];
              double Zi = b_dof[i];
              double mi = lumped_mass_matrix[i];
              double one_over_hiReg = 2*hi/(hi*hi+std::pow(fmax(hi,hEps),2)); // hEps
              double ui = hui*one_over_hiReg;
              double vi = hvi*one_over_hiReg;

              // HIGH ORDER DISSIPATIVE TERMS
              double
                ith_dHij_minus_muHij_times_hStarStates=0.,
                ith_dHij_minus_muHij_times_huStarStates=0.,
                ith_dHij_minus_muHij_times_hvStarStates=0.,
                ith_muHij_times_hStates=0.,
                ith_muHij_times_huStates=0.,
                ith_muHij_times_hvStates=0.;

              // loop over the sparsity pattern of the i-th DOF
              for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
                {
                  int j = csrColumnOffsets_DofLoops[offset];
                  double hj = h_dof_old[j];
                  double huj = hu_dof_old[j];
                  double hvj = hv_dof_old[j];
                  double Zj = b_dof[j];
                  double one_over_hjReg = 2*hj/(hj*hj+std::pow(fmax(hj,hEps),2)); //hEps
                  double uj = huj*one_over_hjReg;
                  double vj = hvj*one_over_hjReg;

                  // COMPUTE STAR SOLUTION // hStar, huStar and hvStar
                  double hStarij  = fmax(0., hi + Zi - fmax(Zi,Zj));
                  double huStarij = hui*hStarij*one_over_hiReg;
                  double hvStarij = hvi*hStarij*one_over_hiReg;

                  double hStarji  = fmax(0., hj + Zj - fmax(Zi,Zj));
                  double huStarji = huj*hStarji*one_over_hjReg;
                  double hvStarji = hvj*hStarji*one_over_hjReg;

                  // Dissipative well balancing term
                  double muLowij = 0., muLij = 0., muHij = 0.;
                  double dLowij = 0., dLij = 0., dHij = 0.;
                  if (i != j) // This is not necessary. See formula for ith_dissipative_terms
                    {
                      ////////////////////////
                      // DISSIPATIVE MATRIX //
                      ////////////////////////
                      if (lstage == 0)
                        dLowij = dLow[ij];
                      else
                        {
                          double cij_norm = sqrt(Cx[ij]*Cx[ij] + Cy[ij]*Cy[ij]);
                          double cji_norm = sqrt(CTx[ij]*CTx[ij] + CTy[ij]*CTy[ij]);
                          double nxij = Cx[ij]/cij_norm, nyij = Cy[ij]/cij_norm;
                          double nxji = CTx[ij]/cji_norm, nyji = CTy[ij]/cji_norm;
                          dLowij = fmax(maxWaveSpeedSharpInitialGuess(g,nxij,nyij,
                                                                      hi,hui,hvi,
                                                                      hj,huj,hvj,
                                                                      hEps,hEps,false)*cij_norm,
                                        maxWaveSpeedSharpInitialGuess(g,nxji,nyji,
                                                                      hj,huj,hvj,
                                                                      hi,hui,hvi,
                                                                      hEps,hEps,false)*cji_norm);
                        }
                      dLij = dLowij;//*fmax(psi[i],psi[j]); // enhance the order to 2nd order. No EV
		      dLow[ij]=dLij; // save dLow for limiting step

                      ///////////////////////////////////////
                      // WELL BALANCING DISSIPATIVE MATRIX //
                      ///////////////////////////////////////
                      muLowij = fmax(fmax(0.,-(ui*Cx[ij] + vi*Cy[ij])),
				     fmax(0,(uj*Cx[ij] + vj*Cy[ij])));
                      muLij = muLowij;//*fmax(psi[i],psi[j]); // enhance the order to 2nd order.

		      ////////////////////////
		      // COMPUTE BAR STATES //
		      ////////////////////////
		      double hBar_ij = 0, hTilde_ij = 0,
			huBar_ij = 0, huTilde_ij = 0,
			hvBar_ij = 0, hvTilde_ij = 0;
		      if (dLij != 0)
			{
			  // h component
			  hBar_ij = -1./(2*dLij)*((huj-hui)*Cx[ij] + (hvj-hvi)*Cy[ij])+0.5*(hj+hi);
			  hTilde_ij = (dLij-muLij)/(2*dLij)*((hStarji-hj)-(hStarij-hi));
			  // hu component
			  huBar_ij = -1./(2*dLij)*((uj*huj-ui*hui)*Cx[ij] +
						   (uj*hvj-ui*hvi)*Cy[ij]) +0.5*(huj+hui);
			  huTilde_ij = (dLij-muLij)/(2*dLij)*((huStarji-huj)-(huStarij-hui));
			  // hv component
			  hvBar_ij = -1./(2*dLij)*((vj*huj-vi*hui)*Cx[ij] +
						   (vj*hvj-vi*hvi)*Cy[ij]) +0.5*(hvj+hvi);
			  hvTilde_ij = (dLij-muLij)/(2*dLij)*((hvStarji-hvj)-(hvStarij-hvi));

			}
		      hBT[ij] = hBar_ij + hTilde_ij;
		      huBT[ij] = huBar_ij + huTilde_ij;
		      hvBT[ij] = hvBar_ij + hvTilde_ij;

                      ///////////////////////
                      // ENTROPY VISCOSITY //
                      ///////////////////////
                      double dEVij = cE*fmax(global_entropy_residual[i],global_entropy_residual[j]);
		      dHij  = fmin(dLowij,dEVij);
		      muHij = fmin(muLowij,dEVij);

                      // compute dij_minus_muij times star solution terms
                      ith_dHij_minus_muHij_times_hStarStates  += (dHij - muHij)*(hStarji-hStarij);
                      ith_dHij_minus_muHij_times_huStarStates += (dHij - muHij)*(huStarji-huStarij);
                      ith_dHij_minus_muHij_times_hvStarStates += (dHij - muHij)*(hvStarji-hvStarij);

                      // compute muij times solution terms
                      ith_muHij_times_hStates  += muHij*(hj-hi);
                      ith_muHij_times_huStates += muHij*(huj-hui);
                      ith_muHij_times_hvStates += muHij*(hvj-hvi);

                      // compute dH_minus_dL
                      dH_minus_dL[ij] = dHij - dLij;
                      muH_minus_muL[ij] = muHij - muLij;
                    }
                  else // i==j
                    {
                      dH_minus_dL[ij]=0.; //Not true but the prod of this times Uj-Ui will be zero
                      muH_minus_muL[ij]=0.; //Not true but the prod of this times Uj-Ui will be zero
                    }
                  // update ij
                  ij+=1;
                }
	      if (LUMPED_MASS_MATRIX==1)
		{
		  globalResidual[offset_h+stride_h*i]
		    = hi - dt/mi*(hyp_flux_h[i]
				  - ith_dHij_minus_muHij_times_hStarStates
				  - ith_muHij_times_hStates);
		  globalResidual[offset_hu+stride_hu*i]
		    = hui - dt/mi*(hyp_flux_hu[i]
				   + extendedSourceTerm_hu[i]
				   - ith_dHij_minus_muHij_times_huStarStates
				   - ith_muHij_times_huStates);
		  globalResidual[offset_hv+stride_hv*i]
		    = hvi - dt/mi*(hyp_flux_hv[i]
				   + extendedSourceTerm_hu[i]
				   - ith_dHij_minus_muHij_times_hvStarStates
				   - ith_muHij_times_hvStates);
		  // clean up potential negative water height due to machine precision
		  if (globalResidual[offset_h+stride_h*i] >= -hEps && globalResidual[offset_h+stride_h*i] < hEps)
		    globalResidual[offset_h+stride_h*i] = 0;
		}
	      else
		{
		  // Distribute residual
		  // NOTE: MASS MATRIX IS CONSISTENT
		  globalResidual[offset_h+stride_h*i]
		    += dt*(hyp_flux_h[i]
			   - ith_dHij_minus_muHij_times_hStarStates
			   - ith_muHij_times_hStates
			   );
		  globalResidual[offset_hu+stride_hu*i]
		    += dt*(hyp_flux_hu[i]
			   + extendedSourceTerm_hu[i]
			   - ith_dHij_minus_muHij_times_huStarStates
			   - ith_muHij_times_huStates
			   );
		  globalResidual[offset_hv+stride_hv*i]
		    += dt*(hyp_flux_hv[i]
			   + extendedSourceTerm_hu[i]
			   - ith_dHij_minus_muHij_times_hvStarStates
			   - ith_muHij_times_hvStates
			   );
		}
	    }
          // ********** END OF LOOP IN DOFs ********** //
        }

      // ********** COMPUTE NORMALS ********** //
      if (COMPUTE_NORMALS==1)
        {
          // This is to identify the normals and create a vector of normal components
          for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++)
            {
              register int ebN = exteriorElementBoundariesArray[ebNE],
                eN  = elementBoundaryElementsArray[ebN*2+0],
                ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+0];
              register double normal[3];
              { // "Loop" in quad points
                int kb = 0; // NOTE: I need to consider just one quad point since the element is not curved so the normal is constant per element
                register int ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb;
                register double
                  jac_ext[nSpace*nSpace],jacDet_ext,jacInv_ext[nSpace*nSpace],boundaryJac[nSpace*(nSpace-1)],
                  metricTensor[(nSpace-1)*(nSpace-1)],metricTensorDetSqrt,x_ext,y_ext;
                /* compute information about mapping from reference element to physical element */
                ck.calculateMapping_elementBoundary(eN,
                                                    ebN_local,
                                                    kb,
                                                    ebN_local_kb,
                                                    mesh_dof,
                                                    mesh_l2g,
                                                    mesh_trial_trace_ref,
                                                    mesh_grad_trial_trace_ref,
                                                    boundaryJac_ref,
                                                    jac_ext,
                                                    jacDet_ext,
                                                    jacInv_ext,
                                                    boundaryJac,
                                                    metricTensor,
                                                    metricTensorDetSqrt,
                                                    normal_ref,
                                                    normal,
                                                    x_ext,y_ext);
              }
              // distribute the normal vectors
              for (int i=0;i<nDOF_test_element;i++)
                {
                  int eN_i = eN*nDOF_test_element+i;
                  int gi = h_l2g[eN_i];
                  normalx[gi] += 0.5*normal[0]*(i==ebN_local ? 0. : 1.);
                  normaly[gi] += 0.5*normal[1]*(i==ebN_local ? 0. : 1.);
                }
            }
          // normalize
          for (int gi=0; gi<numDOFsPerEqn; gi++)
            {
              double norm_factor = sqrt(std::pow(normalx[gi],2) + std::pow(normaly[gi],2));
              if (norm_factor != 0)
                {
                  normalx[gi] /= norm_factor;
                  normaly[gi] /= norm_factor;
                }
            }
        }
      // ********** END OF COMPUTING NORMALS ********** //
    }

    void calculateMassMatrix(//element
                             double* mesh_trial_ref,
                             double* mesh_grad_trial_ref,
                             double* mesh_dof,
                             double* mesh_velocity_dof,
                             double MOVING_DOMAIN,
                             int* mesh_l2g,
                             double* dV_ref,
                             double* h_trial_ref,
                             double* h_grad_trial_ref,
                             double* h_test_ref,
                             double* h_grad_test_ref,
                             double* vel_trial_ref,
                             double* vel_grad_trial_ref,
                             double* vel_test_ref,
                             double* vel_grad_test_ref,
                             //element boundary
                             double* mesh_trial_trace_ref,
                             double* mesh_grad_trial_trace_ref,
                             double* dS_ref,
                             double* h_trial_trace_ref,
                             double* h_grad_trial_trace_ref,
                             double* h_test_trace_ref,
                             double* h_grad_test_trace_ref,
                             double* vel_trial_trace_ref,
                             double* vel_grad_trial_trace_ref,
                             double* vel_test_trace_ref,
                             double* vel_grad_test_trace_ref,
                             double* normal_ref,
                             double* boundaryJac_ref,
                             //physics
                             double* elementDiameter,
                             int nElements_global,
                             double useRBLES,
                             double useMetrics,
                             double alphaBDF,
                             double nu,
                             double g,
                             int* h_l2g,
                             int* vel_l2g,
                             double* b_dof,
                             double* h_dof,
                             double* hu_dof,
                             double* hv_dof,
                             double* h_dof_sge,
                             double* hu_dof_sge,
                             double* hv_dof_sge,
                             double* q_mass_acc_beta_bdf,
                             double* q_mom_hu_acc_beta_bdf,
                             double* q_mom_hv_acc_beta_bdf,
                             double* q_cfl,
                             int* sdInfo_hu_hu_rowptr,
                             int* sdInfo_hu_hu_colind,
                             int* sdInfo_hu_hv_rowptr,
                             int* sdInfo_hu_hv_colind,
                             int* sdInfo_hv_hv_rowptr,
                             int* sdInfo_hv_hv_colind,
                             int* sdInfo_hv_hu_rowptr,
                             int* sdInfo_hv_hu_colind,
                             int* csrRowIndeces_h_h,
                             int* csrColumnOffsets_h_h,
                             int* csrRowIndeces_h_hu,
                             int* csrColumnOffsets_h_hu,
                             int* csrRowIndeces_h_hv,
                             int* csrColumnOffsets_h_hv,
                             int* csrRowIndeces_hu_h,
                             int* csrColumnOffsets_hu_h,
                             int* csrRowIndeces_hu_hu,
                             int* csrColumnOffsets_hu_hu,
                             int* csrRowIndeces_hu_hv,
                             int* csrColumnOffsets_hu_hv,
                             int* csrRowIndeces_hv_h,
                             int* csrColumnOffsets_hv_h,
                             int* csrRowIndeces_hv_hu,
                             int* csrColumnOffsets_hv_hu,
                             int* csrRowIndeces_hv_hv,
                             int* csrColumnOffsets_hv_hv,
                             double* globalJacobian,
                             int nExteriorElementBoundaries_global,
                             int* exteriorElementBoundariesArray,
                             int* elementBoundaryElementsArray,
                             int* elementBoundaryLocalElementBoundariesArray,
                             int* isDOFBoundary_h,
                             int* isDOFBoundary_hu,
                             int* isDOFBoundary_hv,
                             int* isAdvectiveFluxBoundary_h,
                             int* isAdvectiveFluxBoundary_hu,
                             int* isAdvectiveFluxBoundary_hv,
                             int* isDiffusiveFluxBoundary_hu,
                             int* isDiffusiveFluxBoundary_hv,
                             double* ebqe_bc_h_ext,
                             double* ebqe_bc_flux_mass_ext,
                             double* ebqe_bc_flux_mom_hu_adv_ext,
                             double* ebqe_bc_flux_mom_hv_adv_ext,
                             double* ebqe_bc_hu_ext,
                             double* ebqe_bc_flux_hu_diff_ext,
                             double* ebqe_penalty_ext,
                             double* ebqe_bc_hv_ext,
                             double* ebqe_bc_flux_hv_diff_ext,
                             int* csrColumnOffsets_eb_h_h,
                             int* csrColumnOffsets_eb_h_hu,
                             int* csrColumnOffsets_eb_h_hv,
                             int* csrColumnOffsets_eb_hu_h,
                             int* csrColumnOffsets_eb_hu_hu,
                             int* csrColumnOffsets_eb_hu_hv,
                             int* csrColumnOffsets_eb_hv_h,
                             int* csrColumnOffsets_eb_hv_hu,
                             int* csrColumnOffsets_eb_hv_hv,
                             double dt)
    {
      //
      //loop over elements to compute volume integrals and load them into the element Jacobians and global Jacobian
      //
      for(int eN=0;eN<nElements_global;eN++)
        {
          register double
            elementJacobian_h_h[nDOF_test_element][nDOF_trial_element],
            elementJacobian_hu_hu[nDOF_test_element][nDOF_trial_element],
            elementJacobian_hv_hv[nDOF_test_element][nDOF_trial_element];
          for (int i=0;i<nDOF_test_element;i++)
            for (int j=0;j<nDOF_trial_element;j++)
              {
                elementJacobian_h_h[i][j]=0.0;
                elementJacobian_hu_hu[i][j]=0.0;
                elementJacobian_hv_hv[i][j]=0.0;
              }
          for  (int k=0;k<nQuadraturePoints_element;k++)
            {
              int eN_k = eN*nQuadraturePoints_element+k, //index to a scalar at a quadrature point
                eN_k_nSpace = eN_k*nSpace,
                eN_nDOF_trial_element = eN*nDOF_trial_element; //index to a vector at a quadrature point

              //declare local storage
              register double
                jac[nSpace*nSpace],
                jacDet,
                jacInv[nSpace*nSpace],
                dV,
                h_test_dV[nDOF_test_element],
                vel_test_dV[nDOF_test_element],
                x,y,xt,yt;
              //get jacobian, etc for mapping reference element
              ck.calculateMapping_element(eN,
                                          k,
                                          mesh_dof,
                                          mesh_l2g,
                                          mesh_trial_ref,
                                          mesh_grad_trial_ref,
                                          jac,
                                          jacDet,
                                          jacInv,
                                          x,y);
              //get the physical integration weight
              dV = fabs(jacDet)*dV_ref[k];
              //precalculate test function products with integration weights
              for (int j=0;j<nDOF_trial_element;j++)
                {
                  h_test_dV[j] = h_test_ref[k*nDOF_trial_element+j]*dV;
                  vel_test_dV[j] = vel_test_ref[k*nDOF_trial_element+j]*dV;
                }
              for(int i=0;i<nDOF_test_element;i++)
                {
                  register int i_nSpace = i*nSpace;
                  for(int j=0;j<nDOF_trial_element;j++)
                    {
                      register int j_nSpace = j*nSpace;
                      elementJacobian_h_h[i][j] += h_trial_ref[k*nDOF_trial_element+j]*h_test_dV[i];
                      elementJacobian_hu_hu[i][j] += vel_trial_ref[k*nDOF_trial_element+j]*vel_test_dV[i];
                      elementJacobian_hv_hv[i][j] += vel_trial_ref[k*nDOF_trial_element+j]*vel_test_dV[i];
                    }//j
                }//i
            }//k
          //
          //load into element Jacobian into global Jacobian
          //
          for (int i=0;i<nDOF_test_element;i++)
            {
              register int eN_i = eN*nDOF_test_element+i;
              for (int j=0;j<nDOF_trial_element;j++)
                {
                  register int eN_i_j = eN_i*nDOF_trial_element+j;
                  globalJacobian[csrRowIndeces_h_h[eN_i] +
				 csrColumnOffsets_h_h[eN_i_j]] += elementJacobian_h_h[i][j];
                  globalJacobian[csrRowIndeces_hu_hu[eN_i] +
				 csrColumnOffsets_hu_hu[eN_i_j]] += elementJacobian_hu_hu[i][j];
                  globalJacobian[csrRowIndeces_hv_hv[eN_i] +
				 csrColumnOffsets_hv_hv[eN_i_j]] += elementJacobian_hv_hv[i][j];
                }//j
            }//i
        }//elements

    }

    void calculateLumpedMassMatrix(//element
                                   double* mesh_trial_ref,
                                   double* mesh_grad_trial_ref,
                                   double* mesh_dof,
                                   double* mesh_velocity_dof,
                                   double MOVING_DOMAIN,
                                   int* mesh_l2g,
                                   double* dV_ref,
                                   double* h_trial_ref,
                                   double* h_grad_trial_ref,
                                   double* h_test_ref,
                                   double* h_grad_test_ref,
                                   double* vel_trial_ref,
                                   double* vel_grad_trial_ref,
                                   double* vel_test_ref,
                                   double* vel_grad_test_ref,
                                   //element boundary
                                   double* mesh_trial_trace_ref,
                                   double* mesh_grad_trial_trace_ref,
                                   double* dS_ref,
                                   double* h_trial_trace_ref,
                                   double* h_grad_trial_trace_ref,
                                   double* h_test_trace_ref,
                                   double* h_grad_test_trace_ref,
                                   double* vel_trial_trace_ref,
                                   double* vel_grad_trial_trace_ref,
                                   double* vel_test_trace_ref,
                                   double* vel_grad_test_trace_ref,
                                   double* normal_ref,
                                   double* boundaryJac_ref,
                                   //physics
                                   double* elementDiameter,
                                   int nElements_global,
                                   double useRBLES,
                                   double useMetrics,
                                   double alphaBDF,
                                   double nu,
                                   double g,
                                   int* h_l2g,
                                   int* vel_l2g,
                                   double* b_dof,
                                   double* h_dof,
                                   double* hu_dof,
                                   double* hv_dof,
                                   double* h_dof_sge,
                                   double* hu_dof_sge,
                                   double* hv_dof_sge,
                                   double* q_mass_acc_beta_bdf,
                                   double* q_mom_hu_acc_beta_bdf,
                                   double* q_mom_hv_acc_beta_bdf,
                                   double* q_cfl,
                                   int* sdInfo_hu_hu_rowptr,
                                   int* sdInfo_hu_hu_colind,
                                   int* sdInfo_hu_hv_rowptr,
                                   int* sdInfo_hu_hv_colind,
                                   int* sdInfo_hv_hv_rowptr,
                                   int* sdInfo_hv_hv_colind,
                                   int* sdInfo_hv_hu_rowptr,
                                   int* sdInfo_hv_hu_colind,
                                   int* csrRowIndeces_h_h,
                                   int* csrColumnOffsets_h_h,
                                   int* csrRowIndeces_h_hu,
                                   int* csrColumnOffsets_h_hu,
                                   int* csrRowIndeces_h_hv,
                                   int* csrColumnOffsets_h_hv,
                                   int* csrRowIndeces_hu_h,
                                   int* csrColumnOffsets_hu_h,
                                   int* csrRowIndeces_hu_hu,
                                   int* csrColumnOffsets_hu_hu,
                                   int* csrRowIndeces_hu_hv,
                                   int* csrColumnOffsets_hu_hv,
                                   int* csrRowIndeces_hv_h,
                                   int* csrColumnOffsets_hv_h,
                                   int* csrRowIndeces_hv_hu,
                                   int* csrColumnOffsets_hv_hu,
                                   int* csrRowIndeces_hv_hv,
                                   int* csrColumnOffsets_hv_hv,
                                   double* globalJacobian,
                                   int nExteriorElementBoundaries_global,
                                   int* exteriorElementBoundariesArray,
                                   int* elementBoundaryElementsArray,
                                   int* elementBoundaryLocalElementBoundariesArray,
                                   int* isDOFBoundary_h,
                                   int* isDOFBoundary_hu,
                                   int* isDOFBoundary_hv,
                                   int* isAdvectiveFluxBoundary_h,
                                   int* isAdvectiveFluxBoundary_hu,
                                   int* isAdvectiveFluxBoundary_hv,
                                   int* isDiffusiveFluxBoundary_hu,
                                   int* isDiffusiveFluxBoundary_hv,
                                   double* ebqe_bc_h_ext,
                                   double* ebqe_bc_flux_mass_ext,
                                   double* ebqe_bc_flux_mom_hu_adv_ext,
                                   double* ebqe_bc_flux_mom_hv_adv_ext,
                                   double* ebqe_bc_hu_ext,
                                   double* ebqe_bc_flux_hu_diff_ext,
                                   double* ebqe_penalty_ext,
                                   double* ebqe_bc_hv_ext,
                                   double* ebqe_bc_flux_hv_diff_ext,
                                   int* csrColumnOffsets_eb_h_h,
                                   int* csrColumnOffsets_eb_h_hu,
                                   int* csrColumnOffsets_eb_h_hv,
                                   int* csrColumnOffsets_eb_hu_h,
                                   int* csrColumnOffsets_eb_hu_hu,
                                   int* csrColumnOffsets_eb_hu_hv,
                                   int* csrColumnOffsets_eb_hv_h,
                                   int* csrColumnOffsets_eb_hv_hu,
                                   int* csrColumnOffsets_eb_hv_hv,
                                   double dt)
    {
      //
      //loop over elements to compute volume integrals and load them into the element Jacobians and global Jacobian
      //
      for(int eN=0;eN<nElements_global;eN++)
        {
          register double
	    elementJacobian_h_h[nDOF_test_element][nDOF_trial_element],
            elementJacobian_hu_hu[nDOF_test_element][nDOF_trial_element],
            elementJacobian_hv_hv[nDOF_test_element][nDOF_trial_element];
          for (int i=0;i<nDOF_test_element;i++)
            for (int j=0;j<nDOF_trial_element;j++)
              {
                elementJacobian_h_h[i][j]=0.0;
                elementJacobian_hu_hu[i][j]=0.0;
                elementJacobian_hv_hv[i][j]=0.0;
              }
          for  (int k=0;k<nQuadraturePoints_element;k++)
            {
              int eN_k = eN*nQuadraturePoints_element+k, //index to a scalar at a quadrature point
                eN_k_nSpace = eN_k*nSpace,
                eN_nDOF_trial_element = eN*nDOF_trial_element; //index to a vector at a quadrature point

              //declare local storage
              register double
                jac[nSpace*nSpace],
                jacDet,
                jacInv[nSpace*nSpace],
                dV,
                h_test_dV[nDOF_test_element],
                vel_test_dV[nDOF_test_element],
                x,y,xt,yt;
              //get jacobian, etc for mapping reference element
              ck.calculateMapping_element(eN,
                                          k,
                                          mesh_dof,
                                          mesh_l2g,
                                          mesh_trial_ref,
                                          mesh_grad_trial_ref,
                                          jac,
                                          jacDet,
                                          jacInv,
                                          x,y);
              //get the physical integration weight
              dV = fabs(jacDet)*dV_ref[k];
              //precalculate test function products with integration weights
              for (int j=0;j<nDOF_trial_element;j++)
                {
                  h_test_dV[j] = h_test_ref[k*nDOF_trial_element+j]*dV;
                  vel_test_dV[j] = vel_test_ref[k*nDOF_trial_element+j]*dV;
                }

              for(int i=0;i<nDOF_test_element;i++)
                {
                  register int i_nSpace = i*nSpace;
                  for(int j=0;j<nDOF_trial_element;j++)
                    {
                      register int j_nSpace = j*nSpace;
                      elementJacobian_h_h[i][j] += (i==j ? 1.0 : 0.0)*h_test_dV[i];
                      elementJacobian_hu_hu[i][j] += (i==j ? 1.0 : 0.0)*vel_test_dV[i];
                      elementJacobian_hv_hv[i][j] += (i==j ? 1.0 : 0.0)*vel_test_dV[i];
                    }//j
                }//i
            }//k
          //
          //load into element Jacobian into global Jacobian
          //
          for (int i=0;i<nDOF_test_element;i++)
            {
              register int eN_i = eN*nDOF_test_element+i;
              for (int j=0;j<nDOF_trial_element;j++)
                {
                  register int eN_i_j = eN_i*nDOF_trial_element+j;
                  globalJacobian[csrRowIndeces_h_h[eN_i] + csrColumnOffsets_h_h[eN_i_j]] += elementJacobian_h_h[i][j];
                  globalJacobian[csrRowIndeces_hu_hu[eN_i] + csrColumnOffsets_hu_hu[eN_i_j]] += elementJacobian_hu_hu[i][j];
                  globalJacobian[csrRowIndeces_hv_hv[eN_i] + csrColumnOffsets_hv_hv[eN_i_j]] += elementJacobian_hv_hv[i][j];
                }//j
            }//i
        }//elements
    }
  };//SW2DCV

  inline SW2DCV_base* newSW2DCV(int nSpaceIn,
                            int nQuadraturePoints_elementIn,
                            int nDOF_mesh_trial_elementIn,
                            int nDOF_trial_elementIn,
                            int nDOF_test_elementIn,
                            int nQuadraturePoints_elementBoundaryIn,
                            int CompKernelFlag)
  {
    return proteus::chooseAndAllocateDiscretization2D<SW2DCV_base,SW2DCV,CompKernel>(nSpaceIn,
                                                                                 nQuadraturePoints_elementIn,
                                                                                 nDOF_mesh_trial_elementIn,
                                                                                 nDOF_trial_elementIn,
                                                                                 nDOF_test_elementIn,
                                                                                 nQuadraturePoints_elementBoundaryIn,
                                                                                 CompKernelFlag);
  }
}//proteus

#endif

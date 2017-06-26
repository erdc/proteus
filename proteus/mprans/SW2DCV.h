#ifndef SW2DCV_H
#define SW2DCV_H
#include <cmath>
#include <iostream>
#include "CompKernel.h"
#include "ModelFactory.h"
#include <assert.h>

//cek todo
//2. Get stabilization right
//3. Add Riemann solvers for external flux
//4. Add Riemann solvers for internal flux and DG terms 
//5. Try other choices of variables h,hu,hv, Bova-Carey symmetrization?

#define cMax 0.25
#define IMPLICIT 0
#define POWER_SMOOTHNESS_INDICATOR 2
#define BETAij 1
#define LINEAR_FRICTION 0

// FOR CELL BASED ENTROPY VISCOSITY 
#define ENTROPY(g,h,hu,hv,one_over_hReg) 0.5*(g*h*h+hu*hu*one_over_hReg+hv*hv*one_over_hReg)
#define DENTROPY_DH(g,h,hu,hv,one_over_hReg) g*h - 0.5*(hu*hu+hv*hv)*std::pow(one_over_hReg,2) 
#define DENTROPY_DHU(g,h,hu,hv,one_over_hReg) hu*one_over_hReg
#define DENTROPY_DHV(g,h,hu,hv,one_over_hReg) hv*one_over_hReg
#define D_ENTROPY(g,h,hu,hv,hx,hux,hvx,one_over_hnReg) g*h*hx + hu*one_over_hnReg*(hux-0.5*hx*(hu*one_over_hnReg)) + hv*one_over_hnReg*(hvx-0.5*hx*(hv*one_over_hnReg))

#define ENTROPY_FLUX1(g,h,hu,hv,one_over_hReg) (ENTROPY(g,h,hu,hv,one_over_hReg)+0.5*g*h*h)*hu*one_over_hReg
#define ENTROPY_FLUX2(g,h,hu,hv,one_over_hReg) (ENTROPY(g,h,hu,hv,one_over_hReg)+0.5*g*h*h)*hv*one_over_hReg

  // FOR INVARIANT DOMAIN PRESERVING 
#define f(g,h,hZ) ( (h <= hZ) ? 2.*(sqrt(g*h)-sqrt(g*hZ)) : (h-hZ)*sqrt(0.5*g*(h+hZ)/h/hZ) )
#define phi(g,h,hL,hR,uL,uR) ( f(g,h,hL) + f(g,h,hR) + uR - uL )

#define fp(g,h,hZ) ( (h <= hZ) ? sqrt(g/h) : g*(2*h*h+h*hZ+hZ*hZ)/(2*sqrt(2*g)*h*h*hZ*sqrt(1/h+1/hZ)) )
#define phip(g,h,hL,hR) ( fp(g,h,hL) + fp(g,h,hR) )

#define nu1(g,hStar,hL,uL) ( uL - sqrt(g*hL)*sqrt( (1+fmax((hStar-hL)/2/hL,0.0)) * (1+fmax((hStar-hL)/hL,0.)) ) )
#define nu3(g,hStar,hR,uR) ( uR + sqrt(g*hR)*sqrt( (1+fmax((hStar-hR)/2/hR,0.0)) * (1+fmax((hStar-hR)/hR,0.)) ) )

#define phiDiff(g,h1k,h2k,hL,hR,uL,uR)   ( (phi(g,h2k,hL,hR,uL,uR) - phi(g,h1k,hL,hR,uL,uR))/(h2k-h1k)    )
#define phiDDiff1(g,h1k,h2k,hL,hR,uL,uR) ( (phiDiff(g,h1k,h2k,hL,hR,uL,uR) - phip(g,h1k,hL,hR))/(h2k-h1k) )
#define phiDDiff2(g,h1k,h2k,hL,hR,uL,uR) ( (phip(g,h2k,hL,hR) - phiDiff(g,h1k,h2k,hL,hR,uL,uR))/(h2k-h1k) )

#define hStarLFromQuadPhiFromAbove(g,hStarL,hStarR,hL,hR,uL,uR) ( hStarL-2*phi(g,hStarL,hL,hR,uL,uR)/(phip(g,hStarL,hL,hR)+sqrt(std::pow(phip(g,hStarL,hL,hR),2)-4*phi(g,hStarL,hL,hR,uL,uR)*phiDDiff1(g,hStarL,hStarR,hL,hR,uL,uR))) )
#define hStarRFromQuadPhiFromBelow(g,hStarL,hStarR,hL,hR,uL,uR) ( hStarR-2*phi(g,hStarR,hL,hR,uL,uR)/(phip(g,hStarR,hL,hR)+sqrt(std::pow(phip(g,hStarR,hL,hR),2)-4*phi(g,hStarR,hL,hR,uL,uR)*phiDDiff2(g,hStarL,hStarR,hL,hR,uL,uR))) )

namespace proteus
{
  class SW2DCV_base
  {
  public:
    virtual ~SW2DCV_base(){}
    virtual void FCTStep(double dt, 
			 int NNZ, //number on non-zero entries on sparsity pattern
			 int numDOFs, //number of DOFs
			 double* lumped_mass_matrix, //lumped mass matrix (as vector)
			 double* h_lstage, //DOFs of solution at last stage
			 double* hu_lstage, 
			 double* hv_lstage, 
			 double* b_dof,
			 double* high_order_hnp1, //DOFs of high order solution at tnp1
			 double* high_order_hunp1, 
			 double* high_order_hvnp1, 
			 double* low_order_hnp1, //operators to construct low order solution
			 double* low_order_hunp1,
			 double* low_order_hvnp1,
			 double* limited_hnp1,
			 double* limited_hunp1,
			 double* limited_hvnp1,
			 int* csrRowIndeces_DofLoops, //csr row indeces 
			 int* csrColumnOffsets_DofLoops, //csr column offsets 
			 double* MassMatrix, //mass matrix
			 double* dEV_minus_dL,
			 double hEps, 
			 int LUMPED_MASS_MATRIX
			 )=0;
    virtual void calculateResidual_SUPG(//element
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
				   double* h_dof_old_old,
 				   double* hu_dof_old_old, 
				   double* hv_dof_old_old,
				   double* h_dof_lstage,
 				   double* hu_dof_lstage, 
				   double* hv_dof_lstage,
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
				   double* q_velocity_sge,
				   double* q_cfl,
				   double* q_numDiff_h,
				   double* q_numDiff_hu, 
				   double* q_numDiff_hv,
				   double* q_numDiff_h_last, 
				   double* q_numDiff_hu_last, 
				   double* q_numDiff_hv_last,
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
				   // FOR EDGE BASED METHODS
				   double* lumped_mass_matrix, 
				   double* edge_based_cfl, 
				   double cfl_run,
				   double hEps, 
				   int recompute_lumped_mass_matrix,
				   // SAVE SOLUTION (mql)
				   double* hnp1_at_quad_point,
				   double* hunp1_at_quad_point,
				   double* hvnp1_at_quad_point,
				   // GALERKIN SOLUTION 
				   double* h_dof_galerkin,
				   double* hu_dof_galerkin, 
				   double* hv_dof_galerkin,
				   // TO COMPUTE LOW ORDER 
				   double* low_order_hnp1,
				   double* low_order_hunp1,
				   double* low_order_hvnp1,
				   double* dEV_minus_dL,
				   double cE, 
				   int LUMPED_MASS_MATRIX, 
				   int USE_EV_BASED_ON_GALERKIN,
				   double dt,
				   double mannings
				   )=0;
    virtual void calculateResidual_cell_based_entropy_viscosity(//element
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
								double* h_dof_old_old, 
								double* hu_dof_old_old, 
								double* hv_dof_old_old,
								double* h_dof_lstage, 
								double* hu_dof_lstage, 
								double* hv_dof_lstage,
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
								double* q_velocity_sge,
								double* q_cfl,
								double* q_numDiff_h,
								double* q_numDiff_hu, 
								double* q_numDiff_hv,
								double* q_numDiff_h_last, 
								double* q_numDiff_hu_last, 
								double* q_numDiff_hv_last,
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
								double* edge_based_cfl,
								double cfl_run,
								double hEps,
								int recompute_lumped_mass_matrix,
								// SAVE SOLUTION (mql)
								double* hnp1_at_quad_point,
								double* hunp1_at_quad_point,
								double* hvnp1_at_quad_point,
								// GALERKIN SOLUTION 
								double* h_dof_galerkin,
								double* hu_dof_galerkin, 
								double* hv_dof_galerkin,
								// TO COMPUTE LOW ORDER 
								double* low_order_hnp1,
								double* low_order_hunp1,
								double* low_order_hvnp1,
								double* dEV_minus_dL,
								double cE,
								int LUMPED_MASS_MATRIX,
								int USE_EV_BASED_ON_GALERKIN,
								double dt,
								double mannings
								)=0;
    virtual void calculateResidual_first_order_flatB_GP(//element
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
					    double* h_dof_old_old, 
					    double* hu_dof_old_old, 
					    double* hv_dof_old_old,
					    double* h_dof_lstage, 
					    double* hu_dof_lstage, 
					    double* hv_dof_lstage,
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
					    double* q_velocity_sge,
					    double* q_cfl,
					    double* q_numDiff_h,
					    double* q_numDiff_hu, 
					    double* q_numDiff_hv,
					    double* q_numDiff_h_last, 
					    double* q_numDiff_hu_last, 
					    double* q_numDiff_hv_last,
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
					    double* edge_based_cfl,
					    double cfl_run,
					    double hEps,
					    int recompute_lumped_mass_matrix,
					    // SAVE SOLUTION (mql)
					    double* hnp1_at_quad_point,
					    double* hunp1_at_quad_point,
					    double* hvnp1_at_quad_point,
					    // GALERKIN SOLUTION 
					    double* h_dof_galerkin,
					    double* hu_dof_galerkin, 
					    double* hv_dof_galerkin,
					    // TO COMPUTE LOW ORDER 
					    double* low_order_hnp1,
					    double* low_order_hunp1,
					    double* low_order_hvnp1,
					    double* dEV_minus_dL,
					    double cE,
					    int LUMPED_MASS_MATRIX,
					    int USE_EV_BASED_ON_GALERKIN,
					    double dt,
					    double mannings
							)=0;
    virtual void calculateResidual_second_order_flatB_GP(//element
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
					       double* h_dof_old_old, 
					       double* hu_dof_old_old, 
					       double* hv_dof_old_old,
					       double* h_dof_lstage, 
					       double* hu_dof_lstage, 
					       double* hv_dof_lstage,
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
					       double* q_velocity_sge,
					       double* q_cfl,
					       double* q_numDiff_h,
					       double* q_numDiff_hu, 
					       double* q_numDiff_hv,
					       double* q_numDiff_h_last, 
					       double* q_numDiff_hu_last, 
					       double* q_numDiff_hv_last,
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
					       double* edge_based_cfl, 
					       double cfl_run,
					       double hEps,
					       int recompute_lumped_mass_matrix,
					       // SAVE SOLUTION (mql)
					       double* hnp1_at_quad_point,
					       double* hunp1_at_quad_point,
					       double* hvnp1_at_quad_point,
					       // GALERKIN SOLUTION 
					       double* h_dof_galerkin,
					       double* hu_dof_galerkin, 
					       double* hv_dof_galerkin,
					       // TO COMPUTE LOW ORDER 
					       double* low_order_hnp1,
					       double* low_order_hunp1,
					       double* low_order_hvnp1,
					       double* dEV_minus_dL,
					       double cE,
					       int LUMPED_MASS_MATRIX,
					       int USE_EV_BASED_ON_GALERKIN,
					       double dt,
					       double mannings
							 )=0;
    virtual void calculateResidual_second_order_NonFlatB_GP(//element
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
					       double* h_dof_old_old, 
					       double* hu_dof_old_old, 
					       double* hv_dof_old_old,
					       double* h_dof_lstage, 
					       double* hu_dof_lstage, 
					       double* hv_dof_lstage,
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
					       double* q_velocity_sge,
					       double* q_cfl,
					       double* q_numDiff_h,
					       double* q_numDiff_hu, 
					       double* q_numDiff_hv,
					       double* q_numDiff_h_last, 
					       double* q_numDiff_hu_last, 
					       double* q_numDiff_hv_last,
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
					       double* edge_based_cfl, 
					       double cfl_run,
					       double hEps,
					       int recompute_lumped_mass_matrix,
					       // SAVE SOLUTION (mql)
					       double* hnp1_at_quad_point,
					       double* hunp1_at_quad_point,
					       double* hvnp1_at_quad_point,
					       // GALERKIN SOLUTION 
					       double* h_dof_galerkin,
					       double* hu_dof_galerkin, 
					       double* hv_dof_galerkin,
					       // TO COMPUTE LOW ORDER 
					       double* low_order_hnp1,
					       double* low_order_hunp1,
					       double* low_order_hvnp1,
					       double* dEV_minus_dL,
					       double cE,
					       int LUMPED_MASS_MATRIX,
					       int USE_EV_BASED_ON_GALERKIN,
					       double dt,
					       double mannings
							    )=0;
    virtual void calculateResidual_second_order_NonFlatB_with_EV(//element
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
					       double* h_dof_old_old, 
					       double* hu_dof_old_old, 
					       double* hv_dof_old_old,
					       double* h_dof_lstage, 
					       double* hu_dof_lstage, 
					       double* hv_dof_lstage,
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
					       double* q_velocity_sge,
					       double* q_cfl,
					       double* q_numDiff_h,
					       double* q_numDiff_hu, 
					       double* q_numDiff_hv,
					       double* q_numDiff_h_last, 
					       double* q_numDiff_hu_last, 
					       double* q_numDiff_hv_last,
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
					       double* edge_based_cfl, 
					       double cfl_run,
					       double hEps,
					       int recompute_lumped_mass_matrix,
					       // SAVE SOLUTION (mql)
					       double* hnp1_at_quad_point,
					       double* hunp1_at_quad_point,
					       double* hvnp1_at_quad_point,
					       // GALERKIN SOLUTION 
					       double* h_dof_galerkin,
					       double* hu_dof_galerkin, 
					       double* hv_dof_galerkin,
					       // TO COMPUTE LOW ORDER 
					       double* low_order_hnp1,
					       double* low_order_hunp1,
					       double* low_order_hvnp1,
					       double* dEV_minus_dL,
					       double cE,
					       int LUMPED_MASS_MATRIX,
					       int USE_EV_BASED_ON_GALERKIN,
					       double dt,
					       double mannings
								 )=0;
    virtual void calculateResidual_galerkin(//element
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
					    double* h_dof_old_old, 
					    double* hu_dof_old_old, 
					    double* hv_dof_old_old,
					    double* h_dof_lstage, 
					    double* hu_dof_lstage, 
					    double* hv_dof_lstage,
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
					    double* q_velocity_sge,
					    double* q_cfl,
					    double* q_numDiff_h,
					    double* q_numDiff_hu, 
					    double* q_numDiff_hv,
					    double* q_numDiff_h_last, 
					    double* q_numDiff_hu_last, 
					    double* q_numDiff_hv_last,
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
					    double* edge_based_cfl, 
					    double cfl_run,
					    double hEps,
					    int recompute_lumped_mass_matrix,
					    // SAVE SOLUTION (mql)
					    double* hnp1_at_quad_point,
					    double* hunp1_at_quad_point,
					    double* hvnp1_at_quad_point,
					    // GALERKIN SOLUTION 
					    double* h_dof_galerkin,
					    double* hu_dof_galerkin, 
					    double* hv_dof_galerkin,
					    // TO COMPUTE LOW ORDER 
					    double* low_order_hnp1,
					    double* low_order_hunp1,
					    double* low_order_hvnp1,
					    double* dEV_minus_dL,
					    double cE,
					    int LUMPED_MASS_MATRIX,
					    int USE_EV_BASED_ON_GALERKIN,
					    double dt,
					    double mannings
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
				   double* h_dof_old_old, 
				   double* hu_dof_old_old, 
				   double* hv_dof_old_old,
				   double* h_dof_lstage, 
				   double* hu_dof_lstage, 
				   double* hv_dof_lstage,
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
				   double* q_velocity_sge,
				   double* q_cfl,
				   double* q_numDiff_h,
				   double* q_numDiff_hu, 
				   double* q_numDiff_hv,
				   double* q_numDiff_h_last, 
				   double* q_numDiff_hu_last, 
				   double* q_numDiff_hv_last,
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
				   double* edge_based_cfl, 
				   double cfl_run,
				   double hEps,
				   int recompute_lumped_mass_matrix,
				   // SAVE SOLUTION (mql)
				   double* hnp1_at_quad_point,
				   double* hunp1_at_quad_point,
				   double* hvnp1_at_quad_point,
				   // GALERKIN SOLUTION 
				   double* h_dof_galerkin,
				   double* hu_dof_galerkin, 
				   double* hv_dof_galerkin,
				   // TO COMPUTE LOW ORDER 
				   double* low_order_hnp1,
				   double* low_order_hunp1,
				   double* low_order_hvnp1,
				   double* dEV_minus_dL,
				   double cE,
				   int LUMPED_MASS_MATRIX,
				   int USE_EV_BASED_ON_GALERKIN,
				   double dt,
				   double mannings,
				   // Quant of interests
				   double* quantDOFs,
				   // Lumped mass matrix
				   double* ML, 
				   int SECOND_CALL_CALCULATE_RESIDUAL
				   )=0;
    virtual void calculateJacobian(//element
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
				   double* q_velocity_sge,
				   double* q_cfl,
				   double* q_numDiff_h_last,
				   double* q_numDiff_hu_last, 
				   double* q_numDiff_hv_last,
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
    virtual void calculateJacobian_cell_based_entropy_viscosity(//element
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
								double* q_velocity_sge,
								double* q_cfl,
								double* q_numDiff_h_last,
								double* q_numDiff_hu_last, 
								double* q_numDiff_hv_last,
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
				     double* q_velocity_sge,
				     double* q_cfl,
				     double* q_numDiff_h_last,
				     double* q_numDiff_hu_last, 
				     double* q_numDiff_hv_last,
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
					 double* q_velocity_sge,
					 double* q_cfl,
					 double* q_numDiff_h_last,
					 double* q_numDiff_hu_last, 
					 double* q_numDiff_hv_last,
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
      void evaluateCoefficients(const double nu,
				const double g,
				const double grad_b[nSpace],
				const double& h,
				const double& hu,
				const double& hv,
				double& mass_acc,
				double& dmass_acc_h,
				double& mom_hu_acc,
				double& dmom_hu_acc_h,
				double& dmom_hu_acc_hu,
				double& mom_hv_acc,
				double& dmom_hv_acc_h,
				double& dmom_hv_acc_hv,
				double mass_adv[nSpace],
				double dmass_adv_h[nSpace],
				double dmass_adv_hu[nSpace],
				double dmass_adv_hv[nSpace],
				double mom_hu_adv[nSpace],
				double dmom_hu_adv_h[nSpace],
				double dmom_hu_adv_hu[nSpace],
				double dmom_hu_adv_hv[nSpace],
				double mom_hv_adv[nSpace],
				double dmom_hv_adv_h[nSpace],
				double dmom_hv_adv_hu[nSpace],
				double dmom_hv_adv_hv[nSpace],
				double mom_hu_diff_ten[nSpace],
				double mom_hv_diff_ten[nSpace],
				double mom_huhv_diff_ten[1],
				double mom_hvhu_diff_ten[1],
				double& mom_hu_source,
				double& dmom_hu_source_h,
				double& mom_hv_source,
				double& dmom_hv_source_h)
    {
      double hStar = fmax(1.0e-8,h);
      //mass accumulation
      mass_acc = h;
      dmass_acc_h = 1.0;
      
      //u momentum accumulation
      mom_hu_acc=hu;
      dmom_hu_acc_h=0.0;
      dmom_hu_acc_hu=1.0;
  
      //v momentum accumulation
      mom_hv_acc=hv;
      dmom_hv_acc_h=0.0;
      dmom_hv_acc_hv=1.0;
  
      //mass advective flux
      mass_adv[0]=hu;
      mass_adv[1]=hv;
  
      dmass_adv_h[0]=0.0;
      dmass_adv_h[1]=0.0;

      dmass_adv_hu[0]=1.0;
      dmass_adv_hu[1]=0.0;

      dmass_adv_hv[0]=0.0;
      dmass_adv_hv[1]=1.0;

      //u momentum advective flux
      mom_hu_adv[0]=hu*hu/hStar  + 0.5*g*h*h;
      mom_hu_adv[1]=hu*hv/hStar;
      
      dmom_hu_adv_h[0]=-hu*hu/(hStar*hStar) + g*h;
      dmom_hu_adv_h[1]=-hu*hv/(hStar*hStar);
  
      dmom_hu_adv_hu[0]=2.0*hu/hStar;
      dmom_hu_adv_hu[1]=hv/hStar;
  
      dmom_hu_adv_hv[0]=0.0;
      dmom_hu_adv_hv[1]=hu/hStar;
  
      //v momentum advective_flux
      mom_hv_adv[0]=hv*hu/hStar;
      mom_hv_adv[1]=hv*hv/hStar + 0.5*g*h*h;
  
      dmom_hv_adv_h[0]=-hv*hu/(hStar*hStar);
      dmom_hv_adv_h[1]=-hv*hv/(hStar*hStar) + g*h;
  
      dmom_hv_adv_hu[0]=hv/hStar;
      dmom_hv_adv_hu[1]=0.0;
  
      dmom_hv_adv_hv[0]=hu/hStar;
      dmom_hv_adv_hv[1]=2.0*hv/hStar;

      //u momentum diffusion tensor
      mom_hu_diff_ten[0] = 2.0*nu;
      mom_hu_diff_ten[1] = nu;
  
      mom_huhv_diff_ten[0]=nu;
  
      //v momentum diffusion tensor
      mom_hv_diff_ten[0] = nu;
      mom_hv_diff_ten[1] = 2.0*nu;
  
      mom_hvhu_diff_ten[0]=nu;
  
      //momentum sources
      mom_hu_source = g*h*grad_b[0];
      dmom_hu_source_h = g*grad_b[0];

      mom_hv_source = g*h*grad_b[1];
      dmom_hv_source_h = g*grad_b[1];
    }

    inline
      void evaluateCoefficientsForResidual(const double g,
					   const double grad_b[nSpace],
					   const double& h,
					   const double& hu,
					   const double& hv,
					   double& mass_acc,
					   double& mom_hu_acc,
					   double& mom_hv_acc,
					   double mass_adv[nSpace],
					   double mom_hu_adv[nSpace],
					   double mom_hv_adv[nSpace],
					   double& mom_hu_source,
					   double& mom_hv_source)					   
    {
      double hStar = fmax(1.0e-8,h);
      //mass accumulation
      mass_acc = h;      
      //u momentum accumulation
      mom_hu_acc=hu;
      //v momentum accumulation
      mom_hv_acc=hv;
      //mass advective flux
      mass_adv[0]=hu;
      mass_adv[1]=hv;
      //u momentum advective flux
      mom_hu_adv[0]=hu*hu/hStar  + 0.5*g*h*h;
      mom_hu_adv[1]=hu*hv/hStar;
      //v momentum advective_flux
      mom_hv_adv[0]=hv*hu/hStar;
      mom_hv_adv[1]=hv*hv/hStar + 0.5*g*h*h;
      //momentum sources
      mom_hu_source = g*h*grad_b[0];
      mom_hv_source = g*h*grad_b[1];
    }

    inline
      void evaluateCoefficientsForJacobian(const double g,
					   const double grad_b[nSpace],
					   const double& h,
					   const double& hu,
					   const double& hv,
					   double& mass_acc,
					   double& dmass_acc_h,
					   double& mom_hu_acc,
					   double& dmom_hu_acc_h,
					   double& dmom_hu_acc_hu,
					   double& mom_hv_acc,
					   double& dmom_hv_acc_h,
					   double& dmom_hv_acc_hv,
					   double mass_adv[nSpace],
					   double dmass_adv_h[nSpace],
					   double dmass_adv_hu[nSpace],
					   double dmass_adv_hv[nSpace],
					   double mom_hu_adv[nSpace],
					   double dmom_hu_adv_h[nSpace],
					   double dmom_hu_adv_hu[nSpace],
					   double dmom_hu_adv_hv[nSpace],
					   double mom_hv_adv[nSpace],
					   double dmom_hv_adv_h[nSpace],
					   double dmom_hv_adv_hu[nSpace],
					   double dmom_hv_adv_hv[nSpace],
					   double& mom_hu_source,
					   double& dmom_hu_source_h,
					   double& mom_hv_source,
					   double& dmom_hv_source_h)
    {
      double hStar = fmax(1.0e-8,h);
      //mass accumulation
      mass_acc = h;
      dmass_acc_h = 1.0;
      
      //u momentum accumulation
      mom_hu_acc=hu;
      dmom_hu_acc_h=0.0;
      dmom_hu_acc_hu=1.0;
  
      //v momentum accumulation
      mom_hv_acc=hv;
      dmom_hv_acc_h=0.0;
      dmom_hv_acc_hv=1.0;
  
      //mass advective flux
      mass_adv[0]=hu;
      mass_adv[1]=hv;
  
      dmass_adv_h[0]=0.0;
      dmass_adv_h[1]=0.0;

      dmass_adv_hu[0]=1.0;
      dmass_adv_hu[1]=0.0;

      dmass_adv_hv[0]=0.0;
      dmass_adv_hv[1]=1.0;

      //u momentum advective flux
      mom_hu_adv[0]=hu*hu/hStar  + 0.5*g*h*h;
      mom_hu_adv[1]=hu*hv/hStar;
      
      dmom_hu_adv_h[0]=-hu*hu/(hStar*hStar) + g*h;
      dmom_hu_adv_h[1]=-hu*hv/(hStar*hStar);
  
      dmom_hu_adv_hu[0]=2.0*hu/hStar;
      dmom_hu_adv_hu[1]=hv/hStar;
  
      dmom_hu_adv_hv[0]=0.0;
      dmom_hu_adv_hv[1]=hu/hStar;
  
      //v momentum advective_flux
      mom_hv_adv[0]=hv*hu/hStar;
      mom_hv_adv[1]=hv*hv/hStar + 0.5*g*h*h;
  
      dmom_hv_adv_h[0]=-hv*hu/(hStar*hStar);
      dmom_hv_adv_h[1]=-hv*hv/(hStar*hStar) + g*h;
  
      dmom_hv_adv_hu[0]=hv/hStar;
      dmom_hv_adv_hu[1]=0.0;
  
      dmom_hv_adv_hv[0]=hu/hStar;
      dmom_hv_adv_hv[1]=2.0*hv/hStar;

      //momentum sources
      mom_hu_source = g*h*grad_b[0];
      dmom_hu_source_h = g*grad_b[0];

      mom_hv_source = g*h*grad_b[1];
      dmom_hv_source_h = g*grad_b[1];
    }

    inline
    void calculateSubgridError_tau(const double& elementDiameter,
    				   const double& nu,
    				   const double& g,
    				   const double& h,
    				   const double& hu,
    				   const double& hv,
    				   double tau[9],
    				   double& cfl)
    {
      double rx[9],ry[9],rxInv[9],ryInv[9],rn=0.0,c=sqrt(fmax(g*1.0e-8,g*h)),lambdax[3],lambday[3],tauxHat[3],tauyHat[3],taux[9],tauy[9],tauc[9],cflx,cfly,L[9],hStar=fmax(1.0e-8,h),u,v;
      u = hu/hStar;
      v = hv/hStar;
      //eigenvalues and eigenvectors for conservation variables h,hu,hv
      lambdax[0] = u - c;
      lambdax[1] = u;
      lambdax[2] = u + c;
      if (u > 0.0)
	cflx = (u+c)/elementDiameter;
      else
	cflx = fabs(u-c)/elementDiameter;

      rn = sqrt(1.0 + (u-c)*(u-c) + v*v);
      //rn = 1.0;
      rx[0*3+0] = 1.0/rn;
      rx[1*3+0] = (u - c)/rn;
      rx[2*3+0] = v/rn;

      rx[0*3+1] = 0.0;
      rx[1*3+1] = 0.0;
      rx[2*3+1] = 1.0;

      rn = sqrt(1.0 + (u+c)*(u+c) + v*v);
      //rn = 1.0;
      rx[0*3+2] = 1.0/rn;
      rx[1*3+2] = (u + c)/rn;
      rx[2*3+2] = v/rn;

      rxInv[0*3+0] = 1.0 - (c - u)/(2.0*c);
      rxInv[1*3+0] = -v;
      rxInv[2*3+0] = (c-u)/(2.0*c);

      rxInv[0*3+1] = -1.0/(2.0*c);
      rxInv[1*3+1] =  0.0;
      rxInv[2*3+1] =  1.0/(2.0*c);

      rxInv[0*3+2] = 0.0;
      rxInv[1*3+2] = 1.0;
      rxInv[2*3+2] = 0.0;

      /* //cek hack assume no viscosity for now */
      /* num = 0.5*fabs(lambdax[0])*elementDiameter + 1.0e-8; */
      /* den = nu + num*1.0e-8; */
      /* pe = num/den; */
      /* tauxHat[0] = 0.5*h*(1.0/tanh(pe) - 1.0/pe)/(Vlin+1.0e-8); */

      tauxHat[0] = 0.5*elementDiameter/(fabs(lambdax[0])+1.0e-8);
      tauxHat[1] = 0.5*elementDiameter/(fabs(lambdax[1])+1.0e-8);
      tauxHat[2] = 0.5*elementDiameter/(fabs(lambdax[2])+1.0e-8);

      lambday[0] = v - c;
      lambday[1] = v;
      lambday[2] = v + c;

      rn=sqrt(1.0 + u*u+(v-c)*(v-c));
      //rn = 1.0;
      ry[0*3+0] = 1.0/rn;
      ry[1*3+0] = u/rn;
      ry[2*3+0] = (v - c)/rn;
      
      ry[0*3+1] =  0.0;
      ry[1*3+1] = -1.0;
      ry[2*3+1] =  0.0;

      rn = sqrt(1.0 + u*u + (v+c)*(v+c));
      //rn = 1.0;
      ry[0*3+2] = 1.0/rn;
      ry[1*3+2] = u/rn;
      ry[2*3+2] = (v + c)/rn;

      ryInv[0*3+0] = 1.0 - (c - v)/(2*c);
      ryInv[1*3+0] = u;
      ryInv[2*3+0] = (c-v)/(2*c);

      ryInv[0*3+1] =  0.0;
      ryInv[1*3+1] = -1.0;
      ryInv[2*3+1] =  0.0;

      ryInv[0*3+2] = -1.0/(2*c);
      ryInv[1*3+2] =  0.0;
      ryInv[2*3+2] =  1.0/(2*c);

      //cek hack assume no viscosity for now
      tauyHat[0] = 0.5*elementDiameter/(fabs(lambday[0])+1.0e-8);
      tauyHat[1] = 0.5*elementDiameter/(fabs(lambday[1])+1.0e-8);
      tauyHat[2] = 0.5*elementDiameter/(fabs(lambday[2])+1.0e-8);
      if (v > 0.0)
	cfly = (v+c)/elementDiameter;
      else
	cfly = fabs(v-c)/elementDiameter;
      cfl = sqrt(cflx*cflx+cfly*cfly);//hack, conservative estimate

      //project back to solution variables
      //initialize, hack do with memset
      double tmpx[9],tmpy[9];
      for (int i=0;i<3;i++)
	for (int j=0;j<3;j++)
	  {
	    taux[i*3+j] = 0.0;
	    tauy[i*3+j] = 0.0;
	    tauc[i*3+j] = 0.0;
	    tau[i*3+j] = 0.0;
	    tmpx[i*3+j] = 0.0;
	    tmpy[i*3+j] = 0.0;
	    L[i*3+j] = 0.0;
	    /* double Ix=0,Iy=0.0; */
	    /* for (int k=0;k<3;k++) */
	    /*   { */
	    /* 	Ix += rx[i*3+k]*rx[j*3+k]; */
	    /* 	Iy += ry[i*3+k]*ry[j*3+k]; */
	    /*   } */
	    /* std::cout<<i<<'\t'<<j<<'\t'<<Ix<<'\t'<<Iy<<std::endl; */
	  }
      //transform from characteristic variables to conservation variables
      for (int i=0;i<3;i++)
	for (int j=0;j<3;j++)
	  for (int m=0;m<3;m++)
	  {
	    if (m==j)
	      {
		tmpx[i*3+m] += rx[i*3+m]*tauxHat[m];
		tmpy[i*3+m] += ry[i*3+m]*tauyHat[m];
	      }
	  }
      for (int i=0;i<3;i++)
	for (int j=0;j<3;j++)
	  for (int m=0;m<3;m++)
	  {
	    taux[i*3+j] += tmpx[i*3+m]*rx[j*3 + m];
	    tauy[i*3+j] += tmpy[i*3+m]*ry[j*3 + m];
	  }
      //matrix norm
      for (int i=0;i<3;i++)
	for (int j=0;j<3;j++)
	  {
	    tauc[i*3+j] = sqrt(taux[i*3+j]*taux[i*3+j] + tauy[i*3+j]*tauy[i*3+j]);
	    tau[i*3+j] = tauc[i*3+j];
	    tau[i*3+j] = 0.0;//hack
	  }
      /* std::cout<<"tau"<<std::endl; */
      /* for (int i=0;i<3;i++) */
      /* 	{ */
      /* 	  for (int j=0;j<3;j++) */
      /* 	    { */
      /* 	      std::cout<<tau[i*3+j]<<'\t'; */
      /* 	    } */
      /* 	  std::cout<<std::endl; */
      /* 	} */
    }


    inline 
      double maxWaveSpeedTwoRarefactions(double g, double nx, double ny,
					 double hL, double huL, double hvL, 
					 double hR, double huR, double hvR) 
    {	
      //1-eigenvalue: uL-sqrt(g*hL)
      //3-eigenvalue: uR+sqrt(g*hR) 

      double hVelL = nx*huL + ny*hvL;
      double hVelR = nx*huR + ny*hvR;
      double velL = hVelL/hL;
      double velR = hVelR/hR;

      // Start computing lambda1 and lambda3 as if we have a 1- and 3-rarefactions 
      double lambda1 = velL - sqrt(g*hL);
      double lambda3 = velR + sqrt(g*hR);
      
      return fmax(fabs(lambda1),fabs(lambda3));
    }    

    inline 
      double maxWaveSpeedSharpInitialGuess(double g, double nx, double ny,
					   double hL, double huL, double hvL, 
					   double hR, double huR, double hvR, 
					   double hEps, bool debugging) 
    {
      double lambda1, lambda3;
      //1-eigenvalue: uL-sqrt(g*hL)
      //3-eigenvalue: uR+sqrt(g*hR) 
      
      double hVelL = nx*huL + ny*hvL;
      double hVelR = nx*huR + ny*hvR;
      double velL = 2*hL/(hL*hL+std::pow(fmax(hL,hEps),2))*hVelL;
      double velR = 2*hR/(hR*hR+std::pow(fmax(hR,hEps),2))*hVelR;

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
      return fmax(fmax(0.,-lambda1), fmax(0,lambda3));
    }
    
    inline 
      double maxWaveSpeedIterativeProcess(double g, double nx, double ny,
					  double hL, double huL, double hvL, 
					  double hR, double huR, double hvR, 
					  double hEps, bool verbose) 
    {
      double tol = 1E-15;
      //1-eigenvalue: uL-sqrt(g*hL)
      //3-eigenvalue: uR+sqrt(g*hR) 
      
      double hVelL = nx*huL + ny*hvL;
      double hVelR = nx*huR + ny*hvR;
      double velL = 2*hL/(hL*hL+std::pow(fmax(hL,hEps),2))*hVelL;
      double velR = 2*hR/(hR*hR+std::pow(fmax(hR,hEps),2))*hVelR;

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

    inline
      void exteriorNumericalAdvectiveFlux(const int& isDOFBoundary_h,
					  const int& isDOFBoundary_hu,
					  const int& isDOFBoundary_hv,
					  const int& isFluxBoundary_h,
					  const int& isFluxBoundary_hu,
					  const int& isFluxBoundary_hv,
					  const double n[nSpace],
					  const double& bc_h,
					  const double bc_f_mass[nSpace],
					  const double bc_f_humom[nSpace],
					  const double bc_f_hvmom[nSpace],
					  const double& bc_flux_mass,
					  const double& bc_flux_humom,
					  const double& bc_flux_hvmom,
					  const double& h,
					  const double f_mass[nSpace],
					  const double f_humom[nSpace],
					  const double f_hvmom[nSpace],
					  const double df_mass_dh[nSpace],
					  const double df_mass_du[nSpace],
					  const double df_mass_dv[nSpace],
					  const double df_humom_dh[nSpace],
					  const double df_humom_du[nSpace],
					  const double df_humom_dv[nSpace],
					  const double df_hvmom_dh[nSpace],
					  const double df_hvmom_du[nSpace],
					  const double df_hvmom_dv[nSpace],
					  double& flux_mass,
					  double& flux_humom,
					  double& flux_hvmom,
					  double* velocity)
    {
      //cek todo, need to do the Riemann solve
      /* double flowDirection; */
      flux_mass = 0.0;
      flux_humom = 0.0;
      flux_hvmom = 0.0;
      /* flowDirection=n[0]*f_mass[0]+n[1]*f_mass[1]; */
      /* if (isDOFBoundary_hu != 1) */
      /* 	{ */
      /* 	  flux_mass += n[0]*f_mass[0]; */
      /* 	  velocity[0] = f_mass[0]; */
      /* 	  if (flowDirection >= 0.0) */
      /* 	    { */
      /* 	      flux_humom += n[0]*f_humom[0]; */
      /* 	      flux_hvmom += n[0]*f_hvmom[0]; */
      /* 	    } */
      /* 	} */
      /* else */
      /* 	{ */
      /* 	  flux_mass += n[0]*bc_f_mass[0]; */
      /* 	  velocity[0] = bc_f_mass[0]; */
      /* 	  //cek still upwind the advection for Dirichlet? */
      /* 	  if (flowDirection >= 0.0) */
      /* 	    { */
      /* 	      flux_humom += n[0]*f_humom[0]; */
      /* 	      flux_hvmom += n[0]*f_hvmom[0]; */
      /* 	    } */
      /* 	  else */
      /* 	    { */
      /* 	      flux_humom+=n[0]*bc_f_humom[0]; */
      /* 	      flux_hvmom+=n[0]*bc_f_hvmom[0]; */
      /* 	    } */
      /* 	} */
      /* if (isDOFBoundary_hv != 1) */
      /* 	{ */
      /* 	  flux_mass+=n[1]*f_mass[1]; */
      /* 	  velocity[1] = f_mass[1]; */
      /* 	  if (flowDirection >= 0.0) */
      /* 	    { */
      /* 	      flux_humom+=n[1]*f_humom[1]; */
      /* 	      flux_hvmom+=n[1]*f_hvmom[1]; */
      /* 	    } */
      /* 	} */
      /* else */
      /* 	{ */
      /* 	  flux_mass+=n[1]*bc_f_mass[1]; */
      /* 	  velocity[1] = bc_f_mass[1]; */
      /* 	  //cek still upwind the advection for Dirichlet? */
      /* 	  if (flowDirection >= 0.0) */
      /* 	    { */
      /* 	      flux_humom+=n[1]*f_humom[1]; */
      /* 	      flux_hvmom+=n[1]*f_hvmom[1]; */
      /* 	    } */
      /* 	  else */
      /* 	    { */
      /* 	      flux_humom+=n[1]*bc_f_humom[1]; */
      /* 	      flux_hvmom+=n[1]*bc_f_hvmom[1]; */
      /* 	    } */
      /* 	} */
      /* else */
      /* 	{ */
      /* 	  flux_mass +=n[2]*bc_f_mass[2]; */
      /* 	  velocity[2] = bc_f_mass[2]; */
      /* 	  //cek still upwind the advection for Dirichlet? */
      /* 	  if (flowDirection >= 0.0) */
      /* 	    { */
      /* 	      flux_humom+=n[2]*f_humom[2]; */
      /* 	      flux_hvmom+=n[2]*f_hvmom[2]; */
      /* 	    } */
      /* 	  else */
      /* 	    { */
      /* 	      flux_humom+=n[2]*bc_f_humom[2]; */
      /* 	      flux_hvmom+=n[2]*bc_f_hvmom[2]; */
      /* 	    } */
      /* 	} */
      /* if (isDOFBoundary_h == 1) */
      /* 	{ */
      /* 	  flux_humom+= n[0]*(bc_h-p)*oneByRho; */
      /* 	  flux_hvmom+= n[1]*(bc_h-p)*oneByRho; */
      /* 	} */
      if (isFluxBoundary_h == 1)
	{
	  //cek todo, not sure if we'll need this for SW2
	  //velocity[0] += (bc_flux_mass - flux_mass)*n[0];
	  //velocity[1] += (bc_flux_mass - flux_mass)*n[1];
	  //velocity[2] += (bc_flux_mass - flux_mass)*n[2];
	  flux_mass = bc_flux_mass;
	}
      if (isFluxBoundary_hu == 1)
	{
	  flux_humom = bc_flux_humom;
	}
      if (isFluxBoundary_hv == 1)
	{
	  flux_hvmom = bc_flux_hvmom;
	}
    }

    inline
    void exteriorNumericalAdvectiveFluxDerivatives(const int& isDOFBoundary_h,
						   const int& isDOFBoundary_hu,
						   const int& isDOFBoundary_hv,
						   const int& isFluxBoundary_h,
						   const int& isFluxBoundary_hu,
						   const int& isFluxBoundary_hv,
						   const double n[nSpace],
						   const double& bc_h,
						   const double bc_f_mass[nSpace],
						   const double bc_f_humom[nSpace],
						   const double bc_f_hvmom[nSpace],
						   const double& bc_flux_mass,
						   const double& bc_flux_humom,
						   const double& bc_flux_hvmom,
						   const double& h,
						   const double f_mass[nSpace],
						   const double f_humom[nSpace],
						   const double f_hvmom[nSpace],
						   const double df_mass_du[nSpace],
						   const double df_mass_dv[nSpace],
						   const double df_humom_dh[nSpace],
						   const double df_humom_du[nSpace],
						   const double df_humom_dv[nSpace],
						   const double df_hvmom_dh[nSpace],
						   const double df_hvmom_du[nSpace],
						   const double df_hvmom_dv[nSpace],
						   double& dflux_mass_dh,
						   double& dflux_mass_du,
						   double& dflux_mass_dv,
						   double& dflux_humom_dh,
						   double& dflux_humom_du,
						   double& dflux_humom_dv,
						   double& dflux_hvmom_dh,
						   double& dflux_hvmom_du,
						   double& dflux_hvmom_dv)
    {
      double flowDirection;
      dflux_mass_dh = 0.0;
      dflux_mass_du = 0.0;
      dflux_mass_dv = 0.0;
  
      dflux_humom_dh = 0.0;
      dflux_humom_du = 0.0;
      dflux_humom_dv = 0.0;
  
      dflux_hvmom_dh = 0.0;
      dflux_hvmom_du = 0.0;
      dflux_hvmom_dv = 0.0;
  
      flowDirection=n[0]*f_mass[0]+n[1]*f_mass[1];
      /* if (isDOFBoundary_hu != 1) */
      /* 	{ */
      /* 	  dflux_mass_du += n[0]*df_mass_du[0]; */
      /* 	  if (flowDirection >= 0.0) */
      /* 	    { */
      /* 	      dflux_humom_du += n[0]*df_humom_du[0]; */
      /* 	      dflux_hvmom_du += n[0]*df_hvmom_du[0]; */
      /* 	      dflux_hvmom_dv += n[0]*df_hvmom_dv[0]; */
      /* 	    } */
      /* 	} */
      /* else */
      /* 	{ */
      /* 	  //cek still upwind the advection for Dirichlet? */
      /* 	  if (flowDirection >= 0.0) */
      /* 	    { */
      /* 	      dflux_humom_du += n[0]*df_humom_du[0]; */
      /* 	      dflux_hvmom_du += n[0]*df_hvmom_du[0]; */
      /* 	      dflux_hvmom_dv += n[0]*df_hvmom_dv[0]; */
      /* 	    } */
      /* 	  else */
      /* 	    { */
      /* 	      if (isDOFBoundary_hv != 1) */
      /* 		dflux_hvmom_dv += n[0]*df_hvmom_dv[0]; */
      /* 	    } */
      /* 	} */
      /* if (isDOFBoundary_hv != 1) */
      /* 	{ */
      /* 	  dflux_mass_dv += n[1]*df_mass_dv[1]; */
      /* 	  if (flowDirection >= 0.0) */
      /* 	    { */
      /* 	      dflux_humom_du += n[1]*df_humom_du[1]; */
      /* 	      dflux_humom_dv += n[1]*df_humom_dv[1]; */
      /* 	      dflux_hvmom_dv += n[1]*df_hvmom_dv[1]; */
      /* 	    } */
      /* 	} */
      /* else */
      /* 	{ */
      /* 	  //cek still upwind the advection for Dirichlet? */
      /* 	  if (flowDirection >= 0.0) */
      /* 	    { */
      /* 	      dflux_humom_du += n[1]*df_humom_du[1]; */
      /* 	      dflux_humom_dv += n[1]*df_humom_dv[1]; */
      /* 	      dflux_hvmom_dv += n[1]*df_hvmom_dv[1]; */
      /* 	    } */
      /* 	  else */
      /* 	    { */
      /* 	      if (isDOFBoundary_hu != 1) */
      /* 		dflux_humom_du += n[1]*df_humom_du[1]; */
      /* 	    } */
      /* 	} */
      /* else */
      /* 	{ */
      /* 	  //cek still upwind the advection for Dirichlet? */
      /* 	  if (flowDirection >= 0.0) */
      /* 	    { */
      /* 	      dflux_humom_du += n[2]*df_humom_du[2]; */
      /* 	      dflux_humom_dw += n[2]*df_humom_dw[2]; */
      /* 	      dflux_hvmom_dv += n[2]*df_hvmom_dv[2]; */
      /* 	    } */
      /* 	  else */
      /* 	    { */
      /* 	      if (isDOFBoundary_hu != 1) */
      /* 		dflux_humom_du += n[2]*df_humom_du[2]; */
      /* 	      if (isDOFBoundary_hv != 1) */
      /* 		dflux_hvmom_dv += n[2]*df_hvmom_dv[2]; */
      /* 	    } */
      /* 	} */
      /* if (isDOFBoundary_h == 1) */
      /* 	{ */
      /* 	  dflux_humom_dp= -n[0]*oneByRho; */
      /* 	  dflux_hvmom_dp= -n[1]*oneByRho; */
      /* 	} */
      if (isFluxBoundary_h == 1)
	{
	  dflux_mass_dh = 0.0;
	  dflux_mass_du = 0.0;
	  dflux_mass_dv = 0.0;
	}
      if (isFluxBoundary_hu == 1)
	{
	  dflux_humom_dh = 0.0;
	  dflux_humom_du = 0.0;
	  dflux_humom_dv = 0.0;
	}
      if (isFluxBoundary_hv == 1)
	{
	  dflux_hvmom_dh = 0.0;
	  dflux_hvmom_du = 0.0;
	  dflux_hvmom_dv = 0.0;
	}
    }

    /* inline */
    /* void exteriorNumericalDiffusiveFlux(const double& eps, */
    /* 					int* rowptr, */
    /* 					int* colind, */
    /* 					const int& isDOFBoundary, */
    /* 					const int& isFluxBoundary, */
    /* 					const double n[nSpace], */
    /* 					double* bc_a, */
    /* 					const double& bc_hu, */
    /* 					const double& bc_flux, */
    /* 					double* a, */
    /* 					const double grad_phi[nSpace], */
    /* 					const double& u, */
    /* 					const double& penalty, */
    /* 					double& flux) */
    /* { */
    /*   double diffusiveVelocityComponent_I,penaltyFlux,max_a; */
    /*   if(isDOFBoundary == 1) */
    /* 	{ */
    /* 	  flux = 0.0; */
    /* 	  max_a=0.0; */
    /* 	  for(int I=0;I<nSpace;I++) */
    /* 	    { */
    /* 	      diffusiveVelocityComponent_I=0.0; */
    /* 	      for(int m=rowptr[I];m<rowptr[I+1];m++) */
    /* 		{ */
    /* 		  diffusiveVelocityComponent_I -= a[m]*grad_phi[colind[m]]; */
    /* 		  max_a = fmax(max_a,a[m]); */
    /* 		} */
    /* 	      flux+= diffusiveVelocityComponent_I*n[I]; */
    /* 	    } */
    /* 	  penaltyFlux = max_a*penalty*(u-bc_hu); */
    /* 	  flux += penaltyFlux; */
    /* 	} */
    /*   else if(isFluxBoundary == 1) */
    /* 	{ */
    /* 	  flux = bc_flux; */
    /* 	} */
    /*   else */
    /* 	{ */
    /* 	  std::cerr<<"warning, diffusion term with no boundary condition set, setting diffusive flux to 0.0"<<std::endl; */
    /* 	  flux = 0.0; */
    /* 	} */
    /* } */

    /* inline */
    /* double ExteriorNumericalDiffusiveFluxJacobian(const double& eps, */
    /* 						  int* rowptr, */
    /* 						  int* colind, */
    /* 						  const int& isDOFBoundary, */
    /* 						  const double n[nSpace], */
    /* 						  double* a, */
    /* 						  const double& v, */
    /* 						  const double grad_hv[nSpace], */
    /* 						  const double& penalty) */
    /* { */
    /*   double dvel_I,tmp=0.0,max_a=0.0; */
    /*   if(isDOFBoundary >= 1) */
    /* 	{ */
    /* 	  for(int I=0;I<nSpace;I++) */
    /* 	    { */
    /* 	      dvel_I=0.0; */
    /* 	      for(int m=rowptr[I];m<rowptr[I+1];m++) */
    /* 		{ */
    /* 		  dvel_I -= a[m]*grad_hv[colind[m]]; */
    /* 		  max_a = fmax(max_a,a[m]); */
    /* 		} */
    /* 	      tmp += dvel_I*n[I]; */
    /* 	    } */
    /* 	  tmp +=max_a*penalty*v; */
    /* 	} */
    /*   return tmp; */
    /* } */

    void FCTStep(double dt, 
		 int NNZ, //number on non-zero entries on sparsity pattern
		 int numDOFs, //number of DOFs
		 double* lumped_mass_matrix, //lumped mass matrix (as vector))
		 double* h_lstage, //DOFs of solution at last stage
		 double* hu_lstage, 
		 double* hv_lstage,
		 double* b_dof,
		 double* high_order_hnp1, //DOFs of high order solution at tnp1
		 double* high_order_hunp1, 
		 double* high_order_hvnp1, 
		 double* low_order_hnp1, //operators to construct low order solution
		 double* low_order_hunp1,
		 double* low_order_hvnp1,
		 double* limited_hnp1,
		 double* limited_hunp1,
		 double* limited_hvnp1,
		 int* csrRowIndeces_DofLoops, //csr row indeces 
		 int* csrColumnOffsets_DofLoops, //csr column offsets 
		 double* MassMatrix, //mass matrix
		 double* dEV_minus_dL,
		 double hEps, 
		 int LUMPED_MASS_MATRIX
		 )
    {
      double h_threshold = 1.0E-3;
      register double Rneg[numDOFs], Rpos[numDOFs];
      //////////////////
      // LOOP in DOFs //
      //////////////////
      int ij=0;
      for (int i=0; i<numDOFs; i++)
	{
	  //read some vectors 
	  double high_order_hnp1i  = high_order_hnp1[i];
	  double hni = h_lstage[i];
	  double Zi = b_dof[i];
	  double mi = lumped_mass_matrix[i];

	  double minH=0.;
	  double Pnegi=0.;
	  // LOOP OVER THE SPARSITY PATTERN (j-LOOP)//
	  for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
	    {
	      int j = csrColumnOffsets_DofLoops[offset];
	      // read some vectors
	      double hnj = h_lstage[j];
	      double Zj = b_dof[j];

	      // COMPUTE STAR SOLUTION // hStar, huStar and hvStar
	      double hStarij  = fmax(0., hni + Zi - fmax(Zi,Zj));
	      double hStarji  = fmax(0., hnj + Zj - fmax(Zi,Zj));
	      
	      // i-th row of flux correction matrix 
	      double ML_minus_MC = (LUMPED_MASS_MATRIX == 1 ? 0. : 
				    (i==j ? 1. : 0.)*mi - MassMatrix[ij]);
	      double FluxCorrectionMatrix1 
		= ML_minus_MC*(high_order_hnp1[j]-hnj - (high_order_hnp1i-hni)) 
		+ dt*dEV_minus_dL[ij]*(hStarji-hStarij);

	      // COMPUTE P VECTORS //
	      Pnegi += FluxCorrectionMatrix1*((FluxCorrectionMatrix1 < 0) ? 1. : 0.);	      

	      //update ij 
	      ij+=1;
	    }
	  ///////////////////////
	  // COMPUTE Q VECTORS //
	  ///////////////////////
	  double Qnegi = mi*(minH-low_order_hnp1[i]);

	  ///////////////////////
	  // COMPUTE R VECTORS //
	  ///////////////////////
	  Rpos[i] = 1.;
	  if (high_order_hnp1[i] < h_threshold)
	    Rneg[i] = 0.;
	  else
	    Rneg[i] = ((Pnegi==0) ? 1. : std::min(1.0,Qnegi/Pnegi));
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
	  double hni = h_lstage[i];
	  double huni = hu_lstage[i];
	  double hvni = hv_lstage[i];
	  double Zi = b_dof[i];
	  double mi = lumped_mass_matrix[i];

	  double ith_Limiter_times_FluxCorrectionMatrix1 = 0.;
	  double ith_Limiter_times_FluxCorrectionMatrix2 = 0.;
	  double ith_Limiter_times_FluxCorrectionMatrix3 = 0.;
	  double Rposi = Rpos[i], Rnegi = Rneg[i];
	  // LOOP OVER THE SPARSITY PATTERN (j-LOOP)//
	  for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
	    {
	      int j = csrColumnOffsets_DofLoops[offset];
	      // read some vectors
	      double hnj = h_lstage[j];
	      double hunj = hu_lstage[j];
	      double hvnj = hv_lstage[j];
	      double Zj = b_dof[j];

	      // COMPUTE STAR SOLUTION // hStar, huStar and hvStar
	      double hStarij  = fmax(0., hni + Zi - fmax(Zi,Zj));
	      double huStarij = (hni <= hEps ? 0. : huni*hStarij/hni);
	      double hvStarij = (hni <= hEps ? 0. : hvni*hStarij/hni);
	      
	      double hStarji  = fmax(0., hnj + Zj - fmax(Zi,Zj));
	      double huStarji = (hnj <= hEps ? 0. : hunj*hStarji/hnj);
	      double hvStarji = (hnj <= hEps ? 0. : hvnj*hStarji/hnj);

	      // COMPUTE FLUX CORRECTION MATRICES
	      double ML_minus_MC = (LUMPED_MASS_MATRIX == 1 ? 0. : 
				    (i==j ? 1. : 0.)*mi - MassMatrix[ij]);
	      double FluxCorrectionMatrix1 = 
		ML_minus_MC*(high_order_hnp1[j]-hnj - (high_order_hnp1i-hni)) 
		+ dt*dEV_minus_dL[ij]*(hStarji-hStarij);

	      double FluxCorrectionMatrix2 = 
		ML_minus_MC*(high_order_hunp1[j]-hunj - (high_order_hunp1i-huni)) 
		+ dt*dEV_minus_dL[ij]*(huStarji-huStarij);

	      double FluxCorrectionMatrix3 = 
		ML_minus_MC*(high_order_hvnp1[j]-hvnj - (high_order_hvnp1i-hvni)) 
		+ dt*dEV_minus_dL[ij]*(hvStarji-hvStarij);

	      // compute limiter based on water height
	      double Lij = (FluxCorrectionMatrix1 > 0. ? std::min(Rposi,Rneg[j]) : std::min(Rnegi,Rpos[j]));

	      ith_Limiter_times_FluxCorrectionMatrix1 += Lij*FluxCorrectionMatrix1;
	      ith_Limiter_times_FluxCorrectionMatrix2 += Lij*FluxCorrectionMatrix2;
	      ith_Limiter_times_FluxCorrectionMatrix3 += Lij*FluxCorrectionMatrix3;
	      //update ij
	      ij+=1;
	    }

	  double one_over_mi = 1.0/lumped_mass_matrix[i];
	  limited_hnp1[i]  = low_order_hnp1[i]  + 0*one_over_mi*ith_Limiter_times_FluxCorrectionMatrix1;
	  limited_hunp1[i] = low_order_hunp1[i] + 0*one_over_mi*ith_Limiter_times_FluxCorrectionMatrix2;
	  limited_hvnp1[i] = low_order_hvnp1[i] + 0*one_over_mi*ith_Limiter_times_FluxCorrectionMatrix3;

	  if (limited_hnp1[i] < -1E-14)
	    {
	      std::cout << "Water height: "
			<<  limited_hnp1[i] 
			<< " ... aborting!" << std::endl;
	      abort();
	    }
	  else 
	    limited_hnp1[i] = fmax(0.,limited_hnp1[i]);

	}
    }

    void calculateResidual_SUPG(//element
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
			   double* h_dof_old_old, 
			   double* hu_dof_old_old, 
			   double* hv_dof_old_old, 
			   double* h_dof_lstage, 
			   double* hu_dof_lstage, 
			   double* hv_dof_lstage, 
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
			   double* q_velocity_sge,
			   double* q_cfl,
			   double* q_numDiff_h, 
			   double* q_numDiff_hu, 
			   double* q_numDiff_hv, 
			   double* q_numDiff_h_last,
			   double* q_numDiff_hu_last, 
			   double* q_numDiff_hv_last,
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
			   double* edge_based_cfl,
			   double cfl_run,
			   double hEps,
			   int recompute_lumped_mass_matrix,
			   // SAVE SOLUTION (mql)
			   double* hnp1_at_quad_point,
			   double* hunp1_at_quad_point,
			   double* hvnp1_at_quad_point,
			   // GALERKIN SOLUTION 
			   double* h_dof_galerkin,
			   double* hu_dof_galerkin,
			   double* hv_dof_galerkin,
			   // TO COMPUTE LOW ORDER 
			   double* low_order_hnp1,
			   double* low_order_hunp1,
			   double* low_order_hvnp1,
			   double* dEV_minus_dL,
			   double cE,
			   int LUMPED_MASS_MATRIX,
			   int USE_EV_BASED_ON_GALERKIN,
			   double dt, 
			   double mannings)
			   
    {
      //
      //loop over elements to compute volume integrals and load them into element and global residual
      //
      double globalConservationError=0.0,tauSum=0.0;
      for(int eN=0;eN<nElements_global;eN++)
      	{
      	  //declare local storage for element residual and initialize
      	  register double elementResidual_h[nDOF_test_element],
      	    elementResidual_hu[nDOF_test_element],
      	    elementResidual_hv[nDOF_test_element];
      	  for (int i=0;i<nDOF_test_element;i++)
      	    {
      	      int eN_i = eN*nDOF_test_element+i;
      	      elementResidual_h_save[eN_i]=0.0;
      	      elementResidual_h[i]=0.0;
      	      elementResidual_hu[i]=0.0;
      	      elementResidual_hv[i]=0.0;
      	    }//i
      	  //
      	  //loop over quadrature points and compute integrands
      	  //
      	  for(int k=0;k<nQuadraturePoints_element;k++)
      	    {
      	      //compute indices and declare local storage
      	      register int eN_k = eN*nQuadraturePoints_element+k,
      		eN_k_nSpace = eN_k*nSpace,
      		eN_nDOF_trial_element = eN*nDOF_trial_element;
      	      register double b=0.0,h=0.0,hu=0.0,hv=0.0,h_sge=0.0,hu_sge=0.0,hv_sge=0.0,
      		grad_b[nSpace],grad_h[nSpace],grad_hu[nSpace],grad_hv[nSpace],
      		mass_acc=0.0,
      		dmass_acc_h=0.0,
      		mom_hu_acc=0.0,
      		dmom_hu_acc_h=0.0,
      		dmom_hu_acc_hu=0.0,
      		mom_hv_acc=0.0,
      		dmom_hv_acc_h=0.0,
      		dmom_hv_acc_hv=0.0,
      		mass_adv[nSpace],
      		dmass_adv_h[nSpace],
      		dmass_adv_hu[nSpace],
      		dmass_adv_hv[nSpace],
      		dmass_adv_h_sge[nSpace],
      		dmass_adv_hu_sge[nSpace],
      		dmass_adv_hv_sge[nSpace],
      		mom_hu_adv[nSpace],
      		dmom_hu_adv_h[nSpace],
      		dmom_hu_adv_hu[nSpace],
      		dmom_hu_adv_hv[nSpace],
      		dmom_hu_adv_h_sge[nSpace],
      		dmom_hu_adv_hu_sge[nSpace],
      		dmom_hu_adv_hv_sge[nSpace],
      		mom_hv_adv[nSpace],
      		dmom_hv_adv_h[nSpace],
      		dmom_hv_adv_hu[nSpace],
      		dmom_hv_adv_hv[nSpace],
      		dmom_hv_adv_h_sge[nSpace],
      		dmom_hv_adv_hu_sge[nSpace],
      		dmom_hv_adv_hv_sge[nSpace],
      		mom_hu_diff_ten[nSpace],
      		mom_hv_diff_ten[nSpace],
      		mom_huhv_diff_ten[1],
      		mom_hvhu_diff_ten[1],
      		mom_hu_source=0.0,
      		dmom_hu_source_h=0.0,
      		mom_hv_source=0.0,
      		dmom_hv_source_h=0.0,
      		mass_acc_t=0.0,
      		dmass_acc_h_t=0.0,
      		mom_hu_acc_t=0.0,
      		dmom_hu_acc_h_t=0.0,
      		dmom_hu_acc_hu_t=0.0,
      		mom_hv_acc_t=0.0,
      		dmom_hv_acc_h_t=0.0,
      		dmom_hv_acc_hv_t=0.0,
		tau[9],
      		pdeResidual_h=0.0,
      		pdeResidual_hu=0.0,
      		pdeResidual_hv=0.0,
      		Lstar_h_h[nDOF_test_element],
      		Lstar_hu_h[nDOF_test_element],
      		Lstar_hv_h[nDOF_test_element],
      		Lstar_h_hu[nDOF_test_element],
      		Lstar_hu_hu[nDOF_test_element],
      		Lstar_hv_hu[nDOF_test_element],
      		Lstar_h_hv[nDOF_test_element],
      		Lstar_hu_hv[nDOF_test_element],
      		Lstar_hv_hv[nDOF_test_element],
      		subgridError_h=0.0,
      		subgridError_hu=0.0,
      		subgridError_hv=0.0,
      		jac[nSpace*nSpace],
      		jacDet,
      		jacInv[nSpace*nSpace],
      		h_grad_trial[nDOF_trial_element*nSpace],vel_grad_trial[nDOF_trial_element*nSpace],
      		h_test_dV[nDOF_trial_element],vel_test_dV[nDOF_trial_element],
      		h_grad_test_dV[nDOF_test_element*nSpace],vel_grad_test_dV[nDOF_test_element*nSpace],
      		dV,x,y,xt,yt,
      		G[nSpace*nSpace],G_dd_G,tr_G,norm_Rv, dmom_adv_star[nSpace],dmom_adv_sge[nSpace];
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
      	      ck.calculateMappingVelocity_element(eN,
      						  k,
      						  mesh_velocity_dof,
      						  mesh_l2g,
      						  mesh_trial_ref,
      						  xt,yt);
      	      //xt=0.0;yt=0.0;
      	      //std::cout<<"xt "<<xt<<'\t'<<yt<<'\t'<<std::endl;
      	      //get the physical integration weight
      	      dV = fabs(jacDet)*dV_ref[k];
      	      //get the trial function gradients
      	      ck.gradTrialFromRef(&h_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,h_grad_trial);
      	      ck.gradTrialFromRef(&vel_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,vel_grad_trial);
      	      //get the solution
      	      ck.valFromDOF(b_dof,&h_l2g[eN_nDOF_trial_element],&h_trial_ref[k*nDOF_trial_element],b);
      	      ck.valFromDOF(h_dof,&h_l2g[eN_nDOF_trial_element],&h_trial_ref[k*nDOF_trial_element],h);
      	      ck.valFromDOF(hu_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hu);
      	      ck.valFromDOF(hv_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hv);
      	      ck.valFromDOF(h_dof_sge,&h_l2g[eN_nDOF_trial_element],&h_trial_ref[k*nDOF_trial_element],h_sge);
      	      ck.valFromDOF(hu_dof_sge,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hu_sge);
      	      ck.valFromDOF(hv_dof_sge,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hv_sge);
      	      //get the solution gradients
      	      ck.gradFromDOF(b_dof,&h_l2g[eN_nDOF_trial_element],h_grad_trial,grad_b);
      	      ck.gradFromDOF(h_dof,&h_l2g[eN_nDOF_trial_element],h_grad_trial,grad_h);
      	      ck.gradFromDOF(hu_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_hu);
      	      ck.gradFromDOF(hv_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_hv);
      	      //precalculate test function products with integration weights
      	      for (int j=0;j<nDOF_trial_element;j++)
      		{
      		  h_test_dV[j] = h_test_ref[k*nDOF_trial_element+j]*dV;
      		  vel_test_dV[j] = vel_test_ref[k*nDOF_trial_element+j]*dV;
      		  for (int I=0;I<nSpace;I++)
      		    {
      		      h_grad_test_dV[j*nSpace+I]   = h_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
      		      vel_grad_test_dV[j*nSpace+I] = vel_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
      		    }
      		}
      	      //save velocity at quadrature points for other models to use
      	      q_velocity[eN_k_nSpace+0]=hu/h;
      	      q_velocity[eN_k_nSpace+1]=hv/h;
      	      //
      	      //calculate pde coefficients at quadrature points
      	      //
      	      evaluateCoefficients(nu,
      				   g,
      				   grad_b,
      				   h,
      				   hu,
      				   hv,
      				   mass_acc,
      				   dmass_acc_h,
      				   mom_hu_acc,
      				   dmom_hu_acc_h,
      				   dmom_hu_acc_hu,
      				   mom_hv_acc,
      				   dmom_hv_acc_h,
      				   dmom_hv_acc_hv,
      				   mass_adv,
      				   dmass_adv_h,
      				   dmass_adv_hu,
      				   dmass_adv_hv,
      				   mom_hu_adv,
      				   dmom_hu_adv_h,
      				   dmom_hu_adv_hu,
      				   dmom_hu_adv_hv,
      				   mom_hv_adv,
      				   dmom_hv_adv_h,
      				   dmom_hv_adv_hu,
      				   dmom_hv_adv_hv,
      				   mom_hu_diff_ten,
      				   mom_hv_diff_ten,
      				   mom_huhv_diff_ten,
      				   mom_hvhu_diff_ten,
      				   mom_hu_source,
      				   dmom_hu_source_h,
      				   mom_hv_source,
      				   dmom_hv_source_h);
      	      //
      	      //save momentum for time history and velocity for subgrid error
      	      //
      	      q_mass_acc[eN_k] = mass_acc;
      	      q_mom_hu_acc[eN_k] = mom_hu_acc;
      	      q_mom_hv_acc[eN_k] = mom_hv_acc;
      	      //subgrid error uses grid scale discharge
      	      q_mass_adv[eN_k_nSpace+0] = hu;
      	      q_mass_adv[eN_k_nSpace+1] = hv;
      	      //
      	      //moving mesh
      	      //
      	      //transform the continuity equation as if the accumulation term was  d(1)/dt
      	      /* mass_adv[0] -= MOVING_DOMAIN*mass_acc*xt; */
      	      /* mass_adv[1] -= MOVING_DOMAIN*mass_acc*yt; */
      	      /* dmass_adv_h[0] -= MOVING_DOMAIN*dmass_acc_h*xt; */
      	      /* dmass_adv_h[1] -= MOVING_DOMAIN*dmass_acc_h*yt; */

      	      /* mom_hu_adv[0] -= MOVING_DOMAIN*mom_hu_acc*xt; */
      	      /* mom_hu_adv[1] -= MOVING_DOMAIN*mom_hu_acc*yt; */
      	      /* dmom_hu_adv_hu[0] -= MOVING_DOMAIN*dmom_hu_acc_hu*xt; */
      	      /* dmom_hu_adv_hu[1] -= MOVING_DOMAIN*dmom_hu_acc_hu*yt; */

      	      /* mom_hv_adv[0] -= MOVING_DOMAIN*mom_hv_acc*xt; */
      	      /* mom_hv_adv[1] -= MOVING_DOMAIN*mom_hv_acc*yt; */
      	      /* dmom_hv_adv_hv[0] -= MOVING_DOMAIN*dmom_hv_acc_hv*xt; */
      	      /* dmom_hv_adv_hv[1] -= MOVING_DOMAIN*dmom_hv_acc_hv*yt; */

      	      //
      	      //calculate time derivative at quadrature points
      	      //
      	      ck.bdf(alphaBDF,
      		     q_mass_acc_beta_bdf[eN_k],
      		     mass_acc,
      		     dmass_acc_h,
      		     mass_acc_t,
      		     dmass_acc_h_t);
      	      ck.bdf(alphaBDF,
		     q_mom_hu_acc_beta_bdf[eN_k],
		     mom_hu_acc,
		     dmom_hu_acc_hu,
		     mom_hu_acc_t,
		     dmom_hu_acc_hu_t);
      	      ck.bdf(alphaBDF,
		     q_mom_hv_acc_beta_bdf[eN_k],
		     mom_hv_acc,
		     dmom_hv_acc_hv,
		     mom_hv_acc_t,
		     dmom_hv_acc_hv_t);
      	      //
      	      //calculate subgrid error (strong residual and adjoint)
      	      //
	      double hStar_sge = fmax(1.0e-8,h_sge);

      	      //calculate strong residual
	      /* dmass_adv_h_sge[0]  = dmass_adv_h[0]; */
	      /* dmass_adv_h_sge[1]  = dmass_adv_h[1]; */
	      /* dmass_adv_hu_sge[0]  = dmass_adv_hu[0]; */
	      /* dmass_adv_hu_sge[1]  = dmass_adv_hu[1]; */
	      /* dmass_adv_hv_sge[0]  = dmass_adv_hv[0]; */
	      /* dmass_adv_hv_sge[1]  = dmass_adv_hv[1]; */
	      /* dmom_hu_adv_h_sge[0] = dmom_hu_adv_h[0]; */
	      /* dmom_hu_adv_h_sge[1] = dmom_hu_adv_h[1]; */
	      /* dmom_hu_adv_hu_sge[0] = dmom_hu_adv_hu[0]; */
	      /* dmom_hu_adv_hu_sge[1] = dmom_hu_adv_hu[1]; */
	      /* dmom_hu_adv_hv_sge[0] = dmom_hu_adv_hv[0]; */
	      /* dmom_hu_adv_hv_sge[1] = dmom_hu_adv_hv[1]; */
	      /* dmom_hv_adv_h_sge[0] = dmom_hv_adv_h[0]; */
	      /* dmom_hv_adv_h_sge[1] = dmom_hv_adv_h[1]; */
	      /* dmom_hv_adv_hu_sge[0] = dmom_hv_adv_hu[0]; */
	      /* dmom_hv_adv_hu_sge[1] = dmom_hv_adv_hu[1]; */
	      /* dmom_hv_adv_hv_sge[0] = dmom_hv_adv_hv[0]; */
	      /* dmom_hv_adv_hv_sge[1] = dmom_hv_adv_hv[1]; */

	      //lagged strong residual coefficients
	      //mass advective flux
	      dmass_adv_h_sge[0]=0.0;
	      dmass_adv_h_sge[1]=0.0;
	      
	      dmass_adv_hu_sge[0]=1.0;
	      dmass_adv_hu_sge[1]=0.0;
	      
	      dmass_adv_hv_sge[0]=0.0;
	      dmass_adv_hv_sge[1]=1.0;
	      
	      //u momentum advective flux
	      dmom_hu_adv_h_sge[0]=-hu_sge*hu_sge/(hStar_sge*hStar_sge) + g*h_sge;
	      dmom_hu_adv_h_sge[1]=-hu_sge*hv_sge/(hStar_sge*hStar_sge);
	      
	      dmom_hu_adv_hu_sge[0]=2.0*hu_sge/hStar_sge;
	      dmom_hu_adv_hu_sge[1]=hv_sge/hStar_sge;
	      
	      dmom_hu_adv_hv_sge[0]=0.0;
	      dmom_hu_adv_hv_sge[1]=hu_sge/hStar_sge;
	      
	      //v momentum advective_flux
	      dmom_hv_adv_h_sge[0]=-hv_sge*hu_sge/(hStar_sge*hStar_sge);
	      dmom_hv_adv_h_sge[1]=-hv_sge*hv_sge/(hStar_sge*hStar_sge) + g*h_sge;
	      
	      dmom_hv_adv_hu_sge[0]=hv_sge/hStar_sge;
	      dmom_hv_adv_hu_sge[1]=0.0;
	      
	      dmom_hv_adv_hv_sge[0]=hu_sge/hStar_sge;
	      dmom_hv_adv_hv_sge[1]=2.0*hv_sge/hStar_sge;
	      
	      //approximate linearization

	      /* //mass advective flux */
	      /* dmass_adv_h_sge[0]=0.0; */
	      /* dmass_adv_h_sge[1]=0.0; */
	      
	      /* dmass_adv_hu_sge[0]=1.0; */
	      /* dmass_adv_hu_sge[1]=0.0; */
	      
	      /* dmass_adv_hv_sge[0]=0.0; */
	      /* dmass_adv_hv_sge[1]=1.0; */
	      
	      /* //u momentum advective flux */
	      /* dmom_hu_adv_h_sge[0]= g*h_sge; */
	      /* dmom_hu_adv_h_sge[1]= 0.0; */
	      
	      /* dmom_hu_adv_hu_sge[0]= hu_sge/hStar_sge; */
	      /* dmom_hu_adv_hu_sge[1]= hv_sge/hStar_sge; */
	      
	      /* dmom_hu_adv_hv_sge[0]= 0.0; */
	      /* dmom_hu_adv_hv_sge[1]= 0.0; */
  
	      /* //v momentum advective_flux */
	      /* dmom_hv_adv_h_sge[0]= 0.0; */
	      /* dmom_hv_adv_h_sge[1]= g*h_sge; */
	      
	      /* dmom_hv_adv_hu_sge[0]= 0.0; */
	      /* dmom_hv_adv_hu_sge[1]= 0.0; */
	      
	      /* dmom_hv_adv_hv_sge[0]= hu_sge/hStar_sge; */
	      /* dmom_hv_adv_hv_sge[1]= hv_sge/hStar_sge; */
	      
	      //
      	      pdeResidual_h = ck.Mass_strong(mass_acc_t) +
      	      	ck.Advection_strong(dmass_adv_hu_sge,grad_hu) +
      	      	ck.Advection_strong(dmass_adv_hv_sge,grad_hv);
	  
      	      pdeResidual_hu = ck.Mass_strong(mom_hu_acc_t) +
      	      	ck.Advection_strong(dmom_hu_adv_h_sge,grad_h) +
      	      	ck.Advection_strong(dmom_hu_adv_hu_sge,grad_hu) +
      	      	ck.Advection_strong(dmom_hu_adv_hv_sge,grad_hv) +
      	      	ck.Reaction_strong(mom_hu_source);
	  
      	      pdeResidual_hv = ck.Mass_strong(mom_hv_acc_t) +
		ck.Advection_strong(dmom_hv_adv_h_sge,grad_h) +
      	      	ck.Advection_strong(dmom_hv_adv_hu_sge,grad_hu) +
      	      	ck.Advection_strong(dmom_hv_adv_hv_sge,grad_hv) +
      	      	ck.Reaction_strong(mom_hv_source);
	  
      	      calculateSubgridError_tau(elementDiameter[eN],
					nu,
					g,
					h_sge,
					hu_sge,
					hv_sge,
					tau,
					q_cfl[eN_k]);
	      for (int i=0;i<9;i++)
		tauSum += tau[i];

	      subgridError_h  = - tau[0*3+0]*pdeResidual_h - tau[0*3+1]*pdeResidual_hu - tau[0*3+2]*pdeResidual_hv;
	      subgridError_hu = - tau[1*3+0]*pdeResidual_h - tau[1*3+1]*pdeResidual_hu - tau[1*3+2]*pdeResidual_hv;
	      subgridError_hv = - tau[2*3+0]*pdeResidual_h - tau[2*3+1]*pdeResidual_hu - tau[2*3+2]*pdeResidual_hv;
         
      	      //adjoint times the test functions
      	      for (int i=0;i<nDOF_test_element;i++)
      	      	{
      	      	  register int i_nSpace = i*nSpace;
      	      	  Lstar_h_h[i]=0.0;
      	      	  Lstar_hu_h[i]=ck.Advection_adjoint(dmass_adv_hu_sge,&h_grad_test_dV[i_nSpace]);
      	      	  Lstar_hv_h[i]=ck.Advection_adjoint(dmass_adv_hv_sge,&h_grad_test_dV[i_nSpace]);

      	      	  Lstar_h_hu[i]=ck.Advection_adjoint(dmom_hu_adv_h_sge,&vel_grad_test_dV[i_nSpace])  +
		    ck.Reaction_adjoint(dmom_hu_source_h,vel_test_dV[i]);
      	      	  Lstar_hu_hu[i]=ck.Advection_adjoint(dmom_hu_adv_hu_sge,&vel_grad_test_dV[i_nSpace]);
      	      	  Lstar_hv_hu[i]=ck.Advection_adjoint(dmom_hu_adv_hv_sge,&vel_grad_test_dV[i_nSpace]);
		  
      	      	  Lstar_h_hv[i]=ck.Advection_adjoint(dmom_hv_adv_h_sge,&vel_grad_test_dV[i_nSpace]) +
		    ck.Reaction_adjoint(dmom_hv_source_h,vel_test_dV[i]);
      	      	  Lstar_hu_hv[i]=ck.Advection_adjoint(dmom_hv_adv_hu_sge,&vel_grad_test_dV[i_nSpace]);
      	      	  Lstar_hv_hv[i]=ck.Advection_adjoint(dmom_hv_adv_hv_sge,&vel_grad_test_dV[i_nSpace]);

      	      	}

      	      norm_Rv = sqrt(pdeResidual_hu*pdeResidual_hu + pdeResidual_hv*pdeResidual_hv);
	      /* double */
	      double norm_grad = 1.0;
	      q_numDiff_hu[eN_k] = 0.5*elementDiameter[eN]*norm_Rv/(norm_grad+1.0e-8);
	      q_numDiff_hv[eN_k] = q_numDiff_hu[eN_k];

	      //std::cout << q_numDiff_hu[eN_k] << "\t" << q_numDiff_hv[eN_k] << std::endl;
	      /* ck.calculateNumericalDiffusion(1.0, */
	      /* 				     elementDiameter[eN], */
	      /* 				     pdeResidual_h, */
	      /* 				     grad_h, */
	      /* 				     q_numDiff_h[eN_k]); */

      	      //update element residual

      	      for(int i=0;i<nDOF_test_element;i++)
      		{
      		  register int i_nSpace=i*nSpace;

      		  elementResidual_h[i] += ck.Mass_weak(mass_acc_t,h_test_dV[i]) +
      		    ck.Advection_weak(mass_adv,&h_grad_test_dV[i_nSpace]) +
      		    ck.SubgridError(subgridError_h,Lstar_h_h[i]) +
      		    ck.SubgridError(subgridError_hu,Lstar_hu_h[i]) +
      		    ck.SubgridError(subgridError_hv,Lstar_hv_h[i]) +
      		    ck.NumericalDiffusion(q_numDiff_h_last[eN_k],grad_h,&h_grad_test_dV[i_nSpace]);
		  
      		  elementResidual_hu[i] += ck.Mass_weak(mom_hu_acc_t,vel_test_dV[i]) +
      		    ck.Advection_weak(mom_hu_adv,&vel_grad_test_dV[i_nSpace]) +
      		    ck.Diffusion_weak(sdInfo_hu_hu_rowptr,sdInfo_hu_hu_colind,mom_hu_diff_ten,grad_hu,&vel_grad_test_dV[i_nSpace]) +
      		    ck.Diffusion_weak(sdInfo_hu_hv_rowptr,sdInfo_hu_hv_colind,mom_huhv_diff_ten,grad_hv,&vel_grad_test_dV[i_nSpace]) +
		    ck.Reaction_weak(mom_hu_source,vel_test_dV[i]) +
      		    ck.SubgridError(subgridError_h,Lstar_h_hu[i]) +
      		    ck.SubgridError(subgridError_hu,Lstar_hu_hu[i]) +
		    ck.SubgridError(subgridError_hv,Lstar_hv_hu[i]) +
      		    ck.NumericalDiffusion(q_numDiff_hu_last[eN_k],grad_hu,&vel_grad_test_dV[i_nSpace]);
		 
      		  elementResidual_hv[i] += ck.Mass_weak(mom_hv_acc_t,vel_test_dV[i]) +
      		    ck.Advection_weak(mom_hv_adv,&vel_grad_test_dV[i_nSpace]) +
      		    ck.Diffusion_weak(sdInfo_hv_hv_rowptr,sdInfo_hv_hv_colind,mom_hv_diff_ten,grad_hv,&vel_grad_test_dV[i_nSpace]) +
      		    ck.Diffusion_weak(sdInfo_hv_hu_rowptr,sdInfo_hv_hu_colind,mom_hvhu_diff_ten,grad_hu,&vel_grad_test_dV[i_nSpace]) +
      		    ck.Reaction_weak(mom_hv_source,vel_test_dV[i]) +
      		    ck.SubgridError(subgridError_h,Lstar_h_hv[i]) +
      		    ck.SubgridError(subgridError_hu,Lstar_hu_hv[i]) +
      		    ck.SubgridError(subgridError_hv,Lstar_hv_hv[i]) +
      		    ck.NumericalDiffusion(q_numDiff_hv_last[eN_k],grad_hv,&vel_grad_test_dV[i_nSpace]);
      		}
      	    }
      	  
      	  //load element into global residual and save element residual
	    
      	  for(int i=0;i<nDOF_test_element;i++)
      	    {
      	      register int eN_i=eN*nDOF_test_element+i;

      	      elementResidual_h_save[eN_i] +=  elementResidual_h[i];
	  
      	      globalResidual[offset_h+stride_h*h_l2g[eN_i]]+=elementResidual_h[i];
      	      globalResidual[offset_hu+stride_hu*vel_l2g[eN_i]]+=elementResidual_hu[i];
      	      globalResidual[offset_hv+stride_hv*vel_l2g[eN_i]]+=elementResidual_hv[i];
      	    }
      	}
      /* */
      /* loop over exterior element boundaries to calculate surface integrals and load into element and global residuals */
      /* */
      /* ebNE is the Exterior element boundary INdex */
      /* ebN is the element boundary INdex */
      /* eN is the element index */
      /* for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++)  */
      /* 	{  */
      /* 	  register int ebN = exteriorElementBoundariesArray[ebNE],  */
      /* 	    eN  = elementBoundaryElementsArray[ebN*2+0], */
      /* 	    ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+0], */
      /* 	    eN_nDOF_trial_element = eN*nDOF_trial_element; */
      /* 	  register double elementResidual_h[nDOF_test_element], */
      /* 	    elementResidual_hu[nDOF_test_element], */
      /* 	    elementResidual_hv[nDOF_test_element], */
      /* 	    eps_rho,eps_mu; */
      /* 	  for (int i=0;i<nDOF_test_element;i++) */
      /* 	    { */
      /* 	      elementResidual_h[i]=0.0; */
      /* 	      elementResidual_hu[i]=0.0; */
      /* 	      elementResidual_hv[i]=0.0; */
      /* 	    } */
      /* 	  for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++)  */
      /* 	    {  */
      /* 	      register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb, */
      /* 		ebNE_kb_nSpace = ebNE_kb*nSpace, */
      /* 		ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb, */
      /* 		ebN_local_kb_nSpace = ebN_local_kb*nSpace; */
      /* 	      register double h_ext=0.0, */
      /* 		u_ext=0.0, */
      /* 		v_ext=0.0, */
      /* 		grad_h_ext[nSpace], */
      /* 		grad_hu_ext[nSpace], */
      /* 		grad_hv_ext[nSpace], */
      /* 		mom_hu_acc_ext=0.0, */
      /* 		dmom_hu_acc_hu_ext=0.0, */
      /* 		mom_hv_acc_ext=0.0, */
      /* 		dmom_hv_acc_hv_ext=0.0, */
      /* 		mass_adv_ext[nSpace], */
      /* 		dmass_adv_hu_ext[nSpace], */
      /* 		dmass_adv_hv_ext[nSpace], */
      /* 		mom_hu_adv_ext[nSpace], */
      /* 		dmom_hu_adv_hu_ext[nSpace], */
      /* 		dmom_hu_adv_hv_ext[nSpace], */
      /* 		mom_hv_adv_ext[nSpace], */
      /* 		dmom_hv_adv_hu_ext[nSpace], */
      /* 		dmom_hv_adv_hv_ext[nSpace], */
      /* 		mom_hu_diff_ten_ext[nSpace], */
      /* 		mom_hv_diff_ten_ext[nSpace], */
      /* 		mom_huhv_diff_ten_ext[1], */
      /* 		mom_hvhu_diff_ten_ext[1], */
      /* 		mom_hu_source_ext=0.0, */
      /* 		mom_hv_source_ext=0.0, */
      /* 		mom_hu_ham_ext=0.0, */
      /* 		dmom_hu_ham_grad_h_ext[nSpace], */
      /* 		mom_hv_ham_ext=0.0, */
      /* 		dmom_hv_ham_grad_h_ext[nSpace], */
      /* 		dmom_hu_adv_h_ext[nSpace], */
      /* 		dmom_hv_adv_h_ext[nSpace], */
      /* 		flux_mass_ext=0.0, */
      /* 		flux_mom_hu_adv_ext=0.0, */
      /* 		flux_mom_hv_adv_ext=0.0, */
      /* 		flux_mom_hu_diff_ext=0.0, */
      /* 		flux_mom_hv_diff_ext=0.0, */
      /* 		bc_h_ext=0.0, */
      /* 		bc_hu_ext=0.0, */
      /* 		bc_hv_ext=0.0, */
      /* 		bc_mom_hu_acc_ext=0.0, */
      /* 		bc_dmom_hu_acc_hu_ext=0.0, */
      /* 		bc_mom_hv_acc_ext=0.0, */
      /* 		bc_dmom_hv_acc_hv_ext=0.0, */
      /* 		bc_mass_adv_ext[nSpace], */
      /* 		bc_dmass_adv_hu_ext[nSpace], */
      /* 		bc_dmass_adv_hv_ext[nSpace], */
      /* 		bc_mom_hu_adv_ext[nSpace], */
      /* 		bc_dmom_hu_adv_hu_ext[nSpace], */
      /* 		bc_dmom_hu_adv_hv_ext[nSpace], */
      /* 		bc_mom_hv_adv_ext[nSpace], */
      /* 		bc_dmom_hv_adv_hu_ext[nSpace], */
      /* 		bc_dmom_hv_adv_hv_ext[nSpace], */
      /* 		bc_mom_hu_diff_ten_ext[nSpace], */
      /* 		bc_mom_hv_diff_ten_ext[nSpace], */
      /* 		bc_mom_huhv_diff_ten_ext[1], */
      /* 		bc_mom_hvhu_diff_ten_ext[1], */
      /* 		bc_mom_hu_source_ext=0.0, */
      /* 		bc_mom_hv_source_ext=0.0, */
      /* 		bc_mom_hu_ham_ext=0.0, */
      /* 		bc_dmom_hu_ham_grad_h_ext[nSpace], */
      /* 		bc_mom_hv_ham_ext=0.0, */
      /* 		bc_dmom_hv_ham_grad_h_ext[nSpace], */
      /* 		jac_ext[nSpace*nSpace], */
      /* 		jacDet_ext, */
      /* 		jacInv_ext[nSpace*nSpace], */
      /* 		boundaryJac[nSpace*(nSpace-1)], */
      /* 		metricTensor[(nSpace-1)*(nSpace-1)], */
      /* 		metricTensorDetSqrt, */
      /* 		dS,h_test_dS[nDOF_test_element],vel_test_dS[nDOF_test_element], */
      /* 		h_grad_trial_trace[nDOF_trial_element*nSpace],vel_grad_trial_trace[nDOF_trial_element*nSpace], */
      /* 		normal[3],x_ext,y_ext,xt_ext,yt_ext,integralScaling, */
      /* 		G[nSpace*nSpace],G_dd_G,tr_G,h_penalty; */
      /* 	      compute information about mapping from reference element to physical element */
      /* 	      ck.calculateMapping_elementBoundary(eN, */
      /* 						  ebN_local, */
      /* 						  kb, */
      /* 						  ebN_local_kb, */
      /* 						  mesh_dof, */
      /* 						  mesh_l2g, */
      /* 						  mesh_trial_trace_ref, */
      /* 						  mesh_grad_trial_trace_ref, */
      /* 						  boundaryJac_ref, */
      /* 						  jac_ext, */
      /* 						  jacDet_ext, */
      /* 						  jacInv_ext, */
      /* 						  boundaryJac, */
      /* 						  metricTensor, */
      /* 						  metricTensorDetSqrt, */
      /* 						  normal_ref, */
      /* 						  normal, */
      /* 						  x_ext,y_ext); */
      /* 	      ck.calculateMappingVelocity_elementBoundary(eN, */
      /* 							  ebN_local, */
      /* 							  kb, */
      /* 							  ebN_local_kb, */
      /* 							  mesh_velocity_dof, */
      /* 							  mesh_l2g, */
      /* 							  mesh_trial_trace_ref, */
      /* 							  xt_ext,yt_ext, */
      /* 							  normal, */
      /* 							  boundaryJac, */
      /* 							  metricTensor, */
      /* 							  integralScaling); */
      /* 	      xt_ext=0.0;yt_ext=0.0; */
      /* 	      std::cout<<"xt_ext "<<xt_ext<<'\t'<<yt_ext<<'\t'<<std::endl; */
      /* 	      std::cout<<"integralScaling - metricTensorDetSrt ==============================="<<integralScaling-metricTensorDetSqrt<<std::endl; */
      /* 	      dS = ((1.0-MOVING_DOMAIN)*metricTensorDetSqrt + MOVING_DOMAIN*integralScaling)*dS_ref[kb]; */
      /* 	      compute shape and solution information */
      /* 	      shape */
      /* 	      ck.gradTrialFromRef(&h_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,h_grad_trial_trace); */
      /* 	      ck.gradTrialFromRef(&vel_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,vel_grad_trial_trace); */
      /* 	      solution and gradients	 */
      /* 	      ck.valFromDOF(h_dof,&h_l2g[eN_nDOF_trial_element],&h_trial_trace_ref[ebN_local_kb*nDOF_test_element],h_ext); */
      /* 	      ck.valFromDOF(hu_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_trace_ref[ebN_local_kb*nDOF_test_element],u_ext); */
      /* 	      ck.valFromDOF(hv_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_trace_ref[ebN_local_kb*nDOF_test_element],v_ext); */
      /* 	      ck.gradFromDOF(h_dof,&h_l2g[eN_nDOF_trial_element],h_grad_trial_trace,grad_h_ext); */
      /* 	      ck.gradFromDOF(hu_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial_trace,grad_hu_ext); */
      /* 	      ck.gradFromDOF(hv_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial_trace,grad_hv_ext); */
      /* 	      precalculate test function products with integration weights */
      /* 	      for (int j=0;j<nDOF_trial_element;j++) */
      /* 		{ */
      /* 		  h_test_dS[j] = h_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS; */
      /* 		  vel_test_dS[j] = vel_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS; */
      /* 		} */
      /* 	      // */
      /* 	      //debugging section for finite element calculations on exterior */
      /* 	      // */
      /* 	      std::cout<<"ebNE = "<<ebNE<<" kb = "<<kb<<std::endl; */
      /* 	      for (int j=0;j<nDOF_trial_element;j++) */
      /* 	        { */
      /* 	          std::cout<<"h_trial_trace["<<j<<"] "<<h_trial_trace_ref[ebN_local_kb*nDOF_trial_element+j]<<std::endl; */
      /* 	          std::cout<<"vel_trial_trace["<<j<<"] "<<vel_trial_trace_ref[ebN_local_kb*nDOF_trial_element+j]<<std::endl; */
      /* 	          std::cout<<"h_test_dS["<<j<<"] "<<h_test_dS[j]<<std::endl; */
      /* 	          std::cout<<"vel_test_dS["<<j<<"] "<<vel_test_dS[j]<<std::endl; */
      /* 	          for (int I=0;I<nSpace;I++) */
      /* 	      	{ */
      /* 	      	  std::cout<<"h_grad_trial_trace["<<j<<","<<I<<"] "<<h_grad_trial_trace[j*nSpace+I]<<std::endl; */
      /* 	      	  std::cout<<"vel_grad_trial_trace["<<j<<","<<I<<"] "<<vel_grad_trial_trace[j*nSpace+I]<<std::endl; */
      /* 	      	} */
      /* 	        } */
      /* 	      std::cout<<"h_ext "<<h_ext<<std::endl; */
      /* 	      std::cout<<"u_ext "<<u_ext<<std::endl; */
      /* 	      std::cout<<"v_ext "<<v_ext<<std::endl; */
      /* 	      for(int I=0;I<nSpace;I++) */
      /* 	        { */
      /* 	          std::cout<<"grad_h_ext["<<I<<"] "<<grad_h_ext[I]<<std::endl; */
      /* 	          std::cout<<"grad_hu_ext["<<I<<"] "<<grad_hu_ext[I]<<std::endl; */
      /* 	          std::cout<<"grad_hv_ext["<<I<<"] "<<grad_hv_ext[I]<<std::endl; */
      /* 	        } */
      /* 	      */
      /* 	      load the boundary values */
      /* 	      */
      /* 	      bc_h_ext = isDOFBoundary_h[ebNE_kb]*ebqe_bc_h_ext[ebNE_kb]+(1-isDOFBoundary_h[ebNE_kb])*h_ext; */
      /* 	      bc_hu_ext = isDOFBoundary_hu[ebNE_kb]*ebqe_bc_hu_ext[ebNE_kb]+(1-isDOFBoundary_hu[ebNE_kb])*u_ext; */
      /* 	      bc_hv_ext = isDOFBoundary_hv[ebNE_kb]*ebqe_bc_hv_ext[ebNE_kb]+(1-isDOFBoundary_hv[ebNE_kb])*v_ext; */
      /* 	       */
      /* 	      calculate the pde coefficients using the solution and the boundary values for the solution  */
      /* 	       */
      /* 	      cek debug */
      /* 	      eps_rho=0.1; */
      /* 	      eps_mu=0.1; */
      /* 	      evaluateCoefficients(eps_rho, */
      /* 				   eps_mu, */
      /* 				   sigma, */
      /* 				   rho_0, */
      /* 				   nu_0, */
      /* 				   rho_1, */
      /* 				   nu_1, */
      /* 				   g, */
      /* 				   h_ext, */
      /* 				   grad_h_ext, */
      /* 				   u_ext, */
      /* 				   v_ext, */
      /* 				   mom_hu_acc_ext, */
      /* 				   dmom_hu_acc_hu_ext, */
      /* 				   mom_hv_acc_ext, */
      /* 				   dmom_hv_acc_hv_ext, */
      /* 				   mass_adv_ext, */
      /* 				   dmass_adv_hu_ext, */
      /* 				   dmass_adv_hv_ext, */
      /* 				   mom_hu_adv_ext, */
      /* 				   dmom_hu_adv_hu_ext, */
      /* 				   dmom_hu_adv_hv_ext, */
      /* 				   mom_hv_adv_ext, */
      /* 				   dmom_hv_adv_hu_ext, */
      /* 				   dmom_hv_adv_hv_ext, */
      /* 				   mom_hu_diff_ten_ext, */
      /* 				   mom_hv_diff_ten_ext, */
      /* 				   mom_huhv_diff_ten_ext, */
      /* 				   mom_hvhu_diff_ten_ext, */
      /* 				   mom_hu_source_ext, */
      /* 				   mom_hv_source_ext, */
      /* 				   mom_hu_ham_ext, */
      /* 				   dmom_hu_ham_grad_h_ext, */
      /* 				   mom_hv_ham_ext, */
      /* 				   dmom_hv_ham_grad_h_ext);           */
      /* 	      evaluateCoefficients(eps_rho, */
      /* 				   eps_mu, */
      /* 				   sigma, */
      /* 				   rho_0, */
      /* 				   nu_0, */
      /* 				   rho_1, */
      /* 				   nu_1, */
      /* 				   g, */
      /* 				   bc_h_ext, */
      /* 				   grad_h_ext,cek should't be used */
      /* 				   bc_hu_ext, */
      /* 				   bc_hv_ext, */
      /* 				   bc_mom_hu_acc_ext, */
      /* 				   bc_dmom_hu_acc_hu_ext, */
      /* 				   bc_mom_hv_acc_ext, */
      /* 				   bc_dmom_hv_acc_hv_ext, */
      /* 				   bc_mass_adv_ext, */
      /* 				   bc_dmass_adv_hu_ext, */
      /* 				   bc_dmass_adv_hv_ext, */
      /* 				   bc_mom_hu_adv_ext, */
      /* 				   bc_dmom_hu_adv_hu_ext, */
      /* 				   bc_dmom_hu_adv_hv_ext, */
      /* 				   bc_mom_hv_adv_ext, */
      /* 				   bc_dmom_hv_adv_hu_ext, */
      /* 				   bc_dmom_hv_adv_hv_ext, */
      /* 				   bc_mom_hu_diff_ten_ext, */
      /* 				   bc_mom_hv_diff_ten_ext, */
      /* 				   bc_mom_huhv_diff_ten_ext, */
      /* 				   bc_mom_hvhu_diff_ten_ext, */
      /* 				   bc_mom_hu_source_ext, */
      /* 				   bc_mom_hv_source_ext, */
      /* 				   bc_mom_hu_ham_ext, */
      /* 				   bc_dmom_hu_ham_grad_h_ext, */
      /* 				   bc_mom_hv_ham_ext, */
      /* 				   bc_dmom_hv_ham_grad_h_ext);           */
      /* 	      */
      /* 	      moving domain */
      /* 	      */
      /* 	      mass_adv_ext[0] -= MOVING_DOMAIN*xt_ext; */
      /* 	      mass_adv_ext[1] -= MOVING_DOMAIN*yt_ext; */

      /* 	      mom_hu_adv_ext[0] -= MOVING_DOMAIN*mom_hu_acc_ext*xt_ext; */
      /* 	      mom_hu_adv_ext[1] -= MOVING_DOMAIN*mom_hu_acc_ext*yt_ext; */
      /* 	      dmom_hu_adv_hu_ext[0] -= MOVING_DOMAIN*dmom_hu_acc_hu_ext*xt_ext; */
      /* 	      dmom_hu_adv_hu_ext[1] -= MOVING_DOMAIN*dmom_hu_acc_hu_ext*yt_ext; */

      /* 	      mom_hv_adv_ext[0] -= MOVING_DOMAIN*mom_hv_acc_ext*xt_ext; */
      /* 	      mom_hv_adv_ext[1] -= MOVING_DOMAIN*mom_hv_acc_ext*yt_ext; */
      /* 	      dmom_hv_adv_hv_ext[0] -= MOVING_DOMAIN*dmom_hv_acc_hv_ext*xt_ext; */
      /* 	      dmom_hv_adv_hv_ext[1] -= MOVING_DOMAIN*dmom_hv_acc_hv_ext*yt_ext; */

      /* 	      bc's */
      /* 	      bc_mom_hu_adv_ext[0] -= MOVING_DOMAIN*bc_mom_hu_acc_ext*xt_ext; */
      /* 	      bc_mom_hu_adv_ext[1] -= MOVING_DOMAIN*bc_mom_hu_acc_ext*yt_ext; */

      /* 	      bc_mom_hv_adv_ext[0] -= MOVING_DOMAIN*bc_mom_hv_acc_ext*xt_ext; */
      /* 	      bc_mom_hv_adv_ext[1] -= MOVING_DOMAIN*bc_mom_hv_acc_ext*yt_ext; */

      /* 	       */
      /* 	      calculate the numerical fluxes  */
      /* 	       */
      /* 	      cek debug */
      /* 	      ebqe_penalty_ext[ebNE_kb] = 10.0; */
      /* 	      */
      /* 	      ck.calculateGScale(G,normal,h_penalty); */
      /* 	      h_penalty = 10.0/h_penalty; */
      /* 	      cek debug, do it the old way */
      /* 	      h_penalty = 100.0/elementDiameter[eN]; */
      /* 	      exteriorNumericalAdvectiveFlux(isDOFBoundary_h[ebNE_kb], */
      /* 					     isDOFBoundary_hu[ebNE_kb], */
      /* 					     isDOFBoundary_hv[ebNE_kb], */
      /* 					     isAdvectiveFluxBoundary_h[ebNE_kb], */
      /* 					     isAdvectiveFluxBoundary_hu[ebNE_kb], */
      /* 					     isAdvectiveFluxBoundary_hv[ebNE_kb], */
      /* 					     dmom_hu_ham_grad_h_ext[0],=1/rho, */
      /* 					     normal, */
      /* 					     bc_h_ext, */
      /* 					     bc_mass_adv_ext, */
      /* 					     bc_mom_hu_adv_ext, */
      /* 					     bc_mom_hv_adv_ext, */
      /* 					     ebqe_bc_flux_mass_ext[ebNE_kb], */
      /* 					     ebqe_bc_flux_mom_hu_adv_ext[ebNE_kb], */
      /* 					     ebqe_bc_flux_mom_hv_adv_ext[ebNE_kb], */
      /* 					     h_ext, */
      /* 					     mass_adv_ext, */
      /* 					     mom_hu_adv_ext, */
      /* 					     mom_hv_adv_ext, */
      /* 					     dmass_adv_hu_ext, */
      /* 					     dmass_adv_hv_ext, */
      /* 					     dmom_hu_adv_h_ext, */
      /* 					     dmom_hu_adv_hu_ext, */
      /* 					     dmom_hu_adv_hv_ext, */
      /* 					     dmom_hv_adv_h_ext, */
      /* 					     dmom_hv_adv_hu_ext, */
      /* 					     dmom_hv_adv_hv_ext, */
      /* 					     flux_mass_ext, */
      /* 					     flux_mom_hu_adv_ext, */
      /* 					     flux_mom_hv_adv_ext, */
      /* 					     &ebqe_velocity[ebNE_kb_nSpace]); */
      /* 	      cek todo need to switch to full stress and add adjoint consistency */
      /* 	      exteriorNumericalDiffusiveFlux(eps_rho, */
      /* 					     sdInfo_hu_hu_rowptr, */
      /* 					     sdInfo_hu_hu_colind, */
      /* 					     isDOFBoundary_hu[ebNE_kb], */
      /* 					     isDiffusiveFluxBoundary_hu[ebNE_kb], */
      /* 					     normal, */
      /* 					     bc_mom_hu_diff_ten_ext, */
      /* 					     bc_hu_ext, */
      /* 					     ebqe_bc_flux_hu_diff_ext[ebNE_kb], */
      /* 					     mom_hu_diff_ten_ext, */
      /* 					     grad_hu_ext, */
      /* 					     u_ext, */
      /* 					     h_penalty,ebqe_penalty_ext[ebNE_kb], */
      /* 					     flux_mom_hu_diff_ext); */
      /* 	      exteriorNumericalDiffusiveFlux(eps_rho, */
      /* 					     sdInfo_hv_hv_rowptr, */
      /* 					     sdInfo_hv_hv_colind, */
      /* 					     isDOFBoundary_hv[ebNE_kb], */
      /* 					     isDiffusiveFluxBoundary_hv[ebNE_kb], */
      /* 					     normal, */
      /* 					     bc_mom_hv_diff_ten_ext, */
      /* 					     bc_hv_ext, */
      /* 					     ebqe_bc_flux_hv_diff_ext[ebNE_kb], */
      /* 					     mom_hv_diff_ten_ext, */
      /* 					     grad_hv_ext, */
      /* 					     v_ext, */
      /* 					     h_penalty,ebqe_penalty_ext[ebNE_kb], */
      /* 					     flux_mom_hv_diff_ext); */
      /* 	      flux[ebN*nQuadraturePoints_elementBoundary+kb] = flux_mass_ext; */
      /* 	      flux[ebN*nQuadraturePoints_elementBoundary+kb] = 0.0;//cek debug */
      /* 	      */
      /* 	      update residuals */
      /* 	      */
      /* 	      for (int i=0;i<nDOF_test_element;i++) */
      /* 		{ */
      /* 		  elementResidual_h[i] += ck.ExteriorElementBoundaryFlux(flux_mass_ext,h_test_dS[i]); */
      /* 		  globalConservationError += ck.ExteriorElementBoundaryFlux(flux_mass_ext,h_test_dS[i]); */
		  
      /* 		  elementResidual_hu[i] += ck.ExteriorElementBoundaryFlux(flux_mom_hu_adv_ext,vel_test_dS[i])+ */
      /* 		    ck.ExteriorElementBoundaryFlux(flux_mom_hu_diff_ext,vel_test_dS[i]);  */

      /* 		  elementResidual_hv[i] += ck.ExteriorElementBoundaryFlux(flux_mom_hv_adv_ext,vel_test_dS[i]) + */
      /* 		    ck.ExteriorElementBoundaryFlux(flux_mom_hv_diff_ext,vel_test_dS[i]);  */
	       
      /* 		}i */
      /* 	    }kb */
      /* 	  */
      /* 	  update the element and global residual storage */
      /* 	  */
      /* 	  for (int i=0;i<nDOF_test_element;i++) */
      /* 	    { */
      /* 	      int eN_i = eN*nDOF_test_element+i; */
	      
      /* 	      elementResidual_h_save[eN_i] +=  elementResidual_h[i]; */
		  
      /* 	      globalResidual[offset_h+stride_h*h_l2g[eN_i]]+=elementResidual_h[i]; */
      /* 	      globalResidual[offset_hu+stride_hu*vel_l2g[eN_i]]+=elementResidual_hu[i]; */
      /* 	      globalResidual[offset_hv+stride_hv*vel_l2g[eN_i]]+=elementResidual_hv[i]; */
      /* 	    }i */
      /* 	  // */
      /* 	  //debug */
      /* 	  // */
      /* 	  for(int i=0;i<nDOF_test_element;i++)  */
      /* 	    {  */
      /* 	  	  std::cout<<"ebNE "<<ebNE<<" i "<<i<<std::endl; */
      /* 	  	  std::cout<<"r_h"<<elementResidual_h[i]<<std::endl; */
      /* 	  	  std::cout<<"r_hu"<<elementResidual_hu[i]<<std::endl; */
      /* 	  	  std::cout<<"r_hv"<<elementResidual_hv[i]<<std::endl; */
      /* 	  	} */

      /* 	}ebNE */
    }

    void calculateResidual_cell_based_entropy_viscosity(//element
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
							double* h_dof_old_old, 
							double* hu_dof_old_old, 
							double* hv_dof_old_old, 
							double* h_dof_lstage, 
							double* hu_dof_lstage, 
							double* hv_dof_lstage, 
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
							double* q_velocity_sge,
							double* q_cfl,
							double* q_numDiff_h, 
							double* q_numDiff_hu, 
							double* q_numDiff_hv, 
							double* q_numDiff_h_last,
							double* q_numDiff_hu_last, 
							double* q_numDiff_hv_last,
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
							double* edge_based_cfl,
							double cfl_run,
							double hEps,
							int recompute_lumped_mass_matrix, 
							// SAVE SOLUTION (mql)
							double* hnp1_at_quad_point,
							double* hunp1_at_quad_point,
							double* hvnp1_at_quad_point,
							// GALERKIN SOLUTION 
							double* h_dof_galerkin,
							double* hu_dof_galerkin,
							double* hv_dof_galerkin,
							// TO COMPUTE LOW ORDER 
							double* low_order_hnp1,
							double* low_order_hunp1,
							double* low_order_hvnp1,
							double* dEV_minus_dL,
							double cE,
							int LUMPED_MASS_MATRIX,
							int USE_EV_BASED_ON_GALERKIN,
							double dt,
							double mannings)
    {
      // ** COMPUTE QUANTITIES PER CELL (MQL) ** //
      // for linear viscosity //
      double max_speed_per_cell[nElements_global];
      double max_speed = 0, cell_max_speed;
      // for entropy viscosity //
      double entropy_max=-1.E10, entropy_min=1.E10, cell_entropy_mean, entropy_mean=0; 
      double cell_volume, volume=0;
      double cell_entropy_residual, entropy_residual[nElements_global];
      double entropy_normalization_factor=1.0;

      // loop over cells
      for(int eN=0;eN<nElements_global;eN++)
	{
	  cell_max_speed = 0;
	  cell_volume = 0;
	  cell_entropy_mean = 0;
	  cell_entropy_residual = 0;
	  // loop over quadrature points
	  for(int k=0;k<nQuadraturePoints_element;k++)
	    {
	      //get the physical integration weight
	      register double dV,x,y,jac[nSpace*nSpace],jacDet,jacInv[nSpace*nSpace],
		h_grad_trial[nDOF_trial_element*nSpace],vel_grad_trial[nDOF_trial_element*nSpace];
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
      	      dV = fabs(jacDet)*dV_ref[k];
	      //get the trial function gradients
      	      ck.gradTrialFromRef(&h_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,h_grad_trial);
      	      ck.gradTrialFromRef(&vel_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,vel_grad_trial);
	      // SOLUTION AT QUADRATURE POINTS
	      register double hn=0.0, hun=0.0, hvn=0.0, hnm1=0.0, hunm1=0.0, hvnm1=0.0, 
		grad_hn[nSpace],grad_hun[nSpace],grad_hvn[nSpace];
	      register int eN_nDOF_trial_element = eN*nDOF_trial_element;
	      // calculate solution at tn at quadrature points
	      ck.valFromDOF(h_dof_lstage,&h_l2g[eN_nDOF_trial_element],&h_trial_ref[k*nDOF_trial_element],hn);
	      ck.valFromDOF(hu_dof_lstage,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hun);
	      ck.valFromDOF(hv_dof_lstage,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hvn);
	      // calculate solution at tnm1 at quadrature points
      	      ck.valFromDOF(h_dof_old_old,&h_l2g[eN_nDOF_trial_element],&h_trial_ref[k*nDOF_trial_element],hnm1);
      	      ck.valFromDOF(hu_dof_old_old,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hunm1);
      	      ck.valFromDOF(hv_dof_old_old,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hvnm1);
	      // calculate grad of solution at tn at quadrature points
      	      ck.gradFromDOF(h_dof_lstage,&h_l2g[eN_nDOF_trial_element],h_grad_trial,grad_hn);
      	      ck.gradFromDOF(hu_dof_lstage,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_hun);
      	      ck.gradFromDOF(hv_dof_lstage,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_hvn);
	      ///////////////
	      // MAX SPEED //
	      ///////////////
	      double un = hun/hn;
	      double vn = hvn/hn;
	      cell_max_speed = std::max(cell_max_speed,
					std::max(std::abs(un)+std::sqrt(g*hn),std::abs(vn)+std::sqrt(g*hn)));
	      // entropy residual and entropy min and max
	      entropy_max = std::max(entropy_max,ENTROPY(g,hn,un,vn,1./hn));
	      entropy_min = std::min(entropy_min,ENTROPY(g,hn,un,vn,1./hn));
	      cell_entropy_mean += ENTROPY(g,hn,un,vn,1./hn)*dV;
	      cell_volume += dV;
	      double gradX_Entropy = D_ENTROPY(g,hn,hun,hvn,grad_hn[0],grad_hun[0],grad_hvn[0],1./hn);
	      double gradY_Entropy = D_ENTROPY(g,hn,hun,hvn,grad_hn[1],grad_hun[1],grad_hvn[1],1./hn);
	      cell_entropy_residual 
		= std::max(cell_entropy_residual,
			   std::abs(
				    (ENTROPY(g,hn,hun,hvn,1./hn) - ENTROPY(g,hnm1,hunm1,hvnm1,1./hnm1))/dt
				    +ENTROPY(g,hn,hun,hvn,1./hn)/hn/hn*(hn*(grad_hun[0]+grad_hvn[1])-(hun*grad_hn[0]+hvn*grad_hn[1]))
				    +1./hn*(gradX_Entropy*hun+gradY_Entropy*hvn)
				    +0.5*g*hn*(grad_hun[0]+grad_hvn[1])
				    +0.5*g*(hun*grad_hn[0]+hvn*grad_hn[1])));
	    }
	  max_speed_per_cell[eN] = cell_max_speed;
	  max_speed = std::max(max_speed,cell_max_speed);
	  volume += cell_volume;
	  entropy_mean += cell_entropy_mean;
	  entropy_residual[eN] = cell_entropy_residual;
	}
      entropy_mean /= volume;
      entropy_normalization_factor = std::max(std::abs(entropy_max-entropy_mean),
      				      std::abs(entropy_min-entropy_mean));
      //entropy_normalization_factor = entropy_max - entropy_min;
      //
      //loop over elements to compute volume integrals and load them into element and global residual
      //
      double globalConservationError=0.0,tauSum=0.0;
      for(int eN=0;eN<nElements_global;eN++)
      	{
      	  //declare local storage for element residual and initialize
      	  register double elementResidual_h[nDOF_test_element],
      	    elementResidual_hu[nDOF_test_element],
      	    elementResidual_hv[nDOF_test_element];
      	  for (int i=0;i<nDOF_test_element;i++)
      	    {
      	      int eN_i = eN*nDOF_test_element+i;
      	      elementResidual_h_save[eN_i]=0.0;
      	      elementResidual_h[i]=0.0;
      	      elementResidual_hu[i]=0.0;
      	      elementResidual_hv[i]=0.0;
      	    }//i
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
		b=0.0,h=0.0,hu=0.0,hv=0.0, // solution at current time
		h_tn=0.0, hu_tn=0.0, hv_tn=0.0, // solution at tn
		h_star=0.0, hu_star=0.0, hv_star=0.0, // solution at t star
      		grad_b[nSpace],grad_h[nSpace],grad_hu[nSpace],grad_hv[nSpace], //grad at current time
		grad_h_tn[nSpace],grad_hu_tn[nSpace],grad_hv_tn[nSpace], //grad at tn
		grad_h_star[nSpace],grad_hu_star[nSpace],grad_hv_star[nSpace], //grad at t star
      		mass_acc=0.0,mom_hu_acc=0.0,mom_hv_acc=0.0, //accumulation variables 
      		dmass_acc_h=0.0, dmom_hu_acc_hu=0.0, dmom_hv_acc_hv=0.0,
		mass_adv[nSpace], mom_hu_adv[nSpace], mom_hv_adv[nSpace], //adv terms at current time
		mom_hu_source=0.0, mom_hv_source=0.0, //source terms at current time
		mass_acc_t=0.0, dmass_acc_h_t=0.0, //dt of mass accumulation
      		mom_hu_acc_t=0.0, dmom_hu_acc_h_t=0.0, dmom_hu_acc_hu_t=0.0, //dt of x-mom accumulation 
      		mom_hv_acc_t=0.0, dmom_hv_acc_h_t=0.0, dmom_hv_acc_hv_t=0.0, //dt of y-mom accumulation 
		jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace],
      		h_grad_trial[nDOF_trial_element*nSpace],vel_grad_trial[nDOF_trial_element*nSpace],
      		h_test_dV[nDOF_trial_element],vel_test_dV[nDOF_trial_element],
      		h_grad_test_dV[nDOF_test_element*nSpace],vel_grad_test_dV[nDOF_test_element*nSpace],
      		dV,x,y,xt,yt;

	      // FOR EXPLICIT TIME INTEGRATION 
	      register double 
		mass_acc_star, mom_hu_acc_star, mom_hv_acc_star,
		mass_adv_star[nSpace], mom_hu_adv_star[nSpace], mom_hv_adv_star[nSpace], 
      		mom_hu_source_star=0.0, mom_hv_source_star=0.0; 

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
      	      //get the trial function gradients
      	      ck.gradTrialFromRef(&h_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,h_grad_trial);
      	      ck.gradTrialFromRef(&vel_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,vel_grad_trial);
      	      //get the solution at current time
      	      ck.valFromDOF(b_dof,&h_l2g[eN_nDOF_trial_element],&h_trial_ref[k*nDOF_trial_element],b);
      	      ck.valFromDOF(h_dof,&h_l2g[eN_nDOF_trial_element],&h_trial_ref[k*nDOF_trial_element],h);
      	      ck.valFromDOF(hu_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hu);
      	      ck.valFromDOF(hv_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hv);
	      //get the solution at time tn (old time)
	      ck.valFromDOF(h_dof_lstage,&h_l2g[eN_nDOF_trial_element],&h_trial_ref[k*nDOF_trial_element],h_tn);
      	      ck.valFromDOF(hu_dof_lstage,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hu_tn);
      	      ck.valFromDOF(hv_dof_lstage,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hv_tn);
	      //get the solution gradients at current time
      	      ck.gradFromDOF(b_dof,&h_l2g[eN_nDOF_trial_element],h_grad_trial,grad_b);
      	      ck.gradFromDOF(h_dof,&h_l2g[eN_nDOF_trial_element],h_grad_trial,grad_h);
      	      ck.gradFromDOF(hu_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_hu);
      	      ck.gradFromDOF(hv_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_hv);
	      //get the solution gradients at tn (old time)
      	      ck.gradFromDOF(h_dof_lstage,&h_l2g[eN_nDOF_trial_element],h_grad_trial,grad_h_tn);
      	      ck.gradFromDOF(hu_dof_lstage,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_hu_tn);
      	      ck.gradFromDOF(hv_dof_lstage,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_hv_tn);
      	      //precalculate test function products with integration weights
      	      for (int j=0;j<nDOF_trial_element;j++)
      		{
      		  h_test_dV[j] = h_test_ref[k*nDOF_trial_element+j]*dV;
      		  vel_test_dV[j] = vel_test_ref[k*nDOF_trial_element+j]*dV;
      		  for (int I=0;I<nSpace;I++)
      		    {
      		      h_grad_test_dV[j*nSpace+I]   = h_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
      		      vel_grad_test_dV[j*nSpace+I] = vel_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
      		    }
      		}
	      // COMPUTE solution "star" to allow quick change between implicit or explicit time integration
	      h_star = IMPLICIT*h+(1-IMPLICIT)*h_tn;
	      hu_star = IMPLICIT*hu+(1-IMPLICIT)*hu_tn;
	      hv_star = IMPLICIT*hv+(1-IMPLICIT)*hv_tn;
	      for (int I=0; I<nSpace; I++)
		{
		  grad_h_star[I] = IMPLICIT*grad_h[I]+(1-IMPLICIT)*grad_h_tn[I];
		  grad_hu_star[I] = IMPLICIT*grad_hu[I]+(1-IMPLICIT)*grad_hu_tn[I];
		  grad_hv_star[I] = IMPLICIT*grad_hv[I]+(1-IMPLICIT)*grad_hv_tn[I];
		}
      	      //save velocity at quadrature points for other models to use
      	      q_velocity[eN_k_nSpace+0]=hu/h;
      	      q_velocity[eN_k_nSpace+1]=hv/h;
      	      //
      	      //calculate pde coefficients at quadrature points
      	      //
      	      
	      evaluateCoefficientsForResidual( // WITH "CURRENT" SOLUTION
					      // ********** INPUT ********** //
					      g, // gravity
					      grad_b, // grad of bathymetry
					      h,
					      hu,
					      hv,
					      // ********** OUTPUT ********** //
					      mass_acc, 
					      mom_hu_acc,
					      mom_hv_acc,
					      mass_adv, // [h*u, h*v] 
					      mom_hu_adv, // [h*u^2+0.5*g*h^2, h*u*v]
					      mom_hv_adv, // [h*u*v, h*v^2+0.5*g*h^2]
					      mom_hu_source, // x-momentum source
					      mom_hv_source); // y-momentum source
	      evaluateCoefficientsForResidual( // WITH "STAR" SOLUTION
					      // ********** INPUT ********** //
					      g, // gravity
					      grad_b, // grad of bathymetry
					      h_star,
					      hu_star,
					      hv_star,
					      // ********** OUTPUT ********** //
					      mass_acc_star, //dummy
					      mom_hu_acc_star, //dummy
					      mom_hv_acc_star, //dummy
					      mass_adv_star, // [h*u, h*v] 
					      mom_hu_adv_star, // [h*u^2+0.5*g*h^2, h*u*v]
					      mom_hv_adv_star, // [h*u*v, h*v^2+0.5*g*h^2]
					      mom_hu_source_star, // x-momentum source
					      mom_hv_source_star); // y-momentum source
      	      //
      	      //save momentum for time history and velocity for subgrid error
      	      //
      	      q_mass_acc[eN_k] = mass_acc;
      	      q_mom_hu_acc[eN_k] = mom_hu_acc;
      	      q_mom_hv_acc[eN_k] = mom_hv_acc;
      	      //subgrid error uses grid scale discharge
      	      q_mass_adv[eN_k_nSpace+0] = hu;
      	      q_mass_adv[eN_k_nSpace+1] = hv;
      	      //
      	      //moving mesh (TODO)
      	      //
      	      //
      	      //calculate time derivative at quadrature points
      	      //
      	      ck.bdf(alphaBDF,
      		     q_mass_acc_beta_bdf[eN_k],
      		     mass_acc,
      		     dmass_acc_h,
      		     mass_acc_t,
      		     dmass_acc_h_t);
      	      ck.bdf(alphaBDF,
		     q_mom_hu_acc_beta_bdf[eN_k],
		     mom_hu_acc,
		     dmom_hu_acc_hu,
		     mom_hu_acc_t,
		     dmom_hu_acc_hu_t);
      	      ck.bdf(alphaBDF,
		     q_mom_hv_acc_beta_bdf[eN_k],
		     mom_hv_acc,
		     dmom_hv_acc_hv,
		     mom_hv_acc_t,
		     dmom_hv_acc_hv_t);

	      calculateCFL(elementDiameter[eN],
			   g,
			   h_tn,
			   hu_tn,
			   hv_tn,
			   hEps,
			   q_cfl[eN_k]);

	      /////////////////////////////////
	      // COMPUTE NUMERICAL DIFFUSION //
	      /////////////////////////////////
	      // LINEAR VISCOSITY //
	      double linear_viscosity = cMax*elementDiameter[eN]*max_speed_per_cell[eN];
	      q_numDiff_h[eN_k] = linear_viscosity;
	      q_numDiff_hu[eN_k] = linear_viscosity;
	      q_numDiff_hv[eN_k] = linear_viscosity;
	      // ENTROPY VISCOSITY //
	      double entropy_viscosity = cE*std::pow(elementDiameter[eN],2)*
		entropy_residual[eN]/entropy_normalization_factor;
	      // NUMERICAL VISCOSITY //
	      q_numDiff_h[eN_k] = std::min(linear_viscosity,entropy_viscosity);
	      q_numDiff_hu[eN_k] = std::min(linear_viscosity,entropy_viscosity);
	      q_numDiff_hv[eN_k] = std::min(linear_viscosity,entropy_viscosity);

      	      //update element residual
      	      for(int i=0;i<nDOF_test_element;i++)
      		{
      		  register int i_nSpace=i*nSpace;
		  int eN_i=eN*nDOF_test_element+i;
		  int h_gi = h_l2g[eN_i]; //global i-th index for h variable
		  int vel_gi = vel_l2g[eN_i]; //global i-th index for velocity variables

      		  elementResidual_h[i] += 
		    dt*ck.Mass_weak(mass_acc_t,h_test_dV[i]) + // Mass matrix is NOT lumped
      		    dt*ck.Advection_weak(mass_adv_star,&h_grad_test_dV[i_nSpace]) +
		    dt*ck.NumericalDiffusion(q_numDiff_h_last[eN_k],grad_h_star,&h_grad_test_dV[i_nSpace]);
		  
      		  elementResidual_hu[i] += 
		    dt*ck.Mass_weak(mom_hu_acc_t,vel_test_dV[i]) + // Mass matrix is NOT lumped
      		    dt*ck.Advection_weak(mom_hu_adv_star,&vel_grad_test_dV[i_nSpace]) +
		    //dt*ck.Reaction_weak(mom_hu_source_star,vel_test_dV[i]) +
		    dt*ck.NumericalDiffusion(q_numDiff_hu_last[eN_k],grad_hu_star,&vel_grad_test_dV[i_nSpace]);
		 
      		  elementResidual_hv[i] += 
		    dt*ck.Mass_weak(mom_hv_acc_t,vel_test_dV[i]) + // Mass matrix is NOT lumped
      		    dt*ck.Advection_weak(mom_hv_adv_star,&vel_grad_test_dV[i_nSpace]) +
		    //dt*ck.Reaction_weak(mom_hv_source_star,vel_test_dV[i]) +
		    dt*ck.NumericalDiffusion(q_numDiff_hv_last[eN_k],grad_hv_star,&vel_grad_test_dV[i_nSpace]);
      		}
      	    }
      	  
      	  //load element into global residual and save element residual
	    
      	  for(int i=0;i<nDOF_test_element;i++)
      	    {
      	      register int eN_i=eN*nDOF_test_element+i;
	      int h_gi = h_l2g[eN_i]; //global i-th index for h
	      int vel_gi = vel_l2g[eN_i]; //global i-th index for velocities 
		
	      elementResidual_h_save[eN_i] +=  elementResidual_h[i];//* (h_dof[h_gi] - h_dof_lstage[h_gi]);
	        	      
      	      globalResidual[offset_h+stride_h*h_gi]  += elementResidual_h[i];
      	      globalResidual[offset_hu+stride_hu*vel_gi] += elementResidual_hu[i];
      	      globalResidual[offset_hv+stride_hv*vel_gi] += elementResidual_hv[i];
      	    }
      	}
    }
    
    void calculateResidual_first_order_flatB_GP(//element
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
				    double* h_dof_old_old, 
				    double* hu_dof_old_old, 
				    double* hv_dof_old_old, 
				    double* h_dof_lstage, 
				    double* hu_dof_lstage, 
				    double* hv_dof_lstage, 
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
				    double* q_velocity_sge,
				    double* q_cfl,
				    double* q_numDiff_h, 
				    double* q_numDiff_hu, 
				    double* q_numDiff_hv, 
				    double* q_numDiff_h_last,
				    double* q_numDiff_hu_last, 
				    double* q_numDiff_hv_last,
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
				    double* edge_based_cfl,
				    double cfl_run,
				    double hEps,
				    int recompute_lumped_mass_matrix,
				    // SAVE SOLUTION (mql)
				    double* hnp1_at_quad_point,
				    double* hunp1_at_quad_point,
				    double* hvnp1_at_quad_point,
				    // GALERKIN SOLUTION 
				    double* h_dof_galerkin,
				    double* hu_dof_galerkin,
				    double* hv_dof_galerkin,
				    // TO COMPUTE LOW ORDER 
				    double* low_order_hnp1,
				    double* low_order_hunp1,
				    double* low_order_hvnp1,
				    double* dEV_minus_dL,
				    double cE,
				    int LUMPED_MASS_MATRIX,
				    int USE_EV_BASED_ON_GALERKIN,
				    double dt, 
				    double mannings)
    {
      ////////////////
      // CELL LOOPS //
      ////////////////
      // To compute: 
      //      * lumped_mass_matrix
      //      * Cell based CFL
      //      * velocity at quad points for other models 
      // init lumped mass matrix to zero
      if (recompute_lumped_mass_matrix==1)
	{
	  for (int i=0; i<numDOFsPerEqn; i++)
	    lumped_mass_matrix[i] = 0;
	  for(int eN=0;eN<nElements_global;eN++)
	    {
	      //declare local storage for element residual and initialize
	      register double element_lumped_mass_matrix[nDOF_test_element];
	      for (int i=0;i<nDOF_test_element;i++)
		element_lumped_mass_matrix[i]=0.0;
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
		    b=0.0,h=0.0,hu=0.0,hv=0.0, // solution at current time
		    h_tn=0.0, hu_tn=0.0, hv_tn=0.0, // solution at tn
		    jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace],
		    h_test_dV[nDOF_trial_element],vel_test_dV[nDOF_trial_element],
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
		  //get the solution at current time. This is to compute velocity for other models
		  ck.valFromDOF(b_dof,&h_l2g[eN_nDOF_trial_element],&h_trial_ref[k*nDOF_trial_element],b);
		  ck.valFromDOF(h_dof,&h_l2g[eN_nDOF_trial_element],&h_trial_ref[k*nDOF_trial_element],h);
		  ck.valFromDOF(hu_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hu);
		  ck.valFromDOF(hv_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hv);
		  //get the solution at time tn (old time). This is needed to compute the CFL
		  ck.valFromDOF(h_dof_lstage,&h_l2g[eN_nDOF_trial_element],&h_trial_ref[k*nDOF_trial_element],h_tn);
		  ck.valFromDOF(hu_dof_lstage,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hu_tn);
		  ck.valFromDOF(hv_dof_lstage,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hv_tn);
		  //precalculate test function products with integration weights
		  for (int j=0;j<nDOF_trial_element;j++)
		    {
		      h_test_dV[j] = h_test_ref[k*nDOF_trial_element+j]*dV;
		      vel_test_dV[j] = vel_test_ref[k*nDOF_trial_element+j]*dV;
		    }
		  //save velocity at quadrature points for other models to use
		  q_velocity[eN_k_nSpace+0] = 2*h/(h*h+std::pow(fmax(h,hEps),2))*hu;
		  q_velocity[eN_k_nSpace+1] = 2*h/(h*h+std::pow(fmax(h,hEps),2))*hv;
		  hnp1_at_quad_point[eN_k] = h;
		  // calculatte CFL
		  calculateCFL(elementDiameter[eN],
			       g,
			       h_tn,
			       hu_tn,
			       hv_tn,
			       hEps,
			       q_cfl[eN_k]);
		  //update element residual. Part about the lumped mass matrix 
		  for(int i=0;i<nDOF_test_element;i++)
		    element_lumped_mass_matrix[i] += h_test_dV[i];
		}
	      // distribute
	      for(int i=0;i<nDOF_test_element;i++)
		{
		  register int eN_i=eN*nDOF_test_element+i;
		  int h_gi = h_l2g[eN_i]; //global i-th index for h
		  lumped_mass_matrix[h_gi]  += element_lumped_mass_matrix[i];	  
		}
	    }
	}

      //////////////////
      // Loop on DOFs //
      //////////////////
      int ij = 0;
      for (int i=0; i<numDOFsPerEqn; i++)
	{
	  double hi = h_dof_lstage[i];
	  double hui = hu_dof_lstage[i];
	  double hvi = hv_dof_lstage[i];
	  // regularization of 1/hi 
	  double one_over_hiReg = 2*hi/(hi*hi+std::pow(fmax(hi,hEps),2));

	  double ith_flux_term1=0., ith_flux_term2=0., ith_flux_term3=0.;
	  double ith_dissipative_term1=0., ith_dissipative_term2=0., ith_dissipative_term3=0.;

	  double dLii = 0;
	  // loop over the sparsity pattern of the i-th DOF
	  for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
	    {
	      int j = csrColumnOffsets_DofLoops[offset];
	      double hj = h_dof_lstage[j];
	      double huj = hu_dof_lstage[j];
	      double hvj = hv_dof_lstage[j];

	      // regularization of 1/hj
	      double one_over_hjReg = 2*hj/(hj*hj+std::pow(fmax(hj,hEps),2));

	      // Nodal projection of fluxes
	      ith_flux_term1 += huj*Cx[ij] + hvj*Cy[ij]; // f1*C
	      ith_flux_term2 += huj*huj*one_over_hjReg*Cx[ij] + huj*hvj*one_over_hjReg*Cy[ij] + 0.5*g*hj*hj*Cx[ij];
	      ith_flux_term3 += huj*hvj*one_over_hjReg*Cx[ij] + hvj*hvj*one_over_hjReg*Cy[ij] + 0.5*g*hj*hj*Cy[ij];

	      // Dissipative term
	      double dLij = 0;
	      if (i != j) // This is not necessary. See formula for ith_dissipative_terms
		{
		  // norm of the C and C transpose matrices
		  double cij_norm = sqrt(Cx[ij]*Cx[ij] + Cy[ij]*Cy[ij]);
		  double cji_norm = sqrt(CTx[ij]*CTx[ij] + CTy[ij]*CTy[ij]);

		  double nxij = Cx[ij]/cij_norm, nyij = Cy[ij]/cij_norm;
		  double nxji = CTx[ij]/cji_norm, nyji = CTy[ij]/cji_norm;

		  //dLij =  fmax(maxWaveSpeedIterativeProcess(g,nxij,nyij,
		  //		    hi,hui,hvi,hj,huj,hvj)*cij_norm,
		  //   maxWaveSpeedIterativeProcess(g,nxji,nyji,
		  //	    hj,huj,hvj,hi,hui,hvi)*cji_norm);
		  dLij =  fmax(maxWaveSpeedSharpInitialGuess(g,nxij,nyij,
							     hi,hui,hvi,hj,huj,hvj,hEps,false)*cij_norm,
			       maxWaveSpeedSharpInitialGuess(g,nxji,nyji,
							     hj,huj,hvj,hi,hui,hvi,hEps,false)*cji_norm);
		 
		  ith_dissipative_term1 += dLij*(hj-hi);
		  ith_dissipative_term2 += dLij*(huj-hui);
		  ith_dissipative_term3 += dLij*(hvj-hvi);
		  // compute dLii (for debugging and selecting time step) 
		  dLii -= dLij;
		}
	      // update ij
	      ij+=1;
	    }
	  double mi = lumped_mass_matrix[i];
	  // compute residual
	  globalResidual[offset_h+stride_h*i]   = mi*(h_dof[i] - hi) + dt*(ith_flux_term1 - ith_dissipative_term1);
	  globalResidual[offset_hu+stride_hu*i] = mi*(hu_dof[i] - hui) + dt*(ith_flux_term2 - ith_dissipative_term2);
	  globalResidual[offset_hv+stride_hv*i] = mi*(hv_dof[i] - hvi) + dt*(ith_flux_term3 - ith_dissipative_term3);

	  // calculate edge based CFL
	  edge_based_cfl[i] = 4*fabs(dLii)/mi;
	}
    }
 
    void calculateResidual_second_order_flatB_GP(//element
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
				       double* h_dof_old_old, 
				       double* hu_dof_old_old, 
				       double* hv_dof_old_old, 
				       double* h_dof_lstage, 
				       double* hu_dof_lstage, 
				       double* hv_dof_lstage, 
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
				       double* q_velocity_sge,
				       double* q_cfl,
				       double* q_numDiff_h, 
				       double* q_numDiff_hu, 
				       double* q_numDiff_hv, 
				       double* q_numDiff_h_last,
				       double* q_numDiff_hu_last, 
				       double* q_numDiff_hv_last,
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
				       double* edge_based_cfl,
				       double cfl_run,
				       double hEps,
				       int recompute_lumped_mass_matrix, 
				       // SAVE SOLUTION (mql)
				       double* hnp1_at_quad_point,
				       double* hunp1_at_quad_point,
				       double* hvnp1_at_quad_point,
				       // GALERKIN SOLUTION 
				       double* h_dof_galerkin,
				       double* hu_dof_galerkin,
				       double* hv_dof_galerkin,
				       // TO COMPUTE LOW ORDER 
				       double* low_order_hnp1,
				       double* low_order_hunp1,
				       double* low_order_hvnp1,
				       double* dEV_minus_dL,
				       double cE, 
				       int LUMPED_MASS_MATRIX,
				       int USE_EV_BASED_ON_GALERKIN,
				       double dt,
				       double mannings)
    {
      ////////////////
      // CELL LOOPS //
      ////////////////
      // To compute: 
      //      * lumped_mass_matrix
      //      * Cell based CFL
      //      * velocity at quad points for other models 
      // init lumped mass matrix to zero
      if (recompute_lumped_mass_matrix==1)
	{
	  for (int i=0; i<numDOFsPerEqn; i++)
	    lumped_mass_matrix[i] = 0;
	  for(int eN=0;eN<nElements_global;eN++)
	    {
	      //declare local storage for element residual and initialize
	      register double element_lumped_mass_matrix[nDOF_test_element];
	      for (int i=0;i<nDOF_test_element;i++)
		element_lumped_mass_matrix[i]=0.0;
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
		    b=0.0,h=0.0,hu=0.0,hv=0.0, // solution at current time
		    h_tn=0.0, hu_tn=0.0, hv_tn=0.0, // solution at tn
		    jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace],
		    h_test_dV[nDOF_trial_element],vel_test_dV[nDOF_trial_element],
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
		  //get the solution at current time. This is to compute velocity for other models
		  ck.valFromDOF(b_dof,&h_l2g[eN_nDOF_trial_element],&h_trial_ref[k*nDOF_trial_element],b);
		  ck.valFromDOF(h_dof,&h_l2g[eN_nDOF_trial_element],&h_trial_ref[k*nDOF_trial_element],h);
		  ck.valFromDOF(hu_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hu);
		  ck.valFromDOF(hv_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hv);
		  //get the solution at time tn (old time). This is needed to compute the CFL
		  ck.valFromDOF(h_dof_lstage,&h_l2g[eN_nDOF_trial_element],&h_trial_ref[k*nDOF_trial_element],h_tn);
		  ck.valFromDOF(hu_dof_lstage,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hu_tn);
		  ck.valFromDOF(hv_dof_lstage,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hv_tn);
		  //precalculate test function products with integration weights
		  for (int j=0;j<nDOF_trial_element;j++)
		    {
		      h_test_dV[j] = h_test_ref[k*nDOF_trial_element+j]*dV;
		      vel_test_dV[j] = vel_test_ref[k*nDOF_trial_element+j]*dV;
		    }
		  //save velocity at quadrature points for other models to use
		  q_velocity[eN_k_nSpace+0] = 2*h/(h*h+std::pow(fmax(h,hEps),2))*hu;
		  q_velocity[eN_k_nSpace+1] = 2*h/(h*h+std::pow(fmax(h,hEps),2))*hv;
		  hnp1_at_quad_point[eN_k] = h;
		  hunp1_at_quad_point[eN_k] = hu;
		  hvnp1_at_quad_point[eN_k] = hv;
		  // calculatte CFL
		  calculateCFL(elementDiameter[eN],
			       g,
			       h_tn,
			       hu_tn,
			       hv_tn,
			       hEps,
			       q_cfl[eN_k]);
		  //update element residual. Part about the lumped mass matrix 
		  for(int i=0;i<nDOF_test_element;i++)
		    element_lumped_mass_matrix[i] += h_test_dV[i];
		}
	      // distribute
	      for(int i=0;i<nDOF_test_element;i++)
		{
		  register int eN_i=eN*nDOF_test_element+i;
		  int h_gi = h_l2g[eN_i]; //global i-th index for h
		  lumped_mass_matrix[h_gi]  += element_lumped_mass_matrix[i];	  
		}
	    }
	} // end of computing lumped mass matrix

      //////////////////////////////////
      // COMPUTE SMOOTHNESS INDICATOR // and // COMPUTE LOCAL MAX OF ENT RESIDUALS //
      //////////////////////////////////
      // Smoothness indicator is based on the solution. psi_i = psi_i(alpha_i); 
      // alpha_i = |sum(uj-ui)|/sum|uj-ui|
      register double psi[numDOFsPerEqn];
      for (int i=0; i<numDOFsPerEqn; i++)
	{
	  double alphai, alphai_numerator=0, alphai_denominator=0; // smoothness indicator of solution
	  double hi = h_dof_lstage[i]; // solution at time tn for the ith DOF
	  for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
	    { //loop in j (sparsity pattern)
	      int j = csrColumnOffsets_DofLoops[offset];
	      double hj = h_dof_lstage[j]; // solution at time tn for the jth DOF
	      alphai_numerator += (hj-hi);
	      alphai_denominator += std::abs(hj-hi);
	    }
	  if (hi < hEps)
	    alphai = 1.;
	  else
	    {
	      if (alphai_denominator==0)
		alphai = 1.;
	      else 
		alphai = std::abs(alphai_numerator)/alphai_denominator;
	    }
	  if (POWER_SMOOTHNESS_INDICATOR==0)
	    psi[i] = 1.0;
	  else
	    psi[i] = std::pow(alphai,POWER_SMOOTHNESS_INDICATOR); //NOTE: they use alpha^2 in the paper
	}

      //////////////////
      // Loop on DOFs //
      //////////////////
      int ij = 0;
      for (int i=0; i<numDOFsPerEqn; i++)
	{
	  double hi = h_dof_lstage[i];
	  double hui = hu_dof_lstage[i];
	  double hvi = hv_dof_lstage[i];
	  // regularization of 1/hi 
	  double one_over_hiReg = 2*hi/(hi*hi+std::pow(fmax(hi,hEps),2));

	  double ith_flux_term1=0., ith_flux_term2=0., ith_flux_term3=0.;
	  double ith_dissipative_term1=0., ith_dissipative_term2=0., ith_dissipative_term3=0.;

	  double dLii = 0;
	  // loop over the sparsity pattern of the i-th DOF
	  for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
	    {
	      int j = csrColumnOffsets_DofLoops[offset];
	      double hj = h_dof_lstage[j];
	      double huj = hu_dof_lstage[j];
	      double hvj = hv_dof_lstage[j];

	      // regularization of 1/hj
	      double one_over_hjReg = 2*hj/(hj*hj+std::pow(fmax(hj,hEps),2));

	      // Nodal projection of fluxes
	      ith_flux_term1 += huj*Cx[ij] + hvj*Cy[ij]; // f1*C
	      ith_flux_term2 += huj*huj*one_over_hjReg*Cx[ij] + huj*hvj*one_over_hjReg*Cy[ij] + g*hi*hj*Cx[ij];
	      ith_flux_term3 += huj*hvj*one_over_hjReg*Cx[ij] + hvj*hvj*one_over_hjReg*Cy[ij] + g*hi*hj*Cy[ij];

	      // Dissipative term
	      double dLij = 0;
	      if (i != j) // This is not necessary. See formula for ith_dissipative_terms
		{
		  // norm of the C and C transpose matrices
		  double cij_norm = sqrt(Cx[ij]*Cx[ij] + Cy[ij]*Cy[ij]);
		  double cji_norm = sqrt(CTx[ij]*CTx[ij] + CTy[ij]*CTy[ij]);

		  double nxij = Cx[ij]/cij_norm, nyij = Cy[ij]/cij_norm;
		  double nxji = CTx[ij]/cji_norm, nyji = CTy[ij]/cji_norm;

		  //dLij =  fmax(maxWaveSpeedIterativeProcess(g,nxij,nyij,
		  //				    hi,hui,hvi,hj,huj,hvj)*cij_norm,
		  //       maxWaveSpeedIterativeProcess(g,nxji,nyji,
		  //				    hj,huj,hvj,hi,hui,hvi)*cji_norm);
		  dLij =  fmax(maxWaveSpeedSharpInitialGuess(g,nxij,nyij,
							     hi,hui,hvi,hj,huj,hvj,hEps,false)*cij_norm,
			       maxWaveSpeedSharpInitialGuess(g,nxji,nyji,
							     hj,huj,hvj,hi,hui,hvi,hEps,false)*cji_norm);
		  dLij *= std::max(psi[i],psi[j]); 
		  
		  ith_dissipative_term1 += dLij*(hj-hi);
		  ith_dissipative_term2 += dLij*(huj-hui);
		  ith_dissipative_term3 += dLij*(hvj-hvi);
		  // compute dLii (for debugging and selecting time step) 
		  dLii -= dLij;
		}
	      // update ij
	      ij+=1;
	    }
	  double mi = lumped_mass_matrix[i];
	  // compute residualo
	  globalResidual[offset_h+stride_h*i]   = mi*(h_dof[i] - hi) + dt*(ith_flux_term1 - ith_dissipative_term1);
	  globalResidual[offset_hu+stride_hu*i] = mi*(hu_dof[i] - hui) + dt*(ith_flux_term2 - ith_dissipative_term2);
	  globalResidual[offset_hv+stride_hv*i] = mi*(hv_dof[i] - hvi) + dt*(ith_flux_term3 - ith_dissipative_term3);

	  // calculate edge based CFL
	  edge_based_cfl[i] = 4*fabs(dLii)/mi;
	}
    }

    void calculateResidual_second_order_NonFlatB_GP(//element
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
				       double* h_dof_old_old, 
				       double* hu_dof_old_old, 
				       double* hv_dof_old_old, 
				       double* h_dof_lstage, 
				       double* hu_dof_lstage, 
				       double* hv_dof_lstage, 
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
				       double* q_velocity_sge,
				       double* q_cfl,
				       double* q_numDiff_h, 
				       double* q_numDiff_hu, 
				       double* q_numDiff_hv, 
				       double* q_numDiff_h_last,
				       double* q_numDiff_hu_last, 
				       double* q_numDiff_hv_last,
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
				       double* edge_based_cfl,
				       double cfl_run,
				       double hEps,
				       int recompute_lumped_mass_matrix, 
				       // SAVE SOLUTION (mql)
				       double* hnp1_at_quad_point,
				       double* hunp1_at_quad_point,
				       double* hvnp1_at_quad_point,
				       // GALERKIN SOLUTION 
				       double* h_dof_galerkin,
				       double* hu_dof_galerkin,
				       double* hv_dof_galerkin,
				       // TO COMPUTE LOW ORDER 
				       double* low_order_hnp1,
				       double* low_order_hunp1,
				       double* low_order_hvnp1,
				       double* dEV_minus_dL,
				       double cE,
				       int LUMPED_MASS_MATRIX,
				       int USE_EV_BASED_ON_GALERKIN,
				       double dt,
				       double mannings)
    {
      //TMP FOR FRICTION//
      double n2 = std::pow(mannings,2);
      double gamma=4./3;
      double xi=2.;

      double min_h = 1E10;
      ////////////////
      // CELL LOOPS //
      ////////////////
      // To compute: 
      //      * lumped_mass_matrix
      //      * Cell based CFL
      //      * velocity at quad points for other models 
      // init lumped mass matrix to zero
      if (recompute_lumped_mass_matrix==1)
	{
	  for (int i=0; i<numDOFsPerEqn; i++)
	    lumped_mass_matrix[i] = 0;
	  for(int eN=0;eN<nElements_global;eN++)
	    {
	      //declare local storage for element residual and initialize
	      register double element_lumped_mass_matrix[nDOF_test_element];
	      for (int i=0;i<nDOF_test_element;i++)
		element_lumped_mass_matrix[i]=0.0;
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
		  //get the solution at current time. This is to compute velocity for other models
		  ck.valFromDOF(h_dof,&h_l2g[eN_nDOF_trial_element],&h_trial_ref[k*nDOF_trial_element],h);
		  ck.valFromDOF(hu_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hu);
		  ck.valFromDOF(hv_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hv);
		  //precalculate test function products with integration weights
		  for (int j=0;j<nDOF_trial_element;j++)
		      h_test_dV[j] = h_test_ref[k*nDOF_trial_element+j]*dV;
		  //save velocity at quadrature points for other models to use
		  q_velocity[eN_k_nSpace+0] = 2*h/(h*h+std::pow(fmax(h,hEps),2))*hu;
		  q_velocity[eN_k_nSpace+1] = 2*h/(h*h+std::pow(fmax(h,hEps),2))*hv;
		  hnp1_at_quad_point[eN_k] = h;
		  hunp1_at_quad_point[eN_k] = hu;
		  hvnp1_at_quad_point[eN_k] = hv;
		  //update element residual. Part about the lumped mass matrix 
		  for(int i=0;i<nDOF_test_element;i++)
		    element_lumped_mass_matrix[i] += h_test_dV[i];
		}
	      // distribute
	      for(int i=0;i<nDOF_test_element;i++)
		{
		  register int eN_i=eN*nDOF_test_element+i;
		  int h_gi = h_l2g[eN_i]; //global i-th index for h
		  lumped_mass_matrix[h_gi]  += element_lumped_mass_matrix[i];	  
		}
	      // COMPUTE MIN h //
	      min_h = fmin(min_h,elementDiameter[eN]);
	    }
	} // end of computing lumped mass matrix
      
      ///////////////////////////////////////////
      // COMPUTE SMOOTHNESS INDICATOR and dLij //
      ///////////////////////////////////////////
      // Smoothness indicator is based on the solution. psi_i = psi_i(alpha_i); 
      // alpha_i = |sum(uj-ui)|/sum|uj-ui|
      int ij = 0;
      double max_cfl = 0.;
      register double psi[numDOFsPerEqn], dL[NNZ];
      for (int i=0; i<numDOFsPerEqn; i++)
	{
	  double alphai, alphai_numerator=0, alphai_denominator=0; // smoothness indicator of solution
	  double hi = h_dof_lstage[i]; // solution at time tn for the ith DOF
	  double hui = hu_dof_lstage[i]; 
	  double hvi = hv_dof_lstage[i]; 
	  double Zi = b_dof[i];
	  double dLii = 0.;

	  for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
	    { //loop in j (sparsity pattern)
	      int j = csrColumnOffsets_DofLoops[offset];
	      double hj = h_dof_lstage[j]; // solution at time tn for the jth DOF
	      double huj = hu_dof_lstage[j];
	      double hvj = hv_dof_lstage[j];
	      double Zj = b_dof[j];

	      //////////////////////////////////
	      // COMPUTE SMOOTHNESS INDICATOR //
	      //////////////////////////////////
	      alphai_numerator += (hj-hi);
	      alphai_denominator += std::abs(hj-hi);

	      ///////////////////////////
	      // COMPUTE STAR SOLUTION // hStar, huStar and hvStar
	      ///////////////////////////
	      double hStarij  = fmax(0., hi + Zi - fmax(Zi,Zj));
	      double huStarij = (hi <= hEps ? 0. : hui*hStarij/hi);
	      double hvStarij = (hi <= hEps ? 0. : hvi*hStarij/hi);
	      
	      double hStarji  = fmax(0., hj + Zj - fmax(Zi,Zj));
	      double huStarji = (hj <= hEps ? 0. : huj*hStarji/hj);
	      double hvStarji = (hj <= hEps ? 0. : hvj*hStarji/hj);

	      //////////////////
	      // COMPUTE dLij // 
	      //////////////////
	      // NOTE: I don't have anything in the entry dLij when i=j
	      if (i != j)
		{		 
		  // compute dLij
		  double cij_norm = sqrt(Cx[ij]*Cx[ij] + Cy[ij]*Cy[ij]);
		  double cji_norm = sqrt(CTx[ij]*CTx[ij] + CTy[ij]*CTy[ij]);
		  
		  double nxij = Cx[ij]/cij_norm, nyij = Cy[ij]/cij_norm;
		  double nxji = CTx[ij]/cji_norm, nyji = CTy[ij]/cji_norm;
		  
		  //dissipative terms for flux
		  dL[ij] = fmax(fmax(maxWaveSpeedSharpInitialGuess(g,nxij,nyij, 
								   hi,hui,hvi,
								   hStarji,huStarji,hvStarji,
								   hEps,false),
				     maxWaveSpeedSharpInitialGuess(g,nxij,nyij, 
								   hi,hui,hvi, 
								   hStarij,huStarij,hvStarij,
								   hEps,false))*cij_norm,
				fmax(maxWaveSpeedSharpInitialGuess(g,nxji,nyji, 
								   hj,huj,hvj, 
								   hStarij,huStarij,hvStarij,
								   hEps,false),
				     maxWaveSpeedSharpInitialGuess(g,nxji,nyji, 
								   hj,huj,hvj, 
								   hStarji,huStarji,hvStarji,
								   hEps,false))*cji_norm);  
		  /*
		  // for debugging
		  if (dL[ij] > 100)
		    {		      
		      std::cout << "dLij: " << dL[ij] << std::endl;
		      std::cout << "hi, hui, hvi: " 
				<< hi << ", " 
				<< hui << ", " 
				<< hvi << std::endl;
		      std::cout << "hj, huj, hvj: " 
				<< hj << ", " 
				<< huj << ", " 
				<< hvj << std::endl;
		      std::cout << "hStarij, huStarij, hvStarij: " 
				<< hStarij << ", " 
				<< huStarij << ", " 
				<< hvStarij << std::endl;
		      std::cout << "hStarji, huStarji, hvStarji: " 
				<< hStarji << ", " 
				<< huStarji << ", " 
				<< hvStarji << std::endl;
		      std::cout << "*************************" << std::endl;
		      std::cout << maxWaveSpeedSharpInitialGuess(g,nxij,nyij, 
								 hi,hui,hvi,
								 hStarji,huStarji,hvStarji,
								 hEps,true)
				<< std::endl;
		      std::cout << "*************************" << std::endl;
		      std::cout << maxWaveSpeedSharpInitialGuess(g,nxij,nyij, 
								 hi,hui,hvi, 
								 hStarij,huStarij,hvStarij,
								 hEps,true)
				<< std::endl;
		      std::cout << "*************************" << std::endl;
		      std::cout << maxWaveSpeedSharpInitialGuess(g,nxji,nyji, 
								 hj,huj,hvj, 
								 hStarij,huStarij,hvStarij,
								 hEps,true)
				<< std::endl;
		      std::cout << "*************************" << std::endl;
		      std::cout << maxWaveSpeedSharpInitialGuess(g,nxji,nyji, 
								 hj,huj,hvj, 
								 hStarji,huStarji,hvStarji,
								 hEps,true)
				<< std::endl;
		      std::cout << "*************************" << std::endl;
		      abort();
		      }
		  */
		  dLii -= dL[ij];
		}
	      //update ij
	      ij+=1;
	    }
	  double mi = lumped_mass_matrix[i];
	  // calculate edge based CFL
	  edge_based_cfl[i] = 2*fabs(dLii)/mi;
	  max_cfl = fmax(max_cfl,edge_based_cfl[i]);
	  if (hi < hEps)
	    alphai = 1.;
	  else
	    {
	      if (alphai_denominator==0)
		alphai = 1.;
	      else 
		alphai = std::abs(alphai_numerator)/alphai_denominator;
	    }
	  if (POWER_SMOOTHNESS_INDICATOR==0)
	    psi[i] = 1.0;
	  else
	    psi[i] = std::pow(alphai,POWER_SMOOTHNESS_INDICATOR); //NOTE: they use alpha^2 in the paper
	}
      //////////////////////
      // CHOOSE TIME STEP //
      //////////////////////
      //double dt = cfl_run/max_cfl;
      //double dt = 0.1*min_h;

      //////////////////
      // Loop on DOFs // to compute flux and dissipative terms
      //////////////////
      ij = 0;
      for (int i=0; i<numDOFsPerEqn; i++)
	{
	  double hi = h_dof_lstage[i];
	  double hui = hu_dof_lstage[i];
	  double hvi = hv_dof_lstage[i];	  
	  double Zi = b_dof[i];
	  
	  double ui = 2*hi/(hi*hi+std::pow(fmax(hi,hEps),2))*hui;
	  double vi = 2*hi/(hi*hi+std::pow(fmax(hi,hEps),2))*hvi;

	  // regularization of 1/hi 
	  double one_over_hiReg = 2*hi/(hi*hi+std::pow(fmax(hi,hEps),2));

	  double ith_flux_term1=0., ith_flux_term2=0., ith_flux_term3=0.;
	  double ith_dissipative_term1=0., ith_dissipative_term2=0., ith_dissipative_term3=0.;
	  double ith_well_balancing_term1=0., ith_well_balancing_term2=0., ith_well_balancing_term3=0.;
	  double aux1_to_compute_hnp1=0., aux2_to_compute_hnp1=0.;

	  // friction term lumped
	  double veli_norm = std::sqrt(ui*ui+vi*vi);
	  double mi = lumped_mass_matrix[i];
	  
	  double hi_to_the_gamma = std::pow(hi,gamma);
	  double ith_friction_term2 =  
	    veli_norm==0 ? 0. : 2*g*n2*hui*veli_norm*mi/(hi_to_the_gamma+fmax(hi_to_the_gamma,xi*g*n2*dt*veli_norm));
	  double ith_friction_term3 =  
	    veli_norm==0 ? 0. : 2*g*n2*hvi*veli_norm*mi/(hi_to_the_gamma+fmax(hi_to_the_gamma,xi*g*n2*dt*veli_norm));


	  // loop over the sparsity pattern of the i-th DOF
	  for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
	    {
	      int j = csrColumnOffsets_DofLoops[offset];
	      double hj = h_dof_lstage[j];
	      double huj = hu_dof_lstage[j];
	      double hvj = hv_dof_lstage[j];
	      double Zj = b_dof[j];
	      double uj = 2*hj/(hj*hj+std::pow(fmax(hj,hEps),2))*huj;
	      double vj = 2*hj/(hj*hj+std::pow(fmax(hj,hEps),2))*hvj;

	      // regularization of 1/hj
	      double one_over_hjReg = 2*hj/(hj*hj+std::pow(fmax(hj,hEps),2));

	      // Nodal projection of fluxes
	      ith_flux_term1 += huj*Cx[ij] + hvj*Cy[ij]; // f1*C
	      ith_flux_term2 += ( huj*huj*one_over_hjReg*Cx[ij] + huj*hvj*one_over_hjReg*Cy[ij] 
				  + g*hi*(hj+Zj)*Cx[ij] );
	      ith_flux_term3 += ( huj*hvj*one_over_hjReg*Cx[ij] + hvj*hvj*one_over_hjReg*Cy[ij] 
				  + g*hi*(hj+Zj)*Cy[ij] );

	      // Dissipative term
	      double muij = 0.;
	      if (i != j) // This is not necessary. See formula for ith_dissipative_terms
		{
		  ///////////////////////////
		  // COMPUTE STAR SOLUTION // hStar, huStar and hvStar
		  ///////////////////////////
		  double hStarij  = fmax(0., hi + Zi - fmax(Zi,Zj));
		  double huStarij = (hi <= hEps ? 0. : hui*hStarij/hi);
		  double hvStarij = (hi <= hEps ? 0. : hvi*hStarij/hi);
		  
		  double hStarji  = fmax(0., hj + Zj - fmax(Zi,Zj));
		  double huStarji = (hj <= hEps ? 0. : huj*hStarji/hj);
		  double hvStarji = (hj <= hEps ? 0. : hvj*hStarji/hj);
		  
		  //dissipative terms for well balancing
		  muij = fmax(fmax(0.,-(ui*Cx[ij] + vi*Cy[ij])),fmax(0,(uj*Cx[ij] + vj*Cy[ij])));
		  muij *= std::max(psi[i],psi[j]);
		  // compute dissipative terms for well balancing for second and third equations
		  ith_well_balancing_term1 += muij*(hj-hStarji-(hi-hStarij));
		  ith_well_balancing_term2 += muij*(huj-huStarji-(hui-huStarij));
		  ith_well_balancing_term3 += muij*(hvj-hvStarji-(hvi-hvStarij));

		  // compute dissipative terms for second and third equations
		  dL[ij] *= std::max(psi[i],psi[j]);
		  ith_dissipative_term1 += dL[ij]*(hStarji-hStarij);
		  ith_dissipative_term2 += dL[ij]*(huStarji-huStarij);
		  ith_dissipative_term3 += dL[ij]*(hvStarji-hvStarij);

		  // compute aux quantities for first equation 
		  // aux1 = Veli*Cii + sum_j[ muij + (dLij-muij)*hStarij/hi ]
		  // aux2 = sum_j[ (muij-Velj*Cij)*hj + (dLij-muij)*hStarji ] >= 0
		  aux1_to_compute_hnp1 += muij + (hi == 0. ? 0. : (dL[ij] - muij)*hStarij/hi);
		  aux2_to_compute_hnp1 += (muij - (uj*Cx[ij] + vj*Cy[ij]))*hj + (dL[ij] - muij)*hStarji;
		}
	      else // i==j
		{
		  aux1_to_compute_hnp1 += (ui*Cx[ij] + vi*Cy[ij])*hi;
		}
	      // update ij
	      ij+=1;
	    }
	  //if (1 - dt/mi*aux1_to_compute_hnp1 < 0)
	  //{
	  //  std::cout << "aux1 < 0: " << 1-dt/mi*aux1_to_compute_hnp1 << std::endl;
	  //  std::cout << "dt, cfl: " << dt << "\t" << edge_based_cfl[i] << std::endl;
	  //  abort();
	  //}

	  // compute residual
	  double hinp1 = hi*(1 - dt/mi*aux1_to_compute_hnp1) + dt/mi*aux2_to_compute_hnp1;
	  globalResidual[offset_h+stride_h*i] = fmax(0.,hinp1);	 
	  //globalResidual[offset_h+stride_h*i] = hinp1;	 
	  //globalResidual[offset_h+stride_h*i] = hi - dt/mi*(ith_flux_term1 - ith_dissipative_term1 - ith_well_balancing_term1);
	  globalResidual[offset_hu+stride_hu*i] = hui - dt/mi*(ith_flux_term2 - ith_dissipative_term2 - ith_well_balancing_term2 + ith_friction_term2);
	  globalResidual[offset_hv+stride_hv*i] = hvi - dt/mi*(ith_flux_term3 - ith_dissipative_term3 - ith_well_balancing_term3 + ith_friction_term3);
	}
    }

    void calculateResidual_second_order_NonFlatB_with_EV(//element
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
							 double* h_dof_old_old, 
							 double* hu_dof_old_old, 
							 double* hv_dof_old_old, 
							 double* h_dof_lstage, 
							 double* hu_dof_lstage, 
							 double* hv_dof_lstage, 
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
							 double* q_velocity_sge,
							 double* q_cfl,
							 double* q_numDiff_h, 
							 double* q_numDiff_hu, 
							 double* q_numDiff_hv, 
							 double* q_numDiff_h_last,
							 double* q_numDiff_hu_last, 
							 double* q_numDiff_hv_last,
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
							 double* edge_based_cfl,
							 double cfl_run,
							 double hEps,
							 int recompute_lumped_mass_matrix, 
							 // SAVE SOLUTION (mql)
							 double* hnp1_at_quad_point,
							 double* hunp1_at_quad_point,
							 double* hvnp1_at_quad_point,
							 // GALERKIN SOLUTION 
							 double* h_dof_galerkin,
							 double* hu_dof_galerkin,
							 double* hv_dof_galerkin, 
							 // TO COMPUTE LOW ORDER 
							 double* low_order_hnp1,
							 double* low_order_hunp1,
							 double* low_order_hvnp1,
							 double* dEV_minus_dL,
							 double cE,
							 int LUMPED_MASS_MATRIX,
							 int USE_EV_BASED_ON_GALERKIN,
							 double dt,
							 double mannings)
    {
      //FOR FRICTION//
      double n2 = std::pow(mannings,2);
      double gamma=4./3;
      double xi=2.;

      ////////////////
      // CELL LOOPS //
      ////////////////
      // To compute: 
      //      * Entropy residual
      //      * lumped_mass_matrix
      //      * Cell based CFL
      //      * velocity at quad points for other models 
      // create global entropy residual vectors for each equation
      register double global_entropy_residual[numDOFsPerEqn]; 
      
      // init lumped mass matrix and ent residual vectors to zero
      for (int i=0; i<numDOFsPerEqn; i++)
	{
	  lumped_mass_matrix[i] = 0;
	  global_entropy_residual[i]=0;
	}
      for(int eN=0;eN<nElements_global;eN++)
	{
	  //declare local storage for element residual and initialize
	  register double 
	    element_lumped_mass_matrix[nDOF_test_element],
	    elementResidual_h[nDOF_test_element],
	    elementResidual_hu[nDOF_test_element],
	    elementResidual_hv[nDOF_test_element],
	    element_entropy_residual[nDOF_test_element]; 
	    
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      element_lumped_mass_matrix[i]=0.0;
	      elementResidual_h[i]=0.0;
	      elementResidual_hu[i]=0.0;
	      elementResidual_hv[i]=0.0;
	      element_entropy_residual[i] =0;
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
		hnm1=0.0,hunm1=0.0,hvnm1=0.0, // solution at tnm1
		hn=0.0,hun=0.0,hvn=0.0, // solution at tn
		h_lstage=0.0,hu_lstage=0.0,hv_lstage=0.0, // solution at lstage
		hG=0.0,huG=0.0,hvG=0.0, // galerkin solution
		grad_hn[nSpace],grad_hun[nSpace],grad_hvn[nSpace],grad_z[nSpace],		
		jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace],
		h_test_dV[nDOF_trial_element],
		h_grad_trial[nDOF_trial_element*nSpace],vel_grad_trial[nDOF_trial_element*nSpace],
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
	      //get the trial function gradients
	      ck.gradTrialFromRef(&h_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,h_grad_trial);
	      ck.gradTrialFromRef(&vel_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,vel_grad_trial);
	      //get the solution at current time. This is to compute velocity for other models
	      ck.valFromDOF(h_dof,&h_l2g[eN_nDOF_trial_element],&h_trial_ref[k*nDOF_trial_element],h);
	      ck.valFromDOF(hu_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hu);
	      ck.valFromDOF(hv_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hv);
	      // calculate the solution at time tnm1
	      ck.valFromDOF(h_dof_old_old,&h_l2g[eN_nDOF_trial_element],&h_trial_ref[k*nDOF_trial_element],hnm1);
	      ck.valFromDOF(hu_dof_old_old,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hunm1);
	      ck.valFromDOF(hv_dof_old_old,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hvnm1);
	      // get the solution at time tn for ENTROPY VISCOSITY
	      ck.valFromDOF(h_dof_old,&h_l2g[eN_nDOF_trial_element],&h_trial_ref[k*nDOF_trial_element],hn);
	      ck.valFromDOF(hu_dof_old,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hun);
	      ck.valFromDOF(hv_dof_old,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hvn);
	      // get the solution at the lstage
	      ck.valFromDOF(h_dof_lstage,&h_l2g[eN_nDOF_trial_element],&h_trial_ref[k*nDOF_trial_element],h_lstage);
	      ck.valFromDOF(hu_dof_lstage,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hu_lstage);
	      ck.valFromDOF(hv_dof_lstage,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hv_lstage);
	      // get the galerkin solution
	      if (USE_EV_BASED_ON_GALERKIN==1)
		{
		  ck.valFromDOF(h_dof_galerkin,&h_l2g[eN_nDOF_trial_element],&h_trial_ref[k*nDOF_trial_element],hG);
		  ck.valFromDOF(hu_dof_galerkin,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],huG);
		  ck.valFromDOF(hv_dof_galerkin,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hvG);
		}
	      // calculate grad of solution at tn for ENTROPY VISCOSITY
	      ck.gradFromDOF(h_dof_old,&h_l2g[eN_nDOF_trial_element],h_grad_trial,grad_hn);
	      ck.gradFromDOF(hu_dof_old,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_hun);
	      ck.gradFromDOF(hv_dof_old,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_hvn);
	      ck.gradFromDOF(b_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_z);
	      // calculate cell based CFL
	      calculateCFL(elementDiameter[eN],
			   g,
			   hn,
			   hun,
			   hvn,
			   hEps,
			   q_cfl[eN_k]);
	      //precalculate test function products with integration weights
	      for (int j=0;j<nDOF_trial_element;j++)
		h_test_dV[j] = h_test_ref[k*nDOF_trial_element+j]*dV;
	      //save velocity at quadrature points for other models to use
	      q_velocity[eN_k_nSpace+0] = 2*h/(h*h+std::pow(fmax(h,hEps),2))*hu;
	      q_velocity[eN_k_nSpace+1] = 2*h/(h*h+std::pow(fmax(h,hEps),2))*hv;
	      hnp1_at_quad_point[eN_k] = h;
	      hunp1_at_quad_point[eN_k] = hu;
	      hvnp1_at_quad_point[eN_k] = hv;
	      //update element residual. Part about the lumped mass matrix 
	      //double one_over_hnm1Reg = 2*hnm1/(hnm1*hnm1+std::pow(fmax(hnm1,hEps),2));
	      double one_over_hnReg = 2*hn/(hn*hn+std::pow(fmax(hn,hEps),2));
	      //double one_over_hGReg = 2*hG/(hG*hG+std::pow(fmax(hG,hEps),2));
	      //double dxf1 = grad_hun[0];
	      //double dyf1 = grad_hvn[1];
	      //double dxf2 = (2*hn*hun*grad_hun[0] - hun*hun*grad_hn[0])*std::pow(one_over_hnReg,2)+g*hn*grad_hn[0];
	      //double dyf2 = (hn*(hun*grad_hvn[1]+hvn*grad_hun[1])-hun*hvn*grad_hn[1])*std::pow(one_over_hnReg,2);
	      //double dxf3 = (hn*(hun*grad_hvn[0]+hvn*grad_hun[0])-hun*hvn*grad_hn[0])*std::pow(one_over_hnReg,2);
	      //double dyf3 = (2*hn*hvn*grad_hvn[1] - hvn*hvn*grad_hn[1])*std::pow(one_over_hnReg,2)+g*hn*grad_hn[1];

	      // bathymetry source terms
	      double b_dot_gradx_z = g*hn*grad_z[0];
	      double b_dot_grady_z = g*hn*grad_z[1];
	      // fprime_dot_grad(u)
	      double un = 2*hn/(hn*hn+std::pow(fmax(hn,hEps),2))*hun;
	      double vn = 2*hn/(hn*hn+std::pow(fmax(hn,hEps),2))*hvn;

	      double fp_dot_grad_u1 = grad_hun[0] + grad_hvn[1];
	      double fp_dot_grad_u2 = grad_hn[0]*(-un*un+g*hn)+2*un*grad_hun[0]-un*vn*grad_hn[1]+vn*grad_hun[1]+un*grad_hvn[1];
	      double fp_dot_grad_u3 = -un*vn*grad_hn[0]+vn*grad_hun[0]+un*grad_hvn[0]+grad_hn[1]*(-vn*vn+g*hn)+2*vn*grad_hvn[1];

	      double entropy_residual1 = ((hn-hnm1)/dt   + (fp_dot_grad_u1))                 * DENTROPY_DH (g,hn,hun,hvn,one_over_hnReg);
	      double entropy_residual2 = ((hun-hunm1)/dt + (fp_dot_grad_u2 + b_dot_gradx_z)) * DENTROPY_DHU(g,hn,hun,hvn,one_over_hnReg);
	      double entropy_residual3 = ((hvn-hvnm1)/dt + (fp_dot_grad_u3 + b_dot_grady_z)) * DENTROPY_DHV(g,hn,hun,hvn,one_over_hnReg);
	      if (USE_EV_BASED_ON_GALERKIN==1)
		{
		  entropy_residual1 = ((hG-hn)/dt   + (fp_dot_grad_u1))                 * DENTROPY_DH (g,hn,hun,hvn,one_over_hnReg);
		  entropy_residual2 = ((huG-hun)/dt + (fp_dot_grad_u2 + b_dot_gradx_z)) * DENTROPY_DHU(g,hn,hun,hvn,one_over_hnReg);
		  entropy_residual3 = ((hvG-hvn)/dt + (fp_dot_grad_u3 + b_dot_grady_z)) * DENTROPY_DHV(g,hn,hun,hvn,one_over_hnReg);
		}

	      for(int i=0;i<nDOF_test_element;i++)
		{
		  // lumped mass matrix
		  element_lumped_mass_matrix[i] += h_test_dV[i];
		  element_entropy_residual[i] += 
		    (entropy_residual1+entropy_residual2+entropy_residual3)*h_test_dV[i];

		  // compute time derivative part of global residual. NOTE: no lumping
		  elementResidual_h[i]  += (h-h_lstage)*h_test_dV[i];
		  elementResidual_hu[i] += (hu-hu_lstage)*h_test_dV[i];
		  elementResidual_hv[i] += (hv-hv_lstage)*h_test_dV[i];
		}
	    }
	  // distribute
	  for(int i=0;i<nDOF_test_element;i++)
	    {
	      register int eN_i=eN*nDOF_test_element+i;
	      int h_gi = h_l2g[eN_i]; //global i-th index for h
	      int vel_gi = vel_l2g[eN_i]; //global i-th index for velocities 
	      
	      lumped_mass_matrix[h_gi]  += element_lumped_mass_matrix[i];
	      // distribute element entropy residual to global vectors
	      global_entropy_residual[h_gi] += element_entropy_residual[i];
	      
	      // distribute time derivative to global residual
	      globalResidual[offset_h+stride_h*h_gi]  += elementResidual_h[i];
	      globalResidual[offset_hu+stride_hu*vel_gi] += elementResidual_hu[i];
	      globalResidual[offset_hv+stride_hv*vel_gi] += elementResidual_hv[i];
	    }
	}
     

      ////////////////////////////////////
      // COMPUTE FAKE GALERKIN SOLUTION //
      ////////////////////////////////////
      int ij = 0;
      for (int i=0; i<numDOFsPerEqn; i++)
	{
	  double hni = h_dof_lstage[i];
	  double huni = hu_dof_lstage[i];
	  double hvni = hv_dof_lstage[i];
	  double ith_flux_term1=0., ith_flux_term2=0., ith_flux_term3=0.;
	  
	  for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
	    {
	      int j = csrColumnOffsets_DofLoops[offset];
	      double hnj = h_dof_lstage[j];
	      double hunj = hu_dof_lstage[j];
	      double hvnj = hv_dof_lstage[j];
	      double Zj = b_dof[j];
	      
	      // regularization of 1/hj
	      double one_over_hnjReg = 2*hnj/(hnj*hnj+std::pow(fmax(hnj,hEps),2));

	      // Nodal projection of fluxes
	      ith_flux_term1 += hunj*Cx[ij] + hvnj*Cy[ij]; // f1*C
	      ith_flux_term2 += ( hunj*hunj*one_over_hnjReg*Cx[ij] + hunj*hvnj*one_over_hnjReg*Cy[ij] 
				  + g*hni*(hnj+Zj)*Cx[ij] );
	      ith_flux_term3 += ( hunj*hvnj*one_over_hnjReg*Cx[ij] + hvnj*hvnj*one_over_hnjReg*Cy[ij] 
				  + g*hni*(hnj+Zj)*Cy[ij] );
	    }
	  double mi = lumped_mass_matrix[i];
	  h_dof_galerkin[i] = hni - dt/mi*ith_flux_term1;
	  hu_dof_galerkin[i] = huni - dt/mi*ith_flux_term2;
	  hv_dof_galerkin[i] = hvni - dt/mi*ith_flux_term3;
	}
 
      ////////////////////////////////
      // COMPUTE ENTROPY AT ith DOF //
      ////////////////////////////////
      // compute entropy (defined as eta) corresponding to ith node
      register double eta[numDOFsPerEqn];
      register double etaG[numDOFsPerEqn];
      for (int i=0; i<numDOFsPerEqn; i++)
	{
	  // COMPUTE ENTROPY BASED ON OLD STAGE
	  //double hin = h_dof_lstage[i];
	  //double one_over_hinReg = 2*hin/(hin*hin+std::pow(fmax(hin,hEps),2));
	  //eta[i] = ENTROPY(g,hin,hu_dof_lstage[i],hv_dof_lstage[i],one_over_hinReg);

	  // COMPUTE ENTROPY BASED ON OLD SOLUTION 
	  double hni = h_dof_lstage[i]; 
	  double huni = h_dof_lstage[i]; 
	  double hvni = h_dof_lstage[i]; 
	  double one_over_hniReg = 2*hni/(hni*hni+std::pow(fmax(hni,hEps),2));
	  eta[i] = ENTROPY(g,hni,huni,hvni,one_over_hniReg);

	  // COMPUTE ENTROPY BASED ON FAKE GALERKIN
	  double hGi = h_dof_galerkin[i]; 
	  double huGi = hu_dof_galerkin[i]; 
	  double hvGi = hv_dof_galerkin[i]; 
	  double one_over_hGiReg = 2*hGi/(hGi*hGi+std::pow(fmax(hGi,hEps),2));
	  etaG[i] = ENTROPY(g,hGi,huGi,hvGi,one_over_hGiReg);

	  //////////////////////////
	  // NEW ENTROPY RESIDUAL //
	  //////////////////////////
	  double ith_flux_term1=0., ith_flux_term2=0., ith_flux_term3=0.;
	  double eta_prime1 = DENTROPY_DH  (g,hni,huni,hvni,one_over_hniReg); 
	  double eta_prime2 = DENTROPY_DHU (g,hni,huni,hvni,one_over_hniReg); 
	  double eta_prime3 = DENTROPY_DHV (g,hni,huni,hvni,one_over_hniReg); 

	  // loop over the sparsity pattern of the i-th DOF
	  for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
	    {
	      int j = csrColumnOffsets_DofLoops[offset];
	      double hnj = h_dof_lstage[j];
	      double hunj = hu_dof_lstage[j];
	      double hvnj = hv_dof_lstage[j];
	      double Zj = b_dof[j];

	      // regularization of 1/hj
	      double one_over_hnjReg = 2*hnj/(hnj*hnj+std::pow(fmax(hnj,hEps),2));

	      // Nodal projection of fluxes
	      ith_flux_term1 += hunj*Cx[ij] + hvnj*Cy[ij]; // f1*C
	      ith_flux_term2 += ( hunj*hunj*one_over_hnjReg*Cx[ij] + hunj*hvnj*one_over_hnjReg*Cy[ij] 
				  + g*hni*(hnj+Zj)*Cx[ij] );
	      ith_flux_term3 += ( hunj*hvnj*one_over_hnjReg*Cx[ij] + hvnj*hvnj*one_over_hnjReg*Cy[ij] 
				  + g*hni*(hnj+Zj)*Cy[ij] );	      
	    }
	  double mi = lumped_mass_matrix[i];
	  // option 1
	  //global_entropy_residual[i] = 
	  //(etaG[i] - eta[i])*mi/dt 
	  //+ ith_flux_term1*eta_prime1 + ith_flux_term2*eta_prime2 + ith_flux_term3*eta_prime3;
	  // option 2
	  global_entropy_residual[i] = ( Cx[ij]*ENTROPY_FLUX1(g,hni,huni,hvni,one_over_hniReg) + 
					 Cy[ij]*ENTROPY_FLUX2(g,hni,huni,hvni,one_over_hniReg) )
	    -(ith_flux_term1*eta_prime1 + ith_flux_term2*eta_prime2 + ith_flux_term3*eta_prime3);

	  //double aux = ( Cx[ij]*ENTROPY_FLUX1(g,hni,huni,hvni,one_over_hniReg) + 
	  //	 Cy[ij]*ENTROPY_FLUX2(g,hni,huni,hvni,one_over_hniReg) )
	  //-(ith_flux_term1*eta_prime1 + ith_flux_term2*eta_prime2 + ith_flux_term3*eta_prime3);

	  //std::cout << global_entropy_residual[i] - aux << std::endl;
	}

      //////////////////////////////////////////////////////////////
      // COMPUTE SMOOTHNESS INDICATOR, dLij and etaMin and etaMax //
      //////////////////////////////////////////////////////////////
      // Smoothness indicator is based on the solution. psi_i = psi_i(alpha_i); 
      // alpha_i = |sum(uj-ui)|/sum|uj-ui|
      ij = 0;
      double max_cfl = 0.;
      register double psi[numDOFsPerEqn], dL[NNZ], etaMax[numDOFsPerEqn], etaMin[numDOFsPerEqn];
      for (int i=0; i<numDOFsPerEqn; i++)
	{
	  double alphai, alphai_numerator=0, alphai_denominator=0; // smoothness indicator of solution
	  double hi = h_dof_lstage[i]; // solution at time tn for the ith DOF
	  double hui = hu_dof_lstage[i]; 
	  double hvi = hv_dof_lstage[i]; 
	  double Zi = b_dof[i];
	  double dLii = 0.;

	  etaMax[i] = std::abs(eta[i]);
	  etaMin[i] = std::abs(eta[i]);
	  for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
	    { //loop in j (sparsity pattern)
	      int j = csrColumnOffsets_DofLoops[offset];
	      double hj = h_dof_lstage[j]; // solution at time tn for the jth DOF
	      double huj = hu_dof_lstage[j];
	      double hvj = hv_dof_lstage[j];
	      double Zj = b_dof[j];

	      //////////////////////////////////
	      // COMPUTE SMOOTHNESS INDICATOR //
	      //////////////////////////////////
	      alphai_numerator += (hj-hi);
	      alphai_denominator += std::abs(hj-hi);

	      ///////////////////////////
	      // COMPUTE STAR SOLUTION // hStar, huStar and hvStar
	      ///////////////////////////
	      double hStarij  = fmax(0., hi + Zi - fmax(Zi,Zj));
	      double huStarij = (hi <= hEps ? 0. : hui*hStarij/hi);
	      double hvStarij = (hi <= hEps ? 0. : hvi*hStarij/hi);
	      
	      double hStarji  = fmax(0., hj + Zj - fmax(Zi,Zj));
	      double huStarji = (hj <= hEps ? 0. : huj*hStarji/hj);
	      double hvStarji = (hj <= hEps ? 0. : hvj*hStarji/hj);

	      //////////////////
	      // COMPUTE dLij // 
	      //////////////////
	      // NOTE: I don't have anything in the entry dLij when i=j
	      if (i != j)
		{		 
		  // compute dLij
		  double cij_norm = sqrt(Cx[ij]*Cx[ij] + Cy[ij]*Cy[ij]);
		  double cji_norm = sqrt(CTx[ij]*CTx[ij] + CTy[ij]*CTy[ij]);
		  
		  double nxij = Cx[ij]/cij_norm, nyij = Cy[ij]/cij_norm;
		  double nxji = CTx[ij]/cji_norm, nyji = CTy[ij]/cji_norm;
		  
		  //dissipative terms for flux
		  dL[ij] = fmax(fmax(maxWaveSpeedSharpInitialGuess(g,nxij,nyij, 
								   hi,hui,hvi,
								   hStarji,huStarji,hvStarji,
								   hEps,false),
				     maxWaveSpeedSharpInitialGuess(g,nxij,nyij, 
								   hi,hui,hvi, 
								   hStarij,huStarij,hvStarij,
								   hEps,false))*cij_norm,
				fmax(maxWaveSpeedSharpInitialGuess(g,nxji,nyji, 
								   hj,huj,hvj, 
								   hStarij,huStarij,hvStarij,
								   hEps,false),
				     maxWaveSpeedSharpInitialGuess(g,nxji,nyji, 
								   hj,huj,hvj, 
								   hStarji,huStarji,hvStarji,
								   hEps,false))*cji_norm);  
		  /*
		  if (dL[ij] > 100) //TMP
		    {
		      std::cout << dL[ij] << std::endl;

		      std::cout << maxWaveSpeedSharpInitialGuess(g,nxij,nyij, 
								 hi,hui,hvi,
								 hStarji,huStarji,hvStarji,
								 hEps,false)
				<< std::endl;
		      std::cout << maxWaveSpeedSharpInitialGuess(g,nxij,nyij, 
								 hi,hui,hvi, 
								 hStarij,huStarij,hvStarij,
								 hEps,false)
				<< std::endl;
		      std::cout << maxWaveSpeedSharpInitialGuess(g,nxji,nyji, 
								 hj,huj,hvj, 
								 hStarij,huStarij,hvStarij,
								 hEps,false)
				<< std::endl;
		      
		      std::cout << hj << "\t"
				<< huj << "\t" 
				<< hvj
				<< std::endl;
		      std::cout << hStarji << "\t"
				<< huStarji << "\t" 
				<< hvStarji
				<< std::endl;
		      std::cout << maxWaveSpeedSharpInitialGuess(g,nxji,nyji, 
								 hj,huj,hvj, 
								 hStarji,huStarji,hvStarji,
								 hEps,true)
				<< std::endl;

	      

		      abort();
		    }
		  */
		  dLii -= dL[ij];

		  // compute entropy max and min
		  etaMax[i] = fmax(etaMax[i],std::abs(eta[j]));
    		  etaMin[i] = fmin(etaMin[i],std::abs(eta[j]));
		}
	      //update ij
	      ij+=1;
	    }
	  double mi = lumped_mass_matrix[i];
	  // calculate edge based CFL
	  edge_based_cfl[i] = 2*fabs(dLii)/mi;
	  max_cfl = fmax(max_cfl,edge_based_cfl[i]);
	  if (hi < hEps)
	    alphai = 1.;
	  else
	    {
	      if (alphai_denominator==0)
		alphai = 1.;
	      else 
		alphai = std::abs(alphai_numerator)/alphai_denominator;
	    }
	  if (POWER_SMOOTHNESS_INDICATOR==0)
	    psi[i] = 1.0;
	  else
	    psi[i] = std::pow(alphai,POWER_SMOOTHNESS_INDICATOR); //NOTE: they use alpha^2 in the paper
	}
      //////////////////////
      // CHOOSE TIME STEP //
      //////////////////////
      //double dt = cfl_run/max_cfl;
      //double dt = 0.1*min_h;

      //////////////////
      // Loop on DOFs // to compute flux and dissipative terms
      //////////////////
      ij = 0;
      for (int i=0; i<numDOFsPerEqn; i++)
	{
	  double hi = h_dof_lstage[i];
	  double hui = hu_dof_lstage[i];
	  double hvi = hv_dof_lstage[i];	  
	  double Zi = b_dof[i];
	  
	  double ui = 2*hi/(hi*hi+std::pow(fmax(hi,hEps),2))*hui;
	  double vi = 2*hi/(hi*hi+std::pow(fmax(hi,hEps),2))*hvi;


	  // entropy normalization factor 
	  double one_over_entNormFactori = etaMax[i] == etaMin[i] ? 0. : 1./(etaMax[i]-etaMin[i]);

	  // regularization of 1/hi 
	  double one_over_hiReg = 2*hi/(hi*hi+std::pow(fmax(hi,hEps),2));

	  double ith_flux_term1=0., ith_flux_term2=0., ith_flux_term3=0.;
	  double ith_dissipative_low_order_term1=0., ith_dissipative_low_order_term2=0., ith_dissipative_low_order_term3=0.,
	    ith_dissipative_term1=0., ith_dissipative_term2=0., ith_dissipative_term3=0.;
	  double ith_well_balancing_term1=0., ith_well_balancing_term2=0., ith_well_balancing_term3=0.;
	  double aux1_to_compute_hnp1=0., aux2_to_compute_hnp1=0.;


	  // friction term lumped
	  double veli_norm = std::sqrt(ui*ui+vi*vi);
	  double mi = lumped_mass_matrix[i];
	  
	  double hi_to_the_gamma = std::pow(hi,gamma);
	  double ith_friction_term2 =  
	    veli_norm==0 ? 0. : 2*g*n2*hui*veli_norm*mi/(hi_to_the_gamma+fmax(hi_to_the_gamma,xi*g*n2*dt*veli_norm));
	  double ith_friction_term3 =  
	    veli_norm==0 ? 0. : 2*g*n2*hvi*veli_norm*mi/(hi_to_the_gamma+fmax(hi_to_the_gamma,xi*g*n2*dt*veli_norm));
	  
	  // loop over the sparsity pattern of the i-th DOF
	  for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
	    {
	      int j = csrColumnOffsets_DofLoops[offset];
	      double hj = h_dof_lstage[j];
	      double huj = hu_dof_lstage[j];
	      double hvj = hv_dof_lstage[j];
	      double Zj = b_dof[j];
	      double uj = 2*hj/(hj*hj+std::pow(fmax(hj,hEps),2))*huj;
	      double vj = 2*hj/(hj*hj+std::pow(fmax(hj,hEps),2))*hvj;

	      // regularization of 1/hj
	      double one_over_hjReg = 2*hj/(hj*hj+std::pow(fmax(hj,hEps),2));

	      // Nodal projection of fluxes
	      ith_flux_term1 += huj*Cx[ij] + hvj*Cy[ij]; // f1*C
	      ith_flux_term2 += ( huj*huj*one_over_hjReg*Cx[ij] + huj*hvj*one_over_hjReg*Cy[ij] 
				  + g*hi*(hj+Zj)*Cx[ij] );
	      ith_flux_term3 += ( huj*hvj*one_over_hjReg*Cx[ij] + hvj*hvj*one_over_hjReg*Cy[ij] 
				  + g*hi*(hj+Zj)*Cy[ij] );

	      // Dissipative term
	      double muij = 0.;

	      if (i != j) // This is not necessary. See formula for ith_dissipative_terms
		{
		  ///////////////////////////
		  // COMPUTE STAR SOLUTION // hStar, huStar and hvStar
		  ///////////////////////////
		  double hStarij  = fmax(0., hi + Zi - fmax(Zi,Zj));
		  double huStarij = (hi <= hEps ? 0. : hui*hStarij/hi);
		  double hvStarij = (hi <= hEps ? 0. : hvi*hStarij/hi);
		  
		  double hStarji  = fmax(0., hj + Zj - fmax(Zi,Zj));
		  double huStarji = (hj <= hEps ? 0. : huj*hStarji/hj);
		  double hvStarji = (hj <= hEps ? 0. : hvj*hStarji/hj);
		  
		  //dissipative terms for well balancing
		  muij = fmax(fmax(0.,-(ui*Cx[ij] + vi*Cy[ij])),fmax(0,(uj*Cx[ij] + vj*Cy[ij])));
		  muij *= std::max(psi[i],psi[j]);
		  // compute dissipative terms for well balancing for second and third equations
		  ith_well_balancing_term1 += muij*(hj-hStarji-(hi-hStarij));
		  ith_well_balancing_term2 += muij*(huj-huStarji-(hui-huStarij));
		  ith_well_balancing_term3 += muij*(hvj-hvStarji-(hvi-hvStarij));

		  // compute dissipative terms for second and third equations
		  dL[ij] *= std::max(psi[i],psi[j]);
		  // compute dLij*(uStarji-uStarij)
		  double ith_dissipative_low_order_term1_ij = dL[ij]*(hStarji-hStarij);
		  double ith_dissipative_low_order_term2_ij = dL[ij]*(huStarji-huStarij);
		  double ith_dissipative_low_order_term3_ij = dL[ij]*(hvStarji-hvStarij);
		  // sum dissipative matrix contribution to global vectors
		  ith_dissipative_low_order_term1 += ith_dissipative_low_order_term1_ij;
		  ith_dissipative_low_order_term2 += ith_dissipative_low_order_term2_ij;
		  ith_dissipative_low_order_term3 += ith_dissipative_low_order_term3_ij;

		  // compute entropy residual for each component 
		  double one_over_entNormFactorj = etaMax[j] == etaMin[j] ? 0. : 1./(etaMax[j]-etaMin[j]);
		  double dEVij = cE*fmax(std::abs(global_entropy_residual[i])*one_over_entNormFactori,
					 std::abs(global_entropy_residual[j])*one_over_entNormFactorj);
		  double ith_dissipative_term1_ij = fmin(dL[ij],dEVij)*(hStarji-hStarij);
		  double ith_dissipative_term2_ij = fmin(dL[ij],dEVij)*(huStarji-huStarij);
		  double ith_dissipative_term3_ij = fmin(dL[ij],dEVij)*(hvStarji-hvStarij);

		  if (cE >= 1000) //HACK to switch between EV or not
		    {
		      ith_dissipative_term1_ij = dL[ij]*(hStarji-hStarij);
		      ith_dissipative_term2_ij = dL[ij]*(huStarji-huStarij);
		      ith_dissipative_term3_ij = dL[ij]*(hvStarji-hvStarij);
		    }
		  // sum high order dissipative matrix contribution to global vectors
		  ith_dissipative_term1 += ith_dissipative_term1_ij;
		  ith_dissipative_term2 += ith_dissipative_term2_ij;
		  ith_dissipative_term3 += ith_dissipative_term3_ij;

		  // compute dEV_minus_dL
		  dEV_minus_dL[ij] = fmin(dL[ij],dEVij) - dL[ij];

		  // compute aux quantities for first equation 
		  // aux1 = Veli*Cii + sum_j[ muij + (dLij-muij)*hStarij/hi ]
		  // aux2 = sum_j[ (muij-Velj*Cij)*hj + (dLij-muij)*hStarji ] >= 0
		  aux1_to_compute_hnp1 += muij + (hi == 0. ? 0. : (dL[ij] - muij)*hStarij/hi);
		  aux2_to_compute_hnp1 += (muij - (uj*Cx[ij] + vj*Cy[ij]))*hj + (dL[ij] - muij)*hStarji;
		}
	      else // i==j
		{
		  aux1_to_compute_hnp1 += (ui*Cx[ij] + vi*Cy[ij])*hi;
		  dEV_minus_dL[ij]=0.; //Not true but the prod of this times Uj-Ui will be zero
		}
	      // update ij
	      ij+=1;
	    }
	  //double mi = lumped_mass_matrix[i];

	  // Compute low order solution: lumped mass matrix and low order dissipative matrix
	  low_order_hnp1[i]  = hi*(1-dt/mi*aux1_to_compute_hnp1) + dt/mi*aux2_to_compute_hnp1;
	  //low_order_hnp1[i]  = hi  - dt/mi*(ith_flux_term1 - ith_dissipative_low_order_term1);
	  low_order_hunp1[i] = hui - dt/mi*(ith_flux_term2 - ith_dissipative_low_order_term2 - ith_well_balancing_term2 + ith_friction_term2);
	  low_order_hvnp1[i] = hvi - dt/mi*(ith_flux_term3 - ith_dissipative_low_order_term3 - ith_well_balancing_term3 + ith_friction_term3); 

	  if (low_order_hnp1[i] < 0.)
	    {
	      std::cout << "low order hnp1: " << low_order_hnp1[i] << std::endl;
	      std::cout << hi << "\t"
			<< 1-dt/mi*aux1_to_compute_hnp1 << "\t" 
			<< dt/mi*aux2_to_compute_hnp1 
			<< std::endl;
	    }

	  // distribute local to global flux and dissipative contributions
	  if (LUMPED_MASS_MATRIX==1)
	    {
	      globalResidual[offset_h+stride_h*i]   = hi  - dt/mi*(ith_flux_term1 - ith_dissipative_term1 - ith_well_balancing_term1);
	      globalResidual[offset_hu+stride_hu*i] = hui - dt/mi*(ith_flux_term2 - ith_dissipative_term2 - ith_well_balancing_term2 + ith_friction_term2);
	      globalResidual[offset_hv+stride_hv*i] = hvi - dt/mi*(ith_flux_term3 - ith_dissipative_term3 - ith_well_balancing_term3 + ith_friction_term3);
	    }
	  else
	    {
	      globalResidual[offset_h+stride_h*i]   += dt*(ith_flux_term1 - ith_dissipative_term1 - ith_well_balancing_term1);
	      globalResidual[offset_hu+stride_hu*i] += dt*(ith_flux_term2 - ith_dissipative_term2 - ith_well_balancing_term2 + ith_friction_term2);
	      globalResidual[offset_hv+stride_hv*i] += dt*(ith_flux_term3 - ith_dissipative_term3 - ith_well_balancing_term3 + ith_friction_term3);
	    }

	}
    }

    void calculateResidual_galerkin(//element
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
				    double* h_dof_old_old, 
				    double* hu_dof_old_old, 
				    double* hv_dof_old_old, 
				    double* h_dof_lstage, 
				    double* hu_dof_lstage, 
				    double* hv_dof_lstage, 
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
				    double* q_velocity_sge,
				    double* q_cfl,
				    double* q_numDiff_h, 
				    double* q_numDiff_hu, 
				    double* q_numDiff_hv, 
				    double* q_numDiff_h_last,
				    double* q_numDiff_hu_last, 
				    double* q_numDiff_hv_last,
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
				    double* edge_based_cfl,
				    double cfl_run,
				    double hEps,
				    int recompute_lumped_mass_matrix, 
				    // SAVE SOLUTION (mql)
				    double* hnp1_at_quad_point,
				    double* hunp1_at_quad_point,
				    double* hvnp1_at_quad_point,
				    // GALERKIN SOLUTION 
				    double* h_dof_galerkin,
				    double* hu_dof_galerkin,
				    double* hv_dof_galerkin,
				    // TO COMPUTE LOW ORDER 
				    double* low_order_hnp1,
				    double* low_order_hunp1,
				    double* low_order_hvnp1,
				    double* dEV_minus_dL,
				    double cE,
				    int LUMPED_MASS_MATRIX,
				    int USE_EV_BASED_ON_GALERKIN,
				    double dt,
				    double mannings)
    {
      ////////////////
      // CELL LOOPS //
      ////////////////
      // To compute residual of galerkin solution
      for(int eN=0;eN<nElements_global;eN++)
	{
	  //declare local storage for element residual and initialize
	  register double elementResidual_h[nDOF_test_element], 
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
		hn=0.0,hun=0.0,hvn=0.0, // solution at current time
		grad_hn[nSpace],grad_hun[nSpace],grad_hvn[nSpace],grad_z[nSpace],    
		jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace],
		h_test_dV[nDOF_trial_element],
		h_grad_trial[nDOF_trial_element*nSpace],vel_grad_trial[nDOF_trial_element*nSpace],
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
	      //get the trial function gradients
	      ck.gradTrialFromRef(&h_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,h_grad_trial);
	      ck.gradTrialFromRef(&vel_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,vel_grad_trial);
	      //get the solution at current time. This is to compute velocity for other models
	      ck.valFromDOF(h_dof,&h_l2g[eN_nDOF_trial_element],&h_trial_ref[k*nDOF_trial_element],h);
	      ck.valFromDOF(hu_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hu);
	      ck.valFromDOF(hv_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hv);
	      // get the solution at time tn
	      ck.valFromDOF(h_dof_lstage,&h_l2g[eN_nDOF_trial_element],&h_trial_ref[k*nDOF_trial_element],hn);
	      ck.valFromDOF(hu_dof_lstage,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hun);
	      ck.valFromDOF(hv_dof_lstage,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hvn);
	      // calculate grad of solution at tn at quadrature points
	      ck.gradFromDOF(h_dof_lstage,&h_l2g[eN_nDOF_trial_element],h_grad_trial,grad_hn);
	      ck.gradFromDOF(hu_dof_lstage,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_hun);
	      ck.gradFromDOF(hv_dof_lstage,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_hvn);
	      ck.gradFromDOF(b_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_z);
	      // calculate cell based CFL
	      calculateCFL(elementDiameter[eN],
			   g,
			   hn,
			   hun,
			   hvn,
			   hEps,
			   q_cfl[eN_k]);
	      //precalculate test function products with integration weights
	      for (int j=0;j<nDOF_trial_element;j++)
		h_test_dV[j] = h_test_ref[k*nDOF_trial_element+j]*dV;
	      //save velocity at quadrature points for other models to use
	      q_velocity[eN_k_nSpace+0] = 2*h/(h*h+std::pow(fmax(h,hEps),2))*hu;
	      q_velocity[eN_k_nSpace+1] = 2*h/(h*h+std::pow(fmax(h,hEps),2))*hv;
	      hnp1_at_quad_point[eN_k] = h;
	      hunp1_at_quad_point[eN_k] = hu;
	      hvnp1_at_quad_point[eN_k] = hv;
	      //update element residual. Part about the lumped mass matrix 
	      double one_over_hnReg = 2*hn/(hn*hn+std::pow(fmax(hn,hEps),2));
	      double dxf1 = grad_hun[0];
	      double dyf1 = grad_hvn[1];
	      double dxf2 = (2*hn*hun*grad_hun[0] - hun*hun*grad_hn[0])*std::pow(one_over_hnReg,2)+g*hn*grad_hn[0];
	      double dyf2 = (hn*(hun*grad_hvn[1]+hvn*grad_hun[1])-hun*hvn*grad_hn[1])*std::pow(one_over_hnReg,2);
	      double dxf3 = (hn*(hun*grad_hvn[0]+hvn*grad_hun[0])-hun*hvn*grad_hn[0])*std::pow(one_over_hnReg,2);
	      double dyf3 = (2*hn*hvn*grad_hvn[1] - hvn*hvn*grad_hn[1])*std::pow(one_over_hnReg,2)+g*hn*grad_hn[1];
	      // bathymetry source terms
	      double b_dot_gradx_z = g*hn*grad_z[0];
	      double b_dot_grady_z = g*hn*grad_z[1];
		
	      for(int i=0;i<nDOF_test_element;i++)
		{
		  elementResidual_h[i]  += ((h-hn)  + dt*(dxf1 + dyf1))*h_test_dV[i];
		  elementResidual_hu[i] += ((hu-hun) + dt*(dxf2 + dyf2) + dt*b_dot_gradx_z)*h_test_dV[i];
		  elementResidual_hv[i] += ((hv-hvn) + dt*(dxf3 + dyf3) + dt*b_dot_grady_z)*h_test_dV[i];
		}
	    }
	  // distribute
	  for(int i=0;i<nDOF_test_element;i++)
	    {
	      register int eN_i=eN*nDOF_test_element+i;
	      int h_gi = h_l2g[eN_i]; //global i-th index for h
	      int vel_gi = vel_l2g[eN_i]; //global i-th index for velocities 

	      // distribute element residual to global residual
	      globalResidual[offset_h+stride_h*h_gi]  += elementResidual_h[i];
	      globalResidual[offset_hu+stride_hu*vel_gi] += elementResidual_hu[i];
	      globalResidual[offset_hv+stride_hv*vel_gi] += elementResidual_hv[i];
	    }
	}
      
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
			   double* h_dof_old_old, 
			   double* hu_dof_old_old, 
			   double* hv_dof_old_old, 
			   double* h_dof_lstage, 
			   double* hu_dof_lstage, 
			   double* hv_dof_lstage, 
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
			   double* q_velocity_sge,
			   double* q_cfl,
			   double* q_numDiff_h, 
			   double* q_numDiff_hu, 
			   double* q_numDiff_hv, 
			   double* q_numDiff_h_last,
			   double* q_numDiff_hu_last, 
			   double* q_numDiff_hv_last,
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
			   double* edge_based_cfl,
			   double cfl_run,
			   double hEps,
			   int recompute_lumped_mass_matrix, 
			   // SAVE SOLUTION (mql)
			   double* hnp1_at_quad_point,
			   double* hunp1_at_quad_point,
			   double* hvnp1_at_quad_point,
			   // GALERKIN SOLUTION 
			   double* h_dof_galerkin,
			   double* hu_dof_galerkin,
			   double* hv_dof_galerkin, 
			   // TO COMPUTE LOW ORDER 
			   double* low_order_hnp1,
			   double* low_order_hunp1,
			   double* low_order_hvnp1,
			   double* dEV_minus_dL,
			   double cE,
			   int LUMPED_MASS_MATRIX,
			   int USE_EV_BASED_ON_GALERKIN,
			   double dt,
			   double mannings,
			   // Quant of interests
			   double* quantDOFs,
			   // Lumped mass matrix
			   double* ML, 
			   int SECOND_CALL_CALCULATE_RESIDUAL)
    {
      //FOR FRICTION//
      double n2 = std::pow(mannings,2.);
      double gamma=4./3;
      double xi=2.;

      ////////////////
      // CELL LOOPS //
      ////////////////
      // To compute: 
      //      * lumped_mass_matrix
      //      * Time derivative term
      //      * Cell based CFL
      //      * Velocity at quad points for other models       
      // init lumped mass matrix and ent residual vectors to zero
      for (int i=0; i<numDOFsPerEqn; i++)
	lumped_mass_matrix[i] = 0;
      for(int eN=0;eN<nElements_global;eN++)
	{
	  //declare local storage for element residual and initialize
	  register double 
	    element_lumped_mass_matrix[nDOF_test_element],
	    elementResidual_h[nDOF_test_element],
	    elementResidual_hu[nDOF_test_element],
	    elementResidual_hv[nDOF_test_element];
	    
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      element_lumped_mass_matrix[i]=0.0;
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
		h_lstage=0.0,hu_lstage=0.0,hv_lstage=0.0, // solution at lstage
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
	      ck.valFromDOF(h_dof,&h_l2g[eN_nDOF_trial_element],&h_trial_ref[k*nDOF_trial_element],h);
	      ck.valFromDOF(hu_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hu);
	      ck.valFromDOF(hv_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hv);
	      // get the solution at the lstage
	      ck.valFromDOF(h_dof_lstage,&h_l2g[eN_nDOF_trial_element],&h_trial_ref[k*nDOF_trial_element],h_lstage);
	      ck.valFromDOF(hu_dof_lstage,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hu_lstage);
	      ck.valFromDOF(hv_dof_lstage,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hv_lstage);
	      // calculate cell based CFL to keep a reference
	      calculateCFL(elementDiameter[eN],
			   g,
			   h_lstage,
			   hu_lstage,
			   hv_lstage,
			   hEps,
			   q_cfl[eN_k]);
	      //precalculate test function products with integration weights
	      for (int j=0;j<nDOF_trial_element;j++)
		h_test_dV[j] = h_test_ref[k*nDOF_trial_element+j]*dV;
	      //save velocity at quadrature points for other models to use
	      q_velocity[eN_k_nSpace+0] = 2*h/(h*h+std::pow(fmax(h,hEps),2))*hu;
	      q_velocity[eN_k_nSpace+1] = 2*h/(h*h+std::pow(fmax(h,hEps),2))*hv;
	      hnp1_at_quad_point[eN_k] = h;
	      hunp1_at_quad_point[eN_k] = hu;
	      hvnp1_at_quad_point[eN_k] = hv;

	      for(int i=0;i<nDOF_test_element;i++)
		{
		  // lumped mass matrix
		  element_lumped_mass_matrix[i] += h_test_dV[i];
		  // compute time derivative part of global residual. NOTE: no lumping
		  elementResidual_h[i]  += (h  - h_lstage)*h_test_dV[i];
		  elementResidual_hu[i] += (hu - hu_lstage)*h_test_dV[i];
		  elementResidual_hv[i] += (hv - hv_lstage)*h_test_dV[i];
		}
	    }
	  // distribute
	  for(int i=0;i<nDOF_test_element;i++)
	    {
	      register int eN_i=eN*nDOF_test_element+i;
	      int h_gi = h_l2g[eN_i]; //global i-th index for h
	      int vel_gi = vel_l2g[eN_i]; //global i-th index for velocities 

	      // distribute lumped mass matrix
	      lumped_mass_matrix[h_gi]  += element_lumped_mass_matrix[i];
	      // distribute time derivative to global residual
	      globalResidual[offset_h+stride_h*h_gi]  += elementResidual_h[i];
	      globalResidual[offset_hu+stride_hu*vel_gi] += elementResidual_hu[i];
	      globalResidual[offset_hv+stride_hv*vel_gi] += elementResidual_hv[i];
	    }
	}
      
      //////////////////////////////////
      // COMPUTE ENTROPY and g VECTOR //
      //////////////////////////////////
      // compute entropy (defined as eta) corresponding to ith node
      register double eta[numDOFsPerEqn], gx[numDOFsPerEqn], gy[numDOFsPerEqn], 
	alpha_numerator_pos[numDOFsPerEqn], alpha_numerator_neg[numDOFsPerEqn],
	alpha_denominator_pos[numDOFsPerEqn], alpha_denominator_neg[numDOFsPerEqn];
      int ij=0;
      for (int i=0; i<numDOFsPerEqn; i++)
	{
	  // COMPUTE ENTROPY BASED ON OLD STAGE
	  double hni = h_dof_lstage[i]; 
	  double one_over_hniReg = 2*hni/(hni*hni+std::pow(fmax(hni,hEps),2));
	  eta[i] = ENTROPY(g,hni,hu_dof_lstage[i],hv_dof_lstage[i],one_over_hniReg);
	  // for copmuting g vector 
	  gx[i]=0.;
	  gy[i]=0.;
	  // for smoothness indicator 
	  alpha_numerator_pos[i] = 0.;
	  alpha_numerator_neg[i] = 0.;
	  alpha_denominator_pos[i] = 0.;
	  alpha_denominator_neg[i] = 0.;

	  for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
	    {
	      int j = csrColumnOffsets_DofLoops[offset];
	      double hnj = h_dof_lstage[j];
	      // for gi vector
	      gx[i] += Cx[ij]*hnj;
	      gy[i] += Cy[ij]*hnj;

	      double alpha_num = hnj - hni;
	      alpha_numerator_pos[i] += alpha_num > 0 ? alpha_num : 0.;
	      alpha_numerator_neg[i] += alpha_num < 0 ? alpha_num : 0.;
	      alpha_denominator_pos[i] += alpha_num > 0 ? alpha_num : 0.;
	      alpha_denominator_neg[i] += alpha_num < 0 ? fabs(alpha_num) : 0.;

	      //update ij
	      ij+=1;
	    }
	  gx[i] /= ML[i];
	  gy[i] /= ML[i];
	}

      /////////////////////////////////////////////////////////////////////////////////////////////////
      // COMPUTE SMOOTHNESS INDICATOR, dLij, etaMin, etaMax, edge_based_cfl, global_entropy_residual //
      /////////////////////////////////////////////////////////////////////////////////////////////////
      // Smoothness indicator is based on the solution. psi_i = psi_i(alpha_i); 
      // alpha_i = |sum(uj-ui)|/sum|uj-ui|
      ij = 0;
      register double global_entropy_residual[numDOFsPerEqn]; 
      register double psi[numDOFsPerEqn], dL[NNZ], etaMax[numDOFsPerEqn], etaMin[numDOFsPerEqn];
      register double SumPos[numDOFsPerEqn], SumNeg[numDOFsPerEqn];
      for (int i=0; i<numDOFsPerEqn; i++)
	{
	  double alphai; // smoothness indicator of solution
	  // get location of i-nodes (for linearity preserving) smoothness indicator
	  double xi = mesh_dof[i*3+0];
	  double yi = mesh_dof[i*3+1];

	  double hi = h_dof_lstage[i]; // solution at time tn for the ith DOF
	  double hui = hu_dof_lstage[i]; 
	  double hvi = hv_dof_lstage[i]; 
	  double Zi = b_dof[i];
	  double dLii = 0.;

	  // For eta min and max
	  etaMax[i] = std::abs(eta[i]);
	  etaMin[i] = std::abs(eta[i]);

	  // FOR ENTROPY RESIDUAL
	  double ith_flux_term1=0., ith_flux_term2=0., ith_flux_term3=0.;
	  double entropy_flux=0.;
	  double one_over_hiReg = 2*hi/(hi*hi+std::pow(fmax(hi,hEps),2));
	  double eta_prime1 = DENTROPY_DH(g,hi,hui,hvi,one_over_hiReg); 
	  double eta_prime2 = DENTROPY_DHU(g,hi,hui,hvi,one_over_hiReg); 
	  double eta_prime3 = DENTROPY_DHV(g,hi,hui,hvi,one_over_hiReg); 

	  // init SumPos and SumNeg
	  SumPos[i] = 0.;
	  SumNeg[i] = 0.;
	  for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
	    { //loop in j (sparsity pattern)
	      int j = csrColumnOffsets_DofLoops[offset];
	      double hj = h_dof_lstage[j]; // solution at time tn for the jth DOF
	      double huj = hu_dof_lstage[j];
	      double hvj = hv_dof_lstage[j];
	      double Zj = b_dof[j];

	      // get location of j-nodes (for linearity preserving) smoothness indicator
	      double xj = mesh_dof[j*3+0];
	      double yj = mesh_dof[j*3+1];

	      // FOR ENTROPY RESIDUAL //
	      double one_over_hjReg = 2*hj/(hj*hj+std::pow(fmax(hj,hEps),2));
	      ith_flux_term1 += huj*Cx[ij] + hvj*Cy[ij]; 
	      ith_flux_term2 += ( huj*huj*one_over_hjReg*Cx[ij] + huj*hvj*one_over_hjReg*Cy[ij] 
				  + g*hi*(hj+Zj)*Cx[ij] );
	      ith_flux_term3 += ( huj*hvj*one_over_hjReg*Cx[ij] + hvj*hvj*one_over_hjReg*Cy[ij] 
				  + g*hi*(hj+Zj)*Cy[ij] );	      
	      
	      entropy_flux += ( Cx[ij]*ENTROPY_FLUX1(g,hj,huj,hvj,one_over_hjReg) + 
				Cy[ij]*ENTROPY_FLUX2(g,hj,huj,hvj,one_over_hjReg) );

	      // COMPUTE STAR SOLUTION // hStar, huStar and hvStar
	      double hStarij  = fmax(0., hi + Zi - fmax(Zi,Zj));
	      double huStarij = (hi <= hEps ? 0. : hui*hStarij/hi);
	      double hvStarij = (hi <= hEps ? 0. : hvi*hStarij/hi);
	      
	      double hStarji  = fmax(0., hj + Zj - fmax(Zi,Zj));
	      double huStarji = (hj <= hEps ? 0. : huj*hStarji/hj);
	      double hvStarji = (hj <= hEps ? 0. : hvj*hStarji/hj);
	   
	      //////////////////////////////
	      // FOR SMOOTHNESS INDICATOR //
	      //////////////////////////////
	      double gj_times_x = gx[j]*(xj-xi) + gy[j]*(yj-yi);
	      SumPos[i] += gj_times_x > 0 ? gj_times_x : 0;
	      SumNeg[i] += gj_times_x < 0 ? gj_times_x : 0;

	      //////////////////
	      // COMPUTE dLij // 
	      //////////////////
	      // NOTE: I don't have anything in the entry dLij when i=j
	      if (i != j)
		{		 
		  // compute dLij
		  double cij_norm = sqrt(Cx[ij]*Cx[ij] + Cy[ij]*Cy[ij]);
		  double cji_norm = sqrt(CTx[ij]*CTx[ij] + CTy[ij]*CTy[ij]);
		  
		  double nxij = Cx[ij]/cij_norm, nyij = Cy[ij]/cij_norm;
		  double nxji = CTx[ij]/cji_norm, nyji = CTy[ij]/cji_norm;
		  
		  //dissipative terms for flux
		  dL[ij] = fmax(fmax(maxWaveSpeedSharpInitialGuess(g,nxij,nyij, 
								   hi,hui,hvi,
								   hStarji,huStarji,hvStarji,
								   hEps,false),
				     maxWaveSpeedSharpInitialGuess(g,nxij,nyij, 
								   hi,hui,hvi, 
								   hStarij,huStarij,hvStarij,
								   hEps,false))*cij_norm,
				fmax(maxWaveSpeedSharpInitialGuess(g,nxji,nyji, 
								   hj,huj,hvj, 
								   hStarij,huStarij,hvStarij,
								   hEps,false),
				     maxWaveSpeedSharpInitialGuess(g,nxji,nyji, 
								   hj,huj,hvj, 
								   hStarji,huStarji,hvStarji,
								   hEps,false))*cji_norm);  
		  dLii -= dL[ij];
		  
		  /////////////////////////////////
		  // COMPUTE ETA MIN AND ETA MAX // 
		  /////////////////////////////////
		  etaMax[i] = fmax(etaMax[i],std::abs(eta[j]));
    		  etaMin[i] = fmin(etaMin[i],std::abs(eta[j]));
		}
	      //update ij
	      ij+=1;
	    }
	  double mi = lumped_mass_matrix[i];

	  //////////////////////////////
	  // CALCULATE EDGE BASED CFL //
	  //////////////////////////////
	  edge_based_cfl[i] = 2*fabs(dLii)/mi;

	  /////////////////////////////////////
	  // COMPUTE GLOBAL ENTROPY RESIDUAL //
	  /////////////////////////////////////
	  global_entropy_residual[i] = entropy_flux 
	    -(ith_flux_term1*eta_prime1 + ith_flux_term2*eta_prime2 + ith_flux_term3*eta_prime3);

	  //////////////////////////////////
	  // COMPUTE SMOOTHNESS INDICATOR //
	  //////////////////////////////////
	  //if (hi < hEps)
	  //alphai = 1.;
	  //else
	    {
	      // Compute sigmaPos and sigmaNeg
	      double Sum = SumPos[i] + SumNeg[i];
	      double sigmaPosi = Sum < 0 ? 1. : (-SumNeg[i])/(SumPos[i]+1E-15);
	      double sigmaNegi = Sum > 0 ? 1. : (SumPos[i])/(-SumNeg[i]+1E-15);

	      double alpha_numi = fabs(sigmaPosi*alpha_numerator_pos[i] + sigmaNegi * alpha_numerator_neg[i]);
	      double alpha_deni = sigmaPosi*alpha_denominator_pos[i] + sigmaNegi * alpha_denominator_neg[i];

	      if (BETAij == 1)
		{
		  alpha_numi = fabs(alpha_numerator_pos[i] + alpha_numerator_neg[i]);
		  alpha_deni = alpha_denominator_pos[i] + alpha_denominator_neg[i];
		}
	      if (fabs(alpha_numerator_pos[i] + alpha_numerator_neg[i]) == 0 ) //constante state
		alphai = 0.;
	      else 
		alphai = alpha_numi/(alpha_deni+1E-15);
	    }

	  quantDOFs[i] = alphai;
	  if (POWER_SMOOTHNESS_INDICATOR==0)
	    psi[i] = 1.0;
	  else
	    psi[i] = std::pow(alphai,POWER_SMOOTHNESS_INDICATOR); //NOTE: they use alpha^2 in the paper
	}

      //////////////////
      // Loop on DOFs // to compute flux and dissipative terms
      //////////////////
      ij = 0;
      for (int i=0; i<numDOFsPerEqn; i++)
	{
	  double hi = h_dof_lstage[i];
	  double hui = hu_dof_lstage[i];
	  double hvi = hv_dof_lstage[i];	  
	  double Zi = b_dof[i];
	  
	  double ui = 2*hi/(hi*hi+std::pow(fmax(hi,hEps),2))*hui;
	  double vi = 2*hi/(hi*hi+std::pow(fmax(hi,hEps),2))*hvi;

	  // entropy normalization factor 
	  double one_over_entNormFactori = etaMax[i] == etaMin[i] ? 0. : 1./(etaMax[i]-etaMin[i]);
	  // regularization of 1/hi 
	  double one_over_hiReg = 2*hi/(hi*hi+std::pow(fmax(hi,hEps),2));

	  double ith_flux_term1=0., ith_flux_term2=0., ith_flux_term3=0.;
	  double 
	    ith_dLij_minus_muij_times_hStarStates=0.,
	    ith_dLij_minus_muij_times_huStarStates=0.,
	    ith_dLij_minus_muij_times_hvStarStates=0.,
	    ith_dHij_minus_muij_times_hStarStates=0.,
	    ith_dHij_minus_muij_times_huStarStates=0.,
	    ith_dHij_minus_muij_times_hvStarStates=0.,
	    ith_muij_times_hStates=0.,
	    ith_muij_times_huStates=0.,
	    ith_muij_times_hvStates=0.;
	  double aux1_to_compute_hnp1=0., aux2_to_compute_hnp1=0.;	    

	  ///////////////////
	  // FRICTION TERM //
	  ///////////////////
	  double veli_norm = std::sqrt(ui*ui+vi*vi);
	  double mi = lumped_mass_matrix[i];	  
	  double hi_to_the_gamma = std::pow(hi,gamma);
	  double ith_friction_term2 =  
	    veli_norm==0 ? 0. : 2*g*n2*hui*veli_norm*mi/(hi_to_the_gamma+fmax(hi_to_the_gamma,xi*g*n2*dt*veli_norm));
	  double ith_friction_term3 =  
	    veli_norm==0 ? 0. : 2*g*n2*hvi*veli_norm*mi/(hi_to_the_gamma+fmax(hi_to_the_gamma,xi*g*n2*dt*veli_norm));
	  
	  if (LINEAR_FRICTION==1)
	    {
	      ith_friction_term2 = mannings*hui*mi;
	      ith_friction_term3 = mannings*hvi*mi;
	    }
	  // loop over the sparsity pattern of the i-th DOF
	  for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
	    {
	      int j = csrColumnOffsets_DofLoops[offset];
	      double hj = h_dof_lstage[j];
	      double huj = hu_dof_lstage[j];
	      double hvj = hv_dof_lstage[j];
	      double Zj = b_dof[j];
	      double uj = 2*hj/(hj*hj+std::pow(fmax(hj,hEps),2))*huj;
	      double vj = 2*hj/(hj*hj+std::pow(fmax(hj,hEps),2))*hvj;

	      // regularization of 1/hj
	      double one_over_hjReg = 2*hj/(hj*hj+std::pow(fmax(hj,hEps),2));

	      // Nodal projection of fluxes
	      ith_flux_term1 += huj*Cx[ij] + hvj*Cy[ij]; // f1*C
	      ith_flux_term2 += ( huj*huj*one_over_hjReg*Cx[ij] + huj*hvj*one_over_hjReg*Cy[ij] 
				  + g*hi*(hj+Zj)*Cx[ij] );
	      ith_flux_term3 += ( huj*hvj*one_over_hjReg*Cx[ij] + hvj*hvj*one_over_hjReg*Cy[ij] 
				  + g*hi*(hj+Zj)*Cy[ij] );

	      // Dissipative well balancing term
	      double muij = 0.;
	      if (i != j) // This is not necessary. See formula for ith_dissipative_terms
		{
		  ///////////////////////////
		  // COMPUTE STAR SOLUTION // hStar, huStar and hvStar
		  ///////////////////////////
		  double hStarij  = fmax(0., hi + Zi - fmax(Zi,Zj));
		  double huStarij = (hi <= hEps ? 0. : hui*hStarij/hi);
		  double hvStarij = (hi <= hEps ? 0. : hvi*hStarij/hi);
		  
		  double hStarji  = fmax(0., hj + Zj - fmax(Zi,Zj));
		  double huStarji = (hj <= hEps ? 0. : huj*hStarji/hj);
		  double hvStarji = (hj <= hEps ? 0. : hvj*hStarji/hj);
		  
		  ///////////////////////////////////////
		  // WELL BALANCING DISSIPATIVE MATRIX //
		  ///////////////////////////////////////
		  muij = fmax(fmax(0.,-(ui*Cx[ij] + vi*Cy[ij])),fmax(0,(uj*Cx[ij] + vj*Cy[ij])));
		  muij *= std::max(psi[i],psi[j]); // enhance its order
		  
		  ////////////////////////
		  // DISSIPATIVE MATRIX //
		  ////////////////////////
		  // compute dissipative terms for second and third equations
		  dL[ij] *= std::max(psi[i],psi[j]); // enhance its order
		  // compute entropy residual for each component 
		  double one_over_entNormFactorj = etaMax[j] == etaMin[j] ? 0. : 1./(etaMax[j]-etaMin[j]);
		  double dEVij = cE*fmax(std::abs(global_entropy_residual[i])*one_over_entNormFactori,
					 std::abs(global_entropy_residual[j])*one_over_entNormFactorj);
		  double dHij = fmin(dL[ij],dEVij);
		  if (cE >= 1000)
		    dHij = dL[ij];

		  // compute dij_minus_muij times star solution terms
		  ith_dHij_minus_muij_times_hStarStates  += (dHij - muij)*(hStarji-hStarij); 
		  ith_dHij_minus_muij_times_huStarStates += (dHij - muij)*(huStarji-huStarij); 
		  ith_dHij_minus_muij_times_hvStarStates += (dHij - muij)*(hvStarji-hvStarij); 

		  ith_dLij_minus_muij_times_hStarStates  += (dL[ij] - muij)*(hStarji-hStarij); 
		  ith_dLij_minus_muij_times_huStarStates += (dL[ij] - muij)*(huStarji-huStarij); 
		  ith_dLij_minus_muij_times_hvStarStates += (dL[ij] - muij)*(hvStarji-hvStarij); 

		  // compute muij times solution terms
		  ith_muij_times_hStates  += muij*(hj-hi);
		  ith_muij_times_huStates += muij*(huj-hui);
		  ith_muij_times_hvStates += muij*(hvj-hvi);

		  // compute dEV_minus_dL
		  dEV_minus_dL[ij] = dHij - dL[ij];

		  // compute aux quantities for first equation 
		  // aux1 = Veli*Cii + sum_j[ muij + (dLij-muij)*hStarij/hi ]
		  // aux2 = sum_j[ (muij-Velj*Cij)*hj + (dLij-muij)*hStarji ] >= 0
		  aux1_to_compute_hnp1 += muij + (hi == 0. ? 0. : (dL[ij] - muij)*hStarij/hi);
		  aux2_to_compute_hnp1 += (muij - (uj*Cx[ij] + vj*Cy[ij]))*hj + (dL[ij] - muij)*hStarji;
		}
	      else // i==j
		{
		  aux1_to_compute_hnp1 += (ui*Cx[ij] + vi*Cy[ij])*hi;
		  dEV_minus_dL[ij]=0.; //Not true but the prod of this times Uj-Ui will be zero
		}
	      // update ij
	      ij+=1;
	    }
	  
	  // Compute low order solution: lumped mass matrix and low order dissipative matrix
	  //low_order_hnp1[i]  = hi*(1-dt/mi*aux1_to_compute_hnp1) + dt/mi*aux2_to_compute_hnp1;
	  low_order_hnp1[i]  = hi  - dt/mi*(ith_flux_term1
					    - ith_dLij_minus_muij_times_hStarStates
					    - ith_muij_times_hStates);
	  low_order_hunp1[i] = hui - dt/mi*(ith_flux_term2 
					    - ith_dLij_minus_muij_times_huStarStates
					    - ith_muij_times_huStates
					    + ith_friction_term2);
	  low_order_hvnp1[i] = hvi - dt/mi*(ith_flux_term3 
					    - ith_dLij_minus_muij_times_hvStarStates
					    - ith_muij_times_hvStates
					    + ith_friction_term3); 
	  // FIX LOW ORDER SOLUTION //
	  low_order_hnp1[i] = fmax(low_order_hnp1[i],0.);
	  double aux = fmax(low_order_hnp1[i],hEps);
	  low_order_hunp1[i] *= 2*std::pow(low_order_hnp1[i],2.)/(std::pow(low_order_hnp1[i],2.)+std::pow(aux,2.));
	  low_order_hvnp1[i] *= 2*std::pow(low_order_hnp1[i],2.)/(std::pow(low_order_hnp1[i],2.)+std::pow(aux,2.));

	  int LOW_ORDER_SOLUTION=0;
	  if (LOW_ORDER_SOLUTION==1)
	    {
	      globalResidual[offset_h+stride_h*i]   = low_order_hnp1[i];
	      globalResidual[offset_hu+stride_hu*i] = low_order_hunp1[i];
	      globalResidual[offset_hv+stride_hv*i] = low_order_hvnp1[i];
	    }	    
	  else
	    {
	      if (LUMPED_MASS_MATRIX==1)
		{
		  globalResidual[offset_h+stride_h*i]   
		    = hi - dt/mi*(ith_flux_term1 
				  - ith_dHij_minus_muij_times_hStarStates
				  - ith_muij_times_hStates);
		  globalResidual[offset_hu+stride_hu*i] 
		    = hui - dt/mi*(ith_flux_term2 
				   - ith_dHij_minus_muij_times_huStarStates
				   - ith_muij_times_huStates
				   + ith_friction_term2);
		  globalResidual[offset_hv+stride_hv*i] 
		    = hvi - dt/mi*(ith_flux_term3 
				   - ith_dHij_minus_muij_times_hvStarStates
				   - ith_muij_times_hvStates
				   + ith_friction_term3); 
		}
	      else
		{
		  // Distribute residual
		  // NOTE: MASS MATRIX IS CONSISTENT
		  globalResidual[offset_h+stride_h*i]   
		    += dt*(ith_flux_term1 
			   - ith_dHij_minus_muij_times_hStarStates
			   - ith_muij_times_hStates);   
		  globalResidual[offset_hu+stride_hu*i] 
		    += dt*(ith_flux_term2 
			   - ith_dHij_minus_muij_times_huStarStates
			   - ith_muij_times_huStates
			   + ith_friction_term2);
		  globalResidual[offset_hv+stride_hv*i] 
		    += dt*(ith_flux_term3 
			   - ith_dHij_minus_muij_times_hvStarStates
			   - ith_muij_times_hvStates
			   + ith_friction_term3); 
		}
	    }
	}
    }

    
    void calculateJacobian(//element
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
			   double* q_velocity_sge,
			   double* q_cfl,
			   double* q_numDiff_h_last,
			   double* q_numDiff_hu_last, 
			   double* q_numDiff_hv_last, 
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
	  register double  elementJacobian_h_h[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_h_hu[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_h_hv[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_hu_h[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_hu_hu[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_hu_hv[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_hv_h[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_hv_hu[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_hv_hv[nDOF_test_element][nDOF_trial_element];
	  for (int i=0;i<nDOF_test_element;i++)
	    for (int j=0;j<nDOF_trial_element;j++)
	      {
		elementJacobian_h_h[i][j]=0.0;
		elementJacobian_h_hu[i][j]=0.0;
		elementJacobian_h_hv[i][j]=0.0;
		elementJacobian_hu_h[i][j]=0.0;
		elementJacobian_hu_hu[i][j]=0.0;
		elementJacobian_hu_hv[i][j]=0.0;
		elementJacobian_hv_h[i][j]=0.0;
		elementJacobian_hv_hu[i][j]=0.0;
		elementJacobian_hv_hv[i][j]=0.0;
	      }
	  for  (int k=0;k<nQuadraturePoints_element;k++)
	    {
	      int eN_k = eN*nQuadraturePoints_element+k, //index to a scalar at a quadrature point
		eN_k_nSpace = eN_k*nSpace,
		eN_nDOF_trial_element = eN*nDOF_trial_element; //index to a vector at a quadrature point

	      //declare local storage
	      register double b=0.0,
		h=0.0,
		hu=0.0,
		hv=0.0,
		h_sge=0.0,
		hu_sge=0.0,
		hv_sge=0.0,
		grad_b[nSpace],
		grad_h[nSpace],
		grad_hu[nSpace],
		grad_hv[nSpace],
		mass_acc=0.0,
		dmass_acc_h=0.0,
		mom_hu_acc=0.0,
		dmom_hu_acc_h=0.0,
		dmom_hu_acc_hu=0.0,
		mom_hv_acc=0.0,
		dmom_hv_acc_h=0.0,
		dmom_hv_acc_hv=0.0,
		mass_adv[nSpace],
		dmass_adv_h[nSpace],
		dmass_adv_hu[nSpace],
		dmass_adv_hv[nSpace],
		dmass_adv_h_sge[nSpace],
		dmass_adv_hu_sge[nSpace],
		dmass_adv_hv_sge[nSpace],
		mom_hu_adv[nSpace],
		dmom_hu_adv_h[nSpace],
		dmom_hu_adv_hu[nSpace],
		dmom_hu_adv_hv[nSpace],
		dmom_hu_adv_h_sge[nSpace],
		dmom_hu_adv_hu_sge[nSpace],
		dmom_hu_adv_hv_sge[nSpace],
		mom_hv_adv[nSpace],
		dmom_hv_adv_h[nSpace],
		dmom_hv_adv_hu[nSpace],
		dmom_hv_adv_hv[nSpace],
		dmom_hv_adv_h_sge[nSpace],
		dmom_hv_adv_hu_sge[nSpace],
		dmom_hv_adv_hv_sge[nSpace],
		mom_hu_diff_ten[nSpace],
		mom_hv_diff_ten[nSpace],
		mom_huhv_diff_ten[1],
		mom_hvhu_diff_ten[1],
		mom_hu_source=0.0,
		dmom_hu_source_h=0.0,
		mom_hv_source=0.0,
		dmom_hv_source_h=0.0,
		mass_acc_t=0.0,
		dmass_acc_h_t=0.0,
		mom_hu_acc_t=0.0,
		dmom_hu_acc_h_t=0.0,
		dmom_hu_acc_hu_t=0.0,
		mom_hv_acc_t=0.0,
		dmom_hv_acc_h_t=0.0,
		dmom_hv_acc_hv_t=0.0,
		tau[9],
		pdeResidual_h=0.0,
		pdeResidual_hu=0.0,
		pdeResidual_hv=0.0,
		dpdeResidual_h_h[nDOF_trial_element],
		dpdeResidual_h_hu[nDOF_trial_element],
		dpdeResidual_h_hv[nDOF_trial_element],
		dpdeResidual_hu_h[nDOF_trial_element],
		dpdeResidual_hu_hu[nDOF_trial_element],
		dpdeResidual_hu_hv[nDOF_trial_element],
		dpdeResidual_hv_h[nDOF_trial_element],
		dpdeResidual_hv_hu[nDOF_trial_element],
		dpdeResidual_hv_hv[nDOF_trial_element],
		Lstar_h_h[nDOF_test_element],
		Lstar_hu_h[nDOF_test_element],
		Lstar_hv_h[nDOF_test_element],
		Lstar_h_hu[nDOF_test_element],
		Lstar_hu_hu[nDOF_test_element],
		Lstar_hv_hu[nDOF_test_element],
		Lstar_h_hv[nDOF_test_element],
		Lstar_hu_hv[nDOF_test_element],
		Lstar_hv_hv[nDOF_test_element],
		subgridError_h=0.0,
		subgridError_hu=0.0,
		subgridError_hv=0.0,
		dsubgridError_h_h[nDOF_trial_element],
		dsubgridError_h_hu[nDOF_trial_element],
		dsubgridError_h_hv[nDOF_trial_element],
		dsubgridError_hu_h[nDOF_trial_element],
		dsubgridError_hu_hu[nDOF_trial_element],
		dsubgridError_hu_hv[nDOF_trial_element],
		dsubgridError_hv_h[nDOF_trial_element],
		dsubgridError_hv_hu[nDOF_trial_element],
		dsubgridError_hv_hv[nDOF_trial_element],
		tau_h=0.0,tau_h0=0.0,tau_h1=0.0,
		tau_hv=0.0,tau_hv0=0.0,tau_hv1=0.0,
		jac[nSpace*nSpace],
		jacDet,
		jacInv[nSpace*nSpace],
		h_grad_trial[nDOF_trial_element*nSpace],
		vel_grad_trial[nDOF_trial_element*nSpace],
		dV,
		h_test_dV[nDOF_test_element],
		vel_test_dV[nDOF_test_element],
		h_grad_test_dV[nDOF_test_element*nSpace],
		vel_grad_test_dV[nDOF_test_element*nSpace],
		x,y,xt,yt,
		G[nSpace*nSpace],G_dd_G,tr_G, dmom_adv_star[nSpace], dmom_adv_sge[nSpace];
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
	      ck.calculateMappingVelocity_element(eN,
						  k,
						  mesh_velocity_dof,
						  mesh_l2g,
						  mesh_trial_ref,
						  xt,yt);
	      //get the physical integration weight
	      dV = fabs(jacDet)*dV_ref[k];

	      //get the trial function gradients
	      ck.gradTrialFromRef(&h_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,h_grad_trial);
	      ck.gradTrialFromRef(&vel_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,vel_grad_trial);
	      //get the solution 	
	      ck.valFromDOF(b_dof,&h_l2g[eN_nDOF_trial_element],&h_trial_ref[k*nDOF_trial_element],b);
	      ck.valFromDOF(h_dof,&h_l2g[eN_nDOF_trial_element],&h_trial_ref[k*nDOF_trial_element],h);
	      ck.valFromDOF(hu_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hu);
	      ck.valFromDOF(hv_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hv);
	      ck.valFromDOF(h_dof_sge,&h_l2g[eN_nDOF_trial_element],&h_trial_ref[k*nDOF_trial_element],h_sge);
	      ck.valFromDOF(hu_dof_sge,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hu_sge);
	      ck.valFromDOF(hv_dof_sge,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hv_sge);
	      //get the solution gradients
	      ck.gradFromDOF(b_dof,&h_l2g[eN_nDOF_trial_element],h_grad_trial,grad_b);
	      ck.gradFromDOF(h_dof,&h_l2g[eN_nDOF_trial_element],h_grad_trial,grad_h);
	      ck.gradFromDOF(hu_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_hu);
	      ck.gradFromDOF(hv_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_hv);
	      //precalculate test function products with integration weights
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  h_test_dV[j] = h_test_ref[k*nDOF_trial_element+j]*dV;
		  vel_test_dV[j] = vel_test_ref[k*nDOF_trial_element+j]*dV;
		  for (int I=0;I<nSpace;I++)
		    {
		      h_grad_test_dV[j*nSpace+I]   = h_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
		      vel_grad_test_dV[j*nSpace+I] = vel_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin}
		    }
		}
	      evaluateCoefficients(nu,
				   g,
				   grad_b,
				   h,
				   hu,
				   hv,
				   mass_acc,
				   dmass_acc_h,
				   mom_hu_acc,
				   dmom_hu_acc_h,
				   dmom_hu_acc_hu,
				   mom_hv_acc,
				   dmom_hv_acc_h,
				   dmom_hv_acc_hv,
				   mass_adv,
				   dmass_adv_h,
				   dmass_adv_hu,
				   dmass_adv_hv,
				   mom_hu_adv,
				   dmom_hu_adv_h,
				   dmom_hu_adv_hu,
				   dmom_hu_adv_hv,
				   mom_hv_adv,
				   dmom_hv_adv_h,
				   dmom_hv_adv_hu,
				   dmom_hv_adv_hv,
				   mom_hu_diff_ten,
				   mom_hv_diff_ten,
				   mom_huhv_diff_ten,
				   mom_hvhu_diff_ten,
				   mom_hu_source,
				   dmom_hu_source_h,
				   mom_hv_source,
				   dmom_hv_source_h);
	      //
	      //moving mesh
	      //
	      /* mass_adv[0] -= MOVING_DOMAIN*mass_acc*xt; */
	      /* mass_adv[1] -= MOVING_DOMAIN*mass_acc*yt; */

	      /* dmass_adv_h[0] -= MOVING_DOMAIN*dmass_acc_h*xt; */
	      /* dmass_adv_h[1] -= MOVING_DOMAIN*dmass_acc_h*yt; */

	      /* mom_hu_adv[0] -= MOVING_DOMAIN*mom_hu_acc*xt; */
	      /* mom_hu_adv[1] -= MOVING_DOMAIN*mom_hu_acc*yt; */

	      /* dmom_hu_adv_h[0] -= MOVING_DOMAIN*dmom_hu_acc_h*xt; */
	      /* dmom_hu_adv_h[1] -= MOVING_DOMAIN*dmom_hu_acc_h*yt; */

	      /* dmom_hu_adv_hu[0] -= MOVING_DOMAIN*dmom_hu_acc_hu*xt; */
	      /* dmom_hu_adv_hu[1] -= MOVING_DOMAIN*dmom_hu_acc_hu*yt; */

	      /* mom_hv_adv[0] -= MOVING_DOMAIN*mom_hv_acc*xt; */
	      /* mom_hv_adv[1] -= MOVING_DOMAIN*mom_hv_acc*yt; */

	      /* dmom_hv_adv_hv[0] -= MOVING_DOMAIN*dmom_hv_acc_h*xt; */
	      /* dmom_hv_adv_hv[1] -= MOVING_DOMAIN*dmom_hv_acc_h*yt; */

	      /* dmom_hv_adv_hv[0] -= MOVING_DOMAIN*dmom_hv_acc_hv*xt; */
	      /* dmom_hv_adv_hv[1] -= MOVING_DOMAIN*dmom_hv_acc_hv*yt; */
	      //
	      //calculate time derivatives
	      //
	      ck.bdf(alphaBDF,
		     q_mass_acc_beta_bdf[eN_k],
		     mass_acc,
		     dmass_acc_h,
		     mass_acc_t,
		     dmass_acc_h_t);
	      ck.bdf(alphaBDF,
		     q_mom_hu_acc_beta_bdf[eN_k],
		     mom_hu_acc,
		     dmom_hu_acc_hu,
		     mom_hu_acc_t,
		     dmom_hu_acc_hu_t);
	      ck.bdf(alphaBDF,
		     q_mom_hv_acc_beta_bdf[eN_k],
		     mom_hv_acc,
		     dmom_hv_acc_hv,
		     mom_hv_acc_t,
		     dmom_hv_acc_hv_t);
	      //
	      //calculate subgrid error contribution to the Jacobian (strong residual, adjoint, jacobian of strong residual)
	      //
	      //
	      //calculate strong residual
	      //
	      double hStar_sge = fmax(1.0e-8,h_sge);

	      /* dmass_adv_h_sge[0]  = dmass_adv_h[0]; */
	      /* dmass_adv_h_sge[1]  = dmass_adv_h[1]; */
	      /* dmass_adv_hu_sge[0]  = dmass_adv_hu[0]; */
	      /* dmass_adv_hu_sge[1]  = dmass_adv_hu[1]; */
	      /* dmass_adv_hv_sge[0]  = dmass_adv_hv[0]; */
	      /* dmass_adv_hv_sge[1]  = dmass_adv_hv[1]; */
	      /* dmom_hu_adv_h_sge[0] = dmom_hu_adv_h[0]; */
	      /* dmom_hu_adv_h_sge[1] = dmom_hu_adv_h[1]; */
	      /* dmom_hu_adv_hu_sge[0] = dmom_hu_adv_hu[0]; */
	      /* dmom_hu_adv_hu_sge[1] = dmom_hu_adv_hu[1]; */
	      /* dmom_hu_adv_hv_sge[0] = dmom_hu_adv_hv[0]; */
	      /* dmom_hu_adv_hv_sge[1] = dmom_hu_adv_hv[1]; */
	      /* dmom_hv_adv_h_sge[0] = dmom_hv_adv_h[0]; */
	      /* dmom_hv_adv_h_sge[1] = dmom_hv_adv_h[1]; */
	      /* dmom_hv_adv_hu_sge[0] = dmom_hv_adv_hu[0]; */
	      /* dmom_hv_adv_hu_sge[1] = dmom_hv_adv_hu[1]; */
	      /* dmom_hv_adv_hv_sge[0] = dmom_hv_adv_hv[0]; */
	      /* dmom_hv_adv_hv_sge[1] = dmom_hv_adv_hv[1]; */

	      //lagged strong residual coefficients
	      //mass advective flux
	      dmass_adv_h_sge[0]=0.0;
	      dmass_adv_h_sge[1]=0.0;
	      
	      dmass_adv_hu_sge[0]=1.0;
	      dmass_adv_hu_sge[1]=0.0;
	      
	      dmass_adv_hv_sge[0]=0.0;
	      dmass_adv_hv_sge[1]=1.0;
	      
	      //u momentum advective flux
	      dmom_hu_adv_h_sge[0]=-hu_sge*hu_sge/(hStar_sge*hStar_sge) + g*h_sge;
	      dmom_hu_adv_h_sge[1]=-hu_sge*hv_sge/(hStar_sge*hStar_sge);
	      
	      dmom_hu_adv_hu_sge[0]=2.0*hu_sge/hStar_sge;
	      dmom_hu_adv_hu_sge[1]=hv_sge/hStar_sge;
	      
	      dmom_hu_adv_hv_sge[0]=0.0;
	      dmom_hu_adv_hv_sge[1]=hu_sge/hStar_sge;
	      
	      //v momentum advective_flux
	      dmom_hv_adv_h_sge[0]=-hv_sge*hu_sge/(hStar_sge*hStar_sge);
	      dmom_hv_adv_h_sge[1]=-hv_sge*hv_sge/(hStar_sge*hStar_sge) + g*h_sge;
	      
	      dmom_hv_adv_hu_sge[0]=hv_sge/hStar_sge;
	      dmom_hv_adv_hu_sge[1]=0.0;
	      
	      dmom_hv_adv_hv_sge[0]=hu_sge/hStar_sge;
	      dmom_hv_adv_hv_sge[1]=2.0*hv_sge/hStar_sge;

	      //approximate linearization

	      /* //mass advective flux */
	      /* dmass_adv_h_sge[0]=0.0; */
	      /* dmass_adv_h_sge[1]=0.0; */
	      
	      /* dmass_adv_hu_sge[0]=1.0; */
	      /* dmass_adv_hu_sge[1]=0.0; */
	      
	      /* dmass_adv_hv_sge[0]=0.0; */
	      /* dmass_adv_hv_sge[1]=1.0; */
	      
	      /* //u momentum advective flux */
	      /* dmom_hu_adv_h_sge[0]= g*h_sge; */
	      /* dmom_hu_adv_h_sge[1]= 0.0; */
	      
	      /* dmom_hu_adv_hu_sge[0]= hu_sge/hStar_sge; */
	      /* dmom_hu_adv_hu_sge[1]= hv_sge/hStar_sge; */
	      
	      /* dmom_hu_adv_hv_sge[0]= 0.0; */
	      /* dmom_hu_adv_hv_sge[1]= 0.0; */
  
	      /* //v momentum advective_flux */
	      /* dmom_hv_adv_h_sge[0]= 0.0; */
	      /* dmom_hv_adv_h_sge[1]= g*h_sge; */
	      
	      /* dmom_hv_adv_hu_sge[0]= 0.0; */
	      /* dmom_hv_adv_hu_sge[1]= 0.0; */
	      
	      /* dmom_hv_adv_hv_sge[0]= hu_sge/hStar_sge; */
	      /* dmom_hv_adv_hv_sge[1]= hv_sge/hStar_sge; */

	      //calculate the Jacobian of strong residual
	      for (int j=0;j<nDOF_trial_element;j++)
	      	{
	      	  register int j_nSpace = j*nSpace;

	      	  dpdeResidual_h_h[j]= ck.MassJacobian_strong(dmass_acc_h_t,h_trial_ref[k*nDOF_trial_element+j]);

	      	  dpdeResidual_h_hu[j]=ck.AdvectionJacobian_strong(dmass_adv_hu_sge,&vel_grad_trial[j_nSpace]);

	      	  dpdeResidual_h_hv[j]=ck.AdvectionJacobian_strong(dmass_adv_hv_sge,&vel_grad_trial[j_nSpace]);

	      	  dpdeResidual_hu_h[j]= ck.AdvectionJacobian_strong(dmom_hu_adv_h_sge,&h_grad_trial[j_nSpace]) +
		    ck.ReactionJacobian_strong(dmom_hu_source_h,h_trial_ref[k*nDOF_trial_element+j]);

	      	  dpdeResidual_hu_hu[j]= ck.MassJacobian_strong(dmom_hu_acc_hu_t,vel_trial_ref[k*nDOF_trial_element+j])+
		    ck.AdvectionJacobian_strong(dmom_hu_adv_hu_sge,&vel_grad_trial[j_nSpace]);

	      	  dpdeResidual_hu_hv[j]=ck.AdvectionJacobian_strong(dmom_hu_adv_hv_sge,&vel_grad_trial[j_nSpace]);

	      	  dpdeResidual_hv_h[j]= ck.AdvectionJacobian_strong(dmom_hv_adv_h_sge,&h_grad_trial[j_nSpace])+
		    ck.ReactionJacobian_strong(dmom_hv_source_h,h_trial_ref[k*nDOF_trial_element+j]);

	      	  dpdeResidual_hv_hu[j]=ck.AdvectionJacobian_strong(dmom_hv_adv_hu_sge,&vel_grad_trial[j_nSpace]);

	      	  dpdeResidual_hv_hv[j]= ck.MassJacobian_strong(dmom_hv_acc_hv_t,vel_trial_ref[k*nDOF_trial_element+j])+
		    ck.AdvectionJacobian_strong(dmom_hv_adv_hv_sge,&vel_grad_trial[j_nSpace]);

	      	}
	      //calculate tau and tau*Res
	      calculateSubgridError_tau(elementDiameter[eN],
	      				nu,
					g,
					h_sge,
					hu_sge,
					hv_sge,
					tau,
	      				q_cfl[eN_k]);
					
	      for (int j=0;j<nDOF_trial_element;j++)
	      	{
	      	  dsubgridError_h_h[j]  = - tau[0*3+0]*dpdeResidual_h_h[j]  - tau[0*3+1]*dpdeResidual_hu_h[j]  - tau[0*3+2]*dpdeResidual_hv_h[j];
	      	  dsubgridError_h_hu[j] = - tau[0*3+0]*dpdeResidual_h_hu[j] - tau[0*3+1]*dpdeResidual_hu_hu[j] - tau[0*3+2]*dpdeResidual_hv_hu[j];
	      	  dsubgridError_h_hv[j] = - tau[0*3+0]*dpdeResidual_h_hv[j] - tau[0*3+1]*dpdeResidual_hu_hv[j] - tau[0*3+2]*dpdeResidual_hv_hv[j];

	      	  dsubgridError_hu_h[j]  = - tau[1*3+0]*dpdeResidual_h_h[j]  - tau[1*3+1]*dpdeResidual_hu_h[j]  - tau[1*3+2]*dpdeResidual_hv_h[j];
	      	  dsubgridError_hu_hu[j] = - tau[1*3+0]*dpdeResidual_h_hu[j] - tau[1*3+1]*dpdeResidual_hu_hu[j] - tau[1*3+2]*dpdeResidual_hv_hu[j];
	      	  dsubgridError_hu_hv[j] = - tau[1*3+0]*dpdeResidual_h_hv[j] - tau[1*3+1]*dpdeResidual_hu_hv[j] - tau[1*3+2]*dpdeResidual_hv_hv[j];

	      	  dsubgridError_hv_h[j]  = - tau[2*3+0]*dpdeResidual_h_h[j]  - tau[2*3+1]*dpdeResidual_hu_h[j]  - tau[2*3+2]*dpdeResidual_hv_h[j];
	      	  dsubgridError_hv_hu[j] = - tau[2*3+0]*dpdeResidual_h_hu[j] - tau[2*3+1]*dpdeResidual_hu_hu[j] - tau[2*3+2]*dpdeResidual_hv_hu[j];
	      	  dsubgridError_hv_hv[j] = - tau[2*3+0]*dpdeResidual_h_hv[j] - tau[2*3+1]*dpdeResidual_hu_hv[j] - tau[2*3+2]*dpdeResidual_hv_hv[j];
	      	}
      	      //adjoint times the test functions
      	      for (int i=0;i<nDOF_test_element;i++)
      	      	{
      	      	  register int i_nSpace = i*nSpace;
      	      	  Lstar_h_h[i]=0.0;
      	      	  Lstar_hu_h[i]=ck.Advection_adjoint(dmass_adv_hu_sge,&h_grad_test_dV[i_nSpace]);
      	      	  Lstar_hv_h[i]=ck.Advection_adjoint(dmass_adv_hv_sge,&h_grad_test_dV[i_nSpace]);

      	      	  Lstar_h_hu[i]=ck.Advection_adjoint(dmom_hu_adv_h_sge,&vel_grad_test_dV[i_nSpace]) +
		    ck.Reaction_adjoint(dmom_hu_source_h,vel_test_dV[i]);
      	      	  Lstar_hu_hu[i]=ck.Advection_adjoint(dmom_hu_adv_hu_sge,&vel_grad_test_dV[i_nSpace]);
      	      	  Lstar_hv_hu[i]=ck.Advection_adjoint(dmom_hu_adv_hv_sge,&vel_grad_test_dV[i_nSpace]);

      	      	  Lstar_h_hv[i]=ck.Advection_adjoint(dmom_hv_adv_h_sge,&vel_grad_test_dV[i_nSpace])+
		    ck.Reaction_adjoint(dmom_hv_source_h,vel_test_dV[i]);
      	      	  Lstar_hu_hv[i]=ck.Advection_adjoint(dmom_hv_adv_hu_sge,&vel_grad_test_dV[i_nSpace]);
      	      	  Lstar_hv_hv[i]=ck.Advection_adjoint(dmom_hv_adv_hv_sge,&vel_grad_test_dV[i_nSpace]);
      	      	}

	      for(int i=0;i<nDOF_test_element;i++)
		{
		  register int i_nSpace = i*nSpace;
		  for(int j=0;j<nDOF_trial_element;j++) 
		    { 
		      register int j_nSpace = j*nSpace;
		      //h
		      elementJacobian_h_h[i][j] += ck.MassJacobian_weak(dmass_acc_h_t,h_trial_ref[k*nDOF_trial_element+j],h_test_dV[i]) + 
		      	ck.SubgridErrorJacobian(dsubgridError_h_h[j],Lstar_h_h[i]) +
		      	ck.SubgridErrorJacobian(dsubgridError_hu_h[j],Lstar_hu_h[i]) +
		      	ck.SubgridErrorJacobian(dsubgridError_hv_h[j],Lstar_hv_h[i]) +
			ck.NumericalDiffusionJacobian(q_numDiff_h_last[eN_k],&h_grad_trial[j_nSpace],&h_grad_test_dV[i_nSpace]);
		      
		      elementJacobian_h_hu[i][j] += ck.AdvectionJacobian_weak(dmass_adv_hu,vel_trial_ref[k*nDOF_trial_element+j],&h_grad_test_dV[i_nSpace]) +
		      	ck.SubgridErrorJacobian(dsubgridError_h_hu[j],Lstar_h_h[i]) +
		      	ck.SubgridErrorJacobian(dsubgridError_hu_hu[j],Lstar_hu_h[i]) +
		      	ck.SubgridErrorJacobian(dsubgridError_hv_hu[j],Lstar_hv_h[i]);

		      elementJacobian_h_hv[i][j] += ck.AdvectionJacobian_weak(dmass_adv_hv,vel_trial_ref[k*nDOF_trial_element+j],&h_grad_test_dV[i_nSpace]) +
		      	ck.SubgridErrorJacobian(dsubgridError_h_hv[j],Lstar_h_h[i]) +
		      	ck.SubgridErrorJacobian(dsubgridError_hu_hv[j],Lstar_hu_h[i]) +
		      	ck.SubgridErrorJacobian(dsubgridError_hv_hv[j],Lstar_hv_h[i]);

		      //u
		      elementJacobian_hu_h[i][j] += ck.AdvectionJacobian_weak(dmom_hu_adv_h,h_trial_ref[k*nDOF_trial_element+j],&vel_grad_test_dV[i_nSpace]) +
			ck.ReactionJacobian_weak(dmom_hu_source_h,h_trial_ref[k*nDOF_trial_element+j],vel_test_dV[i]) +
		      	ck.SubgridErrorJacobian(dsubgridError_h_h[j],Lstar_h_hu[i]) +
		      	ck.SubgridErrorJacobian(dsubgridError_hu_h[j],Lstar_hu_hu[i]) +
		      	ck.SubgridErrorJacobian(dsubgridError_hv_h[j],Lstar_hv_hu[i]);

		      elementJacobian_hu_hu[i][j] += ck.MassJacobian_weak(dmom_hu_acc_hu_t,vel_trial_ref[k*nDOF_trial_element+j],vel_test_dV[i]) + 
			ck.AdvectionJacobian_weak(dmom_hu_adv_hu,vel_trial_ref[k*nDOF_trial_element+j],&vel_grad_test_dV[i_nSpace]) +
			ck.SimpleDiffusionJacobian_weak(sdInfo_hu_hu_rowptr,sdInfo_hu_hu_colind,mom_hu_diff_ten,&vel_grad_trial[j_nSpace],&vel_grad_test_dV[i_nSpace]) +
		      	ck.SubgridErrorJacobian(dsubgridError_h_hu[j],Lstar_h_hu[i]) +
		      	ck.SubgridErrorJacobian(dsubgridError_hu_hu[j],Lstar_hu_hu[i]) +
		      	ck.SubgridErrorJacobian(dsubgridError_hv_hu[j],Lstar_hv_hu[i]) +
			ck.NumericalDiffusionJacobian(q_numDiff_hu_last[eN_k],&vel_grad_trial[j_nSpace],&vel_grad_test_dV[i_nSpace]);

		      elementJacobian_hu_hv[i][j] += ck.AdvectionJacobian_weak(dmom_hu_adv_hv,vel_trial_ref[k*nDOF_trial_element+j],&vel_grad_test_dV[i_nSpace]) +  
		       	ck.SimpleDiffusionJacobian_weak(sdInfo_hu_hv_rowptr,sdInfo_hu_hv_colind,mom_huhv_diff_ten,&vel_grad_trial[j_nSpace],&vel_grad_test_dV[i_nSpace]) +
		      	ck.SubgridErrorJacobian(dsubgridError_h_hv[j],Lstar_h_hu[i]) +
		      	ck.SubgridErrorJacobian(dsubgridError_hu_hv[j],Lstar_hu_hu[i]) +
		      	ck.SubgridErrorJacobian(dsubgridError_hv_hv[j],Lstar_hv_hu[i]);
		      
		      //v
		      elementJacobian_hv_h[i][j] += ck.AdvectionJacobian_weak(dmom_hv_adv_h,h_trial_ref[k*nDOF_trial_element+j],&vel_grad_test_dV[i_nSpace]) +
			ck.ReactionJacobian_weak(dmom_hv_source_h,h_trial_ref[k*nDOF_trial_element+j],vel_test_dV[i]) +
		      	ck.SubgridErrorJacobian(dsubgridError_h_h[j],Lstar_h_hv[i]) +
		      	ck.SubgridErrorJacobian(dsubgridError_hu_h[j],Lstar_hu_hv[i]) +
		      	ck.SubgridErrorJacobian(dsubgridError_hv_h[j],Lstar_hv_hv[i]);

		      elementJacobian_hv_hu[i][j] += ck.AdvectionJacobian_weak(dmom_hv_adv_hu,vel_trial_ref[k*nDOF_trial_element+j],&vel_grad_test_dV[i_nSpace]) + 
			ck.SimpleDiffusionJacobian_weak(sdInfo_hv_hu_rowptr,sdInfo_hv_hu_colind,mom_hvhu_diff_ten,&vel_grad_trial[j_nSpace],&vel_grad_test_dV[i_nSpace]) +
		      	ck.SubgridErrorJacobian(dsubgridError_h_hu[j],Lstar_h_hv[i]) +
		      	ck.SubgridErrorJacobian(dsubgridError_hu_hu[j],Lstar_hu_hv[i]) +
		      	ck.SubgridErrorJacobian(dsubgridError_hv_hu[j],Lstar_hv_hv[i]);

		      elementJacobian_hv_hv[i][j] += ck.MassJacobian_weak(dmom_hv_acc_hv_t,vel_trial_ref[k*nDOF_trial_element+j],vel_test_dV[i]) + 
			ck.AdvectionJacobian_weak(dmom_hv_adv_hv,vel_trial_ref[k*nDOF_trial_element+j],&vel_grad_test_dV[i_nSpace]) +
			ck.SimpleDiffusionJacobian_weak(sdInfo_hv_hv_rowptr,sdInfo_hv_hv_colind,mom_hv_diff_ten,&vel_grad_trial[j_nSpace],&vel_grad_test_dV[i_nSpace]) +
		      	ck.SubgridErrorJacobian(dsubgridError_h_hv[j],Lstar_h_hv[i]) +
		      	ck.SubgridErrorJacobian(dsubgridError_hu_hv[j],Lstar_hu_hv[i]) +
		      	ck.SubgridErrorJacobian(dsubgridError_hv_hv[j],Lstar_hv_hv[i]) +
		      	ck.NumericalDiffusionJacobian(q_numDiff_hv_last[eN_k],&vel_grad_trial[j_nSpace],&vel_grad_test_dV[i_nSpace]);
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
		  globalJacobian[csrRowIndeces_h_hu[eN_i] + csrColumnOffsets_h_hu[eN_i_j]] += elementJacobian_h_hu[i][j];
		  globalJacobian[csrRowIndeces_h_hv[eN_i] + csrColumnOffsets_h_hv[eN_i_j]] += elementJacobian_h_hv[i][j];

		  globalJacobian[csrRowIndeces_hu_h[eN_i] + csrColumnOffsets_hu_h[eN_i_j]] += elementJacobian_hu_h[i][j];
		  globalJacobian[csrRowIndeces_hu_hu[eN_i] + csrColumnOffsets_hu_hu[eN_i_j]] += elementJacobian_hu_hu[i][j];
		  globalJacobian[csrRowIndeces_hu_hv[eN_i] + csrColumnOffsets_hu_hv[eN_i_j]] += elementJacobian_hu_hv[i][j];

		  globalJacobian[csrRowIndeces_hv_h[eN_i] + csrColumnOffsets_hv_h[eN_i_j]] += elementJacobian_hv_h[i][j];
		  globalJacobian[csrRowIndeces_hv_hu[eN_i] + csrColumnOffsets_hv_hu[eN_i_j]] += elementJacobian_hv_hu[i][j];
		  globalJacobian[csrRowIndeces_hv_hv[eN_i] + csrColumnOffsets_hv_hv[eN_i_j]] += elementJacobian_hv_hv[i][j];
		}//j
	    }//i
	}//elements
      //
      //loop over exterior element boundaries to compute the surface integrals and load them into the global Jacobian
      //
      /* for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++)  */
      /* 	{  */
      /* 	  register int ebN = exteriorElementBoundariesArray[ebNE], */
      /* 	    eN  = elementBoundaryElementsArray[ebN*2+0], */
      /* 	    eN_nDOF_trial_element = eN*nDOF_trial_element, */
      /* 	    ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+0]; */
      /* 	  register double eps_rho,eps_mu; */
      /* 	  for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++)  */
      /* 	    {  */
      /* 	      register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb, */
      /* 		ebNE_kb_nSpace = ebNE_kb*nSpace, */
      /* 		ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb, */
      /* 		ebN_local_kb_nSpace = ebN_local_kb*nSpace; */

      /* 	      register double h_ext=0.0, */
      /* 		u_ext=0.0, */
      /* 		v_ext=0.0, */
      /* 		grad_h_ext[nSpace], */
      /* 		grad_hu_ext[nSpace], */
      /* 		grad_hv_ext[nSpace], */
      /* 		mom_hu_acc_ext=0.0, */
      /* 		dmom_hu_acc_hu_ext=0.0, */
      /* 		mom_hv_acc_ext=0.0, */
      /* 		dmom_hv_acc_hv_ext=0.0, */
      /* 		mass_adv_ext[nSpace], */
      /* 		dmass_adv_hu_ext[nSpace], */
      /* 		dmass_adv_hv_ext[nSpace], */
      /* 		mom_hu_adv_ext[nSpace], */
      /* 		dmom_hu_adv_hu_ext[nSpace], */
      /* 		dmom_hu_adv_hv_ext[nSpace], */
      /* 		mom_hv_adv_ext[nSpace], */
      /* 		dmom_hv_adv_hu_ext[nSpace], */
      /* 		dmom_hv_adv_hv_ext[nSpace], */
      /* 		mom_hu_diff_ten_ext[nSpace], */
      /* 		mom_hv_diff_ten_ext[nSpace], */
      /* 		mom_huhv_diff_ten_ext[1], */
      /* 		mom_hvhu_diff_ten_ext[1], */
      /* 		mom_hu_source_ext=0.0, */
      /* 		mom_hv_source_ext=0.0, */
      /* 		mom_hu_ham_ext=0.0, */
      /* 		dmom_hu_ham_grad_h_ext[nSpace], */
      /* 		mom_hv_ham_ext=0.0, */
      /* 		dmom_hv_ham_grad_h_ext[nSpace], */
      /* 		dmom_hu_adv_h_ext[nSpace], */
      /* 		dmom_hv_adv_h_ext[nSpace], */
      /* 		dflux_mass_hu_ext=0.0, */
      /* 		dflux_mass_hv_ext=0.0, */
      /* 		dflux_mom_hu_adv_h_ext=0.0, */
      /* 		dflux_mom_hu_adv_hu_ext=0.0, */
      /* 		dflux_mom_hu_adv_hv_ext=0.0, */
      /* 		dflux_mom_hv_adv_h_ext=0.0, */
      /* 		dflux_mom_hv_adv_hu_ext=0.0, */
      /* 		dflux_mom_hv_adv_hv_ext=0.0, */
      /* 		bc_h_ext=0.0, */
      /* 		bc_hu_ext=0.0, */
      /* 		bc_hv_ext=0.0, */
      /* 		bc_mom_hu_acc_ext=0.0, */
      /* 		bc_dmom_hu_acc_hu_ext=0.0, */
      /* 		bc_mom_hv_acc_ext=0.0, */
      /* 		bc_dmom_hv_acc_hv_ext=0.0, */
      /* 		bc_mass_adv_ext[nSpace], */
      /* 		bc_dmass_adv_hu_ext[nSpace], */
      /* 		bc_dmass_adv_hv_ext[nSpace], */
      /* 		bc_mom_hu_adv_ext[nSpace], */
      /* 		bc_dmom_hu_adv_hu_ext[nSpace], */
      /* 		bc_dmom_hu_adv_hv_ext[nSpace], */
      /* 		bc_mom_hv_adv_ext[nSpace], */
      /* 		bc_dmom_hv_adv_hu_ext[nSpace], */
      /* 		bc_dmom_hv_adv_hv_ext[nSpace], */
      /* 		bc_mom_hu_diff_ten_ext[nSpace], */
      /* 		bc_mom_hv_diff_ten_ext[nSpace], */
      /* 		bc_mom_huhv_diff_ten_ext[1], */
      /* 		bc_mom_hvhu_diff_ten_ext[1], */
      /* 		bc_mom_hu_source_ext=0.0, */
      /* 		bc_mom_hv_source_ext=0.0, */
      /* 		bc_mom_hu_ham_ext=0.0, */
      /* 		bc_dmom_hu_ham_grad_h_ext[nSpace], */
      /* 		bc_mom_hv_ham_ext=0.0, */
      /* 		bc_dmom_hv_ham_grad_h_ext[nSpace], */
      /* 		fluxJacobian_h_h[nDOF_trial_element], */
      /* 		fluxJacobian_h_hu[nDOF_trial_element], */
      /* 		fluxJacobian_h_hv[nDOF_trial_element], */
      /* 		fluxJacobian_hu_h[nDOF_trial_element], */
      /* 		fluxJacobian_hu_hu[nDOF_trial_element], */
      /* 		fluxJacobian_hu_hv[nDOF_trial_element], */
      /* 		fluxJacobian_hv_h[nDOF_trial_element], */
      /* 		fluxJacobian_hv_hu[nDOF_trial_element], */
      /* 		fluxJacobian_hv_hv[nDOF_trial_element], */
      /* 		jac_ext[nSpace*nSpace], */
      /* 		jacDet_ext, */
      /* 		jacInv_ext[nSpace*nSpace], */
      /* 		boundaryJac[nSpace*(nSpace-1)], */
      /* 		metricTensor[(nSpace-1)*(nSpace-1)], */
      /* 		metricTensorDetSqrt, */
      /* 		h_grad_trial_trace[nDOF_trial_element*nSpace], */
      /* 		vel_grad_trial_trace[nDOF_trial_element*nSpace], */
      /* 		dS, */
      /* 		h_test_dS[nDOF_test_element], */
      /* 		vel_test_dS[nDOF_test_element], */
      /* 		normal[3], */
      /* 		x_ext,y_ext,xt_ext,yt_ext,integralScaling, */
      /* 		G[nSpace*nSpace],G_dd_G,tr_G,h_hhi,h_henalty; */
      /* 	      ck.calculateMapping_elementBoundary(eN, */
      /* 						  ebN_local, */
      /* 						  kb, */
      /* 						  ebN_local_kb, */
      /* 						  mesh_dof, */
      /* 						  mesh_l2g, */
      /* 						  mesh_trial_trace_ref, */
      /* 						  mesh_grad_trial_trace_ref, */
      /* 						  boundaryJac_ref, */
      /* 						  jac_ext, */
      /* 						  jacDet_ext, */
      /* 						  jacInv_ext, */
      /* 						  boundaryJac, */
      /* 						  metricTensor, */
      /* 						  metricTensorDetSqrt, */
      /* 						  normal_ref, */
      /* 						  normal, */
      /* 						  x_ext,y_ext); */
      /* 	      ck.calculateMappingVelocity_elementBoundary(eN, */
      /* 							  ebN_local, */
      /* 							  kb, */
      /* 							  ebN_local_kb, */
      /* 							  mesh_velocity_dof, */
      /* 							  mesh_l2g, */
      /* 							  mesh_trial_trace_ref, */
      /* 							  xt_ext,yt_ext, */
      /* 							  normal, */
      /* 							  boundaryJac, */
      /* 							  metricTensor, */
      /* 							  integralScaling); */
      /* 	      //xt_ext=0.0;yt_ext=0.0; */
      /* 	      //std::cout<<"xt_ext "<<xt_ext<<'\t'<<yt_ext<<'\t'<<std::endl; */
      /* 	      dS = ((1.0-MOVING_DOMAIN)*metricTensorDetSqrt + MOVING_DOMAIN*integralScaling)*dS_ref[kb]; */
      /* 	      ck.calculateG(jacInv_ext,G,G_dd_G,tr_G); */
      /* 	      ck.calculateGScale(G,&ebqe_normal_hhi_ext[ebNE_kb_nSpace],h_hhi); */

      /* 	      eps_rho = epsFact_rho*(useMetrics*h_hhi+(1.0-useMetrics)*elementDiameter[eN]); */
      /* 	      eps_mu  = epsFact_mu *(useMetrics*h_hhi+(1.0-useMetrics)*elementDiameter[eN]); */

      /* 	      //compute shape and solution information */
      /* 	      //shape */
      /* 	      ck.gradTrialFromRef(&h_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,h_grad_trial_trace); */
      /* 	      ck.gradTrialFromRef(&vel_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,vel_grad_trial_trace); */
      /* 	      //solution and gradients	 */
      /* 	      ck.valFromDOF(h_dof,&h_l2g[eN_nDOF_trial_element],&h_trial_trace_ref[ebN_local_kb*nDOF_test_element],h_ext); */
      /* 	      ck.valFromDOF(hu_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_trace_ref[ebN_local_kb*nDOF_test_element],u_ext); */
      /* 	      ck.valFromDOF(hv_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_trace_ref[ebN_local_kb*nDOF_test_element],v_ext); */
      /* 	      ck.gradFromDOF(h_dof,&h_l2g[eN_nDOF_trial_element],h_grad_trial_trace,grad_h_ext); */
      /* 	      ck.gradFromDOF(hu_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial_trace,grad_hu_ext); */
      /* 	      ck.gradFromDOF(hv_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial_trace,grad_hv_ext); */
      /* 	      //precalculate test function products with integration weights */
      /* 	      for (int j=0;j<nDOF_trial_element;j++) */
      /* 		{ */
      /* 		  h_test_dS[j] = h_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS; */
      /* 		  vel_test_dS[j] = vel_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS; */
      /* 		} */
      /* 	      // // */
      /* 	      // //debugging section for finite element calculations on exterior */
      /* 	      // // */
      /* 	      // std::cout<<"ebNE = "<<ebNE<<" kb = "<<kb<<std::endl; */
      /* 	      // for (int j=0;j<nDOF_trial_element;j++) */
      /* 	      //   { */
      /* 	      //     std::cout<<"h_trial_trace["<<j<<"] "<<h_trial_trace_ref[ebN_local_kb*nDOF_trial_element+j]<<std::endl; */
      /* 	      //     std::cout<<"vel_trial_trace["<<j<<"] "<<vel_trial_trace_ref[ebN_local_kb*nDOF_trial_element+j]<<std::endl; */
      /* 	      //     std::cout<<"h_test_dS["<<j<<"] "<<h_test_dS[j]<<std::endl; */
      /* 	      //     std::cout<<"vel_test_dS["<<j<<"] "<<vel_test_dS[j]<<std::endl; */
      /* 	      //     for (int I=0;I<nSpace;I++) */
      /* 	      // 	{ */
      /* 	      // 	  std::cout<<"h_grad_trial_trace["<<j<<","<<I<<"] "<<h_grad_trial_trace[j*nSpace+I]<<std::endl; */
      /* 	      // 	  std::cout<<"vel_grad_trial_trace["<<j<<","<<I<<"] "<<vel_grad_trial_trace[j*nSpace+I]<<std::endl; */
      /* 	      // 	} */
      /* 	      //   } */
      /* 	      // std::cout<<"h_ext "<<h_ext<<std::endl; */
      /* 	      // std::cout<<"u_ext "<<u_ext<<std::endl; */
      /* 	      // std::cout<<"v_ext "<<v_ext<<std::endl; */
      /* 	      // for(int I=0;I<nSpace;I++) */
      /* 	      //   { */
      /* 	      //     std::cout<<"grad_h_ext["<<I<<"] "<<grad_h_ext[I]<<std::endl; */
      /* 	      //     std::cout<<"grad_hu_ext["<<I<<"] "<<grad_hu_ext[I]<<std::endl; */
      /* 	      //     std::cout<<"grad_hv_ext["<<I<<"] "<<grad_hv_ext[I]<<std::endl; */
      /* 	      //   } */
      /* 	      // */
      /* 	      //load the boundary values */
      /* 	      // */
      /* 	      bc_h_ext = isDOFBoundary_h[ebNE_kb]*ebqe_bc_h_ext[ebNE_kb]+(1-isDOFBoundary_h[ebNE_kb])*h_ext; */
      /* 	      bc_hu_ext = isDOFBoundary_hu[ebNE_kb]*ebqe_bc_hu_ext[ebNE_kb]+(1-isDOFBoundary_hu[ebNE_kb])*u_ext; */
      /* 	      bc_hv_ext = isDOFBoundary_hv[ebNE_kb]*ebqe_bc_hv_ext[ebNE_kb]+(1-isDOFBoundary_hv[ebNE_kb])*v_ext; */
      /* 	      //  */
      /* 	      //calculate the internal and external trace of the pde coefficients  */
      /* 	      //  */
      /* 	      //cek debug */
      /* 	      //eps_rho=0.1; */
      /* 	      //eps_mu=0.1; */
      /* 	      evaluateCoefficients(eps_rho, */
      /* 				   eps_mu, */
      /* 				   sigma, */
      /* 				   rho_0, */
      /* 				   nu_0, */
      /* 				   rho_1, */
      /* 				   nu_1, */
      /* 				   g, */
      /* 				   ebqe_hhi_ext[ebNE_kb], */
      /* 				   &ebqe_normal_hhi_ext[ebNE_kb_nSpace], */
      /* 				   ebqe_kappa_hhi_ext[ebNE_kb], */
      /* 				   h_ext, */
      /* 				   grad_h_ext, */
      /* 				   u_ext, */
      /* 				   v_ext, */
      /* 				   mom_hu_acc_ext, */
      /* 				   dmom_hu_acc_hu_ext, */
      /* 				   mom_hv_acc_ext, */
      /* 				   dmom_hv_acc_hv_ext, */
      /* 				   mass_adv_ext, */
      /* 				   dmass_adv_hu_ext, */
      /* 				   dmass_adv_hv_ext, */
      /* 				   mom_hu_adv_ext, */
      /* 				   dmom_hu_adv_hu_ext, */
      /* 				   dmom_hu_adv_hv_ext, */
      /* 				   mom_hv_adv_ext, */
      /* 				   dmom_hv_adv_hu_ext, */
      /* 				   dmom_hv_adv_hv_ext, */
      /* 				   mom_hu_diff_ten_ext, */
      /* 				   mom_hv_diff_ten_ext, */
      /* 				   mom_huhv_diff_ten_ext, */
      /* 				   mom_hvhu_diff_ten_ext, */
      /* 				   mom_hu_source_ext, */
      /* 				   mom_hv_source_ext, */
      /* 				   mom_hu_ham_ext, */
      /* 				   dmom_hu_ham_grad_h_ext, */
      /* 				   mom_hv_ham_ext, */
      /* 				   dmom_hv_ham_grad_h_ext);           */
      /* 	      evaluateCoefficients(eps_rho, */
      /* 				   eps_mu, */
      /* 				   sigma, */
      /* 				   rho_0, */
      /* 				   nu_0, */
      /* 				   rho_1, */
      /* 				   nu_1, */
      /* 				   g, */
      /* 				   ebqe_hhi_ext[ebNE_kb], */
      /* 				   &ebqe_normal_hhi_ext[ebNE_kb_nSpace], */
      /* 				   ebqe_kappa_hhi_ext[ebNE_kb], */
      /* 				   bc_h_ext, */
      /* 				   grad_h_ext, //cek shouldn't be used */
      /* 				   bc_hu_ext, */
      /* 				   bc_hv_ext, */
      /* 				   bc_mom_hu_acc_ext, */
      /* 				   bc_dmom_hu_acc_hu_ext, */
      /* 				   bc_mom_hv_acc_ext, */
      /* 				   bc_dmom_hv_acc_hv_ext, */
      /* 				   bc_mass_adv_ext, */
      /* 				   bc_dmass_adv_hu_ext, */
      /* 				   bc_dmass_adv_hv_ext, */
      /* 				   bc_mom_hu_adv_ext, */
      /* 				   bc_dmom_hu_adv_hu_ext, */
      /* 				   bc_dmom_hu_adv_hv_ext, */
      /* 				   bc_mom_hv_adv_ext, */
      /* 				   bc_dmom_hv_adv_hu_ext, */
      /* 				   bc_dmom_hv_adv_hv_ext, */
      /* 				   bc_mom_hu_diff_ten_ext, */
      /* 				   bc_mom_hv_diff_ten_ext, */
      /* 				   bc_mom_huhv_diff_ten_ext, */
      /* 				   bc_mom_hvhu_diff_ten_ext, */
      /* 				   bc_mom_hu_source_ext, */
      /* 				   bc_mom_hv_source_ext, */
      /* 				   bc_mom_hu_ham_ext, */
      /* 				   bc_dmom_hu_ham_grad_h_ext, */
      /* 				   bc_mom_hv_ham_ext, */
      /* 				   bc_dmom_hv_ham_grad_h_ext);           */
      /* 	      // */
      /* 	      //moving domain */
      /* 	      // */
      /* 	      mass_adv_ext[0] -= MOVING_DOMAIN*xt_ext; */
      /* 	      mass_adv_ext[1] -= MOVING_DOMAIN*yt_ext; */

      /* 	      mom_hu_adv_ext[0] -= MOVING_DOMAIN*mom_hu_acc_ext*xt_ext; */
      /* 	      mom_hu_adv_ext[1] -= MOVING_DOMAIN*mom_hu_acc_ext*yt_ext; */
      /* 	      dmom_hu_adv_hu_ext[0] -= MOVING_DOMAIN*dmom_hu_acc_hu_ext*xt_ext; */
      /* 	      dmom_hu_adv_hu_ext[1] -= MOVING_DOMAIN*dmom_hu_acc_hu_ext*yt_ext; */
	      
      /* 	      mom_hv_adv_ext[0] -= MOVING_DOMAIN*mom_hv_acc_ext*xt_ext; */
      /* 	      mom_hv_adv_ext[1] -= MOVING_DOMAIN*mom_hv_acc_ext*yt_ext; */
      /* 	      dmom_hv_adv_hv_ext[0] -= MOVING_DOMAIN*dmom_hv_acc_hv_ext*xt_ext; */
      /* 	      dmom_hv_adv_hv_ext[1] -= MOVING_DOMAIN*dmom_hv_acc_hv_ext*yt_ext; */
	      
      /* 	      //moving domain bc's */
      /* 	      bc_mom_hu_adv_ext[0] -= MOVING_DOMAIN*bc_mom_hu_acc_ext*xt_ext; */
      /* 	      bc_mom_hu_adv_ext[1] -= MOVING_DOMAIN*bc_mom_hu_acc_ext*yt_ext; */
	      
      /* 	      bc_mom_hv_adv_ext[0] -= MOVING_DOMAIN*bc_mom_hv_acc_ext*xt_ext; */
      /* 	      bc_mom_hv_adv_ext[1] -= MOVING_DOMAIN*bc_mom_hv_acc_ext*yt_ext; */

      /* 	      //  */
      /* 	      //calculate the numerical fluxes  */
      /* 	      //  */
      /* 	      exteriorNumericalAdvectiveFluxDerivatives(isDOFBoundary_h[ebNE_kb], */
      /* 							isDOFBoundary_hu[ebNE_kb], */
      /* 							isDOFBoundary_hv[ebNE_kb], */
      /* 							isAdvectiveFluxBoundary_h[ebNE_kb], */
      /* 							isAdvectiveFluxBoundary_hu[ebNE_kb], */
      /* 							isAdvectiveFluxBoundary_hv[ebNE_kb], */
      /* 							dmom_hu_ham_grad_h_ext[0],//=1/rho */
      /* 							normal, */
      /* 							bc_h_ext, */
      /* 							bc_mass_adv_ext, */
      /* 							bc_mom_hu_adv_ext, */
      /* 							bc_mom_hv_adv_ext, */
      /* 							ebqe_bc_flux_mass_ext[ebNE_kb], */
      /* 							ebqe_bc_flux_mom_hu_adv_ext[ebNE_kb], */
      /* 							ebqe_bc_flux_mom_hv_adv_ext[ebNE_kb], */
      /* 							h_ext, */
      /* 							mass_adv_ext, */
      /* 							mom_hu_adv_ext, */
      /* 							mom_hv_adv_ext, */
      /* 							dmass_adv_hu_ext, */
      /* 							dmass_adv_hv_ext, */
      /* 							dmom_hu_adv_h_ext, */
      /* 							dmom_hu_adv_hu_ext, */
      /* 							dmom_hu_adv_hv_ext, */
      /* 							dmom_hv_adv_h_ext, */
      /* 							dmom_hv_adv_hu_ext, */
      /* 							dmom_hv_adv_hv_ext, */
      /* 							dflux_mass_hu_ext, */
      /* 							dflux_mass_hv_ext, */
      /* 							dflux_mom_hu_adv_h_ext, */
      /* 							dflux_mom_hu_adv_hu_ext, */
      /* 							dflux_mom_hu_adv_hv_ext, */
      /* 							dflux_mom_hv_adv_h_ext, */
      /* 							dflux_mom_hv_adv_hu_ext, */
      /* 							dflux_mom_hv_adv_hv_ext); */
      /* 	      // */
      /* 	      //calculate the flux jacobian */
      /* 	      // */
      /* 	      ck.calculateGScale(G,normal,h_henalty); */
      /* 	      h_henalty = 10.0/h_henalty; */
      /* 	      //cek debug, do it the old way */
      /* 	      h_henalty = 100.0/elementDiameter[eN]; */
      /* 	      for (int j=0;j<nDOF_trial_element;j++) */
      /* 		{ */
      /* 		  register int j_nSpace = j*nSpace,ebN_local_kb_j=ebN_local_kb*nDOF_trial_element+j; */
      /* 		  //cek debug */
      /* 		  //ebqe_henalty_ext[ebNE_kb] = 10.0; */
      /* 		  // */
      /* 		  //cek todo add full stress on boundaries */

      /* 		  fluxJacobian_h_h[j]=0.0; */
      /* 		  fluxJacobian_h_hu[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mass_hu_ext,vel_trial_trace_ref[ebN_local_kb_j]); */
      /* 		  fluxJacobian_h_hv[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mass_hv_ext,vel_trial_trace_ref[ebN_local_kb_j]); */

      /* 		  fluxJacobian_hu_h[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_hu_adv_h_ext,h_trial_trace_ref[ebN_local_kb_j]); */
      /* 		  fluxJacobian_hu_hu[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_hu_adv_hu_ext,vel_trial_trace_ref[ebN_local_kb_j]) + */
      /* 		    ExteriorNumericalDiffusiveFluxJacobian(eps_rho, */
      /* 							   ebqe_hhi_ext[ebNE_kb], */
      /* 							   sdInfo_hu_hu_rowptr, */
      /* 							   sdInfo_hu_hu_colind, */
      /* 							   isDOFBoundary_hu[ebNE_kb], */
      /* 							   normal, */
      /* 							   mom_hu_diff_ten_ext, */
      /* 							   vel_trial_trace_ref[ebN_local_kb_j], */
      /* 							   &vel_grad_trial_trace[j_nSpace], */
      /* 							   h_henalty);//ebqe_henalty_ext[ebNE_kb]); */
      /* 		  fluxJacobian_hu_hv[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_hu_adv_hv_ext,vel_trial_trace_ref[ebN_local_kb_j]); */

      /* 		  fluxJacobian_hv_h[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_hv_adv_h_ext,h_trial_trace_ref[ebN_local_kb_j]); */
      /* 		  fluxJacobian_hv_hu[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_hv_adv_hu_ext,vel_trial_trace_ref[ebN_local_kb_j]); */
      /* 		  fluxJacobian_hv_hv[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_hv_adv_hv_ext,vel_trial_trace_ref[ebN_local_kb_j]) + */
      /* 		    ExteriorNumericalDiffusiveFluxJacobian(eps_rho, */
      /* 							   ebqe_hhi_ext[ebNE_kb], */
      /* 							   sdInfo_hv_hv_rowptr, */
      /* 							   sdInfo_hv_hv_colind, */
      /* 							   isDOFBoundary_hv[ebNE_kb], */
      /* 							   normal, */
      /* 							   mom_hv_diff_ten_ext, */
      /* 							   vel_trial_trace_ref[ebN_local_kb_j], */
      /* 							   &vel_grad_trial_trace[j_nSpace], */
      /* 							   h_henalty);//ebqe_henalty_ext[ebNE_kb]); */

      /* 		  fluxJacobian_h_h[j]=0.0; */
      /* 		  fluxJacobian_h_hu[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mass_hu_ext,vel_trial_trace_ref[ebN_local_kb_j]); */
      /* 		  fluxJacobian_h_hv[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mass_hv_ext,vel_trial_trace_ref[ebN_local_kb_j]); */

      /* 		  fluxJacobian_hu_h[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_hu_adv_h_ext,h_trial_trace_ref[ebN_local_kb_j]); */
      /* 		  fluxJacobian_hu_hu[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_hu_adv_hu_ext,vel_trial_trace_ref[ebN_local_kb_j]) + */
      /* 		    ExteriorNumericalDiffusiveFluxJacobian(eps_rho, */
      /* 							   ebqe_hhi_ext[ebNE_kb], */
      /* 							   sdInfo_hu_hu_rowptr, */
      /* 							   sdInfo_hu_hu_colind, */
      /* 							   isDOFBoundary_hu[ebNE_kb], */
      /* 							   normal, */
      /* 							   mom_hu_diff_ten_ext, */
      /* 							   vel_trial_trace_ref[ebN_local_kb_j], */
      /* 							   &vel_grad_trial_trace[j_nSpace], */
      /* 							   h_henalty);//ebqe_henalty_ext[ebNE_kb]); */
      /* 		  fluxJacobian_hu_hv[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_hu_adv_hv_ext,vel_trial_trace_ref[ebN_local_kb_j]); */

      /* 		  fluxJacobian_hv_h[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_hv_adv_h_ext,h_trial_trace_ref[ebN_local_kb_j]); */
      /* 		  fluxJacobian_hv_hu[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_hv_adv_hu_ext,vel_trial_trace_ref[ebN_local_kb_j]); */
      /* 		  fluxJacobian_hv_hv[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_hv_adv_hv_ext,vel_trial_trace_ref[ebN_local_kb_j]) + */
      /* 		    ExteriorNumericalDiffusiveFluxJacobian(eps_rho, */
      /* 							   ebqe_hhi_ext[ebNE_kb], */
      /* 							   sdInfo_hv_hv_rowptr, */
      /* 							   sdInfo_hv_hv_colind, */
      /* 							   isDOFBoundary_hv[ebNE_kb], */
      /* 							   normal, */
      /* 							   mom_hv_diff_ten_ext, */
      /* 							   vel_trial_trace_ref[ebN_local_kb_j], */
      /* 							   &vel_grad_trial_trace[j_nSpace], */
      /* 							   h_henalty);//ebqe_henalty_ext[ebNE_kb]); */
      /* 		  // //cek debug */
      /* 		  // fluxJacobian_h_h[j]=0.0; */
      /* 		  // fluxJacobian_h_hu[j]=0.0; */
      /* 		  // fluxJacobian_h_hv[j]=0.0; */

      /* 		  // fluxJacobian_hu_h[j]=0.0; */
      /* 		  // fluxJacobian_hu_hu[j]=0.0; */
      /* 		  // fluxJacobian_hu_hv[j]=0.0; */

      /* 		  // fluxJacobian_hv_h[j]=0.0; */
      /* 		  // fluxJacobian_hv_hu[j]=0.0; */
      /* 		  // fluxJacobian_hv_hv[j]=0.0; */

      /* 		  // //cek debug */
      /* 		}//j */
      /* 	      // */
      /* 	      //update the global Jacobian from the flux Jacobian */
      /* 	      // */
      /* 	      for (int i=0;i<nDOF_test_element;i++) */
      /* 		{ */
      /* 		  register int eN_i = eN*nDOF_test_element+i; */
      /* 		  for (int j=0;j<nDOF_trial_element;j++) */
      /* 		    { */
      /* 		      register int ebN_i_j = ebN*4*nDOF_test_X_trial_element + i*nDOF_trial_element + j; */
		  
      /* 		      globalJacobian[csrRowIndeces_h_h[eN_i] + csrColumnOffsets_eb_h_h[ebN_i_j]] += fluxJacobian_h_h[j]*h_test_dS[i]; */
      /* 		      globalJacobian[csrRowIndeces_h_hu[eN_i] + csrColumnOffsets_eb_h_hu[ebN_i_j]] += fluxJacobian_h_hu[j]*h_test_dS[i]; */
      /* 		      globalJacobian[csrRowIndeces_h_hv[eN_i] + csrColumnOffsets_eb_h_hv[ebN_i_j]] += fluxJacobian_h_hv[j]*h_test_dS[i]; */
		   
      /* 		      globalJacobian[csrRowIndeces_hu_h[eN_i] + csrColumnOffsets_eb_hu_h[ebN_i_j]] += fluxJacobian_hu_h[j]*vel_test_dS[i]; */
      /* 		      globalJacobian[csrRowIndeces_hu_hu[eN_i] + csrColumnOffsets_eb_hu_hu[ebN_i_j]] += fluxJacobian_hu_hu[j]*vel_test_dS[i]; */
      /* 		      globalJacobian[csrRowIndeces_hu_hv[eN_i] + csrColumnOffsets_eb_hu_hv[ebN_i_j]] += fluxJacobian_hu_hv[j]*vel_test_dS[i]; */
		   
      /* 		      globalJacobian[csrRowIndeces_hv_h[eN_i] + csrColumnOffsets_eb_hv_h[ebN_i_j]] += fluxJacobian_hv_h[j]*vel_test_dS[i]; */
      /* 		      globalJacobian[csrRowIndeces_hv_hu[eN_i] + csrColumnOffsets_eb_hv_hu[ebN_i_j]] += fluxJacobian_hv_hu[j]*vel_test_dS[i]; */
      /* 		      globalJacobian[csrRowIndeces_hv_hv[eN_i] + csrColumnOffsets_eb_hv_hv[ebN_i_j]] += fluxJacobian_hv_hv[j]*vel_test_dS[i]; */
		   
      /* 		    }//j */
      /* 		}//i */
      /* 	      // //debug */
      /* 	      // std::cout<<"flux jacobian ebNE "<<ebNE<<" kb "<<kb<<std::endl; */
      /* 	      // for (int i=0;i<nDOF_test_element;i++) */
      /* 	      //   { */
      /* 	      //     for (int j=0;j<nDOF_trial_element;j++) */
      /* 	      // 	{ */
      /* 	      // 	  std::cout<< fluxJacobian_h_h[j]*h_test_dS[i]<<std::endl; */
      /* 	      // 	  std::cout<< fluxJacobian_h_hu[j]*h_test_dS[i]<<std::endl; */
      /* 	      // 	  std::cout<< fluxJacobian_h_hv[j]*h_test_dS[i]<<std::endl; */
		  
      /* 	      // 	  std::cout<< fluxJacobian_hu_h[j]*vel_test_dS[i]<<std::endl; */
      /* 	      // 	  std::cout<< fluxJacobian_hu_hu[j]*vel_test_dS[i]<<std::endl; */
      /* 	      // 	  std::cout<< fluxJacobian_hu_hv[j]*vel_test_dS[i]<<std::endl; */
		  
      /* 	      // 	  std::cout<< fluxJacobian_hv_h[j]*vel_test_dS[i]<<std::endl; */
      /* 	      // 	  std::cout<< fluxJacobian_hv_hu[j]*vel_test_dS[i]<<std::endl; */
      /* 	      // 	  std::cout<< fluxJacobian_hv_hv[j]*vel_test_dS[i]<<std::endl; */
      /* 	      // 	}//j */
      /* 	      //   }//i */
      /* 	    }//kb */
      /* 	}//ebNE */
    }//computeJacobian

    /* void calculateVelocityAverage(int nExteriorElementBoundaries_global, */
    /* 				  int* exteriorElementBoundariesArray, */
    /* 				  int nInteriorElementBoundaries_global, */
    /* 				  int* interiorElementBoundariesArray, */
    /* 				  int* elementBoundaryElementsArray, */
    /* 				  int* elementBoundaryLocalElementBoundariesArray, */
    /* 				  double* mesh_dof, */
    /* 				  int* mesh_l2g, */
    /* 				  double* mesh_trial_trace_ref, */
    /* 				  double* mesh_grad_trial_trace_ref, */
    /* 				  double* normal_ref, */
    /* 				  double* boundaryJac_ref, */
    /* 				  int* vel_l2g, */
    /* 				  double* hu_dof, */
    /* 				  double* hv_dof, */
    /* 				  double* vel_trial_trace_ref, */
    /* 				  double* ebqe_velocity, */
    /* 				  double* velocityAverage) */
    /* { */
    /*   int permutations[nQuadraturePoints_elementBoundary]; */
    /*   double xArray_left[nQuadraturePoints_elementBoundary*3], */
    /* 	xArray_right[nQuadraturePoints_elementBoundary*3]; */
    /*   for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) */
    /* 	{ */
    /* 	  register int ebN = exteriorElementBoundariesArray[ebNE]; */
    /* 	  for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) */
    /* 	    { */
    /* 	      register int ebN_kb_nSpace = ebN*nQuadraturePoints_elementBoundary*nSpace+kb*nSpace, */
    /* 		ebNE_kb_nSpace = ebNE*nQuadraturePoints_elementBoundary*nSpace+kb*nSpace; */
    /* 	      velocityAverage[ebN_kb_nSpace+0]=ebqe_velocity[ebNE_kb_nSpace+0]; */
    /* 	      velocityAverage[ebN_kb_nSpace+1]=ebqe_velocity[ebNE_kb_nSpace+1]; */
    /* 	      velocityAverage[ebN_kb_nSpace+2]=ebqe_velocity[ebNE_kb_nSpace+2]; */
    /* 	    }//ebNE */
    /* 	} */
    /*   for (int ebNI = 0; ebNI < nInteriorElementBoundaries_global; ebNI++) */
    /* 	{ */
    /* 	  register int ebN = interiorElementBoundariesArray[ebNI], */
    /* 	    left_eN_global   = elementBoundaryElementsArray[ebN*2+0], */
    /* 	    left_ebN_element  = elementBoundaryLocalElementBoundariesArray[ebN*2+0], */
    /* 	    right_eN_global  = elementBoundaryElementsArray[ebN*2+1], */
    /* 	    right_ebN_element = elementBoundaryLocalElementBoundariesArray[ebN*2+1], */
    /* 	    left_eN_nDOF_trial_element = left_eN_global*nDOF_trial_element, */
    /* 	    right_eN_nDOF_trial_element = right_eN_global*nDOF_trial_element; */
    /* 	  double jac[nSpace*nSpace], */
    /* 	    jacDet, */
    /* 	    jacInv[nSpace*nSpace], */
    /* 	    boundaryJac[nSpace*(nSpace-1)], */
    /* 	    metricTensor[(nSpace-1)*(nSpace-1)], */
    /* 	    metricTensorDetSqrt, */
    /* 	    normal[3], */
    /* 	    x,y; */
    /* 	  //double G[nSpace*nSpace],G_dd_G,tr_G,h_hhi,h_henalty; */
	  
    /* 	  for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) */
    /* 	    { */
    /* 	      ck.calculateMapping_elementBoundary(left_eN_global, */
    /* 						  left_ebN_element, */
    /* 						  kb, */
    /* 						  left_ebN_element*kb, */
    /* 						  mesh_dof, */
    /* 						  mesh_l2g, */
    /* 						  mesh_trial_trace_ref, */
    /* 						  mesh_grad_trial_trace_ref, */
    /* 						  boundaryJac_ref, */
    /* 						  jac, */
    /* 						  jacDet, */
    /* 						  jacInv, */
    /* 						  boundaryJac, */
    /* 						  metricTensor, */
    /* 						  metricTensorDetSqrt, */
    /* 						  normal_ref, */
    /* 						  normal, */
    /* 						  x,y); */
    /* 	      xArray_left[kb*3+0] = x; */
    /* 	      xArray_left[kb*3+1] = y; */
    /* 	      ck.calculateMapping_elementBoundary(right_eN_global, */
    /* 						  right_ebN_element, */
    /* 						  kb, */
    /* 						  right_ebN_element*kb, */
    /* 						  mesh_dof, */
    /* 						  mesh_l2g, */
    /* 						  mesh_trial_trace_ref, */
    /* 						  mesh_grad_trial_trace_ref, */
    /* 						  boundaryJac_ref, */
    /* 						  jac, */
    /* 						  jacDet, */
    /* 						  jacInv, */
    /* 						  boundaryJac, */
    /* 						  metricTensor, */
    /* 						  metricTensorDetSqrt, */
    /* 						  normal_ref, */
    /* 						  normal, */
    /* 						  x,y); */
    /* 	      xArray_right[kb*3+0] = x; */
    /* 	      xArray_right[kb*3+1] = y; */
    /* 	    } */
    /* 	  for  (int kb_left=0;kb_left<nQuadraturePoints_elementBoundary;kb_left++) */
    /* 	    { */
    /* 	      double errorNormMin = 1.0; */
    /* 	      for  (int kb_right=0;kb_right<nQuadraturePoints_elementBoundary;kb_right++) */
    /* 		{ */
    /* 		  double errorNorm=0.0; */
    /* 		  for (int I=0;I<nSpace;I++) */
    /* 		    { */
    /* 		      errorNorm += fabs(xArray_left[kb_left*3+I] */
    /* 					- */
    /* 					xArray_right[kb_right*3+I]); */
    /* 		    } */
    /* 		  if (errorNorm < errorNormMin) */
    /* 		    { */
    /* 		      permutations[kb_right] = kb_left; */
    /* 		      errorNormMin = errorNorm; */
    /* 		    } */
    /* 		} */
    /* 	    } */
    /* 	  for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) */
    /* 	    { */
    /* 	      register int ebN_kb_nSpace = ebN*nQuadraturePoints_elementBoundary*nSpace+kb*nSpace; */
    /* 	      register double u_left=0.0, */
    /* 		v_left=0.0, */
    /* 		u_right=0.0, */
    /* 		v_right=0.0; */
    /* 	      register int left_kb = kb, */
    /* 		right_kb = permutations[kb], */
    /* 		left_ebN_element_kb_nDOF_test_element=left_ebN_element*left_kb*nDOF_test_element, */
    /* 		right_ebN_element_kb_nDOF_test_element=right_ebN_element*right_kb*nDOF_test_element; */
    /* 	      // */
    /* 	      //calculate the velocity solution at quadrature points on left and right */
    /* 	      // */
    /* 	      ck.valFromDOF(hu_dof,&vel_l2g[left_eN_nDOF_trial_element],&vel_trial_trace_ref[left_ebN_element_kb_nDOF_test_element],u_left); */
    /* 	      ck.valFromDOF(hv_dof,&vel_l2g[left_eN_nDOF_trial_element],&vel_trial_trace_ref[left_ebN_element_kb_nDOF_test_element],v_left); */
    /* 	      // */
    /* 	      ck.valFromDOF(hu_dof,&vel_l2g[right_eN_nDOF_trial_element],&vel_trial_trace_ref[right_ebN_element_kb_nDOF_test_element],u_right); */
    /* 	      ck.valFromDOF(hv_dof,&vel_l2g[right_eN_nDOF_trial_element],&vel_trial_trace_ref[right_ebN_element_kb_nDOF_test_element],v_right); */
    /* 	      // */
    /* 	      velocityAverage[ebN_kb_nSpace+0]=0.5*(u_left + u_right); */
    /* 	      velocityAverage[ebN_kb_nSpace+1]=0.5*(v_left + v_right); */
    /* 	    }//ebNI */
    /* 	} */
    /* } */

    void calculateJacobian_cell_based_entropy_viscosity(//element
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
							double* q_velocity_sge,
							double* q_cfl,
							double* q_numDiff_h_last,
							double* q_numDiff_hu_last, 
							double* q_numDiff_hv_last, 
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
	  register double  elementJacobian_h_h[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_h_hu[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_h_hv[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_hu_h[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_hu_hu[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_hu_hv[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_hv_h[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_hv_hu[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_hv_hv[nDOF_test_element][nDOF_trial_element];
	  for (int i=0;i<nDOF_test_element;i++)
	    for (int j=0;j<nDOF_trial_element;j++)
	      {
		elementJacobian_h_h[i][j]=0.0;
		elementJacobian_h_hu[i][j]=0.0;
		elementJacobian_h_hv[i][j]=0.0;
		elementJacobian_hu_h[i][j]=0.0;
		elementJacobian_hu_hu[i][j]=0.0;
		elementJacobian_hu_hv[i][j]=0.0;
		elementJacobian_hv_h[i][j]=0.0;
		elementJacobian_hv_hu[i][j]=0.0;
		elementJacobian_hv_hv[i][j]=0.0;
	      }
	  for  (int k=0;k<nQuadraturePoints_element;k++)
	    {
	      int eN_k = eN*nQuadraturePoints_element+k, //index to a scalar at a quadrature point
		eN_k_nSpace = eN_k*nSpace,
		eN_nDOF_trial_element = eN*nDOF_trial_element; //index to a vector at a quadrature point

	      //declare local storage
	      register double b=0.0,
		h=0.0,
		hu=0.0,
		hv=0.0,
		grad_b[nSpace],
		grad_h[nSpace],
		grad_hu[nSpace],
		grad_hv[nSpace],
		mass_acc=0.0,
		dmass_acc_h=0.0,
		mom_hu_acc=0.0,
		dmom_hu_acc_h=0.0,
		dmom_hu_acc_hu=0.0,
		mom_hv_acc=0.0,
		dmom_hv_acc_h=0.0,
		dmom_hv_acc_hv=0.0,
		mass_adv[nSpace],
		dmass_adv_h[nSpace],
		dmass_adv_hu[nSpace],
		dmass_adv_hv[nSpace],
		mom_hu_adv[nSpace],
		dmom_hu_adv_h[nSpace],
		dmom_hu_adv_hu[nSpace],
		dmom_hu_adv_hv[nSpace],
		mom_hv_adv[nSpace],
		dmom_hv_adv_h[nSpace],
		dmom_hv_adv_hu[nSpace],
		dmom_hv_adv_hv[nSpace],
		mom_hu_source=0.0,
		dmom_hu_source_h=0.0,
		mom_hv_source=0.0,
		dmom_hv_source_h=0.0,
		mass_acc_t=0.0,
		dmass_acc_h_t=0.0,
		mom_hu_acc_t=0.0,
		dmom_hu_acc_h_t=0.0,
		dmom_hu_acc_hu_t=0.0,
		mom_hv_acc_t=0.0,
		dmom_hv_acc_h_t=0.0,
		dmom_hv_acc_hv_t=0.0,
		jac[nSpace*nSpace],
		jacDet,
		jacInv[nSpace*nSpace],
		h_grad_trial[nDOF_trial_element*nSpace],
		vel_grad_trial[nDOF_trial_element*nSpace],
		dV,
		h_test_dV[nDOF_test_element],
		vel_test_dV[nDOF_test_element],
		h_grad_test_dV[nDOF_test_element*nSpace],
		vel_grad_test_dV[nDOF_test_element*nSpace],
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
	      //ck.calculateMappingVelocity_element(eN,
	      //				  k,
	      //				  mesh_velocity_dof,
	      //				  mesh_l2g,
	      //				  mesh_trial_ref,
	      //				  xt,yt);
	      //get the physical integration weight
	      dV = fabs(jacDet)*dV_ref[k];
	      //get the trial function gradients
	      ck.gradTrialFromRef(&h_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,h_grad_trial);
	      ck.gradTrialFromRef(&vel_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,vel_grad_trial);
	      //get the solution 	
	      ck.valFromDOF(b_dof,&h_l2g[eN_nDOF_trial_element],&h_trial_ref[k*nDOF_trial_element],b);
	      ck.valFromDOF(h_dof,&h_l2g[eN_nDOF_trial_element],&h_trial_ref[k*nDOF_trial_element],h);
	      ck.valFromDOF(hu_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hu);
	      ck.valFromDOF(hv_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hv);
	      //get the solution gradients
	      ck.gradFromDOF(b_dof,&h_l2g[eN_nDOF_trial_element],h_grad_trial,grad_b);
	      ck.gradFromDOF(h_dof,&h_l2g[eN_nDOF_trial_element],h_grad_trial,grad_h);
	      ck.gradFromDOF(hu_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_hu);
	      ck.gradFromDOF(hv_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_hv);
	      //precalculate test function products with integration weights
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  h_test_dV[j] = h_test_ref[k*nDOF_trial_element+j]*dV;
		  vel_test_dV[j] = vel_test_ref[k*nDOF_trial_element+j]*dV;
		  for (int I=0;I<nSpace;I++)
		    {
		      h_grad_test_dV[j*nSpace+I]   = h_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
		      vel_grad_test_dV[j*nSpace+I] = vel_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin}
		    }
		}
	      evaluateCoefficientsForJacobian(g,
					      grad_b,
					      h,
					      hu,
					      hv,
					      mass_acc,
					      dmass_acc_h,
					      mom_hu_acc,
					      dmom_hu_acc_h,
					      dmom_hu_acc_hu,
					      mom_hv_acc,
					      dmom_hv_acc_h,
					      dmom_hv_acc_hv,
					      mass_adv,
					      dmass_adv_h,
					      dmass_adv_hu,
					      dmass_adv_hv,
					      mom_hu_adv,
					      dmom_hu_adv_h,
					      dmom_hu_adv_hu,
					      dmom_hu_adv_hv,
					      mom_hv_adv,
					      dmom_hv_adv_h,
					      dmom_hv_adv_hu,
					      dmom_hv_adv_hv,
					      mom_hu_source,
					      dmom_hu_source_h,
					      mom_hv_source,
					      dmom_hv_source_h);
	      //
	      //moving mesh (TODO)
	      //
	      //
	      //calculate time derivatives
	      //
	      ck.bdf(alphaBDF,
		     q_mass_acc_beta_bdf[eN_k],
		     mass_acc,
		     dmass_acc_h,
		     mass_acc_t,
		     dmass_acc_h_t);
	      ck.bdf(alphaBDF,
		     q_mom_hu_acc_beta_bdf[eN_k],
		     mom_hu_acc,
		     dmom_hu_acc_hu,
		     mom_hu_acc_t,
		     dmom_hu_acc_hu_t);
	      ck.bdf(alphaBDF,
		     q_mom_hv_acc_beta_bdf[eN_k],
		     mom_hv_acc,
		     dmom_hv_acc_hv,
		     mom_hv_acc_t,
		     dmom_hv_acc_hv_t);

	      for(int i=0;i<nDOF_test_element;i++)
		{
		  register int i_nSpace = i*nSpace;
		  for(int j=0;j<nDOF_trial_element;j++) 
		    { 
		      register int j_nSpace = j*nSpace;
		      //////////////////////
		      // h: h_h, h_u, h_v //
		      //////////////////////
		      if (IMPLICIT==1)
			{
			  elementJacobian_h_h[i][j] += 
			    dt*ck.MassJacobian_weak(dmass_acc_h_t,h_trial_ref[k*nDOF_trial_element+j],h_test_dV[i]) + 
			    dt*ck.NumericalDiffusionJacobian(q_numDiff_h_last[eN_k],&h_grad_trial[j_nSpace],&h_grad_test_dV[i_nSpace]);
			  
			  elementJacobian_h_hu[i][j] += 
			    dt*ck.AdvectionJacobian_weak(dmass_adv_hu,vel_trial_ref[k*nDOF_trial_element+j],&h_grad_test_dV[i_nSpace]);
			  
			  elementJacobian_h_hv[i][j] += 
			    dt*ck.AdvectionJacobian_weak(dmass_adv_hv,vel_trial_ref[k*nDOF_trial_element+j],&h_grad_test_dV[i_nSpace]);
			}
		      else //EXPLICIT
			{ // Mass matrix is NOT lumped
			  elementJacobian_h_h[i][j] += dt*ck.MassJacobian_weak(dmass_acc_h_t,h_trial_ref[k*nDOF_trial_element+j],h_test_dV[i]);
			  elementJacobian_h_hu[i][j] += 0;
			  elementJacobian_h_hv[i][j] += 0;
			}

		      //////////////////////
		      // u: u_h, u_u, u_v //
		      //////////////////////
		      if (IMPLICIT==1)
			{		      
			  elementJacobian_hu_h[i][j] += 
			    dt*ck.AdvectionJacobian_weak(dmom_hu_adv_h,h_trial_ref[k*nDOF_trial_element+j],&vel_grad_test_dV[i_nSpace]) +
			    dt*ck.ReactionJacobian_weak(dmom_hu_source_h,h_trial_ref[k*nDOF_trial_element+j],vel_test_dV[i]);
			  
			  elementJacobian_hu_hu[i][j] += 
			    dt*ck.MassJacobian_weak(dmom_hu_acc_hu_t,vel_trial_ref[k*nDOF_trial_element+j],vel_test_dV[i]) + 
			    dt*ck.AdvectionJacobian_weak(dmom_hu_adv_hu,vel_trial_ref[k*nDOF_trial_element+j],&vel_grad_test_dV[i_nSpace]) +
			    dt*ck.NumericalDiffusionJacobian(q_numDiff_hu_last[eN_k],&vel_grad_trial[j_nSpace],&vel_grad_test_dV[i_nSpace]);
			  
			  elementJacobian_hu_hv[i][j] += 
			    dt*ck.AdvectionJacobian_weak(dmom_hu_adv_hv,vel_trial_ref[k*nDOF_trial_element+j],&vel_grad_test_dV[i_nSpace]);
			}
		      else //EXPLICIT
			{ // Mass matrix is NOT lumped
			  elementJacobian_hu_h[i][j] += 0;			  
			  elementJacobian_hu_hu[i][j] += dt*ck.MassJacobian_weak(dmom_hu_acc_hu_t,vel_trial_ref[k*nDOF_trial_element+j],vel_test_dV[i]);
			  elementJacobian_hu_hv[i][j] += 0;
			}

		      //////////////////////
		      // v: v_h, v_u, v_v //
		      //////////////////////
		      if (IMPLICIT==1)
			{
			  elementJacobian_hv_h[i][j] += 
			    dt*ck.AdvectionJacobian_weak(dmom_hv_adv_h,h_trial_ref[k*nDOF_trial_element+j],&vel_grad_test_dV[i_nSpace]) +
			    dt*ck.ReactionJacobian_weak(dmom_hv_source_h,h_trial_ref[k*nDOF_trial_element+j],vel_test_dV[i]);
			  
			  elementJacobian_hv_hu[i][j] += 
			    dt*ck.AdvectionJacobian_weak(dmom_hv_adv_hu,vel_trial_ref[k*nDOF_trial_element+j],&vel_grad_test_dV[i_nSpace]);
			  
			  elementJacobian_hv_hv[i][j] += 
			    dt*ck.MassJacobian_weak(dmom_hv_acc_hv_t,vel_trial_ref[k*nDOF_trial_element+j],vel_test_dV[i]) + 
			    dt*ck.AdvectionJacobian_weak(dmom_hv_adv_hv,vel_trial_ref[k*nDOF_trial_element+j],&vel_grad_test_dV[i_nSpace]) + 
			    dt*ck.NumericalDiffusionJacobian(q_numDiff_hv_last[eN_k],&vel_grad_trial[j_nSpace],&vel_grad_test_dV[i_nSpace]);
			}
		      else //EXPLICIT
			{ // Mass matrix is NOT lumped
			  elementJacobian_hv_h[i][j] += 0;
			  elementJacobian_hv_hu[i][j] += 0;
			  elementJacobian_hv_hv[i][j] += dt*ck.MassJacobian_weak(dmom_hv_acc_hv_t,vel_trial_ref[k*nDOF_trial_element+j],vel_test_dV[i]);
			}
		      
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
		  globalJacobian[csrRowIndeces_h_hu[eN_i] + csrColumnOffsets_h_hu[eN_i_j]] += elementJacobian_h_hu[i][j];
		  globalJacobian[csrRowIndeces_h_hv[eN_i] + csrColumnOffsets_h_hv[eN_i_j]] += elementJacobian_h_hv[i][j];

		  globalJacobian[csrRowIndeces_hu_h[eN_i] + csrColumnOffsets_hu_h[eN_i_j]] += elementJacobian_hu_h[i][j];
		  globalJacobian[csrRowIndeces_hu_hu[eN_i] + csrColumnOffsets_hu_hu[eN_i_j]] += elementJacobian_hu_hu[i][j];
		  globalJacobian[csrRowIndeces_hu_hv[eN_i] + csrColumnOffsets_hu_hv[eN_i_j]] += elementJacobian_hu_hv[i][j];

		  globalJacobian[csrRowIndeces_hv_h[eN_i] + csrColumnOffsets_hv_h[eN_i_j]] += elementJacobian_hv_h[i][j];
		  globalJacobian[csrRowIndeces_hv_hu[eN_i] + csrColumnOffsets_hv_hu[eN_i_j]] += elementJacobian_hv_hu[i][j];
		  globalJacobian[csrRowIndeces_hv_hv[eN_i] + csrColumnOffsets_hv_hv[eN_i_j]] += elementJacobian_hv_hv[i][j];
		}//j
	    }//i
	}//elements
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
			     double* q_velocity_sge,
			     double* q_cfl,
			     double* q_numDiff_h_last,
			     double* q_numDiff_hu_last, 
			     double* q_numDiff_hv_last, 
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
	    elementJacobian_h_hu[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_h_hv[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_hu_h[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_hu_hu[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_hu_hv[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_hv_h[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_hv_hu[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_hv_hv[nDOF_test_element][nDOF_trial_element];
	  for (int i=0;i<nDOF_test_element;i++)
	    for (int j=0;j<nDOF_trial_element;j++)
	      {
		elementJacobian_h_h[i][j]=0.0;
		elementJacobian_h_hu[i][j]=0.0;
		elementJacobian_h_hv[i][j]=0.0;
		elementJacobian_hu_h[i][j]=0.0;
		elementJacobian_hu_hu[i][j]=0.0;
		elementJacobian_hu_hv[i][j]=0.0;
		elementJacobian_hv_h[i][j]=0.0;
		elementJacobian_hv_hu[i][j]=0.0;
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
		      //////////////////////
		      // h: h_h, h_u, h_v //
		      //////////////////////
		      // EXPLICIT AND LUMPED 
		      elementJacobian_h_h[i][j] += h_trial_ref[k*nDOF_trial_element+j]*h_test_dV[i];
		      elementJacobian_h_hu[i][j] += 0;
		      elementJacobian_h_hv[i][j] += 0;
		      
		      //////////////////////
		      // u: u_h, u_u, u_v //
		      //////////////////////
		      elementJacobian_hu_h[i][j] += 0;			  
		      elementJacobian_hu_hu[i][j] += vel_trial_ref[k*nDOF_trial_element+j]*vel_test_dV[i];
		      elementJacobian_hu_hv[i][j] += 0;
			
		      //////////////////////
		      // v: v_h, v_u, v_v //
		      //////////////////////
		      elementJacobian_hv_h[i][j] += 0;
		      elementJacobian_hv_hu[i][j] += 0;
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
		  globalJacobian[csrRowIndeces_h_h[eN_i] + csrColumnOffsets_h_h[eN_i_j]] += elementJacobian_h_h[i][j];
		  globalJacobian[csrRowIndeces_h_hu[eN_i] + csrColumnOffsets_h_hu[eN_i_j]] += elementJacobian_h_hu[i][j];
		  globalJacobian[csrRowIndeces_h_hv[eN_i] + csrColumnOffsets_h_hv[eN_i_j]] += elementJacobian_h_hv[i][j];

		  globalJacobian[csrRowIndeces_hu_h[eN_i] + csrColumnOffsets_hu_h[eN_i_j]] += elementJacobian_hu_h[i][j];
		  globalJacobian[csrRowIndeces_hu_hu[eN_i] + csrColumnOffsets_hu_hu[eN_i_j]] += elementJacobian_hu_hu[i][j];
		  globalJacobian[csrRowIndeces_hu_hv[eN_i] + csrColumnOffsets_hu_hv[eN_i_j]] += elementJacobian_hu_hv[i][j];

		  globalJacobian[csrRowIndeces_hv_h[eN_i] + csrColumnOffsets_hv_h[eN_i_j]] += elementJacobian_hv_h[i][j];
		  globalJacobian[csrRowIndeces_hv_hu[eN_i] + csrColumnOffsets_hv_hu[eN_i_j]] += elementJacobian_hv_hu[i][j];
		  globalJacobian[csrRowIndeces_hv_hv[eN_i] + csrColumnOffsets_hv_hv[eN_i_j]] += elementJacobian_hv_hv[i][j];
		}//j
	    }//i
	}//elements
      
    }

    ////
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
				   double* q_velocity_sge,
				   double* q_cfl,
				   double* q_numDiff_h_last,
				   double* q_numDiff_hu_last, 
				   double* q_numDiff_hv_last, 
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
	  register double  elementJacobian_h_h[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_h_hu[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_h_hv[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_hu_h[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_hu_hu[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_hu_hv[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_hv_h[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_hv_hu[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_hv_hv[nDOF_test_element][nDOF_trial_element];
	  for (int i=0;i<nDOF_test_element;i++)
	    for (int j=0;j<nDOF_trial_element;j++)
	      {
		elementJacobian_h_h[i][j]=0.0;
		elementJacobian_h_hu[i][j]=0.0;
		elementJacobian_h_hv[i][j]=0.0;
		elementJacobian_hu_h[i][j]=0.0;
		elementJacobian_hu_hu[i][j]=0.0;
		elementJacobian_hu_hv[i][j]=0.0;
		elementJacobian_hv_h[i][j]=0.0;
		elementJacobian_hv_hu[i][j]=0.0;
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
		      //////////////////////
		      // h: h_h, h_u, h_v //
		      //////////////////////
		      // EXPLICIT AND LUMPED 
		      elementJacobian_h_h[i][j] += (i==j ? 1.0 : 0.0)*h_test_dV[i];
		      elementJacobian_h_hu[i][j] += 0;
		      elementJacobian_h_hv[i][j] += 0;
		      
		      //////////////////////
		      // u: u_h, u_u, u_v //
		      //////////////////////
		      elementJacobian_hu_h[i][j] += 0;			  
		      elementJacobian_hu_hu[i][j] += (i==j ? 1.0 : 0.0)*vel_test_dV[i];
		      elementJacobian_hu_hv[i][j] += 0;
			
		      //////////////////////
		      // v: v_h, v_u, v_v //
		      //////////////////////
		      elementJacobian_hv_h[i][j] += 0;
		      elementJacobian_hv_hu[i][j] += 0;
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
		  globalJacobian[csrRowIndeces_h_hu[eN_i] + csrColumnOffsets_h_hu[eN_i_j]] += elementJacobian_h_hu[i][j];
		  globalJacobian[csrRowIndeces_h_hv[eN_i] + csrColumnOffsets_h_hv[eN_i_j]] += elementJacobian_h_hv[i][j];

		  globalJacobian[csrRowIndeces_hu_h[eN_i] + csrColumnOffsets_hu_h[eN_i_j]] += elementJacobian_hu_h[i][j];
		  globalJacobian[csrRowIndeces_hu_hu[eN_i] + csrColumnOffsets_hu_hu[eN_i_j]] += elementJacobian_hu_hu[i][j];
		  globalJacobian[csrRowIndeces_hu_hv[eN_i] + csrColumnOffsets_hu_hv[eN_i_j]] += elementJacobian_hu_hv[i][j];

		  globalJacobian[csrRowIndeces_hv_h[eN_i] + csrColumnOffsets_hv_h[eN_i_j]] += elementJacobian_hv_h[i][j];
		  globalJacobian[csrRowIndeces_hv_hu[eN_i] + csrColumnOffsets_hv_hu[eN_i_j]] += elementJacobian_hv_hu[i][j];
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

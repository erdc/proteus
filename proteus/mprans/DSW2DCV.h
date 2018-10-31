#ifndef DSW2DCV_H
#define DSW2DCV_H
#include <cmath>
#include <iostream>
#include <fstream>
#include "CompKernel.h"
#include "ModelFactory.h"
#include <assert.h>


//cek todo
//2. Get stabilization right
//3. Add Riemann solvers for external flux
//4. Add Riemann solvers for internal flux and DG terms
//5. Try other choices of variables h,hu,hv, Bova-Carey symmetrization?

#define lambda 1.0 // For Dispersive Model
#define WHICH_DISP_MODEL 1
// EJT.
// Recall this is the hyperbolic modification of the Green Naghdi
// equations. See Guermond, Popov, Tovar.

#define POWER_SMOOTHNESS_INDICATOR 2
#define VEL_FIX_POWER 2.
#define REESTIMATE_MAX_EDGE_BASED_CFL 1

// FOR CELL BASED ENTROPY VISCOSITY
#define ENTROPY(g,h,hu,hv,z,one_over_hReg) 0.5*(g*h*h + one_over_hReg*(hu*hu+hv*hv) + 2.*g*h*z)
#define DENTROPY_DH(g,h,hu,hv,z,one_over_hReg) g*h - 0.5*(hu*hu+hv*hv)*std::pow(one_over_hReg,2) + g*z
#define DENTROPY_DHU(g,h,hu,hv,z,one_over_hReg) hu*one_over_hReg
#define DENTROPY_DHV(g,h,hu,hv,z,one_over_hReg) hv*one_over_hReg

#define ENTROPY_FLUX1(g,h,hu,hv,z,one_over_hReg) (ENTROPY(g,h,hu,hv,z,one_over_hReg) + 0.5*g*h*h + g*h*z)*hu*one_over_hReg
#define ENTROPY_FLUX2(g,h,hu,hv,z,one_over_hReg) (ENTROPY(g,h,hu,hv,z,one_over_hReg) + 0.5*g*h*h + g*h*z)*hv*one_over_hReg

// FOR ESTIMATING MAX WAVE SPEEDS
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
  class DSW2DCV_base
  {
  public:
    virtual ~DSW2DCV_base(){}
      virtual void FCTStep(
                         double dt,
                            int NNZ, //number on non-zero entries on sparsity pattern
                            int numDOFs, //number of DOFs
                            double* lumped_mass_matrix, //lumped mass matrix (as vector)
                            double* h_old, //DOFs of solution at last stage
                            double* hu_old,
                            double* hv_old,
			    double* heta_old,
			    double* hw_old,
                            double* b_dof,
                            double* high_order_hnp1, //DOFs of high order solution at tnp1
                            double* high_order_hunp1,
                            double* high_order_hvnp1,
			    double* high_order_hetanp1,
			    double* high_order_hwnp1,
                            double* low_order_hnp1, //operators to construct low order solution
                            double* low_order_hunp1,
                            double* low_order_hvnp1,
			    double* low_order_hetanp1,
			    double* low_order_hwnp1,
                            double* limited_hnp1,
                            double* limited_hunp1,
                            double* limited_hvnp1,
			    double* limited_hetanp1,
			    double* limited_hwnp1,
                            int* csrRowIndeces_DofLoops, //csr row indeces
                            int* csrColumnOffsets_DofLoops, //csr column offsets
                            double* MassMatrix, //mass matrix
                            double* dH_minus_dL,
                            double* muH_minus_muL,
                            double hEps,
                            double* hReg,
                            int LUMPED_MASS_MATRIX
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
                                           double* edge_based_cfl
                                           )=0;
    virtual void calculateResidual_entropy_viscosity(// last EDGE BASED version
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
				   double* heta_dof_old,
				   double* hw_dof_old,
                                   double* b_dof,
                                   double* h_dof,
                                   double* hu_dof,
                                   double* hv_dof,
                                   double* heta_dof,
                        				   double* hw_dof,
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
				                           int offset_heta,
				                           int offset_hw,
                                   int stride_h,
                                   int stride_hu,
                                   int stride_hv,
                        				   int stride_heta,
                        				   int stride_hw,
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
                                   double* low_order_hnp1,
                                   double* low_order_hunp1,
                                   double* low_order_hvnp1,
				   double* low_order_hetanp1,
				   double* low_order_hwnp1,
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
                                     // h
                                                              int* csrRowIndeces_h_h,
                                                              int* csrColumnOffsets_h_h,
                                                              int* csrRowIndeces_h_hu,
                                                              int* csrColumnOffsets_h_hu,
                                                              int* csrRowIndeces_h_hv,
                                                              int* csrColumnOffsets_h_hv,
                                     int* csrRowIndeces_h_heta,
                                                              int* csrColumnOffsets_h_heta,
                                     int* csrRowIndeces_h_hw,
                                                              int* csrColumnOffsets_h_hw,
                                     // hu
                                                              int* csrRowIndeces_hu_h,
                                                              int* csrColumnOffsets_hu_h,
                                                              int* csrRowIndeces_hu_hu,
                                                              int* csrColumnOffsets_hu_hu,
                                                              int* csrRowIndeces_hu_hv,
                                                              int* csrColumnOffsets_hu_hv,
                                     int* csrRowIndeces_hu_heta,
                                                              int* csrColumnOffsets_hu_heta,
                                     int* csrRowIndeces_hu_hw,
                                                              int* csrColumnOffsets_hu_hw,
                                     // hv
                                                              int* csrRowIndeces_hv_h,
                                                              int* csrColumnOffsets_hv_h,
                                                              int* csrRowIndeces_hv_hu,
                                                              int* csrColumnOffsets_hv_hu,
                                                              int* csrRowIndeces_hv_hv,
                                                              int* csrColumnOffsets_hv_hv,
                                     int* csrRowIndeces_hv_heta,
                                                              int* csrColumnOffsets_hv_heta,
                                     int* csrRowIndeces_hv_hw,
                                                              int* csrColumnOffsets_hv_hw,
                                     // heta
                                     int* csrRowIndeces_heta_h,
                                                              int* csrColumnOffsets_heta_h,
                                                              int* csrRowIndeces_heta_hu,
                                                              int* csrColumnOffsets_heta_hu,
                                                              int* csrRowIndeces_heta_hv,
                                                              int* csrColumnOffsets_heta_hv,
                                     int* csrRowIndeces_heta_heta,
                                                              int* csrColumnOffsets_heta_heta,
                                     int* csrRowIndeces_heta_hw,
                                                              int* csrColumnOffsets_heta_hw,
                                     //hw
                                     int* csrRowIndeces_hw_h,
                                                              int* csrColumnOffsets_hw_h,
                                                              int* csrRowIndeces_hw_hu,
                                                              int* csrColumnOffsets_hw_hu,
                                                              int* csrRowIndeces_hw_hv,
                                                              int* csrColumnOffsets_hw_hv,
                                     int* csrRowIndeces_hw_heta,
                                                              int* csrColumnOffsets_hw_heta,
                                     int* csrRowIndeces_hw_hw,
                                                              int* csrColumnOffsets_hw_hw,
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
  class DSW2DCV : public DSW2DCV_base
  {
  public:
    const int nDOF_test_X_trial_element;
    CompKernelType ck;
    DSW2DCV():
      nDOF_test_X_trial_element(nDOF_test_element*nDOF_trial_element),
      ck()
    {
      std::cout<<"Constructing DSW2DCV<CompKernelTemplate<"
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
            /*  Ix += rx[i*3+k]*rx[j*3+k]; */
            /*  Iy += ry[i*3+k]*ry[j*3+k]; */
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
      /*        { */
      /*          for (int j=0;j<3;j++) */
      /*            { */
      /*              std::cout<<tau[i*3+j]<<'\t'; */
      /*            } */
      /*          std::cout<<std::endl; */
      /*        } */
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
        }
      else if (hR==0) // right dry state
        {
          lambda1 = velL-sqrt(g*hL);
          lambda3 = velL+2*sqrt(g*hL);
          if (debugging && false)
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

          if (debugging && false)
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
      if (debugging && false)
        std::cout << "lambda1, lambda3: " << lambda1 << ", " << lambda3 << std::endl;
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
      /*        { */
      /*          flux_mass += n[0]*f_mass[0]; */
      /*          velocity[0] = f_mass[0]; */
      /*          if (flowDirection >= 0.0) */
      /*            { */
      /*              flux_humom += n[0]*f_humom[0]; */
      /*              flux_hvmom += n[0]*f_hvmom[0]; */
      /*            } */
      /*        } */
      /* else */
      /*        { */
      /*          flux_mass += n[0]*bc_f_mass[0]; */
      /*          velocity[0] = bc_f_mass[0]; */
      /*          //cek still upwind the advection for Dirichlet? */
      /*          if (flowDirection >= 0.0) */
      /*            { */
      /*              flux_humom += n[0]*f_humom[0]; */
      /*              flux_hvmom += n[0]*f_hvmom[0]; */
      /*            } */
      /*          else */
      /*            { */
      /*              flux_humom+=n[0]*bc_f_humom[0]; */
      /*              flux_hvmom+=n[0]*bc_f_hvmom[0]; */
      /*            } */
      /*        } */
      /* if (isDOFBoundary_hv != 1) */
      /*        { */
      /*          flux_mass+=n[1]*f_mass[1]; */
      /*          velocity[1] = f_mass[1]; */
      /*          if (flowDirection >= 0.0) */
      /*            { */
      /*              flux_humom+=n[1]*f_humom[1]; */
      /*              flux_hvmom+=n[1]*f_hvmom[1]; */
      /*            } */
      /*        } */
      /* else */
      /*        { */
      /*          flux_mass+=n[1]*bc_f_mass[1]; */
      /*          velocity[1] = bc_f_mass[1]; */
      /*          //cek still upwind the advection for Dirichlet? */
      /*          if (flowDirection >= 0.0) */
      /*            { */
      /*              flux_humom+=n[1]*f_humom[1]; */
      /*              flux_hvmom+=n[1]*f_hvmom[1]; */
      /*            } */
      /*          else */
      /*            { */
      /*              flux_humom+=n[1]*bc_f_humom[1]; */
      /*              flux_hvmom+=n[1]*bc_f_hvmom[1]; */
      /*            } */
      /*        } */
      /* else */
      /*        { */
      /*          flux_mass +=n[2]*bc_f_mass[2]; */
      /*          velocity[2] = bc_f_mass[2]; */
      /*          //cek still upwind the advection for Dirichlet? */
      /*          if (flowDirection >= 0.0) */
      /*            { */
      /*              flux_humom+=n[2]*f_humom[2]; */
      /*              flux_hvmom+=n[2]*f_hvmom[2]; */
      /*            } */
      /*          else */
      /*            { */
      /*              flux_humom+=n[2]*bc_f_humom[2]; */
      /*              flux_hvmom+=n[2]*bc_f_hvmom[2]; */
      /*            } */
      /*        } */
      /* if (isDOFBoundary_h == 1) */
      /*        { */
      /*          flux_humom+= n[0]*(bc_h-p)*oneByRho; */
      /*          flux_hvmom+= n[1]*(bc_h-p)*oneByRho; */
      /*        } */
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
      /*        { */
      /*          dflux_mass_du += n[0]*df_mass_du[0]; */
      /*          if (flowDirection >= 0.0) */
      /*            { */
      /*              dflux_humom_du += n[0]*df_humom_du[0]; */
      /*              dflux_hvmom_du += n[0]*df_hvmom_du[0]; */
      /*              dflux_hvmom_dv += n[0]*df_hvmom_dv[0]; */
      /*            } */
      /*        } */
      /* else */
      /*        { */
      /*          //cek still upwind the advection for Dirichlet? */
      /*          if (flowDirection >= 0.0) */
      /*            { */
      /*              dflux_humom_du += n[0]*df_humom_du[0]; */
      /*              dflux_hvmom_du += n[0]*df_hvmom_du[0]; */
      /*              dflux_hvmom_dv += n[0]*df_hvmom_dv[0]; */
      /*            } */
      /*          else */
      /*            { */
      /*              if (isDOFBoundary_hv != 1) */
      /*                dflux_hvmom_dv += n[0]*df_hvmom_dv[0]; */
      /*            } */
      /*        } */
      /* if (isDOFBoundary_hv != 1) */
      /*        { */
      /*          dflux_mass_dv += n[1]*df_mass_dv[1]; */
      /*          if (flowDirection >= 0.0) */
      /*            { */
      /*              dflux_humom_du += n[1]*df_humom_du[1]; */
      /*              dflux_humom_dv += n[1]*df_humom_dv[1]; */
      /*              dflux_hvmom_dv += n[1]*df_hvmom_dv[1]; */
      /*            } */
      /*        } */
      /* else */
      /*        { */
      /*          //cek still upwind the advection for Dirichlet? */
      /*          if (flowDirection >= 0.0) */
      /*            { */
      /*              dflux_humom_du += n[1]*df_humom_du[1]; */
      /*              dflux_humom_dv += n[1]*df_humom_dv[1]; */
      /*              dflux_hvmom_dv += n[1]*df_hvmom_dv[1]; */
      /*            } */
      /*          else */
      /*            { */
      /*              if (isDOFBoundary_hu != 1) */
      /*                dflux_humom_du += n[1]*df_humom_du[1]; */
      /*            } */
      /*        } */
      /* else */
      /*        { */
      /*          //cek still upwind the advection for Dirichlet? */
      /*          if (flowDirection >= 0.0) */
      /*            { */
      /*              dflux_humom_du += n[2]*df_humom_du[2]; */
      /*              dflux_humom_dw += n[2]*df_humom_dw[2]; */
      /*              dflux_hvmom_dv += n[2]*df_hvmom_dv[2]; */
      /*            } */
      /*          else */
      /*            { */
      /*              if (isDOFBoundary_hu != 1) */
      /*                dflux_humom_du += n[2]*df_humom_du[2]; */
      /*              if (isDOFBoundary_hv != 1) */
      /*                dflux_hvmom_dv += n[2]*df_hvmom_dv[2]; */
      /*            } */
      /*        } */
      /* if (isDOFBoundary_h == 1) */
      /*        { */
      /*          dflux_humom_dp= -n[0]*oneByRho; */
      /*          dflux_hvmom_dp= -n[1]*oneByRho; */
      /*        } */
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
    /*                                  int* rowptr, */
    /*                                  int* colind, */
    /*                                  const int& isDOFBoundary, */
    /*                                  const int& isFluxBoundary, */
    /*                                  const double n[nSpace], */
    /*                                  double* bc_a, */
    /*                                  const double& bc_hu, */
    /*                                  const double& bc_flux, */
    /*                                  double* a, */
    /*                                  const double grad_phi[nSpace], */
    /*                                  const double& u, */
    /*                                  const double& penalty, */
    /*                                  double& flux) */
    /* { */
    /*   double diffusiveVelocityComponent_I,penaltyFlux,max_a; */
    /*   if(isDOFBoundary == 1) */
    /*  { */
    /*    flux = 0.0; */
    /*    max_a=0.0; */
    /*    for(int I=0;I<nSpace;I++) */
    /*      { */
    /*        diffusiveVelocityComponent_I=0.0; */
    /*        for(int m=rowptr[I];m<rowptr[I+1];m++) */
    /*          { */
    /*            diffusiveVelocityComponent_I -= a[m]*grad_phi[colind[m]]; */
    /*            max_a = fmax(max_a,a[m]); */
    /*          } */
    /*        flux+= diffusiveVelocityComponent_I*n[I]; */
    /*      } */
    /*    penaltyFlux = max_a*penalty*(u-bc_hu); */
    /*    flux += penaltyFlux; */
    /*  } */
    /*   else if(isFluxBoundary == 1) */
    /*  { */
    /*    flux = bc_flux; */
    /*  } */
    /*   else */
    /*  { */
    /*    std::cerr<<"warning, diffusion term with no boundary condition set, setting diffusive flux to 0.0"<<std::endl; */
    /*    flux = 0.0; */
    /*  } */
    /* } */

    /* inline */
    /* double ExteriorNumericalDiffusiveFluxJacobian(const double& eps, */
    /*                                            int* rowptr, */
    /*                                            int* colind, */
    /*                                            const int& isDOFBoundary, */
    /*                                            const double n[nSpace], */
    /*                                            double* a, */
    /*                                            const double& v, */
    /*                                            const double grad_hv[nSpace], */
    /*                                            const double& penalty) */
    /* { */
    /*   double dvel_I,tmp=0.0,max_a=0.0; */
    /*   if(isDOFBoundary >= 1) */
    /*  { */
    /*    for(int I=0;I<nSpace;I++) */
    /*      { */
    /*        dvel_I=0.0; */
    /*        for(int m=rowptr[I];m<rowptr[I+1];m++) */
    /*          { */
    /*            dvel_I -= a[m]*grad_hv[colind[m]]; */
    /*            max_a = fmax(max_a,a[m]); */
    /*          } */
    /*        tmp += dvel_I*n[I]; */
    /*      } */
    /*    tmp +=max_a*penalty*v; */
    /*  } */
    /*   return tmp; */
    /* } */

    void FCTStep(double dt,
                    int NNZ, //number on non-zero entries on sparsity pattern
                    int numDOFs, //number of DOFs
                    double* lumped_mass_matrix, //lumped mass matrix (as vector))
                    double* h_old, //DOFs of solution at last stage
                    double* hu_old,
                    double* hv_old,
		 double* heta_old,
		 double* hw_old,
                    double* b_dof,
                    double* high_order_hnp1, //DOFs of high order solution at tnp1
                    double* high_order_hunp1,
                    double* high_order_hvnp1,
		 double* high_order_hetanp1,
		 double* high_order_hwnp1,
                    double* low_order_hnp1, //operators to construct low order solution
                    double* low_order_hunp1,
                    double* low_order_hvnp1,
		 double* low_order_hetanp1,
		 double* low_order_hwnp1,
                    double* limited_hnp1,
                    double* limited_hunp1,
                    double* limited_hvnp1,
		 double* limited_hetanp1,
		 double* limited_hwnp1,
                    int* csrRowIndeces_DofLoops, //csr row indeces
                    int* csrColumnOffsets_DofLoops, //csr column offsets
                    double* MassMatrix, //mass matrix
                    double* dH_minus_dL,
                    double* muH_minus_muL,
                    double hEps,
                    double* hReg,
                    int LUMPED_MASS_MATRIX
                    )
    {
      register double Rneg[numDOFs];
      //////////////////
      // LOOP in DOFs //
      //////////////////
      int ij=0;
      for (int i=0; i<numDOFs; i++)
        {
          //read some vectors
          double high_order_hnp1i  = high_order_hnp1[i];
          double hni = h_old[i];
          double Zi = b_dof[i];
          double mi = lumped_mass_matrix[i];

          double minH=0.;
          double Pnegi=0.;
          // LOOP OVER THE SPARSITY PATTERN (j-LOOP)//
          for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
            {
              int j = csrColumnOffsets_DofLoops[offset];
              // read some vectors
              double hnj = h_old[j];
              double Zj = b_dof[j];

              // COMPUTE STAR SOLUTION // hStar, huStar and hvStar
              double hStarij  = fmax(0., hni + Zi - fmax(Zi,Zj));
              double hStarji  = fmax(0., hnj + Zj - fmax(Zi,Zj));

              // i-th row of flux correction matrix
              double ML_minus_MC = (LUMPED_MASS_MATRIX == 1 ? 0. : (i==j ? 1. : 0.)*mi - MassMatrix[ij]);
              double FluxCorrectionMatrix1 =
                ML_minus_MC*(high_order_hnp1[j]-hnj - (high_order_hnp1i-hni))
                + dt*(dH_minus_dL[ij]-muH_minus_muL[ij])*(hStarji-hStarij)
                + dt*muH_minus_muL[ij]*(hnj-hni);

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
          if (high_order_hnp1[i] <= hReg[i]) //hEps
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
	  double high_order_hetanp1i = high_order_hetanp1[i];
	  double high_order_hwnp1i = high_order_hwnp1[i];
          double hni = h_old[i];
          double huni = hu_old[i];
          double hvni = hv_old[i];
	  double hetani = heta_old[i];
	  double hwni = hw_old[i];
          double Zi = b_dof[i];
          double mi = lumped_mass_matrix[i];
          double one_over_hiReg = 2*hni/(hni*hni+std::pow(fmax(hni,hEps),2)); //hEps

          double ith_Limiter_times_FluxCorrectionMatrix1 = 0.;
          double ith_Limiter_times_FluxCorrectionMatrix2 = 0.;
          double ith_Limiter_times_FluxCorrectionMatrix3 = 0.;
	  double ith_Limiter_times_FluxCorrectionMatrix4 = 0.;
	  double ith_Limiter_times_FluxCorrectionMatrix5 = 0.;
          double Rnegi = Rneg[i];
          // LOOP OVER THE SPARSITY PATTERN (j-LOOP)//
          for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
            {
              int j = csrColumnOffsets_DofLoops[offset];
              // read some vectors
              double hnj = h_old[j];
              double hunj = hu_old[j];
              double hvnj = hv_old[j];
	      double hetanj = heta_old[j];
	      double hwnj = hw_old[j];
              double Zj = b_dof[j];
              double one_over_hjReg = 2*hnj/(hnj*hnj+std::pow(fmax(hnj,hEps),2)); //hEps

              // COMPUTE STAR SOLUTION // hStar, huStar and hvStar
              double hStarij  = fmax(0., hni + Zi - fmax(Zi,Zj));
              double huStarij = huni*hStarij*one_over_hiReg;
              double hvStarij = hvni*hStarij*one_over_hiReg;
	      double hetaStarij = hetani*hStarij*one_over_hiReg;
	      double hwStarij = hwni*hStarij*one_over_hiReg;

              double hStarji  = fmax(0., hnj + Zj - fmax(Zi,Zj));
              double huStarji = hunj*hStarji*one_over_hjReg;
              double hvStarji = hvnj*hStarji*one_over_hjReg;
	      double hetaStarji = hetanj*hStarji*one_over_hjReg;
	      double hwStarji = hwnj*hStarji*one_over_hjReg;

              // COMPUTE FLUX CORRECTION MATRICES
              double ML_minus_MC = (LUMPED_MASS_MATRIX == 1 ? 0. : (i==j ? 1. : 0.)*mi - MassMatrix[ij]);
              double FluxCorrectionMatrix1 =
                ML_minus_MC*(high_order_hnp1[j]-hnj - (high_order_hnp1i-hni))
                + dt*(dH_minus_dL[ij]-muH_minus_muL[ij])*(hStarji-hStarij)
                + dt*muH_minus_muL[ij]*(hnj-hni);

              double FluxCorrectionMatrix2 =
                ML_minus_MC*(high_order_hunp1[j]-hunj - (high_order_hunp1i-huni))
                + dt*(dH_minus_dL[ij]-muH_minus_muL[ij])*(huStarji-huStarij)
                + dt*muH_minus_muL[ij]*(hunj-huni);

              double FluxCorrectionMatrix3 =
                ML_minus_MC*(high_order_hvnp1[j]-hvnj - (high_order_hvnp1i-hvni))
                + dt*(dH_minus_dL[ij]-muH_minus_muL[ij])*(hvStarji-hvStarij)
                + dt*muH_minus_muL[ij]*(hvnj-hvni);

	       double FluxCorrectionMatrix4 =
		 ML_minus_MC*(high_order_hetanp1[j]-hetanj - (high_order_hetanp1i-hetani))
		 + dt*(dH_minus_dL[ij]-muH_minus_muL[ij])*(hetaStarji-hetaStarij)
		  + dt*muH_minus_muL[ij]*(hetanj-hetani);

	        double FluxCorrectionMatrix5 =
		 ML_minus_MC*(high_order_hwnp1[j]-hwnj - (high_order_hwnp1i-hwni))
		 + dt*(dH_minus_dL[ij]-muH_minus_muL[ij])*(hwStarji-hwStarij)
		  + dt*muH_minus_muL[ij]*(hwnj-hwni);

              // compute limiter based on water height
              double Lij = (FluxCorrectionMatrix1 > 0. ? std::min(1.,Rneg[j]) : std::min(Rnegi,1.));

              ith_Limiter_times_FluxCorrectionMatrix1 += Lij*FluxCorrectionMatrix1;
              ith_Limiter_times_FluxCorrectionMatrix2 += Lij*FluxCorrectionMatrix2;
              ith_Limiter_times_FluxCorrectionMatrix3 += Lij*FluxCorrectionMatrix3;
	      ith_Limiter_times_FluxCorrectionMatrix4 += Lij*FluxCorrectionMatrix4;
	      ith_Limiter_times_FluxCorrectionMatrix5 += Lij*FluxCorrectionMatrix5;
              //update ij
              ij+=1;
            }

          double one_over_mi = 1.0/lumped_mass_matrix[i];
          limited_hnp1[i]  = low_order_hnp1[i]  + one_over_mi*ith_Limiter_times_FluxCorrectionMatrix1;
          limited_hunp1[i] = low_order_hunp1[i] + one_over_mi*ith_Limiter_times_FluxCorrectionMatrix2;
          limited_hvnp1[i] = low_order_hvnp1[i] + one_over_mi*ith_Limiter_times_FluxCorrectionMatrix3;
	  limited_hetanp1[i] = low_order_hetanp1[i] + one_over_mi*ith_Limiter_times_FluxCorrectionMatrix4;
	  limited_hwnp1[i] = low_order_hwnp1[i] + one_over_mi*ith_Limiter_times_FluxCorrectionMatrix5;
          if (limited_hnp1[i] < -1E-14 && dt < 1.0)
            {
              std::cout << "Limited water height is negative: "
                        <<  limited_hnp1[i]
                        << " ... aborting!" << std::endl;
              abort();
            }
          else
            {
              limited_hnp1[i] = fmax(limited_hnp1[i],0.);
              //double aux = fmax(limited_hnp1[i],hEps); // hEps
              double aux = fmax(limited_hnp1[i],hReg[i]); // hEps. Using hReg makes the code more robust
              limited_hunp1[i] *= 2*std::pow(limited_hnp1[i],VEL_FIX_POWER)/(std::pow(limited_hnp1[i],VEL_FIX_POWER)+std::pow(aux,VEL_FIX_POWER));
              limited_hvnp1[i] *= 2*std::pow(limited_hnp1[i],VEL_FIX_POWER)/(std::pow(limited_hnp1[i],VEL_FIX_POWER)+std::pow(aux,VEL_FIX_POWER));
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
                                 double* edge_based_cfl)
    {
      register double psi[numDOFsPerEqn];
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
                                                                hEps,hEps,false)*cij_norm, //hEps
                                  maxWaveSpeedSharpInitialGuess(g,nxji,nyji,
                                                                hj,huj,hvj,
                                                                hi,hui,hvi,
                                                                hEps,hEps,false)*cji_norm); //hEps
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
            psi[i] = std::pow(alphai,POWER_SMOOTHNESS_INDICATOR); //NOTE: they use alpha^2 in the paper
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
                      dLij = dLow[ij]*fmax(psi[i],psi[j]); // enhance the order to 2nd order. No EV

                      muLowij = fmax(fmax(0.,-(ui*Cx[ij] + vi*Cy[ij])),fmax(0,(uj*Cx[ij] + vj*Cy[ij])));
                      muLij = muLowij*fmax(psi[i],psi[j]); // enhance the order to 2nd order. No EV
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

    void calculateResidual_entropy_viscosity(// last EDGE BASED version
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
			   double* heta_dof_old,
			   double* hw_dof_old,
                           double* b_dof,
                           double* h_dof,
                           double* hu_dof,
                           double* hv_dof,
                           double* heta_dof,
                           double* hw_dof,
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
			   int offset_heta,
			   int offset_hw,
                           int stride_h,
                           int stride_hu,
                           int stride_hv,
			   int stride_heta,
			   int stride_hw,
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
                           double* low_order_hnp1,
                           double* low_order_hunp1,
                           double* low_order_hvnp1,
			   double* low_order_hetanp1,
			   double* low_order_hwnp1,
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
            elementResidual_hv[nDOF_test_element],
            elementResidual_heta[nDOF_test_element],
      	    elementResidual_hw[nDOF_test_element];

          for (int i=0;i<nDOF_test_element;i++)
            {
              elementResidual_h[i]=0.0;
              elementResidual_hu[i]=0.0;
              elementResidual_hv[i]=0.0;
              elementResidual_heta[i]=0.0;
      	      elementResidual_hw[i]=0.0;
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
                h=0.0,hu=0.0,hv=0.0,heta=0.0,hw=0.0, // solution at current time
                h_old=0.0,hu_old=0.0,hv_old=0.0,heta_old=0.0,hw_old=0.0, // solution at lstage
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
              ck.valFromDOF(heta_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],heta);
              ck.valFromDOF(hw_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hw);
              // get the solution at the lstage
              ck.valFromDOF(h_dof_old,&h_l2g[eN_nDOF_trial_element],&h_trial_ref[k*nDOF_trial_element],h_old);
              ck.valFromDOF(hu_dof_old,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hu_old);
              ck.valFromDOF(hv_dof_old,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hv_old);
              ck.valFromDOF(heta_dof_old,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],heta_old);
      	      ck.valFromDOF(hw_dof_old,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hw_old);
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
              //save velocity at quadrature points for other models to use
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
                  elementResidual_heta[i] += (heta - heta_old)*h_test_dV[i];
            		  elementResidual_hw[i] += (hw - hw_old)*h_test_dV[i];
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
              globalResidual[offset_heta+stride_heta*vel_gi] += elementResidual_heta[i];
      	      globalResidual[offset_hw+stride_hw*vel_gi] += elementResidual_hw[i];
            }
        }
      // ********** END OF CELL LOOPS ********** //

      if (SECOND_CALL_CALCULATE_RESIDUAL==0) // This is to same some time
        {
          ///////////////////////////////////////////
          // ********** COMPUTE ENTROPY ********** //
          ///////////////////////////////////////////
          // compute entropy (defined as eta) corresponding to ith node
          register double eta[numDOFsPerEqn];
          for (int i=0; i<numDOFsPerEqn; i++)
            {
              // COMPUTE ENTROPY. NOTE: WE CONSIDER A FLAT BOTTOM
              double hni = h_dof_old[i];
              double one_over_hniReg = 2*hni/(hni*hni+std::pow(fmax(hni,hEps),2)); //hEps
              eta[i] = ENTROPY(g,hni,hu_dof_old[i],hv_dof_old[i],0.,one_over_hniReg);
            }
          // ********** END OF COMPUTING ENTROPY ********** //

          ////////////////////////////////////////////////////////////////////////////////////////
          // ********** COMPUTE SMOOTHNESS INDICATOR, GLOBAL ENTROPY RESIDUAL and dL ********** //
          ////////////////////////////////////////////////////////////////////////////////////////
          // Smoothness indicator is based on the solution. psi_i = psi_i(alpha_i);
          // alpha_i = |sum(uj-ui)|/sum|uj-ui|
          int ij = 0;
          register double global_entropy_residual[numDOFsPerEqn];
          register double psi[numDOFsPerEqn], etaMax[numDOFsPerEqn], etaMin[numDOFsPerEqn];
          for (int i=0; i<numDOFsPerEqn; i++)
            {
              double hi = h_dof_old[i]; // solution at time tn for the ith DOF
              double hui = hu_dof_old[i];
              double hvi = hv_dof_old[i];
              double one_over_hiReg = 2*hi/(hi*hi+std::pow(fmax(hi,hEps),2)); //hEps

              // For eta min and max
              etaMax[i] = fabs(eta[i]);
              etaMin[i] = fabs(eta[i]);

              // FOR ENTROPY RESIDUAL
              double ith_flux_term1=0., ith_flux_term2=0., ith_flux_term3=0.;
              double entropy_flux=0.;
              double eta_prime1 = DENTROPY_DH(g,hi,hui,hvi,0.,one_over_hiReg); // NOTE: WE CONSIDER A FLAT BOTTOM
              double eta_prime2 = DENTROPY_DHU(g,hi,hui,hvi,0.,one_over_hiReg);
              double eta_prime3 = DENTROPY_DHV(g,hi,hui,hvi,0.,one_over_hiReg);

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
                  // NOTE: WE CONSIDER FLAT BOTTOM. i.e., Zj=0
                  double Zj = 0.;

                  // FOR ENTROPY RESIDUAL //
                  ith_flux_term1 += huj*Cx[ij] + hvj*Cy[ij]; // f1*C
                  ith_flux_term2 += uj*huj*Cx[ij] + uj*hvj*Cy[ij] + g*hi*(hj+Zj)*Cx[ij];
                  ith_flux_term3 += vj*huj*Cx[ij] + vj*hvj*Cy[ij] + g*hi*(hj+Zj)*Cy[ij];

                  // NOTE: WE CONSIDER FLAT BOTTOM
                  entropy_flux += ( Cx[ij]*ENTROPY_FLUX1(g,hj,huj,hvj,0.,one_over_hjReg) +
                                    Cy[ij]*ENTROPY_FLUX2(g,hj,huj,hvj,0.,one_over_hjReg) );

                  /////////////////////////////////
                  // COMPUTE ETA MIN AND ETA MAX //
                  /////////////////////////////////
                  etaMax[i] = fmax(etaMax[i],fabs(eta[j]));
                  etaMin[i] = fmin(etaMin[i],fabs(eta[j]));

                  // FOR SMOOTHNESS INDICATOR //
                  alpha_numerator += hj - hi;
                  alpha_denominator += fabs(hj - hi);
                  //update ij
                  ij+=1;
                }
              /////////////////////////////////////
              // COMPUTE GLOBAL ENTROPY RESIDUAL //
              /////////////////////////////////////
              double one_over_entNormFactori = 2./(etaMax[i]-etaMin[i]+1E-15);
              global_entropy_residual[i] = one_over_entNormFactori*
                fabs(entropy_flux -(ith_flux_term1*eta_prime1 + ith_flux_term2*eta_prime2 + ith_flux_term3*eta_prime3));

              //////////////////////////////////
              // COMPUTE SMOOTHNESS INDICATOR //
              //////////////////////////////////
              if (hi <= hReg[i]) //hEps, hReg makes the method more robust
                {
                  alphai = 1.;
                  global_entropy_residual[i] = 1E10;
                }
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
                psi[i] = std::pow(alphai,POWER_SMOOTHNESS_INDICATOR); //NOTE: they use alpha^2 in the paper

            }
          // ********** END OF COMPUTING SMOOTHNESS INDICATOR, and GLOBAL ENTROPY RESIDUAL ********** //

          //double max_ui = 0.; //TMP
          //double max_vi = 0.;
          //double max_hui = 0.; //TMP
          //double max_hvi = 0.;
          //double max_hxui = 0.;
          //double max_hxvi = 0.;

          ////////////////////////////////////////
          // ********** Loop on DOFs ********** // to compute flux and dissipative terms
          ////////////////////////////////////////
          ij = 0;

          for (int i=0; i<numDOFsPerEqn; i++)
            {
              double hi = h_dof_old[i];
              double hui = hu_dof_old[i];
              double hvi = hv_dof_old[i];
	            double hetai = heta_dof_old[i];
	            double hwi = hw_dof_old[i];

              double Zi = b_dof[i];
              double mi = lumped_mass_matrix[i];
	            double cut_offi = 1.0/std::pow(10.0,8.0);
              double one_over_hiReg = 2*hi/(hi*hi+std::pow(fmax(hi,cut_offi),2)); // hEps
              double ui = hui*one_over_hiReg;
              double vi = hvi*one_over_hiReg;

              double ith_flux_term1=0., ith_flux_term2=0., ith_flux_term3=0., ith_flux_term4=0., ith_flux_term5=0.;
              // LOW ORDER DISSIPATIVE TERMS
              double
                ith_dLij_minus_muLij_times_hStarStates=0.,
                ith_dLij_minus_muLij_times_huStarStates=0.,
                ith_dLij_minus_muLij_times_hvStarStates=0.,
		            ith_dLij_minus_muLij_times_hetaStarStates=0.,
		            ith_dLij_minus_muLij_times_hwStarStates=0.,
                ith_muLij_times_hStates=0.,
                ith_muLij_times_huStates=0.,
                ith_muLij_times_hvStates=0.,
		            ith_muLij_times_hetaStates=0.,
	              ith_muLij_times_hwStates=0.;

              // HIGH ORDER DISSIPATIVE TERMS
              double
                ith_dHij_minus_muHij_times_hStarStates=0.,
                ith_dHij_minus_muHij_times_huStarStates=0.,
                ith_dHij_minus_muHij_times_hvStarStates=0.,
		            ith_dHij_minus_muHij_times_hetaStarStates=0.,
	            	ith_dHij_minus_muHij_times_hwStarStates=0.,
                ith_muHij_times_hStates=0.,
                ith_muHij_times_huStates=0.,
                ith_muHij_times_hvStates=0.,
		            ith_muHij_times_hetaStates=0.,
		            ith_muHij_times_hwStates=0.;

              ///////////////////
              // FRICTION TERM //
              ///////////////////
              double veli_norm = std::sqrt(ui*ui+vi*vi);
              double hi_to_the_gamma = std::pow(hi,gamma);
              double friction_aux =
                veli_norm == 0. ? 0.: (2*g*n2*veli_norm*mi/
                                       (hi_to_the_gamma+fmax(hi_to_the_gamma,xi*g*n2*dt*veli_norm)));
              double ith_friction_term2 = friction_aux*hui;
              double ith_friction_term3 = friction_aux*hvi;

              if (LINEAR_FRICTION==1)
                {
                  ith_friction_term2 = mannings*hui*mi;
                  ith_friction_term3 = mannings*hvi*mi;
                }

        // Define things here to make life easier.
        double meshi = std::sqrt(mi);
        double alphai = lambda / (3.0 * meshi);
        double constanti = lambda * g / meshi;
        double etai = hetai*one_over_hiReg;
        double eta_over_h_i = etai * one_over_hiReg;
        double ratio_i =  2.0 * hetai/(std::pow(etai,2.0)+std::pow(hi,2.0)+cut_offi);
        double x0_i =  fmin(2.0, std::sqrt(1.0+1.0/(2*alphai*(fmax(etai,0.0)+cut_offi))));

        /////////////////////
	      // For Force Terms //
        /////////////////////
        double hSqd_GammaPi = 4.0*(x0_i-1.0)*hetai+(1-std::pow(x0_i,2.0))*std::pow(hi,2.0);
        if (etai < 0.0)
        {
          hSqd_GammaPi = 0.0;
        }
        else if (etai <= x0_i*hi)
        {
          hSqd_GammaPi = 3.0*std::pow(etai,2.0)+std::pow(hi,2.0)-4.0*hetai;
        }
	      double force_term_hetai = hwi*mi*std::pow(ratio_i,3.0);
	      double force_term_hwi = -constanti*hSqd_GammaPi*std::pow(ratio_i,3.0)*mi;



              // loop over the sparsity pattern of the i-th DOF
              for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
                {
                  int j = csrColumnOffsets_DofLoops[offset];
                  double hj = h_dof_old[j];
                  double huj = hu_dof_old[j];
                  double hvj = hv_dof_old[j];
            	    double hetaj = heta_dof_old[j];
            	    double hwj = hw_dof_old[j];
                  double Zj = b_dof[j];
		              double cut_off = 1.0 / std::pow(10.0,8.0);
                  double one_over_hjReg = 2*hj/(hj*hj+std::pow(fmax(hj,cut_off),2)); //hEps
                  double uj = huj*one_over_hjReg;
		              double vj = hvj*one_over_hjReg;
            		  // EJT. Corrected terms to match our paper.
                  double mj = lumped_mass_matrix[j];
                  double meshj = std::sqrt(mj); // local mesh size in 2d
                  double alphaj = lambda / (3.0 * meshj);
                  double etaj = hetaj*one_over_hjReg;
                  double eta_over_h_j = etaj * one_over_hjReg;
                  double ratio_j =  2.0 * hetaj/(std::pow(etaj,2.0)+std::pow(hj,2.0)+cut_off);
                  double x0_j =  fmin(2.0, std::sqrt(1.0+1.0/(2*alphaj*(fmax(etaj,0.0)+cut_off))));

                  /////////////////////
          	      // For Flux Terms //
                  /////////////////////
                  double pTildej = -alphaj*g*(std::pow(x0_j,2.0)-1.0)*hetaj*hj;
                  if (etaj < 0.0)
                  {
                    pTildej = 0.0;
                  }
                  else if (etaj <= x0_j*hj)
                  {
                    pTildej = -alphaj*g*etaj*(std::pow(etaj,2.0)-std::pow(hj,2.0));
                  }




                  // Nodal projection of fluxes
                  ith_flux_term1 += huj*Cx[ij] + hvj*Cy[ij]; // f1*C
	                ith_flux_term2 += uj*huj*Cx[ij] + uj*hvj*Cy[ij] + (g*hi*(hj+Zj) + pTildej )*Cx[ij];
                  ith_flux_term3 += vj*huj*Cx[ij] + vj*hvj*Cy[ij] + (g*hi*(hj+Zj) + pTildej )*Cy[ij];
		              ith_flux_term4 += hetaj*uj*Cx[ij] + hetaj*vj*Cy[ij];
	                ith_flux_term5 += hwj*uj*Cx[ij] + hwj*vj*Cy[ij];

                  // COMPUTE STAR SOLUTION // hStar, huStar and hvStar
                  double hStarij  = fmax(0., hi + Zi - fmax(Zi,Zj));
                  double huStarij = hui*hStarij*one_over_hiReg;
                  double hvStarij = hvi*hStarij*one_over_hiReg;
            	    double hetaStarij = hetai*hStarij*one_over_hiReg;
            	    double hwStarij = hwi*hStarij*one_over_hiReg;

                  double hStarji  = fmax(0., hj + Zj - fmax(Zi,Zj));
                  double huStarji = huj*hStarji*one_over_hjReg;
                  double hvStarji = hvj*hStarji*one_over_hjReg;
            	    double hetaStarji = hetaj*hStarji*one_over_hjReg;
            	    double hwStarji = hwj*hStarji*one_over_hjReg;

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
                                                                      hEps,hEps,false)*cij_norm, //hEps
                                        maxWaveSpeedSharpInitialGuess(g,nxji,nyji,
                                                                      hj,huj,hvj,
                                                                      hi,hui,hvi,
                                                                      hEps,hEps,false)*cji_norm); //hEps
                        }
                      dLij = dLowij*fmax(psi[i],psi[j]); // enhance the order to 2nd order. No EV

                      ///////////////////////////////////////
                      // WELL BALANCING DISSIPATIVE MATRIX //
                      ///////////////////////////////////////
                      muLowij = fmax(fmax(0.,-(ui*Cx[ij] + vi*Cy[ij])),fmax(0,(uj*Cx[ij] + vj*Cy[ij])));
                      muLij = muLowij*fmax(psi[i],psi[j]); // enhance the order to 2nd order. No EV

                      ///////////////////////
                      // ENTROPY VISCOSITY //
                      ///////////////////////
                      double dEVij = cE*fmax(global_entropy_residual[i],global_entropy_residual[j]);

                      if (cE < 1000) // Hack to quickly deactivate EV
                        {
                          dHij  = fmin(dLowij,dEVij);
                          muHij = fmin(muLowij,dEVij);
                        }
                      else
                        {
                          dHij  = dLij;
                          muHij = muLij;
                        }

                      // compute dij_minus_muij times star solution terms
                      ith_dHij_minus_muHij_times_hStarStates  += (dHij - muHij)*(hStarji-hStarij);
                      ith_dHij_minus_muHij_times_huStarStates += (dHij - muHij)*(huStarji-huStarij);
                      ith_dHij_minus_muHij_times_hvStarStates += (dHij - muHij)*(hvStarji-hvStarij);
            	      ith_dHij_minus_muHij_times_hetaStarStates += (dHij - muHij)*(hetaStarji-hetaStarij);
            	      ith_dHij_minus_muHij_times_hwStarStates += (dHij - muHij)*(hwStarji-hwStarij);

                      ith_dLij_minus_muLij_times_hStarStates  += (dLij - muLij)*(hStarji-hStarij);
                      ith_dLij_minus_muLij_times_huStarStates += (dLij - muLij)*(huStarji-huStarij);
                      ith_dLij_minus_muLij_times_hvStarStates += (dLij - muLij)*(hvStarji-hvStarij);
            	      ith_dLij_minus_muLij_times_hetaStarStates += (dLij - muLij)*(hetaStarji-hetaStarij);
            	      ith_dLij_minus_muLij_times_hwStarStates += (dLij - muLij)*(hwStarji-hwStarij);

                      // compute muij times solution terms
                      ith_muHij_times_hStates  += muHij*(hj-hi);
                      ith_muHij_times_huStates += muHij*(huj-hui);
                      ith_muHij_times_hvStates += muHij*(hvj-hvi);
            	      ith_muHij_times_hetaStates += muHij*(hetaj-hetai);
            	      ith_muHij_times_hwStates += muHij*(hwj-hwi);

                      ith_muLij_times_hStates  += muLij*(hj-hi);
		      ith_muLij_times_huStates += muLij*(huj-hui);
                      ith_muLij_times_hvStates += muLij*(hvj-hvi);
            	      ith_muLij_times_hetaStates += muLij*(hetaj-hetai);
            	      ith_muLij_times_hwStates += muLij*(hwj-hwi);

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

              ////////////////////////
              // LOW ORDER SOLUTION //: lumped mass matrix and low order dissipative matrix
              ////////////////////////


      	      low_order_hnp1[i]  = hi  - dt/mi*(ith_flux_term1
          						- ith_dLij_minus_muLij_times_hStarStates
          						- ith_muLij_times_hStates);

              low_order_hunp1[i] = hui - dt/mi*(ith_flux_term2
          						- ith_dLij_minus_muLij_times_huStarStates
						- ith_muLij_times_huStates//);
						+ ith_friction_term2);
	      // EJT. Set v = 0 for 1d setting.
              low_order_hvnp1[i] = 0.0*( hvi - dt/mi*(ith_flux_term3
						          - ith_dLij_minus_muLij_times_hvStarStates
						          - ith_muLij_times_hvStates));

	      low_order_hetanp1[i] = hetai - dt/mi*(ith_flux_term4
      		             			          - ith_dLij_minus_muLij_times_hetaStarStates
      						          - ith_muLij_times_hetaStates) + dt/mi*force_term_hetai;

              low_order_hwnp1[i] = hwi - dt/mi*(ith_flux_term5
          						- ith_dLij_minus_muLij_times_hwStarStates
          						- ith_muLij_times_hwStates) + dt/mi*force_term_hwi;
              // FIX LOW ORDER SOLUTION //
              if (low_order_hnp1[i] < -1E-12 && dt < 1.0)
                {
                  std::cout << "dt taken: " << dt
                            << std::endl;
                  std::cout << "********.... "
                            << "Low order water height is negative: "
                            << hi << "\t"
                            << low_order_hnp1[i]
                            << " ... aborting!" << std::endl;
                  abort();
                }
              else
                {
                  low_order_hnp1[i] = fmax(low_order_hnp1[i],0.);
                  //double aux = fmax(low_order_hnp1[i],hEps); //hEps
                  double aux = fmax(low_order_hnp1[i],hReg[i]); //hEps. Using hReg makes the code more robust
                  low_order_hunp1[i] *= 2*std::pow(low_order_hnp1[i],VEL_FIX_POWER)/(std::pow(low_order_hnp1[i],VEL_FIX_POWER)+std::pow(aux,VEL_FIX_POWER));
                  low_order_hvnp1[i] *= 2*std::pow(low_order_hnp1[i],VEL_FIX_POWER)/(std::pow(low_order_hnp1[i],VEL_FIX_POWER)+std::pow(aux,VEL_FIX_POWER));
		  //low_order_hetanp1[i] = fmax(std::pow(low_order_hnp1[i],2.0),0.);
		  //  low_order_hwnp1[i] *= 2*std::pow(low_order_hnp1[i],VEL_FIX_POWER)/(std::pow(low_order_hnp1[i],VEL_FIX_POWER)+std::pow(aux,VEL_FIX_POWER));
                }
              int LOW_ORDER_SOLUTION=0; // FOR DEBUGGING
              if (LOW_ORDER_SOLUTION==1)
                {
                  globalResidual[offset_h+stride_h*i]   = low_order_hnp1[i];
                  globalResidual[offset_hu+stride_hu*i] = low_order_hunp1[i];
                  globalResidual[offset_hv+stride_hv*i] = low_order_hvnp1[i];
            	  globalResidual[offset_heta+stride_heta*i] = low_order_hetanp1[i];
            	  globalResidual[offset_hw+stride_hw*i] = low_order_hwnp1[i];
                }
              else
                {
                  if (LUMPED_MASS_MATRIX==1)
                    {
                      globalResidual[offset_h+stride_h*i]
			= hi - dt/mi*(ith_flux_term1
				      - ith_dHij_minus_muHij_times_hStarStates
				      - ith_muHij_times_hStates);
                      globalResidual[offset_hu+stride_hu*i]
			= hui - dt/mi*(ith_flux_term2
				       - ith_dHij_minus_muHij_times_huStarStates
				       - ith_muHij_times_huStates//);
		      + ith_friction_term2);
                      globalResidual[offset_hv+stride_hv*i]
			= 0.0* ( hvi - dt/mi*(ith_flux_term3
				       - ith_dHij_minus_muHij_times_hvStarStates
					      - ith_muHij_times_hvStates)// );
		      + ith_friction_term3);
		      globalResidual[offset_heta+stride_heta*i]
			= hetai - dt/mi*(ith_flux_term4
				       - ith_dHij_minus_muHij_times_hetaStarStates
				       - ith_muHij_times_hetaStates)
				       + dt/mi * force_term_hetai;
	              globalResidual[offset_hw+stride_hw*i]
			= hwi - dt/mi*(ith_flux_term5
				       - ith_dHij_minus_muHij_times_hwStarStates
				       - ith_muHij_times_hwStates)
				       + dt/mi * force_term_hwi;
                      // clean up potential negative water height due to machine precision
                      if (globalResidual[offset_h+stride_h*i] >= -1E-14)
                        globalResidual[offset_h+stride_h*i] = fmax(globalResidual[offset_h+stride_h*i],0.);
                    }
                  else
                    {
                      // Distribute residual
                      // NOTE: MASS MATRIX IS CONSISTENT
                      globalResidual[offset_h+stride_h*i]
                        += dt*(ith_flux_term1
                               - ith_dHij_minus_muHij_times_hStarStates
                               - ith_muHij_times_hStates);
                      globalResidual[offset_hu+stride_hu*i]
                        += dt*(ith_flux_term2
                               - ith_dHij_minus_muHij_times_huStarStates
                               - ith_muHij_times_huStates//);
		      + ith_friction_term2);
                      globalResidual[offset_hv+stride_hv*i]
                        += dt*(ith_flux_term3
                               - ith_dHij_minus_muHij_times_hvStarStates
                               - ith_muHij_times_hvStates//);
		      + ith_friction_term3);
		      globalResidual[offset_heta+stride_heta*i]
			+=dt*(ith_flux_term4
			      - ith_dHij_minus_muHij_times_hetaStarStates
			      - ith_muHij_times_hetaStates
			      + force_term_hetai);
		      globalResidual[offset_hw+stride_hw*i]
			+=dt*(ith_flux_term5
			      - ith_dHij_minus_muHij_times_hwStarStates
			      - ith_muHij_times_hwStates
			      + force_term_hwi);
                    }
                }
            }
          //std::cout << "Max of vel: "
          //        << max_ui << "\t"
          //        << max_vi << std::endl;
          //std::cout << "Max of hu: "
          //        << max_hui << "\t"
          //        << max_hvi << std::endl;
          //std::cout << "Max of hxu: "
          //        << max_hxui << "\t"
          //        << max_hxvi << std::endl;

          // ********** END OF LOOP IN DOFs ********** //
        }

      ///////////////////////////////////////////
      // ********** COMPUTE NORMALS ********** //
      ///////////////////////////////////////////
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
      ////////////////////////////////////////////////////
      // ********** END OF COMPUTING NORMALS ********** //
      ////////////////////////////////////////////////////
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
                             // h
                            int* csrRowIndeces_h_h,
                            int* csrColumnOffsets_h_h,
                            int* csrRowIndeces_h_hu,
                            int* csrColumnOffsets_h_hu,
                            int* csrRowIndeces_h_hv,
                            int* csrColumnOffsets_h_hv,
                            int* csrRowIndeces_h_heta,
                            int* csrColumnOffsets_h_heta,
                            int* csrRowIndeces_h_hw,
                            int* csrColumnOffsets_h_hw,
                            // hu
                            int* csrRowIndeces_hu_h,
                            int* csrColumnOffsets_hu_h,
                            int* csrRowIndeces_hu_hu,
                            int* csrColumnOffsets_hu_hu,
                            int* csrRowIndeces_hu_hv,
                            int* csrColumnOffsets_hu_hv,
                            int* csrRowIndeces_hu_heta,
                            int* csrColumnOffsets_hu_heta,
                            int* csrRowIndeces_hu_hw,
                            int* csrColumnOffsets_hu_hw,
                            // hv
                            int* csrRowIndeces_hv_h,
                            int* csrColumnOffsets_hv_h,
                            int* csrRowIndeces_hv_hu,
                            int* csrColumnOffsets_hv_hu,
                            int* csrRowIndeces_hv_hv,
                            int* csrColumnOffsets_hv_hv,
                            int* csrRowIndeces_hv_heta,
                            int* csrColumnOffsets_hv_heta,
                            int* csrRowIndeces_hv_hw,
                            int* csrColumnOffsets_hv_hw,
                            // heta
                            int* csrRowIndeces_heta_h,
                            int* csrColumnOffsets_heta_h,
                            int* csrRowIndeces_heta_hu,
                            int* csrColumnOffsets_heta_hu,
                            int* csrRowIndeces_heta_hv,
                            int* csrColumnOffsets_heta_hv,
                            int* csrRowIndeces_heta_heta,
                            int* csrColumnOffsets_heta_heta,
                            int* csrRowIndeces_heta_hw,
                            int* csrColumnOffsets_heta_hw,
                            //hw
                            int* csrRowIndeces_hw_h,
                            int* csrColumnOffsets_hw_h,
                            int* csrRowIndeces_hw_hu,
                            int* csrColumnOffsets_hw_hu,
                            int* csrRowIndeces_hw_hv,
                            int* csrColumnOffsets_hw_hv,
                            int* csrRowIndeces_hw_heta,
                            int* csrColumnOffsets_hw_heta,
                            int* csrRowIndeces_hw_hw,
                            int* csrColumnOffsets_hw_hw,
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
          //h
                elementJacobian_h_h[nDOF_test_element][nDOF_trial_element],
                elementJacobian_h_hu[nDOF_test_element][nDOF_trial_element],
                elementJacobian_h_hv[nDOF_test_element][nDOF_trial_element],
          elementJacobian_h_heta[nDOF_test_element][nDOF_trial_element],
          elementJacobian_h_hw[nDOF_test_element][nDOF_trial_element],
          //hu
                elementJacobian_hu_h[nDOF_test_element][nDOF_trial_element],
                elementJacobian_hu_hu[nDOF_test_element][nDOF_trial_element],
                elementJacobian_hu_hv[nDOF_test_element][nDOF_trial_element],
          elementJacobian_hu_heta[nDOF_test_element][nDOF_trial_element],
          elementJacobian_hu_hw[nDOF_test_element][nDOF_trial_element],
          //hv
                elementJacobian_hv_h[nDOF_test_element][nDOF_trial_element],
                elementJacobian_hv_hu[nDOF_test_element][nDOF_trial_element],
          elementJacobian_hv_hv[nDOF_test_element][nDOF_trial_element],
          elementJacobian_hv_heta[nDOF_test_element][nDOF_trial_element],
                elementJacobian_hv_hw[nDOF_test_element][nDOF_trial_element],
          //heta
                elementJacobian_heta_h[nDOF_test_element][nDOF_trial_element],
                elementJacobian_heta_hu[nDOF_test_element][nDOF_trial_element],
          elementJacobian_heta_hv[nDOF_test_element][nDOF_trial_element],
          elementJacobian_heta_heta[nDOF_test_element][nDOF_trial_element],
                elementJacobian_heta_hw[nDOF_test_element][nDOF_trial_element],
          //hw
                elementJacobian_hw_h[nDOF_test_element][nDOF_trial_element],
                elementJacobian_hw_hu[nDOF_test_element][nDOF_trial_element],
          elementJacobian_hw_hv[nDOF_test_element][nDOF_trial_element],
          elementJacobian_hw_heta[nDOF_test_element][nDOF_trial_element],
                elementJacobian_hw_hw[nDOF_test_element][nDOF_trial_element];

          for (int i=0;i<nDOF_test_element;i++)
            for (int j=0;j<nDOF_trial_element;j++)
              {
                // h
                            elementJacobian_h_h[i][j]=0.0;
                            elementJacobian_h_hu[i][j]=0.0;
                            elementJacobian_h_hv[i][j]=0.0;
                elementJacobian_h_heta[i][j]=0.0;
                elementJacobian_h_hw[i][j]=0.0;
                // hu
                            elementJacobian_hu_h[i][j]=0.0;
                            elementJacobian_hu_hu[i][j]=0.0;
                            elementJacobian_hu_hv[i][j]=0.0;
                elementJacobian_hu_heta[i][j]=0.0;
                elementJacobian_hu_hw[i][j]=0.0;
                // hv
                            elementJacobian_hv_h[i][j]=0.0;
                            elementJacobian_hv_hu[i][j]=0.0;
                            elementJacobian_hv_hv[i][j]=0.0;
                elementJacobian_hv_heta[i][j]=0.0;
                elementJacobian_hv_hw[i][j]=0.0;
                // heta
                            elementJacobian_heta_h[i][j]=0.0;
                            elementJacobian_heta_hu[i][j]=0.0;
                            elementJacobian_heta_hv[i][j]=0.0;
                elementJacobian_heta_heta[i][j]=0.0;
                elementJacobian_heta_hw[i][j]=0.0;
                // hw
                            elementJacobian_hw_h[i][j]=0.0;
                            elementJacobian_hw_hu[i][j]=0.0;
                            elementJacobian_hw_hv[i][j]=0.0;
                elementJacobian_hw_heta[i][j]=0.0;
                elementJacobian_hw_hw[i][j]=0.0;
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
                      elementJacobian_h_heta[i][j] += 0;
            		      elementJacobian_h_hw[i][j] += 0;

                      //////////////////////
                      // u: u_h, u_u, u_v //
                      //////////////////////
                      elementJacobian_hu_h[i][j] += 0;
                      elementJacobian_hu_hu[i][j] += vel_trial_ref[k*nDOF_trial_element+j]*vel_test_dV[i];
                      elementJacobian_hu_hv[i][j] += 0;
                      elementJacobian_hu_heta[i][j] += 0;
            		      elementJacobian_hu_hw[i][j] += 0;

                      //////////////////////
                      // v: v_h, v_u, v_v //
                      //////////////////////
                      elementJacobian_hv_h[i][j] += 0;
                      elementJacobian_hv_hu[i][j] += 0;
                      elementJacobian_hv_hv[i][j] += vel_trial_ref[k*nDOF_trial_element+j]*vel_test_dV[i];
                      elementJacobian_hv_heta[i][j] += 0;
                      elementJacobian_hv_hw[i][j] += 0;

                                  //////////////////////
                                  // heta: v_h, v_u, v_v //
                                  //////////////////////
                                  elementJacobian_heta_h[i][j] += 0;
                                  elementJacobian_heta_hu[i][j] += 0;
                                  elementJacobian_heta_hv[i][j] += 0;
                      elementJacobian_heta_heta[i][j] += vel_trial_ref[k*nDOF_trial_element+j]*vel_test_dV[i];
                      elementJacobian_heta_hw[i][j] += 0;

                                  //////////////////////
                                  // hw: v_h, v_u, v_v //
                                  //////////////////////
                                  elementJacobian_hw_h[i][j] += 0;
                                  elementJacobian_hw_hu[i][j] += 0;
                                  elementJacobian_hw_hv[i][j] += 0;
                      elementJacobian_hw_heta[i][j] += 0;
                      elementJacobian_hw_hw[i][j] += vel_trial_ref[k*nDOF_trial_element+j]*vel_test_dV[i];

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
                  globalJacobian[csrRowIndeces_h_heta[eN_i] + csrColumnOffsets_h_heta[eN_i_j]] += elementJacobian_h_heta[i][j];
                  globalJacobian[csrRowIndeces_h_hw[eN_i] + csrColumnOffsets_h_hw[eN_i_j]] += elementJacobian_h_hw[i][j];

                  globalJacobian[csrRowIndeces_hu_h[eN_i] + csrColumnOffsets_hu_h[eN_i_j]] += elementJacobian_hu_h[i][j];
                  globalJacobian[csrRowIndeces_hu_hu[eN_i] + csrColumnOffsets_hu_hu[eN_i_j]] += elementJacobian_hu_hu[i][j];
                  globalJacobian[csrRowIndeces_hu_hv[eN_i] + csrColumnOffsets_hu_hv[eN_i_j]] += elementJacobian_hu_hv[i][j];
                  globalJacobian[csrRowIndeces_hu_heta[eN_i] + csrColumnOffsets_hu_heta[eN_i_j]] += elementJacobian_hu_heta[i][j];
            		  globalJacobian[csrRowIndeces_hu_hw[eN_i] + csrColumnOffsets_hu_hw[eN_i_j]] += elementJacobian_hu_hw[i][j];

                  globalJacobian[csrRowIndeces_hv_h[eN_i] + csrColumnOffsets_hv_h[eN_i_j]] += elementJacobian_hv_h[i][j];
                  globalJacobian[csrRowIndeces_hv_hu[eN_i] + csrColumnOffsets_hv_hu[eN_i_j]] += elementJacobian_hv_hu[i][j];
                  globalJacobian[csrRowIndeces_hv_hv[eN_i] + csrColumnOffsets_hv_hv[eN_i_j]] += elementJacobian_hv_hv[i][j];
                  globalJacobian[csrRowIndeces_hv_heta[eN_i] + csrColumnOffsets_hv_heta[eN_i_j]] += elementJacobian_hv_heta[i][j];
            		  globalJacobian[csrRowIndeces_hv_hw[eN_i] + csrColumnOffsets_hv_hw[eN_i_j]] += elementJacobian_hv_hw[i][j];

                  globalJacobian[csrRowIndeces_heta_h[eN_i] + csrColumnOffsets_heta_h[eN_i_j]] += elementJacobian_heta_h[i][j];
                  globalJacobian[csrRowIndeces_heta_hu[eN_i] + csrColumnOffsets_heta_hu[eN_i_j]] += elementJacobian_heta_hu[i][j];
                  globalJacobian[csrRowIndeces_heta_hv[eN_i] + csrColumnOffsets_heta_hv[eN_i_j]] += elementJacobian_heta_hv[i][j];
            		  globalJacobian[csrRowIndeces_heta_heta[eN_i] + csrColumnOffsets_heta_heta[eN_i_j]] += elementJacobian_heta_heta[i][j];
            		  globalJacobian[csrRowIndeces_heta_hw[eN_i] + csrColumnOffsets_heta_hw[eN_i_j]] += elementJacobian_heta_hw[i][j];

                  globalJacobian[csrRowIndeces_hw_h[eN_i] + csrColumnOffsets_hw_h[eN_i_j]] += elementJacobian_hw_h[i][j];
                  globalJacobian[csrRowIndeces_hw_hu[eN_i] + csrColumnOffsets_hw_hu[eN_i_j]] += elementJacobian_hw_hu[i][j];
                  globalJacobian[csrRowIndeces_hw_hv[eN_i] + csrColumnOffsets_hw_hv[eN_i_j]] += elementJacobian_hw_hv[i][j];
            		  globalJacobian[csrRowIndeces_hw_heta[eN_i] + csrColumnOffsets_hw_heta[eN_i_j]] += elementJacobian_hw_heta[i][j];
            		  globalJacobian[csrRowIndeces_hw_hw[eN_i] + csrColumnOffsets_hw_hw[eN_i_j]] += elementJacobian_hw_hw[i][j];
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
  };//DSW2DCV

  inline DSW2DCV_base* newDSW2DCV(int nSpaceIn,
				  int nQuadraturePoints_elementIn,
				  int nDOF_mesh_trial_elementIn,
				  int nDOF_trial_elementIn,
				  int nDOF_test_elementIn,
				  int nQuadraturePoints_elementBoundaryIn,
				  int CompKernelFlag)
  {
    return proteus::chooseAndAllocateDiscretization2D<DSW2DCV_base,DSW2DCV,CompKernel>(nSpaceIn,
										       nQuadraturePoints_elementIn,
										       nDOF_mesh_trial_elementIn,
										       nDOF_trial_elementIn,
										       nDOF_test_elementIn,
										       nQuadraturePoints_elementBoundaryIn,
										       CompKernelFlag);
  }
}//proteus

#endif

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
  class SW2DCV_base
  {
  public:
    virtual ~SW2DCV_base(){}
      virtual void FCTStep(
                         double dt,
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
                            double* low_order_hnp1, //operators to construct low order solution
                            double* low_order_hunp1,
                            double* low_order_hvnp1,
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
                                   // FOR FCT
                                   double* dH_minus_dL,
                                   double* muH_minus_muL,
                                   double cE,
                                   int LUMPED_MASS_MATRIX,
                                   double dt,
                                   int LINEAR_FRICTION,
                                   double mannings,
                                   double* quantDOFs,
                                   int SECOND_CALL_CALCULATE_RESIDUAL,
                                   // NORMAL COMPONENTS
                                   int COMPUTE_NORMALS,
                                   double* normalx,
                                   double* normaly,
                                   // DISSIPATIVE LOW ORDER MATRIX
                                   double* dLow,
                                   int lstage
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
    virtual void calculateJacobian_SUPG(//element
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
          double hni = h_old[i];
          double huni = hu_old[i];
          double hvni = hv_old[i];
          double Zi = b_dof[i];
          double mi = lumped_mass_matrix[i];
          double one_over_hiReg = 2*hni/(hni*hni+std::pow(fmax(hni,hEps),2)); //hEps

          double ith_Limiter_times_FluxCorrectionMatrix1 = 0.;
          double ith_Limiter_times_FluxCorrectionMatrix2 = 0.;
          double ith_Limiter_times_FluxCorrectionMatrix3 = 0.;
          double Rnegi = Rneg[i];
          // LOOP OVER THE SPARSITY PATTERN (j-LOOP)//
          for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
            {
              int j = csrColumnOffsets_DofLoops[offset];
              // read some vectors
              double hnj = h_old[j];
              double hunj = hu_old[j];
              double hvnj = hv_old[j];
              double Zj = b_dof[j];
              double one_over_hjReg = 2*hnj/(hnj*hnj+std::pow(fmax(hnj,hEps),2)); //hEps

              // COMPUTE STAR SOLUTION // hStar, huStar and hvStar
              double hStarij  = fmax(0., hni + Zi - fmax(Zi,Zj));
              double huStarij = huni*hStarij*one_over_hiReg;
              double hvStarij = hvni*hStarij*one_over_hiReg;

              double hStarji  = fmax(0., hnj + Zj - fmax(Zi,Zj));
              double huStarji = hunj*hStarji*one_over_hjReg;
              double hvStarji = hvnj*hStarji*one_over_hjReg;

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

              // compute limiter based on water height
              double Lij = (FluxCorrectionMatrix1 > 0. ? std::min(1.,Rneg[j]) : std::min(Rnegi,1.));

              ith_Limiter_times_FluxCorrectionMatrix1 += Lij*FluxCorrectionMatrix1;
              ith_Limiter_times_FluxCorrectionMatrix2 += Lij*FluxCorrectionMatrix2;
              ith_Limiter_times_FluxCorrectionMatrix3 += Lij*FluxCorrectionMatrix3;
              //update ij
              ij+=1;
            }

          double one_over_mi = 1.0/lumped_mass_matrix[i];
          limited_hnp1[i]  = low_order_hnp1[i]  + one_over_mi*ith_Limiter_times_FluxCorrectionMatrix1;
          limited_hunp1[i] = low_order_hunp1[i] + one_over_mi*ith_Limiter_times_FluxCorrectionMatrix2;
          limited_hvnp1[i] = low_order_hvnp1[i] + one_over_mi*ith_Limiter_times_FluxCorrectionMatrix3;

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
              /*                                     elementDiameter[eN], */
              /*                                     pdeResidual_h, */
              /*                                     grad_h, */
              /*                                     q_numDiff_h[eN_k]); */

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
      /*        {  */
      /*          register int ebN = exteriorElementBoundariesArray[ebNE],  */
      /*            eN  = elementBoundaryElementsArray[ebN*2+0], */
      /*            ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+0], */
      /*            eN_nDOF_trial_element = eN*nDOF_trial_element; */
      /*          register double elementResidual_h[nDOF_test_element], */
      /*            elementResidual_hu[nDOF_test_element], */
      /*            elementResidual_hv[nDOF_test_element], */
      /*            eps_rho,eps_mu; */
      /*          for (int i=0;i<nDOF_test_element;i++) */
      /*            { */
      /*              elementResidual_h[i]=0.0; */
      /*              elementResidual_hu[i]=0.0; */
      /*              elementResidual_hv[i]=0.0; */
      /*            } */
      /*          for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++)  */
      /*            {  */
      /*              register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb, */
      /*                ebNE_kb_nSpace = ebNE_kb*nSpace, */
      /*                ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb, */
      /*                ebN_local_kb_nSpace = ebN_local_kb*nSpace; */
      /*              register double h_ext=0.0, */
      /*                u_ext=0.0, */
      /*                v_ext=0.0, */
      /*                grad_h_ext[nSpace], */
      /*                grad_hu_ext[nSpace], */
      /*                grad_hv_ext[nSpace], */
      /*                mom_hu_acc_ext=0.0, */
      /*                dmom_hu_acc_hu_ext=0.0, */
      /*                mom_hv_acc_ext=0.0, */
      /*                dmom_hv_acc_hv_ext=0.0, */
      /*                mass_adv_ext[nSpace], */
      /*                dmass_adv_hu_ext[nSpace], */
      /*                dmass_adv_hv_ext[nSpace], */
      /*                mom_hu_adv_ext[nSpace], */
      /*                dmom_hu_adv_hu_ext[nSpace], */
      /*                dmom_hu_adv_hv_ext[nSpace], */
      /*                mom_hv_adv_ext[nSpace], */
      /*                dmom_hv_adv_hu_ext[nSpace], */
      /*                dmom_hv_adv_hv_ext[nSpace], */
      /*                mom_hu_diff_ten_ext[nSpace], */
      /*                mom_hv_diff_ten_ext[nSpace], */
      /*                mom_huhv_diff_ten_ext[1], */
      /*                mom_hvhu_diff_ten_ext[1], */
      /*                mom_hu_source_ext=0.0, */
      /*                mom_hv_source_ext=0.0, */
      /*                mom_hu_ham_ext=0.0, */
      /*                dmom_hu_ham_grad_h_ext[nSpace], */
      /*                mom_hv_ham_ext=0.0, */
      /*                dmom_hv_ham_grad_h_ext[nSpace], */
      /*                dmom_hu_adv_h_ext[nSpace], */
      /*                dmom_hv_adv_h_ext[nSpace], */
      /*                flux_mass_ext=0.0, */
      /*                flux_mom_hu_adv_ext=0.0, */
      /*                flux_mom_hv_adv_ext=0.0, */
      /*                flux_mom_hu_diff_ext=0.0, */
      /*                flux_mom_hv_diff_ext=0.0, */
      /*                bc_h_ext=0.0, */
      /*                bc_hu_ext=0.0, */
      /*                bc_hv_ext=0.0, */
      /*                bc_mom_hu_acc_ext=0.0, */
      /*                bc_dmom_hu_acc_hu_ext=0.0, */
      /*                bc_mom_hv_acc_ext=0.0, */
      /*                bc_dmom_hv_acc_hv_ext=0.0, */
      /*                bc_mass_adv_ext[nSpace], */
      /*                bc_dmass_adv_hu_ext[nSpace], */
      /*                bc_dmass_adv_hv_ext[nSpace], */
      /*                bc_mom_hu_adv_ext[nSpace], */
      /*                bc_dmom_hu_adv_hu_ext[nSpace], */
      /*                bc_dmom_hu_adv_hv_ext[nSpace], */
      /*                bc_mom_hv_adv_ext[nSpace], */
      /*                bc_dmom_hv_adv_hu_ext[nSpace], */
      /*                bc_dmom_hv_adv_hv_ext[nSpace], */
      /*                bc_mom_hu_diff_ten_ext[nSpace], */
      /*                bc_mom_hv_diff_ten_ext[nSpace], */
      /*                bc_mom_huhv_diff_ten_ext[1], */
      /*                bc_mom_hvhu_diff_ten_ext[1], */
      /*                bc_mom_hu_source_ext=0.0, */
      /*                bc_mom_hv_source_ext=0.0, */
      /*                bc_mom_hu_ham_ext=0.0, */
      /*                bc_dmom_hu_ham_grad_h_ext[nSpace], */
      /*                bc_mom_hv_ham_ext=0.0, */
      /*                bc_dmom_hv_ham_grad_h_ext[nSpace], */
      /*                jac_ext[nSpace*nSpace], */
      /*                jacDet_ext, */
      /*                jacInv_ext[nSpace*nSpace], */
      /*                boundaryJac[nSpace*(nSpace-1)], */
      /*                metricTensor[(nSpace-1)*(nSpace-1)], */
      /*                metricTensorDetSqrt, */
      /*                dS,h_test_dS[nDOF_test_element],vel_test_dS[nDOF_test_element], */
      /*                h_grad_trial_trace[nDOF_trial_element*nSpace],vel_grad_trial_trace[nDOF_trial_element*nSpace], */
      /*                normal[3],x_ext,y_ext,xt_ext,yt_ext,integralScaling, */
      /*                G[nSpace*nSpace],G_dd_G,tr_G,h_penalty; */
      /*              compute information about mapping from reference element to physical element */
      /*              ck.calculateMapping_elementBoundary(eN, */
      /*                                                  ebN_local, */
      /*                                                  kb, */
      /*                                                  ebN_local_kb, */
      /*                                                  mesh_dof, */
      /*                                                  mesh_l2g, */
      /*                                                  mesh_trial_trace_ref, */
      /*                                                  mesh_grad_trial_trace_ref, */
      /*                                                  boundaryJac_ref, */
      /*                                                  jac_ext, */
      /*                                                  jacDet_ext, */
      /*                                                  jacInv_ext, */
      /*                                                  boundaryJac, */
      /*                                                  metricTensor, */
      /*                                                  metricTensorDetSqrt, */
      /*                                                  normal_ref, */
      /*                                                  normal, */
      /*                                                  x_ext,y_ext); */
      /*              ck.calculateMappingVelocity_elementBoundary(eN, */
      /*                                                          ebN_local, */
      /*                                                          kb, */
      /*                                                          ebN_local_kb, */
      /*                                                          mesh_velocity_dof, */
      /*                                                          mesh_l2g, */
      /*                                                          mesh_trial_trace_ref, */
      /*                                                          xt_ext,yt_ext, */
      /*                                                          normal, */
      /*                                                          boundaryJac, */
      /*                                                          metricTensor, */
      /*                                                          integralScaling); */
      /*              xt_ext=0.0;yt_ext=0.0; */
      /*              std::cout<<"xt_ext "<<xt_ext<<'\t'<<yt_ext<<'\t'<<std::endl; */
      /*              std::cout<<"integralScaling - metricTensorDetSrt ==============================="<<integralScaling-metricTensorDetSqrt<<std::endl; */
      /*              dS = ((1.0-MOVING_DOMAIN)*metricTensorDetSqrt + MOVING_DOMAIN*integralScaling)*dS_ref[kb]; */
      /*              compute shape and solution information */
      /*              shape */
      /*              ck.gradTrialFromRef(&h_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,h_grad_trial_trace); */
      /*              ck.gradTrialFromRef(&vel_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,vel_grad_trial_trace); */
      /*              solution and gradients     */
      /*              ck.valFromDOF(h_dof,&h_l2g[eN_nDOF_trial_element],&h_trial_trace_ref[ebN_local_kb*nDOF_test_element],h_ext); */
      /*              ck.valFromDOF(hu_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_trace_ref[ebN_local_kb*nDOF_test_element],u_ext); */
      /*              ck.valFromDOF(hv_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_trace_ref[ebN_local_kb*nDOF_test_element],v_ext); */
      /*              ck.gradFromDOF(h_dof,&h_l2g[eN_nDOF_trial_element],h_grad_trial_trace,grad_h_ext); */
      /*              ck.gradFromDOF(hu_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial_trace,grad_hu_ext); */
      /*              ck.gradFromDOF(hv_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial_trace,grad_hv_ext); */
      /*              precalculate test function products with integration weights */
      /*              for (int j=0;j<nDOF_trial_element;j++) */
      /*                { */
      /*                  h_test_dS[j] = h_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS; */
      /*                  vel_test_dS[j] = vel_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS; */
      /*                } */
      /*              // */
      /*              //debugging section for finite element calculations on exterior */
      /*              // */
      /*              std::cout<<"ebNE = "<<ebNE<<" kb = "<<kb<<std::endl; */
      /*              for (int j=0;j<nDOF_trial_element;j++) */
      /*                { */
      /*                  std::cout<<"h_trial_trace["<<j<<"] "<<h_trial_trace_ref[ebN_local_kb*nDOF_trial_element+j]<<std::endl; */
      /*                  std::cout<<"vel_trial_trace["<<j<<"] "<<vel_trial_trace_ref[ebN_local_kb*nDOF_trial_element+j]<<std::endl; */
      /*                  std::cout<<"h_test_dS["<<j<<"] "<<h_test_dS[j]<<std::endl; */
      /*                  std::cout<<"vel_test_dS["<<j<<"] "<<vel_test_dS[j]<<std::endl; */
      /*                  for (int I=0;I<nSpace;I++) */
      /*                { */
      /*                  std::cout<<"h_grad_trial_trace["<<j<<","<<I<<"] "<<h_grad_trial_trace[j*nSpace+I]<<std::endl; */
      /*                  std::cout<<"vel_grad_trial_trace["<<j<<","<<I<<"] "<<vel_grad_trial_trace[j*nSpace+I]<<std::endl; */
      /*                } */
      /*                } */
      /*              std::cout<<"h_ext "<<h_ext<<std::endl; */
      /*              std::cout<<"u_ext "<<u_ext<<std::endl; */
      /*              std::cout<<"v_ext "<<v_ext<<std::endl; */
      /*              for(int I=0;I<nSpace;I++) */
      /*                { */
      /*                  std::cout<<"grad_h_ext["<<I<<"] "<<grad_h_ext[I]<<std::endl; */
      /*                  std::cout<<"grad_hu_ext["<<I<<"] "<<grad_hu_ext[I]<<std::endl; */
      /*                  std::cout<<"grad_hv_ext["<<I<<"] "<<grad_hv_ext[I]<<std::endl; */
      /*                } */
      /*              */
      /*              load the boundary values */
      /*              */
      /*              bc_h_ext = isDOFBoundary_h[ebNE_kb]*ebqe_bc_h_ext[ebNE_kb]+(1-isDOFBoundary_h[ebNE_kb])*h_ext; */
      /*              bc_hu_ext = isDOFBoundary_hu[ebNE_kb]*ebqe_bc_hu_ext[ebNE_kb]+(1-isDOFBoundary_hu[ebNE_kb])*u_ext; */
      /*              bc_hv_ext = isDOFBoundary_hv[ebNE_kb]*ebqe_bc_hv_ext[ebNE_kb]+(1-isDOFBoundary_hv[ebNE_kb])*v_ext; */
      /*               */
      /*              calculate the pde coefficients using the solution and the boundary values for the solution  */
      /*               */
      /*              cek debug */
      /*              eps_rho=0.1; */
      /*              eps_mu=0.1; */
      /*              evaluateCoefficients(eps_rho, */
      /*                                   eps_mu, */
      /*                                   sigma, */
      /*                                   rho_0, */
      /*                                   nu_0, */
      /*                                   rho_1, */
      /*                                   nu_1, */
      /*                                   g, */
      /*                                   h_ext, */
      /*                                   grad_h_ext, */
      /*                                   u_ext, */
      /*                                   v_ext, */
      /*                                   mom_hu_acc_ext, */
      /*                                   dmom_hu_acc_hu_ext, */
      /*                                   mom_hv_acc_ext, */
      /*                                   dmom_hv_acc_hv_ext, */
      /*                                   mass_adv_ext, */
      /*                                   dmass_adv_hu_ext, */
      /*                                   dmass_adv_hv_ext, */
      /*                                   mom_hu_adv_ext, */
      /*                                   dmom_hu_adv_hu_ext, */
      /*                                   dmom_hu_adv_hv_ext, */
      /*                                   mom_hv_adv_ext, */
      /*                                   dmom_hv_adv_hu_ext, */
      /*                                   dmom_hv_adv_hv_ext, */
      /*                                   mom_hu_diff_ten_ext, */
      /*                                   mom_hv_diff_ten_ext, */
      /*                                   mom_huhv_diff_ten_ext, */
      /*                                   mom_hvhu_diff_ten_ext, */
      /*                                   mom_hu_source_ext, */
      /*                                   mom_hv_source_ext, */
      /*                                   mom_hu_ham_ext, */
      /*                                   dmom_hu_ham_grad_h_ext, */
      /*                                   mom_hv_ham_ext, */
      /*                                   dmom_hv_ham_grad_h_ext);           */
      /*              evaluateCoefficients(eps_rho, */
      /*                                   eps_mu, */
      /*                                   sigma, */
      /*                                   rho_0, */
      /*                                   nu_0, */
      /*                                   rho_1, */
      /*                                   nu_1, */
      /*                                   g, */
      /*                                   bc_h_ext, */
      /*                                   grad_h_ext,cek should't be used */
      /*                                   bc_hu_ext, */
      /*                                   bc_hv_ext, */
      /*                                   bc_mom_hu_acc_ext, */
      /*                                   bc_dmom_hu_acc_hu_ext, */
      /*                                   bc_mom_hv_acc_ext, */
      /*                                   bc_dmom_hv_acc_hv_ext, */
      /*                                   bc_mass_adv_ext, */
      /*                                   bc_dmass_adv_hu_ext, */
      /*                                   bc_dmass_adv_hv_ext, */
      /*                                   bc_mom_hu_adv_ext, */
      /*                                   bc_dmom_hu_adv_hu_ext, */
      /*                                   bc_dmom_hu_adv_hv_ext, */
      /*                                   bc_mom_hv_adv_ext, */
      /*                                   bc_dmom_hv_adv_hu_ext, */
      /*                                   bc_dmom_hv_adv_hv_ext, */
      /*                                   bc_mom_hu_diff_ten_ext, */
      /*                                   bc_mom_hv_diff_ten_ext, */
      /*                                   bc_mom_huhv_diff_ten_ext, */
      /*                                   bc_mom_hvhu_diff_ten_ext, */
      /*                                   bc_mom_hu_source_ext, */
      /*                                   bc_mom_hv_source_ext, */
      /*                                   bc_mom_hu_ham_ext, */
      /*                                   bc_dmom_hu_ham_grad_h_ext, */
      /*                                   bc_mom_hv_ham_ext, */
      /*                                   bc_dmom_hv_ham_grad_h_ext);           */
      /*              */
      /*              moving domain */
      /*              */
      /*              mass_adv_ext[0] -= MOVING_DOMAIN*xt_ext; */
      /*              mass_adv_ext[1] -= MOVING_DOMAIN*yt_ext; */

      /*              mom_hu_adv_ext[0] -= MOVING_DOMAIN*mom_hu_acc_ext*xt_ext; */
      /*              mom_hu_adv_ext[1] -= MOVING_DOMAIN*mom_hu_acc_ext*yt_ext; */
      /*              dmom_hu_adv_hu_ext[0] -= MOVING_DOMAIN*dmom_hu_acc_hu_ext*xt_ext; */
      /*              dmom_hu_adv_hu_ext[1] -= MOVING_DOMAIN*dmom_hu_acc_hu_ext*yt_ext; */

      /*              mom_hv_adv_ext[0] -= MOVING_DOMAIN*mom_hv_acc_ext*xt_ext; */
      /*              mom_hv_adv_ext[1] -= MOVING_DOMAIN*mom_hv_acc_ext*yt_ext; */
      /*              dmom_hv_adv_hv_ext[0] -= MOVING_DOMAIN*dmom_hv_acc_hv_ext*xt_ext; */
      /*              dmom_hv_adv_hv_ext[1] -= MOVING_DOMAIN*dmom_hv_acc_hv_ext*yt_ext; */

      /*              bc's */
      /*              bc_mom_hu_adv_ext[0] -= MOVING_DOMAIN*bc_mom_hu_acc_ext*xt_ext; */
      /*              bc_mom_hu_adv_ext[1] -= MOVING_DOMAIN*bc_mom_hu_acc_ext*yt_ext; */

      /*              bc_mom_hv_adv_ext[0] -= MOVING_DOMAIN*bc_mom_hv_acc_ext*xt_ext; */
      /*              bc_mom_hv_adv_ext[1] -= MOVING_DOMAIN*bc_mom_hv_acc_ext*yt_ext; */

      /*               */
      /*              calculate the numerical fluxes  */
      /*               */
      /*              cek debug */
      /*              ebqe_penalty_ext[ebNE_kb] = 10.0; */
      /*              */
      /*              ck.calculateGScale(G,normal,h_penalty); */
      /*              h_penalty = 10.0/h_penalty; */
      /*              cek debug, do it the old way */
      /*              h_penalty = 100.0/elementDiameter[eN]; */
      /*              exteriorNumericalAdvectiveFlux(isDOFBoundary_h[ebNE_kb], */
      /*                                             isDOFBoundary_hu[ebNE_kb], */
      /*                                             isDOFBoundary_hv[ebNE_kb], */
      /*                                             isAdvectiveFluxBoundary_h[ebNE_kb], */
      /*                                             isAdvectiveFluxBoundary_hu[ebNE_kb], */
      /*                                             isAdvectiveFluxBoundary_hv[ebNE_kb], */
      /*                                             dmom_hu_ham_grad_h_ext[0],=1/rho, */
      /*                                             normal, */
      /*                                             bc_h_ext, */
      /*                                             bc_mass_adv_ext, */
      /*                                             bc_mom_hu_adv_ext, */
      /*                                             bc_mom_hv_adv_ext, */
      /*                                             ebqe_bc_flux_mass_ext[ebNE_kb], */
      /*                                             ebqe_bc_flux_mom_hu_adv_ext[ebNE_kb], */
      /*                                             ebqe_bc_flux_mom_hv_adv_ext[ebNE_kb], */
      /*                                             h_ext, */
      /*                                             mass_adv_ext, */
      /*                                             mom_hu_adv_ext, */
      /*                                             mom_hv_adv_ext, */
      /*                                             dmass_adv_hu_ext, */
      /*                                             dmass_adv_hv_ext, */
      /*                                             dmom_hu_adv_h_ext, */
      /*                                             dmom_hu_adv_hu_ext, */
      /*                                             dmom_hu_adv_hv_ext, */
      /*                                             dmom_hv_adv_h_ext, */
      /*                                             dmom_hv_adv_hu_ext, */
      /*                                             dmom_hv_adv_hv_ext, */
      /*                                             flux_mass_ext, */
      /*                                             flux_mom_hu_adv_ext, */
      /*                                             flux_mom_hv_adv_ext, */
      /*                                             &ebqe_velocity[ebNE_kb_nSpace]); */
      /*              cek todo need to switch to full stress and add adjoint consistency */
      /*              exteriorNumericalDiffusiveFlux(eps_rho, */
      /*                                             sdInfo_hu_hu_rowptr, */
      /*                                             sdInfo_hu_hu_colind, */
      /*                                             isDOFBoundary_hu[ebNE_kb], */
      /*                                             isDiffusiveFluxBoundary_hu[ebNE_kb], */
      /*                                             normal, */
      /*                                             bc_mom_hu_diff_ten_ext, */
      /*                                             bc_hu_ext, */
      /*                                             ebqe_bc_flux_hu_diff_ext[ebNE_kb], */
      /*                                             mom_hu_diff_ten_ext, */
      /*                                             grad_hu_ext, */
      /*                                             u_ext, */
      /*                                             h_penalty,ebqe_penalty_ext[ebNE_kb], */
      /*                                             flux_mom_hu_diff_ext); */
      /*              exteriorNumericalDiffusiveFlux(eps_rho, */
      /*                                             sdInfo_hv_hv_rowptr, */
      /*                                             sdInfo_hv_hv_colind, */
      /*                                             isDOFBoundary_hv[ebNE_kb], */
      /*                                             isDiffusiveFluxBoundary_hv[ebNE_kb], */
      /*                                             normal, */
      /*                                             bc_mom_hv_diff_ten_ext, */
      /*                                             bc_hv_ext, */
      /*                                             ebqe_bc_flux_hv_diff_ext[ebNE_kb], */
      /*                                             mom_hv_diff_ten_ext, */
      /*                                             grad_hv_ext, */
      /*                                             v_ext, */
      /*                                             h_penalty,ebqe_penalty_ext[ebNE_kb], */
      /*                                             flux_mom_hv_diff_ext); */
      /*              flux[ebN*nQuadraturePoints_elementBoundary+kb] = flux_mass_ext; */
      /*              flux[ebN*nQuadraturePoints_elementBoundary+kb] = 0.0;//cek debug */
      /*              */
      /*              update residuals */
      /*              */
      /*              for (int i=0;i<nDOF_test_element;i++) */
      /*                { */
      /*                  elementResidual_h[i] += ck.ExteriorElementBoundaryFlux(flux_mass_ext,h_test_dS[i]); */
      /*                  globalConservationError += ck.ExteriorElementBoundaryFlux(flux_mass_ext,h_test_dS[i]); */

      /*                  elementResidual_hu[i] += ck.ExteriorElementBoundaryFlux(flux_mom_hu_adv_ext,vel_test_dS[i])+ */
      /*                    ck.ExteriorElementBoundaryFlux(flux_mom_hu_diff_ext,vel_test_dS[i]);  */

      /*                  elementResidual_hv[i] += ck.ExteriorElementBoundaryFlux(flux_mom_hv_adv_ext,vel_test_dS[i]) + */
      /*                    ck.ExteriorElementBoundaryFlux(flux_mom_hv_diff_ext,vel_test_dS[i]);  */

      /*                }i */
      /*            }kb */
      /*          */
      /*          update the element and global residual storage */
      /*          */
      /*          for (int i=0;i<nDOF_test_element;i++) */
      /*            { */
      /*              int eN_i = eN*nDOF_test_element+i; */

      /*              elementResidual_h_save[eN_i] +=  elementResidual_h[i]; */

      /*              globalResidual[offset_h+stride_h*h_l2g[eN_i]]+=elementResidual_h[i]; */
      /*              globalResidual[offset_hu+stride_hu*vel_l2g[eN_i]]+=elementResidual_hu[i]; */
      /*              globalResidual[offset_hv+stride_hv*vel_l2g[eN_i]]+=elementResidual_hv[i]; */
      /*            }i */
      /*          // */
      /*          //debug */
      /*          // */
      /*          for(int i=0;i<nDOF_test_element;i++)  */
      /*            {  */
      /*                  std::cout<<"ebNE "<<ebNE<<" i "<<i<<std::endl; */
      /*                  std::cout<<"r_h"<<elementResidual_h[i]<<std::endl; */
      /*                  std::cout<<"r_hu"<<elementResidual_hu[i]<<std::endl; */
      /*                  std::cout<<"r_hv"<<elementResidual_hv[i]<<std::endl; */
      /*                } */

      /*        }ebNE */
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
              ck.valFromDOF(h_dof,&h_l2g[eN_nDOF_trial_element],&h_trial_ref[k*nDOF_trial_element],h);
              ck.valFromDOF(hu_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hu);
              ck.valFromDOF(hv_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hv);
              // get the solution at the lstage
              ck.valFromDOF(h_dof_old,&h_l2g[eN_nDOF_trial_element],&h_trial_ref[k*nDOF_trial_element],h_old);
              ck.valFromDOF(hu_dof_old,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hu_old);
              ck.valFromDOF(hv_dof_old,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],hv_old);
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
              double Zi = b_dof[i];
              double mi = lumped_mass_matrix[i];
              double one_over_hiReg = 2*hi/(hi*hi+std::pow(fmax(hi,hEps),2)); // hEps
              double ui = hui*one_over_hiReg;
              double vi = hvi*one_over_hiReg;

              //max_ui = fmax(max_ui,fabs(ui));
              //max_vi = fmax(max_vi,fabs(vi));
              //max_hui = fmax(max_hui,fabs(hui));
              //max_hvi = fmax(max_hvi,fabs(hvi));
              //max_hxui = fmax(max_hxui,fabs(hi*ui));
              //max_hxvi = fmax(max_hxvi,fabs(hi*vi));

              double ith_flux_term1=0., ith_flux_term2=0., ith_flux_term3=0.;
              // LOW ORDER DISSIPATIVE TERMS
              double
                ith_dLij_minus_muLij_times_hStarStates=0.,
                ith_dLij_minus_muLij_times_huStarStates=0.,
                ith_dLij_minus_muLij_times_hvStarStates=0.,
                ith_muLij_times_hStates=0.,
                ith_muLij_times_huStates=0.,
                ith_muLij_times_hvStates=0.;
              // HIGH ORDER DISSIPATIVE TERMS
              double
                ith_dHij_minus_muHij_times_hStarStates=0.,
                ith_dHij_minus_muHij_times_huStarStates=0.,
                ith_dHij_minus_muHij_times_hvStarStates=0.,
                ith_muHij_times_hStates=0.,
                ith_muHij_times_huStates=0.,
                ith_muHij_times_hvStates=0.;

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

                  // Nodal projection of fluxes
                  ith_flux_term1 += huj*Cx[ij] + hvj*Cy[ij]; // f1*C
                  //ith_flux_term1 += hj*(uj*Cx[ij] + vj*Cy[ij]); // f1*C
                  ith_flux_term2 += uj*huj*Cx[ij] + uj*hvj*Cy[ij] + g*hi*(hj+Zj)*Cx[ij];
                  ith_flux_term3 += vj*huj*Cx[ij] + vj*hvj*Cy[ij] + g*hi*(hj+Zj)*Cy[ij];

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

                      ith_dLij_minus_muLij_times_hStarStates  += (dLij - muLij)*(hStarji-hStarij);
                      ith_dLij_minus_muLij_times_huStarStates += (dLij - muLij)*(huStarji-huStarij);
                      ith_dLij_minus_muLij_times_hvStarStates += (dLij - muLij)*(hvStarji-hvStarij);

                      // compute muij times solution terms
                      ith_muHij_times_hStates  += muHij*(hj-hi);
                      ith_muHij_times_huStates += muHij*(huj-hui);
                      ith_muHij_times_hvStates += muHij*(hvj-hvi);

                      ith_muLij_times_hStates  += muLij*(hj-hi);
                      ith_muLij_times_huStates += muLij*(huj-hui);
                      ith_muLij_times_hvStates += muLij*(hvj-hvi);

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
                                                - ith_muLij_times_huStates
                                                + ith_friction_term2);
              low_order_hvnp1[i] = hvi - dt/mi*(ith_flux_term3
                                                - ith_dLij_minus_muLij_times_hvStarStates
                                                - ith_muLij_times_hvStates
                                                + ith_friction_term3);
              // FIX LOW ORDER SOLUTION //
              if (low_order_hnp1[i] < -1E-14 && dt < 1.0)
                {
                  std::cout << "dt taken: " << dt
                            << std::endl;
                  std::cout << "********.... "
                            << "Low order water height is negative: "
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
                }
              int LOW_ORDER_SOLUTION=0; // FOR DEBUGGING
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
                                      - ith_dHij_minus_muHij_times_hStarStates
                                      - ith_muHij_times_hStates);
                      globalResidual[offset_hu+stride_hu*i]
                        = hui - dt/mi*(ith_flux_term2
                                       - ith_dHij_minus_muHij_times_huStarStates
                                       - ith_muHij_times_huStates
                                       + ith_friction_term2);
                      globalResidual[offset_hv+stride_hv*i]
                        = hvi - dt/mi*(ith_flux_term3
                                       - ith_dHij_minus_muHij_times_hvStarStates
                                       - ith_muHij_times_hvStates
                                       + ith_friction_term3);
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
                               - ith_muHij_times_huStates
                               + ith_friction_term2);
                      globalResidual[offset_hv+stride_hv*i]
                        += dt*(ith_flux_term3
                               - ith_dHij_minus_muHij_times_hvStarStates
                               - ith_muHij_times_hvStates
                               + ith_friction_term3);
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

    void calculateJacobian_SUPG(//element
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
      /*        {  */
      /*          register int ebN = exteriorElementBoundariesArray[ebNE], */
      /*            eN  = elementBoundaryElementsArray[ebN*2+0], */
      /*            eN_nDOF_trial_element = eN*nDOF_trial_element, */
      /*            ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+0]; */
      /*          register double eps_rho,eps_mu; */
      /*          for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++)  */
      /*            {  */
      /*              register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb, */
      /*                ebNE_kb_nSpace = ebNE_kb*nSpace, */
      /*                ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb, */
      /*                ebN_local_kb_nSpace = ebN_local_kb*nSpace; */

      /*              register double h_ext=0.0, */
      /*                u_ext=0.0, */
      /*                v_ext=0.0, */
      /*                grad_h_ext[nSpace], */
      /*                grad_hu_ext[nSpace], */
      /*                grad_hv_ext[nSpace], */
      /*                mom_hu_acc_ext=0.0, */
      /*                dmom_hu_acc_hu_ext=0.0, */
      /*                mom_hv_acc_ext=0.0, */
      /*                dmom_hv_acc_hv_ext=0.0, */
      /*                mass_adv_ext[nSpace], */
      /*                dmass_adv_hu_ext[nSpace], */
      /*                dmass_adv_hv_ext[nSpace], */
      /*                mom_hu_adv_ext[nSpace], */
      /*                dmom_hu_adv_hu_ext[nSpace], */
      /*                dmom_hu_adv_hv_ext[nSpace], */
      /*                mom_hv_adv_ext[nSpace], */
      /*                dmom_hv_adv_hu_ext[nSpace], */
      /*                dmom_hv_adv_hv_ext[nSpace], */
      /*                mom_hu_diff_ten_ext[nSpace], */
      /*                mom_hv_diff_ten_ext[nSpace], */
      /*                mom_huhv_diff_ten_ext[1], */
      /*                mom_hvhu_diff_ten_ext[1], */
      /*                mom_hu_source_ext=0.0, */
      /*                mom_hv_source_ext=0.0, */
      /*                mom_hu_ham_ext=0.0, */
      /*                dmom_hu_ham_grad_h_ext[nSpace], */
      /*                mom_hv_ham_ext=0.0, */
      /*                dmom_hv_ham_grad_h_ext[nSpace], */
      /*                dmom_hu_adv_h_ext[nSpace], */
      /*                dmom_hv_adv_h_ext[nSpace], */
      /*                dflux_mass_hu_ext=0.0, */
      /*                dflux_mass_hv_ext=0.0, */
      /*                dflux_mom_hu_adv_h_ext=0.0, */
      /*                dflux_mom_hu_adv_hu_ext=0.0, */
      /*                dflux_mom_hu_adv_hv_ext=0.0, */
      /*                dflux_mom_hv_adv_h_ext=0.0, */
      /*                dflux_mom_hv_adv_hu_ext=0.0, */
      /*                dflux_mom_hv_adv_hv_ext=0.0, */
      /*                bc_h_ext=0.0, */
      /*                bc_hu_ext=0.0, */
      /*                bc_hv_ext=0.0, */
      /*                bc_mom_hu_acc_ext=0.0, */
      /*                bc_dmom_hu_acc_hu_ext=0.0, */
      /*                bc_mom_hv_acc_ext=0.0, */
      /*                bc_dmom_hv_acc_hv_ext=0.0, */
      /*                bc_mass_adv_ext[nSpace], */
      /*                bc_dmass_adv_hu_ext[nSpace], */
      /*                bc_dmass_adv_hv_ext[nSpace], */
      /*                bc_mom_hu_adv_ext[nSpace], */
      /*                bc_dmom_hu_adv_hu_ext[nSpace], */
      /*                bc_dmom_hu_adv_hv_ext[nSpace], */
      /*                bc_mom_hv_adv_ext[nSpace], */
      /*                bc_dmom_hv_adv_hu_ext[nSpace], */
      /*                bc_dmom_hv_adv_hv_ext[nSpace], */
      /*                bc_mom_hu_diff_ten_ext[nSpace], */
      /*                bc_mom_hv_diff_ten_ext[nSpace], */
      /*                bc_mom_huhv_diff_ten_ext[1], */
      /*                bc_mom_hvhu_diff_ten_ext[1], */
      /*                bc_mom_hu_source_ext=0.0, */
      /*                bc_mom_hv_source_ext=0.0, */
      /*                bc_mom_hu_ham_ext=0.0, */
      /*                bc_dmom_hu_ham_grad_h_ext[nSpace], */
      /*                bc_mom_hv_ham_ext=0.0, */
      /*                bc_dmom_hv_ham_grad_h_ext[nSpace], */
      /*                fluxJacobian_h_h[nDOF_trial_element], */
      /*                fluxJacobian_h_hu[nDOF_trial_element], */
      /*                fluxJacobian_h_hv[nDOF_trial_element], */
      /*                fluxJacobian_hu_h[nDOF_trial_element], */
      /*                fluxJacobian_hu_hu[nDOF_trial_element], */
      /*                fluxJacobian_hu_hv[nDOF_trial_element], */
      /*                fluxJacobian_hv_h[nDOF_trial_element], */
      /*                fluxJacobian_hv_hu[nDOF_trial_element], */
      /*                fluxJacobian_hv_hv[nDOF_trial_element], */
      /*                jac_ext[nSpace*nSpace], */
      /*                jacDet_ext, */
      /*                jacInv_ext[nSpace*nSpace], */
      /*                boundaryJac[nSpace*(nSpace-1)], */
      /*                metricTensor[(nSpace-1)*(nSpace-1)], */
      /*                metricTensorDetSqrt, */
      /*                h_grad_trial_trace[nDOF_trial_element*nSpace], */
      /*                vel_grad_trial_trace[nDOF_trial_element*nSpace], */
      /*                dS, */
      /*                h_test_dS[nDOF_test_element], */
      /*                vel_test_dS[nDOF_test_element], */
      /*                normal[3], */
      /*                x_ext,y_ext,xt_ext,yt_ext,integralScaling, */
      /*                G[nSpace*nSpace],G_dd_G,tr_G,h_hhi,h_henalty; */
      /*              ck.calculateMapping_elementBoundary(eN, */
      /*                                                  ebN_local, */
      /*                                                  kb, */
      /*                                                  ebN_local_kb, */
      /*                                                  mesh_dof, */
      /*                                                  mesh_l2g, */
      /*                                                  mesh_trial_trace_ref, */
      /*                                                  mesh_grad_trial_trace_ref, */
      /*                                                  boundaryJac_ref, */
      /*                                                  jac_ext, */
      /*                                                  jacDet_ext, */
      /*                                                  jacInv_ext, */
      /*                                                  boundaryJac, */
      /*                                                  metricTensor, */
      /*                                                  metricTensorDetSqrt, */
      /*                                                  normal_ref, */
      /*                                                  normal, */
      /*                                                  x_ext,y_ext); */
      /*              ck.calculateMappingVelocity_elementBoundary(eN, */
      /*                                                          ebN_local, */
      /*                                                          kb, */
      /*                                                          ebN_local_kb, */
      /*                                                          mesh_velocity_dof, */
      /*                                                          mesh_l2g, */
      /*                                                          mesh_trial_trace_ref, */
      /*                                                          xt_ext,yt_ext, */
      /*                                                          normal, */
      /*                                                          boundaryJac, */
      /*                                                          metricTensor, */
      /*                                                          integralScaling); */
      /*              //xt_ext=0.0;yt_ext=0.0; */
      /*              //std::cout<<"xt_ext "<<xt_ext<<'\t'<<yt_ext<<'\t'<<std::endl; */
      /*              dS = ((1.0-MOVING_DOMAIN)*metricTensorDetSqrt + MOVING_DOMAIN*integralScaling)*dS_ref[kb]; */
      /*              ck.calculateG(jacInv_ext,G,G_dd_G,tr_G); */
      /*              ck.calculateGScale(G,&ebqe_normal_hhi_ext[ebNE_kb_nSpace],h_hhi); */

      /*              eps_rho = epsFact_rho*(useMetrics*h_hhi+(1.0-useMetrics)*elementDiameter[eN]); */
      /*              eps_mu  = epsFact_mu *(useMetrics*h_hhi+(1.0-useMetrics)*elementDiameter[eN]); */

      /*              //compute shape and solution information */
      /*              //shape */
      /*              ck.gradTrialFromRef(&h_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,h_grad_trial_trace); */
      /*              ck.gradTrialFromRef(&vel_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,vel_grad_trial_trace); */
      /*              //solution and gradients   */
      /*              ck.valFromDOF(h_dof,&h_l2g[eN_nDOF_trial_element],&h_trial_trace_ref[ebN_local_kb*nDOF_test_element],h_ext); */
      /*              ck.valFromDOF(hu_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_trace_ref[ebN_local_kb*nDOF_test_element],u_ext); */
      /*              ck.valFromDOF(hv_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_trace_ref[ebN_local_kb*nDOF_test_element],v_ext); */
      /*              ck.gradFromDOF(h_dof,&h_l2g[eN_nDOF_trial_element],h_grad_trial_trace,grad_h_ext); */
      /*              ck.gradFromDOF(hu_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial_trace,grad_hu_ext); */
      /*              ck.gradFromDOF(hv_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial_trace,grad_hv_ext); */
      /*              //precalculate test function products with integration weights */
      /*              for (int j=0;j<nDOF_trial_element;j++) */
      /*                { */
      /*                  h_test_dS[j] = h_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS; */
      /*                  vel_test_dS[j] = vel_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS; */
      /*                } */
      /*              // // */
      /*              // //debugging section for finite element calculations on exterior */
      /*              // // */
      /*              // std::cout<<"ebNE = "<<ebNE<<" kb = "<<kb<<std::endl; */
      /*              // for (int j=0;j<nDOF_trial_element;j++) */
      /*              //   { */
      /*              //     std::cout<<"h_trial_trace["<<j<<"] "<<h_trial_trace_ref[ebN_local_kb*nDOF_trial_element+j]<<std::endl; */
      /*              //     std::cout<<"vel_trial_trace["<<j<<"] "<<vel_trial_trace_ref[ebN_local_kb*nDOF_trial_element+j]<<std::endl; */
      /*              //     std::cout<<"h_test_dS["<<j<<"] "<<h_test_dS[j]<<std::endl; */
      /*              //     std::cout<<"vel_test_dS["<<j<<"] "<<vel_test_dS[j]<<std::endl; */
      /*              //     for (int I=0;I<nSpace;I++) */
      /*              //        { */
      /*              //          std::cout<<"h_grad_trial_trace["<<j<<","<<I<<"] "<<h_grad_trial_trace[j*nSpace+I]<<std::endl; */
      /*              //          std::cout<<"vel_grad_trial_trace["<<j<<","<<I<<"] "<<vel_grad_trial_trace[j*nSpace+I]<<std::endl; */
      /*              //        } */
      /*              //   } */
      /*              // std::cout<<"h_ext "<<h_ext<<std::endl; */
      /*              // std::cout<<"u_ext "<<u_ext<<std::endl; */
      /*              // std::cout<<"v_ext "<<v_ext<<std::endl; */
      /*              // for(int I=0;I<nSpace;I++) */
      /*              //   { */
      /*              //     std::cout<<"grad_h_ext["<<I<<"] "<<grad_h_ext[I]<<std::endl; */
      /*              //     std::cout<<"grad_hu_ext["<<I<<"] "<<grad_hu_ext[I]<<std::endl; */
      /*              //     std::cout<<"grad_hv_ext["<<I<<"] "<<grad_hv_ext[I]<<std::endl; */
      /*              //   } */
      /*              // */
      /*              //load the boundary values */
      /*              // */
      /*              bc_h_ext = isDOFBoundary_h[ebNE_kb]*ebqe_bc_h_ext[ebNE_kb]+(1-isDOFBoundary_h[ebNE_kb])*h_ext; */
      /*              bc_hu_ext = isDOFBoundary_hu[ebNE_kb]*ebqe_bc_hu_ext[ebNE_kb]+(1-isDOFBoundary_hu[ebNE_kb])*u_ext; */
      /*              bc_hv_ext = isDOFBoundary_hv[ebNE_kb]*ebqe_bc_hv_ext[ebNE_kb]+(1-isDOFBoundary_hv[ebNE_kb])*v_ext; */
      /*              //  */
      /*              //calculate the internal and external trace of the pde coefficients  */
      /*              //  */
      /*              //cek debug */
      /*              //eps_rho=0.1; */
      /*              //eps_mu=0.1; */
      /*              evaluateCoefficients(eps_rho, */
      /*                                   eps_mu, */
      /*                                   sigma, */
      /*                                   rho_0, */
      /*                                   nu_0, */
      /*                                   rho_1, */
      /*                                   nu_1, */
      /*                                   g, */
      /*                                   ebqe_hhi_ext[ebNE_kb], */
      /*                                   &ebqe_normal_hhi_ext[ebNE_kb_nSpace], */
      /*                                   ebqe_kappa_hhi_ext[ebNE_kb], */
      /*                                   h_ext, */
      /*                                   grad_h_ext, */
      /*                                   u_ext, */
      /*                                   v_ext, */
      /*                                   mom_hu_acc_ext, */
      /*                                   dmom_hu_acc_hu_ext, */
      /*                                   mom_hv_acc_ext, */
      /*                                   dmom_hv_acc_hv_ext, */
      /*                                   mass_adv_ext, */
      /*                                   dmass_adv_hu_ext, */
      /*                                   dmass_adv_hv_ext, */
      /*                                   mom_hu_adv_ext, */
      /*                                   dmom_hu_adv_hu_ext, */
      /*                                   dmom_hu_adv_hv_ext, */
      /*                                   mom_hv_adv_ext, */
      /*                                   dmom_hv_adv_hu_ext, */
      /*                                   dmom_hv_adv_hv_ext, */
      /*                                   mom_hu_diff_ten_ext, */
      /*                                   mom_hv_diff_ten_ext, */
      /*                                   mom_huhv_diff_ten_ext, */
      /*                                   mom_hvhu_diff_ten_ext, */
      /*                                   mom_hu_source_ext, */
      /*                                   mom_hv_source_ext, */
      /*                                   mom_hu_ham_ext, */
      /*                                   dmom_hu_ham_grad_h_ext, */
      /*                                   mom_hv_ham_ext, */
      /*                                   dmom_hv_ham_grad_h_ext);           */
      /*              evaluateCoefficients(eps_rho, */
      /*                                   eps_mu, */
      /*                                   sigma, */
      /*                                   rho_0, */
      /*                                   nu_0, */
      /*                                   rho_1, */
      /*                                   nu_1, */
      /*                                   g, */
      /*                                   ebqe_hhi_ext[ebNE_kb], */
      /*                                   &ebqe_normal_hhi_ext[ebNE_kb_nSpace], */
      /*                                   ebqe_kappa_hhi_ext[ebNE_kb], */
      /*                                   bc_h_ext, */
      /*                                   grad_h_ext, //cek shouldn't be used */
      /*                                   bc_hu_ext, */
      /*                                   bc_hv_ext, */
      /*                                   bc_mom_hu_acc_ext, */
      /*                                   bc_dmom_hu_acc_hu_ext, */
      /*                                   bc_mom_hv_acc_ext, */
      /*                                   bc_dmom_hv_acc_hv_ext, */
      /*                                   bc_mass_adv_ext, */
      /*                                   bc_dmass_adv_hu_ext, */
      /*                                   bc_dmass_adv_hv_ext, */
      /*                                   bc_mom_hu_adv_ext, */
      /*                                   bc_dmom_hu_adv_hu_ext, */
      /*                                   bc_dmom_hu_adv_hv_ext, */
      /*                                   bc_mom_hv_adv_ext, */
      /*                                   bc_dmom_hv_adv_hu_ext, */
      /*                                   bc_dmom_hv_adv_hv_ext, */
      /*                                   bc_mom_hu_diff_ten_ext, */
      /*                                   bc_mom_hv_diff_ten_ext, */
      /*                                   bc_mom_huhv_diff_ten_ext, */
      /*                                   bc_mom_hvhu_diff_ten_ext, */
      /*                                   bc_mom_hu_source_ext, */
      /*                                   bc_mom_hv_source_ext, */
      /*                                   bc_mom_hu_ham_ext, */
      /*                                   bc_dmom_hu_ham_grad_h_ext, */
      /*                                   bc_mom_hv_ham_ext, */
      /*                                   bc_dmom_hv_ham_grad_h_ext);           */
      /*              // */
      /*              //moving domain */
      /*              // */
      /*              mass_adv_ext[0] -= MOVING_DOMAIN*xt_ext; */
      /*              mass_adv_ext[1] -= MOVING_DOMAIN*yt_ext; */

      /*              mom_hu_adv_ext[0] -= MOVING_DOMAIN*mom_hu_acc_ext*xt_ext; */
      /*              mom_hu_adv_ext[1] -= MOVING_DOMAIN*mom_hu_acc_ext*yt_ext; */
      /*              dmom_hu_adv_hu_ext[0] -= MOVING_DOMAIN*dmom_hu_acc_hu_ext*xt_ext; */
      /*              dmom_hu_adv_hu_ext[1] -= MOVING_DOMAIN*dmom_hu_acc_hu_ext*yt_ext; */

      /*              mom_hv_adv_ext[0] -= MOVING_DOMAIN*mom_hv_acc_ext*xt_ext; */
      /*              mom_hv_adv_ext[1] -= MOVING_DOMAIN*mom_hv_acc_ext*yt_ext; */
      /*              dmom_hv_adv_hv_ext[0] -= MOVING_DOMAIN*dmom_hv_acc_hv_ext*xt_ext; */
      /*              dmom_hv_adv_hv_ext[1] -= MOVING_DOMAIN*dmom_hv_acc_hv_ext*yt_ext; */

      /*              //moving domain bc's */
      /*              bc_mom_hu_adv_ext[0] -= MOVING_DOMAIN*bc_mom_hu_acc_ext*xt_ext; */
      /*              bc_mom_hu_adv_ext[1] -= MOVING_DOMAIN*bc_mom_hu_acc_ext*yt_ext; */

      /*              bc_mom_hv_adv_ext[0] -= MOVING_DOMAIN*bc_mom_hv_acc_ext*xt_ext; */
      /*              bc_mom_hv_adv_ext[1] -= MOVING_DOMAIN*bc_mom_hv_acc_ext*yt_ext; */

      /*              //  */
      /*              //calculate the numerical fluxes  */
      /*              //  */
      /*              exteriorNumericalAdvectiveFluxDerivatives(isDOFBoundary_h[ebNE_kb], */
      /*                                                        isDOFBoundary_hu[ebNE_kb], */
      /*                                                        isDOFBoundary_hv[ebNE_kb], */
      /*                                                        isAdvectiveFluxBoundary_h[ebNE_kb], */
      /*                                                        isAdvectiveFluxBoundary_hu[ebNE_kb], */
      /*                                                        isAdvectiveFluxBoundary_hv[ebNE_kb], */
      /*                                                        dmom_hu_ham_grad_h_ext[0],//=1/rho */
      /*                                                        normal, */
      /*                                                        bc_h_ext, */
      /*                                                        bc_mass_adv_ext, */
      /*                                                        bc_mom_hu_adv_ext, */
      /*                                                        bc_mom_hv_adv_ext, */
      /*                                                        ebqe_bc_flux_mass_ext[ebNE_kb], */
      /*                                                        ebqe_bc_flux_mom_hu_adv_ext[ebNE_kb], */
      /*                                                        ebqe_bc_flux_mom_hv_adv_ext[ebNE_kb], */
      /*                                                        h_ext, */
      /*                                                        mass_adv_ext, */
      /*                                                        mom_hu_adv_ext, */
      /*                                                        mom_hv_adv_ext, */
      /*                                                        dmass_adv_hu_ext, */
      /*                                                        dmass_adv_hv_ext, */
      /*                                                        dmom_hu_adv_h_ext, */
      /*                                                        dmom_hu_adv_hu_ext, */
      /*                                                        dmom_hu_adv_hv_ext, */
      /*                                                        dmom_hv_adv_h_ext, */
      /*                                                        dmom_hv_adv_hu_ext, */
      /*                                                        dmom_hv_adv_hv_ext, */
      /*                                                        dflux_mass_hu_ext, */
      /*                                                        dflux_mass_hv_ext, */
      /*                                                        dflux_mom_hu_adv_h_ext, */
      /*                                                        dflux_mom_hu_adv_hu_ext, */
      /*                                                        dflux_mom_hu_adv_hv_ext, */
      /*                                                        dflux_mom_hv_adv_h_ext, */
      /*                                                        dflux_mom_hv_adv_hu_ext, */
      /*                                                        dflux_mom_hv_adv_hv_ext); */
      /*              // */
      /*              //calculate the flux jacobian */
      /*              // */
      /*              ck.calculateGScale(G,normal,h_henalty); */
      /*              h_henalty = 10.0/h_henalty; */
      /*              //cek debug, do it the old way */
      /*              h_henalty = 100.0/elementDiameter[eN]; */
      /*              for (int j=0;j<nDOF_trial_element;j++) */
      /*                { */
      /*                  register int j_nSpace = j*nSpace,ebN_local_kb_j=ebN_local_kb*nDOF_trial_element+j; */
      /*                  //cek debug */
      /*                  //ebqe_henalty_ext[ebNE_kb] = 10.0; */
      /*                  // */
      /*                  //cek todo add full stress on boundaries */

      /*                  fluxJacobian_h_h[j]=0.0; */
      /*                  fluxJacobian_h_hu[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mass_hu_ext,vel_trial_trace_ref[ebN_local_kb_j]); */
      /*                  fluxJacobian_h_hv[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mass_hv_ext,vel_trial_trace_ref[ebN_local_kb_j]); */

      /*                  fluxJacobian_hu_h[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_hu_adv_h_ext,h_trial_trace_ref[ebN_local_kb_j]); */
      /*                  fluxJacobian_hu_hu[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_hu_adv_hu_ext,vel_trial_trace_ref[ebN_local_kb_j]) + */
      /*                    ExteriorNumericalDiffusiveFluxJacobian(eps_rho, */
      /*                                                           ebqe_hhi_ext[ebNE_kb], */
      /*                                                           sdInfo_hu_hu_rowptr, */
      /*                                                           sdInfo_hu_hu_colind, */
      /*                                                           isDOFBoundary_hu[ebNE_kb], */
      /*                                                           normal, */
      /*                                                           mom_hu_diff_ten_ext, */
      /*                                                           vel_trial_trace_ref[ebN_local_kb_j], */
      /*                                                           &vel_grad_trial_trace[j_nSpace], */
      /*                                                           h_henalty);//ebqe_henalty_ext[ebNE_kb]); */
      /*                  fluxJacobian_hu_hv[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_hu_adv_hv_ext,vel_trial_trace_ref[ebN_local_kb_j]); */

      /*                  fluxJacobian_hv_h[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_hv_adv_h_ext,h_trial_trace_ref[ebN_local_kb_j]); */
      /*                  fluxJacobian_hv_hu[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_hv_adv_hu_ext,vel_trial_trace_ref[ebN_local_kb_j]); */
      /*                  fluxJacobian_hv_hv[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_hv_adv_hv_ext,vel_trial_trace_ref[ebN_local_kb_j]) + */
      /*                    ExteriorNumericalDiffusiveFluxJacobian(eps_rho, */
      /*                                                           ebqe_hhi_ext[ebNE_kb], */
      /*                                                           sdInfo_hv_hv_rowptr, */
      /*                                                           sdInfo_hv_hv_colind, */
      /*                                                           isDOFBoundary_hv[ebNE_kb], */
      /*                                                           normal, */
      /*                                                           mom_hv_diff_ten_ext, */
      /*                                                           vel_trial_trace_ref[ebN_local_kb_j], */
      /*                                                           &vel_grad_trial_trace[j_nSpace], */
      /*                                                           h_henalty);//ebqe_henalty_ext[ebNE_kb]); */

      /*                  fluxJacobian_h_h[j]=0.0; */
      /*                  fluxJacobian_h_hu[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mass_hu_ext,vel_trial_trace_ref[ebN_local_kb_j]); */
      /*                  fluxJacobian_h_hv[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mass_hv_ext,vel_trial_trace_ref[ebN_local_kb_j]); */

      /*                  fluxJacobian_hu_h[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_hu_adv_h_ext,h_trial_trace_ref[ebN_local_kb_j]); */
      /*                  fluxJacobian_hu_hu[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_hu_adv_hu_ext,vel_trial_trace_ref[ebN_local_kb_j]) + */
      /*                    ExteriorNumericalDiffusiveFluxJacobian(eps_rho, */
      /*                                                           ebqe_hhi_ext[ebNE_kb], */
      /*                                                           sdInfo_hu_hu_rowptr, */
      /*                                                           sdInfo_hu_hu_colind, */
      /*                                                           isDOFBoundary_hu[ebNE_kb], */
      /*                                                           normal, */
      /*                                                           mom_hu_diff_ten_ext, */
      /*                                                           vel_trial_trace_ref[ebN_local_kb_j], */
      /*                                                           &vel_grad_trial_trace[j_nSpace], */
      /*                                                           h_henalty);//ebqe_henalty_ext[ebNE_kb]); */
      /*                  fluxJacobian_hu_hv[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_hu_adv_hv_ext,vel_trial_trace_ref[ebN_local_kb_j]); */

      /*                  fluxJacobian_hv_h[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_hv_adv_h_ext,h_trial_trace_ref[ebN_local_kb_j]); */
      /*                  fluxJacobian_hv_hu[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_hv_adv_hu_ext,vel_trial_trace_ref[ebN_local_kb_j]); */
      /*                  fluxJacobian_hv_hv[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_hv_adv_hv_ext,vel_trial_trace_ref[ebN_local_kb_j]) + */
      /*                    ExteriorNumericalDiffusiveFluxJacobian(eps_rho, */
      /*                                                           ebqe_hhi_ext[ebNE_kb], */
      /*                                                           sdInfo_hv_hv_rowptr, */
      /*                                                           sdInfo_hv_hv_colind, */
      /*                                                           isDOFBoundary_hv[ebNE_kb], */
      /*                                                           normal, */
      /*                                                           mom_hv_diff_ten_ext, */
      /*                                                           vel_trial_trace_ref[ebN_local_kb_j], */
      /*                                                           &vel_grad_trial_trace[j_nSpace], */
      /*                                                           h_henalty);//ebqe_henalty_ext[ebNE_kb]); */
      /*                  // //cek debug */
      /*                  // fluxJacobian_h_h[j]=0.0; */
      /*                  // fluxJacobian_h_hu[j]=0.0; */
      /*                  // fluxJacobian_h_hv[j]=0.0; */

      /*                  // fluxJacobian_hu_h[j]=0.0; */
      /*                  // fluxJacobian_hu_hu[j]=0.0; */
      /*                  // fluxJacobian_hu_hv[j]=0.0; */

      /*                  // fluxJacobian_hv_h[j]=0.0; */
      /*                  // fluxJacobian_hv_hu[j]=0.0; */
      /*                  // fluxJacobian_hv_hv[j]=0.0; */

      /*                  // //cek debug */
      /*                }//j */
      /*              // */
      /*              //update the global Jacobian from the flux Jacobian */
      /*              // */
      /*              for (int i=0;i<nDOF_test_element;i++) */
      /*                { */
      /*                  register int eN_i = eN*nDOF_test_element+i; */
      /*                  for (int j=0;j<nDOF_trial_element;j++) */
      /*                    { */
      /*                      register int ebN_i_j = ebN*4*nDOF_test_X_trial_element + i*nDOF_trial_element + j; */

      /*                      globalJacobian[csrRowIndeces_h_h[eN_i] + csrColumnOffsets_eb_h_h[ebN_i_j]] += fluxJacobian_h_h[j]*h_test_dS[i]; */
      /*                      globalJacobian[csrRowIndeces_h_hu[eN_i] + csrColumnOffsets_eb_h_hu[ebN_i_j]] += fluxJacobian_h_hu[j]*h_test_dS[i]; */
      /*                      globalJacobian[csrRowIndeces_h_hv[eN_i] + csrColumnOffsets_eb_h_hv[ebN_i_j]] += fluxJacobian_h_hv[j]*h_test_dS[i]; */

      /*                      globalJacobian[csrRowIndeces_hu_h[eN_i] + csrColumnOffsets_eb_hu_h[ebN_i_j]] += fluxJacobian_hu_h[j]*vel_test_dS[i]; */
      /*                      globalJacobian[csrRowIndeces_hu_hu[eN_i] + csrColumnOffsets_eb_hu_hu[ebN_i_j]] += fluxJacobian_hu_hu[j]*vel_test_dS[i]; */
      /*                      globalJacobian[csrRowIndeces_hu_hv[eN_i] + csrColumnOffsets_eb_hu_hv[ebN_i_j]] += fluxJacobian_hu_hv[j]*vel_test_dS[i]; */

      /*                      globalJacobian[csrRowIndeces_hv_h[eN_i] + csrColumnOffsets_eb_hv_h[ebN_i_j]] += fluxJacobian_hv_h[j]*vel_test_dS[i]; */
      /*                      globalJacobian[csrRowIndeces_hv_hu[eN_i] + csrColumnOffsets_eb_hv_hu[ebN_i_j]] += fluxJacobian_hv_hu[j]*vel_test_dS[i]; */
      /*                      globalJacobian[csrRowIndeces_hv_hv[eN_i] + csrColumnOffsets_eb_hv_hv[ebN_i_j]] += fluxJacobian_hv_hv[j]*vel_test_dS[i]; */

      /*                    }//j */
      /*                }//i */
      /*              // //debug */
      /*              // std::cout<<"flux jacobian ebNE "<<ebNE<<" kb "<<kb<<std::endl; */
      /*              // for (int i=0;i<nDOF_test_element;i++) */
      /*              //   { */
      /*              //     for (int j=0;j<nDOF_trial_element;j++) */
      /*              //        { */
      /*              //          std::cout<< fluxJacobian_h_h[j]*h_test_dS[i]<<std::endl; */
      /*              //          std::cout<< fluxJacobian_h_hu[j]*h_test_dS[i]<<std::endl; */
      /*              //          std::cout<< fluxJacobian_h_hv[j]*h_test_dS[i]<<std::endl; */

      /*              //          std::cout<< fluxJacobian_hu_h[j]*vel_test_dS[i]<<std::endl; */
      /*              //          std::cout<< fluxJacobian_hu_hu[j]*vel_test_dS[i]<<std::endl; */
      /*              //          std::cout<< fluxJacobian_hu_hv[j]*vel_test_dS[i]<<std::endl; */

      /*              //          std::cout<< fluxJacobian_hv_h[j]*vel_test_dS[i]<<std::endl; */
      /*              //          std::cout<< fluxJacobian_hv_hu[j]*vel_test_dS[i]<<std::endl; */
      /*              //          std::cout<< fluxJacobian_hv_hv[j]*vel_test_dS[i]<<std::endl; */
      /*              //        }//j */
      /*              //   }//i */
      /*            }//kb */
      /*        }//ebNE */
    }//computeJacobian

    /* void calculateVelocityAverage(int nExteriorElementBoundaries_global, */
    /*                            int* exteriorElementBoundariesArray, */
    /*                            int nInteriorElementBoundaries_global, */
    /*                            int* interiorElementBoundariesArray, */
    /*                            int* elementBoundaryElementsArray, */
    /*                            int* elementBoundaryLocalElementBoundariesArray, */
    /*                            double* mesh_dof, */
    /*                            int* mesh_l2g, */
    /*                            double* mesh_trial_trace_ref, */
    /*                            double* mesh_grad_trial_trace_ref, */
    /*                            double* normal_ref, */
    /*                            double* boundaryJac_ref, */
    /*                            int* vel_l2g, */
    /*                            double* hu_dof, */
    /*                            double* hv_dof, */
    /*                            double* vel_trial_trace_ref, */
    /*                            double* ebqe_velocity, */
    /*                            double* velocityAverage) */
    /* { */
    /*   int permutations[nQuadraturePoints_elementBoundary]; */
    /*   double xArray_left[nQuadraturePoints_elementBoundary*3], */
    /*  xArray_right[nQuadraturePoints_elementBoundary*3]; */
    /*   for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) */
    /*  { */
    /*    register int ebN = exteriorElementBoundariesArray[ebNE]; */
    /*    for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) */
    /*      { */
    /*        register int ebN_kb_nSpace = ebN*nQuadraturePoints_elementBoundary*nSpace+kb*nSpace, */
    /*          ebNE_kb_nSpace = ebNE*nQuadraturePoints_elementBoundary*nSpace+kb*nSpace; */
    /*        velocityAverage[ebN_kb_nSpace+0]=ebqe_velocity[ebNE_kb_nSpace+0]; */
    /*        velocityAverage[ebN_kb_nSpace+1]=ebqe_velocity[ebNE_kb_nSpace+1]; */
    /*        velocityAverage[ebN_kb_nSpace+2]=ebqe_velocity[ebNE_kb_nSpace+2]; */
    /*      }//ebNE */
    /*  } */
    /*   for (int ebNI = 0; ebNI < nInteriorElementBoundaries_global; ebNI++) */
    /*  { */
    /*    register int ebN = interiorElementBoundariesArray[ebNI], */
    /*      left_eN_global   = elementBoundaryElementsArray[ebN*2+0], */
    /*      left_ebN_element  = elementBoundaryLocalElementBoundariesArray[ebN*2+0], */
    /*      right_eN_global  = elementBoundaryElementsArray[ebN*2+1], */
    /*      right_ebN_element = elementBoundaryLocalElementBoundariesArray[ebN*2+1], */
    /*      left_eN_nDOF_trial_element = left_eN_global*nDOF_trial_element, */
    /*      right_eN_nDOF_trial_element = right_eN_global*nDOF_trial_element; */
    /*    double jac[nSpace*nSpace], */
    /*      jacDet, */
    /*      jacInv[nSpace*nSpace], */
    /*      boundaryJac[nSpace*(nSpace-1)], */
    /*      metricTensor[(nSpace-1)*(nSpace-1)], */
    /*      metricTensorDetSqrt, */
    /*      normal[3], */
    /*      x,y; */
    /*    //double G[nSpace*nSpace],G_dd_G,tr_G,h_hhi,h_henalty; */

    /*    for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) */
    /*      { */
    /*        ck.calculateMapping_elementBoundary(left_eN_global, */
    /*                                            left_ebN_element, */
    /*                                            kb, */
    /*                                            left_ebN_element*kb, */
    /*                                            mesh_dof, */
    /*                                            mesh_l2g, */
    /*                                            mesh_trial_trace_ref, */
    /*                                            mesh_grad_trial_trace_ref, */
    /*                                            boundaryJac_ref, */
    /*                                            jac, */
    /*                                            jacDet, */
    /*                                            jacInv, */
    /*                                            boundaryJac, */
    /*                                            metricTensor, */
    /*                                            metricTensorDetSqrt, */
    /*                                            normal_ref, */
    /*                                            normal, */
    /*                                            x,y); */
    /*        xArray_left[kb*3+0] = x; */
    /*        xArray_left[kb*3+1] = y; */
    /*        ck.calculateMapping_elementBoundary(right_eN_global, */
    /*                                            right_ebN_element, */
    /*                                            kb, */
    /*                                            right_ebN_element*kb, */
    /*                                            mesh_dof, */
    /*                                            mesh_l2g, */
    /*                                            mesh_trial_trace_ref, */
    /*                                            mesh_grad_trial_trace_ref, */
    /*                                            boundaryJac_ref, */
    /*                                            jac, */
    /*                                            jacDet, */
    /*                                            jacInv, */
    /*                                            boundaryJac, */
    /*                                            metricTensor, */
    /*                                            metricTensorDetSqrt, */
    /*                                            normal_ref, */
    /*                                            normal, */
    /*                                            x,y); */
    /*        xArray_right[kb*3+0] = x; */
    /*        xArray_right[kb*3+1] = y; */
    /*      } */
    /*    for  (int kb_left=0;kb_left<nQuadraturePoints_elementBoundary;kb_left++) */
    /*      { */
    /*        double errorNormMin = 1.0; */
    /*        for  (int kb_right=0;kb_right<nQuadraturePoints_elementBoundary;kb_right++) */
    /*          { */
    /*            double errorNorm=0.0; */
    /*            for (int I=0;I<nSpace;I++) */
    /*              { */
    /*                errorNorm += fabs(xArray_left[kb_left*3+I] */
    /*                                  - */
    /*                                  xArray_right[kb_right*3+I]); */
    /*              } */
    /*            if (errorNorm < errorNormMin) */
    /*              { */
    /*                permutations[kb_right] = kb_left; */
    /*                errorNormMin = errorNorm; */
    /*              } */
    /*          } */
    /*      } */
    /*    for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) */
    /*      { */
    /*        register int ebN_kb_nSpace = ebN*nQuadraturePoints_elementBoundary*nSpace+kb*nSpace; */
    /*        register double u_left=0.0, */
    /*          v_left=0.0, */
    /*          u_right=0.0, */
    /*          v_right=0.0; */
    /*        register int left_kb = kb, */
    /*          right_kb = permutations[kb], */
    /*          left_ebN_element_kb_nDOF_test_element=left_ebN_element*left_kb*nDOF_test_element, */
    /*          right_ebN_element_kb_nDOF_test_element=right_ebN_element*right_kb*nDOF_test_element; */
    /*        // */
    /*        //calculate the velocity solution at quadrature points on left and right */
    /*        // */
    /*        ck.valFromDOF(hu_dof,&vel_l2g[left_eN_nDOF_trial_element],&vel_trial_trace_ref[left_ebN_element_kb_nDOF_test_element],u_left); */
    /*        ck.valFromDOF(hv_dof,&vel_l2g[left_eN_nDOF_trial_element],&vel_trial_trace_ref[left_ebN_element_kb_nDOF_test_element],v_left); */
    /*        // */
    /*        ck.valFromDOF(hu_dof,&vel_l2g[right_eN_nDOF_trial_element],&vel_trial_trace_ref[right_ebN_element_kb_nDOF_test_element],u_right); */
    /*        ck.valFromDOF(hv_dof,&vel_l2g[right_eN_nDOF_trial_element],&vel_trial_trace_ref[right_ebN_element_kb_nDOF_test_element],v_right); */
    /*        // */
    /*        velocityAverage[ebN_kb_nSpace+0]=0.5*(u_left + u_right); */
    /*        velocityAverage[ebN_kb_nSpace+1]=0.5*(v_left + v_right); */
    /*      }//ebNI */
    /*  } */
    /* } */

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

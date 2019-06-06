#ifndef RANS3PF2D_H
#define RANS3PF2D_H
#include <cmath>
#include <valarray>
#include <iostream>
#include <vector>
#include <set>
#include <cstring>
#include "CompKernel.h"
#include "ModelFactory.h"
#include "SedClosure.h"
#define DRAG_FAC 1.0
#define TURB_FORCE_FAC 0.0
#define CUT_CELL_INTEGRATION 0
double sgn(double val)
{
  return double((0.0 < val) - (val < 0.0));
}
//////////////////////
// ***** TODO ***** //
//////////////////////
// *fix the following w.r.t. not dividing momentum eqn by rho
//      * updateSolidParticleTerms
// *Double check the following w.r.t. not dividing momentum eqn by rho
//      * updateDarcyForchheimerTerms_Ergun
//      * updateTurbulenceClosure
//      * check pdeResidual_p. In particular check the term with q_dvos_dt
//      * double check exteriorNumericalAdvectiveFlux. I multiply from outside porosity*rho
//      * MOVING MESH. Double check.
//      * Turbulence: double check eddy_viscosity within evaluateCoefficients
// ***** END OF TODO *****

#define CELL_BASED_EV_COEFF 1
#define POWER_SMOOTHNESS_INDICATOR 2
#define EPS_FOR_GAMMA_INDICATOR 1E-10
#define C_FOR_GAMMA_INDICATOR 0.25 // increase gamma to make the indicator more agressive (less dissipative)
#define USE_GAMMA_INDICATOR 0
#define ANISOTROPIC_DIFFUSION 0

inline void baryCoords(const double r0[2],
                       const double r1[2],
                       const double r2[2],
                       const double r[2],
                       double *lambda)
{
  double detT = (r1[1] - r2[1]) * (r0[0] - r2[0]) + (r2[0] - r1[0]) * (r0[1] - r2[1]);
  lambda[0] = ((r1[1] - r2[1]) * (r[0] - r2[0]) + (r2[0] - r1[0]) * (r[1] - r2[1])) / detT;
  lambda[1] = ((r2[1] - r0[1]) * (r[0] - r2[0]) + (r0[0] - r2[0]) * (r[1] - r2[1])) / detT;
  lambda[2] = 1.0 - lambda[0] - lambda[1];
}

namespace proteus
{
class cppRANS3PF2D_base
{

public:
  std::valarray<double> TransportMatrix, TransposeTransportMatrix;
  std::valarray<double> uStar_psi, vStar_psi, wStar_psi;
  std::valarray<double> uStar_hi, vStar_hi, wStar_hi, den_hi;
  std::valarray<double> uStar_min_hiHe, vStar_min_hiHe, wStar_min_hiHe;
  std::valarray<double> uStar_gamma, vStar_gamma, wStar_gamma;
  virtual ~cppRANS3PF2D_base() {}
  virtual void setSedClosure(double aDarcy,
                             double betaForch,
                             double grain,
                             double packFraction,
                             double packMargin,
                             double maxFraction,
                             double frFraction,
                             double sigmaC,
                             double C3e,
                             double C4e,
                             double eR,
                             double fContact,
                             double mContact,
                             double nContact,
                             double angFriction, double vos_limiter, double mu_fr_limiter) {}
  virtual void calculateResidual(double *mesh_trial_ref,
                                 double *mesh_grad_trial_ref,
                                 double *mesh_dof,
                                 double *mesh_velocity_dof,
                                 double MOVING_DOMAIN, //0 or 1
                                 double PSTAB,
                                 int *mesh_l2g,
                                 double *dV_ref,
                                 int nDOF_per_element_pressure,
                                 double *p_trial_ref,
                                 double *p_grad_trial_ref,
                                 double *p_test_ref,
                                 double *p_grad_test_ref,
                                 double *q_ptt_old,/////////////data for pressure
                                 double *q_pt_old,
                                 double *q_p_old,
                                 double *q_ptt,
                                 double *q_pt,
                                 double *q_p,
                                 double *utt_u_dof_old,/////////data for velocity
                                 double *utt_v_dof_old,
                                 double *ut_u_dof_old,
                                 double *ut_v_dof_old,
                                 double *u_u_dof_old,
                                 double *u_v_dof_old,
                                 double *utt_u_dof,
                                 double *utt_v_dof,
                                 double *ut_u_dof,
                                 double *ut_v_dof,
                                 double *u_u_dof,
                                 double *u_v_dof,
                                 double *vel_trial_ref,
                                 double *vel_grad_trial_ref,
                                 double *vel_hess_trial_ref,
                                 double *vel_test_ref,
                                 double *vel_grad_test_ref,
                                 double *mesh_trial_trace_ref,
                                 double *mesh_grad_trial_trace_ref,
                                 double *dS_ref,
                                 double *p_trial_trace_ref,
                                 double *p_grad_trial_trace_ref,
                                 double *p_test_trace_ref,
                                 double *p_grad_test_trace_ref,
                                 double *vel_trial_trace_ref,
                                 double *vel_grad_trial_trace_ref,
                                 double *vel_test_trace_ref,
                                 double *vel_grad_test_trace_ref,
                                 double *normal_ref,
                                 double *boundaryJac_ref,
                                 double eb_adjoint_sigma,
                                 double *elementDiameter,
                                 double *nodeDiametersArray,
                                 double hFactor,
                                 int nElements_global,
                                 int nElements_owned,
                                 int nElementBoundaries_global,
                                 int nElementBoundaries_owned,
                                 int nNodes_owned,
                                 double useRBLES,
                                 double useMetrics,
                                 double alphaBDF,
                                 double epsFact_rho,
                                 double epsFact_mu,
                                 double sigma,
                                 double rho_0,
                                 double nu_0,
                                 double rho_1,
                                 double nu_1,
                                 double smagorinskyConstant,
                                 int turbulenceClosureModel,
                                 double Ct_sge,
                                 double Cd_sge,
                                 double C_dc,
                                 double C_b,
                                 const double *eps_solid,
                                 const double *ebq_global_phi_solid,
                                 const double *ebq_global_grad_phi_solid,
                                 const double *ebq_particle_velocity_solid,
                                 double *phi_solid_nodes,
                                 double *phi_solid,
                                 const double *q_velocity_solid,
                                 const double *q_velocityStar_solid,
                                 const double *q_vos,
                                 const double *q_dvos_dt,
                                 const double *q_grad_vos,
                                 const double *q_dragAlpha,
                                 const double *q_dragBeta,
                                 const double *q_mass_source,
                                 const double *q_turb_var_0,
                                 const double *q_turb_var_1,
                                 const double *q_turb_var_grad_0,
                                 double *q_eddy_viscosity,
                                 int *p_l2g,
                                 int *vel_l2g,
                                 double *p_dof,
                                 double *u_dof,
                                 double *v_dof,
                                 double *w_dof,
                                 double *u_dof_old,
                                 double *v_dof_old,
                                 double *w_dof_old,
                                 double *u_dof_old_old,
                                 double *v_dof_old_old,
                                 double *w_dof_old_old,
                                 double *uStar_dof,
                                 double *vStar_dof,
                                 double *wStar_dof,
                                 double *g,
                                 const double useVF,
                                 double *vf,
                                 double *phi,
                                 double *normal_phi,
                                 double *kappa_phi,
                                 double *q_mom_u_acc,
                                 double *q_mom_v_acc,
                                 double *q_mom_w_acc,
                                 double *q_mass_adv,
                                 double *q_mom_u_acc_beta_bdf,
                                 double *q_mom_v_acc_beta_bdf,
                                 double *q_mom_w_acc_beta_bdf,
                                 double *q_dV,
                                 double *q_dV_last,
                                 double *q_velocity_sge,
                                 double *ebqe_velocity_star,
                                 double *q_cfl,
                                 double *q_numDiff_u,
                                 double *q_numDiff_v,
                                 double *q_numDiff_w,
                                 double *q_numDiff_u_last,
                                 double *q_numDiff_v_last,
                                 double *q_numDiff_w_last,
                                 int *sdInfo_u_u_rowptr,
                                 int *sdInfo_u_u_colind,
                                 int *sdInfo_u_v_rowptr,
                                 int *sdInfo_u_v_colind,
                                 int *sdInfo_u_w_rowptr,
                                 int *sdInfo_u_w_colind,
                                 int *sdInfo_v_v_rowptr,
                                 int *sdInfo_v_v_colind,
                                 int *sdInfo_v_u_rowptr,
                                 int *sdInfo_v_u_colind,
                                 int *sdInfo_v_w_rowptr,
                                 int *sdInfo_v_w_colind,
                                 int *sdInfo_w_w_rowptr,
                                 int *sdInfo_w_w_colind,
                                 int *sdInfo_w_u_rowptr,
                                 int *sdInfo_w_u_colind,
                                 int *sdInfo_w_v_rowptr,
                                 int *sdInfo_w_v_colind,
                                 int offset_p,
                                 int offset_u,
                                 int offset_v,
                                 int offset_w,
                                 int stride_p,
                                 int stride_u,
                                 int stride_v,
                                 int stride_w,
                                 double *globalResidual,
                                 int nExteriorElementBoundaries_global,
                                 int *exteriorElementBoundariesArray,
                                 int *elementBoundariesArray,
                                 int *elementBoundaryElementsArray,
                                 int *elementBoundaryLocalElementBoundariesArray,
                                 double *ebqe_vf_ext,
                                 double *bc_ebqe_vf_ext,
                                 double *ebqe_phi_ext,
                                 double *bc_ebqe_phi_ext,
                                 double *ebqe_normal_phi_ext,
                                 double *ebqe_kappa_phi_ext,
                                 const double *ebqe_vos_ext,
                                 const double *ebqe_turb_var_0,
                                 const double *ebqe_turb_var_1,
                                 int *isDOFBoundary_p,
                                 int *isDOFBoundary_u,
                                 int *isDOFBoundary_v,
                                 int *isDOFBoundary_w,
                                 int *isAdvectiveFluxBoundary_p,
                                 int *isAdvectiveFluxBoundary_u,
                                 int *isAdvectiveFluxBoundary_v,
                                 int *isAdvectiveFluxBoundary_w,
                                 int *isDiffusiveFluxBoundary_u,
                                 int *isDiffusiveFluxBoundary_v,
                                 int *isDiffusiveFluxBoundary_w,
                                 double *ebqe_bc_p_ext,
                                 double *ebqe_bc_flux_mass_ext,
                                 double *ebqe_bc_flux_mom_u_adv_ext,
                                 double *ebqe_bc_flux_mom_v_adv_ext,
                                 double *ebqe_bc_flux_mom_w_adv_ext,
                                 double *ebqe_bc_u_ext,
                                 double *ebqe_bc_flux_u_diff_ext,
                                 double *ebqe_penalty_ext,
                                 double *ebqe_bc_v_ext,
                                 double *ebqe_bc_flux_v_diff_ext,
                                 double *ebqe_bc_w_ext,
                                 double *ebqe_bc_flux_w_diff_ext,
                                 double *q_x,
                                 double *q_velocity,
                                 double *ebqe_velocity,
                                 double *q_grad_u,
                                 double *q_grad_v,
                                 double *q_grad_w,
                                 double *q_divU,
                                 double *ebqe_grad_u,
                                 double *ebqe_grad_v,
                                 double *ebqe_grad_w,
                                 double *flux,
                                 double *elementResidual_p,
                                 int *elementFlags,
                                 int *boundaryFlags,
                                 double *barycenters,
                                 double *wettedAreas,
                                 double *netForces_p,
                                 double *netForces_v,
                                 double *netMoments,
                                 double *q_rho,
                                 double *ebqe_rho,
                                 double *q_nu,
                                 double *ebqe_nu,
                                 int nParticles,
                                 double particle_epsFact,
                                 double particle_alpha,
                                 double particle_beta,
                                 double particle_penalty_constant,
                                 double *particle_signed_distances,
                                 double *particle_signed_distance_normals,
                                 double *particle_velocities,
                                 double *particle_centroids,
                                 double *particle_netForces,
                                 double *particle_netMoments,
                                 double *particle_surfaceArea,
                                 double particle_nitsche,
                                 int use_ball_as_particle,
                                 double *ball_center,
                                 double *ball_radius,
                                 double *ball_velocity,
                                 double *ball_angular_velocity,
                                 double *phisError,
                                 double *phisErrorNodal,
                                 int USE_SUPG,
                                 int ARTIFICIAL_VISCOSITY,
                                 double cMax,
                                 double cE,
                                 int MULTIPLY_EXTERNAL_FORCE_BY_DENSITY,
                                 double *forcex,
                                 double *forcey,
                                 double *forcez,
                                 int KILL_PRESSURE_TERM,
                                 double dt,
                                 double *quantDOFs,
                                 int MATERIAL_PARAMETERS_AS_FUNCTION,
                                 double *density_as_function,
                                 double *dynamic_viscosity_as_function,
                                 double *ebqe_density_as_function,
                                 double *ebqe_dynamic_viscosity_as_function,
                                 double order_polynomial,
                                 double *isActiveDOF,
                                 int USE_SBM,
                                 double *ncDrag,
                                 double *betaDrag,
                                 double *vos_vel_nodes,
                                 // For edge based dissipation
                                 double *entropyResidualPerNode,
                                 double *laggedEntropyResidualPerNode,
                                 double *uStar_dMatrix,
                                 double *vStar_dMatrix,
                                 double *wStar_dMatrix,
                                 int numDOFs_1D,
                                 int NNZ_1D,
                                 int *csrRowIndeces_1D, int *csrColumnOffsets_1D,
                                 int *rowptr_1D, int *colind_1D,
                                 double *isBoundary_1D,
                                 // int by parts pressure
                                 int INT_BY_PARTS_PRESSURE) = 0;
  virtual void calculateJacobian( //element
      double *mesh_trial_ref,
      double *mesh_grad_trial_ref,
      double *mesh_dof,
      double *mesh_velocity_dof,
      double MOVING_DOMAIN,
      double PSTAB,
      int *mesh_l2g,
      double *dV_ref,
      double *p_trial_ref,
      double *p_grad_trial_ref,
      double *p_test_ref,
      double *p_grad_test_ref,
      double *q_p,
      double *q_grad_p,
      double *ebqe_p,
      double *ebqe_grad_p,
      double *vel_trial_ref,
      double *vel_grad_trial_ref,
      double *vel_hess_trial_ref,
      double *vel_test_ref,
      double *vel_grad_test_ref,
      //element boundary
      double *mesh_trial_trace_ref,
      double *mesh_grad_trial_trace_ref,
      double *dS_ref,
      double *p_trial_trace_ref,
      double *p_grad_trial_trace_ref,
      double *p_test_trace_ref,
      double *p_grad_test_trace_ref,
      double *vel_trial_trace_ref,
      double *vel_grad_trial_trace_ref,
      double *vel_test_trace_ref,
      double *vel_grad_test_trace_ref,
      double *normal_ref,
      double *boundaryJac_ref,
      //physics
      double eb_adjoint_sigma,
      double *elementDiameter,
      double *nodeDiametersArray,
      double hFactor,
      int nElements_global,
      int nElements_owned,
      int nElementBoundaries_global,
      int nElementBoundaries_owned,
      int nNodes_owned,
      double useRBLES,
      double useMetrics,
      double alphaBDF,
      double epsFact_rho,
      double epsFact_mu,
      double sigma,
      double rho_0,
      double nu_0,
      double rho_1,
      double nu_1,
      double smagorinskyConstant,
      int turbulenceClosureModel,
      double Ct_sge,
      double Cd_sge,
      double C_dg,
      double C_b,
      //VRANS
      const double *eps_solid,
      const double *ebq_global_phi_solid,
      const double *ebq_global_grad_phi_solid,
      const double *ebq_particle_velocity_solid,
      double *phi_solid_nodes,
      const double *phi_solid,
      const double *q_velocity_solid,
      const double *q_velocityStar_solid,
      const double *q_vos,
      const double *q_dvos_dt,
      const double *q_grad_vos,
      const double *q_dragAlpha,
      const double *q_dragBeta,
      const double *q_mass_source,
      const double *q_turb_var_0,
      const double *q_turb_var_1,
      const double *q_turb_var_grad_0,
      int *p_l2g,
      int *vel_l2g,
      double *p_dof, double *u_dof, double *v_dof, double *w_dof,
      double *g,
      const double useVF,
      double *vf,
      double *phi,
      double *normal_phi,
      double *kappa_phi,
      double *q_mom_u_acc_beta_bdf, double *q_mom_v_acc_beta_bdf, double *q_mom_w_acc_beta_bdf,
      double *q_dV,
      double *q_dV_last,
      double *q_velocity_sge,
      double *ebqe_velocity_star,
      double *q_cfl,
      double *q_numDiff_u_last, double *q_numDiff_v_last, double *q_numDiff_w_last,
      int *sdInfo_u_u_rowptr, int *sdInfo_u_u_colind,
      int *sdInfo_u_v_rowptr, int *sdInfo_u_v_colind,
      int *sdInfo_u_w_rowptr, int *sdInfo_u_w_colind,
      int *sdInfo_v_v_rowptr, int *sdInfo_v_v_colind,
      int *sdInfo_v_u_rowptr, int *sdInfo_v_u_colind,
      int *sdInfo_v_w_rowptr, int *sdInfo_v_w_colind,
      int *sdInfo_w_w_rowptr, int *sdInfo_w_w_colind,
      int *sdInfo_w_u_rowptr, int *sdInfo_w_u_colind,
      int *sdInfo_w_v_rowptr, int *sdInfo_w_v_colind,
      int *csrRowIndeces_p_p, int *csrColumnOffsets_p_p,
      int *csrRowIndeces_p_u, int *csrColumnOffsets_p_u,
      int *csrRowIndeces_p_v, int *csrColumnOffsets_p_v,
      int *csrRowIndeces_p_w, int *csrColumnOffsets_p_w,
      int *csrRowIndeces_u_p, int *csrColumnOffsets_u_p,
      int *csrRowIndeces_u_u, int *csrColumnOffsets_u_u,
      int *csrRowIndeces_u_v, int *csrColumnOffsets_u_v,
      int *csrRowIndeces_u_w, int *csrColumnOffsets_u_w,
      int *csrRowIndeces_v_p, int *csrColumnOffsets_v_p,
      int *csrRowIndeces_v_u, int *csrColumnOffsets_v_u,
      int *csrRowIndeces_v_v, int *csrColumnOffsets_v_v,
      int *csrRowIndeces_v_w, int *csrColumnOffsets_v_w,
      int *csrRowIndeces_w_p, int *csrColumnOffsets_w_p,
      int *csrRowIndeces_w_u, int *csrColumnOffsets_w_u,
      int *csrRowIndeces_w_v, int *csrColumnOffsets_w_v,
      int *csrRowIndeces_w_w, int *csrColumnOffsets_w_w,
      double *globalJacobian,
      int nExteriorElementBoundaries_global,
      int *exteriorElementBoundariesArray,
      int *elementBoundariesArray,
      int *elementBoundaryElementsArray,
      int *elementBoundaryLocalElementBoundariesArray,
      double *ebqe_vf_ext,
      double *bc_ebqe_vf_ext,
      double *ebqe_phi_ext,
      double *bc_ebqe_phi_ext,
      double *ebqe_normal_phi_ext,
      double *ebqe_kappa_phi_ext,
      //VRANS
      const double *ebqe_vos_ext,
      const double *ebqe_turb_var_0,
      const double *ebqe_turb_var_1,
      //VRANS end
      int *isDOFBoundary_p,
      int *isDOFBoundary_u,
      int *isDOFBoundary_v,
      int *isDOFBoundary_w,
      int *isAdvectiveFluxBoundary_p,
      int *isAdvectiveFluxBoundary_u,
      int *isAdvectiveFluxBoundary_v,
      int *isAdvectiveFluxBoundary_w,
      int *isDiffusiveFluxBoundary_u,
      int *isDiffusiveFluxBoundary_v,
      int *isDiffusiveFluxBoundary_w,
      double *ebqe_bc_p_ext,
      double *ebqe_bc_flux_mass_ext,
      double *ebqe_bc_flux_mom_u_adv_ext,
      double *ebqe_bc_flux_mom_v_adv_ext,
      double *ebqe_bc_flux_mom_w_adv_ext,
      double *ebqe_bc_u_ext,
      double *ebqe_bc_flux_u_diff_ext,
      double *ebqe_penalty_ext,
      double *ebqe_bc_v_ext,
      double *ebqe_bc_flux_v_diff_ext,
      double *ebqe_bc_w_ext,
      double *ebqe_bc_flux_w_diff_ext,
      int *csrColumnOffsets_eb_p_p,
      int *csrColumnOffsets_eb_p_u,
      int *csrColumnOffsets_eb_p_v,
      int *csrColumnOffsets_eb_p_w,
      int *csrColumnOffsets_eb_u_p,
      int *csrColumnOffsets_eb_u_u,
      int *csrColumnOffsets_eb_u_v,
      int *csrColumnOffsets_eb_u_w,
      int *csrColumnOffsets_eb_v_p,
      int *csrColumnOffsets_eb_v_u,
      int *csrColumnOffsets_eb_v_v,
      int *csrColumnOffsets_eb_v_w,
      int *csrColumnOffsets_eb_w_p,
      int *csrColumnOffsets_eb_w_u,
      int *csrColumnOffsets_eb_w_v,
      int *csrColumnOffsets_eb_w_w,
      int *elementFlags,
      int nParticles,
      double particle_epsFact,
      double particle_alpha,
      double particle_beta,
      double particle_penalty_constant,
      double *particle_signed_distances,
      double *particle_signed_distance_normals,
      double *particle_velocities,
      double *particle_centroids,
      double particle_nitsche,
      int use_ball_as_particle,
      double *ball_center,
      double *ball_radius,
      double *ball_velocity,
      double *ball_angular_velocity,
      int USE_SUPG,
      int KILL_PRESSURE_TERM,
      double dt,
      int MATERIAL_PARAMETERS_AS_FUNCTION,
      double *density_as_function,
      double *dynamic_viscosity_as_function,
      double *ebqe_density_as_function,
      double *ebqe_dynamic_viscosity_as_function,
      int USE_SBM,
      // For edge based dissipation
      int ARTIFICIAL_VISCOSITY,
      double *uStar_dMatrix,
      double *vStar_dMatrix,
      double *wStar_dMatrix,
      int numDOFs_1D,
      int offset_u, int offset_v, int offset_w,
      int stride_u, int stride_v, int stride_w,
      int *rowptr_1D, int *colind_1D,
      int *rowptr, int *colind,
      int INT_BY_PARTS_PRESSURE) = 0;
  virtual void calculateVelocityAverage(int nExteriorElementBoundaries_global,
                                        int *exteriorElementBoundariesArray,
                                        int nInteriorElementBoundaries_global,
                                        int *interiorElementBoundariesArray,
                                        int *elementBoundaryElementsArray,
                                        int *elementBoundaryLocalElementBoundariesArray,
                                        double *mesh_dof,
                                        double *mesh_velocity_dof,
                                        double MOVING_DOMAIN, //0 or 1
                                        int *mesh_l2g,
                                        double *mesh_trial_trace_ref,
                                        double *mesh_grad_trial_trace_ref,
                                        double *normal_ref,
                                        double *boundaryJac_ref,
                                        int *vel_l2g,
                                        double *u_dof,
                                        double *v_dof,
                                        double *w_dof,
                                        double *vos_dof,
                                        double *vel_trial_trace_ref,
                                        double *ebqe_velocity,
                                        double *velocityAverage) = 0;
  virtual void getBoundaryDOFs(double *mesh_dof,
                               int *mesh_l2g,
                               double *mesh_trial_trace_ref,
                               double *mesh_grad_trial_trace_ref,
                               double *dS_ref,
                               double *vel_test_trace_ref,
                               double *normal_ref,
                               double *boundaryJac_ref,
                               int *vel_l2g,
                               int nExteriorElementBoundaries_global,
                               int *exteriorElementBoundariesArray,
                               int *elementBoundaryElementsArray,
                               int *elementBoundaryLocalElementBoundariesArray,
                               double *isBoundary_1D) = 0;
};

template <class CompKernelType,
          int nSpace,
          int nQuadraturePoints_element,
          int nDOF_mesh_trial_element,
          int nDOF_trial_element,
          int nDOF_test_element,
          int nQuadraturePoints_elementBoundary>
class cppRANS3PF2D : public cppRANS3PF2D_base
{
public:
  std::vector<int> surrogate_boundaries, surrogate_boundary_elements, surrogate_boundary_particle;
  std::valarray<double> TransportMatrix, TransposeTransportMatrix, psi;
  double C_sbm, beta_sbm;
  cppHsuSedStress<2> closure;
  const int nDOF_test_X_trial_element,
      nSpace2;
  CompKernelType ck;
  cppRANS3PF2D() : nSpace2(4),
                   closure(150.0,
                           0.0,
                           0.0102,
                           0.2,
                           0.01,
                           0.635,
                           0.57,
                           1.1,
                           1.2,
                           1.0,
                           0.8,
                           0.02,
                           2.0,
                           5.0,
                           M_PI / 6., 0.05, 1.00),
                   nDOF_test_X_trial_element(nDOF_test_element * nDOF_trial_element),
                   ck(),
                   C_sbm(10.0),
                   beta_sbm(0.0)
  { /*        std::cout<<"Constructing cppRANS3PF2D<CompKernelTemplate<"
                     <<0<<","
                     <<0<<","
                     <<0<<","
                     <<0<<">,"*/
    /*  <<nSpaceIn<<","
                <<nQuadraturePoints_elementIn<<","
                <<nDOF_mesh_trial_elementIn<<","
                <<nDOF_trial_elementIn<<","
                <<nDOF_test_elementIn<<","
                <<nQuadraturePoints_elementBoundaryIn<<">());"*/
    /*  <<std::endl<<std::flush; */
  }

  void setSedClosure(double aDarcy,
                     double betaForch,
                     double grain,
                     double packFraction,
                     double packMargin,
                     double maxFraction,
                     double frFraction,
                     double sigmaC,
                     double C3e,
                     double C4e,
                     double eR,
                     double fContact,
                     double mContact,
                     double nContact,
                     double angFriction, double vos_limiter, double mu_fr_limiter)
  {
    closure = cppHsuSedStress<2>(aDarcy,
                                 betaForch,
                                 grain,
                                 packFraction,
                                 packMargin,
                                 maxFraction,
                                 frFraction,
                                 sigmaC,
                                 C3e,
                                 C4e,
                                 eR,
                                 fContact,
                                 mContact,
                                 nContact,
                                 angFriction, vos_limiter, mu_fr_limiter);
  }

  inline double Dot(const double vec1[nSpace],
                    const double vec2[nSpace])
  {
    double dot = 0;
    for (int I = 0; I < nSpace; I++)
      dot += vec1[I] * vec2[I];
    return dot;
  }

  inline void calculateTangentialGradient(const double normal[nSpace],
                                          const double vel_grad[nSpace],
                                          double vel_tgrad[nSpace])
  {
    double normal_dot_vel_grad = Dot(normal, vel_grad);
    for (int I = 0; I < nSpace; I++)
      vel_tgrad[I] = vel_grad[I] - normal_dot_vel_grad * normal[I];
  }

  inline double smoothedHeaviside(double eps, double phi)
  {
    double H;
    if (phi > eps)
      H = 1.0;
    else if (phi < -eps)
      H = 0.0;
    else if (phi == 0.0)
      H = 0.5;
    else
      H = 0.5 * (1.0 + phi / eps + sin(M_PI * phi / eps) / M_PI);
    return H;
  }

  inline double smoothedHeaviside_integral(double eps, double phi)
  {
    double HI;
    if (phi > eps)
    {
      HI = phi - eps + 0.5 * (eps + 0.5 * eps * eps / eps - eps * cos(M_PI * eps / eps) / (M_PI * M_PI)) - 0.5 * ((-eps) + 0.5 * (-eps) * (-eps) / eps - eps * cos(M_PI * (-eps) / eps) / (M_PI * M_PI));
    }
    else if (phi < -eps)
    {
      HI = 0.0;
    }
    else
    {
      HI = 0.5 * (phi + 0.5 * phi * phi / eps - eps * cos(M_PI * phi / eps) / (M_PI * M_PI)) - 0.5 * ((-eps) + 0.5 * (-eps) * (-eps) / eps - eps * cos(M_PI * (-eps) / eps) / (M_PI * M_PI));
    }
    return HI;
  }

  inline double smoothedDirac(double eps, double phi)
  {
    double d;
    if (phi > eps)
      d = 0.0;
    else if (phi < -eps)
      d = 0.0;
    else
      d = 0.5 * (1.0 + cos(M_PI * phi / eps)) / eps;
    return d;
  }
  inline void evaluateCoefficients(const double eps_rho,
                                   const double eps_mu,
                                   const double eps_s,
                                   const double sigma,
                                   const double rho_0,
                                   double nu_0,
                                   const double rho_1,
                                   double nu_1,
                                   const double h_e,
                                   const double smagorinskyConstant,
                                   const int turbulenceClosureModel,
                                   const double g[nSpace],
                                   const double useVF,
                                   const double &vf,
                                   const double &phi,
                                   const double n[nSpace],
                                   const double distance_to_omega_solid,
                                   const double &kappa,
                                   const double porosity, //VRANS specific
                                   const double &p,
                                   const double grad_p[nSpace],
                                   const double grad_u[nSpace],
                                   const double grad_v[nSpace],
                                   const double grad_w[nSpace],
                                   const double &u,
                                   const double &v,
                                   const double &w,
                                   const double &uStar,
                                   const double &vStar,
                                   const double &wStar,
                                   double &eddy_viscosity,
                                   double &mom_u_acc,
                                   double &dmom_u_acc_u,
                                   double &mom_v_acc,
                                   double &dmom_v_acc_v,
                                   double &mom_w_acc,
                                   double &dmom_w_acc_w,
                                   double mass_adv[nSpace],
                                   double dmass_adv_u[nSpace],
                                   double dmass_adv_v[nSpace],
                                   double dmass_adv_w[nSpace],
                                   double mom_u_adv[nSpace],
                                   double dmom_u_adv_u[nSpace],
                                   double dmom_u_adv_v[nSpace],
                                   double dmom_u_adv_w[nSpace],
                                   double mom_v_adv[nSpace],
                                   double dmom_v_adv_u[nSpace],
                                   double dmom_v_adv_v[nSpace],
                                   double dmom_v_adv_w[nSpace],
                                   double mom_w_adv[nSpace],
                                   double dmom_w_adv_u[nSpace],
                                   double dmom_w_adv_v[nSpace],
                                   double dmom_w_adv_w[nSpace],
                                   double mom_uu_diff_ten[nSpace],
                                   double mom_vv_diff_ten[nSpace],
                                   double mom_ww_diff_ten[nSpace],
                                   double mom_uv_diff_ten[1],
                                   double mom_uw_diff_ten[1],
                                   double mom_vu_diff_ten[1],
                                   double mom_vw_diff_ten[1],
                                   double mom_wu_diff_ten[1],
                                   double mom_wv_diff_ten[1],
                                   double &mom_u_source,
                                   double &mom_v_source,
                                   double &mom_w_source,
                                   double &mom_u_ham,
                                   double dmom_u_ham_grad_p[nSpace],
                                   double dmom_u_ham_grad_u[nSpace],
                                   double &mom_v_ham,
                                   double dmom_v_ham_grad_p[nSpace],
                                   double dmom_v_ham_grad_v[nSpace],
                                   double &mom_w_ham,
                                   double dmom_w_ham_grad_p[nSpace],
                                   double dmom_w_ham_grad_w[nSpace],
                                   double &rhoSave,
                                   double &nuSave,
                                   int KILL_PRESSURE_TERM,
                                   int MULTIPLY_EXTERNAL_FORCE_BY_DENSITY,
                                   double forcex,
                                   double forcey,
                                   double forcez,
                                   int MATERIAL_PARAMETERS_AS_FUNCTION,
                                   double density_as_function,
                                   double dynamic_viscosity_as_function,
                                   int USE_SBM,
                                   double x, double y, double z,
                                   int use_ball_as_particle,
                                   double *ball_center,
                                   double *ball_radius,
                                   double *ball_velocity,
                                   double *ball_angular_velocity,
                                   // int by parts pressure
                                   int INT_BY_PARTS_PRESSURE)
  {
    double rho, nu, mu, H_rho, d_rho, H_mu, d_mu, norm_n, nu_t0 = 0.0, nu_t1 = 0.0, nu_t;
    H_rho = (1.0 - useVF) * smoothedHeaviside(eps_rho, phi) + useVF * fmin(1.0, fmax(0.0, vf));
    d_rho = (1.0 - useVF) * smoothedDirac(eps_rho, phi);
    H_mu = (1.0 - useVF) * smoothedHeaviside(eps_mu, phi) + useVF * fmin(1.0, fmax(0.0, vf));
    d_mu = (1.0 - useVF) * smoothedDirac(eps_mu, phi);

    rho = rho_0;
    nu  = nu_0;
    mu = nu*rho;

    rhoSave = rho;
    nuSave = nu;

    //u momentum accumulation
    mom_u_acc = u; //trick for non-conservative form
    dmom_u_acc_u = rho;

    //v momentum accumulation
    mom_v_acc = v;
    dmom_v_acc_v = rho;

    //mass advective flux
    mass_adv[0] = u;
    mass_adv[1] = v;
    
    dmass_adv_u[0] = 1.0;
    dmass_adv_u[1] = 0.0;

    dmass_adv_v[0] = 0.0;
    dmass_adv_v[1] = 1.0;
    
    //u momentum advective flux
    mom_u_adv[0] = 0.0;
    mom_u_adv[1] = 0.0;
    
    dmom_u_adv_u[0] = 0.0;
    dmom_u_adv_u[1] = 0.0;

    dmom_u_adv_v[0] = 0.0;
    dmom_u_adv_v[1] = 0.0;

    //v momentum advective_flux
    mom_v_adv[0] = 0.0;
    mom_v_adv[1] = 0.0;

    dmom_v_adv_u[0] = 0.0;
    dmom_v_adv_u[1] = 0.0;

    dmom_v_adv_v[0] = 0.0;
    dmom_v_adv_v[1] = 0.0;

    if(1){///for symmetric gradient
      //u momentum diffusion tensor
      mom_uu_diff_ten[0] = 2.0 * mu;
      mom_uu_diff_ten[1] = mu;
      mom_uv_diff_ten[0] = mu;

      //v momentum diffusion tensor
      mom_vv_diff_ten[0] = mu;
      mom_vv_diff_ten[1] = 2.0 * mu;
      mom_vu_diff_ten[0] = mu;
    }
    else
    {
      //u momentum diffusion tensor
      mom_uu_diff_ten[0] = mu;
      mom_uu_diff_ten[1] = mu;
      
      //v momentum diffusion tensor
      mom_vv_diff_ten[0] = mu;
      mom_vv_diff_ten[1] = mu;
    }
    //momentum sources
    mom_u_source = -rho * g[0];
    mom_v_source = -rho * g[1];

    //add general force term
    mom_u_source -= forcex;
    mom_v_source -= forcey;
    
    if(0)
    {
    // u momentum Hamiltonian (pressure)
    double aux_pressure = 0.0;
    mom_u_ham = grad_p[0] * aux_pressure;
    dmom_u_ham_grad_p[0] = aux_pressure;
    dmom_u_ham_grad_p[1] = 0.0;
    //v momentum Hamiltonian (pressure)
    mom_v_ham = grad_p[1] * aux_pressure;
    dmom_v_ham_grad_p[0] = 0.0;
    dmom_v_ham_grad_p[1] = aux_pressure;

    //u momentum Hamiltonian (advection)
    mom_u_ham += rho * (uStar * grad_u[0] + vStar * grad_u[1]);
    dmom_u_ham_grad_u[0] = rho * uStar;
    dmom_u_ham_grad_u[1] = rho * vStar;
    /* dmom_u_ham_grad_u[2]=porosity*rho*wStar; */

    //v momentum Hamiltonian (advection)
    mom_v_ham += rho * (uStar * grad_v[0] + vStar * grad_v[1]);
    dmom_v_ham_grad_v[0] = rho * uStar;
    dmom_v_ham_grad_v[1] = rho * vStar;
  }
  else
  {
    mom_u_ham = 0.0;
    dmom_u_ham_grad_p[0] = 0.0;
    dmom_u_ham_grad_p[1] = 0.0;
    mom_v_ham = 0.0;
    dmom_v_ham_grad_p[0] = 0.0;
    dmom_v_ham_grad_p[1] = 0.0;
    dmom_u_ham_grad_u[0] = 0.0;
    dmom_u_ham_grad_u[1] = 0.0;
    dmom_v_ham_grad_v[0] = 0.0;
    dmom_v_ham_grad_v[1] = 0.0;
  }
  }

  //VRANS specific
  inline void updateDarcyForchheimerTerms_Ergun(/* const double linearDragFactor, */
                                                /* const double nonlinearDragFactor, */
                                                /* const double porosity, */
                                                /* const double meanGrainSize, */
                                                const double alpha,
                                                const double beta,
                                                const double eps_rho,
                                                const double eps_mu,
                                                const double rho_0,
                                                const double nu_0,
                                                const double rho_1,
                                                const double nu_1,
                                                double nu_t,
                                                const double useVF,
                                                const double vf,
                                                const double phi,
                                                const double u,
                                                const double v,
                                                const double w,
                                                const double uStar,
                                                const double vStar,
                                                const double wStar,
                                                const double eps_s,
                                                const double phi_s,
                                                const double u_s,
                                                const double v_s,
                                                const double w_s,
                                                const double uStar_s,
                                                const double vStar_s,
                                                const double wStar_s,
                                                double &mom_u_source,
                                                double &mom_v_source,
                                                double &mom_w_source,
                                                double dmom_u_source[nSpace],
                                                double dmom_v_source[nSpace],
                                                double dmom_w_source[nSpace],
                                                double gradC_x,
                                                double gradC_y,
                                                double gradC_z)
  {
  }

  inline void updateSolidParticleTerms(bool element_owned,
                                       const double particle_nitsche,
                                       const double dV,
                                       const int nParticles,
                                       const int sd_offset,
                                       double *particle_signed_distances,
                                       double *particle_signed_distance_normals,
                                       double *particle_velocities,
                                       double *particle_centroids,
                                       int use_ball_as_particle,
                                       double *ball_center,
                                       double *ball_radius,
                                       double *ball_velocity,
                                       double *ball_angular_velocity,
                                       const double porosity, //VRANS specific
                                       const double penalty,
                                       const double alpha,
                                       const double beta,
                                       const double eps_rho,
                                       const double eps_mu,
                                       const double rho_0,
                                       const double nu_0,
                                       const double rho_1,
                                       const double nu_1,
                                       const double useVF,
                                       const double vf,
                                       const double phi,
                                       const double x,
                                       const double y,
                                       const double z,
                                       const double p,
                                       const double u,
                                       const double v,
                                       const double w,
                                       const double uStar,
                                       const double vStar,
                                       const double wStar,
                                       const double eps_s,
                                       const double grad_u[nSpace],
                                       const double grad_v[nSpace],
                                       const double grad_w[nSpace],
                                       double &mom_u_source,
                                       double &mom_v_source,
                                       double &mom_w_source,
                                       double dmom_u_source[nSpace],
                                       double dmom_v_source[nSpace],
                                       double dmom_w_source[nSpace],
                                       double mom_u_adv[nSpace],
                                       double mom_v_adv[nSpace],
                                       double mom_w_adv[nSpace],
                                       double dmom_u_adv_u[nSpace],
                                       double dmom_v_adv_v[nSpace],
                                       double dmom_w_adv_w[nSpace],
                                       double &mom_u_ham,
                                       double dmom_u_ham_grad_u[nSpace],
                                       double &mom_v_ham,
                                       double dmom_v_ham_grad_v[nSpace],
                                       double &mom_w_ham,
                                       double dmom_w_ham_grad_w[nSpace],
                                       double *particle_netForces,
                                       double *particle_netMoments,
                                       double *particle_surfaceArea)
  {

  }
  inline void compute_force_around_solid(bool element_owned,
                                         const double dV,
                                         const int nParticles,
                                         const int sd_offset,
                                         double *particle_signed_distances,
                                         double *particle_signed_distance_normals,
                                         double *particle_velocities,
                                         double *particle_centroids,
                                         int use_ball_as_particle,
                                         double *ball_center,
                                         double *ball_radius,
                                         double *ball_velocity,
                                         double *ball_angular_velocity,
                                         const double penalty,
                                         const double alpha,
                                         const double beta,
                                         const double eps_rho,
                                         const double eps_mu,
                                         const double rho_0,
                                         const double nu_0,
                                         const double rho_1,
                                         const double nu_1,
                                         const double useVF,
                                         const double vf,
                                         const double phi,
                                         const double x,
                                         const double y,
                                         const double z,
                                         const double p,
                                         const double u,
                                         const double v,
                                         const double w,
                                         const double uStar,
                                         const double vStar,
                                         const double wStar,
                                         const double eps_s,
                                         const double grad_u[nSpace],
                                         const double grad_v[nSpace],
                                         const double grad_w[nSpace],
                                         double *particle_netForces,
                                         double *particle_netMoments)
  {
    double C, rho, mu, nu, H_mu, uc, duc_du, duc_dv, duc_dw, H_s, D_s, phi_s, u_s, v_s, w_s, force_x, force_y, r_x, r_y;
    double phi_s_normal[2];
    double fluid_outward_normal[2];
    double vel[2];
    double center[2];
    H_mu = (1.0 - useVF) * smoothedHeaviside(eps_mu, phi) + useVF * fmin(1.0, fmax(0.0, vf));
    nu = nu_0 * (1.0 - H_mu) + nu_1 * H_mu;
    rho = rho_0 * (1.0 - H_mu) + rho_1 * H_mu;
    mu = rho_0 * nu_0 * (1.0 - H_mu) + rho_1 * nu_1 * H_mu;
    C = 0.0;
    for (int i = 0; i < nParticles; i++)
    {
      if (use_ball_as_particle == 1)
      {
        get_distance_to_ith_ball(nParticles, ball_center, ball_radius, i, x, y, z, phi_s);
        get_normal_to_ith_ball(nParticles, ball_center, ball_radius, i, x, y, z, phi_s_normal[0], phi_s_normal[1]);
        get_velocity_to_ith_ball(nParticles, ball_center, ball_radius,
                                 ball_velocity, ball_angular_velocity,
                                 i, x, y, z,
                                 vel[0], vel[1]);
        center[0] = ball_center[3 * i + 0];
        center[1] = ball_center[3 * i + 1];
      }
      else
      {
        phi_s = particle_signed_distances[i * sd_offset];
        phi_s_normal[0] = particle_signed_distance_normals[i * sd_offset * 3 + 0];
        phi_s_normal[1] = particle_signed_distance_normals[i * sd_offset * 3 + 1];
        vel[0] = particle_velocities[i * sd_offset * 3 + 0];
        vel[1] = particle_velocities[i * sd_offset * 3 + 1];
        center[0] = particle_centroids[3 * i + 0];
        center[1] = particle_centroids[3 * i + 1];
      }
      fluid_outward_normal[0] = -phi_s_normal[0];
      fluid_outward_normal[1] = -phi_s_normal[1];
      u_s = vel[0];
      v_s = vel[1];
      w_s = 0;
      H_s = smoothedHeaviside(eps_s, phi_s);
      D_s = smoothedDirac(eps_s, phi_s);
      double rel_vel_norm = sqrt((uStar - u_s) * (uStar - u_s) + (vStar - v_s) * (vStar - v_s) + (wStar - w_s) * (wStar - w_s));
      double C_surf = (phi_s > 0.0) ? 0.0 : nu * penalty;
      double C_vol = (phi_s > 0.0) ? 0.0 : (alpha + beta * rel_vel_norm);

      C = (D_s * C_surf + (1.0 - H_s) * C_vol);
      force_x = dV * D_s * (p * fluid_outward_normal[0] - mu * (fluid_outward_normal[0] * 2 * grad_u[0] + fluid_outward_normal[1] * (grad_u[1] + grad_v[0])));
      //+dV*D_s*C_surf*rel_vel_norm*(u-u_s)*rho
      //+dV * (1.0 - H_s) * C_vol * (u - u_s) * rho;
      force_y = dV * D_s * (p * fluid_outward_normal[1] - mu * (fluid_outward_normal[0] * (grad_u[1] + grad_v[0]) + fluid_outward_normal[1] * 2 * grad_v[1]));
      //+dV*D_s*C_surf*rel_vel_norm*(v-v_s)*rho
      //+dV * (1.0 - H_s) * C_vol * (v - v_s) * rho;

      //always 3D for particle centroids
      r_x = x - center[0];
      r_y = y - center[1];

      if (element_owned)
      {
        particle_netForces[i * 3 + 0] += force_x;
        particle_netForces[i * 3 + 1] += force_y;
        particle_netMoments[i * 3 + 2] += (r_x * force_y - r_y * force_x);
      }
    }
  }
  inline void calculateCFL(const double &hFactor,
                           const double &elementDiameter,
                           const double &dm,
                           const double df[nSpace],
                           double &cfl)
  {
  }

  inline void updateTurbulenceClosure(const int turbulenceClosureModel,
                                      const double eps_rho,
                                      const double eps_mu,
                                      const double rho_0,
                                      const double nu_0,
                                      const double rho_1,
                                      const double nu_1,
                                      const double useVF,
                                      const double vf,
                                      const double phi,
                                      const double porosity,
                                      const double eddy_visc_coef_0,
                                      const double turb_var_0,          //k for k-eps or k-omega
                                      const double turb_var_1,          //epsilon for k-epsilon, omega for k-omega
                                      const double turb_grad_0[nSpace], //grad k for k-eps,k-omega
                                      double &eddy_viscosity,
                                      double mom_uu_diff_ten[nSpace],
                                      double mom_vv_diff_ten[nSpace],
                                      double mom_ww_diff_ten[nSpace],
                                      double mom_uv_diff_ten[1],
                                      double mom_uw_diff_ten[1],
                                      double mom_vu_diff_ten[1],
                                      double mom_vw_diff_ten[1],
                                      double mom_wu_diff_ten[1],
                                      double mom_wv_diff_ten[1],
                                      double &mom_u_source,
                                      double &mom_v_source,
                                      double &mom_w_source)
  {
  }

  inline void calculateSubgridError_tau(const double &hFactor,
                                        const double &elementDiameter,
                                        const double &dmt,
                                        const double &dm,
                                        const double df[nSpace],
                                        const double &a,
                                        const double &pfac,
                                        double &tau_v,
                                        double &tau_p,
                                        double &cfl)
  {
  }

  inline void calculateSubgridError_tau(const double &Ct_sge,
                                        const double &Cd_sge,
                                        const double G[nSpace * nSpace],
                                        const double &G_dd_G,
                                        const double &tr_G,
                                        const double &A0,
                                        const double Ai[nSpace],
                                        const double &Kij,
                                        const double &pfac,
                                        double &tau_v,
                                        double &tau_p,
                                        double &q_cfl)
  {
  }

  inline void calculateSubgridError_tauRes(const double &tau_p,
                                           const double &tau_v,
                                           const double &pdeResidualP,
                                           const double &pdeResidualU,
                                           const double &pdeResidualV,
                                           const double &pdeResidualW,
                                           double &subgridErrorP,
                                           double &subgridErrorU,
                                           double &subgridErrorV,
                                           double &subgridErrorW)
  {
  }

  inline void calculateSubgridErrorDerivatives_tauRes(const double &tau_p,
                                                      const double &tau_v,
                                                      const double dpdeResidualP_du[nDOF_trial_element],
                                                      const double dpdeResidualP_dv[nDOF_trial_element],
                                                      const double dpdeResidualP_dw[nDOF_trial_element],
                                                      const double dpdeResidualU_dp[nDOF_trial_element],
                                                      const double dpdeResidualU_du[nDOF_trial_element],
                                                      const double dpdeResidualV_dp[nDOF_trial_element],
                                                      const double dpdeResidualV_dv[nDOF_trial_element],
                                                      const double dpdeResidualW_dp[nDOF_trial_element],
                                                      const double dpdeResidualW_dw[nDOF_trial_element],
                                                      double dsubgridErrorP_du[nDOF_trial_element],
                                                      double dsubgridErrorP_dv[nDOF_trial_element],
                                                      double dsubgridErrorP_dw[nDOF_trial_element],
                                                      double dsubgridErrorU_dp[nDOF_trial_element],
                                                      double dsubgridErrorU_du[nDOF_trial_element],
                                                      double dsubgridErrorV_dp[nDOF_trial_element],
                                                      double dsubgridErrorV_dv[nDOF_trial_element],
                                                      double dsubgridErrorW_dp[nDOF_trial_element],
                                                      double dsubgridErrorW_dw[nDOF_trial_element])
  {
  }

  inline void exteriorNumericalAdvectiveFlux(const int &isDOFBoundary_p,
                                             const int &isDOFBoundary_u,
                                             const int &isDOFBoundary_v,
                                             const int &isDOFBoundary_w,
                                             const int &isFluxBoundary_p,
                                             const int &isFluxBoundary_u,
                                             const int &isFluxBoundary_v,
                                             const int &isFluxBoundary_w,
                                             const double &oneByRho,
                                             const double &bc_oneByRho,
                                             const double n[nSpace],
                                             const double &porosity, //mql. CHECK. Multiply by rho outside
                                             const double &bc_p,
                                             const double &bc_u,
                                             const double &bc_v,
                                             const double &bc_w,
                                             const double bc_f_mass[nSpace],
                                             const double bc_f_umom[nSpace],
                                             const double bc_f_vmom[nSpace],
                                             const double bc_f_wmom[nSpace],
                                             const double &bc_flux_mass,
                                             const double &bc_flux_umom,
                                             const double &bc_flux_vmom,
                                             const double &bc_flux_wmom,
                                             const double &p,
                                             const double &u,
                                             const double &v,
                                             const double &w,
                                             const double f_mass[nSpace],
                                             const double f_umom[nSpace],
                                             const double f_vmom[nSpace],
                                             const double f_wmom[nSpace],
                                             const double df_mass_du[nSpace],
                                             const double df_mass_dv[nSpace],
                                             const double df_mass_dw[nSpace],
                                             const double df_umom_dp[nSpace],
                                             const double df_umom_du[nSpace],
                                             const double df_umom_dv[nSpace],
                                             const double df_umom_dw[nSpace],
                                             const double df_vmom_dp[nSpace],
                                             const double df_vmom_du[nSpace],
                                             const double df_vmom_dv[nSpace],
                                             const double df_vmom_dw[nSpace],
                                             const double df_wmom_dp[nSpace],
                                             const double df_wmom_du[nSpace],
                                             const double df_wmom_dv[nSpace],
                                             const double df_wmom_dw[nSpace],
                                             double &flux_mass,
                                             double &flux_umom,
                                             double &flux_vmom,
                                             double &flux_wmom,
                                             double *velocity_star,
                                             double *velocity)
  {
  }

  inline void exteriorNumericalAdvectiveFluxDerivatives(const int &isDOFBoundary_p,
                                                        const int &isDOFBoundary_u,
                                                        const int &isDOFBoundary_v,
                                                        const int &isDOFBoundary_w,
                                                        const int &isFluxBoundary_p,
                                                        const int &isFluxBoundary_u,
                                                        const int &isFluxBoundary_v,
                                                        const int &isFluxBoundary_w,
                                                        const double &oneByRho,
                                                        const double n[nSpace],
                                                        const double &porosity, //mql. CHECK. Multiply by rho outside
                                                        const double &bc_p,
                                                        const double &bc_u,
                                                        const double &bc_v,
                                                        const double &bc_w,
                                                        const double bc_f_mass[nSpace],
                                                        const double bc_f_umom[nSpace],
                                                        const double bc_f_vmom[nSpace],
                                                        const double bc_f_wmom[nSpace],
                                                        const double &bc_flux_mass,
                                                        const double &bc_flux_umom,
                                                        const double &bc_flux_vmom,
                                                        const double &bc_flux_wmom,
                                                        const double &p,
                                                        const double &u,
                                                        const double &v,
                                                        const double &w,
                                                        const double f_mass[nSpace],
                                                        const double f_umom[nSpace],
                                                        const double f_vmom[nSpace],
                                                        const double f_wmom[nSpace],
                                                        const double df_mass_du[nSpace],
                                                        const double df_mass_dv[nSpace],
                                                        const double df_mass_dw[nSpace],
                                                        const double df_umom_dp[nSpace],
                                                        const double df_umom_du[nSpace],
                                                        const double df_umom_dv[nSpace],
                                                        const double df_umom_dw[nSpace],
                                                        const double df_vmom_dp[nSpace],
                                                        const double df_vmom_du[nSpace],
                                                        const double df_vmom_dv[nSpace],
                                                        const double df_vmom_dw[nSpace],
                                                        const double df_wmom_dp[nSpace],
                                                        const double df_wmom_du[nSpace],
                                                        const double df_wmom_dv[nSpace],
                                                        const double df_wmom_dw[nSpace],
                                                        double &dflux_mass_du,
                                                        double &dflux_mass_dv,
                                                        double &dflux_mass_dw,
                                                        double &dflux_umom_dp,
                                                        double &dflux_umom_du,
                                                        double &dflux_umom_dv,
                                                        double &dflux_umom_dw,
                                                        double &dflux_vmom_dp,
                                                        double &dflux_vmom_du,
                                                        double &dflux_vmom_dv,
                                                        double &dflux_vmom_dw,
                                                        double &dflux_wmom_dp,
                                                        double &dflux_wmom_du,
                                                        double &dflux_wmom_dv,
                                                        double &dflux_wmom_dw,
                                                        double *velocity_star)
  {
  }

  inline void exteriorNumericalDiffusiveFlux(const double &eps,
                                             const double &phi,
                                             int *rowptr,
                                             int *colind,
                                             const int &isDOFBoundary,
                                             const int &isFluxBoundary,
                                             const double n[nSpace],
                                             double *bc_a,
                                             const double &bc_u,
                                             const double &bc_flux,
                                             double *a,
                                             const double grad_potential[nSpace],
                                             const double &u,
                                             const double &penalty,
                                             double &flux)
  {
  }

  inline double ExteriorNumericalDiffusiveFluxJacobian(const double &eps,
                                                       const double &phi,
                                                       int *rowptr,
                                                       int *colind,
                                                       const int &isDOFBoundary,
                                                       const int &isFluxBoundary,
                                                       const double n[nSpace],
                                                       double *a,
                                                       const double &v,
                                                       const double grad_v[nSpace],
                                                       const double &penalty)
  {
  }
  void get_symmetric_gradient_dot_vec(const double *grad_u, const double *grad_v, const double *n, double res[2])
  {
    //          res[0] =         2.0*grad_u[0]*n[0]+(grad_u[1]+grad_v[0])*n[1];
    //          res[1] = (grad_v[0]+grad_u[1])*n[0]+          2*grad_v[1]*n[1];
    res[0] = grad_u[0] * n[0] + grad_u[1] * n[1];
    res[1] = grad_v[0] * n[0] + grad_v[1] * n[1];
  }
  double get_cross_product(const double *u, const double *v)
  {
    return u[0] * v[1] - u[1] * v[0];
  }
  double get_dot_product(const double *u, const double *v)
  {
    return u[0] * v[0] + u[1] * v[1];
  }
  int get_distance_to_ball(int n_balls, double *ball_center, double *ball_radius, double x, double y, double z, double &distance)
  {
    distance = 1e10;
    int index = -1;
    double d_ball_i;
    for (int i = 0; i < n_balls; ++i)
    {
      d_ball_i = std::sqrt((ball_center[i * 3 + 0] - x) * (ball_center[i * 3 + 0] - x) + (ball_center[i * 3 + 1] - y) * (ball_center[i * 3 + 1] - y)
                           //                                  +(ball_center[i*3+2]-z)*(ball_center[i*3+2]-z)
                           ) -
                 ball_radius[i];
      if (d_ball_i < distance)
      {
        distance = d_ball_i;
        index = i;
      }
    }
    return index;
  }
  void get_distance_to_ith_ball(int n_balls, double *ball_center, double *ball_radius,
                                int I,
                                double x, double y, double z,
                                double &distance)
  {
    distance = std::sqrt((ball_center[I * 3 + 0] - x) * (ball_center[I * 3 + 0] - x) + (ball_center[I * 3 + 1] - y) * (ball_center[I * 3 + 1] - y)
                         //                                  + (ball_center[I*3+2]-z)*(ball_center[I*3+2]-z)
                         ) -
               ball_radius[I];
  }
  void get_normal_to_ith_ball(int n_balls, double *ball_center, double *ball_radius,
                              int I,
                              double x, double y, double z,
                              double &nx, double &ny)
  {
    double distance = std::sqrt((ball_center[I * 3 + 0] - x) * (ball_center[I * 3 + 0] - x) + (ball_center[I * 3 + 1] - y) * (ball_center[I * 3 + 1] - y)
                                //                                  + (ball_center[I*3+2]-z)*(ball_center[I*3+2]-z)
    );
    nx = (x - ball_center[I * 3 + 0]) / (distance + 1e-10);
    ny = (y - ball_center[I * 3 + 1]) / (distance + 1e-10);
  }
  void get_velocity_to_ith_ball(int n_balls, double *ball_center, double *ball_radius,
                                double *ball_velocity, double *ball_angular_velocity,
                                int I,
                                double x, double y, double z,
                                double &vx, double &vy)
  {
    vx = ball_velocity[3 * I + 0] - ball_angular_velocity[3 * I + 2] * (y - ball_center[3 * I + 1]);
    vy = ball_velocity[3 * I + 1] + ball_angular_velocity[3 * I + 2] * (x - ball_center[3 * I + 0]);
  }
  // void get_curl_cross(double u, double v, double *grad_u, double *grad_v, double *curl_cross)
  // {
  //   curl_cross[0] = -v*(grad_v[0]-grad_u[1]);
  //   curl_cross[1] =  u*(grad_v[0]-grad_u[1]);
  // }
  void get_curl_cross(double u, double v, double *grad_u, double *grad_v, double *curl_cross)
  {
    curl_cross[0] = u*grad_u[0]+v*grad_u[1];
    curl_cross[1] = u*grad_v[0]+v*grad_v[1];
  }
  void calculateResidual( //element
      double *mesh_trial_ref,
      double *mesh_grad_trial_ref,
      double *mesh_dof,
      double *mesh_velocity_dof,
      double MOVING_DOMAIN,
      double PSTAB,
      int *mesh_l2g,
      double *dV_ref,
      int nDOF_per_element_pressure,
      double *p_trial_ref,
      double *p_grad_trial_ref,
      double *p_test_ref,
      double *p_grad_test_ref,
      double *q_ptt_old,/////////////data for pressure
      double *q_pt_old,
      double *q_p_old,
      double *q_ptt,
      double *q_pt,
      double *q_p,
      double *utt_u_dof_old,/////////data for velocity
      double *utt_v_dof_old,
      double *ut_u_dof_old,
      double *ut_v_dof_old,
      double *u_u_dof_old,
      double *u_v_dof_old,
      double *utt_u_dof,
      double *utt_v_dof,
      double *ut_u_dof,
      double *ut_v_dof,
      double *u_u_dof,
      double *u_v_dof,
      double *vel_trial_ref,
      double *vel_grad_trial_ref,
      double *vel_hess_trial_ref,
      double *vel_test_ref,
      double *vel_grad_test_ref,
      //element boundary
      double *mesh_trial_trace_ref,
      double *mesh_grad_trial_trace_ref,
      double *dS_ref,
      double *p_trial_trace_ref,
      double *p_grad_trial_trace_ref,
      double *p_test_trace_ref,
      double *p_grad_test_trace_ref,
      double *vel_trial_trace_ref,
      double *vel_grad_trial_trace_ref,
      double *vel_test_trace_ref,
      double *vel_grad_test_trace_ref,
      double *normal_ref,
      double *boundaryJac_ref,
      //physics
      double eb_adjoint_sigma,
      double *elementDiameter,
      double *nodeDiametersArray,
      double hFactor,
      int nElements_global,
      int nElements_owned,
      int nElementBoundaries_global,
      int nElementBoundaries_owned,
      int nNodes_owned,
      double useRBLES,
      double useMetrics,
      double alphaBDF,
      double epsFact_rho,
      double epsFact_mu,
      double sigma,
      double rho_0,
      double nu_0,
      double rho_1,
      double nu_1,
      double smagorinskyConstant,
      int turbulenceClosureModel,
      double Ct_sge,
      double Cd_sge,
      double C_dc,
      double C_b,
      //VRANS
      const double *eps_solid,
      const double *ebq_global_phi_solid,
      const double *ebq_global_grad_phi_solid,
      const double *ebq_particle_velocity_solid,
      double *phi_solid_nodes,
      double *phi_solid,
      const double *q_velocity_solid,
      const double *q_velocityStar_solid,
      const double *q_vos,
      const double *q_dvos_dt,
      const double *q_grad_vos,
      const double *q_dragAlpha,
      const double *q_dragBeta,
      const double *q_mass_source,
      const double *q_turb_var_0,
      const double *q_turb_var_1,
      const double *q_turb_var_grad_0,
      double *q_eddy_viscosity,
      //
      int *p_l2g,
      int *vel_l2g,
      double *p_dof,
      double *u_dof,
      double *v_dof,
      double *w_dof,
      double *u_dof_old,
      double *v_dof_old,
      double *w_dof_old,
      double *u_dof_old_old,
      double *v_dof_old_old,
      double *w_dof_old_old,
      double *uStar_dof,
      double *vStar_dof,
      double *wStar_dof,
      double *g,
      const double useVF,
      double *vf,
      double *phi,
      double *normal_phi,
      double *kappa_phi,
      double *q_mom_u_acc,
      double *q_mom_v_acc,
      double *q_mom_w_acc,
      double *q_mass_adv,
      double *q_mom_u_acc_beta_bdf, double *q_mom_v_acc_beta_bdf, double *q_mom_w_acc_beta_bdf,
      double *q_dV,
      double *q_dV_last,
      double *q_velocity_sge,
      double *ebqe_velocity_star,
      double *q_cfl,
      double *q_numDiff_u, double *q_numDiff_v, double *q_numDiff_w,
      double *q_numDiff_u_last, double *q_numDiff_v_last, double *q_numDiff_w_last,
      int *sdInfo_u_u_rowptr, int *sdInfo_u_u_colind,
      int *sdInfo_u_v_rowptr, int *sdInfo_u_v_colind,
      int *sdInfo_u_w_rowptr, int *sdInfo_u_w_colind,
      int *sdInfo_v_v_rowptr, int *sdInfo_v_v_colind,
      int *sdInfo_v_u_rowptr, int *sdInfo_v_u_colind,
      int *sdInfo_v_w_rowptr, int *sdInfo_v_w_colind,
      int *sdInfo_w_w_rowptr, int *sdInfo_w_w_colind,
      int *sdInfo_w_u_rowptr, int *sdInfo_w_u_colind,
      int *sdInfo_w_v_rowptr, int *sdInfo_w_v_colind,
      int offset_p, int offset_u, int offset_v, int offset_w,
      int stride_p, int stride_u, int stride_v, int stride_w,
      double *globalResidual,
      int nExteriorElementBoundaries_global,
      int *exteriorElementBoundariesArray,
      int *elementBoundariesArray,
      int *elementBoundaryElementsArray,
      int *elementBoundaryLocalElementBoundariesArray,
      double *ebqe_vf_ext,
      double *bc_ebqe_vf_ext,
      double *ebqe_phi_ext,
      double *bc_ebqe_phi_ext,
      double *ebqe_normal_phi_ext,
      double *ebqe_kappa_phi_ext,
      //VRANS
      const double *ebqe_vos_ext,
      const double *ebqe_turb_var_0,
      const double *ebqe_turb_var_1,
      //VRANS end
      int *isDOFBoundary_p,
      int *isDOFBoundary_u,
      int *isDOFBoundary_v,
      int *isDOFBoundary_w,
      int *isAdvectiveFluxBoundary_p,
      int *isAdvectiveFluxBoundary_u,
      int *isAdvectiveFluxBoundary_v,
      int *isAdvectiveFluxBoundary_w,
      int *isDiffusiveFluxBoundary_u,
      int *isDiffusiveFluxBoundary_v,
      int *isDiffusiveFluxBoundary_w,
      double *ebqe_bc_p_ext,
      double *ebqe_bc_flux_mass_ext,
      double *ebqe_bc_flux_mom_u_adv_ext,
      double *ebqe_bc_flux_mom_v_adv_ext,
      double *ebqe_bc_flux_mom_w_adv_ext,
      double *ebqe_bc_u_ext,
      double *ebqe_bc_flux_u_diff_ext,
      double *ebqe_penalty_ext,
      double *ebqe_bc_v_ext,
      double *ebqe_bc_flux_v_diff_ext,
      double *ebqe_bc_w_ext,
      double *ebqe_bc_flux_w_diff_ext,
      double *q_x,
      double *q_velocity,
      double *ebqe_velocity,
      double *q_grad_u,
      double *q_grad_v,
      double *q_grad_w,
      double *q_divU,
      double *ebqe_grad_u,
      double *ebqe_grad_v,
      double *ebqe_grad_w,
      double *flux,
      double *elementResidual_p_save,
      int *elementFlags,
      int *boundaryFlags,
      double *barycenters,
      double *wettedAreas,
      double *netForces_p,
      double *netForces_v,
      double *netMoments,
      double *q_rho,
      double *ebqe_rho,
      double *q_nu,
      double *ebqe_nu,
      int nParticles,
      double particle_epsFact,
      double particle_alpha,
      double particle_beta,
      double particle_penalty_constant,
      double *particle_signed_distances,
      double *particle_signed_distance_normals,
      double *particle_velocities,
      double *particle_centroids,
      double *particle_netForces,
      double *particle_netMoments,
      double *particle_surfaceArea,
      double particle_nitsche,
      int use_ball_as_particle,
      double *ball_center,
      double *ball_radius,
      double *ball_velocity,
      double *ball_angular_velocity,
      double *phisError,
      double *phisErrorNodal,
      int USE_SUPG,
      int ARTIFICIAL_VISCOSITY,
      double cMax,
      double cE,
      int MULTIPLY_EXTERNAL_FORCE_BY_DENSITY,
      double *forcex,
      double *forcey,
      double *forcez,
      int KILL_PRESSURE_TERM,
      double dt,
      double *quantDOFs,
      int MATERIAL_PARAMETERS_AS_FUNCTION,
      double *density_as_function,
      double *dynamic_viscosity_as_function,
      double *ebqe_density_as_function,
      double *ebqe_dynamic_viscosity_as_function,
      double order_polynomial,
      double *isActiveDOF,
      int USE_SBM,
      double *ncDrag,
      double *betaDrag,
      double *vos_vel_nodes,
      // For edge based discretization
      double *entropyResidualPerNode,
      double *laggedEntropyResidualPerNode,
      double *uStar_dMatrix,
      double *vStar_dMatrix,
      double *wStar_dMatrix,
      int numDOFs_1D,
      int NNZ_1D,
      int *csrRowIndeces_1D, int *csrColumnOffsets_1D,
      int *rowptr_1D, int *colind_1D,
      double *isBoundary_1D,
      // int by parts pressure
      int INT_BY_PARTS_PRESSURE)
  {

    double cut_cell_boundary_length = 0.0, p_force_x = 0.0, p_force_y = 0.0;
    register double element_uStar_He[nElements_global], element_vStar_He[nElements_global];

    //
    //Loop over elements to compute volume integrals and load them into element and global residual
    //
    double mesh_volume_conservation = 0.0,
           mesh_volume_conservation_weak = 0.0,
           mesh_volume_conservation_err_max = 0.0,
           mesh_volume_conservation_err_max_weak = 0.0;
    double globalConservationError = 0.0;
    const int nQuadraturePoints_global(nElements_global * nQuadraturePoints_element);

    //std::set<int> active_velocity_dof;
    for (int eN = 0; eN < nElements_global; eN++)
    {
      register double elementTransport[nDOF_test_element][nDOF_trial_element];
      register double elementTransposeTransport[nDOF_test_element][nDOF_trial_element];
      //declare local storage for element residual and initialize
      register double elementResidual_p[nDOF_test_element], elementResidual_mesh[nDOF_test_element],
          elementResidual_u[nDOF_test_element],
          elementResidual_v[nDOF_test_element],
          mom_u_source_i[nDOF_test_element],
          mom_v_source_i[nDOF_test_element],
          betaDrag_i[nDOF_test_element],
          vos_i[nDOF_test_element],
          phisErrorElement[nDOF_test_element],
          //elementResidual_w[nDOF_test_element],
          elementEntropyResidual[nDOF_test_element],
          eps_rho, eps_mu;
      //const double* elementResidual_w(NULL);
      double element_active = 1.0; //use 1 since by default it is ibm
      double mesh_volume_conservation_element = 0.0,
             mesh_volume_conservation_element_weak = 0.0;
      // for entropy viscosity
      double linVisc_eN = 0, nlinVisc_eN_num = 0, nlinVisc_eN_den = 0;
      // for hessians of uStar
      double det_hess_uStar_Ke = 0.0, det_hess_vStar_Ke = 0.0, area_Ke = 0.0;
      for (int i = 0; i < nDOF_test_element; i++)
      {
        int eN_i = eN * nDOF_test_element + i;
        elementResidual_p_save[eN_i] = 0.0;
        elementResidual_mesh[i] = 0.0;
        elementResidual_p[i] = 0.0;
        elementResidual_u[i] = 0.0;
        elementResidual_v[i] = 0.0;
        mom_u_source_i[i] = 0.0;
        mom_v_source_i[i] = 0.0;
        betaDrag_i[i] = 0.0;
        vos_i[i] = 0.0;
        phisErrorElement[i] = 0.0;
        /* elementResidual_w[i]=0.0; */
      } //i
      //
      //loop over quadrature points and compute integrands
      //
      for (int k = 0; k < nQuadraturePoints_element; k++)
      {
        //compute indices and declare local storage
        register int eN_k = eN * nQuadraturePoints_element + k,
                     eN_k_nSpace = eN_k * nSpace,
                     eN_k_3d = eN_k * 3,
                     eN_nDOF_trial_element = eN * nDOF_trial_element;
        register double p = 0.0, u = 0.0, v = 0.0, w = 0.0, un = 0.0, vn = 0.0, wn = 0.0,
                        grad_p[nSpace], grad_u[nSpace], grad_v[nSpace], grad_w[nSpace],
                        grad_un[nSpace],grad_vn[nSpace],
                        
                        u2_n = 0.0, u1_n = 0.0, u0_n = 0.0,
                        u2   = 0.0, u1   = 0.0, u0   = 0.0,
                        v2_n = 0.0, v1_n = 0.0, v0_n = 0.0,
                        v2   = 0.0, v1   = 0.0, v0   = 0.0,
                        p2_n = 0.0, p1_n = 0.0, p0_n = 0.0,
                        p2   = 0.0, p1   = 0.0, p0   = 0.0,
                        grad_u2_n[nSpace],grad_v2_n[nSpace],
                        grad_u1_n[nSpace],grad_v1_n[nSpace],
                        grad_u0_n[nSpace],grad_v0_n[nSpace],
                        grad_u2[nSpace],grad_v2[nSpace],
                        grad_u1[nSpace],grad_v1[nSpace],
                        grad_u0[nSpace],grad_v0[nSpace],

                        hess_u[nSpace2], hess_v[nSpace2],
                        mom_u_acc = 0.0,
                        dmom_u_acc_u = 0.0,
                        mom_v_acc = 0.0,
                        dmom_v_acc_v = 0.0,
                        mom_w_acc = 0.0,
                        dmom_w_acc_w = 0.0,
                        mass_adv[nSpace],
                        dmass_adv_u[nSpace],
                        dmass_adv_v[nSpace],
                        dmass_adv_w[nSpace],
                        mom_u_adv[nSpace],
                        dmom_u_adv_u[nSpace],
                        dmom_u_adv_v[nSpace],
                        dmom_u_adv_w[nSpace],
                        mom_v_adv[nSpace],
                        dmom_v_adv_u[nSpace],
                        dmom_v_adv_v[nSpace],
                        dmom_v_adv_w[nSpace],
                        mom_w_adv[nSpace],
                        dmom_w_adv_u[nSpace],
                        dmom_w_adv_v[nSpace],
                        dmom_w_adv_w[nSpace],
                        mom_uu_diff_ten[nSpace],
                        mom_vv_diff_ten[nSpace],
                        mom_ww_diff_ten[nSpace],
                        mom_uv_diff_ten[1],
                        mom_uw_diff_ten[1],
                        mom_vu_diff_ten[1],
                        mom_vw_diff_ten[1],
                        mom_wu_diff_ten[1],
                        mom_wv_diff_ten[1],
                        mom_u_source = 0.0,
                        mom_v_source = 0.0,
                        mom_w_source = 0.0,
                        mom_u_ham = 0.0,
                        dmom_u_ham_grad_p[nSpace],
                        dmom_u_ham_grad_u[nSpace],
                        mom_v_ham = 0.0,
                        dmom_v_ham_grad_p[nSpace],
                        dmom_v_ham_grad_v[nSpace],
                        mom_w_ham = 0.0,
                        dmom_w_ham_grad_p[nSpace],
                        dmom_w_ham_grad_w[nSpace],
                        mom_u_acc_t = 0.0,
                        dmom_u_acc_u_t = 0.0,
                        mom_v_acc_t = 0.0,
                        dmom_v_acc_v_t = 0.0,
                        mom_w_acc_t = 0.0,
                        dmom_w_acc_w_t = 0.0,
                        pdeResidual_p = 0.0,
                        pdeResidual_u = 0.0,
                        pdeResidual_v = 0.0,
                        pdeResidual_w = 0.0,
                        Lstar_u_p[nDOF_test_element],
                        Lstar_v_p[nDOF_test_element],
                        Lstar_w_p[nDOF_test_element],
                        Lstar_u_u[nDOF_test_element],
                        Lstar_v_v[nDOF_test_element],
                        Lstar_w_w[nDOF_test_element],
                        Lstar_p_u[nDOF_test_element],
                        Lstar_p_v[nDOF_test_element],
                        Lstar_p_w[nDOF_test_element],
                        subgridError_p = 0.0,
                        subgridError_u = 0.0,
                        subgridError_v = 0.0,
                        subgridError_w = 0.0,
                        tau_p = 0.0, tau_p0 = 0.0, tau_p1 = 0.0,
                        tau_v = 0.0, tau_v0 = 0.0, tau_v1 = 0.0,
                        jac[nSpace * nSpace],
                        jacDet,
                        jacInv[nSpace * nSpace],
                        p_grad_trial[nDOF_trial_element * nSpace], vel_grad_trial[nDOF_trial_element * nSpace],
                        vel_hess_trial[nDOF_trial_element * nSpace2],
                        p_test_dV[nDOF_trial_element], vel_test_dV[nDOF_trial_element],
                        p_grad_test_dV[nDOF_test_element * nSpace], vel_grad_test_dV[nDOF_test_element * nSpace],
                        u_times_vel_grad_test_dV[nDOF_test_element * nSpace], // For entropy residual
            v_times_vel_grad_test_dV[nDOF_test_element * nSpace],             // For entropy residual
            dV, x, y, z, xt, yt, zt,
                        //
            porosity,
                        //meanGrainSize,
            mass_source,
                        dmom_u_source[nSpace],
                        dmom_v_source[nSpace],
                        dmom_w_source[nSpace],
                        //
            velStar[nSpace], hess_uStar[nSpace2], hess_vStar[nSpace2],
                        //
            G[nSpace * nSpace], G_dd_G, tr_G, norm_Rv, h_phi, dmom_adv_star[nSpace], dmom_adv_sge[nSpace];
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
                                    x, y, z);
        ck.calculateH_element(eN,
                              k,
                              nodeDiametersArray,
                              mesh_l2g,
                              mesh_trial_ref,
                              h_phi);
        ck.calculateMappingVelocity_element(eN,
                                            k,
                                            mesh_velocity_dof,
                                            mesh_l2g,
                                            mesh_trial_ref,
                                            xt, yt, zt);
        //xt=0.0;yt=0.0;zt=0.0;
        //std::cout<<"xt "<<xt<<'\t'<<yt<<'\t'<<zt<<std::endl;
        //get the physical integration weight
        dV = fabs(jacDet) * dV_ref[k];
        ck.calculateG(jacInv, G, G_dd_G, tr_G);
        //ck.calculateGScale(G,&normal_phi[eN_k_nSpace],h_phi);

        eps_rho = epsFact_rho * (useMetrics * h_phi + (1.0 - useMetrics) * elementDiameter[eN]);
        eps_mu = epsFact_mu * (useMetrics * h_phi + (1.0 - useMetrics) * elementDiameter[eN]);
        double particle_eps = particle_epsFact * (useMetrics * h_phi + (1.0 - useMetrics) * elementDiameter[eN]);

        //get the trial function gradients
        /* ck.gradTrialFromRef(&p_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,p_grad_trial); */
        ck.gradTrialFromRef(&vel_grad_trial_ref[k * nDOF_trial_element * nSpace], jacInv, vel_grad_trial);
        // ck.hessTrialFromRef(&vel_hess_trial_ref[k * nDOF_trial_element * nSpace2], jacInv, vel_hess_trial);
        //get the solution
        /* ck.valFromDOF(p_dof,&p_l2g[eN_nDOF_trial_element],&p_trial_ref[k*nDOF_trial_element],p); */
        p = q_p[eN_k];
        // get solution at quad points
        ck.valFromDOF(u_dof, &vel_l2g[eN_nDOF_trial_element], &vel_trial_ref[k * nDOF_trial_element], u);
        ck.valFromDOF(v_dof, &vel_l2g[eN_nDOF_trial_element], &vel_trial_ref[k * nDOF_trial_element], v);
        /* ck.valFromDOF(w_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],w); */
        // get old solution at quad points
        ck.valFromDOF(u_dof_old, &vel_l2g[eN_nDOF_trial_element], &vel_trial_ref[k * nDOF_trial_element], un);
        ck.valFromDOF(v_dof_old, &vel_l2g[eN_nDOF_trial_element], &vel_trial_ref[k * nDOF_trial_element], vn);
        /* ck.valFromDOF(w_dof_old,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],wn); */
        //get the solution gradients
        /* ck.gradFromDOF(p_dof,&p_l2g[eN_nDOF_trial_element],p_grad_trial,grad_p); */
        // for (int I = 0; I < nSpace; I++)
        //   grad_p[I] = q_grad_p[eN_k_nSpace + I];
        ck.gradFromDOF(u_dof, &vel_l2g[eN_nDOF_trial_element], vel_grad_trial, grad_u);
        ck.gradFromDOF(v_dof, &vel_l2g[eN_nDOF_trial_element], vel_grad_trial, grad_v);
        ck.gradFromDOF(u_dof_old, &vel_l2g[eN_nDOF_trial_element], vel_grad_trial, grad_un);
        ck.gradFromDOF(v_dof_old, &vel_l2g[eN_nDOF_trial_element], vel_grad_trial, grad_vn);

        //
        // COMPUTE VALUES
        //
        p2_n = q_ptt_old[eN_k];
        p1_n = q_pt_old[eN_k];
        p0_n = q_p_old[eN_k];
        p2   = q_ptt[eN_k];
        p1   = q_pt[eN_k];
        p0   = q_p[eN_k];

        ck.valFromDOF(utt_u_dof_old, &vel_l2g[eN_nDOF_trial_element], &vel_trial_ref[k * nDOF_trial_element], u2_n);
        ck.valFromDOF(utt_v_dof_old, &vel_l2g[eN_nDOF_trial_element], &vel_trial_ref[k * nDOF_trial_element], v2_n);
        ck.valFromDOF(ut_u_dof_old,  &vel_l2g[eN_nDOF_trial_element], &vel_trial_ref[k * nDOF_trial_element], u1_n);
        ck.valFromDOF(ut_v_dof_old,  &vel_l2g[eN_nDOF_trial_element], &vel_trial_ref[k * nDOF_trial_element], v1_n);
        ck.valFromDOF(u_u_dof_old,   &vel_l2g[eN_nDOF_trial_element], &vel_trial_ref[k * nDOF_trial_element], u0_n);
        ck.valFromDOF(u_v_dof_old,   &vel_l2g[eN_nDOF_trial_element], &vel_trial_ref[k * nDOF_trial_element], v0_n);
        ck.valFromDOF(utt_u_dof, &vel_l2g[eN_nDOF_trial_element], &vel_trial_ref[k * nDOF_trial_element], u2);
        ck.valFromDOF(utt_v_dof, &vel_l2g[eN_nDOF_trial_element], &vel_trial_ref[k * nDOF_trial_element], v2);
        ck.valFromDOF(ut_u_dof,  &vel_l2g[eN_nDOF_trial_element], &vel_trial_ref[k * nDOF_trial_element], u1);
        ck.valFromDOF(ut_v_dof,  &vel_l2g[eN_nDOF_trial_element], &vel_trial_ref[k * nDOF_trial_element], v1);
        ck.valFromDOF(u_u_dof,   &vel_l2g[eN_nDOF_trial_element], &vel_trial_ref[k * nDOF_trial_element], u0);
        ck.valFromDOF(u_v_dof,   &vel_l2g[eN_nDOF_trial_element], &vel_trial_ref[k * nDOF_trial_element], v0);
        
        ck.gradFromDOF(utt_u_dof_old, &vel_l2g[eN_nDOF_trial_element], vel_grad_trial, grad_u2_n);
        ck.gradFromDOF(utt_v_dof_old, &vel_l2g[eN_nDOF_trial_element], vel_grad_trial, grad_v2_n);
        ck.gradFromDOF(ut_u_dof_old,  &vel_l2g[eN_nDOF_trial_element], vel_grad_trial, grad_u1_n);
        ck.gradFromDOF(ut_v_dof_old,  &vel_l2g[eN_nDOF_trial_element], vel_grad_trial, grad_v1_n);
        ck.gradFromDOF(u_u_dof_old,   &vel_l2g[eN_nDOF_trial_element], vel_grad_trial, grad_u0_n);
        ck.gradFromDOF(u_v_dof_old,   &vel_l2g[eN_nDOF_trial_element], vel_grad_trial, grad_v0_n);
        ck.gradFromDOF(utt_u_dof, &vel_l2g[eN_nDOF_trial_element], vel_grad_trial, grad_u2);
        ck.gradFromDOF(utt_v_dof, &vel_l2g[eN_nDOF_trial_element], vel_grad_trial, grad_v2);
        ck.gradFromDOF(ut_u_dof,  &vel_l2g[eN_nDOF_trial_element], vel_grad_trial, grad_u1);
        ck.gradFromDOF(ut_v_dof,  &vel_l2g[eN_nDOF_trial_element], vel_grad_trial, grad_v1);
        ck.gradFromDOF(u_u_dof,   &vel_l2g[eN_nDOF_trial_element], vel_grad_trial, grad_u0);
        ck.gradFromDOF(u_v_dof,   &vel_l2g[eN_nDOF_trial_element], vel_grad_trial, grad_v0);
        

        double curl_cross_old[2];
        get_curl_cross(un, vn, grad_un, grad_vn, curl_cross_old);

        /* ck.gradFromDOF(w_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_w); */
        //precalculate test function products with integration weights
        for (int j = 0; j < nDOF_trial_element; j++)
        {
          /* p_test_dV[j] = p_test_ref[k*nDOF_trial_element+j]*dV; */
          vel_test_dV[j] = vel_test_ref[k * nDOF_trial_element + j] * dV;
          for (int I = 0; I < nSpace; I++)
          {
            /* p_grad_test_dV[j*nSpace+I]   = p_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin */
            vel_grad_test_dV[j * nSpace + I] = vel_grad_trial[j * nSpace + I] * dV; //cek warning won't work for Petrov-Galerkin
          }
        }

        //cek hack
        double div_mesh_velocity = 0.0;
        int NDOF_MESH_TRIAL_ELEMENT = 3;
        for (int j = 0; j < NDOF_MESH_TRIAL_ELEMENT; j++)
        {
          int eN_j = eN * NDOF_MESH_TRIAL_ELEMENT + j;
          div_mesh_velocity +=
              mesh_velocity_dof[mesh_l2g[eN_j] * 3 + 0] * vel_grad_trial[j * nSpace + 0] +
              mesh_velocity_dof[mesh_l2g[eN_j] * 3 + 1] * vel_grad_trial[j * nSpace + 1];
        }
        mesh_volume_conservation_element += (alphaBDF * (dV - q_dV_last[eN_k]) / dV - div_mesh_velocity) * dV;
        div_mesh_velocity = DM3 * div_mesh_velocity + (1.0 - DM3) * alphaBDF * (dV - q_dV_last[eN_k]) / dV;

        porosity = 1.0;
        //
        q_x[eN_k_3d + 0] = x;
        q_x[eN_k_3d + 1] = y;
        
        double distance_to_omega_solid = 1e10;
        phi_solid[eN_k] = distance_to_omega_solid; //save it

        evaluateCoefficients(eps_rho,
                             eps_mu,
                             particle_eps,
                             sigma,
                             rho_0,
                             nu_0,
                             rho_1,
                             nu_1,
                             elementDiameter[eN],
                             smagorinskyConstant,
                             turbulenceClosureModel,
                             g,
                             useVF,
                             vf[eN_k],
                             phi[eN_k],
                             &normal_phi[eN_k_nSpace],
                             distance_to_omega_solid,
                             kappa_phi[eN_k],
                             porosity,
                             p,
                             grad_p,
                             grad_u,
                             grad_v,
                             grad_w,
                             u,
                             v,
                             w,
                             q_velocity_sge[eN_k_nSpace + 0],
                             q_velocity_sge[eN_k_nSpace + 1],
                             q_velocity_sge[eN_k_nSpace + 1], //hack, shouldn't  be used
                             q_eddy_viscosity[eN_k],
                             mom_u_acc,
                             dmom_u_acc_u,
                             mom_v_acc,
                             dmom_v_acc_v,
                             mom_w_acc,
                             dmom_w_acc_w,
                             mass_adv,
                             dmass_adv_u,
                             dmass_adv_v,
                             dmass_adv_w,
                             mom_u_adv,
                             dmom_u_adv_u,
                             dmom_u_adv_v,
                             dmom_u_adv_w,
                             mom_v_adv,
                             dmom_v_adv_u,
                             dmom_v_adv_v,
                             dmom_v_adv_w,
                             mom_w_adv,
                             dmom_w_adv_u,
                             dmom_w_adv_v,
                             dmom_w_adv_w,
                             mom_uu_diff_ten,
                             mom_vv_diff_ten,
                             mom_ww_diff_ten,
                             mom_uv_diff_ten,
                             mom_uw_diff_ten,
                             mom_vu_diff_ten,
                             mom_vw_diff_ten,
                             mom_wu_diff_ten,
                             mom_wv_diff_ten,
                             mom_u_source,
                             mom_v_source,
                             mom_w_source,
                             mom_u_ham,
                             dmom_u_ham_grad_p,
                             dmom_u_ham_grad_u,
                             mom_v_ham,
                             dmom_v_ham_grad_p,
                             dmom_v_ham_grad_v,
                             mom_w_ham,
                             dmom_w_ham_grad_p,
                             dmom_w_ham_grad_w,
                             q_rho[eN_k],
                             q_nu[eN_k],
                             KILL_PRESSURE_TERM,
                             MULTIPLY_EXTERNAL_FORCE_BY_DENSITY,
                             forcex[eN_k],
                             forcey[eN_k],
                             forcez[eN_k],
                             MATERIAL_PARAMETERS_AS_FUNCTION,
                             density_as_function[eN_k],
                             dynamic_viscosity_as_function[eN_k],
                             USE_SBM,
                             x, y, z,
                             use_ball_as_particle,
                             ball_center,
                             ball_radius,
                             ball_velocity,
                             ball_angular_velocity,
                             INT_BY_PARTS_PRESSURE);

        //VRANS
        mass_source = q_mass_source[eN_k];
        for (int I = 0; I < nSpace; I++)
        {
          dmom_u_source[I] = 0.0;
          dmom_v_source[I] = 0.0;
          dmom_w_source[I] = 0.0;
        }

        double C_particles = 0.0;

        //
        //save momentum for time history and velocity for subgrid error
        //
        q_mom_u_acc[eN_k] = mom_u_acc;
        q_mom_v_acc[eN_k] = mom_v_acc;
        /* q_mom_w_acc[eN_k] = mom_w_acc; */
        //subgrid error uses grid scale velocity
        q_mass_adv[eN_k_nSpace + 0] = u;
        q_mass_adv[eN_k_nSpace + 1] = v;
        /* q_mass_adv[eN_k_nSpace+2] = w; */
        //
        //moving mesh
        //
        mom_u_adv[0] -= MOVING_DOMAIN * dmom_u_acc_u * mom_u_acc * xt; // multiply by rho*porosity. mql. CHECK.
        mom_u_adv[1] -= MOVING_DOMAIN * dmom_u_acc_u * mom_u_acc * yt;
        /* mom_u_adv[2] -= MOVING_DOMAIN*dmom_u_acc_u*mom_u_acc*zt; */
        dmom_u_adv_u[0] -= MOVING_DOMAIN * dmom_u_acc_u * xt;
        dmom_u_adv_u[1] -= MOVING_DOMAIN * dmom_u_acc_u * yt;
        /* dmom_u_adv_u[2] -= MOVING_DOMAIN*dmom_u_acc_u*zt; */

        mom_v_adv[0] -= MOVING_DOMAIN * dmom_v_acc_v * mom_v_acc * xt;
        mom_v_adv[1] -= MOVING_DOMAIN * dmom_v_acc_v * mom_v_acc * yt;
        /* mom_v_adv[2] -= MOVING_DOMAIN*dmom_v_acc_v*mom_v_acc*zt; */
        dmom_v_adv_v[0] -= MOVING_DOMAIN * dmom_v_acc_v * xt;
        dmom_v_adv_v[1] -= MOVING_DOMAIN * dmom_v_acc_v * yt;
        /* dmom_v_adv_v[2] -= MOVING_DOMAIN*dmom_v_acc_v*zt; */

        /* mom_w_adv[0] -= MOVING_DOMAIN*dmom_w_acc_w*mom_w_acc*xt; */
        /* mom_w_adv[1] -= MOVING_DOMAIN*dmom_w_acc_w*mom_w_acc*yt; */
        /* mom_w_adv[2] -= MOVING_DOMAIN*dmom_w_acc_w*mom_w_acc*zt; */
        /* dmom_w_adv_w[0] -= MOVING_DOMAIN*dmom_w_acc_w*xt; */
        /* dmom_w_adv_w[1] -= MOVING_DOMAIN*dmom_w_acc_w*yt; */
        /* dmom_w_adv_w[2] -= MOVING_DOMAIN*dmom_w_acc_w*zt; */
        //
        //calculate time derivative at quadrature points
        //
        if (q_dV_last[eN_k] <= -100)
          q_dV_last[eN_k] = dV;
        q_dV[eN_k] = dV;
        ck.bdf(alphaBDF,
               q_mom_u_acc_beta_bdf[eN_k] * q_dV_last[eN_k] / dV,
               mom_u_acc,
               dmom_u_acc_u,
               mom_u_acc_t,
               dmom_u_acc_u_t);
        ck.bdf(alphaBDF,
               q_mom_v_acc_beta_bdf[eN_k] * q_dV_last[eN_k] / dV,
               mom_v_acc,
               dmom_v_acc_v,
               mom_v_acc_t,
               dmom_v_acc_v_t);

        /* ck.bdf(alphaBDF, */
        /*           q_mom_w_acc_beta_bdf[eN_k]*q_dV_last[eN_k]/dV, */
        /*           mom_w_acc, */
        /*           dmom_w_acc_w, */
        /*           mom_w_acc_t, */
        /*           dmom_w_acc_w_t); */
        /* // */

        mom_u_acc_t *= dmom_u_acc_u; //multiply by rho*porosity. mql. CHECK.
        mom_v_acc_t *= dmom_v_acc_v;

        //
        //update element residual
        //
        double mesh_vel[2];
        // Save velocity and its gradient (to be used in other models and to compute errors)
        q_velocity[eN_k_nSpace + 0] = u;
        q_velocity[eN_k_nSpace + 1] = v;
        /* q_velocity[eN_k_nSpace+2]=w; */
        for (int I = 0; I < nSpace; I++)
        {
          q_grad_u[eN_k_nSpace + I] = grad_u[I];
          q_grad_v[eN_k_nSpace + I] = grad_v[I];
          /* q_grad_w[eN_k_nSpace+I] = grad_w[I]; */
        }
        // save divergence of velocity
        q_divU[eN_k] = q_grad_u[eN_k_nSpace + 0] + q_grad_v[eN_k_nSpace + 1];
        const double lambda = 1.0;
        double p_prediction = 0.0;
        double  grad_u_np2[nSpace], grad_v_np2[nSpace],curl_cross_np2[nSpace],
                grad_u_np1[nSpace], grad_v_np1[nSpace],curl_cross_np1[nSpace],
                grad_u_np0[nSpace], grad_v_np0[nSpace],curl_cross_np0[nSpace],
                u_np2,v_np2,
                u_np1,v_np1,
                u_np0,v_np0;
        u_np2 = u0_n + 2.0*u1_n/alphaBDF + 2.0*u2_n/alphaBDF/alphaBDF;
        v_np2 = v0_n + 2.0*v1_n/alphaBDF + 2.0*v2_n/alphaBDF/alphaBDF;
        u_np1 = u0_n + u1_n/alphaBDF + 0.5*u2_n/alphaBDF/alphaBDF;
        v_np1 = v0_n + v1_n/alphaBDF + 0.5*v2_n/alphaBDF/alphaBDF;
        u_np0 = u0_n;
        v_np0 = v0_n;
        for (int I = 0; I < nSpace; I++)
        {
          grad_u_np2[I] = grad_u0_n[I] + 2.0*grad_u1_n[I]/alphaBDF + 2.0*grad_u2_n[I]/alphaBDF/alphaBDF;
          grad_v_np2[I] = grad_v0_n[I] + 2.0*grad_v1_n[I]/alphaBDF + 2.0*grad_v2_n[I]/alphaBDF/alphaBDF;
          grad_u_np1[I] = grad_u0_n[I] +     grad_u1_n[I]/alphaBDF + 0.5*grad_u2_n[I]/alphaBDF/alphaBDF;
          grad_v_np1[I] = grad_v0_n[I] +     grad_v1_n[I]/alphaBDF + 0.5*grad_v2_n[I]/alphaBDF/alphaBDF;
          grad_u_np0[I] = grad_u0_n[I];
          grad_v_np0[I] = grad_v0_n[I];
        }
        get_curl_cross(u_np2, v_np2, grad_u_np2, grad_v_np2, curl_cross_np2);
        get_curl_cross(u_np1, v_np1, grad_u_np1, grad_v_np1, curl_cross_np1);
        get_curl_cross(u_np0, v_np0, grad_u_np0, grad_v_np0, curl_cross_np0);

        if(USE_SBM==-1)////stage-1 or (3.12)
        {
          p_prediction = p2_n;
          for (int I = 0; I < nSpace; I++)
          {
            curl_cross_old[I] = (curl_cross_np2[I]-2.0*curl_cross_np1[I]+curl_cross_np0[I])/alphaBDF/alphaBDF;
          }
        }
        else if(USE_SBM==-2)////stage-2 or (3.13)
        {
          p_prediction = p1_n + p2 / alphaBDF;
          mom_u_source = mom_u_source + 0.5 * (u2 - u2_n);//mom_u_source contains minus sign since it is on lhs
          mom_v_source = mom_v_source + 0.5 * (v2 - v2_n);
          for (int I = 0; I < nSpace; I++)
          {
            curl_cross_old[I] = (curl_cross_np2[I]-curl_cross_np0[I])*0.5/alphaBDF;
          }
        }
        else if(USE_SBM==-3)////stage-3 or (3.14)
        {
          p_prediction = p0_n + p1 / alphaBDF - 0.5 * p2 / alphaBDF / alphaBDF;
          mom_u_source = mom_u_source + 0.5 * (u1 - u1_n) + (u2 - u2_n)/12.0/alphaBDF;//mom_u_source contains minus sign since it is on lhs
          mom_v_source = mom_v_source + 0.5 * (v1 - v1_n) + (v2 - v2_n)/12.0/alphaBDF;
          for (int I = 0; I < nSpace; I++)
          {
            curl_cross_old[I] = curl_cross_np0[I];
          }
        }
        
        for (int i = 0; i < nDOF_test_element; i++)
        {
          register int i_nSpace = i * nSpace;

          elementResidual_u[i] += // mql. CHECK.
              ck.Mass_weak(mom_u_acc_t, vel_test_dV[i])
              + ck.Advection_weak(mom_u_adv, &vel_grad_test_dV[i_nSpace])
              // + mom_uu_diff_ten[0] * grad_u[0] * vel_grad_test_dV[i_nSpace + 0]
              // + mom_uu_diff_ten[1] * grad_u[1] * vel_grad_test_dV[i_nSpace + 1]
              - lambda * q_divU[eN_k] * (vel_grad_test_dV[i_nSpace + 0] + vel_grad_test_dV[i_nSpace + 1])
              + ck.Diffusion_weak(sdInfo_u_u_rowptr, sdInfo_u_u_colind, mom_uu_diff_ten, grad_u, &vel_grad_test_dV[i_nSpace])
              + ck.Diffusion_weak(sdInfo_u_v_rowptr, sdInfo_u_v_colind, mom_uv_diff_ten, grad_v, &vel_grad_test_dV[i_nSpace])
              + ck.Reaction_weak(mom_u_source, vel_test_dV[i])
              + ck.Hamiltonian_weak(mom_u_ham, vel_test_dV[i])
              + curl_cross_old[0] * vel_test_dV[i]
              - p_prediction * vel_grad_test_dV[i_nSpace + 0]
              //ck.SubgridError(subgridError_p,Lstar_p_u[i]) +
              // USE_SUPG * ck.SubgridError(subgridError_u, Lstar_u_u[i]) +
              // ck.NumericalDiffusion(q_numDiff_u_last[eN_k], grad_u, &vel_grad_test_dV[i_nSpace])
              ; //imp.

          elementResidual_v[i] +=
              ck.Mass_weak(mom_v_acc_t, vel_test_dV[i])
              + ck.Advection_weak(mom_v_adv, &vel_grad_test_dV[i_nSpace])
              // + mom_vv_diff_ten[0] * grad_v[0] * vel_grad_test_dV[i_nSpace + 0]
              // + mom_vv_diff_ten[1] * grad_v[1] * vel_grad_test_dV[i_nSpace + 1]
              - lambda * q_divU[eN_k] * (vel_grad_test_dV[i_nSpace + 0] + vel_grad_test_dV[i_nSpace + 1])
              + ck.Diffusion_weak(sdInfo_v_u_rowptr, sdInfo_v_u_colind, mom_vu_diff_ten, grad_u, &vel_grad_test_dV[i_nSpace])
              + ck.Diffusion_weak(sdInfo_v_v_rowptr, sdInfo_v_v_colind, mom_vv_diff_ten, grad_v, &vel_grad_test_dV[i_nSpace])
              + ck.Reaction_weak(mom_v_source, vel_test_dV[i])
              + curl_cross_old[1] * vel_test_dV[i]
              + ck.Hamiltonian_weak(mom_v_ham, vel_test_dV[i])
              - p_prediction * vel_grad_test_dV[i_nSpace + 1]
              //ck.SubgridError(subgridError_p,Lstar_p_v[i]) +
              // USE_SUPG * ck.SubgridError(subgridError_v, Lstar_v_v[i]) +
              // ck.NumericalDiffusion(q_numDiff_v_last[eN_k], grad_v, &vel_grad_test_dV[i_nSpace])
              ; //imp.

        } //i
      }
      element_uStar_He[eN] = det_hess_uStar_Ke / area_Ke;
      element_vStar_He[eN] = det_hess_vStar_Ke / area_Ke;

      //
      //load element into global residual and save element residual
      //
      for (int i = 0; i < nDOF_test_element; i++)
      {
        register int eN_i = eN * nDOF_test_element + i;
        /* elementResidual_p_save[eN_i] +=  elementResidual_p[i]; */
        /* mesh_volume_conservation_element_weak += elementResidual_mesh[i]; */
        /* globalResidual[offset_p+stride_p*p_l2g[eN_i]]+=elementResidual_p[i]; */
        globalResidual[offset_u + stride_u * vel_l2g[eN_i]] += element_active * elementResidual_u[i];
        globalResidual[offset_v + stride_v * vel_l2g[eN_i]] += element_active * elementResidual_v[i];
        /* globalResidual[offset_w+stride_w*vel_l2g[eN_i]]+=elementResidual_w[i]; */

      } //i
      /* mesh_volume_conservation += mesh_volume_conservation_element; */
      /* mesh_volume_conservation_weak += mesh_volume_conservation_element_weak; */
      /* mesh_volume_conservation_err_max=fmax(mesh_volume_conservation_err_max,fabs(mesh_volume_conservation_element)); */
      /* mesh_volume_conservation_err_max_weak=fmax(mesh_volume_conservation_err_max_weak,fabs(mesh_volume_conservation_element_weak)); */
    } //elements

  }

  void calculateJacobian( //element
      double *mesh_trial_ref,
      double *mesh_grad_trial_ref,
      double *mesh_dof,
      double *mesh_velocity_dof,
      double MOVING_DOMAIN,
      double PSTAB,
      int *mesh_l2g,
      double *dV_ref,
      double *p_trial_ref,
      double *p_grad_trial_ref,
      double *p_test_ref,
      double *p_grad_test_ref,
      double *q_p,
      double *q_grad_p,
      double *ebqe_p,
      double *ebqe_grad_p,
      double *vel_trial_ref,
      double *vel_grad_trial_ref,
      double *vel_hess_trial_ref,
      double *vel_test_ref,
      double *vel_grad_test_ref,
      //element boundary
      double *mesh_trial_trace_ref,
      double *mesh_grad_trial_trace_ref,
      double *dS_ref,
      double *p_trial_trace_ref,
      double *p_grad_trial_trace_ref,
      double *p_test_trace_ref,
      double *p_grad_test_trace_ref,
      double *vel_trial_trace_ref,
      double *vel_grad_trial_trace_ref,
      double *vel_test_trace_ref,
      double *vel_grad_test_trace_ref,
      double *normal_ref,
      double *boundaryJac_ref,
      //physics
      double eb_adjoint_sigma,
      double *elementDiameter,
      double *nodeDiametersArray,
      double hFactor,
      int nElements_global,
      int nElements_owned,
      int nElementBoundaries_global,
      int nElementBoundaries_owned,
      int nNodes_owned,
      double useRBLES,
      double useMetrics,
      double alphaBDF,
      double epsFact_rho,
      double epsFact_mu,
      double sigma,
      double rho_0,
      double nu_0,
      double rho_1,
      double nu_1,
      double smagorinskyConstant,
      int turbulenceClosureModel,
      double Ct_sge,
      double Cd_sge,
      double C_dg,
      double C_b,
      //VRANS
      const double *eps_solid,
      const double *ebq_global_phi_solid,
      const double *ebq_global_grad_phi_solid,
      const double *ebq_particle_velocity_solid,
      double *phi_solid_nodes,
      const double *phi_solid,
      const double *q_velocity_solid,
      const double *q_velocityStar_solid,
      const double *q_vos,
      const double *q_dvos_dt,
      const double *q_grad_vos,
      const double *q_dragAlpha,
      const double *q_dragBeta,
      const double *q_mass_source,
      const double *q_turb_var_0,
      const double *q_turb_var_1,
      const double *q_turb_var_grad_0,
      //
      int *p_l2g,
      int *vel_l2g,
      double *p_dof, double *u_dof, double *v_dof, double *w_dof,
      double *g,
      const double useVF,
      double *vf,
      double *phi,
      double *normal_phi,
      double *kappa_phi,
      double *q_mom_u_acc_beta_bdf, double *q_mom_v_acc_beta_bdf, double *q_mom_w_acc_beta_bdf,
      double *q_dV,
      double *q_dV_last,
      double *q_velocity_sge,
      double *ebqe_velocity_star,
      double *q_cfl,
      double *q_numDiff_u_last, double *q_numDiff_v_last, double *q_numDiff_w_last,
      int *sdInfo_u_u_rowptr, int *sdInfo_u_u_colind,
      int *sdInfo_u_v_rowptr, int *sdInfo_u_v_colind,
      int *sdInfo_u_w_rowptr, int *sdInfo_u_w_colind,
      int *sdInfo_v_v_rowptr, int *sdInfo_v_v_colind,
      int *sdInfo_v_u_rowptr, int *sdInfo_v_u_colind,
      int *sdInfo_v_w_rowptr, int *sdInfo_v_w_colind,
      int *sdInfo_w_w_rowptr, int *sdInfo_w_w_colind,
      int *sdInfo_w_u_rowptr, int *sdInfo_w_u_colind,
      int *sdInfo_w_v_rowptr, int *sdInfo_w_v_colind,
      int *csrRowIndeces_p_p, int *csrColumnOffsets_p_p,
      int *csrRowIndeces_p_u, int *csrColumnOffsets_p_u,
      int *csrRowIndeces_p_v, int *csrColumnOffsets_p_v,
      int *csrRowIndeces_p_w, int *csrColumnOffsets_p_w,
      int *csrRowIndeces_u_p, int *csrColumnOffsets_u_p,
      int *csrRowIndeces_u_u, int *csrColumnOffsets_u_u,
      int *csrRowIndeces_u_v, int *csrColumnOffsets_u_v,
      int *csrRowIndeces_u_w, int *csrColumnOffsets_u_w,
      int *csrRowIndeces_v_p, int *csrColumnOffsets_v_p,
      int *csrRowIndeces_v_u, int *csrColumnOffsets_v_u,
      int *csrRowIndeces_v_v, int *csrColumnOffsets_v_v,
      int *csrRowIndeces_v_w, int *csrColumnOffsets_v_w,
      int *csrRowIndeces_w_p, int *csrColumnOffsets_w_p,
      int *csrRowIndeces_w_u, int *csrColumnOffsets_w_u,
      int *csrRowIndeces_w_v, int *csrColumnOffsets_w_v,
      int *csrRowIndeces_w_w, int *csrColumnOffsets_w_w,
      double *globalJacobian,
      int nExteriorElementBoundaries_global,
      int *exteriorElementBoundariesArray,
      int *elementBoundariesArray,
      int *elementBoundaryElementsArray,
      int *elementBoundaryLocalElementBoundariesArray,
      double *ebqe_vf_ext,
      double *bc_ebqe_vf_ext,
      double *ebqe_phi_ext,
      double *bc_ebqe_phi_ext,
      double *ebqe_normal_phi_ext,
      double *ebqe_kappa_phi_ext,
      //VRANS
      const double *ebqe_vos_ext,
      const double *ebqe_turb_var_0,
      const double *ebqe_turb_var_1,
      //
      int *isDOFBoundary_p,
      int *isDOFBoundary_u,
      int *isDOFBoundary_v,
      int *isDOFBoundary_w,
      int *isAdvectiveFluxBoundary_p,
      int *isAdvectiveFluxBoundary_u,
      int *isAdvectiveFluxBoundary_v,
      int *isAdvectiveFluxBoundary_w,
      int *isDiffusiveFluxBoundary_u,
      int *isDiffusiveFluxBoundary_v,
      int *isDiffusiveFluxBoundary_w,
      double *ebqe_bc_p_ext,
      double *ebqe_bc_flux_mass_ext,
      double *ebqe_bc_flux_mom_u_adv_ext,
      double *ebqe_bc_flux_mom_v_adv_ext,
      double *ebqe_bc_flux_mom_w_adv_ext,
      double *ebqe_bc_u_ext,
      double *ebqe_bc_flux_u_diff_ext,
      double *ebqe_penalty_ext,
      double *ebqe_bc_v_ext,
      double *ebqe_bc_flux_v_diff_ext,
      double *ebqe_bc_w_ext,
      double *ebqe_bc_flux_w_diff_ext,
      int *csrColumnOffsets_eb_p_p,
      int *csrColumnOffsets_eb_p_u,
      int *csrColumnOffsets_eb_p_v,
      int *csrColumnOffsets_eb_p_w,
      int *csrColumnOffsets_eb_u_p,
      int *csrColumnOffsets_eb_u_u,
      int *csrColumnOffsets_eb_u_v,
      int *csrColumnOffsets_eb_u_w,
      int *csrColumnOffsets_eb_v_p,
      int *csrColumnOffsets_eb_v_u,
      int *csrColumnOffsets_eb_v_v,
      int *csrColumnOffsets_eb_v_w,
      int *csrColumnOffsets_eb_w_p,
      int *csrColumnOffsets_eb_w_u,
      int *csrColumnOffsets_eb_w_v,
      int *csrColumnOffsets_eb_w_w,
      int *elementFlags,
      int nParticles,
      double particle_epsFact,
      double particle_alpha,
      double particle_beta,
      double particle_penalty_constant,
      double *particle_signed_distances,
      double *particle_signed_distance_normals,
      double *particle_velocities,
      double *particle_centroids,
      double particle_nitsche,
      int use_ball_as_particle,
      double *ball_center,
      double *ball_radius,
      double *ball_velocity,
      double *ball_angular_velocity,
      int USE_SUPG,
      int KILL_PRESSURE_TERM,
      double dt,
      int MATERIAL_PARAMETERS_AS_FUNCTION,
      double *density_as_function,
      double *dynamic_viscosity_as_function,
      double *ebqe_density_as_function,
      double *ebqe_dynamic_viscosity_as_function,
      int USE_SBM,
      // For edge based dissipation
      int ARTIFICIAL_VISCOSITY,
      double *uStar_dMatrix,
      double *vStar_dMatrix,
      double *wStar_dMatrix,
      int numDOFs_1D,
      int offset_u, int offset_v, int offset_w,
      int stride_u, int stride_v, int stride_w,
      int *rowptr_1D, int *colind_1D,
      int *rowptr, int *colind,
      int INT_BY_PARTS_PRESSURE)
  {
    //
    //loop over elements to compute volume integrals and load them into the element Jacobians and global Jacobian
    //
    std::valarray<double> particle_surfaceArea(nParticles), particle_netForces(nParticles * 3 * 3), particle_netMoments(nParticles * 3);
    const int nQuadraturePoints_global(nElements_global * nQuadraturePoints_element);
    //std::set<int> active_velocity_dof;

    for (int eN = 0; eN < nElements_global; eN++)
    {
      register double eps_rho, eps_mu;
      double element_active = 1.0; //value 1 is because it is ibm by default

      register double elementJacobian_p_p[nDOF_test_element][nDOF_trial_element],
          elementJacobian_p_u[nDOF_test_element][nDOF_trial_element],
          elementJacobian_p_v[nDOF_test_element][nDOF_trial_element],
          elementJacobian_p_w[nDOF_test_element][nDOF_trial_element],
          elementJacobian_u_p[nDOF_test_element][nDOF_trial_element],
          elementJacobian_u_u[nDOF_test_element][nDOF_trial_element],
          elementJacobian_u_v[nDOF_test_element][nDOF_trial_element],
          elementJacobian_u_w[nDOF_test_element][nDOF_trial_element],
          elementJacobian_v_p[nDOF_test_element][nDOF_trial_element],
          elementJacobian_v_u[nDOF_test_element][nDOF_trial_element],
          elementJacobian_v_v[nDOF_test_element][nDOF_trial_element],
          elementJacobian_v_w[nDOF_test_element][nDOF_trial_element],
          elementJacobian_w_p[nDOF_test_element][nDOF_trial_element],
          elementJacobian_w_u[nDOF_test_element][nDOF_trial_element],
          elementJacobian_w_v[nDOF_test_element][nDOF_trial_element],
          elementJacobian_w_w[nDOF_test_element][nDOF_trial_element];
      for (int i = 0; i < nDOF_test_element; i++)
        for (int j = 0; j < nDOF_trial_element; j++)
        {
          elementJacobian_p_p[i][j] = 0.0;
          elementJacobian_p_u[i][j] = 0.0;
          elementJacobian_p_v[i][j] = 0.0;
          elementJacobian_p_w[i][j] = 0.0;
          elementJacobian_u_p[i][j] = 0.0;
          elementJacobian_u_u[i][j] = 0.0;
          elementJacobian_u_v[i][j] = 0.0;
          elementJacobian_u_w[i][j] = 0.0;
          elementJacobian_v_p[i][j] = 0.0;
          elementJacobian_v_u[i][j] = 0.0;
          elementJacobian_v_v[i][j] = 0.0;
          elementJacobian_v_w[i][j] = 0.0;
          elementJacobian_w_p[i][j] = 0.0;
          elementJacobian_w_u[i][j] = 0.0;
          elementJacobian_w_v[i][j] = 0.0;
          elementJacobian_w_w[i][j] = 0.0;
        }

      for (int k = 0; k < nQuadraturePoints_element; k++)
      {
        int eN_k = eN * nQuadraturePoints_element + k, //index to a scalar at a quadrature point
            eN_k_nSpace = eN_k * nSpace,
            eN_k_3d = eN_k * 3,
            eN_nDOF_trial_element = eN * nDOF_trial_element; //index to a vector at a quadrature point

        //declare local storage
        register double p = 0.0, u = 0.0, v = 0.0, w = 0.0,
                        grad_p[nSpace], grad_u[nSpace], grad_v[nSpace], grad_w[nSpace],
                        hess_u[nSpace2], hess_v[nSpace2],
                        mom_u_acc = 0.0,
                        dmom_u_acc_u = 0.0,
                        mom_v_acc = 0.0,
                        dmom_v_acc_v = 0.0,
                        mom_w_acc = 0.0,
                        dmom_w_acc_w = 0.0,
                        mass_adv[nSpace],
                        dmass_adv_u[nSpace],
                        dmass_adv_v[nSpace],
                        dmass_adv_w[nSpace],
                        mom_u_adv[nSpace],
                        dmom_u_adv_u[nSpace],
                        dmom_u_adv_v[nSpace],
                        dmom_u_adv_w[nSpace],
                        mom_v_adv[nSpace],
                        dmom_v_adv_u[nSpace],
                        dmom_v_adv_v[nSpace],
                        dmom_v_adv_w[nSpace],
                        mom_w_adv[nSpace],
                        dmom_w_adv_u[nSpace],
                        dmom_w_adv_v[nSpace],
                        dmom_w_adv_w[nSpace],
                        mom_uu_diff_ten[nSpace],
                        mom_vv_diff_ten[nSpace],
                        mom_ww_diff_ten[nSpace],
                        mom_uv_diff_ten[1],
                        mom_uw_diff_ten[1],
                        mom_vu_diff_ten[1],
                        mom_vw_diff_ten[1],
                        mom_wu_diff_ten[1],
                        mom_wv_diff_ten[1],
                        mom_u_source = 0.0,
                        mom_v_source = 0.0,
                        mom_w_source = 0.0,
                        mom_u_ham = 0.0,
                        dmom_u_ham_grad_p[nSpace],
                        dmom_u_ham_grad_u[nSpace],
                        mom_v_ham = 0.0,
                        dmom_v_ham_grad_p[nSpace],
                        dmom_v_ham_grad_v[nSpace],
                        mom_w_ham = 0.0,
                        dmom_w_ham_grad_p[nSpace],
                        dmom_w_ham_grad_w[nSpace],
                        mom_u_acc_t = 0.0,
                        dmom_u_acc_u_t = 0.0,
                        mom_v_acc_t = 0.0,
                        dmom_v_acc_v_t = 0.0,
                        mom_w_acc_t = 0.0,
                        dmom_w_acc_w_t = 0.0,
                        pdeResidual_p = 0.0,
                        pdeResidual_u = 0.0,
                        pdeResidual_v = 0.0,
                        pdeResidual_w = 0.0,
                        dpdeResidual_p_u[nDOF_trial_element], dpdeResidual_p_v[nDOF_trial_element], dpdeResidual_p_w[nDOF_trial_element],
                        dpdeResidual_u_p[nDOF_trial_element], dpdeResidual_u_u[nDOF_trial_element],
                        dpdeResidual_v_p[nDOF_trial_element], dpdeResidual_v_v[nDOF_trial_element],
                        dpdeResidual_w_p[nDOF_trial_element], dpdeResidual_w_w[nDOF_trial_element],
                        Lstar_u_p[nDOF_test_element],
                        Lstar_v_p[nDOF_test_element],
                        Lstar_w_p[nDOF_test_element],
                        Lstar_u_u[nDOF_test_element],
                        Lstar_v_v[nDOF_test_element],
                        Lstar_w_w[nDOF_test_element],
                        Lstar_p_u[nDOF_test_element],
                        Lstar_p_v[nDOF_test_element],
                        Lstar_p_w[nDOF_test_element],
                        subgridError_p = 0.0,
                        subgridError_u = 0.0,
                        subgridError_v = 0.0,
                        subgridError_w = 0.0,
                        dsubgridError_p_u[nDOF_trial_element],
                        dsubgridError_p_v[nDOF_trial_element],
                        dsubgridError_p_w[nDOF_trial_element],
                        dsubgridError_u_p[nDOF_trial_element],
                        dsubgridError_u_u[nDOF_trial_element],
                        dsubgridError_v_p[nDOF_trial_element],
                        dsubgridError_v_v[nDOF_trial_element],
                        dsubgridError_w_p[nDOF_trial_element],
                        dsubgridError_w_w[nDOF_trial_element],
                        tau_p = 0.0, tau_p0 = 0.0, tau_p1 = 0.0,
                        tau_v = 0.0, tau_v0 = 0.0, tau_v1 = 0.0,
                        jac[nSpace * nSpace],
                        jacDet,
                        jacInv[nSpace * nSpace],
                        p_grad_trial[nDOF_trial_element * nSpace], vel_grad_trial[nDOF_trial_element * nSpace],
                        vel_hess_trial[nDOF_trial_element * nSpace2],
                        dV,
                        p_test_dV[nDOF_test_element], vel_test_dV[nDOF_test_element],
                        p_grad_test_dV[nDOF_test_element * nSpace], vel_grad_test_dV[nDOF_test_element * nSpace],
                        x, y, z, xt, yt, zt,
                        //VRANS
            porosity,
                        //meanGrainSize,
            dmom_u_source[nSpace],
                        dmom_v_source[nSpace],
                        dmom_w_source[nSpace],
                        mass_source,
                        //
            G[nSpace * nSpace], G_dd_G, tr_G, h_phi, dmom_adv_star[nSpace], dmom_adv_sge[nSpace];
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
                                    x, y, z);
        ck.calculateH_element(eN,
                              k,
                              nodeDiametersArray,
                              mesh_l2g,
                              mesh_trial_ref,
                              h_phi);
        ck.calculateMappingVelocity_element(eN,
                                            k,
                                            mesh_velocity_dof,
                                            mesh_l2g,
                                            mesh_trial_ref,
                                            xt, yt, zt);
        //xt=0.0;yt=0.0;zt=0.0;
        //std::cout<<"xt "<<xt<<'\t'<<yt<<'\t'<<zt<<std::endl;
        //get the physical integration weight
        dV = fabs(jacDet) * dV_ref[k];
        ck.calculateG(jacInv, G, G_dd_G, tr_G);
        //ck.calculateGScale(G,&normal_phi[eN_k_nSpace],h_phi);

        eps_rho = epsFact_rho * (useMetrics * h_phi + (1.0 - useMetrics) * elementDiameter[eN]);
        eps_mu = epsFact_mu * (useMetrics * h_phi + (1.0 - useMetrics) * elementDiameter[eN]);
        const double particle_eps = particle_epsFact * (useMetrics * h_phi + (1.0 - useMetrics) * elementDiameter[eN]);

        //get the trial function gradients
        /* ck.gradTrialFromRef(&p_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,p_grad_trial); */
        ck.gradTrialFromRef(&vel_grad_trial_ref[k * nDOF_trial_element * nSpace], jacInv, vel_grad_trial);
        // ck.hessTrialFromRef(&vel_hess_trial_ref[k * nDOF_trial_element * nSpace2], jacInv, vel_hess_trial);
        //get the solution
        /* ck.valFromDOF(p_dof,&p_l2g[eN_nDOF_trial_element],&p_trial_ref[k*nDOF_trial_element],p); */
        p = q_p[eN_k];
        ck.valFromDOF(u_dof, &vel_l2g[eN_nDOF_trial_element], &vel_trial_ref[k * nDOF_trial_element], u);
        ck.valFromDOF(v_dof, &vel_l2g[eN_nDOF_trial_element], &vel_trial_ref[k * nDOF_trial_element], v);
        /* ck.valFromDOF(w_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],w); */
        //get the solution gradients
        /* ck.gradFromDOF(p_dof,&p_l2g[eN_nDOF_trial_element],p_grad_trial,grad_p); */
        for (int I = 0; I < nSpace; I++)
          grad_p[I] = q_grad_p[eN_k_nSpace + I];
        ck.gradFromDOF(u_dof, &vel_l2g[eN_nDOF_trial_element], vel_grad_trial, grad_u);
        ck.gradFromDOF(v_dof, &vel_l2g[eN_nDOF_trial_element], vel_grad_trial, grad_v);
        /* ck.gradFromDOF(w_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_w); */
        //precalculate test function products with integration weights
        for (int j = 0; j < nDOF_trial_element; j++)
        {
          /* p_test_dV[j] = p_test_ref[k*nDOF_trial_element+j]*dV; */
          vel_test_dV[j] = vel_test_ref[k * nDOF_trial_element + j] * dV;
          for (int I = 0; I < nSpace; I++)
          {
            /* p_grad_test_dV[j*nSpace+I]   = p_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin */
            vel_grad_test_dV[j * nSpace + I] = vel_grad_trial[j * nSpace + I] * dV; //cek warning won't work for Petrov-Galerkin}
          }
        }
        //cek hack
        double div_mesh_velocity = 0.0;
        //
        //VRANS
        porosity = 1.0 - q_vos[eN_k];
        //
        //
        //calculate pde coefficients and derivatives at quadrature points
        //
        double distance_to_omega_solid = phi_solid[eN_k]; //computed in getResidual
        double eddy_viscosity(0.), rhoSave, nuSave;       //not really interested in saving eddy_viscosity in jacobian
        evaluateCoefficients(eps_rho,
                             eps_mu,
                             particle_eps,
                             sigma,
                             rho_0,
                             nu_0,
                             rho_1,
                             nu_1,
                             elementDiameter[eN],
                             smagorinskyConstant,
                             turbulenceClosureModel,
                             g,
                             useVF,
                             vf[eN_k],
                             phi[eN_k],
                             &normal_phi[eN_k_nSpace],
                             distance_to_omega_solid,
                             kappa_phi[eN_k],
                             //VRANS
                             porosity,
                             //
                             p,
                             grad_p,
                             grad_u,
                             grad_v,
                             grad_w,
                             u,
                             v,
                             w,
                             q_velocity_sge[eN_k_nSpace + 0],
                             q_velocity_sge[eN_k_nSpace + 1],
                             q_velocity_sge[eN_k_nSpace + 1], //hack, shouldn't be used
                             eddy_viscosity,
                             mom_u_acc,
                             dmom_u_acc_u,
                             mom_v_acc,
                             dmom_v_acc_v,
                             mom_w_acc,
                             dmom_w_acc_w,
                             mass_adv,
                             dmass_adv_u,
                             dmass_adv_v,
                             dmass_adv_w,
                             mom_u_adv,
                             dmom_u_adv_u,
                             dmom_u_adv_v,
                             dmom_u_adv_w,
                             mom_v_adv,
                             dmom_v_adv_u,
                             dmom_v_adv_v,
                             dmom_v_adv_w,
                             mom_w_adv,
                             dmom_w_adv_u,
                             dmom_w_adv_v,
                             dmom_w_adv_w,
                             mom_uu_diff_ten,
                             mom_vv_diff_ten,
                             mom_ww_diff_ten,
                             mom_uv_diff_ten,
                             mom_uw_diff_ten,
                             mom_vu_diff_ten,
                             mom_vw_diff_ten,
                             mom_wu_diff_ten,
                             mom_wv_diff_ten,
                             mom_u_source,
                             mom_v_source,
                             mom_w_source,
                             mom_u_ham,
                             dmom_u_ham_grad_p,
                             dmom_u_ham_grad_u,
                             mom_v_ham,
                             dmom_v_ham_grad_p,
                             dmom_v_ham_grad_v,
                             mom_w_ham,
                             dmom_w_ham_grad_p,
                             dmom_w_ham_grad_w,
                             rhoSave,
                             nuSave,
                             KILL_PRESSURE_TERM,
                             0,
                             0., // mql: the force term doesn't play a role in the Jacobian
                             0.,
                             0.,
                             MATERIAL_PARAMETERS_AS_FUNCTION,
                             density_as_function[eN_k],
                             dynamic_viscosity_as_function[eN_k],
                             USE_SBM,
                             x, y, z,
                             use_ball_as_particle,
                             ball_center,
                             ball_radius,
                             ball_velocity,
                             ball_angular_velocity,
                             INT_BY_PARTS_PRESSURE);
        //VRANS
        mass_source = q_mass_source[eN_k];
        for (int I = 0; I < nSpace; I++)
        {
          dmom_u_source[I] = 0.0;
          dmom_v_source[I] = 0.0;
          dmom_w_source[I] = 0.0;
        }

        double C_particles = 0.0;

        //
        //calculate time derivatives
        //
        ck.bdf(alphaBDF,
               q_mom_u_acc_beta_bdf[eN_k] * q_dV_last[eN_k] / dV,
               mom_u_acc,
               dmom_u_acc_u,
               mom_u_acc_t,
               dmom_u_acc_u_t);
        ck.bdf(alphaBDF,
               q_mom_v_acc_beta_bdf[eN_k] * q_dV_last[eN_k] / dV,
               mom_v_acc,
               dmom_v_acc_v,
               mom_v_acc_t,
               dmom_v_acc_v_t);
        /* ck.bdf(alphaBDF, */
        /*           q_mom_w_acc_beta_bdf[eN_k]*q_dV_last[eN_k]/dV, */
        /*           mom_w_acc, */
        /*           dmom_w_acc_w, */
        /*           mom_w_acc_t, */
        /*           dmom_w_acc_w_t); */
        //
        //calculate subgrid error contribution to the Jacobian (strong residual, adjoint, jacobian of strong residual)

        mom_u_acc_t *= dmom_u_acc_u; //multiply by porosity*rho. mql. CHECK.
        mom_v_acc_t *= dmom_v_acc_v;

        //
        dmom_adv_sge[0] = dmom_u_acc_u * (q_velocity_sge[eN_k_nSpace + 0] - MOVING_DOMAIN * xt);
        dmom_adv_sge[1] = dmom_u_acc_u * (q_velocity_sge[eN_k_nSpace + 1] - MOVING_DOMAIN * yt);
        /* dmom_adv_sge[2] = dmom_u_acc_u*(q_velocity_sge[eN_k_nSpace+2] - MOVING_DOMAIN*zt); */
        //cek todo add RBLES terms consistent to residual modifications or ignore the partials w.r.t the additional RBLES terms
        const double lambda = 1.0;
        for (int i = 0; i < nDOF_test_element; i++)
        {
          register int i_nSpace = i * nSpace;
          for (int j = 0; j < nDOF_trial_element; j++)
          {
            register int j_nSpace = j * nSpace;

            /* elementJacobian_p_p[i][j] += ck.SubgridErrorJacobian(dsubgridError_u_p[j],Lstar_u_p[i]) +  */
            /*   ck.SubgridErrorJacobian(dsubgridError_v_p[j],Lstar_v_p[i]);// +  */
            /*   /\* ck.SubgridErrorJacobian(dsubgridError_w_p[j],Lstar_w_p[i]);  *\/ */

            /* elementJacobian_p_u[i][j] += ck.AdvectionJacobian_weak(dmass_adv_u,vel_trial_ref[k*nDOF_trial_element+j],&p_grad_test_dV[i_nSpace]) +  */
            /*   ck.SubgridErrorJacobian(dsubgridError_u_u[j],Lstar_u_p[i]);  */
            /* elementJacobian_p_v[i][j] += ck.AdvectionJacobian_weak(dmass_adv_v,vel_trial_ref[k*nDOF_trial_element+j],&p_grad_test_dV[i_nSpace]) +  */
            /*   ck.SubgridErrorJacobian(dsubgridError_v_v[j],Lstar_v_p[i]);  */
            /* elementJacobian_p_w[i][j] += ck.AdvectionJacobian_weak(dmass_adv_w,vel_trial_ref[k*nDOF_trial_element+j],&p_grad_test_dV[i_nSpace]) +  */
            /*      ck.SubgridErrorJacobian(dsubgridError_w_w[j],Lstar_w_p[i]);  */

            /* elementJacobian_u_p[i][j] += ck.HamiltonianJacobian_weak(dmom_u_ham_grad_p,&p_grad_trial[j_nSpace],vel_test_dV[i]) +  */
            /*   ck.SubgridErrorJacobian(dsubgridError_u_p[j],Lstar_u_u[i]);  */
            elementJacobian_u_u[i][j] +=
                ck.MassJacobian_weak(dmom_u_acc_u_t, vel_trial_ref[k * nDOF_trial_element + j], vel_test_dV[i])
                + ck.HamiltonianJacobian_weak(dmom_u_ham_grad_u, &vel_grad_trial[j_nSpace], vel_test_dV[i])
                + ck.AdvectionJacobian_weak(dmom_u_adv_u, vel_trial_ref[k * nDOF_trial_element + j], &vel_grad_test_dV[i_nSpace])
                + ck.SimpleDiffusionJacobian_weak(sdInfo_u_u_rowptr, sdInfo_u_u_colind, mom_uu_diff_ten, &vel_grad_trial[j_nSpace], &vel_grad_test_dV[i_nSpace])
                // + mom_uu_diff_ten[0] * vel_grad_trial[j_nSpace + 0] * vel_grad_test_dV[i_nSpace + 0]
                // + mom_uu_diff_ten[1] * vel_grad_trial[j_nSpace + 1] * vel_grad_test_dV[i_nSpace + 1]
                - lambda * vel_grad_trial[j_nSpace + 0] * (vel_grad_test_dV[i_nSpace + 0] + vel_grad_test_dV[i_nSpace + 1])
                //VRANS
                // + ck.ReactionJacobian_weak(dmom_u_source[0], vel_trial_ref[k * nDOF_trial_element + j], vel_test_dV[i])
                //
                //ck.SubgridErrorJacobian(dsubgridError_p_u[j],Lstar_p_u[i]) +
                // USE_SUPG * ck.SubgridErrorJacobian(dsubgridError_u_u[j], Lstar_u_u[i]) +
                // ck.NumericalDiffusionJacobian(q_numDiff_u_last[eN_k], &vel_grad_trial[j_nSpace], &vel_grad_test_dV[i_nSpace])
                ;

            elementJacobian_u_v[i][j] +=
                ck.AdvectionJacobian_weak(dmom_u_adv_v, vel_trial_ref[k * nDOF_trial_element + j], &vel_grad_test_dV[i_nSpace])
                + ck.SimpleDiffusionJacobian_weak(sdInfo_u_v_rowptr, sdInfo_u_v_colind, mom_uv_diff_ten, &vel_grad_trial[j_nSpace], &vel_grad_test_dV[i_nSpace])
                - lambda * vel_grad_trial[j_nSpace + 1] * (vel_grad_test_dV[i_nSpace + 0] + vel_grad_test_dV[i_nSpace + 1])
                //VRANS
                // + ck.ReactionJacobian_weak(dmom_u_source[1], vel_trial_ref[k * nDOF_trial_element + j], vel_test_dV[i])
                ;
            elementJacobian_v_u[i][j] +=
                ck.AdvectionJacobian_weak(dmom_v_adv_u, vel_trial_ref[k * nDOF_trial_element + j], &vel_grad_test_dV[i_nSpace])
                + ck.SimpleDiffusionJacobian_weak(sdInfo_v_u_rowptr, sdInfo_v_u_colind, mom_vu_diff_ten, &vel_grad_trial[j_nSpace], &vel_grad_test_dV[i_nSpace])
                - lambda * vel_grad_trial[j_nSpace + 0] * (vel_grad_test_dV[i_nSpace + 0] + vel_grad_test_dV[i_nSpace + 1])
                //VRANS
                // + ck.ReactionJacobian_weak(dmom_v_source[0], vel_trial_ref[k * nDOF_trial_element + j], vel_test_dV[i])
                ;
            elementJacobian_v_v[i][j] +=
                ck.MassJacobian_weak(dmom_v_acc_v_t, vel_trial_ref[k * nDOF_trial_element + j], vel_test_dV[i])
                + ck.HamiltonianJacobian_weak(dmom_v_ham_grad_v, &vel_grad_trial[j_nSpace], vel_test_dV[i])
                + ck.AdvectionJacobian_weak(dmom_v_adv_v, vel_trial_ref[k * nDOF_trial_element + j], &vel_grad_test_dV[i_nSpace])
                + ck.SimpleDiffusionJacobian_weak(sdInfo_v_v_rowptr, sdInfo_v_v_colind, mom_vv_diff_ten, &vel_grad_trial[j_nSpace], &vel_grad_test_dV[i_nSpace])
                // + mom_vv_diff_ten[0] * vel_grad_trial[j_nSpace + 0] * vel_grad_test_dV[i_nSpace + 0]
                // + mom_vv_diff_ten[1] * vel_grad_trial[j_nSpace + 1] * vel_grad_test_dV[i_nSpace + 1]
                - lambda * vel_grad_trial[j_nSpace + 1] * (vel_grad_test_dV[i_nSpace + 0] + vel_grad_test_dV[i_nSpace + 1])
                // + ck.ReactionJacobian_weak(dmom_v_source[1], vel_trial_ref[k * nDOF_trial_element + j], vel_test_dV[i])
                // USE_SUPG * ck.SubgridErrorJacobian(dsubgridError_v_v[j], Lstar_v_v[i]) +
                // ck.NumericalDiffusionJacobian(q_numDiff_v_last[eN_k], &vel_grad_trial[j_nSpace], &vel_grad_test_dV[i_nSpace])
                ;

          } //j
        }   //i
      }     //k
      //
      //load into element Jacobian into global Jacobian
      //
      for (int i = 0; i < nDOF_test_element; i++)
      {
        register int eN_i = eN * nDOF_test_element + i;
        for (int j = 0; j < nDOF_trial_element; j++)
        {
          register int eN_i_j = eN_i * nDOF_trial_element + j;
          /* globalJacobian[csrRowIndeces_p_p[eN_i] + csrColumnOffsets_p_p[eN_i_j]] += elementJacobian_p_p[i][j]; */
          /* globalJacobian[csrRowIndeces_p_u[eN_i] + csrColumnOffsets_p_u[eN_i_j]] += elementJacobian_p_u[i][j]; */
          /* globalJacobian[csrRowIndeces_p_v[eN_i] + csrColumnOffsets_p_v[eN_i_j]] += elementJacobian_p_v[i][j]; */
          /* globalJacobian[csrRowIndeces_p_w[eN_i] + csrColumnOffsets_p_w[eN_i_j]] += elementJacobian_p_w[i][j]; */

          /* globalJacobian[csrRowIndeces_u_p[eN_i] + csrColumnOffsets_u_p[eN_i_j]] += elementJacobian_u_p[i][j]; */
          globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_u_u[eN_i_j]] += element_active * elementJacobian_u_u[i][j];
          globalJacobian[csrRowIndeces_u_v[eN_i] + csrColumnOffsets_u_v[eN_i_j]] += element_active * elementJacobian_u_v[i][j];
          /* globalJacobian[csrRowIndeces_u_w[eN_i] + csrColumnOffsets_u_w[eN_i_j]] += elementJacobian_u_w[i][j]; */

          /* globalJacobian[csrRowIndeces_v_p[eN_i] + csrColumnOffsets_v_p[eN_i_j]] += elementJacobian_v_p[i][j]; */
          globalJacobian[csrRowIndeces_v_u[eN_i] + csrColumnOffsets_v_u[eN_i_j]] += element_active * elementJacobian_v_u[i][j];
          globalJacobian[csrRowIndeces_v_v[eN_i] + csrColumnOffsets_v_v[eN_i_j]] += element_active * elementJacobian_v_v[i][j];
          /* globalJacobian[csrRowIndeces_v_w[eN_i] + csrColumnOffsets_v_w[eN_i_j]] += elementJacobian_v_w[i][j]; */

          /* globalJacobian[csrRowIndeces_w_p[eN_i] + csrColumnOffsets_w_p[eN_i_j]] += elementJacobian_w_p[i][j]; */
          /* globalJacobian[csrRowIndeces_w_u[eN_i] + csrColumnOffsets_w_u[eN_i_j]] += elementJacobian_w_u[i][j]; */
          /* globalJacobian[csrRowIndeces_w_v[eN_i] + csrColumnOffsets_w_v[eN_i_j]] += elementJacobian_w_v[i][j]; */
          /* globalJacobian[csrRowIndeces_w_w[eN_i] + csrColumnOffsets_w_w[eN_i_j]] += elementJacobian_w_w[i][j]; */
        } //j
      }   //i
    }     //elements

  }           //computeJacobian

  void calculateVelocityAverage(int nExteriorElementBoundaries_global,
                                int *exteriorElementBoundariesArray,
                                int nInteriorElementBoundaries_global,
                                int *interiorElementBoundariesArray,
                                int *elementBoundaryElementsArray,
                                int *elementBoundaryLocalElementBoundariesArray,
                                double *mesh_dof,
                                double *mesh_velocity_dof,
                                double MOVING_DOMAIN, //0 or 1
                                int *mesh_l2g,
                                double *mesh_trial_trace_ref,
                                double *mesh_grad_trial_trace_ref,
                                double *normal_ref,
                                double *boundaryJac_ref,
                                int *vel_l2g,
                                double *u_dof,
                                double *v_dof,
                                double *w_dof,
                                double *vos_dof,
                                double *vel_trial_trace_ref,
                                double *ebqe_velocity,
                                double *velocityAverage)
  {
  }

  void getBoundaryDOFs( //element
      double *mesh_dof,
      int *mesh_l2g,
      double *mesh_trial_trace_ref,
      double *mesh_grad_trial_trace_ref,
      double *dS_ref,
      double *vel_test_trace_ref,
      double *normal_ref,
      double *boundaryJac_ref,
      int *vel_l2g,
      int nExteriorElementBoundaries_global,
      int *exteriorElementBoundariesArray,
      int *elementBoundaryElementsArray,
      int *elementBoundaryLocalElementBoundariesArray,
      double *isBoundary_1D)
  {
  }
}; //RANS3PF2D

inline cppRANS3PF2D_base *newRANS3PF2D(int nSpaceIn,
                                       int nQuadraturePoints_elementIn,
                                       int nDOF_mesh_trial_elementIn,
                                       int nDOF_trial_elementIn,
                                       int nDOF_test_elementIn,
                                       int nQuadraturePoints_elementBoundaryIn,
                                       int CompKernelFlag,
                                       double aDarcy,
                                       double betaForch,
                                       double grain,
                                       double packFraction,
                                       double packMargin,
                                       double maxFraction,
                                       double frFraction,
                                       double sigmaC,
                                       double C3e,
                                       double C4e,
                                       double eR,
                                       double fContact,
                                       double mContact,
                                       double nContact,
                                       double angFriction, double vos_limiter, double mu_fr_limiter)
{
  cppRANS3PF2D_base *rvalue = proteus::chooseAndAllocateDiscretization2D<cppRANS3PF2D_base, cppRANS3PF2D, CompKernel>(nSpaceIn,
                                                                                                                      nQuadraturePoints_elementIn,
                                                                                                                      nDOF_mesh_trial_elementIn,
                                                                                                                      nDOF_trial_elementIn,
                                                                                                                      nDOF_test_elementIn,
                                                                                                                      nQuadraturePoints_elementBoundaryIn,
                                                                                                                      CompKernelFlag);
  rvalue->setSedClosure(aDarcy,
                        betaForch,
                        grain,
                        packFraction,
                        packMargin,
                        maxFraction,
                        frFraction,
                        sigmaC,
                        C3e,
                        C4e,
                        eR,
                        fContact,
                        mContact,
                        nContact,
                        angFriction, vos_limiter, mu_fr_limiter);
  return rvalue;
}
} // namespace proteus

#endif

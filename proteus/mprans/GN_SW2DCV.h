#ifndef GN_SW2DCV_H
#define GN_SW2DCV_H
#include "CompKernel.h"
#include "ModelFactory.h"
#include <assert.h>
#include <cmath>
#include <iostream>
#include <valarray>
#include <vector>

// cek todo
// 2. Get stabilization right
// 3. Add Riemann solvers for external flux
// 4. Add Riemann solvers for internal flux and DG terms
// 5. Try other choices of variables h,hu,hv, Bova-Carey symmetrization?

#define GLOBAL_FCT 0
#define POWER_SMOOTHNESS_INDICATOR 0
#define VEL_FIX_POWER 2.
#define REESTIMATE_MAX_EDGE_BASED_CFL 0
#define LAMBDA_MGN 1
#define IF_BOTH_GAMMA_BRANCHES 0

namespace proteus {
/* FOR mGN, maybe we don't need a lot of these inline functions for now after
 * first two that I defined */

// for mGN stuff "max" wave speeds. See section 4.4 of first mGN paper by
// Guermond, Popov, Tovar, Kees for formulas
inline double GN_nu1(const double &g, const double &hL, const double &uL,
                     const double &etaL, const double &meshSizeL) {

  double augL = LAMBDA_MGN / (3. * meshSizeL) * (6. * hL + 12. * (hL - etaL));

  if (etaL >= hL) {
    augL = LAMBDA_MGN / (3. * meshSizeL) * (6. * hL);
  }
  augL = augL * std::pow(meshSizeL / fmax(meshSizeL, hL), 2.0);

  return uL - sqrt(g * hL) * sqrt(1. + augL);
}
inline double GN_nu3(const double &g, const double &hR, const double &uR,
                     const double &etaR, const double &meshSizeR) {
  // See section 4.4 of first mGN paper by Guermond, Popov, Tovar, Kees for
  // formulas
  double augR = LAMBDA_MGN / (3. * meshSizeR) * (6. * hR + 12. * (hR - etaR));

  if (etaR >= hR) {
    augR = LAMBDA_MGN / (3. * meshSizeR) * (6. * hR);
  }
  augR = augR * std::pow(meshSizeR / fmax(meshSizeR, hR), 2.0);

  return uR + sqrt(g * hR) * sqrt(1 + augR);
}

// FOR CELL BASED ENTROPY VISCOSITY
inline double ENTROPY(const double &g, const double &h, const double &hu,
                      const double &hv, const double &z,
                      const double &one_over_hReg) {
  return 0.5 *
         (g * h * h + one_over_hReg * (hu * hu + hv * hv) + 2. * g * h * z);
}
inline double DENTROPY_DH(const double &g, const double &h, const double &hu,
                          const double &hv, const double &z,
                          const double &one_over_hReg) {
  return g * h - 0.5 * (hu * hu + hv * hv) * std::pow(one_over_hReg, 2) + g * z;
}
inline double DENTROPY_DHU(const double &g, const double &h, const double &hu,
                           const double &hv, const double &z,
                           const double &one_over_hReg) {
  return hu * one_over_hReg;
}
inline double DENTROPY_DHV(const double &g, const double &h, const double &hu,
                           const double &hv, const double &z,
                           const double &one_over_hReg) {
  return hv * one_over_hReg;
}
inline double ENTROPY_FLUX1(const double &g, const double &h, const double &hu,
                            const double &hv, const double &z,
                            const double &one_over_hReg) {
  return (ENTROPY(g, h, hu, hv, z, one_over_hReg) + 0.5 * g * h * h +
          g * h * z) *
         hu * one_over_hReg;
}
inline double ENTROPY_FLUX2(const double &g, const double &h, const double &hu,
                            const double &hv, const double &z,
                            const double &one_over_hReg) {
  return (ENTROPY(g, h, hu, hv, z, one_over_hReg) + 0.5 * g * h * h +
          g * h * z) *
         hv * one_over_hReg;
}
// FOR ESTIMATING MAX WAVE SPEEDS, GIVE REFERENCE
inline double f(const double &g, const double &h, const double &hZ) {
  return ((h <= hZ) ? 2. * (sqrt(g * h) - sqrt(g * hZ))
                    : (h - hZ) * sqrt(0.5 * g * (h + hZ) / h / hZ));
}
inline double phi(const double &g, const double &h, const double &hL,
                  const double &hR, const double &uL, const double &uR) {
  return (f(g, h, hL) + f(g, h, hR) + uR - uL);
}
inline double fp(const double &g, const double &h, const double &hZ) {
  return ((h <= hZ)
              ? sqrt(g / h)
              : g * (2 * h * h + h * hZ + hZ * hZ) /
                    (2 * sqrt(2 * g) * h * h * hZ * sqrt(1 / h + 1 / hZ)));
}
inline double phip(const double &g, const double &h, const double &hL,
                   const double &hR) {
  return (fp(g, h, hL) + fp(g, h, hR));
}
inline double nu1(const double &g, const double &hStar, const double &hL,
                  const double &uL) {
  return (uL - sqrt(g * hL) * sqrt((1 + fmax((hStar - hL) / 2 / hL, 0.0)) *
                                   (1 + fmax((hStar - hL) / hL, 0.))));
}
inline double nu3(const double &g, const double &hStar, const double &hR,
                  const double &uR) {
  return (uR + sqrt(g * hR) * sqrt((1 + fmax((hStar - hR) / 2 / hR, 0.0)) *
                                   (1 + fmax((hStar - hR) / hR, 0.))));
}
inline double phiDiff(const double &g, const double &h1k, const double &h2k,
                      const double &hL, const double &hR, const double &uL,
                      const double &uR) {
  return ((phi(g, h2k, hL, hR, uL, uR) - phi(g, h1k, hL, hR, uL, uR)) /
          (h2k - h1k));
}
inline double phiDDiff1(const double &g, const double &h1k, const double &h2k,
                        const double &hL, const double &hR, const double &uL,
                        const double &uR) {
  return ((phiDiff(g, h1k, h2k, hL, hR, uL, uR) - phip(g, h1k, hL, hR)) /
          (h2k - h1k));
}
inline double phiDDiff2(const double &g, const double &h1k, const double &h2k,
                        const double &hL, const double &hR, const double &uL,
                        const double &uR) {
  return ((phip(g, h2k, hL, hR) - phiDiff(g, h1k, h2k, hL, hR, uL, uR)) /
          (h2k - h1k));
}
inline double hStarLFromQuadPhiFromAbove(const double &g, const double &hStarL,
                                         const double &hStarR, const double &hL,
                                         const double &hR, const double &uL,
                                         const double &uR) {
  return (hStarL -
          2 * phi(g, hStarL, hL, hR, uL, uR) /
              (phip(g, hStarL, hL, hR) +
               sqrt(std::pow(phip(g, hStarL, hL, hR), 2) -
                    4 * phi(g, hStarL, hL, hR, uL, uR) *
                        phiDDiff1(g, hStarL, hStarR, hL, hR, uL, uR))));
}
inline double hStarRFromQuadPhiFromBelow(const double &g, const double &hStarL,
                                         const double &hStarR, const double &hL,
                                         const double &hR, const double &uL,
                                         const double &uR) {
  return (hStarR -
          2 * phi(g, hStarR, hL, hR, uL, uR) /
              (phip(g, hStarR, hL, hR) +
               sqrt(std::pow(phip(g, hStarR, hL, hR), 2) -
                    4 * phi(g, hStarR, hL, hR, uL, uR) *
                        phiDDiff2(g, hStarL, hStarR, hL, hR, uL, uR))));
}
} // namespace proteus

namespace proteus {
class GN_SW2DCV_base {
public:
  std::valarray<double> Rneg, Rpos, Rneg_heta, Rpos_heta, hLow, huLow, hvLow,
      hetaLow, hwLow, Kmax;
  virtual ~GN_SW2DCV_base() {}
  virtual void
  FCTStep(double dt,
          int NNZ,     // number on non-zero entries on sparsity pattern
          int numDOFs, // number of DOFs
          double *lumped_mass_matrix, // lumped mass matrix (as vector)
          double *h_old,              // DOFs of solution at last stage
          double *hu_old, double *hv_old, double *heta_old, double *hw_old,
          double *b_dof,
          double *high_order_hnp1, // DOFs of high order solution at tnp1
          double *high_order_hunp1, double *high_order_hvnp1,
          double *high_order_hetanp1, double *high_order_hwnp1,
          double *extendedSourceTerm_hu, double *extendedSourceTerm_hv,
          double *extendedSourceTerm_heta, double *extendedSourceTerm_hw,
          double *limited_hnp1, double *limited_hunp1, double *limited_hvnp1,
          double *limited_hetanp1, double *limited_hwnp1,
          int *csrRowIndeces_DofLoops,    // csr row indeces
          int *csrColumnOffsets_DofLoops, // csr column offsets
          double *MassMatrix,             // mass matrix
          double *dH_minus_dL, double *muH_minus_muL, double hEps, double *hReg,
          int LUMPED_MASS_MATRIX,
          // LOCAL LIMITING
          double *dLow, double *hBT, double *huBT, double *hvBT, double *hetaBT,
          double *hwBT) = 0;
  virtual void convexLimiting(
      double dt,
      int NNZ,     // number on non-zero entries on sparsity pattern
      int numDOFs, // number of DOFs
      double *lumped_mass_matrix, // lumped mass matrix (as vector)
      double *h_old,              // DOFs of solution at last stage
      double *hu_old, double *hv_old, double *heta_old, double *hw_old,
      double *b_dof,
      double *high_order_hnp1, // DOFs of high order solution at tnp1
      double *high_order_hunp1, double *high_order_hvnp1,
      double *high_order_hetanp1, double *high_order_hwnp1,
      double *extendedSourceTerm_hu, double *extendedSourceTerm_hv,
      double *extendedSourceTerm_heta, double *extendedSourceTerm_hw,
      double *limited_hnp1, double *limited_hunp1, double *limited_hvnp1,
      double *limited_hetanp1, double *limited_hwnp1,
      int *csrRowIndeces_DofLoops,    // csr row indeces
      int *csrColumnOffsets_DofLoops, // csr column offsets
      double *MassMatrix,             // mass matrix
      double *dH_minus_dL, double *muH_minus_muL, double hEps, double *hReg,
      int LUMPED_MASS_MATRIX,
      // LOCAL LIMITING
      double *dLow, double *hBT, double *huBT, double *hvBT, double *hetaBT,
      double *hwBT, double *new_SourceTerm_hu, double *new_SourceTerm_hv,
      double *new_SourceTerm_heta, double *new_SourceTerm_hw) = 0;
  virtual double calculateEdgeBasedCFL(
      double g,
      int numDOFsPerEqn,          // number of DOFs
      double *lumped_mass_matrix, // lumped mass matrix (as vector)
      double *h_old,              // DOFs of solution at last stage
      double *hu_old, double *hv_old, double *heta_old, double *b_dof,
      int *csrRowIndeces_DofLoops,    // csr row indeces
      int *csrColumnOffsets_DofLoops, // csr column offsets
      double hEps, double *hReg, double *Cx, double *Cy, double *CTx,
      double *CTy, double *dLow, double run_cfl, double *edge_based_cfl,
      int debug) = 0;
  virtual void calculateResidual( // last EDGE BASED version
      double *mesh_trial_ref, double *mesh_grad_trial_ref, double *mesh_dof,
      double *mesh_velocity_dof,
      double MOVING_DOMAIN, // 0 or 1
      int *mesh_l2g, double *dV_ref, double *h_trial_ref,
      double *h_grad_trial_ref, double *h_test_ref, double *h_grad_test_ref,
      double *vel_trial_ref, double *vel_grad_trial_ref, double *vel_test_ref,
      double *vel_grad_test_ref,
      // element boundary
      double *mesh_trial_trace_ref, double *mesh_grad_trial_trace_ref,
      double *dS_ref, double *h_trial_trace_ref, double *h_grad_trial_trace_ref,
      double *h_test_trace_ref, double *h_grad_test_trace_ref,
      double *vel_trial_trace_ref, double *vel_grad_trial_trace_ref,
      double *vel_test_trace_ref, double *vel_grad_test_trace_ref,
      double *normal_ref, double *boundaryJac_ref,
      // physics
      double *elementDiameter, int nElements_global, double useRBLES,
      double useMetrics, double alphaBDF, double nu, double g, int *h_l2g,
      int *vel_l2g, double *h_dof_old, double *hu_dof_old, double *hv_dof_old,
      double *heta_dof_old, double *hw_dof_old, double *b_dof, double *h_dof,
      double *hu_dof, double *hv_dof, double *heta_dof, double *hw_dof,
      double *h_dof_sge, double *hu_dof_sge, double *hv_dof_sge,
      double *heta_dof_sge, double *hw_dof_sge, double *q_mass_acc,
      double *q_mom_hu_acc, double *q_mom_hv_acc, double *q_mass_adv,
      double *q_mass_acc_beta_bdf, double *q_mom_hu_acc_beta_bdf,
      double *q_mom_hv_acc_beta_bdf, double *q_cfl, int *sdInfo_hu_hu_rowptr,
      int *sdInfo_hu_hu_colind, int *sdInfo_hu_hv_rowptr,
      int *sdInfo_hu_hv_colind, int *sdInfo_hv_hv_rowptr,
      int *sdInfo_hv_hv_colind, int *sdInfo_hv_hu_rowptr,
      int *sdInfo_hv_hu_colind, int offset_h, int offset_hu, int offset_hv,
      int offset_heta, int offset_hw, int stride_h, int stride_hu,
      int stride_hv, int stride_heta, int stride_hw, double *globalResidual,
      int nExteriorElementBoundaries_global,
      int *exteriorElementBoundariesArray, int *elementBoundaryElementsArray,
      int *elementBoundaryLocalElementBoundariesArray, int *isDOFBoundary_h,
      int *isDOFBoundary_hu, int *isDOFBoundary_hv,
      int *isAdvectiveFluxBoundary_h, int *isAdvectiveFluxBoundary_hu,
      int *isAdvectiveFluxBoundary_hv, int *isDiffusiveFluxBoundary_hu,
      int *isDiffusiveFluxBoundary_hv, double *ebqe_bc_h_ext,
      double *ebqe_bc_flux_mass_ext, double *ebqe_bc_flux_mom_hu_adv_ext,
      double *ebqe_bc_flux_mom_hv_adv_ext, double *ebqe_bc_hu_ext,
      double *ebqe_bc_flux_hu_diff_ext, double *ebqe_penalty_ext,
      double *ebqe_bc_hv_ext, double *ebqe_bc_flux_hv_diff_ext,
      double *q_velocity, double *ebqe_velocity, double *flux,
      double *elementResidual_h,
      // C matrices
      double *Cx, double *Cy, double *CTx, double *CTy,
      // PARAMETERS FOR EDGE BASED STABILIZATION
      int numDOFsPerEqn, int NNZ, int *csrRowIndeces_DofLoops,
      int *csrColumnOffsets_DofLoops,
      // LUMPED MASS MATRIX
      double *lumped_mass_matrix, double cfl_run, double hEps, double *hReg,
      // SAVE SOLUTION (mql)
      double *hnp1_at_quad_point, double *hunp1_at_quad_point,
      double *hvnp1_at_quad_point, double *hetanp1_at_quad_point,
      double *hwnp1_at_quad_point,
      // TO COMPUTE LOW ORDER
      double *extendedSourceTerm_hu, double *extendedSourceTerm_hv,
      double *extendedSourceTerm_heta, double *extendedSourceTerm_hw,
      // FOR FCT
      double *dH_minus_dL, double *muH_minus_muL, double cE,
      int LUMPED_MASS_MATRIX, double dt, int LINEAR_FRICTION, double mannings,
      // Quant of interests
      double *quantDOFs, int SECOND_CALL_CALCULATE_RESIDUAL,
      // NORMAL COMPONENTS
      int COMPUTE_NORMALS, double *normalx, double *normaly, double *dLow,
      double *hBT, double *huBT, double *hvBT, double *hetaBT, double *hwBT,
      int lstage, double *new_SourceTerm_hu, double *new_SourceTerm_hv,
      double *new_SourceTerm_heta, double *new_SourceTerm_hw) = 0;
  virtual void calculateMassMatrix( // element
      double *mesh_trial_ref, double *mesh_grad_trial_ref, double *mesh_dof,
      double *mesh_velocity_dof, double MOVING_DOMAIN, int *mesh_l2g,
      double *dV_ref, double *h_trial_ref, double *h_grad_trial_ref,
      double *h_test_ref, double *h_grad_test_ref, double *vel_trial_ref,
      double *vel_grad_trial_ref, double *vel_test_ref,
      double *vel_grad_test_ref,
      // element boundary
      double *mesh_trial_trace_ref, double *mesh_grad_trial_trace_ref,
      double *dS_ref, double *h_trial_trace_ref, double *h_grad_trial_trace_ref,
      double *h_test_trace_ref, double *h_grad_test_trace_ref,
      double *vel_trial_trace_ref, double *vel_grad_trial_trace_ref,
      double *vel_test_trace_ref, double *vel_grad_test_trace_ref,
      double *normal_ref, double *boundaryJac_ref,
      // physics
      double *elementDiameter, int nElements_global, double useRBLES,
      double useMetrics, double alphaBDF, double nu, double g, int *h_l2g,
      int *vel_l2g, double *b_dof, double *h_dof, double *hu_dof,
      double *hv_dof, double *heta_dof, double *hw_dof, double *h_dof_sge,
      double *hu_dof_sge, double *hv_dof_sge, double *heta_dof_sge,
      double *hw_dof_sge, double *q_mass_acc_beta_bdf,
      double *q_mom_hu_acc_beta_bdf, double *q_mom_hv_acc_beta_bdf,
      double *q_cfl, int *sdInfo_hu_hu_rowptr, int *sdInfo_hu_hu_colind,
      int *sdInfo_hu_hv_rowptr, int *sdInfo_hu_hv_colind,
      int *sdInfo_hv_hv_rowptr, int *sdInfo_hv_hv_colind,
      int *sdInfo_hv_hu_rowptr, int *sdInfo_hv_hu_colind,
      // h
      int *csrRowIndeces_h_h, int *csrColumnOffsets_h_h,
      int *csrRowIndeces_h_hu, int *csrColumnOffsets_h_hu,
      int *csrRowIndeces_h_hv, int *csrColumnOffsets_h_hv,
      int *csrRowIndeces_h_heta, int *csrColumnOffsets_h_heta,
      int *csrRowIndeces_h_hw, int *csrColumnOffsets_h_hw,
      // hu
      int *csrRowIndeces_hu_h, int *csrColumnOffsets_hu_h,
      int *csrRowIndeces_hu_hu, int *csrColumnOffsets_hu_hu,
      int *csrRowIndeces_hu_hv, int *csrColumnOffsets_hu_hv,
      int *csrRowIndeces_hu_heta, int *csrColumnOffsets_hu_heta,
      int *csrRowIndeces_hu_hw, int *csrColumnOffsets_hu_hw,
      // hv
      int *csrRowIndeces_hv_h, int *csrColumnOffsets_hv_h,
      int *csrRowIndeces_hv_hu, int *csrColumnOffsets_hv_hu,
      int *csrRowIndeces_hv_hv, int *csrColumnOffsets_hv_hv,
      int *csrRowIndeces_hv_heta, int *csrColumnOffsets_hv_heta,
      int *csrRowIndeces_hv_hw, int *csrColumnOffsets_hv_hw,
      // heta
      int *csrRowIndeces_heta_h, int *csrColumnOffsets_heta_h,
      int *csrRowIndeces_heta_hu, int *csrColumnOffsets_heta_hu,
      int *csrRowIndeces_heta_hv, int *csrColumnOffsets_heta_hv,
      int *csrRowIndeces_heta_heta, int *csrColumnOffsets_heta_heta,
      int *csrRowIndeces_heta_hw, int *csrColumnOffsets_heta_hw,
      // hw
      int *csrRowIndeces_hw_h, int *csrColumnOffsets_hw_h,
      int *csrRowIndeces_hw_hu, int *csrColumnOffsets_hw_hu,
      int *csrRowIndeces_hw_hv, int *csrColumnOffsets_hw_hv,
      int *csrRowIndeces_hw_heta, int *csrColumnOffsets_hw_heta,
      int *csrRowIndeces_hw_hw, int *csrColumnOffsets_hw_hw,
      //
      double *globalJacobian, int nExteriorElementBoundaries_global,
      int *exteriorElementBoundariesArray, int *elementBoundaryElementsArray,
      int *elementBoundaryLocalElementBoundariesArray, int *isDOFBoundary_h,
      int *isDOFBoundary_hu, int *isDOFBoundary_hv,
      int *isAdvectiveFluxBoundary_h, int *isAdvectiveFluxBoundary_hu,
      int *isAdvectiveFluxBoundary_hv, int *isDiffusiveFluxBoundary_hu,
      int *isDiffusiveFluxBoundary_hv, double *ebqe_bc_h_ext,
      double *ebqe_bc_flux_mass_ext, double *ebqe_bc_flux_mom_hu_adv_ext,
      double *ebqe_bc_flux_mom_hv_adv_ext, double *ebqe_bc_hu_ext,
      double *ebqe_bc_flux_hu_diff_ext, double *ebqe_penalty_ext,
      double *ebqe_bc_hv_ext, double *ebqe_bc_flux_hv_diff_ext,
      int *csrColumnOffsets_eb_h_h, int *csrColumnOffsets_eb_h_hu,
      int *csrColumnOffsets_eb_h_hv, int *csrColumnOffsets_eb_hu_h,
      int *csrColumnOffsets_eb_hu_hu, int *csrColumnOffsets_eb_hu_hv,
      int *csrColumnOffsets_eb_hv_h, int *csrColumnOffsets_eb_hv_hu,
      int *csrColumnOffsets_eb_hv_hv, double dt) = 0;
  virtual void calculateLumpedMassMatrix( // element
      double *mesh_trial_ref, double *mesh_grad_trial_ref, double *mesh_dof,
      double *mesh_velocity_dof, double MOVING_DOMAIN, int *mesh_l2g,
      double *dV_ref, double *h_trial_ref, double *h_grad_trial_ref,
      double *h_test_ref, double *h_grad_test_ref, double *vel_trial_ref,
      double *vel_grad_trial_ref, double *vel_test_ref,
      double *vel_grad_test_ref,
      // element boundary
      double *mesh_trial_trace_ref, double *mesh_grad_trial_trace_ref,
      double *dS_ref, double *h_trial_trace_ref, double *h_grad_trial_trace_ref,
      double *h_test_trace_ref, double *h_grad_test_trace_ref,
      double *vel_trial_trace_ref, double *vel_grad_trial_trace_ref,
      double *vel_test_trace_ref, double *vel_grad_test_trace_ref,
      double *normal_ref, double *boundaryJac_ref,
      // physics
      double *elementDiameter, int nElements_global, double useRBLES,
      double useMetrics, double alphaBDF, double nu, double g, int *h_l2g,
      int *vel_l2g, double *b_dof, double *h_dof, double *hu_dof,
      double *hv_dof, double *h_dof_sge, double *hu_dof_sge, double *hv_dof_sge,
      double *q_mass_acc_beta_bdf, double *q_mom_hu_acc_beta_bdf,
      double *q_mom_hv_acc_beta_bdf, double *q_cfl, int *sdInfo_hu_hu_rowptr,
      int *sdInfo_hu_hu_colind, int *sdInfo_hu_hv_rowptr,
      int *sdInfo_hu_hv_colind, int *sdInfo_hv_hv_rowptr,
      int *sdInfo_hv_hv_colind, int *sdInfo_hv_hu_rowptr,
      int *sdInfo_hv_hu_colind,
      // h
      int *csrRowIndeces_h_h, int *csrColumnOffsets_h_h,
      int *csrRowIndeces_h_hu, int *csrColumnOffsets_h_hu,
      int *csrRowIndeces_h_hv, int *csrColumnOffsets_h_hv,
      int *csrRowIndeces_h_heta, int *csrColumnOffsets_h_heta,
      int *csrRowIndeces_h_hw, int *csrColumnOffsets_h_hw,
      // hu
      int *csrRowIndeces_hu_h, int *csrColumnOffsets_hu_h,
      int *csrRowIndeces_hu_hu, int *csrColumnOffsets_hu_hu,
      int *csrRowIndeces_hu_hv, int *csrColumnOffsets_hu_hv,
      int *csrRowIndeces_hu_heta, int *csrColumnOffsets_hu_heta,
      int *csrRowIndeces_hu_hw, int *csrColumnOffsets_hu_hw,
      // hv
      int *csrRowIndeces_hv_h, int *csrColumnOffsets_hv_h,
      int *csrRowIndeces_hv_hu, int *csrColumnOffsets_hv_hu,
      int *csrRowIndeces_hv_hv, int *csrColumnOffsets_hv_hv,
      int *csrRowIndeces_hv_heta, int *csrColumnOffsets_hv_heta,
      int *csrRowIndeces_hv_hw, int *csrColumnOffsets_hv_hw,
      // heta
      int *csrRowIndeces_heta_h, int *csrColumnOffsets_heta_h,
      int *csrRowIndeces_heta_hu, int *csrColumnOffsets_heta_hu,
      int *csrRowIndeces_heta_hv, int *csrColumnOffsets_heta_hv,
      int *csrRowIndeces_heta_heta, int *csrColumnOffsets_heta_heta,
      int *csrRowIndeces_heta_hw, int *csrColumnOffsets_heta_hw,
      // hw
      int *csrRowIndeces_hw_h, int *csrColumnOffsets_hw_h,
      int *csrRowIndeces_hw_hu, int *csrColumnOffsets_hw_hu,
      int *csrRowIndeces_hw_hv, int *csrColumnOffsets_hw_hv,
      int *csrRowIndeces_hw_heta, int *csrColumnOffsets_hw_heta,
      int *csrRowIndeces_hw_hw, int *csrColumnOffsets_hw_hw,
      double *globalJacobian, int nExteriorElementBoundaries_global,
      int *exteriorElementBoundariesArray, int *elementBoundaryElementsArray,
      int *elementBoundaryLocalElementBoundariesArray, int *isDOFBoundary_h,
      int *isDOFBoundary_hu, int *isDOFBoundary_hv,
      int *isAdvectiveFluxBoundary_h, int *isAdvectiveFluxBoundary_hu,
      int *isAdvectiveFluxBoundary_hv, int *isDiffusiveFluxBoundary_hu,
      int *isDiffusiveFluxBoundary_hv, double *ebqe_bc_h_ext,
      double *ebqe_bc_flux_mass_ext, double *ebqe_bc_flux_mom_hu_adv_ext,
      double *ebqe_bc_flux_mom_hv_adv_ext, double *ebqe_bc_hu_ext,
      double *ebqe_bc_flux_hu_diff_ext, double *ebqe_penalty_ext,
      double *ebqe_bc_hv_ext, double *ebqe_bc_flux_hv_diff_ext,
      int *csrColumnOffsets_eb_h_h, int *csrColumnOffsets_eb_h_hu,
      int *csrColumnOffsets_eb_h_hv, int *csrColumnOffsets_eb_hu_h,
      int *csrColumnOffsets_eb_hu_hu, int *csrColumnOffsets_eb_hu_hv,
      int *csrColumnOffsets_eb_hv_h, int *csrColumnOffsets_eb_hv_hu,
      int *csrColumnOffsets_eb_hv_hv, double dt) = 0;
};

template <class CompKernelType, int nSpace, int nQuadraturePoints_element,
          int nDOF_mesh_trial_element, int nDOF_trial_element,
          int nDOF_test_element, int nQuadraturePoints_elementBoundary>
class GN_SW2DCV : public GN_SW2DCV_base {
public:
  const int nDOF_test_X_trial_element;
  CompKernelType ck;
  GN_SW2DCV()
      : nDOF_test_X_trial_element(nDOF_test_element * nDOF_trial_element),
        ck() {
    std::cout << "Constructing GN_SW2DCV<CompKernelTemplate<" << nSpace << ","
              << nQuadraturePoints_element << "," << nDOF_mesh_trial_element
              << "," << nDOF_trial_element << "," << nDOF_test_element << ","
              << nQuadraturePoints_elementBoundary << ">());" << std::endl
              << std::flush;
  }

  inline double maxWaveSpeedSharpInitialGuess(double g, double nx, double ny,
                                              double hL, double huL, double hvL,
                                              double hetaL, double lumpedL,
                                              double hR, double huR, double hvR,
                                              double hetaR, double lumpedR,
                                              double hEpsL, double hEpsR,
                                              bool debugging) {
    double lambda1, lambda3;

    double hVelL = nx * huL + ny * hvL;
    double hVelR = nx * huR + ny * hvR;
    double velL = 2 * hL / (hL * hL + std::pow(fmax(hL, hEpsL), 2)) * hVelL;
    double velR = 2 * hR / (hR * hR + std::pow(fmax(hR, hEpsR), 2)) * hVelR;
    double etaL = 2 * hL / (hL * hL + std::pow(fmax(hL, hEpsL), 2)) * hetaL;
    double etaR = 2 * hR / (hR * hR + std::pow(fmax(hR, hEpsR), 2)) * hetaR;
    double meshSizeL = sqrt(lumpedL);
    double meshSizeR = sqrt(lumpedR);

    /* See equation 4.12 from mGN paper:
      1-eigenvalue: uL-sqrt(g*hL)*sqrt(1 + augL)
      3-eigenvalue: uR+sqrt(g*hR)*sqrt(1 + augR)
    */
    lambda1 = GN_nu1(g, hL, velL, etaL, meshSizeL);
    lambda3 = GN_nu3(g, hR, velR, etaR, meshSizeR);
    return fmax(fabs(lambda1), fabs(lambda3));
  }

// Don't need this right now but keeping it just in case we need it later
#if 0
  inline double maxWaveSpeedIterativeProcess(double g, double nx, double ny,
                                             double hL, double huL, double
                                             hvL, double hR, double huR,
                                             double hvR, double hEpsL, double
                                             hEpsR, bool verbose) {
    double tol = 1E-15;
    // 1-eigenvalue: uL-sqrt(g*hL)
    // 3-eigenvalue: uR+sqrt(g*hR)

    double hVelL = nx * huL + ny * hvL;
    double hVelR = nx * huR + ny * hvR;
    double velL = 2 * hL / (hL * hL + std::pow(fmax(hL, hEpsL), 2)) * hVelL;
    double velR = 2 * hR / (hR * hR + std::pow(fmax(hR, hEpsR), 2)) * hVelR;

    double lambda1, lambda3;

    // CHECK IF BOTH STATES ARE DRY:
    if (hL == 0 && hR == 0) {
      lambda1 = 0.;
      lambda3 = 0.;
      return 0.;
    } else if (hL == 0) // left dry state
    {
      lambda1 = velR - 2 * sqrt(g * hR);
      lambda3 = velR + sqrt(g * hR);
      return fmax(fabs(lambda1), fabs(lambda3));
    } else if (hR == 0) // right dry state
    {
      lambda1 = velL - sqrt(g * hL);
      lambda3 = velL + 2 * sqrt(g * hL);
      return fmax(fabs(lambda1), fabs(lambda3));
    } else {
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
      double hMin = fmin(hL, hR);
      double hMax = fmin(hL, hR);
      double phiMin = phi(g, hMin, hL, hR, velL, velR);
      double phiMax = phi(g, hMax, hL, hR, velL, velR);
      if (0 <= phiMin) {
        // This is a 1- and 3-rarefactions situation. We know the solution in
        // this case
        lambda1 = velL - sqrt(g * hL);
        lambda3 = velR + sqrt(g * hR);

        std::cout << "lambda Min, lambda Max: " << lambda1 << ", " << lambda3
                  << std::endl;

        return fmax(fabs(lambda1), fabs(lambda3));
      }
      if (phiMax == 0) // if hMax "hits" hStar (very unlikely)
      {
        hStar = hMax;
        lambda1 = nu1(g, hStar, hL, velL);
        lambda3 = nu3(g, hStar, hR, velR);
        return fmax(fabs(lambda1), fabs(lambda3));
      }
      double hStarTwoRarefactions =
          std::pow(velL - velR + 2 * sqrt(g) * (sqrt(hL) + sqrt(hR)), 2) / 16
          / g;
      if (phiMax < 0) // This is a 1- and 3-shock situation
      {
        hStarL = hMax;
        hStarR = hStarTwoRarefactions;
      } else // Here we have one shock and one rarefaction
      {
        hStarL = hMin;
        hStarR = fmin(hMax, hStarTwoRarefactions);
      }

      // improve estimate from below via one newton step (not required)
      hStarL = fmax(hStarL, hStarR - phi(g, hStarR, hL, hR, velL, velR) /
                                         phip(g, hStarR, hL, hR));
      // COMPUTE lambdaMin0 and lambdaMax0
      double nu11 = nu1(g, hStarR, hL, velL);
      double nu12 = nu1(g, hStarL, hL, velL);
      double nu31 = nu3(g, hStarL, hR, velR);
      double nu32 = nu3(g, hStarR, hR, velR);

      double lambdaMin = fmax(fmax(nu31, 0), fmax(-nu12, 0));
      double lambdaMax = fmax(fmax(nu32, 0), fmax(-nu11, 0));

      if (verbose) {
        std::cout << "hStarL, hStarR: " << hStarL << ", " << hStarR << "\t"
                  << "lambda Min, lambda Max: " << lambdaMin << ", "
                  << lambdaMax << std::endl;
      }
      // CHECK IF TOL IS SATISFIED. O.W. GOT TO ALGORITHM 2 //
      if (lambdaMin > 0 && lambdaMax / lambdaMin - 1 <= tol)
        return lambdaMax;
      else // Proceed to algorithm 2
      {
        ///////////////////////////////////////////
        // ALGORITHM 2: ESTIMATION OF LAMBDA MAX //
        ///////////////////////////////////////////
        // Requires: hStarL, hStarR
        // Ensures: lambdaMax
        int aux_counter = 0;
        while (true) {
          aux_counter++;
          // Start having lambdaMin and lambdaMax
          // Check if current lambdaMin and lambdaMax satisfy the tolerance
          if (verbose) {
            std::cout << lambdaMin << ", " << lambdaMax << std::endl;
          }
          if (lambdaMin > 0 && lambdaMax / lambdaMin - 1 <= tol)
            return lambdaMax;
          // Check for round off error
          if (phi(g, hStarL, hL, hR, velL, velR) > 0 ||
              phi(g, hStarR, hL, hR, velL, velR) < 0)
            return lambdaMax;

          // Save old estimates of hStar
          double hStarL_old = hStarL;
          double hStarR_old = hStarR;
          // Compute new estimates on hStarL and hStarR
          // NOTE (MQL): hStarL and hStarR must be computed using the old
          values hStarL = hStarLFromQuadPhiFromAbove(g, hStarL_old,
          hStarR_old, hL, hR,
                                              velL, velR);
          hStarR = hStarRFromQuadPhiFromBelow(g, hStarL_old, hStarR_old, hL,
          hR,
                                              velL, velR);

          // Compute lambdaMax and lambdaMin
          nu11 = nu1(g, hStarR, hL, velL);
          nu12 = nu1(g, hStarL, hL, velL);
          nu31 = nu3(g, hStarL, hR, velR);
          nu32 = nu3(g, hStarR, hR, velR);

          lambdaMin = fmax(fmax(nu31, 0), fmax(-nu12, 0));
          lambdaMax = fmax(fmax(nu32, 0), fmax(-nu11, 0));

          if (aux_counter > 1000) // TMP
          {
            std::cout << "**** AUX COUNTER > 1000... aborting!" << std::endl;
            std::cout << "**** Initial guess hStar: " << hStar0 << std::endl;

            hStar = hStar0;
            lambda1 = nu1(g, hStar, hL, velL);
            lambda3 = nu3(g, hStar, hR, velR);
            std::cout << "**** Initial estimate of max wave speed: "
                      << fmax(fabs(lambda1), fabs(lambda3)) << std::endl;

            abort();
          }
          // else
          //{
          //  std::cout << "*****... AUX COUNTER: " << aux_counter <<
          std::endl;
          //  //TMP
          //}
        }
      }
    }
  }
#endif

  /* I don't think this does anything of importance. -EJT */
  inline void calculateCFL(const double &elementDiameter, const double &g,
                           const double &h, const double &hu, const double &hv,
                           const double hEps, double &cfl) {
    double cflx, cfly, c = sqrt(fmax(g * hEps, g * h));
    double u = 2 * h / (h * h + std::pow(fmax(h, hEps), 2)) * hu;
    double v = 2 * h / (h * h + std::pow(fmax(h, hEps), 2)) * hv;

    if (u > 0.0)
      cflx = (u + c) / elementDiameter;
    else
      cflx = fabs(u - c) / elementDiameter;

    if (v > 0.0)
      cfly = (v + c) / elementDiameter;
    else
      cfly = fabs(v - c) / elementDiameter;
    cfl = sqrt(cflx * cflx + cfly * cfly); // hack, conservative estimate
  }

  void FCTStep(double dt,
               int NNZ,     // number on non-zero entries on sparsity pattern
               int numDOFs, // number of DOFs
               double *lumped_mass_matrix, // lumped mass matrix (as vector))
               double *h_old,              // DOFs of solution at last stage
               double *hu_old, double *hv_old, double *heta_old, double *hw_old,
               double *b_dof,
               double *high_order_hnp1, // DOFs of high order solution at tnp1
               double *high_order_hunp1, double *high_order_hvnp1,
               double *high_order_hetanp1, double *high_order_hwnp1,
               double *extendedSourceTerm_hu, double *extendedSourceTerm_hv,
               double *extendedSourceTerm_heta, double *extendedSourceTerm_hw,
               double *limited_hnp1, double *limited_hunp1,
               double *limited_hvnp1, double *limited_hetanp1,
               double *limited_hwnp1,
               int *csrRowIndeces_DofLoops,    // csr row indeces
               int *csrColumnOffsets_DofLoops, // csr column offsets
               double *MassMatrix,             // mass matrix
               double *dH_minus_dL, double *muH_minus_muL, double hEps,
               double *hReg, int LUMPED_MASS_MATRIX, double *dLow, double *hBT,
               double *huBT, double *hvBT, double *hetaBT, double *hwBT) {
    Rneg.resize(numDOFs, 0.0);
    Rpos.resize(numDOFs, 0.0);
    Rneg_heta.resize(numDOFs, 0.0);
    Rpos_heta.resize(numDOFs, 0.0);
    hLow.resize(numDOFs, 0.0);
    huLow.resize(numDOFs, 0.0);
    hvLow.resize(numDOFs, 0.0);
    hetaLow.resize(numDOFs, 0.0);
    hwLow.resize(numDOFs, 0.0);
    std::cout << "FCT TEST?" << std::endl;

    /* For reference, see section 6.3 of friction shallow water paper
    by Guermond, Quezada De Luna, Popov, Kees, Farthing */

    ////////////////////////
    // FIRST LOOP in DOFs //
    ////////////////////////
    int ij = 0;
    for (int i = 0; i < numDOFs; i++) {
      // read some vectors
      double high_order_hnp1i = high_order_hnp1[i];
      double hi = h_old[i];
      double hui = hu_old[i];
      double hvi = hv_old[i];
      double high_order_hetanp1i = high_order_hetanp1[i];
      double hetai = heta_old[i];
      double hwi = hw_old[i];
      double Zi = b_dof[i];
      double mi = lumped_mass_matrix[i];

      double hiMin = hi;
      double hiMax = hi;
      double Pnegi = 0., Pposi = 0.;

      /* Actually it's probably better to put this all
      in a function with inputs to be able to limit any of the conserved
      variables - EJT */

      // For limiting heta
      double hetai_Min = hetai;
      double hetai_Max = hetai;
      double Pnegi_heta = 0., Pposi_heta = 0.;

      // LOW ORDER SOLUTION without extended source term. Eqn 6.23
      hLow[i] = hi;
      huLow[i] = hui;
      hvLow[i] = hvi;
      hetaLow[i] = hetai;
      hwLow[i] = hwi;

      // LOOP OVER THE SPARSITY PATTERN (j-LOOP)//
      for (int offset = csrRowIndeces_DofLoops[i];
           offset < csrRowIndeces_DofLoops[i + 1]; offset++) {
        int j = csrColumnOffsets_DofLoops[offset];

        // read some vectors
        double hj = h_old[j];
        double hetaj = heta_old[j];
        double Zj = b_dof[j];
        // COMPUTE STAR SOLUTION // hStar, huStar and hvStar
        double hStarij = fmax(0., hi + Zi - fmax(Zi, Zj));
        double hStarji = fmax(0., hj + Zj - fmax(Zi, Zj));

        // i-th row of flux correction matrix
        double ML_minus_MC = (LUMPED_MASS_MATRIX == 1
                                  ? 0.
                                  : (i == j ? 1. : 0.) * mi - MassMatrix[ij]);
        double FluxCorrectionMatrix1 =
            ML_minus_MC * (high_order_hnp1[j] - hj - (high_order_hnp1i - hi)) +
            dt * (dH_minus_dL[ij] - muH_minus_muL[ij]) * (hStarji - hStarij) +
            dt * muH_minus_muL[ij] * (hj - hi);

        double FluxCorrectionMatrix_heta =
            ML_minus_MC *
                (high_order_hetanp1[j] - hetaj - (high_order_hetanp1i - hi)) +
            dt * (dH_minus_dL[ij] - muH_minus_muL[ij]) * (hStarji - hStarij) +
            dt * muH_minus_muL[ij] * (hetaj - hetai);

        // COMPUTE P VECTORS //
        Pnegi +=
            FluxCorrectionMatrix1 * ((FluxCorrectionMatrix1 < 0) ? 1. : 0.);
        Pposi +=
            FluxCorrectionMatrix1 * ((FluxCorrectionMatrix1 > 0) ? 1. : 0.);

        Pnegi_heta += FluxCorrectionMatrix_heta *
                      ((FluxCorrectionMatrix_heta < 0) ? 1. : 0.);
        Pposi_heta += FluxCorrectionMatrix_heta *
                      ((FluxCorrectionMatrix_heta > 0) ? 1. : 0.);

        // COMPUTE LOCAL BOUNDS //
        hiMin = std::min(hiMin, hBT[ij]);
        hiMax = std::max(hiMax, hBT[ij]);
        hetai_Min = std::min(hetai_Min, hetaBT[ij]);
        hetai_Max = std::max(hetai_Max, hetaBT[ij]);

        // COMPUTE LOW ORDER SOLUTION (WITHOUT EXTENDED SOURCE) //
        // SEE EQUATION 6.23 FROM SW FRICTION PAPER -EJT //
        if (i != j) {
          hLow[i] += hi * (-dt / mi * 2 * dLow[ij]) +
                     dt / mi * (2 * dLow[ij] * hBT[ij]);
          huLow[i] += hui * (-dt / mi * 2 * dLow[ij]) +
                      dt / mi * (2 * dLow[ij] * huBT[ij]);
          hvLow[i] += hvi * (-dt / mi * 2 * dLow[ij]) +
                      dt / mi * (2 * dLow[ij] * hvBT[ij]);
          hetaLow[i] += hetai * (-dt / mi * 2 * dLow[ij]) +
                        dt / mi * (2 * dLow[ij] * hetaBT[ij]);
          hwLow[i] += hwi * (-dt / mi * 2 * dLow[ij]) +
                      dt / mi * (2 * dLow[ij] * hwBT[ij]);
        }
        // update ij
        ij += 1;
      }
      // clean up hLow from round off error
      if (hLow[i] < hEps)
        hLow[i] = 0;
      ///////////////////////
      // COMPUTE Q VECTORS //
      ///////////////////////
      if (GLOBAL_FCT == 1)
        hiMin = 0;

      double Qnegi = mi * (hiMin - hLow[i]);
      double Qposi = mi * (hiMax - hLow[i]);

      double Qnegi_heta =
          mi * (hetai_Min - (hetaLow[i] - extendedSourceTerm_heta[i]));
      double Qposi_heta =
          mi * (hetai_Max - (hetaLow[i] - extendedSourceTerm_heta[i]));

      ///////////////////////
      // COMPUTE R VECTORS //
      ///////////////////////
      if (high_order_hnp1[i] <= hReg[i]) // hEps
      {
        Rneg[i] = 0.;
        Rpos[i] = 0.;
        // for heta - EJT
        Rneg_heta[i] = 0.;
        Rpos_heta[i] = 0.;
      } else {
        Rneg[i] = ((Pnegi == 0) ? 1. : std::min(1.0, Qnegi / Pnegi));
        Rpos[i] = ((Pposi == 0) ? 1. : std::min(1.0, Qposi / Pposi));
        // for heta
        Rneg_heta[i] =
            ((Pnegi_heta == 0) ? 1. : std::min(1.0, Qnegi_heta / Pnegi_heta));
        Rpos_heta[i] =
            ((Pposi_heta == 0) ? 1. : std::min(1.0, Qposi_heta / Pposi_heta));
      }
    } // i DOFs

    //////////////////////
    // COMPUTE LIMITERS //
    //////////////////////
    ij = 0;
    for (int i = 0; i < numDOFs; i++) {
      /* read some vectors */
      // high order
      double high_order_hnp1i = high_order_hnp1[i];
      double high_order_hunp1i = high_order_hunp1[i];
      double high_order_hvnp1i = high_order_hvnp1[i];
      double high_order_hetanp1i = high_order_hetanp1[i];
      double high_order_hwnp1i = high_order_hwnp1[i];
      // old
      double hi = h_old[i];
      double huni = hu_old[i];
      double hvni = hv_old[i];
      double hetani = heta_old[i];
      double hwni = hw_old[i];
      double Zi = b_dof[i];
      double mi = lumped_mass_matrix[i];
      double one_over_hiReg =
          2 * hi / (hi * hi + std::pow(fmax(hi, hEps), 2)); // hEps

      double ith_Limiter_times_FluxCorrectionMatrix1 = 0.;
      double ith_Limiter_times_FluxCorrectionMatrix2 = 0.;
      double ith_Limiter_times_FluxCorrectionMatrix3 = 0.;
      double ith_Limiter_times_FluxCorrectionMatrix4 = 0.;
      double ith_Limiter_times_FluxCorrectionMatrix5 = 0.;

      // LOOP OVER THE SPARSITY PATTERN (j-LOOP)//
      for (int offset = csrRowIndeces_DofLoops[i];
           offset < csrRowIndeces_DofLoops[i + 1]; offset++) {
        int j = csrColumnOffsets_DofLoops[offset];
        // read some vectors
        double hj = h_old[j];
        double hunj = hu_old[j];
        double hvnj = hv_old[j];
        double hetanj = heta_old[j];
        double hwnj = hw_old[j];
        double Zj = b_dof[j];
        double one_over_hjReg =
            2 * hj / (hj * hj + std::pow(fmax(hj, hEps), 2)); // hEps

        /*  COMPUTE STAR SOLUTIONS:
                hStar, huStar,  hvStar, hetaStar, and hwStar
         */
        double hStarij = fmax(0., hi + Zi - fmax(Zi, Zj));
        double huStarij = huni * hStarij * one_over_hiReg;
        double hvStarij = hvni * hStarij * one_over_hiReg;
        double hetaStarij = hetani * std::pow(hStarij * one_over_hiReg, 2);
        double hwStarij = hwni * hStarij * one_over_hiReg;

        double hStarji = fmax(0., hj + Zj - fmax(Zi, Zj));
        double huStarji = hunj * hStarji * one_over_hjReg;
        double hvStarji = hvnj * hStarji * one_over_hjReg;
        double hetaStarji = hetanj * std::pow(hStarji * one_over_hjReg, 2);
        double hwStarji = hwnj * hStarji * one_over_hjReg;

        // COMPUTE FLUX CORRECTION MATRICES
        double ML_minus_MC = (LUMPED_MASS_MATRIX == 1
                                  ? 0.
                                  : (i == j ? 1. : 0.) * mi - MassMatrix[ij]);
        double FluxCorrectionMatrix1 =
            ML_minus_MC * (high_order_hnp1[j] - hj - (high_order_hnp1i - hi)) +
            dt * (dH_minus_dL[ij] - muH_minus_muL[ij]) * (hStarji - hStarij) +
            dt * muH_minus_muL[ij] * (hj - hi);

        double FluxCorrectionMatrix2 =
            ML_minus_MC *
                (high_order_hunp1[j] - hunj - (high_order_hunp1i - huni)) +
            dt * (dH_minus_dL[ij] - muH_minus_muL[ij]) * (huStarji - huStarij) +
            dt * muH_minus_muL[ij] * (hunj - huni);

        double FluxCorrectionMatrix3 =
            ML_minus_MC *
                (high_order_hvnp1[j] - hvnj - (high_order_hvnp1i - hvni)) +
            dt * (dH_minus_dL[ij] - muH_minus_muL[ij]) * (hvStarji - hvStarij) +
            dt * muH_minus_muL[ij] * (hvnj - hvni);

        double FluxCorrectionMatrix4 =
            ML_minus_MC * (high_order_hetanp1[j] - hetanj -
                           (high_order_hetanp1i - hetani)) +
            dt * (dH_minus_dL[ij] - muH_minus_muL[ij]) *
                (hetaStarji - hetaStarij) +
            dt * muH_minus_muL[ij] * (hetanj - hetani);

        double FluxCorrectionMatrix5 =
            ML_minus_MC *
                (high_order_hwnp1[j] - hwnj - (high_order_hwnp1i - hwni)) +
            dt * (dH_minus_dL[ij] - muH_minus_muL[ij]) * (hwStarji - hwStarij) +
            dt * muH_minus_muL[ij] * (hwnj - hwni);

        // COMPUTE LIMITER based on water height
        double Lij = 0.;
        if (FluxCorrectionMatrix1 >= 0)
          Lij = std::min(Rpos[i], Rneg[j]);
        else
          Lij = std::min(Rneg[i], Rpos[j]);

        // COMPUTE LIMITER based on heta -EJT
        // Note that we set lij = min(lij_h,lij_heta)
        if (FluxCorrectionMatrix4 >= 0)
          Lij = std::min(Lij, std::min(Rpos_heta[i], Rneg_heta[j]));
        else
          Lij = std::min(Lij, std::min(Rneg_heta[i], Rpos_heta[j]));

        // take off all instances of GLOBAL_FCT -EJT
        if (GLOBAL_FCT == 1)
          Lij = (FluxCorrectionMatrix1 >= 0. ? std::min(1., Rneg[j])
                                             : std::min(Rneg[i], 1.));

        // COMPUTE LIMITED FLUX //
        ith_Limiter_times_FluxCorrectionMatrix1 += Lij * FluxCorrectionMatrix1;
        ith_Limiter_times_FluxCorrectionMatrix2 += Lij * FluxCorrectionMatrix2;
        ith_Limiter_times_FluxCorrectionMatrix3 += Lij * FluxCorrectionMatrix3;
        ith_Limiter_times_FluxCorrectionMatrix4 += Lij * FluxCorrectionMatrix4;
        ith_Limiter_times_FluxCorrectionMatrix5 += Lij * FluxCorrectionMatrix5;

        // update ij
        ij += 1;
      }

      double one_over_mi = 1.0 / lumped_mass_matrix[i];
      limited_hnp1[i] =
          hLow[i] + one_over_mi * ith_Limiter_times_FluxCorrectionMatrix1;
      limited_hunp1[i] =
          ((huLow[i] -
            dt / mi * extendedSourceTerm_hu[i]) // low_order_hunp1+...
           + one_over_mi * ith_Limiter_times_FluxCorrectionMatrix2);
      limited_hvnp1[i] =
          ((hvLow[i] -
            dt / mi * extendedSourceTerm_hv[i]) // low_order_hvnp1+...
           + one_over_mi * ith_Limiter_times_FluxCorrectionMatrix3);
      limited_hetanp1[i] =
          ((hetaLow[i] -
            dt / mi * extendedSourceTerm_heta[i]) // low_order_hetanp1+...
           + one_over_mi * ith_Limiter_times_FluxCorrectionMatrix4);
      limited_hwnp1[i] =
          ((hwLow[i] -
            dt / mi * extendedSourceTerm_hw[i]) // low_order_hwnp1+...
           + one_over_mi * ith_Limiter_times_FluxCorrectionMatrix5);

      if (limited_hnp1[i] < -hEps && dt < 1.0) {
        std::cout << "Limited water height is negative: "
                  << "hLow: " << hLow[i] << "\t"
                  << "hHigh: " << limited_hnp1[i] << "\t"
                  << " ... aborting!" << std::endl;
        abort();
      } else {
        // clean up uHigh from round off error
        if (limited_hnp1[i] < hEps)
          limited_hnp1[i] = 0;
        // double aux = fmax(limited_hnp1[i],hEps); // hEps
        double aux =
            fmax(limited_hnp1[i], hReg[i]); // hReg makes the code more robust
        limited_hunp1[i] *= 2 * std::pow(limited_hnp1[i], VEL_FIX_POWER) /
                            (std::pow(limited_hnp1[i], VEL_FIX_POWER) +
                             std::pow(aux, VEL_FIX_POWER));
        limited_hvnp1[i] *= 2 * std::pow(limited_hnp1[i], VEL_FIX_POWER) /
                            (std::pow(limited_hnp1[i], VEL_FIX_POWER) +
                             std::pow(aux, VEL_FIX_POWER));
        limited_hetanp1[i] *= 2 * std::pow(limited_hnp1[i], VEL_FIX_POWER) /
                              (std::pow(limited_hnp1[i], VEL_FIX_POWER) +
                               std::pow(aux, VEL_FIX_POWER));
        limited_hwnp1[i] *= 2 * std::pow(limited_hnp1[i], VEL_FIX_POWER) /
                            (std::pow(limited_hnp1[i], VEL_FIX_POWER) +
                             std::pow(aux, VEL_FIX_POWER));
      }
    }
  }

  void convexLimiting(
      double dt,
      int NNZ,     // number on non-zero entries on sparsity pattern
      int numDOFs, // number of DOFs
      double *lumped_mass_matrix, // lumped mass matrix (as vector))
      double *h_old,              // DOFs of solution at last stage
      double *hu_old, double *hv_old, double *heta_old, double *hw_old,
      double *b_dof,
      double *high_order_hnp1, // DOFs of high order solution at tnp1
      double *high_order_hunp1, double *high_order_hvnp1,
      double *high_order_hetanp1, double *high_order_hwnp1,
      double *extendedSourceTerm_hu, double *extendedSourceTerm_hv,
      double *extendedSourceTerm_heta, double *extendedSourceTerm_hw,
      double *limited_hnp1, double *limited_hunp1, double *limited_hvnp1,
      double *limited_hetanp1, double *limited_hwnp1,
      int *csrRowIndeces_DofLoops,    // csr row indeces
      int *csrColumnOffsets_DofLoops, // csr column offsets
      double *MassMatrix,             // mass matrix
      double *dH_minus_dL, double *muH_minus_muL, double hEps, double *hReg,
      int LUMPED_MASS_MATRIX, double *dLow, double *hBT, double *huBT,
      double *hvBT, double *hetaBT, double *hwBT, double *new_SourceTerm_hu,
      double *new_SourceTerm_hv, double *new_SourceTerm_heta,
      double *new_SourceTerm_hw) {

    Rneg.resize(numDOFs, 0.0);
    Rpos.resize(numDOFs, 0.0);
    Rneg_heta.resize(numDOFs, 0.0);
    Rpos_heta.resize(numDOFs, 0.0);
    hLow.resize(numDOFs, 0.0);
    huLow.resize(numDOFs, 0.0);
    hvLow.resize(numDOFs, 0.0);
    hetaLow.resize(numDOFs, 0.0);
    hwLow.resize(numDOFs, 0.0);
    Kmax.resize(numDOFs, 0.0);

#if 0
    //////////////////////////////
    // For relaxation of bounds //
    //////////////////////////////

    // Change all of these to valarray

    std::vector<double> urelax;
    urelax.resize(numDOFs, 0.0);
    std::vector<double> drelax;
    drelax.resize(numDOFs, 0.0);

    std::vector<double> delta_Sqd_h;
    delta_Sqd_h.resize(numDOFs, 0.0);
    std::vector<double> bar_deltaSqd_h;
    bar_deltaSqd_h.resize(numDOFs, 0.0);
    // for heta
    std::vector<double> delta_Sqd_heta;
    delta_Sqd_heta.resize(numDOFs, 0.0);
    std::vector<double> bar_deltaSqd_heta;
    bar_deltaSqd_heta.resize(numDOFs, 0.0);

    double size_of_domain = 0.0;
    for (int i = 0; i < numDOFs; i++) {
      size_of_domain += lumped_mass_matrix[i];
    }

    // First loop
    int ij = 0;
    for (int i = 0; i < numDOFs; i++) {
      urelax[i] =
          1.0 +
          2.0 * std::pow(sqrt(sqrt(lumped_mass_matrix[i] / size_of_domain)), 3);
      drelax[i] =
          1.0 -
          2.0 * std::pow(sqrt(sqrt(lumped_mass_matrix[i] / size_of_domain)), 3);
      for (int offset = csrRowIndeces_DofLoops[i];
           offset < csrRowIndeces_DofLoops[i + 1]; offset++) {
        int j = csrColumnOffsets_DofLoops[offset];
        if (i != j) {
          delta_Sqd_h[i] += h_old[i] - h_old[j];
          delta_Sqd_heta[i] += heta_old[i] - heta_old[j];
        }
        // UPDATE ij //
        ij += 1;
      } // j loop ends here
    }   // i loops ends here

    // Second loop
    ij = 0;
    for (int i = 0; i < numDOFs; i++) {
      for (int offset = csrRowIndeces_DofLoops[i];
           offset < csrRowIndeces_DofLoops[i + 1]; offset++) {
        int j = csrColumnOffsets_DofLoops[offset];
        if (i != j) {
          bar_deltaSqd_h[i] += delta_Sqd_h[j] + delta_Sqd_h[i];
          bar_deltaSqd_heta[i] += delta_Sqd_heta[j] + delta_Sqd_heta[i];
        }
        // UPDATE ij //
        ij += 1;
      } // j loop ends here
    }   // i loops ends here
    for (int i = 0; i < numDOFs; i++) {
      bar_deltaSqd_h[i] =
          bar_deltaSqd_h[i] /
          (csrRowIndeces_DofLoops[i + 1] - csrRowIndeces_DofLoops[i]) / 2.0;
      bar_deltaSqd_heta[i] =
          bar_deltaSqd_heta[i] /
          (csrRowIndeces_DofLoops[i + 1] - csrRowIndeces_DofLoops[i]) / 2.0;
    }
#endif

    ////////////////////////
    // FIRST LOOP in DOFs //
    ////////////////////////
    int ij = 0;
    for (int i = 0; i < numDOFs; i++) {
      // read some vectors
      double high_order_hnp1i = high_order_hnp1[i];
      double hi = h_old[i];
      double hui = hu_old[i];
      double hvi = hv_old[i];
      double high_order_hetanp1i = high_order_hetanp1[i];
      double hetai = heta_old[i];
      double hwi = hw_old[i];
      double Zi = b_dof[i];
      double mi = lumped_mass_matrix[i];
      double one_over_hiReg =
          2 * hi / (hi * hi + std::pow(fmax(hi, hEps), 2)); // hEps

      // For limiting h, initializing hMin and hMax
      double hiMin = hi;
      double hiMax = hi;
      double Pnegi = 0., Pposi = 0.;

      // For limiting heta, initializing heta_Min, heta_max
      double hetai_Min = hetai;
      double hetai_Max = hetai;
      double Pnegi_heta = 0., Pposi_heta = 0.;

      /* LOW ORDER SOLUTION */
      hLow[i] = hi;
      huLow[i] = hui;
      hvLow[i] = hvi;
      hetaLow[i] = hetai;
      hwLow[i] = hwi;
      Kmax[i] = 0;

      // LOOP OVER THE SPARSITY PATTERN (j-LOOP)//
      for (int offset = csrRowIndeces_DofLoops[i];
           offset < csrRowIndeces_DofLoops[i + 1]; offset++) {
        int j = csrColumnOffsets_DofLoops[offset];
        double psi_ij = 0;
        double one_over_hBT =
            2.0 * hBT[ij] /
            (hBT[ij] * hBT[ij] + std::pow(fmax(hBT[ij], hEps), 2));
        psi_ij = one_over_hBT * (huBT[ij] * huBT[ij] + huBT[ij] * huBT[ij]) /
                 2.0; // Eqn (6.31)
        Kmax[i] = fmax(psi_ij, Kmax[i]);

        // read some vectors
        double hj = h_old[j];
        double one_over_hjReg =
            2 * hj / (hj * hj + std::pow(fmax(hj, hEps), 2)); // hEps
        double hetaj = heta_old[j];
        double Zj = b_dof[j];

        // COMPUTE STAR SOLUTION // hStar, huStar and hvStar
        double hStarij = fmax(0., hi + Zi - fmax(Zi, Zj));
        double hStarji = fmax(0., hj + Zj - fmax(Zi, Zj));
        double hetaStarji = hetaj * std::pow(hStarji * one_over_hjReg, 2);
        double hetaStarij = hetai * std::pow(hStarij * one_over_hiReg, 2);

        // i-th row of flux correction matrix
        double ML_minus_MC = (LUMPED_MASS_MATRIX == 1
                                  ? 0.
                                  : (i == j ? 1. : 0.) * mi - MassMatrix[ij]);
        double FluxCorrectionMatrix1 =
            ML_minus_MC * (high_order_hnp1[j] - hj - (high_order_hnp1i - hi)) +
            dt * (dH_minus_dL[ij] - muH_minus_muL[ij]) * (hStarji - hStarij) +
            dt * muH_minus_muL[ij] * (hj - hi);

        double FluxCorrectionMatrix_heta =
            ML_minus_MC * (high_order_hetanp1[j] - hetaj -
                           (high_order_hetanp1i - hetai)) +
            dt * (dH_minus_dL[ij] - muH_minus_muL[ij]) *
                (hetaStarji - hetaStarij) +
            dt * muH_minus_muL[ij] * (hetaj - hetai);

        // COMPUTE P VECTORS //
        Pnegi +=
            FluxCorrectionMatrix1 * ((FluxCorrectionMatrix1 < 0) ? 1. : 0.);
        Pposi +=
            FluxCorrectionMatrix1 * ((FluxCorrectionMatrix1 > 0) ? 1. : 0.);

        Pnegi_heta += FluxCorrectionMatrix_heta *
                      ((FluxCorrectionMatrix_heta < 0) ? 1. : 0.);
        Pposi_heta += FluxCorrectionMatrix_heta *
                      ((FluxCorrectionMatrix_heta > 0) ? 1. : 0.);

        // COMPUTE LOCAL BOUNDS //
        hiMin = std::min(hiMin, hBT[ij]);
        hiMax = std::max(hiMax, hBT[ij]);

        hetai_Min = std::min(hetai_Min, hetaBT[ij]);
        hetai_Max = std::max(hetai_Max, hetaBT[ij]);

// Then do relaxation of bounds here. If confused, see (4.12) of Euler
// convex limiting paper. Turn this off for now
#if 0
        hiMin = std::max(drelax[i] * hiMin,
                         hiMin - std::abs(bar_deltaSqd_h[i]) / 2.0);
        hetai_Min = std::max(drelax[i] * hetai_Min,
                             hetai_Min - std::abs(bar_deltaSqd_heta[i])
                             / 2.0);

        hiMax = std::min(urelax[i] * hiMax,
                         hiMax + std::abs(bar_deltaSqd_h[i]) / 2.0);
        hetai_Max = std::min(urelax[i] * hetai_Max,
                             hetai_Max + std::abs(bar_deltaSqd_heta[i])
                             / 2.0);
#endif

        /* COMPUTE LOW ORDER SOLUTION. See EQN 6.23 in
         * SW friction paper */
        if (i != j) {
          hLow[i] += hi * (-dt / mi * 2 * dLow[ij]) +
                     dt / mi * (2 * dLow[ij] * hBT[ij]);
          huLow[i] += hui * (-dt / mi * 2 * dLow[ij]) +
                      dt / mi * (2 * dLow[ij] * huBT[ij]);
          hvLow[i] += hvi * (-dt / mi * 2 * dLow[ij]) +
                      dt / mi * (2 * dLow[ij] * hvBT[ij]);
          hetaLow[i] += hetai * (-dt / mi * 2 * dLow[ij]) +
                        dt / mi * (2 * dLow[ij] * hetaBT[ij]);
          hwLow[i] += hwi * (-dt / mi * 2 * dLow[ij]) +
                      dt / mi * (2 * dLow[ij] * hwBT[ij]);
        }
        // UPDATE ij //
        ij += 1;
      } // j loop ends here

      // clean up hLow from round off error
      if (hLow[i] < hEps)
        hLow[i] = 0.0;

      ///////////////////////
      // COMPUTE Q VECTORS //
      ///////////////////////
      double Qnegi = std::min(mi * (hiMin - hLow[i]), 0.0);
      double Qposi = std::max(mi * (hiMax - hLow[i]), 0.0);
      // for heta, note that we subtract the sources term from the low order
      // solution
      double Qnegi_heta = std::min(
          mi * (hetai_Min - (hetaLow[i] - new_SourceTerm_heta[i] * dt / mi)),
          0.0);
      double Qposi_heta = std::max(
          mi * (hetai_Max - (hetaLow[i] - new_SourceTerm_heta[i] * dt / mi)),
          0.0);

      ///////////////////////
      // COMPUTE R VECTORS //
      ///////////////////////
      if (high_order_hnp1[i] <= hEps) // hEps
      {
        Rneg[i] = 0.;
        Rpos[i] = 0.;
        Rneg_heta[i] = 0.;
        Rpos_heta[i] = 0.;
      } else {
        Rneg[i] = ((Pnegi == 0) ? 1. : std::min(1.0, Qnegi / Pnegi));
        Rpos[i] = ((Pposi == 0) ? 1. : std::min(1.0, Qposi / Pposi));
        Rneg_heta[i] =
            ((Pnegi_heta == 0) ? 1. : std::min(1.0, Qnegi_heta / Pnegi_heta));
        Rpos_heta[i] =
            ((Pposi_heta == 0) ? 1. : std::min(1.0, Qposi_heta / Pposi_heta));
      }
    } // i DOFs

    //////////////////////
    // COMPUTE LIMITERS //
    //////////////////////
    ij = 0;
    for (int i = 0; i < numDOFs; i++) {
      // read some vectors
      double high_order_hnp1i = high_order_hnp1[i];
      double high_order_hunp1i = high_order_hunp1[i];
      double high_order_hvnp1i = high_order_hvnp1[i];
      double high_order_hetanp1i = high_order_hetanp1[i];
      double high_order_hwnp1i = high_order_hwnp1[i];
      double hi = h_old[i];
      double huni = hu_old[i];
      double hvni = hv_old[i];
      double hetani = heta_old[i];
      double hwni = hw_old[i];
      double Zi = b_dof[i];
      double mi = lumped_mass_matrix[i];
      double one_over_hiReg =
          2 * hi / (hi * hi + std::pow(fmax(hi, hEps), 2)); // hEps

      double ith_Limiter_times_FluxCorrectionMatrix1 = 0.;
      double ith_Limiter_times_FluxCorrectionMatrix2 = 0.;
      double ith_Limiter_times_FluxCorrectionMatrix3 = 0.;
      double ith_Limiter_times_FluxCorrectionMatrix4 = 0.;
      double ith_Limiter_times_FluxCorrectionMatrix5 = 0.;

      double huLow_minus_source_i = (huLow[i] - new_SourceTerm_hu[i] * dt / mi);
      double hvLow_minus_source_i = (hvLow[i] - new_SourceTerm_hv[i] * dt / mi);

      double ci =
          Kmax[i] * hLow[i] -
          0.5 * (huLow_minus_source_i * huLow_minus_source_i +
                 hvLow_minus_source_i * hvLow_minus_source_i); // for KE lim.

      // LOOP OVER THE SPARSITY PATTERN (j-LOOP)//
      for (int offset = csrRowIndeces_DofLoops[i];
           offset < csrRowIndeces_DofLoops[i + 1]; offset++) {
        int j = csrColumnOffsets_DofLoops[offset];
        // read some vectors
        double hj = h_old[j];
        double hunj = hu_old[j];
        double hvnj = hv_old[j];
        double hetanj = heta_old[j];
        double hwnj = hw_old[j];
        double Zj = b_dof[j];
        double one_over_hjReg =
            2 * hj / (hj * hj + std::pow(fmax(hj, hEps), 2)); // hEps

        // COMPUTE STAR SOLUTION // hStar, huStar, hvStar, hetaStar, and
        // hwStar
        double hStarij = fmax(0., hi + Zi - fmax(Zi, Zj));
        double huStarij = huni * hStarij * one_over_hiReg;
        double hvStarij = hvni * hStarij * one_over_hiReg;
        double hetaStarij = hetani * std::pow(hStarij * one_over_hiReg, 2);
        double hwStarij = hwni * std::pow(hStarij * one_over_hiReg, 2);

        double hStarji = fmax(0., hj + Zj - fmax(Zi, Zj));
        double huStarji = hunj * hStarji * one_over_hjReg;
        double hvStarji = hvnj * hStarji * one_over_hjReg;
        double hetaStarji = hetanj * std::pow(hStarji * one_over_hjReg, 2);
        double hwStarji = hwnj * std::pow(hStarij * one_over_hiReg, 2);

        // COMPUTE FLUX CORRECTION MATRICES
        double ML_minus_MC = (LUMPED_MASS_MATRIX == 1
                                  ? 0.
                                  : (i == j ? 1. : 0.) * mi - MassMatrix[ij]);
        double FluxCorrectionMatrix1 =
            ML_minus_MC * (high_order_hnp1[j] - hj - (high_order_hnp1i - hi)) +
            dt * (dH_minus_dL[ij] - muH_minus_muL[ij]) * (hStarji - hStarij) +
            dt * muH_minus_muL[ij] * (hj - hi);

        double FluxCorrectionMatrix2 =
            ML_minus_MC *
                (high_order_hunp1[j] - hunj - (high_order_hunp1i - huni)) +
            dt * (dH_minus_dL[ij] - muH_minus_muL[ij]) * (huStarji - huStarij) +
            dt * muH_minus_muL[ij] * (hunj - huni);

        double FluxCorrectionMatrix3 =
            ML_minus_MC *
                (high_order_hvnp1[j] - hvnj - (high_order_hvnp1i - hvni)) +
            dt * (dH_minus_dL[ij] - muH_minus_muL[ij]) * (hvStarji - hvStarij) +
            dt * muH_minus_muL[ij] * (hvnj - hvni);
        double FluxCorrectionMatrix4 =
            ML_minus_MC * (high_order_hetanp1[j] - hetanj -
                           (high_order_hetanp1i - hetani)) +
            dt * (dH_minus_dL[ij] - muH_minus_muL[ij]) *
                (hetaStarji - hetaStarij) +
            dt * muH_minus_muL[ij] * (hetanj - hetani);

        double FluxCorrectionMatrix5 =
            ML_minus_MC *
                (high_order_hwnp1[j] - hwnj - (high_order_hwnp1i - hwni)) +
            dt * (dH_minus_dL[ij] - muH_minus_muL[ij]) * (hwStarji - hwStarij) +
            dt * muH_minus_muL[ij] * (hwnj - hwni);

        // compute limiter based on water height
        double Lij = 1.;
        if (FluxCorrectionMatrix1 >= 0)
          Lij = fmin(Lij, std::min(Rneg[j], Rpos[i]));
        else
          Lij = fmin(Lij, std::min(Rneg[i], Rpos[j]));

        // COMPUTE LIMITER based on heta -EJT
        // Note that we set lij = min(lij_h,lij_heta)
        if (FluxCorrectionMatrix4 >= 0)
          Lij = fmin(Lij, std::min(Rneg_heta[j], Rpos_heta[i]));
        else
          Lij = fmin(Lij, std::min(Rneg_heta[i], Rpos_heta[j]));

        // CONVEX LIMITING // for kinetic energy
        // root of ith-DOF
        double huLow_minus_source_j =
            (huLow[j] - new_SourceTerm_hu[j] * dt / mi);
        double hvLow_minus_source_j =
            (hvLow[j] - new_SourceTerm_hv[j] * dt / mi);

        double lambdaj =
            csrRowIndeces_DofLoops[i + 1] - csrRowIndeces_DofLoops[i] - 1;
        double Ph_ij = FluxCorrectionMatrix1 / mi / lambdaj;
        double Phu_ij = FluxCorrectionMatrix2 / mi / lambdaj;
        double Phv_ij = FluxCorrectionMatrix3 / mi / lambdaj;

        double ai = -0.5 * (Phu_ij * Phu_ij + Phv_ij * Phv_ij);
        double bi = Kmax[i] * Ph_ij - (huLow_minus_source_i * Phu_ij +
                                       hvLow_minus_source_i * Phv_ij);

        double r1 = ai == 0
                        ? (bi == 0 ? 1. : -ci / bi)
                        : (-bi + std::sqrt(bi * bi - 4 * ai * ci)) / 2. / ai;
        double r2 = ai == 0
                        ? (bi == 0 ? 1. : -ci / bi)
                        : (-bi - std::sqrt(bi * bi - 4 * ai * ci)) / 2. / ai;
        if (r1 < 0 && r2 < 0) {
          r1 = 1.;
          r2 = 1.;
        }
        double ri = fabs(fmax(r1, r2));

        // root of jth-DOF (To compute transpose component)
        double lambdai =
            csrRowIndeces_DofLoops[j + 1] - csrRowIndeces_DofLoops[j] - 1;
        double mj = lumped_mass_matrix[j];
        double cj = Kmax[j] * hLow[j] -
                    0.5 * (huLow_minus_source_j * huLow_minus_source_j +
                           hvLow_minus_source_j * hvLow_minus_source_j);
        double Ph_ji = -FluxCorrectionMatrix1 / mj / lambdai; // Aij=-Aji
        double Phu_ji = -FluxCorrectionMatrix2 / mj / lambdai;
        double Phv_ji = -FluxCorrectionMatrix3 / mj / lambdai;
        double aj = -0.5 * (Phu_ji * Phu_ji + Phv_ji * Phv_ji);
        double bj = Kmax[j] * Ph_ji - (huLow_minus_source_j * Phu_ji +
                                       hvLow_minus_source_j * Phv_ji);

        r1 = aj == 0 ? (bj == 0 ? 1. : -cj / bj)
                     : (-bj + std::sqrt(bj * bj - 4 * aj * cj)) / 2. / aj;
        r2 = aj == 0 ? (bj == 0 ? 1. : -cj / bj)
                     : (-bj - std::sqrt(bj * bj - 4 * aj * cj)) / 2. / aj;
        if (r1 < 0 && r2 < 0) {
          r1 = 1.;
          r2 = 1.;
        }
        double rj = fabs(fmax(r1, r2));

        // COMPUTE LIMITER //
        Lij = fmin(fmin(ri, Lij), fmin(rj, Lij)); // Lij=Lji

        // COMPUTE LIMITED FLUX //
        ith_Limiter_times_FluxCorrectionMatrix1 += Lij * FluxCorrectionMatrix1;
        ith_Limiter_times_FluxCorrectionMatrix2 += Lij * FluxCorrectionMatrix2;
        ith_Limiter_times_FluxCorrectionMatrix3 += Lij * FluxCorrectionMatrix3;
        ith_Limiter_times_FluxCorrectionMatrix4 += Lij * FluxCorrectionMatrix4;
        ith_Limiter_times_FluxCorrectionMatrix5 += Lij * FluxCorrectionMatrix5;

        // update ij
        ij += 1;
      }
      double one_over_mi = 1.0 / lumped_mass_matrix[i];
      limited_hnp1[i] =
          hLow[i] + one_over_mi * ith_Limiter_times_FluxCorrectionMatrix1;
      limited_hunp1[i] =
          ((huLow[i] + dt / mi * new_SourceTerm_hu[i]) +
           one_over_mi * ith_Limiter_times_FluxCorrectionMatrix2);
      limited_hvnp1[i] =
          ((hvLow[i] + dt / mi * new_SourceTerm_hv[i]) +
           one_over_mi * ith_Limiter_times_FluxCorrectionMatrix3);
      limited_hetanp1[i] =
          ((hetaLow[i] + dt / mi * new_SourceTerm_heta[i]) +
           one_over_mi * ith_Limiter_times_FluxCorrectionMatrix4);
      limited_hwnp1[i] =
          ((hwLow[i] + dt / mi * new_SourceTerm_hw[i]) +
           one_over_mi * ith_Limiter_times_FluxCorrectionMatrix5);

      if (limited_hnp1[i] < -hEps && dt < 1.0) {
        std::cout << "Limited water height is negative: "
                  << "hLow: " << hLow[i] << "\t"
                  << "hHigh: " << limited_hnp1[i] << "\t"
                  << "hEps:  " << hEps << "\t"
                  << " ... aborting!" << std::endl;
        abort();
      } else {
        // clean up uHigh from round off error
        if (limited_hnp1[i] < hEps)
          limited_hnp1[i] = 0.0;
        // double aux = fmax(limited_hnp1[i], hEps); // hEps
        double aux =
            fmax(limited_hnp1[i], hEps); // hReg makes the code more robust
        limited_hunp1[i] *= 2 * std::pow(limited_hnp1[i], VEL_FIX_POWER) /
                            (std::pow(limited_hnp1[i], VEL_FIX_POWER) +
                             std::pow(aux, VEL_FIX_POWER));
        limited_hvnp1[i] *= 2 * std::pow(limited_hnp1[i], VEL_FIX_POWER) /
                            (std::pow(limited_hnp1[i], VEL_FIX_POWER) +
                             std::pow(aux, VEL_FIX_POWER));
        limited_hetanp1[i] *= 2 * std::pow(limited_hnp1[i], VEL_FIX_POWER) /
                              (std::pow(limited_hnp1[i], VEL_FIX_POWER) +
                               std::pow(aux, VEL_FIX_POWER));
        limited_hwnp1[i] *= 2 * std::pow(limited_hnp1[i], VEL_FIX_POWER) /
                            (std::pow(limited_hnp1[i], VEL_FIX_POWER) +
                             std::pow(aux, VEL_FIX_POWER));
      }
    }
  }

  double calculateEdgeBasedCFL(
      double g,
      int numDOFsPerEqn,          // number of DOFs
      double *lumped_mass_matrix, // lumped mass matrix (as vector))
      double *h_dof_old,          // DOFs of solution at last stage
      double *hu_dof_old, double *hv_dof_old, double *heta_dof_old,
      double *b_dof,
      int *csrRowIndeces_DofLoops,    // csr row indeces
      int *csrColumnOffsets_DofLoops, // csr column offsets
      double hEps, double *hReg, double *Cx, double *Cy, double *CTx,
      double *CTy, double *dLow, double run_cfl, double *edge_based_cfl,
      int debug) {

    /* note that for the CFL condition, we use only the values of dij and
     * don't do the dij = Max(dij,muij) thing */

    std::valarray<double> psi(numDOFsPerEqn);
    double max_edge_based_cfl = 0.;
    int ij = 0;
    for (int i = 0; i < numDOFsPerEqn; i++) {
      // solution at time tn for the ith DOF
      double hi = h_dof_old[i];
      double hui = hu_dof_old[i];
      double hvi = hv_dof_old[i];
      double hetai = heta_dof_old[i];
      double mi = lumped_mass_matrix[i];
      double dLowii = 0.;

      for (int offset = csrRowIndeces_DofLoops[i];
           offset < csrRowIndeces_DofLoops[i + 1]; offset++) {
        // loop in j (sparsity pattern)
        // solution at time tn for the jth DOF
        int j = csrColumnOffsets_DofLoops[offset];
        double hj = h_dof_old[j];
        double huj = hu_dof_old[j];
        double hvj = hv_dof_old[j];
        double hetaj = heta_dof_old[j];
        double mj = lumped_mass_matrix[j];

        if (i != j) {
          ////////////////////////
          // DISSIPATIVE MATRIX //
          ////////////////////////
          double cij_norm = sqrt(Cx[ij] * Cx[ij] + Cy[ij] * Cy[ij]);
          double cji_norm = sqrt(CTx[ij] * CTx[ij] + CTy[ij] * CTy[ij]);
          double nxij = Cx[ij] / cij_norm, nyij = Cy[ij] / cij_norm;
          double nxji = CTx[ij] / cji_norm, nyji = CTy[ij] / cji_norm;

          dLow[ij] =
              fmax(maxWaveSpeedSharpInitialGuess(g, nxij, nyij, hi, hui, hvi,
                                                 hetai, mi, hj, huj, hvj, hetaj,
                                                 mj, hEps, hEps, debug) *
                       cij_norm, // hEps
                   maxWaveSpeedSharpInitialGuess(g, nxji, nyji, hj, huj, hvj,
                                                 hetaj, mj, hi, hui, hvi, hetai,
                                                 mi, hEps, hEps, debug) *
                       cji_norm); // hEps
          dLowii -= dLow[ij];

        } else
          dLow[ij] = 0.;
        // muLow[ij] = 0.;
        // update ij
        ij += 1;
      }
      //////////////////////////////
      // CALCULATE EDGE BASED CFL //
      //////////////////////////////
      mi = lumped_mass_matrix[i];
      edge_based_cfl[i] = 1.0 * fabs(dLowii) / mi;
      max_edge_based_cfl = fmax(max_edge_based_cfl, edge_based_cfl[i]);

      // At the moment, we don't need to compute the smoothness indicator for
      // edge CFL
      //////////////////////////////////
      // COMPUTE SMOOTHNESS INDICATOR //
      //////////////////////////////////
    }

    return max_edge_based_cfl;
  } // End calculateEdgeBasedCFL

  void calculateResidual( // last EDGE BASED version
      double *mesh_trial_ref, double *mesh_grad_trial_ref, double *mesh_dof,
      double *mesh_velocity_dof, double MOVING_DOMAIN, int *mesh_l2g,
      double *dV_ref, double *h_trial_ref, double *h_grad_trial_ref,
      double *h_test_ref, double *h_grad_test_ref, double *vel_trial_ref,
      double *vel_grad_trial_ref, double *vel_test_ref,
      double *vel_grad_test_ref,
      // element boundary
      double *mesh_trial_trace_ref, double *mesh_grad_trial_trace_ref,
      double *dS_ref, double *h_trial_trace_ref, double *h_grad_trial_trace_ref,
      double *h_test_trace_ref, double *h_grad_test_trace_ref,
      double *vel_trial_trace_ref, double *vel_grad_trial_trace_ref,
      double *vel_test_trace_ref, double *vel_grad_test_trace_ref,
      double *normal_ref, double *boundaryJac_ref,
      // physics
      double *elementDiameter, int nElements_global, double useRBLES,
      double useMetrics, double alphaBDF, double nu, double g, int *h_l2g,
      int *vel_l2g, double *h_dof_old, double *hu_dof_old, double *hv_dof_old,
      double *heta_dof_old, double *hw_dof_old, double *b_dof, double *h_dof,
      double *hu_dof, double *hv_dof, double *heta_dof, double *hw_dof,
      double *h_dof_sge, double *hu_dof_sge, double *hv_dof_sge,
      double *heta_dof_sge, double *hw_dof_sge, double *q_mass_acc,
      double *q_mom_hu_acc, double *q_mom_hv_acc, double *q_mass_adv,
      double *q_mass_acc_beta_bdf, double *q_mom_hu_acc_beta_bdf,
      double *q_mom_hv_acc_beta_bdf, double *q_cfl, int *sdInfo_hu_hu_rowptr,
      int *sdInfo_hu_hu_colind, int *sdInfo_hu_hv_rowptr,
      int *sdInfo_hu_hv_colind, int *sdInfo_hv_hv_rowptr,
      int *sdInfo_hv_hv_colind, int *sdInfo_hv_hu_rowptr,
      int *sdInfo_hv_hu_colind, int offset_h, int offset_hu, int offset_hv,
      int offset_heta, int offset_hw, int stride_h, int stride_hu,
      int stride_hv, int stride_heta, int stride_hw, double *globalResidual,
      int nExteriorElementBoundaries_global,
      int *exteriorElementBoundariesArray, int *elementBoundaryElementsArray,
      int *elementBoundaryLocalElementBoundariesArray, int *isDOFBoundary_h,
      int *isDOFBoundary_hu, int *isDOFBoundary_hv,
      int *isAdvectiveFluxBoundary_h, int *isAdvectiveFluxBoundary_hu,
      int *isAdvectiveFluxBoundary_hv, int *isDiffusiveFluxBoundary_hu,
      int *isDiffusiveFluxBoundary_hv, double *ebqe_bc_h_ext,
      double *ebqe_bc_flux_mass_ext, double *ebqe_bc_flux_mom_hu_adv_ext,
      double *ebqe_bc_flux_mom_hv_adv_ext, double *ebqe_bc_hu_ext,
      double *ebqe_bc_flux_hu_diff_ext, double *ebqe_penalty_ext,
      double *ebqe_bc_hv_ext, double *ebqe_bc_flux_hv_diff_ext,
      double *q_velocity, double *ebqe_velocity, double *flux,
      double *elementResidual_h_save,
      // C matrices
      double *Cx, double *Cy, double *CTx, double *CTy,
      // PARAMETERS FOR EDGE BASED STABILIZATION
      int numDOFsPerEqn, int NNZ, int *csrRowIndeces_DofLoops,
      int *csrColumnOffsets_DofLoops,
      // LUMPED MASS MATRIX
      double *lumped_mass_matrix, double cfl_run, double hEps, double *hReg,
      // SAVE SOLUTION (mql)
      double *hnp1_at_quad_point, double *hunp1_at_quad_point,
      double *hvnp1_at_quad_point, double *hetanp1_at_quad_point,
      double *hwnp1_at_quad_point,
      // TO COMPUTE LOW ORDER
      double *extendedSourceTerm_hu, double *extendedSourceTerm_hv,
      double *extendedSourceTerm_heta, double *extendedSourceTerm_hw,
      // FOR FCT
      double *dH_minus_dL, double *muH_minus_muL, double cE,
      int LUMPED_MASS_MATRIX, double dt, int LINEAR_FRICTION, double mannings,
      // Quant of interests
      double *quantDOFs, int SECOND_CALL_CALCULATE_RESIDUAL,
      // NORMAL COMPONENTS
      int COMPUTE_NORMALS, double *normalx, double *normaly, double *dLow,
      // LOCAL LIMITING
      double *hBT, double *huBT, double *hvBT, double *hetaBT, double *hwBT,
      int lstage, double *new_SourceTerm_hu, double *new_SourceTerm_hv,
      double *new_SourceTerm_heta, double *new_SourceTerm_hw) {

    // FOR FRICTION//
    double n2 = std::pow(mannings, 2.);
    double gamma = 4. / 3;
    // mql. This parameter relaxes the cfl restriction.
    // It is now conservative. I might change it after the local bounds are
    // implemented
    double xi = 10.;

    //////////////////////////////////////
    // ********** CELL LOOPS ********** //
    //////////////////////////////////////
    // To compute:
    //      * Time derivative term
    //      * Cell based CFL
    //      * Velocity and soln at quad points (for other models)
    for (int eN = 0; eN < nElements_global; eN++) {
      // declare local storage for element residual and initialize
      register double elementResidual_h[nDOF_test_element],
          elementResidual_hu[nDOF_test_element],
          elementResidual_hv[nDOF_test_element],
          elementResidual_heta[nDOF_test_element],
          elementResidual_hw[nDOF_test_element];

      for (int i = 0; i < nDOF_test_element; i++) {
        elementResidual_h[i] = 0.0;
        elementResidual_hu[i] = 0.0;
        elementResidual_hv[i] = 0.0;
        elementResidual_heta[i] = 0.0;
        elementResidual_hw[i] = 0.0;
      }
      //
      // loop over quadrature points and compute integrands
      //
      for (int k = 0; k < nQuadraturePoints_element; k++) {
        // compute indices and declare local storage
        register int eN_k = eN * nQuadraturePoints_element + k,
                     eN_k_nSpace = eN_k * nSpace,
                     eN_nDOF_trial_element = eN * nDOF_trial_element;
        register double h = 0.0, hu = 0.0, hv = 0.0, heta = 0.0,
                        hw = 0.0, // solution at current time
            h_old = 0.0, hu_old = 0.0, hv_old = 0.0, heta_old = 0.0,
                        hw_old = 0.0, // solution at lstage
            jac[nSpace * nSpace], jacDet, jacInv[nSpace * nSpace],
                        h_test_dV[nDOF_trial_element], dV, x, y, xt, yt;
        // get jacobian, etc for mapping reference element
        ck.calculateMapping_element(eN, k, mesh_dof, mesh_l2g, mesh_trial_ref,
                                    mesh_grad_trial_ref, jac, jacDet, jacInv, x,
                                    y);
        // get the physical integration weight
        dV = fabs(jacDet) * dV_ref[k];
        // get the solution at current time
        ck.valFromDOF(h_dof, &h_l2g[eN_nDOF_trial_element],
                      &h_trial_ref[k * nDOF_trial_element], h);
        ck.valFromDOF(hu_dof, &vel_l2g[eN_nDOF_trial_element],
                      &vel_trial_ref[k * nDOF_trial_element], hu);
        ck.valFromDOF(hv_dof, &vel_l2g[eN_nDOF_trial_element],
                      &vel_trial_ref[k * nDOF_trial_element], hv);
        ck.valFromDOF(heta_dof, &vel_l2g[eN_nDOF_trial_element],
                      &vel_trial_ref[k * nDOF_trial_element], heta);
        ck.valFromDOF(hw_dof, &vel_l2g[eN_nDOF_trial_element],
                      &vel_trial_ref[k * nDOF_trial_element], hw);
        // get the solution at the lstage
        ck.valFromDOF(h_dof_old, &h_l2g[eN_nDOF_trial_element],
                      &h_trial_ref[k * nDOF_trial_element], h_old);
        ck.valFromDOF(hu_dof_old, &vel_l2g[eN_nDOF_trial_element],
                      &vel_trial_ref[k * nDOF_trial_element], hu_old);
        ck.valFromDOF(hv_dof_old, &vel_l2g[eN_nDOF_trial_element],
                      &vel_trial_ref[k * nDOF_trial_element], hv_old);
        ck.valFromDOF(heta_dof_old, &vel_l2g[eN_nDOF_trial_element],
                      &vel_trial_ref[k * nDOF_trial_element], heta_old);
        ck.valFromDOF(hw_dof_old, &vel_l2g[eN_nDOF_trial_element],
                      &vel_trial_ref[k * nDOF_trial_element], hw_old);
        // calculate cell based CFL to keep a reference
        calculateCFL(elementDiameter[eN], g, h_old, hu_old, hv_old, hEps,
                     q_cfl[eN_k]);
        // precalculate test function products with integration weights
        for (int j = 0; j < nDOF_trial_element; j++)
          h_test_dV[j] = h_test_ref[k * nDOF_trial_element + j] * dV;

        // SAVE VELOCITY // at quadrature points for other models to use
        q_velocity[eN_k_nSpace + 0] =
            2 * h / (h * h + std::pow(fmax(h, hEps), 2)) * hu;
        q_velocity[eN_k_nSpace + 1] =
            2 * h / (h * h + std::pow(fmax(h, hEps), 2)) * hv;
        hnp1_at_quad_point[eN_k] = h;
        hunp1_at_quad_point[eN_k] = hu;
        hvnp1_at_quad_point[eN_k] = hv;
        hetanp1_at_quad_point[eN_k] = heta;
        hwnp1_at_quad_point[eN_k] = hw;

        for (int i = 0; i < nDOF_test_element; i++) {
          // compute time derivative part of global residual. NOTE: no lumping
          elementResidual_h[i] += (h - h_old) * h_test_dV[i];
          elementResidual_hu[i] += (hu - hu_old) * h_test_dV[i];
          elementResidual_hv[i] += (hv - hv_old) * h_test_dV[i];
          elementResidual_heta[i] += (heta - heta_old) * h_test_dV[i];
          elementResidual_hw[i] += (hw - hw_old) * h_test_dV[i];
        }
      }
      // distribute
      for (int i = 0; i < nDOF_test_element; i++) {
        register int eN_i = eN * nDOF_test_element + i;
        int h_gi = h_l2g[eN_i];     // global i-th index for h
        int vel_gi = vel_l2g[eN_i]; // global i-th index for velocities

        // distribute time derivative to global residual
        globalResidual[offset_h + stride_h * h_gi] += elementResidual_h[i];
        globalResidual[offset_hu + stride_hu * vel_gi] += elementResidual_hu[i];
        globalResidual[offset_hv + stride_hv * vel_gi] += elementResidual_hv[i];
        globalResidual[offset_heta + stride_heta * vel_gi] +=
            elementResidual_heta[i];
        globalResidual[offset_hw + stride_hw * vel_gi] += elementResidual_hw[i];
      }
    }
    // ********** END OF CELL LOOPS ********** //

    if (SECOND_CALL_CALCULATE_RESIDUAL == 0) // This is to save some time
    {
      //////////////////////////////////////////////
      // ********** FIRST LOOP ON DOFs ********** //
      //////////////////////////////////////////////

      // To compute:
      //     * Entropy at i-th node
      std::valarray<double> eta(numDOFsPerEqn);
      for (int i = 0; i < numDOFsPerEqn; i++) {
        // COMPUTE ENTROPY. NOTE: WE CONSIDER A FLAT BOTTOM
        double hi = h_dof_old[i];
        double one_over_hiReg =
            2 * hi / (hi * hi + std::pow(fmax(hi, hEps), 2)); // hEps
        eta[i] =
            ENTROPY(g, hi, hu_dof_old[i], hv_dof_old[i], 0., one_over_hiReg);
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
      std::valarray<double> hyp_flux_h(numDOFsPerEqn),
          hyp_flux_hu(numDOFsPerEqn), hyp_flux_hv(numDOFsPerEqn),
          hyp_flux_heta(numDOFsPerEqn), hyp_flux_hw(numDOFsPerEqn),
          global_entropy_residual(numDOFsPerEqn), psi(numDOFsPerEqn),
          etaMax(numDOFsPerEqn), etaMin(numDOFsPerEqn);

      for (int i = 0; i < numDOFsPerEqn; i++) {

        // solution at time tn for the ith DOF
        double hi = h_dof_old[i];
        double hui = hu_dof_old[i];
        double hvi = hv_dof_old[i];
        double hetai = heta_dof_old[i];
        double hwi = hw_dof_old[i];
        double Zi = b_dof[i];
        // Define some things using above
        double one_over_hiReg =
            2 * hi / (hi * hi + std::pow(fmax(hi, hEps), 2)); // hEps
        double ui = hui * one_over_hiReg;
        double vi = hvi * one_over_hiReg;
        double etai = hetai * one_over_hiReg;
        double mi = lumped_mass_matrix[i];
        double meshSizei = std::sqrt(mi);

        // For eta min and max
        etaMax[i] = fabs(eta[i]);
        etaMin[i] = fabs(eta[i]);

        /* COMPUTE EXTENDED SOURCE TERMS:
         * Friction terms
         * Source for heta equation
         * Source for hw equation
         * NOTE: Be careful with sign of source terms.
         */

        // Friction
        if (LINEAR_FRICTION == 1) {
          extendedSourceTerm_hu[i] = mannings * hui * mi;
          extendedSourceTerm_hv[i] = mannings * hvi * mi;
          // For use in the convex limiting function -EJT
          new_SourceTerm_hu[i] = -mannings * hui * mi;
          new_SourceTerm_hv[i] = -mannings * hvi * mi;
        } else {
          double veli_norm = std::sqrt(ui * ui + vi * vi);
          double hi_to_the_gamma = std::pow(fmax(hi, hEps), gamma);
          double friction_aux =
              veli_norm == 0.
                  ? 0.
                  : (2 * g * n2 * veli_norm * mi /
                     (hi_to_the_gamma +
                      fmax(hi_to_the_gamma, xi * g * n2 * dt * veli_norm)));
          extendedSourceTerm_hu[i] = friction_aux * hui;
          extendedSourceTerm_hv[i] = friction_aux * hvi;
          // For use in the convex limiting functio -EJT
          new_SourceTerm_hu[i] = -friction_aux * hui;
          new_SourceTerm_hv[i] = -friction_aux * hvi;
        }

        // Define some things for heta and hw sources
        double ratioi = (2.0 * hetai) / (etai * etai + hi * hi + hEps);
        double diff_over_h_i = (hetai - hi * hi) * one_over_hiReg;
        double hSqd_GammaPi = 6.0 * (hetai - hi * hi);
        if (IF_BOTH_GAMMA_BRANCHES) {
          if (hetai > std::pow(hi, 2.0)) {
            hSqd_GammaPi = 6.0 * etai * diff_over_h_i;
          }
        }

        // heta source
        extendedSourceTerm_heta[i] = -hwi * mi * ratioi;
        new_SourceTerm_heta[i] = hwi * mi * ratioi;

        // hw source
        extendedSourceTerm_hw[i] =
            (LAMBDA_MGN * g / meshSizei) * hSqd_GammaPi * mi * ratioi;
        new_SourceTerm_hw[i] =
            -(LAMBDA_MGN * g / meshSizei) * hSqd_GammaPi * mi * ratioi;

        /* HYPERBOLIC FLUXES */
        hyp_flux_h[i] = 0;
        hyp_flux_hu[i] = 0;
        hyp_flux_hv[i] = 0;
        hyp_flux_heta[i] = 0;
        hyp_flux_hw[i] = 0;

        // FOR ENTROPY RESIDUAL //
        double ith_flux_term1 = 0., ith_flux_term2 = 0., ith_flux_term3 = 0.;
        double ith_flux_term4 = 0., ith_flux_term5 = 0.;
        double entropy_flux = 0.;

        // FOR SMOOTHNESS INDICATOR //
        double alphai; // smoothness indicator of solution
        double alpha_numerator = 0;
        double alpha_denominator = 0;
        double alpha_zero = 0.75;
        double alpha_factor = (1.0 / (1.0 - alpha_zero));

        // loop in j (sparsity pattern)
        for (int offset = csrRowIndeces_DofLoops[i];
             offset < csrRowIndeces_DofLoops[i + 1]; offset++) {

          int j = csrColumnOffsets_DofLoops[offset];

          // solution at time tn for the jth DOF
          double hj = h_dof_old[j];
          double huj = hu_dof_old[j];
          double hvj = hv_dof_old[j];
          double hetaj = heta_dof_old[j];
          double hwj = hw_dof_old[j];
          double Zj = b_dof[j];

          // Then define some things here using above
          double one_over_hjReg =
              2 * hj / (hj * hj + std::pow(fmax(hj, hEps), 2)); // hEps
          double uj = huj * one_over_hjReg;
          double vj = hvj * one_over_hjReg;
          double etaj = hetaj * one_over_hjReg;
          double meshSizej =
              std::sqrt(lumped_mass_matrix[j]); // local mesh size in 2d

          // pTilde at jth node gets defined here
          double diff_over_h_j = (hetaj - hj * hj) * one_over_hjReg;
          double pTildej = -(LAMBDA_MGN * g / (3.0 * meshSizej)) * 6.0 * hj *
                           (hetaj - hj * hj);

          if (IF_BOTH_GAMMA_BRANCHES) {
            if (hetaj > std::pow(hj, 2.0)) {
              pTildej = -(LAMBDA_MGN * g / (3.0 * meshSizej)) * 2.0 *
                        diff_over_h_j * (etaj * etaj + etaj * hj + hj * hj);
            }
          }

          // auxiliary functions to compute fluxes
          double aux_h =
              (uj * hj - ui * hi) * Cx[ij] + (vj * hi - vi * hi) * Cy[ij];
          double aux_hu =
              (uj * huj - ui * hui) * Cx[ij] + (vj * huj - vi * hui) * Cy[ij];
          double aux_hv =
              (uj * hvj - ui * hvi) * Cx[ij] + (vj * hvj - vi * hvi) * Cy[ij];
          double aux_heta = (uj * hetaj - ui * hetai) * Cx[ij] +
                            (vj * hetaj - vi * hetai) * Cy[ij];
          double aux_hw =
              (uj * hwj - ui * hwi) * Cx[ij] + (vj * hwj - vi * hwi) * Cy[ij];

          /* HYPERBOLIC FLUX */
          hyp_flux_h[i] += aux_h;
          hyp_flux_hu[i] += aux_hu + pTildej * Cx[ij];
          hyp_flux_hv[i] += aux_hv + pTildej * Cy[ij];
          hyp_flux_heta[i] += aux_heta;
          hyp_flux_hw[i] += aux_hw;

          // EXTENDED SOURCE, USING 6.13 //
          extendedSourceTerm_hu[i] += g * hi * (hj + Zj) * Cx[ij];
          extendedSourceTerm_hv[i] += g * hi * (hj + Zj) * Cy[ij];

          new_SourceTerm_hu[i] +=
              g * (-hi * (Zj - Zi) + 0.5 * std::pow(hj - hi, 2)) * Cx[ij];
          new_SourceTerm_hv[i] +=
              g * (-hi * (Zj - Zi) + 0.5 * std::pow(hj - hi, 2)) * Cy[ij];

          // flux for entropy
          ith_flux_term1 += aux_h;
          ith_flux_term2 += aux_hu + g * hi * (hj + 0.) * Cx[ij]; // NOTE: Zj=0
          ith_flux_term3 += aux_hv + g * hi * (hj + 0.) * Cy[ij]; // NOTE: Zj=0

          // NOTE: WE CONSIDER FLAT BOTTOM
          entropy_flux +=
              (Cx[ij] * ENTROPY_FLUX1(g, hj, huj, hvj, 0., one_over_hjReg) +
               Cy[ij] * ENTROPY_FLUX2(g, hj, huj, hvj, 0., one_over_hjReg));

          // COMPUTE ETA MIN AND ETA MAX //
          etaMax[i] = fmax(etaMax[i], fabs(eta[j]));
          etaMin[i] = fmin(etaMin[i], fabs(eta[j]));

          // FOR SMOOTHNESS INDICATOR //
          alpha_numerator += hj - hi;
          alpha_denominator += fabs(hj - hi);

          // update ij
          ij += 1;
        }

        // COMPUTE ENTROPY RESIDUAL //
        double one_over_entNormFactori = 2. / (etaMax[i] - etaMin[i] + 1E-15);
        double eta_prime1 = DENTROPY_DH(g, hi, hui, hvi, 0.,
                                        one_over_hiReg); // NOTE: FLAT BOTTOM
        double eta_prime2 = DENTROPY_DHU(g, hi, hui, hvi, 0., one_over_hiReg);
        double eta_prime3 = DENTROPY_DHV(g, hi, hui, hvi, 0., one_over_hiReg);
        global_entropy_residual[i] =
            one_over_entNormFactori *
            fabs(entropy_flux -
                 (ith_flux_term1 * eta_prime1 + ith_flux_term2 * eta_prime2 +
                  ith_flux_term3 * eta_prime3));

        // COMPUTE SMOOTHNESS INDICATOR //
        if (hi <= hEps) {
          alphai = 1.;
          global_entropy_residual[i] = 1E10;
        } else {
          // Force alphai=0 in constant states
          if (fabs(alpha_numerator) <= hEps) {
            alphai = 0.;
          } else {
            alphai =
                (fabs(alpha_numerator) - hEps) / fabs(alpha_denominator - hEps);
          }
          alphai = fmax(alphai - alpha_zero, 0.0) * alpha_factor;
        }
        if (POWER_SMOOTHNESS_INDICATOR == 0)
          psi[i] = 1.0;
        else
          psi[i] = std::pow(
              alphai, POWER_SMOOTHNESS_INDICATOR); // NOTE: alpha^4 for mGN
      }
      // ********** END OF 2nd LOOP ON DOFS ********** //

      /////////////////////////////////////////////
      // ********** MAIN LOOP ON DOFs **********
      // To compute:
      //      * dissipative terms
      //      * bar states
      /////////////////////////////////////////////

      ij = 0;
      for (int i = 0; i < numDOFsPerEqn; i++) {

        double hi = h_dof_old[i];
        double hui = hu_dof_old[i];
        double hvi = hv_dof_old[i];
        double hetai = heta_dof_old[i];
        double hwi = hw_dof_old[i];
        double Zi = b_dof[i];
        double mi = lumped_mass_matrix[i];

        double one_over_hiReg =
            2 * hi / (hi * hi + std::pow(fmax(hi, hEps), 2));
        double ui = hui * one_over_hiReg;
        double vi = hvi * one_over_hiReg;
        double etai = hetai * one_over_hiReg;
        double meshSizei = std::sqrt(mi); // local mesh size in 2d

        // We define pTilde at ith node here
        double diff_over_h_i = (hetai - hi * hi) * one_over_hiReg;
        double pTildei = -(LAMBDA_MGN * g / (3.0 * meshSizei)) * 6.0 * hi *
                         (hetai - hi * hi);

        if (IF_BOTH_GAMMA_BRANCHES) {
          if (hetai > std::pow(hi, 2.0)) {
            pTildei = -(LAMBDA_MGN * g / (3.0 * meshSizei)) * 2.0 *
                      diff_over_h_i * (etai * etai + etai * hi + hi * hi);
          }
        }

        // Define full pressure at ith node for definition of bar states below
        double pressure_i = 0.5 * g * hi * hi + pTildei;

        // HIGH ORDER DISSIPATIVE TERMS, for Aij matrix
        double ith_dHij_minus_muHij_times_hStarStates = 0.,
               ith_dHij_minus_muHij_times_huStarStates = 0.,
               ith_dHij_minus_muHij_times_hvStarStates = 0.,
               ith_dHij_minus_muHij_times_hetaStarStates = 0.,
               ith_dHij_minus_muHij_times_hwStarStates = 0.,
               ith_muHij_times_hStates = 0., ith_muHij_times_huStates = 0.,
               ith_muHij_times_hvStates = 0., ith_muHij_times_hetaStates = 0.,
               ith_muHij_times_hwStates = 0.;

        // loop over the sparsity pattern of the i-th DOF
        for (int offset = csrRowIndeces_DofLoops[i];
             offset < csrRowIndeces_DofLoops[i + 1]; offset++) {

          int j = csrColumnOffsets_DofLoops[offset];

          double hj = h_dof_old[j];
          double huj = hu_dof_old[j];
          double hvj = hv_dof_old[j];
          double hetaj = heta_dof_old[j];
          double hwj = hw_dof_old[j];
          double Zj = b_dof[j];

          double one_over_hjReg =
              2 * hj / (hj * hj + std::pow(fmax(hj, hEps), 2));
          double uj = huj * one_over_hjReg;
          double vj = hvj * one_over_hjReg;
          double etaj = hetaj * one_over_hjReg;
          double wj = hwj * one_over_hjReg;
          double mj = lumped_mass_matrix[j];
          double meshSizej = std::sqrt(mj); // local mesh size in 2d

          // Here we define pTilde at jth node
          double diff_over_h_j = (hetaj - hj * hj) * one_over_hjReg;
          double pTildej = -(LAMBDA_MGN * g / (3.0 * meshSizej)) * 6.0 * hj *
                           (hetaj - hj * hj);

          if (IF_BOTH_GAMMA_BRANCHES) {
            if (hetaj > std::pow(hj, 2.0)) {
              pTildej = -(LAMBDA_MGN * g / (3.0 * meshSizej)) * 2.0 *
                        diff_over_h_j * (etaj * etaj + etaj * hj + hj * hj);
            }
          }

          // define pressure at jth node
          double pressure_j = 0.5 * g * hj * hj + pTildej;

          // COMPUTE STAR SOLUTION // hStar, huStar, hvStar, hetaStar, hwStar
          double hStarij = fmax(0., hi + Zi - fmax(Zi, Zj));
          double huStarij = hui * hStarij * one_over_hiReg;
          double hvStarij = hvi * hStarij * one_over_hiReg;
          double hetaStarij = hetai * std::pow(hStarij * one_over_hiReg, 2);
          double hwStarij = hwi * std::pow(hStarij * one_over_hiReg, 2);

          double hStarji = fmax(0., hj + Zj - fmax(Zi, Zj));
          double huStarji = huj * hStarji * one_over_hjReg;
          double hvStarji = hvj * hStarji * one_over_hjReg;
          double hetaStarji = hetaj * std::pow(hStarji * one_over_hjReg, 2);
          double hwStarji = hwj * std::pow(hStarji * one_over_hiReg, 2);

          // Dissipative well balancing term
          double muLowij = 0., muLij = 0., muHij = 0.;
          double dLowij = 0., dLij = 0., dHij = 0.;
          if (i != j) // This is not necessary. See formula for
                      // ith_dissipative_terms
          {
            ////////////////////////
            // DISSIPATIVE MATRIX //
            ////////////////////////
            if (lstage == 0)
              dLowij = dLow[ij];
            else {
              double cij_norm = sqrt(Cx[ij] * Cx[ij] + Cy[ij] * Cy[ij]);
              double cji_norm = sqrt(CTx[ij] * CTx[ij] + CTy[ij] * CTy[ij]);
              double nxij = Cx[ij] / cij_norm, nyij = Cy[ij] / cij_norm;
              double nxji = CTx[ij] / cji_norm, nyji = CTy[ij] / cji_norm;
              dLowij = fmax(maxWaveSpeedSharpInitialGuess(
                                g, nxij, nyij, hi, hui, hvi, hetai, mi, hj, huj,
                                hvj, hetaj, mj, hEps, hEps, false) *
                                cij_norm,
                            maxWaveSpeedSharpInitialGuess(
                                g, nxji, nyji, hj, huj, hvj, hetaj, mj, hi, hui,
                                hvi, hetai, mi, hEps, hEps, false) *
                                cji_norm);
            }

            ///////////////////////////////////////
            // WELL BALANCING DISSIPATIVE MATRIX //
            ///////////////////////////////////////
            muLowij = fmax(fmax(0., -(ui * Cx[ij] + vi * Cy[ij])),
                           fmax(0, (uj * Cx[ij] + vj * Cy[ij])));

            // Define dLij as low order dijs
            muLij = muLowij;
            dLij = fmax(dLowij, muLij);

            // Then save dLow for limiting step, maybe a bit confusing
            dLow[ij] = fmax(dLij, muLij);

            ////////////////////////
            // COMPUTE BAR STATES //
            ////////////////////////
            double hBar_ij = 0, hTilde_ij = 0, huBar_ij = 0, huTilde_ij = 0,
                   hvBar_ij = 0, hvTilde_ij = 0, hetaBar_ij = 0,
                   hetaTilde_ij = 0, hwBar_ij = 0, hwTilde_ij = 0;

            if (dLij != 0) {
              // h component
              hBar_ij = -1. / (2 * dLij) *
                            ((uj * hj - ui * hi) * Cx[ij] +
                             (vj * hj - vi * hi) * Cy[ij]) +
                        0.5 * (hj + hi);
              hTilde_ij = (dLij - muLij) / (2 * dLij) *
                          ((hStarji - hj) - (hStarij - hi));
              // hu component
              huBar_ij = -1. / (2 * dLij) *
                             ((uj * huj - ui * hui + pressure_j - pressure_i) *
                                  Cx[ij] +
                              (vj * huj - vi * hui) * Cy[ij]) +
                         0.5 * (huj + hui);
              huTilde_ij = (dLij - muLij) / (2 * dLij) *
                           ((huStarji - huj) - (huStarij - hui));
              // hv component
              hvBar_ij = -1. / (2 * dLij) *
                             ((uj * hvj - ui * hvi) * Cx[ij] +
                              (vj * hvj - vi * hvi + pressure_j - pressure_i) *
                                  Cy[ij]) +
                         0.5 * (hvj + hvi);
              hvTilde_ij = (dLij - muLij) / (2 * dLij) *
                           ((hvStarji - hvj) - (hvStarij - hvi));
              // heta component
              hetaBar_ij = -1. / (2 * dLij) *
                               ((uj * hetaj - ui * hetai) * Cx[ij] +
                                (vj * hetaj - vi * hetai) * Cy[ij]) +
                           0.5 * (hetaj + hetai);
              hetaTilde_ij = (dLij - muLij) / (2 * dLij) *
                             ((hetaStarji - hetaj) - (hetaStarij - hetai));
              // hw component
              hwBar_ij = -1. / (2 * dLij) *
                             ((uj * hwj - ui * hwi) * Cx[ij] +
                              (vj * hwj - vi * hwi) * Cy[ij]) +
                         0.5 * (hwj + hwi);
              hwTilde_ij = (dLij - muLij) / (2 * dLij) *
                           ((hwStarji - hwj) - (hwStarij - hwi));
            }
            // Here we define uBar + uTilde
            hBT[ij] = hBar_ij + hTilde_ij;
            huBT[ij] = huBar_ij + huTilde_ij;
            hvBT[ij] = hvBar_ij + hvTilde_ij;
            hetaBT[ij] = hetaBar_ij + hetaTilde_ij;
            hwBT[ij] = hwBar_ij + hwTilde_ij;

            ///////////////////////
            // ENTROPY VISCOSITY //
            ///////////////////////
            double dEVij = cE * fmax(global_entropy_residual[i],
                                     global_entropy_residual[j]);
            // dHij = fmin(dLowij, dEVij);
            // muHij = fmin(muLowij, dEVij);

            // Assume no EV for now and just use the alpha limiting for
            // higher order method. -EJT
            dHij = fmax(psi[i], psi[j]) * dLij;
            muHij = fmax(psi[i], psi[j]) * muLij;

            // compute dij_minus_muij times star solution terms
            // see: eqn (6.13)
            ith_dHij_minus_muHij_times_hStarStates +=
                (dHij - muHij) * (hStarji - hStarij);
            ith_dHij_minus_muHij_times_huStarStates +=
                (dHij - muHij) * (huStarji - huStarij);
            ith_dHij_minus_muHij_times_hvStarStates +=
                (dHij - muHij) * (hvStarji - hvStarij);
            ith_dHij_minus_muHij_times_hetaStarStates +=
                (dHij - muHij) * (hetaStarji - hetaStarij);
            ith_dHij_minus_muHij_times_hwStarStates +=
                (dHij - muHij) * (hwStarji - hwStarij);

            // compute muij times solution terms
            ith_muHij_times_hStates += muHij * (hj - hi);
            ith_muHij_times_huStates += muHij * (huj - hui);
            ith_muHij_times_hvStates += muHij * (hvj - hvi);
            ith_muHij_times_hetaStates += muHij * (hetaj - hetai);
            ith_muHij_times_hwStates += muHij * (hwj - hwi);

            // compute dH_minus_dL
            dH_minus_dL[ij] = dHij - dLij;
            muH_minus_muL[ij] = muHij - muLij;
          } else // i==j
          {
            dH_minus_dL[ij] = 0.;   // Not true but the prod of this times
                                    // Uj-Ui will be zero
            muH_minus_muL[ij] = 0.; // Not true but the prod of this times
                                    // Uj-Ui will be zero
          }

          // update ij
          ij += 1;

        } /*j loop ends here*/

        /* Define global residual */
        if (LUMPED_MASS_MATRIX == 1) {
          globalResidual[offset_h + stride_h * i] =
              hi - dt / mi *
                       (hyp_flux_h[i] - ith_dHij_minus_muHij_times_hStarStates -
                        ith_muHij_times_hStates);
          globalResidual[offset_hu + stride_hu * i] =
              hui - dt / mi *
                        ((hyp_flux_hu[i] + extendedSourceTerm_hu[i]) -
                         ith_dHij_minus_muHij_times_huStarStates -
                         ith_muHij_times_huStates);
          globalResidual[offset_hv + stride_hv * i] =
              hvi - dt / mi *
                        ((hyp_flux_hv[i] + extendedSourceTerm_hv[i]) -
                         ith_dHij_minus_muHij_times_hvStarStates -
                         ith_muHij_times_hvStates);
          globalResidual[offset_heta + stride_heta * i] =
              hetai - dt / mi *
                          ((hyp_flux_heta[i] + extendedSourceTerm_heta[i]) -
                           ith_dHij_minus_muHij_times_hetaStarStates -
                           ith_muHij_times_hetaStates);
          globalResidual[offset_hw + stride_hw * i] =
              hwi - dt / mi *
                        ((hyp_flux_hw[i] + extendedSourceTerm_hw[i]) -
                         ith_dHij_minus_muHij_times_hwStarStates -
                         ith_muHij_times_hwStates);

          // clean up potential negative water height due to machine
          // precision
          if (globalResidual[offset_h + stride_h * i] >= -hEps &&
              globalResidual[offset_h + stride_h * i] < hEps)
            globalResidual[offset_h + stride_h * i] = 0;

        } else {
          // Distribute residual
          // NOTE: MASS MATRIX IS CONSISTENT
          globalResidual[offset_h + stride_h * i] +=
              dt * (hyp_flux_h[i] - ith_dHij_minus_muHij_times_hStarStates -
                    ith_muHij_times_hStates);
          globalResidual[offset_hu + stride_hu * i] +=
              dt * (hyp_flux_hu[i] + extendedSourceTerm_hu[i] -
                    ith_dHij_minus_muHij_times_huStarStates -
                    ith_muHij_times_huStates);
          globalResidual[offset_hv + stride_hv * i] +=
              dt * (hyp_flux_hv[i] + extendedSourceTerm_hv[i] -
                    ith_dHij_minus_muHij_times_hvStarStates -
                    ith_muHij_times_hvStates);
          globalResidual[offset_heta + stride_heta * i] +=
              dt * (hyp_flux_heta[i] + extendedSourceTerm_heta[i] -
                    ith_dHij_minus_muHij_times_hetaStarStates -
                    ith_muHij_times_hetaStates);
          globalResidual[offset_hw + stride_hw * i] +=
              dt * (hyp_flux_hw[i] + extendedSourceTerm_hw[i] -
                    ith_dHij_minus_muHij_times_hwStarStates -
                    ith_muHij_times_hwStates);
        }
      }
      // ********** END OF LOOP IN DOFs ********** //
    }

    // ********** COMPUTE NORMALS ********** //
    if (COMPUTE_NORMALS == 1) {
      // This is to identify the normals and create a vector of normal
      // components
      for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) {
        register int
            ebN = exteriorElementBoundariesArray[ebNE],
            eN = elementBoundaryElementsArray[ebN * 2 + 0],
            ebN_local = elementBoundaryLocalElementBoundariesArray[ebN * 2 + 0];
        register double normal[3];
        {             // "Loop" in quad points
          int kb = 0; // NOTE: I need to consider just one quad point since
                      // the element is not curved so the normal is constant
                      // per element
          register int ebN_local_kb =
              ebN_local * nQuadraturePoints_elementBoundary + kb;
          register double jac_ext[nSpace * nSpace], jacDet_ext,
              jacInv_ext[nSpace * nSpace], boundaryJac[nSpace * (nSpace - 1)],
              metricTensor[(nSpace - 1) * (nSpace - 1)], metricTensorDetSqrt,
              x_ext, y_ext;
          /* compute information about mapping from reference element to
           * physical element */
          ck.calculateMapping_elementBoundary(
              eN, ebN_local, kb, ebN_local_kb, mesh_dof, mesh_l2g,
              mesh_trial_trace_ref, mesh_grad_trial_trace_ref, boundaryJac_ref,
              jac_ext, jacDet_ext, jacInv_ext, boundaryJac, metricTensor,
              metricTensorDetSqrt, normal_ref, normal, x_ext, y_ext);
        }
        // distribute the normal vectors
        for (int i = 0; i < nDOF_test_element; i++) {
          int eN_i = eN * nDOF_test_element + i;
          int gi = h_l2g[eN_i];
          normalx[gi] += 0.5 * normal[0] * (i == ebN_local ? 0. : 1.);
          normaly[gi] += 0.5 * normal[1] * (i == ebN_local ? 0. : 1.);
        }
      }
      // normalize
      for (int gi = 0; gi < numDOFsPerEqn; gi++) {
        double norm_factor =
            sqrt(std::pow(normalx[gi], 2) + std::pow(normaly[gi], 2));
        if (norm_factor != 0) {
          normalx[gi] /= norm_factor;
          normaly[gi] /= norm_factor;
        }
      }
    }
    // ********** END OF COMPUTING NORMALS ********** //
  } // namespace proteus

  void calculateMassMatrix( // element
      double *mesh_trial_ref, double *mesh_grad_trial_ref, double *mesh_dof,
      double *mesh_velocity_dof, double MOVING_DOMAIN, int *mesh_l2g,
      double *dV_ref, double *h_trial_ref, double *h_grad_trial_ref,
      double *h_test_ref, double *h_grad_test_ref, double *vel_trial_ref,
      double *vel_grad_trial_ref, double *vel_test_ref,
      double *vel_grad_test_ref,
      // element boundary
      double *mesh_trial_trace_ref, double *mesh_grad_trial_trace_ref,
      double *dS_ref, double *h_trial_trace_ref, double *h_grad_trial_trace_ref,
      double *h_test_trace_ref, double *h_grad_test_trace_ref,
      double *vel_trial_trace_ref, double *vel_grad_trial_trace_ref,
      double *vel_test_trace_ref, double *vel_grad_test_trace_ref,
      double *normal_ref, double *boundaryJac_ref,
      // physics
      double *elementDiameter, int nElements_global, double useRBLES,
      double useMetrics, double alphaBDF, double nu, double g, int *h_l2g,
      int *vel_l2g, double *b_dof, double *h_dof, double *hu_dof,
      double *hv_dof, double *heta_dof, double *hw_dof, double *h_dof_sge,
      double *hu_dof_sge, double *hv_dof_sge, double *heta_dof_sge,
      double *hw_dof_sge, double *q_mass_acc_beta_bdf,
      double *q_mom_hu_acc_beta_bdf, double *q_mom_hv_acc_beta_bdf,
      double *q_cfl, int *sdInfo_hu_hu_rowptr, int *sdInfo_hu_hu_colind,
      int *sdInfo_hu_hv_rowptr, int *sdInfo_hu_hv_colind,
      int *sdInfo_hv_hv_rowptr, int *sdInfo_hv_hv_colind,
      int *sdInfo_hv_hu_rowptr, int *sdInfo_hv_hu_colind,
      // h
      int *csrRowIndeces_h_h, int *csrColumnOffsets_h_h,
      int *csrRowIndeces_h_hu, int *csrColumnOffsets_h_hu,
      int *csrRowIndeces_h_hv, int *csrColumnOffsets_h_hv,
      int *csrRowIndeces_h_heta, int *csrColumnOffsets_h_heta,
      int *csrRowIndeces_h_hw, int *csrColumnOffsets_h_hw,
      // hu
      int *csrRowIndeces_hu_h, int *csrColumnOffsets_hu_h,
      int *csrRowIndeces_hu_hu, int *csrColumnOffsets_hu_hu,
      int *csrRowIndeces_hu_hv, int *csrColumnOffsets_hu_hv,
      int *csrRowIndeces_hu_heta, int *csrColumnOffsets_hu_heta,
      int *csrRowIndeces_hu_hw, int *csrColumnOffsets_hu_hw,
      // hv
      int *csrRowIndeces_hv_h, int *csrColumnOffsets_hv_h,
      int *csrRowIndeces_hv_hu, int *csrColumnOffsets_hv_hu,
      int *csrRowIndeces_hv_hv, int *csrColumnOffsets_hv_hv,
      int *csrRowIndeces_hv_heta, int *csrColumnOffsets_hv_heta,
      int *csrRowIndeces_hv_hw, int *csrColumnOffsets_hv_hw,
      // heta
      int *csrRowIndeces_heta_h, int *csrColumnOffsets_heta_h,
      int *csrRowIndeces_heta_hu, int *csrColumnOffsets_heta_hu,
      int *csrRowIndeces_heta_hv, int *csrColumnOffsets_heta_hv,
      int *csrRowIndeces_heta_heta, int *csrColumnOffsets_heta_heta,
      int *csrRowIndeces_heta_hw, int *csrColumnOffsets_heta_hw,
      // hw
      int *csrRowIndeces_hw_h, int *csrColumnOffsets_hw_h,
      int *csrRowIndeces_hw_hu, int *csrColumnOffsets_hw_hu,
      int *csrRowIndeces_hw_hv, int *csrColumnOffsets_hw_hv,
      int *csrRowIndeces_hw_heta, int *csrColumnOffsets_hw_heta,
      int *csrRowIndeces_hw_hw, int *csrColumnOffsets_hw_hw,
      //
      double *globalJacobian, int nExteriorElementBoundaries_global,
      int *exteriorElementBoundariesArray, int *elementBoundaryElementsArray,
      int *elementBoundaryLocalElementBoundariesArray, int *isDOFBoundary_h,
      int *isDOFBoundary_hu, int *isDOFBoundary_hv,
      int *isAdvectiveFluxBoundary_h, int *isAdvectiveFluxBoundary_hu,
      int *isAdvectiveFluxBoundary_hv, int *isDiffusiveFluxBoundary_hu,
      int *isDiffusiveFluxBoundary_hv, double *ebqe_bc_h_ext,
      double *ebqe_bc_flux_mass_ext, double *ebqe_bc_flux_mom_hu_adv_ext,
      double *ebqe_bc_flux_mom_hv_adv_ext, double *ebqe_bc_hu_ext,
      double *ebqe_bc_flux_hu_diff_ext, double *ebqe_penalty_ext,
      double *ebqe_bc_hv_ext, double *ebqe_bc_flux_hv_diff_ext,
      int *csrColumnOffsets_eb_h_h, int *csrColumnOffsets_eb_h_hu,
      int *csrColumnOffsets_eb_h_hv, int *csrColumnOffsets_eb_hu_h,
      int *csrColumnOffsets_eb_hu_hu, int *csrColumnOffsets_eb_hu_hv,
      int *csrColumnOffsets_eb_hv_h, int *csrColumnOffsets_eb_hv_hu,
      int *csrColumnOffsets_eb_hv_hv, double dt) {
    //
    // loop over elements to compute volume integrals and load them into the
    // element Jacobians and global Jacobian
    //
    for (int eN = 0; eN < nElements_global; eN++) {
      register double elementJacobian_h_h[nDOF_test_element]
                                         [nDOF_trial_element],
          elementJacobian_hu_hu[nDOF_test_element][nDOF_trial_element],
          elementJacobian_hv_hv[nDOF_test_element][nDOF_trial_element],
          elementJacobian_heta_heta[nDOF_test_element][nDOF_trial_element],
          elementJacobian_hw_hw[nDOF_test_element][nDOF_trial_element];
      for (int i = 0; i < nDOF_test_element; i++)
        for (int j = 0; j < nDOF_trial_element; j++) {
          elementJacobian_h_h[i][j] = 0.0;
          elementJacobian_hu_hu[i][j] = 0.0;
          elementJacobian_hv_hv[i][j] = 0.0;
          elementJacobian_heta_heta[i][j] = 0.0;
          elementJacobian_hw_hw[i][j] = 0.0;
        }
      for (int k = 0; k < nQuadraturePoints_element; k++) {
        int eN_k = eN * nQuadraturePoints_element +
                   k, // index to a scalar at a quadrature point
            eN_k_nSpace = eN_k * nSpace,
            eN_nDOF_trial_element =
                eN * nDOF_trial_element; // index to a vector at a
                                         // quadrature point

        // declare local storage
        register double jac[nSpace * nSpace], jacDet, jacInv[nSpace * nSpace],
            dV, h_test_dV[nDOF_test_element], vel_test_dV[nDOF_test_element], x,
            y, xt, yt;
        // get jacobian, etc for mapping reference element
        ck.calculateMapping_element(eN, k, mesh_dof, mesh_l2g, mesh_trial_ref,
                                    mesh_grad_trial_ref, jac, jacDet, jacInv, x,
                                    y);
        // get the physical integration weight
        dV = fabs(jacDet) * dV_ref[k];
        // precalculate test function products with integration weights
        for (int j = 0; j < nDOF_trial_element; j++) {
          h_test_dV[j] = h_test_ref[k * nDOF_trial_element + j] * dV;
          vel_test_dV[j] = vel_test_ref[k * nDOF_trial_element + j] * dV;
        }
        for (int i = 0; i < nDOF_test_element; i++) {
          register int i_nSpace = i * nSpace;
          for (int j = 0; j < nDOF_trial_element; j++) {
            register int j_nSpace = j * nSpace;
            elementJacobian_h_h[i][j] +=
                h_trial_ref[k * nDOF_trial_element + j] * h_test_dV[i];
            elementJacobian_hu_hu[i][j] +=
                vel_trial_ref[k * nDOF_trial_element + j] * vel_test_dV[i];
            elementJacobian_hv_hv[i][j] +=
                vel_trial_ref[k * nDOF_trial_element + j] * vel_test_dV[i];
            elementJacobian_heta_heta[i][j] +=
                vel_trial_ref[k * nDOF_trial_element + j] * vel_test_dV[i];
            elementJacobian_hw_hw[i][j] +=
                vel_trial_ref[k * nDOF_trial_element + j] * vel_test_dV[i];
          } // j
        }   // i
      }     // k
      //
      // load into element Jacobian into global Jacobian
      //
      for (int i = 0; i < nDOF_test_element; i++) {
        register int eN_i = eN * nDOF_test_element + i;
        for (int j = 0; j < nDOF_trial_element; j++) {
          register int eN_i_j = eN_i * nDOF_trial_element + j;
          globalJacobian[csrRowIndeces_h_h[eN_i] +
                         csrColumnOffsets_h_h[eN_i_j]] +=
              elementJacobian_h_h[i][j];
          globalJacobian[csrRowIndeces_hu_hu[eN_i] +
                         csrColumnOffsets_hu_hu[eN_i_j]] +=
              elementJacobian_hu_hu[i][j];
          globalJacobian[csrRowIndeces_hv_hv[eN_i] +
                         csrColumnOffsets_hv_hv[eN_i_j]] +=
              elementJacobian_hv_hv[i][j];
          globalJacobian[csrRowIndeces_heta_heta[eN_i] +
                         csrColumnOffsets_heta_heta[eN_i_j]] +=
              elementJacobian_heta_heta[i][j];
          globalJacobian[csrRowIndeces_hw_hw[eN_i] +
                         csrColumnOffsets_hw_hw[eN_i_j]] +=
              elementJacobian_hw_hw[i][j];
        } // j
      }   // i
    }     // elements
  }

  void calculateLumpedMassMatrix( // element
      double *mesh_trial_ref, double *mesh_grad_trial_ref, double *mesh_dof,
      double *mesh_velocity_dof, double MOVING_DOMAIN, int *mesh_l2g,
      double *dV_ref, double *h_trial_ref, double *h_grad_trial_ref,
      double *h_test_ref, double *h_grad_test_ref, double *vel_trial_ref,
      double *vel_grad_trial_ref, double *vel_test_ref,
      double *vel_grad_test_ref,
      // element boundary
      double *mesh_trial_trace_ref, double *mesh_grad_trial_trace_ref,
      double *dS_ref, double *h_trial_trace_ref, double *h_grad_trial_trace_ref,
      double *h_test_trace_ref, double *h_grad_test_trace_ref,
      double *vel_trial_trace_ref, double *vel_grad_trial_trace_ref,
      double *vel_test_trace_ref, double *vel_grad_test_trace_ref,
      double *normal_ref, double *boundaryJac_ref,
      // physics
      double *elementDiameter, int nElements_global, double useRBLES,
      double useMetrics, double alphaBDF, double nu, double g, int *h_l2g,
      int *vel_l2g, double *b_dof, double *h_dof, double *hu_dof,
      double *hv_dof, double *h_dof_sge, double *hu_dof_sge, double *hv_dof_sge,
      double *q_mass_acc_beta_bdf, double *q_mom_hu_acc_beta_bdf,
      double *q_mom_hv_acc_beta_bdf, double *q_cfl, int *sdInfo_hu_hu_rowptr,
      int *sdInfo_hu_hu_colind, int *sdInfo_hu_hv_rowptr,
      int *sdInfo_hu_hv_colind, int *sdInfo_hv_hv_rowptr,
      int *sdInfo_hv_hv_colind, int *sdInfo_hv_hu_rowptr,
      int *sdInfo_hv_hu_colind,
      // h
      int *csrRowIndeces_h_h, int *csrColumnOffsets_h_h,
      int *csrRowIndeces_h_hu, int *csrColumnOffsets_h_hu,
      int *csrRowIndeces_h_hv, int *csrColumnOffsets_h_hv,
      int *csrRowIndeces_h_heta, int *csrColumnOffsets_h_heta,
      int *csrRowIndeces_h_hw, int *csrColumnOffsets_h_hw,
      // hu
      int *csrRowIndeces_hu_h, int *csrColumnOffsets_hu_h,
      int *csrRowIndeces_hu_hu, int *csrColumnOffsets_hu_hu,
      int *csrRowIndeces_hu_hv, int *csrColumnOffsets_hu_hv,
      int *csrRowIndeces_hu_heta, int *csrColumnOffsets_hu_heta,
      int *csrRowIndeces_hu_hw, int *csrColumnOffsets_hu_hw,
      // hv
      int *csrRowIndeces_hv_h, int *csrColumnOffsets_hv_h,
      int *csrRowIndeces_hv_hu, int *csrColumnOffsets_hv_hu,
      int *csrRowIndeces_hv_hv, int *csrColumnOffsets_hv_hv,
      int *csrRowIndeces_hv_heta, int *csrColumnOffsets_hv_heta,
      int *csrRowIndeces_hv_hw, int *csrColumnOffsets_hv_hw,
      // heta
      int *csrRowIndeces_heta_h, int *csrColumnOffsets_heta_h,
      int *csrRowIndeces_heta_hu, int *csrColumnOffsets_heta_hu,
      int *csrRowIndeces_heta_hv, int *csrColumnOffsets_heta_hv,
      int *csrRowIndeces_heta_heta, int *csrColumnOffsets_heta_heta,
      int *csrRowIndeces_heta_hw, int *csrColumnOffsets_heta_hw,
      // hw
      int *csrRowIndeces_hw_h, int *csrColumnOffsets_hw_h,
      int *csrRowIndeces_hw_hu, int *csrColumnOffsets_hw_hu,
      int *csrRowIndeces_hw_hv, int *csrColumnOffsets_hw_hv,
      int *csrRowIndeces_hw_heta, int *csrColumnOffsets_hw_heta,
      int *csrRowIndeces_hw_hw, int *csrColumnOffsets_hw_hw,
      //
      double *globalJacobian, int nExteriorElementBoundaries_global,
      int *exteriorElementBoundariesArray, int *elementBoundaryElementsArray,
      int *elementBoundaryLocalElementBoundariesArray, int *isDOFBoundary_h,
      int *isDOFBoundary_hu, int *isDOFBoundary_hv,
      int *isAdvectiveFluxBoundary_h, int *isAdvectiveFluxBoundary_hu,
      int *isAdvectiveFluxBoundary_hv, int *isDiffusiveFluxBoundary_hu,
      int *isDiffusiveFluxBoundary_hv, double *ebqe_bc_h_ext,
      double *ebqe_bc_flux_mass_ext, double *ebqe_bc_flux_mom_hu_adv_ext,
      double *ebqe_bc_flux_mom_hv_adv_ext, double *ebqe_bc_hu_ext,
      double *ebqe_bc_flux_hu_diff_ext, double *ebqe_penalty_ext,
      double *ebqe_bc_hv_ext, double *ebqe_bc_flux_hv_diff_ext,
      int *csrColumnOffsets_eb_h_h, int *csrColumnOffsets_eb_h_hu,
      int *csrColumnOffsets_eb_h_hv, int *csrColumnOffsets_eb_hu_h,
      int *csrColumnOffsets_eb_hu_hu, int *csrColumnOffsets_eb_hu_hv,
      int *csrColumnOffsets_eb_hv_h, int *csrColumnOffsets_eb_hv_hu,
      int *csrColumnOffsets_eb_hv_hv, double dt) {
    //
    // loop over elements to compute volume integrals and load them into the
    // element Jacobians and global Jacobian
    //
    for (int eN = 0; eN < nElements_global; eN++) {
      register double elementJacobian_h_h[nDOF_test_element]
                                         [nDOF_trial_element],
          elementJacobian_hu_hu[nDOF_test_element][nDOF_trial_element],
          elementJacobian_hv_hv[nDOF_test_element][nDOF_trial_element],
          elementJacobian_heta_heta[nDOF_test_element][nDOF_trial_element],
          elementJacobian_hw_hw[nDOF_test_element][nDOF_trial_element];
      for (int i = 0; i < nDOF_test_element; i++)
        for (int j = 0; j < nDOF_trial_element; j++) {
          elementJacobian_h_h[i][j] = 0.0;
          elementJacobian_hu_hu[i][j] = 0.0;
          elementJacobian_hv_hv[i][j] = 0.0;
          elementJacobian_heta_heta[i][j] = 0.0;
          elementJacobian_heta_heta[i][j] = 0.0;
        }
      for (int k = 0; k < nQuadraturePoints_element; k++) {
        int eN_k = eN * nQuadraturePoints_element +
                   k, // index to a scalar at a quadrature point
            eN_k_nSpace = eN_k * nSpace,
            eN_nDOF_trial_element =
                eN * nDOF_trial_element; // index to a vector at a
                                         // quadrature point

        // declare local storage
        register double jac[nSpace * nSpace], jacDet, jacInv[nSpace * nSpace],
            dV, h_test_dV[nDOF_test_element], vel_test_dV[nDOF_test_element], x,
            y, xt, yt;
        // get jacobian, etc for mapping reference element
        ck.calculateMapping_element(eN, k, mesh_dof, mesh_l2g, mesh_trial_ref,
                                    mesh_grad_trial_ref, jac, jacDet, jacInv, x,
                                    y);
        // get the physical integration weight
        dV = fabs(jacDet) * dV_ref[k];
        // precalculate test function products with integration weights
        for (int j = 0; j < nDOF_trial_element; j++) {
          h_test_dV[j] = h_test_ref[k * nDOF_trial_element + j] * dV;
          vel_test_dV[j] = vel_test_ref[k * nDOF_trial_element + j] * dV;
        }

        for (int i = 0; i < nDOF_test_element; i++) {
          register int i_nSpace = i * nSpace;
          for (int j = 0; j < nDOF_trial_element; j++) {
            register int j_nSpace = j * nSpace;
            elementJacobian_h_h[i][j] += (i == j ? 1.0 : 0.0) * h_test_dV[i];
            elementJacobian_hu_hu[i][j] +=
                (i == j ? 1.0 : 0.0) * vel_test_dV[i];
            elementJacobian_hv_hv[i][j] +=
                (i == j ? 1.0 : 0.0) * vel_test_dV[i];
            elementJacobian_heta_heta[i][j] +=
                (i == j ? 1.0 : 0.0) * vel_test_dV[i];
            elementJacobian_hw_hw[i][j] +=
                (i == j ? 1.0 : 0.0) * vel_test_dV[i];
          } // j
        }   // i
      }     // k
      //
      // load into element Jacobian into global Jacobian
      //
      for (int i = 0; i < nDOF_test_element; i++) {
        register int eN_i = eN * nDOF_test_element + i;
        for (int j = 0; j < nDOF_trial_element; j++) {
          register int eN_i_j = eN_i * nDOF_trial_element + j;
          globalJacobian[csrRowIndeces_h_h[eN_i] +
                         csrColumnOffsets_h_h[eN_i_j]] +=
              elementJacobian_h_h[i][j];
          globalJacobian[csrRowIndeces_hu_hu[eN_i] +
                         csrColumnOffsets_hu_hu[eN_i_j]] +=
              elementJacobian_hu_hu[i][j];
          globalJacobian[csrRowIndeces_hv_hv[eN_i] +
                         csrColumnOffsets_hv_hv[eN_i_j]] +=
              elementJacobian_hv_hv[i][j];
          globalJacobian[csrRowIndeces_heta_heta[eN_i] +
                         csrColumnOffsets_heta_heta[eN_i_j]] +=
              elementJacobian_heta_heta[i][j];
          globalJacobian[csrRowIndeces_hw_hw[eN_i] +
                         csrColumnOffsets_hw_hw[eN_i_j]] +=
              elementJacobian_hw_hw[i][j];
        } // j
      }   // i
    }     // elements
  }
}; // GN_SW2DCV

inline GN_SW2DCV_base *
newGN_SW2DCV(int nSpaceIn, int nQuadraturePoints_elementIn,
             int nDOF_mesh_trial_elementIn, int nDOF_trial_elementIn,
             int nDOF_test_elementIn, int nQuadraturePoints_elementBoundaryIn,
             int CompKernelFlag) {
  return proteus::chooseAndAllocateDiscretization2D<GN_SW2DCV_base, GN_SW2DCV,
                                                    CompKernel>(
      nSpaceIn, nQuadraturePoints_elementIn, nDOF_mesh_trial_elementIn,
      nDOF_trial_elementIn, nDOF_test_elementIn,
      nQuadraturePoints_elementBoundaryIn, CompKernelFlag);
}
} // namespace proteus

#endif

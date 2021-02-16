#ifndef GN_SW2DCV_H
#define GN_SW2DCV_H
#include "ArgumentsDict.h"
#include "CompKernel.h"
#include "ModelFactory.h"
#include "xtensor-python/pyarray.hpp"
#include <assert.h>
#include <cmath>
#include <iostream>
#include <valarray>

namespace py = pybind11;

#define POWER_SMOOTHNESS_INDICATOR 2
#define VEL_FIX_POWER 2.
#define LAMBDA_MGN 1
#define IF_BOTH_GAMMA_BRANCHES 0
#define LIMITING_ITERATION 2

/* inline functions */
namespace proteus {

// for Serre-Green-Naghdi "max" wave speeds. See section 4.4 of  JCP 2019
// Guermond, Popov, Tovar, Kees for formulas
inline double GN_nu1(const double &g, const double &hL, const double &uL,
                     const double &etaL, const double &meshSizeL) {

  double augL = LAMBDA_MGN / (3. * meshSizeL) * (6. * hL + 12. * (hL - etaL));

  if (etaL >= hL) {
    augL = LAMBDA_MGN / (3. * meshSizeL) * (6. * hL);
  }
  augL = augL * std::pow(meshSizeL / fmax(meshSizeL, hL), 2);

  return uL - sqrt(g * hL) * sqrt(1. + augL);
}
inline double GN_nu3(const double &g, const double &hR, const double &uR,
                     const double &etaR, const double &meshSizeR) {
  double augR = LAMBDA_MGN / (3. * meshSizeR) * (6. * hR + 12. * (hR - etaR));

  if (etaR >= hR) {
    augR = LAMBDA_MGN / (3. * meshSizeR) * (6. * hR);
  }
  augR = augR * std::pow(meshSizeR / fmax(meshSizeR, hR), 2);
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
} // namespace proteus

namespace proteus {

class GN_SW2DCV_base {
public:
  virtual ~GN_SW2DCV_base() {}
  virtual void convexLimiting(arguments_dict &args) = 0;
  virtual double calculateEdgeBasedCFL(arguments_dict &args) = 0;
  virtual void calculateEV(arguments_dict &args) = 0;
  virtual void calculateResidual(arguments_dict &args) = 0;
  virtual void calculateMassMatrix(arguments_dict &args) = 0;
  virtual void calculateLumpedMassMatrix(arguments_dict &args) = 0;
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
                                              double hEps) {
    double lambda1, lambda3;
    double one_over_hL = 2.0 * hL / (hL * hL + std::pow(fmax(hL, hEps), 2.0));
    double one_over_hR = 2.0 * hR / (hR * hR + std::pow(fmax(hR, hEps), 2.0));
    double hVelL = nx * huL + ny * hvL;
    double hVelR = nx * huR + ny * hvR;
    double velL = one_over_hL * hVelL;
    double velR = one_over_hR * hVelR;
    double etaL = 2.0 * hL / (hL * hL + std::pow(fmax(hL, hEps), 2)) * hetaL;
    double etaR = 2.0 * hR / (hR * hR + std::pow(fmax(hR, hEps), 2)) * hetaR;
    double meshSizeL = sqrt(lumpedL);
    double meshSizeR = sqrt(lumpedR);

    /* See equation 4.12 JCP 2019 paper:
      1-eigenvalue: uL-sqrt(g*hL)*sqrt(1 + augL)
      3-eigenvalue: uR+sqrt(g*hR)*sqrt(1 + augR)
    */
    lambda1 = GN_nu1(g, hL, velL, etaL, meshSizeL);
    lambda3 = GN_nu3(g, hR, velR, etaR, meshSizeR);

    return fmax(fabs(lambda1), fabs(lambda3));
  }

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

  void convexLimiting(arguments_dict &args) {
    double dt = args.scalar<double>("dt");
    int NNZ = args.scalar<int>("NNZ");
    int numDOFs = args.scalar<int>("numDOFs");
    xt::pyarray<double> &lumped_mass_matrix =
        args.array<double>("lumped_mass_matrix");
    xt::pyarray<double> &h_old = args.array<double>("h_old");
    xt::pyarray<double> &hu_old = args.array<double>("hu_old");
    xt::pyarray<double> &hv_old = args.array<double>("hv_old");
    xt::pyarray<double> &heta_old = args.array<double>("heta_old");
    xt::pyarray<double> &hw_old = args.array<double>("hw_old");
    xt::pyarray<double> &hbeta_old = args.array<double>("hbeta_old");
    xt::pyarray<double> &b_dof = args.array<double>("b_dof");
    xt::pyarray<double> &high_order_hnp1 =
        args.array<double>("high_order_hnp1");
    xt::pyarray<double> &high_order_hunp1 =
        args.array<double>("high_order_hunp1");
    xt::pyarray<double> &high_order_hvnp1 =
        args.array<double>("high_order_hvnp1");
    xt::pyarray<double> &high_order_hetanp1 =
        args.array<double>("high_order_hetanp1");
    xt::pyarray<double> &high_order_hwnp1 =
        args.array<double>("high_order_hwnp1");
    xt::pyarray<double> &high_order_hbetanp1 =
        args.array<double>("high_order_hbetanp1");
    xt::pyarray<double> &extendedSourceTerm_hu =
        args.array<double>("extendedSourceTerm_hu");
    xt::pyarray<double> &extendedSourceTerm_hv =
        args.array<double>("extendedSourceTerm_hv");
    xt::pyarray<double> &extendedSourceTerm_heta =
        args.array<double>("extendedSourceTerm_heta");
    xt::pyarray<double> &extendedSourceTerm_hw =
        args.array<double>("extendedSourceTerm_hw");
    xt::pyarray<double> &extendedSourceTerm_hbeta =
        args.array<double>("extendedSourceTerm_hbeta");
    xt::pyarray<double> &limited_hnp1 = args.array<double>("limited_hnp1");
    xt::pyarray<double> &limited_hunp1 = args.array<double>("limited_hunp1");
    xt::pyarray<double> &limited_hvnp1 = args.array<double>("limited_hvnp1");
    xt::pyarray<double> &limited_hetanp1 =
        args.array<double>("limited_hetanp1");
    xt::pyarray<double> &limited_hwnp1 = args.array<double>("limited_hwnp1");
    xt::pyarray<double> &limited_hbetanp1 =
        args.array<double>("limited_hbetanp1");
    xt::pyarray<int> &csrRowIndeces_DofLoops =
        args.array<int>("csrRowIndeces_DofLoops");
    xt::pyarray<int> &csrColumnOffsets_DofLoops =
        args.array<int>("csrColumnOffsets_DofLoops");
    xt::pyarray<double> &MassMatrix = args.array<double>("MassMatrix");
    xt::pyarray<double> &dH_minus_dL = args.array<double>("dH_minus_dL");
    xt::pyarray<double> &muH_minus_muL = args.array<double>("muH_minus_muL");
    double hEps = args.scalar<double>("hEps");
    int LUMPED_MASS_MATRIX = args.scalar<int>("LUMPED_MASS_MATRIX");
    xt::pyarray<double> &dLow = args.array<double>("dLow");
    xt::pyarray<double> &new_SourceTerm_hu =
        args.array<double>("new_SourceTerm_hu");
    xt::pyarray<double> &new_SourceTerm_hv =
        args.array<double>("new_SourceTerm_hv");
    xt::pyarray<double> &new_SourceTerm_heta =
        args.array<double>("new_SourceTerm_heta");
    xt::pyarray<double> &new_SourceTerm_hw =
        args.array<double>("new_SourceTerm_hw");
    xt::pyarray<double> &hLow = args.array<double>("hLow");
    xt::pyarray<double> &huLow = args.array<double>("huLow");
    xt::pyarray<double> &hvLow = args.array<double>("hvLow");
    xt::pyarray<double> &hetaLow = args.array<double>("hetaLow");
    xt::pyarray<double> &hwLow = args.array<double>("hwLow");
    xt::pyarray<double> &hbetaLow = args.array<double>("hbetaLow");
    xt::pyarray<double> &h_min = args.array<double>("h_min");
    xt::pyarray<double> &h_max = args.array<double>("h_max");
    xt::pyarray<double> &heta_min = args.array<double>("heta_min");
    xt::pyarray<double> &heta_max = args.array<double>("heta_max");
    xt::pyarray<double> &kin_max = args.array<double>("kin_max");

    // Declare stuff for limiting on h and h*heta
    std::valarray<double> Rneg(0.0, numDOFs), Rpos(0.0, numDOFs),
        Rneg_heta(0.0, numDOFs), Rpos_heta(0.0, numDOFs);

    // Create FCT component matrices in vector form
    std::valarray<double> FCT_h(0.0, dH_minus_dL.size()),
        FCT_hu(0.0, dH_minus_dL.size()), FCT_hv(0.0, dH_minus_dL.size()),
        FCT_heta(0.0, dH_minus_dL.size()), FCT_hw(0.0, dH_minus_dL.size()),
        FCT_hbeta(0.0, dH_minus_dL.size());

    ////////////////////////////////////////////////////
    // Loop to define FCT matrices for each component //
    ////////////////////////////////////////////////////
    int ij = 0;
    for (int i = 0; i < numDOFs; i++) {
      // Read some vectors
      double high_order_hnp1i = high_order_hnp1[i];
      double high_order_hunp1i = high_order_hunp1[i];
      double high_order_hvnp1i = high_order_hvnp1[i];
      double high_order_hetanp1i = high_order_hetanp1[i];
      double high_order_hwnp1i = high_order_hwnp1[i];
      double high_order_hbetanp1i = high_order_hbetanp1[i];
      double hi = h_old[i];
      double huni = hu_old[i];
      double hvni = hv_old[i];
      double hetani = heta_old[i];
      double hwni = hw_old[i];
      double hbetani = hbeta_old[i];
      double Zi = b_dof[i];
      double mi = lumped_mass_matrix[i];
      double one_over_hiReg =
          2 * hi / (hi * hi + std::pow(fmax(hi, hEps), 2)); // hEps

      // LOOP OVER THE SPARSITY PATTERN (j-LOOP)//
      for (int offset = csrRowIndeces_DofLoops[i];
           offset < csrRowIndeces_DofLoops[i + 1]; offset++) {

        int j = csrColumnOffsets_DofLoops[offset];

        if (i != j) {
          // Read some vectors
          double hj = h_old[j];
          double hunj = hu_old[j];
          double hvnj = hv_old[j];
          double hetanj = heta_old[j];
          double hwnj = hw_old[j];
          double hbetanj = hbeta_old[j];
          double Zj = b_dof[j];
          double one_over_hjReg =
              2. * hj / (hj * hj + std::pow(fmax(hj, hEps), 2));

          // Compute star states
          double hStarij = fmax(0., hi + Zi - fmax(Zi, Zj));
          double huStarij = huni * hStarij * one_over_hiReg;
          double hvStarij = hvni * hStarij * one_over_hiReg;
          double hetaStarij = hetani * std::pow(hStarij * one_over_hiReg, 2);
          double hwStarij = hwni * hStarij * one_over_hiReg;
          double hbetaStarij = hbetani * hStarij * one_over_hiReg;

          double hStarji = fmax(0., hj + Zj - fmax(Zi, Zj));
          double huStarji = hunj * hStarji * one_over_hjReg;
          double hvStarji = hvnj * hStarji * one_over_hjReg;
          double hetaStarji = hetanj * std::pow(hStarji * one_over_hjReg, 2);
          double hwStarji = hwnj * hStarji * one_over_hjReg;
          double hbetaStarji = hbetanj * hStarji * one_over_hjReg;

          // i-th row of flux correction matrix
          double ML_minus_MC = (LUMPED_MASS_MATRIX == 1
                                    ? 0.
                                    : (i == j ? 1. : 0.) * mi - MassMatrix[ij]);

          FCT_h[ij] =
              ML_minus_MC *
                  (high_order_hnp1[j] - hj - (high_order_hnp1i - hi)) +
              dt * (dH_minus_dL[ij] - muH_minus_muL[ij]) * (hStarji - hStarij) +
              dt * muH_minus_muL[ij] * (hj - hi);

          FCT_hu[ij] = ML_minus_MC * (high_order_hunp1[j] - hunj -
                                      (high_order_hunp1i - huni)) +
                       dt * (dH_minus_dL[ij] - muH_minus_muL[ij]) *
                           (huStarji - huStarij) +
                       dt * muH_minus_muL[ij] * (hunj - huni);

          FCT_hv[ij] = ML_minus_MC * (high_order_hvnp1[j] - hvnj -
                                      (high_order_hvnp1i - hvni)) +
                       dt * (dH_minus_dL[ij] - muH_minus_muL[ij]) *
                           (hvStarji - hvStarij) +
                       dt * muH_minus_muL[ij] * (hvnj - hvni);

          FCT_heta[ij] = ML_minus_MC * (high_order_hetanp1[j] - hetanj -
                                        (high_order_hetanp1i - hetani)) +
                         dt * (dH_minus_dL[ij] - muH_minus_muL[ij]) *
                             (hetaStarji - hetaStarij) +
                         dt * muH_minus_muL[ij] * (hetanj - hetani);

          FCT_hw[ij] = ML_minus_MC * (high_order_hwnp1[j] - hwnj -
                                      (high_order_hwnp1i - hwni)) +
                       dt * (dH_minus_dL[ij] - muH_minus_muL[ij]) *
                           (hwStarji - hwStarij) +
                       dt * muH_minus_muL[ij] * (hwnj - hwni);

          FCT_hbeta[ij] = ML_minus_MC * (high_order_hbetanp1[j] - hbetanj -
                                         (high_order_hbetanp1i - hbetani)) +
                          dt * (dH_minus_dL[ij] - muH_minus_muL[ij]) *
                              (hbetaStarji - hbetaStarij) +
                          dt * muH_minus_muL[ij] * (hbetanj - hbetani);
        } else {
          FCT_h[ij] = 0.0;
          FCT_hu[ij] = 0.0;
          FCT_hv[ij] = 0.0;
          FCT_heta[ij] = 0.0;
          FCT_hw[ij] = 0.0;
          FCT_hbeta[ij] = 0.0;
        }

        // UPDATE ij //
        ij += 1;
      } // j loop ends here
    }   // i loop ends here

    ////////////////////////////////////////////////////////////////////
    // Main loop to define limiters and computed limited solution //////
    ////////////////////////////////////////////////////////////////////

    // Create Lij_array for initialization
    std::valarray<double> Lij_array(1.0, dH_minus_dL.size());

    /* Loop over limiting iterations */
    for (int limit_iter = 0; limit_iter < LIMITING_ITERATION; limit_iter++) {

      /* Loop to define FCT Rpos and Rneg values */
      ij = 0;
      for (int i = 0; i < numDOFs; i++) {

        double hi = h_old[i];
        double mi = lumped_mass_matrix[i];
        // Initialize Pneg and Ppos quantities at ith node
        double Pnegi = 0., Pposi = 0.;
        double Pnegi_heta = 0., Pposi_heta = 0.;

        // LOOP OVER THE SPARSITY PATTERN (j-LOOP)
        for (int offset = csrRowIndeces_DofLoops[i];
             offset < csrRowIndeces_DofLoops[i + 1]; offset++) {

          int j = csrColumnOffsets_DofLoops[offset];

          // COMPUTE P VECTORS
          Pnegi += FCT_h[ij] * ((FCT_h[ij] < 0.) ? 1. : 0.);
          Pposi += FCT_h[ij] * ((FCT_h[ij] > 0.) ? 1. : 0.);
          Pnegi_heta += FCT_heta[ij] * ((FCT_heta[ij] < 0.) ? 1. : 0.);
          Pposi_heta += FCT_heta[ij] * ((FCT_heta[ij] > 0.) ? 1. : 0.);

          // UPDATE ij
          ij += 1;
        } // j loop ends here

        double psmall_h = 1E-14 * fmax(fabs(Pnegi), fabs(Pposi));
        double psmall_heta = 1E-14 * fmax(fabs(Pnegi_heta), fabs(Pposi_heta));

        ///////////////////////
        // COMPUTE Q VECTORS //
        ///////////////////////
        double Qnegi = std::min(mi * (h_min[i] - hLow[i]), 0.0);
        double Qposi = std::max(mi * (h_max[i] - hLow[i]), 0.0);
        double Qnegi_heta = std::min(mi * (heta_min[i] - hetaLow[i]), 0.0);
        double Qposi_heta = std::max(mi * (heta_max[i] - hetaLow[i]), 0.0);

        ///////////////////////
        // COMPUTE R VECTORS //
        ///////////////////////
        if (hi <= hEps) {
          Rneg[i] = 0.;
          Rpos[i] = 0.;
          Rneg_heta[i] = 0.;
          Rpos_heta[i] = 0.;
        } else {
          // for h
          if (Pnegi >= -psmall_h) {
            Rneg[i] = 1.0;
          } else {
            Rneg[i] = std::min(1.0, Qnegi / Pnegi);
          }
          if (Pposi <= psmall_h) {
            Rpos[i] = 1.0;
          } else {
            Rpos[i] = std::min(1.0, Qposi / Pposi);
          }

          // for heta
          if (Pnegi_heta >= -psmall_heta) {
            Rneg_heta[i] = 1.0;
          } else {
            Rneg_heta[i] = std::min(1.0, Qnegi_heta / Pnegi_heta);
          }
          if (Pposi_heta <= psmall_heta) {
            Rpos_heta[i] = 1.0;
          } else {
            Rpos_heta[i] = std::min(1.0, Qposi_heta / Pposi_heta);
          }
        }
      } // i loop ends here

      // for limiting kinetic energy
      auto KE_tiny = xt::amax(kin_max)[0] * hEps;

      /* Here we compute the limiters Lij_array */
      ij = 0;
      for (int i = 0; i < numDOFs; i++) {

        double mi = lumped_mass_matrix[i];
        double ci =
            kin_max[i] * hLow[i] -
            0.5 * (huLow[i] * huLow[i] + hvLow[i] * hvLow[i]); // for KE lim.

        // LOOP OVER THE SPARSITY PATTERN (j-LOOP)//
        for (int offset = csrRowIndeces_DofLoops[i];
             offset < csrRowIndeces_DofLoops[i + 1]; offset++) {

          int j = csrColumnOffsets_DofLoops[offset];

          if (j != i) {

            // Compute limiter based on water height
            if (FCT_h[ij] >= 0.) {
              Lij_array[ij] = fmin(Lij_array[ij], std::min(Rneg[j], Rpos[i]));
            } else {
              Lij_array[ij] = fmin(Lij_array[ij], std::min(Rneg[i], Rpos[j]));
            }

            // COMPUTE LIMITER based on heta
            if (FCT_heta[ij] >= 0.) {
              Lij_array[ij] =
                  fmin(Lij_array[ij], std::min(Rneg_heta[j], Rpos_heta[i]));
            } else {
              Lij_array[ij] =
                  fmin(Lij_array[ij], std::min(Rneg_heta[i], Rpos_heta[j]));
            }

            /*======================================================*/
            /*            KE limiting                               */
            double lambdaj =
                csrRowIndeces_DofLoops[i + 1] - csrRowIndeces_DofLoops[i] - 1;
            double Ph_ij = FCT_h[ij] / mi / lambdaj;
            double Phu_ij = FCT_hu[ij] / mi / lambdaj;
            double Phv_ij = FCT_hv[ij] / mi / lambdaj;

            // Here we have some conditions for the KE limiting

            double KE_limiter = 1;

            // We first check if local kinetic energy > 0
            if (kin_max[i] * hLow[i] <= KE_tiny) {
              KE_limiter = 0.;
            }

            // We then check if KE bound is already satisfied
            double hi_with_lijPij = hLow[i] + Lij_array[ij] * Ph_ij;
            double hui_with_lijPij = huLow[i] + Lij_array[ij] * Phu_ij;
            double hvi_with_lijPij = hvLow[i] + Lij_array[ij] * Phv_ij;
            double psi = kin_max[i] * hi_with_lijPij -
                         0.5 * (hui_with_lijPij * hui_with_lijPij +
                                hvi_with_lijPij * hvi_with_lijPij);
            if (psi > -KE_tiny) {
              KE_limiter = fmin(KE_limiter, Lij_array[ij]);
            }

            /*======================================================*/

            double ai = -0.5 * (Phu_ij * Phu_ij + Phv_ij * Phv_ij);
            double bi =
                kin_max[i] * Ph_ij - (huLow[i] * Phu_ij + hvLow[i] * Phv_ij);
            double delta_i = bi * bi - 4. * ai * ci;
            double ri = 1.;

            if (delta_i < 0. || ai >= -0.) {
              KE_limiter = fmin(KE_limiter, Lij_array[ij]);
            } else {
              double neg_root_i = (-bi - std::sqrt(delta_i)) / 2. / ai;
              KE_limiter = fmin(KE_limiter, neg_root_i);
            }

            // root of jth-DOF (To compute transpose component)
            double lambdai =
                csrRowIndeces_DofLoops[j + 1] - csrRowIndeces_DofLoops[j] - 1;
            double mj = lumped_mass_matrix[j];
            double cj = kin_max[j] * hLow[j] -
                        0.5 * (huLow[j] * huLow[j] + hvLow[j] * hvLow[j]);
            double Ph_ji = -FCT_h[ij] / mj / lambdai; // Aij=-Aji
            double Phu_ji = -FCT_hu[ij] / mj / lambdai;
            double Phv_ji = -FCT_hv[ij] / mj / lambdai;
            double aj = -0.5 * (Phu_ji * Phu_ji + Phv_ji * Phv_ji);
            double bj =
                kin_max[j] * Ph_ji - (huLow[j] * Phu_ji + hvLow[j] * Phv_ji);

            double delta_j = bj * bj - 4. * aj * cj;
            double rj = 1.;

            if (delta_j < 0. || aj >= -0.) {
              KE_limiter = fmin(KE_limiter, Lij_array[ij]);
            } else {
              double neg_root_j = (-bj - std::sqrt(delta_j)) / 2. / aj;
              KE_limiter = fmin(KE_limiter, neg_root_j);
            }

            // Compute KE limiter
            Lij_array[ij] = fmin(KE_limiter, Lij_array[ij]);
            // Lij_array[ij] = fmin(fmin(fabs(ri), fabs(rj)), Lij_array[ij]);
          } else {
            // if i = j then lij = 0
            Lij_array[ij] = 0.;
          }

          // update ij
          ij += 1;
        } // end j loop
      }   // end i loop

      /* Final loop to apply limiting and then define limited solution */
      ij = 0;
      for (int i = 0; i < numDOFs; i++) {

        double one_over_mi = 1.0 / lumped_mass_matrix[i];
        double ith_Limiter_times_FluxCorrectionMatrix1 = 0.;
        double ith_Limiter_times_FluxCorrectionMatrix2 = 0.;
        double ith_Limiter_times_FluxCorrectionMatrix3 = 0.;
        double ith_Limiter_times_FluxCorrectionMatrix4 = 0.;
        double ith_Limiter_times_FluxCorrectionMatrix5 = 0.;
        double ith_Limiter_times_FluxCorrectionMatrix6 = 0.;

        // LOOP OVER THE SPARSITY PATTERN (j-LOOP)//
        for (int offset = csrRowIndeces_DofLoops[i];
             offset < csrRowIndeces_DofLoops[i + 1]; offset++) {

          int j = csrColumnOffsets_DofLoops[offset];

          // COMPUTE LIMITED FLUX //
          ith_Limiter_times_FluxCorrectionMatrix1 += Lij_array[ij] * FCT_h[ij];
          ith_Limiter_times_FluxCorrectionMatrix2 += Lij_array[ij] * FCT_hu[ij];
          ith_Limiter_times_FluxCorrectionMatrix3 += Lij_array[ij] * FCT_hv[ij];
          ith_Limiter_times_FluxCorrectionMatrix4 +=
              Lij_array[ij] * FCT_heta[ij];
          ith_Limiter_times_FluxCorrectionMatrix5 += Lij_array[ij] * FCT_hw[ij];
          ith_Limiter_times_FluxCorrectionMatrix6 +=
              Lij_array[ij] * FCT_hbeta[ij];

          // update ij
          ij += 1;
        } // end j loop

        // then we add lij*Aij to uLow
        hLow[i] += one_over_mi * ith_Limiter_times_FluxCorrectionMatrix1;
        huLow[i] += one_over_mi * ith_Limiter_times_FluxCorrectionMatrix2;
        hvLow[i] += one_over_mi * ith_Limiter_times_FluxCorrectionMatrix3;
        hetaLow[i] += one_over_mi * ith_Limiter_times_FluxCorrectionMatrix4;
        hwLow[i] += one_over_mi * ith_Limiter_times_FluxCorrectionMatrix5;
        hbetaLow[i] += one_over_mi * ith_Limiter_times_FluxCorrectionMatrix6;

        // Finally define the limited solution
        limited_hnp1[i] = hLow[i];
        limited_hunp1[i] = huLow[i] + dt * one_over_mi * new_SourceTerm_hu[i];
        limited_hvnp1[i] = hvLow[i] + dt * one_over_mi * new_SourceTerm_hv[i];
        limited_hetanp1[i] =
            hetaLow[i] - dt * one_over_mi * extendedSourceTerm_heta[i];
        limited_hwnp1[i] =
            hwLow[i] - dt * one_over_mi * extendedSourceTerm_hw[i];
        limited_hbetanp1[i] =
            hbetaLow[i] - dt * one_over_mi * extendedSourceTerm_hbeta[i];

        if (limited_hnp1[i] < -hEps && dt < 1.0) {
          std::cout << "Limited water height is negative: \n "
                    << "hLow: " << hLow[i] << "\n"
                    << "hHigh: " << limited_hnp1[i] << "\n"
                    << "hEps:  " << hEps << "\n"
                    << " ... aborting!" << std::endl;
          abort();
        } else {
          // clean up uHigh from round off error
          if (limited_hnp1[i] < hEps) {
            limited_hnp1[i] = 0.0;
            limited_hetanp1[i] = 0.0;
          }
          //
          double aux = fmax(limited_hnp1[i], hEps);
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
          limited_hbetanp1[i] *= 2 * std::pow(limited_hnp1[i], VEL_FIX_POWER) /
                                 (std::pow(limited_hnp1[i], VEL_FIX_POWER) +
                                  std::pow(aux, VEL_FIX_POWER));
        }
      } // end i loop

      // update FCT matrices as Fct = (1 - Lij)*Fct
      FCT_h = (1.0 - Lij_array) * FCT_h;
      FCT_hu = (1.0 - Lij_array) * FCT_hu;
      FCT_hv = (1.0 - Lij_array) * FCT_hv;
      FCT_heta = (1.0 - Lij_array) * FCT_heta;
      FCT_hw = (1.0 - Lij_array) * FCT_hw;
      FCT_hbeta = (1.0 - Lij_array) * FCT_hbeta;
    } // end loop for limiting iteration
  }   // end convex limiting function

  double calculateEdgeBasedCFL(arguments_dict &args) {
    double g = args.scalar<double>("g");
    int numDOFsPerEqn = args.scalar<int>("numDOFsPerEqn");
    xt::pyarray<double> &lumped_mass_matrix =
        args.array<double>("lumped_mass_matrix");
    xt::pyarray<double> &h_dof_old = args.array<double>("h_dof_old");
    xt::pyarray<double> &hu_dof_old = args.array<double>("hu_dof_old");
    xt::pyarray<double> &hv_dof_old = args.array<double>("hv_dof_old");
    xt::pyarray<double> &heta_dof_old = args.array<double>("heta_dof_old");
    xt::pyarray<double> &b_dof = args.array<double>("b_dof");
    xt::pyarray<int> &csrRowIndeces_DofLoops =
        args.array<int>("csrRowIndeces_DofLoops");
    xt::pyarray<int> &csrColumnOffsets_DofLoops =
        args.array<int>("csrColumnOffsets_DofLoops");
    double hEps = args.scalar<double>("hEps");
    xt::pyarray<double> &Cx = args.array<double>("Cx");
    xt::pyarray<double> &Cy = args.array<double>("Cy");
    xt::pyarray<double> &CTx = args.array<double>("CTx");
    xt::pyarray<double> &CTy = args.array<double>("CTy");
    xt::pyarray<double> &dLow = args.array<double>("dLow");
    double run_cfl = args.scalar<double>("run_cfl");
    xt::pyarray<double> &edge_based_cfl = args.array<double>("edge_based_cfl");

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
        int j = csrColumnOffsets_DofLoops[offset];

        // solution at time tn for the jth DOF
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
          dLow[ij] = fmax(
              maxWaveSpeedSharpInitialGuess(g, nxij, nyij, hi, hui, hvi, hetai,
                                            mi, hj, huj, hvj, hetaj, mj, hEps) *
                  cij_norm, // hEps
              maxWaveSpeedSharpInitialGuess(g, nxji, nyji, hj, huj, hvj, hetaj,
                                            mj, hi, hui, hvi, hetai, mi, hEps) *
                  cji_norm); // hEps
          dLowii -= dLow[ij];

        } else
          dLow[ij] = 0.;
        // update ij
        ij += 1;
      }
      //////////////////////////////
      // CALCULATE EDGE BASED CFL //
      //////////////////////////////
      edge_based_cfl[i] = 1.0 * fabs(dLowii) / mi;
      max_edge_based_cfl = fmax(max_edge_based_cfl, edge_based_cfl[i]);
    }
    return max_edge_based_cfl;
  } // End calculateEdgeBasedCFL

  void calculateEV(arguments_dict &args) {
    double g = args.scalar<double>("g");
    xt::pyarray<double> &h_dof_old = args.array<double>("h_dof_old");
    xt::pyarray<double> &hu_dof_old = args.array<double>("hu_dof_old");
    xt::pyarray<double> &hv_dof_old = args.array<double>("hv_dof_old");
    xt::pyarray<double> &b_dof = args.array<double>("b_dof");
    xt::pyarray<double> &Cx = args.array<double>("Cx");
    xt::pyarray<double> &Cy = args.array<double>("Cy");
    xt::pyarray<double> &CTx = args.array<double>("CTx");
    xt::pyarray<double> &CTy = args.array<double>("CTy");
    int numDOFsPerEqn = args.scalar<int>("numDOFsPerEqn");
    xt::pyarray<int> &csrRowIndeces_DofLoops =
        args.array<int>("csrRowIndeces_DofLoops");
    xt::pyarray<int> &csrColumnOffsets_DofLoops =
        args.array<int>("csrColumnOffsets_DofLoops");
    xt::pyarray<double> &lumped_mass_matrix =
        args.array<double>("lumped_mass_matrix");
    double eps = args.scalar<double>("eps");
    double hEps = args.scalar<double>("hEps");
    xt::pyarray<double> &global_entropy_residual =
        args.array<double>("global_entropy_residual");
    double &dij_small = args.scalar<double>("dij_small");

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
      eta[i] = ENTROPY(g, hi, hu_dof_old[i], hv_dof_old[i], 0., one_over_hiReg);
    }

    // ********** END OF FIRST LOOP ON DOFs ********** //

    ///////////////////////////////////////////////
    // ********** SECOND LOOP ON DOFs ********** //
    ///////////////////////////////////////////////
    // To compute:
    //     * global entropy residual
    //     * dij_small to avoid division by 0

    int ij = 0;
    std::valarray<double> etaMax(numDOFsPerEqn), etaMin(numDOFsPerEqn);

    // speed = sqrt(g max(h_0)), I divide by h_epsilon to get max(h_0)
    double speed = std::sqrt(g * hEps / eps);
    dij_small = 0.0;

    for (int i = 0; i < numDOFsPerEqn; i++) {

      // solution at time tn for the ith DOF
      double hi = h_dof_old[i];
      double hui = hu_dof_old[i];
      double hvi = hv_dof_old[i];
      double Zi = b_dof[i];

      // Define some things using above
      double one_over_hiReg =
          2 * hi / (hi * hi + std::pow(fmax(hi, hEps), 2)); // hEps
      double ui = hui * one_over_hiReg;
      double vi = hvi * one_over_hiReg;
      double mi = lumped_mass_matrix[i];

      // initialize etaMax and etaMin
      etaMax[i] = fabs(eta[i]);
      etaMin[i] = fabs(eta[i]);

      // FOR ENTROPY RESIDUAL, NOTE: FLAT BOTTOM //
      double ith_flux_term1 = 0., ith_flux_term2 = 0., ith_flux_term3 = 0.;
      double entropy_flux = 0.;
      double sum_entprime_flux = 0.;
      double eta_prime1 = DENTROPY_DH(g, hi, hui, hvi, 0., one_over_hiReg);
      double eta_prime2 = DENTROPY_DHU(g, hi, hui, hvi, 0., one_over_hiReg);
      double eta_prime3 = DENTROPY_DHV(g, hi, hui, hvi, 0., one_over_hiReg);

      // loop in j (sparsity pattern)
      for (int offset = csrRowIndeces_DofLoops[i];
           offset < csrRowIndeces_DofLoops[i + 1]; offset++) {

        int j = csrColumnOffsets_DofLoops[offset];

        // solution at time tn for the jth DOF
        double hj = h_dof_old[j];
        double huj = hu_dof_old[j];
        double hvj = hv_dof_old[j];
        double Zj = b_dof[j];

        // Then define some things here using above
        double one_over_hjReg =
            2.0 * hj / (hj * hj + std::pow(fmax(hj, hEps), 2));
        double uj = huj * one_over_hjReg;
        double vj = hvj * one_over_hjReg;

        // auxiliary functions to compute fluxes
        double aux_h =
            (uj * hj - ui * hi) * Cx[ij] + (vj * hj - vi * hi) * Cy[ij];
        double aux_hu =
            (uj * huj - ui * hui) * Cx[ij] + (vj * huj - vi * hui) * Cy[ij];
        double aux_hv =
            (uj * hvj - ui * hvi) * Cx[ij] + (vj * hvj - vi * hvi) * Cy[ij];

        // flux for entropy
        ith_flux_term1 += aux_h;
        ith_flux_term2 += aux_hu + 0.5 * g * hj * hj * Cx[ij];
        ith_flux_term3 += aux_hv + 0.5 * g * hj * hj * Cy[ij];

        // NOTE: WE CONSIDER FLAT BOTTOM
        entropy_flux +=
            (Cx[ij] * ENTROPY_FLUX1(g, hj, huj, hvj, 0., one_over_hjReg) +
             Cy[ij] * ENTROPY_FLUX2(g, hj, huj, hvj, 0., one_over_hjReg));

        // COMPUTE ETA MIN AND ETA MAX //
        etaMax[i] = fmax(etaMax[i], fabs(eta[j]));
        etaMin[i] = fmin(etaMin[i], fabs(eta[j]));

        // define dij_small in j loop
        double x = fabs(Cx[ij]) + fabs(Cy[ij]);
        dij_small = fmax(dij_small, x * speed);

        // update ij
        ij += 1;
      } // end j loop

      // define sum of entprime*flux
      sum_entprime_flux =
          (ith_flux_term1 * eta_prime1 + ith_flux_term2 * eta_prime2 +
           ith_flux_term3 * eta_prime3);

      // define rescale for normalization
      double small_rescale = g * hEps * hEps / eps;
      double rescale = fmax(fabs(etaMax[i] - etaMin[i]) / 2., small_rescale);

      // COMPUTE ENTROPY RESIDUAL //
      double one_over_entNormFactori = 1.0 / rescale;
      global_entropy_residual[i] =
          one_over_entNormFactori * fabs(entropy_flux - sum_entprime_flux);
      if (hi <= hEps) {
        global_entropy_residual[i] = 1.0;
      }
    } // end i loop

    // Finally dij_small here
    dij_small = 1E-14 * dij_small;
    // ********** END OF LOOP IN DOFs ********** //
  } // end calculateEV

  void calculateResidual(arguments_dict &args) {
    xt::pyarray<double> &mesh_trial_ref = args.array<double>("mesh_trial_ref");
    xt::pyarray<double> &mesh_grad_trial_ref =
        args.array<double>("mesh_grad_trial_ref");
    xt::pyarray<double> &mesh_dof = args.array<double>("mesh_dof");
    xt::pyarray<int> &mesh_l2g = args.array<int>("mesh_l2g");
    xt::pyarray<double> &dV_ref = args.array<double>("dV_ref");
    xt::pyarray<double> &h_trial_ref = args.array<double>("h_trial_ref");
    xt::pyarray<double> &h_grad_trial_ref =
        args.array<double>("h_grad_trial_ref");
    xt::pyarray<double> &h_test_ref = args.array<double>("h_test_ref");
    xt::pyarray<double> &h_grad_test_ref =
        args.array<double>("h_grad_test_ref");
    xt::pyarray<double> &vel_trial_ref = args.array<double>("vel_trial_ref");
    xt::pyarray<double> &vel_grad_trial_ref =
        args.array<double>("vel_grad_trial_ref");
    xt::pyarray<double> &vel_test_ref = args.array<double>("vel_test_ref");
    xt::pyarray<double> &vel_grad_test_ref =
        args.array<double>("vel_grad_test_ref");
    xt::pyarray<double> &mesh_trial_trace_ref =
        args.array<double>("mesh_trial_trace_ref");
    xt::pyarray<double> &mesh_grad_trial_trace_ref =
        args.array<double>("mesh_grad_trial_trace_ref");
    xt::pyarray<double> &h_trial_trace_ref =
        args.array<double>("h_trial_trace_ref");
    xt::pyarray<double> &h_grad_trial_trace_ref =
        args.array<double>("h_grad_trial_trace_ref");
    xt::pyarray<double> &h_test_trace_ref =
        args.array<double>("h_test_trace_ref");
    xt::pyarray<double> &h_grad_test_trace_ref =
        args.array<double>("h_grad_test_trace_ref");
    xt::pyarray<double> &vel_trial_trace_ref =
        args.array<double>("vel_trial_trace_ref");
    xt::pyarray<double> &vel_grad_trial_trace_ref =
        args.array<double>("vel_grad_trial_trace_ref");
    xt::pyarray<double> &vel_test_trace_ref =
        args.array<double>("vel_test_trace_ref");
    xt::pyarray<double> &vel_grad_test_trace_ref =
        args.array<double>("vel_grad_test_trace_ref");
    xt::pyarray<double> &normal_ref = args.array<double>("normal_ref");
    xt::pyarray<double> &boundaryJac_ref =
        args.array<double>("boundaryJac_ref");
    xt::pyarray<double> &elementDiameter =
        args.array<double>("elementDiameter");
    int nElements_global = args.scalar<int>("nElements_global");
    double g = args.scalar<double>("g");
    xt::pyarray<int> &h_l2g = args.array<int>("h_l2g");
    xt::pyarray<int> &vel_l2g = args.array<int>("vel_l2g");
    xt::pyarray<double> &h_dof_old = args.array<double>("h_dof_old");
    xt::pyarray<double> &hu_dof_old = args.array<double>("hu_dof_old");
    xt::pyarray<double> &hv_dof_old = args.array<double>("hv_dof_old");
    xt::pyarray<double> &heta_dof_old = args.array<double>("heta_dof_old");
    xt::pyarray<double> &hw_dof_old = args.array<double>("hw_dof_old");
    xt::pyarray<double> &hbeta_dof_old = args.array<double>("hbeta_dof_old");
    xt::pyarray<double> &b_dof = args.array<double>("b_dof");
    xt::pyarray<double> &h_dof = args.array<double>("h_dof");
    xt::pyarray<double> &hu_dof = args.array<double>("hu_dof");
    xt::pyarray<double> &hv_dof = args.array<double>("hv_dof");
    xt::pyarray<double> &heta_dof = args.array<double>("heta_dof");
    xt::pyarray<double> &hw_dof = args.array<double>("hw_dof");
    xt::pyarray<double> &hbeta_dof = args.array<double>("hbeta_dof");
    xt::pyarray<double> &q_cfl = args.array<double>("q_cfl");
    xt::pyarray<int> &sdInfo_hu_hu_rowptr =
        args.array<int>("sdInfo_hu_hu_rowptr");
    xt::pyarray<int> &sdInfo_hu_hu_colind =
        args.array<int>("sdInfo_hu_hu_colind");
    xt::pyarray<int> &sdInfo_hu_hv_rowptr =
        args.array<int>("sdInfo_hu_hv_rowptr");
    xt::pyarray<int> &sdInfo_hu_hv_colind =
        args.array<int>("sdInfo_hu_hv_colind");
    xt::pyarray<int> &sdInfo_hv_hv_rowptr =
        args.array<int>("sdInfo_hv_hv_rowptr");
    xt::pyarray<int> &sdInfo_hv_hv_colind =
        args.array<int>("sdInfo_hv_hv_colind");
    xt::pyarray<int> &sdInfo_hv_hu_rowptr =
        args.array<int>("sdInfo_hv_hu_rowptr");
    xt::pyarray<int> &sdInfo_hv_hu_colind =
        args.array<int>("sdInfo_hv_hu_colind");
    int offset_h = args.scalar<int>("offset_h");
    int offset_hu = args.scalar<int>("offset_hu");
    int offset_hv = args.scalar<int>("offset_hv");
    int offset_heta = args.scalar<int>("offset_heta");
    int offset_hw = args.scalar<int>("offset_hw");
    int offset_hbeta = args.scalar<int>("offset_hbeta");
    int stride_h = args.scalar<int>("stride_h");
    int stride_hu = args.scalar<int>("stride_hu");
    int stride_hv = args.scalar<int>("stride_hv");
    int stride_heta = args.scalar<int>("stride_heta");
    int stride_hw = args.scalar<int>("stride_hw");
    int stride_hbeta = args.scalar<int>("stride_hbeta");
    xt::pyarray<double> &globalResidual = args.array<double>("globalResidual");
    int nExteriorElementBoundaries_global =
        args.scalar<int>("nExteriorElementBoundaries_global");
    xt::pyarray<int> &exteriorElementBoundariesArray =
        args.array<int>("exteriorElementBoundariesArray");
    xt::pyarray<int> &elementBoundaryElementsArray =
        args.array<int>("elementBoundaryElementsArray");
    xt::pyarray<int> &elementBoundaryLocalElementBoundariesArray =
        args.array<int>("elementBoundaryLocalElementBoundariesArray");
    xt::pyarray<int> &isDOFBoundary_h = args.array<int>("isDOFBoundary_h");
    xt::pyarray<int> &isDOFBoundary_hu = args.array<int>("isDOFBoundary_hu");
    xt::pyarray<int> &isDOFBoundary_hv = args.array<int>("isDOFBoundary_hv");
    xt::pyarray<int> &isAdvectiveFluxBoundary_h =
        args.array<int>("isAdvectiveFluxBoundary_h");
    xt::pyarray<int> &isAdvectiveFluxBoundary_hu =
        args.array<int>("isAdvectiveFluxBoundary_hu");
    xt::pyarray<int> &isAdvectiveFluxBoundary_hv =
        args.array<int>("isAdvectiveFluxBoundary_hv");
    xt::pyarray<int> &isDiffusiveFluxBoundary_hu =
        args.array<int>("isDiffusiveFluxBoundary_hu");
    xt::pyarray<int> &isDiffusiveFluxBoundary_hv =
        args.array<int>("isDiffusiveFluxBoundary_hv");
    xt::pyarray<double> &ebqe_bc_h_ext = args.array<double>("ebqe_bc_h_ext");
    xt::pyarray<double> &ebqe_bc_flux_mass_ext =
        args.array<double>("ebqe_bc_flux_mass_ext");
    xt::pyarray<double> &ebqe_bc_flux_mom_hu_adv_ext =
        args.array<double>("ebqe_bc_flux_mom_hu_adv_ext");
    xt::pyarray<double> &ebqe_bc_flux_mom_hv_adv_ext =
        args.array<double>("ebqe_bc_flux_mom_hv_adv_ext");
    xt::pyarray<double> &ebqe_bc_hu_ext = args.array<double>("ebqe_bc_hu_ext");
    xt::pyarray<double> &ebqe_bc_flux_hu_diff_ext =
        args.array<double>("ebqe_bc_flux_hu_diff_ext");
    xt::pyarray<double> &ebqe_penalty_ext =
        args.array<double>("ebqe_penalty_ext");
    xt::pyarray<double> &ebqe_bc_hv_ext = args.array<double>("ebqe_bc_hv_ext");
    xt::pyarray<double> &ebqe_bc_flux_hv_diff_ext =
        args.array<double>("ebqe_bc_flux_hv_diff_ext");
    xt::pyarray<double> &q_velocity = args.array<double>("q_velocity");
    xt::pyarray<double> &ebqe_velocity = args.array<double>("ebqe_velocity");
    xt::pyarray<double> &flux = args.array<double>("flux");
    xt::pyarray<double> &elementResidual_h_save =
        args.array<double>("elementResidual_h_save");
    xt::pyarray<double> &Cx = args.array<double>("Cx");
    xt::pyarray<double> &Cy = args.array<double>("Cy");
    xt::pyarray<double> &CTx = args.array<double>("CTx");
    xt::pyarray<double> &CTy = args.array<double>("CTy");
    int numDOFsPerEqn = args.scalar<int>("numDOFsPerEqn");
    int NNZ = args.scalar<int>("NNZ");
    xt::pyarray<int> &csrRowIndeces_DofLoops =
        args.array<int>("csrRowIndeces_DofLoops");
    xt::pyarray<int> &csrColumnOffsets_DofLoops =
        args.array<int>("csrColumnOffsets_DofLoops");
    xt::pyarray<double> &lumped_mass_matrix =
        args.array<double>("lumped_mass_matrix");
    double cfl_run = args.scalar<double>("cfl_run");
    double eps = args.scalar<double>("eps");
    double hEps = args.scalar<double>("hEps");
    xt::pyarray<double> &hnp1_at_quad_point =
        args.array<double>("hnp1_at_quad_point");
    xt::pyarray<double> &hunp1_at_quad_point =
        args.array<double>("hunp1_at_quad_point");
    xt::pyarray<double> &hvnp1_at_quad_point =
        args.array<double>("hvnp1_at_quad_point");
    xt::pyarray<double> &hetanp1_at_quad_point =
        args.array<double>("hetanp1_at_quad_point");
    xt::pyarray<double> &hwnp1_at_quad_point =
        args.array<double>("hwnp1_at_quad_point");
    xt::pyarray<double> &hbetanp1_at_quad_point =
        args.array<double>("hbetanp1_at_quad_point");
    xt::pyarray<double> &extendedSourceTerm_hu =
        args.array<double>("extendedSourceTerm_hu");
    xt::pyarray<double> &extendedSourceTerm_hv =
        args.array<double>("extendedSourceTerm_hv");
    xt::pyarray<double> &extendedSourceTerm_heta =
        args.array<double>("extendedSourceTerm_heta");
    xt::pyarray<double> &extendedSourceTerm_hw =
        args.array<double>("extendedSourceTerm_hw");
    xt::pyarray<double> &extendedSourceTerm_hbeta =
        args.array<double>("extendedSourceTerm_hbeta");
    xt::pyarray<double> &dH_minus_dL = args.array<double>("dH_minus_dL");
    xt::pyarray<double> &muH_minus_muL = args.array<double>("muH_minus_muL");
    double cE = args.scalar<double>("cE");
    int LUMPED_MASS_MATRIX = args.scalar<int>("LUMPED_MASS_MATRIX");
    double dt = args.scalar<double>("dt");
    int LINEAR_FRICTION = args.scalar<int>("LINEAR_FRICTION");
    double mannings = args.scalar<double>("mannings");
    xt::pyarray<double> &quantDOFs = args.array<double>("quantDOFs");
    int SECOND_CALL_CALCULATE_RESIDUAL =
        args.scalar<int>("SECOND_CALL_CALCULATE_RESIDUAL");
    int COMPUTE_NORMALS = args.scalar<int>("COMPUTE_NORMALS");
    xt::pyarray<double> &normalx = args.array<double>("normalx");
    xt::pyarray<double> &normaly = args.array<double>("normaly");
    xt::pyarray<double> &dLow = args.array<double>("dLow");
    int lstage = args.scalar<int>("lstage");
    xt::pyarray<double> &new_SourceTerm_hu =
        args.array<double>("new_SourceTerm_hu");
    xt::pyarray<double> &new_SourceTerm_hv =
        args.array<double>("new_SourceTerm_hv");
    xt::pyarray<double> &new_SourceTerm_heta =
        args.array<double>("new_SourceTerm_heta");
    xt::pyarray<double> &new_SourceTerm_hw =
        args.array<double>("new_SourceTerm_hw");
    xt::pyarray<double> &global_entropy_residual =
        args.array<double>("global_entropy_residual");
    double dij_small = args.scalar<double>("dij_small");
    xt::pyarray<double> &hLow = args.array<double>("hLow");
    xt::pyarray<double> &huLow = args.array<double>("huLow");
    xt::pyarray<double> &hvLow = args.array<double>("hvLow");
    xt::pyarray<double> &hetaLow = args.array<double>("hetaLow");
    xt::pyarray<double> &hwLow = args.array<double>("hwLow");
    xt::pyarray<double> &hbetaLow = args.array<double>("hbetaLow");
    xt::pyarray<double> &h_min = args.array<double>("h_min");
    xt::pyarray<double> &h_max = args.array<double>("h_max");
    xt::pyarray<double> &heta_min = args.array<double>("heta_min");
    xt::pyarray<double> &heta_max = args.array<double>("heta_max");
    xt::pyarray<double> &kin_max = args.array<double>("kin_max");
    double size_of_domain = args.scalar<double>("size_of_domain");
    xt::pyarray<double> &urelax = args.array<double>("urelax");
    xt::pyarray<double> &drelax = args.array<double>("drelax");
    // FOR FRICTION//
    double n2 = std::pow(mannings, 2.);
    double gamma = 4. / 3;
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
          elementResidual_hw[nDOF_test_element],
          elementResidual_hbeta[nDOF_test_element];

      for (int i = 0; i < nDOF_test_element; i++) {
        elementResidual_h[i] = 0.0;
        elementResidual_hu[i] = 0.0;
        elementResidual_hv[i] = 0.0;
        elementResidual_heta[i] = 0.0;
        elementResidual_hw[i] = 0.0;
        elementResidual_hbeta[i] = 0.0;
      }
      //
      // loop over quadrature points and compute integrands
      //
      for (int k = 0; k < nQuadraturePoints_element; k++) {
        // compute indices and declare local storage
        register int eN_k = eN * nQuadraturePoints_element + k,
                     eN_k_nSpace = eN_k * nSpace,
                     eN_nDOF_trial_element = eN * nDOF_trial_element;
        register double h = 0.0, hu = 0.0, hv = 0.0, heta = 0.0, hw = 0.0,
                        hbeta = 0.0, // solution at current time
            h_old = 0.0, hu_old = 0.0, hv_old = 0.0, heta_old = 0.0,
                        hw_old = 0.0,
                        hbeta_old = 0.0, // solution at lstage
            jac[nSpace * nSpace], jacDet, jacInv[nSpace * nSpace],
                        h_test_dV[nDOF_trial_element], dV, x, y, xt, yt;
        // get jacobian, etc for mapping reference element
        ck.calculateMapping_element(
            eN, k, mesh_dof.data(), mesh_l2g.data(), mesh_trial_ref.data(),
            mesh_grad_trial_ref.data(), jac, jacDet, jacInv, x, y);
        // get the physical integration weight
        dV = fabs(jacDet) * dV_ref[k];
        // get the solution at current time
        ck.valFromDOF(h_dof.data(), &h_l2g.data()[eN_nDOF_trial_element],
                      &h_trial_ref.data()[k * nDOF_trial_element], h);
        ck.valFromDOF(hu_dof.data(), &vel_l2g[eN_nDOF_trial_element],
                      &vel_trial_ref.data()[k * nDOF_trial_element], hu);
        ck.valFromDOF(hv_dof.data(), &vel_l2g[eN_nDOF_trial_element],
                      &vel_trial_ref.data()[k * nDOF_trial_element], hv);
        ck.valFromDOF(heta_dof.data(), &vel_l2g.data()[eN_nDOF_trial_element],
                      &vel_trial_ref.data()[k * nDOF_trial_element], heta);
        ck.valFromDOF(hw_dof.data(), &vel_l2g.data()[eN_nDOF_trial_element],
                      &vel_trial_ref.data()[k * nDOF_trial_element], hw);
        ck.valFromDOF(hbeta_dof.data(), &vel_l2g.data()[eN_nDOF_trial_element],
                      &vel_trial_ref.data()[k * nDOF_trial_element], hbeta);
        // get the solution at the lstage
        ck.valFromDOF(h_dof_old.data(), &h_l2g.data()[eN_nDOF_trial_element],
                      &h_trial_ref.data()[k * nDOF_trial_element], h_old);
        ck.valFromDOF(hu_dof_old.data(), &vel_l2g.data()[eN_nDOF_trial_element],
                      &vel_trial_ref.data()[k * nDOF_trial_element], hu_old);
        ck.valFromDOF(hv_dof_old.data(), &vel_l2g.data()[eN_nDOF_trial_element],
                      &vel_trial_ref.data()[k * nDOF_trial_element], hv_old);
        ck.valFromDOF(heta_dof_old.data(), &vel_l2g[eN_nDOF_trial_element],
                      &vel_trial_ref.data()[k * nDOF_trial_element], heta_old);
        ck.valFromDOF(hw_dof_old.data(), &vel_l2g.data()[eN_nDOF_trial_element],
                      &vel_trial_ref.data()[k * nDOF_trial_element], hw_old);
        ck.valFromDOF(hbeta_dof_old.data(),
                      &vel_l2g.data()[eN_nDOF_trial_element],
                      &vel_trial_ref.data()[k * nDOF_trial_element], hbeta_old);
        // calculate cell based CFL to keep a reference
        calculateCFL(elementDiameter.data()[eN], g, h_old, hu_old, hv_old, hEps,
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
        hbetanp1_at_quad_point[eN_k] = hbeta;

        for (int i = 0; i < nDOF_test_element; i++) {
          // compute time derivative part of global residual. NOTE: no lumping
          elementResidual_h[i] += (h - h_old) * h_test_dV[i];
          elementResidual_hu[i] += (hu - hu_old) * h_test_dV[i];
          elementResidual_hv[i] += (hv - hv_old) * h_test_dV[i];
          elementResidual_heta[i] += (heta - heta_old) * h_test_dV[i];
          elementResidual_hw[i] += (hw - hw_old) * h_test_dV[i];
          elementResidual_hbeta[i] += (hbeta - hbeta_old) * h_test_dV[i];
        }
      }
      // distribute
      for (int i = 0; i < nDOF_test_element; i++) {
        register int eN_i = eN * nDOF_test_element + i;

        // global i-th index for h (this is same for vel_l2g)
        int h_gi = h_l2g[eN_i];

        // distribute time derivative to global residual
        globalResidual[offset_h + stride_h * h_gi] += elementResidual_h[i];
        globalResidual[offset_hu + stride_hu * h_gi] += elementResidual_hu[i];
        globalResidual[offset_hv + stride_hv * h_gi] += elementResidual_hv[i];
        globalResidual[offset_heta + stride_heta * h_gi] +=
            elementResidual_heta[i];
        globalResidual[offset_hw + stride_hw * h_gi] += elementResidual_hw[i];
        globalResidual[offset_hbeta + stride_hbeta * h_gi] +=
            elementResidual_hbeta[i];
      }
    }
    // ********** END OF CELL LOOPS ********** //

    if (SECOND_CALL_CALCULATE_RESIDUAL == 0) // This is to save some time
    {
      /////////////////////////////////////////////////////
      // ********** FIRST SET OF LOOP ON DOFs ********** //
      /////////////////////////////////////////////////////
      // To compute:
      //     * Local bounds for limiting
      //     * Low order solution (in terms of bar states)

      // Here we declare some arrays for local bounds
      std::valarray<double> delta_Sqd_h(0.0, numDOFsPerEqn),
          bar_deltaSqd_h(0.0, numDOFsPerEqn),
          delta_Sqd_heta(0.0, numDOFsPerEqn),
          bar_deltaSqd_heta(0.0, numDOFsPerEqn),
          delta_Sqd_kin(0.0, numDOFsPerEqn),
          bar_deltaSqd_kin(0.0, numDOFsPerEqn);
      xt::pyarray<double> kin(numDOFsPerEqn), max_of_h_and_hEps(numDOFsPerEqn);

      // Define kinetic energy, kin = 1/2 q^2 / h
      max_of_h_and_hEps = xt::where(h_dof_old > hEps, h_dof_old, hEps);
      kin = 0.5 * (hu_dof_old * hu_dof_old + hv_dof_old * hv_dof_old);
      kin *= 2.0 * h_dof_old /
             (h_dof_old * h_dof_old + max_of_h_and_hEps * max_of_h_and_hEps);

      /* First loop to define: delta_Sqd_h, delta_Sqd_heta, delta_Sqd_kin */
      for (int i = 0; i < numDOFsPerEqn; i++) {

        for (int offset = csrRowIndeces_DofLoops[i];
             offset < csrRowIndeces_DofLoops[i + 1]; offset++) {

          int j = csrColumnOffsets_DofLoops[offset];

          if (i != j) {
            delta_Sqd_h[i] += h_dof_old[i] - h_dof_old[j];
            delta_Sqd_heta[i] += heta_dof_old[i] - heta_dof_old[j];
            delta_Sqd_kin[i] += kin[i] - kin[j];
          }
        } // j loop ends here
      }   // i loops ends here

      // Stuff for bar states here (BT = BarTilde)
      std::valarray<double> hBT(0.0, dH_minus_dL.size()),
          huBT(0.0, dH_minus_dL.size()), hvBT(0.0, dH_minus_dL.size()),
          hetaBT(0.0, dH_minus_dL.size()), hwBT(0.0, dH_minus_dL.size()),
          hbetaBT(0.0, dH_minus_dL.size());

      /* Second loop to compute bar states (variable)BT */
      int ij = 0;
      for (int i = 0; i < numDOFsPerEqn; i++) {

        // define things at ith node
        double hi = h_dof_old[i];
        double hui = hu_dof_old[i];
        double hvi = hv_dof_old[i];
        double hetai = heta_dof_old[i];
        double hwi = hw_dof_old[i];
        double hbetai = hbeta_dof_old[i];
        double Zi = b_dof[i];
        double mi = lumped_mass_matrix[i];
        double one_over_hiReg =
            2.0 * hi / (hi * hi + std::pow(fmax(hi, hEps), 2));
        double ui = hui * one_over_hiReg;
        double vi = hvi * one_over_hiReg;
        double etai = hetai * one_over_hiReg;
        double meshSizei = std::sqrt(mi);

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

        // define full pressure at ith node for definition of bar states
        double pressure_i = 0.5 * g * hi * hi + pTildei;

        // loop over the sparsity pattern of the i-th DOF
        for (int offset = csrRowIndeces_DofLoops[i];
             offset < csrRowIndeces_DofLoops[i + 1]; offset++) {
          int j = csrColumnOffsets_DofLoops[offset];

          // define things at jth node
          double hj = h_dof_old[j];
          double huj = hu_dof_old[j];
          double hvj = hv_dof_old[j];
          double hetaj = heta_dof_old[j];
          double hwj = hw_dof_old[j];
          double hbetaj = hbeta_dof_old[j];
          double Zj = b_dof[j];
          double one_over_hjReg =
              2.0 * hj / (hj * hj + std::pow(fmax(hj, hEps), 2));
          double uj = huj * one_over_hjReg;
          double vj = hvj * one_over_hjReg;
          double etaj = hetaj * one_over_hjReg;
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

          // define pressure at jth node for bar states
          double pressure_j = 0.5 * g * hj * hj + pTildej;

          // Compute star states
          double hStarij = fmax(0., hi + Zi - fmax(Zi, Zj));
          double huStarij = hui * hStarij * one_over_hiReg;
          double hvStarij = hvi * hStarij * one_over_hiReg;
          double hetaStarij = hetai * std::pow(hStarij * one_over_hiReg, 2);
          double hwStarij = hwi * hStarij * one_over_hiReg;
          double hbetaStarij = hbetai * hStarij * one_over_hiReg;

          double hStarji = fmax(0., hj + Zj - fmax(Zi, Zj));
          double huStarji = huj * hStarji * one_over_hjReg;
          double hvStarji = hvj * hStarji * one_over_hjReg;
          double hetaStarji = hetaj * std::pow(hStarji * one_over_hjReg, 2);
          double hwStarji = hwj * hStarji * one_over_hjReg;
          double hbetaStarji = hbetaj * hStarji * one_over_hjReg;

          // for bar_deltaSqd_h, bar_deltaSqd_heta, bar_deltaSqd_kin
          double muLowij = 0., muLij = 0., dLowij = 0., dLij = 0.;
          if (i != j) {
            // Put these computations first  before it gets messy
            bar_deltaSqd_h[i] += delta_Sqd_h[j] + delta_Sqd_h[i];
            bar_deltaSqd_heta[i] += delta_Sqd_heta[j] + delta_Sqd_heta[i];
            bar_deltaSqd_kin[i] += delta_Sqd_kin[j] + delta_Sqd_kin[i];

            if (lstage == 0)
              dLowij = dLow[ij];
            else {
              double cij_norm = sqrt(Cx[ij] * Cx[ij] + Cy[ij] * Cy[ij]);
              double cji_norm = sqrt(CTx[ij] * CTx[ij] + CTy[ij] * CTy[ij]);
              double nxij = Cx[ij] / cij_norm, nyij = Cy[ij] / cij_norm;
              double nxji = CTx[ij] / cji_norm, nyji = CTy[ij] / cji_norm;
              dLowij =
                  fmax(maxWaveSpeedSharpInitialGuess(g, nxij, nyij, hi, hui,
                                                     hvi, hetai, mi, hj, huj,
                                                     hvj, hetaj, mj, hEps) *
                           cij_norm,
                       maxWaveSpeedSharpInitialGuess(g, nxji, nyji, hj, huj,
                                                     hvj, hetaj, mj, hi, hui,
                                                     hvi, hetai, mi, hEps) *
                           cji_norm);
            }
            // save dLij
            dLij = dLowij;

            // compute muij
            muLij = fmax(fmax(0., -(ui * Cx[ij] + vi * Cy[ij])),
                         fmax(0., (uj * Cx[ij] + vj * Cy[ij])));

            // Define dLij as max of dLij and muLij
            dLij = fmax(dLowij, muLij);
            dLow[ij] = fmax(dLij, muLij);

            ////////////////////////
            // COMPUTE BAR STATES //
            ////////////////////////
            double hBar_ij = 0, hTilde_ij = 0, huBar_ij = 0, huTilde_ij = 0,
                   hvBar_ij = 0, hvTilde_ij = 0, hetaBar_ij = 0,
                   hetaTilde_ij = 0, hwBar_ij = 0, hwTilde_ij = 0,
                   hbetaBar_ij = 0, hbetaTilde_ij = 0;

            // h component
            hBar_ij = -1. / (fmax(2.0 * dLij, dij_small)) *
                          ((uj * hj - ui * hi) * Cx[ij] +
                           (vj * hj - vi * hi) * Cy[ij]) +
                      0.5 * (hj + hi);
            hTilde_ij = (dLij - muLij) / (fmax(2.0 * dLij, dij_small)) *
                        ((hStarji - hj) - (hStarij - hi));
            // hu component
            huBar_ij =
                -1. / (fmax(2.0 * dLij, dij_small)) *
                    ((uj * huj - ui * hui + pressure_j - pressure_i) * Cx[ij] +
                     (vj * huj - vi * hui) * Cy[ij]) +
                0.5 * (huj + hui);
            huTilde_ij = (dLij - muLij) / (fmax(2.0 * dLij, dij_small)) *
                         ((huStarji - huj) - (huStarij - hui));
            // hv component
            hvBar_ij =
                -1. / (fmax(2.0 * dLij, dij_small)) *
                    ((uj * hvj - ui * hvi) * Cx[ij] +
                     (vj * hvj - vi * hvi + pressure_j - pressure_i) * Cy[ij]) +
                0.5 * (hvj + hvi);
            hvTilde_ij = (dLij - muLij) / (fmax(2.0 * dLij, dij_small)) *
                         ((hvStarji - hvj) - (hvStarij - hvi));
            // heta component
            hetaBar_ij = -1. / (fmax(2.0 * dLij, dij_small)) *
                             ((uj * hetaj - ui * hetai) * Cx[ij] +
                              (vj * hetaj - vi * hetai) * Cy[ij]) +
                         0.5 * (hetaj + hetai);
            hetaTilde_ij = (dLij - muLij) / (fmax(2.0 * dLij, dij_small)) *
                           ((hetaStarji - hetaj) - (hetaStarij - hetai));
            // hw component
            hwBar_ij = -1. / (fmax(2.0 * dLij, dij_small)) *
                           ((uj * hwj - ui * hwi) * Cx[ij] +
                            (vj * hwj - vi * hwi) * Cy[ij]) +
                       0.5 * (hwj + hwi);
            hwTilde_ij = (dLij - muLij) / (fmax(2.0 * dLij, dij_small)) *
                         ((hwStarji - hwj) - (hwStarij - hwi));
            // hbeta component
            hbetaBar_ij = -1. / (fmax(2.0 * dLij, dij_small)) *
                              ((uj * hbetaj - ui * hbetai) * Cx[ij] +
                               (vj * hbetaj - vi * hbetai) * Cy[ij]) +
                          0.5 * (hbetaj + hbetai);
            hbetaTilde_ij = (dLij - muLij) / (fmax(2.0 * dLij, dij_small)) *
                            ((hbetaStarji - hbetaj) - (hbetaStarij - hbetai));

            // Here we define uBar + uTilde
            hBT[ij] = hBar_ij + hTilde_ij;
            huBT[ij] = huBar_ij + huTilde_ij;
            hvBT[ij] = hvBar_ij + hvTilde_ij;
            hetaBT[ij] = hetaBar_ij + hetaTilde_ij;
            hwBT[ij] = hwBar_ij + hwTilde_ij;
            hbetaBT[ij] = hbetaBar_ij + hbetaTilde_ij;
          } else {
            // i==j
            // Bar states by definition satisfy Utilde_ii + Ubar_ii = U_i
            hBT[ij] = hi;
            huBT[ij] = hui;
            hvBT[ij] = hvi;
            hetaBT[ij] = hetai;
            hwBT[ij] = hwi;
            hbetaBT[ij] = hbetai;
          }

          // UPDATE ij //
          ij += 1;
        } // j loop ends here

        // for bar_deltaSqd_h, bar_deltaSqd_heta, bar_deltaSqd_kin
        bar_deltaSqd_h[i] =
            bar_deltaSqd_h[i] /
            (csrRowIndeces_DofLoops[i + 1] - csrRowIndeces_DofLoops[i] - 1) /
            2.0;
        bar_deltaSqd_heta[i] =
            bar_deltaSqd_heta[i] /
            (csrRowIndeces_DofLoops[i + 1] - csrRowIndeces_DofLoops[i] - 1) /
            2.0;
        bar_deltaSqd_kin[i] =
            bar_deltaSqd_kin[i] /
            (csrRowIndeces_DofLoops[i + 1] - csrRowIndeces_DofLoops[i] - 1) /
            2.0;

      } // i loops ends here

      /* Then final loop of first set to get local bounds */
      ij = 0;
      for (int i = 0; i < numDOFsPerEqn; i++) {

        // define m_i
        double mi = lumped_mass_matrix[i];

        /* Initialize hmin, hmax, heta_min, heta_max */
        h_min[i] = h_dof_old[i];
        h_max[i] = h_dof_old[i];
        heta_min[i] = heta_dof_old[i];
        heta_max[i] = heta_dof_old[i];

        /* Initialize low order solution */
        hLow[i] = h_dof_old[i];
        huLow[i] = hu_dof_old[i];
        hvLow[i] = hv_dof_old[i];
        hetaLow[i] = heta_dof_old[i];
        hwLow[i] = hw_dof_old[i];
        hbetaLow[i] = hbeta_dof_old[i];
        kin_max[i] = kin[i];

        // loop in j (sparsity pattern)
        for (int offset = csrRowIndeces_DofLoops[i];
             offset < csrRowIndeces_DofLoops[i + 1]; offset++) {

          int j = csrColumnOffsets_DofLoops[offset];

          double one_over_hBT =
              2.0 * hBT[ij] /
              (hBT[ij] * hBT[ij] + std::pow(fmax(hBT[ij], hEps), 2));
          double psi_ij =
              one_over_hBT * 0.5 *
              (huBT[ij] * huBT[ij] + hvBT[ij] * hvBT[ij]); // Eqn (6.31)

          // COMPUTE LOCAL BOUNDS //
          kin_max[i] = fmax(psi_ij, kin_max[i]);
          h_min[i] = std::min(h_min[i], hBT[ij]);
          h_max[i] = std::max(h_max[i], hBT[ij]);
          heta_min[i] = std::min(heta_min[i], hetaBT[ij]);
          heta_max[i] = std::max(heta_max[i], hetaBT[ij]);

          // Then do relaxation of bounds here. If confused, see convex
          // limiting paper
          kin_max[i] =
              std::min(urelax[i] * kin_max[i],
                       kin_max[i] + std::abs(bar_deltaSqd_kin[i]) / 2.0);
          h_min[i] = std::max(drelax[i] * h_min[i],
                              h_min[i] - std::abs(bar_deltaSqd_h[i]) / 2.0);
          h_max[i] = std::min(urelax[i] * h_max[i],
                              h_max[i] + std::abs(bar_deltaSqd_h[i]) / 2.0);
          heta_min[i] =
              std::max(drelax[i] * heta_min[i],
                       heta_min[i] - std::abs(bar_deltaSqd_heta[i]) / 2.0);
          heta_max[i] =
              std::min(urelax[i] * heta_max[i],
                       heta_max[i] + std::abs(bar_deltaSqd_heta[i]) / 2.0);

          /* COMPUTE LOW ORDER SOLUTION. See EQN 6.23 in SW friction paper */
          // This is low order solution WITHOUT sources
          if (i != j) {
            hLow[i] += h_dof_old[i] * (-dt / mi * 2 * dLow[ij]) +
                       dt / mi * (2 * dLow[ij] * hBT[ij]);
            huLow[i] += hu_dof_old[i] * (-dt / mi * 2 * dLow[ij]) +
                        dt / mi * (2 * dLow[ij] * huBT[ij]);
            hvLow[i] += hv_dof_old[i] * (-dt / mi * 2 * dLow[ij]) +
                        dt / mi * (2 * dLow[ij] * hvBT[ij]);
            hetaLow[i] += heta_dof_old[i] * (-dt / mi * 2 * dLow[ij]) +
                          dt / mi * (2 * dLow[ij] * hetaBT[ij]);
            hwLow[i] += hw_dof_old[i] * (-dt / mi * 2 * dLow[ij]) +
                        dt / mi * (2 * dLow[ij] * hwBT[ij]);
            hbetaLow[i] += hbeta_dof_old[i] * (-dt / mi * 2 * dLow[ij]) +
                           dt / mi * (2 * dLow[ij] * hbetaBT[ij]);
          }
          // UPDATE ij //
          ij += 1;
        } // j loop ends here

        // clean up hLow from round off error
        if (hLow[i] < hEps) {
          hLow[i] = 0.0;
          hetaLow[i] = 0.0;
        }
      } // i loop ends here

      ////////////////////////////////////////////////////////
      // ********** Second set of loops on dofs ********** //
      ///////////////////////////////////////////////////////
      // To compute:
      //     * Hyperbolic part of the flux
      //     * Extended source terms
      //     * Smoothness indicator

      ij = 0;
      std::valarray<double> hyp_flux_h(numDOFsPerEqn),
          hyp_flux_hu(numDOFsPerEqn), hyp_flux_hv(numDOFsPerEqn),
          hyp_flux_heta(numDOFsPerEqn), hyp_flux_hw(numDOFsPerEqn),
          hyp_flux_hbeta(numDOFsPerEqn), psi(numDOFsPerEqn),
          etaMax(numDOFsPerEqn), etaMin(numDOFsPerEqn);

      for (int i = 0; i < numDOFsPerEqn; i++) {
        // solution at time tn for the ith DOF
        double hi = h_dof_old[i];
        double hui = hu_dof_old[i];
        double hvi = hv_dof_old[i];
        double hetai = heta_dof_old[i];
        double hwi = hw_dof_old[i];
        double hbetai = hbeta_dof_old[i];
        double Zi = b_dof[i];
        // Define some things using above
        double one_over_hiReg =
            2 * hi / (hi * hi + std::pow(fmax(hi, hEps), 2)); // hEps
        double ui = hui * one_over_hiReg;
        double vi = hvi * one_over_hiReg;
        double etai = hetai * one_over_hiReg;
        double mi = lumped_mass_matrix[i];
        double meshSizei = std::sqrt(mi);

        /* COMPUTE EXTENDED SOURCE TERMS for all equations:
         * Friction terms
         * NOTE: Be careful with sign of source terms. Extended sources are on
         * left side of equations. "new_SourceTerm" variables are on right
         */

        // FRICTION
        if (LINEAR_FRICTION == 1) {
          extendedSourceTerm_hu[i] = mannings * hui * mi;
          extendedSourceTerm_hv[i] = mannings * hvi * mi;
          // For use in the convex limiting function
          // actually didn't need to do this but it helps with signs
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
          // For use in the convex limiting function
          new_SourceTerm_hu[i] = -friction_aux * hui;
          new_SourceTerm_hv[i] = -friction_aux * hvi;
        }

        // Define some things for heta, hw and hbeta sources
        double ratio_i = 1.0; //(2.0 * hetai) / (etai * etai + hi * hi + hEps);
        double diff_over_h_i = (hetai - hi * hi) * one_over_hiReg;
        double hSqd_GammaPi = 6.0 * (hetai - hi * hi);
        double s_i = LAMBDA_MGN * g / meshSizei * hSqd_GammaPi;
        double mu_bar_i = LAMBDA_MGN * std::sqrt(g * hEps / eps) / meshSizei;
        double grad_Z_x_i = 0.0;
        double grad_Z_y_i = 0.0;
        double q_dot_gradZ = 0.0;

        if (IF_BOTH_GAMMA_BRANCHES) {
          if (hetai > std::pow(hi, 2.0)) {
            hSqd_GammaPi = 6.0 * etai * diff_over_h_i;
          }
        }

        // heta source first part
        extendedSourceTerm_heta[i] = -hwi * mi * ratio_i;
        new_SourceTerm_heta[i] = hwi * mi * ratio_i;

        // hw source
        extendedSourceTerm_hw[i] =
            (LAMBDA_MGN * g / meshSizei) * hSqd_GammaPi * mi * ratio_i;
        new_SourceTerm_hw[i] =
            -(LAMBDA_MGN * g / meshSizei) * hSqd_GammaPi * mi * ratio_i;

        /* HYPERBOLIC FLUXES */
        hyp_flux_h[i] = 0;
        hyp_flux_hu[i] = 0;
        hyp_flux_hv[i] = 0;
        hyp_flux_heta[i] = 0;
        hyp_flux_hw[i] = 0;
        hyp_flux_hbeta[i] = 0;

        // FOR SMOOTHNESS INDICATOR //
        double alphai; // smoothness indicator of solution
        double alpha_numerator = 0;
        double alpha_denominator = 0;
        double alpha_zero = 0.5; // if only want smoothness
        double alpha_factor = 1.0 / (1.0 - alpha_zero);

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
          double hbetaj = hbeta_dof_old[j];
          double Zj = b_dof[j];
          // Then define some things here using above
          double one_over_hjReg =
              2.0 * hj / (hj * hj + std::pow(fmax(hj, hEps), 2));
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
              (uj * hj - ui * hi) * Cx[ij] + (vj * hj - vi * hi) * Cy[ij];
          double aux_hu =
              (uj * huj - ui * hui) * Cx[ij] + (vj * huj - vi * hui) * Cy[ij];
          double aux_hv =
              (uj * hvj - ui * hvi) * Cx[ij] + (vj * hvj - vi * hvi) * Cy[ij];
          double aux_heta = (uj * hetaj - ui * hetai) * Cx[ij] +
                            (vj * hetaj - vi * hetai) * Cy[ij];
          double aux_hw =
              (uj * hwj - ui * hwi) * Cx[ij] + (vj * hwj - vi * hwi) * Cy[ij];
          double aux_hbeta = (uj * hbetaj - ui * hbetai) * Cx[ij] +
                             (vj * hbetaj - vi * hbetai) * Cy[ij];

          // Define grad_Z_x and grad_Z_y
          grad_Z_x_i += Zj * Cx[ij];
          grad_Z_y_i += Zj * Cy[ij];

          /* HYPERBOLIC FLUX */
          hyp_flux_h[i] += aux_h;
          hyp_flux_hu[i] += aux_hu + pTildej * Cx[ij];
          hyp_flux_hv[i] += aux_hv + pTildej * Cy[ij];
          hyp_flux_heta[i] += aux_heta;
          hyp_flux_hw[i] += aux_hw;
          hyp_flux_hbeta[i] += aux_hbeta;

          // EXTENDED SOURCE, USING 6.13 //
          extendedSourceTerm_hu[i] += g * hi * (hj + Zj) * Cx[ij];
          extendedSourceTerm_hv[i] += g * hi * (hj + Zj) * Cy[ij];

          new_SourceTerm_hu[i] +=
              g * (-hi * (Zj - Zi) + 0.5 * std::pow(hj - hi, 2)) * Cx[ij];
          new_SourceTerm_hv[i] +=
              g * (-hi * (Zj - Zi) + 0.5 * std::pow(hj - hi, 2)) * Cy[ij];

          // FOR SMOOTHNESS INDICATOR //
          alpha_numerator += hj - hi;
          alpha_denominator += fabs(hj - hi);

          // update ij
          ij += 1;
        } // end j loop

        // update sources with gradZ contribution
        q_dot_gradZ = hui * grad_Z_x_i + hvi * grad_Z_y_i;
        double Spi = mu_bar_i * (q_dot_gradZ / mi - hbetai);
        double r_i = -s_i / 2.0 + Spi / 4.0;
        // q1
        extendedSourceTerm_heta[i] += 1.5 * q_dot_gradZ;
        new_SourceTerm_heta -= 1.5 * q_dot_gradZ;
        // q3
        extendedSourceTerm_hbeta[i] = -mi * Spi;
        // momentum equations
        extendedSourceTerm_hu[i] -= -r_i * grad_Z_x_i;
        new_SourceTerm_hu[i] += -r_i * grad_Z_x_i;
        extendedSourceTerm_hv[i] -= -r_i * grad_Z_y_i;
        new_SourceTerm_hv[i] += -r_i * grad_Z_y_i;

        // COMPUTE SMOOTHNESS INDICATOR //
        if (hi <= hEps) {
          alphai = 1.0;
          psi[i] = 1.0;
        } else {
          // Force alphai=0 in constant states
          if (fabs(alpha_numerator) <= hEps) {
            alphai = 0.;
          } else {
            alphai =
                (fabs(alpha_numerator) - hEps) / fabs(alpha_denominator - hEps);
          }
          alphai = fmax(alphai - alpha_zero, 0.0) * alpha_factor;
          if (POWER_SMOOTHNESS_INDICATOR == 0)
            psi[i] = 1.0;
          else
            psi[i] = std::pow(
                alphai, POWER_SMOOTHNESS_INDICATOR); // NOTE: alpha^4 for mGN
        }
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
        double hbetai = hbeta_dof_old[i];
        double Zi = b_dof[i];
        double mi = lumped_mass_matrix[i];

        double one_over_hiReg =
            2.0 * hi / (hi * hi + std::pow(fmax(hi, hEps), 2));
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
               ith_dHij_minus_muHij_times_hbetaStarStates = 0.,
               ith_muHij_times_hStates = 0., ith_muHij_times_huStates = 0.,
               ith_muHij_times_hvStates = 0., ith_muHij_times_hetaStates = 0.,
               ith_muHij_times_hwStates = 0., ith_muHij_times_hbetaStates = 0.;

        // loop over the sparsity pattern of the i-th DOF
        for (int offset = csrRowIndeces_DofLoops[i];
             offset < csrRowIndeces_DofLoops[i + 1]; offset++) {

          int j = csrColumnOffsets_DofLoops[offset];

          double hj = h_dof_old[j];
          double huj = hu_dof_old[j];
          double hvj = hv_dof_old[j];
          double hetaj = heta_dof_old[j];
          double hwj = hw_dof_old[j];
          double hbetaj = hbeta_dof_old[j];
          double Zj = b_dof[j];

          double one_over_hjReg =
              2.0 * hj / (hj * hj + std::pow(fmax(hj, hEps), 2));
          double uj = huj * one_over_hjReg;
          double vj = hvj * one_over_hjReg;
          double etaj = hetaj * one_over_hjReg;
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

          // COMPUTE STAR SOLUTION // hStar, huStar, hvStar, hetaStar, hwStar,
          // hbetaStar
          double hStarij = fmax(0., hi + Zi - fmax(Zi, Zj));
          double huStarij = hui * hStarij * one_over_hiReg;
          double hvStarij = hvi * hStarij * one_over_hiReg;
          double hetaStarij = hetai * std::pow(hStarij * one_over_hiReg, 2);
          double hwStarij = hwi * hStarij * one_over_hiReg;
          double hbetaStarij = hbetai * hStarij * one_over_hiReg;

          double hStarji = fmax(0., hj + Zj - fmax(Zi, Zj));
          double huStarji = huj * hStarji * one_over_hjReg;
          double hvStarji = hvj * hStarji * one_over_hjReg;
          double hetaStarji = hetaj * std::pow(hStarji * one_over_hjReg, 2);
          double hwStarji = hwj * hStarji * one_over_hjReg;
          double hbetaStarji = hbetaj * hStarji * one_over_hjReg;

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
              dLowij =
                  fmax(maxWaveSpeedSharpInitialGuess(g, nxij, nyij, hi, hui,
                                                     hvi, hetai, mi, hj, huj,
                                                     hvj, hetaj, mj, hEps) *
                           cij_norm,
                       maxWaveSpeedSharpInitialGuess(g, nxji, nyji, hj, huj,
                                                     hvj, hetaj, mj, hi, hui,
                                                     hvi, hetai, mi, hEps) *
                           cji_norm);
            }
            // this is standard low-order dij, can you use dLij =
            // dLowij*fmax(psi[i],psi[j]) if want smoothness based as low
            // order
            dLij = dLowij;

            ///////////////////////////////////////
            // WELL BALANCING DISSIPATIVE MATRIX //
            ///////////////////////////////////////
            muLowij = fmax(fmax(0., -(ui * Cx[ij] + vi * Cy[ij])),
                           fmax(0., (uj * Cx[ij] + vj * Cy[ij])));
            muLij = muLowij;

            // Define dLij as low order dijs
            muLij = muLowij;
            dLij = fmax(dLowij, muLij);

            // Then save dLow for limiting step, maybe a bit confusing
            dLow[ij] = fmax(dLij, muLij);

            ///////////////////////
            // ENTROPY VISCOSITY //
            ///////////////////////
            double dEVij = cE * fmax(global_entropy_residual[i],
                                     global_entropy_residual[j]);
            dHij = fmin(dLowij, dEVij);
            muHij = fmin(muLowij, dEVij);

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
            ith_dHij_minus_muHij_times_hbetaStarStates +=
                (dHij - muHij) * (hbetaStarji - hbetaStarij);

            // compute muij times solution terms
            ith_muHij_times_hStates += muHij * (hj - hi);
            ith_muHij_times_huStates += muHij * (huj - hui);
            ith_muHij_times_hvStates += muHij * (hvj - hvi);
            ith_muHij_times_hetaStates += muHij * (hetaj - hetai);
            ith_muHij_times_hwStates += muHij * (hwj - hwi);
            ith_muHij_times_hbetaStates += muHij * (hbetaj - hbetai);

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
        } // j loop ends here

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
          globalResidual[offset_hbeta + stride_hbeta * i] =
              hbetai - dt / mi *
                           ((hyp_flux_hbeta[i] + extendedSourceTerm_hbeta[i]) -
                            ith_dHij_minus_muHij_times_hbetaStarStates -
                            ith_muHij_times_hbetaStates);

          // clean up potential negative water height due to machine
          // precision
          if (globalResidual[offset_h + stride_h * i] >= -hEps &&
              globalResidual[offset_h + stride_h * i] < hEps) {
            globalResidual[offset_h + stride_h * i] = 0.0;
            globalResidual[offset_heta + stride_heta * i] = 0.0;
          }

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
          globalResidual[offset_hbeta + stride_hbeta * i] +=
              dt * (hyp_flux_hbeta[i] + extendedSourceTerm_hbeta[i] -
                    ith_dHij_minus_muHij_times_hbetaStarStates -
                    ith_muHij_times_hbetaStates);
        }
      }
      // ********** END OF LOOP IN DOFs ********** //
    } // end SECOND_CALL_CALCULATE_RESIDUAL

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
              eN, ebN_local, kb, ebN_local_kb, mesh_dof.data(), mesh_l2g.data(),
              mesh_trial_trace_ref.data(), mesh_grad_trial_trace_ref.data(),
              boundaryJac_ref.data(), jac_ext, jacDet_ext, jacInv_ext,
              boundaryJac, metricTensor, metricTensorDetSqrt, normal_ref.data(),
              normal, x_ext, y_ext);
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
  } // end calculateResidual

  void calculateMassMatrix(arguments_dict &args) {
    xt::pyarray<double> &mesh_trial_ref = args.array<double>("mesh_trial_ref");
    xt::pyarray<double> &mesh_grad_trial_ref =
        args.array<double>("mesh_grad_trial_ref");
    xt::pyarray<double> &mesh_dof = args.array<double>("mesh_dof");
    xt::pyarray<double> &mesh_velocity_dof =
        args.array<double>("mesh_velocity_dof");
    double MOVING_DOMAIN = args.scalar<double>("MOVING_DOMAIN");
    xt::pyarray<int> &mesh_l2g = args.array<int>("mesh_l2g");
    xt::pyarray<double> &dV_ref = args.array<double>("dV_ref");
    xt::pyarray<double> &h_trial_ref = args.array<double>("h_trial_ref");
    xt::pyarray<double> &h_grad_trial_ref =
        args.array<double>("h_grad_trial_ref");
    xt::pyarray<double> &h_test_ref = args.array<double>("h_test_ref");
    xt::pyarray<double> &h_grad_test_ref =
        args.array<double>("h_grad_test_ref");
    xt::pyarray<double> &vel_trial_ref = args.array<double>("vel_trial_ref");
    xt::pyarray<double> &vel_grad_trial_ref =
        args.array<double>("vel_grad_trial_ref");
    xt::pyarray<double> &vel_test_ref = args.array<double>("vel_test_ref");
    xt::pyarray<double> &vel_grad_test_ref =
        args.array<double>("vel_grad_test_ref");
    xt::pyarray<double> &mesh_trial_trace_ref =
        args.array<double>("mesh_trial_trace_ref");
    xt::pyarray<double> &mesh_grad_trial_trace_ref =
        args.array<double>("mesh_grad_trial_trace_ref");
    xt::pyarray<double> &dS_ref = args.array<double>("dS_ref");
    xt::pyarray<double> &h_trial_trace_ref =
        args.array<double>("h_trial_trace_ref");
    xt::pyarray<double> &h_grad_trial_trace_ref =
        args.array<double>("h_grad_trial_trace_ref");
    xt::pyarray<double> &h_test_trace_ref =
        args.array<double>("h_test_trace_ref");
    xt::pyarray<double> &h_grad_test_trace_ref =
        args.array<double>("h_grad_test_trace_ref");
    xt::pyarray<double> &vel_trial_trace_ref =
        args.array<double>("vel_trial_trace_ref");
    xt::pyarray<double> &vel_grad_trial_trace_ref =
        args.array<double>("vel_grad_trial_trace_ref");
    xt::pyarray<double> &vel_test_trace_ref =
        args.array<double>("vel_test_trace_ref");
    xt::pyarray<double> &vel_grad_test_trace_ref =
        args.array<double>("vel_grad_test_trace_ref");
    xt::pyarray<double> &normal_ref = args.array<double>("normal_ref");
    xt::pyarray<double> &boundaryJac_ref =
        args.array<double>("boundaryJac_ref");
    xt::pyarray<double> &elementDiameter =
        args.array<double>("elementDiameter");
    int nElements_global = args.scalar<int>("nElements_global");
    double g = args.scalar<double>("g");
    xt::pyarray<int> &h_l2g = args.array<int>("h_l2g");
    xt::pyarray<int> &vel_l2g = args.array<int>("vel_l2g");
    xt::pyarray<double> &b_dof = args.array<double>("b_dof");
    xt::pyarray<double> &h_dof = args.array<double>("h_dof");
    xt::pyarray<double> &hu_dof = args.array<double>("hu_dof");
    xt::pyarray<double> &hv_dof = args.array<double>("hv_dof");
    xt::pyarray<double> &heta_dof = args.array<double>("heta_dof");
    xt::pyarray<double> &hw_dof = args.array<double>("hw_dof");
    xt::pyarray<double> &hbeta_dof = args.array<double>("hbeta_dof");
    xt::pyarray<double> &q_cfl = args.array<double>("q_cfl");
    xt::pyarray<int> &sdInfo_hu_hu_rowptr =
        args.array<int>("sdInfo_hu_hu_rowptr");
    xt::pyarray<int> &sdInfo_hu_hu_colind =
        args.array<int>("sdInfo_hu_hu_colind");
    xt::pyarray<int> &sdInfo_hu_hv_rowptr =
        args.array<int>("sdInfo_hu_hv_rowptr");
    xt::pyarray<int> &sdInfo_hu_hv_colind =
        args.array<int>("sdInfo_hu_hv_colind");
    xt::pyarray<int> &sdInfo_hv_hv_rowptr =
        args.array<int>("sdInfo_hv_hv_rowptr");
    xt::pyarray<int> &sdInfo_hv_hv_colind =
        args.array<int>("sdInfo_hv_hv_colind");
    xt::pyarray<int> &sdInfo_hv_hu_rowptr =
        args.array<int>("sdInfo_hv_hu_rowptr");
    xt::pyarray<int> &sdInfo_hv_hu_colind =
        args.array<int>("sdInfo_hv_hu_colind");
    xt::pyarray<int> &csrRowIndeces_h_h = args.array<int>("csrRowIndeces_h_h");
    xt::pyarray<int> &csrColumnOffsets_h_h =
        args.array<int>("csrColumnOffsets_h_h");
    xt::pyarray<int> &csrRowIndeces_h_hu =
        args.array<int>("csrRowIndeces_h_hu");
    xt::pyarray<int> &csrColumnOffsets_h_hu =
        args.array<int>("csrColumnOffsets_h_hu");
    xt::pyarray<int> &csrRowIndeces_h_hv =
        args.array<int>("csrRowIndeces_h_hv");
    xt::pyarray<int> &csrColumnOffsets_h_hv =
        args.array<int>("csrColumnOffsets_h_hv");
    xt::pyarray<int> &csrRowIndeces_h_heta =
        args.array<int>("csrRowIndeces_h_heta");
    xt::pyarray<int> &csrColumnOffsets_h_heta =
        args.array<int>("csrColumnOffsets_h_heta");
    xt::pyarray<int> &csrRowIndeces_h_hw =
        args.array<int>("csrRowIndeces_h_hw");
    xt::pyarray<int> &csrColumnOffsets_h_hw =
        args.array<int>("csrColumnOffsets_h_hw");
    xt::pyarray<int> &csrRowIndeces_h_hbeta =
        args.array<int>("csrRowIndeces_h_hbeta");
    xt::pyarray<int> &csrColumnOffsets_h_hbeta =
        args.array<int>("csrColumnOffsets_h_hbeta");
    xt::pyarray<int> &csrRowIndeces_hu_h =
        args.array<int>("csrRowIndeces_hu_h");
    xt::pyarray<int> &csrColumnOffsets_hu_h =
        args.array<int>("csrColumnOffsets_hu_h");
    xt::pyarray<int> &csrRowIndeces_hu_hu =
        args.array<int>("csrRowIndeces_hu_hu");
    xt::pyarray<int> &csrColumnOffsets_hu_hu =
        args.array<int>("csrColumnOffsets_hu_hu");
    xt::pyarray<int> &csrRowIndeces_hu_hv =
        args.array<int>("csrRowIndeces_hu_hv");
    xt::pyarray<int> &csrColumnOffsets_hu_hv =
        args.array<int>("csrColumnOffsets_hu_hv");
    xt::pyarray<int> &csrRowIndeces_hu_heta =
        args.array<int>("csrRowIndeces_hu_heta");
    xt::pyarray<int> &csrColumnOffsets_hu_heta =
        args.array<int>("csrColumnOffsets_hu_heta");
    xt::pyarray<int> &csrRowIndeces_hu_hw =
        args.array<int>("csrRowIndeces_hu_hw");
    xt::pyarray<int> &csrColumnOffsets_hu_hw =
        args.array<int>("csrColumnOffsets_hu_hw");
    xt::pyarray<int> &csrRowIndeces_hu_hbeta =
        args.array<int>("csrRowIndeces_hu_hbeta");
    xt::pyarray<int> &csrColumnOffsets_hu_hbeta =
        args.array<int>("csrColumnOffsets_hu_hbeta");
    xt::pyarray<int> &csrRowIndeces_hv_h =
        args.array<int>("csrRowIndeces_hv_h");
    xt::pyarray<int> &csrColumnOffsets_hv_h =
        args.array<int>("csrColumnOffsets_hv_h");
    xt::pyarray<int> &csrRowIndeces_hv_hu =
        args.array<int>("csrRowIndeces_hv_hu");
    xt::pyarray<int> &csrColumnOffsets_hv_hu =
        args.array<int>("csrColumnOffsets_hv_hu");
    xt::pyarray<int> &csrRowIndeces_hv_hv =
        args.array<int>("csrRowIndeces_hv_hv");
    xt::pyarray<int> &csrColumnOffsets_hv_hv =
        args.array<int>("csrColumnOffsets_hv_hv");
    xt::pyarray<int> &csrRowIndeces_hv_heta =
        args.array<int>("csrRowIndeces_hv_heta");
    xt::pyarray<int> &csrColumnOffsets_hv_heta =
        args.array<int>("csrColumnOffsets_hv_heta");
    xt::pyarray<int> &csrRowIndeces_hv_hw =
        args.array<int>("csrRowIndeces_hv_hw");
    xt::pyarray<int> &csrColumnOffsets_hv_hw =
        args.array<int>("csrColumnOffsets_hv_hw");
    xt::pyarray<int> &csrRowIndeces_hv_hbeta =
        args.array<int>("csrRowIndeces_hv_hbeta");
    xt::pyarray<int> &csrColumnOffsets_hv_hbeta =
        args.array<int>("csrColumnOffsets_hv_hbeta");
    xt::pyarray<int> &csrRowIndeces_heta_h =
        args.array<int>("csrRowIndeces_heta_h");
    xt::pyarray<int> &csrColumnOffsets_heta_h =
        args.array<int>("csrColumnOffsets_heta_h");
    xt::pyarray<int> &csrRowIndeces_heta_hu =
        args.array<int>("csrRowIndeces_heta_hu");
    xt::pyarray<int> &csrColumnOffsets_heta_hu =
        args.array<int>("csrColumnOffsets_heta_hu");
    xt::pyarray<int> &csrRowIndeces_heta_hv =
        args.array<int>("csrRowIndeces_heta_hv");
    xt::pyarray<int> &csrColumnOffsets_heta_hv =
        args.array<int>("csrColumnOffsets_heta_hv");
    xt::pyarray<int> &csrRowIndeces_heta_heta =
        args.array<int>("csrRowIndeces_heta_heta");
    xt::pyarray<int> &csrColumnOffsets_heta_heta =
        args.array<int>("csrColumnOffsets_heta_heta");
    xt::pyarray<int> &csrRowIndeces_heta_hw =
        args.array<int>("csrRowIndeces_heta_hw");
    xt::pyarray<int> &csrColumnOffsets_heta_hw =
        args.array<int>("csrColumnOffsets_heta_hw");
    xt::pyarray<int> &csrRowIndeces_heta_hbeta =
        args.array<int>("csrRowIndeces_heta_hbeta");
    xt::pyarray<int> &csrColumnOffsets_heta_hbeta =
        args.array<int>("csrColumnOffsets_heta_hbeta");
    xt::pyarray<int> &csrRowIndeces_hw_h =
        args.array<int>("csrRowIndeces_hw_h");
    xt::pyarray<int> &csrColumnOffsets_hw_h =
        args.array<int>("csrColumnOffsets_hw_h");
    xt::pyarray<int> &csrRowIndeces_hw_hu =
        args.array<int>("csrRowIndeces_hw_hu");
    xt::pyarray<int> &csrColumnOffsets_hw_hu =
        args.array<int>("csrColumnOffsets_hw_hu");
    xt::pyarray<int> &csrRowIndeces_hw_hv =
        args.array<int>("csrRowIndeces_hw_hv");
    xt::pyarray<int> &csrColumnOffsets_hw_hv =
        args.array<int>("csrColumnOffsets_hw_hv");
    xt::pyarray<int> &csrRowIndeces_hw_heta =
        args.array<int>("csrRowIndeces_hw_heta");
    xt::pyarray<int> &csrColumnOffsets_hw_heta =
        args.array<int>("csrColumnOffsets_hw_heta");
    xt::pyarray<int> &csrRowIndeces_hw_hw =
        args.array<int>("csrRowIndeces_hw_hw");
    xt::pyarray<int> &csrColumnOffsets_hw_hw =
        args.array<int>("csrColumnOffsets_hw_hw");
    xt::pyarray<int> &csrRowIndeces_hw_hbeta =
        args.array<int>("csrRowIndeces_hw_hbeta");
    xt::pyarray<int> &csrColumnOffsets_hw_hbeta =
        args.array<int>("csrColumnOffsets_hw_hbeta");
    //
    xt::pyarray<int> &csrRowIndeces_hbeta_h =
        args.array<int>("csrRowIndeces_hbeta_h");
    xt::pyarray<int> &csrColumnOffsets_hbeta_h =
        args.array<int>("csrColumnOffsets_hbeta_h");
    xt::pyarray<int> &csrRowIndeces_hbeta_hu =
        args.array<int>("csrRowIndeces_hbeta_hu");
    xt::pyarray<int> &csrColumnOffsets_hbeta_hu =
        args.array<int>("csrColumnOffsets_hbeta_hu");
    xt::pyarray<int> &csrRowIndeces_hbeta_hv =
        args.array<int>("csrRowIndeces_hbeta_hv");
    xt::pyarray<int> &csrColumnOffsets_hbeta_hv =
        args.array<int>("csrColumnOffsets_hbeta_hv");
    xt::pyarray<int> &csrRowIndeces_hbeta_heta =
        args.array<int>("csrRowIndeces_hbeta_heta");
    xt::pyarray<int> &csrColumnOffsets_hbeta_heta =
        args.array<int>("csrColumnOffsets_hbeta_heta");
    xt::pyarray<int> &csrRowIndeces_hbeta_hw =
        args.array<int>("csrRowIndeces_hbeta_hw");
    xt::pyarray<int> &csrColumnOffsets_hbeta_hw =
        args.array<int>("csrColumnOffsets_hbeta_hw");
    xt::pyarray<int> &csrRowIndeces_hbeta_hbeta =
        args.array<int>("csrRowIndeces_hbeta_hbeta");
    xt::pyarray<int> &csrColumnOffsets_hbeta_hbeta =
        args.array<int>("csrColumnOffsets_hbeta_hbeta");
    xt::pyarray<double> &globalJacobian = args.array<double>("globalJacobian");
    int nExteriorElementBoundaries_global =
        args.scalar<int>("nExteriorElementBoundaries_global");
    xt::pyarray<int> &exteriorElementBoundariesArray =
        args.array<int>("exteriorElementBoundariesArray");
    xt::pyarray<int> &elementBoundaryElementsArray =
        args.array<int>("elementBoundaryElementsArray");
    xt::pyarray<int> &elementBoundaryLocalElementBoundariesArray =
        args.array<int>("elementBoundaryLocalElementBoundariesArray");
    xt::pyarray<int> &isDOFBoundary_h = args.array<int>("isDOFBoundary_h");
    xt::pyarray<int> &isDOFBoundary_hu = args.array<int>("isDOFBoundary_hu");
    xt::pyarray<int> &isDOFBoundary_hv = args.array<int>("isDOFBoundary_hv");
    xt::pyarray<int> &isAdvectiveFluxBoundary_h =
        args.array<int>("isAdvectiveFluxBoundary_h");
    xt::pyarray<int> &isAdvectiveFluxBoundary_hu =
        args.array<int>("isAdvectiveFluxBoundary_hu");
    xt::pyarray<int> &isAdvectiveFluxBoundary_hv =
        args.array<int>("isAdvectiveFluxBoundary_hv");
    xt::pyarray<int> &isDiffusiveFluxBoundary_hu =
        args.array<int>("isDiffusiveFluxBoundary_hu");
    xt::pyarray<int> &isDiffusiveFluxBoundary_hv =
        args.array<int>("isDiffusiveFluxBoundary_hv");
    xt::pyarray<double> &ebqe_bc_h_ext = args.array<double>("ebqe_bc_h_ext");
    xt::pyarray<double> &ebqe_bc_flux_mass_ext =
        args.array<double>("ebqe_bc_flux_mass_ext");
    xt::pyarray<double> &ebqe_bc_flux_mom_hu_adv_ext =
        args.array<double>("ebqe_bc_flux_mom_hu_adv_ext");
    xt::pyarray<double> &ebqe_bc_flux_mom_hv_adv_ext =
        args.array<double>("ebqe_bc_flux_mom_hv_adv_ext");
    xt::pyarray<double> &ebqe_bc_hu_ext = args.array<double>("ebqe_bc_hu_ext");
    xt::pyarray<double> &ebqe_bc_flux_hu_diff_ext =
        args.array<double>("ebqe_bc_flux_hu_diff_ext");
    xt::pyarray<double> &ebqe_penalty_ext =
        args.array<double>("ebqe_penalty_ext");
    xt::pyarray<double> &ebqe_bc_hv_ext = args.array<double>("ebqe_bc_hv_ext");
    xt::pyarray<double> &ebqe_bc_flux_hv_diff_ext =
        args.array<double>("ebqe_bc_flux_hv_diff_ext");
    xt::pyarray<int> &csrColumnOffsets_eb_h_h =
        args.array<int>("csrColumnOffsets_eb_h_h");
    xt::pyarray<int> &csrColumnOffsets_eb_h_hu =
        args.array<int>("csrColumnOffsets_eb_h_hu");
    xt::pyarray<int> &csrColumnOffsets_eb_h_hv =
        args.array<int>("csrColumnOffsets_eb_h_hv");
    xt::pyarray<int> &csrColumnOffsets_eb_hu_h =
        args.array<int>("csrColumnOffsets_eb_hu_h");
    xt::pyarray<int> &csrColumnOffsets_eb_hu_hu =
        args.array<int>("csrColumnOffsets_eb_hu_hu");
    xt::pyarray<int> &csrColumnOffsets_eb_hu_hv =
        args.array<int>("csrColumnOffsets_eb_hu_hv");
    xt::pyarray<int> &csrColumnOffsets_eb_hv_h =
        args.array<int>("csrColumnOffsets_eb_hv_h");
    xt::pyarray<int> &csrColumnOffsets_eb_hv_hu =
        args.array<int>("csrColumnOffsets_eb_hv_hu");
    xt::pyarray<int> &csrColumnOffsets_eb_hv_hv =
        args.array<int>("csrColumnOffsets_eb_hv_hv");
    double dt = args.scalar<double>("dt");
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
          elementJacobian_hw_hw[nDOF_test_element][nDOF_trial_element],
          elementJacobian_hbeta_hbeta[nDOF_test_element][nDOF_trial_element];
      for (int i = 0; i < nDOF_test_element; i++)
        for (int j = 0; j < nDOF_trial_element; j++) {
          elementJacobian_h_h[i][j] = 0.0;
          elementJacobian_hu_hu[i][j] = 0.0;
          elementJacobian_hv_hv[i][j] = 0.0;
          elementJacobian_heta_heta[i][j] = 0.0;
          elementJacobian_hw_hw[i][j] = 0.0;
          elementJacobian_hbeta_hbeta[i][j] = 0.0;
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
        ck.calculateMapping_element(
            eN, k, mesh_dof.data(), mesh_l2g.data(), mesh_trial_ref.data(),
            mesh_grad_trial_ref.data(), jac, jacDet, jacInv, x, y);
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
            elementJacobian_hbeta_hbeta[i][j] +=
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
          globalJacobian[csrRowIndeces_hbeta_hbeta[eN_i] +
                         csrColumnOffsets_hbeta_hbeta[eN_i_j]] +=
              elementJacobian_hbeta_hbeta[i][j];
        } // j
      }   // i
    }     // elements
  }

  void calculateLumpedMassMatrix(arguments_dict &args) {
    xt::pyarray<double> &mesh_trial_ref = args.array<double>("mesh_trial_ref");
    xt::pyarray<double> &mesh_grad_trial_ref =
        args.array<double>("mesh_grad_trial_ref");
    xt::pyarray<double> &mesh_dof = args.array<double>("mesh_dof");
    xt::pyarray<double> &mesh_velocity_dof =
        args.array<double>("mesh_velocity_dof");
    double MOVING_DOMAIN = args.scalar<double>("MOVING_DOMAIN");
    xt::pyarray<int> &mesh_l2g = args.array<int>("mesh_l2g");
    xt::pyarray<double> &dV_ref = args.array<double>("dV_ref");
    xt::pyarray<double> &h_trial_ref = args.array<double>("h_trial_ref");
    xt::pyarray<double> &h_grad_trial_ref =
        args.array<double>("h_grad_trial_ref");
    xt::pyarray<double> &h_test_ref = args.array<double>("h_test_ref");
    xt::pyarray<double> &h_grad_test_ref =
        args.array<double>("h_grad_test_ref");
    xt::pyarray<double> &vel_trial_ref = args.array<double>("vel_trial_ref");
    xt::pyarray<double> &vel_grad_trial_ref =
        args.array<double>("vel_grad_trial_ref");
    xt::pyarray<double> &vel_test_ref = args.array<double>("vel_test_ref");
    xt::pyarray<double> &vel_grad_test_ref =
        args.array<double>("vel_grad_test_ref");
    xt::pyarray<double> &mesh_trial_trace_ref =
        args.array<double>("mesh_trial_trace_ref");
    xt::pyarray<double> &mesh_grad_trial_trace_ref =
        args.array<double>("mesh_grad_trial_trace_ref");
    xt::pyarray<double> &dS_ref = args.array<double>("dS_ref");
    xt::pyarray<double> &h_trial_trace_ref =
        args.array<double>("h_trial_trace_ref");
    xt::pyarray<double> &h_grad_trial_trace_ref =
        args.array<double>("h_grad_trial_trace_ref");
    xt::pyarray<double> &h_test_trace_ref =
        args.array<double>("h_test_trace_ref");
    xt::pyarray<double> &h_grad_test_trace_ref =
        args.array<double>("h_grad_test_trace_ref");
    xt::pyarray<double> &vel_trial_trace_ref =
        args.array<double>("vel_trial_trace_ref");
    xt::pyarray<double> &vel_grad_trial_trace_ref =
        args.array<double>("vel_grad_trial_trace_ref");
    xt::pyarray<double> &vel_test_trace_ref =
        args.array<double>("vel_test_trace_ref");
    xt::pyarray<double> &vel_grad_test_trace_ref =
        args.array<double>("vel_grad_test_trace_ref");
    xt::pyarray<double> &normal_ref = args.array<double>("normal_ref");
    xt::pyarray<double> &boundaryJac_ref =
        args.array<double>("boundaryJac_ref");
    xt::pyarray<double> &elementDiameter =
        args.array<double>("elementDiameter");
    int nElements_global = args.scalar<int>("nElements_global");
    double g = args.scalar<double>("g");
    xt::pyarray<int> &h_l2g = args.array<int>("h_l2g");
    xt::pyarray<int> &vel_l2g = args.array<int>("vel_l2g");
    xt::pyarray<double> &b_dof = args.array<double>("b_dof");
    xt::pyarray<double> &h_dof = args.array<double>("h_dof");
    xt::pyarray<double> &hu_dof = args.array<double>("hu_dof");
    xt::pyarray<double> &hv_dof = args.array<double>("hv_dof");
    xt::pyarray<double> &q_cfl = args.array<double>("q_cfl");
    xt::pyarray<int> &sdInfo_hu_hu_rowptr =
        args.array<int>("sdInfo_hu_hu_rowptr");
    xt::pyarray<int> &sdInfo_hu_hu_colind =
        args.array<int>("sdInfo_hu_hu_colind");
    xt::pyarray<int> &sdInfo_hu_hv_rowptr =
        args.array<int>("sdInfo_hu_hv_rowptr");
    xt::pyarray<int> &sdInfo_hu_hv_colind =
        args.array<int>("sdInfo_hu_hv_colind");
    xt::pyarray<int> &sdInfo_hv_hv_rowptr =
        args.array<int>("sdInfo_hv_hv_rowptr");
    xt::pyarray<int> &sdInfo_hv_hv_colind =
        args.array<int>("sdInfo_hv_hv_colind");
    xt::pyarray<int> &sdInfo_hv_hu_rowptr =
        args.array<int>("sdInfo_hv_hu_rowptr");
    xt::pyarray<int> &sdInfo_hv_hu_colind =
        args.array<int>("sdInfo_hv_hu_colind");
    // h
    xt::pyarray<int> &csrRowIndeces_h_h = args.array<int>("csrRowIndeces_h_h");
    xt::pyarray<int> &csrColumnOffsets_h_h =
        args.array<int>("csrColumnOffsets_h_h");
    xt::pyarray<int> &csrRowIndeces_h_hu =
        args.array<int>("csrRowIndeces_h_hu");
    xt::pyarray<int> &csrColumnOffsets_h_hu =
        args.array<int>("csrColumnOffsets_h_hu");
    xt::pyarray<int> &csrRowIndeces_h_hv =
        args.array<int>("csrRowIndeces_h_hv");
    xt::pyarray<int> &csrColumnOffsets_h_hv =
        args.array<int>("csrColumnOffsets_h_hv");
    xt::pyarray<int> &csrRowIndeces_h_heta =
        args.array<int>("csrRowIndeces_h_heta");
    xt::pyarray<int> &csrColumnOffsets_h_heta =
        args.array<int>("csrColumnOffsets_h_heta");
    xt::pyarray<int> &csrRowIndeces_h_hw =
        args.array<int>("csrRowIndeces_h_hw");
    xt::pyarray<int> &csrColumnOffsets_h_hw =
        args.array<int>("csrColumnOffsets_h_hw");
    // hu
    xt::pyarray<int> &csrRowIndeces_hu_h =
        args.array<int>("csrRowIndeces_hu_h");
    xt::pyarray<int> &csrColumnOffsets_hu_h =
        args.array<int>("csrColumnOffsets_hu_h");
    xt::pyarray<int> &csrRowIndeces_hu_hu =
        args.array<int>("csrRowIndeces_hu_hu");
    xt::pyarray<int> &csrColumnOffsets_hu_hu =
        args.array<int>("csrColumnOffsets_hu_hu");
    xt::pyarray<int> &csrRowIndeces_hu_hv =
        args.array<int>("csrRowIndeces_hu_hv");
    xt::pyarray<int> &csrColumnOffsets_hu_hv =
        args.array<int>("csrColumnOffsets_hu_hv");
    xt::pyarray<int> &csrRowIndeces_hu_heta =
        args.array<int>("csrRowIndeces_hu_heta");
    xt::pyarray<int> &csrColumnOffsets_hu_heta =
        args.array<int>("csrColumnOffsets_hu_heta");
    xt::pyarray<int> &csrRowIndeces_hu_hw =
        args.array<int>("csrRowIndeces_hu_hw");
    xt::pyarray<int> &csrColumnOffsets_hu_hw =
        args.array<int>("csrColumnOffsets_hu_hw");
    // hv
    xt::pyarray<int> &csrRowIndeces_hv_h =
        args.array<int>("csrRowIndeces_hv_h");
    xt::pyarray<int> &csrColumnOffsets_hv_h =
        args.array<int>("csrColumnOffsets_hv_h");
    xt::pyarray<int> &csrRowIndeces_hv_hu =
        args.array<int>("csrRowIndeces_hv_hu");
    xt::pyarray<int> &csrColumnOffsets_hv_hu =
        args.array<int>("csrColumnOffsets_hv_hu");
    xt::pyarray<int> &csrRowIndeces_hv_hv =
        args.array<int>("csrRowIndeces_hv_hv");
    xt::pyarray<int> &csrColumnOffsets_hv_hv =
        args.array<int>("csrColumnOffsets_hv_hv");
    xt::pyarray<int> &csrRowIndeces_hv_heta =
        args.array<int>("csrRowIndeces_hv_heta");
    xt::pyarray<int> &csrColumnOffsets_hv_heta =
        args.array<int>("csrColumnOffsets_hv_heta");
    xt::pyarray<int> &csrRowIndeces_hv_hw =
        args.array<int>("csrRowIndeces_hv_hw");
    xt::pyarray<int> &csrColumnOffsets_hv_hw =
        args.array<int>("csrColumnOffsets_hv_hw");
    // heta
    xt::pyarray<int> &csrRowIndeces_heta_h =
        args.array<int>("csrRowIndeces_heta_h");
    xt::pyarray<int> &csrColumnOffsets_heta_h =
        args.array<int>("csrColumnOffsets_heta_h");
    xt::pyarray<int> &csrRowIndeces_heta_hu =
        args.array<int>("csrRowIndeces_heta_hu");
    xt::pyarray<int> &csrColumnOffsets_heta_hu =
        args.array<int>("csrColumnOffsets_heta_hu");
    xt::pyarray<int> &csrRowIndeces_heta_hv =
        args.array<int>("csrRowIndeces_heta_hv");
    xt::pyarray<int> &csrColumnOffsets_heta_hv =
        args.array<int>("csrColumnOffsets_heta_hv");
    xt::pyarray<int> &csrRowIndeces_heta_heta =
        args.array<int>("csrRowIndeces_heta_heta");
    xt::pyarray<int> &csrColumnOffsets_heta_heta =
        args.array<int>("csrColumnOffsets_heta_heta");
    xt::pyarray<int> &csrRowIndeces_heta_hw =
        args.array<int>("csrRowIndeces_heta_hw");
    xt::pyarray<int> &csrColumnOffsets_heta_hw =
        args.array<int>("csrColumnOffsets_heta_hw");
    // hw
    xt::pyarray<int> &csrRowIndeces_hw_h =
        args.array<int>("csrRowIndeces_hw_h");
    xt::pyarray<int> &csrColumnOffsets_hw_h =
        args.array<int>("csrColumnOffsets_hw_h");
    xt::pyarray<int> &csrRowIndeces_hw_hu =
        args.array<int>("csrRowIndeces_hw_hu");
    xt::pyarray<int> &csrColumnOffsets_hw_hu =
        args.array<int>("csrColumnOffsets_hw_hu");
    xt::pyarray<int> &csrRowIndeces_hw_hv =
        args.array<int>("csrRowIndeces_hw_hv");
    xt::pyarray<int> &csrColumnOffsets_hw_hv =
        args.array<int>("csrColumnOffsets_hw_hv");
    xt::pyarray<int> &csrRowIndeces_hw_heta =
        args.array<int>("csrRowIndeces_hw_heta");
    xt::pyarray<int> &csrColumnOffsets_hw_heta =
        args.array<int>("csrColumnOffsets_hw_heta");
    xt::pyarray<int> &csrRowIndeces_hw_hw =
        args.array<int>("csrRowIndeces_hw_hw");
    xt::pyarray<int> &csrColumnOffsets_hw_hw =
        args.array<int>("csrColumnOffsets_hw_hw");
    xt::pyarray<int> &csrRowIndeces_hbeta_hbeta =
        args.array<int>("csrRowIndeces_hbeta_hbeta");
    xt::pyarray<int> &csrColumnOffsets_hbeta_hbeta =
        args.array<int>("csrColumnOffsets_hbeta_hbeta");
    xt::pyarray<double> &globalJacobian = args.array<double>("globalJacobian");
    int nExteriorElementBoundaries_global =
        args.scalar<int>("nExteriorElementBoundaries_global");
    xt::pyarray<int> &exteriorElementBoundariesArray =
        args.array<int>("exteriorElementBoundariesArray");
    xt::pyarray<int> &elementBoundaryElementsArray =
        args.array<int>("elementBoundaryElementsArray");
    xt::pyarray<int> &elementBoundaryLocalElementBoundariesArray =
        args.array<int>("elementBoundaryLocalElementBoundariesArray");
    xt::pyarray<int> &isDOFBoundary_h = args.array<int>("isDOFBoundary_h");
    xt::pyarray<int> &isDOFBoundary_hu = args.array<int>("isDOFBoundary_hu");
    xt::pyarray<int> &isDOFBoundary_hv = args.array<int>("isDOFBoundary_hv");
    xt::pyarray<int> &isAdvectiveFluxBoundary_h =
        args.array<int>("isAdvectiveFluxBoundary_h");
    xt::pyarray<int> &isAdvectiveFluxBoundary_hu =
        args.array<int>("isAdvectiveFluxBoundary_hu");
    xt::pyarray<int> &isAdvectiveFluxBoundary_hv =
        args.array<int>("isAdvectiveFluxBoundary_hv");
    xt::pyarray<int> &isDiffusiveFluxBoundary_hu =
        args.array<int>("isDiffusiveFluxBoundary_hu");
    xt::pyarray<int> &isDiffusiveFluxBoundary_hv =
        args.array<int>("isDiffusiveFluxBoundary_hv");
    xt::pyarray<double> &ebqe_bc_h_ext = args.array<double>("ebqe_bc_h_ext");
    xt::pyarray<double> &ebqe_bc_flux_mass_ext =
        args.array<double>("ebqe_bc_flux_mass_ext");
    xt::pyarray<double> &ebqe_bc_flux_mom_hu_adv_ext =
        args.array<double>("ebqe_bc_flux_mom_hu_adv_ext");
    xt::pyarray<double> &ebqe_bc_flux_mom_hv_adv_ext =
        args.array<double>("ebqe_bc_flux_mom_hv_adv_ext");
    xt::pyarray<double> &ebqe_bc_hu_ext = args.array<double>("ebqe_bc_hu_ext");
    xt::pyarray<double> &ebqe_bc_flux_hu_diff_ext =
        args.array<double>("ebqe_bc_flux_hu_diff_ext");
    xt::pyarray<double> &ebqe_penalty_ext =
        args.array<double>("ebqe_penalty_ext");
    xt::pyarray<double> &ebqe_bc_hv_ext = args.array<double>("ebqe_bc_hv_ext");
    xt::pyarray<double> &ebqe_bc_flux_hv_diff_ext =
        args.array<double>("ebqe_bc_flux_hv_diff_ext");
    xt::pyarray<int> &csrColumnOffsets_eb_h_h =
        args.array<int>("csrColumnOffsets_eb_h_h");
    xt::pyarray<int> &csrColumnOffsets_eb_h_hu =
        args.array<int>("csrColumnOffsets_eb_h_hu");
    xt::pyarray<int> &csrColumnOffsets_eb_h_hv =
        args.array<int>("csrColumnOffsets_eb_h_hv");
    xt::pyarray<int> &csrColumnOffsets_eb_hu_h =
        args.array<int>("csrColumnOffsets_eb_hu_h");
    xt::pyarray<int> &csrColumnOffsets_eb_hu_hu =
        args.array<int>("csrColumnOffsets_eb_hu_hu");
    xt::pyarray<int> &csrColumnOffsets_eb_hu_hv =
        args.array<int>("csrColumnOffsets_eb_hu_hv");
    xt::pyarray<int> &csrColumnOffsets_eb_hv_h =
        args.array<int>("csrColumnOffsets_eb_hv_h");
    xt::pyarray<int> &csrColumnOffsets_eb_hv_hu =
        args.array<int>("csrColumnOffsets_eb_hv_hu");
    xt::pyarray<int> &csrColumnOffsets_eb_hv_hv =
        args.array<int>("csrColumnOffsets_eb_hv_hv");
    double dt = args.scalar<double>("dt");
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
          elementJacobian_hw_hw[nDOF_test_element][nDOF_trial_element],
          elementJacobian_hbeta_hbeta[nDOF_test_element][nDOF_trial_element];
      for (int i = 0; i < nDOF_test_element; i++)
        for (int j = 0; j < nDOF_trial_element; j++) {
          elementJacobian_h_h[i][j] = 0.0;
          elementJacobian_hu_hu[i][j] = 0.0;
          elementJacobian_hv_hv[i][j] = 0.0;
          elementJacobian_heta_heta[i][j] = 0.0;
          elementJacobian_hw_hw[i][j] = 0.0;
          elementJacobian_hbeta_hbeta[i][j] = 0.0;
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
        ck.calculateMapping_element(
            eN, k, mesh_dof.data(), mesh_l2g.data(), mesh_trial_ref.data(),
            mesh_grad_trial_ref.data(), jac, jacDet, jacInv, x, y);
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
            elementJacobian_hbeta_hbeta[i][j] +=
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
          globalJacobian[csrRowIndeces_hbeta_hbeta[eN_i] +
                         csrColumnOffsets_hbeta_hbeta[eN_i_j]] +=
              elementJacobian_hbeta_hbeta[i][j];
        } // j
      }   // i
    }     // elements
  }
}; // namespace proteus

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
} // end namespace proteus

#endif

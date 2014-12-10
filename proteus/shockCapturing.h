#ifndef SHOCKCAPTURING_H
#define SHOCKCAPTURING_H
/*!
 \file shockCapturing.h
 \brief C implementations of shock capturing diffusion calculations.
*/

/**
   \defgroup shockCapturing shockCapturing
 \brief C implementations of shock capturing diffusion calculations
 @{
*/
extern void calculateNumericalDiffusionResGrad(int nElements_global,
                                               int nQuadraturePoints_element,
                                               int nSpace,
                                               double shockCapturingDiffusion,
                                               double* elementDiameter,
                                               double* strong_residual,
                                               double* grad_u,
                                               double* numDiff);
extern void calculateNumericalDiffusionResGradQuad(int nElements_global,
                                                   int nQuadraturePoints_element,
                                                   int nSpace,
                                                   double shockCapturingDiffusion,
                                                   double *elementDiameter,
                                                   double *strong_residual,
                                                   double *grad_u,
                                                   double *numDiff);
extern void calculateNumericalDiffusionHJ(int nElements_global,
                                          int nQuadraturePoints_element,
                                          char shockCapturing,
                                          double shockCapturingDiffusion,
                                          double* elementDiameter,
                                          double* strong_residual,
                                          double* mt,
                                          double* H,
                                          double* numDiff);
extern void calculateNumericalDiffusionHJV2(int nElements_global,
                                            int nQuadraturePoints_element,
                                            char shockCapturing,
                                            double shockCapturingDiffusion,
                                            double* elementDiameter,
                                            double* strong_residual,
                                            double* mt,
                                            double* H,
                                            double* numDiff);
extern void calculateNumericalDiffusion_A_1(int nElements_global,
                                            int nQuadraturePoints_element,
                                            int nSpace,
                                            double shockCapturingFactor,
                                            double* elementDiameter,
                                            double* strong_residual,
                                            double* mt,
                                            double* df,
                                            double* numDiff);
extern void calculateNumericalDiffusionResGradJuanes(int nElements_global,
                                                     int nQuadraturePoints_element,
                                                     int nSpace,
                                                     double shockCapturingDiffusion,
                                                     double uSC,
                                                     double *elementDiameter,
                                                     double *strong_residual,
                                                     double *grad_u,
                                                     double *numDiff);
/** @} */
#endif

#ifndef FEMINTEGRALS_H
#define FEMINTEGRALS_H

/*!
 \file femIntegrals.h
 \brief C implementation of fem integral calculations
*/

/**
   \defgroup femIntegrals femIntegrals
 \brief C implementation of fem integral calculations
   @{ 
*/

extern void parametricFiniteElementSpace_getHessianValues(int nElements_global, int nQuadraturePoints_element, int nDOF_element, int nSpace_global, double *Hessian_psi, double *inverseJacobianArray, double *Hessian_vArray);
extern void     updateDiffusion2_strong(int nElements_global, int nQuadraturePoints_element, int nSpace, double *a, double *Hess_phi, double *strong_residual);
extern void updateDiffusion2_strong_sd(int nElements_global, int nQuadraturePoints_element,int nSpace, int* rowptr, int* colind, double* a, double* Hess_phi, double* strong_residual);
extern void     updateDiffusionJacobian2_strong(int nElements_global, int nQuadraturePoints_element, int nDOF_trial_element, int nSpace, int *l2g, double *a, double *da, double *v, double *Hess_phi, double *dphi, double *Hess_v, double *dstrong_residual);
extern void updateDiffusionJacobian2_strong_sd(int nElements_global, 
					       int nQuadraturePoints_element,
					       int nDOF_trial_element,
					       int nSpace,
					       int* rowptr,
					       int* colind,
					       int* l2g,
					       double* a,
					       double* da,
					       double* v,
					       double* Hess_phi,
					       double* dphi,
					       double* Hess_v,
					       double* dstrong_residual);
extern void     updateDiffusion2_adjoint(int nElements_global, int nQuadraturePoints_element, int nDOF_test_element, int nSpace, double *a, double *Hess_w_dV, double *Lstar_w_dV);
void updateDiffusion2_adjoint_sd(int nElements_global, int nQuadraturePoints_element, int nDOF_test_element, int nSpace, int* rowptr, int* colind, double* a, double* Hess_w_dV, double* Lstar_w_dV);
extern void     calculateWeightedShapeHessians(int nElements_global, int nQuadraturePoints_element, int nDOF_test_element, int nSpace, double *dVR, double *abs_det_jac, double *Hess_w, double *Hess_w_dV);
extern void     calculateFiniteElementFunctionHessianValues(int nElements_global, int nQuadraturePoints_element, int nDOF_trial_element, int nComponents, int nSpace, int *l2g, double *dof, double *Hessian_v, double *Hessian_u);
extern void     updateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_2sided_CSR(int nInteriorElementBoundaries_global, int nElementBoundaries_element, int nQuadraturePoints_elementBoundary, int nDOF_test_element, int nDOF_trial_element, int *interiorElementBoundaries, int *elementBoundaryElements, int *elementBoundaryLocalElementBoundaries, int *nFreeDOF_element_r, int *freeLocal_r, int *nFreeDOF_element_u, int *freeLocal_u, int *csrRowIndeces_ru, int *csrColumnOffsets_eb_ru, double *elementBoundaryFluxJacobian_2sided, double *w_dS, double *jac);
extern void     updateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_2sided_dense(int nInteriorElementBoundaries_global, int nElementBoundaries_element, int nQuadraturePoints_elementBoundary, int nDOF_test_element, int nDOF_trial_element, int offset_r, int stride_r, int offset_u, int stride_u, int nFreeVDOF_global, int *interiorElementBoundaries, int *elementBoundaryElements, int *elementBoundaryLocalElementBoundaries, int *nFreeDOF_element_r, int *freeLocal_r, int *freeGlobal_r, int *nFreeDOF_element_u, int *freeLocal_u, int *freeGlobal_u, double *elementBoundaryFluxJacobian_2sided, double *w_dS, double *jac);
extern void     updateInteriorTwoSidedElementBoundaryFlux(int nInteriorElementBoundaries_global, int nElementBoundaries_element, int nQuadraturePoints_elementBoundary, int nDOF_test_element, int *interiorElementBoundaries, int *elementBoundaryElements, int *elementBoundaryLocalElementBoundaries, double *flux, double *w_dS, double *residual);
extern void     calculateCFLADR2speeds(int nElements_global, int nQuadraturePoints_element, int nSpace, double *elementDiameter, double *dm, double *df1, double *df2, double *cfl);
extern int      checkElementBoundaryAndExteriorElementBoundaryArraysSame(int nElementBoundaries_element, int nExteriorElementBoundaries_global, int nQuadraturePoints_elementBoundary, int nValuesPerQuadraturePoint, double tolerance, const int *exteriorElementBoundariesArray, const int *elementBoundaryElementsArray, const int *elementBoundaryLocalElementBoundariesArray, const double *ebq_val, const double *ebqe_val, int *firstBadIndex);
extern int      checkGlobalElementBoundaryAndExteriorElementBoundaryArraysSame(int nExteriorElementBoundaries_global, int nQuadraturePoints_elementBoundary, int nValuesPerQuadraturePoint, double tolerance, const int *exteriorElementBoundariesArray, const int *elementBoundaryElementsArray, const int *elementBoundaryLocalElementBoundariesArray, const double *ebq_global_val, const double *ebqe_val, int *firstBadIndex);
extern void     calculateExteriorElementBoundaryStress3D(int nExteriorElementBoundaries_global, int nQuadraturePoints_elementBoundary, int *elementBoundaryMaterialTypes, int *exteriorElementBoundaries, int *elementBoundaryElements, int *elementBoundaryLocalElementBoundaries, double *p, double *mom_flux_vec_u, double *mom_flux_vec_v, double *mom_flux_vec_w, double *dS, double *n, double *F);
extern void     calculateExteriorElementBoundaryStress2D(int nExteriorElementBoundaries_global, int nQuadraturePoints_elementBoundary, int *elementBoundaryMaterialTypes, int *exteriorElementBoundaries, int *elementBoundaryElements, int *elementBoundaryLocalElementBoundaries, double *p, double *mom_flux_vec_u, double *mom_flux_vec_v, double *dS, double *n, double *F);

extern void copyLeftElementBoundaryInfo(
  int nElementBoundaries_element,
  int nElementBoundaryQuadraturePoints_elementBoundary,
  int nSpace_global,
  int nExteriorElementBoundaries_global,
  int nInteriorElementBoundaries_global,
  int *elementBoundaryElementsArray,
  int *elementBoundaryLocalElementBoundariesArray,
  int *exteriorElementBoundariesArray,
  int *interiorElementBoundariesArray,
  double *x,
  double *n,
  double *xg,
  double *ng
);
extern void parametricFiniteElementSpace_getValues(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_element,
  double *psi,
  double *vArray
);
extern void parametricFiniteElementSpace_getValuesTrace(
  int nElements_global,
  int nElementBoundaries_element,
  int nElementBoundaryQuadraturePoints_elementBoundary,
  int nDOF_element,
  double *psi,
  int *permutations,
  double *vArray
);
extern void parametricFiniteElementSpace_getGradientValues(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_element,
  int nSpace_global,
  double *grad_psi,
  double *inverseJacobianArray,
  double *grad_vArray
);
extern void parametricFiniteElementSpace_getGradientValuesTrace(
  int nElements_global,
  int nElementBoundaries_element,
  int nElementBoundaryQuadraturePoints_elementBoundary,
  int nDOF_element,
  int nSpace_global,
  double *grad_psi,
  int *permutations,
  double *inverseJacobianArray,
  double *grad_vArray
);
extern void parametricMaps_getPermutations(
  int nElements_global,
  int nElementBoundaries_element,
  int nElementBoundaryQuadraturePoints_elementBoundary,
  int nSpace_global,
  double *xiArray,
  int *permutations
);
extern void parametricMaps_getPermutationsGlobalExterior(int nElementBoundaryQuadraturePoints_elementBoundary,
							 int nSpace_global,
							 int nExteriorElementBoundaries_global,
							 const int * exteriorElementBoundariesArray,
							 const int * elementBoundaryElementsArray,
							 const int * elementBoundaryLocalElementBoundariesArray,
							 double* xiArray,
							 int* permutations);
extern void getPermutationsGlobal(int nElementBoundaries_global,
				  int nElementBoundaryQuadraturePoints_elementBoundary,
				  double* xArray,
				  double* xArrayNew,
				  int* permutations);
extern void parametricMaps_getValues(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_element,
  int nSpace_global,
  double *psi,
  int *l2g,
  double *nodeArray,
  double *xArray
);
extern void parametricMaps_getValuesTrace(
  int nElements_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_element,
  int nDOF_element,
  int nSpace_global,
  double *psi,
  int *l2g,
  double *nodeArray,
  double *xArray
);
extern void parametricMaps_getInverseValues(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_element,
  int nSpace_global,
  double *inverseJacobian,
  int *l2g,
  double *nodeArray,
  double *xArray,
  double *xiArray
);
extern void parametricMaps_getInverseValuesTrace(
  int nElements_global,
  int nElementBoundaries_element,
  int nElementBoundaryQuadraturePoints_elementBoundary,
  int nDOF_element,
  int nSpace_global,
  double *inverseJacobian,
  int *l2g,
  double *nodeArray,
  double *xArray,
  double *xiArray
);
extern void parametricMaps_getJacobianValues3D(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_element,
  double *grad_psi,
  int *l2g,
  double *nodeArray,
  double *jacobianArray,
  double *jacobianDeterminantArray,
  double *jacobianInverseArray
);
extern void parametricMaps_getJacobianValues2D(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_element,
  double *grad_psi,
  int *l2g,
  double *nodeArray,
  double *jacobianArray,
  double *jacobianDeterminantArray,
  double *jacobianInverseArray
);
extern void parametricMaps_getJacobianValues1D(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_element,
  double *grad_psi,
  int *l2g,
  double *nodeArray,
  double *jacobianArray,
  double *jacobianDeterminantArray,
  double *jacobianInverseArray
);
extern void parametricMaps_getJacobianValuesTrace3D(
  int nElements_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_element,
  int nDOF_element,
  double *grad_psi,
  double *boundaryNormals,
  double *boundaryJacobians,
  int *l2g,
  double *nodeArray,
  double *jacobianInverseArray,
  double *metricTensorArray,
  double *metricTensorDeterminantSqrtArray,
  double *unitNormalArray
);
extern void parametricMaps_getJacobianValuesTrace2D(
  int nElements_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_element,
  int nDOF_element,
  double *grad_psi,
  double *boundaryNormals,
  double *boundaryJacobians,
  int *l2g,
  double *nodeArray,
  double *jacobianInverseArray,
  double *metricTensorArray,
  double *metricTensorDeterminantSqrtArray,
  double *unitNormalArray
);
extern void parametricMaps_getJacobianValuesTrace1D(
  int nElements_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_element,
  int nDOF_element,
  double *grad_psi,
  double *boundaryNormals,
  double *boundaryJacobians,
  int *l2g,
  double *nodeArray,
  double *jacobianInverseArray,
  double *metricTensorArray,
  double *metricTensorDeterminantSqrtArray,
  double *unitNormalArray
);
extern void updateMass_weak(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_test_element,
  double *mt,
  double *w_dV,
  double *weak_residual
);
extern void updateMassJacobian_weak(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_trial_element,
  int nDOF_test_element,
  double *dmt,
  double *v_X_w_dV,
  double *jacobian_weak_residual
);
extern void updateMassJacobian_weak_lowmem(int nElements_global,
                                           int nQuadraturePoints_element,
                                           int nDOF_trial_element,
                                           int nDOF_test_element,
                                           double* dmt,
                                           double* v,
                                           double* w_dV,
                                           double* jacobian_weak_residual);
extern void updateMass_strong(
  int nElements_global,
  int nQuadraturePoints_element,
  double *mt,
  double *strong_residual
);
extern void updateMassJacobian_strong(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_trial_element,
  double *dmt,
  double *v,
  double *dstrong_residual
);
extern void updateMass_adjoint(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_test_element,
  double *dmt,
  double *w_dV,
  double *Lstar_w_dV
);
extern void updateAdvection_weak(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_test_element,
  int nSpace,
  double *f,
  double *grad_w_dV,
  double *weak_residual
);
extern void updateAdvectionJacobian_weak(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_trial_element,
  int nDOF_test_element,
  int nSpace,
  double *df,
  double *v_X_grad_w_dV,
  double *jacobian_weak_residual
);
extern void updateAdvectionJacobian_weak_lowmem(int nElements_global,
                                                int nQuadraturePoints_element,
                                                int nDOF_trial_element,
                                                int nDOF_test_element,
                                                int nSpace,
                                                double* df,
                                                double* v,
                                                double* grad_w_dV,
                                                double* jacobian_weak_residual);
extern void updateAdvection_strong(
  int nElements_global,
  int nQuadraturePoints_element,
  int nSpace,
  double *df,
  double *grad_u,
  double *strong_residual
);
extern void updateAdvectionJacobian_strong(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_trial_element,
  int nSpace,
  double *df,
  double *grad_v,
  double *dstrong_residual
);
extern void updateAdvection_adjoint(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_test_element,
  int nSpace,
  double *df,
  double *grad_w_dV,
  double *Lstar_w_dV
);
extern void updateHamiltonian_weak(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_test_element,
  double *H,
  double *w_dV,
  double *weak_residual
);
extern void updateHamiltonianJacobian_weak(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_trial_element,
  int nDOF_test_element,
  int nSpace,
  double *dH,
  double *grad_v_X_w_dV,
  double *jacobian_weak_residual
);
extern void updateHamiltonianJacobian_weak_lowmem(int nElements_global,
                                                  int nQuadraturePoints_element,
                                                  int nDOF_trial_element,
                                                  int nDOF_test_element,
                                                  int nSpace,
                                                  double* dH,
                                                  double* grad_v,
                                                  double* w_dV,
                                                  double* jacobian_weak_residual);
extern void updateHamiltonian_strong(
  int nElements_global,
  int nQuadraturePoints_element,
  int nSpace,
  double *dH,
  double *grad_u,
  double *strong_residual
);
extern void updateHamiltonianJacobian_strong(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_trial_element,
  int nSpace,
  double *dH,
  double *grad_v,
  double *dstrong_residual
);
extern void updateHamiltonian_adjoint(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_test_element,
  int nSpace,
  double *dH,
  double *grad_w_dV,
  double *Lstar_w_dV
);
extern void updateDiffusion_weak(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_test_element,
  int nSpace,
  double *a,
  double *grad_phi_X_grad_w_dV,
  double *weak_residual
);
extern void updateDiffusion_weak_lowmem(int nElements_global,
                                 int nQuadraturePoints_element,
                                 int nDOF_test_element,
                                 int nSpace,
                                 double* a,
                                 double* grad_phi,
                                 double* grad_w_dV,
                                 double* weak_residual);
extern void updateDiffusion_weak_sd(int nElements_global,
				    int nQuadraturePoints_element,
				    int nDOF_test_element,
				    int nSpace,
				    int* rowptr,
				    int* colind,
				    double* a,
				    double* grad_phi,
				    double* grad_w_dV,
				    double* weak_residual);
extern void updateDiffusionJacobian_weak(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_trial_element,
  int nDOF_test_element,
  int nSpace,
  int *l2g,
  double *a,
  double *da,
  double *grad_phi_X_grad_w_dV,
  double *dphi,
  double *v,
  double *grad_v_X_grad_w_dV,
  double *jacobian_weak_residual
);
extern void updateDiffusionJacobian_weak_lowmem(int nElements_global,
                                                int nQuadraturePoints_element,
                                                int nDOF_trial_element,
                                                int nDOF_test_element,
                                                int nSpace,
                                                int* l2g,
                                                double* a,
                                                double* da,
                                                double* grad_phi,
                                                double* grad_w_dV,
                                                double* dphi,
                                                double* v,
                                                double* grad_v,
                                                double* jacobian_weak_residual);
extern void updateDiffusionJacobian_weak_sd(int nElements_global,
					    int nQuadraturePoints_element,
					    int nDOF_trial_element,
					    int nDOF_test_element,
					    int nSpace,
					    int* rowptr,
					    int* colind,
					    int* l2g,
					    double* a,
					    double* da,
					    double* grad_phi,
					    double* grad_w_dV,
					    double* dphi,
					    double* v,
					    double* grad_v,
					    double* jacobian_weak_residual);
extern void updateDiffusion_strong(
  int nElements_global,
  int nQuadraturePoints_element,
  int nSpace,
  double *da,
  double *grad_phi,
  double *grad_u,
  double *strong_residual
);
extern void updateDiffusion_strong_sd(int nElements_global,
				      int nQuadraturePoints_element,
				      int nSpace,
				      int* rowptr,
				      int* colind,
				      double* da,
				      double* grad_phi,
				      double* grad_u,
				      double* strong_residual);
extern void updateDiffusionJacobian_strong(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_trial_element,
  int nSpace,
  int *l2g,
  double *da,
  double *dphi,
  double *grad_phi,
  double *grad_u,
  double *grad_v,
  double *dstrong_residual
);
extern void updateDiffusionJacobian_strong_sd(int nElements_global,
					      int nQuadraturePoints_element,
					      int nDOF_trial_element,
					      int nSpace,
					      int* rowptr,
					      int* colind,
					      int* l2g,
					      double* da,
					      double* dphi,
					      double* grad_phi,
					      double* grad_u,
					      double* grad_v,
					      double* dstrong_residual);
extern void updateDiffusion_adjoint(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_test_element,
  int nSpace,
  double *da,
  double *grad_phi,
  double *grad_w_dV,
  double *Lstar_w_dV
);
extern void updateDiffusion_adjoint_sd(int nElements_global,
				       int nQuadraturePoints_element,
				       int nDOF_test_element,
				       int nSpace,
				       int* rowptr,
				       int* colind,
				       double* da,
				       double* grad_phi,
				       double* grad_w_dV,
				       double* Lstar_w_dV);
extern void updateReaction_weak(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_test_element,
  double *r,
  double *w_dV,
  double *weak_residual
);
extern void updateReactionJacobian_weak(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_trial_element,
  int nDOF_test_element,
  double *dr,
  double *v_X_w_dV,
  double *jacobian_weak_residual
);
extern void updateReactionJacobian_weak_lowmem(int nElements_global,
                                               int nQuadraturePoints_element,
                                               int nDOF_trial_element,
                                               int nDOF_test_element,
                                               double* dr,
                                               double* v,
                                               double* w_dV,
                                               double* jacobian_weak_residual);
extern void updateReaction_strong(
  int nElements_global,
  int nQuadraturePoints_element,
  double *r,
  double *strong_residual
);
extern void updateReactionJacobian_strong(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_trial_element,
  double *dr,
  double *v,
  double *dstrong_residual
);
extern void updateReaction_adjoint(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_test_element,
  double *dr,
  double *w_dV,
  double *Lstar_w_dV
);
extern void updateSubgridError(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_test_element,
  double *error,
  double *Lstar_w_dV,
  double *weak_residual
);
extern void updateSubgridErrorJacobian(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_trial_element,
  int nDOF_test_element,
  double *derror,
  double *Lstar_w_dV,
  double *jacobian_weak_residual
);
extern void updateNumericalDiffusion(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_test_element,
  int nSpace,
  double *numDiff,
  double *grad_u_X_grad_w_dV,
  double *weak_residual
);
extern void updateNumericalDiffusion_lowmem(int nElements_global,
                                     int nQuadraturePoints_element,
                                     int nDOF_test_element,
                                     int nSpace,
                                     double* numDiff,
                                     double* grad_u,
                                     double* grad_w_dV,
                                     double* weak_residual);
extern void updateNumericalDiffusionJacobian(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_trial_element,
  int nDOF_test_element,
  int nSpace,
  double *numDiff,
  double *grad_v_X_grad_w_dV,
  double *jacobian_weak_residual
);
extern void updateNumericalDiffusionJacobian_lowmem(int nElements_global,
                                                    int nQuadraturePoints_element,
                                                    int nDOF_trial_element,
                                                    int nDOF_test_element,
                                                    int nSpace,
                                                    double* numDiff,
                                                    double* grad_v,
                                                    double* grad_w_dV,
                                                    double* jacobian_weak_residual);
extern void calculateScalarScalarProduct(
  int nElements_global,
  int nQuadraturePoints_element,
  double *s1,
  double *s2,
  double *sResult
);
extern void calculateVectorScalarProduct(
  int nElements_global,
  int nQuadraturePoints_element,
  int nSpace,
  double *v,
  double *s,
  double *vResult
);
extern void calculateTensorScalarProduct(
  int nElements_global,
  int nQuadraturePoints_element,
  int nSpace,
  double *t,
  double *s,
  double *tResult
);
extern void updateInteriorElementBoundaryFlux(
  int nInteriorElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nDOF_test_element,
  int *interiorElementBoundaries,
  int *elementBoundaryElements,
  int *elementBoundaryLocalElementBoundaries,
  double *flux,
  double *w_dS,
  double *residual
);
extern void updateExteriorElementBoundaryFlux(
  int nExteriorElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nDOF_test_element,
  int *exteriorElementBoundaries,
  int *elementBoundaryElements,
  int *elementBoundaryLocalElementBoundaries,
  double *flux,
  double *w_dS,
  double *residual
);
extern void updateGlobalResidualFromElementResidual(
  int nElements_global,
  int nDOF_test_element,
  int offset_r,
  int stride_r,
  int *nFreeDOF_element_r,
  int *freeLocal_r,
  int *freeGlobal_r,
  double *elementResidual,
  double *globalResidual
);
extern void updateGlobalJacobianFromElementJacobian_dense(
  int nElements_global,
  int nDOF_test_element,
  int nDOF_trial_element,
  int offset_r,
  int stride_r,
  int offset_u,
  int stride_u,
  int nFreeVDOF_global,
  int *nFreeDOF_element_r,
  int *freeLocal_r,
  int *freeGlobal_r,
  int *nFreeDOF_element_u,
  int *freeLocal_u,
  int *freeGlobal_u,
  double *elementJacobian,
  double *globalJacobian
);
extern void updateGlobalJacobianFromElementJacobian_eb_dense(
  int *elementNeighbors,
  int nElements_global,
  int nElementBoundaries_element,
  int nDOF_test_element,
  int nDOF_trial_element,
  int offset_r,
  int stride_r,
  int offset_u,
  int stride_u,
  int nFreeVDOF_global,
  int *nFreeDOF_element_r,
  int *freeLocal_r,
  int *freeGlobal_r,
  int *nFreeDOF_element_u,
  int *freeLocal_u,
  int *freeGlobal_u,
  double *elementJacobian_eb,
  double *globalJacobian
);
extern void updateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_dense(
  int nInteriorElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nDOF_test_element,
  int nDOF_trial_element,
  int offset_r,
  int stride_r,
  int offset_u,
  int stride_u,
  int nFreeVDOF_global,
  int *interiorElementBoundaries,
  int *elementBoundaryElements,
  int *elementBoundaryLocalElementBoundaries,
  int *nFreeDOF_element_r,
  int *freeLocal_r,
  int *freeGlobal_r,
  int *nFreeDOF_element_u,
  int *freeLocal_u,
  int *freeGlobal_u,
  double *elementBoundaryFluxJacobian,
  double *w_dS,
  double *jac
);
extern void
updateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_eb_dense(
  int *elementNeighbors,
  int nElements_global,
  int nInteriorElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nDOF_test_element,
  int nDOF_trial_element,
  int offset_r,
  int stride_r,
  int offset_u,
  int stride_u,
  int nFreeVDOF_global,
  int *interiorElementBoundaries,
  int *elementBoundaryElements,
  int *elementBoundaryLocalElementBoundaries,
  int *nFreeDOF_element_r,
  int *freeLocal_r,
  int *freeGlobal_r,
  int *nFreeDOF_element_u,
  int *freeLocal_u,
  int *freeGlobal_u,
  double *elementBoundaryFluxJacobian_eb,
  double *w_dS,
  double *jac
);
extern void updateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_dense(
  int nExteriorElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nDOF_test_element,
  int nDOF_trial_element,
  int offset_r,
  int stride_r,
  int offset_u,
  int stride_u,
  int nFreeVDOF_global,
  int *exteriorElementBoundaries,
  int *elementBoundaryElements,
  int *elementBoundaryLocalElementBoundaries,
  int *nFreeDOF_element_r,
  int *freeLocal_r,
  int *freeGlobal_r,
  int *nFreeDOF_element_u,
  int *freeLocal_u,
  int *freeGlobal_u,
  double *elementBoundaryFluxJacobian,
  double *w_dS,
  double *jac
);
extern void
updateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_eb_dense(
  int *elementNeighbors,
  int nElements_global,
  int nExteriorElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nDOF_test_element,
  int nDOF_trial_element,
  int offset_r,
  int stride_r,
  int offset_u,
  int stride_u,
  int nFreeVDOF_global,
  int *exteriorElementBoundaries,
  int *elementBoundaryElements,
  int *elementBoundaryLocalElementBoundaries,
  int *nFreeDOF_element_r,
  int *freeLocal_r,
  int *freeGlobal_r,
  int *nFreeDOF_element_u,
  int *freeLocal_u,
  int *freeGlobal_u,
  double *elementBoundaryFluxJacobian_eb,
  double *w_dS,
  double *jac
);
extern void updateGlobalJacobianFromElementJacobian_CSR(
  int nElements_global,
  int nDOF_test_element,
  int nDOF_trial_element,
  int *nFreeDOF_element_r,
  int *freeLocal_r,
  int *nFreeDOF_element_u,
  int *freeLocal_u,
  int *csrRowIndeces_ru,
  int *csrColumnOffsets_ru,
  double *elementJacobian,
  double *globalJacobian
);
extern void updateGlobalJacobianFromElementJacobian_eb_CSR(
  int *elementNeighbors,
  int nElements_global,
  int nElementBoundaries_element,
  int nDOF_test_element,
  int nDOF_trial_element,
  int *nFreeDOF_element_r,
  int *freeLocal_r,
  int *nFreeDOF_element_u,
  int *freeLocal_u,
  int *csrRowIndeces_ru,
  int *csrColumnOffsets_eb_ru,
  double *elementJacobian_eb,
  double *globalJacobian
);
extern void updateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_CSR(
  int nInteriorElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nDOF_test_element,
  int nDOF_trial_element,
  int *interiorElementBoundaries,
  int *elementBoundaryElements,
  int *elementBoundaryLocalElementBoundaries,
  int *nFreeDOF_element_r,
  int *freeLocal_r,
  int *nFreeDOF_element_u,
  int *freeLocal_u,
  int *csrRowIndeces_ru,
  int *csrColumnOffsets_eb_ru,
  double *elementBoundaryFluxJacobian,
  double *w_dS,
  double *jac
);
extern void updateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_CSR(
  int nExteriorElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nDOF_test_element,
  int nDOF_trial_element,
  int *exteriorElementBoundaries,
  int *elementBoundaryElements,
  int *elementBoundaryLocalElementBoundaries,
  int *nFreeDOF_element_r,
  int *freeLocal_r,
  int *nFreeDOF_element_u,
  int *freeLocal_u,
  int *csrRowIndeces_ru,
  int *csrColumnOffsets_eb_ru,
  double *elementBoundaryFluxJacobian,
  double *w_dS,
  double *jac
);
extern void updateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_eb_CSR(
  int *elementNeighbors,
  int nInteriorElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nDOF_test_element,
  int nDOF_trial_element,
  int *interiorElementBoundaries,
  int *elementBoundaryElements,
  int *elementBoundaryLocalElementBoundaries,
  int *nFreeDOF_element_r,
  int *freeLocal_r,
  int *nFreeDOF_element_u,
  int *freeLocal_u,
  int *csrRowIndeces_ru,
  int *csrColumnOffsets_eb_eNebN_ru,
  double *elementBoundaryFluxJacobian_eb,
  double *w_dS,
  double *jac
);
extern void updateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_eb_CSR(
  int *elementNeighbors,
  int nExteriorElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nDOF_test_element,
  int nDOF_trial_element,
  int *exteriorElementBoundaries,
  int *elementBoundaryElements,
  int *elementBoundaryLocalElementBoundaries,
  int *nFreeDOF_element_r,
  int *freeLocal_r,
  int *nFreeDOF_element_u,
  int *freeLocal_u,
  int *csrRowIndeces_ru,
  int *csrColumnOffsets_eb_eNebN_ru,
  double *elementBoundaryFluxJacobian_eb,
  double *w_dS,
  double *jac
);
extern void calculateWeightedShape(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_test_element,
  double *dVR,
  double *abs_det_jac,
  double *w,
  double *w_dV
);
extern void calculateWeightedShapeGradients(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_test_element,
  int nSpace,
  double *dVR,
  double *abs_det_jac,
  double *grad_w,
  double *grad_w_dV
);
extern void calculateShape_X_weightedShape(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_trial_element,
  int nDOF_test_element,
  double *v,
  double *w_dV,
  double *v_X_w_dV
);
extern void calculateShape_X_weightedGradShape(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_trial_element,
  int nDOF_test_element,
  int nSpace,
  double *v,
  double *grad_w_dV,
  double *v_X_grad_w_dV
);
extern void calculateGradShape_X_weightedShape(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_trial_element,
  int nDOF_test_element,
  int nSpace,
  double *grad_v,
  double *w_dV,
  double *grad_v_X_w_dV
);
extern void calculateGradShape_X_weightedGradShape(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_trial_element,
  int nDOF_test_element,
  int nSpace,
  double *grad_v,
  double *grad_w_dV,
  double *grad_v_X_grad_w_dV
);
extern void calculateWeightedShapeTrace(
  int nElements_global,
  int nElementBoundaries_element,
  int nElementBoundaryQuadraturePoints_elementBoundary,
  int nDOF_test_element,
  double *dSR,
  double *sqrt_det_g,
  double *w,
  double *w_dS
);
extern void calculateShape_X_weightedShapeTrace(
  int nElements_global,
  int nElementBoundaries_element,
  int nElementBoundaryQuadraturePoints_elementBoundary,
  int nDOF_trial_element,
  int nDOF_test_element,
  double *v,
  double *w_dS,
  double *v_X_w_dS
);
extern void calculateGradShape_X_weightedShapeTrace(
  int nElements_global,
  int nElementBoundaries_element,
  int nElementBoundaryQuadraturePoints_elementBoundary,
  int nDOF_trial_element,
  int nDOF_test_element,
  int nSpace,
  double *grad_v,
  double *w_dS,
  double *grad_v_X_w_dS
);
extern void calculateIntegrationWeights(
  int nElements_global,
  int nQuadraturePoints_element,
  double *abs_det_J,
  double *referenceWeights,
  double *weights
);
extern void calculateElementBoundaryIntegrationWeights(
  int nElements_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  double *sqrt_det_g,
  double *referenceWeights,
  double *weights
);
extern void calculateFiniteElementFunctionValues(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_trial_element,
  int nComponents,
  int *l2g,
  double *dof,
  double *v,
  double *u
);
extern void calculateFiniteElementFunctionGradientValues(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_trial_element,
  int nComponents,
  int nSpace,
  int *l2g,
  double *dof,
  double *grad_v,
  double *grad_u
);
extern void calculateFiniteElementFunctionGradientTensorValues(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_trial_element,
  int nDOF_test_element,
  int nComponents,
  int nSpace,
  int *l2g,
  double *dof,
  double *grad_v_X_grad_w_dV,
  double *grad_u_X_grad_w_dV
);
extern void calculateFiniteElementFunctionValuesTrace(
  int nElements_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nDOF_trial_element,
  int nComponents,
  int *l2g,
  double *dof,
  double *v,
  double *u
);
extern void calculateFiniteElementFunctionGradientValuesTrace(
  int nElements_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nDOF_trial_element,
  int nComponents,
  int nSpace,
  int *l2g,
  double *dof,
  double *grad_v,
  double *grad_u
);
extern void calculateFlowVelocity(
  int nElements_global,
  int nQuadraturePoints_element,
  int nSpace,
  double *f,
  double *a,
  double *grad_phi,
  double *v
);
extern void updateAddJacobian_CSR(
  int jacIndex,
  double val,
  double *jac
);
extern void zeroJacobian_CSR(
  int nNonzeros,
  double *jac
);
extern void setInflowFlux(
  int nExteriorElementBoundaries_global,
  int nQuadraturePoints_elementBoundary,
  int *exteriorElementBoundaries,
  double *inflowFlux,
  double *flux
);
extern void calculateInteriorElementBoundaryVelocities(
  int nInteriorElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nSpace,
  int *interiorElementBoundaries,
  int *elementBoundaryElements,
  int *elementBoundaryLocalElementBoundaries,
  double *m,
  double *a,
  double *grad_phi,
  double *f,
  double *vAverage,
  double *vJump,
  double *mAverage,
  double *mJump
);
extern void calculateExteriorElementBoundaryVelocities(
  int nExteriorElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nSpace,
  int *exteriorElementBoundaries,
  int *elementBoundaryElements,
  int *elementBoundaryLocalElementBoundaries,
  double *m,
  double *a,
  double *grad_phi,
  double *f,
  double *vAverage,
  double *vJump,
  double *mAverage,
  double *mJump
);
extern void calculateConservationResidualPWL(
  int nElements_global,
  int nInteriorElementBoundaries_global,
  int nExteriorElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nNodes_element,
  int nSpace,
  int *interiorElementBoundaries,
  int *exteriorElementBoundaries,
  int *elementBoundaryElements,
  int *elementBoundaryLocalElementBoundaries,
  int *elementNodes,
  int *nodeStarElements,
  int *nodeStarElementNeighbors,
  int *nodeStarOffsets,
  int *nElements_node,
  double *elementResidual,
  double *vAverage,
  double *starU,
  double *dX,
  double *w,
  double *normal,
  double *conservationResidual,
  double *starR,
  double *vConservative,
  double *vConservative_element
);
extern void calculateConservationJacobianPWL(
  int nNodes_global,
  int nNodes_internal,
  int nElements_global,
  int nInteriorElementBoundaries_global,
  int nExteriorElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nNodes_element,
  int nSpace,
  int *interiorElementBoundaries,
  int *exteriorElementBoundaries,
  int *elementBoundaryElements,
  int *elementBoundaryLocalElementBoundaries,
  int *elementNodes,
  int *nodeStarElements,
  int *nodeStarElementNeighbors,
  int *nodeStarOffsets,
  int *nodeStarJacobianOffsets,
  int *nElements_node,
  int *internalNodes,
  double *w,
  double *normal,
  double *starJacobian
);
extern void calculateConservationFluxPWL(
  int nNodes_global,
  int nNodes_internal,
  int *nElements_node,
  int *nodeStarOffsets,
  int *nodeStarJacobianOffsets,
  int *internalNodes,
  double *starR,
  double *starJ,
  double *starU
);
extern void setExteriorGlobalElementBoundaryVelocityValues(
  int updateFluxValues,
  int nExteriorElementBoundaries_global,
  int nQuadraturePoints_elementBoundary,
  int nSpace,
  int *exteriorElementBoundaries,
  int *elementBoundaryElements,
  int *elementBoundaryLocalElementBoundaries,
  double *n,
  double *vn_in,
  double *v_out
);
extern void calculateDimensionlessNumbersADR(
  int nElements_global,
  int nQuadraturePoints_element,
  int nSpace,
  int computeDiffusiveTimeStepLimit,
  double *elementDiameter,
  double *df,
  double *a,
  double *dphi,
  double *dr,
  double *dmt,
  double *pe,
  double *cfl
);
extern void calculateCFLADR(
  int nElements_global,
  int nQuadraturePoints_element,
  int nSpace,
  double *elementDiameter,
  double *dm,
  double *df,
  double *cfl
);
extern void updateInteriorElementBoundaryDiffusiveVelocity(
  int nInteriorElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nSpace,
  int *interiorElementBoundaries,
  int *elementBoundaryElements,
  int *elementBoundaryLocalElementBoundaries,
  double *a,
  double *grad_phi,
  double *velocity
);
extern void updateInteriorElementBoundaryDiffusiveVelocity_sd(int nInteriorElementBoundaries_global,
							      int nElementBoundaries_element,
							      int nQuadraturePoints_elementBoundary,
							      int nSpace,
							      int* rowptr,
							      int* colind,
							      int* interiorElementBoundaries,
							      int* elementBoundaryElements,
							      int* elementBoundaryLocalElementBoundaries,
							      double* a,
							      double* grad_phi,
							      double* velocity);
extern void updateExteriorElementBoundaryDiffusiveVelocity(
  int nExteriorElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nSpace,
  int *exteriorElementBoundaries,
  int *elementBoundaryElements,
  int *elementBoundaryLocalElementBoundaries,
  double *a,
  double *grad_phi,
  double *velocity
);
extern void updateExteriorElementBoundaryDiffusiveVelocity_sd(int nExteriorElementBoundaries_global,
							      int nElementBoundaries_element,
							      int nQuadraturePoints_elementBoundary,
							      int nSpace,
							      int* rowptr,
							      int* colind,
							      int* exteriorElementBoundaries,
							      int* elementBoundaryElements,
							      int* elementBoundaryLocalElementBoundaries,
							      double* a,
							      double* grad_phi,
							      double* velocity);
extern void updateInteriorElementBoundaryAdvectiveVelocity(
  int nInteriorElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nSpace,
  int *interiorElementBoundaries,
  int *elementBoundaryElements,
  int *elementBoundaryLocalElementBoundaries,
  double *f,
  double *velocity
);
extern void updateExteriorElementBoundaryAdvectiveVelocity(
  int nExteriorElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nSpace,
  int *exteriorElementBoundaries,
  int *elementBoundaryElements,
  int *elementBoundaryLocalElementBoundaries,
  double *f,
  double *velocity
);
extern void updateInteriorElementBoundaryShockCapturingVelocity(
  int nInteriorElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nQuadraturePoints_element,
  int nSpace,
  int *interiorElementBoundaries,
  int *elementBoundaryElements,
  int *elementBoundaryLocalElementBoundaries,
  double *numDiff,
  double *grad_u,
  double *velocity
);
extern void updateExteriorElementBoundaryShockCapturingVelocity(
  int nExteriorElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nQuadraturePoints_element,
  int nSpace,
  int *exteriorElementBoundaries,
  int *elementBoundaryElements,
  int *elementBoundaryLocalElementBoundaries,
  double *numDiff,
  double *grad_u,
  double *velocity
);
extern void calculateInteriorElementBoundaryAverageVelocity(
  int nInteriorElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nSpace,
  int *interiorElementBoundaries,
  int *elementBoundaryElements,
  int *elementBoundaryLocalElementBoundaries,
  double *v,
  double *vAverage
);
extern void calculateExteriorElementBoundaryAverageVelocity(
  int nExteriorElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nSpace,
  int *exteriorElementBoundaries,
  int *elementBoundaryElements,
  int *elementBoundaryLocalElementBoundaries,
  double *v,
  double *vAverage
);
extern void calculateConservationResidualDG(
  int nElements_global,
  int nDOF_test_element,
  double *elementResidual,
  double *conservationResidual
);
extern void calculateConservationResidual(
  int nElements_global,
  int nDOF_test_element,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nSpace,
  double *n,
  double *dS_u,
  double *elementResidual,
  double *velocity,
  double *conservationResidual
);
extern void calculateConservationResidualGlobalBoundaries(
  int nElements_global,
  int nInteriorElementBoundaries_global,
  int nExteriorElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nNodes_element,
  int nSpace,
  int *interiorElementBoundaries,
  int *exteriorElementBoundaries,
  int *elementBoundaryElements,
  int *elementBoundaryLocalElementBoundaries,
  double *dS,
  double *normal,
  double *elementResidual,
  double *velocity,
  double *conservationResidual
);
extern void copyGlobalElementBoundaryVelocityToElementBoundary(
  int nElements_global,
  int nInteriorElementBoundaries_global,
  int nExteriorElementBoundaries_global,
  int nElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nSpace,
  int *interiorElementBoundaries,
  int *exteriorElementBoundaries,
  int *elementBoundaryElementsArray,
  int *elementBoundaryLocalElementBoundariesArray,
  double *velocityBoundary_global,
  double *velocityBoundary_element
);
extern void loadBoundaryFluxIntoGlobalElementBoundaryVelocity(
  int nExteriorElementBoundaries_global,
  int nQuadraturePoints_elementBoundary,
  int nSpace,
  int *exteriorElementBoundaries,
  int *fluxElementBoundaries,
  double *normal,
  double *flux,
  double updateCoef,
  double *velocity
);
extern void calculateInteriorNumericalTrace_Potential(
  int nInteriorElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int *interiorElementBoundaries,
  int *elementBoundaryElements,
  int *elementBoundaryLocalElementBoundaries,
  double *phi,
  double *dphi,
  double *phi_trace,
  double *dphi_trace_left,
  double *dphi_trace_right
);
extern void calculateExteriorNumericalTrace_Potential(
  int *isDOFBoundary,
  int nExteriorElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int *exteriorElementBoundaries,
  int *elementBoundaryElements,
  int *elementBoundaryLocalElementBoundaries,
  double *phi_bc,
  double *phi,
  double *dphi,
  double *phi_trace,
  double *dphi_trace_left
);
extern void updateInteriorElementBoundary_MixedForm_weak(
  int nInteriorElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nDOF_test_element,
  int nSpace,
  int *interiorElementBoundaries,
  int *elementBoundaryElements,
  int *elementBoundaryLocalElementBoundaries,
  double *n,
  double *phi_trace,
  double *w_dS,
  double *b
);
extern void updateInteriorElementBoundary_MixedForm_weakJacobian(
  int nInteriorElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nDOF_test_element,
  int nSpace,
  int *interiorElementBoundaries,
  int *elementBoundaryElements,
  int *elementBoundaryLocalElementBoundaries,
  double *n,
  double *dphi_trace_left,
  double *dphi_trace_right,
  double *v,
  double *w_dS,
  double *db,
  double *db_eb
);
extern void updateExteriorElementBoundary_MixedForm_weak(
  int nExteriorElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nDOF_test_element,
  int nSpace,
  int *exteriorElementBoundaries,
  int *elementBoundaryElements,
  int *elementBoundaryLocalElementBoundaries,
  double *n,
  double *phi_trace,
  double *w_dS,
  double *b
);
extern void updateExteriorElementBoundary_MixedForm_weakJacobian(
  int nExteriorElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nDOF_test_element,
  int nSpace,
  int *exteriorElementBoundaries,
  int *elementBoundaryElements,
  int *elementBoundaryLocalElementBoundaries,
  double *n,
  double *dphi_trace_left,
  double *v,
  double *w_dS,
  double *db,
  double *db_eb
);
extern void updatePotential_MixedForm_weak(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_test_element,
  int nSpace,
  double *phi,
  double *grad_w_dV,
  double *b
);
extern void updatePotential_MixedForm_weakJacobian(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_test_element,
  int nSpace,
  double *dphi,
  double *v,
  double *grad_w_dV,
  double *db
);
extern void calculateVelocityQuadrature_MixedForm(
  int nElements_global,
  int nElementBoundaries_element,
  int nElementBoundaryQuadraturePoints_elementBoundary,
  int nDOF_element,
  int nSpace,
  int nQuadraturePoints_element,
  double *A_inv,
  double *b,
  double *v,
  double *V,
  double *qv,
  double *qV
);
extern void calculateVelocityQuadrature_MixedForm_Jacobian(
  int nElements_global,
  int nElementBoundaries_element,
  int nElementBoundaryQuadraturePoints_elementBoundary,
  int nDOF_element,
  int nSpace,
  int nQuadraturePoints_element,
  double *A_inv,
  double *db,
  double *db_eb,
  double *v,
  double *DV,
  double *DV_eb,
  double *qv,
  double *qDV,
  double *qDV_eb
);
extern void calculateVelocityProjectionMatrixLDG(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_element,
  double *vXw_dV,
  double *A_inv
);
extern void updateDiffusion_MixedForm_weak(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_test_element,
  int nSpace,
  double *a,
  double *qV,
  double *grad_w_dV,
  double *weak_residual
);
extern void updateDiffusionJacobian_MixedForm_weak(
  int nElements_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_element,
  int nDOF_trial_element,
  int nDOF_test_element,
  int nSpace,
  double *a,
  double *da,
  double *qV,
  double *qDV,
  double *qDV_eb,
  double *grad_w_dV,
  double *v,
  double *jacobian_weak_residual,
  double *jacobian_weak_residual_eb
);
extern void estimate_mt(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_element,
  double *v,
  double *vXw_dV,
  double *elementSpatialResidual,
  double *mt
  );
extern void estimate_mt_lowmem(int nElements_global,
                               int nQuadraturePoints_element,
                               int nDOF_element,
                               double* v,
                               double* w_dV,
                               double* elementSpatialResidual,
                               double* mt);
extern double scalarDomainIntegral(int nElements_global,
                                   int nQuadraturePoints_element,
                                   double* dV,
                                   double* nValueArray);
extern double scalarHeavisideDomainIntegral(int nElements_global,
                                            int nQuadraturePoints_element,
                                            double* dV,
                                            double* nValueArray);
extern double scalarSmoothedHeavisideDomainIntegral(int nElements_global,
						    int nQuadraturePoints_element,
						    double epsFact,
						    double* elementDiameter,
						    double* dV,
						    double* nValueArray);
extern double fluxDomainBoundaryIntegral(int nExteriorElementBoundaries,
					 int nElementBoundaries_owned,
                                         int nQuadraturePoints_elementBoundary,
                                         int* flag,
                                         int* exteriorElementBoundariesArray,
                                         double* dS,
                                         double* nValueArray);

extern double fluxDomainBoundaryIntegralFromVector(int nExteriorElementBoundaries,
						   int nElementBoundaries_owned,
						   int nQuadraturePoints_elementBoundary,
						   int nSpace,
						   int* flag,
						   int* exteriorElementBoundaries,
						   double* dS,
						   double* nValueArray,
						   double* normal);
extern
void copyExteriorElementBoundaryValuesFromElementBoundaryValues(int nExteriorElementBoundaries_global,
								int nElements_global,
								int nElementBoundaries_element,
								int nQuadraturePoints_elementBoundary,
								int nValuesPerQuadraturePoint,
								const int * exteriorElementBoundaries,
								const int* elementBoundaryElements,
								const int * elementBoundaryLocalElementBoundaries,
								const double * ebq_val,
								double * ebqe_val);

extern
void copyExteriorElementBoundaryValuesToElementBoundaryValues(int nExteriorElementBoundaries_global,
							      int nElements_global,
							      int nElementBoundaries_element,
							      int nQuadraturePoints_elementBoundary,
							      int nValuesPerQuadraturePoint,
							      const int * exteriorElementBoundaries,
							      const int* elementBoundaryElements,
							      const int * elementBoundaryLocalElementBoundaries,
							      const double * ebqe_val,
							      double * ebq_val);

extern
void copyExteriorElementBoundaryValuesToGlobalElementBoundaryValues(int nExteriorElementBoundaries_global,
								    int nQuadraturePoints_elementBoundary,
								    int nValuesPerQuadraturePoint,
								    const int * exteriorElementBoundaries,
								    const int* elementBoundaryElements,
								    const int * elementBoundaryLocalElementBoundaries,
								    const double * ebqe_val,
								    double * ebq_global_val);



extern
void copyExteriorElementBoundaryValuesFromGlobalElementBoundaryValues(int nExteriorElementBoundaries_global,
								      int nQuadraturePoints_elementBoundary,
								      int nValuesPerQuadraturePoint,
								      const int * exteriorElementBoundaries,
								      const int* elementBoundaryElements,
								      const int * elementBoundaryLocalElementBoundaries,
								      const double * ebq_global_val,
								      double * ebqe_val);




extern
void computeC0P1InterpolantDGP0(int nElements_global,
				int nNodes_global,
				int nNodes_element,
				int nDOF_element,
				int dim_dof,
				const int* elementNodesArray,
				const int* nodeElementOffsets,
				const int* nodeElementsArray,
				const int* l2g,
				const double * dof,
				double* nodalAverage);

extern
void computeC0P1InterpolantNCP1(int nElements_global,
				int nNodes_global,
				int nNodes_element,
				int nDOF_element,
				int dim_dof,
				const int* elementNodesArray,
				const int* nodeElementOffsets,
				const int* nodeElementsArray,
				const int* l2g,
				const double * dof,
				double* nodalAverage);

extern
void computeC0P1InterpolantDGP12(int nElements_global,
				 int nNodes_global,
				 int nNodes_element,
				 int nDOF_element,
				 int dim_dof,
				 const int* elementNodesArray,
				 const int* nodeElementOffsets,
				 const int* nodeElementsArray,
				 const int* l2g,
				 const double * dof,
				 double* nodalAverage);




extern
void parametricFiniteElementSpace_getValuesGlobalExteriorTrace(int nElementBoundaries_element,
							       int nElementBoundaryQuadraturePoints_elementBoundary,
							       int nDOF_element,
							       int nExteriorElementBoundaries_global,
							       const int* exteriorElementBoundariesArray,
							       const int* elementBoundaryElementsArray,
							       const int* elementBoundaryLocalElementBoundariesArray,
							       double* psi,
							       double* vArray);


extern
void parametricFiniteElementSpace_getGradientValuesGlobalExteriorTrace(int nElementBoundaries_element,
								       int nElementBoundaryQuadraturePoints_elementBoundary,
								       int nDOF_element,
								       int nSpace_global,
								       int nExteriorElementBoundaries_global,
								       const int *exteriorElementBoundariesArray,
								       const int *elementBoundaryElementsArray,
								       const int *elementBoundaryLocalElementBoundariesArray,
								       double* grad_psi,
								       double* inverseJacobianArray,
								       double* grad_vArray);


extern
void parametricMaps_getValuesGlobalExteriorTrace(int nQuadraturePoints_elementBoundary,
						 int nDOF_element,
						 int nSpace_global,
						 int nExteriorElementBoundaries_global,
						 const int* exteriorElementBoundariesArray,
						 const int* elementBoundaryElementsArray,
						 const int* elementBoundaryLocalElementBoundariesArray,
						 double* psi,
						 int* l2g,
						 double* nodeArray,
						 double* xArray);


extern
void parametricMaps_getInverseValuesGlobalExteriorTrace(int nElementBoundaryQuadraturePoints_elementBoundary,
							int nDOF_element,
							int nSpace_global,
							int nExteriorElementBoundaries_global,
							const int* exteriorElementBoundariesArray,
							const int* elementBoundaryElementsArray,
							const int* elementBoundaryLocalElementBoundariesArray,
							double* inverseJacobian,
							int* l2g,
							double* nodeArray,
							double* xArray,
							double* xiArray);



extern
void parametricMaps_getJacobianValuesGlobalExteriorTrace1D(int nQuadraturePoints_element,
							   int nDOF_element,
							   int nExteriorElementBoundaries_global,
							   const int * exteriorElementBoundariesArray,
							   const int * elementBoundaryElementsArray,
							   const int * elementBoundaryLocalElementBoundariesArray,
							   double* grad_psi,
							   double* boundaryNormals,
							   double* boundaryJacobians,
							   int* l2g,
							   double* nodeArray,
							   double* jacobianInverseArray,
							   double* metricTensorArray,
							   double* metricTensorDeterminantSqrtArray,
							   double* unitNormalArray);





extern
void parametricMaps_getJacobianValuesGlobalExteriorTrace2D(int nQuadraturePoints_element,
							   int nDOF_element,
							   int nExteriorElementBoundaries_global,
							   const int * exteriorElementBoundariesArray,
							   const int * elementBoundaryElementsArray,
							   const int * elementBoundaryLocalElementBoundariesArray,
							   double* grad_psi,
							   double* boundaryNormals,
							   double* boundaryJacobians,
							   int* l2g,
							   double* nodeArray,
							   double* jacobianInverseArray,
							   double* metricTensorArray,
							   double* metricTensorDeterminantSqrtArray,
							   double* unitNormalArray);
extern
void parametricMaps_getJacobianValuesGlobalExteriorTrace2D_movingDomain(int nQuadraturePoints_element,
									int nDOF_element,
									int nExteriorElementBoundaries_global,
									const int * exteriorElementBoundariesArray,
									const int * elementBoundaryElementsArray,
									const int * elementBoundaryLocalElementBoundariesArray,
									double* xtArray,
									double* grad_psi,
									double* boundaryNormals,
									double* boundaryJacobians,
									int* l2g,
									double* nodeArray,
									double* jacobianInverseArray,
									double* metricTensorArray,
									double* metricTensorDeterminantSqrtArray,
									double* unitNormalArray);
extern
void parametricMaps_getJacobianValuesGlobalExteriorTrace3D(int nQuadraturePoints_element,
							   int nDOF_element,
							   int nExteriorElementBoundaries_global,
							   const int * exteriorElementBoundariesArray,
							   const int * elementBoundaryElementsArray,
							   const int * elementBoundaryLocalElementBoundariesArray,
							   double* grad_psi,
							   double* boundaryNormals,
							   double* boundaryJacobians,
							   int* l2g,
							   double* nodeArray,
							   double* jacobianInverseArray,
							   double* metricTensorArray,
							   double* metricTensorDeterminantSqrtArray,
							   double* unitNormalArray);


extern
void updateGlobalExteriorElementBoundaryFlux(int nExteriorElementBoundaries_global,
					     int nQuadraturePoints_elementBoundary,
					     int nDOF_test_element,
					     int* exteriorElementBoundaries,
					     int* elementBoundaryElements,
					     int* elementBoundaryLocalElementBoundaries,
					     double* flux,
					     double* w_dS,
					     double* residual);



extern
void updateGlobalJacobianFromGlobalExteriorElementBoundaryFluxJacobian_eb_dense(int* elementNeighbors,
										int nElements_global,
										int nExteriorElementBoundaries_global,
										int nElementBoundaries_element,
										int nQuadraturePoints_elementBoundary,
										int nDOF_test_element,
										int nDOF_trial_element,
										int offset_r,
										int stride_r,
										int offset_u,
										int stride_u,
										int nFreeVDOF_global,
										int* exteriorElementBoundaries,
										int* elementBoundaryElements,
										int* elementBoundaryLocalElementBoundaries,
										int* nFreeDOF_element_r,
										int* freeLocal_r,
										int* freeGlobal_r,
										int* nFreeDOF_element_u,
										int* freeLocal_u,
										int* freeGlobal_u,
										double* elementBoundaryFluxJacobian_eb,
										double* w_dS,
										double* jac);

extern
void updateGlobalJacobianFromGlobalExteriorElementBoundaryFluxJacobian_dense(int nExteriorElementBoundaries_global,
									     int nQuadraturePoints_elementBoundary,
									     int nDOF_test_element,
									     int nDOF_trial_element,
									     int offset_r,
									     int stride_r,
									     int offset_u,
									     int stride_u,
									     int nFreeVDOF_global,
									     int* exteriorElementBoundaries,
									     int* elementBoundaryElements,
									     int* elementBoundaryLocalElementBoundaries,
									     int* nFreeDOF_element_r,
									     int* freeLocal_r,
									     int* freeGlobal_r,
									     int* nFreeDOF_element_u,
									     int* freeLocal_u,
									     int* freeGlobal_u,
									     double* elementBoundaryFluxJacobian,
									     double* w_dS,
									     double* jac);

extern
void updateGlobalJacobianFromGlobalExteriorElementBoundaryFluxJacobian_CSR(int nExteriorElementBoundaries_global,
									   int nQuadraturePoints_elementBoundary,
									   int nDOF_test_element,
									   int nDOF_trial_element,
									   int* exteriorElementBoundaries,
									   int* elementBoundaryElements,
									   int* elementBoundaryLocalElementBoundaries,
									   int* nFreeDOF_element_r,
									   int* freeLocal_r,
									   int* nFreeDOF_element_u,
									   int* freeLocal_u,
									   int* csrRowIndeces_ru,
									   int* csrColumnOffsets_eb_ru,
									   double* elementBoundaryFluxJacobian,
									   double* w_dS,
									   double* jac);



extern
void updateGlobalJacobianFromGlobalExteriorElementBoundaryFluxJacobian_eb_CSR(int* elementNeighbors,
									      int nExteriorElementBoundaries_global,
									      int nElementBoundaries_element,
									      int nQuadraturePoints_elementBoundary,
									      int nDOF_test_element,
									      int nDOF_trial_element,
									      int* exteriorElementBoundaries,
									      int* elementBoundaryElements,
									      int* elementBoundaryLocalElementBoundaries,
									      int* nFreeDOF_element_r,
									      int* freeLocal_r,
									      int* nFreeDOF_element_u,
									      int* freeLocal_u,
									      int* csrRowIndeces_ru,
									      int* csrColumnOffsets_eb_eNebN_ru,
									      double* elementBoundaryFluxJacobian_eb,
									      double* w_dS,
									      double* jac);

extern
void calculateWeightedShapeGlobalExteriorTrace(int nElementBoundaryQuadraturePoints_elementBoundary,
					       int nDOF_test_element,
					       int nExteriorElementBoundaries_global,
					       const int* exteriorElementBoundariesArray,
					       const int* elementBoundaryElementsArray,
					       const int* elementBoundaryLocalElementBoundariesArray,
					       double* dSR,
					       double* sqrt_det_g,
					       double* w,
					       double* w_dS);


extern
void calculateShape_X_weightedShapeGlobalExteriorTrace(int nElementBoundaryQuadraturePoints_elementBoundary,
						       int nDOF_trial_element,
						       int nDOF_test_element,
						       int nExteriorElementBoundaries_global,
						       const int* exteriorElementBoundariesArray,
						       const int* elementBoundaryElementsArray,
						       const int* elementBoundaryLocalElementBoundariesArray,
						       double* v,
						       double* w_dS,
						       double* v_X_w_dS);

extern
void calculateGradShape_X_weightedShapeGlobalExteriorTrace(int nElementBoundaryQuadraturePoints_elementBoundary,
							   int nDOF_trial_element,
							   int nDOF_test_element,
							   int nSpace,
							   int nExteriorElementBoundaries_global,
							   const int* exteriorElementBoundariesArray,
							   const int* elementBoundaryElementsArray,
							   const int* elementBoundaryLocalElementBoundariesArray,
							   double* grad_v,
							   double* w_dS,
							   double* grad_v_X_w_dS);


extern
void calculateGlobalExteriorElementBoundaryIntegrationWeights(int nQuadraturePoints_elementBoundary,
							      int nExteriorElementBoundaries_global,
							      double* sqrt_det_g,
							      double* referenceWeights,
							      double* weights);

extern
void calculateFiniteElementFunctionValuesGlobalExteriorTrace(int nQuadraturePoints_elementBoundary,
							     int nDOF_trial_element,
							     int nComponents,
							     int nExteriorElementBoundaries_global,
							     const int * exteriorElementBoundariesArray,
							     const int * elementBoundaryElementsArray,
							     const int * elementBoundaryLocalElementBoundariesArray,
							     int* l2g,
							     double* dof,
							     double* v,
							     double* u);

extern
void calculateFiniteElementFunctionGradientValuesGlobalExteriorTrace(int nQuadraturePoints_elementBoundary,
								     int nDOF_trial_element,
								     int nComponents,
								     int nSpace,
								     int nExteriorElementBoundaries_global,
								     const int * exteriorElementBoundariesArray,
								     const int * elementBoundaryElementsArray,
								     const int * elementBoundaryLocalElementBoundariesArray,
								     int* l2g,
								     double* dof,
								     double* grad_v,
								     double* grad_u);


extern
void updateGlobalExteriorElementBoundaryDiffusiveVelocity(int nExteriorElementBoundaries_global,
							  int nQuadraturePoints_elementBoundary,
							  int nSpace,
							  int* exteriorElementBoundaries,
							  int* elementBoundaryElements,
							  int* elementBoundaryLocalElementBoundaries,
							  double* a,
							  double* grad_phi,
							  double* velocity);


extern
void updateGlobalExteriorElementBoundaryAdvectiveVelocity(int nExteriorElementBoundaries_global,
							  int nQuadraturePoints_elementBoundary,
							  int nSpace,
							  int* exteriorElementBoundaries,
							  int* elementBoundaryElements,
							  int* elementBoundaryLocalElementBoundaries,
							  double* f,
							  double* velocity);


extern
void updateGlobalExteriorElementBoundaryShockCapturingVelocity(int nExteriorElementBoundaries_global,
							       int nQuadraturePoints_elementBoundary,
							       int nSpace,
							       int* exteriorElementBoundaries,
							       int* elementBoundaryElements,
							       int* elementBoundaryLocalElementBoundaries,
							       double* numDiff,
							       double* grad_u,
							       double* velocity);



extern
void copyFreeUnknownsToGlobalUnknowns(int nDOF2set,
				      int offset,
				      int stride,
				      const int* globalDOFids,
				      const int* freeDOFids,
				      const double * free_u,
				      double * u);
extern
void copyGlobalUnknownsToFreeUnknowns(int nDOF2set,
				      int offset,
				      int stride,
				      const int* globalDOFids,
				      const int* freeDOFids,
				      const double * u,
				      double * free_u);
extern
void updateInteriorElementBoundaryDiffusionAdjoint(int nInteriorElementBoundaries_global,
						   int nElementBoundaries_element,
						   int nQuadraturePoints_elementBoundary,
						   int nDOF_test_element,
						   int nSpace,
						   int* interiorElementBoundaries,
						   int* elementBoundaryElements,
						   int* elementBoundaryLocalElementBoundaries,
						   double sigma,
						   double* u,
						   double* n,
						   double* a,
						   double* grad_w,
						   double* dS,
						   double* residual);
extern
void updateExteriorElementBoundaryDiffusionAdjoint(int nExteriorElementBoundaries_global,
						   int nQuadraturePoints_elementBoundary,
						   int nDOF_test_element,
						   int nSpace,
						   int* isDOFBoundary,
						   int* exteriorElementBoundaries,
						   int* elementBoundaryElements,
						   int* elementBoundaryLocalElementBoundaries,
						   double sigma,
						   double* u,
						   double* ub,
						   double* n,
						   double* a,
						   double* grad_w,
						   double* dS,
						   double* residual);
extern
void updateGlobalJacobianFromInteriorElementBoundaryDiffusionAdjoint_dense(int nInteriorElementBoundaries_global,
									   int nElementBoundaries_element,
									   int nQuadraturePoints_elementBoundary,
									   int nDOF_test_element,
									   int nDOF_trial_element,
									   int nSpace,
									   int offset_r,
									   int stride_r,
									   int offset_u,
									   int stride_u,
									   int nFreeVDOF_global,
									   int* interiorElementBoundaries,
									   int* elementBoundaryElements,
									   int* elementBoundaryLocalElementBoundaries,
									   int* nFreeDOF_element_r,
									   int* freeLocal_r,
									   int* freeGlobal_r,
									   int* nFreeDOF_element_u,
									   int* freeLocal_u,
									   int* freeGlobal_u,
									   double sigma,
									   double* v,
									   double* n,
									   double* a,
									   double* grad_w,
									   double* dS,
									   double* jac);
extern
void updateGlobalJacobianFromExteriorElementBoundaryDiffusionAdjoint_dense(int nExteriorElementBoundaries_global,
									   int nQuadraturePoints_elementBoundary,
									   int nDOF_test_element,
									   int nDOF_trial_element,
									   int nSpace,
									   int offset_r,
									   int stride_r,
									   int offset_u,
									   int stride_u,
									   int nFreeVDOF_global,
									   int* exteriorElementBoundaries,
									   int* elementBoundaryElements,
									   int* elementBoundaryLocalElementBoundaries,
									   int* nFreeDOF_element_r,
									   int* freeLocal_r,
									   int* freeGlobal_r,
									   int* nFreeDOF_element_u,
									   int* freeLocal_u,
									   int* freeGlobal_u,
									   int* isDOFBoundary,
									   double sigma,
									   double* v,
									   double* n,
									   double* a,
									   double* grad_w,
									   double* dS,
									   double* jac);
extern
void updateInteriorElementBoundaryDiffusionAdjoint_sd(int nInteriorElementBoundaries_global,
						      int nElementBoundaries_element,
						      int nQuadraturePoints_elementBoundary,
						      int nDOF_test_element,
						      int nSpace,
						      int* rowptr,
						      int* colind,
						      int* interiorElementBoundaries,
						      int* elementBoundaryElements,
						      int* elementBoundaryLocalElementBoundaries,
						      double sigma,
						      double* u,
						      double* n,
						      double* a,
						      double* grad_w,
						      double* dS,
						      double* residual);
extern
void updateExteriorElementBoundaryDiffusionAdjoint_sd(int nExteriorElementBoundaries_global,
						      int nQuadraturePoints_elementBoundary,
						      int nDOF_test_element,
						      int nSpace,
						      int* rowptr,
						      int* colind,
						      int* isDOFBoundary,
						      int* exteriorElementBoundaries,
						      int* elementBoundaryElements,
						      int* elementBoundaryLocalElementBoundaries,
						      double sigma,
						      double* u,
						      double* ub,
						      double* n,
						      double* a,
						      double* grad_w,
						      double* dS,
						      double* residual);
extern
void updateGlobalJacobianFromInteriorElementBoundaryDiffusionAdjoint_dense_sd(int nInteriorElementBoundaries_global,
									      int nElementBoundaries_element,
									      int nQuadraturePoints_elementBoundary,
									      int nDOF_test_element,
									      int nDOF_trial_element,
									      int nSpace,
									      int* rowptr,
									      int* colind,
									      int offset_r,
									      int stride_r,
									      int offset_u,
									      int stride_u,
									      int nFreeVDOF_global,
									      int* interiorElementBoundaries,
									      int* elementBoundaryElements,
									      int* elementBoundaryLocalElementBoundaries,
									      int* nFreeDOF_element_r,
									      int* freeLocal_r,
									      int* freeGlobal_r,
									      int* nFreeDOF_element_u,
									      int* freeLocal_u,
									      int* freeGlobal_u,
									      double sigma,
									      double* v,
									      double* n,
									      double* a,
									      double* grad_w,
									      double* dS,
									      double* jac);
extern
void updateGlobalJacobianFromExteriorElementBoundaryDiffusionAdjoint_dense_sd(int nExteriorElementBoundaries_global,
									      int nQuadraturePoints_elementBoundary,
									      int nDOF_test_element,
									      int nDOF_trial_element,
									      int nSpace,
									      int* rowptr,
									      int* colind,
									      int offset_r,
									      int stride_r,
									      int offset_u,
									      int stride_u,
									      int nFreeVDOF_global,
									      int* exteriorElementBoundaries,
									      int* elementBoundaryElements,
									      int* elementBoundaryLocalElementBoundaries,
									      int* nFreeDOF_element_r,
									      int* freeLocal_r,
									      int* freeGlobal_r,
									      int* nFreeDOF_element_u,
									      int* freeLocal_u,
									      int* freeGlobal_u,
									      int* isDOFBoundary,
									      double sigma,
									      double* v,
									      double* n,
									      double* a,
									      double* grad_w,
									      double* dS,
									      double* jac);
extern
void updateGlobalJacobianFromInteriorElementBoundaryDiffusionAdjoint_CSR_sd(int nInteriorElementBoundaries_global,
									    int nElementBoundaries_element,
									    int nQuadraturePoints_elementBoundary,
									    int nDOF_test_element,
									    int nDOF_trial_element,
									    int nSpace,
									    int* rowptr,
									    int* colind,
									    int offset_r,
									    int stride_r,
									    int offset_u,
									    int stride_u,
									    int nFreeVDOF_global,
									    int* interiorElementBoundaries,
									    int* elementBoundaryElements,
									    int* elementBoundaryLocalElementBoundaries,
									    int* nFreeDOF_element_r,
									    int* freeLocal_r,
									    int* freeGlobal_r,
									    int* nFreeDOF_element_u,
									    int* freeLocal_u,
									    int* freeGlobal_u,
									    int* csrRowIndeces_ru,
									    int* csrColumnOffsets_eb_ru,
									    double sigma,
									    double* v,
									    double* n,
									    double* a,
									    double* grad_w,
									    double* dS,
									    double* jac);
extern
void updateGlobalJacobianFromExteriorElementBoundaryDiffusionAdjoint_CSR_sd(int nExteriorElementBoundaries_global,
									    int nQuadraturePoints_elementBoundary,
									    int nDOF_test_element,
									    int nDOF_trial_element,
									    int nSpace,
									    int* rowptr,
									    int* colind,
									    int offset_r,
									    int stride_r,
									    int offset_u,
									    int stride_u,
									    int nFreeVDOF_global,
									    int* exteriorElementBoundaries,
									    int* elementBoundaryElements,
									    int* elementBoundaryLocalElementBoundaries,
									    int* nFreeDOF_element_r,
									    int* freeLocal_r,
									    int* freeGlobal_r,
									    int* nFreeDOF_element_u,
									    int* freeLocal_u,
									    int* freeGlobal_u,
									    int* csrRowIndeces_ru,
									    int* csrColumnOffsets_eb_ru,
									    int* isDOFBoundary,
									    double sigma,
									    double* v,
									    double* n,
									    double* a,
									    double* grad_w,
									    double* dS,
									    double* jac);
extern
void update_f_movingDomain_q(int nElements_global,
			     int nQuadraturePoints_element,
			     int nSpace,
			     double* xt,
			     double* m,
			     double* f);
extern
void update_f_movingDomain_constantMass_q(int nElements_global,
					  int nQuadraturePoints_element,
					  int nSpace,
					  double* xt,
					  double* f);
extern
void update_f_movingDomain_ebq(int nElements_global,
			       int nElementBoundaries_element,
			       int nQuadraturePoints_elementBoundary,
			       int nSpace,
			       double* xt,
			       double* m,
			       double* f);
extern
void update_f_movingDomain_constantMass_ebq(int nElements_global,
					    int nElementBoundaries_element,
					    int nQuadraturePoints_elementBoundary,
					    int nSpace,
					    double* xt,
					    double* f);
extern
void updateStress_weak(int nElements_global,
		       int nQuadraturePoints_element,
		       int nDOF_test_element,
		       int nSpace,
		       double* sigma,
		       double* grad_w_dV,
		       double* weak_residual_x,
		       double* weak_residual_y,
		       double* weak_residual_z);
extern
void updateStressJacobian_weak(int nElements_global,
			       int nQuadraturePoints_element,
			       int nDOF_trial_element,
			       int nDOF_test_element,
			       int nSpace,
			       double* dsigma_xx,
			       double* dsigma_xy,
			       double* dsigma_xz,
			       double* dsigma_yx,
			       double* dsigma_yy,
			       double* dsigma_yz,
			       double* dsigma_zx,
			       double* dsigma_zy,
			       double* dsigma_zz,
			       double* grad_v,
			       double* grad_w_dV,
			       double* jacobian_weak_residual_xx,
			       double* jacobian_weak_residual_xy,
			       double* jacobian_weak_residual_xz,
			       double* jacobian_weak_residual_yx,
			       double* jacobian_weak_residual_yy,
			       double* jacobian_weak_residual_yz,
			       double* jacobian_weak_residual_zx,
			       double* jacobian_weak_residual_zy,
			       double* jacobian_weak_residual_zz);
extern
void projectFromNodalInterpolationConditions(int nElements_global,
					     int nDOF_element,
					     int dim_dof,
					     const int * l2g,
					     const int * functional_map_element,
					     const double * interpolationValues,
					     double * dofs);
/** @} */
#endif

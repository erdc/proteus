#ifndef CLSVOF_H
#define CLSVOF_H
#include <cmath>
#include <iostream>
#include "CompKernel.h"
#include "ModelFactory.h"

#define betaNormGrad(gradu2,beta2) std::sqrt(gradu2+beta2)

namespace proteus
{
  class CLSVOF_base
  {
    //The base class defining the interface
  public:
    virtual ~CLSVOF_base(){}
    virtual void calculateResidual_MCorr_with_VOF(//element
						  double dt,
						  double* mesh_trial_ref,
						  double* mesh_grad_trial_ref,
						  double* mesh_dof,
						  double* mesh_velocity_dof,
						  double MOVING_DOMAIN,
						  int* mesh_l2g,
						  double* dV_ref,
						  double* u_trial_ref,
						  double* u_grad_trial_ref,
						  double* u_test_ref,
						  double* u_grad_test_ref,
						  //element boundary
						  double* mesh_trial_trace_ref,
						  double* mesh_grad_trial_trace_ref,
						  double* dS_ref,
						  double* u_trial_trace_ref,
						  double* u_grad_trial_trace_ref,
						  double* u_test_trace_ref,
						  double* u_grad_test_trace_ref,
						  double* normal_ref,
						  double* boundaryJac_ref,
						  //physics
						  int nElements_global,
						  double useMetrics, 
						  double alphaBDF,
						  int lag_shockCapturing,
						  double shockCapturingDiffusion,
						  double sc_uref, 
						  double sc_alpha,
						  //VRANS
						  const double* q_porosity,
						  const double* porosity_dof,
						  //
						  int* u_l2g, 
						  double* elementDiameter,
						  int degree_polynomial,
						  double* u_dof,
						  double* u_dof_old,
						  double* uStar_dof,
						  double* velocity,
						  double* q_m,
						  double* q_u,
						  double* q_m_betaBDF,
						  double* q_dV,
						  double* q_dV_last,
						  double* cfl,
						  double* edge_based_cfl,
						  double* q_numDiff_u, 
						  double* q_numDiff_u_last, 
						  int offset_u, int stride_u, 
						  double* globalResidual,
						  int nExteriorElementBoundaries_global,
						  int* exteriorElementBoundariesArray,
						  int* elementBoundaryElementsArray,
						  int* elementBoundaryLocalElementBoundariesArray,
						  double* ebqe_velocity_ext,
						  //VRANS
						  const double* ebqe_porosity_ext,
						  //
						  int* isDOFBoundary_u,
						  double* ebqe_bc_u_ext,
						  int* isFluxBoundary_u,
						  double* ebqe_bc_flux_u_ext,
						  double* ebqe_phi,double epsFact,
						  double* ebqe_u,
						  double* ebqe_flux,
						  // PARAMETERS FOR EDGE BASED STABILIZATION
						  double cE,
						  double cK,
						  // PARAMETERS FOR LOG BASED ENTROPY FUNCTION 
						  double uL, 
						  double uR, 
						  // PARAMETERS FOR EDGE VISCOSITY 
						  int numDOFs,
						  int NNZ,
						  int* csrRowIndeces_DofLoops,
						  int* csrColumnOffsets_DofLoops,
						  int* csrRowIndeces_CellLoops,
						  int* csrColumnOffsets_CellLoops,
						  int* csrColumnOffsets_eb_CellLoops,
						  // C matrices
						  double* Cx, 
						  double* Cy,
						  double* Cz,
						  double* CTx,
						  double* CTy,
						  double* CTz,
						  double* ML,
						  // PARAMETERS FOR 1st or 2nd ORDER MPP METHOD
						  int LUMPED_MASS_MATRIX, 
						  int STABILIZATION_TYPE,
						  int ENTROPY_TYPE,
						  // FOR FCT
						  double* low_order_solution,
						  double* dt_times_dH_minus_dL,
						  double* min_u_bc,
						  double* max_u_bc,
						  // FOR NONLINEAR VOF; i.e., MCorr with VOF
						  int useFullNewton,
						  double epsFactHeaviside,
						  double epsFactDirac,
						  double epsFactDiffusion,
						  double* phin_dof,
						  double* phiHat_dof,
						  double* lumped_wx,
						  double* lumped_wy,
						  // AUX QUANTITIES OF INTEREST
						  double* quantDOFs)=0;
    virtual void calculateResidual_MCorr_with_VOF2(//element
						  double dt,
						  double* mesh_trial_ref,
						  double* mesh_grad_trial_ref,
						  double* mesh_dof,
						  double* mesh_velocity_dof,
						  double MOVING_DOMAIN,
						  int* mesh_l2g,
						  double* dV_ref,
						  double* u_trial_ref,
						  double* u_grad_trial_ref,
						  double* u_test_ref,
						  double* u_grad_test_ref,
						  //element boundary
						  double* mesh_trial_trace_ref,
						  double* mesh_grad_trial_trace_ref,
						  double* dS_ref,
						  double* u_trial_trace_ref,
						  double* u_grad_trial_trace_ref,
						  double* u_test_trace_ref,
						  double* u_grad_test_trace_ref,
						  double* normal_ref,
						  double* boundaryJac_ref,
						  //physics
						  int nElements_global,
						  double useMetrics, 
						  double alphaBDF,
						  int lag_shockCapturing,
						  double shockCapturingDiffusion,
						  double sc_uref, 
						  double sc_alpha,
						  //VRANS
						  const double* q_porosity,
						  const double* porosity_dof,
						  //
						  int* u_l2g, 
						  double* elementDiameter,
						  int degree_polynomial,
						  double* u_dof,
						  double* u_dof_old,
						  double* uStar_dof,
						  double* velocity,
						  double* q_m,
						  double* q_u,
						  double* q_m_betaBDF,
						  double* q_dV,
						  double* q_dV_last,
						  double* cfl,
						  double* edge_based_cfl,
						  double* q_numDiff_u, 
						  double* q_numDiff_u_last, 
						  int offset_u, int stride_u, 
						  double* globalResidual,
						  int nExteriorElementBoundaries_global,
						  int* exteriorElementBoundariesArray,
						  int* elementBoundaryElementsArray,
						  int* elementBoundaryLocalElementBoundariesArray,
						  double* ebqe_velocity_ext,
						  //VRANS
						  const double* ebqe_porosity_ext,
						  //
						  int* isDOFBoundary_u,
						  double* ebqe_bc_u_ext,
						  int* isFluxBoundary_u,
						  double* ebqe_bc_flux_u_ext,
						  double* ebqe_phi,double epsFact,
						  double* ebqe_u,
						  double* ebqe_flux,
						  // PARAMETERS FOR EDGE BASED STABILIZATION
						  double cE,
						  double cK,
						  // PARAMETERS FOR LOG BASED ENTROPY FUNCTION 
						  double uL, 
						  double uR, 
						  // PARAMETERS FOR EDGE VISCOSITY 
						  int numDOFs,
						  int NNZ,
						  int* csrRowIndeces_DofLoops,
						  int* csrColumnOffsets_DofLoops,
						  int* csrRowIndeces_CellLoops,
						  int* csrColumnOffsets_CellLoops,
						  int* csrColumnOffsets_eb_CellLoops,
						  // C matrices
						  double* Cx, 
						  double* Cy,
						  double* Cz,
						  double* CTx,
						  double* CTy,
						  double* CTz,
						  double* ML,
						  // PARAMETERS FOR 1st or 2nd ORDER MPP METHOD
						  int LUMPED_MASS_MATRIX, 
						  int STABILIZATION_TYPE,
						  int ENTROPY_TYPE,
						  // FOR FCT
						  double* low_order_solution,
						  double* dt_times_dH_minus_dL,
						  double* min_u_bc,
						  double* max_u_bc,
						  // FOR NONLINEAR VOF; i.e., MCorr with VOF
						  int useFullNewton,
						  double epsFactHeaviside,
						  double epsFactDirac,
						  double epsFactDiffusion,
						  double* phin_dof,
						  double* phiHat_dof,
						  double* lumped_wx,
						  double* lumped_wy,
						  // AUX QUANTITIES OF INTEREST
						  double* quantDOFs)=0;
    virtual void calculateResidual_MCorr_with_VOF3(//element ... working 1st version of consLS with projected reg/pen
						  double dt,
						  double* mesh_trial_ref,
						  double* mesh_grad_trial_ref,
						  double* mesh_dof,
						  double* mesh_velocity_dof,
						  double MOVING_DOMAIN,
						  int* mesh_l2g,
						  double* dV_ref,
						  double* u_trial_ref,
						  double* u_grad_trial_ref,
						  double* u_test_ref,
						  double* u_grad_test_ref,
						  //element boundary
						  double* mesh_trial_trace_ref,
						  double* mesh_grad_trial_trace_ref,
						  double* dS_ref,
						  double* u_trial_trace_ref,
						  double* u_grad_trial_trace_ref,
						  double* u_test_trace_ref,
						  double* u_grad_test_trace_ref,
						  double* normal_ref,
						  double* boundaryJac_ref,
						  //physics
						  int nElements_global,
						  double useMetrics, 
						  double alphaBDF,
						  int lag_shockCapturing,
						  double shockCapturingDiffusion,
						  double sc_uref, 
						  double sc_alpha,
						  //VRANS
						  const double* q_porosity,
						  const double* porosity_dof,
						  //
						  int* u_l2g, 
						  double* elementDiameter,
						  int degree_polynomial,
						  double* u_dof,
						  double* u_dof_old,
						  double* uStar_dof,
						  double* velocity,
						  double* q_m,
						  double* q_u,
						  double* q_m_betaBDF,
						  double* q_dV,
						  double* q_dV_last,
						  double* cfl,
						  double* edge_based_cfl,
						  double* q_numDiff_u, 
						  double* q_numDiff_u_last, 
						  int offset_u, int stride_u, 
						  double* globalResidual,
						  int nExteriorElementBoundaries_global,
						  int* exteriorElementBoundariesArray,
						  int* elementBoundaryElementsArray,
						  int* elementBoundaryLocalElementBoundariesArray,
						  double* ebqe_velocity_ext,
						  //VRANS
						  const double* ebqe_porosity_ext,
						  //
						  int* isDOFBoundary_u,
						  double* ebqe_bc_u_ext,
						  int* isFluxBoundary_u,
						  double* ebqe_bc_flux_u_ext,
						  double* ebqe_phi,double epsFact,
						  double* ebqe_u,
						  double* ebqe_flux,
						  // PARAMETERS FOR EDGE BASED STABILIZATION
						  double cE,
						  double cK,
						  // PARAMETERS FOR LOG BASED ENTROPY FUNCTION 
						  double uL, 
						  double uR, 
						  // PARAMETERS FOR EDGE VISCOSITY 
						  int numDOFs,
						  int NNZ,
						  int* csrRowIndeces_DofLoops,
						  int* csrColumnOffsets_DofLoops,
						  int* csrRowIndeces_CellLoops,
						  int* csrColumnOffsets_CellLoops,
						  int* csrColumnOffsets_eb_CellLoops,
						  // C matrices
						  double* Cx, 
						  double* Cy,
						  double* Cz,
						  double* CTx,
						  double* CTy,
						  double* CTz,
						  double* ML,
						  // PARAMETERS FOR 1st or 2nd ORDER MPP METHOD
						  int LUMPED_MASS_MATRIX, 
						  int STABILIZATION_TYPE,
						  int ENTROPY_TYPE,
						  // FOR FCT
						  double* low_order_solution,
						  double* dt_times_dH_minus_dL,
						  double* min_u_bc,
						  double* max_u_bc,
						  // FOR NONLINEAR VOF; i.e., MCorr with VOF
						  int useFullNewton,
						  double epsFactHeaviside,
						  double epsFactDirac,
						  double epsFactDiffusion,
						  double* phin_dof,
						  double* phiHat_dof,
						  double* lumped_wx,
						  double* lumped_wy,
						  // AUX QUANTITIES OF INTEREST
						  double* quantDOFs)=0;
    virtual void calculateResidual_MCorr_with_VOF4(//element... checking anisotropic regularization/penalization
						  double dt,
						  double* mesh_trial_ref,
						  double* mesh_grad_trial_ref,
						  double* mesh_dof,
						  double* mesh_velocity_dof,
						  double MOVING_DOMAIN,
						  int* mesh_l2g,
						  double* dV_ref,
						  double* u_trial_ref,
						  double* u_grad_trial_ref,
						  double* u_test_ref,
						  double* u_grad_test_ref,
						  //element boundary
						  double* mesh_trial_trace_ref,
						  double* mesh_grad_trial_trace_ref,
						  double* dS_ref,
						  double* u_trial_trace_ref,
						  double* u_grad_trial_trace_ref,
						  double* u_test_trace_ref,
						  double* u_grad_test_trace_ref,
						  double* normal_ref,
						  double* boundaryJac_ref,
						  //physics
						  int nElements_global,
						  double useMetrics, 
						  double alphaBDF,
						  int lag_shockCapturing,
						  double shockCapturingDiffusion,
						  double sc_uref, 
						  double sc_alpha,
						  //VRANS
						  const double* q_porosity,
						  const double* porosity_dof,
						  //
						  int* u_l2g, 
						  double* elementDiameter,
						  int degree_polynomial,
						  double* u_dof,
						  double* u_dof_old,
						  double* uStar_dof,
						  double* velocity,
						  double* q_m,
						  double* q_u,
						  double* q_m_betaBDF,
						  double* q_dV,
						  double* q_dV_last,
						  double* cfl,
						  double* edge_based_cfl,
						  double* q_numDiff_u, 
						  double* q_numDiff_u_last, 
						  int offset_u, int stride_u, 
						  double* globalResidual,
						  int nExteriorElementBoundaries_global,
						  int* exteriorElementBoundariesArray,
						  int* elementBoundaryElementsArray,
						  int* elementBoundaryLocalElementBoundariesArray,
						  double* ebqe_velocity_ext,
						  //VRANS
						  const double* ebqe_porosity_ext,
						  //
						  int* isDOFBoundary_u,
						  double* ebqe_bc_u_ext,
						  int* isFluxBoundary_u,
						  double* ebqe_bc_flux_u_ext,
						  double* ebqe_phi,double epsFact,
						  double* ebqe_u,
						  double* ebqe_flux,
						  // PARAMETERS FOR EDGE BASED STABILIZATION
						  double cE,
						  double cK,
						  // PARAMETERS FOR LOG BASED ENTROPY FUNCTION 
						  double uL, 
						  double uR, 
						  // PARAMETERS FOR EDGE VISCOSITY 
						  int numDOFs,
						  int NNZ,
						  int* csrRowIndeces_DofLoops,
						  int* csrColumnOffsets_DofLoops,
						  int* csrRowIndeces_CellLoops,
						  int* csrColumnOffsets_CellLoops,
						  int* csrColumnOffsets_eb_CellLoops,
						  // C matrices
						  double* Cx, 
						  double* Cy,
						  double* Cz,
						  double* CTx,
						  double* CTy,
						  double* CTz,
						  double* ML,
						  // PARAMETERS FOR 1st or 2nd ORDER MPP METHOD
						  int LUMPED_MASS_MATRIX, 
						  int STABILIZATION_TYPE,
						  int ENTROPY_TYPE,
						  // FOR FCT
						  double* low_order_solution,
						  double* dt_times_dH_minus_dL,
						  double* min_u_bc,
						  double* max_u_bc,
						  // FOR NONLINEAR VOF; i.e., MCorr with VOF
						  int useFullNewton,
						  double epsFactHeaviside,
						  double epsFactDirac,
						  double epsFactDiffusion,
						  double* phin_dof,
						  double* phiHat_dof,
						  double* lumped_wx,
						  double* lumped_wy,
						  // AUX QUANTITIES OF INTEREST
						  double* quantDOFs)=0;    
    virtual void calculateJacobian_MCorr_with_VOF(//element
				     double dt,
				     double* mesh_trial_ref,
				     double* mesh_grad_trial_ref,
				     double* mesh_dof,
				     double* mesh_velocity_dof,
				     double MOVING_DOMAIN,
				     int* mesh_l2g,
				     double* dV_ref,
				     double* u_trial_ref,
				     double* u_grad_trial_ref,
				     double* u_test_ref,
				     double* u_grad_test_ref,
				     //element boundary
				     double* mesh_trial_trace_ref,
				     double* mesh_grad_trial_trace_ref,
				     double* dS_ref,
				     double* u_trial_trace_ref,
				     double* u_grad_trial_trace_ref,
				     double* u_test_trace_ref,
				     double* u_grad_test_trace_ref,
				     double* normal_ref,
				     double* boundaryJac_ref,
				     //physics
				     int nElements_global,
				     double useMetrics, 
				     double alphaBDF,
				     int lag_shockCapturing,/*mwf not used yet*/
				     double shockCapturingDiffusion,
				     //VRANS
				     const double* q_porosity,
				     //
				     int* u_l2g,
				     double* elementDiameter,
				     int degree_polynomial,
				     double* u_dof, 
				     double* velocity,
				     double* q_m_betaBDF, 
				     double* cfl,
				     double* q_numDiff_u_last, 
				     int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
				     double* globalJacobian,
				     int nExteriorElementBoundaries_global,
				     int* exteriorElementBoundariesArray,
				     int* elementBoundaryElementsArray,
				     int* elementBoundaryLocalElementBoundariesArray,
				     double* ebqe_velocity_ext,
				     //VRANS
				     const double* ebqe_porosity_ext,
				     //
				     int* isDOFBoundary_u,
				     double* ebqe_bc_u_ext,
				     int* isFluxBoundary_u,
				     double* ebqe_bc_flux_u_ext,
				     int* csrColumnOffsets_eb_u_u,
				     int LUMPED_MASS_MATRIX,
				     // FOR NONLINEAR VOF; i.e., MCorr with VOF
				     int useFullNewton,
				     double epsFactHeaviside,
				     double epsFactDirac,
				     double epsFactDiffusion,
				     double cK,
				     double uL,
				     double uR,
				     double* phin_dof,
				     double* phiHat_dof,
				     double* lumped_wx,
				     double* lumped_wy
						  )=0;
    virtual void calculateJacobian_MCorr_with_VOF2(//element
				     double dt,
				     double* mesh_trial_ref,
				     double* mesh_grad_trial_ref,
				     double* mesh_dof,
				     double* mesh_velocity_dof,
				     double MOVING_DOMAIN,
				     int* mesh_l2g,
				     double* dV_ref,
				     double* u_trial_ref,
				     double* u_grad_trial_ref,
				     double* u_test_ref,
				     double* u_grad_test_ref,
				     //element boundary
				     double* mesh_trial_trace_ref,
				     double* mesh_grad_trial_trace_ref,
				     double* dS_ref,
				     double* u_trial_trace_ref,
				     double* u_grad_trial_trace_ref,
				     double* u_test_trace_ref,
				     double* u_grad_test_trace_ref,
				     double* normal_ref,
				     double* boundaryJac_ref,
				     //physics
				     int nElements_global,
				     double useMetrics, 
				     double alphaBDF,
				     int lag_shockCapturing,/*mwf not used yet*/
				     double shockCapturingDiffusion,
				     //VRANS
				     const double* q_porosity,
				     //
				     int* u_l2g,
				     double* elementDiameter,
				     int degree_polynomial,
				     double* u_dof, 
				     double* velocity,
				     double* q_m_betaBDF, 
				     double* cfl,
				     double* q_numDiff_u_last, 
				     int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
				     double* globalJacobian,
				     int nExteriorElementBoundaries_global,
				     int* exteriorElementBoundariesArray,
				     int* elementBoundaryElementsArray,
				     int* elementBoundaryLocalElementBoundariesArray,
				     double* ebqe_velocity_ext,
				     //VRANS
				     const double* ebqe_porosity_ext,
				     //
				     int* isDOFBoundary_u,
				     double* ebqe_bc_u_ext,
				     int* isFluxBoundary_u,
				     double* ebqe_bc_flux_u_ext,
				     int* csrColumnOffsets_eb_u_u,
				     int LUMPED_MASS_MATRIX,
				     // FOR NONLINEAR VOF; i.e., MCorr with VOF
				     int useFullNewton,
				     double epsFactHeaviside,
				     double epsFactDirac,
				     double epsFactDiffusion,
				     double cK,
				     double uL,
				     double uR,
				     double* phin_dof,
				     double* phiHat_dof,
				     double* lumped_wx,
				     double* lumped_wy)=0;
    virtual void calculateJacobian_MCorr_with_VOF4(//element
				     double dt,
				     double* mesh_trial_ref,
				     double* mesh_grad_trial_ref,
				     double* mesh_dof,
				     double* mesh_velocity_dof,
				     double MOVING_DOMAIN,
				     int* mesh_l2g,
				     double* dV_ref,
				     double* u_trial_ref,
				     double* u_grad_trial_ref,
				     double* u_test_ref,
				     double* u_grad_test_ref,
				     //element boundary
				     double* mesh_trial_trace_ref,
				     double* mesh_grad_trial_trace_ref,
				     double* dS_ref,
				     double* u_trial_trace_ref,
				     double* u_grad_trial_trace_ref,
				     double* u_test_trace_ref,
				     double* u_grad_test_trace_ref,
				     double* normal_ref,
				     double* boundaryJac_ref,
				     //physics
				     int nElements_global,
				     double useMetrics, 
				     double alphaBDF,
				     int lag_shockCapturing,/*mwf not used yet*/
				     double shockCapturingDiffusion,
				     //VRANS
				     const double* q_porosity,
				     //
				     int* u_l2g,
				     double* elementDiameter,
				     int degree_polynomial,
				     double* u_dof, 
				     double* velocity,
				     double* q_m_betaBDF, 
				     double* cfl,
				     double* q_numDiff_u_last, 
				     int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
				     double* globalJacobian,
				     int nExteriorElementBoundaries_global,
				     int* exteriorElementBoundariesArray,
				     int* elementBoundaryElementsArray,
				     int* elementBoundaryLocalElementBoundariesArray,
				     double* ebqe_velocity_ext,
				     //VRANS
				     const double* ebqe_porosity_ext,
				     //
				     int* isDOFBoundary_u,
				     double* ebqe_bc_u_ext,
				     int* isFluxBoundary_u,
				     double* ebqe_bc_flux_u_ext,
				     int* csrColumnOffsets_eb_u_u,
				     int LUMPED_MASS_MATRIX,
				     // FOR NONLINEAR VOF; i.e., MCorr with VOF
				     int useFullNewton,
				     double epsFactHeaviside,
				     double epsFactDirac,
				     double epsFactDiffusion,
				     double cK,
				     double uL,
				     double uR,
				     double* phin_dof,
				     double* phiHat_dof,
				     double* lumped_wx,
				     double* lumped_wy)=0;
    virtual void calculateRhsL2p(
				 double* mesh_trial_ref,
				 double* mesh_grad_trial_ref,
				 double* mesh_dof,
				 int* mesh_l2g,
				 double* dV_ref,
				 double* u_trial_ref,
				 double* u_grad_trial_ref,
				 double* u_test_ref,
				 //physics
				 int nElements_global,
				 int* u_l2g, 
				 double* elementDiameter,
				 //double* nodeDiametersArray,
				 double* u_dof,
				 double* phiHat_dof,
				 double* phiExact_dof,
				 int offset_u, int stride_u, 
				 double* globalResidual,
				 double* global_mass_error,
				 double* global_L2_interface,
				 double* global_H1_interface,
				 double* global_L2_Hinterface,
				 double* global_H1_Hinterface,
				 double* global_L2_u,
				 double* global_H1_u)=0;				 
    virtual void normalReconstruction(
				      double* mesh_trial_ref,
				      double* mesh_grad_trial_ref,
				      double* mesh_dof,
				      int* mesh_l2g,
				      double* dV_ref,
				      double* u_trial_ref,
				      double* u_grad_trial_ref,
				      double* u_test_ref,
				      int nElements_global,
				      int* u_l2g, 
				      double* elementDiameter,
				      double* phi_dof,
				      int offset_u, int stride_u, 
				      int numDOFs,
				      double* lumped_wx,
				      double* lumped_wy)=0;
				      //double* rhs_mass_correction,
				      //double* lumped_L2p,
				      //double* lumped_mass_matrix,
				      // FOR NONLINEAR VOF; i.e., MCorr with VOF
				      //double epsFactHeaviside,
				      //double epsFactDiffusion,
				      //double* phiHat_dof)=0;           
  };

  template<class CompKernelType,
	   int nSpace,
	   int nQuadraturePoints_element,
	   int nDOF_mesh_trial_element,
	   int nDOF_trial_element,
	   int nDOF_test_element,
	   int nQuadraturePoints_elementBoundary>
  class CLSVOF : public CLSVOF_base
  {
  public:
    const int nDOF_test_X_trial_element;
    CompKernelType ck;
    CLSVOF():
      nDOF_test_X_trial_element(nDOF_test_element*nDOF_trial_element),
      ck()
    {}
    inline
    void evaluateCoefficients(const double v[nSpace],
			      const double& u,
			      const double& porosity, //VRANS specific
			      double& m,
			      double& dm,
			      double f[nSpace],
			      double df[nSpace])
    {
    m = porosity*u;
    dm= porosity;
    for (int I=0; I < nSpace; I++)
      {
	f[I] = v[I]*porosity*u;
	df[I] = v[I]*porosity;
      }
    }

    inline
    void calculateCFL(const double& elementDiameter,
		      const double df[nSpace],
		      double& cfl)
    {
      double h,nrm_v;
      h = elementDiameter;
      nrm_v=0.0;
      for(int I=0;I<nSpace;I++)
	nrm_v+=df[I]*df[I];
      nrm_v = sqrt(nrm_v);
      cfl = nrm_v/h;
    }

    inline
    void calculateSubgridError_tau(const double& elementDiameter,
				   const double& dmt,
				   const double dH[nSpace],
				   double& cfl,
				   double& tau)
    {
      double h,nrm_v,oneByAbsdt;
      h = elementDiameter;
      nrm_v=0.0;
      for(int I=0;I<nSpace;I++)
	nrm_v+=dH[I]*dH[I];
      nrm_v = sqrt(nrm_v);
      cfl = nrm_v/h;
      oneByAbsdt =  fabs(dmt);
      tau = 1.0/(2.0*nrm_v/h + oneByAbsdt + 1.0e-8);
    }

 
    inline
    void calculateSubgridError_tau(     const double&  Ct_sge,
                                        const double   G[nSpace*nSpace],
					const double&  A0,
					const double   Ai[nSpace],
					double& tau_v,
					double& cfl)	
    {
      double v_d_Gv=0.0; 
      for(int I=0;I<nSpace;I++) 
         for (int J=0;J<nSpace;J++) 
           v_d_Gv += Ai[I]*G[I*nSpace+J]*Ai[J];     
    
      tau_v = 1.0/sqrt(Ct_sge*A0*A0 + v_d_Gv + 1.0e-8);    
    } 
 
 

    inline 
    void calculateNumericalDiffusion(const double& shockCapturingDiffusion,
				     const double& elementDiameter,
				     const double& strong_residual,
				     const double grad_u[nSpace],
				     double& numDiff)
    {
      double h,
	num,
	den,
	n_grad_u;
      h = elementDiameter;
      n_grad_u = 0.0;
      for (int I=0;I<nSpace;I++)
	n_grad_u += grad_u[I]*grad_u[I];
      num = shockCapturingDiffusion*0.5*h*fabs(strong_residual);
      den = sqrt(n_grad_u) + 1.0e-8;
      numDiff = num/den;
    }

    inline
    void exteriorNumericalAdvectiveFlux(const int& isDOFBoundary_u,
					const int& isFluxBoundary_u,
					const double n[nSpace],
					const double& bc_u,
					const double& bc_flux_u,
					const double& u,
					const double velocity[nSpace],
					double& flux)
    {

      double flow=0.0;
      for (int I=0; I < nSpace; I++)
	flow += n[I]*velocity[I];
      //std::cout<<" isDOFBoundary_u= "<<isDOFBoundary_u<<" flow= "<<flow<<std::endl;
      if (isDOFBoundary_u == 1)
	{
	  //std::cout<<"Dirichlet boundary u and bc_u "<<u<<'\t'<<bc_u<<std::endl;
	  if (flow >= 0.0)
	    {
	      flux = u*flow;
	      //flux = flow;
	    }
	  else
	    {
	      flux = bc_u*flow;
	      //flux = flow;
	    }
	}
      else if (isFluxBoundary_u == 1)
	{
	  flux = bc_flux_u;
	  //std::cout<<"Flux boundary flux and flow"<<flux<<'\t'<<flow<<std::endl;
	}
      else
	{
	  //std::cout<<"No BC boundary flux and flow"<<flux<<'\t'<<flow<<std::endl;
	  if (flow >= 0.0)
	    {
	      flux = u*flow;
	    }
	  else
	    {
	      std::cout<<"warning: VOF open boundary with no external trace, setting to zero for inflow"<<std::endl;
	      flux = 0.0;
	    }

	}
      //flux = flow;
      //std::cout<<"flux error "<<flux-flow<<std::endl;
      //std::cout<<"flux in computationa"<<flux<<std::endl;
    }

    inline
    void exteriorNumericalAdvectiveFluxDerivative(const int& isDOFBoundary_u,
						  const int& isFluxBoundary_u,
						  const double n[nSpace],
						  const double velocity[nSpace],
						  double& dflux)
    {
      double flow=0.0;
      for (int I=0; I < nSpace; I++)
	flow += n[I]*velocity[I];
      //double flow=n[0]*velocity[0]+n[1]*velocity[1]+n[2]*velocity[2];
      dflux=0.0;//default to no flux
      if (isDOFBoundary_u == 1)
	{
	  if (flow >= 0.0)
	    {
	      dflux = flow;
	    }
	  else
	    {
	      dflux = 0.0;
	    }
	}
      else if (isFluxBoundary_u == 1)
	{
	  dflux = 0.0;
	}
      else
	{
	  if (flow >= 0.0)
	    {
	      dflux = flow;
	    }
	}
    }

    inline void Mult(const double mat[nSpace][nSpace],
		     const double vec[nSpace],
		     double mat_times_vector[nSpace])
    {
      for (int I=0; I<nSpace; I++)
	{
	  mat_times_vector[I] = 0.;
	  for (int J=0; J<nSpace; J++)
	    mat_times_vector[I] += mat[I][J] * vec[J];
	}
    }
    
    inline double smoothedHeaviside(double eps, double phi)
    {
      double H;
      if (phi > eps)
	H=1.0;
      else if (phi < -eps)
	H=0.0;
      else if (phi==0.0)
	H=0.5;
      else
	H = 0.5*(1.0 + phi/eps + sin(M_PI*phi/eps)/M_PI);
      return H; //MQL. USE SIGNED FUNCTION. TMP
      //return 2*H-1; //MQL. USE SIGNED FUNCTION. TMP
    }

    inline double smoothedDirac(double eps, double phi)
    {
      double d;
      if (phi > eps)
	d=0.0;
      else if (phi < -eps)
	d=0.0;
      else
	d = 0.5*(1.0 + cos(M_PI*phi/eps))/eps;
      return d; //MQL. USE SIGNED FUNCTION. TMP
      //return 2*d; //MQL. USE SIGNED FUNCTION. TMP
    }
    
    void calculateResidual_MCorr_with_VOF(//element
					  double dt,
					  double* mesh_trial_ref,
					  double* mesh_grad_trial_ref,
					  double* mesh_dof,
					  double* mesh_velocity_dof,
					  double MOVING_DOMAIN,
					  int* mesh_l2g,
					  double* dV_ref,
					  double* u_trial_ref,
					  double* u_grad_trial_ref,
					  double* u_test_ref,
					  double* u_grad_test_ref,
					  //element boundary
					  double* mesh_trial_trace_ref,
					  double* mesh_grad_trial_trace_ref,
					  double* dS_ref,
					  double* u_trial_trace_ref,
					  double* u_grad_trial_trace_ref,
					  double* u_test_trace_ref,
					  double* u_grad_test_trace_ref,
					  double* normal_ref,
					  double* boundaryJac_ref,
					  //physics
					  int nElements_global,
					  double useMetrics, 
					  double alphaBDF,
					  int lag_shockCapturing, 
					  double shockCapturingDiffusion,
					  double sc_uref, double sc_alpha,
					  //VRANS
					  const double* q_porosity,
					  const double* porosity_dof,
					  //
					  int* u_l2g, 
					  double* elementDiameter,
					  int degree_polynomial,
					  double* u_dof,
					  double* u_dof_old,
					  double* uStar_dof,
					  double* velocity,
					  double* q_m,
					  double* q_u,
					  double* q_m_betaBDF,
					  double* q_dV,
					  double* q_dV_last,
					  double* cfl,
					  double* edge_based_cfl,
					  double* q_numDiff_u, 
					  double* q_numDiff_u_last, 
					  int offset_u, int stride_u, 
					  double* globalResidual,
					  int nExteriorElementBoundaries_global,
					  int* exteriorElementBoundariesArray,
					  int* elementBoundaryElementsArray,
					  int* elementBoundaryLocalElementBoundariesArray,
					  double* ebqe_velocity_ext,
					  //VRANS
					  const double* ebqe_porosity_ext,
					  //
					  int* isDOFBoundary_u,
					  double* ebqe_bc_u_ext,
					  int* isFluxBoundary_u,
					  double* ebqe_bc_flux_u_ext,
					  double* ebqe_phi,double epsFact,
					  double* ebqe_u,
					  double* ebqe_flux,
					  // PARAMETERS FOR EDGE BASED STABILIZATION
					  double cE,
					  double cK,
					  // PARAMETERS FOR LOG BASED ENTROPY FUNCTION 
					  double uL, 
					  double uR,
					  // PARAMETERS FOR EDGE VISCOSITY 
					  int numDOFs,
					  int NNZ,
					  int* csrRowIndeces_DofLoops,
					  int* csrColumnOffsets_DofLoops,
					  int* csrRowIndeces_CellLoops,
					  int* csrColumnOffsets_CellLoops,
					  int* csrColumnOffsets_eb_CellLoops,
					  // C matrices
					  double* Cx, 
					  double* Cy, 
					  double* Cz, 
					  double* CTx,
					  double* CTy, 
					  double* CTz, 
					  double* ML,
					  // PARAMETERS FOR 1st or 2nd ORDER MPP METHOD
					  int LUMPED_MASS_MATRIX, 
					  int STABILIZATION_TYPE,
					  int ENTROPY_TYPE,
					  // FOR FCT
					  double* low_order_solution,
					  double* dt_times_dH_minus_dL,
					  double* min_u_bc,
					  double* max_u_bc,
					  // FOR NONLINEAR VOF; i.e., MCorr with VOF
					  int useFullNewton,
					  double epsFactHeaviside,
					  double epsFactDirac,
					  double epsFactDiffusion,
					  double* phin_dof,
					  double* phiHat_dof,
					  double* lumped_wx,
					  double* lumped_wy,
					  // AUX QUANTITIES OF INTEREST 
					  double* quantDOFs)
    {
      for(int eN=0;eN<nElements_global;eN++)
	{
	  //declare local storage for local contributions and initialize
	  register double elementResidual_u[nDOF_test_element];
	  for (int i=0;i<nDOF_test_element;i++)
	    elementResidual_u[i]=0.0;
	  //loop over quadrature points and compute integrands
	  for  (int k=0;k<nQuadraturePoints_element;k++)
	    {
	      //compute indeces and declare local storage
	      register int eN_k = eN*nQuadraturePoints_element+k,
		eN_k_nSpace = eN_k*nSpace,
		eN_nDOF_trial_element = eN*nDOF_trial_element;
	      register double 
		//for mass matrix contributions
		u, grad_u[nSpace], relative_velocity[nSpace], f[nSpace], //f=velocity*H(phi)
		phiHatnp1, phin,
		u_test_dV[nDOF_trial_element], 
		u_grad_trial[nDOF_trial_element*nSpace], 
		u_grad_test_dV[nDOF_test_element*nSpace],
		//for general use
		jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace],
		dV,x,y,z,xt,yt,zt;
	      //get the physical integration weight
	      ck.calculateMapping_element(eN,
					  k,
					  mesh_dof,
					  mesh_l2g,
					  mesh_trial_ref,
					  mesh_grad_trial_ref,
					  jac,
					  jacDet,
					  jacInv,
					  x,y,z);
	      ck.calculateMappingVelocity_element(eN,
						  k,
						  mesh_velocity_dof,
						  mesh_l2g,
						  mesh_trial_ref,
						  xt,yt,zt);	      
	      dV = fabs(jacDet)*dV_ref[k];
	      //get the solution (of Newton's solver)
	      ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],u);
	      //get the solution gradients at quad points
	      ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,u_grad_trial);
	      ck.gradFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_u);
	      // get phin and phiHatnp1 at quad points
	      ck.valFromDOF(phin_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],phin);
	      ck.valFromDOF(phiHat_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],phiHatnp1);
	      //precalculate test function products with integration weights for mass matrix terms
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
		  for (int I=0;I<nSpace;I++)
		    u_grad_test_dV[j*nSpace+I] = u_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
		}
	      //calculate time derivative at quadrature points
	      if (q_dV_last[eN_k] <= -100)
		q_dV_last[eN_k] = dV;
	      q_dV[eN_k] = dV;
	      //
	      //moving mesh
	      //
	      double mesh_velocity[3];
	      mesh_velocity[0] = xt;
	      mesh_velocity[1] = yt;
	      mesh_velocity[2] = zt;

	      double epsHeaviside = epsFactHeaviside*elementDiameter[eN]; 
	      double epsDiffusion = epsFactDiffusion*elementDiameter[eN]; //kappa = const*h
	      double Hn = smoothedHeaviside(epsHeaviside,phin);
	      //double Hnp1 = smoothedHeaviside(epsHeaviside,phiHatnp1+u);
	      for (int I=0;I<nSpace;I++)
	      {
	        relative_velocity[I] = (velocity[eN_k_nSpace+I]-MOVING_DOMAIN*mesh_velocity[I]);
		f[I] = relative_velocity[I]*Hn;
		//f[I] = relative_velocity[I]*Hnp1; //implicit advection term
	      }
	      //////////////////////////////
	      // CALCULATE CELL BASED CFL //
	      //////////////////////////////
	      calculateCFL(elementDiameter[eN],relative_velocity,cfl[eN_k]); 
	      //double time_derivative_residual = (smoothedHeaviside(epsHeaviside,phiHatnp1+u)-Hn)/dt;
	      double time_derivative_residual = (smoothedHeaviside(epsHeaviside,phiHatnp1+u)-Hn); //TMP
	      //double time_derivative_residual = (smoothedHeaviside(epsHeaviside,
	      //						   phiHatnp1+epsFactDiffusion*u)
	      //				 -Hn)/dt; //TMP

	      //////////////
	      // ith-LOOP //
	      //////////////	      
	      for(int i=0;i<nDOF_test_element;i++) 
		{ 
		  register int i_nSpace=i*nSpace;
		  elementResidual_u[i] += 
		    time_derivative_residual*u_test_dV[i]
		    //+ ck.Advection_weak(f,&u_grad_test_dV[i_nSpace])
		    //+ ck.NumericalDiffusion(epsDiffusion/dt,grad_u,&u_grad_test_dV[i_nSpace]);
		    + ck.NumericalDiffusion(epsDiffusion,grad_u,&u_grad_test_dV[i_nSpace]); //TMP
		  //+ ck.NumericalDiffusion(1.0,grad_u,&u_grad_test_dV[i_nSpace]); //TMP
		}//i
	      //save solution for other models 
	      q_u[eN_k] = u;
	      q_m[eN_k] = u;//porosity*u;

	    }
	  /////////////////
	  // DISTRIBUTE // load cell based element into global residual
	  ////////////////
	  for(int i=0;i<nDOF_test_element;i++) 
	    { 
	      int eN_i=eN*nDOF_test_element+i;
	      int gi = offset_u+stride_u*u_l2g[eN_i]; //global i-th index
	      // distribute global residual for (lumped) mass matrix
	      globalResidual[gi] += elementResidual_u[i];
	    }//i
	}//elements
      //////////////
      // BOUNDARY //
      //////////////
      if(false)
      for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) 
	{ 
	  register int ebN = exteriorElementBoundariesArray[ebNE]; 
	  register int eN  = elementBoundaryElementsArray[ebN*2+0],
	    ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+0],
	    eN_nDOF_trial_element = eN*nDOF_trial_element;
	  register double elementResidual_u[nDOF_test_element];
	  for (int i=0;i<nDOF_test_element;i++)
	    elementResidual_u[i]=0.0;
	  // loop on quad points
	  for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) 
	    { 
	      register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
		ebNE_kb_nSpace = ebNE_kb*nSpace,
		ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb;
	      register double 
		phin=0.0, 
		relative_velocity[nSpace],
		jac_ext[nSpace*nSpace],
		jacDet_ext,
		jacInv_ext[nSpace*nSpace],
		boundaryJac[nSpace*(nSpace-1)],
		metricTensor[(nSpace-1)*(nSpace-1)],
		metricTensorDetSqrt,
		dS,
		u_test_dS[nDOF_test_element],
		normal[nSpace],
		x_ext,y_ext,z_ext,xt_ext,yt_ext,zt_ext,integralScaling,porosity_ext;
	      // calculate mappings 
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
						  x_ext,y_ext,z_ext);
	      ck.calculateMappingVelocity_elementBoundary(eN,
							  ebN_local,
							  kb,
							  ebN_local_kb,
							  mesh_velocity_dof,
							  mesh_l2g,
							  mesh_trial_trace_ref,
							  xt_ext,yt_ext,zt_ext,
							  normal,
							  boundaryJac,
							  metricTensor,
							  integralScaling);
	      dS = ((1.0-MOVING_DOMAIN)*metricTensorDetSqrt +
		    MOVING_DOMAIN*integralScaling)*dS_ref[kb];
	      //compute shape and solution information
	      ck.valFromDOF(phin_dof,&u_l2g[eN_nDOF_trial_element],
			    &u_trial_trace_ref[ebN_local_kb*nDOF_test_element],phin);
	      //precalculate test function products with integration weights
	      for (int j=0;j<nDOF_trial_element;j++)
		u_test_dS[j] = u_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;
	      //
	      //moving mesh
	      //
	      double mesh_velocity[3];
	      mesh_velocity[0] = xt_ext;
	      mesh_velocity[1] = yt_ext;
	      mesh_velocity[2] = zt_ext;
	      //std::cout<<"mesh_velocity ext"<<std::endl;
	      for (int I=0;I<nSpace;I++)
		relative_velocity[I] = (ebqe_velocity_ext[ebNE_kb_nSpace+I]
					- MOVING_DOMAIN*mesh_velocity[I]);
	      double flow = 0.;
	      for (int I=0; I < nSpace; I++)
		flow += normal[I]*relative_velocity[I];
	      
	      double epsHeaviside = epsFactHeaviside*elementDiameter[eN]; 
	      double Hn = smoothedHeaviside(epsHeaviside,phin);
	      for (int i=0;i<nDOF_trial_element;i++)
		elementResidual_u[i] += flow*Hn*u_test_dS[i];
	    }//kb	  
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      int eN_i = eN*nDOF_test_element+i;
	      int gi = offset_u+stride_u*u_l2g[eN_i]; //global i-th index
	      globalResidual[gi] += elementResidual_u[i];
	    }
	}//ebNE
      // END OF BOUNDARY //
    }

    void calculateResidual_MCorr_with_VOF2(//element
					  double dt,
					  double* mesh_trial_ref,
					  double* mesh_grad_trial_ref,
					  double* mesh_dof,
					  double* mesh_velocity_dof,
					  double MOVING_DOMAIN,
					  int* mesh_l2g,
					  double* dV_ref,
					  double* u_trial_ref,
					  double* u_grad_trial_ref,
					  double* u_test_ref,
					  double* u_grad_test_ref,
					  //element boundary
					  double* mesh_trial_trace_ref,
					  double* mesh_grad_trial_trace_ref,
					  double* dS_ref,
					  double* u_trial_trace_ref,
					  double* u_grad_trial_trace_ref,
					  double* u_test_trace_ref,
					  double* u_grad_test_trace_ref,
					  double* normal_ref,
					  double* boundaryJac_ref,
					  //physics
					  int nElements_global,
					  double useMetrics, 
					  double alphaBDF,
					  int lag_shockCapturing, 
					  double shockCapturingDiffusion,
					  double sc_uref, double sc_alpha,
					  //VRANS
					  const double* q_porosity,
					  const double* porosity_dof,
					  //
					  int* u_l2g, 
					  double* elementDiameter,
					  int degree_polynomial,
					  double* u_dof,
					  double* u_dof_old,
					  double* uStar_dof,
					  double* velocity,
					  double* q_m,
					  double* q_u,
					  double* q_m_betaBDF,
					  double* q_dV,
					  double* q_dV_last,
					  double* cfl,
					  double* edge_based_cfl,
					  double* q_numDiff_u, 
					  double* q_numDiff_u_last, 
					  int offset_u, int stride_u, 
					  double* globalResidual,
					  int nExteriorElementBoundaries_global,
					  int* exteriorElementBoundariesArray,
					  int* elementBoundaryElementsArray,
					  int* elementBoundaryLocalElementBoundariesArray,
					  double* ebqe_velocity_ext,
					  //VRANS
					  const double* ebqe_porosity_ext,
					  //
					  int* isDOFBoundary_u,
					  double* ebqe_bc_u_ext,
					  int* isFluxBoundary_u,
					  double* ebqe_bc_flux_u_ext,
					  double* ebqe_phi,double epsFact,
					  double* ebqe_u,
					  double* ebqe_flux,
					  // PARAMETERS FOR EDGE BASED STABILIZATION
					  double cE,
					  double cK,
					  // PARAMETERS FOR LOG BASED ENTROPY FUNCTION 
					  double uL, 
					  double uR,
					  // PARAMETERS FOR EDGE VISCOSITY 
					  int numDOFs,
					  int NNZ,
					  int* csrRowIndeces_DofLoops,
					  int* csrColumnOffsets_DofLoops,
					  int* csrRowIndeces_CellLoops,
					  int* csrColumnOffsets_CellLoops,
					  int* csrColumnOffsets_eb_CellLoops,
					  // C matrices
					  double* Cx, 
					  double* Cy, 
					  double* Cz, 
					  double* CTx,
					  double* CTy, 
					  double* CTz, 
					  double* ML,
					  // PARAMETERS FOR 1st or 2nd ORDER MPP METHOD
					  int LUMPED_MASS_MATRIX, 
					  int STABILIZATION_TYPE,
					  int ENTROPY_TYPE,
					  // FOR FCT
					  double* low_order_solution,
					  double* dt_times_dH_minus_dL,
					  double* min_u_bc,
					  double* max_u_bc,
					  // FOR NONLINEAR VOF; i.e., MCorr with VOF
					  int useFullNewton,
					  double epsFactHeaviside,
					  double epsFactDirac,
					  double epsFactDiffusion,
					  double* phin_dof,
					  double* phiHat_dof,
					  double* lumped_wx,
					  double* lumped_wy,
					  // AUX QUANTITIES OF INTEREST 
					  double* quantDOFs)
    {
      for(int eN=0;eN<nElements_global;eN++)
	{
	  //declare local storage for local contributions and initialize
	  register double elementResidual_u[nDOF_test_element];
	  for (int i=0;i<nDOF_test_element;i++)
	    elementResidual_u[i]=0.0;
	  //loop over quadrature points and compute integrands
	  for  (int k=0;k<nQuadraturePoints_element;k++)
	    {
	      //compute indeces and declare local storage
	      register int eN_k = eN*nQuadraturePoints_element+k,
		eN_k_nSpace = eN_k*nSpace,
		eN_nDOF_trial_element = eN*nDOF_trial_element;
	      register double 
		//for mass matrix contributions
		u, grad_u[nSpace],
		un, grad_un[nSpace], uStar, grad_uStar[nSpace],
		relative_velocity[nSpace], f[nSpace], //f=velocity*H(phi)
		phiHatnp1, phin,
		u_test_dV[nDOF_trial_element], 
		u_grad_trial[nDOF_trial_element*nSpace], 
		u_grad_test_dV[nDOF_test_element*nSpace],
		//for general use
		jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace],
		dV,x,y,z,xt,yt,zt;
	      //get the physical integration weight
	      ck.calculateMapping_element(eN,
					  k,
					  mesh_dof,
					  mesh_l2g,
					  mesh_trial_ref,
					  mesh_grad_trial_ref,
					  jac,
					  jacDet,
					  jacInv,
					  x,y,z);
	      ck.calculateMappingVelocity_element(eN,
						  k,
						  mesh_velocity_dof,
						  mesh_l2g,
						  mesh_trial_ref,
						  xt,yt,zt);	      
	      dV = fabs(jacDet)*dV_ref[k];
	      ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],
				  jacInv,u_grad_trial);
	      
	      // get the solution (of Newton's solver)
	      ck.valFromDOF(u_dof,
			    &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			    u);
	      // get old solutions
	      ck.valFromDOF(u_dof_old,
			    &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			    un);
	      ck.valFromDOF(uStar_dof,
			    &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			    uStar);
	      //get the solution gradients at quad points	      
	      ck.gradFromDOF(u_dof,
			     &u_l2g[eN_nDOF_trial_element],u_grad_trial,
			     grad_u);
	      ck.gradFromDOF(u_dof_old,
			     &u_l2g[eN_nDOF_trial_element],u_grad_trial,
			     grad_un);
	      ck.gradFromDOF(uStar_dof,
			     &u_l2g[eN_nDOF_trial_element],u_grad_trial,
			     grad_uStar);
	      //precalculate test function products with integration weights for mass matrix terms
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
		  for (int I=0;I<nSpace;I++)
		    u_grad_test_dV[j*nSpace+I] = u_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
		}
	      //calculate time derivative at quadrature points
	      if (q_dV_last[eN_k] <= -100)
		q_dV_last[eN_k] = dV;
	      q_dV[eN_k] = dV;
	      //
	      //moving mesh
	      //
	      double mesh_velocity[3];
	      mesh_velocity[0] = xt;
	      mesh_velocity[1] = yt;
	      mesh_velocity[2] = zt;

	      double gradu2 = 0;
	      for(int I=0;I<nSpace;I++)
		gradu2 += grad_u[I]*grad_u[I];
	      double beta_norm_grad_u = betaNormGrad(gradu2,1.E-10);

	      double gradun2 = 0;
	      for(int I=0;I<nSpace;I++)
		gradun2 += grad_un[I]*grad_un[I];
	      double beta_norm_grad_un = betaNormGrad(gradun2,1.E-10);

	      double graduStar2 = 0;
	      for(int I=0;I<nSpace;I++)
		graduStar2 += grad_uStar[I]*grad_uStar[I];
	      double beta_norm_grad_uStar = betaNormGrad(graduStar2,1.E-10);
	      
	      double lambda = epsFactDiffusion;	      
	      //double lambda = epsFactDiffusion*elementDiameter[eN]/dt;

	      //double coeffFullNewton = -1./beta_norm_grad_u; // single potential
	      double coeffFullNewton = 2*std::pow(beta_norm_grad_u,2)-3*beta_norm_grad_u; //2 pot.
	      
	      //double coeffLinNewton = -1./beta_norm_grad_un; // single potential
	      double coeffLinNewton = 2*std::pow(beta_norm_grad_un,2)-3*beta_norm_grad_un; //2 pot.
 	      //-3+3./2*beta_norm_grad_un+0.5/beta_norm_grad_un;

	      double epsHeaviside = epsFactHeaviside*elementDiameter[eN]; 	    
	      double Hn = smoothedHeaviside(epsHeaviside,un);
	      for (int I=0;I<nSpace;I++)
	      {
	        relative_velocity[I] = (velocity[eN_k_nSpace+I]-MOVING_DOMAIN*mesh_velocity[I]);
		f[I] = relative_velocity[I]*Hn;
		//f[I] = relative_velocity[I]*un;
	      }
	      //////////////////////////////
	      // CALCULATE CELL BASED CFL //
	      //////////////////////////////
	      calculateCFL(elementDiameter[eN],relative_velocity,cfl[eN_k]); 
	      double time_derivative_residual = (smoothedHeaviside(epsHeaviside,u)-Hn)/dt;
	      //double time_derivative_residual = (u-un)/dt;
	      //////////////
	      // ith-LOOP //
	      //////////////
	      if (useFullNewton==false)
		for(int i=0;i<nDOF_test_element;i++) 
		  { 
		    register int i_nSpace=i*nSpace;
		    elementResidual_u[i] += 
		      time_derivative_residual*u_test_dV[i]
		      + ck.Advection_weak(f,&u_grad_test_dV[i_nSpace])
		      + lambda*ck.NumericalDiffusion(1.0,grad_u,&u_grad_test_dV[i_nSpace])
		      + lambda*ck.NumericalDiffusion(coeffLinNewton,
						     grad_un,&u_grad_test_dV[i_nSpace]);
		  }//i
	      else
		{
		  for(int i=0;i<nDOF_test_element;i++) 
		    { 
		      register int i_nSpace=i*nSpace;
		      elementResidual_u[i] += 
			time_derivative_residual*u_test_dV[i]
			+ ck.Advection_weak(f,&u_grad_test_dV[i_nSpace])
			+ lambda*ck.NumericalDiffusion(1.0,grad_u,&u_grad_test_dV[i_nSpace])
			+ lambda*ck.NumericalDiffusion(coeffFullNewton,
						       grad_u,&u_grad_test_dV[i_nSpace]);
		    }//i
		}
	      //save solution for other models 
	      q_u[eN_k] = u;
	      q_m[eN_k] = u;//porosity*u;
	    }
	  /////////////////
	  // DISTRIBUTE // load cell based element into global residual
	  ////////////////
	  for(int i=0;i<nDOF_test_element;i++) 
	    { 
	      int eN_i=eN*nDOF_test_element+i;
	      int gi = offset_u+stride_u*u_l2g[eN_i]; //global i-th index
	      // distribute global residual for (lumped) mass matrix
	      globalResidual[gi] += elementResidual_u[i];
	    }//i
	}//elements
      
      //////////////
      // BOUNDARY //
      //////////////
      if(false)
      for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) 
	{ 
	  register int ebN = exteriorElementBoundariesArray[ebNE]; 
	  register int eN  = elementBoundaryElementsArray[ebN*2+0],
	    ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+0],
	    eN_nDOF_trial_element = eN*nDOF_trial_element;
	  register double elementResidual_u[nDOF_test_element];
	  for (int i=0;i<nDOF_test_element;i++)
	    elementResidual_u[i]=0.0;
	  // loop on quad points
	  for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) 
	    { 
	      register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
		ebNE_kb_nSpace = ebNE_kb*nSpace,
		ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb;
	      register double 
		phin=0.0, 
		relative_velocity[nSpace],
		jac_ext[nSpace*nSpace],
		jacDet_ext,
		jacInv_ext[nSpace*nSpace],
		boundaryJac[nSpace*(nSpace-1)],
		metricTensor[(nSpace-1)*(nSpace-1)],
		metricTensorDetSqrt,
		dS,
		u_test_dS[nDOF_test_element],
		normal[nSpace],
		x_ext,y_ext,z_ext,xt_ext,yt_ext,zt_ext,integralScaling,porosity_ext;
	      // calculate mappings 
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
						  x_ext,y_ext,z_ext);
	      ck.calculateMappingVelocity_elementBoundary(eN,
							  ebN_local,
							  kb,
							  ebN_local_kb,
							  mesh_velocity_dof,
							  mesh_l2g,
							  mesh_trial_trace_ref,
							  xt_ext,yt_ext,zt_ext,
							  normal,
							  boundaryJac,
							  metricTensor,
							  integralScaling);
	      dS = ((1.0-MOVING_DOMAIN)*metricTensorDetSqrt +
		    MOVING_DOMAIN*integralScaling)*dS_ref[kb];
	      //compute shape and solution information
	      ck.valFromDOF(phin_dof,&u_l2g[eN_nDOF_trial_element],
			    &u_trial_trace_ref[ebN_local_kb*nDOF_test_element],phin);
	      //precalculate test function products with integration weights
	      for (int j=0;j<nDOF_trial_element;j++)
		u_test_dS[j] = u_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;
	      //
	      //moving mesh
	      //
	      double mesh_velocity[3];
	      mesh_velocity[0] = xt_ext;
	      mesh_velocity[1] = yt_ext;
	      mesh_velocity[2] = zt_ext;
	      //std::cout<<"mesh_velocity ext"<<std::endl;
	      for (int I=0;I<nSpace;I++)
		relative_velocity[I] = (ebqe_velocity_ext[ebNE_kb_nSpace+I]
					- MOVING_DOMAIN*mesh_velocity[I]);
	      double flow = 0.;
	      for (int I=0; I < nSpace; I++)
		flow += normal[I]*relative_velocity[I];
	      
	      double epsHeaviside = epsFactHeaviside*elementDiameter[eN]; 
	      double Hn = smoothedHeaviside(epsHeaviside,phin);
	      for (int i=0;i<nDOF_trial_element;i++)
		elementResidual_u[i] += flow*Hn*u_test_dS[i];
	    }//kb	  
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      int eN_i = eN*nDOF_test_element+i;
	      int gi = offset_u+stride_u*u_l2g[eN_i]; //global i-th index
	      globalResidual[gi] += elementResidual_u[i];
	    }
	}//ebNE
      // END OF BOUNDARY //
    }

    void calculateResidual_MCorr_with_VOF3(//element
					  double dt,
					  double* mesh_trial_ref,
					  double* mesh_grad_trial_ref,
					  double* mesh_dof,
					  double* mesh_velocity_dof,
					  double MOVING_DOMAIN,
					  int* mesh_l2g,
					  double* dV_ref,
					  double* u_trial_ref,
					  double* u_grad_trial_ref,
					  double* u_test_ref,
					  double* u_grad_test_ref,
					  //element boundary
					  double* mesh_trial_trace_ref,
					  double* mesh_grad_trial_trace_ref,
					  double* dS_ref,
					  double* u_trial_trace_ref,
					  double* u_grad_trial_trace_ref,
					  double* u_test_trace_ref,
					  double* u_grad_test_trace_ref,
					  double* normal_ref,
					  double* boundaryJac_ref,
					  //physics
					  int nElements_global,
					  double useMetrics, 
					  double alphaBDF,
					  int lag_shockCapturing, 
					  double shockCapturingDiffusion,
					  double sc_uref, double sc_alpha,
					  //VRANS
					  const double* q_porosity,
					  const double* porosity_dof,
					  //
					  int* u_l2g, 
					  double* elementDiameter,
					  int degree_polynomial,
					  double* u_dof,
					  double* u_dof_old,
					  double* uStar_dof,
					  double* velocity,
					  double* q_m,
					  double* q_u,
					  double* q_m_betaBDF,
					  double* q_dV,
					  double* q_dV_last,
					  double* cfl,
					  double* edge_based_cfl,
					  double* q_numDiff_u, 
					  double* q_numDiff_u_last, 
					  int offset_u, int stride_u, 
					  double* globalResidual,
					  int nExteriorElementBoundaries_global,
					  int* exteriorElementBoundariesArray,
					  int* elementBoundaryElementsArray,
					  int* elementBoundaryLocalElementBoundariesArray,
					  double* ebqe_velocity_ext,
					  //VRANS
					  const double* ebqe_porosity_ext,
					  //
					  int* isDOFBoundary_u,
					  double* ebqe_bc_u_ext,
					  int* isFluxBoundary_u,
					  double* ebqe_bc_flux_u_ext,
					  double* ebqe_phi,double epsFact,
					  double* ebqe_u,
					  double* ebqe_flux,
					  // PARAMETERS FOR EDGE BASED STABILIZATION
					  double cE,
					  double cK,
					  // PARAMETERS FOR LOG BASED ENTROPY FUNCTION 
					  double uL, 
					  double uR,
					  // PARAMETERS FOR EDGE VISCOSITY 
					  int numDOFs,
					  int NNZ,
					  int* csrRowIndeces_DofLoops,
					  int* csrColumnOffsets_DofLoops,
					  int* csrRowIndeces_CellLoops,
					  int* csrColumnOffsets_CellLoops,
					  int* csrColumnOffsets_eb_CellLoops,
					  // C matrices
					  double* Cx, 
					  double* Cy, 
					  double* Cz, 
					  double* CTx,
					  double* CTy, 
					  double* CTz, 
					  double* ML,
					  // PARAMETERS FOR 1st or 2nd ORDER MPP METHOD
					  int LUMPED_MASS_MATRIX, 
					  int STABILIZATION_TYPE,
					  int ENTROPY_TYPE,
					  // FOR FCT
					  double* low_order_solution,
					  double* dt_times_dH_minus_dL,
					  double* min_u_bc,
					  double* max_u_bc,
					  // FOR NONLINEAR VOF; i.e., MCorr with VOF
					  int useFullNewton,
					  double epsFactHeaviside,
					  double epsFactDirac,
					  double epsFactDiffusion,
					  double* phin_dof,
					  double* phiHat_dof,
					  // normal reconstruction
					  double* lumped_wx,
					  double* lumped_wy,
					  // AUX QUANTITIES OF INTEREST 
					  double* quantDOFs)
    {
      for(int eN=0;eN<nElements_global;eN++)
	{
	  //declare local storage for local contributions and initialize
	  register double elementResidual_u[nDOF_test_element];
	  for (int i=0;i<nDOF_test_element;i++)
	    elementResidual_u[i]=0.0;
	  //loop over quadrature points and compute integrands
	  for  (int k=0;k<nQuadraturePoints_element;k++)
	    {
	      //compute indeces and declare local storage
	      register int eN_k = eN*nQuadraturePoints_element+k,
		eN_k_nSpace = eN_k*nSpace,
		eN_nDOF_trial_element = eN*nDOF_trial_element;
	      register double 
		//for mass matrix contributions
		u, grad_u[nSpace],
		normalReconstruction[nSpace], wx, wy,
		un, grad_un[nSpace], uStar, grad_uStar[nSpace],
		relative_velocity[nSpace], f[nSpace], //f=velocity*H(phi)
		phiHatnp1, phin,
		u_test_dV[nDOF_trial_element], 
		u_grad_trial[nDOF_trial_element*nSpace], 
		u_grad_test_dV[nDOF_test_element*nSpace],
		//for general use
		jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace],
		dV,x,y,z,xt,yt,zt;
	      //get the physical integration weight
	      ck.calculateMapping_element(eN,
					  k,
					  mesh_dof,
					  mesh_l2g,
					  mesh_trial_ref,
					  mesh_grad_trial_ref,
					  jac,
					  jacDet,
					  jacInv,
					  x,y,z);
	      ck.calculateMappingVelocity_element(eN,
						  k,
						  mesh_velocity_dof,
						  mesh_l2g,
						  mesh_trial_ref,
						  xt,yt,zt);	      
	      dV = fabs(jacDet)*dV_ref[k];
	      ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],
				  jacInv,u_grad_trial);
	      // get the components of the normal reconstruction
	      ck.valFromDOF(lumped_wx,
			    &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			    wx);
	      ck.valFromDOF(lumped_wy,
			    &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			    wy);
	      // get the solution (of Newton's solver)
	      ck.valFromDOF(u_dof,
			    &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			    u);
	      // get old solutions
	      ck.valFromDOF(u_dof_old,
			    &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			    un);
	      ck.valFromDOF(uStar_dof,
			    &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			    uStar);
	      //get the solution gradients at quad points	      
	      ck.gradFromDOF(u_dof,
			     &u_l2g[eN_nDOF_trial_element],u_grad_trial,
			     grad_u);
	      ck.gradFromDOF(u_dof_old,
			     &u_l2g[eN_nDOF_trial_element],u_grad_trial,
			     grad_un);
	      ck.gradFromDOF(uStar_dof,
			     &u_l2g[eN_nDOF_trial_element],u_grad_trial,
			     grad_uStar);
	      //precalculate test function products with integration weights for mass matrix terms
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
		  for (int I=0;I<nSpace;I++)
		    u_grad_test_dV[j*nSpace+I] = u_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
		}
	      //calculate time derivative at quadrature points
	      if (q_dV_last[eN_k] <= -100)
		q_dV_last[eN_k] = dV;
	      q_dV[eN_k] = dV;
	      //
	      //moving mesh
	      //
	      double mesh_velocity[3];
	      mesh_velocity[0] = xt;
	      mesh_velocity[1] = yt;
	      mesh_velocity[2] = zt;

	      double gradu2 = 0;
	      for(int I=0;I<nSpace;I++)
		gradu2 += grad_u[I]*grad_u[I];
	      double beta_norm_grad_u = betaNormGrad(gradu2,1.E-10);

	      double gradun2 = 0;
	      for(int I=0;I<nSpace;I++)
		gradun2 += grad_un[I]*grad_un[I];
	      double beta_norm_grad_un = betaNormGrad(gradun2,1.E-10);

	      double graduStar2 = 0;
	      for(int I=0;I<nSpace;I++)
		graduStar2 += grad_uStar[I]*grad_uStar[I];
	      double beta_norm_grad_uStar = betaNormGrad(graduStar2,1.E-10);
	      
	      double lambda = epsFactDiffusion;	      
	      //double lambda = epsFactDiffusion*elementDiameter[eN]/dt;

	      //double coeffFullNewton = -1./beta_norm_grad_u; // single potential
	      double coeffFullNewton = 2*std::pow(beta_norm_grad_u,2)-3*beta_norm_grad_u; //2 pot.
	      
	      //double coeffLinNewton = -1./beta_norm_grad_un; // single potential
	      double coeffLinNewton = 2*std::pow(beta_norm_grad_un,2)-3*beta_norm_grad_un; //2 pot.
 	      //-3+3./2*beta_norm_grad_un+0.5/beta_norm_grad_un;

	      double epsHeaviside = epsFactHeaviside*elementDiameter[eN]; 	    
	      double Hn = smoothedHeaviside(epsHeaviside,un);
	      double Hnp1 = smoothedHeaviside(epsHeaviside,u);
	      for (int I=0;I<nSpace;I++)
	      {
	        relative_velocity[I] = (velocity[eN_k_nSpace+I]-MOVING_DOMAIN*mesh_velocity[I]);
		//f[I] = relative_velocity[I]*Hn;
		f[I] = relative_velocity[I]*Hnp1; //implicit advection via BDF
		//f[I] = 0.5*relative_velocity[I]*(Hnp1+Hn); //implicit advection via CN
	      }
	      //////////////////////////////
	      // CALCULATE CELL BASED CFL //
	      //////////////////////////////
	      calculateCFL(elementDiameter[eN],relative_velocity,cfl[eN_k]); 
	      double time_derivative_residual = (smoothedHeaviside(epsHeaviside,u)-Hn)/dt;
	      //double time_derivative_residual = (u-un)/dt;
	      //////////////
	      // ith-LOOP //
	      //////////////
	      normalReconstruction[0] = wx;
	      normalReconstruction[1] = wy;
	      for(int i=0;i<nDOF_test_element;i++) 
		{ 
		  register int i_nSpace=i*nSpace;
		  elementResidual_u[i] += 
		    time_derivative_residual*u_test_dV[i]
		    + ck.Advection_weak(f,&u_grad_test_dV[i_nSpace])
		    + lambda*ck.NumericalDiffusion(1.0,grad_u,&u_grad_test_dV[i_nSpace])
		    + lambda*ck.NumericalDiffusion(-1.0,
						   normalReconstruction,
						   &u_grad_test_dV[i_nSpace]);
		}//i
	      //save solution for other models 
	      q_u[eN_k] = u;
	      q_m[eN_k] = u;//porosity*u;
	    }
	  /////////////////
	  // DISTRIBUTE // load cell based element into global residual
	  ////////////////
	  for(int i=0;i<nDOF_test_element;i++) 
	    { 
	      int eN_i=eN*nDOF_test_element+i;
	      int gi = offset_u+stride_u*u_l2g[eN_i]; //global i-th index
	      // distribute global residual for (lumped) mass matrix
	      globalResidual[gi] += elementResidual_u[i];
	    }//i
	}//elements
      
      //////////////
      // BOUNDARY //
      //////////////
      if(false)
      for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) 
	{ 
	  register int ebN = exteriorElementBoundariesArray[ebNE]; 
	  register int eN  = elementBoundaryElementsArray[ebN*2+0],
	    ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+0],
	    eN_nDOF_trial_element = eN*nDOF_trial_element;
	  register double elementResidual_u[nDOF_test_element];
	  for (int i=0;i<nDOF_test_element;i++)
	    elementResidual_u[i]=0.0;
	  // loop on quad points
	  for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) 
	    { 
	      register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
		ebNE_kb_nSpace = ebNE_kb*nSpace,
		ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb;
	      register double 
		phin=0.0, 
		relative_velocity[nSpace],
		jac_ext[nSpace*nSpace],
		jacDet_ext,
		jacInv_ext[nSpace*nSpace],
		boundaryJac[nSpace*(nSpace-1)],
		metricTensor[(nSpace-1)*(nSpace-1)],
		metricTensorDetSqrt,
		dS,
		u_test_dS[nDOF_test_element],
		normal[nSpace],
		x_ext,y_ext,z_ext,xt_ext,yt_ext,zt_ext,integralScaling,porosity_ext;
	      // calculate mappings 
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
						  x_ext,y_ext,z_ext);
	      ck.calculateMappingVelocity_elementBoundary(eN,
							  ebN_local,
							  kb,
							  ebN_local_kb,
							  mesh_velocity_dof,
							  mesh_l2g,
							  mesh_trial_trace_ref,
							  xt_ext,yt_ext,zt_ext,
							  normal,
							  boundaryJac,
							  metricTensor,
							  integralScaling);
	      dS = ((1.0-MOVING_DOMAIN)*metricTensorDetSqrt +
		    MOVING_DOMAIN*integralScaling)*dS_ref[kb];
	      //compute shape and solution information
	      ck.valFromDOF(phin_dof,&u_l2g[eN_nDOF_trial_element],
			    &u_trial_trace_ref[ebN_local_kb*nDOF_test_element],phin);
	      //precalculate test function products with integration weights
	      for (int j=0;j<nDOF_trial_element;j++)
		u_test_dS[j] = u_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;
	      //
	      //moving mesh
	      //
	      double mesh_velocity[3];
	      mesh_velocity[0] = xt_ext;
	      mesh_velocity[1] = yt_ext;
	      mesh_velocity[2] = zt_ext;
	      //std::cout<<"mesh_velocity ext"<<std::endl;
	      for (int I=0;I<nSpace;I++)
		relative_velocity[I] = (ebqe_velocity_ext[ebNE_kb_nSpace+I]
					- MOVING_DOMAIN*mesh_velocity[I]);
	      double flow = 0.;
	      for (int I=0; I < nSpace; I++)
		flow += normal[I]*relative_velocity[I];
	      
	      double epsHeaviside = epsFactHeaviside*elementDiameter[eN]; 
	      double Hn = smoothedHeaviside(epsHeaviside,phin);
	      for (int i=0;i<nDOF_trial_element;i++)
		elementResidual_u[i] += flow*Hn*u_test_dS[i];
	    }//kb	  
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      int eN_i = eN*nDOF_test_element+i;
	      int gi = offset_u+stride_u*u_l2g[eN_i]; //global i-th index
	      globalResidual[gi] += elementResidual_u[i];
	    }
	}//ebNE
      // END OF BOUNDARY //
    }

    void calculateResidual_MCorr_with_VOF4(//element
					  double dt,
					  double* mesh_trial_ref,
					  double* mesh_grad_trial_ref,
					  double* mesh_dof,
					  double* mesh_velocity_dof,
					  double MOVING_DOMAIN,
					  int* mesh_l2g,
					  double* dV_ref,
					  double* u_trial_ref,
					  double* u_grad_trial_ref,
					  double* u_test_ref,
					  double* u_grad_test_ref,
					  //element boundary
					  double* mesh_trial_trace_ref,
					  double* mesh_grad_trial_trace_ref,
					  double* dS_ref,
					  double* u_trial_trace_ref,
					  double* u_grad_trial_trace_ref,
					  double* u_test_trace_ref,
					  double* u_grad_test_trace_ref,
					  double* normal_ref,
					  double* boundaryJac_ref,
					  //physics
					  int nElements_global,
					  double useMetrics, 
					  double alphaBDF,
					  int lag_shockCapturing, 
					  double shockCapturingDiffusion,
					  double sc_uref, double sc_alpha,
					  //VRANS
					  const double* q_porosity,
					  const double* porosity_dof,
					  //
					  int* u_l2g, 
					  double* elementDiameter,
					  int degree_polynomial,
					  double* u_dof,
					  double* u_dof_old,
					  double* uStar_dof,
					  double* velocity,
					  double* q_m,
					  double* q_u,
					  double* q_m_betaBDF,
					  double* q_dV,
					  double* q_dV_last,
					  double* cfl,
					  double* edge_based_cfl,
					  double* q_numDiff_u, 
					  double* q_numDiff_u_last, 
					  int offset_u, int stride_u, 
					  double* globalResidual,
					  int nExteriorElementBoundaries_global,
					  int* exteriorElementBoundariesArray,
					  int* elementBoundaryElementsArray,
					  int* elementBoundaryLocalElementBoundariesArray,
					  double* ebqe_velocity_ext,
					  //VRANS
					  const double* ebqe_porosity_ext,
					  //
					  int* isDOFBoundary_u,
					  double* ebqe_bc_u_ext,
					  int* isFluxBoundary_u,
					  double* ebqe_bc_flux_u_ext,
					  double* ebqe_phi,double epsFact,
					  double* ebqe_u,
					  double* ebqe_flux,
					  // PARAMETERS FOR EDGE BASED STABILIZATION
					  double cE,
					  double cK,
					  // PARAMETERS FOR LOG BASED ENTROPY FUNCTION 
					  double uL, 
					  double uR,
					  // PARAMETERS FOR EDGE VISCOSITY 
					  int numDOFs,
					  int NNZ,
					  int* csrRowIndeces_DofLoops,
					  int* csrColumnOffsets_DofLoops,
					  int* csrRowIndeces_CellLoops,
					  int* csrColumnOffsets_CellLoops,
					  int* csrColumnOffsets_eb_CellLoops,
					  // C matrices
					  double* Cx, 
					  double* Cy, 
					  double* Cz, 
					  double* CTx,
					  double* CTy, 
					  double* CTz, 
					  double* ML,
					  // PARAMETERS FOR 1st or 2nd ORDER MPP METHOD
					  int LUMPED_MASS_MATRIX, 
					  int STABILIZATION_TYPE,
					  int ENTROPY_TYPE,
					  // FOR FCT
					  double* low_order_solution,
					  double* dt_times_dH_minus_dL,
					  double* min_u_bc,
					  double* max_u_bc,
					  // FOR NONLINEAR VOF; i.e., MCorr with VOF
					  int useFullNewton,
					  double epsFactHeaviside,
					  double epsFactDirac,
					  double epsFactDiffusion,
					  double* phin_dof,
					  double* phiHat_dof,
					  // normal reconstruction
					  double* lumped_wx,
					  double* lumped_wy,
					  // AUX QUANTITIES OF INTEREST 
					  double* quantDOFs)
    {
      for(int eN=0;eN<nElements_global;eN++)
	{
	  //declare local storage for local contributions and initialize
	  register double elementResidual_u[nDOF_test_element];
	  for (int i=0;i<nDOF_test_element;i++)
	    elementResidual_u[i]=0.0;
	  //loop over quadrature points and compute integrands
	  for  (int k=0;k<nQuadraturePoints_element;k++)
	    {
	      //compute indeces and declare local storage
	      register int eN_k = eN*nQuadraturePoints_element+k,
		eN_k_nSpace = eN_k*nSpace,
		eN_nDOF_trial_element = eN*nDOF_trial_element;
	      register double 
		//for mass matrix contributions
		u, grad_u[nSpace], proj_times_grad_u[nSpace],
		normalReconstruction[nSpace], wx, wy,
		un, grad_un[nSpace], uStar, grad_uStar[nSpace],
		relative_velocity[nSpace], f[nSpace], //f=velocity*H(phi)
		phiHatnp1, phin,
		u_test_dV[nDOF_trial_element], 
		u_grad_trial[nDOF_trial_element*nSpace], 
		u_grad_test_dV[nDOF_test_element*nSpace],
		//for general use
		jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace],
		dV,x,y,z,xt,yt,zt;
	      //get the physical integration weight
	      ck.calculateMapping_element(eN,
					  k,
					  mesh_dof,
					  mesh_l2g,
					  mesh_trial_ref,
					  mesh_grad_trial_ref,
					  jac,
					  jacDet,
					  jacInv,
					  x,y,z);
	      ck.calculateMappingVelocity_element(eN,
						  k,
						  mesh_velocity_dof,
						  mesh_l2g,
						  mesh_trial_ref,
						  xt,yt,zt);	      
	      dV = fabs(jacDet)*dV_ref[k];
	      ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],
				  jacInv,u_grad_trial);
	      // get the components of the normal reconstruction
	      ck.valFromDOF(lumped_wx,
			    &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			    wx);
	      ck.valFromDOF(lumped_wy,
			    &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			    wy);
	      // get the solution (of Newton's solver)
	      ck.valFromDOF(u_dof,
			    &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			    u);
	      // get old solutions
	      ck.valFromDOF(u_dof_old,
			    &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			    un);
	      ck.valFromDOF(uStar_dof,
			    &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			    uStar);
	      //get the solution gradients at quad points	      
	      ck.gradFromDOF(u_dof,
			     &u_l2g[eN_nDOF_trial_element],u_grad_trial,
			     grad_u);
	      ck.gradFromDOF(u_dof_old,
			     &u_l2g[eN_nDOF_trial_element],u_grad_trial,
			     grad_un);
	      ck.gradFromDOF(uStar_dof,
			     &u_l2g[eN_nDOF_trial_element],u_grad_trial,
			     grad_uStar);
	      //precalculate test function products with integration weights for mass matrix terms
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
		  for (int I=0;I<nSpace;I++)
		    u_grad_test_dV[j*nSpace+I] = u_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
		}
	      //calculate time derivative at quadrature points
	      if (q_dV_last[eN_k] <= -100)
		q_dV_last[eN_k] = dV;
	      q_dV[eN_k] = dV;
	      //
	      //moving mesh
	      //
	      double mesh_velocity[3];
	      mesh_velocity[0] = xt;
	      mesh_velocity[1] = yt;
	      mesh_velocity[2] = zt;

	      double gradu2 = 0;
	      for(int I=0;I<nSpace;I++)
		gradu2 += grad_u[I]*grad_u[I];
	      double beta_norm_grad_u = betaNormGrad(gradu2,1.E-10);

	      double gradun2 = 0;
	      for(int I=0;I<nSpace;I++)
		gradun2 += grad_un[I]*grad_un[I];
	      double beta_norm_grad_un = betaNormGrad(gradun2,1.E-10);

	      double graduStar2 = 0;
	      for(int I=0;I<nSpace;I++)
		graduStar2 += grad_uStar[I]*grad_uStar[I];
	      double beta_norm_grad_uStar = betaNormGrad(graduStar2,1.E-10);
	      
	      double lambda = epsFactDiffusion;	      
	      //double lambda = epsFactDiffusion*elementDiameter[eN]/dt;

	      //double coeffFullNewton = -1./beta_norm_grad_u; // single potential
	      double coeffFullNewton = 2*std::pow(beta_norm_grad_u,2)-3*beta_norm_grad_u; //2 pot.
	      
	      //double coeffLinNewton = -1./beta_norm_grad_un; // single potential
	      double coeffLinNewton = 2*std::pow(beta_norm_grad_un,2)-3*beta_norm_grad_un; //2 pot.
 	      //-3+3./2*beta_norm_grad_un+0.5/beta_norm_grad_un;

	      double epsHeaviside = epsFactHeaviside*elementDiameter[eN]; 	    
	      double Hn = smoothedHeaviside(epsHeaviside,un);
	      double Hnp1 = smoothedHeaviside(epsHeaviside,u);
	      for (int I=0;I<nSpace;I++)
	      {
	        relative_velocity[I] = (velocity[eN_k_nSpace+I]-MOVING_DOMAIN*mesh_velocity[I]);
		//f[I] = relative_velocity[I]*Hn;
		f[I] = relative_velocity[I]*Hnp1; //implicit advection via BDF
		//f[I] = 0.5*relative_velocity[I]*(Hnp1+Hn); //implicit advection via CN
	      }
	      //////////////////////////////
	      // CALCULATE CELL BASED CFL //
	      //////////////////////////////
	      calculateCFL(elementDiameter[eN],relative_velocity,cfl[eN_k]); 
	      double time_derivative_residual = (smoothedHeaviside(epsHeaviside,u)-Hn)/dt;
	      //double time_derivative_residual = (u-un)/dt;
	      //////////////
	      // ith-LOOP //
	      //////////////
	      normalReconstruction[0] = wx;
	      normalReconstruction[1] = wy;

	      // PROJECTOR TO NORMAL FIELD //
	      double norm2_w = wx*wx + wy*wy;
	      double coeff = 1.0-std::sqrt(norm2_w);
	      norm2_w += 1E-10;
	      
	      double projector[nSpace][nSpace]; //nn^T
	      for (int I = 0; I < nSpace; ++I)
		for (int J = 0; J < nSpace; ++J)
		  projector[I][J] = normalReconstruction[I]*normalReconstruction[J]/norm2_w;
	      Mult(projector,grad_u,proj_times_grad_u);
	      double alpha=0.5;
	      double hK=elementDiameter[eN];
	      for(int i=0;i<nDOF_test_element;i++) 
		{ 
		  register int i_nSpace=i*nSpace;
		  elementResidual_u[i] += 
		    time_derivative_residual*u_test_dV[i]
		    + ck.Advection_weak(f,&u_grad_test_dV[i_nSpace])
		    // option 1
		    +hK*lambda*ck.NumericalDiffusion(std::sqrt(norm2_w),
							    proj_times_grad_u,
							    &u_grad_test_dV[i_nSpace])
		    +hK*lambda*ck.NumericalDiffusion(coeff,grad_u,&u_grad_test_dV[i_nSpace])
		    // end of option 1
		    // option 2
		    //+ lambda*ck.NumericalDiffusion(1.0-alpha,
		    //				   proj_times_grad_u,&u_grad_test_dV[i_nSpace])
		    //+ lambda*ck.NumericalDiffusion(alpha,grad_u,&u_grad_test_dV[i_nSpace])
		    //end of option 2
		    + hK*lambda*ck.NumericalDiffusion(-1.0,
						   normalReconstruction,
						   &u_grad_test_dV[i_nSpace]);
		}//i
	      //save solution for other models 
	      q_u[eN_k] = u;
	      q_m[eN_k] = u;//porosity*u;
	    }
	  /////////////////
	  // DISTRIBUTE // load cell based element into global residual
	  ////////////////
	  for(int i=0;i<nDOF_test_element;i++) 
	    { 
	      int eN_i=eN*nDOF_test_element+i;
	      int gi = offset_u+stride_u*u_l2g[eN_i]; //global i-th index
	      // distribute global residual for (lumped) mass matrix
	      globalResidual[gi] += elementResidual_u[i];
	    }//i
	}//elements
      
      //////////////
      // BOUNDARY //
      //////////////
      if(false)
      for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) 
	{ 
	  register int ebN = exteriorElementBoundariesArray[ebNE]; 
	  register int eN  = elementBoundaryElementsArray[ebN*2+0],
	    ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+0],
	    eN_nDOF_trial_element = eN*nDOF_trial_element;
	  register double elementResidual_u[nDOF_test_element];
	  for (int i=0;i<nDOF_test_element;i++)
	    elementResidual_u[i]=0.0;
	  // loop on quad points
	  for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) 
	    { 
	      register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
		ebNE_kb_nSpace = ebNE_kb*nSpace,
		ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb;
	      register double 
		phin=0.0, 
		relative_velocity[nSpace],
		jac_ext[nSpace*nSpace],
		jacDet_ext,
		jacInv_ext[nSpace*nSpace],
		boundaryJac[nSpace*(nSpace-1)],
		metricTensor[(nSpace-1)*(nSpace-1)],
		metricTensorDetSqrt,
		dS,
		u_test_dS[nDOF_test_element],
		normal[nSpace],
		x_ext,y_ext,z_ext,xt_ext,yt_ext,zt_ext,integralScaling,porosity_ext;
	      // calculate mappings 
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
						  x_ext,y_ext,z_ext);
	      ck.calculateMappingVelocity_elementBoundary(eN,
							  ebN_local,
							  kb,
							  ebN_local_kb,
							  mesh_velocity_dof,
							  mesh_l2g,
							  mesh_trial_trace_ref,
							  xt_ext,yt_ext,zt_ext,
							  normal,
							  boundaryJac,
							  metricTensor,
							  integralScaling);
	      dS = ((1.0-MOVING_DOMAIN)*metricTensorDetSqrt +
		    MOVING_DOMAIN*integralScaling)*dS_ref[kb];
	      //compute shape and solution information
	      ck.valFromDOF(phin_dof,&u_l2g[eN_nDOF_trial_element],
			    &u_trial_trace_ref[ebN_local_kb*nDOF_test_element],phin);
	      //precalculate test function products with integration weights
	      for (int j=0;j<nDOF_trial_element;j++)
		u_test_dS[j] = u_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;
	      //
	      //moving mesh
	      //
	      double mesh_velocity[3];
	      mesh_velocity[0] = xt_ext;
	      mesh_velocity[1] = yt_ext;
	      mesh_velocity[2] = zt_ext;
	      //std::cout<<"mesh_velocity ext"<<std::endl;
	      for (int I=0;I<nSpace;I++)
		relative_velocity[I] = (ebqe_velocity_ext[ebNE_kb_nSpace+I]
					- MOVING_DOMAIN*mesh_velocity[I]);
	      double flow = 0.;
	      for (int I=0; I < nSpace; I++)
		flow += normal[I]*relative_velocity[I];
	      
	      double epsHeaviside = epsFactHeaviside*elementDiameter[eN]; 
	      double Hn = smoothedHeaviside(epsHeaviside,phin);
	      for (int i=0;i<nDOF_trial_element;i++)
		elementResidual_u[i] += flow*Hn*u_test_dS[i];
	    }//kb	  
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      int eN_i = eN*nDOF_test_element+i;
	      int gi = offset_u+stride_u*u_l2g[eN_i]; //global i-th index
	      globalResidual[gi] += elementResidual_u[i];
	    }
	}//ebNE
      // END OF BOUNDARY //
    }
        
    void calculateJacobian_MCorr_with_VOF(//element
			     double dt,
			     double* mesh_trial_ref,
			     double* mesh_grad_trial_ref,
			     double* mesh_dof,
			     double* mesh_velocity_dof,
			     double MOVING_DOMAIN,
			     int* mesh_l2g,
			     double* dV_ref,
			     double* u_trial_ref,
			     double* u_grad_trial_ref,
			     double* u_test_ref,
			     double* u_grad_test_ref,
			     //element boundary
			     double* mesh_trial_trace_ref,
			     double* mesh_grad_trial_trace_ref,
			     double* dS_ref,
			     double* u_trial_trace_ref,
			     double* u_grad_trial_trace_ref,
			     double* u_test_trace_ref,
			     double* u_grad_test_trace_ref,
			     double* normal_ref,
			     double* boundaryJac_ref,
			     //physics
			     int nElements_global,
			     double useMetrics, 
			     double alphaBDF,
			     int lag_shockCapturing,/*mwf not used yet*/
			     double shockCapturingDiffusion,
			     //VRANS
			     const double* q_porosity,
			     //
			     int* u_l2g,
			     double* elementDiameter,
			     int degree_polynomial,
			     double* u_dof, 
			     double* velocity,
			     double* q_m_betaBDF, 
			     double* cfl,
			     double* q_numDiff_u_last, 
			     int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
			     double* globalJacobian,
			     int nExteriorElementBoundaries_global,
			     int* exteriorElementBoundariesArray,
			     int* elementBoundaryElementsArray,
			     int* elementBoundaryLocalElementBoundariesArray,
			     double* ebqe_velocity_ext,
			     //VRANS
			     const double* ebqe_porosity_ext,
			     //
			     int* isDOFBoundary_u,
			     double* ebqe_bc_u_ext,
			     int* isFluxBoundary_u,
			     double* ebqe_bc_flux_u_ext,
			     int* csrColumnOffsets_eb_u_u,
			     int LUMPED_MASS_MATRIX,
			     // FOR NONLINEAR VOF; i.e., MCorr with VOF
			     int useFullNewton,
			     double epsFactHeaviside,
			     double epsFactDirac,
			     double epsFactDiffusion,
			     double cK,
			     double uL,
			     double uR,
			     double* phin_dof,
			     double* phiHat_dof,
			     double* lumped_wx,
			     double* lumped_wy)
    {
      ////////////////////////
      // loop over elements //
      ////////////////////////
      for(int eN=0;eN<nElements_global;eN++)
	{
	  register double  elementJacobian_u_u[nDOF_test_element][nDOF_trial_element];
	  for (int i=0;i<nDOF_test_element;i++)
	    for (int j=0;j<nDOF_trial_element;j++)
	      elementJacobian_u_u[i][j]=0.0;	      
	  for  (int k=0;k<nQuadraturePoints_element;k++)
	    {
	      int eN_k = eN*nQuadraturePoints_element+k, //index to a scalar at a quadrature point
		eN_k_nSpace = eN_k*nSpace,
		eN_nDOF_trial_element = eN*nDOF_trial_element; //index to a vector at a quadrature point
	      //declare local storage
	      register double u, phiHatnp1, u_grad_trial[nDOF_trial_element*nSpace],
		df[nSpace], relative_velocity[nSpace],
		jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace],		
		dV, u_test_dV[nDOF_test_element], u_grad_test_dV[nDOF_test_element*nSpace],
		x,y,z,xt,yt,zt;
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
					  x,y,z);
	      ck.calculateMappingVelocity_element(eN,
						  k,
						  mesh_velocity_dof,
						  mesh_l2g,
						  mesh_trial_ref,
						  xt,yt,zt);	      
	      //get the physical integration weight
	      dV = fabs(jacDet)*dV_ref[k];
	      //get the trial function gradients
	      ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,u_grad_trial);
	      //get the solution 	
	      ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],u);
	      //get phiHat at tnp1
	      ck.valFromDOF(phiHat_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],phiHatnp1);
	      //precalculate test function products with integration weights
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
		  for (int I=0;I<nSpace;I++)
		    u_grad_test_dV[j*nSpace+I]   = u_grad_trial[j*nSpace+I]*dV;
		}
	      double epsDiffusion = epsFactDiffusion*elementDiameter[eN];
	      double epsDirac = epsFactDirac*elementDiameter[eN];
	      //double time_derivative_jacobian = smoothedDirac(epsDirac,phiHatnp1+u)/dt;
	      double time_derivative_jacobian = smoothedDirac(epsDirac,phiHatnp1+u); //TMP
	      //double time_derivative_jacobian=epsFactDiffusion*smoothedDirac(epsDirac,
	      //							     phiHatnp1+epsFactDiffusion*u)/dt; //TMP
	      
	      ////////////////////////////////
	      // FOR IMPLICIT NONLINEAR VOF // (implicit advection term)
	      ////////////////////////////////
	      //double mesh_velocity[3];
	      //mesh_velocity[0] = xt;
	      //mesh_velocity[1] = yt;
	      //mesh_velocity[2] = zt;
	      //double dHnp1 = smoothedDirac(epsDirac,phiHatnp1+u);
	      //for (int I=0;I<nSpace;I++)
	      //{
	      //relative_velocity[I] = (velocity[eN_k_nSpace+I]-MOVING_DOMAIN*mesh_velocity[I]);
	      //df[I] = relative_velocity[I]*dHnp1;
	      //}
	      ////////////////////////////////
  	      for(int i=0;i<nDOF_test_element;i++)
		{
		  for(int j=0;j<nDOF_trial_element;j++)
		    {
		      int j_nSpace = j*nSpace;
		      int i_nSpace = i*nSpace;		
		      elementJacobian_u_u[i][j] +=
			time_derivative_jacobian*u_trial_ref[k*nDOF_trial_element+j]*u_test_dV[i]
			//+ ck.AdvectionJacobian_weak(df,u_trial_ref[k*nDOF_trial_element+j],
			//			    &u_grad_test_dV[i_nSpace])
			//+ ck.NumericalDiffusionJacobian(epsDiffusion/dt,
			//				&u_grad_trial[j_nSpace],
			//				&u_grad_test_dV[i_nSpace]);
			+ ck.NumericalDiffusionJacobian(epsDiffusion,
							&u_grad_trial[j_nSpace],
							&u_grad_test_dV[i_nSpace]); //TMP
			//+ ck.NumericalDiffusionJacobian(1.0,
			//				&u_grad_trial[j_nSpace],
			//				&u_grad_test_dV[i_nSpace]); //TMP
		    }//j
		}//i
	    }//k
	  //
	  //load into element Jacobian into global Jacobian
	  //
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      int eN_i = eN*nDOF_test_element+i;
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  int eN_i_j = eN_i*nDOF_trial_element+j;
		  globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_u_u[eN_i_j]] +=
		    elementJacobian_u_u[i][j];
		}//j
	    }//i
	}//elements
    }//computeJacobian for MCorr with VOF

    void calculateJacobian_MCorr_with_VOF2(//element
			     double dt,
			     double* mesh_trial_ref,
			     double* mesh_grad_trial_ref,
			     double* mesh_dof,
			     double* mesh_velocity_dof,
			     double MOVING_DOMAIN,
			     int* mesh_l2g,
			     double* dV_ref,
			     double* u_trial_ref,
			     double* u_grad_trial_ref,
			     double* u_test_ref,
			     double* u_grad_test_ref,
			     //element boundary
			     double* mesh_trial_trace_ref,
			     double* mesh_grad_trial_trace_ref,
			     double* dS_ref,
			     double* u_trial_trace_ref,
			     double* u_grad_trial_trace_ref,
			     double* u_test_trace_ref,
			     double* u_grad_test_trace_ref,
			     double* normal_ref,
			     double* boundaryJac_ref,
			     //physics
			     int nElements_global,
			     double useMetrics, 
			     double alphaBDF,
			     int lag_shockCapturing,/*mwf not used yet*/
			     double shockCapturingDiffusion,
			     //VRANS
			     const double* q_porosity,
			     //
			     int* u_l2g,
			     double* elementDiameter,
			     int degree_polynomial,
			     double* u_dof, 
			     double* velocity,
			     double* q_m_betaBDF, 
			     double* cfl,
			     double* q_numDiff_u_last, 
			     int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
			     double* globalJacobian,
			     int nExteriorElementBoundaries_global,
			     int* exteriorElementBoundariesArray,
			     int* elementBoundaryElementsArray,
			     int* elementBoundaryLocalElementBoundariesArray,
			     double* ebqe_velocity_ext,
			     //VRANS
			     const double* ebqe_porosity_ext,
			     //
			     int* isDOFBoundary_u,
			     double* ebqe_bc_u_ext,
			     int* isFluxBoundary_u,
			     double* ebqe_bc_flux_u_ext,
			     int* csrColumnOffsets_eb_u_u,
			     int LUMPED_MASS_MATRIX,
			     // FOR NONLINEAR VOF; i.e., MCorr with VOF
			     int useFullNewton,
			     double epsFactHeaviside,
			     double epsFactDirac,
			     double epsFactDiffusion,
			     double cK,
			     double uL,
			     double uR,
			     double* phin_dof,
			     double* phiHat_dof,
			     double* lumped_wx,
			     double* lumped_wy)
    {
      ////////////////////////
      // loop over elements //
      ////////////////////////
      for(int eN=0;eN<nElements_global;eN++)
	{
	  register double  elementJacobian_u_u[nDOF_test_element][nDOF_trial_element];
	  for (int i=0;i<nDOF_test_element;i++)
	    for (int j=0;j<nDOF_trial_element;j++)
	      elementJacobian_u_u[i][j]=0.0;	      
	  for  (int k=0;k<nQuadraturePoints_element;k++)
	    {
	      int eN_k = eN*nQuadraturePoints_element+k, //index to a scalar at a quadrature point
		eN_k_nSpace = eN_k*nSpace,
		eN_nDOF_trial_element = eN*nDOF_trial_element; //index to a vector at a quadrature point
	      //declare local storage
	      register double u, grad_u[nSpace], phiHatnp1, u_grad_trial[nDOF_trial_element*nSpace],
		df[nSpace], relative_velocity[nSpace],
		jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace],		
		dV, u_test_dV[nDOF_test_element], u_grad_test_dV[nDOF_test_element*nSpace],
		x,y,z,xt,yt,zt;
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
					  x,y,z);
	      ck.calculateMappingVelocity_element(eN,
						  k,
						  mesh_velocity_dof,
						  mesh_l2g,
						  mesh_trial_ref,
						  xt,yt,zt);	      
	      //get the physical integration weight
	      dV = fabs(jacDet)*dV_ref[k];
	      //get the trial function gradients
	      ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],
				  jacInv,u_grad_trial);
	      //get the solution 	
	      ck.valFromDOF(u_dof,
			    &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			    u);
	      ck.gradFromDOF(u_dof,
			     &u_l2g[eN_nDOF_trial_element],u_grad_trial,
			     grad_u);
	      //precalculate test function products with integration weights
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
		  for (int I=0;I<nSpace;I++)
		    u_grad_test_dV[j*nSpace+I]   = u_grad_trial[j*nSpace+I]*dV;
		}

	      double gradu2 = 0;
	      for(int I=0;I<nSpace;I++)
		gradu2 += grad_u[I]*grad_u[I];
	      double beta_norm_grad_u = betaNormGrad(gradu2,1E-15);
	      
	      double lambda = epsFactDiffusion;
	      //double lambda = epsFactDiffusion*elementDiameter[eN]/dt;

	      // single potential
	      //double coeff1FullNonlinear = 1.-1./beta_norm_grad_u; 
	      //double coeff2FullNonlinear = 1./std::pow(beta_norm_grad_u,3);
	      
	      // double potential
	      double coeff1FullNonlinear=fmax(1E-10,
					      1+2*std::pow(beta_norm_grad_u,2)-3*beta_norm_grad_u);
	      double coeff2FullNonlinear=fmax(1E-10,
					      4.-3./beta_norm_grad_u);

	      double epsDirac = epsFactDirac*elementDiameter[eN];
	      double time_derivative_jacobian = smoothedDirac(epsDirac,u)/dt;
	      //double time_derivative_jacobian = 1/dt;

	      ////////////////////////////////
	      // FOR IMPLICIT NONLINEAR VOF // (implicit advection term)
	      ////////////////////////////////
	      double mesh_velocity[3];
	      mesh_velocity[0] = xt;
	      mesh_velocity[1] = yt;
	      mesh_velocity[2] = zt;
	      double dHnp1 = smoothedDirac(epsDirac,u);
	      for (int I=0;I<nSpace;I++)
		{
		  relative_velocity[I] = (velocity[eN_k_nSpace+I]-MOVING_DOMAIN*mesh_velocity[I]);
		  df[I] = relative_velocity[I]*dHnp1;
		}
	      /////////////////////////////////
	      if (useFullNewton==false)
		for(int i=0;i<nDOF_test_element;i++)
		  {
		    for(int j=0;j<nDOF_trial_element;j++)
		      {
			int j_nSpace = j*nSpace;
			int i_nSpace = i*nSpace;		
			elementJacobian_u_u[i][j] +=
			  time_derivative_jacobian*u_trial_ref[k*nDOF_trial_element+j]*u_test_dV[i]
			  + ck.AdvectionJacobian_weak(df,u_trial_ref[k*nDOF_trial_element+j],
						      &u_grad_test_dV[i_nSpace])
			  + lambda*ck.NumericalDiffusionJacobian(1.0,
								 &u_grad_trial[j_nSpace],
								 &u_grad_test_dV[i_nSpace]);
		      }//j
		  }//i
	      else //use full newton method
		for(int i=0;i<nDOF_test_element;i++)
		  {
		    for(int j=0;j<nDOF_trial_element;j++)
		      {
			int j_nSpace = j*nSpace;
			int i_nSpace = i*nSpace;		
			elementJacobian_u_u[i][j] +=
			  time_derivative_jacobian*u_trial_ref[k*nDOF_trial_element+j]*u_test_dV[i]
			  + lambda*ck.NumericalDiffusionJacobian(coeff1FullNonlinear,
								 &u_grad_trial[j_nSpace],
								 &u_grad_test_dV[i_nSpace])
			  + lambda*( coeff2FullNonlinear*dV*
				     ck.NumericalDiffusion(1.0,grad_u,&u_grad_trial[i_nSpace])*
				     ck.NumericalDiffusion(1.0,grad_u,&u_grad_trial[j_nSpace]) );
		      }//j
		  }//i
	    }//k
	  //
	  //load into element Jacobian into global Jacobian
	  //
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      int eN_i = eN*nDOF_test_element+i;
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  int eN_i_j = eN_i*nDOF_trial_element+j;
		  globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_u_u[eN_i_j]] +=
		    elementJacobian_u_u[i][j];
		}//j
	    }//i
	}//elements
    }//computeJacobian for MCorr with VOF

    void calculateJacobian_MCorr_with_VOF4(//element
			     double dt,
			     double* mesh_trial_ref,
			     double* mesh_grad_trial_ref,
			     double* mesh_dof,
			     double* mesh_velocity_dof,
			     double MOVING_DOMAIN,
			     int* mesh_l2g,
			     double* dV_ref,
			     double* u_trial_ref,
			     double* u_grad_trial_ref,
			     double* u_test_ref,
			     double* u_grad_test_ref,
			     //element boundary
			     double* mesh_trial_trace_ref,
			     double* mesh_grad_trial_trace_ref,
			     double* dS_ref,
			     double* u_trial_trace_ref,
			     double* u_grad_trial_trace_ref,
			     double* u_test_trace_ref,
			     double* u_grad_test_trace_ref,
			     double* normal_ref,
			     double* boundaryJac_ref,
			     //physics
			     int nElements_global,
			     double useMetrics, 
			     double alphaBDF,
			     int lag_shockCapturing,/*mwf not used yet*/
			     double shockCapturingDiffusion,
			     //VRANS
			     const double* q_porosity,
			     //
			     int* u_l2g,
			     double* elementDiameter,
			     int degree_polynomial,
			     double* u_dof, 
			     double* velocity,
			     double* q_m_betaBDF, 
			     double* cfl,
			     double* q_numDiff_u_last, 
			     int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
			     double* globalJacobian,
			     int nExteriorElementBoundaries_global,
			     int* exteriorElementBoundariesArray,
			     int* elementBoundaryElementsArray,
			     int* elementBoundaryLocalElementBoundariesArray,
			     double* ebqe_velocity_ext,
			     //VRANS
			     const double* ebqe_porosity_ext,
			     //
			     int* isDOFBoundary_u,
			     double* ebqe_bc_u_ext,
			     int* isFluxBoundary_u,
			     double* ebqe_bc_flux_u_ext,
			     int* csrColumnOffsets_eb_u_u,
			     int LUMPED_MASS_MATRIX,
			     // FOR NONLINEAR VOF; i.e., MCorr with VOF
			     int useFullNewton,
			     double epsFactHeaviside,
			     double epsFactDirac,
			     double epsFactDiffusion,
			     double cK,
			     double uL,
			     double uR,
			     double* phin_dof,
			     double* phiHat_dof,
			     double* lumped_wx,
			     double* lumped_wy)
    {
      ////////////////////////
      // loop over elements //
      ////////////////////////
      for(int eN=0;eN<nElements_global;eN++)
	{
	  register double  elementJacobian_u_u[nDOF_test_element][nDOF_trial_element];
	  for (int i=0;i<nDOF_test_element;i++)
	    for (int j=0;j<nDOF_trial_element;j++)
	      elementJacobian_u_u[i][j]=0.0;	      
	  for  (int k=0;k<nQuadraturePoints_element;k++)
	    {
	      int eN_k = eN*nQuadraturePoints_element+k, //index to a scalar at a quadrature point
		eN_k_nSpace = eN_k*nSpace,
		eN_nDOF_trial_element = eN*nDOF_trial_element; //index to a vector at a quadrature point
	      //declare local storage
	      register double u, grad_u[nSpace], phiHatnp1, u_grad_trial[nDOF_trial_element*nSpace],
		normalReconstruction[nSpace], wx, wy, proj_times_grad_trial[nSpace],
		df[nSpace], relative_velocity[nSpace],
		jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace],		
		dV, u_test_dV[nDOF_test_element], u_grad_test_dV[nDOF_test_element*nSpace],
		x,y,z,xt,yt,zt;
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
					  x,y,z);
	      ck.calculateMappingVelocity_element(eN,
						  k,
						  mesh_velocity_dof,
						  mesh_l2g,
						  mesh_trial_ref,
						  xt,yt,zt);	      
	      //get the physical integration weight
	      dV = fabs(jacDet)*dV_ref[k];
	      //get the trial function gradients
	      ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],
				  jacInv,u_grad_trial);
	      // get the components of the normal reconstruction
	      ck.valFromDOF(lumped_wx,
			    &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			    wx);
	      ck.valFromDOF(lumped_wy,
			    &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			    wy);
	      //get the solution 	
	      ck.valFromDOF(u_dof,
			    &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			    u);
	      ck.gradFromDOF(u_dof,
			     &u_l2g[eN_nDOF_trial_element],u_grad_trial,
			     grad_u);
	      //precalculate test function products with integration weights
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
		  for (int I=0;I<nSpace;I++)
		    u_grad_test_dV[j*nSpace+I]   = u_grad_trial[j*nSpace+I]*dV;
		}

	      double gradu2 = 0;
	      for(int I=0;I<nSpace;I++)
		gradu2 += grad_u[I]*grad_u[I];
	      double beta_norm_grad_u = betaNormGrad(gradu2,1E-15);
	      
	      double lambda = epsFactDiffusion;
	      //double lambda = epsFactDiffusion*elementDiameter[eN]/dt;

	      // single potential
	      //double coeff1FullNonlinear = 1.-1./beta_norm_grad_u; 
	      //double coeff2FullNonlinear = 1./std::pow(beta_norm_grad_u,3);
	      
	      // double potential
	      double coeff1FullNonlinear=fmax(1E-10,
					      1+2*std::pow(beta_norm_grad_u,2)-3*beta_norm_grad_u);
	      double coeff2FullNonlinear=fmax(1E-10,
					      4.-3./beta_norm_grad_u);

	      double epsDirac = epsFactDirac*elementDiameter[eN];
	      double time_derivative_jacobian = smoothedDirac(epsDirac,u)/dt;
	      //double time_derivative_jacobian = 1/dt;

	      ////////////////////////////////
	      // FOR IMPLICIT NONLINEAR VOF // (implicit advection term)
	      ////////////////////////////////
	      double mesh_velocity[3];
	      mesh_velocity[0] = xt;
	      mesh_velocity[1] = yt;
	      mesh_velocity[2] = zt;
	      double dHnp1 = smoothedDirac(epsDirac,u);
	      for (int I=0;I<nSpace;I++)
		{
		  relative_velocity[I] = (velocity[eN_k_nSpace+I]-MOVING_DOMAIN*mesh_velocity[I]);
		  df[I] = relative_velocity[I]*dHnp1;
		}

	      normalReconstruction[0] = wx;
	      normalReconstruction[1] = wy;
	      // PROJECTOR TO NORMAL FIELD //
	      double norm2_w = wx*wx + wy*wy;
	      double coeff = 1.0-std::sqrt(norm2_w);
	      norm2_w += 1E-10;
	      
	      double projector[nSpace][nSpace]; //nn^T
	      for (int I = 0; I < nSpace; ++I)
		for (int J = 0; J < nSpace; ++J)
		  projector[I][J] = normalReconstruction[I]*normalReconstruction[J]/norm2_w;

	      double alpha=0.5;
	      double hK=elementDiameter[eN];
	      /////////////////////////////////
	      if (useFullNewton==false)
		for(int i=0;i<nDOF_test_element;i++)
		  {
		    for(int j=0;j<nDOF_trial_element;j++)
		      {
			int j_nSpace = j*nSpace;
			int i_nSpace = i*nSpace;
			
			Mult(projector,&u_grad_trial[j_nSpace],proj_times_grad_trial);
				      
			elementJacobian_u_u[i][j] +=
			  time_derivative_jacobian*u_trial_ref[k*nDOF_trial_element+j]*u_test_dV[i]
			  + ck.AdvectionJacobian_weak(df,u_trial_ref[k*nDOF_trial_element+j],
						      &u_grad_test_dV[i_nSpace])
			  // option 1
			  + hK*lambda*ck.NumericalDiffusionJacobian(std::sqrt(norm2_w),
			  					 proj_times_grad_trial,
			  					 &u_grad_test_dV[i_nSpace])
			  + hK*lambda*ck.NumericalDiffusionJacobian(coeff,
			  					 &u_grad_trial[j_nSpace],
			  					 &u_grad_test_dV[i_nSpace]);
			  // option 2
			  //+ lambda*ck.NumericalDiffusionJacobian(1.0-alpha,
			  //					 proj_times_grad_trial,
			  //					 &u_grad_test_dV[i_nSpace])
			  //+ lambda*ck.NumericalDiffusionJacobian(alpha,
			  //					 &u_grad_trial[j_nSpace],
			  //					 &u_grad_test_dV[i_nSpace]);
			
		      }//j
		  }//i
	      else //use full newton method
		for(int i=0;i<nDOF_test_element;i++)
		  {
		    for(int j=0;j<nDOF_trial_element;j++)
		      {
			int j_nSpace = j*nSpace;
			int i_nSpace = i*nSpace;		
			elementJacobian_u_u[i][j] +=
			  time_derivative_jacobian*u_trial_ref[k*nDOF_trial_element+j]*u_test_dV[i]
			  + lambda*ck.NumericalDiffusionJacobian(coeff1FullNonlinear,
								 &u_grad_trial[j_nSpace],
								 &u_grad_test_dV[i_nSpace])
			  + lambda*( coeff2FullNonlinear*dV*
				     ck.NumericalDiffusion(1.0,grad_u,&u_grad_trial[i_nSpace])*
				     ck.NumericalDiffusion(1.0,grad_u,&u_grad_trial[j_nSpace]) );
		      }//j
		  }//i
	    }//k
	  //
	  //load into element Jacobian into global Jacobian
	  //
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      int eN_i = eN*nDOF_test_element+i;
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  int eN_i_j = eN_i*nDOF_trial_element+j;
		  globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_u_u[eN_i_j]] +=
		    elementJacobian_u_u[i][j];
		}//j
	    }//i
	}//elements
    }//computeJacobian for MCorr with VOF

    void calculateRhsL2p(
			 double* mesh_trial_ref, //
			 double* mesh_grad_trial_ref, //
			 double* mesh_dof, //
			 int* mesh_l2g, //
			 double* dV_ref, //
			 double* u_trial_ref, 
			 double* u_grad_trial_ref,
			 double* u_test_ref, //
			 //physics
			 int nElements_global, //
			 int* u_l2g, //
			 double* elementDiameter,
			 //double* nodeDiametersArray,
			 double* u_dof,
			 double* phiHat_dof,
			 double* phiExact_dof,
			 int offset_u, int stride_u, 
			 double* globalResidual,
			 double* global_mass_error,
			 double* global_L2_interface,
			 double* global_H1_interface,
			 double* global_L2_Hinterface,
			 double* global_H1_Hinterface,
			 double* global_L2_u,
			 double* global_H1_u)			 
    {
      double global_mass_exact = 0.0;
      *global_mass_error = 0.0;
      *global_L2_interface = 0.0;
      *global_H1_interface = 0.0;
      *global_L2_Hinterface = 0.0;
      *global_H1_Hinterface = 0.0;
      *global_L2_u = 0.0;
      *global_H1_u = 0.0;
      //////////////////////////////////////////////
      // ** LOOP IN CELLS FOR CELL BASED TERMS ** //
      //////////////////////////////////////////////
      for(int eN=0;eN<nElements_global;eN++)
	{
	  //declare local storage for local contributions and initialize
	  register double 
	    elementResidual_u[nDOF_test_element];
	  double cell_mass_error = 0., cell_mass_exact = 0.,
	    cell_L2_u = 0., cell_L2_int = 0., cell_L2_Hint = 0.,
	    cell_H1_u = 0., cell_H1_int = 0., cell_H1_Hint = 0.;
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      elementResidual_u[i]=0.0;
	    }
	  
	  //loop over quadrature points and compute integrands
	  for  (int k=0;k<nQuadraturePoints_element;k++)
	    {
	      //compute indeces and declare local storage
	      register int eN_k = eN*nQuadraturePoints_element+k,
		eN_k_nSpace = eN_k*nSpace,
		eN_nDOF_trial_element = eN*nDOF_trial_element;
	      register double
		u, phiHat, phiExact,
		u_grad_trial[nDOF_trial_element*nSpace],
		grad_u[nSpace], grad_phiHat[nSpace], grad_phiExact[nSpace],
		grad_int[nSpace], grad_Hint[nSpace],
		u_test_dV[nDOF_trial_element],
		//for general use
		jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace],
		dV,x,y,z;
	      //get the physical integration weight
	      ck.calculateMapping_element(eN,
					  k,
					  mesh_dof,
					  mesh_l2g,
					  mesh_trial_ref,
					  mesh_grad_trial_ref,
					  jac,
					  jacDet,
					  jacInv,
					  x,y,z);
	      dV = fabs(jacDet)*dV_ref[k];
	      // get functions at quad points
	      ck.valFromDOF(u_dof,
			    &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			    u);
	      ck.valFromDOF(phiHat_dof,
			    &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			    phiHat);
	      ck.valFromDOF(phiExact_dof,
			    &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			    phiExact);
	      // get gradients
	      ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],
				  jacInv,
				  u_grad_trial);
	      ck.gradFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_u);
	      ck.gradFromDOF(phiHat_dof,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_phiHat);
	      ck.gradFromDOF(phiExact_dof,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_phiExact);
		      
	      //precalculate test function products with integration weights for mass matrix terms
	      for (int j=0;j<nDOF_trial_element;j++)
		u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;

	      double epsHeaviside = 1.5*elementDiameter[eN];
	      double epsDirac = 1.5*elementDiameter[eN];

	      // gradients of errors in the interfaces
	      for (int I=0;I<nSpace;I++)
		{
		  grad_int[I] = grad_phiHat[I]+grad_u[I] - grad_phiExact[I];
		  grad_Hint[I] = (smoothedDirac(epsDirac,phiHat+u)*(grad_phiHat[I]+grad_u[I])
				  -smoothedDirac(epsDirac,phiExact)*(grad_phiExact[I]));
		}
	      // cell mass error
	      cell_mass_error += smoothedHeaviside(epsHeaviside,phiHat+u)*dV;
	      cell_mass_exact += smoothedHeaviside(epsHeaviside,phiExact)*dV;
	      // L2 component
	      double L2_u_tmp = u*u*dV;
	      double L2_int_tmp = std::pow(phiHat+u - phiExact,2)*dV; 
	      double L2_Hint_tmp = std::pow(smoothedHeaviside(epsHeaviside,phiHat+u) -
					    smoothedHeaviside(epsHeaviside,phiExact),2)*dV;
	      // H1 Semi norm component
	      double H1Semi_u_tmp = ck.NumericalDiffusion(dV,grad_u,grad_u);
	      double H1Semi_int_tmp = ck.NumericalDiffusion(dV,grad_int,grad_int);
	      double H1Semi_Hint_tmp = ck.NumericalDiffusion(dV,grad_Hint,grad_Hint);
	      // cell L2 norms
	      cell_L2_u    += L2_u_tmp;
	      cell_L2_int  += L2_int_tmp;
	      cell_L2_Hint += L2_Hint_tmp;
	      // cell H1 norms
	      cell_H1_u    += L2_u_tmp    + H1Semi_u_tmp;
	      cell_H1_int  += L2_int_tmp  + H1Semi_int_tmp;
	      cell_H1_Hint += L2_Hint_tmp + H1Semi_Hint_tmp;
	      
	      // ith-LOOP //	      
	      for(int i=0;i<nDOF_test_element;i++)
		elementResidual_u[i] += smoothedHeaviside(epsHeaviside,phiHat+u)*u_test_dV[i];
	    }
	  *global_mass_error += cell_mass_error;
	  global_mass_exact += cell_mass_exact;
	  *global_L2_interface += cell_L2_int;
	  *global_H1_interface += cell_H1_int;
	  *global_L2_Hinterface += cell_L2_Hint;
	  *global_H1_Hinterface += cell_H1_Hint;
	  *global_L2_u += cell_L2_u;
	  *global_H1_u += cell_H1_u;
	  /////////////////
	  // DISTRIBUTE // load cell based element into global residual
	  ////////////////
	  for(int i=0;i<nDOF_test_element;i++) 
	    { 
	      int eN_i=eN*nDOF_test_element+i;
	      int gi = offset_u+stride_u*u_l2g[eN_i]; //global i-th index
	      // distribute global residual
	      globalResidual[gi] += elementResidual_u[i];
	    }//i
	}//elements
      *global_mass_error -= global_mass_exact;
    }
    
    void normalReconstruction(//element
			      double* mesh_trial_ref,//
			      double* mesh_grad_trial_ref,
			      double* mesh_dof, //
			      int* mesh_l2g,//
			      double* dV_ref,//
			      double* u_trial_ref,
			      double* u_grad_trial_ref,
			      double* u_test_ref,
			      //physics
			      int nElements_global,//
			      int* u_l2g, //
			      double* elementDiameter,//
			      double* phi_dof,//
			      int offset_u, int stride_u, 
			      // PARAMETERS FOR EDGE VISCOSITY 
			      int numDOFs,
			      double* lumped_wx,
			      double* lumped_wy)
      //double* rhs_mass_correction,
      //double* lumped_L2p,
      //double* lumped_mass_matrix,
      // FOR NONLINEAR VOF; i.e., MCorr with VOF
      //double epsFactHeaviside,
      //double epsFactDiffusion,
    //double* phiHat_dof)
    {
      register double
	weighted_lumped_mass_matrix[numDOFs],
	rhsx_normal_reconstruction[numDOFs],
	rhsy_normal_reconstruction[numDOFs];
      for (int i=0; i<numDOFs; i++)
	{
	  lumped_wx[i]=0.;
	  lumped_wy[i]=0.;
	  weighted_lumped_mass_matrix[i]=0.;
	  rhsx_normal_reconstruction[i]=0.;
	  rhsy_normal_reconstruction[i]=0.;
	}
      for(int eN=0;eN<nElements_global;eN++)
	{
	  //declare local storage for local contributions and initialize
	  register double
	    element_weighted_lumped_mass_matrix[nDOF_test_element],
	    element_rhsx_normal_reconstruction[nDOF_test_element],
	    element_rhsy_normal_reconstruction[nDOF_test_element];	    
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      element_weighted_lumped_mass_matrix[i]=0.0;
	      element_rhsx_normal_reconstruction[i]=0.0;
	      element_rhsy_normal_reconstruction[i]=0.0;
	    }
	  //loop over quadrature points and compute integrands
	  for  (int k=0;k<nQuadraturePoints_element;k++)
	    {
	      //compute indeces and declare local storage
	      register int eN_k = eN*nQuadraturePoints_element+k,
		eN_k_nSpace = eN_k*nSpace,
		eN_nDOF_trial_element = eN*nDOF_trial_element;
	      register double 
		//for mass matrix contributions
		grad_phi[nSpace],
		u_grad_trial[nDOF_trial_element*nSpace],
		u_test_dV[nDOF_trial_element], 
		//for general use
		jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace],
		dV,x,y,z;
	      //get the physical integration weight
	      ck.calculateMapping_element(eN,
					  k,
					  mesh_dof,
					  mesh_l2g,
					  mesh_trial_ref,
					  mesh_grad_trial_ref,
					  jac,
					  jacDet,
					  jacInv,
					  x,y,z);
	      dV = fabs(jacDet)*dV_ref[k];
	      ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],
				  jacInv,
				  u_grad_trial);
	      ck.gradFromDOF(phi_dof,
			     &u_l2g[eN_nDOF_trial_element],u_grad_trial,
			     grad_phi);	      
	      //precalculate test function products with integration weights for mass matrix terms
	      for (int j=0;j<nDOF_trial_element;j++)
		u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
	      
	      double rhsx = grad_phi[0];
	      double rhsy = grad_phi[1];
	      double grad_phi2 = 0;
	      for (int I=0;I<nSpace; I++)
		grad_phi2 += grad_phi[I]*grad_phi[I];
	      double beta_norm_grad_phi = betaNormGrad(grad_phi2,1E-10);
	      
	      for(int i=0;i<nDOF_test_element;i++)
		{
		  element_weighted_lumped_mass_matrix[i] += beta_norm_grad_phi*u_test_dV[i];
		  element_rhsx_normal_reconstruction[i] += rhsx*u_test_dV[i];
		  element_rhsy_normal_reconstruction[i] += rhsy*u_test_dV[i];
		}
	    } //k
	  // DISTRIBUTE //
	  for(int i=0;i<nDOF_test_element;i++) 
	    { 
	      int eN_i=eN*nDOF_test_element+i;
	      int gi = offset_u+stride_u*u_l2g[eN_i]; //global i-th index	      
	      weighted_lumped_mass_matrix[gi] += element_weighted_lumped_mass_matrix[i];
	      rhsx_normal_reconstruction[gi] += element_rhsx_normal_reconstruction[i];
	      rhsy_normal_reconstruction[gi] += element_rhsy_normal_reconstruction[i];
	    }//i
	}//elements
      // COMPUTE LUMPED L2 PROJECTION
      for (int i=0; i<numDOFs; i++)
	{	 
	  double mi = weighted_lumped_mass_matrix[i];
	  lumped_wx[i] = 1./mi*rhsx_normal_reconstruction[i];
	  lumped_wy[i] = 1./mi*rhsy_normal_reconstruction[i];
	}
    }
    
  };//VOF

  inline CLSVOF_base* newCLSVOF(int nSpaceIn,
				int nQuadraturePoints_elementIn,
				int nDOF_mesh_trial_elementIn,
				int nDOF_trial_elementIn,
				int nDOF_test_elementIn,
				int nQuadraturePoints_elementBoundaryIn,
				int CompKernelFlag)
  {
    if (nSpaceIn == 2)
      return proteus::chooseAndAllocateDiscretization2D<CLSVOF_base,CLSVOF,CompKernel>(nSpaceIn,
										 nQuadraturePoints_elementIn,
										 nDOF_mesh_trial_elementIn,
										 nDOF_trial_elementIn,
										 nDOF_test_elementIn,
										 nQuadraturePoints_elementBoundaryIn,
										 CompKernelFlag);
    else
      return proteus::chooseAndAllocateDiscretization<CLSVOF_base,CLSVOF,CompKernel>(nSpaceIn,
									       nQuadraturePoints_elementIn,
									       nDOF_mesh_trial_elementIn,
									       nDOF_trial_elementIn,
									       nDOF_test_elementIn,
									       nQuadraturePoints_elementBoundaryIn,
									       CompKernelFlag);
  }
}//proteus
#endif

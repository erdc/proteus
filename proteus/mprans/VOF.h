#ifndef VOF_H
#define VOF_H
#include <cmath>
#include <iostream>
#include "CompKernel.h"
#include "ModelFactory.h"

#define POWER_SMOOTHNESS_INDICATOR 2
#define IS_BETAij_ONE 0
#define GLOBAL_FCT 0

/////////////////////
//ENTROPY FUNCTION //
/////////////////////
// Power entropy //
#define entropy_power 2. // phiL and phiR are dummy variables
#define ENTROPY(phi,phiL,phiR) 1./entropy_power*std::pow(fabs(phi),entropy_power)
#define DENTROPY(phi,phiL,phiR) std::pow(fabs(phi),entropy_power-1.)*(phi>=0 ? 1 : -1)
// Log entropy //
// LOG ENTROPY FOR LEVEL SET FROM 0 to 1
#define ENTROPY_LOG(phi,phiL,phiR) std::log(fabs((phi-phiL)*(phiR-phi))+1E-14)
#define DENTROPY_LOG(phi,phiL,phiR) (phiL+phiR-2*phi)*((phi-phiL)*(phiR-phi)>=0 ? 1 : -1)/(fabs((phi-phiL)*(phiR-phi))+1E-14) 

namespace proteus
{
  class VOF_base
  {
    //The base class defining the interface
  public:
    virtual ~VOF_base(){}
    virtual void FCTStepL2p(int NNZ, //number on non-zero entries on sparsity pattern
			    int numDOFs, //number of DOFs
			    double* lumped_mass_matrix, //lumped mass matrix (as vector)
			    double* solH, //DOFs of high order solution at tnp1
			    double* solL,
			    double* limited_solution,
			    int* csrRowIndeces_DofLoops, //csr row indeces 
			    int* csrColumnOffsets_DofLoops, //csr column offsets 
			    double* MassMatrix //mass matrix
			    )=0;
    virtual void FCTStep(int NNZ, //number on non-zero entries on sparsity pattern
			 int numDOFs, //number of DOFs
			 double* lumped_mass_matrix, //lumped mass matrix (as vector)
			 double* soln, //DOFs of solution at time tn
			 double* solH, //DOFs of high order solution at tnp1
			 double* low_order_solution,
			 double* limited_solution,
			 int* csrRowIndeces_DofLoops, //csr row indeces 
			 int* csrColumnOffsets_DofLoops, //csr column offsets 
			 double* MassMatrix, //mass matrix
			 double* dt_times_dH_minus_dL, //low minus high order dissipative matrices
			 double* min_u_bc, //min/max value at BCs. If DOF is not at boundary then min=1E10, max=-1E10
			 double* max_u_bc,
			 int LUMPED_MASS_MATRIX
			 )=0;
    virtual void calculateResidual(//element
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
				   double epsFactHeaviside,
				   double epsFactDirac,
				   double epsFactDiffusion,
				   double* phin_dof,
				   double* phiHat_dof,
				   // AUX QUANTITIES OF INTEREST
				   double* quantDOFs)=0;
    virtual void calculateResidual_entropy_viscosity(//element
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
						     double epsFactHeaviside,
						     double epsFactDirac,
						     double epsFactDiffusion,
						     double* phin_dof,
						     double* phiHat_dof,
						     // AUX QUANTITIES OF INTEREST
						     double* quantDOFs)=0;
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
						  double epsFactHeaviside,
						  double epsFactDirac,
						  double epsFactDiffusion,
						  double* phin_dof,
						  double* phiHat_dof,
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
						  double epsFactHeaviside,
						  double epsFactDirac,
						  double epsFactDiffusion,
						  double* phin_dof,
						  double* phiHat_dof,
						  // AUX QUANTITIES OF INTEREST
						  double* quantDOFs)=0;
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
    virtual double calculateRhsQuadratureMass(//element
					    //double dt,
					    double* mesh_trial_ref,
					    double* mesh_grad_trial_ref,
					    double* mesh_dof,
					    //double* mesh_velocity_dof,
					    //double MOVING_DOMAIN,
					    int* mesh_l2g,
					    double* dV_ref,
					    double* u_trial_ref,
					    //double* u_grad_trial_ref,
					    double* u_test_ref,
					    //double* u_grad_test_ref,
					    //element boundary
					    //double* mesh_trial_trace_ref,
					    //double* mesh_grad_trial_trace_ref,
					    //double* dS_ref,
					    //double* u_trial_trace_ref,
					    //double* u_grad_trial_trace_ref,
					    //double* u_test_trace_ref,
					    //double* u_grad_test_trace_ref,
					    //double* normal_ref,
					    //double* boundaryJac_ref,
					    //physics
					    int nElements_global,
					    //double useMetrics, 
					    //double alphaBDF,
					    //int lag_shockCapturing,
					    //double shockCapturingDiffusion,
					    //double sc_uref, 
					    //double sc_alpha,
					    //VRANS
					    //const double* q_porosity,
					    //const double* porosity_dof,
					    //
					    int* u_l2g, 
					    double* elementDiameter,
					    //int degree_polynomial,
					    double* u_dof,
					    //double* u_dof_old,
					    //double* velocity,
					    //double* q_m,
					    //double* q_u,
					    //double* q_m_betaBDF,
					    //double* q_dV,
					    //double* q_dV_last,
					    //double* cfl,
					    //double* edge_based_cfl,
					    //double* q_numDiff_u, 
					    //double* q_numDiff_u_last, 
					    int offset_u, int stride_u, 
					    //double* globalResidual,
					    //int nExteriorElementBoundaries_global,
					    //int* exteriorElementBoundariesArray,
					    //int* elementBoundaryElementsArray,
					    //int* elementBoundaryLocalElementBoundariesArray,
					    //double* ebqe_velocity_ext,
					    //VRANS
					    //const double* ebqe_porosity_ext,
					    //
					    //int* isDOFBoundary_u,
					    //double* ebqe_bc_u_ext,
					    //int* isFluxBoundary_u,
					    //double* ebqe_bc_flux_u_ext,
					    //double* ebqe_phi,double epsFact,
					    //double* ebqe_u,
					    //double* ebqe_flux,
					    // PARAMETERS FOR EDGE BASED STABILIZATION
					    //double cE,
					    //double cK,
					    // PARAMETERS FOR LOG BASED ENTROPY FUNCTION 
					    //double uL, 
					    //double uR, 
					    // PARAMETERS FOR EDGE VISCOSITY 
					    int numDOFs,
					    //int NNZ,
					    //int* csrRowIndeces_DofLoops,
					    //int* csrColumnOffsets_DofLoops,
					    //int* csrRowIndeces_CellLoops,
					    //int* csrColumnOffsets_CellLoops,
					    //int* csrColumnOffsets_eb_CellLoops,
					    // C matrices
					    //double* Cx, 
					    //double* Cy,
					    //double* Cz,
					    //double* CTx,
					    //double* CTy,
					    //double* CTz,
					    //double* ML,
					    // PARAMETERS FOR 1st or 2nd ORDER MPP METHOD
					    //int LUMPED_MASS_MATRIX, 
					    //int STABILIZATION_TYPE,
					    //int ENTROPY_TYPE,
					    // FOR FCT
					    //double* low_order_solution,
					    //double* dt_times_dH_minus_dL,
					    //double* min_u_bc,
					    //double* max_u_bc,
					    // FOR NONLINEAR VOF; i.e., MCorr with VOF
					    //double epsFactHeaviside,
					    //double epsFactDirac,
					    //double epsFactDiffusion,
					    //double* phin_dof,
					    //double* phiHat_dof,
					    // AUX QUANTITIES OF INTEREST
					    //double* quantDOFs)=0;
					    // FOR FCT
					    double* rhs_mass_correction,
					    double* lumped_L2p,
					    double* lumped_mass_matrix,
					    // FOR NONLINEAR VOF; i.e., MCorr with VOF
					    double epsFactHeaviside,
					    double epsFactDiffusion,
					    double* phiHat_dof)=0;    
    virtual void calculateJacobian(//element
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
				   double epsFactHeaviside,
				   double epsFactDirac,
				   double epsFactDiffusion,
				   double cK,
				   double uL,
				   double uR,
				   double* phin_dof,
				   double* phiHat_dof)=0;
    virtual void calculateMassMatrix(//element
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
				     double epsFactHeaviside,
				     double epsFactDirac,
				     double epsFactDiffusion,
				     double cK,
				     double uL,
				     double uR,
				     double* phin_dof,
				     double* phiHat_dof)=0;
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
				     double epsFactHeaviside,
				     double epsFactDirac,
				     double epsFactDiffusion,
				     double cK,
				     double uL,
				     double uR,
				     double* phin_dof,
				     double* phiHat_dof)=0;
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
				     double epsFactHeaviside,
				     double epsFactDirac,
				     double epsFactDiffusion,
				     double cK,
				     double uL,
				     double uR,
				     double* phin_dof,
				     double* phiHat_dof)=0;
    
  };

  template<class CompKernelType,
	   int nSpace,
	   int nQuadraturePoints_element,
	   int nDOF_mesh_trial_element,
	   int nDOF_trial_element,
	   int nDOF_test_element,
	   int nQuadraturePoints_elementBoundary>
  class VOF : public VOF_base
  {
  public:
    const int nDOF_test_X_trial_element;
    CompKernelType ck;
    VOF():
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

    void FCTStepL2p(int NNZ, //number on non-zero entries on sparsity pattern
		    int numDOFs, //number of DOFs
		    double* lumped_mass_matrix, //lumped mass matrix (as vector)
		    double* solH, //DOFs of high order solution at tnp1
		    double* solL,
		    double* limited_solution,
		    int* csrRowIndeces_DofLoops, //csr row indeces 
		    int* csrColumnOffsets_DofLoops, //csr column offsets 
		    double* MassMatrix //mass matrix
		    )
    {
      register double Rpos[numDOFs], Rneg[numDOFs];
      register double FluxCorrectionMatrix[NNZ];
      //////////////////
      // LOOP in DOFs //
      //////////////////
      int ij=0;
      for (int i=0; i<numDOFs; i++)
	{
	  //read some vectors 
	  double solHi = solH[i];
	  double solLi = solL[i];
	  double mi = lumped_mass_matrix[i];

	  double mini=0.0, maxi=1.0; //MQL. USE SIGNED FUNCTION. TMP
	  double Pposi=0, Pnegi=0;
	  // LOOP OVER THE SPARSITY PATTERN (j-LOOP)//
	  for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
	    {
	      int j = csrColumnOffsets_DofLoops[offset];
	      // i-th row of flux correction matrix 
	      FluxCorrectionMatrix[ij] = ((i==j ? 1. : 0.)*mi - MassMatrix[ij]) * (solH[j]-solHi);

	      ///////////////////////
	      // COMPUTE P VECTORS //
	      ///////////////////////
	      Pposi += FluxCorrectionMatrix[ij]*((FluxCorrectionMatrix[ij] > 0) ? 1. : 0.);
	      Pnegi += FluxCorrectionMatrix[ij]*((FluxCorrectionMatrix[ij] < 0) ? 1. : 0.);

	      //update ij 
	      ij+=1;
	    }
	  ///////////////////////
	  // COMPUTE Q VECTORS //
	  ///////////////////////
	  double Qposi = mi*(maxi-solLi);
	  double Qnegi = mi*(mini-solLi);

	  ///////////////////////
	  // COMPUTE R VECTORS //
	  ///////////////////////
	  Rpos[i] = ((Pposi==0) ? 1. : std::min(1.0,Qposi/Pposi));
	  Rneg[i] = ((Pnegi==0) ? 1. : std::min(1.0,Qnegi/Pnegi));
	} // i DOFs
      
      //////////////////////
      // COMPUTE LIMITERS // 
      //////////////////////
      ij=0;
      for (int i=0; i<numDOFs; i++)
	{
	  double ith_Limiter_times_FluxCorrectionMatrix = 0.;
	  double Rposi = Rpos[i], Rnegi = Rneg[i];
	  // LOOP OVER THE SPARSITY PATTERN (j-LOOP)//
	  for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
	    {
	      int j = csrColumnOffsets_DofLoops[offset];
	      ith_Limiter_times_FluxCorrectionMatrix += 
		((FluxCorrectionMatrix[ij]>0) ? std::min(Rposi,Rneg[j]) : std::min(Rnegi,Rpos[j])) 
		* FluxCorrectionMatrix[ij];
	      //ith_Limiter_times_FluxCorrectionMatrix += FluxCorrectionMatrix[ij];
	      //update ij
	      ij+=1;
	    }
	  //limited_solution[i] = fmax(0.0,solL[i] + 1./lumped_mass_matrix[i]*ith_Limiter_times_FluxCorrectionMatrix);
	  limited_solution[i] = solL[i] + 1./lumped_mass_matrix[i]*ith_Limiter_times_FluxCorrectionMatrix;
	}
    }
    
    void FCTStep(int NNZ, //number on non-zero entries on sparsity pattern
		 int numDOFs, //number of DOFs
		 double* lumped_mass_matrix, //lumped mass matrix (as vector)
		 double* soln, //DOFs of solution at time tn
		 double* solH, //DOFs of high order solution at tnp1
		 double* low_order_solution,
		 double* limited_solution,
		 int* csrRowIndeces_DofLoops, //csr row indeces 
		 int* csrColumnOffsets_DofLoops, //csr column offsets 
		 double* MassMatrix, //mass matrix
		 double* dt_times_dH_minus_dL, //low minus high order dissipative matrices
		 double* min_u_bc, //min/max value at BCs. If DOF is not at boundary then min=1E10, max=-1E10
		 double* max_u_bc, 
		 int LUMPED_MASS_MATRIX
		 )
    {
      register double Rpos[numDOFs], Rneg[numDOFs];
      register double FluxCorrectionMatrix[NNZ];
      register double solL[numDOFs];
      //////////////////
      // LOOP in DOFs //
      //////////////////
      int ij=0;
      for (int i=0; i<numDOFs; i++)
	{
	  //read some vectors 
	  double solHi = solH[i];
	  double solni = soln[i];
	  double mi = lumped_mass_matrix[i];
	  // compute low order solution
	  // mi*(uLi-uni) + dt*sum_j[(Tij+dLij)*unj] = 0
	  solL[i] = low_order_solution[i];

	  double mini=min_u_bc[i], maxi=max_u_bc[i]; // init min/max with value at BCs (NOTE: if no boundary then min=1E10, max=-1E10)
	  if (GLOBAL_FCT==1)
	    {
	      mini = 0.;
	      maxi = 1.;
	    }

	  double Pposi=0, Pnegi=0;
	  // LOOP OVER THE SPARSITY PATTERN (j-LOOP)//
	  for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
	    {
	      int j = csrColumnOffsets_DofLoops[offset];
	      ////////////////////////
	      // COMPUTE THE BOUNDS //
	      ////////////////////////
	      if (GLOBAL_FCT == 0)
		{
		  mini = fmin(mini,soln[j]);
		  maxi = fmax(maxi,soln[j]);
		}      
	      // i-th row of flux correction matrix 
	      double ML_minus_MC = (LUMPED_MASS_MATRIX == 1 ? 0. : (i==j ? 1. : 0.)*mi - MassMatrix[ij]);
	      FluxCorrectionMatrix[ij] = ML_minus_MC * (solH[j]-soln[j] - (solHi-solni)) 
		+ dt_times_dH_minus_dL[ij]*(soln[j]-solni);

	      ///////////////////////
	      // COMPUTE P VECTORS //
	      ///////////////////////
	      Pposi += FluxCorrectionMatrix[ij]*((FluxCorrectionMatrix[ij] > 0) ? 1. : 0.);
	      Pnegi += FluxCorrectionMatrix[ij]*((FluxCorrectionMatrix[ij] < 0) ? 1. : 0.);

	      //update ij 
	      ij+=1;
	    }
	  ///////////////////////
	  // COMPUTE Q VECTORS //
	  ///////////////////////
	  double Qposi = mi*(maxi-solL[i]);
	  double Qnegi = mi*(mini-solL[i]);

	  ///////////////////////
	  // COMPUTE R VECTORS //
	  ///////////////////////
	  Rpos[i] = ((Pposi==0) ? 1. : fmin(1.0,Qposi/Pposi));
	  Rneg[i] = ((Pnegi==0) ? 1. : fmin(1.0,Qnegi/Pnegi));
	} // i DOFs
      
      //////////////////////
      // COMPUTE LIMITERS // 
      //////////////////////
      ij=0;
      for (int i=0; i<numDOFs; i++)
	{
	  double ith_Limiter_times_FluxCorrectionMatrix = 0.;
	  double Rposi = Rpos[i], Rnegi = Rneg[i];
	  // LOOP OVER THE SPARSITY PATTERN (j-LOOP)//
	  for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
	    {
	      int j = csrColumnOffsets_DofLoops[offset];
	      ith_Limiter_times_FluxCorrectionMatrix += 
	      ((FluxCorrectionMatrix[ij]>0) ? fmin(Rposi,Rneg[j]) : fmin(Rnegi,Rpos[j])) 
		* FluxCorrectionMatrix[ij];
	      //ith_Limiter_times_FluxCorrectionMatrix += FluxCorrectionMatrix[ij];
	      //update ij
	      ij+=1;
	    }
	  limited_solution[i] = solL[i] + 1./lumped_mass_matrix[i]*ith_Limiter_times_FluxCorrectionMatrix;	  
	}
    }

    
    void calculateResidual(//element
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
			   int lag_shockCapturing, /*mwf not used yet*/
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
			   double epsFactHeaviside,
			   double epsFactDirac,
			   double epsFactDiffusion,
			   double* phin_dof,
			   double* phiHat_dof,
			   // AUX QUANTITIES OF INTEREST 
			   double* quantDOFs)
    {
      double Ct_sge = 4.0;	  
      //
      //loop over elements to compute volume integrals and load them into element and global residual
      //
      //eN is the element index
      //eN_k is the quadrature point index for a scalar
      //eN_k_nSpace is the quadrature point index for a vector
      //eN_i is the element test function index
      //eN_j is the element trial function index
      //eN_k_j is the quadrature point index for a trial function
      //eN_k_i is the quadrature point index for a trial function
      for(int eN=0;eN<nElements_global;eN++)
	{
	  //declare local storage for element residual and initialize
	  register double elementResidual_u[nDOF_test_element];
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      elementResidual_u[i]=0.0;
	    }//i
	  //loop over quadrature points and compute integrands
	  for  (int k=0;k<nQuadraturePoints_element;k++)
	    {
	      //compute indeces and declare local storage
	      register int eN_k = eN*nQuadraturePoints_element+k,
		eN_k_nSpace = eN_k*nSpace,
		eN_nDOF_trial_element = eN*nDOF_trial_element;
	      register double u=0.0,grad_u[nSpace],grad_u_old[nSpace],
		m=0.0,dm=0.0,
		f[nSpace],df[nSpace],
		m_t=0.0,dm_t=0.0,
		pdeResidual_u=0.0,
		Lstar_u[nDOF_test_element],
		subgridError_u=0.0,
		tau=0.0,tau0=0.0,tau1=0.0,
		numDiff0=0.0,numDiff1=0.0,
		jac[nSpace*nSpace],
		jacDet,
		jacInv[nSpace*nSpace],
		u_grad_trial[nDOF_trial_element*nSpace],
		u_test_dV[nDOF_trial_element],
		u_grad_test_dV[nDOF_test_element*nSpace],
		dV,x,y,z,xt,yt,zt,
		//VRANS
		porosity,
		//
		G[nSpace*nSpace],G_dd_G,tr_G;//norm_Rv;
	      // //
	      // //compute solution and gradients at quadrature points
	      // //
	      // u=0.0;
	      // for (int I=0;I<nSpace;I++)
	      //   {
	      //     grad_u[I]=0.0;
	      //   }
	      // for (int j=0;j<nDOF_trial_element;j++)
	      //   {
	      //     int eN_j=eN*nDOF_trial_element+j;
	      //     int eN_k_j=eN_k*nDOF_trial_element+j;
	      //     int eN_k_j_nSpace = eN_k_j*nSpace;
	      //     u += valFromDOF_c(u_dof[u_l2g[eN_j]],u_trial[eN_k_j]);
	      //     for (int I=0;I<nSpace;I++)
	      //       {
	      //         grad_u[I] += gradFromDOF_c(u_dof[u_l2g[eN_j]],u_grad_trial[eN_k_j_nSpace+I]);
	      //       }
	      //   }
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
	      ck.calculateG(jacInv,G,G_dd_G,tr_G);
	      //get the trial function gradients
	      ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,u_grad_trial);
	      //get the solution
	      ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],u);
	      //get the solution gradients
	      ck.gradFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_u);
	      ck.gradFromDOF(u_dof_old,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_u_old);
	      //precalculate test function products with integration weights
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
		  for (int I=0;I<nSpace;I++)
		    {
		      u_grad_test_dV[j*nSpace+I]   = u_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
		    }
		}
	      //VRANS
	      porosity = q_porosity[eN_k];
	      //
	      //
	      //calculate pde coefficients at quadrature points
	      //
	      evaluateCoefficients(&velocity[eN_k_nSpace],
				   u,
				   //VRANS
				   porosity,
				   //
				   m,
				   dm,
				   f,
				   df);
	      //
	      //moving mesh
	      //
	      double mesh_velocity[3];
	      mesh_velocity[0] = xt;
	      mesh_velocity[1] = yt;
	      mesh_velocity[2] = zt;
	      //std::cout<<"q mesh_velocity"<<std::endl;
	      for (int I=0;I<nSpace;I++)
		{
		  //std::cout<<mesh_velocity[I]<<std::endl;
		  f[I] -= MOVING_DOMAIN*m*mesh_velocity[I];
		  df[I] -= MOVING_DOMAIN*dm*mesh_velocity[I];
		}
	      //
	      //calculate time derivative at quadrature points
	      //
              if (q_dV_last[eN_k] <= -100)
                q_dV_last[eN_k] = dV;
              q_dV[eN_k] = dV;
	      ck.bdf(alphaBDF,
		     q_m_betaBDF[eN_k]*q_dV_last[eN_k]/dV,//ensure prior mass integral is correct for  m_t with BDF1
		     m,
		     dm,
		     m_t,
		     dm_t);
	      //
	      //calculate subgrid error (strong residual and adjoint)
	      //
	      //calculate strong residual
	      pdeResidual_u = ck.Mass_strong(m_t) + ck.Advection_strong(df,grad_u);
	      //calculate adjoint
	      for (int i=0;i<nDOF_test_element;i++)
		{
		  // register int eN_k_i_nSpace = (eN_k*nDOF_trial_element+i)*nSpace;
		  // Lstar_u[i]  = ck.Advection_adjoint(df,&u_grad_test_dV[eN_k_i_nSpace]);
		  register int i_nSpace = i*nSpace;
		  Lstar_u[i]  = ck.Advection_adjoint(df,&u_grad_test_dV[i_nSpace]);
		}
	      //calculate tau and tau*Res
	      calculateSubgridError_tau(elementDiameter[eN],dm_t,df,cfl[eN_k],tau0);
              calculateSubgridError_tau(Ct_sge,
                                        G,
					dm_t,
					df,
					tau1,
				        cfl[eN_k]);
					
              tau = useMetrics*tau1+(1.0-useMetrics)*tau0;

	      subgridError_u = -tau*pdeResidual_u;
	      //
	      //calculate shock capturing diffusion
	      //
	      ck.calculateNumericalDiffusion(shockCapturingDiffusion,elementDiameter[eN],pdeResidual_u,grad_u,numDiff0);	      
	      //ck.calculateNumericalDiffusion(shockCapturingDiffusion,G,pdeResidual_u,grad_u_old,numDiff1);
	      ck.calculateNumericalDiffusion(shockCapturingDiffusion,sc_uref, sc_alpha,G,G_dd_G,pdeResidual_u,grad_u,numDiff1);
	      q_numDiff_u[eN_k] = useMetrics*numDiff1+(1.0-useMetrics)*numDiff0;
              //std::cout<<tau<<"   "<<q_numDiff_u[eN_k]<<'\t'<<numDiff0<<'\t'<<numDiff1<<'\t'<<pdeResidual_u<<std::endl;
	      // 
	      //update element residual 
	      // 
	      /*	      std::cout<<m_t<<'\t'
		       <<f[0]<<'\t'
		       <<f[1]<<'\t'
		       <<df[0]<<'\t'
		       <<df[1]<<'\t'
		       <<subgridError_u<<'\t'
		       <<q_numDiff_u_last[eN_k]<<std::endl;*/
	    
	      for(int i=0;i<nDOF_test_element;i++) 
		{ 
		  //register int eN_k_i=eN_k*nDOF_test_element+i,
		    //eN_k_i_nSpace = eN_k_i*nSpace,
		   register int i_nSpace=i*nSpace;
		   elementResidual_u[i] += ck.Mass_weak(m_t,u_test_dV[i]) + 
		     ck.Advection_weak(f,&u_grad_test_dV[i_nSpace]) + 
		     ck.SubgridError(subgridError_u,Lstar_u[i]) + 
		     ck.NumericalDiffusion(q_numDiff_u_last[eN_k],grad_u,&u_grad_test_dV[i_nSpace]);		   
		}//i
	      //
	      //cek/ido todo, get rid of m, since u=m
	      //save momentum for time history and velocity for subgrid error
	      //save solution for other models 
	      //
	      q_u[eN_k] = u;
	      q_m[eN_k] = m;
	    }
	  //
	  //load element into global residual and save element residual
	  //
	  for(int i=0;i<nDOF_test_element;i++) 
	    { 
	      register int eN_i=eN*nDOF_test_element+i;          
	      globalResidual[offset_u+stride_u*u_l2g[eN_i]] += elementResidual_u[i];
	    }//i
	}//elements
      //
      //loop over exterior element boundaries to calculate surface integrals and load into element and global residuals
      //
      //ebNE is the Exterior element boundary INdex
      //ebN is the element boundary INdex
      //eN is the element index
      for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) 
	{ 
	  register int ebN = exteriorElementBoundariesArray[ebNE], 
	    eN  = elementBoundaryElementsArray[ebN*2+0],
	    ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+0],
	    eN_nDOF_trial_element = eN*nDOF_trial_element;
	  register double elementResidual_u[nDOF_test_element];
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      elementResidual_u[i]=0.0;
	    }
	  for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) 
	    { 
	      register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
		ebNE_kb_nSpace = ebNE_kb*nSpace,
		ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb,
		ebN_local_kb_nSpace = ebN_local_kb*nSpace;
	      register double u_ext=0.0,
		grad_u_ext[nSpace],
		m_ext=0.0,
		dm_ext=0.0,
		f_ext[nSpace],
		df_ext[nSpace],
		flux_ext=0.0,
		bc_u_ext=0.0,
		//bc_grad_u_ext[nSpace],
		bc_m_ext=0.0,
		bc_dm_ext=0.0,
		bc_f_ext[nSpace],
		bc_df_ext[nSpace],
		jac_ext[nSpace*nSpace],
		jacDet_ext,
		jacInv_ext[nSpace*nSpace],
		boundaryJac[nSpace*(nSpace-1)],
		metricTensor[(nSpace-1)*(nSpace-1)],
		metricTensorDetSqrt,
		dS,
		u_test_dS[nDOF_test_element],
		u_grad_trial_trace[nDOF_trial_element*nSpace],
		normal[nSpace],x_ext,y_ext,z_ext,xt_ext,yt_ext,zt_ext,integralScaling,
		//VRANS
		porosity_ext,
		//
		G[nSpace*nSpace],G_dd_G,tr_G;
	      // 
	      //calculate the solution and gradients at quadrature points 
	      // 
	      //compute information about mapping from reference element to physical element
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
	      //std::cout<<"metricTensorDetSqrt "<<metricTensorDetSqrt<<" integralScaling "<<integralScaling<<std::endl;
	      dS = ((1.0-MOVING_DOMAIN)*metricTensorDetSqrt + MOVING_DOMAIN*integralScaling)*dS_ref[kb];
	      //get the metric tensor
	      //cek todo use symmetry
	      ck.calculateG(jacInv_ext,G,G_dd_G,tr_G);
	      //compute shape and solution information
	      //shape
	      ck.gradTrialFromRef(&u_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,u_grad_trial_trace);
	      //solution and gradients	
	      ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_trace_ref[ebN_local_kb*nDOF_test_element],u_ext);
	      ck.gradFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],u_grad_trial_trace,grad_u_ext);
	      //precalculate test function products with integration weights
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  u_test_dS[j] = u_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;
		}
	      //
	      //load the boundary values
	      //
	      bc_u_ext = isDOFBoundary_u[ebNE_kb]*ebqe_bc_u_ext[ebNE_kb]+(1-isDOFBoundary_u[ebNE_kb])*u_ext;
	      //VRANS
	      porosity_ext = ebqe_porosity_ext[ebNE_kb];
	      //
	      // 
	      //calculate the pde coefficients using the solution and the boundary values for the solution 
	      // 
	      evaluateCoefficients(&ebqe_velocity_ext[ebNE_kb_nSpace],
				   u_ext,
				   //VRANS
				   porosity_ext,
				   //
				   m_ext,
				   dm_ext,
				   f_ext,
				   df_ext);
	      evaluateCoefficients(&ebqe_velocity_ext[ebNE_kb_nSpace],
				   bc_u_ext,
				   //VRANS
				   porosity_ext,
				   //
				   bc_m_ext,
				   bc_dm_ext,
				   bc_f_ext,
				   bc_df_ext);    
	      //
	      //moving mesh
	      //
	      double mesh_velocity[3];
	      mesh_velocity[0] = xt_ext;
	      mesh_velocity[1] = yt_ext;
	      mesh_velocity[2] = zt_ext;
	      //std::cout<<"mesh_velocity ext"<<std::endl;
	      for (int I=0;I<nSpace;I++)
		{
		  //std::cout<<mesh_velocity[I]<<std::endl;
		  f_ext[I] -= MOVING_DOMAIN*m_ext*mesh_velocity[I];
		  df_ext[I] -= MOVING_DOMAIN*dm_ext*mesh_velocity[I];
		  bc_f_ext[I] -= MOVING_DOMAIN*bc_m_ext*mesh_velocity[I];
		  bc_df_ext[I] -= MOVING_DOMAIN*bc_dm_ext*mesh_velocity[I];
		}
	      // 
	      //calculate the numerical fluxes 
	      // 
	      exteriorNumericalAdvectiveFlux(isDOFBoundary_u[ebNE_kb],
					     isFluxBoundary_u[ebNE_kb],
					     normal,
					     bc_u_ext,
					     ebqe_bc_flux_u_ext[ebNE_kb],
					     u_ext,//smoothedHeaviside(eps,ebqe_phi[ebNE_kb]),//cek hack
					     df_ext,//VRANS includes porosity
					     flux_ext);
	      ebqe_flux[ebNE_kb] = flux_ext;
	      //save for other models? cek need to be consistent with numerical flux
	      if(flux_ext >=0.0)
		ebqe_u[ebNE_kb] = u_ext;
	      else
		ebqe_u[ebNE_kb] = bc_u_ext;
	      //
	      //update residuals
	      //
	      for (int i=0;i<nDOF_test_element;i++)
		{
		  //int ebNE_kb_i = ebNE_kb*nDOF_test_element+i;

		  elementResidual_u[i] += ck.ExteriorElementBoundaryFlux(flux_ext,u_test_dS[i]);
		}//i
	    }//kb
	  //
	  //update the element and global residual storage
	  //
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      int eN_i = eN*nDOF_test_element+i;

	      globalResidual[offset_u+stride_u*u_l2g[eN_i]] += elementResidual_u[i];
	    }//i
	}//ebNE
    }
    
    void calculateResidual_entropy_viscosity(//element
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
					     double epsFactHeaviside,
					     double epsFactDirac,
					     double epsFactDiffusion,
					     double* phin_dof,
					     double* phiHat_dof,
					     // AUX QUANTITIES OF INTEREST 
					     double* quantDOFs)
    {
      // NOTE: This function follows a different (but equivalent) implementation of the smoothness based indicator than NCLS.h
      // Allocate space for the transport matrices
      // This is used for first order KUZMIN'S METHOD
      register double TransportMatrix[NNZ], TransposeTransportMatrix[NNZ];
      for (int i=0; i<NNZ; i++)
	{
	  TransportMatrix[i] = 0.;
	  TransposeTransportMatrix[i] = 0.;
	}

      // compute entropy and init global_entropy_residual and boundary_integral
      register double psi[numDOFs], eta[numDOFs], global_entropy_residual[numDOFs], boundary_integral[numDOFs];
      for (int i=0; i<numDOFs; i++)
	{
	  // NODAL ENTROPY //
	  if (STABILIZATION_TYPE==1) //EV stab
	    {
	      double porosity_times_solni = porosity_dof[i]*u_dof_old[i];
	      eta[i] = ENTROPY_TYPE == 1 ? ENTROPY(porosity_times_solni,uL,uR) : ENTROPY_LOG(porosity_times_solni,uL,uR);
	      global_entropy_residual[i]=0.;
	    }
	  boundary_integral[i]=0.;
	} 
	  
      //////////////////////////////////////////////
      // ** LOOP IN CELLS FOR CELL BASED TERMS ** //
      //////////////////////////////////////////////
      // HERE WE COMPUTE: 
      //    * Time derivative term. porosity*u_t
      //    * cell based CFL (for reference)
      //    * Entropy residual
      //    * Transport matrices
      for(int eN=0;eN<nElements_global;eN++)
	{
	  //declare local storage for local contributions and initialize
	  register double 
	    elementResidual_u[nDOF_test_element], 
	    element_entropy_residual[nDOF_test_element];
	  register double  elementTransport[nDOF_test_element][nDOF_trial_element];
	  register double  elementTransposeTransport[nDOF_test_element][nDOF_trial_element];
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      elementResidual_u[i]=0.0;
	      element_entropy_residual[i]=0.0;
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  elementTransport[i][j]=0.0;
		  elementTransposeTransport[i][j]=0.0;
		}
	    }
	  //loop over quadrature points and compute integrands
	  for  (int k=0;k<nQuadraturePoints_element;k++)
	    {
	      //compute indeces and declare local storage
	      register int eN_k = eN*nQuadraturePoints_element+k,
		eN_k_nSpace = eN_k*nSpace,
		eN_nDOF_trial_element = eN*nDOF_trial_element;
	      register double 
		// for entropy residual
		aux_entropy_residual=0., DENTROPY_un, DENTROPY_uni,
		//for mass matrix contributions
		u=0.0, un=0.0, grad_un[nSpace], porosity_times_velocity[nSpace], 
		u_test_dV[nDOF_trial_element], 
		u_grad_trial[nDOF_trial_element*nSpace], 
		u_grad_test_dV[nDOF_test_element*nSpace],
		//for general use
		jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace],
		dV,x,y,z,xt,yt,zt,
		//VRANS
		porosity;
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
	      //get the solution (of Newton's solver). To compute time derivative term
	      ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],u);
	      //get the solution at quad point at tn and tnm1 for entropy viscosity
	      ck.valFromDOF(u_dof_old,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],un);
	      //get the solution gradients at tn for entropy viscosity
	      ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,u_grad_trial);
	      ck.gradFromDOF(u_dof_old,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_un);
	      
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
	      //VRANS
	      porosity = q_porosity[eN_k];
	      //
	      //moving mesh
	      //
	      double mesh_velocity[3];
	      mesh_velocity[0] = xt;
	      mesh_velocity[1] = yt;
	      mesh_velocity[2] = zt;
	      //relative velocity at tn
	      for (int I=0;I<nSpace;I++)
		porosity_times_velocity[I] = porosity*(velocity[eN_k_nSpace+I]-MOVING_DOMAIN*mesh_velocity[I]);

	      //////////////////////////////
	      // CALCULATE CELL BASED CFL //
	      //////////////////////////////
	      calculateCFL(elementDiameter[eN]/degree_polynomial,porosity_times_velocity,cfl[eN_k]); 
	      
	      //////////////////////////////////////////////
	      // CALCULATE ENTROPY RESIDUAL AT QUAD POINT //
	      //////////////////////////////////////////////
	      if (STABILIZATION_TYPE==1) // EV stab
		{
		  for (int I=0;I<nSpace;I++)
		    aux_entropy_residual += porosity_times_velocity[I]*grad_un[I];
		  DENTROPY_un = ENTROPY_TYPE==1 ? DENTROPY(porosity*un,uL,uR) : DENTROPY_LOG(porosity*un,uL,uR);
		}
	      //////////////
	      // ith-LOOP //
	      //////////////	      
	      for(int i=0;i<nDOF_test_element;i++) 
		{ 
		  // VECTOR OF ENTROPY RESIDUAL //
		  int eN_i=eN*nDOF_test_element+i;
		  if (STABILIZATION_TYPE==1) // EV stab
		    {
		      int gi = offset_u+stride_u*u_l2g[eN_i]; //global i-th index
		      double porosity_times_uni = porosity_dof[gi]*u_dof_old[gi];
		      DENTROPY_uni = ENTROPY_TYPE == 1 ? DENTROPY(porosity_times_uni,uL,uR) : DENTROPY_LOG(porosity_times_uni,uL,uR);
		      element_entropy_residual[i] += (DENTROPY_un - DENTROPY_uni)*aux_entropy_residual*u_test_dV[i];
		    }
		  elementResidual_u[i] += porosity*(u-un)*u_test_dV[i];		  
		  ///////////////
		  // j-th LOOP // To construct transport matrices
		  ///////////////
		  for(int j=0;j<nDOF_trial_element;j++) 
		    { 
		      int j_nSpace = j*nSpace;
		      int i_nSpace = i*nSpace;
		      elementTransport[i][j] += // -int[(vel.grad_wi)*wj*dx]
			ck.AdvectionJacobian_weak(porosity_times_velocity,
						  u_trial_ref[k*nDOF_trial_element+j],&u_grad_test_dV[i_nSpace]);
		      elementTransposeTransport[i][j] += // -int[(vel.grad_wj)*wi*dx]
			ck.AdvectionJacobian_weak(porosity_times_velocity,
						  u_trial_ref[k*nDOF_trial_element+i],&u_grad_test_dV[j_nSpace]);
		    }
		}//i
	      //save solution for other models 
	      q_u[eN_k] = u;
	      q_m[eN_k] = porosity*u;
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
	      // distribute entropy_residual
	      if (STABILIZATION_TYPE==1) // EV Stab
		global_entropy_residual[gi] += element_entropy_residual[i];

	      // distribute transport matrices
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  int eN_i_j = eN_i*nDOF_trial_element+j;
		  TransportMatrix[csrRowIndeces_CellLoops[eN_i] + csrColumnOffsets_CellLoops[eN_i_j]] 
		    += elementTransport[i][j];
		  TransposeTransportMatrix[csrRowIndeces_CellLoops[eN_i] + csrColumnOffsets_CellLoops[eN_i_j]] 
		    += elementTransposeTransport[i][j];
		}//j
	    }//i
	}//elements

      //////////////////////////////////////////////////////////////////////////////////////////
      // ADD OUTFLOW BOUNDARY TERM TO TRANSPORT MATRICES AND COMPUTE INFLOW BOUNDARY INTEGRAL //
      //////////////////////////////////////////////////////////////////////////////////////////
      //   * Compute outflow boundary integral as a matrix; i.e., int_B[ (vel.normal)*wi*wj*dx]
      for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) 
	{ 
	  double min_u_bc_local = 1E10, max_u_bc_local = -1E10;
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
		ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb,
		ebN_local_kb_nSpace = ebN_local_kb*nSpace;	      
	      register double 
		u_ext=0.0, bc_u_ext=0.0,
		porosity_times_velocity[nSpace],
		flux_ext=0.0, dflux_ext=0.0,
		fluxTransport[nDOF_trial_element],
		jac_ext[nSpace*nSpace],
		jacDet_ext,
		jacInv_ext[nSpace*nSpace],
		boundaryJac[nSpace*(nSpace-1)],
		metricTensor[(nSpace-1)*(nSpace-1)],
		metricTensorDetSqrt,
		dS,
		u_test_dS[nDOF_test_element],
		normal[nSpace],x_ext,y_ext,z_ext,xt_ext,yt_ext,zt_ext,integralScaling,porosity_ext;
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
	      dS = ((1.0-MOVING_DOMAIN)*metricTensorDetSqrt + MOVING_DOMAIN*integralScaling)*dS_ref[kb];
	      //compute shape and solution information
	      ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_trace_ref[ebN_local_kb*nDOF_test_element],u_ext);
	      //precalculate test function products with integration weights
	      for (int j=0;j<nDOF_trial_element;j++)
		u_test_dS[j] = u_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;

	      //VRANS
	      porosity_ext = ebqe_porosity_ext[ebNE_kb];
	      //
	      //moving mesh
	      //
	      double mesh_velocity[3];
	      mesh_velocity[0] = xt_ext;
	      mesh_velocity[1] = yt_ext;
	      mesh_velocity[2] = zt_ext;
	      //std::cout<<"mesh_velocity ext"<<std::endl;
	      for (int I=0;I<nSpace;I++)
		porosity_times_velocity[I] = porosity_ext*(ebqe_velocity_ext[ebNE_kb_nSpace+I] - MOVING_DOMAIN*mesh_velocity[I]);
	      //
	      //calculate the fluxes
	      //
	      double flow = 0.;
	      for (int I=0; I < nSpace; I++)
		flow += normal[I]*porosity_times_velocity[I];

	      if (flow >= 0) //outflow. This is handled via the transport matrices. Then flux_ext=0 and dflux_ext!=0
		{
		  dflux_ext = flow; 
		  flux_ext = 0;
		  // save external u
		  ebqe_u[ebNE_kb] = u_ext;
 		}
	      else // inflow. This is handled via the boundary integral. Then flux_ext!=0 and dflux_ext=0
		{
		  dflux_ext = 0;
		  // save external u
		  ebqe_u[ebNE_kb] = isDOFBoundary_u[ebNE_kb]*ebqe_bc_u_ext[ebNE_kb]+(1-isDOFBoundary_u[ebNE_kb])*u_ext;
		  if (isDOFBoundary_u[ebNE_kb] == 1)
		    flux_ext = ebqe_bc_u_ext[ebNE_kb]*flow;
		  else if (isFluxBoundary_u[ebNE_kb] == 1)
		    flux_ext = ebqe_bc_flux_u_ext[ebNE_kb];
		  else
		    {
		      std::cout<<"warning: VOF open boundary with no external trace, setting to zero for inflow"<<std::endl;
		      flux_ext = 0.0;
		    }
		}
	      
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  // elementResidual. This is to include the inflow boundary integral. 
		  // NOTE: here I assume that we use a Galerkin approach st nDOF_test_element = nDOF_trial_element
		  elementResidual_u[j] += flux_ext*u_test_dS[j];
		  register int ebN_local_kb_j=ebN_local_kb*nDOF_trial_element+j;
		  fluxTransport[j] = dflux_ext*u_trial_trace_ref[ebN_local_kb_j];
		}//j
	      ///////////////////////////////////////////////////////
	      // DISTRIBUTE OUTFLOW BOUNDARY TO TRANSPORT MATRICES //
	      ///////////////////////////////////////////////////////
	      for (int i=0;i<nDOF_test_element;i++)
		{
		  register int eN_i = eN*nDOF_test_element+i;
		  for (int j=0;j<nDOF_trial_element;j++)
		    {
		      register int ebN_i_j = ebN*4*nDOF_test_X_trial_element + i*nDOF_trial_element + j;
		      TransportMatrix[csrRowIndeces_CellLoops[eN_i] + csrColumnOffsets_eb_CellLoops[ebN_i_j]] 
			+= fluxTransport[j]*u_test_dS[i];
		      TransposeTransportMatrix[csrRowIndeces_CellLoops[eN_i] + csrColumnOffsets_eb_CellLoops[ebN_i_j]] 
			+= fluxTransport[i]*u_test_dS[j];
		    }//j
		}//i	      
	      // local min/max at boundary
	      min_u_bc_local = fmin(ebqe_u[ebNE_kb], min_u_bc_local);
	      max_u_bc_local = fmax(ebqe_u[ebNE_kb], max_u_bc_local);
	    }//kb	  
	  // global min/max at boundary 
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      int eN_i = eN*nDOF_test_element+i;
	      int gi = offset_u+stride_u*u_l2g[eN_i]; //global i-th index
	      globalResidual[gi] += dt*elementResidual_u[i];
	      boundary_integral[gi] += elementResidual_u[i];	      
	      min_u_bc[gi] = fmin(min_u_bc_local,min_u_bc[gi]);
	      max_u_bc[gi] = fmax(max_u_bc_local,max_u_bc[gi]);
	    }
	}//ebNE
      // END OF ADDING BOUNDARY TERM TO TRANSPORT MATRICES and COMPUTING BOUNDARY INTEGRAL //

      /////////////////////////////////////////////////////////////////
      // COMPUTE SMOOTHNESS INDICATOR and NORMALIZE ENTROPY RESIDUAL //
      /////////////////////////////////////////////////////////////////
      // NOTE: see NCLS.h for a different but equivalent implementation of this. 
      int ij = 0;
      for (int i=0; i<numDOFs; i++)
	{
	  double gi[nSpace], Cij[nSpace], xi[nSpace], etaMaxi, etaMini;
	  if (STABILIZATION_TYPE==1) //EV Stabilization
	    {
	      // For eta min and max
	      etaMaxi = fabs(eta[i]);
	      etaMini = fabs(eta[i]);	  
	    }
	  double porosity_times_solni = porosity_dof[i]*u_dof_old[i];
	  // initialize gi and compute xi
	  for (int I=0; I < nSpace; I++)
	    {
	      gi[I] = 0.;
	      xi[I] = mesh_dof[i*3+I];
	    }
	  // for smoothness indicator // 
	  double alpha_numerator_pos = 0., alpha_numerator_neg = 0., alpha_denominator_pos = 0., alpha_denominator_neg = 0.;
	  for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
	    { // First loop in j (sparsity pattern)
	      int j = csrColumnOffsets_DofLoops[offset];
	      if (STABILIZATION_TYPE==1) //EV Stabilization
		{
		  // COMPUTE ETA MIN AND ETA MAX // 
		  etaMaxi = fmax(etaMaxi,fabs(eta[j]));
		  etaMini = fmin(etaMini,fabs(eta[j]));
		}
	      double porosity_times_solnj = porosity_dof[j]*u_dof_old[j];
	      // Update Cij matrices
	      Cij[0] = Cx[ij];
	      Cij[1] = Cy[ij];
#if nSpace == 3
	      Cij[2] = Cz[ij];
#endif
	      // COMPUTE gi VECTOR. gi=1/mi*sum_j(Cij*solj)
	      for (int I=0; I < nSpace; I++)
		gi[I] += Cij[I]*porosity_times_solnj;

	      // COMPUTE numerator and denominator of smoothness indicator
	      double alpha_num = porosity_times_solni - porosity_times_solnj;
	      if (alpha_num >= 0.)
		{
		  alpha_numerator_pos += alpha_num;
		  alpha_denominator_pos += alpha_num;
		}
	      else
		{
		  alpha_numerator_neg += alpha_num;
		  alpha_denominator_neg += fabs(alpha_num);
		}
	      //update ij
	      ij+=1;
	    }
	  // scale g vector by lumped mass matrix
	  for (int I=0; I < nSpace; I++)
	    gi[I] /= ML[i];
	  if (STABILIZATION_TYPE==1) //EV Stab
	    {
	      // Normalizae entropy residual 
	      global_entropy_residual[i] *= etaMini == etaMaxi ? 0. : 2*cE/(etaMaxi-etaMini);
	      quantDOFs[i] = fabs(global_entropy_residual[i]); 
	    }

	  // Now that I have the gi vectors, I can use them for the current i-th DOF
	  double SumPos=0., SumNeg=0.;
	  for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
	    { // second loop in j (sparsity pattern)
	      int j = csrColumnOffsets_DofLoops[offset];
	      // compute xj
	      double xj[nSpace];
	      for (int I=0; I < nSpace; I++)
		xj[I] = mesh_dof[j*3+I];
	      // compute gi*(xi-xj)
	      double gi_times_x=0.;
	      for (int I=0; I < nSpace; I++)
		gi_times_x += gi[I]*(xi[I]-xj[I]);
	      // compute the positive and negative part of gi*(xi-xj)
	      SumPos += gi_times_x > 0 ? gi_times_x : 0;
	      SumNeg += gi_times_x < 0 ? gi_times_x : 0;
	    }
	  double sigmaPosi = fmin(1.,(fabs(SumNeg)+1E-15)/(SumPos+1E-15));
	  double sigmaNegi = fmin(1.,(SumPos+1E-15)/(fabs(SumNeg)+1E-15));
	  double alpha_numi = fabs(sigmaPosi*alpha_numerator_pos + sigmaNegi*alpha_numerator_neg);
	  double alpha_deni = sigmaPosi*alpha_denominator_pos + sigmaNegi*alpha_denominator_neg;
	  if (IS_BETAij_ONE == 1)
	    {
	      alpha_numi = fabs(alpha_numerator_pos + alpha_numerator_neg);
	      alpha_deni = alpha_denominator_pos + alpha_denominator_neg;
	    }
	  double alphai = alpha_numi/(alpha_deni+1E-15);
	  quantDOFs[i] = alphai;

	  if (POWER_SMOOTHNESS_INDICATOR==0)
	    psi[i] = 1.0;
	  else
	    psi[i] = std::pow(alphai,POWER_SMOOTHNESS_INDICATOR); //NOTE: they use alpha^2 in the paper	  	  
	}
      /////////////////////////////////////////////
      // ** LOOP IN DOFs FOR EDGE BASED TERMS ** //
      /////////////////////////////////////////////
      ij=0;
      for (int i=0; i<numDOFs; i++)
	{
	  // NOTE: Transport matrices already have the porosity considered. ---> Dissipation matrices as well. 
	  double solni = u_dof_old[i]; // solution at time tn for the ith DOF
	  double porosityi = porosity_dof[i];
	  double ith_dissipative_term = 0;
	  double ith_low_order_dissipative_term = 0;
	  double ith_flux_term = 0;
	  double dLii = 0.;
	  
	  // loop over the sparsity pattern of the i-th DOF
	  for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
	    {
	      int j = csrColumnOffsets_DofLoops[offset];
	      double solnj = u_dof_old[j]; // solution at time tn for the jth DOF
	      double porosityj = porosity_dof[j];
	      double dLowij, dLij, dEVij, dHij;
	      
	      ith_flux_term += TransportMatrix[ij]*solnj;
	      if (i != j) //NOTE: there is really no need to check for i!=j (see formula for ith_dissipative_term)
		{
		  // artificial compression
		  double solij = 0.5*(porosityi*solni+porosityj*solnj);
		  double Compij = cK*fmax(solij*(1.0-solij),0.0)/(fabs(porosityi*solni-porosityj*solnj)+1E-14);
		  // first-order dissipative operator
		  dLowij = fmax(fabs(TransportMatrix[ij]),fabs(TransposeTransportMatrix[ij]));
		  //dLij = fmax(0.,fmax(psi[i]*TransportMatrix[ij], // Approach by S. Badia
		  //		  psi[j]*TransposeTransportMatrix[ij]));
		  dLij = dLowij*fmax(psi[i],psi[j]); // enhance the order to 2nd order. No EV
		  if (STABILIZATION_TYPE==1) //EV Stab
		    {
		      // high-order (entropy viscosity) dissipative operator 		  
		      dEVij = fmax(fabs(global_entropy_residual[i]),fabs(global_entropy_residual[j]));
		      dHij = fmin(dLowij,dEVij) * fmax(1.0-Compij,0.0); // artificial compression
		    }
		  else // smoothness based indicator
		    {
		      dHij = dLij * fmax(1.0-Compij,0.0); // artificial compression
		    }	     	      	     
		  //dissipative terms
		  ith_dissipative_term += dHij*(solnj-solni);
		  ith_low_order_dissipative_term += dLij*(solnj-solni);
		  //dHij - dLij. This matrix is needed during FCT step
		  dt_times_dH_minus_dL[ij] = dt*(dHij - dLij);
		  dLii -= dLij;
		}
	      else //i==j
		{
		  // NOTE: this is incorrect. Indeed, dLii = -sum_{j!=i}(dLij) and similarly for dCii. 
		  // However, it is irrelevant since during the FCT step we do (dL-dC)*(solnj-solni)
		  dt_times_dH_minus_dL[ij]=0;
		}
	      //update ij
	      ij+=1;
	    }
	  double mi = ML[i];
	  // compute edge_based_cfl
	  edge_based_cfl[i] = 2.*fabs(dLii)/mi;
	  low_order_solution[i] = u_dof_old[i] - dt/mi*(ith_flux_term 
							+ boundary_integral[i] 
							- ith_low_order_dissipative_term);
	  // update residual
	  if (LUMPED_MASS_MATRIX==1)
	    globalResidual[i] = u_dof_old[i] - dt/mi*(ith_flux_term + boundary_integral[i] - ith_dissipative_term);
	  else
	    globalResidual[i] += dt*(ith_flux_term - ith_dissipative_term);
	}
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
					  double epsFactHeaviside,
					  double epsFactDirac,
					  double epsFactDiffusion,
					  double* phin_dof,
					  double* phiHat_dof,
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
					  double epsFactHeaviside,
					  double epsFactDirac,
					  double epsFactDiffusion,
					  double* phin_dof,
					  double* phiHat_dof,
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

	      double hK = elementDiameter[eN];
	      double epsHeaviside = epsFactHeaviside*hK;
	      double epsDiffusion = epsFactDiffusion*hK; //kappa = const*h
	      double epsDirac = epsFactDirac*hK;
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
	      calculateCFL(hK,relative_velocity,cfl[eN_k]); 
	      double time_derivative_residual = (smoothedHeaviside(epsHeaviside,phiHatnp1+u)-Hn)/dt;

	      //////////////
	      // ith-LOOP //
	      //////////////
	      double norm_grad_Hn=smoothedDirac(epsDirac,Hn)*std::sqrt(grad_u[0]*grad_u[0]
								       +grad_u[1]*grad_u[1]);
	      double viscosity = 1.0;
	      double diffusion = epsFactDiffusion;
	      //double diffusion = viscosity*fmax(1-cK*(Hn-uL)*(uR-Hn)/(hK*norm_grad_Hn+1E-10),0.0);
	      //double diffusion = viscosity*(1-cK*(Hn-uL)*(uR-Hn))/(hK*norm_grad_Hn+1E-10);
				    
	      for(int i=0;i<nDOF_test_element;i++) 
		{ 
		  register int i_nSpace=i*nSpace;
		  elementResidual_u[i] += 
		    time_derivative_residual*u_test_dV[i]
		    + ck.Advection_weak(f,&u_grad_test_dV[i_nSpace])
		    + ck.NumericalDiffusion(epsDiffusion/dt,grad_u,&u_grad_test_dV[i_nSpace]);
		  //+ ck.NumericalDiffusion(diffusion,grad_u,&u_grad_test_dV[i_nSpace]);
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
		elementResidual_u[i] += smoothedHeaviside(epsHeaviside,phiHat+0*u)*u_test_dV[i];
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
    
    double calculateRhsQuadratureMass(//element
				    double* mesh_trial_ref,//
				    double* mesh_grad_trial_ref,
				    double* mesh_dof, //
				    int* mesh_l2g,//
				    double* dV_ref,//
				    double* u_trial_ref,
				    double* u_test_ref,
				    //physics
				    int nElements_global,//
				    int* u_l2g, //
				    double* elementDiameter,//
				    double* u_dof,//
				    int offset_u, int stride_u, 
				    // PARAMETERS FOR EDGE VISCOSITY 
				    int numDOFs,
				    // FOR FCT
				    double* rhs_mass_correction,
				    double* lumped_L2p,
				    double* lumped_mass_matrix,
				    // FOR NONLINEAR VOF; i.e., MCorr with VOF
				    double epsFactHeaviside,
				    double epsFactDiffusion,
				    double* phiHat_dof)
    {
      double mass = 0;
      for(int eN=0;eN<nElements_global;eN++)
	{
	  double cell_mass = 0.;
	  //declare local storage for local contributions and initialize
	  register double element_rhs_mass_correction[nDOF_test_element];
	  for (int i=0;i<nDOF_test_element;i++)
	    element_rhs_mass_correction[i]=0.0;
	  //loop over quadrature points and compute integrands
	  for  (int k=0;k<nQuadraturePoints_element;k++)
	    {
	      //compute indeces and declare local storage
	      register int eN_k = eN*nQuadraturePoints_element+k,
		eN_k_nSpace = eN_k*nSpace,
		eN_nDOF_trial_element = eN*nDOF_trial_element;
	      register double 
		//for mass matrix contributions
		u, phiHatnp1,
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
	      ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],u);
	      ck.valFromDOF(phiHat_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],phiHatnp1);
	      //precalculate test function products with integration weights for mass matrix terms
	      for (int j=0;j<nDOF_trial_element;j++)
		u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
	      
	      double rhs = smoothedHeaviside(epsFactHeaviside*elementDiameter[eN],
					     phiHatnp1+u);
	      //double rhs = smoothedHeaviside(epsFactHeaviside*elementDiameter[eN],
	      //			     phiHatnp1+epsFactDiffusion*u); //TMP
	      for(int i=0;i<nDOF_test_element;i++) 
		element_rhs_mass_correction[i] += rhs*u_test_dV[i];
	      // compute cell mass
	      cell_mass += rhs*dV;
	    } //k
	  // DISTRIBUTE //
	  mass += cell_mass;
	  for(int i=0;i<nDOF_test_element;i++) 
	    { 
	      int eN_i=eN*nDOF_test_element+i;
	      int gi = offset_u+stride_u*u_l2g[eN_i]; //global i-th index
	      // distribute global residual for (lumped) mass matrix
	      rhs_mass_correction[gi] += element_rhs_mass_correction[i];
	    }//i
	}//elements
      // COMPUTE LUMPED L2 PROJECTION
      for (int i=0; i<numDOFs; i++)
	{
	  double mi = lumped_mass_matrix[i];
	  lumped_L2p[i] = 1./mi*rhs_mass_correction[i];
	}
      return mass;
    }
    
    void calculateJacobian(//element
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
			   double epsFactHeaviside,
			   double epsFactDirac,
			   double epsFactDiffusion,
			   double cK,
			   double uL,
			   double uR,
			   double* phin_dof,
			   double* phiHat_dof)
    {
      //std::cout<<"ndjaco  address "<<q_numDiff_u_last<<std::endl;
      double Ct_sge = 4.0;
      //
      //loop over elements to compute volume integrals and load them into the element Jacobians and global Jacobian
      //
      for(int eN=0;eN<nElements_global;eN++)
	{
	  register double  elementJacobian_u_u[nDOF_test_element][nDOF_trial_element];
	  for (int i=0;i<nDOF_test_element;i++)
	    for (int j=0;j<nDOF_trial_element;j++)
	      {
		elementJacobian_u_u[i][j]=0.0;
	      }
	  for  (int k=0;k<nQuadraturePoints_element;k++)
	    {
	      int eN_k = eN*nQuadraturePoints_element+k, //index to a scalar at a quadrature point
		eN_k_nSpace = eN_k*nSpace,
		eN_nDOF_trial_element = eN*nDOF_trial_element; //index to a vector at a quadrature point

	      //declare local storage
	      register double u=0.0,
		grad_u[nSpace],
		m=0.0,dm=0.0,
		f[nSpace],df[nSpace],
		m_t=0.0,dm_t=0.0,
		dpdeResidual_u_u[nDOF_trial_element],
		Lstar_u[nDOF_test_element],
		dsubgridError_u_u[nDOF_trial_element],
		tau=0.0,tau0=0.0,tau1=0.0,
		jac[nSpace*nSpace],
		jacDet,
		jacInv[nSpace*nSpace],
		u_grad_trial[nDOF_trial_element*nSpace],
		dV,
		u_test_dV[nDOF_test_element],
		u_grad_test_dV[nDOF_test_element*nSpace],
		x,y,z,xt,yt,zt,
		//VRANS
		porosity,
		//
		G[nSpace*nSpace],G_dd_G,tr_G;
	      //
	      //calculate solution and gradients at quadrature points
	      //
	      // u=0.0;
	      // for (int I=0;I<nSpace;I++)
	      //   {
	      //     grad_u[I]=0.0;
	      //   }
	      // for (int j=0;j<nDOF_trial_element;j++)
	      //   {
	      //     int eN_j=eN*nDOF_trial_element+j;
	      //     int eN_k_j=eN_k*nDOF_trial_element+j;
	      //     int eN_k_j_nSpace = eN_k_j*nSpace;
              
	      //     u += valFromDOF_c(u_dof[u_l2g[eN_j]],u_trial[eN_k_j]);
	      //     for (int I=0;I<nSpace;I++)
	      //       {
	      //         grad_u[I] += gradFromDOF_c(u_dof[u_l2g[eN_j]],u_grad_trial[eN_k_j_nSpace+I]);
	      //       }
	      //   }
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
	      ck.calculateG(jacInv,G,G_dd_G,tr_G);
	      //get the trial function gradients
	      ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,u_grad_trial);
	      //get the solution 	
	      ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],u);
	      //get the solution gradients
	      ck.gradFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_u);
	      //precalculate test function products with integration weights
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
		  for (int I=0;I<nSpace;I++)
		    {
		      u_grad_test_dV[j*nSpace+I]   = u_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
		    }
		}
	      //VRANS
	      porosity = q_porosity[eN_k];
	      //
	      //
	      //calculate pde coefficients and derivatives at quadrature points
	      //
	      evaluateCoefficients(&velocity[eN_k_nSpace],
				   u,
				   //VRANS
				   porosity,
				   //
				   m,
				   dm,
				   f,
				   df);
	      //
	      //moving mesh
	      //
	      double mesh_velocity[3];
	      mesh_velocity[0] = xt;
	      mesh_velocity[1] = yt;
	      mesh_velocity[2] = zt;
	      //std::cout<<"qj mesh_velocity"<<std::endl;
	      for(int I=0;I<nSpace;I++)
		{
		  //std::cout<<mesh_velocity[I]<<std::endl;
		  f[I] -= MOVING_DOMAIN*m*mesh_velocity[I];
		  df[I] -= MOVING_DOMAIN*dm*mesh_velocity[I];
		}
	      //
	      //calculate time derivatives
	      //
	      ck.bdf(alphaBDF,
		     q_m_betaBDF[eN_k],//since m_t isn't used, we don't have to correct mass
		     m,
		     dm,
		     m_t,
		     dm_t);
	      //
	      //calculate subgrid error contribution to the Jacobian (strong residual, adjoint, jacobian of strong residual)
	      //
	      //calculate the adjoint times the test functions
	      for (int i=0;i<nDOF_test_element;i++)
		{
		  // int eN_k_i_nSpace = (eN_k*nDOF_trial_element+i)*nSpace;
		  // Lstar_u[i]=ck.Advection_adjoint(df,&u_grad_test_dV[eN_k_i_nSpace]);	      
		  register int i_nSpace = i*nSpace;
		  Lstar_u[i]=ck.Advection_adjoint(df,&u_grad_test_dV[i_nSpace]);	      
		}
	      //calculate the Jacobian of strong residual
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  //int eN_k_j=eN_k*nDOF_trial_element+j;
		  //int eN_k_j_nSpace = eN_k_j*nSpace;
		  int j_nSpace = j*nSpace;
		  dpdeResidual_u_u[j]= ck.MassJacobian_strong(dm_t,u_trial_ref[k*nDOF_trial_element+j]) +
		    ck.AdvectionJacobian_strong(df,&u_grad_trial[j_nSpace]);
		}
	      //tau and tau*Res
	      calculateSubgridError_tau(elementDiameter[eN],
					dm_t,
					df,
					cfl[eN_k],
					tau0);
  
              calculateSubgridError_tau(Ct_sge,
                                        G,
					dm_t,
					df,
					tau1,
				        cfl[eN_k]);
              tau = useMetrics*tau1+(1.0-useMetrics)*tau0;

	      for(int j=0;j<nDOF_trial_element;j++)
		dsubgridError_u_u[j] = -tau*dpdeResidual_u_u[j];
	      //double h=elementDiameter[eN];
	      for(int i=0;i<nDOF_test_element;i++)
		{
		  //int eN_k_i=eN_k*nDOF_test_element+i;
		  //int eN_k_i_nSpace=eN_k_i*nSpace;
		  for(int j=0;j<nDOF_trial_element;j++)
		    {
		      //int eN_k_j=eN_k*nDOF_trial_element+j;
		      //int eN_k_j_nSpace = eN_k_j*nSpace;
		      int j_nSpace = j*nSpace;
		      int i_nSpace = i*nSpace;
		      //std::cout<<"jac "<<'\t'<<q_numDiff_u_last[eN_k]<<'\t'<<dm_t<<'\t'<<df[0]<<df[1]<<'\t'<<dsubgridError_u_u[j]<<std::endl;
		      elementJacobian_u_u[i][j] += 
			ck.MassJacobian_weak(dm_t,u_trial_ref[k*nDOF_trial_element+j],u_test_dV[i]) + 
			ck.AdvectionJacobian_weak(df,u_trial_ref[k*nDOF_trial_element+j],&u_grad_test_dV[i_nSpace]) +
			ck.SubgridErrorJacobian(dsubgridError_u_u[j],Lstar_u[i]) +
			ck.NumericalDiffusionJacobian(q_numDiff_u_last[eN_k],&u_grad_trial[j_nSpace],&u_grad_test_dV[i_nSpace]); //implicit
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
		  globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_u_u[eN_i_j]] += elementJacobian_u_u[i][j];
		}//j
	    }//i
	}//elements
      //
      //loop over exterior element boundaries to compute the surface integrals and load them into the global Jacobian
      //
      for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) 
	{ 
	  register int ebN = exteriorElementBoundariesArray[ebNE]; 
	  register int eN  = elementBoundaryElementsArray[ebN*2+0],
            ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+0],
            eN_nDOF_trial_element = eN*nDOF_trial_element;
	  for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) 
	    { 
	      register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
		ebNE_kb_nSpace = ebNE_kb*nSpace,
		ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb,
		ebN_local_kb_nSpace = ebN_local_kb*nSpace;
	      register double u_ext=0.0,
		grad_u_ext[nSpace],
		m_ext=0.0,
		dm_ext=0.0,
		f_ext[nSpace],
		df_ext[nSpace],
		dflux_u_u_ext=0.0,
		bc_u_ext=0.0,
		//bc_grad_u_ext[nSpace],
		bc_m_ext=0.0,
		bc_dm_ext=0.0,
		bc_f_ext[nSpace],
		bc_df_ext[nSpace],
		fluxJacobian_u_u[nDOF_trial_element],
		jac_ext[nSpace*nSpace],
		jacDet_ext,
		jacInv_ext[nSpace*nSpace],
		boundaryJac[nSpace*(nSpace-1)],
		metricTensor[(nSpace-1)*(nSpace-1)],
		metricTensorDetSqrt,
		dS,
		u_test_dS[nDOF_test_element],
		u_grad_trial_trace[nDOF_trial_element*nSpace],
		normal[nSpace],x_ext,y_ext,z_ext,xt_ext,yt_ext,zt_ext,integralScaling,
		//VRANS
		porosity_ext,
		//
		G[nSpace*nSpace],G_dd_G,tr_G;
	      // 
	      //calculate the solution and gradients at quadrature points 
	      // 
	      // u_ext=0.0;
	      // for (int I=0;I<nSpace;I++)
	      //   {
	      //     grad_u_ext[I] = 0.0;
	      //     bc_grad_u_ext[I] = 0.0;
	      //   }
	      // for (int j=0;j<nDOF_trial_element;j++) 
	      //   { 
	      //     register int eN_j = eN*nDOF_trial_element+j,
	      //       ebNE_kb_j = ebNE_kb*nDOF_trial_element+j,
	      //       ebNE_kb_j_nSpace= ebNE_kb_j*nSpace;
	      //     u_ext += valFromDOF_c(u_dof[u_l2g[eN_j]],u_trial_ext[ebNE_kb_j]); 
	                     
	      //     for (int I=0;I<nSpace;I++)
	      //       {
	      //         grad_u_ext[I] += gradFromDOF_c(u_dof[u_l2g[eN_j]],u_grad_trial_ext[ebNE_kb_j_nSpace+I]); 
	      //       } 
	      //   }
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
	      //std::cout<<"J mtsqrdet "<<metricTensorDetSqrt<<" integralScaling "<<integralScaling<<std::endl;
	      dS = ((1.0-MOVING_DOMAIN)*metricTensorDetSqrt + MOVING_DOMAIN*integralScaling)*dS_ref[kb];
	      //dS = metricTensorDetSqrt*dS_ref[kb];
	      ck.calculateG(jacInv_ext,G,G_dd_G,tr_G);
	      //compute shape and solution information
	      //shape
	      ck.gradTrialFromRef(&u_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,u_grad_trial_trace);
	      //solution and gradients	
	      ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_trace_ref[ebN_local_kb*nDOF_test_element],u_ext);
	      ck.gradFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],u_grad_trial_trace,grad_u_ext);
	      //precalculate test function products with integration weights
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  u_test_dS[j] = u_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;
		}
	      //
	      //load the boundary values
	      //
	      bc_u_ext = isDOFBoundary_u[ebNE_kb]*ebqe_bc_u_ext[ebNE_kb]+(1-isDOFBoundary_u[ebNE_kb])*u_ext;
	      //VRANS
	      porosity_ext = ebqe_porosity_ext[ebNE_kb];
	      //
	      // 
	      //calculate the internal and external trace of the pde coefficients 
	      // 
	      evaluateCoefficients(&ebqe_velocity_ext[ebNE_kb_nSpace],
				   u_ext,
				   //VRANS
				   porosity_ext,
				   //
				   m_ext,
				   dm_ext,
				   f_ext,
				   df_ext);
	      evaluateCoefficients(&ebqe_velocity_ext[ebNE_kb_nSpace],
				   bc_u_ext,
				   //VRANS
				   porosity_ext,
				   //
				   bc_m_ext,
				   bc_dm_ext,
				   bc_f_ext,
				   bc_df_ext);
	      //
	      //moving domain
	      //
	      double mesh_velocity[3];
	      mesh_velocity[0] = xt_ext;
	      mesh_velocity[1] = yt_ext;
	      mesh_velocity[2] = zt_ext;
	      //std::cout<<"ext J mesh_velocity"<<std::endl;
	      for (int I=0;I<nSpace;I++)
		{
		  //std::cout<<mesh_velocity[I]<<std::endl;
		  f_ext[I] -= MOVING_DOMAIN*m_ext*mesh_velocity[I];
		  df_ext[I] -= MOVING_DOMAIN*dm_ext*mesh_velocity[I];
		  bc_f_ext[I] -= MOVING_DOMAIN*bc_m_ext*mesh_velocity[I];
		  bc_df_ext[I] -= MOVING_DOMAIN*bc_dm_ext*mesh_velocity[I];
		}
	      // 
	      //calculate the numerical fluxes 
	      // 
	      exteriorNumericalAdvectiveFluxDerivative(isDOFBoundary_u[ebNE_kb],
						       isFluxBoundary_u[ebNE_kb],
						       normal,
						       df_ext,//VRANS holds porosity
						       dflux_u_u_ext);
	      //
	      //calculate the flux jacobian
	      //
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  //register int ebNE_kb_j = ebNE_kb*nDOF_trial_element+j;
		  register int ebN_local_kb_j=ebN_local_kb*nDOF_trial_element+j;
		  fluxJacobian_u_u[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_u_u_ext,u_trial_trace_ref[ebN_local_kb_j]);
		}//j
	      //
	      //update the global Jacobian from the flux Jacobian
	      //
	      for (int i=0;i<nDOF_test_element;i++)
		{
		  register int eN_i = eN*nDOF_test_element+i;
		  //register int ebNE_kb_i = ebNE_kb*nDOF_test_element+i;
		  for (int j=0;j<nDOF_trial_element;j++)
		    {
		      register int ebN_i_j = ebN*4*nDOF_test_X_trial_element + i*nDOF_trial_element + j;
		      globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] += fluxJacobian_u_u[j]*u_test_dS[i];
		    }//j
		}//i
	    }//kb
	}//ebNE
    }//computeJacobian

    void calculateMassMatrix(//element
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
			     double epsFactHeaviside,
			     double epsFactDirac,
			     double epsFactDiffusion,
			     double cK,
			     double uL,
			     double uR,
			     double* phin_dof,
			     double* phiHat_dof)
    {
      //std::cout<<"ndjaco  address "<<q_numDiff_u_last<<std::endl;
      double Ct_sge = 4.0;
      //
      //loop over elements to compute volume integrals and load them into the element Jacobians and global Jacobian
      //
      for(int eN=0;eN<nElements_global;eN++)
	{
	  register double  elementJacobian_u_u[nDOF_test_element][nDOF_trial_element];
	  for (int i=0;i<nDOF_test_element;i++)
	    for (int j=0;j<nDOF_trial_element;j++)
	      {
		elementJacobian_u_u[i][j]=0.0;
	      }
	  for  (int k=0;k<nQuadraturePoints_element;k++)
	    {
	      int eN_k = eN*nQuadraturePoints_element+k, //index to a scalar at a quadrature point
		eN_k_nSpace = eN_k*nSpace,
		eN_nDOF_trial_element = eN*nDOF_trial_element; //index to a vector at a quadrature point

	      //declare local storage
	      register double u=0.0,
		grad_u[nSpace],
		m=0.0,dm=0.0,
		f[nSpace],df[nSpace],
		m_t=0.0,dm_t=0.0,
		dpdeResidual_u_u[nDOF_trial_element],
		Lstar_u[nDOF_test_element],
		dsubgridError_u_u[nDOF_trial_element],
		tau=0.0,tau0=0.0,tau1=0.0,
		jac[nSpace*nSpace],
		jacDet,
		jacInv[nSpace*nSpace],
		u_grad_trial[nDOF_trial_element*nSpace],
		dV,
		u_test_dV[nDOF_test_element],
		u_grad_test_dV[nDOF_test_element*nSpace],
		x,y,z,xt,yt,zt,
		//VRANS
		porosity,
		//
		G[nSpace*nSpace],G_dd_G,tr_G;

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
	      ck.calculateG(jacInv,G,G_dd_G,tr_G);
	      //get the trial function gradients
	      ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,u_grad_trial);
	      //get the solution 	
	      ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],u);
	      //get the solution gradients
	      ck.gradFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_u);
	      //precalculate test function products with integration weights
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
		  for (int I=0;I<nSpace;I++)
		    {
		      u_grad_test_dV[j*nSpace+I]   = u_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
		    }
		}
	      //VRANS
	      porosity = q_porosity[eN_k];
	      //
	      //
	      //calculate pde coefficients and derivatives at quadrature points
	      //
	      evaluateCoefficients(&velocity[eN_k_nSpace],
				   u,
				   //VRANS
				   porosity,
				   //
				   m,
				   dm,
				   f,
				   df);
	      //
	      //moving mesh
	      //
	      double mesh_velocity[3];
	      mesh_velocity[0] = xt;
	      mesh_velocity[1] = yt;
	      mesh_velocity[2] = zt;
	      //std::cout<<"qj mesh_velocity"<<std::endl;
	      for(int I=0;I<nSpace;I++)
		{
		  //std::cout<<mesh_velocity[I]<<std::endl;
		  f[I] -= MOVING_DOMAIN*m*mesh_velocity[I];
		  df[I] -= MOVING_DOMAIN*dm*mesh_velocity[I];
		}
	      //
	      //calculate time derivatives
	      //
	      ck.bdf(alphaBDF,
		     q_m_betaBDF[eN_k],//since m_t isn't used, we don't have to correct mass
		     m,
		     dm,
		     m_t,
		     dm_t);
	      //
	      //calculate subgrid error contribution to the Jacobian (strong residual, adjoint, jacobian of strong residual)
	      //
	      //calculate the adjoint times the test functions
	      for (int i=0;i<nDOF_test_element;i++)
		{
		  // int eN_k_i_nSpace = (eN_k*nDOF_trial_element+i)*nSpace;
		  // Lstar_u[i]=ck.Advection_adjoint(df,&u_grad_test_dV[eN_k_i_nSpace]);	      
		  register int i_nSpace = i*nSpace;
		  Lstar_u[i]=ck.Advection_adjoint(df,&u_grad_test_dV[i_nSpace]);	      
		}
	      //calculate the Jacobian of strong residual
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  //int eN_k_j=eN_k*nDOF_trial_element+j;
		  //int eN_k_j_nSpace = eN_k_j*nSpace;
		  int j_nSpace = j*nSpace;
		  dpdeResidual_u_u[j]= ck.MassJacobian_strong(dm_t,u_trial_ref[k*nDOF_trial_element+j]) +
		    ck.AdvectionJacobian_strong(df,&u_grad_trial[j_nSpace]);
		}
	      //tau and tau*Res
	      calculateSubgridError_tau(elementDiameter[eN],
					dm_t,
					df,
					cfl[eN_k],
					tau0);
  
              calculateSubgridError_tau(Ct_sge,
                                        G,
					dm_t,
					df,
					tau1,
				        cfl[eN_k]);
              tau = useMetrics*tau1+(1.0-useMetrics)*tau0;

	      for(int j=0;j<nDOF_trial_element;j++)
		dsubgridError_u_u[j] = -tau*dpdeResidual_u_u[j];
	      //double h=elementDiameter[eN];
	      for(int i=0;i<nDOF_test_element;i++)
		{
		  //int eN_k_i=eN_k*nDOF_test_element+i;
		  //int eN_k_i_nSpace=eN_k_i*nSpace;
		  for(int j=0;j<nDOF_trial_element;j++)
		    {
		      if (LUMPED_MASS_MATRIX==1)
			{
			  if (i==j)
			    elementJacobian_u_u[i][j] += u_test_dV[i];
			}
		      else
			{
			  //int eN_k_j=eN_k*nDOF_trial_element+j;
			  //int eN_k_j_nSpace = eN_k_j*nSpace;
			  int j_nSpace = j*nSpace;
			  int i_nSpace = i*nSpace;
			  //std::cout<<"jac "<<'\t'<<q_numDiff_u_last[eN_k]<<'\t'<<dm_t<<'\t'<<df[0]<<df[1]<<'\t'<<dsubgridError_u_u[j]<<std::endl;
			  elementJacobian_u_u[i][j] += 
			    dt*ck.MassJacobian_weak(dm_t,u_trial_ref[k*nDOF_trial_element+j],u_test_dV[i]);
			}
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
		  globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_u_u[eN_i_j]] += elementJacobian_u_u[i][j];
		}//j
	    }//i
	}//elements
    }//computeJacobian

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
			     double epsFactHeaviside,
			     double epsFactDirac,
			     double epsFactDiffusion,
			     double cK,
			     double uL,
			     double uR,
			     double* phin_dof,
			     double* phiHat_dof)
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
			     double epsFactHeaviside,
			     double epsFactDirac,
			     double epsFactDiffusion,
			     double cK,
			     double uL,
			     double uR,
			     double* phin_dof,
			     double* phiHat_dof)
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
	      register double u, u_grad_trial[nDOF_trial_element*nSpace], grad_u[nSpace],
		phin, phiHatnp1, 
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
	      ck.gradFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_u);
	      //get the solution 	
	      ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],u);
	      //get phiHat at tnp1
	      ck.valFromDOF(phiHat_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],phiHatnp1);
	      ck.valFromDOF(phin_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],phin);
	      //precalculate test function products with integration weights
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
		  for (int I=0;I<nSpace;I++)
		    u_grad_test_dV[j*nSpace+I]   = u_grad_trial[j*nSpace+I]*dV;
		}
	      double hK = elementDiameter[eN];
	      double epsHeaviside = epsFactHeaviside*hK;
	      double epsDiffusion = epsFactDiffusion*hK;
	      double epsDirac = epsFactDirac*hK;
	      double Hn = smoothedHeaviside(epsHeaviside,phin);
	      double time_derivative_jacobian = smoothedDirac(epsDirac,phiHatnp1+u)/dt;
	      
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

	      double norm_grad_Hn=smoothedDirac(epsDirac,Hn)*std::sqrt(grad_u[0]*grad_u[0]
								       +grad_u[1]*grad_u[1]);
	      double viscosity = 1.0;
	      double diffusion = epsFactDiffusion;
	      //double diffusion = viscosity*fmax(1-cK*(Hn-uL)*(uR-Hn)/(hK*norm_grad_Hn+1E-10),0.0);
	      //double diffusion = viscosity*(1-cK*(Hn-uL)*(uR-Hn))/(hK*norm_grad_Hn+1E-10);
	      
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
			+ ck.NumericalDiffusionJacobian(epsDiffusion/dt,
							&u_grad_trial[j_nSpace],
							&u_grad_test_dV[i_nSpace]);
			//+ ck.NumericalDiffusionJacobian(diffusion,
			//				&u_grad_trial[j_nSpace],
			//				&u_grad_test_dV[i_nSpace]); 
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
    
  };//VOF

  inline VOF_base* newVOF(int nSpaceIn,
				int nQuadraturePoints_elementIn,
				int nDOF_mesh_trial_elementIn,
				int nDOF_trial_elementIn,
				int nDOF_test_elementIn,
				int nQuadraturePoints_elementBoundaryIn,
				int CompKernelFlag)
  {
    if (nSpaceIn == 2)
      return proteus::chooseAndAllocateDiscretization2D<VOF_base,VOF,CompKernel>(nSpaceIn,
										 nQuadraturePoints_elementIn,
										 nDOF_mesh_trial_elementIn,
										 nDOF_trial_elementIn,
										 nDOF_test_elementIn,
										 nQuadraturePoints_elementBoundaryIn,
										 CompKernelFlag);
    else
      return proteus::chooseAndAllocateDiscretization<VOF_base,VOF,CompKernel>(nSpaceIn,
									       nQuadraturePoints_elementIn,
									       nDOF_mesh_trial_elementIn,
									       nDOF_trial_elementIn,
									       nDOF_test_elementIn,
									       nQuadraturePoints_elementBoundaryIn,
									       CompKernelFlag);
  }
}//proteus
#endif

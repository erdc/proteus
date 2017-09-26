#ifndef VOF_H
#define VOF_H
#include <cmath>
#include <iostream>
#include "CompKernel.h"
#include "ModelFactory.h"

#define KUZMINS_METHOD 1
#define INTEGRATE_BY_PARTS 1
#define QUANTITIES_OF_INTEREST 0
#define FIX_BOUNDARY_KUZMINS 1
#define POWER_SMOOTHNESS_INDICATOR 2
#define BETAij 1
#define GLOBAL_FCT 0

/////////////////////
//ENTROPY FUNCTION //
/////////////////////
// Power entropy //
//#define entropy_power 2. // phiL and phiR are dummy variables
//#define ENTROPY(phi,phiL,phiR) 1./entropy_power*std::pow(std::abs(phi),entropy_power)
//#define DENTROPY(phi,phiL,phiR) std::pow(std::abs(phi),entropy_power-1.)*(phi>=0 ? 1 : -1)
//#define ENTROPY_GRAD(phi,phix,phiL,phiR) std::pow(std::abs(phi),entropy_power-1.)*(phi>=0 ? 1 : -1)*phix
// Log entropy //
// LOG ENTROPY FOR LEVEL SET FROM 0 to 1
#define ENTROPY(phi,phiL,phiR) std::log(std::abs((phi-phiL)*(phiR-phi))+1E-14)
#define DENTROPY(phi,phiL,phiR) (phiL+phiR-2*phi)*((phi-phiL)*(phiR-phi)>=0 ? 1 : -1)/(std::abs((phi-phiL)*(phiR-phi))+1E-14) 
#define ENTROPY_GRAD(phi,phix,phiL,phiR) (phiL+phiR-2*phi)*phix*((phi-phiL)*(phiR-phi)>=0 ? 1 : -1)/(std::abs((phi-phiL)*(phiR-phi))+1E-14) 

namespace proteus
{
  class VOF_base
  {
    //The base class defining the interface
  public:
    virtual ~VOF_base(){}
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
			 double* dt_times_dC_minus_dL, //low minus high order dissipative matrices
			 double* min_u_bc, //min/max value at BCs. If DOF is not at boundary then min=1E10, max=-1E10
			 double* max_u_bc,
			 int LUMPED_MASS_MATRIX
			 )=0;
    virtual void calculateResidual(//element
				   double* mesh_trial_ref,
				   double* mesh_grad_trial_ref,
				   double* mesh_dof,
				   double* meshVelocity_dof,
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
				   //
				   int* u_l2g, 
				   double* elementDiameter,
				   int degree_polynomial,
				   double* u_dof,
				   double* u_dof_lstage,
				   double* u_dof_old,
				   double* u_dof_old_old,
				   double* velocity,
				   double* velx_dof,
				   double* vely_dof,
				   double* velz_dof,
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
				   int EDGE_VISCOSITY, 
				   int ENTROPY_VISCOSITY, 
				   // PARAMETERS FOR ENTROPY_VISCOSITY 
				   double cE,
				   double cMax, 
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
				   double* CTx,
				   double* CTy,
				   // PARAMETERS FOR 1st or 2nd ORDER MPP METHOD
				   int LUMPED_MASS_MATRIX, 
 				   // FOR FCT
				   double* low_order_solution,
				   double* dt_times_dC_minus_dL,
				   double* min_u_bc,
				   double* max_u_bc,
				   // AUX QUANTITIES OF INTEREST
				   double* quantDOFs)=0;
    virtual void calculateResidual_entropy_viscosity(//element
				   double* mesh_trial_ref,
				   double* mesh_grad_trial_ref,
				   double* mesh_dof,
				   double* meshVelocity_dof,
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
				   //
				   int* u_l2g, 
				   double* elementDiameter,
				   int degree_polynomial,
				   double* u_dof,
				   double* u_dof_lstage,
				   double* u_dof_old,
				   double* u_dof_old_old,
				   double* velocity,
				   double* velx_dof,
				   double* vely_dof,
				   double* velz_dof,			     
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
				   int EDGE_VISCOSITY, 
				   int ENTROPY_VISCOSITY, 
				   // PARAMETERS FOR ENTROPY_VISCOSITY 
				   double cE,
				   double cMax, 
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
				   double* CTx,
				   double* CTy,
				   // PARAMETERS FOR 1st or 2nd ORDER MPP METHOD
				   int LUMPED_MASS_MATRIX, 
 				   // FOR FCT
				   double* low_order_solution,
				   double* dt_times_dC_minus_dL,
				   double* min_u_bc,
				   double* max_u_bc,
				   // AUX QUANTITIES OF INTEREST
				   double* quantDOFs)=0;
    virtual void calculateJacobian(//element
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
				   // PARAMETERS FOR EDGE_VISCOSITY
				   int EDGE_VISCOSITY, 
				   int ENTROPY_VISCOSITY,
				   int LUMPED_MASS_MATRIX)=0;
    virtual void calculateMassMatrix(//element
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
				   // PARAMETERS FOR EDGE_VISCOSITY
				   int EDGE_VISCOSITY, 
				   int ENTROPY_VISCOSITY,
				   int LUMPED_MASS_MATRIX)=0;
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
		 double* dt_times_dC_minus_dL, //low minus high order dissipative matrices
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
		  mini = std::min(mini,soln[j]);
		  maxi = std::max(maxi,soln[j]);
		}      
	      // i-th row of flux correction matrix 
	      double ML_minus_MC = (LUMPED_MASS_MATRIX == 1 ? 0. : (i==j ? 1. : 0.)*mi - MassMatrix[ij]);
	      FluxCorrectionMatrix[ij] = ML_minus_MC * (solH[j]-soln[j] - (solHi-solni)) 
		+ dt_times_dC_minus_dL[ij]*(soln[j]-solni);

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
	  limited_solution[i] = solL[i] + 1./lumped_mass_matrix[i]*ith_Limiter_times_FluxCorrectionMatrix;	  
	}
    }

    void calculateResidual(//element
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
			   //
			   int* u_l2g, 
			   double* elementDiameter,
			   int degree_polynomial,
			   double* u_dof,
			   double* u_dof_lstage,
			   double* u_dof_old,
			   double* u_dof_old_old,
			   double* velocity,
			   double* velx_dof,
			   double* vely_dof,
			   double* velz_dof,				   
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
			   int EDGE_VISCOSITY, 
			   int ENTROPY_VISCOSITY, 
			   double cE,
			   double cMax, 
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
			   double* CTx,
			   double* CTy, 
			   // PARAMETERS FOR 1st or 2nd ORDER MPP METHOD
			   int LUMPED_MASS_MATRIX, 
			   // FOR FCT
			   double* low_order_solution,
			   double* dt_times_dC_minus_dL,
			   double* min_u_bc,
			   double* max_u_bc,
			   // AUX QUANTITIES OF INTEREST 
			   double* quantDOFs)
    {
      double dt = 1./alphaBDF; // HACKED to work just for BDF1
      if (EDGE_VISCOSITY==1)
	{
	  // Allocate space for the transport matrices
	  // This is used for first order KUZMIN'S METHOD
	  register double TransportMatrix[NNZ], TransposeTransportMatrix[NNZ];
	  if (KUZMINS_METHOD==1)
	    for (int i=0; i<NNZ; i++)
	      {
		TransportMatrix[i] = 0.;
		TransposeTransportMatrix[i] = 0.;
	      }

	  // Allocate and init to zero the Entropy residual vector
	  register double EntResVector[numDOFs], entropy_residual;
	  if (ENTROPY_VISCOSITY==1)
	    for (int i=0; i<numDOFs; i++)
	      EntResVector[i]=0.;	    
	  
	  //////////////////////////////////////////////
	  // ** LOOP IN CELLS FOR CELL BASED TERMS ** //
	  //////////////////////////////////////////////
	  // HERE WE COMPUTE: 
	  //    * CFL for VOF model
	  //    * Entropy residual (if ENTROPY_VISCOSITY=1)
	  //    * (Lumped) Mass matrix times (soln - soln_tn) 
	  //    * Transport matrices (if KUZMINS=1). 
	  for(int eN=0;eN<nElements_global;eN++)
	    {
	      //declare local storage for local contributions and initialize
	      register double elementResidual_u[nDOF_test_element], elementEntResVector[nDOF_test_element];
	      register double  elementTransport[nDOF_test_element][nDOF_trial_element];
	      register double  elementTransposeTransport[nDOF_test_element][nDOF_trial_element];
	      for (int i=0;i<nDOF_test_element;i++)
		{
		  elementResidual_u[i]=0.0;
		  elementEntResVector[i]=0.0;
		  for (int j=0;j<nDOF_trial_element;j++)
		    {
		      elementTransport[i][j]=0.0;
		      elementTransposeTransport[i][j]=0.0;
		    }
		}
	      //restart cell based quantities 
	      entropy_residual = 0;
	      //loop over quadrature points and compute integrands
	      for  (int k=0;k<nQuadraturePoints_element;k++)
		{
		  //compute indeces and declare local storage
		  register int eN_k = eN*nQuadraturePoints_element+k,
		    eN_k_nSpace = eN_k*nSpace,
		    eN_nDOF_trial_element = eN*nDOF_trial_element;
		  register double 
		    //for mass matrix contributions
		    u=0.0, m=0.0, dm=0.0, f[nSpace], df[nSpace], m_t=0.0,dm_t=0.0, 
		    u_test_dV[nDOF_trial_element], 
		    //for entropy viscosity
		    un=0.0, unm1=0.0, grad_un[nSpace], vn[nSpace], 
		    u_grad_trial[nDOF_trial_element*nSpace], 
		    u_grad_test_dV[nDOF_test_element*nSpace],
		    //for general use
		    jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace],
		    dV,x,y,z,
		    //VRANS
		    porosity;
		  //get the physical integration weight
		  ck.calculateMapping_element(eN,k,mesh_dof,mesh_l2g,mesh_trial_ref,mesh_grad_trial_ref,jac,jacDet,jacInv,x,y,z);
		  dV = fabs(jacDet)*dV_ref[k];
		  //get the solution (of Newton's solver). To compute time derivative term
		  ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],u);
		  //get the solution at quad point at tn and tnm1 for entropy viscosity
		  ck.valFromDOF(u_dof_old,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],un);
		  ck.valFromDOF(u_dof_old_old,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],unm1);
		  //get the solution gradients at tn for entropy viscosity
		  ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,u_grad_trial);
		  ck.gradFromDOF(u_dof_old,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_un);

		  //precalculate test function products with integration weights for mass matrix terms
		  for (int j=0;j<nDOF_trial_element;j++)
		    {
		      u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
		      for (int I=0;I<nSpace;I++)
			{
			  u_grad_test_dV[j*nSpace+I] = u_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
			}
		    }
		  //evaluate coefficients to compute time derivative (for term with mass matrix)
		  porosity = q_porosity[eN_k];
		  evaluateCoefficients(&velocity[eN_k_nSpace],
				       u,
				       //VRANS
				       porosity,
				       //
				       m,
				       dm,
				       f,
				       df);
		  //calculate time derivative at quadrature points
		  if (q_dV_last[eN_k] <= -100)
		    q_dV_last[eN_k] = dV;
		  q_dV[eN_k] = dV;
		  ck.bdf(alphaBDF,
			 q_m_betaBDF[eN_k]*q_dV_last[eN_k]/dV,//ensure prior mass integral is correct for  m_t with BDF1
			 m,
			 dm,
			 m_t,
			 dm_t);
		  // CALCULATE CFL //
		  calculateCFL(elementDiameter[eN]/degree_polynomial,df,cfl[eN_k]); // TODO: ADJUST SPEED IF MESH IS MOVING
		  
		  // CALCULATE ENTROPY RESIDUAL AT QUAD POINT //
		  //velocity at tn for entropy viscosity
		  if (ENTROPY_VISCOSITY==1)
		    {
		      vn[0] = velocity[eN_k_nSpace];
		      vn[1] = velocity[eN_k_nSpace+1];
		      entropy_residual = 
			(ENTROPY(un,uL,uR) - ENTROPY(unm1,uL,uR))/dt // time derivative
			+ vn[0]*ENTROPY_GRAD(un,grad_un[0],uL,uR)+vn[1]*ENTROPY_GRAD(un,grad_un[1],uL,uR) // velocity * grad(entropy)
			+ ENTROPY(un,uL,uR)*(vn[0]+vn[1]); // For non div free velocities
		    }

		  //////////////
		  // ith-LOOP //
		  //////////////
		  for(int i=0;i<nDOF_test_element;i++) 
		    { 
		      // MASS MATRIX TERM //
		      if (LUMPED_MASS_MATRIX==1)
			elementResidual_u[i] += u_test_dV[i]; // LUMPING. We multiply times the DOFs during distribution
		      else 
			elementResidual_u[i] += dt*ck.Mass_weak(m_t,u_test_dV[i]);

		      // VECTOR OF ENTROPY RESIDUAL //
		      if (ENTROPY_VISCOSITY==1)
			elementEntResVector[i] += entropy_residual*u_test_dV[i];

		      ///////////////
		      // j-th LOOP // To construct transport matrices
		      ///////////////
		      if (KUZMINS_METHOD==1)
			for(int j=0;j<nDOF_trial_element;j++) 
			  { 
			    int j_nSpace = j*nSpace;
			    int i_nSpace = i*nSpace;
			    // COMPUTE ELEMENT TRANSPORT MATRIX (MQL)
			    if (INTEGRATE_BY_PARTS==1)
			      {
				elementTransport[i][j] += // -int[(vel.grad_wi)*wj*dx]
				  ck.AdvectionJacobian_weak(df,u_trial_ref[k*nDOF_trial_element+j],&u_grad_test_dV[i_nSpace]);
				elementTransposeTransport[i][j] += // -int[(vel.grad_wj)*wi*dx]
				  ck.AdvectionJacobian_weak(df,u_trial_ref[k*nDOF_trial_element+i],&u_grad_test_dV[j_nSpace]);
			      }
			    else
			      {
				// Transport matrix
				//elementTransport[i][j] += // int[(vel.grad_wj + div(vel)*wj)*wi*dx]
				//-ck.AdvectionJacobian_weak(df,u_trial_ref[k*nDOF_trial_element+i],&u_grad_test_dV[j_nSpace]) 
				//+ div_velocity[eN_k]*u_trial_ref[k*nDOF_trial_element+i]*u_trial_ref[k*nDOF_trial_element+j];
				// Transpose of Transport matrix
				//elementTransposeTransport[i][j] += // int[(vel.grad_wi + div(v)*wi)*wj*dx]
				//-ck.AdvectionJacobian_weak(df,u_trial_ref[k*nDOF_trial_element+j],&u_grad_test_dV[i_nSpace])
				//+ div_velocity[eN_k]*u_trial_ref[k*nDOF_trial_element+i]*u_trial_ref[k*nDOF_trial_element+j];
			      }
			  }
		    }//i
		  //save solution for other models 
		  q_u[eN_k] = u;
		  q_m[eN_k] = m;
		}
	      /////////////////
	      // DISTRIBUTE // load cell based element into global residual
	      ////////////////
	      for(int i=0;i<nDOF_test_element;i++) 
		{ 
		  int eN_i=eN*nDOF_test_element+i;
		  int gi = offset_u+stride_u*u_l2g[eN_i]; //global i-th index
		  // distribute global residual for (lumped) mass matrix 

		  if (LUMPED_MASS_MATRIX==1)
		    globalResidual[gi] += elementResidual_u[i]*(u_dof[gi]-u_dof_lstage[gi]); //LUMPING
		  else
		    globalResidual[gi] += elementResidual_u[i];
		  // distribute EntResVector 
		  if (ENTROPY_VISCOSITY==1)
		    EntResVector[gi] += elementEntResVector[i];
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

	  //////////////////////////////////////////////////////
	  // ADD OUTFLOW BOUNDARY TERM TO TRANSPORT MATRICES //
	  /////////////////////////////////////////////////////
	  //   * Compute outflow boundary integral as a matrix; i.e., int_B[ (vel.normal)*wi*wj*dx]
	  if (KUZMINS_METHOD==1 && INTEGRATE_BY_PARTS==1 && FIX_BOUNDARY_KUZMINS==1)
	    {
	      for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) 
		{ 
		  double min_u_bc_local = 1E10, max_u_bc_local = -1E10;
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
			  
		      register double u_ext=0.0,m_ext=0.0,dm_ext=0.0,f_ext[nSpace],df_ext[nSpace],
			dflux_ext=0.0,
			fluxTransport[nDOF_trial_element],
			jac_ext[nSpace*nSpace],
			jacDet_ext,
			jacInv_ext[nSpace*nSpace],
			boundaryJac[nSpace*(nSpace-1)],
			metricTensor[(nSpace-1)*(nSpace-1)],
			metricTensorDetSqrt,
			dS,
			u_test_dS[nDOF_test_element],
			normal[nSpace],x_ext,y_ext,z_ext,xt_ext,yt_ext,zt_ext,integralScaling,
			//VRANS
			porosity_ext;
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
		      //std::cout<<"J mtsqrdet "<<metricTensorDetSqrt<<" integralScaling "<<integralScaling<<std::endl;
		      dS = ((1.0-MOVING_DOMAIN)*metricTensorDetSqrt + MOVING_DOMAIN*integralScaling)*dS_ref[kb];
		      //compute shape and solution information
		      ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_trace_ref[ebN_local_kb*nDOF_test_element],u_ext);
		      //precalculate test function products with integration weights
		      for (int j=0;j<nDOF_trial_element;j++)
			u_test_dS[j] = u_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;
		      //VRANS
		      porosity_ext = ebqe_porosity_ext[ebNE_kb];
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
		      //
		      //moving domain
		      //
		      double mesh_velocity[3];
		      mesh_velocity[0] = xt_ext;
		      mesh_velocity[1] = yt_ext;
		      mesh_velocity[2] = zt_ext;
		      for (int I=0;I<nSpace;I++)
			df_ext[I] -= MOVING_DOMAIN*dm_ext*mesh_velocity[I];
		      // 
		      //calculate the numerical fluxes 
		      // 
		      exteriorNumericalAdvectiveFluxDerivative(isDOFBoundary_u[ebNE_kb],
							       isFluxBoundary_u[ebNE_kb],
							       normal,
							       df_ext,//VRANS holds porosity
							       dflux_ext);
		      //
		      //calculate the flux jacobian
		      //
		      double flow = 0.;
		      for (int I=0; I < nSpace; I++)
			flow += normal[I]*ebqe_velocity_ext[ebNE_kb_nSpace+I];
		      if (flow >= 0)
			dflux_ext = flow;
		      else 
			dflux_ext = 0;

		      for (int j=0;j<nDOF_trial_element;j++)
			{
			  //register int ebNE_kb_j = ebNE_kb*nDOF_trial_element+j;
			  register int ebN_local_kb_j=ebN_local_kb*nDOF_trial_element+j;
			  fluxTransport[j]=
			    ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_ext,u_trial_trace_ref[ebN_local_kb_j]);
			}//j
		      //
		      //update the global Transport Matrices
		      //
		      for (int i=0;i<nDOF_test_element;i++)
			{
			  register int eN_i = eN*nDOF_test_element+i;
			  //register int ebNE_kb_i = ebNE_kb*nDOF_test_element+i;
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
		      min_u_bc_local = std::min(ebqe_bc_u_ext[ebNE_kb], min_u_bc_local);
		      max_u_bc_local = std::max(ebqe_bc_u_ext[ebNE_kb], max_u_bc_local);
		    }//kb
		  // global min/max at boundary 
		  for (int i=0;i<nDOF_test_element;i++)
		    {
		      int eN_i = eN*nDOF_test_element+i;
		      int gi = offset_u+stride_u*u_l2g[eN_i]; //global i-th index
		      min_u_bc[gi] = std::min(min_u_bc_local,min_u_bc[gi]);
		      max_u_bc[gi] = std::max(max_u_bc_local,max_u_bc[gi]);		      
		    }
		}//ebNE
	      // END OF ADDING BOUNDARY TERM TO TRANSPORT MATRICES //
	    }

	  /////////////////////////////////////////////////////////////
	  // MULTIPLICATION OF TRANSPORT MATRIX TIMES SOLUTION AT tn //
	  /////////////////////////////////////////////////////////////
	  int ij=0;
	  if (KUZMINS_METHOD==1)
	    {	  
	      for (int i=0; i<numDOFs; i++)
		{		  
		  double ith_TransportMatrix_times_solution = 0.;
		  // LOOP OVER THE SPARSITY PATTERN (j-LOOP)//
		  for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
		    {
		      int j = csrColumnOffsets_DofLoops[offset];
		      ith_TransportMatrix_times_solution += TransportMatrix[ij]*u_dof_lstage[j];
		      //update ij
		      ij+=1;
		    }
		  globalResidual[i] += dt*ith_TransportMatrix_times_solution;
		}
	    }

	  ////////////////////////////////////////
	  // COMPUTE (INFLOW) BOUNDARY INTEGRAL //
	  ////////////////////////////////////////
	  if (KUZMINS_METHOD==1 && INTEGRATE_BY_PARTS==1 && FIX_BOUNDARY_KUZMINS==1)
	    {
	      //ebNE is the Exterior element boundary INdex
	      //ebN is the element boundary INdex
	      //eN is the element index
	      for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) 
		{ 
		  register int ebN = exteriorElementBoundariesArray[ebNE], 
		    eN  = elementBoundaryElementsArray[ebN*2+0],
		    ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+0],
		    eN_nDOF_trial_element = eN*nDOF_trial_element;
		  register double elementResidual_u[nDOF_test_element], cell_flow;
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
		      register double
			u_ext=0.0,
			flux_ext=0.0,
			bc_u_ext=0.0,
			jac_ext[nSpace*nSpace],
			jacDet_ext,
			jacInv_ext[nSpace*nSpace],
			boundaryJac[nSpace*(nSpace-1)],
			metricTensor[(nSpace-1)*(nSpace-1)],
			metricTensorDetSqrt,
			dS,
			u_test_dS[nDOF_test_element],
			u_grad_trial_trace[nDOF_trial_element*nSpace],
			normal[nSpace],x_ext,y_ext,z_ext;
		      cell_flow = 0.0;
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
		      dS = metricTensorDetSqrt*dS_ref[kb];
		      //precalculate test function products with integration weights
		      // compute flow = v.n
		      double flow=0.0;
		      for (int I=0; I < nSpace; I++)
			flow += normal[I]*ebqe_velocity_ext[ebNE_kb_nSpace+I];		      		      
		      for (int j=0;j<nDOF_trial_element;j++)
			u_test_dS[j] = u_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;
		      //load the boundary values
		      bc_u_ext = ebqe_bc_u_ext[ebNE_kb];
		      // compute flux = flow*bc_u
		      if (flow < 0)
			flux_ext = bc_u_ext*flow;
		      else 
			flux_ext = 0;
		      
		      ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_trace_ref[ebN_local_kb*nDOF_test_element],u_ext);
		      ebqe_flux[ebNE_kb] = flux_ext;
		      //save for other models? cek need to be consistent with numerical flux
		      if(flux_ext >=0.0)
			ebqe_u[ebNE_kb] = u_ext;
		      else
			ebqe_u[ebNE_kb] = bc_u_ext;

		      for (int i=0;i<nDOF_test_element;i++)
			elementResidual_u[i] += ck.ExteriorElementBoundaryFlux(flux_ext,u_test_dS[i]); 		    
		    }//kb
		  // Distribute
		  for (int i=0;i<nDOF_test_element;i++)
		    {
		      int eN_i = eN*nDOF_test_element+i;
		      int gi = offset_u+stride_u*u_l2g[eN_i];
		      globalResidual[gi] += dt*elementResidual_u[i];
		    }//i
		}//ebNE	  
	      // END OF BOUNDARY INTEGRAL //
	    }

	  //////////////////////////////////////////////////////
	  // COMPUTE (FULL: INFLOW+OUTFLOW) BOUNDARY INTEGRAL //
	  //////////////////////////////////////////////////////
	  if (KUZMINS_METHOD==1 && INTEGRATE_BY_PARTS==1 && FIX_BOUNDARY_KUZMINS==0)
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
			elementResidual_u[i] += dt*ck.ExteriorElementBoundaryFlux(flux_ext,u_test_dS[i]);
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
	  // END OF FULL BOUNDARY INTEGRAL

	  //////////////////////////////////
	  // COMPUTE SMOOTHNESS INDICATOR // and // COMPUTE LOCAL MAX OF ENT RESIDUALS //
	  //////////////////////////////////
	  // Smoothness indicator is based on the solution. psi_i = psi_i(alpha_i); alpha_i = |sum(betaij*(uj-ui))|/sum(betaij*|uj-ui|)
	  // MaxEntResi = max_{j\in Si} (|EntRes_j|)
	  register double MaxEntResVector[numDOFs], psi[numDOFs];
	  for (int i=0; i<numDOFs; i++)
	    {
	      double MaxEntResi=0;
	      double alphai, alphai_numerator=0, alphai_denominator=0; // smoothness indicator of solution
	      int num_columns = csrRowIndeces_DofLoops[i+1]-csrRowIndeces_DofLoops[i];
	      double betaij = 1./(num_columns-1); // weight within smoothness indicator
	      double solni = u_dof_lstage[i]; // solution at time tn for the ith DOF
	      for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
		{ //loop in j (sparsity pattern)
		  int j = csrColumnOffsets_DofLoops[offset];
		  if (ENTROPY_VISCOSITY==1)
		    MaxEntResi = std::max(MaxEntResi,std::abs(EntResVector[j]));
		  double solnj = u_dof_lstage[j]; // solution at time tn for the jth DOF
		  alphai_numerator += betaij*(solnj-solni);
		  alphai_denominator += betaij*std::abs(solnj-solni);
		}
	      if (ENTROPY_VISCOSITY==1)
		MaxEntResVector[i] = MaxEntResi + 1E-14; // tolerance is used to avoid division by zero
	      alphai = std::abs(alphai_numerator)/(alphai_denominator+1E-14);
	      if (POWER_SMOOTHNESS_INDICATOR==0)
		psi[i] = 1.; //This means don't enhance the first order viscosity
	      else
		psi[i] = std::pow(alphai,POWER_SMOOTHNESS_INDICATOR); //NOTE: they use alpha^2 in the paper
	    }

	  /////////////////////////////////////////////
	  // ** LOOP IN DOFs FOR EDGE BASED TERMS ** //
	  /////////////////////////////////////////////
	  ij=0;
	  if (KUZMINS_METHOD==1)
	    {
	      for (int i=0; i<numDOFs; i++)
		{
		  double solni = u_dof_lstage[i]; // solution at time tn for the ith DOF
		  double ith_dissipative_term = 0;
		  double ith_low_order_dissipative_term = 0;
		  
		  // loop over the sparsity pattern of the i-th DOF
		  for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
		    {
		      int j = csrColumnOffsets_DofLoops[offset];
		      double solnj = u_dof_lstage[j]; // solution at time tn for the jth DOF
		      double dLij=0., dCij=0.;
		      
		      if (i != j) //NOTE: there is really no need to check for i!=j (see formula for ith_dissipative_term)
			{
			  // first-order dissipative operator
			  dLij = -std::max(0.0,std::max(TransportMatrix[ij],TransposeTransportMatrix[ij]));
			  // weight low-order dissipative matrix to make it higher order
			  dLij *= std::max(psi[i],psi[j]); 
			  
			  // high-order (entropy viscosity) dissipative operator 
			  double dEij = dLij;
			  if (ENTROPY_VISCOSITY==1)
			    {
			      double alphai = std::abs(EntResVector[i])/MaxEntResVector[i];
			      double alphaj = std::abs(EntResVector[j])/MaxEntResVector[j];
			      dEij *= std::max(alphai,alphaj);
			    }
			  
			  // artificial compression
			  double solij = 0.5*(solni+solnj);
			  double Compij = cK*std::max(solij*(1.0-solij),0.0)/(std::abs(solni-solnj)+1E-14);
			  dCij = dEij*std::max(1.0-Compij,0.0);
			  //dissipative terms
			  ith_dissipative_term += dCij*(solnj-solni);
			  ith_low_order_dissipative_term += dLij*(solnj-solni);
			  //dLij - dCij. This matrix is needed during FCT step
			  dt_times_dC_minus_dL[ij] = dt*(dLij - dCij);
			}
		      else //i==j
			{
			  // NOTE: this is incorrect. Indeed, dLii = -sum_{j!=i}(dLij) and similarly for dCii. 
			  // However, it is irrelevant since during the FCT step we do (dL-dC)*(solnj-solni)
			  dt_times_dC_minus_dL[ij]=0;
			}
		      //update ij
		      ij+=1;
		    }
		  // update residual 
		  globalResidual[i] += dt*ith_dissipative_term;
		  //globalResidual[i] += dt*ith_low_order_dissipative_term; //TMP
		}
	    } //END OF KUZMINS METHOD
	  else //EDGE BASED ON GUERMOND/POPOV'S
	    {
	      /*
	      for (int i=0; i<numDOFs; i++)
		{
		  double vxi = velx_tn_dof[i];
		  double vyi = vely_tn_dof[i]; // velocity at time tn for the ith DOF
		  double solni = u_dof_lstage[i]; // solution at time tn for the ith DOF
		  
		  double ith_flux_term = 0;
		  double ith_dissipative_term = 0;
		  double ith_low_order_dissipative_term = 0;
		  
		  // loop over the sparsity pattern of the i-th DOF
		  for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
		    {
		      int j = csrColumnOffsets_DofLoops[offset];
		      double vxj = velx_tn_dof[j];
		      double vyj = vely_tn_dof[j]; // velocity at time tn for the jth DOF
		      double solnj = u_dof_lstage[j]; // solution at time tn for the jth DOF

		      double dLij=0., dCij=0.;
		      ith_flux_term += solnj*(vxj*Cx[ij] + vyj*Cy[ij]);

		      if (i != j) //NOTE: there is really no need to check for i!=j (see formula for ith_dissipative_term)
			{
			  // first-order dissipative operator
			  dLij = -std::max(std::abs(vxi*Cx[ij] + vyi*Cy[ij]),std::abs(vxj*CTx[ij] + vyj*CTy[ij]));
			  // weight low-order dissipative matrix to make it higher order
			  dLij *= std::max(psi[i],psi[j]);

			  // high-order (entropy viscosity) dissipative operator 
			  double dEij = dLij;
			  if (ENTROPY_VISCOSITY==1)
			    {
			      double alphai = std::abs(EntResVector[i])/MaxEntResVector[i];
			      double alphaj = std::abs(EntResVector[j])/MaxEntResVector[j];
			      dEij *= std::max(alphai,alphaj);
			    }
			  
			  // artificial compression
			  double solij = 0.5*(solni+solnj);
			  double Compij = cK*std::max(solij*(1.0-solij),0.0)/(std::abs(solni-solnj)+1E-14);
			  dCij = dEij*std::max(1.0-Compij,0.0);
			  //dissipative terms
			  ith_dissipative_term += dCij*(solnj-solni);
			  ith_low_order_dissipative_term += dLij*(solnj-solni);
			  //dLij - dCij. This matrix is needed during FCT step
			  dt_times_dC_minus_dL[ij] = dLij - dCij;		      
			}
		      else //i==j
			{
			  // NOTE: this is incorrect. Indeed, dLii = -sum_{j!=i}(dLij) and similarly for dCii. 
			  // However, it is irrelevant since during the FCT step we do (dL-dC)*(solnj-solni)
			  dt_times_dC_minus_dL[ij]=0;
			}
		      //update ij
		      ij+=1;
		    }
		  // update residual 
		  globalResidual[i] += dt*(ith_flux_term + ith_dissipative_term);
		  //globalResidual[i] += dt*(ith_flux_term + ith_low_order_dissipative_term);
		
		  }
	      */
	    }
	}
      else // CELL BASED VISCOSITY/METHODS
	{
	  // ** COMPUTE QUANTITIES PER CELL (MQL) ** //
	  double entropy_max=-1.E10, entropy_min=1.E10, cell_entropy_mean, cell_volume, volume=0, entropy_mean=0;
	  double cell_vel_max, cell_entropy_residual;
	  double entropy_residual[nElements_global], vel_max[nElements_global];
	  double entropy_normalization_factor=1.0;	  
	  if (ENTROPY_VISCOSITY==1)
	    {
	      for(int eN=0;eN<nElements_global;eN++)
		{
		  cell_volume = 0;
		  cell_entropy_mean = 0;
		  cell_vel_max = 0;
		  cell_entropy_residual = 0;
		  //loop over quadrature points and compute integrands
		  for  (int k=0;k<nQuadraturePoints_element;k++)
		    {
		      //compute indeces and declare local storage
		      register int eN_k = eN*nQuadraturePoints_element+k,
			eN_k_nSpace = eN_k*nSpace,
			eN_nDOF_trial_element = eN*nDOF_trial_element;
		      register double un=0.0,unm1=0, grad_un[nSpace], vn[nSpace],
			jac[nSpace*nSpace],jacDet,jacInv[nSpace*nSpace],
			u_grad_trial[nDOF_trial_element*nSpace],
			dV,x,y,z,
			porosity=q_porosity[eN_k];
		      ck.calculateMapping_element(eN,k,mesh_dof,mesh_l2g,mesh_trial_ref,mesh_grad_trial_ref,jac,jacDet,jacInv,x,y,z);
		      //get the physical integration weight
		      dV = fabs(jacDet)*dV_ref[k];
		      //get the trial function gradients
		      ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,u_grad_trial);
		      //get the solution at quad point at tn and tnm1
		      ck.valFromDOF(u_dof_old,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],un);
		      ck.valFromDOF(u_dof_old_old,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],unm1);
		      //get the solution gradients at tn
		      ck.gradFromDOF(u_dof_old,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_un);
		      //velocity at tn
		      vn[0] = velocity[eN_k_nSpace];
		      vn[1] = velocity[eN_k_nSpace+1];
		      // compute entropy min and max
		      entropy_max = std::max(entropy_max,ENTROPY(un,uL,uR));
		      entropy_min = std::min(entropy_min,ENTROPY(un,uL,uR));
		      cell_entropy_mean += ENTROPY(un,uL,uR)*dV;
		      cell_volume += dV;
		      cell_vel_max = std::max(cell_vel_max,std::max(std::abs(vn[0]),std::abs(vn[1])));
		      cell_entropy_residual 
			= std::max(std::abs((ENTROPY(un,uL,uR) - ENTROPY(unm1,uL,uR))/dt
					    + vn[0]*ENTROPY_GRAD(un,grad_un[0],uL,uR)+vn[1]*ENTROPY_GRAD(un,grad_un[1],uL,uR) 
					    + ENTROPY(un,uL,uR)*(vn[0]+vn[1])),cell_entropy_residual);
		    }
		  volume += cell_volume;
		  entropy_mean += cell_entropy_mean;
		  vel_max[eN]=cell_vel_max;
		  entropy_residual[eN] = cell_entropy_residual;
		}//elements
	      entropy_mean /= volume;
	      // ** END OF CELL COMPUTATIONS (MQL) ** //
	      entropy_normalization_factor = std::max(std::abs(entropy_max-entropy_mean),
						      std::abs(entropy_min-entropy_mean));
	    } //ENTROPY_VISCOSITY==1

	  double Ct_sge = 4.0;	  
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
		  register double u=0.0, u_old=0.0, u_star=0.0, grad_u[nSpace], grad_u_old[nSpace], grad_u_star[nSpace],
		    m_star=0.0, dm_star=0.0, m=0.0, dm=0.0, 
		    f_star[nSpace], df_star[nSpace], f[nSpace], df[nSpace],
		    m_t=0.0,dm_t=0.0,
		    pdeResidual_u_star=0.0,
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
		  ck.valFromDOF(u_dof_lstage,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],u_old);
		  //get the solution gradients
		  ck.gradFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_u);
		  ck.gradFromDOF(u_dof_lstage,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_u_old);
		  //precalculate test function products with integration weights
		  for (int j=0;j<nDOF_trial_element;j++)
		    {
		      u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
		      for (int I=0;I<nSpace;I++)
			{
			  u_grad_test_dV[j*nSpace+I] = u_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
			}
		    }
		  //VRANS
		  porosity = q_porosity[eN_k];
		  // COMPUTE u and u_grad star to allow easy change between BACKWARD OR FORWARD EULER (for transport)
		  int IMPLICIT = (ENTROPY_VISCOSITY == 1 ? 0. : 1.);
		  u_star = IMPLICIT*u+(1-IMPLICIT)*u_old;
		  for (int I=0; I<nSpace; I++)
		    grad_u_star[I] = IMPLICIT*grad_u[I]+(1-IMPLICIT)*grad_u_old[I];
		  //
		  //calculate pde coefficients at quadrature points
		  //
		  evaluateCoefficients(&velocity[eN_k_nSpace],
				       u_star,
				       //VRANS
				       porosity,
				       //
				       m_star,
				       dm_star,
				       f_star,
				       df_star);
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
		      f_star[I] -= MOVING_DOMAIN*m_star*mesh_velocity[I];
		      df_star[I] -= MOVING_DOMAIN*dm_star*mesh_velocity[I];
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
		  if (ENTROPY_VISCOSITY==0)
		    {
		      //
		      //calculate subgrid error (strong residual and adjoint)
		      //
		      //calculate strong residual
		      pdeResidual_u_star = ck.Mass_strong(m_t) + ck.Advection_strong(df_star,grad_u_star);
		      //calculate adjoint
		      for (int i=0;i<nDOF_test_element;i++)
			{
			  // register int eN_k_i_nSpace = (eN_k*nDOF_trial_element+i)*nSpace;
			  // Lstar_u[i]  = ck.Advection_adjoint(df,&u_grad_test_dV[eN_k_i_nSpace]);
			  register int i_nSpace = i*nSpace;
			  Lstar_u[i]  = ck.Advection_adjoint(df_star,&u_grad_test_dV[i_nSpace]);
			}
		      //calculate tau and tau*Res
		      calculateSubgridError_tau(elementDiameter[eN],dm_t,df_star,cfl[eN_k],tau0);
		      calculateSubgridError_tau(Ct_sge,
						G,
						dm_t,
						df_star,
						tau1,
						cfl[eN_k]);					
		      tau = useMetrics*tau1+(1.0-useMetrics)*tau0;
		      subgridError_u = -tau*pdeResidual_u_star;
		      //
		      //calculate shock capturing diffusion
		      //
		      ck.calculateNumericalDiffusion(shockCapturingDiffusion,elementDiameter[eN],pdeResidual_u_star,grad_u_star,numDiff0);
		      ck.calculateNumericalDiffusion(shockCapturingDiffusion,sc_uref, sc_alpha,G,G_dd_G,pdeResidual_u_star,grad_u_star,numDiff1);
		      q_numDiff_u[eN_k] = useMetrics*numDiff1+(1.0-useMetrics)*numDiff0;
		    } //ENTROPY_VISCOSITY=0
		  else
		    {
		      double h=elementDiameter[eN]/degree_polynomial;
		      // CALCULATE CFL //
		      calculateCFL(h,df_star,cfl[eN_k]); // TODO: ADJUST SPEED IF MESH IS MOVING
		      // ** LINEAR DIFFUSION (MQL) ** //
		      // calculate linear viscosity 
		      double linear_viscosity = cMax*h*vel_max[eN]; // Cell based
		      
		      // ** ENTROPY VISCOSITY (MQL) ** //
		      double entropy_viscosity = cE*h*h*entropy_residual[eN]/entropy_normalization_factor;
		      q_numDiff_u[eN_k] = std::min(linear_viscosity,entropy_viscosity);
		      
		      // ** ARTIFICIAL COMPRESSION (MQL) ** //
		      double n_grad_u=0.0; 
		      for (int I=0; I<nSpace; I++)
			n_grad_u += grad_u_old[I]*grad_u_old[I];
		      n_grad_u = sqrt(n_grad_u);
		      double compression_factor = fmax(1-cK*fmax(u*(1.0-u),0.)/(h*n_grad_u+1.0e-8),0.);
		      q_numDiff_u[eN_k] *= compression_factor;

		      //q_numDiff_u[eN_k] = linear_viscosity;
		    }
		  // 
		  //update element residual 
		  // 		      		  
		  for(int i=0;i<nDOF_test_element;i++) 
		    { 
		      //register int eN_k_i=eN_k*nDOF_test_element+i,
		      //eN_k_i_nSpace = eN_k_i*nSpace,
		      register int i_nSpace=i*nSpace;
		      elementResidual_u[i] += 
			dt*ck.Mass_weak(m_t,u_test_dV[i]) + 
			dt*ck.Advection_weak(f_star,&u_grad_test_dV[i_nSpace]) + 
			dt*(ENTROPY_VISCOSITY==1 ? 0. : 1.)*ck.SubgridError(subgridError_u,Lstar_u[i]) + 
			dt*ck.NumericalDiffusion(q_numDiff_u_last[eN_k],grad_u_star,&u_grad_test_dV[i_nSpace]); 
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
		  //precalculate test function products with integration weights
		  for (int j=0;j<nDOF_trial_element;j++)
		    {
		      u_test_dS[j] = u_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;
		    }
		  //solution and gradients	
		  ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_trace_ref[ebN_local_kb*nDOF_test_element],u_ext);
		  ck.gradFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],u_grad_trial_trace,grad_u_ext);
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
		      elementResidual_u[i] += dt*ck.ExteriorElementBoundaryFlux(flux_ext,u_test_dS[i]);
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
      ////////////////////////////////////
      // COMPUTE QUANTITIES OF INTEREST //
      ////////////////////////////////////
      if (QUANTITIES_OF_INTEREST==1)
	{
	  // INTERIOR BOUNDARY // 
	  for(int eN=0;eN<nElements_global;eN++)
	    {
	      //declare local storage for local contributions and initialize
	      register double elementQuantDOFs[nDOF_test_element];
	      for (int i=0;i<nDOF_test_element;i++)
		elementQuantDOFs[i]=0.0;
	      //loop over quadrature points and compute integrands
	      for  (int k=0;k<nQuadraturePoints_element;k++)
		{
		  //compute indeces and declare local storage
		  register int eN_k = eN*nQuadraturePoints_element+k, eN_k_nSpace = eN_k*nSpace;
		  register double 
		    //for mass matrix contributions
		    u_test_dV[nDOF_trial_element], u_grad_test_dV[nDOF_test_element*nSpace], u_grad_trial[nDOF_trial_element*nSpace],
		    //for integration weight
		    jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace],dV,x,y,z;
		  //get the physical integration weight
		  ck.calculateMapping_element(eN,k,mesh_dof,mesh_l2g,mesh_trial_ref,mesh_grad_trial_ref,jac,jacDet,jacInv,x,y,z);
		  dV = fabs(jacDet)*dV_ref[k];
		  ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,u_grad_trial);
		  //precalculate test function products with integration weights for mass matrix terms
		  for (int j=0;j<nDOF_trial_element;j++)
		    {
		      u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
		      for (int I=0;I<nSpace;I++)
			u_grad_test_dV[j*nSpace+I] = u_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
		    }
		  // ith-LOOP //
		  for(int i=0;i<nDOF_test_element;i++) 
		    {
		      elementQuantDOFs[i] += // -int[vel*grad(wi)*dV] 
			- (velocity[eN_k_nSpace]*u_grad_test_dV[i*nSpace] + velocity[eN_k_nSpace+1]*u_grad_test_dV[i*nSpace+1]);
		    }
		}
	      // DISTRIBUTE // load cell based element into global vectors
	      for(int i=0;i<nDOF_test_element;i++) 
		{ 
		  int eN_i=eN*nDOF_test_element+i;
		  int gi = offset_u+stride_u*u_l2g[eN_i]; //global i-th index
		  quantDOFs[gi] += elementQuantDOFs[i];
		}//i
	    }//elements
	  // END OF INTERIOR BOUNDARY //
	  // EXTERIOR BOUNDARY //
	  for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) 
	    { 
	      register int ebN = exteriorElementBoundariesArray[ebNE], 
		eN  = elementBoundaryElementsArray[ebN*2+0],
		ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+0],
		eN_nDOF_trial_element = eN*nDOF_trial_element;
	      register double elementQuantDOFs[nDOF_test_element];
	      for (int i=0;i<nDOF_test_element;i++)
		elementQuantDOFs[i]=0.0;
	      for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) 
		{ 
		  register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
		    ebNE_kb_nSpace = ebNE_kb*nSpace,
		    ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb,
		    ebN_local_kb_nSpace = ebN_local_kb*nSpace;		      
		  register double
		    jac_ext[nSpace*nSpace],
		    jacDet_ext,
		    jacInv_ext[nSpace*nSpace],
		    boundaryJac[nSpace*(nSpace-1)],
		    metricTensor[(nSpace-1)*(nSpace-1)],
		    metricTensorDetSqrt,
		    dS,
		    u_test_dS[nDOF_test_element],		    
		    normal[nSpace],x_ext,y_ext,z_ext;
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
		  dS = metricTensorDetSqrt*dS_ref[kb];
		  // compute flow = v.n
		  double flow=0.0;
		  for (int I=0; I < nSpace; I++)
		    flow += normal[I]*ebqe_velocity_ext[ebNE_kb_nSpace+I];
		  //precalculate test function products with integration weights
		  for (int j=0;j<nDOF_trial_element;j++)
		    u_test_dS[j] = u_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;
		  // compute int[(vel*normal)*wi*dS]
		  for (int i=0;i<nDOF_test_element;i++)
		    elementQuantDOFs[i] += flow*u_test_dS[i];
		}//kb
	      // Distribute
	      for (int i=0;i<nDOF_test_element;i++)
		{
		  int eN_i = eN*nDOF_test_element+i;
		  int gi = offset_u+stride_u*u_l2g[eN_i]; // global i-th index
		  quantDOFs[gi] += elementQuantDOFs[i];
		}//i
	    }//ebNE
	  // END OF EXTERIOR BOUNDARY //
	}
      /////////////////////////////////////////////////
      // END OF COMPUTING AUX QUANTITIES OF INTEREST //
      /////////////////////////////////////////////////
    }
    
    void calculateResidual_entropy_viscosity(//element
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
			   //
			   int* u_l2g, 
			   double* elementDiameter,
			   int degree_polynomial,
			   double* u_dof,
			   double* u_dof_lstage,
			   double* u_dof_old,
			   double* u_dof_old_old,
			   double* velocity,
			   double* velx_dof,
			   double* vely_dof,
			   double* velz_dof,				   
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
			   int EDGE_VISCOSITY, 
			   int ENTROPY_VISCOSITY, 
			   double cE,
			   double cMax, 
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
			   double* CTx,
			   double* CTy, 
			   // PARAMETERS FOR 1st or 2nd ORDER MPP METHOD
			   int LUMPED_MASS_MATRIX, 
			   // FOR FCT
			   double* low_order_solution,
			   double* dt_times_dC_minus_dL,
			   double* min_u_bc,
			   double* max_u_bc,			   
			   // AUX QUANTITIES OF INTEREST 
			   double* quantDOFs)
    {
      double dt = 1./alphaBDF; // HACKED to work just for BDF1
      // Allocate space for the transport matrices
      // This is used for first order KUZMIN'S METHOD
      register double TransportMatrix[NNZ], TransposeTransportMatrix[NNZ];
      for (int i=0; i<NNZ; i++)
	{
	  TransportMatrix[i] = 0.;
	  TransposeTransportMatrix[i] = 0.;
	}

      // Allocate and init to zero the Entropy residual vector
      register double global_entropy_residual[numDOFs], lumped_mass_matrix[numDOFs], boundary_integral[numDOFs];
      for (int i=0; i<numDOFs; i++)
	{
	  lumped_mass_matrix[i]=0.;
	  global_entropy_residual[i]=0.;	   
	  boundary_integral[i]=0.;
	} 
	  
      //////////////////////////////////////////////
      // ** LOOP IN CELLS FOR CELL BASED TERMS ** //
      //////////////////////////////////////////////
      // HERE WE COMPUTE: 
      //    * lumped mass matrix
      //    * Time derivative term 
      //    * cell based CFL (for reference)
      //    * Entropy residual
      //    * Transport matrices
      for(int eN=0;eN<nElements_global;eN++)
	{
	  //declare local storage for local contributions and initialize
	  register double 
	    elementResidual_u[nDOF_test_element], 
	    elementEntResVector[nDOF_test_element],
	    element_lumped_mass_matrix[nDOF_test_element];
	  register double  elementTransport[nDOF_test_element][nDOF_trial_element];
	  register double  elementTransposeTransport[nDOF_test_element][nDOF_trial_element];
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      elementResidual_u[i]=0.0;
	      elementEntResVector[i]=0.0;
	      element_lumped_mass_matrix[i]=0.;
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
		//for mass matrix contributions
		u=0.0, un=0.0, unm1=0.0, grad_un[nSpace], vn[nSpace], 
		u_test_dV[nDOF_trial_element], 
		u_grad_trial[nDOF_trial_element*nSpace], 
		u_grad_test_dV[nDOF_test_element*nSpace],
		//for general use
		jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace],
		dV,x,y,z;
	      //get the physical integration weight
	      ck.calculateMapping_element(eN,k,mesh_dof,mesh_l2g,mesh_trial_ref,mesh_grad_trial_ref,jac,jacDet,jacInv,x,y,z);
	      dV = fabs(jacDet)*dV_ref[k];
	      //get the solution (of Newton's solver). To compute time derivative term
	      ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],u);
	      //get the solution at quad point at tn and tnm1 for entropy viscosity
	      ck.valFromDOF(u_dof_old,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],un);
	      ck.valFromDOF(u_dof_old_old,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],unm1);
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

	      //velocity at tn
	      vn[0] = velocity[eN_k_nSpace];
	      vn[1] = velocity[eN_k_nSpace+1];
	      //////////////////////////////
	      // CALCULATE CELL BASED CFL //
	      //////////////////////////////
	      calculateCFL(elementDiameter[eN]/degree_polynomial,vn,cfl[eN_k]); // TODO: ADJUST SPEED IF MESH IS MOVING
	      
	      //////////////////////////////////////////////
	      // CALCULATE ENTROPY RESIDUAL AT QUAD POINT //
	      //////////////////////////////////////////////
	      //double entropy_residual = 
	      //((un-unm1)/dt + vn[0]*grad_un[0] + vn[1]*grad_un[1] + un*(vn[0]+vn[1]))*DENTROPY(un,uL,uR);
	      double entropy_residual = (vn[0]*grad_un[0] + vn[1]*grad_un[1] + (vn[0]+vn[1])*un);
	      double DENTROPY_un = DENTROPY(un,uL,uR);
	      //double m_t = alphaBDF*u + q_m_betaBDF[eN_k];
	      //////////////
	      // ith-LOOP //
	      //////////////	      
	      for(int i=0;i<nDOF_test_element;i++) 
		{ 
		  // VECTOR OF ENTROPY RESIDUAL //
		  int eN_i=eN*nDOF_test_element+i;
		  int gi = offset_u+stride_u*u_l2g[eN_i]; //global i-th index
		  double uni = u_dof_old[gi];

		  elementEntResVector[i] += (DENTROPY_un - DENTROPY(uni,uL,uR))*entropy_residual*u_test_dV[i];
		  //elementEntResVector[i] += entropy_residual*u_test_dV[i];

		  element_lumped_mass_matrix[i] += u_test_dV[i];
		  elementResidual_u[i] += (u-un)*u_test_dV[i];		  
		  //elementResidual_u[i] += dt*m_t*u_test_dV[i];
		  
		  ///////////////
		  // j-th LOOP // To construct transport matrices
		  ///////////////
		  for(int j=0;j<nDOF_trial_element;j++) 
		    { 
		      int j_nSpace = j*nSpace;
		      int i_nSpace = i*nSpace;
		      elementTransport[i][j] += // -int[(vel.grad_wi)*wj*dx]
			ck.AdvectionJacobian_weak(vn,u_trial_ref[k*nDOF_trial_element+j],&u_grad_test_dV[i_nSpace]);
		      elementTransposeTransport[i][j] += // -int[(vel.grad_wj)*wi*dx]
			ck.AdvectionJacobian_weak(vn,u_trial_ref[k*nDOF_trial_element+i],&u_grad_test_dV[j_nSpace]);
		    }
		}//i
	      //save solution for other models 
	      q_u[eN_k] = u;
	      q_m[eN_k] = u*q_porosity[eN_k];
	    }
	  /////////////////
	  // DISTRIBUTE // load cell based element into global residual
	  ////////////////
	  for(int i=0;i<nDOF_test_element;i++) 
	    { 
	      int eN_i=eN*nDOF_test_element+i;
	      int gi = offset_u+stride_u*u_l2g[eN_i]; //global i-th index

	      // distribute lumped mass matrix
	      lumped_mass_matrix[gi]  += element_lumped_mass_matrix[i];	      
	      // distribute global residual for (lumped) mass matrix
	      globalResidual[gi] += elementResidual_u[i];
	      // distribute EntResVector 
	      global_entropy_residual[gi] += elementEntResVector[i];

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

      //////////////////////////////////////////////////////
      // ADD OUTFLOW BOUNDARY TERM TO TRANSPORT MATRICES //
      /////////////////////////////////////////////////////
      //   * Compute outflow boundary integral as a matrix; i.e., int_B[ (vel.normal)*wi*wj*dx]
      for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) 
	{ 
	  double min_u_bc_local = 1E10, max_u_bc_local = -1E10;
	  register int ebN = exteriorElementBoundariesArray[ebNE]; 
	  register int eN  = elementBoundaryElementsArray[ebN*2+0],
	    ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+0],
	    eN_nDOF_trial_element = eN*nDOF_trial_element;
	  register double elementResidual_u[nDOF_test_element], cell_flow;
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
		normal[nSpace],x_ext,y_ext,z_ext;
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
	      dS = metricTensorDetSqrt*dS_ref[kb];
	      //compute shape and solution information
	      ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_trace_ref[ebN_local_kb*nDOF_test_element],u_ext);
	      //precalculate test function products with integration weights
	      for (int j=0;j<nDOF_trial_element;j++)
		u_test_dS[j] = u_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;

	      //load the boundary values
	      bc_u_ext = ebqe_bc_u_ext[ebNE_kb];
	      //
	      //calculate the fluxes
	      //
	      double flow = 0.;
	      for (int I=0; I < nSpace; I++)
		flow += normal[I]*ebqe_velocity_ext[ebNE_kb_nSpace+I];

	      if (flow >= 0)
		{
		  dflux_ext = flow;
		  flux_ext = 0;
		}
	      else 
		{
		  dflux_ext = 0;
		  flux_ext = bc_u_ext*flow;
		}

	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  // elementResidual. This is to include the inflow boundary integral. 
		  // NOTE: here I assume that we use a Galerkin approach st nDOF_test_element = nDOF_trial_element
		  elementResidual_u[j] += flux_ext*u_test_dS[j];
		  register int ebN_local_kb_j=ebN_local_kb*nDOF_trial_element+j;
		  fluxTransport[j] = dflux_ext*u_trial_trace_ref[ebN_local_kb_j];
		}//j
	      ////////////////
	      // DISTRIBUTE //
	      ////////////////
	      for (int i=0;i<nDOF_test_element;i++)
		{
		  register int eN_i = eN*nDOF_test_element+i;
		  int gi = offset_u+stride_u*u_l2g[eN_i];
		  globalResidual[gi] += dt*elementResidual_u[i];
		  boundary_integral[gi] += elementResidual_u[i];

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
	      min_u_bc_local = std::min(ebqe_bc_u_ext[ebNE_kb], min_u_bc_local);
	      max_u_bc_local = std::max(ebqe_bc_u_ext[ebNE_kb], max_u_bc_local);
	    }//kb
	  // global min/max at boundary 
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      int eN_i = eN*nDOF_test_element+i;
	      int gi = offset_u+stride_u*u_l2g[eN_i]; //global i-th index
	      min_u_bc[gi] = std::min(min_u_bc_local,min_u_bc[gi]);
	      max_u_bc[gi] = std::max(max_u_bc_local,max_u_bc[gi]);		      
	    }
	}//ebNE
      // END OF ADDING BOUNDARY TERM TO TRANSPORT MATRICES and COMPUTING BOUNDARY INTEGRAL //
	  
      //////////////////////////////////////////
      // COMPUTE g vector and ENTROPY AT DOFs //
      //////////////////////////////////////////
      // gi=1/mi*sum_j(Cij*uj)
      register double gx[numDOFs], gy[numDOFs], eta[numDOFs], 
	alpha_numerator_pos[numDOFs], alpha_numerator_neg[numDOFs],
	alpha_denominator_pos[numDOFs], alpha_denominator_neg[numDOFs];
      int ij = 0;
      for (int i=0; i<numDOFs; i++)
	{
	  double solni = u_dof_old[i];
	  gx[i]=0.;
	  gy[i]=0.;
	  // entropy // 
	  eta[i] = ENTROPY(u_dof_old[i],uL,uR);
	  // for smoothness indicator // 
	  alpha_numerator_pos[i] = 0.;
	  alpha_numerator_neg[i] = 0.;
	  alpha_denominator_pos[i] = 0.;
	  alpha_denominator_neg[i] = 0.;

	  for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
	    {
	      int j = csrColumnOffsets_DofLoops[offset];
	      double solnj = u_dof_old[j];
	      // for gi vector
	      gx[i] += Cx[ij]*solnj;
	      gy[i] += Cy[ij]*solnj;

	      double alpha_num = solnj - solni;
	      alpha_numerator_pos[i] += alpha_num > 0. ? alpha_num : 0.;
	      alpha_numerator_neg[i] += alpha_num < 0. ? alpha_num : 0.;
	      alpha_denominator_pos[i] += alpha_num > 0. ? alpha_num : 0.;
	      alpha_denominator_neg[i] += alpha_num < 0. ? fabs(alpha_num) : 0.;
	      //update ij
	      ij+=1;
	    }	  
	  gx[i] /= lumped_mass_matrix[i];
	  gy[i] /= lumped_mass_matrix[i];
	}

      //////////////////////////////////////////////////////////
      // COMPUTE SMOOTHNESS INDICATOR and ENTROPY MIN and MAX //
      //////////////////////////////////////////////////////////
      register double psi[numDOFs], etaMax[numDOFs], etaMin[numDOFs],
	SumPos[numDOFs], SumNeg[numDOFs];
      for (int i=0; i<numDOFs; i++)
	{
	  double alphai;
	  double xi = mesh_dof[i*3+0];
	  double yi = mesh_dof[i*3+1];
	  // For eta min and max
	  etaMax[i] = fabs(eta[i]);
	  etaMin[i] = fabs(eta[i]);	  

	  SumPos[i] = 0.;
	  SumNeg[i] = 0.;

	  for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
	    { //loop in j (sparsity pattern)
	      int j = csrColumnOffsets_DofLoops[offset];
	      double xj = mesh_dof[j*3+0];
	      double yj = mesh_dof[j*3+1];
	      /////////////////////////////////
	      // COMPUTE ETA MIN AND ETA MAX // 
	      /////////////////////////////////
	      etaMax[i] = fmax(etaMax[i],std::abs(eta[j]));
	      etaMin[i] = fmin(etaMin[i],std::abs(eta[j]));
	      //////////////////////////////
	      // FOR SMOOTHNESS INDICATOR //
	      //////////////////////////////
	      double gj_times_x = gx[j]*(xj-xi) + gy[j]*(yj-yi);
	      SumPos[i] += gj_times_x > 0 ? gj_times_x : 0;
	      SumNeg[i] += gj_times_x < 0 ? gj_times_x : 0;
	    }
	  // Compute sigmaPos and sigmaNeg
	  double Sum = SumPos[i] + SumNeg[i];
	  double sigmaPosi = fmin(1.,(-SumNeg[i]+1E-15)/(SumPos[i]+1E-15));
	  double sigmaNegi = fmin(1.,(SumPos[i]+1E-15)/(-SumNeg[i]+1E-15));

	  double alpha_numi = fabs(sigmaPosi*alpha_numerator_pos[i] + sigmaNegi * alpha_numerator_neg[i]);
	  double alpha_deni = sigmaPosi*alpha_denominator_pos[i] + sigmaNegi * alpha_denominator_neg[i];

	  double alpha_num_threshold = 0.1;
	  if (BETAij == 1)
	    {
	      alpha_numi = fabs(alpha_numerator_pos[i] + alpha_numerator_neg[i]);
	      alpha_deni = alpha_denominator_pos[i] + alpha_denominator_neg[i];
	    }
	  //double g = std::sqrt(std::pow(gx[i],2)+std::pow(gy[i],2));
	  //if (g <= alpha_num_threshold) // Extrema and const state
	  //alphai = 1.;
	  //else 
	  alphai = alpha_numi/(alpha_deni+1E-15);
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
	  double solni = u_dof_old[i]; // solution at time tn for the ith DOF
	  double ith_dissipative_term = 0;
	  double ith_low_order_dissipative_term = 0;
	  double ith_flux_term = 0;
	  double dLii = 0.;

	  double one_over_entNormFactori = 1./(etaMax[i]-etaMin[i]+1E-15);
	  // loop over the sparsity pattern of the i-th DOF
	  for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
	    {
	      int j = csrColumnOffsets_DofLoops[offset];
	      double solnj = u_dof_old[j]; // solution at time tn for the jth DOF
	      double dLij=0., dCij=0.;

	      ith_flux_term += TransportMatrix[ij]*solnj;
	      if (i != j) //NOTE: there is really no need to check for i!=j (see formula for ith_dissipative_term)
		{
		  if (cE>=1000) // NO ENTROPY VISCOSITY
		    dLij = fmax(0.,fmax(psi[i]*TransportMatrix[ij],
					psi[j]*TransposeTransportMatrix[ij]));
		  else
		    dLij = fmax(0.,fmax(TransportMatrix[ij],
					TransposeTransportMatrix[ij]));
		  // high-order (entropy viscosity) dissipative operator 
		  double one_over_entNormFactorj = 1./(etaMax[j]-etaMin[j]+1E-15);
		  double dEVij = fmax(fabs(global_entropy_residual[i])*one_over_entNormFactori,
				      fabs(global_entropy_residual[j])*one_over_entNormFactorj);
		  // artificial compression
		  double solij = 0.5*(solni+solnj);
		  double Compij = cK*std::max(solij*(1.0-solij),0.0)/(std::abs(solni-solnj)+1E-14);
		  if (cE>=1000) // NO ENTROPY VISCOSITY 
		    dCij = dLij * fmax(1.0-Compij,0.0);
		  else
		    dCij = fmin(dLij,cE*dEVij) * fmax(1.0-Compij,0.0);
		  //dissipative terms
		  ith_dissipative_term += dCij*(solnj-solni);
		  ith_low_order_dissipative_term += dLij*(solnj-solni);
		  //dLij - dCij. This matrix is needed during FCT step
		  dt_times_dC_minus_dL[ij] = dt*(dCij - dLij);
		  dLii -= dCij;
		}
	      else //i==j
		{
		  // NOTE: this is incorrect. Indeed, dLii = -sum_{j!=i}(dLij) and similarly for dCii. 
		  // However, it is irrelevant since during the FCT step we do (dL-dC)*(solnj-solni)
		  dt_times_dC_minus_dL[ij]=0;
		}
	      //update ij
	      ij+=1;
	    }
	  double mi = lumped_mass_matrix[i];
	  // compute edge_based_cfl
	  edge_based_cfl[i] = 2.*fabs(dLii)/mi;

	  low_order_solution[i] = u_dof_old[i] - dt/mi*(ith_flux_term 
							+ boundary_integral[i] 
							- ith_low_order_dissipative_term);

	  // update residual
	  if (LUMPED_MASS_MATRIX==1)
	    globalResidual[i] = mi*(u_dof[i] - u_dof_old[i]) + dt*(ith_flux_term 
								   + boundary_integral[i]
								   - ith_dissipative_term);
	  else
	    globalResidual[i] += dt*(ith_flux_term - ith_dissipative_term);
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
			   int EDGE_VISCOSITY,
			   int ENTROPY_VISCOSITY,
			   int LUMPED_MASS_MATRIX)
    {
      double dt = 1./alphaBDF; // valid just for forward/backward euler
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
	      int IMPLICIT = (EDGE_VISCOSITY==1 ? 0. : 1.)*(ENTROPY_VISCOSITY==1 ? 0. : 1.); 	
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
			    dt*ck.MassJacobian_weak(dm_t,u_trial_ref[k*nDOF_trial_element+j],u_test_dV[i]) + 
			    dt*IMPLICIT*ck.AdvectionJacobian_weak(df,u_trial_ref[k*nDOF_trial_element+j],&u_grad_test_dV[i_nSpace]) +
			    dt*IMPLICIT*ck.SubgridErrorJacobian(dsubgridError_u_u[j],Lstar_u[i]) +
			    dt*IMPLICIT*ck.NumericalDiffusionJacobian(q_numDiff_u_last[eN_k],&u_grad_trial[j_nSpace],&u_grad_test_dV[i_nSpace]); //implicit
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
		  if (EDGE_VISCOSITY==0)
		    fluxJacobian_u_u[j]=dt*ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_u_u_ext,u_trial_trace_ref[ebN_local_kb_j]);
		  else 
		    fluxJacobian_u_u[j]=
		      (FIX_BOUNDARY_KUZMINS==0 ? 1. : 0.)*dt*ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_u_u_ext,u_trial_trace_ref[ebN_local_kb_j]);		  
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
			   int EDGE_VISCOSITY,
			   int ENTROPY_VISCOSITY,
			   int LUMPED_MASS_MATRIX)
    {
      double dt = 1./alphaBDF; // valid just for forward/backward euler
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

#ifndef NCLS3P_H
#define NCLS3P_H
#include <cmath>
#include <iostream>
#include "CompKernel.h"
#include "ModelFactory.h"

#define POWER_SMOOTHNESS_INDICATOR 1
#define LUMPED_MASS_MATRIX 0
#define FIX_BOUNDARY_KUZMINS 0

/////////////////////
//ENTROPY FUNCTION //
/////////////////////
// Power entropy //
#define entropy_power 2. // phiL and phiR are dummy variables
#define ENTROPY(phi,phiL,phiR) 1./entropy_power*std::pow(std::abs(phi),entropy_power)
#define ENTROPY_GRAD(phi,phix,phiL,phiR) std::pow(std::abs(phi),entropy_power-1.)*(phi>=0 ? 1 : -1)*phix
// Log entropy //
// LOG ENTROPY FOR LEVEL SET FROM 0 to 1
//#define ENTROPY(phi,phiL,phiR) std::log(std::abs((phi-phiL)*(phiR-phi))+1E-14)
//#define ENTROPY_GRAD(phi,phix,phiL,phiR) (phiL+phiR-2*phi)*phix*((phi-phiL)*(phiR-phi)>=0 ? 1 : -1)/(std::abs((phi-phiL)*(phiR-phi))+1E-14) 

namespace proteus
{
  class cppNCLS3P_base
  {
    //The base class defining the interface
  public:
    virtual ~cppNCLS3P_base(){}
    virtual void FCTStep(double dt, 
			 int NNZ, //number on non-zero entries on sparsity pattern
			 int numDOFs, //number of DOFs
			 double* lumped_mass_matrix, //lumped mass matrix (as vector)
			 double* soln, //DOFs of solution at time tn
			 double* solH, //DOFs of high order solution at tnp1
			 double* flux_plus_dLij_times_soln, //operators to construct low order solution
			 int* csrRowIndeces_DofLoops, //csr row indeces 
			 int* csrColumnOffsets_DofLoops, //csr column offsets 
			 double* MassMatrix, //mass matrix
			 double* dL_minus_dC, //low minus high order dissipative matrices
			 double* min_u_bc, //min/max value at BCs. If DOF is not at boundary then min=1E10, max=-1E10
			 double* max_u_bc
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
				   int lag_shockCapturing, /*mwf not used yet*/
				   double shockCapturingDiffusion,
		                   double sc_uref, double sc_alpha,
				   int* u_l2g, 
				   double* elementDiameter,
				   double* u_dof,
				   double* u_dof_old,	
				   double* u_dof_old_old,	
				   double* velocity,
				   double* q_m,
				   double* q_u,				   
				   double* q_n,
				   double* q_dH,
				   double* q_m_betaBDF,
                                   double* q_dV,
                                   double* q_dV_last,
				   double* cfl,
				   double* q_numDiff_u, 
				   double* q_numDiff_u_last, 
				   int offset_u, int stride_u, 
				   double* globalResidual,
				   int nExteriorElementBoundaries_global,
				   int* exteriorElementBoundariesArray,
				   int* elementBoundaryElementsArray,
				   int* elementBoundaryLocalElementBoundariesArray,
				   double* ebqe_velocity_ext,
				   int* isDOFBoundary_u,
				   double* ebqe_rd_u_ext,				   
				   double* ebqe_bc_u_ext,				   
				   double* ebqe_u,
				   // PARAMETERS FOR EDGE BASED STABILIZATION
				   int EDGE_VISCOSITY, 
				   int ENTROPY_VISCOSITY,
				   int numDOFs,
				   int NNZ,
				   int* csrRowIndeces_DofLoops,
				   int* csrColumnOffsets_DofLoops,
				   int* csrRowIndeces_CellLoops,
				   int* csrColumnOffsets_CellLoops,
				   int* csrColumnOffsets_eb_CellLoops,
				   // FOR FCT
				   double* flux_plus_dLij_times_soln,
				   double* dL_minus_dE, 
				   double* min_u_bc,
				   double* max_u_bc)=0;
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
				   int* u_l2g,
				   double* elementDiameter,
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
				   int* isDOFBoundary_u,
				   double* ebqe_rd_u_ext,
				   double* ebqe_bc_u_ext,
				   int* csrColumnOffsets_eb_u_u,
				   int EDGE_VISCOSITY)=0;
    virtual void calculateWaterline(//element
                                   int* wlc,
	                           double* waterline,
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
				   int lag_shockCapturing, /*mwf not used yet*/
				   double shockCapturingDiffusion,
		                   double sc_uref, double sc_alpha,
				   int* u_l2g, 
				   double* elementDiameter,
				   double* u_dof,double* u_dof_old,	
				   double* velocity,
				   double* q_m,
				   double* q_u,				   
				   double* q_n,
				   double* q_dH,
				   double* q_m_betaBDF,
				   double* cfl,
				   double* q_numDiff_u, 
				   double* q_numDiff_u_last, 
				   int offset_u, int stride_u, 
				   int nExteriorElementBoundaries_global,
				   int* exteriorElementBoundariesArray,
				   int* elementBoundaryElementsArray,
				   int* elementBoundaryLocalElementBoundariesArray,
				   int* elementBoundaryMaterialTypes,
				   double* ebqe_velocity_ext,
				   int* isDOFBoundary_u,
				   double* ebqe_bc_u_ext,				   
				   double* ebqe_u)=0;	
  };
  
  template<class CompKernelType,
	   int nSpace,
	   int nQuadraturePoints_element,
	   int nDOF_mesh_trial_element,
	   int nDOF_trial_element,
	   int nDOF_test_element,
	   int nQuadraturePoints_elementBoundary>
  class cppNCLS3P : public cppNCLS3P_base
  {
  public:
    const int nDOF_test_X_trial_element;
    CompKernelType ck;
  cppNCLS3P():      
    nDOF_test_X_trial_element(nDOF_test_element*nDOF_trial_element),
    ck()
    {}
    
    inline void evaluateCoefficients(const double v[nSpace],
				     const double& u,
				     const double grad_u[nSpace],
				     double& m,
				     double& dm,
				     double& H,
				     double dH[nSpace])
    {
      m = u;
      dm=1.0;
      H = 0.0;
      for (int I=0; I < nSpace; I++)
	{
	  H += v[I]*grad_u[I];
	  dH[I] = v[I];
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
 
    void exteriorNumericalFlux(const double n[nSpace],
			       const double& bc_u,
			       const double& u,
			       const double velocity[nSpace],
			       const double velocity_movingDomain[nSpace],
			       double& flux)
    {
      double flow_total=0.0,flow_fluid=0.0,flow_movingDomain=0.0;
      for (int I=0; I < nSpace; I++)
	{
	  flow_fluid += n[I]*velocity[I];
	  //flow_movingDomain -= n[I]*velocity_movingDomain[I];
	}
      flow_total = flow_fluid+flow_movingDomain;
      if (flow_total > 0.0)
	{
	  flux = u*flow_movingDomain;
	}
      else
	{
	  flux = bc_u*flow_movingDomain - flow_fluid*(u-bc_u);
	}
    }
    
    inline
    void exteriorNumericalFluxDerivative(const double n[nSpace],
					 const double velocity[nSpace],
					 const double velocity_movingDomain[nSpace],
					 double& dflux)
    {
      double flow_total=0.0,flow_fluid=0.0,flow_movingDomain=0.0;
      for (int I=0; I < nSpace; I++)
	{
	  flow_fluid += n[I]*velocity[I];
	  //flow_movingDomain -= n[I]*velocity_movingDomain[I];
	}
      flow_total=flow_fluid+flow_movingDomain;
      if (flow_total > 0.0)
	{
	  dflux = flow_movingDomain;
	}
      else
	{
	  dflux = -flow_fluid;
	}
    }

    void FCTStep(double dt, 
		 int NNZ, //number on non-zero entries on sparsity pattern
		 int numDOFs, //number of DOFs
		 double* lumped_mass_matrix, //lumped mass matrix (as vector)
		 double* soln, //DOFs of solution at time tn
		 double* solH, //DOFs of high order solution at tnp1
		 double* flux_plus_dLij_times_soln, //operators to construct low order solution
		 int* csrRowIndeces_DofLoops, //csr row indeces 
		 int* csrColumnOffsets_DofLoops, //csr column offsets 
		 double* MassMatrix, //mass matrix
		 double* dL_minus_dC, //low minus high order dissipative matrices
		 double* min_u_bc, //min/max value at BCs. If DOF is not at boundary then min=1E10, max=-1E10
		 double* max_u_bc)
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
	  solL[i] = solni-dt/mi*flux_plus_dLij_times_soln[i];

	  double mini=min_u_bc[i], maxi=max_u_bc[i]; // init min/max with value at BCs (NOTE: if no boundary then min=1E10, max=-1E10)
	  //double mini=1E10, maxi=-1E-10;
	  double Pposi=0, Pnegi=0;

	  // LOOP OVER THE SPARSITY PATTERN (j-LOOP)//
	  for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
	    {
	      int j = csrColumnOffsets_DofLoops[offset];
	      ////////////////////////
	      // COMPUTE THE BOUNDS //
	      ////////////////////////
	      mini = std::min(mini,soln[j]);
	      maxi = std::max(maxi,soln[j]);
	      
	      // i-th row of flux correction matrix 
	      FluxCorrectionMatrix[ij] = (((i==j) ? 1 : 0)*mi - MassMatrix[ij])*(solH[j]-soln[j] - (solHi-solni)) + dt*dL_minus_dC[ij]*(soln[j]-solni);

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
		((FluxCorrectionMatrix[ij]>0) ? std::min(Rposi,Rneg[j]) : std::min(Rnegi,Rpos[j])) * FluxCorrectionMatrix[ij];
	      //ith_Limiter_times_FluxCorrectionMatrix += FluxCorrectionMatrix[ij]; //TMP
	      //update ij
	      ij+=1;
	    }
	  solH[i] = solL[i] + 1./lumped_mass_matrix[i]*ith_Limiter_times_FluxCorrectionMatrix;
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
			   int* u_l2g, 
			   double* elementDiameter,
			   double* u_dof,
			   double* u_dof_old,			   
			   double* u_dof_old_old,
			   double* velocity,
			   double* q_m,
			   double* q_u,				   
			   double* q_n,
			   double* q_dH,
			   double* q_m_betaBDF,
                           double* q_dV,
                           double* q_dV_last,
			   double* cfl,
			   double* q_numDiff_u, 
			   double* q_numDiff_u_last, 
			   int offset_u, int stride_u, 
			   double* globalResidual,
			   int nExteriorElementBoundaries_global,
			   int* exteriorElementBoundariesArray,
			   int* elementBoundaryElementsArray,
			   int* elementBoundaryLocalElementBoundariesArray,
			   double* ebqe_velocity_ext,
			   int* isDOFBoundary_u,
			   double* ebqe_rd_u_ext,
			   double* ebqe_bc_u_ext,
			   double* ebqe_u,
			   // PARAMETERS FOR EDGE BASED STABILIZATION
			   int EDGE_VISCOSITY, 
			   int ENTROPY_VISCOSITY,
			   int numDOFs,
			   int NNZ,
			   int* csrRowIndeces_DofLoops,
			   int* csrColumnOffsets_DofLoops,
			   int* csrRowIndeces_CellLoops,
			   int* csrColumnOffsets_CellLoops,
			   int* csrColumnOffsets_eb_CellLoops,
			   // FOR FCT
			   double* flux_plus_dLij_times_soln,
			   double* dL_minus_dE, 
			   double* min_u_bc,
			   double* max_u_bc)
    {
      double dt = 1./alphaBDF; // HACKED to work just for BDF1
      if (EDGE_VISCOSITY==1)
	{
	  // Allocate space for the transport matrices
	  // This is used for first order KUZMIN'S METHOD
	  register double TransportMatrix[NNZ], TransposeTransportMatrix[NNZ];
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
	  //    * Transport matrices 
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
		    u=0.0, m=0.0, dm=0.0, H=0.0, dH[nSpace], m_t=0.0,dm_t=0.0, 
		    u_test_dV[nDOF_trial_element], 
		    //for entropy viscosity
		    un=0.0, unm1=0.0, grad_u[nSpace], grad_un[nSpace], vn[nSpace], 
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
		  // solution and grads at old times at quad points
		  ck.valFromDOF(u_dof_old,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],un);
		  ck.valFromDOF(u_dof_old_old,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],unm1);
		  //get the solution gradients at tn for entropy viscosity
		  ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,u_grad_trial);
		  ck.gradFromDOF(u_dof_old,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_un);
		  ck.gradFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_u);
		      
		  //precalculate test function products with integration weights for mass matrix terms
		  for (int j=0;j<nDOF_trial_element;j++)
		    {
		      u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
		      for (int I=0;I<nSpace;I++)
			u_grad_test_dV[j*nSpace+I] = u_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
		    }
		  //evaluate coefficients to compute time derivative (for term with mass matrix)
		  porosity = 1.0;// - q_vos[eN_k]; //TMP
		  evaluateCoefficients(&velocity[eN_k_nSpace],
				       u,
				       grad_u,
				       m,
				       dm,
				       H,
				       dH);
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
		  calculateCFL(elementDiameter[eN],dH,cfl[eN_k]); // TODO: ADJUST SPEED IF MESH IS MOVING
		      
		  // CALCULATE ENTROPY RESIDUAL AT QUAD POINT //
		  if (ENTROPY_VISCOSITY==1)
		    {
		      //velocity at tn for entropy viscosity
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
		      for(int j=0;j<nDOF_trial_element;j++) 
			{ 
			  int j_nSpace = j*nSpace;
			  int i_nSpace = i*nSpace;
			  // COMPUTE ELEMENT TRANSPORT MATRIX (MQL)
			  //elementTransport[i][j] += // int[(vel.grad_wj)*wi*dx]
			  //-ck.AdvectionJacobian_weak(dH,u_trial_ref[k*nDOF_trial_element+i],&u_grad_test_dV[j_nSpace]);
			  //elementTransposeTransport[i][j] += // -int[(vel.grad_wj)*wi*dx]
			  //-ck.AdvectionJacobian_weak(dH,u_trial_ref[k*nDOF_trial_element+j],&u_grad_test_dV[i_nSpace]);
			  elementTransport[i][j] += // int[(vel.grad_wj)*wi*dx]
			    ck.HamiltonianJacobian_weak(dH,&u_grad_trial[j_nSpace],u_test_dV[i]);
			  //-ck.AdvectionJacobian_weak(dH,u_trial_ref[k*nDOF_trial_element+i],&u_grad_test_dV[j_nSpace]);
			  elementTransposeTransport[i][j] += // -int[(vel.grad_wj)*wi*dx]
			    ck.HamiltonianJacobian_weak(dH,&u_grad_trial[i_nSpace],u_test_dV[j]);
			    //-ck.AdvectionJacobian_weak(dH,u_trial_ref[k*nDOF_trial_element+j],&u_grad_test_dV[i_nSpace]);
			}
		    }//i
		  //save solution for other models 
		  q_u[eN_k] = u;
		  q_m[eN_k] = m;
		}
	      /////////////////
	      // DISTRIBUTE // load cell based element into global residual
	      ////////////////
	      double h=elementDiameter[eN];
	      for(int i=0;i<nDOF_test_element;i++) 
		{ 
		  int eN_i=eN*nDOF_test_element+i;
		  int gi = offset_u+stride_u*u_l2g[eN_i]; //global i-th index
		  // distribute global residual for (lumped) mass matrix 
		      
		  if (LUMPED_MASS_MATRIX==1)
		    globalResidual[gi] += elementResidual_u[i]*(u_dof[gi]-u_dof_old[gi]);
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
	  if (FIX_BOUNDARY_KUZMINS==1)
	    {
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
		  for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) 
		    { 
		      register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
			ebNE_kb_nSpace = ebNE_kb*nSpace,
			ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb,
			ebN_local_kb_nSpace = ebN_local_kb*nSpace;
			  
		      register double u_ext=0.0,m_ext=0.0,dm_ext=0.0,f_ext[nSpace],df_ext[nSpace],
			flux_ext=0.0, dflux_ext=0.0, bc_u_ext=0.0,
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
		      //
		      //load the boundary values
		      //
		      bc_u_ext = isDOFBoundary_u[ebNE_kb]*ebqe_bc_u_ext[ebNE_kb]+(1-isDOFBoundary_u[ebNE_kb])*ebqe_rd_u_ext[ebNE_kb];
		      // TODO (MQL). Consider porosity and consider the possibility of having moving meshes 
		      //VRANS
		      porosity_ext = 1.0; // - ebqe_vos_ext[ebNE_kb]; //TMP
		      
		      // 
		      //calculate the numerical fluxes 
		      // 
		      double flow = 0.;
		      for (int I=0; I < nSpace; I++)
			flow += normal[I]*ebqe_velocity_ext[ebNE_kb_nSpace+I];
		      if (flow < 0)
			{
			  dflux_ext = -flow;
			  flux_ext = flow*bc_u_ext;
			}
		      else 
			dflux_ext = 0;
			  
		      for (int j=0;j<nDOF_trial_element;j++)
			{
			  // boundary integral related to data: int[(normal*velocity)*ubc*wi*dS]
			  elementResidual_u[j] += dt*ck.ExteriorElementBoundaryFlux(flux_ext,u_test_dS[j]);
			  // boundary integral related to solution: -int[(normal*velocity)*un*wi*dS]
			  register int ebN_local_kb_j=ebN_local_kb*nDOF_trial_element+j;
			  fluxTransport[j]= // -(normal*vel)*wj
			    ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_ext,u_trial_trace_ref[ebN_local_kb_j]);
			}//j
		      //
		      //update the global Transport Matrices and global residual
		      //
		      for (int i=0;i<nDOF_test_element;i++)
			{
			  register int eN_i = eN*nDOF_test_element+i;
			  // distribute boundary element to global residual
			  globalResidual[offset_u+stride_u*u_l2g[eN_i]] += elementResidual_u[i];	  
	
			  // distribute local transport matrix to global transport matrix
			  for (int j=0;j<nDOF_trial_element;j++)
			    {
			      register int ebN_i_j = ebN*4*nDOF_test_X_trial_element + i*nDOF_trial_element + j;
			      TransportMatrix[csrRowIndeces_CellLoops[eN_i] + csrColumnOffsets_eb_CellLoops[ebN_i_j]] // -(normal*vel)*wj*wi
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
	  else // FIX_BOUNDARY_KUZMINS=0
	    {
	      ////////////////////////////////////////////////////////
	      // COMPUTE (FULL) BOUNDARY INTEGRAL WITHOUT TREATMENT //
	      ////////////////////////////////////////////////////////
	      for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) 
		{ 
		  register int ebN = exteriorElementBoundariesArray[ebNE], 
		    eN  = elementBoundaryElementsArray[ebN*2+0],
		    ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+0],
		    eN_nDOF_trial_element = eN*nDOF_trial_element;
		  register double elementResidual_u[nDOF_test_element];
		  for (int i=0;i<nDOF_test_element;i++)
		    elementResidual_u[i]=0.0;
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
			H_ext=0.0,
			dH_ext[nSpace],
			bc_u_ext=0.0,
			bc_grad_u_ext[nSpace],
			bc_m_ext=0.0,
			bc_dm_ext=0.0,
			bc_H_ext=0.0,
			bc_dH_ext[nSpace],
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
			G[nSpace*nSpace],G_dd_G,tr_G,flux_ext;
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
			u_test_dS[j] = u_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;
		      //
		      //load the boundary values
		      //
		      bc_u_ext = isDOFBoundary_u[ebNE_kb]*ebqe_bc_u_ext[ebNE_kb]+(1-isDOFBoundary_u[ebNE_kb])*ebqe_rd_u_ext[ebNE_kb];		    
		      // 
		      //calculate the pde coefficients using the solution and the boundary values for the solution 
		      // 
		      evaluateCoefficients(&ebqe_velocity_ext[ebNE_kb_nSpace],
					   u_ext,
					   grad_u_ext,
					   m_ext,
					   dm_ext,
					   H_ext,
					   dH_ext);
		      evaluateCoefficients(&ebqe_velocity_ext[ebNE_kb_nSpace],
					   bc_u_ext,
					   bc_grad_u_ext,
					   bc_m_ext,
					   bc_dm_ext,
					   bc_H_ext,
					   bc_dH_ext);
		      //
		      //moving mesh
		      //
		      double velocity_ext[nSpace];
		      double mesh_velocity[3];
		      mesh_velocity[0] = xt_ext;
		      mesh_velocity[1] = yt_ext;
		      mesh_velocity[2] = zt_ext;
		      for (int I=0;I<nSpace;I++)
			velocity_ext[I] = - MOVING_DOMAIN*mesh_velocity[I];
		      // 
		      //calculate the numerical fluxes 
		      // 
		      exteriorNumericalFlux(normal,
					    bc_u_ext,
					    u_ext,
					    dH_ext,
					    velocity_ext,
					    flux_ext);
		      ebqe_u[ebNE_kb] = u_ext;	      
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
		      //globalResidual[offset_u+stride_u*u_l2g[eN_i]] += MOVING_DOMAIN*elementResidual_u[i];
		      globalResidual[offset_u+stride_u*u_l2g[eN_i]] += elementResidual_u[i];	  		      
		    }//i
		}//ebNE
	      // END OF COMPUTING (FULL) BOUNDARY INTEGRAL WITHOUT TREATMENT //
	    }

	  /////////////////////////////////////////////////////////////
	  // MULTIPLICATION OF TRANSPORT MATRIX TIMES SOLUTION AT tn //
	  /////////////////////////////////////////////////////////////
	  int ij=0;
	  for (int i=0; i<numDOFs; i++)
	    {		  
	      double ith_TransportMatrix_times_solution = 0.;
	      // LOOP OVER THE SPARSITY PATTERN (j-LOOP)//
	      for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
		{
		  int j = csrColumnOffsets_DofLoops[offset];
		  ith_TransportMatrix_times_solution += TransportMatrix[ij]*u_dof_old[j];
		  //update ij
		  ij+=1;
		}
	      globalResidual[i] += dt*ith_TransportMatrix_times_solution;
	      flux_plus_dLij_times_soln[i] = ith_TransportMatrix_times_solution; 
	    }
	  // END OF MULTIPLICATION OF TRANSPORT MATRIX TIMES SOLUTION AT tn

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
	      double solni = u_dof_old[i]; // solution at time tn for the ith DOF
	      for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
		{ //loop in j (sparsity pattern)
		  int j = csrColumnOffsets_DofLoops[offset];
		  if (ENTROPY_VISCOSITY==1)
		    MaxEntResi = std::max(MaxEntResi,std::abs(EntResVector[j]));
		  double solnj = u_dof_old[j]; // solution at time tn for the jth DOF
		  alphai_numerator += betaij*(solnj-solni);
		  alphai_denominator += betaij*std::abs(solnj-solni);
		}
	      if (ENTROPY_VISCOSITY==1)
		MaxEntResVector[i] = MaxEntResi + 1E-14; // tolerance is used to avoid division by zero
	      alphai = std::abs(alphai_numerator)/(alphai_denominator+1E-14);
	      psi[i] = std::pow(alphai,POWER_SMOOTHNESS_INDICATOR); //NOTE: they use alpha^2 in the paper
	    }	  
	  // END OF COMPUTING SMOOTHNESS INDICATOR 

	  /////////////////////////////////////////////
	  // ** LOOP IN DOFs FOR EDGE BASED TERMS ** //
	  /////////////////////////////////////////////
	  ij=0;
	  for (int i=0; i<numDOFs; i++)
	    {
	      double solni = u_dof_old[i]; // solution at time tn for the ith DOF
	      double ith_dissipative_term = 0;
	      double ith_low_order_dissipative_term = 0;
	      
	      // loop over the sparsity pattern of the i-th DOF
	      for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
		{
		  int j = csrColumnOffsets_DofLoops[offset];
		  double solnj = u_dof_old[j]; // solution at time tn for the jth DOF
		  double dLij=0.;
		  
		  if (i != j) //NOTE: there is really no need to check for i!=j (see formula for ith_dissipative_term)
		    {
		      // first-order dissipative operator
		      dLij = -std::max(0.0,std::max(TransportMatrix[ij],TransposeTransportMatrix[ij]));
		      // weight low-order dissipative matrix to make it higher order
		      //dLij *= std::max(psi[i],psi[j]); //TMP
		      
		      // high-order (entropy viscosity) dissipative operator 
		      double dEij = dLij*std::max(psi[i],psi[j]);
		      if (ENTROPY_VISCOSITY==1)
			{
			  double alphai = std::abs(EntResVector[i])/MaxEntResVector[i];
			  double alphaj = std::abs(EntResVector[j])/MaxEntResVector[j];
			  dEij *= std::max(alphai,alphaj);
			}
		      //dissipative terms
		      ith_dissipative_term += dEij*(solnj-solni);
		      ith_low_order_dissipative_term += dLij*(solnj-solni);
		      //dLij - dCij. This matrix is needed during FCT step
		      dL_minus_dE[ij] = dLij - dEij;
		    }
		  else //i==j
		    {
		      // NOTE: this is incorrect. Indeed, dLii = -sum_{j!=i}(dLij) and similarly for dCii. 
		      // However, it is irrelevant since during the FCT step we do (dL-dC)*(solnj-solni)
		      dL_minus_dE[ij]=0;
		    }
		  //update ij
		  ij+=1;
		}
	      // update residual 
	      globalResidual[i] += dt*ith_dissipative_term;
	      //globalResidual[i] += dt*ith_low_order_dissipative_term; //TMP
	      flux_plus_dLij_times_soln[i] += ith_low_order_dissipative_term;
	    }
	  //END OF KUZMINS METHOD
	}
      else // cell based stabilization 
	{
	  //cek should this be read in?
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
		  register double u=0.0,grad_u[nSpace],grad_u_old[nSpace],
		    m=0.0,dm=0.0,
		    H=0.0,dH[nSpace],
		    f[nSpace],df[nSpace],//for MOVING_DOMAIN
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
		    G[nSpace*nSpace],G_dd_G,tr_G;//,norm_Rv;
		  //
		  //compute solution and gradients at quadrature points
		  //
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
		  //save solution at quadrature points for other models to use
		  q_u[eN_k]=u;
		  for (int I=0;I<nSpace;I++)
		    q_n[eN_k_nSpace+I]  = grad_u[I];
		  //
		  //calculate pde coefficients at quadrature points
		  //
		  evaluateCoefficients(&velocity[eN_k_nSpace],
				       u,
				       grad_u,
				       m,
				       dm,
				       H,
				       dH);
		  //
		  //moving mesh
		  //
		  double mesh_velocity[3];
		  mesh_velocity[0] = xt;
		  mesh_velocity[1] = yt;
		  mesh_velocity[2] = zt;
		  for (int I=0;I<nSpace;I++)
		    {
		      f[I] = -MOVING_DOMAIN*m*mesh_velocity[I];
		      df[I] = -MOVING_DOMAIN*dm*mesh_velocity[I];
		    }
		  //
		  //calculate time derivative at quadrature points
		  //
		  if (q_dV_last[eN_k] <= -100)
		    q_dV_last[eN_k] = dV;
		  q_dV[eN_k] = dV;
		  ck.bdf(alphaBDF,
			 q_m_betaBDF[eN_k],
			 m,
			 dm,
			 m_t,
			 dm_t);
		  //
		  //calculate subgrid error (strong residual and adjoint)
		  //
		  //calculate strong residual
		  pdeResidual_u = ck.Mass_strong(m_t) +
		    ck.Hamiltonian_strong(dH,grad_u)+
		    MOVING_DOMAIN*ck.Advection_strong(df,grad_u);//cek I don't think all mesh motion will be divergence free so we may need to go add the divergence
		  
		  //calculate adjoint
		  for (int i=0;i<nDOF_test_element;i++)
		    {
		      //register int eN_k_i_nSpace = (eN_k*nDOF_trial_element+i)*nSpace;
		      //Lstar_u[i]  = ck.Hamiltonian_adjoint(dH,&u_grad_test_dV[eN_k_i_nSpace]);
		      register int i_nSpace = i*nSpace;
		      Lstar_u[i]  = ck.Hamiltonian_adjoint(dH,&u_grad_test_dV[i_nSpace]) + MOVING_DOMAIN*ck.Advection_adjoint(df,&u_grad_test_dV[i_nSpace]);
		    }
		  //calculate tau and tau*Res
		  double subgridErrorVelocity[nSpace];
		  for (int I=0;I<nSpace;I++)
		    subgridErrorVelocity[I] = dH[I] - MOVING_DOMAIN*df[I];

		  calculateSubgridError_tau(elementDiameter[eN],
					    dm_t,
					    subgridErrorVelocity,//dH,
					    cfl[eN_k],
					    tau0);
		  
		  calculateSubgridError_tau(Ct_sge,
					    G,
					    dm_t,
					    subgridErrorVelocity,//dH,
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
		  // 
		  //update element residual 
		  // 

		  for(int i=0;i<nDOF_test_element;i++) 
		    { 
		      //register int eN_k_i=eN_k*nDOF_test_element+i,
		      // eN_k_i_nSpace = eN_k_i*nSpace;
		      register int  i_nSpace=i*nSpace;
		      
		      elementResidual_u[i] += 
			dt*ck.Mass_weak(m_t,u_test_dV[i]) + 
			dt*ck.Hamiltonian_weak(H,u_test_dV[i]) + 
			dt*MOVING_DOMAIN*ck.Advection_weak(f,&u_grad_test_dV[i_nSpace])+
			dt*ck.SubgridError(subgridError_u,Lstar_u[i]) + 
			dt*ck.NumericalDiffusion(q_numDiff_u_last[eN_k],grad_u,&u_grad_test_dV[i_nSpace]); 
		      
		    }//i
		  //
		  //cek/ido todo, get rid of m, since u=m
		  //save momentum for time history and velocity for subgrid error
		  //save solution for other models 
		  //
		  q_m[eN_k] = m;
		  q_u[eN_k] = u;
		  for (int I=0;I<nSpace;I++)
		    {
		      int eN_k_I = eN_k*nSpace+I;
		      q_dH[eN_k_I] = dH[I];                     
		    }
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
		    H_ext=0.0,
		    dH_ext[nSpace],
		    //f_ext[nSpace],//MOVING_DOMAIN
		    //df_ext[nSpace],//MOVING_DOMAIN
		    //flux_ext=0.0,
		    bc_u_ext=0.0,
		    bc_grad_u_ext[nSpace],
		    bc_m_ext=0.0,
		    bc_dm_ext=0.0,
		    bc_H_ext=0.0,
		    bc_dH_ext[nSpace],
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
		    G[nSpace*nSpace],G_dd_G,tr_G,flux_ext;
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
		  bc_u_ext = isDOFBoundary_u[ebNE_kb]*ebqe_bc_u_ext[ebNE_kb]+(1-isDOFBoundary_u[ebNE_kb])*ebqe_rd_u_ext[ebNE_kb];

		  // 
		  //calculate the pde coefficients using the solution and the boundary values for the solution 
		  // 
		  evaluateCoefficients(&ebqe_velocity_ext[ebNE_kb_nSpace],
				       u_ext,
				       grad_u_ext,
				       m_ext,
				       dm_ext,
				       H_ext,
				       dH_ext);
		  evaluateCoefficients(&ebqe_velocity_ext[ebNE_kb_nSpace],
				       bc_u_ext,
				       bc_grad_u_ext,
				       bc_m_ext,
				       bc_dm_ext,
				       bc_H_ext,
				       bc_dH_ext);
		  //
		  //moving mesh
		  //
		  double velocity_ext[nSpace];
		  double mesh_velocity[3];
		  mesh_velocity[0] = xt_ext;
		  mesh_velocity[1] = yt_ext;
		  mesh_velocity[2] = zt_ext;
		  for (int I=0;I<nSpace;I++)
		    velocity_ext[I] = - MOVING_DOMAIN*mesh_velocity[I];
		  // 
		  //calculate the numerical fluxes 
		  // 
		  exteriorNumericalFlux(normal,
					bc_u_ext,
					u_ext,
					dH_ext,
					velocity_ext,
					flux_ext);
		  ebqe_u[ebNE_kb] = u_ext;
	      
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

		  //globalResidual[offset_u+stride_u*u_l2g[eN_i]] += MOVING_DOMAIN*elementResidual_u[i];
		  globalResidual[offset_u+stride_u*u_l2g[eN_i]] += elementResidual_u[i];	  
	      
		}//i
	    }//ebNE
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
			   int* u_l2g,
			   double* elementDiameter,
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
			   int* isDOFBoundary_u,
			   double* ebqe_rd_u_ext,
			   double* ebqe_bc_u_ext,
			   int* csrColumnOffsets_eb_u_u,
			   int EDGE_VISCOSITY)
    {
      double dt = 1./alphaBDF; // HACKED to work just for BDF1

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
		H=0.0,dH[nSpace],
		f[nSpace],df[nSpace],//MOVING_MESH
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
		G[nSpace*nSpace],G_dd_G,tr_G;
	      //
	      //calculate solution and gradients at quadrature points
	      //
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
	      //
	      //calculate pde coefficients and derivatives at quadrature points
	      //
	      evaluateCoefficients(&velocity[eN_k_nSpace],
				   u,
				   grad_u,
				   m,
				   dm,
				   H,
				   dH);
	      //
	      //moving mesh
	      //
	      double mesh_velocity[3];
	      mesh_velocity[0] = xt;
	      mesh_velocity[1] = yt;
	      mesh_velocity[2] = zt;
	      for (int I=0;I<nSpace;I++)
		{
		  f[I] = -MOVING_DOMAIN*m*mesh_velocity[I];
		  df[I] = -MOVING_DOMAIN*dm*mesh_velocity[I];
		}
	      //
	      //calculate time derivatives
	      //
	      ck.bdf(alphaBDF,
		     q_m_betaBDF[eN_k],
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
		  //int eN_k_i_nSpace = (eN_k*nDOF_trial_element+i)*nSpace;
		  //Lstar_u[i]=ck.Hamiltonian_adjoint(dH,&u_grad_test_dV[eN_k_i_nSpace]);
		  register int i_nSpace = i*nSpace;
		  Lstar_u[i]=ck.Hamiltonian_adjoint(dH,&u_grad_test_dV[i_nSpace]) + MOVING_DOMAIN*ck.Advection_adjoint(df,&u_grad_test_dV[i_nSpace]);
	      
		}
	      //calculate the Jacobian of strong residual
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  //int eN_k_j=eN_k*nDOF_trial_element+j;
		  //int eN_k_j_nSpace = eN_k_j*nSpace;
		  int j_nSpace = j*nSpace;
		  dpdeResidual_u_u[j]=ck.MassJacobian_strong(dm_t,u_trial_ref[k*nDOF_trial_element+j]) +
		    ck.HamiltonianJacobian_strong(dH,&u_grad_trial[j_nSpace]) +
		    MOVING_DOMAIN*ck.AdvectionJacobian_strong(df,&u_grad_trial[j_nSpace]);

		}
	      //tau and tau*Res
	      double subgridErrorVelocity[nSpace];
	      for (int I=0;I<nSpace;I++)
		subgridErrorVelocity[I] = dH[I] - MOVING_DOMAIN*df[I];

	      calculateSubgridError_tau(elementDiameter[eN],
					dm_t,
					subgridErrorVelocity,//dH,
	   	         		cfl[eN_k],
					tau0);
  
              calculateSubgridError_tau(Ct_sge,
                                        G,
					dm_t,
					subgridErrorVelocity,//dH,
					tau1,
				        cfl[eN_k]);
              tau = useMetrics*tau1+(1.0-useMetrics)*tau0;

	      for(int j=0;j<nDOF_trial_element;j++)
		dsubgridError_u_u[j] = -tau*dpdeResidual_u_u[j];
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
		      if (LUMPED_MASS_MATRIX==1)
			{
			  if (i==j)
			    elementJacobian_u_u[i][j] += u_test_dV[i];
			}
		      else
			{
			  elementJacobian_u_u[i][j] += 
			    dt*ck.MassJacobian_weak(dm_t,u_trial_ref[k*nDOF_trial_element+j],u_test_dV[i]) + 
			    dt*(EDGE_VISCOSITY==1 ? 0. : 1.)*ck.HamiltonianJacobian_weak(dH,&u_grad_trial[j_nSpace],u_test_dV[i]) +
			    dt*(EDGE_VISCOSITY==1 ? 0. : 1.)*MOVING_DOMAIN*ck.AdvectionJacobian_weak(df,u_trial_ref[k*nDOF_trial_element+j],&u_grad_test_dV[i_nSpace]) +
			    dt*(EDGE_VISCOSITY==1 ? 0. : 1.)*ck.SubgridErrorJacobian(dsubgridError_u_u[j],Lstar_u[i]) + 
			    dt*(EDGE_VISCOSITY==1 ? 0. : 1.)*ck.NumericalDiffusionJacobian(q_numDiff_u_last[eN_k],&u_grad_trial[j_nSpace],&u_grad_test_dV[i_nSpace]); 
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
		H_ext=0.0,
		dH_ext[nSpace],
		bc_grad_u_ext[nSpace],
		bc_H_ext=0.0,
		bc_dH_ext[nSpace],
		//f_ext[nSpace],
		//df_ext[nSpace],
		dflux_u_u_ext=0.0,
		bc_u_ext=0.0,
		//bc_grad_u_ext[nSpace],
		bc_m_ext=0.0,
		bc_dm_ext=0.0,
		//bc_f_ext[nSpace],
		//bc_df_ext[nSpace],
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
	      bc_u_ext = isDOFBoundary_u[ebNE_kb]*ebqe_bc_u_ext[ebNE_kb]+(1-isDOFBoundary_u[ebNE_kb])*ebqe_rd_u_ext[ebNE_kb];
	      // 
	      //calculate the internal and external trace of the pde coefficients 
	      // 
	      evaluateCoefficients(&ebqe_velocity_ext[ebNE_kb_nSpace],
				   u_ext,
				   grad_u_ext,
				   m_ext,
				   dm_ext,
				   H_ext,
				   dH_ext);
	      evaluateCoefficients(&ebqe_velocity_ext[ebNE_kb_nSpace],
				   bc_u_ext,
				   bc_grad_u_ext,
				   bc_m_ext,
				   bc_dm_ext,
				   bc_H_ext,
				   bc_dH_ext);
	      //
	      //moving domain
	      //
	      double velocity_ext[nSpace];
	      double mesh_velocity[3];
	      mesh_velocity[0] = xt_ext;
	      mesh_velocity[1] = yt_ext;
	      mesh_velocity[2] = zt_ext;
	      for (int I=0;I<nSpace;I++)
		{
		  velocity_ext[I] = - MOVING_DOMAIN*mesh_velocity[I];
		}
	      // 
	      //calculate the numerical fluxes 
	      // 
	      exteriorNumericalFluxDerivative(normal,
					      dH_ext,
					      velocity_ext,
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
		      //globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] += MOVING_DOMAIN*fluxJacobian_u_u[j]*u_test_dS[i];
		    }//j
		}//i
	    }//kb
	}//ebNE
    }//computeJacobian

    void calculateWaterline(//element
                           int* wlc,
	                   double* waterline,
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
			   int* u_l2g, 
			   double* elementDiameter,
			   double* u_dof,double* u_dof_old,			   
			   double* velocity,
			   double* q_m,
			   double* q_u,				   
			   double* q_n,
			   double* q_dH,
			   double* q_m_betaBDF,
			   double* cfl,
			   double* q_numDiff_u, 
			   double* q_numDiff_u_last, 
			   int offset_u, int stride_u, 
			   int nExteriorElementBoundaries_global,
			   int* exteriorElementBoundariesArray,
			   int* elementBoundaryElementsArray,
			   int* elementBoundaryLocalElementBoundariesArray,
                           int* elementBoundaryMaterialTypes,
			   double* ebqe_velocity_ext,
			   int* isDOFBoundary_u,
			   double* ebqe_bc_u_ext,
			   double* ebqe_u)
    {

      //  Tetrehedral elements specific extraction routine for waterline extraction
      //  Loops over boundaries and checks if boundary is infact a hull (hardwired check if mattype > 6)
      //  Extracts the nodal values of boundary triangle (4th point is dropped = hardwired assumption we are dealing with linear tet)
      //  Then computes an average value and position for both negative and positive values
      //  If both positive and negative values re found, and we are actually in a triangle containing the interface
      //  a linear interpolation of negative and positive average is reported as interface location (exact in case of linear tets) 

      for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) 
	{ 
          register int ebN = exteriorElementBoundariesArray[ebNE]; 
	  register int eN  = elementBoundaryElementsArray[ebN*2+0];
	  register int bN  = elementBoundaryLocalElementBoundariesArray[ebN*2+0];
	    
	  if (elementBoundaryMaterialTypes[ebN] >6)
	    {
	    	    
	    double val,x,y,z;	    
	    int pos=0, neg=0;
	    double xn=0.0, yn=0.0, zn=0.0, vn=0.0; 
	    double xp=0.0, yp=0.0, zp=0.0, vp=0.0;  
	    
            for (int j=0;j<nDOF_trial_element;j++) 
	    {
	      if (j != bN) {  
	       int eN_nDOF_trial_element_j      = eN*nDOF_trial_element + j;
	       int eN_nDOF_mesh_trial_element_j = eN*nDOF_mesh_trial_element + j;
               val = u_dof[u_l2g[eN_nDOF_trial_element_j]];
	       x = mesh_dof[mesh_l2g[eN_nDOF_mesh_trial_element_j]*3+0];
	       y = mesh_dof[mesh_l2g[eN_nDOF_mesh_trial_element_j]*3+1];
	       z = mesh_dof[mesh_l2g[eN_nDOF_mesh_trial_element_j]*3+2];
	       
               if (val < 0.0)
	       {
	         neg++;
		 vn+=val;
		 xn+=x;
		 yn+=y;
		 zn+=z;
	       }
	       else
	       {
	       	 pos++;
		 vp+=val;
		 xp+=x;
		 yp+=y;
		 zp+=z;
	       }
	      }	            
            } // trail for
	    

            if ((pos > 0) && (neg > 0) )
	    {
	      vp /= pos;
	      vn /= neg;
	  
	      double alpha = vp/(vp -vn);

	      waterline[wlc[0]*3 + 0] =  alpha*(xn/neg) + (1.0-alpha)*(xp/pos);
	      waterline[wlc[0]*3 + 1] =  alpha*(yn/neg) + (1.0-alpha)*(yp/pos);
	      waterline[wlc[0]*3 + 2] =  alpha*(zn/neg) + (1.0-alpha)*(zp/pos);
	      wlc[0]++;
	      	      
	    } // end value if	 
	  
	  } // end bnd mat check
	     
	}//ebNE

    } // calcWaterline

  };//cppNCLS3P

  inline cppNCLS3P_base* newNCLS3P(int nSpaceIn,
                                   int nQuadraturePoints_elementIn,
                                   int nDOF_mesh_trial_elementIn,
                                   int nDOF_trial_elementIn,
                                   int nDOF_test_elementIn,
                                   int nQuadraturePoints_elementBoundaryIn,
                                   int CompKernelFlag)
  {
    if (nSpaceIn == 2)
      return proteus::chooseAndAllocateDiscretization2D<cppNCLS3P_base,cppNCLS3P,CompKernel>(nSpaceIn,
										   nQuadraturePoints_elementIn,
										   nDOF_mesh_trial_elementIn,
										   nDOF_trial_elementIn,
										   nDOF_test_elementIn,
										   nQuadraturePoints_elementBoundaryIn,
										   CompKernelFlag);
    else
      return proteus::chooseAndAllocateDiscretization<cppNCLS3P_base,cppNCLS3P,CompKernel>(nSpaceIn,
										 nQuadraturePoints_elementIn,
										 nDOF_mesh_trial_elementIn,
										 nDOF_trial_elementIn,
										 nDOF_test_elementIn,
										 nQuadraturePoints_elementBoundaryIn,
										 CompKernelFlag);
    /* return proteus::chooseAndAllocateDiscretization<cppNCLS3P_base,cppNCLS3P>(nSpaceIn, */
    /* 									       nQuadraturePoints_elementIn, */
    /* 									       nDOF_mesh_trial_elementIn, */
    /* 									       nDOF_trial_elementIn, */
    /* 									       nDOF_test_elementIn, */
    /* 									       nQuadraturePoints_elementBoundaryIn, */
    /* 									       CompKernelFlag); */
  }
}//proteus
#endif

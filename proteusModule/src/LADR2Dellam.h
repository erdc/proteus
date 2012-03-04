#ifndef LADR2DELLAM_H
#define LADR2DELLAM_H
#include <cmath>
#include <iostream>

#define nSpace 2
#define nDOF_trial_element 3
#define nDOF_test_element 3
#define nDOF_test_X_trial_element 9


namespace LADR2Dellam
{
inline
void evaluateTestFunctionsOnElement(const double x[3], 
				    const double x0[3], const double x1[3], const double x2[3],
				    const double n0[2], const double n1[2], const double n2[3],
				    double w[nDOF_test_element])
{
  //lambda_i = 1 - (x-a_i).n_i/(a_j-a_i).n_i
  double lambda_0 = 1.0 -
    ((x[0] - x0[0])*n0[0] + (x[1]-x0[1])*n0[1])/((x1[0]-x0[0])*n0[0] + (x1[1]-x0[1])*n0[1]);
  double lambda_1 = 1.0 -
    ((x[0] - x1[0])*n1[0] + (x[1]-x1[1])*n1[1])/((x0[0]-x1[0])*n1[0] + (x0[1]-x1[1])*n1[1]);

  w[0] = lambda_0; w[1] = lambda_1; w[2] = 1.0 - lambda_0 - lambda_1;
}

inline
double valFromDOF_c(const double& dof,const double& v)
{
  return dof*v;
}


inline
double gradFromDOF_c(const double& dof,const double& grad_v_I)
{
  return dof*grad_v_I;
}



inline
void evaluateCoefficients_c(const double theta,
			    const double q[nSpace],
			    const double aL,
			    const double aT,
			    const double Dm,
			    const double& u,
			    double& m,
			    double& dm,
			    double f[nSpace],
			    double df[nSpace],
			    double a[nSpace*nSpace])
{
  m = u*theta;
  dm = theta;
  double normq = 0.0;
  for (int I=0; I < nSpace; I++)
    {
      f[I] = q[I]*u;
      df[I] = q[I];
      normq += q[I]*q[I];
    }
  normq = sqrt(normq); 
  if (normq > 0.0)
    {
      for (int I = 0; I < nSpace; I++)
	{
	  a[I*nSpace + I]  = aT*normq + (aL - aT)*q[I]*q[I]/normq + theta*Dm;
	  for (int J=I+1; J < nSpace; J++)
	    {
	      a[I*nSpace + J]  = aT*normq  + (aL - aT)*q[I]*q[J]/normq;
	      a[J*nSpace + I]  = a[I*nSpace*nSpace + I*nSpace + J];
	    }
	}
    }
  else
    for (int I = 0; I < nSpace*nSpace; I++)
      a[I] = 0.0;
}

  
inline
void backwardEuler_c(const double& dt, const double& m_old, const double& m, const double& dm, double& mt, double& dmt)
{  
  mt =(m-m_old)/dt;
  dmt = dm/dt;
}


inline
double Mass_weak_c(const double& mt,
		   const double& w_dV)
{
  return mt*w_dV;
}

inline
double MassJacobian_weak_c(const double& dmt,
			   const double& v,
			   const double& w_dV)
{
  return dmt*v*w_dV;
}


inline
double Mass_strong_c(const double& mt)
{
  return mt; 
}


inline
double MassJacobian_strong_c(const double& dmt,
			     const double& v)
{
  return dmt*v;
}


inline
double Mass_adjoint_c(const double& dmt,
		      const double& w_dV)
{
  return dmt*w_dV;
}




inline
double Advection_weak_c(const double  f[nSpace],
			const double grad_w_dV[nSpace])
{
  double tmp=0.0;
  for(int I=0;I<nSpace;I++)
    tmp -= f[I]*grad_w_dV[I];
  return tmp;
}


inline
double AdvectionJacobian_weak_c(const double df[nSpace],
				const double& v,
				const double grad_w_dV[nSpace])
{
  double tmp=0.0;
  for(int I=0;I<nSpace;I++)
    tmp -= df[I]*v*grad_w_dV[I];
  return tmp;
}


inline
double Advection_strong_c(const double df[nSpace],
			  const double grad_u[nSpace])
{
  double tmp=0.0;
  for(int I=0;I<nSpace;I++)
    tmp += df[I]*grad_u[I];
  return tmp;
}


inline
double AdvectionJacobian_strong_c(const double df[nSpace],
				  const double grad_v[nSpace])
{
  double tmp=0.0;
  for(int I=0;I<nSpace;I++)
    tmp += df[I]*grad_v[I];
  return tmp;
}

inline
double Advection_adjoint_c(const double df[nSpace],
			   const double grad_w_dV[nSpace])
{
  double tmp=0.0;
  for(int I=0;I<nSpace;I++)
    tmp -= df[I]*grad_w_dV[I];
  return tmp;
}



inline
double Diffusion_weak_c(double dt,
			int* rowptr,
			int* colind,
			double* a,
			const double grad_phi[nSpace],
			const double grad_w_dV[nSpace])
{
  double tmp=0.0;
  for(int I=0;I<nSpace;I++)
    for (int m=rowptr[I];m<rowptr[I+1];m++)
      tmp += dt*a[m]*grad_phi[colind[m]]*grad_w_dV[I];
  return tmp;
}

inline
double SimpleDiffusionJacobian_weak_c(double dt,
				      int* rowptr,
				      int* colind,
				      double* a,
				      const double grad_v[nSpace],
				      const double grad_w_dV[nSpace])
{
  double tmp=0.0;
  for (int I=0;I<nSpace;I++)
    for(int m=rowptr[I];m<rowptr[I+1];m++)
      tmp += dt*a[m]*grad_v[colind[m]]*grad_w_dV[I];
  return tmp;
}

inline
double Reaction_weak_c(const double& r,
		       const double& w_dV)
{
  return r*w_dV;
}

inline
double ReactionJacobian_weak_c(const double& dr,
			       const double& v,
			       const double& w_dV)
{
  return dr*v*w_dV;
}


inline
double Reaction_strong_c(const double& r)
{
  return r; 
}

inline
double ReactionJacobian_strong_c(const double& dr,
				     const double& v)
{
  return dr*v;
}

inline
double Reaction_adjoint_c(const double& dr,
			  const double& w_dV)
{
  return dr*w_dV;
}
inline 
void calculateCourantNumber_c(const double& elementDiameter,
			      const double& dm,
			      const double df[nSpace],
			      double& cfl)
{
  double h,nrm_v,oneByAbsdt;
  h = elementDiameter;
  nrm_v=0.0;
  for(int I=0;I<nSpace;I++)
    nrm_v+=df[I]*df[I];
  nrm_v = sqrt(nrm_v)/dm;
  cfl = nrm_v/h;

}
//need advection diffusion formula 
inline
void calculateSubgridError_tau_c(const double& elementDiameter,
				 const double& dmt,
				 const double df[nSpace],
				 double& cfl,
				 double& tau)
{
  double h,nrm_v,oneByAbsdt;
  h = elementDiameter;
  nrm_v=0.0;
  for(int I=0;I<nSpace;I++)
    nrm_v+=df[I]*df[I];
  nrm_v = sqrt(nrm_v);
  cfl = nrm_v/h;
  oneByAbsdt =  fabs(dmt);
  tau = 1.0/(2.0*nrm_v/h + oneByAbsdt + 1.0e-8);
}

inline
void calculateSubgridError_tauRes_c(const double& tau,
				    const double& pdeResidual,
				    double& subgridError)
{
  subgridError = -tau*pdeResidual;
}

inline
void calculateSubgridErrorDerivatives_tauRes_c(const double& tau,
					       const double dpdeResidual_du[nDOF_trial_element],
					       double dsubgridError_du[nDOF_trial_element])
{
  for (int j=0;j<nDOF_trial_element;j++)
    {
      dsubgridError_du[j] = -tau*dpdeResidual_du[j];
     }
}


inline
void calculateNumericalDiffusion_c(const double& shockCapturingDiffusion,
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
double SubgridError_c(const double& error,
		      const double& Lstar_w_dV)
{
  return -error*Lstar_w_dV;
}


inline
double SubgridErrorJacobian_c(const double& derror,
				  const double& Lstar_w_dV)
{
  return -derror*Lstar_w_dV;
}


inline
double NumericalDiffusion_c(const double& numDiff,
			    const double grad_u[nSpace],
			    const double grad_w_dV[nSpace])
{
  double tmp=0.0;
  for (int I=0;I<nSpace;I++)
    tmp +=  numDiff*grad_u[I]*grad_w_dV[I];
  return tmp;
}


inline
double NumericalDiffusionJacobian_c(const double& numDiff,
					const double grad_v[nSpace],
					const double grad_w_dV[nSpace])
{
  double tmp=0.0;
  for (int I=0;I<nSpace;I++)
    tmp += numDiff*grad_v[I]*grad_w_dV[I];
  return tmp;
}



inline
void exteriorNumericalAdvectiveFlux_c(const int& isDOFBoundary_u,
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
      //std::cout<<"flux boundary u and bc_u "<<u<<'\t'<<bc_u<<std::endl;
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
      //std::cout<<"flux boundary flux and flow"<<flux<<'\t'<<flow<<std::endl;
    }
  else
    {
      if (flow >= 0.0)
	{
	  flux = u*flow;
	}
      else
	{
	  std::cout<<"warning: open boundary with no external trace, setting to zero for inflow"<<std::endl;
	  flux = 0.0;
	}

    }
}

inline
void exteriorNumericalAdvectiveFluxDerivative_c(const int& isDOFBoundary_u,
						const int& isFluxBoundary_u,
						const double n[nSpace],
						const double velocity[nSpace],
						double& dflux)
{
  double flow=0.0;
  for (int I=0; I < nSpace; I++)
    flow += n[I]*velocity[I];
  //double flow=n[0]*velocity[0]+n[1]*velocity[1]+n[2]*velocity[2];
  if (flow >= 0.0)
    dflux=flow;//default to outflow
  else
    dflux = 0.0;
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
  if (isFluxBoundary_u == 1)
    {
      dflux = 0.0;
    }
}

inline
double ExteriorElementBoundaryFlux_c(const double& flux,
				     const double& w_dS)
{
  return flux*w_dS;
}


inline
double ExteriorNumericalAdvectiveFluxJacobian_c(const double& dflux_left,
						const double& v)
{
  return dflux_left*v;
}

inline
void exteriorOutflowFlux_c(const double n[nSpace],
			   const double& u,
			   const double velocity[nSpace],
			   double& flux)
{

  double flow=0.0;
  for (int I=0; I < nSpace; I++)
    flow += n[I]*velocity[I];

  if (flow >= 0.0)
    {
      flux = u*flow;
    }
}

inline
void exteriorOutflowFluxDerivative_c(const double n[nSpace],
				     const double velocity[nSpace],
				     double& dflux)
{
  double flow=0.0;
  for (int I=0; I < nSpace; I++)
    flow += n[I]*velocity[I];
  if (flow >= 0.0)
    dflux=flow;
  else
    dflux = 0.0;
}
inline
void tagInflowPointForTracking_c(const double& tn,
				 const double& tnp1,
				 const double& t,
				 const double n[nSpace],
				 const double velocity_n[nSpace],
				 const double velocity_np1[nSpace],
				 int& inflow)
{
  double flow=0.0;
  double dt = tnp1-tn;
  for (int I=0; I < nSpace; I++)
    flow += n[I]*((t-tn)/dt*velocity_np1[I] + (tnp1-t)/dt*velocity_n[I]);
  inflow = -1;
  if (flow < 0.0)
    {
      inflow = 1;
    }
}

}//namespace
extern "C"
{
 

  void calculateResidual_LADR2Dellam (double theta,
				      double aL,
				      double aT,
				      double Dm,
				      double dtnp1,
				      int nElements_global,//mesh representation
				      int nNodes_global,
				      int nNodes_element,
				      int nElementBoundaries_element,
				      int nQuadraturePoints_element,
				      const double * nodeArray,
				      const int * elementNodesArray,
				      const int * elementNeighborsArray, //local boundary id is associated with node across from boundary 
				      const double * elementBoundaryOuterNormals, //local element boundary outer normal constant on face
				      const double* dV,
				      const double* x_track,
				      const double* t_track,
				      const int* element_track,
				      const int* flag_track,
				      double* x_dt,
				      int* u_l2g, 
				      double* elementDiameter,
				      double* u_dof,
				      double* u_trial, 
				      double* u_grad_trial, 
				      double* u_test_dV, 
				      double* u_grad_test_dV, 
				      double* velocity,
				      double* q_u,
				      double* q_m,
				      double* q_m_last,
				      double* q_cfl,
				      double* q_elementResidual_u,
				      int* sdInfo_u_rowptr, int* sdInfo_u_colind,
				      int offset_u, int stride_u, 
				      double* globalResidual,
				      int nExteriorElementBoundaries_global,
				      int nQuadraturePoints_elementBoundary,
				      int* exteriorElementBoundariesArray,
				      int* elementBoundaryElementsArray,
				      int* elementBoundaryLocalElementBoundariesArray,
				      double* u_trial_ext,
				      double* u_grad_trial_ext,
				      double* ebqe_velocity_ext,
				      double* ebqe_n_ext,
				      int* isDOFBoundary_u,
				      double* ebqe_bc_u_ext,
				      int* isFluxBoundary_u,
				      double* ebqe_bc_flux_u_ext,
				      double* ebqe_outflow_flux_last,
				      double* u_test_dS_ext,
				      double* ebqe_u,
				      double* ebqe_outflow_flux);
  void calculateJacobian_LADR2Dellam(double theta,
				     double aL,
				     double aT,
				     double Dm,
				     double dtnp1,
				     int nElements_global,
				     int nQuadraturePoints_element,
				     double* x_dt,
				     int* u_l2g,
				     double* elementDiameter,
				     double* u_dof, 
				     double* u_trial, 
				     double* u_grad_trial, 
				     double* u_test_dV, 
				     double* u_grad_test_dV, 
				     double* velocity,
				     double* q_m_last, 
				     double* q_cfl,
				     int * sdInfo_u_rowptr, int * sdInfo_u_colind,
				     int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
				     double* globalJacobian,
				     int nExteriorElementBoundaries_global,
				     int nQuadraturePoints_elementBoundary,
				     int* exteriorElementBoundariesArray,
				     int* elementBoundaryElementsArray,
				     int* elementBoundaryLocalElementBoundariesArray,
				     double* u_trial_ext,
				     double* u_grad_trial_ext,
				     double* ebqe_velocity_ext,
				     double* ebqe_n,
				     int* isDOFBoundary_u,
				     double* ebqe_bc_u_ext,
				     int* isFluxBoundary_u,
				     double* ebqe_bc_flux_u_ext,
				     double* u_test_dS_ext,
				     int* csrColumnOffsets_eb_u_u);
  void markInflowBoundaryPoints2d(double tn, 
				double tnp1,
				double t,
				int nExteriorElementBoundaries_global,
				int nQuadraturePoints_elementBoundary,
				const int* exteriorElementBoundariesArray,
				const int* elementBoundaryElementsArray,
				const int* elementBoundaryLocalElementBoundariesArray,
				const double* ebqe_x,
				const double* ebqe_n_ext,
				const double* ebqe_velocity_ext_last,
				const double* ebqe_velocity_ext,
				const int* isDOFBoundary_u,
				const int* isFluxBoundary_u,
				int* element_track,//which element is point in
				int* flag_track);  //>=0 track, -1 don't

  void accumulateInflowFluxInGlobalResidual2d(int nElements_global,//mesh representation
					    int nNodes_global,
					    int nNodes_element,
					    int nElementBoundaries_element,
					    int nExteriorElementBoundaries_global,
					    int nQuadraturePoints_elementBoundary,
					    const double * nodeArray,
					    const int * elementNodesArray,
					    const int * elementNeighborsArray, //local boundary id is associated with node across from boundary 
					    const int* exteriorElementBoundariesArray,
					    const int* elementBoundaryElementsArray,
					    const int* elementBoundaryLocalElementBoundariesArray,
					    const double * elementBoundaryOuterNormals, //local element boundary outer normal constant on face
					    double tp, //time level for integration points
					    double timeWeight, //temporal integration weight (uniform for now)
					    const double* dS,  //spatial integration weights
					    const double* x_track_ext, //location of forward tracked integration points, assumed for now 
					    //nExteriorElementBoundaries_global x nQuadraturePoints_elementBoundary
					    const double* t_track_ext,     //time forward tracked points stopped (t^{n+1} or earlier if exited domain)
					    const int* element_track_ext,  //element each forward tracked point ended up in
					    const int* flag_track_ext,     //id for each point, 0 -- interior, 1 exited domain, -1 didn't track for some reason
					    const int* u_l2g, 
					    const double* u_dof,
					    double* q_elementResidual_u, 
					    const int* sdInfo_u_rowptr, 
					    const int* sdInfo_u_colind,
					    int offset_u, int stride_u, 
					    const int* isFluxBoundary_u,
					    const double* ebqe_bc_flux_u_ext,
					    double* globalResidual);




}//extern "C"
#endif

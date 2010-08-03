#ifndef RBLES2P_H
#define RBLES2P_H
#include <cmath>
#include <iostream>

#define nSpace 3
#define nQuadraturePoints_element 5
#define nDOF_trial_element 4
#define nDOF_test_element 4
#define nDOF_test_X_trial_element 16
#define nQuadraturePoints_elementBoundary 4
#define nElementBoundaries_element 4
namespace RBLES2P
{
inline
double smoothedHeaviside(double eps, double phi)
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
  return H;
}

inline
double smoothedHeaviside_integral(double eps, double phi)
{
  double HI;
  if (phi > eps)
    {
      HI= phi - eps + 	0.5*(eps + 0.5*eps*eps/eps - eps*cos(M_PI*eps/eps)/(M_PI*M_PI)) - 0.5*((-eps) + 0.5*(-eps)*(-eps)/eps - eps*cos(M_PI*(-eps)/eps)/(M_PI*M_PI));
    }
  else if (phi < -eps)
    {
      HI=0.0;
    }
  else
    {
      HI = 0.5*(phi + 0.5*phi*phi/eps - eps*cos(M_PI*phi/eps)/(M_PI*M_PI)) - 0.5*((-eps) + 0.5*(-eps)*(-eps)/eps - eps*cos(M_PI*(-eps)/eps)/(M_PI*M_PI));
    }
  return HI;
}

inline
double smoothedDirac(double eps, double phi)
{
  double d;
  if (phi > eps)
    d=0.0;
  else if (phi < -eps)
    d=0.0;
  else
    d = 0.5*(1.0 + cos(M_PI*phi/eps))/eps;
  return d;
}

inline
double valFromDOF_c(const double& dof,const double& v)
{
  return dof*v;
}

inline
void valFromDOF(int eN, int k, int j, double* u,double* dof,int* l2g,double* v)
{
  u[eN*nQuadraturePoints_element+
    k]
    += valFromDOF_c(dof[l2g[eN*nDOF_trial_element+
			    j]],
		    v[eN*nQuadraturePoints_element*nDOF_trial_element+
		      k*nDOF_trial_element+
		      j]);
}


inline
void valFromDOF_ext(int ebNE, int eN, int k, int j, double* u, double* dof, int* l2g, double* v)
{
  u[ebNE*nQuadraturePoints_elementBoundary+
    k] 
    += valFromDOF_c(dof[l2g[eN*nDOF_trial_element+j]],
		    v[ebNE*nQuadraturePoints_elementBoundary*nDOF_trial_element+
		      k*nDOF_trial_element+
		      j]);
}


inline
double gradFromDOF_c(const double& dof,const double& grad_v_I)
{
  return dof*grad_v_I;
}


inline
void gradFromDOF(int eN, int k, int j, double* grad_u,double* dof,int* l2g,double* grad_v)
{
  for (int I=0;I<nSpace;I++)
    grad_u[eN*nQuadraturePoints_element*nSpace+
	   k*nSpace+
	   I] += gradFromDOF_c(dof[l2g[eN*nDOF_trial_element+
				       j]],
			       grad_v[eN*nQuadraturePoints_element*nDOF_trial_element*nSpace+
				      k*nDOF_trial_element*nSpace+
				      j*nSpace+
				      I]);
}

inline
void gradFromDOF_ext(int ebNE, int eN, int k, int j, double* grad_u, double* dof, int* l2g, double* grad_v)
{
  for (int I=0;I<nSpace;I++)
    grad_u[ebNE*nQuadraturePoints_elementBoundary*nSpace+
	   k*nSpace+I]
      += gradFromDOF_c(dof[l2g[eN*nDOF_trial_element+j]],
		       grad_v[ebNE*nQuadraturePoints_elementBoundary*nDOF_trial_element*nSpace+
			      k*nDOF_trial_element*nSpace+
			      j*nSpace+I]);
}
  
inline
void backwardEuler_c(const double& dt, const double& m_old, const double& m, const double& dm, double& mt, double& dmt)
{  
  mt =(m-m_old)/dt;
  dmt = dm/dt;
}

inline
void backwardEuler(int eN, int k, double dt, double* m_old, double* m, double* dm, double* mt, double* dmt)
{  
  backwardEuler_c(dt,
		  m_old[eN*nQuadraturePoints_element+
			k],
		  m[eN*nQuadraturePoints_element+
		    k],
		  dm[eN*nQuadraturePoints_element+
		     k],
		  mt[eN*nQuadraturePoints_element+
		     k], 
		  dmt[eN*nQuadraturePoints_element+
		      k]);
}
inline
void bdf_c(const double& alpha, const double& beta, const double& m, const double& dm, double& mt, double& dmt)
{  
  mt = alpha*m + beta;
  dmt = alpha*dm;
}
inline
void bdf(int eN, int k, double alpha, double* beta, double* m, double* dm, double* mt, double* dmt)
{  
  bdf_c(alpha,
	beta[eN*nQuadraturePoints_element+
	     k],
	m[eN*nQuadraturePoints_element+
	  k],
	dm[eN*nQuadraturePoints_element+
	   k],
	mt[eN*nQuadraturePoints_element+
	   k], 
	dmt[eN*nQuadraturePoints_element+
	    k]);
}


inline
void calculateSubgridError_tau_c(const double&  hFactor,
				 const double& elementDiameter,
				 const double& odt,
				 const double& rho,
				 const double& mu, 
				 const double& u,
				 const double& v,
				 const double& w,				 				 
				 double& tau0,
				 double& tau1,
				 double& cfl)
{
  double h,nrm_v,taut,tauc,taud;
  h = hFactor*elementDiameter;
  nrm_v = sqrt(u*u+v*v+w*w);
  cfl = nrm_v/h;

  taut = 2.0*rho*odt;
  tauc = 2.0*rho*nrm_v/h;
  taud = 12.0*mu/(h*h);

  tau0 = 1.0/sqrt(taut*taut + tauc*tauc + taud*taud );  
  tau1 = h*h/tau0;
  
}

inline
void calculateInteriorPenalty   (const double&  hFactor,
				 const double& elementDiameter,
				 const double& mu, 				 				 
				 double& gamma)
{
   double h = hFactor*elementDiameter;

  gamma = 1114.0*mu/h; 
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
void calculateNumericalDiffusion(int eN,
				 int k,
				 double shockCapturingDiffusion,
				 double* elementDiameter,
				 double* strong_residual,
				 double* grad_u,
				 double* numDiff)
{
  calculateNumericalDiffusion_c(shockCapturingDiffusion,
				elementDiameter[eN],
				strong_residual[eN*nQuadraturePoints_element+k],
				&grad_u[eN*nQuadraturePoints_element*nSpace+k*nSpace],
				numDiff[eN*nQuadraturePoints_element+k]);
}

inline
double SubgridError_c(const double& error,
		      const double& Lstar_w_dV)
{
  return -error*Lstar_w_dV;
}

inline
void updateSubgridError(int eN,
			int k,
			int i,
			double* error,
			double* Lstar_w_dV,
			double* weak_residual)
{
  weak_residual[eN*nDOF_test_element + i] +=  SubgridError_c(error[eN*nQuadraturePoints_element + 
								   k],
							     Lstar_w_dV[eN*nQuadraturePoints_element*nDOF_test_element + 
									k*nDOF_test_element + 
									i]);
}

inline
double SubgridErrorJacobian_c(const double& derror,
				  const double& Lstar_w_dV)
{
  return -derror*Lstar_w_dV;
}

inline
void updateSubgridErrorJacobian(int eN,
				int k,
				int j,
				int i,
				double* derror,
				double* Lstar_w_dV,
				double* jacobian_weak_residual)
{
  jacobian_weak_residual[eN*nDOF_test_X_trial_element + 
			 i*nDOF_trial_element + 
			 j] += SubgridErrorJacobian_c(derror[eN*nQuadraturePoints_element*nDOF_trial_element+
							     k*nDOF_trial_element+
							     j],
						      Lstar_w_dV[eN*nQuadraturePoints_element*nDOF_test_element + 
								 k*nDOF_test_element + 
								 i]);
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
void updateNumericalDiffusion(int eN,
			      int k,
			      int i,
			      double* numDiff,
			      double* grad_u,
			      double* grad_w_dV,
			      double* weak_residual)
{
  weak_residual[eN*nDOF_test_element + i] += NumericalDiffusion_c(numDiff[eN*nQuadraturePoints_element + 
									  k],
								  &grad_u[eN*nQuadraturePoints_element*nSpace + 
									  k*nSpace],
								  &grad_w_dV[eN*nQuadraturePoints_element*nDOF_test_element*nSpace + 
									     k*nDOF_test_element*nSpace + 
									     i*nSpace]);
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
void updateNumericalDiffusionJacobian(int eN,
				      int k,
				      int j,
				      int i,
				      double* numDiff,
				      double* grad_v,
				      double* grad_w_dV,
				      double* jacobian_weak_residual)
{
  jacobian_weak_residual[eN*nDOF_test_X_trial_element + 
			 i*nDOF_trial_element + 
			 j] += NumericalDiffusionJacobian_c(numDiff[eN*nQuadraturePoints_element + 
								    k],
							    &grad_v[eN*nQuadraturePoints_element*nDOF_trial_element*nSpace + 
								    k*nDOF_trial_element*nSpace + 
								    j*nSpace],
							    &grad_w_dV[eN*nQuadraturePoints_element*nDOF_test_element*nSpace + 
								       k*nDOF_test_element*nSpace + 
								       i*nSpace]);
}


}//RBLES2P
extern "C"
{
  void calculateResidual_RBLES2P(int nElements_global,
				double alpha_bdf,
				double eps_rho,
				double eps_mu,
				double sigma,
				double rho_0,
				double nu_0,
				double rho_1,
				double nu_1,
				double hFactor,
				double shockCapturingDiffusion,
				int* p_l2g, int* vel_l2g,
				double* elementDiameter,
				double* p_dof, double* u_dof, double* v_dof, double* w_dof,
				double* p_trial, double* vel_trial,
				double* p_grad_trial, double* vel_grad_trial,
				double* p_test_dV, double* vel_test_dV,
				double* p_grad_test_dV, double* vel_grad_test_dV,
				double* vel_Hess_trial,double* vel_Hess_test_dV,
				double* g,
				double* phi,
				double* n,
				double* kappa,
				double* q_mom_u_acc,
				double* q_mom_v_acc,
				double* q_mom_w_acc,
				double* q_mass_adv,
				double* q_mom_u_acc_beta_bdf, double* q_mom_v_acc_beta_bdf, double* q_mom_w_acc_beta_bdf,
				double* q_velocity_last,
				double* q_cfl,
				double* q_numDiff_u, double* q_numDiff_v, double* q_numDiff_w,
				double* q_numDiff_u_last, double* q_numDiff_v_last, double* q_numDiff_w_last,
				double* q_elementResidual_p, double* q_elementResidual_u, double* q_elementResidual_v, double* q_elementResidual_w,
				int* sdInfo_u_u_rowptr,int* sdInfo_u_u_colind,			      
				int* sdInfo_u_v_rowptr,int* sdInfo_u_v_colind,
				int* sdInfo_u_w_rowptr,int* sdInfo_u_w_colind,
				int* sdInfo_v_v_rowptr,int* sdInfo_v_v_colind,
				int* sdInfo_v_u_rowptr,int* sdInfo_v_u_colind,
				int* sdInfo_v_w_rowptr,int* sdInfo_v_w_colind,
				int* sdInfo_w_w_rowptr,int* sdInfo_w_w_colind,
				int* sdInfo_w_u_rowptr,int* sdInfo_w_u_colind,
				int* sdInfo_w_v_rowptr,int* sdInfo_w_v_colind,
				int offset_p, int offset_u, int offset_v, int offset_w, int stride_p, int stride_u, int stride_v, int stride_w, double* globalResidual,
				int nExteriorElementBoundaries_global,
				int* exteriorElementBoundariesArray,
				int* elementBoundaryElementsArray,
				int* elementBoundaryLocalElementBoundariesArray,
				double* p_trial_ext,
				double* vel_trial_ext,
				double* p_grad_trial_ext,
				double* vel_grad_trial_ext,
				double* ebqe_phi_ext,
				double* ebqe_n_ext,
				double* ebqe_kappa_ext,
				int* isDOFBoundary_p,
				int* isDOFBoundary_u,
				int* isDOFBoundary_v,
				int* isDOFBoundary_w,
				int* isAdvectiveFluxBoundary_p,
				int* isAdvectiveFluxBoundary_u,
				int* isAdvectiveFluxBoundary_v,
				int* isAdvectiveFluxBoundary_w,
				int* isDiffusiveFluxBoundary_u,
				int* isDiffusiveFluxBoundary_v,
				int* isDiffusiveFluxBoundary_w,
				double* ebqe_bc_p_ext,
				double* ebqe_bc_flux_mass_ext,
				double* ebqe_bc_flux_mom_u_adv_ext,
				double* ebqe_bc_flux_mom_v_adv_ext,
				double* ebqe_bc_flux_mom_w_adv_ext,
				double* ebqe_bc_u_ext,
				double* ebqe_bc_flux_u_diff_ext,
				double* ebqe_penalty_ext,
				double* ebqe_bc_v_ext,
				double* ebqe_bc_flux_v_diff_ext,
				double* ebqe_bc_w_ext,
				double* ebqe_bc_flux_w_diff_ext,
				double* p_test_dS_ext,
				double* vel_test_dS_ext,
				double* q_velocity,
				double* ebqe_velocity_ext,
				double* flux);

  void calculateJacobian_RBLES2P(int nElements_global,
				double alpha_bdf,
				double eps_rho,
				double eps_mu,
				double sigma,
				double rho_0,
				double nu_0,
				double rho_1,
				double nu_1,
				double hFactor,
				double shockCapturingDiffusion,
				int* p_l2g, int* vel_l2g,
				double* elementDiameter,
				double* p_dof, double* u_dof, double* v_dof, double* w_dof,
				double* p_trial, double* vel_trial,
				double* p_grad_trial, double* vel_grad_trial,
				double* p_test_dV, double* vel_test_dV,
				double* p_grad_test_dV, double* vel_grad_test_dV,
				double* vel_Hess_trial,double* vel_Hess_test_dV,
				double* g,
				double* phi,
				double* n,
				double* kappa,
				double* q_mom_u_acc_beta_bdf, double* q_mom_v_acc_beta_bdf, double* q_mom_w_acc_beta_bdf,
				double* q_velocity_last,
				double* q_cfl,
				double* q_numDiff_u_last, double* q_numDiff_v_last, double* q_numDiff_w_last,
				int* sdInfo_u_u_rowptr,int* sdInfo_u_u_colind,			      
				int* sdInfo_u_v_rowptr,int* sdInfo_u_v_colind,
				int* sdInfo_u_w_rowptr,int* sdInfo_u_w_colind,
				int* sdInfo_v_v_rowptr,int* sdInfo_v_v_colind,
				int* sdInfo_v_u_rowptr,int* sdInfo_v_u_colind,
				int* sdInfo_v_w_rowptr,int* sdInfo_v_w_colind,
				int* sdInfo_w_w_rowptr,int* sdInfo_w_w_colind,
				int* sdInfo_w_u_rowptr,int* sdInfo_w_u_colind,
				int* sdInfo_w_v_rowptr,int* sdInfo_w_v_colind,
				int* csrRowIndeces_p_p,int* csrColumnOffsets_p_p,
				int* csrRowIndeces_p_u,int* csrColumnOffsets_p_u,
				int* csrRowIndeces_p_v,int* csrColumnOffsets_p_v,
				int* csrRowIndeces_p_w,int* csrColumnOffsets_p_w,
				int* csrRowIndeces_u_p,int* csrColumnOffsets_u_p,
				int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
				int* csrRowIndeces_u_v,int* csrColumnOffsets_u_v,
				int* csrRowIndeces_u_w,int* csrColumnOffsets_u_w,
				int* csrRowIndeces_v_p,int* csrColumnOffsets_v_p,
				int* csrRowIndeces_v_u,int* csrColumnOffsets_v_u,
				int* csrRowIndeces_v_v,int* csrColumnOffsets_v_v,
				int* csrRowIndeces_v_w,int* csrColumnOffsets_v_w,
				int* csrRowIndeces_w_p,int* csrColumnOffsets_w_p,
				int* csrRowIndeces_w_u,int* csrColumnOffsets_w_u,
				int* csrRowIndeces_w_v,int* csrColumnOffsets_w_v,
				int* csrRowIndeces_w_w,int* csrColumnOffsets_w_w,
				double* globalJacobian,
				int nExteriorElementBoundaries_global,
				int* exteriorElementBoundariesArray,
				int* elementBoundaryElementsArray,
				int* elementBoundaryLocalElementBoundariesArray,
				double* p_trial_ext,
				double* vel_trial_ext,
				double* p_grad_trial_ext,
				double* vel_grad_trial_ext,
				double* ebqe_phi_ext,
				double* ebqe_n_ext,
				double* ebqe_kappa_ext,
				int* isDOFBoundary_p,
				int* isDOFBoundary_u,
				int* isDOFBoundary_v,
				int* isDOFBoundary_w,
				int* isAdvectiveFluxBoundary_p,
				int* isAdvectiveFluxBoundary_u,
				int* isAdvectiveFluxBoundary_v,
				int* isAdvectiveFluxBoundary_w,
				int* isDiffusiveFluxBoundary_u,
				int* isDiffusiveFluxBoundary_v,
				int* isDiffusiveFluxBoundary_w,
				double* ebqe_bc_p_ext,
				double* ebqe_bc_flux_mass_ext,
				double* ebqe_bc_flux_mom_u_adv_ext,
				double* ebqe_bc_flux_mom_v_adv_ext,
				double* ebqe_bc_flux_mom_w_adv_ext,
				double* ebqe_bc_u_ext,
				double* ebqe_bc_flux_u_diff_ext,
				double* ebqe_penalty_ext,
				double* ebqe_bc_v_ext,
				double* ebqe_bc_flux_v_diff_ext,
				double* ebqe_bc_w_ext,
				double* ebqe_bc_flux_w_diff_ext,
				double* p_test_dS_ext,
				double* vel_test_dS_ext,
				int* csrColumnOffsets_eb_p_p,
				int* csrColumnOffsets_eb_p_u,
				int* csrColumnOffsets_eb_p_v,
				int* csrColumnOffsets_eb_p_w,
				int* csrColumnOffsets_eb_u_p,
				int* csrColumnOffsets_eb_u_u,
				int* csrColumnOffsets_eb_u_v,
				int* csrColumnOffsets_eb_u_w,
				int* csrColumnOffsets_eb_v_p,
				int* csrColumnOffsets_eb_v_u,
				int* csrColumnOffsets_eb_v_v,
				int* csrColumnOffsets_eb_v_w,
				int* csrColumnOffsets_eb_w_p,
				int* csrColumnOffsets_eb_w_u,
				int* csrColumnOffsets_eb_w_v,
				int* csrColumnOffsets_eb_w_w);
  void calculateVelocityAverage_RBLES2P(int nExteriorElementBoundaries_global,
				       int* exteriorElementBoundariesArray,
				       int nInteriorElementBoundaries_global,
				       int* interiorElementBoundariesArray,
				       int* elementBoundaryElementsArray,
				       int* elementBoundaryLocalElementBoundariesArray,
				       int* vel_l2g, 
				       double* u_dof, double* v_dof, double* w_dof,
				       double* vel_trial,
				       double* ebqe_velocity,
				       double* velocityAverage);
}//extern "C"
#endif

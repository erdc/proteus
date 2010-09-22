#ifndef VOLUME_AVERAGED_VOF_H
#define VOLUME_AVERAGED_VOF_H
#include <cmath>
#include <iostream>

#define VOL_VOF_SPACE_DIM 3
#if (VOL_VOF_SPACE_DIM == 3)
 #define nSpace 3
 #define nQuadraturePoints_element 5
 #define nDOF_trial_element 4
 #define nDOF_test_element 4
 #define nDOF_test_X_trial_element 16
 #define nQuadraturePoints_elementBoundary 4
 #define VOL_VOF_NAME VolumeAveragedVOF3D
 #define VOL_VOF_RES calculateResidual_VolumeAveragedVOF3D
 #define VOL_VOF_JAC calculateJacobian_VolumeAveragedVOF3D
#else
 #if (VOL_VOF_SPACE_DIM == 2)
  #define nSpace 2
  #define nQuadraturePoints_element 6
  #define nDOF_trial_element 3
  #define nDOF_test_element 3
  #define nDOF_test_X_trial_element 9
  #define nQuadraturePoints_elementBoundary 4
  #define VOL_VOF_NAME VolumeAveragedVOF2D
  #define VOL_VOF_RES calculateResidual_VolumeAveragedVOF2D
  #define VOL_VOF_JAC calculateJacobian_VolumeAveragedVOF2D
 #else 
  #define nSpace 1
  #define nQuadraturePoints_element 3
  #define nDOF_trial_element 2
  #define nDOF_test_element 2
  #define nDOF_test_X_trial_element 4
  #define nQuadraturePoints_elementBoundary 1
  #define VOL_VOF_NAME VolumeAveragedVOF1D
  #define VOL_VOF_RES calculateResidual_VolumeAveragedVOF1D
  #define VOL_VOF_JAC calculateJacobian_VolumeAveragedVOF1D
 #endif
#endif

namespace VOL_VOF_NAME
{
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
    return H;
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
void evaluateCoefficients_c(const double v[nSpace],
			    const double& porosity,
			    const double& u,
			    double& m,
			    double& dm,
			    double f[nSpace],
			    double df[nSpace])
{
  m = porosity*u;
  dm = porosity;
  for (int I=0; I < nSpace; I++)
    {
      f[I] = porosity*v[I]*u;
      df[I] = porosity*v[I];
    }
}

  
inline
void evaluateCoefficients(int K,
			  const double *v,
			  const double *porosity,
			  const double *u,
			  double* m,
			  double* dm,
			  double* f,
			  double* df)
{
  evaluateCoefficients_c(&v[K*nSpace],
			 porosity[K],
                         u[K],
                         m[K],
                         dm[K],
                         &f[K*nSpace],
                         &df[K*nSpace]);
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
   mt =alpha*m + beta;
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
double Mass_weak_c(const double& mt,
		   const double& w_dV)
{
  return mt*w_dV;
}

inline
void updateMass_weak(int eN,
                     int k,
                     int i,
                     double* mt,
                     double* w_dV,
                     double* weak_residual)
{
  weak_residual[eN*nDOF_test_element + 
		i]  += Mass_weak_c(mt[eN*nQuadraturePoints_element + 
				      k],
				   w_dV[eN*nQuadraturePoints_element*nDOF_test_element + 
					k*nDOF_test_element + 
					i]);
}

inline
double MassJacobian_weak_c(const double& dmt,
			   const double& v,
			   const double& w_dV)
{
  return dmt*v*w_dV;
}

inline
void updateMassJacobian_weak(int eN,
			     int k,
			     int j,
			     int i,
			     double* dmt,
			     double* v,
			     double* w_dV,
			     double* jacobian_weak_residual)
{
  jacobian_weak_residual[eN*nDOF_test_X_trial_element + 
			 i*nDOF_trial_element+
			 j] += MassJacobian_weak_c(dmt[eN*nQuadraturePoints_element + 
						       k],
						   v[eN*nQuadraturePoints_element*nDOF_trial_element + 
						     k*nDOF_trial_element + 
						     j],
						   w_dV[eN*nQuadraturePoints_element*nDOF_test_element + 
							k*nDOF_test_element + 
							i]);
}

inline
double Mass_strong_c(const double& mt)
{
  return mt; 
}

inline
void updateMass_strong(int eN,
		       int k,
		       double* mt,
		       double* strong_residual)
{
  strong_residual[eN*nQuadraturePoints_element+
		  k] += Mass_strong_c(mt[eN*nQuadraturePoints_element+
					 k]);
}

inline
double MassJacobian_strong_c(const double& dmt,
			     const double& v)
{
  return dmt*v;
}

inline
void updateMassJacobian_strong(int eN,
                               int k,
                               int j,
                               double* dmt,
                               double* v,
                               double* dstrong_residual)
{
  dstrong_residual[eN*nQuadraturePoints_element*nDOF_trial_element+
		   k*nDOF_trial_element + 
		   j] += MassJacobian_strong_c(dmt[eN*nQuadraturePoints_element+
						   k],
					       v[eN*nQuadraturePoints_element*nDOF_trial_element+
						 k*nDOF_trial_element + 
						 j]);
}

inline
double Mass_adjoint_c(const double& dmt,
		      const double& w_dV)
{
  return dmt*w_dV;
}

inline
void updateMass_adjoint(int eN,
			int k,
			int i,
			double* dmt,
			double* w_dV,
			double* Lstar_w_dV)
{
  Lstar_w_dV[eN*nQuadraturePoints_element*nDOF_test_element + 
	     k*nDOF_test_element + 
	     i] += Mass_adjoint_c(dmt[eN*nQuadraturePoints_element + 
				      k],
				  w_dV[eN*nQuadraturePoints_element*nDOF_test_element + 
				       k*nDOF_test_element + 
				       i]);
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
void updateAdvection_weak(int eN,
                          int k,
                          int i,
                          double* f,
                          double* grad_w_dV,
                          double* weak_residual)
{
  weak_residual[eN*nDOF_test_element + 
		i] += Advection_weak_c(&f[eN*nQuadraturePoints_element*nSpace + 
					  k*nSpace],
				       &grad_w_dV[eN*nQuadraturePoints_element*nDOF_test_element*nSpace + 
						  k*nDOF_test_element*nSpace + 
						  i*nSpace]);
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
void updateAdvectionJacobian_weak(int eN,
				  int k,
				  int j,
				  int i,
				  double* df,
				  double* v,
				  double* grad_w_dV,
				  double* jacobian_weak_residual)
{
  jacobian_weak_residual[eN*nDOF_test_X_trial_element + 
			 i*nDOF_trial_element+
			 j] += AdvectionJacobian_weak_c(&df[eN*nQuadraturePoints_element*nSpace + 
							    k*nSpace],
							v[eN*nQuadraturePoints_element*nDOF_trial_element + 
							  k*nDOF_trial_element + 
							  j],
							&grad_w_dV[eN*nQuadraturePoints_element*nDOF_test_element*nSpace + 
								   k*nDOF_test_element*nSpace + 
								   i*nSpace]);
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
void updateAdvection_strong(int eN,
			    int k,
			    double* df,
			    double* grad_u,
			    double* strong_residual)
{
  strong_residual[eN*nQuadraturePoints_element+
		  k] += Advection_strong_c(&df[eN*nQuadraturePoints_element*nSpace + 
					       k*nSpace],
					   &grad_u[eN*nQuadraturePoints_element*nSpace + 
						   k*nSpace]);
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
void updateAdvectionJacobian_strong(int eN,
                                    int k,
                                    int j,
                                    double* df,
                                    double* grad_v,
                                    double* dstrong_residual)
{
  dstrong_residual[eN*nQuadraturePoints_element*nDOF_trial_element+
		   k*nDOF_trial_element + 
		   j] += AdvectionJacobian_strong_c(&df[eN*nQuadraturePoints_element*nSpace + 
							k*nSpace],
						    &grad_v[eN*nQuadraturePoints_element*nSpace*nDOF_trial_element + 
							    k*nSpace*nDOF_trial_element +
							    j*nSpace]);
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
void updateAdvection_adjoint(int eN,
			     int k,
			     int i,
			     double* df,
			     double* grad_w_dV,
			     double* Lstar_w_dV)
{
  Lstar_w_dV[eN*nQuadraturePoints_element*nDOF_test_element + 
	     k*nDOF_test_element + 
	     i] += Advection_adjoint_c(&df[eN*nQuadraturePoints_element*nSpace + 
					   k*nSpace],
				       &grad_w_dV[eN*nQuadraturePoints_element*nDOF_test_element*nSpace + 
						  k*nDOF_test_element*nSpace + 
						  i*nSpace]);
}

inline
void calculateSubgridError_tau_c(const double& elementDiameter,
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
void calculateSubgridError_tau(int eN,
                               int k,
                               double* elementDiameter,
                               double* dmt,
                               double* dH,
                               double* cfl,
                               double* tau)
{
  calculateSubgridError_tau_c(elementDiameter[eN],
			      dmt[eN*nQuadraturePoints_element+k],
			      &dH[eN*nQuadraturePoints_element*nSpace+
				k*nSpace],
			      cfl[eN*nQuadraturePoints_element+k],
			      tau[eN*nQuadraturePoints_element+k]);
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
void calculateSubgridError_tauRes(int eN,
                                  int k,
                                  double* tau,
                                  double* pdeResidual,
                                  double* dpdeResidual_du,
				  double* subgridError,
				  double* dsubgridError_du)
{
  calculateSubgridError_tauRes_c(tau[eN*nQuadraturePoints_element+k],
				 pdeResidual[eN*nQuadraturePoints_element+k],
				 subgridError[eN*nQuadraturePoints_element+k]);
  calculateSubgridErrorDerivatives_tauRes_c(tau[eN*nQuadraturePoints_element+k],
					    &dpdeResidual_du[eN*nQuadraturePoints_element*nDOF_trial_element+k*nDOF_trial_element],
					    &dsubgridError_du[eN*nQuadraturePoints_element*nDOF_trial_element+k*nDOF_trial_element]);
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
      std::cout<<"warning: open boundary with no external trace, setting to zero if inflow"<<std::endl;
      if (flow >= 0.0)
	{
	  flux = u*flow;
	}
      else
	{
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
void updateExteriorElementBoundaryFlux(int ebNE,
				       int eN_global,
				       int k,
                                       int i,
                                       double* flux,
                                       double* w_dS,
                                       double* residual)
{
  residual[eN_global*nDOF_test_element+
	   i] +=  ExteriorElementBoundaryFlux_c(flux[ebNE*nQuadraturePoints_elementBoundary+
						     k],
						w_dS[ebNE*nQuadraturePoints_elementBoundary*nDOF_test_element+
						     k*nDOF_test_element+
						     i]);
}

inline
double ExteriorNumericalAdvectiveFluxJacobian_c(const double& dflux_left,
						const double& v)
{
  return dflux_left*v;
}

inline
void updateExteriorNumericalAdvectiveFluxJacobian(int ebNE,
						  int k,
						  int j,
						  double* dflux_left,
						  double* v,
						  double* fluxJacobian)
{
  fluxJacobian[ebNE*nQuadraturePoints_elementBoundary*nDOF_trial_element+k*nDOF_trial_element+j] 
    += 
    ExteriorNumericalAdvectiveFluxJacobian_c(dflux_left[ebNE*nQuadraturePoints_elementBoundary+k],
					     v[ebNE*nQuadraturePoints_elementBoundary*nDOF_trial_element+k*nDOF_trial_element+j]);
}

}//VOF namespace
extern "C"
{
  //void calculateResidual_VolumeAveragedVOF(int nElements_global,
  void VOL_VOF_RES(double epsFact,
		   int nElements_global,
		   double alphaBDF,
		   int lag_shockCapturing, /*mwf not used yet*/
		   double shockCapturingDiffusion,
		   int* u_l2g, 
		   double* elementDiameter,
		   double* u_dof,
		   double* u_trial, 
		   double* u_grad_trial, 
		   double* u_test_dV, 
		   double* u_grad_test_dV, 
		   double* velocity,
		   double* porosity,
		   double* q_m,
		   double* q_u,
		   double* q_m_betaBDF,
		   double* q_cfl,
		   double* q_numDiff_u, 
		   double* q_numDiff_u_last, 
		   double* q_elementResidual_u, 
		   int offset_u, int stride_u, 
		   double* globalResidual,
		   int nExteriorElementBoundaries_global,
		   int* exteriorElementBoundariesArray,
		   int* elementBoundaryElementsArray,
		   int* elementBoundaryLocalElementBoundariesArray,
		   double* u_trial_ext,
		   double* u_grad_trial_ext,
		   double* ebqe_velocity_ext,
		   double* ebqe_porosity_ext,
		   double* ebqe_n_ext,
		   int* isDOFBoundary_u,
		   double* ebqe_bc_u_ext,
		   int* isFluxBoundary_u,
		   double* ebqe_bc_flux_u_ext,
		   double* u_test_dS_ext,
		   double* ebqe_phi,
		   double* ebqe_u,
		   double* ebqe_flux);
  //  void calculateJacobian_VolumeAveragedVOF(int nElements_global,
  void VOL_VOF_JAC(int nElements_global,
		   double alphaBDF,
		   int lag_shockCapturing,/*mwf not used yet*/
		   double shockCapturingDiffusion,
		   int* u_l2g,
		   double* elementDiameter,
		   double* u_dof, 
		   double* u_trial, 
		   double* u_grad_trial, 
		   double* u_test_dV, 
		   double* u_grad_test_dV, 
		   double* velocity,
		   double* porosity,
		   double* q_m_betaBDF, 
		   double* q_cfl,
		   double* q_numDiff_u_last, 
		   int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
		   double* globalJacobian,
		   int nExteriorElementBoundaries_global,
		   int* exteriorElementBoundariesArray,
		   int* elementBoundaryElementsArray,
		   int* elementBoundaryLocalElementBoundariesArray,
		   double* u_trial_ext,
		   double* u_grad_trial_ext,
		   double* ebqe_velocity_ext,
		   double* ebqe_porosity_ext,
		   double* ebqe_n,
		   int* isDOFBoundary_u,
		   double* ebqe_bc_u_ext,
		   int* isFluxBoundary_u,
		   double* ebqe_bc_flux_u_ext,
		   double* u_test_dS_ext,
		   int* csrColumnOffsets_eb_u_u);
}//extern "C"
#endif

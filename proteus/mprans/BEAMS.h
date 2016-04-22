 #ifndef BEAMS_H
#define BEAMS_H
#include <cmath>
#include <iostream>
#include "CompKernel.h"
#include "ModelFactory.h"
//#define COMPRESSIBLE_FORM
namespace proteus
{
  class BEAMS_base
  {
  public:
    virtual ~BEAMS_base(){}

    virtual void calculateBeams(//element
				double* mesh_trial_ref,
				double* mesh_grad_trial_ref,
				double* mesh_dof,
				double* mesh_velocity_dof,
				double MOVING_DOMAIN,//0 or 1
				int* mesh_l2g,
				double* dV_ref,
				double* p_trial_ref,
				double* p_grad_trial_ref,
				double* p_test_ref,
				double* p_grad_test_ref,
				double* vel_trial_ref,
				double* vel_grad_trial_ref,
				double* vel_test_ref,
				double* vel_grad_test_ref,
				//element boundary
				double* mesh_trial_trace_ref,
				double* mesh_grad_trial_trace_ref,
				double* dS_ref,
				double* p_trial_trace_ref,
				double* p_grad_trial_trace_ref,
				double* p_test_trace_ref,
				double* p_grad_test_trace_ref,
				double* vel_trial_trace_ref,
				double* vel_grad_trial_trace_ref,
				double* vel_test_trace_ref,
				double* vel_grad_test_trace_ref,					 
				double* normal_ref,
				double* boundaryJac_ref,
				//physics
				double eb_adjoint_sigma,
				double* elementDiameter,
				double* nodeDiametersArray,
				double hFactor,
				int nElements_global,
				int nElementBoundaries_owned,
				double useRBLES,
				double useMetrics, 
				double alphaBDF,
				double epsFact_rho,
				double epsFact_mu, 
				double sigma,
				double rho_0,
				double nu_0,
				double rho_1,
				double nu_1,
				double smagorinskyConstant,
				int turbulenceClosureModel,
				double Ct_sge,
				double Cd_sge,
				double C_dc,
				double C_b,
				//VRANS
				const double* eps_solid,
				const double* phi_solid,
				const double* q_velocity_solid,
				const double* q_porosity,
				const double* q_dragAlpha,
				const double* q_dragBeta,
				const double* q_mass_source,
				const double* q_turb_var_0,
				const double* q_turb_var_1,
				const double* q_turb_var_grad_0,
				int* p_l2g, 
				int* vel_l2g, 
				double* p_dof, 
				double* u_dof, 
				double* v_dof, 
				double* w_dof,
				double* g,
				const double useVF,
				double* vf,
				double* phi,
				double* normal_phi,
				double* kappa_phi,
				double* q_mom_u_acc,
				double* q_mom_v_acc,
				double* q_mom_w_acc,
				double* q_mass_adv,
				double* q_mom_u_acc_beta_bdf, double* q_mom_v_acc_beta_bdf, double* q_mom_w_acc_beta_bdf,
				double* q_velocity_sge,
				double* q_cfl,
				double* q_numDiff_u, double* q_numDiff_v, double* q_numDiff_w,
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
				int offset_p, int offset_u, int offset_v, int offset_w, 
				int stride_p, int stride_u, int stride_v, int stride_w, 
				double* globalResidual,
				int nExteriorElementBoundaries_global,
				int* exteriorElementBoundariesArray,
				int* elementBoundaryElementsArray,
				int* elementBoundaryLocalElementBoundariesArray,
				double* ebqe_vf_ext,
				double* bc_ebqe_vf_ext,
				double* ebqe_phi_ext,
				double* bc_ebqe_phi_ext,
				double* ebqe_normal_phi_ext,
				double* ebqe_kappa_phi_ext,
				//VRANS
				const double* ebqe_porosity_ext,
				const double* ebqe_turb_var_0,
				const double* ebqe_turb_var_1,
				//VRANS end
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
				double* q_x,
				double* q_velocity,
				double* ebqe_velocity,
				double* flux,
				double* elementResidual_p,
				int* boundaryFlags,
				double* barycenters,
				double* wettedAreas,
				double* netForces_p,
				double* netForces_v,
				double* netMoments,
				double* q_dragBeam1,
				double* q_dragBeam2,
				double* q_dragBeam3,
				double* ebqe_dragBeam1,
				double* ebqe_dragBeam2,
				double* ebqe_dragBeam3,
				int nBeams,
				int nBeamElements,
				int beam_quadOrder,
				double beam_Cd,
				double* beamRadius,
				double* xq,
				double* yq,
				double* zq,
				double* Beam_h,
				double* dV_beam,
				double* q1,
				double* q2,
				double* q3,
				double* vel_avg,
				double* netBeamDrag,
				int* beamIsLocal)=0;
     
  };
  
  template<class CompKernelType,
	   int nSpace,
	   int nQuadraturePoints_element,
	   int nDOF_mesh_trial_element,
	   int nDOF_trial_element,
	   int nDOF_test_element,
	   int nQuadraturePoints_elementBoundary>
  class BEAMS : public BEAMS_base
  {
  public:
    const int nDOF_test_X_trial_element;
    CompKernelType ck;
    BEAMS():
      nDOF_test_X_trial_element(nDOF_test_element*nDOF_trial_element),
      ck()
	{/*	     std::cout<<"Constructing RANS2P<CompKernelTemplate<"
		      <<0<<","
		      <<0<<","
		      <<0<<","
		      <<0<<">,"*/
		    /*  <<nSpaceIn<<","
		      <<nQuadraturePoints_elementIn<<","
		      <<nDOF_mesh_trial_elementIn<<","
		      <<nDOF_trial_elementIn<<","
		      <<nDOF_test_elementIn<<","
		      <<nQuadraturePoints_elementBoundaryIn<<">());"*/
	  /*  <<std::endl<<std::flush; */
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
      return H;
    }
    
    inline double smoothedHeaviside_integral(double eps, double phi)
    {
      double HI;
      if (phi > eps)
	{
	  HI= phi - eps							\
	    + 0.5*(eps + 0.5*eps*eps/eps - eps*cos(M_PI*eps/eps)/(M_PI*M_PI)) \
	    - 0.5*((-eps) + 0.5*(-eps)*(-eps)/eps - eps*cos(M_PI*(-eps)/eps)/(M_PI*M_PI));
	}
      else if (phi < -eps)
	{
	  HI=0.0;
	}
      else
	{
	  HI = 0.5*(phi + 0.5*phi*phi/eps - eps*cos(M_PI*phi/eps)/(M_PI*M_PI)) \
	    - 0.5*((-eps) + 0.5*(-eps)*(-eps)/eps - eps*cos(M_PI*(-eps)/eps)/(M_PI*M_PI));
	}
      return HI;
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
      return d;
    }
    inline
      double delta_h(const double r)
    {
      double delta;
      if (r <= -2.0)
	delta= 0.0;
      else if (r <= -1.0)
	delta= 1.0/8.0*(5.0+2.0*r-sqrt(-7.0-12.0*r-4*r*r));
      else if  (r <= 0.0)
	delta= 1.0/8.0*(3.0+2.0*r+sqrt(1.0-4.0*r-4.0*r*r));
      else if (r <= 1.0)
	delta= 1.0/8.0*(3.0-2.0*r+sqrt(1.0+4.0*r-4.0*r*r));
      else if (r <= 2.0)
	delta= 1.0/8.0*(5.0-2.0*r-sqrt(-7.0+12.0*r-4.0*r*r));
      else
	delta= 0.0;
      return delta;
    }

    inline
      void calculateBeamSinks(const int nBeams,
			      const int nBeamElements,
			      const int beam_quadOrder,
			      const double beam_Cd,
			      const double* beamRadius,
			      const double x,
			      const double y,
			      const double z,
			      const double *xq,
			      const double *yq,
			      const double *zq,
			      const double *Beam_h,
			      const double *dV_beam,
			      const double u,
			      const double v,
			      const double w,
			      const double eps_rho,
			      const double& phi,
			      const double rho_0,
			      const double rho_1,
			      double& mom_u_source,
			      double& mom_v_source,
			      double& mom_w_source,
			      const double dV,
			      double* netBeamDrag,
			      int* beamIsLocal)
    {
      double rho, H_rho, delt, vel;
      H_rho = smoothedHeaviside(eps_rho,phi);
      rho = rho_0*(1.0-H_rho)+rho_1*H_rho;
      vel = sqrt(u*u+v*v+w*w);
      mom_u_source = 0.0;
      mom_v_source = 0.0;
      mom_w_source = 0.0;
      for(int I=0; I<nBeams; I++)
	{
	  if (beamIsLocal[I])
	    {
	  for(int k=0;k<nBeamElements; k++)
	    {
	      for(int l=0;l<beam_quadOrder; l++)
		{
		  //delt= delta_h((x-xq[I*nBeamElements*beam_quadOrder+beam_quadOrder*k+l])/Beam_h[I*nBeamElements+k])*delta_h((y-yq[I*nBeamElements*beam_quadOrder+beam_quadOrder*k+l])/Beam_h[I*nBeamElements+k])*delta_h((z-zq[I*nBeamElements*beam_quadOrder+beam_quadOrder*k+l])/Beam_h[I*nBeamElements+k])/(Beam_h[I*nBeamElements+k]*Beam_h[I*nBeamElements+k]*Beam_h[I*nBeamElements+k]);
		  delt = delta_h((x-xq[I*nBeamElements*beam_quadOrder+beam_quadOrder*k+l])/Beam_h[I*nBeamElements+k]);
		     if (delt>0.0)
		       {
			 delt *= delta_h((y-yq[I*nBeamElements*beam_quadOrder+beam_quadOrder*k+l])/Beam_h[I*nBeamElements+k]);
			 if (delt > 0.0)
			   {
			     delt *= delta_h((z-zq[I*nBeamElements*beam_quadOrder+beam_quadOrder*k+l])/Beam_h[I*nBeamElements+k]);
			     
			     if (delt > 0.0)
			       {
				 delt = delt/(Beam_h[I*nBeamElements+k]*Beam_h[I*nBeamElements+k]*Beam_h[I*nBeamElements+k]);
				 mom_u_source += beam_Cd*rho*beamRadius[I]*u*vel*delt*dV_beam[I*nBeamElements*beam_quadOrder+beam_quadOrder*k+l];
				 mom_v_source += beam_Cd*rho*beamRadius[I]*v*vel*delt*dV_beam[I*nBeamElements*beam_quadOrder+beam_quadOrder*k+l];
				 mom_w_source += beam_Cd*rho*beamRadius[I]*w*vel*delt*dV_beam[I*nBeamElements*beam_quadOrder+beam_quadOrder*k+l];
			       }
			   }
		       }
		  /* q1[I*nBeamElements*beam_quadOrder+beam_quadOrder*k+l] += beam_Cd*rho*beamRadius[I]*u*vel*delt*dV; */
		  /* q2[I*nBeamElements*beam_quadOrder+beam_quadOrder*k+l] += beam_Cd*rho*beamRadius[I]*v*vel*delt*dV; */
		  /* q3[I*nBeamElements*beam_quadOrder+beam_quadOrder*k+l] += beam_Cd*rho*beamRadius[I]*w*vel*delt*dV; */
		}
	    }
	    }
	}
      if ((y>.1) && ( y< 0.4) && (x >= 1.2) && (x <= 11.0))
	{
	  netBeamDrag[0]+= dV*sqrt(mom_u_source*mom_u_source+mom_v_source*mom_v_source+mom_w_source*mom_w_source);
	}
    }
inline
  void calculateBeamLoads(const int nBeams,
			  const int nBeamElements,
			  const int beam_quadOrder,
			  const double beam_Cd,
			  const double* beamRadius,
			  const double x,
			  const double y,
			  const double z,
			  const double *xq,
			  const double *yq,
			  const double *zq,
			  const double *Beam_h,
			  const double u,
			  const double v,
			  const double w,
			  const double eps_rho,
			  const double& phi,
			  const double rho_0,
			  const double rho_1,
			  double *q1,
			  double *q2,
			  double *q3,
			  const double dV,
			  double *vel_avg,
			  int *beamIsLocal)
    {
      double rho, H_rho, delt, vel,h_save, buoy;
      H_rho = smoothedHeaviside(eps_rho,phi);
      rho = rho_0*(1.0-H_rho)+rho_1*H_rho;
      vel = sqrt(u*u+v*v+w*w);
      if ((x >= 1.2) && (x <= 11.0) && (y >= 0.1) && (y<= 0.4))
	{
	  vel_avg[0] += rho*u*dV;
	  vel_avg[1] += rho*v*dV;
	  vel_avg[2] += rho*w*dV;
	}
      for(int I=0; I<nBeams; I++)
	{
	  if (beamIsLocal[I])
	    {
	  buoy=(rho - 1350.0)*3.14159*beamRadius[I]*beamRadius[I];
	  for(int k=0;k<nBeamElements; k++)
	    {
	      for(int l=0;l<beam_quadOrder; l++)
		{
		  if (k==0)
		    h_save=0.25*Beam_h[I*nBeamElements+k];
		  else if (k==1)
		    h_save = 0.5*Beam_h[I*nBeamElements+k];
		  else
		    h_save = Beam_h[I*nBeamElements+k];
		  
		  delt= delta_h((x-xq[I*nBeamElements*beam_quadOrder+beam_quadOrder*k+l])/h_save)*delta_h((y-yq[I*nBeamElements*beam_quadOrder+beam_quadOrder*k+l])/h_save)*delta_h((z-zq[I*nBeamElements*beam_quadOrder+beam_quadOrder*k+l])/h_save)/(h_save*h_save*h_save);
		  /* mom_u_source += beam_Cd*rho*beamRadius[I]*u*vel*delt*dV_beam[I*nBeamElements*beam_quadOrder+beam_quadOrder*k+l]; */
		  /* mom_v_source += beam_Cd*rho*beamRadius[I]*v*vel*delt*dV_beam[I*nBeamElements*beam_quadOrder+beam_quadOrder*k+l]; */
		  /* mom_w_source += beam_Cd*rho*beamRadius[I]*w*vel*delt*dV_beam[I*nBeamElements*beam_quadOrder+beam_quadOrder*k+l]; */
		  if (delt > 0.0)
		    {
		      q1[I*nBeamElements*beam_quadOrder+beam_quadOrder*k+l] += beam_Cd*rho*beamRadius[I]*u*vel*delt*dV;
		      q2[I*nBeamElements*beam_quadOrder+beam_quadOrder*k+l] += beam_Cd*rho*beamRadius[I]*v*vel*delt*dV;
		      q3[I*nBeamElements*beam_quadOrder+beam_quadOrder*k+l] += (beam_Cd*rho*beamRadius[I]*w*vel+buoy)*delt*dV;
		    }
		}
	    }}
	}
    }
   


    void calculateBeams(//element
		       double* mesh_trial_ref,
		       double* mesh_grad_trial_ref,
		       double* mesh_dof,
		       double* mesh_velocity_dof,
		       double MOVING_DOMAIN,
		       int* mesh_l2g,
		       double* dV_ref,
		       double* p_trial_ref,
		       double* p_grad_trial_ref,
		       double* p_test_ref,
		       double* p_grad_test_ref,
		       double* vel_trial_ref,
		       double* vel_grad_trial_ref,
		       double* vel_test_ref,
		       double* vel_grad_test_ref,
		       //element boundary
		       double* mesh_trial_trace_ref,
		       double* mesh_grad_trial_trace_ref,
		       double* dS_ref,
		       double* p_trial_trace_ref,
		       double* p_grad_trial_trace_ref,
		       double* p_test_trace_ref,
		       double* p_grad_test_trace_ref,
		       double* vel_trial_trace_ref,
		       double* vel_grad_trial_trace_ref,
		       double* vel_test_trace_ref,
		       double* vel_grad_test_trace_ref,					 
		       double* normal_ref,
		       double* boundaryJac_ref,
		       //physics
		       double eb_adjoint_sigma,
		       double* elementDiameter,
		       double* nodeDiametersArray,
		       double hFactor,
		       int nElements_global,
		       int nElementBoundaries_owned,
		       double useRBLES,
		       double useMetrics, 
		       double alphaBDF,
		       double epsFact_rho,
		       double epsFact_mu, 
		       double sigma,
		       double rho_0,
		       double nu_0,
		       double rho_1,
		       double nu_1,
		       double smagorinskyConstant,
		       int turbulenceClosureModel,
		       double Ct_sge,
		       double Cd_sge,
		       double C_dc,
		       double C_b,
		       //VRANS
		       const double* eps_solid,
		       const double* phi_solid,
		       const double* q_velocity_solid,
		       const double* q_porosity,
		       const double* q_dragAlpha,
		       const double* q_dragBeta,
		       const double* q_mass_source,
		       const double* q_turb_var_0,
		       const double* q_turb_var_1,
		       const double* q_turb_var_grad_0,
		       //
		       int* p_l2g, 
		       int* vel_l2g, 
		       double* p_dof, 
		       double* u_dof, 
		       double* v_dof, 
		       double* w_dof,
		       double* g,
		       const double useVF,
		       double* vf,
		       double* phi,
		       double* normal_phi,
		       double* kappa_phi,
		       double* q_mom_u_acc,
		       double* q_mom_v_acc,
		       double* q_mom_w_acc,
		       double* q_mass_adv,
		       double* q_mom_u_acc_beta_bdf, double* q_mom_v_acc_beta_bdf, double* q_mom_w_acc_beta_bdf,
		       double* q_velocity_sge,
		       double* q_cfl,
		       double* q_numDiff_u, double* q_numDiff_v, double* q_numDiff_w,
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
		       int offset_p, int offset_u, int offset_v, int offset_w, 
		       int stride_p, int stride_u, int stride_v, int stride_w, 
		       double* globalResidual,
		       int nExteriorElementBoundaries_global,
		       int* exteriorElementBoundariesArray,
		       int* elementBoundaryElementsArray,
		       int* elementBoundaryLocalElementBoundariesArray,
		       double* ebqe_vf_ext,
		       double* bc_ebqe_vf_ext,
		       double* ebqe_phi_ext,
		       double* bc_ebqe_phi_ext,
		       double* ebqe_normal_phi_ext,
		       double* ebqe_kappa_phi_ext,
		       //VRANS
		       const double* ebqe_porosity_ext,
		       const double* ebqe_turb_var_0,
		       const double* ebqe_turb_var_1,
		       //VRANS end
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
		       double* q_x,
		       double* q_velocity,
		       double* ebqe_velocity,
		       double* flux,
		       double* elementResidual_p_save,
		       int* boundaryFlags,
		       double* barycenters,
		       double* wettedAreas,
		       double* netForces_p,
		       double* netForces_v,
		       double* netMoments,
		       double* q_dragBeam1,
		       double* q_dragBeam2,
		       double* q_dragBeam3,
		       double* ebqe_dragBeam1,
		       double* ebqe_dragBeam2,
		       double* ebqe_dragBeam3,
		       int nBeams,
		       int nBeamElements,
		       int beam_quadOrder,
		       double beam_Cd,
		       double* beamRadius,
		       double* xq,
		       double* yq,
		       double* zq,
		       double* Beam_h,
		       double* dV_beam,
		       double* q1,
		       double* q2,
		       double* q3,
		       double* vel_avg,
		       double* netBeamDrag,
		       int* beamIsLocal)
    {
      //
      //loop over elements to compute volume integrals and load them into element and global residual
      //
      double globalConservationError=0.0;
      for(int eN=0;eN<nElements_global;eN++)
	{
	  register double eps_rho,eps_mu;
	  //
	  //loop over quadrature points and compute integrands
	  //
	  for(int k=0;k<nQuadraturePoints_element;k++)
	    {
	      //compute indices and declare local storage
	      register int eN_k = eN*nQuadraturePoints_element+k,
		eN_k_nSpace = eN_k*nSpace,
		eN_nDOF_trial_element = eN*nDOF_trial_element;
	      register double p=0.0,u=0.0,v=0.0,w=0.0,
		grad_p[nSpace],grad_u[nSpace],grad_v[nSpace],grad_w[nSpace],
		mom_u_acc=0.0,
		dmom_u_acc_u=0.0,
		mom_v_acc=0.0,
		dmom_v_acc_v=0.0,
		mom_w_acc=0.0,
		dmom_w_acc_w=0.0,
		mass_adv[nSpace],
		dmass_adv_u[nSpace],
		dmass_adv_v[nSpace],
		dmass_adv_w[nSpace],
		mom_u_adv[nSpace],
		dmom_u_adv_u[nSpace],
		dmom_u_adv_v[nSpace],
		dmom_u_adv_w[nSpace],
		mom_v_adv[nSpace],
		dmom_v_adv_u[nSpace],
		dmom_v_adv_v[nSpace],
		dmom_v_adv_w[nSpace],
		mom_w_adv[nSpace],
		dmom_w_adv_u[nSpace],
		dmom_w_adv_v[nSpace],
		dmom_w_adv_w[nSpace],
		mom_uu_diff_ten[nSpace],
		mom_vv_diff_ten[nSpace],
		mom_ww_diff_ten[nSpace],
		mom_uv_diff_ten[1],
		mom_uw_diff_ten[1],
		mom_vu_diff_ten[1],
		mom_vw_diff_ten[1],
		mom_wu_diff_ten[1],
		mom_wv_diff_ten[1],
		mom_u_source=0.0,
		mom_v_source=0.0,
		mom_w_source=0.0,
		mom_u_ham=0.0,
		dmom_u_ham_grad_p[nSpace],
		mom_v_ham=0.0,
		dmom_v_ham_grad_p[nSpace],
		mom_w_ham=0.0,
		dmom_w_ham_grad_p[nSpace],
		mom_u_acc_t=0.0,
		dmom_u_acc_u_t=0.0,
		mom_v_acc_t=0.0,
		dmom_v_acc_v_t=0.0,
		mom_w_acc_t=0.0,
		dmom_w_acc_w_t=0.0,
		pdeResidual_p=0.0,
		pdeResidual_u=0.0,
		pdeResidual_v=0.0,
		pdeResidual_w=0.0,
		Lstar_u_p[nDOF_test_element],
		Lstar_v_p[nDOF_test_element],
		Lstar_w_p[nDOF_test_element],
		Lstar_u_u[nDOF_test_element],
		Lstar_v_v[nDOF_test_element],
		Lstar_w_w[nDOF_test_element],
		Lstar_p_u[nDOF_test_element],
		Lstar_p_v[nDOF_test_element],
		Lstar_p_w[nDOF_test_element],
		subgridError_p=0.0,
		subgridError_u=0.0,
		subgridError_v=0.0,
		subgridError_w=0.0,
		tau_p=0.0,tau_p0=0.0,tau_p1=0.0,
		tau_v=0.0,tau_v0=0.0,tau_v1=0.0,
		jac[nSpace*nSpace],
		jacDet,
		jacInv[nSpace*nSpace],
		p_grad_trial[nDOF_trial_element*nSpace],vel_grad_trial[nDOF_trial_element*nSpace],
		p_test_dV[nDOF_trial_element],vel_test_dV[nDOF_trial_element],
		p_grad_test_dV[nDOF_test_element*nSpace],vel_grad_test_dV[nDOF_test_element*nSpace],
		dV,x,y,z,xt,yt,zt,
		//
		porosity,
		//meanGrainSize,
		mass_source,
		dmom_u_source[nSpace],
		dmom_v_source[nSpace],
		dmom_w_source[nSpace],
		//
		G[nSpace*nSpace],G_dd_G,tr_G,norm_Rv,h_phi, dmom_adv_star[nSpace],dmom_adv_sge[nSpace];
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
	      ck.calculateH_element(eN,
				    k,
				    nodeDiametersArray,
				    mesh_l2g,
				    mesh_trial_ref,
				    h_phi);
	      
	      ck.calculateMappingVelocity_element(eN,
						  k,
						  mesh_velocity_dof,
						  mesh_l2g,
						  mesh_trial_ref,
						  xt,yt,zt);
	      //xt=0.0;yt=0.0;zt=0.0;
	      //std::cout<<"xt "<<xt<<'\t'<<yt<<'\t'<<zt<<std::endl;
	      //get the physical integration weight
	      dV = fabs(jacDet)*dV_ref[k];
	      ck.calculateG(jacInv,G,G_dd_G,tr_G);
	      //ck.calculateGScale(G,&normal_phi[eN_k_nSpace],h_phi);
	      
	      eps_rho = epsFact_rho*(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter[eN]);
	      eps_mu  = epsFact_mu *(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter[eN]);
	     
	      //get the trial function gradients
	      ck.gradTrialFromRef(&vel_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,vel_grad_trial);
	      //get the solution
	      ck.valFromDOF(u_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],u);
	      ck.valFromDOF(v_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],v);
	      ck.valFromDOF(w_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],w);

	      //precalculate test function products with integration weights
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  vel_test_dV[j] = vel_test_ref[k*nDOF_trial_element+j]*dV;
		}
	      //VRANS
	      porosity      = q_porosity[eN_k];
	      //meanGrainSize = q_meanGrain[eN_k]; 
	      //
	      //save velocity at quadrature points for other models to use
	      q_velocity[eN_k_nSpace+0]=u;
	      q_velocity[eN_k_nSpace+1]=v;
	      q_velocity[eN_k_nSpace+2]=w;
	      q_x[eN_k_nSpace+0]=x;
	      q_x[eN_k_nSpace+1]=y;
	      q_x[eN_k_nSpace+2]=z;
	      
	      

	      //
	      //moving mesh
	      //
	     
	      //
	      //calculate time derivative at quadrature points
	      //
	      ck.bdf(alphaBDF,
		     q_mom_u_acc_beta_bdf[eN_k],
		     mom_u_acc,
		     dmom_u_acc_u,
		     mom_u_acc_t,
		     dmom_u_acc_u_t);
	      ck.bdf(alphaBDF,
		     q_mom_v_acc_beta_bdf[eN_k],
		     mom_v_acc,
		     dmom_v_acc_v,
		     mom_v_acc_t,
		     dmom_v_acc_v_t);
	      ck.bdf(alphaBDF,
		     q_mom_w_acc_beta_bdf[eN_k],
		     mom_w_acc,
		     dmom_w_acc_w,
		     mom_w_acc_t,
		     dmom_w_acc_w_t);
	      //
	      calculateBeamSinks(nBeams,
				 nBeamElements,
				 beam_quadOrder,
				 beam_Cd,
				 beamRadius,
				 x,
				 y,
				 z,
				 xq,
				 yq,
				 zq,
				 Beam_h,
				 dV_beam,
				 u,
				 v,
				 w,
				 eps_rho,
				 phi[eN_k],
				 rho_0,
				 rho_1,
				 q_dragBeam1[eN_k],
				 q_dragBeam2[eN_k],
				 q_dragBeam3[eN_k],
				 dV,
				 netBeamDrag,
				 beamIsLocal);
	      calculateBeamLoads(nBeams,
				 nBeamElements,
				 beam_quadOrder,
				 beam_Cd,
				 beamRadius,
				 x,
				 y,
				 z,
				 xq,
				 yq,
				 zq,
				 Beam_h,
				 u,
				 v,
				 w,
				 eps_rho,
				 phi[eN_k],
				 rho_0,
				 rho_1,
				 q1,
				 q2,
				 q3,
				 dV,
				 vel_avg,
				 beamIsLocal);
	      
	      // 
	      
	    }
	  //
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
	  register double eps_rho, eps_mu;
	  for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) 
	    { 
	      register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
		ebNE_kb_nSpace = ebNE_kb*nSpace,
		ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb,
		ebN_local_kb_nSpace = ebN_local_kb*nSpace;
	      register double p_ext=0.0,
		u_ext=0.0,
		v_ext=0.0,
		w_ext=0.0,
		grad_p_ext[nSpace],
		grad_u_ext[nSpace],
		grad_v_ext[nSpace],
		grad_w_ext[nSpace],
		mom_u_acc_ext=0.0,
		dmom_u_acc_u_ext=0.0,
		mom_v_acc_ext=0.0,
		dmom_v_acc_v_ext=0.0,
		mom_w_acc_ext=0.0,
		dmom_w_acc_w_ext=0.0,
		mass_adv_ext[nSpace],
		dmass_adv_u_ext[nSpace],
		dmass_adv_v_ext[nSpace],
		dmass_adv_w_ext[nSpace],
		mom_u_adv_ext[nSpace],
		dmom_u_adv_u_ext[nSpace],
		dmom_u_adv_v_ext[nSpace],
		dmom_u_adv_w_ext[nSpace],
		mom_v_adv_ext[nSpace],
		dmom_v_adv_u_ext[nSpace],
		dmom_v_adv_v_ext[nSpace],
		dmom_v_adv_w_ext[nSpace],
		mom_w_adv_ext[nSpace],
		dmom_w_adv_u_ext[nSpace],
		dmom_w_adv_v_ext[nSpace],
		dmom_w_adv_w_ext[nSpace],
		mom_uu_diff_ten_ext[nSpace],
		mom_vv_diff_ten_ext[nSpace],
		mom_ww_diff_ten_ext[nSpace],
		mom_uv_diff_ten_ext[1],
		mom_uw_diff_ten_ext[1],
		mom_vu_diff_ten_ext[1],
		mom_vw_diff_ten_ext[1],
		mom_wu_diff_ten_ext[1],
		mom_wv_diff_ten_ext[1],
		mom_u_source_ext=0.0,
		mom_v_source_ext=0.0,
		mom_w_source_ext=0.0,
		mom_u_ham_ext=0.0,
		dmom_u_ham_grad_p_ext[nSpace],
		mom_v_ham_ext=0.0,
		dmom_v_ham_grad_p_ext[nSpace],
		mom_w_ham_ext=0.0,
		dmom_w_ham_grad_p_ext[nSpace],
		dmom_u_adv_p_ext[nSpace],
		dmom_v_adv_p_ext[nSpace],
		dmom_w_adv_p_ext[nSpace],
		flux_mass_ext=0.0,
		flux_mom_u_adv_ext=0.0,
		flux_mom_v_adv_ext=0.0,
		flux_mom_w_adv_ext=0.0,
		flux_mom_uu_diff_ext=0.0,
		flux_mom_uv_diff_ext=0.0,
		flux_mom_uw_diff_ext=0.0,
		flux_mom_vu_diff_ext=0.0,
		flux_mom_vv_diff_ext=0.0,
		flux_mom_vw_diff_ext=0.0,
		flux_mom_wu_diff_ext=0.0,
		flux_mom_wv_diff_ext=0.0,
		flux_mom_ww_diff_ext=0.0,
		bc_p_ext=0.0,
		bc_u_ext=0.0,
		bc_v_ext=0.0,
		bc_w_ext=0.0,
		bc_mom_u_acc_ext=0.0,
		bc_dmom_u_acc_u_ext=0.0,
		bc_mom_v_acc_ext=0.0,
		bc_dmom_v_acc_v_ext=0.0,
		bc_mom_w_acc_ext=0.0,
		bc_dmom_w_acc_w_ext=0.0,
		bc_mass_adv_ext[nSpace],
		bc_dmass_adv_u_ext[nSpace],
		bc_dmass_adv_v_ext[nSpace],
		bc_dmass_adv_w_ext[nSpace],
		bc_mom_u_adv_ext[nSpace],
		bc_dmom_u_adv_u_ext[nSpace],
		bc_dmom_u_adv_v_ext[nSpace],
		bc_dmom_u_adv_w_ext[nSpace],
		bc_mom_v_adv_ext[nSpace],
		bc_dmom_v_adv_u_ext[nSpace],
		bc_dmom_v_adv_v_ext[nSpace],
		bc_dmom_v_adv_w_ext[nSpace],
		bc_mom_w_adv_ext[nSpace],
		bc_dmom_w_adv_u_ext[nSpace],
		bc_dmom_w_adv_v_ext[nSpace],
		bc_dmom_w_adv_w_ext[nSpace],
		bc_mom_uu_diff_ten_ext[nSpace],
		bc_mom_vv_diff_ten_ext[nSpace],
		bc_mom_ww_diff_ten_ext[nSpace],
		bc_mom_uv_diff_ten_ext[1],
		bc_mom_uw_diff_ten_ext[1],
		bc_mom_vu_diff_ten_ext[1],
		bc_mom_vw_diff_ten_ext[1],
		bc_mom_wu_diff_ten_ext[1],
		bc_mom_wv_diff_ten_ext[1],
		bc_mom_u_source_ext=0.0,
		bc_mom_v_source_ext=0.0,
		bc_mom_w_source_ext=0.0,
		bc_mom_u_ham_ext=0.0,
		bc_dmom_u_ham_grad_p_ext[nSpace],
		bc_mom_v_ham_ext=0.0,
		bc_dmom_v_ham_grad_p_ext[nSpace],
		bc_mom_w_ham_ext=0.0,
		bc_dmom_w_ham_grad_p_ext[nSpace],
		jac_ext[nSpace*nSpace],
		jacDet_ext,
		jacInv_ext[nSpace*nSpace],
		boundaryJac[nSpace*(nSpace-1)],
		metricTensor[(nSpace-1)*(nSpace-1)],
		metricTensorDetSqrt,
		dS,p_test_dS[nDOF_test_element],vel_test_dS[nDOF_test_element],
		p_grad_trial_trace[nDOF_trial_element*nSpace],vel_grad_trial_trace[nDOF_trial_element*nSpace],
		vel_grad_test_dS[nDOF_trial_element*nSpace],
		normal[3],x_ext,y_ext,z_ext,xt_ext,yt_ext,zt_ext,integralScaling,
		//VRANS
		porosity_ext,
		//
		G[nSpace*nSpace],G_dd_G,tr_G,h_phi,h_penalty,penalty,
		force_x,force_y,force_z,force_p_x,force_p_y,force_p_z,force_v_x,force_v_y,force_v_z,r_x,r_y,r_z, junk[nSpace];
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
	      //xt_ext=0.0;yt_ext=0.0;zt_ext=0.0;
	      //std::cout<<"xt_ext "<<xt_ext<<'\t'<<yt_ext<<'\t'<<zt_ext<<std::endl;
	      //std::cout<<"integralScaling - metricTensorDetSrt ==============================="<<integralScaling-metricTensorDetSqrt<<std::endl;
	      /* std::cout<<"metricTensorDetSqrt "<<metricTensorDetSqrt */
	      /* 	       <<"dS_ref[kb]"<<dS_ref[kb]<<std::endl; */
	      dS = ((1.0-MOVING_DOMAIN)*metricTensorDetSqrt + MOVING_DOMAIN*integralScaling)*dS_ref[kb];
	      //get the metric tensor
	      //cek todo use symmetry
	      ck.calculateG(jacInv_ext,G,G_dd_G,tr_G);
	      ck.calculateGScale(G,&ebqe_normal_phi_ext[ebNE_kb_nSpace],h_phi);
	      
	      eps_rho = epsFact_rho*(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter[eN]);
	      eps_mu  = epsFact_mu *(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter[eN]);
	      
	      //compute shape and solution information
	      //shape
	      ck.gradTrialFromRef(&vel_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,vel_grad_trial_trace);
	      //cek hack use trial ck.gradTrialFromRef(&vel_grad_test_trace_ref[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,vel_grad_test_trace);
	      //solution and gradients	
	      ck.valFromDOF(u_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_trace_ref[ebN_local_kb*nDOF_test_element],u_ext);
	      ck.valFromDOF(v_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_trace_ref[ebN_local_kb*nDOF_test_element],v_ext);
	      ck.valFromDOF(w_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_trace_ref[ebN_local_kb*nDOF_test_element],w_ext);
	      
	      //precalculate test function products with integration weights
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  vel_test_dS[j] = vel_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;
		}

	      bc_u_ext = isDOFBoundary_u[ebNE_kb]*ebqe_bc_u_ext[ebNE_kb]+(1-isDOFBoundary_u[ebNE_kb])*u_ext;
	      bc_v_ext = isDOFBoundary_v[ebNE_kb]*ebqe_bc_v_ext[ebNE_kb]+(1-isDOFBoundary_v[ebNE_kb])*v_ext;
	      bc_w_ext = isDOFBoundary_w[ebNE_kb]*ebqe_bc_w_ext[ebNE_kb]+(1-isDOFBoundary_w[ebNE_kb])*w_ext;
	      //VRANS
	      porosity_ext = ebqe_porosity_ext[ebNE_kb];
	      calculateBeamSinks(nBeams,
				 nBeamElements,
				 beam_quadOrder,
				 beam_Cd,
				 beamRadius,
				 x_ext,
				 y_ext,
				 z_ext,
				 xq,
				 yq,
				 zq,
				 Beam_h,
				 dV_beam,
				 u_ext,
				 v_ext,
				 w_ext,
				 eps_rho,
				 phi[ebNE_kb],
				 rho_0,
				 rho_1,
				 ebqe_dragBeam1[ebNE_kb],
				 ebqe_dragBeam2[ebNE_kb],
				 ebqe_dragBeam3[ebNE_kb],
				 dS,
				 junk,
				 beamIsLocal);
	      //
	      
	    }//kb
	  //
	  //update the element and global residual storage
	  //
	 
	}//ebNE
    }



  };//RANS2P
  
  inline BEAMS_base* newBEAMS(int nSpaceIn,
				int nQuadraturePoints_elementIn,
				int nDOF_mesh_trial_elementIn,
				int nDOF_trial_elementIn,
				int nDOF_test_elementIn,
				int nQuadraturePoints_elementBoundaryIn,
				int CompKernelFlag)
  {
    return proteus::chooseAndAllocateDiscretization<BEAMS_base,BEAMS,CompKernel>(nSpaceIn,
										   nQuadraturePoints_elementIn,
										   nDOF_mesh_trial_elementIn,
										   nDOF_trial_elementIn,
										   nDOF_test_elementIn,
										   nQuadraturePoints_elementBoundaryIn,
										   CompKernelFlag);
  }
}//proteus

#endif

#include <cmath>
#include <iostream>
namespace proteus
{  
  // template<int nSpace>
  const int nSpace(2);
  class cppHsuSedStress
{
public:
 cppHsuSedStress(
		 double aDarcy, // darcy parameter for drag term. Default value from Ergun (1952) is 150
		 double betaForch, // forchheimer parameter for drag term. Default value from Ergun (1952) is 1.75
		 double grain, // Grain size, default assumed as d50
		 double packFraction, //Critical volume fraction for switching the drag relation 0.2 by default, see Chen and Hsu 2014
		 double packMargin, // For packFraction +- packmargin, the drag coefficient is calculated by taking a weighted combination of the two relation (for packed and non-packed sediment
		 double sigmaC,
		 double C3e,
		 double C4e,
		 double eR
): 
 
  aDarcy_(aDarcy), 
    betaForch_(betaForch), 
    grain_(grain),
    packFraction_(packFraction),
    packMargin_(packMargin),
    sigmaC_(sigmaC),  
    C3e_(C3e),
    C4e_(C4e),
    eR_(eR)


         
    {}
  
    inline double  betaCoeff(     
			      double sedF, // Sediment fraction
			      double uFluid[nSpace], //Fluid velocity
			      double uSolid[nSpace], //Sediment velocity
			      double nu //Kinematic viscosity
			      )
  {
    double du2 = 0.;
    for (int ii=0; ii<nSpace;  ii++)
      {
	du2+= (uFluid[ii] - uSolid[ii])*(uFluid[ii] - uSolid[ii]);
      }
    double du = sqrt(du2);
    double weight = 1.;
    double gDrag1 = (aDarcy_*sedF*nu/((1. - sedF)*grain_*grain_) + betaForch_*du/grain_); // Darcy forchheimer term     
    double Cd = 0.44;
    double Re = (1-sedF)*du*grain_/nu;
    if (Re < 1000) 
      {
	    Cd = 24*(1. + 0.15*pow(Re, 0.687))/Re;                       //Particle cloud resistance Wen and Yu 1966
      }
    double gDrag2 = 0.75 * Cd * du * pow(1. - sedF, -1.65)/grain_;
	  
    
    if(sedF < packFraction_ + packMargin_) 
      { 
	if (sedF > packFraction_ - packMargin_)
	  {	  
	    weight =  0.5 + 0.5* (sedF - packFraction_) /packMargin_;
	  }
	else
	  {
	    weight =  0.;
	  }
      }
    

    
    return weight*gDrag1 + (1.-weight)*gDrag2;
    }

    inline double gs0(double sedF)
    {
      double g0(0.);
      if(sedF< 0.635)
	{
	  if(sedF< 0.49)
	    {
	      g0 = 0.5*(2.-sedF)/((1-sedF)*(1-sedF)*(1-sedF));
	    }
	  else 
	    {
	      g0= 0.853744035/(0.64-sedF);
	    }
	}
      else g0 = 170.74880702;
      return g0;
    }

       
		      
    inline double kappa_sed(
		      double sedF, // Sediment fraction
		      double rhoFluid,
		      double rhoSolid,
		      double uFluid[nSpace], //Fluid velocity
		      double uSolid[nSpace], //Sediment velocity
		      double gradC[nSpace], //Sediment velocity
		      double nu, //Kinematic viscosity
		      double theta_n,
		      double kappa_n,
		      double epsilon_n,
		      double nuT_n)
			   
    {		   
      double small = 1e-30;
      double beta = betaCoeff(sedF,uFluid,uSolid,nu)+small;
      double gs = gs0(sedF)+small;
      double l_c = sqrt(M_PI)*grain_/(24.*(sedF+small)*gs);
      double t_p = rhoSolid/beta;
      double t_c = l_c/(sqrt(theta_n) + small);
      double t_l = 0.165*kappa_n/(epsilon_n + small);
      double t_cl = std::min(t_c,t_l);
      double alpha= t_cl/(t_cl + t_p);
      double term = beta/(rhoFluid*(1.-sedF));
      double es_1 = 2.*term*rhoSolid*(1-alpha)*sedF*kappa_n;
      double U_gradC = 0.;
      for (int ii=0; ii<nSpace;  ii++)
	{
	  U_gradC+= (uFluid[ii] - uSolid[ii])*gradC[ii];
	}


      double es_2 = term *rhoFluid*nuT_n*U_gradC ;
      
      return -es_1 + es_2;

    }

    inline double eps_sed(
		      double sedF, // Sediment fraction
                     double rhoFluid,
		     double rhoSolid,
		      double uFluid[nSpace], //Fluid velocity
		      double uSolid[nSpace], //Sediment velocity
		      double gradC[nSpace], //Sediment velocity
		      double nu, //Kinematic viscosity
		      double theta_n,
		      double kappa_n,
		     double epsilon_n,
		     double nuT_n)
			   
    {		   
      double small = 1e-30;
      double beta = betaCoeff(sedF,uFluid,uSolid,nu)+small;
      double gs = gs0(sedF)+small;
      double l_c = sqrt(M_PI)*grain_/(24.*(sedF+small)*gs);
      double t_p = rhoSolid/beta;
      double t_c = l_c/(sqrt(theta_n) + small);
      double t_l = 0.165*kappa_n/(epsilon_n + small);
      double t_cl = std::min(t_c,t_l);
      double alpha= t_cl/(t_cl + t_p);
      double term = beta/(rhoFluid*(1.-sedF));
      double es_1 = 2.*term*rhoSolid*(1-alpha)*sedF*kappa_n;
      double U_gradC = 0.;
      for (int ii=0; ii<nSpace;  ii++)
	{
	  U_gradC+= (uFluid[ii] - uSolid[ii])*gradC[ii];
	}


      double es_2 = term *rhoFluid*nuT_n*U_gradC ;
      
      return -C3e_ * es_1 * epsilon_n/kappa_n +C4e_ * es_2 * epsilon_n/kappa_n;

    }

    inline double psc(		      double sedF, 
				      double rhoSolid,
				      double theta )
      
    {
      double gs_0 = gs0(sedF);
      double eRp1 = 1.+eR_;
      double psc = rhoSolid * sedF * (1.+2.*eRp1*sedF*gs_0)*theta;
      return psc;
    }

    inline double psc_term(		      double sedF, 
					      double rhoSolid,
					      double theta_np1,
					      double du_dx,
					      double dv_dy,
					      double dw_dz)
      
    {
      
      return -2.*psc(sedF,rhoSolid,theta_np1)*(du_dx + dv_dy + dw_dz)/(3.*rhoSolid * sedF);
    }

    inline double dpsc_term_dtheta(		      double sedF, 
					      double rhoSolid,
					      double du_dx,
					      double dv_dy,
					      double dw_dz)
      
    {
      
      double gs_0 = gs0(sedF);
      double eRp1 = 1.+eR_;
      double dpsc = rhoSolid * sedF * (1.+2.*eRp1*sedF*gs_0);
      return -2.*dpsc*(du_dx + dv_dy + dw_dz)/(3.*rhoSolid * sedF);
    }

    /*
    inline double mu_sc(		      double sedF, 
					      double rhoSolid, 
					      double theta_n )

    {
      double gs_0 = gs0(sedF);
      double sq_pi = (sqrt(M_PI));
      double sedF2 = sedF * sedF;
      double eRp1 = 1.+eR_;
      double sq_theta = sqrt(theta_n);

      double mu_sc = rhoSolid*grain_*sq_theta*( 0.8*sedF2*gs_0*eRp1/(sq_pi) + (1./15.)*sq_pi*sedF2*gs_0*eRp1 + (1./6.)*sq_pi*sedF + (5./48)*sq_pi/(gs_0*eRp1) );

      return mu_sc;
    }

    inline double lamda(		      double sedF, 
					      double rhoSolid, 
					      double theta_n )

    {
      double gs_0 = gs0(sedF);
      double sq_pi = (sqrt(M_PI));
      double sedF2 = sedF * sedF;
      double eRp1 = 1.+eR_;
      double sq_theta = sqrt(theta_n);
      double lambda = (4./3.)*sedF2*rhoSolid*grain_*gs_0*eRp1*sq_theta/sq_pi;

      return lambda;
    }



    inline double psc_term (      double sedF, 
				  double theta_n,
				  double du_dx,
				  double dv_dy,
				  double dw_dz)
    {


      return -(2./3.)*psc_head(sedF,theta_n)*(du_dx + dv_dy + dw_dz);
      
    }

    inline double d_psc_term_dtheta (      double sedF, 
				  double du_dx,
				  double dv_dy
				  double dw_dz)
    {

      double gs_0 = gs0(sedF);
      double sq_pi = (sqrt(M_PI));
      double sedF2 = sedF * sedF;
      double eRp1 = 1.+eR_;
      double dpsc_dtheta = 1.+2.*eRp1*sedF*gs_0;

      return -(2./3.)* dpsc_dtheta * (du_dx + dv_dy + dw_dz);
      
    }


    inline double tausc_term(
		      double sedF, // IN
		      double rhoSolid, //IN
		      double theta_n, //IN
		      double du_dx,
		      double du_dy,
		      double du_dz,
		      double dv_dx,
		      double dv_dy,
		      double dv_dz,
		      double dw_dx,
		      double dw_dy,
		      double dw_dz)
		     
			   
    {		   

      double gs_0 = gs0(sedF);
      double sq_pi = (sqrt(M_PI));
      double sedF2 = sedF * sedF;
      double eRp1 = 1.+eR_;
      double sq_theta = sqrt(theta_n);


	
      double la = lambda(sedF, rhoSolid, theta_n);
      double mu = mu_sc(sedF, rhoSolid, theta_n);

      
      double divU = du_dx + dv_dy + dw_dz;
      double t11 = 2*mu*du_dx + (la - (2./3.)*mu)*divU ;
      double t22 = 2*mu*dv_dy + (la - (2./3.)*mu)*divU ;
      double t33 = 2*mu*dw_dz + (la - (2./3.)*mu)*divU ;
      double t13 = mu*(du_dz + dw_dx);
      double t23 = mu*(dv_dz + dw_dy);
      double t12 = mu*(du_dy + dv_dx);
      // No need to define t31, t32 and t21 as tensor is symmetric

      //Taking the half of the product, as it is multiplied with 0.5 in the equations
      double term = t11 * du_dx + t22 * dv_dy + t33 * dw_dz;
      term += t12*du_dy + t12*dv_dx;
      term += t13*du_dz + t13*dw_dx;
      term += t23*dv_dz + t23*dw_dy;						   
      
      return term* (2./3.) / (rhoSolid*sedF);

    }

    inline double k_diff(		      double sedF, 
					      double rhoSolid, 
					      double theta_n )

    {
      double gs_0 = gs0(sedF);
      double sq_pi = (sqrt(M_PI));
      double sedF2 = sedF * sedF;
      double eRp1 = 1.+eR_;
      double sq_theta = sqrt(theta_n);
      double k_diff = rhoSolid*grain_*sq_theta*( 2.*sedF2*gs_0*eRp1/(sq_pi) + (0.5625)*sq_pi*sedF2*gs_0*eRp1 + (0.9375)*sq_pi*sedF + (0.390625)*sq_pi/(gs_0*eRp1) );


      return k_diff;
    }
    

    inline double  diffusion_theta_rhs(  double sedF,
					 double rhoSolid,
				     double theta_n,
				     double dtheta_dx,
				     double dtheta_dy,
				     double dtheta_dz
			      )
    {
      
      double gs_0 = gs0(sedF);
      double sq_pi = (sqrt(M_PI));
      double sedF2 = sedF * sedF;
      double eRp1 = 1.+eR_;
      double sq_theta = sqrt(theta_n);

      double divTheta = dtheta_dx + dtheta_dy + dtheta_dz; 

      double k_diff = k_diff(rhoSolid);

      return kappa*divTheta;
      
    }


    inline double  gamma_s(  double sedF,
			     double rhoSolid,
			     double theta_n,
			     double theta_np1,
			     double du_dx,
			     double dv_dy,
			     double dw_dz)

			      )
    {
      
      double gs_0 = gs0(sedF);
      double sq_pi = (sqrt(M_PI));
      double sedF2 = sedF * sedF;
      double eRp1 = 1.+eR_;
      double sq_theta = sqrt(theta_n);
      double divU = du_dx + dv_dy + du_dz;
      double gamma_s = 3*(1-eR_*eR_)*sedF*gs0(sedF)*(4.*sq_theta/(sq_pi*grain_ - divU))*theta_np1
      return kappa*divTheta;
      
    }

 
    inline double  Jint1(  double sedF,
			   double uFluid_n,
			   double uSolid_n,
			   double rhoSolid,
			   double kappa_n,
			   double theta_n,
			   double theta_np1,
			   double du_dx,
			   double dv_dy,
			   double dw_dz,
			   double nu)

			      )
    {
      
      double gs_0 = gs0(sedF);
      double beta = betaCoeff(sedF,uFluid_n,uSolid_n,nu);
      double sedF2 = sedF * sedF;
      double eRp1 = 1.+eR_;
      double sq_theta = sqrt(theta_n);
      double divU = du_dx + dv_dy + du_dz;
      double Jint1 = 2* sedF * beta *
      return kappa*divTheta;
      
      }*/



    

    inline double*  mInt(  double sedF,
			      double uFluid_np1[nSpace], //Fluid velocity
			      double uSolid_np1[nSpace], //Sediment velocity
			      double uFluid_n[nSpace], //Fluid velocity
			      double uSolid_n[nSpace], //Sediment velocity
			      double nu, //Kinematic viscosity
			      double nuT, //Turbulent viscosity
			      double gradc[nSpace]
			      )
    {

      double beta = betaCoeff(sedF,uFluid_n,uSolid_n,nu);
      double* mint2;
      mint2 = new double[nSpace];
      for  (int ii=0; ii<nSpace;  ii++)
	{
	  mint2[ii] = -sedF*beta*(uFluid_np1[ii]-uSolid_np1[ii]) - sedF*beta*nuT*gradc[ii]/sigmaC_;
	    }
      return  mint2;
      
      }




   
    inline double  dmInt_duFluid
                            (  double sedF,
			      double uFluid_n[nSpace], //Fluid velocity
			      double uSolid_n[nSpace], //Sediment velocity
			      double nu //Kinematic viscosity

			      )
    {
      return -sedF*betaCoeff(sedF,uFluid_n,uSolid_n,nu);

    }


       inline double  dmInt_duSolid
                            (  double sedF,
			      double uFluid_n[nSpace], //Fluid velocity
			      double uSolid_n[nSpace], //Sediment velocity
			      double nu //Kinematic viscosity

			      )
    {

      return +sedF*betaCoeff(sedF,uFluid_n,uSolid_n,nu);
    }

    

    
    

  double aDarcy_;
  double betaForch_;
  double grain_; 
  double packFraction_;
  double packMargin_;
  double sigmaC_;
  double C3e_;
  double C4e_;
  double eR_;

 

};
}

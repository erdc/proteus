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
		 double sigmaC
): 
 
  aDarcy_(aDarcy), 
      betaForch_(betaForch), 
      grain_(grain),
      packFraction_(packFraction),
      packMargin_(packMargin),
      sigmaC_(sigmaC)  


         
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

    inline double alphaResp( double rhoSolid,
		      double sedF, // Sediment fraction
		      double uFluid[nSpace], //Fluid velocity
		      double uSolid[nSpace], //Sediment velocity
		      double nu, //Kinematic viscosity
		      double theta_n,
		      double kappa_n,
		      double epsilon_n)
    {
      double small = 1e-30;
      double beta = betaCoeff(sedF,uFluid,uSolid,nu)+small;
      double gs = gs0(sedF)+small;
      double l_c = sqrt(M_PI)*grain_/(24.*(sedF+small)*gs);
      double t_p = rhoSolid/beta;
      double t_c = l_c/(sqrt(theta_n) + small);
      double t_l = 0.165*kappa_n/(epsilon_n + small);
      double t_cl = std::min(t_c,t_l);
      return t_cl/(t_cl + t_p);
    }
    
		      /*
    inline double tke_es1(
			   double sedF, // Sediment fraction
			   double uFluid[nSpace], //Fluid velocity
			   double uSolid[nSpace], //Sediment velocity
			   double nu, //Kinematic viscositydouble 
			   double rhoFluid,
			   double rhoSolid,
			   double gradc[nSpace],
			   double kappa_n
			   )
			   

      double beta = betaCoeff(sedF,uFluid_n,uSolid_n,nu);
      return double 2 * beta * rhoSolid 
		      */


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

 

};
}

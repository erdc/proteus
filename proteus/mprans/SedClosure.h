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

    inline double  turbSusp(  double sedF,
			      double uFluid[nSpace], //Fluid velocity
			      double uSolid[nSpace], //Sediment velocity
			      double nu, //Kinematic viscosity
			      double nuT
			      )
    {

      return betaCoeff(sedF,uFluid,uSolid,nu)*nuT/sigmaC_;
    }
    

    
    

  double aDarcy_;
  double betaForch_;
  double grain_; 
  double packFraction_;
  double packMargin_;
  double sigmaC_;

 

};
}

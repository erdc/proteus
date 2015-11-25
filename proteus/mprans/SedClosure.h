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
		 double parameterIn,  // Dummy parameter for first function
		 double aDarcy, // darcy parameter for drag term. Default value from Ergun (1952) is 150
		 double betaForch, // forchheimer parameter for drag term. Default value from Ergun (1952) is 1.75
		 double grain, // Grain size, default assumed as d50
		 double packFraction //Critical volume fraction for switching the drag relation 0.2 by default, see Chen and Hsu 2014

): 
 
  parameter_(parameterIn),
    aDarcy_(aDarcy), 
    betaForch_(betaForch), 
    grain_(grain),
    packFraction_(packFraction)

         
    {}
  inline double  dummy(double porosity, double pf)
  {
    double M_sf_x = parameter_*(porosity*porosity + pf);
    return M_sf_x; 
  }
  
  inline double  granularDrag(
			      double sedF, // Sediment fraction
			      double uFluid[nSpace], //Fluid velocity
			      double uSolid[nSpace], //Sediment velocity
			      double nu //Kinematic viscosity
			      )
  {
    double gDrag(0.);
    double du2 = 0.;// (uFluid - uSolid);
    for (int ii=0; ii<nSpace;  ii++)
      {
	du2+= (uFluid[ii] - uSolid[ii])*(uFluid[ii] - uSolid[ii]);
      }
    double du = sqrt(du2);
    if(sedF > packFraction_)
      {
	gDrag = aDarcy_*sedF*nu/((1. - sedF)*grain_*grain_) + betaForch_*du/grain_;  // Darcy forchheimer term
      }
    else
      {
	double Cd = 0.44;
	double Re = (1-sedF)*du*grain_/nu;
	if (Re < 1000) 
	  {
	    Cd = 24*(1. + 0.15*pow(Re, 0.687))/Re;                       //Particle cloud resistance Wen and Yu 1966
	  }
	gDrag = 0.75 * Cd * du * pow(1. - sedF, -1.65)/grain_;
	  
      }
      
    return gDrag;
    }

   double parameter_;
  double aDarcy_;
  double betaForch_;
  double grain_; 
  double packFraction_;

 

};
}

#ifndef WAVETOOLS_H
#define WAVETOOLS_H

#include <cmath>
#include <iostream>
namespace proteus
{
 const int nDim(3);
 class cppWaveGen 
 {
 public:
 cppWaveGen()
   {}
  

   inline double eta_mode(double x[nDim], double t, double kDir[nDim], double omega, double phi, double amplitude)
  {

    double phase = x[0]*kDir[0]+x[1]*kDir[1]+x[2]*kDir[2] - omega*t  + phi;
    double eta = amplitude*cos(phase);
    return eta;

  }


      inline double* vel_mode(double x[nDim], double t, double kDir[nDim], 
			   double kAbs, double omega, double phi, double amplitude,double mwl, double depth, double waveDir[nDim], double vDir[nDim])
   {

     double phase = x[0]*kDir[0]+x[1]*kDir[1]+x[2]*kDir[2] - omega*t  + phi;
      double Z =  (vDir[0]*x[0] + vDir[1]*x[1]+ vDir[2]*x[2]) - mwl;
      double UH=amplitude*omega*cosh(kAbs*(Z + depth))*cos( phase )/sinh(kAbs*depth);
      double UV=amplitude*omega*sinh(kAbs*(Z + depth))*sin( phase )/sinh(kAbs*depth);
     //Setting wave direction
      double* V = waveDir;
      for(int ii=0; ii<nDim ; ii++)
	  {
	    waveDir[ii] = waveDir[ii]/kAbs;
	    V[ii] = UH*waveDir[ii] + UV*vDir[ii] ;
	  }
      //Setting wave velocities
      return V;
      }


 
//=======================================================LINEAR - USE above directly================================================================


//---------------------------------------------------------NONLINEAR FENTON-------------------------------------------------------------------------

 };


}
#endif

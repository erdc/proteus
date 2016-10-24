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
	    V[ii] = UH*waveDir[ii] + UV*vDir[ii] ;
	  }
      //Setting wave velocities
      return V;
      }


 
//=======================================================LINEAR - USE above directly================================================================


//---------------------------------------------------------NONLINEAR FENTON-------------------------------------------------------------------------

      inline double etaFenton(double x[nDim], double t, double kDir[nDim], double kAbs, double omega, 
			      double phi0, double amplitude, int Nf, double* Ycoeff)


      {

        int ii =0;
	double HH = 0.;
	double om = 0.;
	double kw[3] = {0.,0.,0.};
	double phi = 0.;

        for (int nn=0; nn<Nf; nn++)
	  {
            ii+=1;
	    om = ii*omega;
	    kw = {ii*kDir[0], ii*kDir[1], ii*kDir[2]};
	    phi = ii*phi0;
	    HH= HH + eta_mode(x,t,kw,om,phi,Ycoeff[nn]);
	  }
        return HH/kAbs;
      }

      inline double* uFenton(double x[nDim],double t,double kDir[nDim],double kAbs,double omega,double phi0,double amplitude,
			    double mwl, double depth, double gAbs, int Nf, double* Bcoeff ,double mV[nDim], double waveDir[nDim], double vDir[nDim] )


      {

	int ii =0;
	double om = 0.;
	double kdir[3] = {0.,0.,0.};
	double phi = 0.;
	double kmode = 0.;
	double amp = 0.;
	double Ufenton[3] = {0.,0.,0.};
	double* Uf;
	  
        for (int nn=0; nn<Nf; nn++)
	  {
	    ii+=1;
	    om = ii*omega; 
	    kdir = {ii*kDir[0], ii*kDir[1], ii*kDir[2]};
	    kmode = ii*kAbs;
	    phi = ii*phi0;
            amp = tanh(kmode*depth)*sqrt(gAbs/kAbs)*Bcoeff[nn]/omega;
	    Uf  =  vel_mode(x, t ,kdir, kmode, om, phi, amp, mwl, depth, waveDir, vDir);
	    for ( int nn = 0; nn<3; nn++)
	      {
		Ufenton[nn] += Uf[nn];
	      }

	  }
	for ( int nn = 0; nn<3; nn++)
	  {
	    Ufenton[nn] = Ufenton[nn]+mV[nn];
	  }
        return Ufenton;
	  
      }


      
 };


}
#endif

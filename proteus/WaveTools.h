#ifndef WAVETOOLS_H
#define WAVETOOLS_H

#include <cmath>
#include <iostream>
namespace proteus
{
 const int nDim(3);
  

 inline double __cpp_eta_mode(double x[nDim], double t, double kDir[nDim], double omega, double phi, double amplitude)
  {

    double phase = x[0]*kDir[0]+x[1]*kDir[1]+x[2]*kDir[2] - omega*t  + phi;
    double eta = amplitude*cos(phase);
    return eta;

  }


 inline double* __cpp_vel_mode(double x[nDim], double t, double kDir[nDim],double kAbs, double omega, double phi, double amplitude,double mwl, double depth, double waveDir[nDim], double vDir[nDim])
   {

     double phase = x[0]*kDir[0]+x[1]*kDir[1]+x[2]*kDir[2] - omega*t  + phi;
      double Z =  (vDir[0]*x[0] + vDir[1]*x[1]+ vDir[2]*x[2]) - mwl;
      double UH=amplitude*omega*cosh(kAbs*(Z + depth))*cos( phase )/sinh(kAbs*depth);
      double UV=amplitude*omega*sinh(kAbs*(Z + depth))*sin( phase )/sinh(kAbs*depth);
     //Setting wave direction
      double* VV;
      VV = new double[nDim];
      for(int ii=0; ii<nDim ; ii++)
	  {
	    VV[ii] = UH*waveDir[ii] + UV*vDir[ii];
	  }
      return VV;
      delete [] VV;
      }


 
//=======================================================LINEAR - USE above directly================================================================


//---------------------------------------------------------NONLINEAR FENTON-------------------------------------------------------------------------

 inline double __cpp_etaFenton(double x[nDim], double t, double kDir[nDim], double kAbs, double omega, 
			      double phi0, double amplitude, int Nf, double* Ycoeff)


      {

        int ii =0;
	double HH = 0.;
	double om = 0.;
	double kw[3] = {0.,0.,0.};
	double phi = 0.;

        for (int nn=0; nn<Nf; nn++)
	  {
            ii= nn + 1;
	    om = ii*omega;
	    kw[0] = ii*kDir[0];
	    kw[1] = ii*kDir[1];
	    kw[2] = ii*kDir[2];
	    phi = ii*phi0;
	    HH= HH + __cpp_eta_mode(x,t,kw,om,phi,Ycoeff[nn]);
	  }
        return HH/kAbs;
      }

 inline double* __cpp_uFenton(double x[nDim],double t,double kDir[nDim],double kAbs,double omega,double phi0,double amplitude,
			    double mwl, double depth, double gAbs, int Nf, double* Bcoeff ,double* mV, double waveDir[nDim], double vDir[nDim] )


      {

	int ii =0;
	double om = 0.;
	double kw[nDim] = {0.,0.,0.};
	double phi = 0.;
	double kmode = 0.;
	double amp = 0.;
	double* Ufenton;
	double* Uf;
	Uf = new double[nDim];
	Uf[0] = 0.;
	Uf[1] = 0.;
	Uf[2] = 0.;



	
	
        for ( int nn=0; nn<Nf; nn++)
	  {
	    ii=nn+1;
	    om = ii*omega; 
	    kmode = ii*kAbs;
	    kw[0] = ii*kDir[0];
	    kw[1] = ii*kDir[1];
	    kw[2] = ii*kDir[2];
	    phi = ii*phi0;
            amp = tanh(kmode*depth)*sqrt(gAbs/kAbs)*Bcoeff[nn]/omega;
	    Ufenton = __cpp_vel_mode(x, t ,kw, kmode, om, phi, amp, mwl, depth, waveDir, vDir); 
	    Uf[0] = Uf[0]+ *(Ufenton);//[0];
	    Uf[1] = Uf[1]+ *(Ufenton+1);//[1];
	    Uf[2] = Uf[2]+ *(Ufenton+2);//[2];

	  }
	
	for ( int kk = 0; kk<3; kk++)
	  {
	    Uf[kk] = Uf[kk]+mV[kk];
	    }

        return Uf;

	delete [] Ufenton;      
	delete [] Uf;


      
 };


}
#endif

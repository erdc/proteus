#ifndef WAVETOOLS_H
#define WAVETOOLS_H

#include <cmath>
#include <iostream>
namespace proteus
{
 const int nDim(3);
 const double PI_ = M_PI;
 const double Pi2_ = (2.*PI_);
 const double 	Pi2inv_ = (1./Pi2_);
 const double 	Pihalf_ = (0.5*PI_);
 const double 	Pihalfinv_ = (1./Pihalf_);
 const double 	Pi03_ = (0.3*PI_);
 const double 	Pi07_ =  (0.7*PI_);
 const double 	Pi17_ =  (1.7*PI_);


 inline double fastcosh(double k, double Z ,double d, bool cosh)
 {

       double Kd = k * (Z+d);
       double Kd2 = Kd *  Kd *0.5;
       double Kd3 = Kd2 * Kd * 3.3333333333E-01;
       double Kd4 = Kd3 * Kd* 2.5000000000E-01;
       double Kd5 = Kd4 * Kd *2.0000000000E-01;
       double Kd6 = Kd5 * Kd*1.6666666667E-01;
       double Kd7 = Kd6 * Kd*1.4285714286E-01;
       double  Kd8 = Kd7 * Kd*1.2500000000E-01; 
       double Kd9 = Kd8 * Kd*1.1111111111E-01;
       double  Kd10 =Kd9 * Kd*0.1;
       double hype = 0.;
       if(cosh)
	 {
	   hype = 1. + Kd2  + Kd4  + Kd6   + Kd8   + Kd10;
	 }
       else
	 {
	   hype =      Kd   + Kd3  + Kd5   + Kd7   + Kd9;
	 }
       
       return hype;
     
 }



  

 inline double __cpp_eta_mode(double x[nDim], double t, double kDir[nDim], double omega, double phi, double amplitude)
  {

    double phase = x[0]*kDir[0]+x[1]*kDir[1]+x[2]*kDir[2] - omega*t  + phi;
    double eta = amplitude*cos(phase);
    return eta;

  }


 inline double* __cpp_vel_mode(double x[nDim], double t, double kDir[nDim],double kAbs, double omega, double phi, double amplitude,double mwl, double depth, double waveDir[nDim], double vDir[nDim], double sinhkd)
   {

     double phase = x[0]*kDir[0]+x[1]*kDir[1]+x[2]*kDir[2] - omega*t  + phi;
      double Z =  (vDir[0]*x[0] + vDir[1]*x[1]+ vDir[2]*x[2]) - mwl;
      

      double fcosh = fastcosh(kAbs, Z, depth, true); 
      double fsinh = fastcosh(kAbs, Z, depth, false); 

      double UH=amplitude*omega*fcosh*cos( phase )/sinhkd;//sinh(kAbs*depth);
      double UV=amplitude*omega*fsinh*sin( phase )/sinhkd;//sinh(kAbs*depth);
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
			      double mwl, double depth, double gAbs, int Nf, double* Bcoeff ,double* mV, double waveDir[nDim], double vDir[nDim], double* sinhF , double* tanhF)


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

	double sqrtAbs(sqrt(gAbs/kAbs));

	
	
        for ( int nn=0; nn<Nf; nn++)
	  {
	    ii=nn+1;
	    om = ii*omega; 
	    kmode = ii*kAbs;
	    kw[0] = ii*kDir[0];
	    kw[1] = ii*kDir[1];
	    kw[2] = ii*kDir[2];
	    phi = ii*phi0;
            amp = tanhF[nn]*sqrtAbs*Bcoeff[nn]/omega;
	    Ufenton = __cpp_vel_mode(x, t ,kw, kmode, om, phi, amp, mwl, depth, waveDir, vDir, sinhF[nn]); 
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

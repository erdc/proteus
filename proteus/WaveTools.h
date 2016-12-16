#ifndef WAVETOOLS_H
#define WAVETOOLS_H

#include <cmath>
#include <iostream>
#include <valarray>
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


 inline double fastcosh(double k, double Z , bool cosh)
 {

       double Kd = k * Z;
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

 inline double fastcos(double phi)
 {
// Setting the phase between 0 and 2pi
   double phiphi = phi - Pi2_ * int(phi * Pi2inv_);
   if(phiphi<0.){phiphi += Pi2_;}

   int signcos = 1;
   //Setting the 1.7pi -> 2p branch to -Pi/4 to 0
   if(phiphi > Pi17_){ phiphi = phiphi - Pi2_; }
   //The angle phi must be between -0.3Pi to 0.7Pi for the polynomials
   if(phiphi > Pi07_ ){  phiphi = phiphi - PI_; signcos = -1;}




	  //Calculating approximation. Accuracy <0.4%
   
   double phi2 = phiphi * phiphi *0.5;
   double phi4 = phi2*phi2*0.16666666666666666667;
   double fastc =  1. - phi2 + phi4;
   double fastcos = signcos*fastc;

	  //Choosing the right Quadrant
   if(phiphi >= Pi03_){
     
     double phi1 = phiphi - Pihalf_;	  
     double phi3 = phi1 * phi1 *phi1 * 0.166666666666667;
     double fastc1 = - phi1 + phi3;

     fastcos = signcos*fastc1;}

   return fastcos;
     
 }


  

 inline double __cpp_eta_mode(double x[nDim], double t, double kDir[nDim], double omega, double phi, double amplitude)
  {

    double phase = x[0]*kDir[0]+x[1]*kDir[1]+x[2]*kDir[2] - omega*t  + phi;
    double eta = amplitude*fastcos(phase);
    return eta;

  }


 inline double* __cpp_vel_mode(double x[nDim], double t, double kDir[nDim],double kAbs, double omega, double phi, double amplitude,double mwl, double depth, double waveDir[nDim], double vDir[nDim], double tanhkd)
   {

     double phase = x[0]*kDir[0]+x[1]*kDir[1]+x[2]*kDir[2] - omega*t  + phi;
      double Z =  (vDir[0]*x[0] + vDir[1]*x[1]+ vDir[2]*x[2]) - mwl;
       double* VV;
      VV = new double[nDim];
      double Uhype =0.;
      double Vhype =0.;
      double fcosh =0.;
      double fsinh =0.;

      
      if(kAbs*Z > -PI_)
      {
	fcosh = fastcosh(kAbs, Z,  true); 
	fsinh = fastcosh(kAbs, Z,  false); 
	Uhype = fcosh / tanhkd + fsinh; 
	Vhype = fsinh / tanhkd + fcosh; 
      }
      double fcos = fastcos(phase);
      double fsin = fastcos(Pihalf_ - phase);
	
      double UH=amplitude*omega*Uhype*fcos;
      double UV=amplitude*omega*Vhype*fsin;
     //Setting wave direction
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
			      double mwl, double depth, double gAbs, int Nf, double* Bcoeff ,double* mV, double waveDir[nDim], double vDir[nDim], double* tanhF)


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
	    Ufenton = __cpp_vel_mode(x, t ,kw, kmode, om, phi, amp, mwl, depth, waveDir, vDir, tanhF[nn]); 
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


      
 }


 


//---------------------------------------------------------PLANE RANDOM-------------------------------------------------------------------------

 inline double __cpp_etaRandom(double x[nDim], double t, double* kDir, double* omega, double* phi, double* amplitude, int N)


      {

	double HH = 0.;
	double kw[nDim] = {0.,0.,0.};
	int ii =0;

        for (int nn=0; nn<N; nn++)
	  {
	    ii = 3*nn;
	    kw[0] = kDir[ii];
	    kw[1] = kDir[ii+1];
	    kw[2] = kDir[ii+2];
	    HH= HH + __cpp_eta_mode(x,t,kw,omega[nn],phi[nn],amplitude[nn]);
	  }
        return HH;
      }

 inline double* __cpp_uRandom(double x[nDim],double t,double* kDir,double* kAbs, double* omega, double* phi, double* amplitude, double mwl, double depth, int N, double waveDir[nDim], double vDir[nDim], double* tanhF )


      {

	double kw[nDim] = {0.,0.,0.};
	double* Ufenton;
	double* Uf;
	Uf = new double[nDim];
	Uf[0] = 0.;
	Uf[1] = 0.;
	Uf[2] = 0.;


	int ii =0;

        for (int nn=0; nn<N; nn++)
	  {
	    ii = 3*nn;
	    kw[0] = kDir[ii];
	    kw[1] = kDir[ii+1];
	    kw[2] = kDir[ii+2];
	    Ufenton = __cpp_vel_mode(x, t ,kw, kAbs[nn], omega[nn], phi[nn], amplitude[nn], mwl, depth, waveDir, vDir, tanhF[nn]); 
	    Uf[0] = Uf[0]+ *(Ufenton);//[0];
	    Uf[1] = Uf[1]+ *(Ufenton+1);//[1];
	    Uf[2] = Uf[2]+ *(Ufenton+2);//[2];

	  }
	
        return Uf;

	delete [] Ufenton;      
	delete [] Uf;


      
 }
//---------------------------------------------------------Directional RANDOM / Velocity-------------------------------------------------------------------------

 inline double* __cpp_uDir(double x[nDim],double t,double* kDir,double* kAbs, double* omega, double* phi, double* amplitude, double mwl, double depth, int N, double* waveDir, double vDir[nDim], double* tanhF )


      {

	double kw[nDim] = {0.,0.,0.};
	double wd[nDim] = {0.,0.,0.};
	double* Ufenton;
	double* Uf;
	Uf = new double[nDim];
	Uf[0] = 0.;
	Uf[1] = 0.;
	Uf[2] = 0.;


	int ii =0;

        for (int nn=0; nn<N; nn++)
	  {
	    ii = 3*nn;
	    kw[0] = kDir[ii];
	    kw[1] = kDir[ii+1];
	    kw[2] = kDir[ii+2];
	    wd[0] = waveDir[ii];
	    wd[1] = waveDir[ii+1];
	    wd[2] = waveDir[ii+2];
	    Ufenton = __cpp_vel_mode(x, t ,kw, kAbs[nn], omega[nn], phi[nn], amplitude[nn], mwl, depth, wd, vDir, tanhF[nn]); 
	    Uf[0] = Uf[0]+ *(Ufenton);
	    Uf[1] = Uf[1]+ *(Ufenton+1);
	    Uf[2] = Uf[2]+ *(Ufenton+2);

	  }
	
        return Uf;

	delete [] Ufenton;      
	delete [] Uf;


      }    
 
//---------------------------------------------------------Time series/ Velocity-------------------------------------------------------------------------

 inline int __cpp_findWindow(double t, double handover, double t0, double Twindow, int Nwindows, double* windows_handover)
 {
  double term = 1. - handover;
  int Nw = 0;
    if (t-t0 >= term*Twindow)
    {
      Nw = std::min(int((t-t0 - term*Twindow)/(Twindow - 2. * handover * Twindow)) + 1, Nwindows-1);
      if (t-t0 < windows_handover[Nw-1] - t0)
	{
	  Nw =Nw - 1;
	}
    }
    return Nw;
 }


inline double __cpp_etaDirect(double x[nDim], double x0[nDim], double t, double* kDir, double* omega, double* phi, double* amplitude, int N)


{ 
  x[0] = x[0] - x0[0];
  x[1] = x[1] - x0[1];
  x[2] = x[2] - x0[2];

  return __cpp_etaRandom(x,  t,  kDir,  omega,  phi,  amplitude,  N); 
}





   
 inline double* __cpp_uDirect(double* x, double* x0, double t, double* kDir, double* kAbs, double* omega, double* phi, double* amplitude, double mwl, double depth, int N, double* waveDir, double* vDir, double* tanhKd )



 {
	x[0] = x[0] - x0[0];
	x[1] = x[1] - x0[1];
	x[2] = x[2] - x0[2];

	return __cpp_uRandom(x, t,  kDir,  kAbs,  omega,  phi,  amplitude,  mwl,  depth,  N,  waveDir,  vDir,  tanhKd );

 }


 inline double __cpp_etaWindow(double x[nDim], double x0[nDim], double t, double* t0, double* kDir, double* omega, double* phi, double* amplitude, int N, int Nw)


{ 
  int Is = Nw*N;
  x[0] = x[0] - x0[0];
  x[1] = x[1] - x0[1];
  x[2] = x[2] - x0[2];
  t = t-t0[Nw];
  double HH = 0.;
  double kw[nDim] = {0.,0.,0.};
  int ii =0;

  for (int nn=Is; nn<Is+N; nn++)
	  {
	    ii = 3*nn;
	    kw[0] = kDir[ii];
	    kw[1] = kDir[ii+1];
	    kw[2] = kDir[ii+2];
	    HH= HH + __cpp_eta_mode(x,t,kw,omega[nn],phi[nn],amplitude[nn]);
	  }
  return HH;


}



 inline double* __cpp_uWindow(double* x, double* x0, double t, double* t0, double* kDir, double* kAbs, double* omega, double* phi, double* amplitude, double mwl, double depth, int N,int Nw, double* waveDir, double* vDir, double* tanhF )


 {
  int Is = Nw*N;
  x[0] = x[0] - x0[0];
  x[1] = x[1] - x0[1];
  x[2] = x[2] - x0[2];
  t = t-t0[Nw];

  double kw[nDim] = {0.,0.,0.};
  double* Ufenton;
  double* Uf;
  Uf = new double[nDim];
  Uf[0] = 0.;
  Uf[1] = 0.;
  Uf[2] = 0.;


  int ii =0;

  for (int nn=Is; nn<Is+N; nn++)
    {
      ii = 3*nn;
      kw[0] = kDir[ii];
      kw[1] = kDir[ii+1];
      kw[2] = kDir[ii+2];
      Ufenton = __cpp_vel_mode(x, t ,kw, kAbs[nn], omega[nn], phi[nn], amplitude[nn], mwl, depth, waveDir, vDir, tanhF[nn]); 
      Uf[0] = Uf[0]+ *(Ufenton);//[0];
      Uf[1] = Uf[1]+ *(Ufenton+1);//[1];
      Uf[2] = Uf[2]+ *(Ufenton+2);//[2];

    }
	
        return Uf;

	delete [] Ufenton;      
	delete [] Uf;


      
 }
 //=========================================2nd order correction==============================================

 inline double __cpp_eta2nd(double x[nDim], double t, double* kDir, double* ki, double* omega, double* phi, double* amplitude, int N, double* sinhKd, double* tanhKd)


 {

	

	double HH = 0.;
	double kw[nDim] = {0.,0.,0.};
	int ii =0;
	double ai_2nd=0.;


        for (int nn=0; nn<N; nn++)
	  {
	    ii = 3*nn;
	    kw[0] =2.* kDir[ii];
	    kw[1] =2.* kDir[ii+1];
	    kw[2] =2.* kDir[ii+2];
            ai_2nd = (amplitude[nn]*amplitude[nn] * ki[nn]*(2+3./(sinhKd[nn]*sinhKd[nn]) )/(4.*tanhKd[nn] ));	    	    
	    HH= HH + __cpp_eta_mode(x,t,kw,2.*omega[nn],2.*phi[nn],ai_2nd);
	  }
        return HH;      
      }

 inline double __cpp_eta_short(double x[nDim], double t, double* kDir, double* ki, double* omega, double* phi, double* amplitude, int N, double* sinhKd, double* tanhKd, double gAbs)


 {

	

	double HH = 0.;
	double kw[nDim] = {0.,0.,0.};
	double kw2[nDim] = {0.,0.,0.};
	int ii =0;
	int jj =0;
	double Dp=0.;
	double Bp=0.;
	double tanhSum = 0.;
	double ai =0.;
        for (int i=0; i<N-1; i++)
	  {
	    ii = 3*i;
	    kw[0] = kDir[ii];
	    kw[1] = kDir[ii+1];
	    kw[2] = kDir[ii+2];

	    for (int j=i+1; j<N; j++)
	      {
		jj = 3*j;
		kw2[0] = kDir[jj]+kw[0];
		kw2[1] = kDir[jj+1]+kw[1];
		kw2[2] = kDir[jj+2]+kw[2];
		
		tanhSum = (tanhKd[i]+tanhKd[j])/(1.+tanhKd[i]*tanhKd[j]);

                Dp = pow(omega[i]+omega[j],2) - gAbs*(ki[i]+ki[j])*tanhSum;
                Bp = (pow(omega[i],2)+pow(omega[j],2))/(2*gAbs);
		Bp = Bp -((omega[i]*omega[j])/(2*gAbs))*(1-1./(tanhKd[i]*tanhKd[j]) )*((pow(omega[i]+omega[j],2) + gAbs*(ki[i]+ki[j])*tanhSum)/Dp);
		Bp = Bp + ((omega[i]+omega[j])/(2*gAbs*Dp))*((pow(omega[i],3)/pow(sinhKd[i],2)) + (pow(omega[j],3)/pow(sinhKd[j],2)));
	    
		ai = amplitude[i]*amplitude[j]*Bp;
		HH= HH + __cpp_eta_mode(x,t,kw2,omega[i]+omega[j],phi[i]+phi[j],ai);
	      }
	  }
        return HH;      
      }

 inline double __cpp_eta_long(double x[nDim], double t, double* kDir, double* ki, double* omega, double* phi, double* amplitude, int N, double* sinhKd, double* tanhKd, double gAbs)


 {

	

	double HH = 0.;
	double kw[nDim] = {0.,0.,0.};
	double kw2[nDim] = {0.,0.,0.};
	int ii =0;
	int jj =0;
	double Dp=0.;
	double Bp=0.;
	double tanhSum = 0.;
	double ai =0.;
        for (int i=0; i<N-1; i++)
	  {
	    ii = 3*i;
	    kw[0] = kDir[ii];
	    kw[1] = kDir[ii+1];
	    kw[2] = kDir[ii+2];

	    for (int j=i+1; j<N; j++)
	      {
		jj = 3*j;
		kw2[0] = -kDir[jj]+kw[0];
		kw2[1] = -kDir[jj+1]+kw[1];
		kw2[2] = -kDir[jj+2]+kw[2];
		
		tanhSum = (tanhKd[i]-tanhKd[j])/(1.-tanhKd[i]*tanhKd[j]);

                Dp = pow(omega[i]-omega[j],2) - gAbs*(ki[i]-ki[j])*tanhSum;
                Bp = (pow(omega[i],2)+pow(omega[j],2))/(2*gAbs);
		Bp = Bp + ((omega[i]*omega[j])/(2*gAbs))*(1+1./(tanhKd[i]*tanhKd[j]) )*((pow(omega[i]-omega[j],2) + gAbs*(ki[i]-ki[j])*tanhSum)/Dp);
		Bp = Bp + ((omega[i]-omega[j])/(2*gAbs*Dp))*((pow(omega[i],3)/pow(sinhKd[i],2)) - (pow(omega[j],3)/pow(sinhKd[j],2)));
	    
		ai = amplitude[i]*amplitude[j]*Bp;
		HH= HH + __cpp_eta_mode(x,t,kw2,omega[i]-omega[j],phi[i]-phi[j],ai);
	      }
	  }
        return HH;      
      }




};      
       

#endif

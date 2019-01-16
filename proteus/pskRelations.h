#ifndef PSKRELATIONS_H
#define PSKRELATIONS_H
#include <cmath>
#include <iostream>
#include <cassert>
#include "densityRelations.h"
#include "SubsurfaceTransportCoefficients.h"
/** \file pskrelations.h
    \defgroup pskrelations pskrelations
    \brief A library of pressure-saturation-permeability relations.
    @{
*/

/** \todo Work out what needs to happen at Se=0,1 for psk relations */

/* jcc two phase flow, modified by cek and mwf*/ 
using namespace std;

class  PskRelation
{
public:
  double Se,
    dSe_dSw,
    Sw_min,
    Sw_max,
    krw,
    dkrw,
    krn,
    dkrn,
    psic,
    dpsic,
    //in case need psic--> Se form
    dSe_dpsic;
  PskRelation(const double* rwork, const int* iwork = 0):
  Se(1.0),dSe_dSw(1.0),
    Sw_min(rwork[0]),
    Sw_max(rwork[1]),
    krw(1.0),dkrw(0.0),
    krn(0.0),dkrn(0.0),psic(0.0),
    dpsic(0.0),dSe_dpsic(0.0)
  {}
  PskRelation():
    Se(1.0),dSe_dSw(1.0),
    Sw_min(0.0),
    Sw_max(1),
    krw(1.0),dkrw(0.0),
    krn(0.0),dkrn(0.0),psic(0.0),
    dpsic(0.0),dSe_dpsic(0.0)
  {}
  virtual inline void setParams(const double* rwork, const int* iwork = 0)
  {
    Sw_min = rwork[0];
    Sw_max = rwork[1];
  }
  /*for fudge factors aka tolerances in various models*/
  virtual inline void setTolerances(const double* rwork_tol)
  {
  }
  virtual ~PskRelation(){}

  /*linear coefficients*/
  inline void calc(const double& Sw)
  {  
    calc_Se(Sw);
    
    krw  = Se;
    dkrw = dSe_dSw;
    
    krn  = (1.0-Se);
    dkrn = -dSe_dSw;
    
    psic  = Se;
    dpsic = dSe_dSw;   
  }
  
  inline void calc_Se(const double& Sw)
  {
    if (Sw < Sw_min)
      {
        Se = 0.0;
        dSe_dSw = 0.0;
      }
    else if (Sw > Sw_max)
      {
        Se = 1.0;
        dSe_dSw =0.0;
      }
    else
      {
        Se = (Sw - Sw_min)/(Sw_max-Sw_min);
        dSe_dSw = 1.0/(Sw_max - Sw_min);
      }
  }
  //for Sw as a function of capillary pressure head
  //note dkrw,dkrn are set to be wrt psic
  virtual inline void calc_from_psic(const double& psicIn)
  {
    bool implemented = false;
    assert(implemented);
  }

};

/* quadratic kr */
class SimplePSK : public PskRelation
{
public:
  SimplePSK(const double* rwork, const int* iwork = 0):
    PskRelation(rwork)
  {}

  inline void calc(const double& Sw)
  {  
    calc_Se(Sw);
    
    krw  = Se*Se;
    dkrw = 2.0*Se*dSe_dSw;
    
    krn  = (1.0-Se)*(1.0-Se); 
    dkrn = 2.0*(Se-1.0)*dSe_dSw; 
    
    psic  =  (1.0-Se)/Sw_min; /*mwf change to test het -Se;*/
    dpsic = -dSe_dSw;   
  }
  //
  virtual inline void calc_from_psic(const double& psicIn)
  {
    psic = max(0.0,psicIn);
    dpsic= 1.0;
    Se = 1.0-psic*Sw_min;
    Se = max(0.0,min(1.0,Se));
    
    krw  = Se*Se;
    dkrw = 2.0*Se*dSe_dSw;
    
    krn  = (1.0-Se)*(1.0-Se); 
    dkrn = 2.0*(Se-1.0)*dSe_dSw; 
    
  }
};

/* Van Genuchten-Mualem */ 
class VGMorig: public PskRelation
{
public:
  double alpha,
    m,
    n,
    Se_eps;
  double Se_eps_const;
  VGMorig()
  {}
  VGMorig(const double* rwork, const int* iwork = 0):
    PskRelation(rwork),
    alpha(rwork[2]),
    m(rwork[3])
  {
    n = (1.0)/(1.0-m);
    Se_eps_const=1.0e-4;
    //mwf debug
    //std::cout<<"VGMorig ctor rwork[2]= "<<rwork[2]<<" alpha= "<<alpha<<" rwork[3]= "<<rwork[3]<<" m= "<<m<<" n= "<<n<<std::endl;
  }
  /*for fudge factors aka tolerances in various models*/
  virtual inline void setTolerances(const double* rwork_tol)
  {
    Se_eps_const = rwork_tol[0];
  }
  
  inline void setParams(const double* rwork, const int* iwork = 0)
  {
    Sw_min = rwork[0];
    Sw_max = rwork[1];
    alpha = rwork[2];
    m = rwork[3];
    n = (1.0)/(1.0-m);  
    //mwf debug
    //std::cout<<"VGMorig setParams rwork[2]= "<<rwork[2]<<" alpha= "<<alpha<<" rwork[3]= "<<rwork[3]<<" m= "<<m<<" n= "<<n<<std::endl;
  }
  inline void calc_Se_eps(const double& Se)
  {
    /* cek todo, work out when you need to stay away from Se=0 and 1 and what to assign in these cases */
    if(Se <= Se_eps_const)
      Se_eps=Se_eps_const;
    else if(Se >= (1.0-Se_eps_const))
      Se_eps=1.0-Se_eps_const;
    else
      Se_eps=Se;
  }

  inline void calc(const double& Sw)
  {
    double Seovmmo,Seovm,Seovmmh,S,Smmo,S2mmo,Sm,S2m,pcesub1,pcesub2; 
    calc_Se(Sw);
    calc_Se_eps(Se);

    Seovmmo  = pow(Se_eps,((1.0/m)-1.0));
    Seovm    = Seovmmo*Se_eps;
    Seovmmh  = pow(Se_eps,((1.0/m)-0.5)); 
    S       = 1.0 - Seovm;
    Smmo    = pow(S,m-1.0);
    S2mmo   = pow(S,2.0*m-1.0);
    Sm      = Smmo*S;
    S2m     = S2mmo*S; 
    
    pcesub1 = pow(((1.0/Seovm)-1.0),((1.0/n)-1.0));
    pcesub2 = pcesub1*((1.0/Seovm)-1.0);
	
    krw  = sqrt(Se)*(1.0-Sm)*(1.0-Sm);
    dkrw = (0.5*(1.0/sqrt(Se_eps))*(1.0-Sm)*(1.0-Sm) + 2.0*(1-Sm)*Smmo*Seovmmh)*dSe_dSw;

    krn  = sqrt(1.0-Se)*S2m;
    dkrn = (-0.5*(1.0/sqrt(1.0-Se_eps))*S2m - 2.0*sqrt(1.0-Se)*(S2mmo)*Seovmmo)*dSe_dSw; 
    
    psic  = pow((pow(Se_eps,(-1.0/m)) - 1.0),(1.0/n))/alpha;
    dpsic = ((-1.0/(alpha*n*m))*(pow((pow(Se_eps,-1.0/m)-1.0),-m))*(pow(Se_eps,((-1.0/m)-1.0))))*dSe_dSw;  
  }
};  

class VGM : public VGMorig
{
 public:
  double ns_del,eps_small,sqrt_eps_small;
  VGM()
  {}
  VGM(const double* rwork, const int* iwork = 0):
    VGMorig(rwork,iwork),
    ns_del(1.0e-8),
    eps_small(1.0e-16),
    sqrt_eps_small(1.0e-8)
  {}
  VGM(const VGM& r)
  {
    Se=r.Se;
    dSe_dSw=r.dSe_dSw;
    Sw_min=r.Sw_min;
    Sw_max=r.Sw_max;
    krw=r.krw;
    dkrw=r.dkrw;
    krn=r.krn;
    dkrn=r.dkrn;
    psic=r.psic;
    dpsic=r.dpsic;
    dSe_dpsic=r.dSe_dpsic;
    alpha=r.alpha;
    m=r.m;
    n=r.n;
    Se_eps_const=r.Se_eps_const;
    ns_del=r.ns_del;
    eps_small=r.eps_small;
    sqrt_eps_small=r.sqrt_eps_small;
  }
  /*for fudge factors aka tolerances in various models*/
  virtual inline void setTolerances(const double* rwork_tol)
  {
    eps_small = rwork_tol[0]; //mwf 072110 don't tie to eps_small? sqrt_eps_small = sqrt(eps_small);
    ns_del    = rwork_tol[1];
  }

  inline void calc_Se(const double& Sw)
  {
    Se = (Sw - Sw_min)/(Sw_max-Sw_min);
    dSe_dSw = 1.0/(Sw_max - Sw_min);
    Se = max(eps_small,min(Se,1.0-eps_small));
  }

  inline void calc(const double& Sw)
  {
    calc_Se(Sw);
    //taken from MualemVanGenuchten2p in pdetk
    double sBar,psiC,DsBar_DpC,DDsBar_DDpC,DkrW_DpC,DkrN_DpC;
    double vBar,uBar,
      alphaPsiC, alphaPsiC_n, alphaPsiC_nM1, alphaPsiC_nM2,
      onePlus_alphaPsiC_n,
      sqrt_sBar, sqrt_1minusSbar,
      sBarByOnePlus_alphaPsiC_n, sBarBy_onePlus_alphaPsiC_n_2;

    sBar = Se;
    //begin MualemVanGenuchten2p setVFraction
    onePlus_alphaPsiC_n = pow(sBar,1.0/-m);
    alphaPsiC_n = onePlus_alphaPsiC_n - 1.0;
    alphaPsiC = pow(alphaPsiC_n,1.0/n);
    psiC = alphaPsiC/alpha;

    alphaPsiC_nM1 = alphaPsiC_n/alphaPsiC;
    sBarByOnePlus_alphaPsiC_n = sBar/onePlus_alphaPsiC_n;
    sqrt_sBar = sqrt(sBar);
    sqrt_1minusSbar = sqrt(1.0 - sBar);
    
    DsBar_DpC = -alpha*(n-1.0)*alphaPsiC_nM1 
      *sBarByOnePlus_alphaPsiC_n;
    //DthetaW_DpC = thetaSR[i] * DsBar_DpC; 

    vBar = 1.0-alphaPsiC_nM1*sBar;
    uBar = alphaPsiC_nM1*sBar;

    //change names krW--> krw, krN--> krn
    krw = sqrt_sBar*vBar*vBar;
    krn = sqrt_1minusSbar*uBar*uBar;
    psic= psiC;
    if(psiC<=0.0) 
    {
      DsBar_DpC = 0.0;
      //DthetaW_DpC = 0.0;
      krw = 1.0;
      krn = 0.0;
    }

    //begin MualemVanGenuchten2p calculateDerivatives
    alphaPsiC_nM2 =   alphaPsiC_nM1/alphaPsiC;      
  
    sBarBy_onePlus_alphaPsiC_n_2 = sBarByOnePlus_alphaPsiC_n
      /onePlus_alphaPsiC_n;
    DDsBar_DDpC =  alpha*alpha*(n-1.)
      *((2*n-1.)*alphaPsiC_nM1*alphaPsiC_nM1
	*sBarBy_onePlus_alphaPsiC_n_2
      -
	(n-1.)*alphaPsiC_nM2
	*sBarByOnePlus_alphaPsiC_n);

    //DDthetaW_DDpC = thetaSR[i]*DDsBar_DDpC;

    DkrW_DpC = (0.5/sqrt_sBar)*DsBar_DpC*vBar*vBar
      -
      2.0*sqrt_sBar*vBar*
      (alpha*(n-1.0)*alphaPsiC_nM2*sBar
       + alphaPsiC_nM1 * DsBar_DpC);

    //DKW_DpC = KWs[i]*DkrW_DpC;

  
    //recalculate if necessary
    if (sqrt_1minusSbar >= sqrt_eps_small)//SQRT_MACHINE_EPSILON)
      {
	DkrN_DpC = -(0.5/sqrt_1minusSbar)*DsBar_DpC*uBar*uBar
	  +
	  2.0*sqrt_1minusSbar*uBar*
	  (alpha*(n-1.0)*alphaPsiC_nM2*sBar
	   + alphaPsiC_nM1 * DsBar_DpC);
      }
    else
      {
	DkrN_DpC =((1.0 - sBar)/eps_small)*2.0*sqrt_eps_small*uBar*
	  (alpha*(n-1.0)*alphaPsiC_nM2*sBar
	   + alphaPsiC_nM1 * DsBar_DpC)
	  - (DsBar_DpC/eps_small)*sqrt_eps_small*uBar*uBar;
      }
    
    //if we're in the nonsmooth regime
    if (psiC < ns_del && psiC > 0.0 )
      {
	DkrW_DpC = 0.0;
      }

    if (psiC <= 0.0)
      {
	DDsBar_DDpC = 0.0;
	//DDthetaW_DDpC = 0.0;
	DkrW_DpC = 0.0;
	DkrN_DpC = 0.0;
      }
    //end calculateDerivatives
    double DpC_Dse = 0.0;
    if (fabs(DsBar_DpC) > 0.0)
      DpC_Dse = 1.0/DsBar_DpC;
    double DpC_Dsw = DpC_Dse*dSe_dSw;
    dkrw = DkrW_DpC*DpC_Dsw;
    dkrn = DkrN_DpC*DpC_Dsw;
    dpsic= DpC_Dsw;

  }
  //TODO add this to other classes
  virtual inline void calc_from_psic(const double& psicIn)
  {
    //taken from MualemVanGenuchten2p in pdetk
    double sBar,psiC,DsBar_DpC,DDsBar_DDpC,DkrW_DpC,DkrN_DpC;
    double vBar,uBar,
      alphaPsiC, alphaPsiC_n, alphaPsiC_nM1, alphaPsiC_nM2,
      onePlus_alphaPsiC_n,
      sqrt_sBar, sqrt_1minusSbar,
      sBarByOnePlus_alphaPsiC_n, sBarBy_onePlus_alphaPsiC_n_2;
    

    psiC = max(0.0,psicIn);
    alphaPsiC = alpha*psiC;
    alphaPsiC_n = pow(alphaPsiC,n);
    alphaPsiC_nM1 = alphaPsiC_n/alphaPsiC;
    onePlus_alphaPsiC_n = 1.0 + alphaPsiC_n;
    sBar = pow(onePlus_alphaPsiC_n,-m);
    sBarByOnePlus_alphaPsiC_n = sBar/onePlus_alphaPsiC_n;
    sqrt_sBar = sqrt(sBar);
    sqrt_1minusSbar = sqrt(1.0 - sBar);
    //thetaW = thetaSR[i]*sBar + thetaR[i];
    DsBar_DpC = -alpha*(n-1.0)*alphaPsiC_nM1 
      *sBarByOnePlus_alphaPsiC_n;
    //DthetaW_DpC = thetaSR[i] * DsBar_DpC; 
    vBar = 1.0-alphaPsiC_nM1*sBar;
    uBar = alphaPsiC_nM1*sBar;

    //change names krW--> krw, krN--> krn
    krw = sqrt_sBar*vBar*vBar;
    krn = sqrt_1minusSbar*uBar*uBar;
    Se = sBar;
    psic=psiC;
    if(psiC<=0.0) 
      {
	sBar = 1.0;
	Se = sBar;
	//thetaW = thetaS[i];
	DsBar_DpC = 0.0;
	//DthetaW_DpC = 0.0;
	krw = 1.0;
	krn = 0.0;
      }    
    //mwf debug
    //std::cout<<"vgm_calc_from_psic alpha= "<<alpha<<" m= "<<m<<" psic= "<<psic<<" sBar= "<<sBar<<" krw= "<<krw<<" krn= "<<krn
    //     <<std::endl;
    //begin MualemVanGenuchten2p calculateDerivatives
    alphaPsiC_nM2 =   alphaPsiC_nM1/alphaPsiC;      
  
    sBarBy_onePlus_alphaPsiC_n_2 = sBarByOnePlus_alphaPsiC_n
      /onePlus_alphaPsiC_n;
    DDsBar_DDpC =  alpha*alpha*(n-1.)
      *((2*n-1.)*alphaPsiC_nM1*alphaPsiC_nM1
	*sBarBy_onePlus_alphaPsiC_n_2
      -
	(n-1.)*alphaPsiC_nM2
	*sBarByOnePlus_alphaPsiC_n);

    //DDthetaW_DDpC = thetaSR[i]*DDsBar_DDpC;

    DkrW_DpC = (0.5/sqrt_sBar)*DsBar_DpC*vBar*vBar
      -
      2.0*sqrt_sBar*vBar*
      (alpha*(n-1.0)*alphaPsiC_nM2*sBar
       + alphaPsiC_nM1 * DsBar_DpC);

    //DKW_DpC = KWs[i]*DkrW_DpC;

  
    //recalculate if necessary
    if (sqrt_1minusSbar >= sqrt_eps_small)//SQRT_MACHINE_EPSILON)
      {
	DkrN_DpC = -(0.5/sqrt_1minusSbar)*DsBar_DpC*uBar*uBar
	  +
	  2.0*sqrt_1minusSbar*uBar*
	  (alpha*(n-1.0)*alphaPsiC_nM2*sBar
	   + alphaPsiC_nM1 * DsBar_DpC);
      }
    else
      {
	DkrN_DpC =((1.0 - sBar)/eps_small)*2.0*sqrt_eps_small*uBar*
	  (alpha*(n-1.0)*alphaPsiC_nM2*sBar
	   + alphaPsiC_nM1 * DsBar_DpC)
	  - (DsBar_DpC/eps_small)*sqrt_eps_small*uBar*uBar;
      }
    
    //if we're in the nonsmooth regime
    if (psiC < ns_del && psiC > 0.0 )
      {
	DkrW_DpC = 0.0;
      }

    if (psiC <= 0.0)
      {
	DDsBar_DDpC = 0.0;
	//DDthetaW_DDpC = 0.0;
	DkrW_DpC = 0.0;
	DkrN_DpC = 0.0;
      }
    //end calculateDerivatives
    dkrw = DkrW_DpC; //note: \od{k_{rw}}{\psi_c} not \od{k_{rw}}{S_w}
    dkrn = DkrN_DpC;
    dpsic= 1.0;
    dSe_dpsic=DsBar_DpC;
  }
};
/* Van Genuchten-Burdine */
class VGB : public  VGM
{
public:
  VGB(const double* rwork, const int* iwork = 0):
    VGM(rwork,iwork)
  {}

  inline void calc(const double& Sw)
  {    
    double S,Smmo,Sm,alpha,Se1ovMmo; 
    calc_Se(Sw);
    calc_Se_eps(Se);

    Se1ovMmo = pow(Se_eps,((1.0/m)-1.0));
    S = 1.0-Se1ovMmo*Se_eps;
    Smmo = pow(S,(m-1.0)); 
    Sm = Smmo*S; 

    //cek from  matlab, needs optimizing
    
    krw  = Se_eps*Se_eps*(1.0-pow(1.0-pow(Se_eps,1.0/m),1.0/m));
    dkrw = (2.0*Se_eps*(1.0-pow(1.0-pow(Se_eps,1.0/m),1.0/m))+Se_eps*pow(1.0-pow(Se_eps,1.0/m),1.0/m)/(m*m)*pow(Se_eps,1.0/m)/(1.0-pow(Se_eps,1.0/m)))*dSe_dSw;

    krn  = pow(1.0-Se_eps,2.0)*pow(1.0-pow(Se_eps,1.0/m),1.0*m);
    dkrn = (-2.0*(1.0-Se_eps)*pow(1.0-pow(Se_eps,1.0/m),1.0*m)-pow(1.0-Se_eps,2.0)*pow(1.0-pow(Se_eps,1.0/m),1.0*m)*pow(Se_eps,1.0/m)/Se_eps/(1.0-pow(Se_eps,1.0/m)))*dSe_dSw;

//     krw  = (Se*Se)*(1.0-Sm); 
//     dkrw = ((2.0*Se)*(1.0-Sm)+(Se*Se)*(Se1ovMmo)*Smmo)*dSe_dSw;
    
//     krn  = (1.0-Se)*(1.0-Se)*(Sm);
//     dkrn = (2.0*(Se-1.0)*Sm-((1.0-Se)*(1.0-Se))*Smmo*Se1ovMmo)*dSe_dSw;
    psic  = pow((pow(Se_eps,(-1.0/m)) - 1.0),(1.0/n))/alpha;
    dpsic = ((-1.0/(alpha*n*m))*(pow((pow(Se_eps,-1.0/m)-1.0),-m))*(pow(Se_eps,((-1.0/m)-1.0))))*dSe_dSw;
  }
};

/* Brooks-Corey-Mualem */					  
class BCM : public PskRelation
{
 public:
  double pd,lambda;

  BCM(const double* rwork, const int* iwork = 0):
    PskRelation(rwork,iwork),
    pd(rwork[2]),
    lambda(rwork[3])
  {}
  
  inline void setParams(const double* rwork, const int* iwork = 0)
  {
    Sw_min = rwork[0];
    Sw_max = rwork[1];
    pd = rwork[2];
    lambda = rwork[3];
  }
  inline void calc(const double& Sw)
  {  
    double Value,Expon,krwovSe,Oovbclpo,Oovbcl,X,sqrt1mu;
    calc_Se(Sw);
    
    Oovbcl   = 1.0/lambda; 
    Oovbclpo = Oovbcl+1.0; 
    X        = pow(Se,Oovbcl);
    Value    = 1.0-X*Se;
    Expon    = ((4.0+5.0*lambda)/(2.0*lambda)); 
    krwovSe   = pow(Se,(Expon-1.0));
    sqrt1mu  = sqrt(1.0-Se); 
    
    krw  = krwovSe*Se; 
    dkrw = (Expon*krwovSe)*dSe_dSw;
    
    krn  = sqrt1mu*Value*Value;
    dkrn = (-0.5*(1.0/sqrt1mu)*Value*Value - 2.0*sqrt1mu*Value*Oovbclpo*X )*dSe_dSw;
    
    psic  = 1.0/X;
    dpsic = (-Oovbcl/(X*Se))*dSe_dSw;
  }
};

/* Brooks-Corey-Burdine */					  
class BCB : public BCM
{
public:
  BCB(const double* rwork, const int* iwork = 0):BCM(rwork,iwork)
  {}

  inline void calc(const double& Sw)
  {
    double Se2ovL,Se2,Se3,omSe,Expon,Semoovlmo,Se_cutOff; 
    calc_Se(Sw);

    Se2ovL     = pow(Se,(2.0/lambda));
    Se2        = Se*Se; 
    Se3        = Se*Se2;
    omSe       = 1.0-Se;   
    Expon     = (2.0+3.0*lambda)/lambda;		
    Se_cutOff = max(1.0e-4,Se);
    Semoovlmo  = pow(Se_cutOff,((-1.0/lambda)-1.0)); 					
    
    krw  = Se2ovL*Se3; 
    dkrw = Expon*Se2ovL*Se2*dSe_dSw;
  
    krn  = (omSe*omSe)*(1.0-Se2ovL*Se);
    dkrn = (2.*omSe*(1.0-Se2ovL*Se) - (omSe*omSe)*(Expon-2.0)*Se2ovL*dSe_dSw)*dSe_dSw;
    
    /* cek debug */
    krw  = pow(Se,(2.0+3.0*lambda)/lambda);
    dkrw = (((2.0+3.0*lambda)/lambda)*pow(Se,(2.0+3.0*lambda)/lambda - 1.0))*dSe_dSw;
    
    krn  = (1.0-Se)*(1.0-Se)*(1.0-pow(Se,(2.0+lambda)/lambda));
    dkrn  = (-2.0*(1.0-Se)*(1.0-pow(Se,(2.0+lambda)/lambda))
             -((2.0+lambda)/lambda)*(1.0-Se)*(1.0-Se)*pow(Se,(2.0+lambda)/lambda-1.0))*dSe_dSw;
    
    psic  = pd*Semoovlmo*Se;
    dpsic = pd*(-1.0/lambda)*Semoovlmo*dSe_dSw;
  }
};

class FractionalFlowVariables
{
public:
  double muw,
    mun,
    lambdaw,
    dlambdaw,
    lambdan,
    dlambdan,
    lambdat,
    dlambdat,
    fw,
    dfw,
    fn,
    dfn;
  
  FractionalFlowVariables(double muwIn,double munIn):
    muw(muwIn),
    mun(munIn)
  {}

  inline void calc(const PskRelation& psk,
                   const DensityRelation& density_w,
                   const DensityRelation& density_n)
  {		
    lambdaw  = density_w.rho*psk.krw/muw; 
    dlambdaw =(density_w.rho/muw)*psk.dkrw;
    
    lambdan  =(density_n.rho*psk.krn)/mun;
    dlambdan =(density_n.rho/mun)*psk.dkrn;
    
    lambdat  = lambdaw + lambdan;
    dlambdat = dlambdaw + dlambdan;
    
    fw      = lambdaw/lambdat;
    dfw     = (dlambdaw*lambdat - lambdaw*dlambdat)/(lambdat*lambdat);
    
    fn      = lambdan/lambdat;
    dfn     = (dlambdan*lambdat - lambdan*dlambdat)/(lambdat*lambdat);
  }
};

struct CompressibleN_FractionalFlowVariables : public FractionalFlowVariables
{
public:
  CompressibleN_FractionalFlowVariables(double muwIn,double munIn):
    FractionalFlowVariables(muwIn,munIn)
  {}
  
  double drhon,
    dlambdaw_psiw, 
    drhon_psiw,
    dlambdan_psiw,
    dlambdat_psiw,
    dfw_psiw,
    dfn_psiw;
  
  inline void calc(const PskRelation& psk, 
                   const DensityRelation& density_w,
                   const DensityRelation& density_n)
  {  
    lambdaw       = density_w.rho*psk.krw/muw; 
    dlambdaw      = (density_w.rho/muw)*psk.dkrw;
    dlambdaw_psiw = 0.0; 
    
    drhon         = density_n.drho*psk.dpsic;
    drhon_psiw    = density_n.drho;
    
    lambdan       = (density_n.rho*psk.krn)/mun;
    dlambdan      = (1.0/mun)*(psk.dkrn*density_n.rho + drhon*psk.krn);
    dlambdan_psiw = drhon_psiw*(psk.krn/mun);
    
    lambdat       = lambdaw + lambdan;
    dlambdat      = dlambdaw + dlambdan;
    dlambdat_psiw = dlambdaw_psiw + dlambdan_psiw;
    
    fw       = lambdaw/lambdat;
    dfw      = (dlambdaw*lambdat - lambdaw*dlambdat)/(lambdat*lambdat);
    dfw_psiw = (dlambdaw_psiw*lambdat - lambdaw*dlambdat_psiw)/(lambdat*lambdat);
    
    fn       = lambdan/lambdat;
    dfn      = (dlambdan*lambdat - lambdan*dlambdat)/(lambdat*lambdat);
    dfn_psiw = (dlambdan_psiw*lambdat - lambdan*dlambdat_psiw)/(lambdat*lambdat);
  }
};


class PskSpline: public PskRelation
{
  /*************************************************************
     spline psk relations. tables of values are held (externally) 
     in splineArray in the order
     
     u (sw, or psic), u^-1 (either sw or psic), krw,krn

     force dSe_dSw = 1.0, so that will be consistent with analytical evaluations
      where dkrw, etc are wrt to se

   ************************************************************/
public:
  PskSpline(const double* rworkIn, const int* iworkIn=0):
    PskRelation(),
    nknots(2),
    lastIndex(0),
    uinvOffset(1),
    krwOffset(2),
    krnOffset(3),
    splineArray(rworkIn)
  {
    assert(iworkIn);
    nknots = iworkIn[0];
  }
  virtual ~PskSpline() {}
  virtual inline void calc(const double& Sw)
  {
    assert(splineArray);
    //
    Se = Sw; //note force Se = Sw since splines are evaluated directly
    dSe_dSw = 1.0;
    piecewiseLinearTableLookup(Sw,
			       nknots,
			       &lastIndex,
			       &psic,
			       &dpsic,
			       splineArray,
			       splineArray+uinvOffset*nknots);
    
    //
    piecewiseLinearTableLookup(Sw,
			       nknots,
			       &lastIndex,
			       &krw,
			       &dkrw,
			       splineArray,
			       splineArray+krwOffset*nknots);
    //
    piecewiseLinearTableLookup(Sw,
			       nknots,
			       &lastIndex,
			       &krn,
			       &dkrn,
			       splineArray,
			       splineArray+krnOffset*nknots);

  }
  virtual inline void calc_from_psic(const double& psicIn)
  {
    assert(splineArray);
    //
    psic = psicIn;
    piecewiseLinearTableLookup(psic,
			       nknots,
			       &lastIndex,
			       &Se, //same as Sw
			       &dSe_dpsic, //same as dSw_dpsic
			       splineArray,
			       splineArray+uinvOffset*nknots);
    
    //
    piecewiseLinearTableLookup(psic,
			       nknots,
			       &lastIndex,
			       &krw,
			       &dkrw,
			       splineArray,
			       splineArray+krwOffset*nknots);
    //
    piecewiseLinearTableLookup(psic,
			       nknots,
			       &lastIndex,
			       &krn,
			       &dkrn,
			       splineArray,
			       splineArray+krnOffset*nknots);

  }
  virtual inline void setParams(const double* rwork, const int* iwork = 0)
  {
    splineArray = rwork;
    if (iwork)
      nknots = iwork[0];
  }
public:
  int nknots,lastIndex,uinvOffset,krwOffset,krnOffset;
  const double* splineArray; 
};
/** @} */
#endif

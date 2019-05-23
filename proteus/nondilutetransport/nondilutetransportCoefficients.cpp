#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <memory.h>
#include <stdlib.h>
#include "TCAT_ND.h"
#include "nondilutetransportCoefficients.h"

void ConMassFluidEvaluate(const int nPoints,
                          const int nSpace,
                          const double K,
                          const double grav,
                          const double dt,
                          const double poro,
                          const double L,
                          const double *u,
                          double *m,
                          double *dm,
                          double *f,
                          double *df,
                          double *phi,
                          double *dphi,
                          double *a,
                          double *da,
                          double *r,
                          double *x,
                          double *mass_frac,
                          double *mass_frac_old)
{
  int k,I;
  const int nSpace2=nSpace*nSpace;
  double density,d_density,mu,density_old,hold;


  for (k=0;k<nPoints;k++){

      density = den(mass_frac[k]);
      density_old = den(mass_frac_old[k]);
      mu = visc(mass_frac[k]);

      m[k] = 0.0;
      dm[k] = 0.0;
      f[k] = 0.0;
      df[k] = 0.0;

      phi[k] = u[k] - density*grav*(x[k*3] - L);
      dphi[k] = 1.0;
      if (dt < 1.e-12){
		r[k] = 0.0;
	  }
	  else{
      r[k] = poro/dt*(density - density_old);
	 }

      for (I=0;I<nSpace;I++){
          a[k*nSpace2+I*nSpace+I] = K/mu*density;
	  }

   }
}


void NonDiluteDispersionEvaluate(const int nPoints,
                          const int nSpace,
                          const double poro,
                          const double diff,
                          const double alpha_L,
                          const double R,
                          const double theta,
                          const double MW_a,
                          const double MW_b,
                          const double beta1,
                          const double beta2,
                          const double *w,
                          const double *act,
                          double *m,
                          double *dm,
                          double *f,
                          double *df,
                          double *r,
                          double *dr0,
                          double *dr1,
                          double *a00,
                          double *da000,
                          double *a01,
                          double *da010,
                          double *da011,
                          double *velocity,
                          double *pressure,
                          double *w_old,
                          const double *grad_w,
                          const double *grad_act,
                          double *phi0,
                          double *dphi00,
                          double *dphi01,
                          double *phi1,
                          double *dphi10,
                          double *dphi11,
                          double *x)
{
  int k,I;
  const int nSpace2=nSpace*nSpace;
  double act_w,d_act_w,d2_act_w,molal_w,d_molal_w,density,d_density,mu,d_mu;
  double Disp,alpha_T,MW_w,d_MW_w,x_a,d_x_a;
  double arg,d_arg,d_alpha_T;
  double arg2;

  double scale = 1.e-15;
  struct TCAT_var TCAT_v;
  TCAT_v.MW_a = MW_a;
  TCAT_v.MW_b = MW_b;

  for (k=0;k<nPoints;k++){

      molal_w = molal(w[k]);
	    d_molal_w = d_molal(w[k]);
      act_w = activity(w[k],molal_w);
      d_act_w = d_act(w[k],molal_w,d_molal_w,act_w);
      d2_act_w = d2_act(w[k],molal_w,d_molal_w,act_w,d_act_w);

      density = den(w[k]);
      d_density = d_den(w[k]);
      mu = visc(w[k]);
      d_mu = visc(w[k]);
	    x_a = mole_frac(w[k],TCAT_v);
	    d_x_a = d_mole_frac(w[k],TCAT_v);
      MW_w = mol_weight(x_a,TCAT_v);
	    d_MW_w = d_mol_weight(d_x_a,TCAT_v);

      velocity[k] = velocity[k]/den(w_old[k]);
      Disp = poro*alpha_L*velocity[k];

    if(d_act_w > 0.0){
		arg = Disp*density*grad_w[k]*(beta1 + beta2*w[k]/act_w*MW_b/MW_w*d_act_w);
        arg2 = Disp*density*grad_w[k]*beta2;
    	d_arg = d_density*Disp*grad_w[k]*(beta1 + beta2*w[k]/act_w*MW_b/MW_w*d_act_w)
    	        + arg2*1.0/act_w*MW_b/MW_w*d_act_w
   	            - arg2*w[k]*d_act_w/(act_w*act_w)*MW_b/MW_w*d_act_w
  	            - arg2*w[k]/act_w*MW_b*(d_MW_w)/(MW_w*MW_w)*d_act_w
  	            + arg2*w[k]/act_w*MW_b/MW_w*d2_act_w ;
	}
	else{
		arg = Disp*density*grad_w[k]*beta1;
    	d_arg = d_density*Disp*grad_w[k]*beta1;
	}

      if ( (1.0 - arg) < 0.0){
			arg = 0.0;
			d_arg = 0.0;
			alpha_T = alpha_L;
            d_alpha_T = 0.0;
      }
      else{
      		arg = sqrt(1. - arg);
      		d_arg = -d_arg/(2.*arg);

      		alpha_T = 2.*alpha_L/(1. + arg);
      		d_alpha_T = -2.*alpha_L*d_arg/((1.+arg)*(1.+arg));
      }


      m[k] = w[k]*poro*density;
      dm[k] = poro*(w[k]*d_density + density);

      phi0[k] = w[k];
      dphi00[k] = 1.0;
      phi1[k] = act_w;
      dphi10[k] = d_act_w;

      f[k] = 0.0;
      df[k] = 0.0;

      r[k] = scale*(act[k] - act_w);
      dr0[k] = -scale*d_act_w;
      dr1[k] =  scale;

      for (I=0;I<nSpace;I++){
          Disp = diff+alpha_T*velocity[k*nSpace+I];
          a00[k*nSpace2+I*nSpace+I] = density*poro*Disp;
          da000[k*nSpace2+I*nSpace+I] = d_density*poro*Disp + density*poro*velocity[k*nSpace+I]*d_alpha_T;


          a01[k*nSpace2+I*nSpace+I] = poro*density*w[k]/act_w*MW_b/MW_w*Disp;
           da010[k*nSpace2+I*nSpace+I] = poro*d_density*w[k]/act_w*MW_b/MW_w*Disp
                                       + poro*density*1.0/act_w*MW_b/MW_w*Disp
                                       - poro*density*w[k]*d_act_w/(act_w*act_w)*MW_b/MW_w*Disp
                                       - poro*density*w[k]/act_w*MW_b*d_MW_w/(MW_w*MW_w)*Disp
                                       + poro*density*w[k]/act_w*MW_b/MW_w*d_alpha_T*velocity[k*nSpace+I];

	  }
   }
}


void NonDilutePhiDispersionEvaluate(const int nPoints,
                          const int nSpace,
                          const double poro,
                          const double diff,
                          const double alpha_L,
                          const double R,
                          const double theta,
                          const double MW_a,
                          const double MW_b,
                          const double beta1,
                          const double beta2,
                          const double *w,
                          const double *act,
                          double *phi0,
                          double *dphi00,
                          double *dphi01,
                          double *phi1,
                          double *dphi10,
                          double *dphi11,
                          double *f)
{
  int k,I;
  const int nSpace2=nSpace*nSpace;
  double act_w,d_act_w,molal_w,d_molal_w;

  for (k=0;k<nPoints;k++){

	molal_w = molal(w[k]);
	d_molal_w = d_molal(w[k]);
    act_w = activity(w[k],molal_w);
    d_act_w = d_act(w[k],molal_w,d_molal_w,act_w);
    phi0[k] = w[k];
    dphi00[k] = 1.0;
    phi1[k] = act_w;
    dphi10[k] = d_act_w;


   }
}





void NonDiluteEvaluate(const int nPoints,
                          const int nSpace,
                          const double poro,
                          const double diff,
                          const double alpha_L,
                          const double R,
                          const double theta,
                          const double MW_a,
                          const double MW_b,
                          const double beta1,
                          const double beta2,
                          const double *w,
                          const double *act,
                          double *m,
                          double *dm,
                          double *f,
                          double *df,
                          double *r,
                          double *dr0,
                          double *dr1,
                          double *a00,
                          double *da000,
                          double *a01,
                          double *da010,
                          double *da011,
                          double *velocity,
                          double *pressure,
                          const double *grad_w,
                          const double *grad_act,
                          double *phi0,
                          double *dphi00,
                          double *dphi01,
                          double *phi1,
                          double *dphi10,
                          double *dphi11,
                          double *x)
{
  int k,I;
  const int nSpace2=nSpace*nSpace;
  double act_w,d_act_w,d2_act_w,molal_w,d_molal_w,density,d_density,mu,d_mu;
  double Disp,alpha_T,MW_w,d_MW_w,x_a,d_x_a;
  double arg,d_arg,d_alpha_T;

  double scale = 1.e-15;
  struct TCAT_var TCAT_v;
  TCAT_v.MW_a = MW_a;
  TCAT_v.MW_b = MW_b;

  for (k=0;k<nPoints;k++){

	  molal_w = molal(w[k]);
	  d_molal_w = d_molal(w[k]);
    act_w = activity(w[k],molal_w);
    d_act_w = d_act(w[k],molal_w,d_molal_w,act_w);
    d2_act_w = d2_act(w[k],molal_w,d_molal_w,act_w,d_act_w);

    density = den(w[k]);
    d_density = d_den(w[k]);
    mu = visc(w[k]);
    d_mu = visc(w[k]);
	  x_a = mole_frac(w[k],TCAT_v);
	  d_x_a = d_mole_frac(w[k],TCAT_v);
    MW_w = mol_weight(x_a,TCAT_v);
	  d_MW_w = d_mol_weight(d_x_a,TCAT_v);

    velocity[k] = velocity[k]/density;
    Disp = poro*alpha_L*velocity[k];

    if(d_act_w > 0.0){
		  arg = Disp*density*grad_w[k]*(beta1 + beta2*w[k]/act_w*MW_b/MW_w*d_act_w);
    	d_arg = d_density*Disp*grad_w[k]*(beta1 + beta2*w[k]/act_w*MW_b/MW_w*d_act_w)
    	        + Disp*density*grad_w[k]*beta2*1.0/act_w*MW_b/MW_w*d_act_w
   	          - Disp*density*grad_w[k]*beta2*w[k]*d_act_w/(act_w*act_w)*MW_b/MW_w*d_act_w
  	          - Disp*density*grad_w[k]*beta2*w[k]/act_w*MW_b*(d_MW_w)/(MW_w*MW_w)*d_act_w
  	          + Disp*density*grad_w[k]*beta2*w[k]/act_w*MW_b/MW_w*d2_act_w ;
	  }
	  else{
		  arg = Disp*density*grad_w[k]*beta1;
    	d_arg = d_density*Disp*grad_w[k]*beta1;
	  }


      if (arg > 1.0){
			arg = 0.0;
			d_arg = 0.0;
			alpha_T = alpha_L;
            d_alpha_T = 0.0;
      }
      else{
      		arg = sqrt(1. - arg);
      		d_arg = -d_arg/(2.*arg);

      		alpha_T = 2.*alpha_L/(1. + arg);
      		d_alpha_T = -2.*alpha_L*d_arg/((1.+arg)*(1.+arg));
      }


      m[k] = w[k]*poro*density;
      dm[k] = poro*(w[k]*d_density + density);

      phi0[k] = w[k];
      dphi00[k] = 1.0;
      phi1[k] = act_w;
      dphi10[k] = d_act_w;

      f[k] = poro*density*velocity[k]*w[k];
      df[k] = poro*velocity[k]*(w[k]*d_density + density);

      //printf("%e %e \n",velocity[k],w[k]);

      r[k] = scale*(act[k] - act_w);
      dr0[k] = -scale*d_act_w;
      dr1[k] =  scale;

      for (I=0;I<nSpace;I++){
          Disp = diff+alpha_T*velocity[k*nSpace+I];
          a00[k*nSpace2+I*nSpace+I] = density*poro*Disp;
          da000[k*nSpace2+I*nSpace+I] = d_density*poro*Disp + density*poro*velocity[k*nSpace+I]*d_alpha_T;


          a01[k*nSpace2+I*nSpace+I] = poro*density*w[k]/act_w*MW_b/MW_w*Disp;
           da010[k*nSpace2+I*nSpace+I] = poro*d_density*w[k]/act_w*MW_b/MW_w*Disp
                                       + poro*density*1.0/act_w*MW_b/MW_w*Disp
                                       - poro*density*w[k]*d_act_w/(act_w*act_w)*MW_b/MW_w*Disp
                                       - poro*density*w[k]/act_w*MW_b*d_MW_w/(MW_w*MW_w)*Disp
                                       + poro*density*w[k]/act_w*MW_b/MW_w*d_alpha_T*velocity[k*nSpace+I];

	  }
   }
}


void NonDilutePhiEvaluate(const int nPoints,
                          const int nSpace,
                          const double poro,
                          const double diff,
                          const double alpha_L,
                          const double R,
                          const double theta,
                          const double MW_a,
                          const double MW_b,
                          const double beta1,
                          const double beta2,
                          const double *w,
                          const double *act,
                          double *phi0,
                          double *dphi00,
                          double *dphi01,
                          double *phi1,
                          double *dphi10,
                          double *dphi11,
                          double *f)
{
  int k,I;
  const int nSpace2=nSpace*nSpace;
  double act_w,d_act_w,molal_w,d_molal_w,density,d_density,mu,d_mu,Disp,alpha_T,mole_frac,MW_w,MW_w3;
  double arg,d_arg,d_alpha_T;


  for (k=0;k<nPoints;k++){

	  molal_w = molal(w[k]);
	  d_molal_w = d_molal(w[k]);
      act_w = activity(w[k],molal_w);
      d_act_w = d_act(w[k],molal_w,d_molal_w,act_w);

     phi0[k] = w[k];
     dphi00[k] = 1.0;
     phi1[k] = act_w;
     dphi10[k] = d_act_w;


   }
}





void AdvectionEvaluate(const int nPoints,
                          const int nSpace,
                          const double poro,
                          const double diff,
                          const double alpha_L,
                          const double *u,
                          double *m,
                          double *dm,
                          double *f,
                          double *df,
                          double *a,
                          double *da,
                          double *velocity)
{
  int k,I;
  const int nSpace2=nSpace*nSpace;
  double density,d_density;

  for (k=0;k<nPoints;k++){

      density = den(u[k]);
      velocity[k] = velocity[k]/density;
      d_density = d_den(u[k]);
      m[k] = u[k]*density*poro;
      dm[k] = poro*(u[k]*d_density + density);
      f[k] = poro*density*velocity[k]*u[k];
      df[k] = poro*velocity[k]*(u[k]*d_density + density);
   }

}



















void NonDiluteTESTEvaluate(const int nPoints,
                          const int nSpace,
                          const double poro,
                          const double grav,
                          const double K,
                          const double L,
                          const double diff,
                          const double alpha_L,
                          const double R,
                          const double theta,
                          const double MW_a,
                          const double MW_b,
                          const double beta1,
                          const double beta2,
                          const double *press,
                          const double *w,
                          double *m0,
                          double *m1,
                          double *dm01,
                          double *dm11,
                          double *phi0,
                          double *phi1,
                          double *dphi00,
                          double *dphi01,
                          double *dphi10,
                          double *dphi11,
                          double *f1,
                          double *df10,
                          double *df11,
                          double *a00,
                          double *a10,
                          double *a11,
                          double *da001,
                          double *da101,
                          double *da111,
                          const double *x,
                          const double *grad_press)
{
  int k,I;
  const int nSpace2=nSpace*nSpace;
  double act_w,d_act_w,d2_act_w,molal_w,d_molal_w,density,d_density,mu,d_mu;
  double Disp,d_Disp,alpha_T,mole_frac,MW_w,d_MW_w;
  double arg,d_arg,d_alpha_T,velocity,d_velocity;

  double scale = 1.e-4;
  struct TCAT_var TCAT_v;
  TCAT_v.MW_a = MW_a;
  TCAT_v.MW_b = MW_b;

  for (k=0;k<nPoints;k++){

	   molal_w = molal(w[k]);
       d_molal_w = d_molal(w[k]);
       act_w = activity(w[k],molal_w);
       d_act_w = d_act(w[k],molal_w,d_molal_w,act_w);
       d2_act_w = d2_act(w[k],molal_w,d_molal_w,act_w,d_act_w);

       density = den(w[k]);
       d_density = d_den(w[k]);
       mu = visc(w[k]);
       d_mu = visc(w[k]);
       MW_w = mol_weight(w[k],TCAT_v);
       d_MW_w = d_mol_weight(w[k],TCAT_v);

//       Disp = poro*alpha_L*velocity[k];
//       arg = Disp*density*grad_w[k]*(beta1 + beta2*w[k]/act_w*MW_b/MW_w*d_act_w);
//       d_arg = d_density*Disp*grad_w[k]*(beta1 - beta2*w[k]/act_w*MW_b/MW_w*d_act_w)
//             + Disp*density*grad_w[k]*beta2*1.0/act_w*MW_b/MW_w*d_act_w
//             - Disp*density*grad_w[k]*beta2*w[k]*d_act_w/(act_w*act_w)*MW_b/MW_w*d_act_w
//             + Disp*density*grad_w[k]*beta2*w[k]/act_w*MW_b/MW_w*d2_act_w ;
//       if (arg > 1.0){
// 			arg = 0.0;
// 			d_arg = 0.0;
// 			alpha_T = alpha_L;
//             d_alpha_T = 0.0;
//       }
//       else{
//       		arg = sqrt(1. - arg);
//       		d_arg = -d_arg/(2.*arg);
//       		alpha_T = 2.*alpha_L/(1. + arg);
//       		d_alpha_T = -2.*alpha_L*d_arg/((1.+arg)*(1.+arg));
//       }



      // FLOW EQUATION //
      m0[k] = density;
      dm01[k] = d_density;
      phi0[k] = press[k] - density*grav*(x[k*3] - L);
      dphi00[k] = 1.0;
      dphi01[k] = -d_density*grav*(x[k*3] - L);


      for (I=0;I<nSpace;I++){
          a00[k*nSpace2+I*nSpace+I] = K/mu*density;
          da001[k*nSpace2+I*nSpace+I] = K*(mu*d_density + density*d_mu)/(mu*mu);
	  }


      velocity = grad_press[k]/density;//0.00577256254903116;//-K/mu*(grad_press[k] - density*grav);
      d_velocity = 0.0;//K*( d_mu*(grad_press[k] - density*grav) + d_density*grav*mu)/(mu*mu);

      //printf("%e %e %e\n",velocity,d_velocity,w[k]);

      // TRANSPORT EQUATION //
      m1[k] = w[k]*density*poro;
      dm11[k] = poro*(w[k]*d_density + density);
      phi1[k] = w[k];
      dphi11[k] = 1.0;
      //f1[k] = density*velocity*w[k]*poro;
      //df10[k] = 1.0;
      //df11[k] = poro*(d_density*velocity*w[k] + d_velocity*density*w[k] +  density*velocity);

     //printf("%e %e %e %e\n",grad_press[k],density,grad_press[k]/density,w[k]);

       for (I=0;I<nSpace;I++){
           Disp = diff+alpha_L*velocity;
           d_Disp = alpha_L*d_velocity;

           a10[k*nSpace2+I*nSpace+I] = K/mu*density*w[k];
           da101[k*nSpace2+I*nSpace+I] = K*(mu*d_density + density*d_mu)/(mu*mu)*w[k] + K/mu*density ;

           a11[k*nSpace2+I*nSpace+I] = density*Disp;
           da111[k*nSpace2+I*nSpace+I] = d_density*Disp + density*d_Disp;  //+ density*poro*velocity[k*nSpace+I]*d_alpha_T;


//           //a01[k*nSpace2+I*nSpace+I] = poro*density*w[k]/act_w*MW_b/MW_w*Disp;
//           //da010[k*nSpace2+I*nSpace+I] = poro*d_density*w[k]/act_w*MW_b/MW_w*Disp
//           //                            + poro*density*1.0/act_w*MW_b/MW_w*Disp
//           //                            - poro*density*w[k]*d_act_w/(act_w*act_w)*MW_b/MW_w*Disp
//           //                            - poro*density*w[k]/act_w*MW_b*d_MW_w/(MW_w*MW_w)*Disp
//           //                            + poro*density*w[k]/act_w*MW_b/MW_w*d_alpha_T*velocity[k*nSpace+I];
//
 	  }
   }
}


void NonDilutePhiTESTEvaluate(const int nPoints,
                          const int nSpace,
                          const double poro,
                          const double grav,
                          const double L,
                          const double *press,
                          const double *w,
                          double *phi0,
                          double *phi1,
                          double *dphi00,
                          double *dphi01,
                          double *dphi10,
                          double *dphi11,
                          double *f,
                          const double *x)
{
  int k,I;
  const int nSpace2=nSpace*nSpace;
  double density,d_density;


  for (k=0;k<nPoints;k++){

      density   = den(w[k]);
      d_density = d_den(w[k]);

      phi0[k]   = press[k] - density*grav*(x[k*3] - L);
      dphi00[k] = 1.0;
      dphi01[k] = -d_density*grav*(x[k*3] - L);
      phi1[k]   = w[k];
      dphi11[k] = 1.0;


   }
}

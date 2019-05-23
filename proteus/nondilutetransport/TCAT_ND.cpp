#include <math.h>
#include <iostream>

double grav = -980.;
double eps = 1.e-6;

struct TCAT_var
{
	double T = 298.15; // Temperature in Kelvin
	double R = 8.314e+07; //Universal Gas Constant
	double MW_a = 199.88; // MW CaBr_2
	double MW_b = 18.015; // MW Water
	double poro;
	double perm;
	double diff;
	double alpha_L;
	double beta1;
	double beta2;
};


double den(double x){
 double d,x2,x3;
 double a0 = 0.9971;
 double a1 = 0.83898;
 double a2 = 0.48132;
 double a3 = 0.86154;
 if (x < eps){
	d = a0;
 }
 else{
 	x2 = x*x; x3 = x2*x;
 	d = a0 + a1*x + a2*x2 + a3*x3;
 }

 return d;
}


double d_den(double x){
 double dd;
 double a0 = 0.9971;
 double a1 = 0.83898;
 double a2 = 0.48132;
 double a3 = 0.86154;

  dd = a1 + x*(2.0*a2 + 3.0*a3*x);

 return dd;
}

double den_root(double x){
 double r,x2,x3;
 double a0 = 0.9971;
 double a1 = 0.83898;
 double a2 = 0.48132;
 double a3 = 0.86154;
  x2 = x*x; x3 = x2*x;
  r = (-a0 + sqrt(4.*a1*x + a0*a0))/(2.*a1);
 return r;
}


double visc(double x){
 double mu,eval,x2,x3;
 double a0 = -0.1165;
 double a1 = 1.318;
 double a2 = -2.636;
 double a3 = 11.49;
 double conv = 0.01;

 	x2 = x*x; x3 = x2*x;
 	eval = a0 + a1*x + a2*x2 + a3*x3;
 	mu = conv*exp(eval);

 return mu;
}

double d_visc(double x)
{
 double d_mu,eval,d_eval,x2,x3;
 double a0 = -0.1165;
 double a1 = 1.318;
 double a2 = -2.636;
 double a3 = 11.49;
 double conv = 0.01;


 	x2 = x*x; x3 = x2*x;
 	eval = a3*x3 + a2*x2 + a1*x + a0;
 	d_eval = 3.*a3*x2 + 2.*a2*x + a1;
 	d_mu = conv*exp(eval)*d_eval;

 return d_mu;
}


double mole_frac(double x, struct TCAT_var TCAT_v){
 double m_f;
 double MW_a = TCAT_v.MW_a;
 double MW_b = TCAT_v.MW_b;
 m_f = (MW_b*x)/(MW_b*x + MW_a*(1.-x));
 return m_f;
}

double d_mole_frac(double x, struct TCAT_var TCAT_v){
 double d_m_f;
 double MW_a = TCAT_v.MW_a;
 double MW_b = TCAT_v.MW_b;
 d_m_f = MW_a*MW_b/((MW_b*x - MW_a*(x-1.))*(MW_b*x - MW_a*(x-1.)));
 return d_m_f;
}


double mol_weight(double x, struct TCAT_var TCAT_v){
 double m_w;
 double MW_a = TCAT_v.MW_a;
 double MW_b = TCAT_v.MW_b;
 m_w = MW_a*x + MW_b*(1.-x);
 return m_w;
}

double d_mol_weight(double d_x, struct TCAT_var TCAT_v){
 double d_m_w;
 double MW_a = TCAT_v.MW_a;
 double MW_b = TCAT_v.MW_b;
 d_m_w = (MW_a-MW_b)*d_x;
 return d_m_w;
}

double molal(double x)
{
 double MW_a = 199.88; //CaBr_2
 double MW_b = 18.015; //WATER
 double mol;
 mol = 1000.*x/((1.-x)*MW_a);
 return mol;
}

double d_molal(double x)
{
 double MW_a = 199.88; //CaBr_2
 double MW_b = 18.015; //WATER
 double d_mol;
 d_mol = 1000./((1.-x)*(1.-x)*MW_a);
 return d_mol;
}

double d2_molal(double x)
{
 double MW_a = 199.88; //CaBr_2
 double MW_b = 18.015; //WATER
 double d2_mol;
 d2_mol = 1000./MW_a *( (2.*x)/((1.-x)*(1.-x)*(1.-x)) + 2./((1.-x)*(1.-x)));
 return d2_mol;
}


double activity(double x,double m)
{
	double a;
	double I,I_hp,m2,m3,m4,m5,m6;
    double ac0 = 2.3525;
    double ac1 = 1.793589;
    double ac2 = 3.244255e-1;
    double ac3 = 2.086318e-1;
    double ac4 = -5.6619089e-2;
    double ac5 = 1.21564998e-2;
    double ac6 = -1.2930306e-3;
    double ac7 = 4.84942655e-5;
	if(x < eps){
		a = 1.0;
	}
	else{
    I = 3.*m;
    I_hp = sqrt(I);
    m2 = m*m;
    m3 = m2*m;
    m4 = m3*m;
    m5 = m4*m;
    m6 = m5*m;
    a = exp( -ac0*I_hp/(1.+ac1*I_hp) + ac2*m + ac3*m2 + ac4*m3 + ac5*m4 + ac6*m5 + ac7*m6 );
	}
	return a;
}

double d_act(double x,double m, double dm, double act)
{
	double arg,deriv,d_ax;
	double I,dI,I_hp,m2,m3,m4,m5,m6;
    double ac0 = 2.3525;
    double ac1 = 1.793589;
    double ac2 = 3.244255e-1;
    double ac3 = 2.086318e-1;
    double ac4 = -5.6619089e-2;
    double ac5 = 1.21564998e-2;
    double ac6 = -1.2930306e-3;
    double ac7 = 4.84942655e-5;

	if(x < eps){
		d_ax = 0.0;
	}
	else{
    I = 3.*m;
    dI = 3.*dm;
    I_hp = sqrt(I);

    m2 = m*m;
    m3 = m2*m;
    m4 = m3*m;
    m5 = m4*m;
    m6 = m5*m;

    deriv = -ac0*dI/((2.0*I_hp)*(ac1*I_hp + 1.0)*(ac1*I_hp + 1.0)) + (ac2 + 2.0*ac3*m + 3.0*ac4*m2 + 4.0*ac5*m3 + 5.0*ac6*m4 + 6.0*ac7*m5)*dm;
    d_ax = act*deriv;
	}
	return d_ax;
}


double d2_act(double x,double m, double dm, double act, double d_act)
{
	double arg,deriv,deriv_2,d2_ax,dI2,dm2;
    double num,denom;
	double I,dI,I_hp,m2,m3,m4,m5,m6;
    double ac0 = 2.3525;
    double ac1 = 1.793589;
    double ac2 = 3.244255e-1;
    double ac3 = 2.086318e-1;
    double ac4 = -5.6619089e-2;
    double ac5 = 1.21564998e-2;
    double ac6 = -1.2930306e-3;
    double ac7 = 4.84942655e-5;

	if(x < eps){
		d2_ax = 0.0;
	}
	else{
    	I = 3.*m;
    	dI = 3.*dm;
        dm2 = d2_molal(x);
        dI2 = 3.*dm2;
    	I_hp = sqrt(I);

   		m2 = m*m;
    	m3 = m2*m;
    	m4 = m3*m;
    	m5 = m4*m;
    	m6 = m5*m;

        arg = ac2 + 2.0*ac3*m + 3.0*ac4*m2 + 4.0*ac5*m3 + 5.0*ac6*m4 + 6.0*ac7*m5;
   		deriv = -ac0*dI/((2.0*I_hp)*(ac1*I_hp + 1.0)*(ac1*I_hp + 1.0)) + arg*dm;
        num = ac0*((3.*ac1*I_hp + 1.)*dI*dI - 2.*I*(ac1*I_hp + 1.)*dI2);
        denom = 4.*I_hp*I*(ac1*I_hp + 1.)*(ac1*I_hp + 1.)*(ac1*I_hp + 1.);
    	deriv_2 = num/denom + dm2*arg + 2.*dm*dm*( ac3 + 3.*ac4*m + 6.*ac5*m2 + 10.*ac6*m3 + 15.*ac7*m4);
    	d2_ax = act*(deriv_2 + deriv*deriv);
	}
	return d2_ax;
}



double  ENTROPY_TCAT(double w, double grad_w, double grad_p, struct TCAT_var TCAT_v){

	double act_w,d_act_w,d2_act_w,molal_w,d_molal_w,density,d_density,mu,d_mu;
	double Disp,alpha_T,alpha_L,MW_w,d_MW_w,x_a;
	double arg,d_arg,d_alpha_T;
  double T,R,beta1,beta2,diff,D_u,MW_a,MW_b,poro,velocity,perm;
  double ent,term1,term2,term3,term4;

	poro = TCAT_v.poro;
	perm = TCAT_v.perm;
	diff = TCAT_v.diff;
	alpha_L = TCAT_v.alpha_L;
	MW_a = TCAT_v.MW_a;
	MW_b = TCAT_v.MW_b;
	T = TCAT_v.T;
	R = TCAT_v.R;
	beta1 = TCAT_v.beta1;
	beta2 = TCAT_v.beta2;

	molal_w = molal(w);
	d_molal_w = d_molal(w);
	act_w = activity(w,molal_w);
	d_act_w = d_act(w,molal_w,d_molal_w,act_w);

	density = den(w);
	d_density = d_den(w);
	mu = visc(w);
	x_a = mole_frac(w,TCAT_v);
	MW_w = mol_weight(w,TCAT_v);

  velocity = -perm/mu*(grad_p - density*grav)/poro;
	Disp = poro*alpha_L*velocity;

  if(d_act_w > 0.0){
		arg = Disp*density*grad_w*(beta1 + beta2*w/act_w*MW_b/MW_w*d_act_w);
	}
	else{
		arg = Disp*density*grad_w*beta1;
	}

	if (arg > 1.0){
		arg = 0.0;
		alpha_T = alpha_L;
	}
	else{
		arg = sqrt(1. - arg);
		alpha_T = 2.*alpha_L/(1. + arg);
	}

  Disp = diff+alpha_T*velocity;
	term1 = mu/(perm*T)*(poro*poro*velocity*velocity);
  term2 = poro/(R*T*T)*density*w*(1.-w)*(MW_a*MW_b/MW_w)*Disp;

  if (w<1.e-9){
		term3 = 0.0;
    term4 = 0.0;
  }
  else{
    term3 = R*T/(w*(1.-w))*MW_w/(MW_a*MW_b)*grad_w;
    term4 = R*T/(act_w*MW_a*(1.-w))*grad_w*d_act_w;
  }

	ent = term1 + term2*(term3+term4)*(term3+term4);
	return ent;
}

double DENTROPY_TCAT(double w, double grad_w, double grad_p, struct TCAT_var TCAT_v) {

	double act_w,d_act_w,d2_act_w,molal_w,d_molal_w,density,d_density,mu,d_mu;
	double Disp,alpha_T,alpha_L,MW_w,d_MW_w,x_a,d_x_a;
	double arg,d_arg,d_alpha_T;
    double T,R,beta1,beta2,diff,D_u,MW_a,MW_b,poro,velocity,perm;
    double d_ent,term1,term2,term3,term4;
    double d_term1,d_term2,d_term3,d_term4;

	poro = TCAT_v.poro;
	perm = TCAT_v.perm;
	diff = TCAT_v.diff;
	alpha_L = TCAT_v.alpha_L;
	MW_a = TCAT_v.MW_a;
	MW_b = TCAT_v.MW_b;
	T = TCAT_v.T;
	R = TCAT_v.R;
	beta1 = TCAT_v.beta1;
	beta2 = TCAT_v.beta2;

	molal_w = molal(w);
	d_molal_w = d_molal(w);
	act_w = activity(w,molal_w);
	d_act_w = d_act(w,molal_w,d_molal_w,act_w);
    d2_act_w = d2_act(w,molal_w,d_molal_w,act_w,d_act_w);

	density = den(w);
	d_density = d_den(w);
	mu = visc(w);
	d_mu = d_visc(w);
	x_a = mole_frac(w,TCAT_v);
	d_x_a = d_mole_frac(w,TCAT_v);
	MW_w = mol_weight(x_a,TCAT_v);
	d_MW_w = d_mol_weight(d_x_a,TCAT_v);

    velocity = -perm/mu*(grad_p - density*grav)/poro;
	Disp = poro*alpha_L*velocity;

    if(d_act_w > 0.0){
		arg = Disp*density*grad_w*(beta1 + beta2*w/act_w*MW_b/MW_w*d_act_w);
    	d_arg = d_density*Disp*grad_w*(beta1 + beta2*w/act_w*MW_b/MW_w*d_act_w)
    	        + Disp*density*grad_w*beta2*1.0/act_w*MW_b/MW_w*d_act_w
   	            - Disp*density*grad_w*beta2*w*d_act_w/(act_w*act_w)*MW_b/MW_w*d_act_w
  	            - Disp*density*grad_w*beta2*w/act_w*MW_b*(d_MW_w)/(MW_w*MW_w)*d_act_w
  	            + Disp*density*grad_w*beta2*w/act_w*MW_b/MW_w*d2_act_w ;
	}
	else{
		arg = Disp*density*grad_w*beta1;
    	d_arg = d_density*Disp*grad_w*beta1;
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

    Disp = diff+alpha_T*velocity;

	d_term1 = d_mu/(perm*T)*(poro*poro*velocity*velocity);

    term2 = MW_a*MW_b*poro/(R*T*T)*( Disp * density * w*(1.-w) * (1./MW_w) );
    d_term2 = MW_a*MW_b*poro/(R*T*T)*(   d_alpha_T*velocity * density * w*(1.-w) * (1./MW_w)
                                       + Disp * d_density * w*(1.-w) * (1./MW_w)
                                       + Disp * density * (1. - 2*w) * (1./MW_w)
                                       - Disp * density * w*(1.-w) * d_MW_w/(MW_w*MW_w) );

    if (w<1.e-6){
		term3 = 0.0;
        d_term3 = 0.0;
        term4 = 0.0;
        d_term4 = 0.0;
   }
   else{
    term3 = R*T/(MW_a*MW_b)*grad_w*( 1./(w*(1.-w)) * MW_w );
    d_term3 = R*T/(MW_a*MW_b)*grad_w*(  (2.*w-1.)/(w*w*(1.-w)*(1.-w)) * MW_w
                                      + 1./(w*(1.-w)) * d_MW_w);

    term4 = R*T/MW_a*grad_w*( 1./(w*(1.-w)) * (w/act_w) * d_act_w );
    d_term4 = R*T/MW_a*grad_w*( (2.*w-1.)/(w*w*(1.-w)*(1.-w)) * (w/act_w) * d_act_w
                              + 1./(w*(1.-w)) *(act_w-w*d_act_w)/(act_w*act_w) * d_act_w
                              + 1./(w*(1.-w)) * (w/act_w) * d2_act_w );
   }

	d_ent = d_term1 + (term3+term4)*( d_term2*(term3+term4) + 2*term2*(d_term3+d_term4) );


	return d_ent;

}

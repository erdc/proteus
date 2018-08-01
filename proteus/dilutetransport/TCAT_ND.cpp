#include <math.h>
#include <iostream>

double grav = -980.;

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
};


double den(double x){
 double d,x2,x3;
 double a0 = 0.9971;
 double a1 = 0.83898;
 double a2 = 0.48132;
 double a3 = 0.86154;
 x2 = x*x; x3 = x2*x; 
 d = a0 + a1*x + a2*x2 + a3*x3;
 return 1.0;
}


double d_den(double x){
 double d;
 double a0 = 0.9971;
 double a1 = 0.83898;
 double a2 = 0.48132;
 double a3 = 0.86154;
 d = a1 + x*(2.0*a2 + 3.0*a3*x);
 return d;
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
 mu = conv * exp(eval);
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


double  ENTROPY_TCAT(double w, double grad_w, double grad_p, struct TCAT_var TCAT_v){

double x_a,MW_w,ent,D,velocity,MW_a,MW_b;
double poro,perm,mu,rho,diff,alpha_L,T,R;
poro = TCAT_v.poro;
perm = TCAT_v.perm;
diff = TCAT_v.diff;
alpha_L = TCAT_v.alpha_L;
MW_a = TCAT_v.MW_a;
MW_b = TCAT_v.MW_b;
T = TCAT_v.T;
R = TCAT_v.R;

mu = visc(0.0);
rho = den(0.0);
x_a = mole_frac(0.0,TCAT_v);
MW_w = mol_weight(x_a,TCAT_v);

velocity = -perm/mu*(grad_p - rho*grav)/poro; 
D = (diff + alpha_L*velocity);

ent = mu/(perm*T)*(poro*poro*velocity*velocity) + poro*rho*R*D/(w*(1.-w))*MW_w/MW_a/MW_b*grad_w*grad_w;
return ent;

}

double DENTROPY_TCAT(double w, double grad_w, double grad_p, struct TCAT_var TCAT_v) {

double x_a,d_x_a,MW_w,ent,d_ent,MW_a,MW_b,D,d_MW_w,velocity;
double poro,perm,mu,d_mu,rho,d_rho,diff,alpha_L,MW_a_2,MW_b_2,T,R,c;


poro = TCAT_v.poro;
perm = TCAT_v.perm;
diff = TCAT_v.diff;
alpha_L = TCAT_v.alpha_L;
T = TCAT_v.T;
R = TCAT_v.R;
MW_a = TCAT_v.MW_a;
MW_b = TCAT_v.MW_b;

mu = visc(0.0);
d_mu = d_visc(0.0);
rho = den(0.0);
velocity = -perm/mu*(grad_p - rho*grav)/poro; 


if (w<1.e-10){
	d_ent = 0.0;
}
else{
	d_rho = d_den(0.0);
	x_a = mole_frac(0.0,TCAT_v);
	d_x_a = d_mole_frac(0.0,TCAT_v);
	MW_w = mol_weight(0.0,TCAT_v);
	d_x_a = d_mole_frac(0.0,TCAT_v);
	d_MW_w = d_mol_weight(d_x_a,TCAT_v);
	D = (diff + alpha_L*velocity);
	c = poro*rho*D*R/MW_a/MW_b*grad_w*grad_w;
	d_ent =  c*(d_MW_w/(w*(1.-w)) + MW_w*(2.*w-1.)/(w*w*(1.-w)*(1.-w)) ) ;

}
return d_ent;

}


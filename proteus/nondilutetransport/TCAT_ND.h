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

double den(double x);
double d_den(double x);
double visc(double x);
double mole_frac(double x, struct TCAT_var TCAT_v);
double d_mole_frac(double x, struct TCAT_var TCAT_v);
double mol_weight(double x, struct TCAT_var TCAT_v);
double d_mol_weight(double d_x, struct TCAT_var TCAT_v);
double molal(double x);
double d_molal(double x);
double d2_molal(double x);
double activity(double x,double m);
double d_act(double x,double m, double dm, double act);
double d2_act(double x,double m, double dm, double act, double d_act);
double ENTROPY_TCAT(double w, double grad_w, double grad_p,  struct TCAT_var TCAT_v );
double DENTROPY_TCAT(double w, double grad_w, double grad_p, struct TCAT_var TCAT_v) ;
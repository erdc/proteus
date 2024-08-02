#include <math.h>
#include "ncursesw/curses.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define	Subprograms
#define	Int	extern int
#define	Double	extern double
#include	"Headers.h"

 

extern double SU;

void
	Title_block(FILE*), Input_Data_block(FILE*);

/**********************************************************************/
int Read_data(void)
/**********************************************************************/
{
Readtext(Title);
iff(Title, FINISH) return(0);
Read(MaxH,lf);
fscanf(Input1,"%s", Case); Skip;
iff(Case,Wavelength)
	{
	Read(L,lf);
	Height = MaxH/L;
	}
iff(Case,Period)
	{
	Read(T,lf);
	Height = MaxH/(T*T);
	}
Read(Current_criterion,d);
Read(Current,lf);
if(Current_criterion == 1) strcpy(Currentname, Current1);
if(Current_criterion == 2) strcpy(Currentname, Current2);

Read(n,d);
Read(nstep,d);

Input_Data_block(monitor);

if(strcmp(Theory,"Stokes")==0)
	{
	iff(Case,Wavelength)
		if(L > 10.)
			{
			printf("\nThe dimensionless wavelength is greater than 10.");
			printf("\nStokes theory should not be applied. Exiting.");
			getch();
			exit(1);
			}
	iff(Case,Period)
		if(T > 10.)
			{
			printf("\nThe dimensionless period is greater than 10.");
			printf("\nStokes theory should not be applied. Exiting.");
			getch();
			exit(1);
			}
	}

// Convergence criteria

Input2=fopen(Convergence_file,"r");
fgets(dummy,400,Input2);
fscanf(Input2,"%d", &number);fgets(dummy,400,Input2);
fscanf(Input2,"%le", &crit);fgets(dummy,400,Input2);
fclose(Input2);

// Number of data points to present results for

Input2 = fopen(Points_file,"r");
fgets(dummy,400,Input2);
// Number of points on surface profile (clustered quadratically near crest)
fscanf(Input2,"%d", &Surface_points);fgets(dummy,400,Input2);
// Number of vertical profiles
fscanf(Input2,"%d", &Nprofiles);fgets(dummy,400,Input2);
// Number of points in each profile
fscanf(Input2,"%d", &Points);fgets(dummy,400,Input2);

fclose(Input2);

return(1);
}

//	PRINT OUT TITLE BLOCKS

void Input_Data_block(FILE* file)
{
fprintf(file,"# %s", Title);
fprintf(file,"\n\n# Printing input data here to check");
fprintf(file,"\n\n# Height/Depth:%6.3f", MaxH);
iff(Case,Wavelength)
	{
	fprintf(file,"\n# Length/Depth:%7.2f", L);
	}
iff(Case,Period)
	{
	fprintf(file,"\n# Dimensionless Period T*sqrt(g/d):%7.2f", T);
	}
fprintf(file,"\n# Current criterion: %s,  Dimensionless value:%6.3lf", Currentname, Current);

if(strcmp(Theory,"Stokes")==0)
	{
	if(n<=5) sprintf(Method, "\n# Solution by %d-order Stokes theory", n);
	else
		{
		n = 5;
		sprintf(Method, "\n# Solution by %d-order Stokes theory", n);
		printf("\n\n# (A value of N > 5 has been specified for the Stokes theory.");
		printf("\n# I do not have a theory for that. The program has set N = 5)");
		}
	}
if(strcmp(Theory,"Fourier")==0)
	sprintf(Method, "\n# Solution by %d-term Fourier series", n);

fprintf(file,"\n%s\n", Method);
}

void Title_block(FILE* file)
{
// Highest wave - eqn (32) of Fenton (1990)
L = 2*pi/z[1];
Highest = (0.0077829*L*L*L+0.0095721*L*L+0.141063*L)
	/(0.0093407*L*L*L+0.0317567*L*L+0.078834*L+1);
fprintf(file,"# %s", Title);
fprintf(file,"\n%s\n", Method);
fprintf(file,"\n# Height/Depth:%6.3f, %3.0lf\%% of the maximum of H/d =%6.3f for this length:",
	z[2]/z[1],z[2]/z[1]/Highest*100., Highest);
fprintf(file,"\n# Length/Depth:%7.2f", 2*pi/z[1]);
fprintf(file,"\n# Dimensionless Period T*sqrt(g/d):%7.2f", z[3]/sqrt(z[1]));
fprintf(file,"\n# Current criterion: %s,  Dimensionless value:%6.3lf\n", Currentname, Current);
}

void Output(void)
{
int 		i, I;
double 	X, eta, y;
 double	Surface(double);
 void Point(double, double);

fprintf(monitor,"\n\n# Solution summary:\n\n");
Title_block(monitor);

// Print out summary file of solution

Title_block(Solution);

kd = z[1];
L=2*pi/z[1];
H=z[2]/z[1];
T=z[3]/sqrt(z[1]);
c=z[4]/sqrt(z[1]);
ce=z[5]/sqrt(z[1]);
cs=z[6]/sqrt(z[1]);
ubar=z[7]/sqrt(z[1]);
Q=ubar-z[8]/pow(kd,1.5);
R=1+z[9]/z[1];

pulse=z[8]+z[1]*z[5];
ke=0.5*(z[4]*pulse-z[5]*Q*pow(kd,1.5));

// Calculate potential energy, not by computing the mean of 1/2 (eta-d)^2
// but by exploiting orthogonality of the cosine functions to give the sum of 1/4 Y[i]^2
pe = 0;
for(i=1;i<=n;++i)
	pe += 0.25*pow(Y[i],2);

ub2=2.*z[9]-z[4]*z[4];
sxx=4.*ke-3.*pe+ub2*z[1]+2.*z[5]*(z[7]*z[1]-z[8]);
f=z[4]*(3.*ke-2.*pe)+0.5*ub2*(pulse+z[4]*z[1])+z[4]*z[5]*(z[7]*z[1]-z[8]);
q=z[7]*z[1]-z[8];
r=z[9]+z[1];
s=sxx-2.*z[4]*pulse+(z[4]*z[4]+0.5*z[1])*z[1];

fprintf(Solution, "\n# Stokes-Ursell number %7.3f", 0.5*z[2]/pow(z[1],3));
fprintf(Solution, "\n\n# Integral quantities - notation from Fenton (1988)");
fprintf(Solution, "\n# (1) Quantity, (2) symbol, solution non-dimensionalised by (3) g & wavenumber, and (4) g & mean depth\n");
fprintf(Solution, "\n# Water depth                        (d)" LO LO, z[1], 1.);
fprintf(Solution, "\n# Wave length                   (lambda)" LO LO, 2*pi, L);
fprintf(Solution, "\n# Wave height                        (H)" LO LO, z[2], H);
fprintf(Solution, "\n# Wave period                      (tau)" LO LO, z[3], T);
fprintf(Solution, "\n# Wave speed                         (c)" LO LO, z[4], c);
fprintf(Solution, "\n# Eulerian current                 (u1_)" LO LO, z[5], ce);
fprintf(Solution, "\n# Stokes current                   (u2_)" LO LO, z[6], cs);
fprintf(Solution, "\n# Mean fluid speed in frame of wave (U_)" LO LO, z[7], ubar);
fprintf(Solution, "\n# Volume flux due to waves           (q)" LO LO, z[8], z[8]/pow(kd,1.5));
fprintf(Solution, "\n# Bernoulli constant                 (r)" LO LO, z[9], z[9]/kd);
fprintf(Solution, "\n# Volume flux                        (Q)" LO LO, Q*pow(kd,1.5), Q);
fprintf(Solution, "\n# Bernoulli constant                 (R)" LO LO, R*kd, R);
fprintf(Solution, "\n# Momentum flux                      (S)" LO LO, s, s/kd/kd );
fprintf(Solution, "\n# Impulse                            (I)" LO LO, pulse, pulse/pow(kd,1.5));
fprintf(Solution, "\n# Kinetic energy                     (T)" LO LO, ke, ke/kd/kd);
fprintf(Solution, "\n# Potential energy                   (V)" LO LO, pe, pe/kd/kd);
fprintf(Solution, "\n# Mean square of bed velocity     (ub2_)" LO LO, ub2, ub2/kd);
fprintf(Solution, "\n# Radiation stress                 (Sxx)" LO LO, sxx, sxx/kd/kd);
fprintf(Solution, "\n# Wave power                         (F)" LO LO, f, f/pow(kd,2.5));

fprintf(Solution, "\n\n# Dimensionless coefficients in Fourier series" );
fprintf(Solution, "\n# Potential/Streamfn\tSurface elevations" );
fprintf(Solution, "\n# j, B[j], & E[j], j=1..N\n" );
for ( i=1 ; i <= n ; i++ )
fprintf(Solution, "\n%2d\t%15.7e\t%15.7e", i, B[i], Y[i]);
fprintf(Solution, "\n\n" );

// Surface - print out coordinates of points on surface for plotting plus check of pressure on surface

fprintf(Elevation,  "# %s\n", Title);
fprintf(Elevation,  "%s\n", Method);
fprintf(Elevation,  "\n# Surface of wave - trough-crest-trough,");
fprintf(Elevation,  " note quadratic point spacing clustered around crest");
fprintf(Elevation,  "\n# Non-dimensionalised with respect to depth");
fprintf(Elevation,  "\n# X/d, eta/d, & check of surface pressure\n");

for ( i=-Surface_points/2 ; i <= Surface_points/2; i++)
	{
	X = 2 * L * (i * fabs(i)/Surface_points/Surface_points);	//NB Quadratic point spacing, clustered near crest
   eta = Surface(X);
   Point(X,eta);
	fprintf(Elevation,  "\n%8.4lf\t%7.4f\t%7.0e", X, eta, Pressure);
	}
fprintf(Elevation,  "\n\n");

// Surface - print out Velocity and acceleration profiles plus check of Bernoulli

fprintf(Flowfield,  "# %s\n", Title);
fprintf(Flowfield,  "%s\n", Method);
fprintf(Flowfield,  "\n# Velocity and acceleration profiles and Bernoulli checks\n");
fprintf(Flowfield,  "\n# All quantities are dimensionless with respect to g and/or d\n");
fprintf(Flowfield,  "\n#*******************************************************************************");
fprintf(Flowfield,  "\n# y        u       v    dphi/dt   du/dt   dv/dt  du/dx   du/dy Bernoulli check  ");
fprintf(Flowfield,  "\n# -     -------------   -------  ------   -----  ------------- ---------------  ");
fprintf(Flowfield,  "\n# d        sqrt(gd)       gd        g       g       sqrt(g/d)        gd         ");
fprintf(Flowfield,  "\n#*******************************************************************************");

for(I = 0; I <= Nprofiles ; ++I)
	{
	X = 0.5 * L * I/(Nprofiles);
	eta = Surface(X);
	fprintf(Flowfield,  "\n\n# X/d = %8.4f, Phase = %6.1f°\n", X, X/L*360);

	for(i=0 ; i <= Points; ++i)
		{
		y = (i)*eta/(Points);
		Point(X, y);
		fprintf(Flowfield,  "\n%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f",
			y, u, v, dphidt, ut, vt, ux, uy, Bernoulli_check);
		}
	}
fprintf(Flowfield,  "\n\n");

/*
Procedure for recording every run - not activated in distribution versions
If the lines below are not commented out the program will add a line to a
file Catalogue.res, which could have this as a header:

# A continuing record of all runs with Fourier, Cnoidal, or Stokes.
# Any run of those programs adds a line to it.
# This can be edited at any time.
# Columns are: name of theory, N, H/d, L/d, Stokes-Ursell Number,
# wave height as a percentage of the highest possible for that L/d,
# mean horizontal velocity on a vertical line under the crest.

# Theory n   H/d       L/d     S-U Highest% u_crest_mean
*/

// To activate, de-comment these lines
/***************************************************************************
FILE *Output1;
double Velo[Points+1], sum1, sum2, ucm;
Output1 = fopen(Diagname,"a");

X = 0.;
eta = Surface(X);

for(i=0 ; i <= Points; ++i)
	{
	y = (i)*eta/(Points);
	Point(X, y);
	Velo[i] = u;
	}

for(i=1, sum1=0; i <= Points-1; i+=2) sum1 += Velo[i];
for(i=2, sum2=0; i <= Points-2; i+=2) sum2 += Velo[i];
ucm = (Velo[0]+4*sum1+2*sum2+Velo[Points])/3./Points;

I = strlen(Theory)+1;
for(i=I; i <=8 ; ++i) strcat(Theory," ");
fprintf(Output1,"\n%s%2d\t%7.4f\t%8.3f\t%7.3f\t%3.0f\t%7.4f",
			Theory, n, H, L, 0.5*z[2]/pow(z[1],3), z[2]/z[1]/Highest*100., ucm);
*************************************************************************/
}

// Surface elevation

double Surface(double x)
{
int j;
static double kEta;

kEta = kd;
for ( j = 1 ; j < n ; j++ )
	kEta += Y[j] * cos(j*x*kd);
kEta += 0.5*Y[n] * cos(n*x*kd);
return (kEta/kd);
}

// Velocities, accelerations, and pressure at a point

void Point(double X, double y)
{
int j;
double Cosh, Sinh, Sin, Cos;
double coshdelta,sinhdelta;

u = v = ux = vx = phi = psi = 0.;

for ( j = 1 ; j <= n ; j++ )
	{
   Cos  = cos(j*X*kd);
	Sin  = sin(j*X*kd);
   coshdelta = cosh(j*kd*(y-1.));
   sinhdelta = sinh(j*kd*(y-1.));
   Cosh = coshdelta+sinhdelta*Tanh[j];
	Sinh = sinhdelta+coshdelta*Tanh[j];
   phi += B[j] * Cosh * Sin;
	u += j * B[j] * Cosh * Cos;
	v += j * B[j] * Sinh * Sin;
	ux += - j * j * B[j] * Cosh * Sin;
	vx += j * j * B[j] * Sinh * Cos;
	}

// All PHI, PSI, u, v, ux and vx are dimensionless w.r.t. g & k.
// Now convert to dimensionless w.r.t. d.

phi /= pow(kd,1.5);
u /= pow(kd,0.5);
v /= pow(kd,0.5);
ux *= pow(kd,0.5);
vx *= pow(kd,0.5);

u = ce + u;
phi = ce * X + phi;
dphidt = -c * u;

ut = -c * ux;
vt = -c * vx;
uy = vx;
vy = -ux;
dudt = ut + u*ux + v*uy;
dvdt = vt + u*vx + v*vy;
Pressure = R - y - 0.5 * ((u-c)*(u-c)+v*v);
Bernoulli_check = dphidt + Pressure + y + 0.5*(u*u+v*v)-(R-0.5*c*c);

return;
}


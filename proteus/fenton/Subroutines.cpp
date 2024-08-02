
//	SUBROUTINES FOR FOURIER APPROXIMATION METHOD

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "ncursesw/curses.h"

#define	iff(x,y)	if(strcmp(x,#y)==0)
#define	pi			3.14159265358979324

#define	Int	extern int
#define	Double	extern double

double *dvector(long, long);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
void 	free_dvector(double *, long , long );
void 	free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);

extern char		Case[];

#include "Headers.h"

// **************************************************
// CALCULATE INITIAL SOLUTION FROM LINEAR WAVE THEORY
// **************************************************

void init()
{
int i;
double a, b, t;

iff(Case,Period)
	{
	a=4.*pi*pi*height/Hoverd;
	b=a/sqrt(tanh(a));
	t=tanh(b);
	z[1]=b+(a-b*t)/(t+b*(1.-t*t));
	}
else
	z[1]=2.*pi*height/Hoverd;

z[2]=z[1]*Hoverd;
z[4]=sqrt(tanh(z[1]));
z[3]=2.*pi/z[4];
if(Current_criterion==1)
	{
	z[5]=Current*sqrt(z[2]);
	z[6]=0.;
	}
else
	{
	z[6]=Current*sqrt(z[2]);
	z[5]=0.;
	}
z[7]=z[4];
z[8]=0.;
z[9]=0.5*z[7]*z[7];
cosa[0]=1.;
sina[0]=0.;
z[10]=0.5*z[2];
	for( i=1 ; i<=n ; i++ )
		{
		cosa[i]=cos(i*pi/n);
		cosa[i+n]=cos((i+n)*pi/n);
		sina[i]=sin(i*pi/n);
		sina[i+n]=sin((i+n)*pi/n);
		z[n+i+10]=0.;
		z[i+10]=0.5*z[2]*cosa[i];
		}
z[n+11]=0.5*z[2]/z[7];

for( i=1 ; i<=9 ; i++ )
	sol[i][1] = z[i];
for( i=10 ; i<=num ; i++ )
	sol[i][1] = 0.;

return;
}

//	EVALUATION OF EQUATIONS.

double Eqns(double *rhs)
{
int i, j, m, it, nm;
double c, e, s, u, v;

rhs[1]=z[2]-z[1]*Hoverd;

iff(Case,Wavelength)
	rhs[2]=z[2]-2.*pi*height;
else
	rhs[2]=z[2]-height*z[3]*z[3];

rhs[3]=z[4]*z[3]-pi-pi;
rhs[4]=z[5]+z[7]-z[4];
rhs[5]=z[6]+z[7]-z[4];

rhs[5]=rhs[5]-z[8]/z[1];
for (i=1; i<=n; i++ )
	{
	coeff[i]=z[n+i+10];
	Tanh[i] = tanh(i*z[1]);
	}
it=6;
if(Current_criterion==1)it=5;
rhs[6]=z[it]-Current*sqrt(z[1]); // Correction made 20.5.2013, z[2] changed to z[1]
rhs[7]=z[10]+z[n+10];
for (i=1 ; i<= n-1 ; i++ )
	rhs[7]=rhs[7]+z[10+i]+z[10+i];
rhs[8]=z[10]-z[n+10]-z[2];
for ( m=0 ; m <= n ; m++ )
	{
	psi=0.;
	u=0.;
	v=0.;
	for (j=1 ; j <= n ; j++ )
		{
		nm = (m*j) % (n+n);
		e=exp(j*(z[10+m]));
		s=0.5*(e-1./e);
		c=0.5*(e+1./e);
		psi=psi+coeff[j]*(s+c*Tanh[j])*cosa[nm];
		u=u+j*coeff[j]*(c+s*Tanh[j])*cosa[nm];
		v=v+j*coeff[j]*(s+c*Tanh[j])*sina[nm];
		}
	rhs[m+9]=psi-z[8]-z[7]*z[m+10];
	rhs[n+m+10]=0.5*(pow((-z[7]+u),2.)+v*v)+z[m+10]-z[9];
	}

for (j=1, s=0. ; j <= num ; j++ ) s += rhs[j]*rhs[j];
return s;
}

// **************************************************
//	SET UP JACOBIAN MATRIX AND SOLVE MATRIX EQUATION
// **************************************************

double Newton(int count)
{
double	**a, *rhs, *x;
double 	h, sum;
void		Solve(double **, double *, int, int, double *, int, int);

int i, j;

Eqns(rhs1);

if(count>=1)
{
++count;
rhs=dvector(1,num);
x=dvector(1,num);
a=dmatrix(1,num,1,num);
}

for ( i=1 ; i<=num ; i++ )
	{
	h=0.01*z[i];
	if(fabs(z[i]) < 1.e-4) h = 1.e-5;
	z[i]=z[i]+h;
	Eqns(rhs2);
	z[i]=z[i]-h;
	rhs[i] = -rhs1[i];
	for ( j=1 ; j<=num ; j++ )
		a[j][i] = (rhs2[j] - rhs1[j])/h;
	}

// **************************************************
// SOLVE MATRIX EQUATION
// **************************************************

Solve(a, rhs, (int)num, (int)num, x, (int)num, (int)num);

for ( i=1 ; i<=num ; i++ )
	z[i] += x[i];

for ( sum = 0., i=10 ; i<= n+10 ; i++ )
			sum += fabs(x[i]);
sum /= n;

free_dvector(rhs,1,num);
free_dvector(x,1,num);
free_dmatrix(a,1,num,1,num);

return(sum);
}


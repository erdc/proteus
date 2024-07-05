
// Steady wave program - C++ version

#include <math.h>
#include <stdio.h>
#include <sys/types.h> 
#include <string.h>
#include "ncursesw/curses.h"
#include <stdlib.h>
#define	ANSI
#include "Allocation.h"

void
	init(void), Solve(void), Title_block(FILE*),	Output(void);
int
	flushall(void);

#define Main
#define	Int		int
#define	Double	double
#include "Headers.h"
//#define	Diagnostic
//#ifdef Diagnostic
//char Diagname[30], Theory[10];
//#endif
double SU;

void runfourier()
{
int	i, j, iter, m;
int	Read_data(void);

double	Newton(int), dhe, dho, error, **CC;
void 	Powell(double *, double **, int, double, int *, double *,double (*)(double *));

Input1 = fopen("Data.dat","r");

strcpy(Convergence_file,"Convergence.dat");
strcpy(Points_file,"Points.dat");
monitor = stdout;
strcpy(Theory,"Fourier");
strcpy(Diagname,"Catalogue.res");

for ( wave=1 ; wave<2; wave++ )
	{
	if (Read_data() == 0) break;
	num=2*n+10;
	dhe=Height/nstep;
	dho=MaxH/nstep;

	CC = dmatrix(1,num,1,num);
	for ( j=1; j <=num ; ++j)
		{
		for ( i=1; i <=num ; ++i)
			CC[j][i] = 0.;
		CC[j][j] = 1.;
		}
	Y = dvector(0,num);
	z = dvector(1,num);
	rhs1 = dvector(1,num);
	rhs2 = dvector(1,num);
	coeff = dvector(0, n);
	cosa = dvector(0,2*n);
	sina = dvector(0,2*n);
	sol = dmatrix(0,num,1,2);
	B = dvector(1, n);
	Tanh = dvector(1,n);

//	Commence stepping through steps in wave height

	for ( ns = 1 ; ns <= nstep ; ns++ )
		{
		height=ns*dhe;
		Hoverd=ns*dho;
		fprintf(monitor,"\n\nHeight step %2d of %2d\n", ns, nstep);

//	Calculate initial linear solution

	if(ns <= 1) init();

//	Or, extrapolate for next wave height, if necessary

	else
		for ( i=1 ; i <= num ; i++ )
			z[i]=2.*sol[i][2]-sol[i][1];

//	Commence iterative solution

	for (iter=1 ; iter <= number ; iter++ )
		{
		fprintf(monitor,"\nIteration%3d:", iter);

//	Calculate right sides of equations and differentiate numerically
//	to obtain Jacobian matrix, then solve matrix equation

		error = Newton(iter);

//	Convergence criterion satisfied?

		fprintf(stdout," Mean of corrections to free surface: %8.1e", error);
		if(ns == nstep)	criter = 1.e-10 ;
		else			criter = crit;
		if((error < criter * fabs(z[1]))  && iter > 1 ) break;
		if(iter == number)
			{
			fprintf(stdout,"\nNote that the program still had not converged to the degree specified\n");
			}

//	Operations for extrapolations if more than one height step used

		if(ns == 1)
			for ( i=1 ; i<=num ; i++ )
				sol[i][2] = z[i];
		else
			for ( i=1 ; i<=num ; i++ )
				{
				sol[i][1] = sol[i][2];
				sol[i][2] = z[i];
				}
			}

//	Fourier coefficients (for surface elevation by slow Fourier transform)

	for ( Y[0] = 0., j = 1 ; j <= n ; j++ )
		{
	   B[j]=z[j+n+10];
      sum = 0.5*(z[10]+z[n+10]*pow(-1.,(double)j));
		for ( m = 1 ; m <= n-1 ; m++ )
			sum += z[10+m]*cosa[(m*j)%(n+n)];
		Y[j] = 2. * sum / n;
		}
	} // End stepping through wave heights

// Print  results

	Solution=fopen("Solution.res","w");
	Elevation = fopen("Surface.res","w");
	Flowfield = fopen("Flowfield.res","w");
	Output();
    fflush(NULL);
        printf("\nTouch key to continue\n\n"); getch();

            free_dmatrix(CC,1,num,1,num);
                free_dvector(Y,0,num);
                    free_dvector(z, 1,num);
                        free_dvector(rhs1, 1,num);
                            free_dvector(rhs2, 1,num);
                                free_dvector(coeff, 0, n);
                                    free_dvector(cosa, 0,2*n);
                                        free_dvector(sina, 0,2*n);
                                            free_dmatrix(sol, 0,num,1,2);
                                                free_dvector( B, 1, n);
                                                    free_dvector(Tanh, 1,n);
} // End stepping through waves

printf("\nFinished\n");
}

int main(void)
{
    runfourier();
} // End main program

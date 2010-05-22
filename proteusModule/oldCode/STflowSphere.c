/*
 Stokes Flow around Sphere of radius r^s and center x^s, y^s, z^s.
 Velocity is u,v,w and viscosity is mu.
 Array rwork contains all the above variables.
 Input	: t, x
		: x[npoints*3] = {x,y,z} for each point
 Output : either v^x, v^y, v^z or P ==> 4 solution functions
*/
#include "analyticalSolutions.h"

int STflowSphere(int *iwork, double *rwork, int nPoints, double t, double *x, double *u)
{
//	double rwork[] = {u,v,w,rS,xS,yS,zS,mu};
//	double *u = {vX,vY,vZ};

	int i, j;
	const int dim=3;
	double in_xBarE1, in_u, U, vR, vTHETA, theta, x_u, r, in_xBar;
	double tempvar[3], e1[3], e2[3], xBar[3] ;
	double *P = malloc(nPoints*sizeof(double));

	if (P==NULL)
	{
		printf("Out of memory at STflowSphere for P\n");
		exit(1);
	}
// values not depending on nPoints.
	for (in_u=0.0,j=0; j<dim; j++)
	{
		in_u += pow(rwork[j],2);
	}
	U = sqrt(in_u);

	for (j=0; j<dim; j++)
	{
		e1[j] = rwork[j]/U;
	}

 	for (i=0; i<nPoints; i++)
 	{
		for (in_xBar=0.0,x_u=0.0,j=0; j<dim; j++)
		{
			xBar[j] = x[i*3+j] - rwork[j+4];
			x_u += xBar[j] * rwork[j];
			in_xBar += pow(xBar[j],2);
		}
		r = sqrt(in_xBar);
		theta = acos( x_u ) / (r * U);

		for (j=0; j<dim; i++)
		{
			tempvar[j] = xBar[j] - xBar[j]*e1[j];
		}
		in_xBarE1 = sqrt( pow(tempvar[0],2) + pow(tempvar[1],2) + pow(tempvar[2],2) );

		e2[0] = (tempvar[0])/in_xBarE1;
		e2[1] = (tempvar[1])/in_xBarE1;
		e2[2] = (tempvar[2])/in_xBarE1;

		vR = U * cos(theta) * (1 - (3*rwork[3]/(2*r)) + pow(rwork[3],3)/(2*pow(r,3)) ) ;
		vTHETA = -U * sin(theta) * (1 - (3*rwork[3])/(4*r) - pow(rwork[3],3)/(4*pow(r,3)) );

// P is currently on returned.
		P[i] = (-3 * rwork[3]*rwork[7]*U*cos(theta))/(2*pow(r,2));
		
		for (j=0; j<dim; j++)
		{
			u[i*3+j] = vR * cos(vTHETA) * e1[j] + vR * sin(vTHETA) * e2[j] ;
		}
 	}
 
	return 0;
}

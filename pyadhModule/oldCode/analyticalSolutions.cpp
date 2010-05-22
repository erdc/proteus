//
// analytical functions
//
#include "analyticalSolutions.h"

int diffusionSin1D(int om_x, int nPoints, double t, double *x, double *u)
{
	
	for (int i=0; i<nPoints; i++)
	{
		u[i] = sin(2.0*PI*om_x*x[i*2+0]);
	}
	
	return 0;
}
int diffusionSin2D(int om_x, int om_y, int nPoints, double t, double *x, double *u)
{

	cout<<"pi = "<<PI<<" om_x = "<<om_x<<" om_y = "<<om_y<<endl;
	for (int i=0; i<nPoints; i++)
	{
		u[i] = sin(2.0*PI*om_x*x[i*3+0]) * sin(2.0*PI*om_y*x[i*3+1]);
		cout<<"i*3 = "<<i*3<<" i*3+1 = "<<i*3+1<<endl;
		cout<<"\t\tu["<<i<<"] = "<<u[i]<<endl;
	}

	return 0;
}
int diffusionSin3D(int om_x, int om_y, int om_z, int nPoints, double t, double *x, double *u)
{
	
	for (int i=0; i<nPoints; i++)
	{
		u[i] = sin(2.0*PI*om_x*x[i*4+0]) * sin(2.0*PI*om_y*x[i*4+1]) * sin(2.0*PI*om_z*x[i*4+2]);
	}
	
	return 0;
}
int diffusionSin1D_r(int om_x, int nPoints, double t, double *x, double *r)
{

	for (int i=0; i<nPoints; i++)
	{
		r[i] = pow(2.0*PI*om_x*sin(2.0*PI*om_x*x[i*2+0]),2.0);
	}
	
	return 0;
}
int diffusionSin2D_r(int om_x, int om_y, int nPoints, double t, double *x, double *r)
{

	for (int i=0; i<nPoints; i++)
	{
		r[i] = pow(2.0*PI*om_x*sin(2.0*PI*om_x*x[i*3+0]),2.0) + pow(2.0*PI*om_y*sin(2.0*PI*om_y*x[i*3+1]),2.0);
	}
	
	return 0;
}
int diffusionSin3D_r(int om_x, int om_y, int om_z, int nPoints, double t, double *x, double *r)
{

	for (int i=0; i<nPoints; i++)
	{
		r[i] = pow(2.0*PI*om_x*sin(2.0*PI*om_x*x[i*4+0]),2.0) + pow(2.0*PI*om_y*sin(2.0*PI*om_y*x[i*4+1]),2.0) +
			   pow(2.0*PI*om_z*sin(2.0*PI*om_z*x[i*4+2]),2.0);
	}
	
	return 0;
}


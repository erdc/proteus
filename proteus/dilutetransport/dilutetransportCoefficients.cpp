#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <memory.h>
#include <stdlib.h>
#include "TCAT_ND.h"
#include "dilutetransportCoefficients.h"

void ConMassFluidEvaluate(const int nPoints,
                          const int nSpace,
                          const double K,
                          const double grav,
                          const double L,
                          const double *u,
                          double *m,
                          double *dm,
                          double *f,
                          double *df,
                          double *phi,
                          double *dphi,
                          double *a,
                          double *r,
                          double *x,
                          double *mass_frac)
{
  int k,I;
  const int nSpace2=nSpace*nSpace;
  double density,d_density,mu;


  for (k=0;k<nPoints;k++){
      m[k] = 0.0;
      dm[k] = 0.0;
      density = den(0.0);
      mu = visc(0.0);

      f[k] = 0.0;
      df[k] = 0.0;
      phi[k] = u[k];
      dphi[k] = 1.0;
      for (I=0;I<nSpace;I++){
          phi[k] = phi[k] - density*grav*(x[k*3+I] - L);
          a[k*nSpace2+I*nSpace+I] = K/mu*density;
	  }
   }
}


void DiluteDispersionEvaluate(const int nPoints,
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
  double density,d_density,mu;

  for (k=0;k<nPoints;k++){


      density = den(0.0);
      mu = visc(0.0);

      m[k] = u[k]*poro;
      dm[k] = poro;

      for (I=0;I<nSpace;I++){
          a[k*nSpace2+I*nSpace+I] = poro*(diff+alpha_L*velocity[k*nSpace+I]);
	  }
   }
}


void DiluteEvaluate(const int nPoints,
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
  double density,d_density,mu;

  for (k=0;k<nPoints;k++){


      density = den(0.0);
      mu = visc(0.0);

      m[k] = u[k]*poro;
      dm[k] = poro;
      f[k] = poro*velocity[k]*u[k];
      df[k] = poro*velocity[k];

      for (I=0;I<nSpace;I++){
          a[k*nSpace2+I*nSpace+I] = poro*(diff+alpha_L*velocity[k*nSpace+I]);
	  }
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
  double density,d_density,mu;

  for (k=0;k<nPoints;k++){


      density = den(0.0);

      m[k] = u[k]*density*poro;
      dm[k] = poro*density;
      f[k] = poro*density*velocity[k]*u[k];
      df[k] = poro*density*velocity[k];

      for (I=0;I<nSpace;I++){
          a[k*nSpace2+I*nSpace+I] = 0.0;
	  }
   }
}

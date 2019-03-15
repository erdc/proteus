//
// analytical functions
//
#include "analyticalSolutions.h"

/**
   \ingroup analyticalSolutions
   @{
*/
/**
  \brief Couette Flow between two parallel plates. One moving relative to the other with constant seperation (width).
  @param iwork  NOT USED
  @param rwork[0] velocity of moving plate in the x-dir, \f$ vx \f$
  @param rwork[1] width between plates, \f$ h \f$
  @param rwork[2] x-axis offset, \f$ xs \f$
  @param rwork[3] y-axis offset, \f$ ys \f$
  @param rwork[4] z-axis offset, \f$ zs \f$
  @param nPoints total number of points
  @param t NOT USED
  @param x input array 
  @param u output array

   x-momentum equation
  \f[ \nabla^2 u_x = 0 \f]
  \f[ u(y) = \frac{vx}{h} * y\f]
  \f[ u(0) = 0 \f]
  \f[ u(h) = vx  \f]

  
  Axis of origin: at bottom plate\n
  

  \return status code
*/

int PlaneCouetteFlow_u(int *iwork, double *rwork, int nPoints, double t, double *x, double *u)
{
/*
  Couette Flow between two parallel plates. One moving relative to the other with constant seperation (width).
  Flow is axial. u != 0, v = w = 0 [Plates are very wide(z) and very long(x)]
  vx = Velocity of moving plate in the x-dir.
  h  = width between plates in the y-dir.
  Axis of orgin: bottom plate
  xs,ys,zs = axis offset

  \del^2(\vecV) = 0
  
  d^2u/dy^2 = 0
  u(y) =vx/h * Y

  u(0) = 0
  u(h) = vx
*/
	int i;
	double vx=rwork[0], h=rwork[1];
        double ys=rwork[3];
        double Y;

 	for (i=0; i<nPoints; i++)
 	{
             Y = x[i*3+1] - ys;
		
	     u[i] = vx*(Y/h);
 	}
 
	return 0;
}

/**
  \brief Sinusoidal 1D diffusion
  @param iwork[0] angular frequency, \f$ \omega_x \f$
  @param rwork NOT USED
  @param nPoints - total number of points
  @param t NOT USED
  @param x input array
  @param u output array

  \f[ -\Delta u = f \f] 
  \f[ f = -(2 \pi \omega_x)^2 \sin(2 \pi \omega_x x) \f]
  
  \return status code
*/

int diffusionSin1D(int *iwork, double *rwork, int nPoints, double t, double *x, double *u)
{
  int i;
  double omega_x = iwork[0];

  for (i=0; i<nPoints; i++)
    {
      u[i] = -pow(2.0*PI*omega_x,2.0)*sin(2.0*PI*omega_x*x[i*3+0]);
    }
	
    return 0;
}
/**
  \brief Sinusoidal 2D diffusion
  @param iwork[0] x-component angular frequency, \f$ \omega_x \f$
  @param iwork[1] y-component angular frequency, \f$ \omega_y \f$
  @param rwork NOT USED
  @param nPoints - total number of points
  @param t NOT USED
  @param x input array
  @param u output array

  \f[ -\Delta u = f \f] 
  \f[ f = -(2\pi\omega_x)^2\sin(2 \pi \omega_x x) -(2 \pi \omega_y)^2\sin(2 \pi \omega_y y)\f]
  
  \return status code
*/
int diffusionSin2D(int *iwork, double *rwork, int nPoints, double t, double *x, double *u)
{
	int i;
        double omega_x=iwork[0];
        double omega_y=iwork[1];
 
	for (i=0; i<nPoints; i++)
	{
	     u[i] = -pow(2.0*PI*omega_x,2.0)*sin(2.0*PI*omega_x*x[i*3+0]) 
                    -pow(2.0*PI*omega_y,2.0)*sin(2.0*PI*omega_y*x[i*3+1]) ;
	}

	return 0;
}
/**
  \brief Sinusoidal 3D diffusion
  @param iwork[0] x-component angular frequency, \f$ \omega_x \f$
  @param iwork[1] y-component angular frequency, \f$ \omega_y \f$
  @param iwork[2] z-component angular frequency, \f$ \omega_z \f$
  @param rwork NOT USED
  @param nPoints - total number of points
  @param t NOT USED
  @param x input array
  @param u output array

  \f[ -\Delta u = f \f] 
  \f[ f = -(2\pi\omega_x)^2\sin(2\pi\omega_x x) -(2\pi\omega_y)^2\sin(2\pi\omega_y y)-(2\pi\omega_z)^2sin(2 \pi \omega_z z)\f]
  
  \return status code
*/
int diffusionSin3D(int *iwork, double *rwork, int nPoints, double t, double *x, double *u)
{
	int i;
        double omega_x=2.0*PI*iwork[0];
        double omega_y=2.0*PI*iwork[1];
        double omega_z=2.0*PI*iwork[2];
 
	for (i=0; i<nPoints; i++)
	{
	     u[i] = -pow(2.0*PI*omega_x,2.0)*sin(2.0*PI*omega_x*x[i*3+0]) 
                    -pow(2.0*PI*omega_y,2.0)*sin(2.0*PI*omega_y*x[i*3+1])
                    -pow(2.0*PI*omega_z,2.0)*sin(2.0*PI*omega_z*x[i*3+2]) ;
	}

	return 0;
}
/**
  \brief Sinusoidal 1D diffusion (reaction)
  @param iwork[0] x-component angular frequency, \f$ \omega_x \f$
  @param rwork NOT USED
  @param nPoints - total number of points
  @param t NOT USED
  @param x input array
  @param u NOT USED
  @param r output array

  \f[ -\Delta u = -r \f] 
  \f[ r = (2 \pi \omega_x)^2 \sin(2 \pi \omega_x x) \f]  
  
  \return status code
*/
int diffusionSin1D_r(int *iwork, double *rwork, int nPoints, double t, double *x, double *u, double *r)
{
  	int i;
        double omega_x=iwork[0];

	for (i=0; i<nPoints; i++)
	{
	     r[i] = pow(2.0*PI*omega_x,2.0)*sin(2.0*PI*omega_x * x[i*3+0]);
	}
	
	return 0;
}
/**
  \brief Sinusoidal 2D diffusion (reaction)
  @param iwork[0] x-component angular frequency, \f$ omega_x \f$
  @param iwork[1] y-component angular frequency, \f$ omega_y \f$
  @param rwork NOT USED
  @param nPoints - total number of points
  @param t NOT USED
  @param x input array
  @param u NOT USED
  @param r output array

  \f[-\Delta u = -r \f] 
  \f[r= (2\pi\omega_x)^2\sin(2 \pi \omega_x x)+(2 \pi \omega_y)^2\sin(2 \pi \omega_y y)\f]
    
  \return status code
*/
int diffusionSin2D_r(int *iwork, double *rwork, int nPoints, double t, double *x, double *u, double *r)
{
	int i;
        double omega_x=iwork[0];
        double omega_y=iwork[1];

	for (i=0; i<nPoints; i++)
	{
          r[i] = pow(2.0*PI*omega_x,2.0) * sin(2.0*PI*omega_x * x[i*3+0]) +
                 pow(2.0*PI*omega_y,2.0) * sin(2.0*PI*omega_y * x[i*3+1]) ;
	}
	
	return 0;
}
/**
  \brief Sinusoidal 3D diffusion (reaction)
  @param iwork[0] x-component angular frequency, \f$ omega_x \f$
  @param iwork[1] y-component angular frequency, \f$ omega_y \f$
  @param iwork[2] z-component angular frequency, \f$ omega_z \f$
  @param rwork NOT USED
  @param nPoints - total number of points
  @param t NOT USED
  @param x input array
  @param u NOT USED
  @param r output array

  \f[-\Delta u = -r \f]
  \f[r=(2\pi\omega_x)^2\sin(2\pi\omega_x x) +(2\pi\omega_y)^2\sin(2\pi\omega_y y)+(2\pi\omega_z)^2sin(2 \pi \omega_z z)\f]
  
  \return status code
*/
int diffusionSin3D_r(int *iwork, double *rwork, int nPoints, double t, double *x, double *u, double *r)
{
	int i;
        double omega_x=iwork[0];
        double omega_y=iwork[1];
        double omega_z=iwork[2];
 
	for (i=0; i<nPoints; i++)
	{
          r[i] = pow(2.0*PI*omega_x,2.0) * sin(2.0*PI*omega_x * x[i*3+0]) +
                 pow(2.0*PI*omega_y,2.0) * sin(2.0*PI*omega_y * x[i*3+1]) +
		 pow(2.0*PI*omega_z,2.0) * sin(2.0*PI*omega_z * x[i*3+2]) ;
	}
	
	return 0;
}
/**
  \brief Linear Advective Diffusion Dirac Initial Condition
  @param iwork NOT USED
  @param rwork[0:2] \f$ b[0:2] \f$
  @param rwork[3] \f$ n \f$
  @param rwork[4] \f$ a \f$
  @param rwork[5] \f$ tStart \f$
  @param rwork[6] \f$ u_0 \f$
  @param rwork[7:9] \f$ x_0[0:2] \f$
  @param nPoints total number of points
  @param t NOT USED
  @param x input array
  @param u output array

  The exact solution of
    
   \f[ u_t + \nabla \cdot (b u - a \nabla u) = 0\f]
    
    on an infinite domain with dirac initial data
    
   \f[u_0 = \int u_0 \delta(x - x_0)\f]

    also returns advective, diffusive, and total flux

  \return status code
*/
int LinearAD_DiracIC(int *iwork, double *rwork, int nPoints, double T, double *x, double *u)
/*
    The exact solution of
    
    u_t + \deld(b u - a \grad u) = 0
    
    on an infinite domain with dirac initial data
    
    u0 = \int u0 \delta(x - x0)

    also returns advective, diffusive, and total flux
*/
{
	int i, j;
	double b[3]={rwork[0],rwork[1],rwork[2]};
	double n=rwork[3], a=rwork[4], tStart=rwork[5], u0=rwork[6];
	double x0[3]={rwork[7],rwork[8],rwork[9]};
	double y[3], exp_arg, yDoty, t = T + tStart;

    for (i=0; i<nPoints; i++)
	{
		for (yDoty=0.0,j=0; j<3; j++)
		{
			y[j] = x[i*3+j] - x0[j] - b[j]*t;
			yDoty += y[j]*y[j];
		}
		exp_arg = yDoty / (4.0 * a * t) ;
		
		if (exp_arg > 100)
		{
			u[i] = 0.0;
		}
		else
		{
			u[i] = u0*exp(-exp_arg) / pow( (4.0*a*PI*t),(n/2.0) );
		}
	}

    return 0;
}
/**
  \brief Linear Advective Diffusion Dirac Initial Condition(advective velocity)
  @param iwork NOT USED
  @param rwork[0:2] \f$ b[0:2] \f$
  @param rwork[3] \f$ n \f$
  @param rwork[4] \f$ a \f$
  @param rwork[5] \f$ tStart \f$
  @param rwork[6] \f$ u_0 \f$
  @param rwork[7:9] \f$ x_0[0:2] \f$
  @param nPoints total number of points
  @param t NOT USED
  @param x input array
  @param f output array, vector

  The exact solution of
    
   \f[ u_t + \nabla \cdot (b u - a \nabla u) = 0\f]
    
    on an infinite domain with dirac initial data
    
   \f[u_0 = \int u_0 \delta(x - x_0)\f]

  \return status code
*/
int LinearAD_DiracIC_advectiveVelocity(int *iwork, double *rwork, int nPoints, double T, double *x, double *f)
{
  double b[3]={rwork[0],rwork[1],rwork[2]};
  double *u = (double *)malloc(nPoints*sizeof(double));

	int i, j, iret;
	
	iret = LinearAD_DiracIC(iwork, rwork,  nPoints, T, x, u);

	for (i=0; i<nPoints; i++)
	{
	    for (j=0; j<3; j++)
	    {
	         f[i*3+j] = b[j]*u[i];
	    }
	}
	free(u);

	return iret;
}
/**
  \brief Linear Advective Diffusion Dirac Initial Condition(diffusive velocity)
  @param iwork NOT USED
  @param rwork[0:2] \f$ b[0:2] \f$
  @param rwork[3] \f$ n \f$
  @param rwork[4] \f$ a \f$
  @param rwork[5] \f$ tStart \f$
  @param rwork[6] \f$ u_0 \f$
  @param rwork[7:9] \f$ x_0[0:2] \f$
  @param nPoints total number of points
  @param t NOT USED
  @param x input array
  @param f output array, vector

  The exact solution of
    
   \f[ u_t + \nabla \cdot (b u - a \nabla u) = 0\f]
    
    on an infinite domain with dirac initial data
    
   \f[u_0 = \int u_0 \delta(x - x_0)\f]

  \return status code
*/
int LinearAD_DiracIC_diffusiveVelocity(int *iwork, double *rwork, int nPoints, double T, double *x, double *f)
{
	int i, j, iret;
	double a=rwork[4];
	
	iret = LinearAD_DiracIC_du(iwork, rwork, nPoints, T, x, f);

	for (i=0; i<nPoints; i++)
	{
	    for (j=0; j<3; j++)
	    {
		 f[i*3+j] = -a * f[i*3+j];
	    }
	}
	return iret;
}
/**
  \brief Linear Advective Diffusion Dirac Initial Condition (du)
  @param iwork NOT USED
  @param rwork[0:2] \f$ b[0:2] \f$
  @param rwork[3] \f$ n \f$
  @param rwork[4] \f$ a \f$
  @param rwork[5] \f$ tStart \f$
  @param rwork[6] \f$ u_0 \f$
  @param rwork[7:9] \f$ x_0[0:2] \f$
  @param nPoints total number of points
  @param t NOT USED
  @param x input array
  @param f output array, vector

  The exact solution of
    
   \f[ u_t + \nabla \cdot (b u - a \nabla u) = 0\f]
    
    on an infinite domain with dirac initial data
    
   \f[u_0 = \int u_0 \delta(x - x_0)\f]

    also returns advective, diffusive, and total flux

  \return status code
*/
int LinearAD_DiracIC_du(int *iwork, double *rwork, int nPoints, double T, double *x, double *f)
{
	int i, j, iret;
	double b[3]={rwork[0],rwork[1],rwork[2]};
	double a=rwork[4], tStart=rwork[5];
	double x0[3]={rwork[7],rwork[8],rwork[9]};
	double y[3], t=T + tStart;
	double *u = (double *)malloc(nPoints*sizeof(double));

	iret = LinearAD_DiracIC(iwork, rwork, nPoints, T, x, u);

	for (i=0; i<nPoints; i++)
	{
		for (j=0; j<3; j++)
		{
			y[j] = x[i*3+j] - x0[j] - b[j]*t;
			f[i*3+j] = u[i]*2.0*y[j] / (4.0*a*t);
		}
	}
	free(u);

	return iret;
}

/**
  \brief Linear Advective Diffusion Dirac Initial Condition(total velocity)
  @param iwork NOT USED
  @param rwork[0:2] \f$ b[0:2] \f$
  @param rwork[3] \f$ n \f$
  @param rwork[4] \f$ a \f$
  @param rwork[5] \f$ tStart \f$
  @param rwork[6] \f$ u_0 \f$
  @param rwork[7:9] \f$ x_0[0:2] \f$
  @param nPoints total number of points
  @param t NOT USED
  @param x input array
  @param f output array, vector

  The exact solution of
    
   \f[ u_t + \nabla \cdot (b u - a \nabla u) = 0\f]
    
    on an infinite domain with dirac initial data
    
   \f[u_0 = \int u_0 \delta(x - x_0)\f]

  \return status code
*/
int LinearAD_DiracIC_totalVelocity(int *iwork, double *rwork, int nPoints, double T, double *x, double *f)
{
	int i, j, iret;
	double *f_adv = (double *)malloc(nPoints*3*sizeof(double));
	
	iret = LinearAD_DiracIC_advectiveVelocity(iwork, rwork, nPoints, T, x, f_adv);
	iret += LinearAD_DiracIC_diffusiveVelocity(iwork, rwork, nPoints, T, x, f);

	for (i=0; i<nPoints; i++)
	{
	    for (j=0; j<3; j++)
	    {
		 f[i*3+j] += f_adv[i*3+j];
	    }
	}
	free(f_adv);

	return iret;
}
/**
  \brief Linear Advection-Diffusion Steady State
  @param iwork NOT USED
  @param rwork[0] \f$ b \f$
  @param rwork[1] \f$ a \f$
  @param nPoints - total number of points
  @param t NOT USED
  @param x input array
  @param u output array

  \f[(bu - au_x)_x = 0\f]
  \f[ u(0) = 1\f]
  \f[ u(1) = 0\f]

  
  \return status code
*/

int LinearAD_SteadyState(int *iwork, double *rwork, int nPoints, double t, double *X, double *u)
/*
    The exact solution for
    (bu - au_x)_x = 0
    u(0) = 1
    u(1) = 0
*/
{
  int i;
  double  b=rwork[0], a=rwork[1];
  double x, D=0.0, C;

  if (b != 0.0)
  {
	D = (1.0/(exp(b/a)-1.0));
  }
  C = (-D)*exp(b/a);

  for (i=0;i<nPoints;i++)
  {
      x = X[i*3+0];
      if (D != 0.0)
      {
        u[i] = (-D)*exp(b*x/a) - C;
      }
      else
      {
	u[i] = 1.0 - x;
      }
  } 
 
  return 0;
}

/**
  \brief Linear Avection Diffusion Reaction Decay Dirac Initial Condition
  @param iwork NOT USED
  @param rwork[0:2] \f$ b[0:2] \f$
  @param rwork[3] \f$ n \f$
  @param rwork[4] \f$ a \f$
  @param rwork[5] \f$ tStart \f$
  @param rwork[6] \f$ u_0 \f$
  @param rwork[7:9] \f$ x_0[0:2] \f$
  @param rwork[10] \f$ c \f$
  @param nPoints total number of points
  @param t NOT USED
  @param x input array
  @param u output array

The exact solution of
    
    \f[ u_t + \nabla \cdot(bu - a \nabla u) + cu= 0 \f]
    
    on an infinite domain with Dirac initial data.
    Also returns the fluxes (by inheritance).

  \return status code
*/
int LinearADR_Decay_DiracIC(int *iwork, double *rwork, int nPoints, double T, double *x, double *u)
/*
The exact solution of
    
    u_t + \deld(bu - a \grad u) + cu= 0
    
    on an infinite domain with Dirac initial data.
    Also returns the fluxes (by inheritance).
*/
{
	int i, iret;
	double c=rwork[10], tStart=rwork[5]; 
	double t = T + tStart;

	iret = LinearAD_DiracIC(iwork, rwork, nPoints, T, x, u);
	for (i=0; i<nPoints; i++)
	{
	     u[i] *= exp(-c*t);
	}

	return iret;
}
/**
  \brief Linear Avection Diffusion Reaction Decay Dirac Initial Condition (dr)
  @param iwork NOT USED
  @param rwork[0:2] \f$ b[0:2] \f$
  @param rwork[3] \f$ n \f$
  @param rwork[4] \f$ a \f$
  @param rwork[5] \f$ tStart \f$
  @param rwork[6] \f$ u_0 \f$
  @param rwork[7:9] \f$ x_0[0:2] \f$
  @param rwork[10] \f$ c \f$
  @param nPoints total number of points
  @param t NOT USED
  @param x input array
  @param u output array
  @param dr output array

The exact solution of
    
    \f[ u_t + \nabla \cdot(bu - a \nabla u) + cu= 0 \f]
    
    on an infinite domain with Dirac initial data.
    Also returns the fluxes (by inheritance).

  \return status code
*/
int LinearADR_Decay_DiracIC_dr(int *iwork, double *rwork, int nPoints, double T, double *x, double *u, double *dr)
{
  int i, iret;
  double c=rwork[10];
	
  iret = LinearADR_Decay_DiracIC(iwork, rwork, nPoints, T, x, u);

  for (i=0; i<nPoints; i++)
  {
       dr[i] = c;
  }

  return iret;
}
/**
  \brief Linear Avection Diffusion Reaction Decay Dirac Initial Condition (reaction)
  @param iwork NOT USED
  @param rwork[0:2] \f$ b[0:2] \f$
  @param rwork[3] \f$ n \f$
  @param rwork[4] \f$ a \f$
  @param rwork[5] \f$ tStart \f$
  @param rwork[6] \f$ u_0 \f$
  @param rwork[7:9] \f$ x_0[0:2] \f$
  @param rwork[10] \f$ c \f$
  @param nPoints total number of points
  @param t NOT USED
  @param x input array
  @param u output array
  @param r output array

The exact solution of
    
    \f[ u_t + \nabla \cdot(bu - a \nabla u) + cu= 0 \f]
    
    on an infinite domain with Dirac initial data.
    Also returns the fluxes (by inheritance).

  \return status code
*/
int LinearADR_Decay_DiracIC_r(int *iwork, double *rwork, int nPoints, double T, double *x, double *u, double *r)
{
  int i, iret;
  double c=rwork[10];

  iret = LinearADR_Decay_DiracIC(iwork, rwork, nPoints, T, x, u);
	
  for (i=0; i<nPoints; i++)
  {
       r[i] = u[i] * c;
  }

  return iret;
}
/**
  \brief Linear Avection Diffusion Reaction Sine function
  @param iwork NOT USED
  @param rwork[0:2] \f$ \omega[0:2] \f$
  @param rwork[3] \f$ \omega_0 \f$
  @param nPoints total number of points
  @param t NOT USED
  @param x input array
  @param u output array

  An exact solution and source term for
    
    \f[ \nabla \cdot (\vec b u - \ddot a \nabla u) + c u + d = 0 \f]
    
  where

   \f[ u(x)  = sin(\vec \omega \cdot \vec x + \omega_0) = sin(Ax - b)\f]
   \f[r(u,x) = - ((\ddot a \vec \omega) \cdot \vec \omega) u - (\vec b \cdot \vec \omega) cos(Ax - b)\f]
   \f[       = c u + D cos(Ax - b) \f]
   \f[       = cu + d \f]
           
  also returns the advective, diffusive, and total velocity at a point.

  \return status code
*/
int LinearADR_Sine(int *iwork, double *rwork, int nPoints, double t, double *x, double *u)
/*
    An exact solution and source term for
    
    \deld (\vec b u - \ten a \grad u) + c u + d = 0
    
    where
    
    u(x) = sin[\gvec \omega \cdot \vec x + \omega_0) ] = sin(Ax - b)
    r(u,x) = - [(\ten a \gvec \omega] \cdot \gvec \omega) u
             - (\vec b \cdot \omega) cos(Ax - b)
           = c u + D cos(Ax - b)
           = cu + d
           
    also returns the advective, diffusive, and total flux at a point.
*/
{
  int i=0;
  double omega[3]={rwork[0],rwork[1],rwork[2]};
  double omega0=rwork[3];
  double Y;

  for (i=0;i<nPoints;i++)
  {
	  Y = omega[0]*x[i*3+0] + omega[1]*x[i*3+1] + omega[2]*x[i*3+2] + omega0;
	  u[i] = sin(Y);
  }

  return 0;
}
/**
  \brief Linear Avection Diffusion Reaction Sine function (advective velocity)
  @param iwork NOT USED
  @param rwork[0:2] \f$ \omega[0:2] \f$
  @param rwork[3] \f$ \omega_0 \f$
  @param rwork[4:6] \f$ b[0:2] \f$
  @param rwork[7:15] \f$ a[0:3,0:3] \f$
  @param rwork[16] \f$ c \f$
  @param nPoints total number of points
  @param t NOT USED
  @param x input array
  @param f output array (vector)

  An exact solution and source term for
    
    \f[ \nabla \cdot (\vec b u - \ddot a \nabla u) + c u + d = 0 \f]
    
  where

   \f[ u(x)  = sin(\vec \omega \cdot \vec x + \omega_0) = sin(Ax - b)\f]
   \f[r(u,x) = - ((\ddot a \vec \omega) \cdot \vec \omega) u - (\vec b \cdot \vec \omega) cos(Ax - b)\f]
   \f[       = c u + D cos(Ax - b) \f]
   \f[       = cu + d \f]
           
   returns the advective velocity at a point.

  \return status code
*/
int LinearADR_Sine_advectiveVelocity(int *iwork, double *rwork, int nPoints, double t, double *x, double *f)
{
	int i, j, iret;
	double b[3]={rwork[4],rwork[5],rwork[6]};
 	double *u = (double *)malloc(nPoints * sizeof(double));

        /*for (i=0;i<nPoints;i++)
	{
          printf(" %d x = %f, %f, %f\n",i,x[i*3+0],x[i*3+1],x[i*3+2]);
        }
        */
	iret = LinearADR_Sine(iwork, rwork, nPoints, t, x, u);

	for (i=0; i<nPoints; i++)
	{
		for (j=0; j<3; j++)
		{
                  f[i*3+j] = u[i]*b[j];
		}
	}
        free( u );

	return iret;
}
/**
  \brief Linear Avection Diffusion Reaction Sine function (diffusive velocity)
  @param iwork NOT USED
  @param rwork[0:2] \f$ \omega[0:2] \f$
  @param rwork[3] \f$ \omega_0 \f$
  @param rwork[4:6] \f$ b[0:2] \f$
  @param rwork[7:15] \f$ a[0:3,0:3] \f$
  @param rwork[16] \f$ c \f$
  @param nPoints total number of points
  @param t NOT USED
  @param x input array
  @param f output array (vector)

  An exact solution and source term for
    
    \f[ \nabla \cdot (\vec b u - \ddot a \nabla u) + c u + d = 0 \f]
    
  where

   \f[ u(x)  = sin(\vec \omega \cdot \vec x + \omega_0) = sin(Ax - b)\f]
   \f[r(u,x) = - ((\ddot a \vec \omega) \cdot \vec \omega) u - (\vec b \cdot \vec \omega) cos(Ax - b)\f]
   \f[       = c u + D cos(Ax - b) \f]
   \f[       = cu + d \f]
           
   returns the diffusive velocity at a point.

  \return status code
*/

int LinearADR_Sine_diffusiveVelocity(int *iwork, double *rwork, int nPoints, double t, double *x, double *f)
{
	int i, j, iret;
	double a[3][3]={ {rwork[7],rwork[8],rwork[9]},{rwork[10],rwork[11],rwork[12]},{rwork[13],rwork[14],rwork[15]} };
	double *du = (double *)malloc(nPoints*3*sizeof(double));

	iret = LinearADR_Sine_du(iwork, rwork, nPoints, t, x, du);

// matrix multiplication of matrix a(3x3) and a vector x of nPoints.
 	for (i=0; i<nPoints; i++)
	{
	    for (j=0; j<3; j++)
	    {
		f[i*3+j] = -( a[j][0]*du[i*3+j] + a[j][1]*du[i*3+j] + a[j][2]*du[i*3+j] );
	    }
	}
        free(du);

	return iret;
}
/**
  \brief Linear Avection Diffusion Reaction Sine function (dr)
  @param iwork NOT USED
  @param rwork[0:2] \f$ \omega[0:2] \f$
  @param rwork[3] \f$ \omega_0 \f$
  @param rwork[4:6] \f$ b[0:2] \f$
  @param rwork[7:15] \f$ a[0:3,0:3] \f$
  @param rwork[16] \f$ c \f$
  @param nPoints total number of points
  @param t NOT USED
  @param x input array
  @param u output array
  @param dr output array

  An exact solution and source term for
    
    \f[ \nabla \cdot (\vec b u - \ddot a \nabla u) + c u + d = 0 \f]
    
  where

   \f[ u(x)  = sin[\vec \omega \cdot \vec x + \omega_0) ] = sin(Ax - b) \f
   \f[r(u,x) = - ((\ddot a \vec \omega) \cdot \vec \omega) u - (\vec b \cdot \vec \omega) cos(Ax - b)\f]
   \f[       = c u + D cos(Ax - b) \f]
   \f[       = cu + d \f]
           
  also returns the advective, diffusive, and total velocity at a point.

  \return status code
*/
int LinearADR_Sine_dr(int *iwork, double *rwork, int nPoints, double t, double *x, double *u, double *dr)
{
	int i, iret;
        double c=rwork[16];
	
        iret = LinearADR_Sine(iwork, rwork, nPoints, t, x, u);

	for (i=0; i<nPoints; i++)
	{
	   dr[i] = c;
	}

	return iret;
}
/**
  \brief Linear Avection Diffusion Reaction Sine function (du)
  @param iwork NOT USED
  @param rwork[0:2] \f$ \omega[0:2] \f$
  @param rwork[3] \f$ \omega_0 \f$
  @param nPoints total number of points
  @param t NOT USED
  @param x input array
  @param f output array (vector)

  An exact solution and source term for
    
    \f[ \nabla \cdot (\vec b u - \ddot a \nabla u) + c u + d = 0 \f]
    
  where

   \f[ u(x)  = sin(\vec \omega \cdot \vec x + \omega_0) = sin(Ax - b)\f]
   \f[r(u,x) = - ((\ddot a \vec \omega) \cdot \vec \omega) u - (\vec b \cdot \vec \omega) cos(Ax - b)\f]
   \f[       = c u + D cos(Ax - b) \f]
   \f[       = cu + d \f]
           
  also returns the advective, diffusive, and total velocity at a point.

  \return status code
*/
int LinearADR_Sine_du(int *iwork, double *rwork, int nPoints, double t, double *x, double *f)
{
  int i, j;
  double omega[3]={rwork[0],rwork[1],rwork[2]};
  double omega0=rwork[3];
  double Y;

  for (i=0;i<nPoints;i++)
  {
      Y = cos( omega[0]*x[i*3+0] + omega[1]*x[i*3+1] + omega[2]*x[i*3+2] + omega0 );
      for (j=0; j<3; j++)
      {
	  f[i*3+j] = omega[j] * Y;
      }
  }

  return 0;
}
/**
  \brief Linear Avection Diffusion Reaction Sine function (reaction)
  @param iwork NOT USED
  @param rwork[0:2] \f$ \omega[0:2] \f$
  @param rwork[3] \f$ \omega_0 \f$
  @param rwork[4:6] \f$ b[0:2] \f$
  @param rwork[7:15] \f$ a[0:3,0:3] \f$
  @param rwork[16] \f$ c \f$
  @param nPoints total number of points
  @param t NOT USED
  @param x input array
  @param u output array
  @param r output array

  An exact solution and source term for
    
    \f[ \nabla \cdot (\vec b u - \ddot a \nabla u) + c u + d = 0 \f]
    
  where

   \f[ u(x)  = sin(\vec \omega \cdot \vec x + \omega_0) = sin(Ax - b)\f]
   \f[r(u,x) = - ((\ddot a \vec \omega) \cdot \vec \omega) u - (\vec b \cdot \vec \omega) cos(Ax - b)\f]
   \f[       = c u + D cos(Ax - b) \f]
   \f[       = cu + d \f]
           
  also returns the advective, diffusive, and total velocity at a point.

  \return status code
*/

int LinearADR_Sine_r(int *iwork, double *rwork, int nPoints, double t, double *x, double *u, double *r)
{
  int i, j, iret;
  double omega[3]={rwork[0],rwork[1],rwork[2]};
  double omega0=rwork[3];
  double b[3]={rwork[4],rwork[5],rwork[6]};
  double a[3][3]={ {rwork[7],rwork[8],rwork[9]},{rwork[10],rwork[11],rwork[12]},{rwork[13],rwork[14],rwork[15]}};
  double c=rwork[16]; 
  double Y, D, E;
  double MM[3];

  iret = LinearADR_Sine(iwork, rwork, nPoints, t, x, u);

// matrix multiplication of 3x3 matrix and 3x1 vector
  for (j=0; j<3; j++)
  {
	  MM[j] = a[j][0]*omega[0] + a[j][1]*omega[1] + a[j][2]*omega[2];
  }
  E = - ( MM[0]*omega[0] + MM[1]*omega[1] + MM[2]*omega[2] ) - c;

// dot product of vector b and vector omega
  D = -( b[0]*omega[0] + b[1]*omega[1] + b[2]*omega[2] );

  for (i=0; i<nPoints; i++)
  {
	  Y = omega[0]*x[i*3+0] + omega[1]*x[i*3+1] + omega[2]*x[i*3+2] + omega0;
	  r[i] = c*u[i] + D*cos(Y) + E*sin(Y);
  }

  return iret;
}
/**
  \brief Linear Avection Diffusion Reaction Sine function (total velocity)
  @param iwork NOT USED
  @param rwork[0:2] \f$ \omega[0:2] \f$
  @param rwork[3] \f$ \omega_0 \f$
  @param rwork[4:6] \f$ b[0:2] \f$
  @param rwork[7:15] \f$ a[0:3,0:3] \f$
  @param rwork[16] \f$ c \f$
  @param nPoints total number of points
  @param t NOT USED
  @param x input array
  @param f output array (vector)

  An exact solution and source term for
    
    \f[ \nabla \cdot (\vec b u - \ddot a \nabla u) + c u + d = 0 \f]
    
  where

   \f[ u(x)  = sin(\vec \omega \cdot \vec x + \omega_0) = sin(Ax - b)\f]
   \f[r(u,x) = - ((\ddot a \vec \omega) \cdot \vec \omega) u - (\vec b \cdot \vec \omega) cos(Ax - b)\f]
   \f[       = c u + D cos(Ax - b) \f]
   \f[       = cu + d \f]
           
   returns the total velocity at a point.

  \return status code
*/
int LinearADR_Sine_totalVelocity(int *iwork, double *rwork, int nPoints, double t, double *x, double *f)
{
	int i, j, iret;
	double *f_adv = (double *)malloc(nPoints*3*sizeof(double));

	iret = LinearADR_Sine_advectiveVelocity(iwork, rwork, nPoints, t, x, f_adv);
	iret += LinearADR_Sine_diffusiveVelocity(iwork, rwork, nPoints, t, x, f);

	for (i=0;i<nPoints;i++)
	{
		for (j=0; j<3; j++)
		{
			f[i*3+j] += f_adv[i*3+j];
		}
	}
	free(f_adv);

	return iret;
}
/**
  \brief Nonlinear Advection-Diffusion Steady State
  @param iwork NOT USED
  @param rwork[0] \f$ q \f$
  @param rwork[1] \f$ r \f$
  @param rwork[2] \f$ b \f$
  @param rwork[3] \f$ a \f$
  @param nPoints - total number of points
  @param t NOT USED
  @param x input array
  @param u output array

  \f[(bu^q - a(u^r)_x)_x = 0\f]
  \f[ u(0) = 1\f]
  \f[ u(1) = 0\f]

  
  \return status code
*/
/********************************************************************/
int NonlinearAD_SteadyState(int *iwork, double *rwork, int nPoints, double t, double *X, double *u)
/*
    (bu^q - a(u^r)_x)_x = 0
    u(0) = 1
    u(1) = 0
*/
{
	int i, iret=0;
	int q=iwork[0], r=iwork[1];
	double b=rwork[0], a=rwork[1];
	double x, lu, f0, ff, fC;
	double C, rtmC=0.0, Ctmp, dC, nlC=0.0, nlD=0.0;

    if (q==2 && r==1)
	{
		if (b != 0.0)
		{
			printf( "Solving for sqrt(-C) for q=2, r=1\n");
			rtmC = sqrt(1.5);
			while ( fabs( f(rtmC,b,a,q,r) ) > 1.0E-8 )
			{
				rtmC -= ( f(rtmC,b,a,q,r) / df(rtmC,b,a,q,r) );
			}
			printf( "sqrt(-C)= %lf \n",rtmC);
		}
	}
    else if (q==1 && r==2)
	{
		printf( "\nSolving for C in q=1,r=2\n" );
		C = 1.0 + 1.0E-10;
		f0 = f(C,b,a,q,r);
		printf("f0 = %lf\n", f0);
		while ( fabs(f(C,b,a,q,r)) > (1.0E-7*fabs(f0) + 1.0E-7))
		{
			dC = -f(C,b,a,q,r)/df(C,b,a,q,r);
                  printf("dc = %lf\n",dC);
                  Ctmp = C + dC;
                  while ( (fabs(ff=f(Ctmp,b,a,q,r)) > 0.99*fabs(fC=f(C,b,a,q,r))) || (Ctmp <= 1.0) )
			{
				printf("f(%lf) = %lf\n", Ctmp, ff);
				printf("f(%lf) = %lf\n", C, fC);
				printf( "ls\n" );
                                dC *= 0.9 ;
                                Ctmp = C + dC ;
			}
                  printf( "out\n" );
                  printf( "%lf \n",Ctmp);
// are we printing the new values from the functions f & df?
                  printf(" f(%lf) = %lf\n",Ctmp, f(Ctmp,b,a,q,r));
                  printf("df(%lf) = %lf\n",Ctmp, df(Ctmp,b,a,q,r));
                  C=Ctmp;
		}
                printf("C = %lf\n",C);
                nlC = C;
                nlD = 0.5*(2.0*C*log(C*(C-1)) - 4.0*C + 2.0 - b/a);
                printf( "D = %lf\n",nlD);
	}
     else
	{
		printf("q,r not implemented\n");
                return 0;
	}
// Main Loop for uOfX

     if(q==1 && r==2)
     {
       iret = LinearAD_SteadyState(iwork, rwork, nPoints, t, X, u) ;
     }

     for (i=0; i<nPoints; i++)
	{
	     x = X[i*3+0];
		
	     if (q==2 && r==1)
		{
		    if (b != 0.0)
		    {
			u[i] = rtmC * tanh( (-b*rtmC/a)*(x - 1.0) );
		    }
		    else
		    {
			u[i] = 1.0 - x;
		    }
		}
	     else if (q==1 && r==2)
		{
		    lu = u[i];
		    f0 = uOfX_f(a, b, nlC, nlD, x, lu) ;
		    while ( fabs(uOfX_f(a, b, nlC, nlD, x, lu)) > (1.0E-6*fabs(f0) + 1.0E-6) )
		    {
			lu -= uOfX_f(a, b, nlC, nlD, x, lu) / uOfX_df(nlC, lu) ;
		    }
		    u[i] = lu;
		}
	      else
		{	
			printf("q,r not implemented");	
		}
         }
      return iret;
}

/**
  \brief Non Linear Avection Diffusion Reaction Decay Dirac Initial Condition
  @param iwork NOT USED
  @param rwork[0:2] \f$ b[0:2] \f$
  @param rwork[3] \f$ n \f$
  @param rwork[4] \f$ a \f$
  @param rwork[5] \f$ tStart \f$
  @param rwork[6] \f$ u_0 \f$
  @param rwork[7:9] \f$ x_0[0:2] \f$
  @param rwork[10] \f$ c \f$
  @param rwork[11] \f$ d \f$
  @param nPoints total number of points
  @param t NOT USED
  @param x input array
  @param u output array
  @param dr output array

The approximate analytical solution of
    
    \f[ u_t + \nabla \cdot(bu - a \nabla u) + cu^d= 0 \f]
    
    on an infinite domain with Dirac initial data.
    Also returns the fluxes (by inheritance).

  \return status code
*/
int NonlinearADR_Decay_DiracIC(int *iwork, double *rwork, int nPoints, double T, double *x, double *u)
/*
The approximate analytical solution of
    
    u_t + \deld(bu - a \grad u) + cu^d= 0
    
    on an infinite domain with Dirac initial data.
    Also returns the fluxes (by inheritance).
*/
{
    int i, iret;
    double c=rwork[10], d=rwork[11], tStart=rwork[5];
    double t = T + tStart;

    iret = LinearAD_DiracIC(iwork, rwork, nPoints, T, x, u);

    for (i=0; i<nPoints; i++)
    {
	 if (u[i] > 0.0)
	 {
	     u[i] *= exp( -(2.0*c*t*pow(u[i],(d-1.0)))/(d+1.0) );
	 }
    }

    return iret;
}
/**
  \brief Non Linear Avection Diffusion Reaction Decay Dirac Initial Condition (dr)
  @param iwork NOT USED
  @param rwork[0:2] \f$ b[0:2] \f$
  @param rwork[3] \f$ n \f$
  @param rwork[4] \f$ a \f$
  @param rwork[5] \f$ tStart \f$
  @param rwork[6] \f$ u_0 \f$
  @param rwork[7:9] \f$ x_0[0:2] \f$
  @param rwork[10] \f$ c \f$
  @param rwork[11] \f$ d \f$
  @param nPoints total number of points
  @param t NOT USED
  @param x input array
  @param u output array
  @param dr output array

The approximate analytical solution of
    
    \f[ u_t + \nabla \cdot(bu - a \nabla u) + cu^d= 0 \f]
    
    on an infinite domain with Dirac initial data.
    Also returns the fluxes (by inheritance).

  \return status code
*/
int NonlinearADR_Decay_DiracIC_dr(int *iwork, double *rwork, int nPoints, double T, double *x, double *u, double* dr)
{
    int i, iret;
    double c=rwork[10], d=rwork[11];
	
    iret = NonlinearADR_Decay_DiracIC(iwork, rwork, nPoints, T, x, u);

    for (i=0; i<nPoints; i++)
    {
	 dr[i] = d*c*pow(u[i],(d-1.0));
    }

    return iret;
}

/**
  \brief Non Linear Avection Diffusion Reaction Decay Dirac Initial Condition (reaction)
  @param iwork NOT USED
  @param rwork[0:2] \f$ b[0:2] \f$
  @param rwork[3] \f$ n \f$
  @param rwork[4] \f$ a \f$
  @param rwork[5] \f$ tStart \f$
  @param rwork[6] \f$ u_0 \f$
  @param rwork[7:9] \f$ x_0[0:2] \f$
  @param rwork[10] \f$ c \f$
  @param rwork[11] \f$ d \f$
  @param nPoints total number of points
  @param t NOT USED
  @param x input array
  @param u output array
  @param r output array

The approximate analytical solution of
    
    \f[ u_t + \nabla \cdot(bu - a \nabla u) + cu^d= 0 \f]
    
    on an infinite domain with Dirac initial data.
    Also returns the fluxes (by inheritance).

  \return status code
*/
int NonlinearADR_Decay_DiracIC_r(int *iwork, double *rwork, int nPoints, double T, double *x, double *u, double* r)
{
    int i, iret;
    double c=rwork[10], d=rwork[11];

    iret = NonlinearADR_Decay_DiracIC(iwork, rwork, nPoints, T, x, u);
	
    for (i=0; i<nPoints; i++)
    {
	 r[i] = c*pow(u[i],d);
    }

    return iret;
}

/**
  \brief Nonlinear Differential-algebraic equations
  @param iwork NOT USED
  @param rwork[0] \f$ a \f$
  @param rwork[1] \f$ p \f$
  @param nPoints - total number of points
  @param t NOT USED
  @param x input array
  @param u output array

  \f[u_t = - a*max(u,0)^p\f]
  \f[u(0) = 1\f]
  
  \return status code
*/
int NonlinearDAE(int *iwork, double *rwork, int nPoints, double T, double *x, double *u)
/*
    The exact solution of
    
    u_t = - a max(u,0)^p
    u(0) = 1
*/
{
  int i;
  double q, t, a=rwork[0], p=rwork[1];

  for (i=0;i<nPoints;i++)
    {
      t = x[i*3+0];
      
      if (p == 1.0)
	{
		u[i] =  exp( -a*t );
	}
      else
	{
		q = 1.0/(1.0 - p);
		u[i] = pow( max( (1.0 - (1.0 - p)*a*t), 0.0 ), q ) ;
	}
    }

  return 0;
}
/**
  \brief Nonlinear Differential-algebraic equations
  @param iwork NOT USED
  @param rwork[0] \f$ a \f$
  @param rwork[1] \f$ p \f$
  @param nPoints - total number of points
  @param t NOT USED
  @param x input array
  @param f output array, vector

  \f[u_t = - a*max(u,0)^p\f]
  \f[u(0) = 1\f]
  
  \return status code
*/
int NonlinearDAE_f(int *iwork, double *rwork, int nPoints, double t, double *x, double *f)
/*
    The exact solution of
    
    u_t = - a max(u,0)^p
    u(0) = 1
*/
{
	int i, j;
	double a=rwork[0], p=rwork[1];

	for (i=0;i<nPoints;i++)
	{
	    for (j=0; j<3; j++)
	    {
		  f[i*3+j] = -a*pow(max(f[i*3+j],0.0),p) ;
	    }
	}

	return 0;
}
   
/**
  \brief Poiseuille Flow between two parallel fixed plates with constant seperation (width).
  @param iwork NOT USED
  @param rwork[0] width between the plates, \f$ h \f$
  @param rwork[1] viscosity of fluid, \f$ \mu \f$
  @param rwork[2] the pressure gradient, \f$ \nabla p \f$ (neg)
  @param rwork[3] the rate of flow per unit width, \f$ q \f$
  @param rwork[4] x-axis offset, \f$ xs \f$
  @param rwork[5] y-axis offset, \f$ ys \f$
  @param rwork[6] z-axis offset, \f$ zs \f$
  @param nPoints total number of points
  @param t NOT USED
  @param x input array 
  @param u output array

  x-momentum equation
  \f[ (\nabla p + \mu \nabla^2 u)_x = 0 \f]
  \f[ u(y) = 4  u_{max} y (\frac{h-y}{h^2}) \f]
  \f[ u(0) = u(h) = 0 \f]
  \f[ u(h/2) = u_{max}  \f]
  \f[u_{max} = \frac{-\nabla p}{8 \mu} h^2 = \frac{3}{2} \frac{q}{h} \f] 
  
  \f$ u_{max} \f$ = maximum velocity at the centerline\n
  Axis of origin: at bottom plate\n
  
  \return status code
*/

int PlanePoiseuilleFlow_u(int *iwork, double *rwork, int nPoints, double t, double *x, double *u)
{
/*
  Poiseuille Flow between two parallel fixed plates with constant seperation (width)

  h  = width between plates
  GradP = Pressure gradient (neg.)
  mu = viscosity
  q = rate of flow per unit width 
  xs,ys,zs = axis offset
  Axis of origin: bottom plate
  Vmax = maximum velocity at the centerline

  -\delp + mu \del^2(\u_x) = 0

  given either -GradP: u_{max} = -GradP / (2 * mu) * (h/2)^2
         or     q    : u_{max} = 3/2 * (q / h)

  u(y) = 4 * u_{max} * Y * (h-Y) / h^2

  u(0) = u(h) = 0
  u(h/2) = u_{max}  
  u_{max} = -GradP / (2 * mu) * (h/2)<sup>2</sup>\n
  or \n
  u_{max} = 3/2 * (q / h) 
*/
	int i;
	double h=rwork[0], mu=rwork[1];
        double GradP = rwork[2], q = rwork[3];
        double ys=rwork[5];
        double umax;
        double Y;

        if(GradP != 0.0)
          {
             umax = -GradP *(1.0/(2.0*mu)) * (h*h/4.0);
          }
        else if(q != 0.0)
          {
             umax = (3.0/2.0)*q/h;
          }
        else
          {
             printf("Please enter input values for either q(Flow Rate per unit depth) or GradP(Pressure Gradient)\n");
             exit (0);
          }

 	for (i=0; i<nPoints; i++)
 	{
             Y = x[i*3+1] - ys;
	       
	     u[i] = (4.0*umax*Y*(h-Y))/(h*h);
 	}

	return 0;
}
/**
  \brief Poiseuille Flow through a circular pipe.
  @param iwork NOT USED
  @param rwork[0] radius of pipe, \f$ R \f$
  @param rwork[1] viscosity of fluid, \f$ \mu \f$
  @param rwork[2] the pressure gradient, \f$ \nabla p \f$ (neg)
  @param rwork[3] the rate of flow,  \f$ Q \f$
  @param rwork[4] x-axis offset, \f$ xs \f$
  @param rwork[5] y-axis offset, \f$ ys \f$
  @param rwork[6] z-axis offset, \f$ zs \f$
  @param nPoints total number of points
  @param t NOT USED
  @param x input array 
  @param u output array

  x-momentum equation
  \f[ (\nabla p + \mu \nabla^2 u)_x = 0 \f]
  \f[ u(r) = u_{max} (1 - \frac{r^2}{R^2}) \f]
  \f[ u(0) = u_{max} \f]
  \f[ u(R) = 0 \f]
  \f[ u_{max} = -\nabla p \frac{R^2}{4 \mu}= \frac{2Q}{\pi R^2}\f]
  
  \f$ u_{max} \f$  = maximum velocity at the centerline\n
  Axis of origin: at center of pipe\n
  

  \return status code

*/

int PoiseuillePipeFlow(int *iwork, double *rwork, int nPoints, double t, double *x, double *u)
{
/*
  Poiseuille Flow through a circular pipe
  Velocity in the direction of flow (z-dir)

  R  = radius of pipe
  GradP = Pressure gradient (neg.)
  mu = viscosity   
  Q = Rate of Flow 
  Axis of origin: center of pipe
  xs,ys,zs = axis offset

  x-momentum equation in cylindrical coordinates
  Vr = Vtheta = 0
  -\\delp + mu \\del^2(Vx) = 0

  No transformation of coordinates is required for x.
  
  Vmax = -GradP * R^2/(4 * mu)\n
         or  \n  
  Vmax = 2 * Q / Pipe Area\n

  Vx(r) = Vmax *  (1 - r^2/R^2)

  Vx(0) = Vmax
  Vx(R) = 0
  Vmax = -GradP * R<sup>2</sup>/(4 * mu)\n
  or\n
  Vmax = 2 * Q / Pipe Area\n

*/
	int i;
	double r=rwork[0], mu = rwork[1];
        double GradP=rwork[2], Q=rwork[3];
        double ys=rwork[5];
        double umax;
        double Y;

        if(GradP != 0.0)
          {
             umax = -GradP *(r*r/(4.0*mu));
          }
        else if(Q != 0.0)
          {
             umax = 2.0*Q/(PI * r*r);
          }
        else
          {
             printf("Please enter input values for either Q(Flow Rate) or GradP(Pressure Gradient)\n");
             exit (0);
          }
 
 	for (i=0; i<nPoints; i++)
 	{
             Y = x[i*3+1] - ys;
		
	     u[i] = umax*(1.0 - ((Y*Y)/(r*r)));
 	}

	return 0;
}
/**
  \brief Poiseuille Flow through a circular pipe.
  @param iwork NOT USED
  @param rwork[0] radius of pipe, \f$ R \f$
  @param rwork[1] viscosity of fluid, \f$ \mu \f$
  @param rwork[2] the pressure gradient, \f$ \nabla p \f$ (neg)
  @param rwork[3] the rate of flow,  \f$ Q \f$
  @param rwork[3] the length of pipe,  \f$ L \f$
  @param rwork[4] x-axis offset, \f$ xs \f$
  @param rwork[5] y-axis offset, \f$ ys \f$
  @param rwork[6] z-axis offset, \f$ zs \f$
  @param nPoints total number of points
  @param t NOT USED
  @param x input array 
  @param u output array

  \f[ p - p_0 = -\nabla p L(1-x)\f]

  
  \f$ p_0 =p(L)\f$ = 0, the atmospheric pressure at the outflow\n
  Axis of origin: at center of pipe\n
  
  
  \return status code
*/

int PoiseuillePipeFlow_P(int *iwork, double *rwork, int nPoints, double t, double *x, double *u)
{
/*
  Poiseuille Flow through a circular pipe
  Velocity in the direction of flow (z-dir)

  R  = radius of pipe
  GradP = Pressure gradient (neg.)
  mu = viscosity   
  Q = Rate of Flow 
  Axis of origin: center of pipe
  xs,ys,zs = axis offset
  P = pressure

  x-momentum equation in cylindrical coordinates
  Vr = Vtheta = 0
  -GradP + mu \del^2(Vx) = 0
 
  For Q given:\n
  GradP = (8 * Q * mu) / (R<sup>2</sup> * Pipe Area)
  assumuption: P1 = atmospheric pressure = 0.0
*/
	int i;
	double r=rwork[0], mu=rwork[1];
        double GradP=rwork[2], Q=rwork[3], L=rwork[4];
        double xs=rwork[5];
        double X;


        if(GradP == 0.0 && Q != 0.0)
          {
             GradP = - (8.0*Q*mu) / (PI*pow(r,4.0));
          }
        else if(GradP == 0.0)
          {
             printf("Please enter input values for either Q(Flow Rate) OR GradP(Pressure Gradient)\n");
             exit(0);
          }

 	for (i=0; i<nPoints; i++)
 	{
             X = x[i*3+0] - xs;
		
             u[i] = -GradP*L*(1.0-X);
 	}
 
	return 0;
}
/**
  \brief Poisson Exponential Equation 1D
  @param iwork NOT USED
  @param rwork[0] \f$ K \f$
  @param nPoints total number of points
  @param t NOT USED
  @param x input array
  @param u output array

  \f[-u_{xx} - f = 0\f]
  \f[u = K x(1-x)y(1-y)z(1-z)e^{x^2 + y^2 + z^2}\f]
  \f[f = -K \{[y(1-y)z(1-z)][4x^3 - 4x^2 + 6x - 2]+
              [x(1-x)z(1-z)][4y^3 - 4y^2 + 6y - 2]+
              [x(1-x)y(1-y)][4z^3 - 4z^2 + 6z - 2]\}e^{x^2 + y^2 + z^2}\f]

  \return status code
*/

int poissonsEquationExp1D(int *iwork, double *rwork, int nPoints, double t, double *X, double *u)
{
/*
    -u_{xx} - f = 0
    u = K x(1-x)y(1-y)z(1-z)e^{x^2 + y^2 + z^2}
    f = -K {\> [y(1-y)z(1-z)][4x^3 - 4x^2 + 6x - 2]+
           [x(1-x)z(1-z)][4y^3 - 4y^2 + 6y - 2]+
           [x(1-x)y(1-y)][4z^3 - 4z^2 + 6z - 2]}e^{x^2 + y^2 + z^2} 
*/

  int i=0;
  double K=rwork[0];
  double x;

  for (i=0;i<nPoints;i++)
    {
      x = X[i*3+0];      
      u[i] = K*x*(1.0-x)*exp(x*x);
    }

  return 0;
}
/**
  \brief Poisson Exponential Equation 2D
  @param iwork NOT USED
  @param rwork[0] \f$ K \f$
  @param nPoints total number of points
  @param t NOT USED
  @param x input array
  @param u output array

  \f[-u_{xx} - f = 0\f]
  \f[u = K x(1-x)y(1-y)z(1-z)e^{x^2 + y^2 + z^2}\f]
  \f[f = -K \{[y(1-y)z(1-z)][4x^3 - 4x^2 + 6x - 2]+
              [x(1-x)z(1-z)][4y^3 - 4y^2 + 6y - 2]+
              [x(1-x)y(1-y)][4z^3 - 4z^2 + 6z - 2]\}e^{x^2 + y^2 + z^2}\f]

  \return status code
*/
int poissonsEquationExp2D(int *iwork, double *rwork, int nPoints, double t, double *X, double *u)
{
  int i=0;
  double K=rwork[0];
  double x, y;

  for (i=0;i<nPoints;i++)
    {
      x = X[i*3+0];      
      y = X[i*3+1];      
      u[i] = K*x*(1.0-x)*y*(1.0-y)*exp(x*x + y*y);
    }
  return 0;
}
/**
  \brief Poisson Exponential Equation 3D
  @param iwork NOT USED
  @param rwork[0] \f$ K \f$
  @param nPoints total number of points
  @param t NOT USED
  @param x input array
  @param u output array

  \f[-u_{xx} - f = 0\f]
  \f[u = K x(1-x)y(1-y)z(1-z)e^{x^2 + y^2 + z^2}\f]
  \f[f = -K \{[y(1-y)z(1-z)][4x^3 - 4x^2 + 6x - 2]+
              [x(1-x)z(1-z)][4y^3 - 4y^2 + 6y - 2]+
              [x(1-x)y(1-y)][4z^3 - 4z^2 + 6z - 2]\}e^{x^2 + y^2 + z^2}\f]

  \return status code
*/
int poissonsEquationExp3D(int *iwork, double *rwork, int nPoints, double t, double *X, double *u)
{
  int i=0;
  double K=rwork[0];
  double x, y, z;

  for (i=0;i<nPoints;i++)
    {
      x = X[i*3+0];      
      y = X[i*3+1];      
      z = X[i*3+2];
      u[i] = K*x*(1.0-x)*y*(1.0-y)*z*(1.0-z)*exp(x*x + y*y + z*z);
    }

  return 0;
}
/**
  \brief Poisson Exponential Equation 3D (dr)
  @param iwork NOT USED
  @param rwork[0] \f$ K \f$
  @param nPoints total number of points
  @param t NOT USED
  @param x input array
  @param u NOT USED
  @param dr output array

  \f[-u_{xx} - f = 0\f]
  \f[u = K x(1-x)y(1-y)z(1-z)e^{x^2 + y^2 + z^2}\f]
  \f[f = -K \{[y(1-y)z(1-z)][4x^3 - 4x^2 + 6x - 2]+
              [x(1-x)z(1-z)][4y^3 - 4y^2 + 6y - 2]+
              [x(1-x)y(1-y)][4z^3 - 4z^2 + 6z - 2]\}e^{x^2 + y^2 + z^2}\f]

  \return status code
*/
int poissonsEquationExp3D_dr(int *iwork, double *rwork, int nPoints, double t, double *X, double *u, double *dr)
{
	int i;
	
	for (i=0; i<nPoints; i++)
	{
	   dr[i] = 0.0;
	}

	return 0;
}

/**
  \brief Poisson Exponential Equation 1D (reaction)
  @param iwork NOT USED
  @param rwork[0] \f$ K \f$
  @param nPoints total number of points
  @param t NOT USED
  @param x input array
  @param u NOT USED
  @param r output array

  \f[-u_{xx} - f = 0\f]
  \f[u = K x(1-x)y(1-y)z(1-z)e^{x^2 + y^2 + z^2}\f]
  \f[f = -K \{[y(1-y)z(1-z)][4x^3 - 4x^2 + 6x - 2]+
              [x(1-x)z(1-z)][4y^3 - 4y^2 + 6y - 2]+
              [x(1-x)y(1-y)][4z^3 - 4z^2 + 6z - 2]\}e^{x^2 + y^2 + z^2}\f]

  \return status code
*/
int poissonsEquationExp1D_r(int *iwork, double *rwork, int nPoints, double t, double *X, double *u, double *r)
{
  int i=0;
  double K=rwork[0];
  double x;

  for (i=0;i<nPoints;i++)
  {
      x = X[i*3+0];  
    
      r[i] = K*(4.0*(1.0-x)*pow(x,3.0) - 4.0*x*x + 6.0*(1.0-x)*x - 2.0)*exp(x*x) ;
  }

  return 0;
}
/**
  \brief Poisson Exponential Equation 2D (reaction)
  @param iwork NOT USED
  @param rwork[0] \f$ K \f$
  @param nPoints total number of points
  @param t NOT USED
  @param x input array
  @param u NOT USED
  @param r output array

  \f[-u_{xx} - f = 0\f]
  \f[u = K x(1-x)y(1-y)z(1-z)e^{x^2 + y^2 + z^2}\f]
  \f[f = -K \{[y(1-y)z(1-z)][4x^3 - 4x^2 + 6x - 2]+
              [x(1-x)z(1-z)][4y^3 - 4y^2 + 6y - 2]+
              [x(1-x)y(1-y)][4z^3 - 4z^2 + 6z - 2]\}e^{x^2 + y^2 + z^2}\f]

  \return status code
*/
int poissonsEquationExp2D_r(int *iwork, double *rwork, int nPoints, double t, double *X, double *u, double *r)
{
  int i=0;
  double K=rwork[0];
  double x, y;

  for (i=0;i<nPoints;i++)
  {
      x = X[i*3+0];      
      y = X[i*3+1]; 
     
      r[i] = K*(y*(1.0-y)*(4.0*(1.0-x)*pow(x,3) - 4.0*x*x + 6.0*x*(1.0-x) - 2.0) + 
		x*(1.0-x)*(4.0*(1.0-y)*pow(y,3) - 4.0*y*y + 6.0*y*(1.0-y) - 2.0))*exp(x*x + y*y);
  }

  return 0;
}
/**
  \brief Poisson Exponential Equation 3D (reaction)
  @param iwork NOT USED
  @param rwork[0] \f$ K \f$
  @param nPoints total number of points
  @param t NOT USED
  @param x input array
  @param u NOT USED
  @param r output array

  \f[-u_{xx} - f = 0\f]
  \f[u = K x(1-x)y(1-y)z(1-z)e^{x^2 + y^2 + z^2}\f]
  \f[f = -K \{[y(1-y)z(1-z)][4x^3 - 4x^2 + 6x - 2]+
              [x(1-x)z(1-z)][4y^3 - 4y^2 + 6y - 2]+
              [x(1-x)y(1-y)][4z^3 - 4z^2 + 6z - 2]\}e^{x^2 + y^2 + z^2}\f]

  \return status code
*/
int poissonsEquationExp3D_r(int *iwork, double *rwork, int nPoints, double t, double *X, double *u, double *r)
{
  int i=0;
  double K=rwork[0];
  double x,y,z;

  for (i=0; i<nPoints; i++)
  {
      x = X[i*3+0];      
      y = X[i*3+1];      
      z = X[i*3+2];

      r[i] = K*(y*(1.0-y)*z*(1.0-z)*(4.0*(1.0-x)*pow(x,3) - 4.0*x*x + 6.0*x*(1.0-x) - 2.0) +
		x*(1.0-x)*z*(1.0-z)*(4.0*(1.0-y)*pow(y,3) - 4.0*y*y + 6.0*y*(1.0-y) - 2.0) +
		x*(1.0-x)*y*(1.0-y)*(4.0*(1.0-z)*pow(z,3) - 4.0*z*z + 6.0*z*(1.0-z) - 2.0))*exp(x*x + y*y + z*z);
  }

  return 0;
}

/**
  \brief Stokes Flow around moving Sphere.
  @param iwork NOT USED
  @param rwork[0] Sphere's x-component of velocity, \f$ vx \f$
  @param rwork[1] Sphere's y-component of velocity, \f$ vy \f$
  @param rwork[2] Sphere's z-component of velocity, \f$ vz \f$
  @param rwork[3] Sphere's radius, \f$ rs \f$
  @param rwork[4] Sphere's center x-component, \f$ xs \f$
  @param rwork[5] Sphere's center y-component, \f$ ys \f$
  @param rwork[6] Sphere's center z-component, \f$ zs \f$
  @param rwork[7] viscosity of fluid, \f$ mu \f$
  @param nPoints - total number of points
  @param t NOT USED
  @param x input array
  @param u output array

  \f[ \nabla p = \mu \nabla^2 u \f]

  \f[ p = \frac{(-3 * rs * \mu \| \vec v\| cos\theta)}{2r^2}\f]

  \return status code
*/
int STflowSphere_P(int *iwork, double *rwork, int nPoints, double t, double *x, double *u)
{
/*
 Stokes Flow around Sphere of radius r^s and center x^s, y^s, z^s.
 Velocity is u,v,w and viscosity is mu.
 Array rwork contains all the above variables.
 Input	: t, x
		: x[npoints*3] = {x,y,z} for each point
		: rwork[] = {uu,v,w,rS,xS,yS,zS,mu}
 Output : P
*/
  int i;
  double vx=rwork[0], vy=rwork[1], vz=rwork[2];
  double rS=rwork[3], xS=rwork[4], yS=rwork[5], zS=rwork[6];
  double mu=rwork[7];
  double eR[3],eTHETA[3];
  double norm_v,theta,r;

  for (i=0; i<nPoints; i++)
    {
      coords(vx,vy,vz,
             xS,yS,zS,
             &x[i*3],
             &r,
             &theta,
             &norm_v,
             eR,
             eTHETA);
      if (r >= rS)
        u[i] = (-3.0*rS*mu*norm_v*cos(theta))/(2.0*r*r);
      else
        u[i]=0.0;
    }

  return 0;
}

/**
  \brief Stokes Flow around moving Sphere.
  @param rwork[0] Sphere's x-component of velocity, \f$ vx \f$
  @param rwork[1] Sphere's y-component of velocity, \f$ vy \f$
  @param rwork[2] Sphere's z-component of velocity, \f$ vz \f$
  @param rwork[3] Sphere's radius, \f$ rs \f$
  @param rwork[4] Sphere's center x-component, \f$ xs \f$
  @param rwork[5] Sphere's center y-component, \f$ ys \f$
  @param rwork[6] Sphere's center z-component, \f$ zs \f$
  @param rwork[7] viscosity of fluid, \f$ mu \f$
  @param nPoints - total number of points
  @param t NOT USED
  @param x input array
  @param u output array

  \f[ \nabla p = \mu \nabla^2 u \f]
  \f[ v^r = \|\vec v \| cos\theta(1-\frac{3rs}{2r} + \frac{rs^3}{2r^3})\f]
  \f[ v^\theta = -\|\vec v \| sin\theta(1-\frac{3rs}{4r} - \frac{rs^3}{4r^3})\f]
  \f[v^r(r>rs) = 0\f]
  \f[v^\theta(r>rs) = 0\f]
  \f[e^R = \frac{\vec x}{\|\vec x \|} ; \{ x_x = x - xs;x_y = y-ys; x_z = z-zs\} \f]
  \f[e^\theta = (e^R - \frac{\vec v}{\|\vec v \|}) -
                [(e^R - \frac{\vec v}{\|\vec v \|}) \cdot e^R]e^R\f]

  \f[u_x = v^r e^R_x + v^\theta e^\theta_x \f]

  \return status code
*/
int STflowSphere_Vx(int *iwork, double *rwork, int nPoints, double t, double *x, double *u)
{
/*
 Stokes Flow around Sphere of radius r^s and center x^s, y^s, z^s.
 Velocity is u,v,w and viscosity is mu..
 Input	: t, x
		: x[npoints*3] = {x,y,z} for each point
		: rwork[] = {uu,v,w,rS,xS,yS,zS,mu}
 Output : Vx
*/
  int i;
  double vx=rwork[0], vy=rwork[1], vz=rwork[2];
  double rS=rwork[3], xS=rwork[4], yS=rwork[5], zS=rwork[6];
  double eR[3],eTHETA[3];
  double norm_v,theta,r;
  double vR,vTHETA;

  for (i=0; i<nPoints; i++)
    {
      coords(vx,vy,vz,
             xS,yS,zS,
             &x[i*3],
             &r,
             &theta,
             &norm_v,
             eR,
             eTHETA);
      vel(rS,norm_v,r,theta,&vR,&vTHETA);
      u[i] = vR * eR[0]+ vTHETA * eTHETA[0];
    }

  return 0;
}

/**
  \brief Stokes Flow around moving Sphere.
  @param iwork NOT USED
  @param rwork[0] Sphere's x-component of velocity, \f$ vx \f$
  @param rwork[1] Sphere's y-component of velocity, \f$ vy \f$
  @param rwork[2] Sphere's z-component of velocity, \f$ vz \f$
  @param rwork[3] Sphere's radius, \f$ rs \f$
  @param rwork[4] Sphere's center x-component, \f$ xs \f$
  @param rwork[5] Sphere's center y-component, \f$ ys \f$
  @param rwork[6] Sphere's center z-component, \f$ zs \f$
  @param rwork[7] viscosity of fluid, \f$ mu \f$
  @param nPoints - total number of points
  @param t NOT USED
  @param x input array
  @param u output array

  \f[ \nabla p = \mu \nabla^2 u \f]
  \f[ v^r = \|\vec v \| cos\theta(1-\frac{3rs}{2r} + \frac{rs^3}{2r^3})\f]
  \f[ v^\theta = -\|\vec v \| sin\theta(1-\frac{3rs}{4r} - \frac{rs^3}{4r^3})\f]
  \f[v^r(r>rs) = 0\f]
  \f[v^\theta(r>rs) = 0\f]
  \f[e^R = \frac{\vec x}{\|\vec x \|} ; \{ x_x = x - xs;x_y = y-ys; x_z = z-zs\} \f]
  \f[e^\theta = (e^R - \frac{\vec v}{\|\vec v \|}) -
                [(e^R - \frac{\vec v}{\|\vec v \|}) \cdot e^R]e^R\f]

  \f[u_y = v^r e^R_y + v^\theta e^\theta_y \f]

  \return status code
*/
int STflowSphere_Vy(int *iwork, double *rwork, int nPoints, double t, double *x, double *u)
{
/*
 Stokes Flow around Sphere of radius r^s and center x^s, y^s, z^s.
 Velocity is u,v,w and viscosity is mu.
 Array rwork contains all the above variables.
 Input	: t, x
		: x[npoints*3] = {x,y,z} for each point
		: rwork[] = {uu,v,w,rS,xS,yS,zS,mu}
 Output : Vy
*/
  int i;
  double vx=rwork[0], vy=rwork[1], vz=rwork[2];
  double rS=rwork[3], xS=rwork[4], yS=rwork[5], zS=rwork[6];
  double eR[3],eTHETA[3];
  double norm_v,theta,r;
  double vR,vTHETA;

  for (i=0; i<nPoints; i++)
    {
      coords(vx,vy,vz,
             xS,yS,zS,
             &x[i*3],
             &r,
             &theta,
             &norm_v,
             eR,
             eTHETA);
      vel(rS,norm_v,r,theta,&vR,&vTHETA);
      u[i] = vR * eR[1] + vTHETA * eTHETA[1];
    }

  return 0;
}
/**
  \brief Stokes Flow around moving Sphere.  
  @param iwork NOT USED
  @param rwork[0] Sphere's x-component of velocity, \f$ vx \f$
  @param rwork[1] Sphere's y-component of velocity, \f$ vy \f$
  @param rwork[2] Sphere's z-component of velocity, \f$ vz \f$
  @param rwork[3] Sphere's radius, \f$ rs \f$
  @param rwork[4] Sphere's center x-component, \f$ xs \f$
  @param rwork[5] Sphere's center y-component, \f$ ys \f$
  @param rwork[6] Sphere's center z-component, \f$ zs \f$
  @param rwork[7] viscosity of fluid, \f$ mu \f$
  @param nPoints - total number of points
  @param t NOT USED
  @param x input array
  @param u output array
 
  \f[ \nabla p = \mu \nabla^2 u \f]
  \f[ v^r = \|\vec v \| cos\theta(1-\frac{3rs}{2r} + \frac{rs^3}{2r^3})\f]
  \f[ v^\theta = -\|\vec v \| sin\theta(1-\frac{3rs}{4r} - \frac{rs^3}{4r^3})\f]
  \f[v^r(r>rs) = 0\f]
  \f[v^\theta(r>rs) = 0\f]
  \f[e^R = \frac{\vec x}{\|\vec x \|} ; \{ x_x = x - xs;x_y = y-ys; x_z = z-zs\} \f]
  \f[e^\theta = (e^R - \frac{\vec v}{\|\vec v \|}) -
                [(e^R - \frac{\vec v}{\|\vec v \|}) \cdot e^R]e^R\f]

  \f[u_z = v^r e^R_z + v^\theta e^\theta_z \f]

  \return status code
*/
int STflowSphere_Vz(int *iwork, double *rwork, int nPoints, double t, double *x, double *u)
{
/*
 Stokes Flow around Sphere of radius r^s and center x^s, y^s, z^s.
 Velocity is u,v,w and viscosity is mu.
 Array rwork contains all the above variables.
 Input	: t, x
		: x[npoints*3] = {x,y,z} for each point
		: rwork[] = {uu,v,w,rS,xS,yS,zS,mu}
 Output : Vz
*/
  int i;
  double vx=rwork[0], vy=rwork[1], vz=rwork[2];
  double rS=rwork[3], xS=rwork[4], yS=rwork[5], zS=rwork[6];
  double eR[3],eTHETA[3];
  double norm_v,theta,r;
  double vR,vTHETA;

  for (i=0; i<nPoints; i++)
    {
      coords(vx,vy,vz,
             xS,yS,zS,
             &x[i*3],
             &r,
             &theta,
             &norm_v,
             eR,
             eTHETA);
      vel(rS,norm_v,r,theta,&vR,&vTHETA);
      u[i] = vR * eR[2] + vTHETA* eTHETA[2];
    }

  return 0;
}

/********************************************************************/
/** @} */

void coords(double vx, double vy, double vz,
            double xS, double yS, double zS, 
            double* x, 
            double* r, 
            double* theta, 
            double* norm_v,
            double* eR, 
            double* eTHETA)
{
  double xBar[3],
    norm_xBar,
    xBar_dot_v,
    eTHETA_dot_eR,
    norm_eTHETA;
  
  *norm_v = sqrt(vx*vx+vy*vy+vz*vz);
  
  xBar[0]=x[0]-xS;
  xBar[1]=x[1]-yS;
  xBar[2]=x[2]-zS;
  
  norm_xBar = sqrt(xBar[0]*xBar[0] + xBar[1]*xBar[1] + xBar[2]*xBar[2]);
  
  *r = norm_xBar;

  if (norm_xBar != 0.0)
    {
      eR[0]=xBar[0]/norm_xBar;
      eR[1]=xBar[1]/norm_xBar;
      eR[2]=xBar[2]/norm_xBar;
    }
  else
    {
      eR[0]=0.0;
      eR[1]=0.0;
      eR[2]=0.0;
    }
  xBar_dot_v = xBar[0]*vx + xBar[1]*vy + xBar[2]*vz;
  
  if (norm_xBar != 0.0)
    *theta = acos(xBar_dot_v/(norm_xBar*(*norm_v)));
  else
    *theta = 0.0;

  eTHETA[0] = eR[0] - vx/(*norm_v);
  eTHETA[1] = eR[1] - vy/(*norm_v);
  eTHETA[2] = eR[2] - vz/(*norm_v);
  
  eTHETA_dot_eR = eTHETA[0]*eR[0] + eTHETA[1]*eR[1] + eTHETA[2]*eR[2];
  
  eTHETA[0] -= eTHETA_dot_eR*eR[0];
  eTHETA[1] -= eTHETA_dot_eR*eR[1];
  eTHETA[2] -= eTHETA_dot_eR*eR[2];
  
  norm_eTHETA = sqrt(eTHETA[0]*eTHETA[0] + eTHETA[1]*eTHETA[1] + eTHETA[2]*eTHETA[2]);
  
  if (norm_eTHETA > 0.0)
    {
      eTHETA[0]/= norm_eTHETA;
      eTHETA[1]/= norm_eTHETA;
      eTHETA[2]/= norm_eTHETA;
    }
  else
    {
      eTHETA[0]=0.0;
      eTHETA[1]=0.0;
      eTHETA[2]=0.0;
    }
}

void vel(double rS, double norm_v, double r, double theta, double* vR, double* vTHETA)
{
  if(r > rS)
    {
      *vR = norm_v * cos(theta) * ( 1.0 - (3.0*rS/(2.0*r)) + (pow(rS,3.0)/(2.0*pow(r,3.0))) );
      *vTHETA = -norm_v * sin(theta) * ( 1.0 - (3.0*rS/(4.0*r)) - (pow(rS,3.0)/(4.0*pow(r,3.0))) );
    }
  else
    {
      *vR = 0.0;
      *vTHETA = 0.0;
    }
}
/********************************************************************/
double uOfX_df(double nlC, double lu)
{
	return ( (2.0*nlC/(lu-nlC)+2.0) );
}
double uOfX_f(double a, double b, double nlC, double nlD, double x, double lu)
{
	 return ( (2.0 * (nlC * log(nlC-lu) - (nlC-lu)) - nlD) - b*x/a );
}
double f(double C, double b, double a, int q, int r)
{
	if (q==2 && r==1)
	{
		if (b != 0.0)
		{
			return ( C*tanh(b*C/a) - 1.0 );
		}
		else
		{
			printf("function f(q=2,r=1,b=0) returns -99999.99\n");
			return -99999.99;
		}
	}
	else if (q==1 && r==2)
	{
		return ( 2.0*C*(log(C-1.0) - log(C)) + 2.0 + b/a );
	}
	else 
	{
		printf("q,r not implemented: function f returns -99999.99\n");
		return -99999.99;
	}
}
double df(double C, double b, double a, int q, int r)
{
	if (q==2 && r==1)
	{
		if (b != 0.0)
		{
			return ( C*(1.0/pow(cosh(b*C/a),2))*(b/a) + tanh(b*C/a) );
		}
		else
		{
			printf("function df(q=2,r=1,b=0) returns -99999.99\n");
			return -99999.99;
		}
	}	
	else if (q==1 && r==2)
	{
		return ( 2.0*(log(C-1.0) - log(C)) + 2.0*C*(1.0/(C-1.0) - 1.0/C) ) ;

	}
	else 
	{
		printf("q,r not implemented: function df returns -99999.99\n");
		return -99999.99;
	}

}
/********************************************************************/


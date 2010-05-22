#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <strings.h>
#include <mach/mach_time.h>
/*#include <Accelerate/Accelerate.h>*/

//****************************************************
#pragma mark -
#pragma mark * local ( static ) function prototypes *
//----------------------------------------------------

static double CurrentTime(void);

//****************************************************
#pragma mark -
#pragma mark * exported function implementations *
//----------------------------------------------------

int main(int argc, char* argv[])
{
  const int vlen=atoi(argv[1]),stride=1;
  register int i,m,ntimes;
  register double a=atof(argv[2]),b=atof(argv[3]),c=atof(argv[4]),xi,yi;
  double *x,*y,*z,*u,*v,*w,start,stop,diff;
  ntimes =0;
  while (ntimes< 25)
    {
      ntimes+=1;
  if(argc < 5)
    {
      printf("usage: zerotest vlen da db dc");
      exit(1);
    }
  start = CurrentTime();
  x = malloc(sizeof(double)*vlen);
  y = malloc(sizeof(double)*vlen);
  z = malloc(sizeof(double)*vlen);
  u = malloc(sizeof(double)*vlen);
  v = malloc(sizeof(double)*vlen);
  w = malloc(sizeof(double)*vlen);
  stop = CurrentTime();
  printf("%10.2e s %10.2e seconds per double; malloc\n",(stop-start),(stop-start)/(3.0*vlen));
  
  start = CurrentTime();
  memset(x,0,vlen*sizeof(double));
  memset(y,0,vlen*sizeof(double));
  memset(z,0,vlen*sizeof(double));
  stop = CurrentTime();
  printf("%10.2e s %10.2e seconds per double; memset\n",(stop-start),(stop-start)/(3.0*vlen));

  start = CurrentTime();
  for (i=0;i<vlen;i++)
    {
      u[i]=c;
      v[i]=b;
      w[i]=a;
    }
  stop = CurrentTime();
  printf("%10.2e s %10.2e seconds per double; x[i]=c \n",(stop-start),(stop-start)/(3.0*vlen));
  
  start = CurrentTime();
  for (i=0;i<vlen;i++)
    {
      y[i]=x[i];
    }
  stop = CurrentTime();
  printf("%10.2e s %10.2e seconds per double; y[i]=x[i]\n",(stop-start),(stop-start)/(float)(vlen));
  
  start = CurrentTime();
  cblas_dcopy(vlen,x,stride,y,stride);
  stop = CurrentTime();
  printf("%10.2e s %10.2e seconds per double; dcopy y[i]=x[i]\n",(stop-start),(stop-start)/(float)(vlen));
  
  start = CurrentTime();
  for (i=0;i<vlen;i++)
    {
      x[i]*=a;
    }
  stop = CurrentTime();
  printf("%10.2e s %10.2f MFLOPS; y[i]=a*x[i]\n",(stop-start),((stop-start)/(float)(vlen))/1.0e6);
  
  start = CurrentTime();
  cblas_dscal(vlen,a,x,stride);
  stop = CurrentTime();
  printf("%10.2e s %10.2f MFLOPS; dscal y[i]=a*x[i]\n",(stop-start),((stop-start)/(float)(vlen))/1.0e6);
  
  start = CurrentTime();
  for(i=0;i<vlen;i++)
    y[i] += a*x[i];
  stop = CurrentTime();
  printf("%10.2e s %10.2f MFLOPS; y[i] += a*x[i]\n",(stop-start),(vlen*2/(stop-start))/1.0e6);

  start = CurrentTime();
  cblas_daxpy(vlen,a,x,stride,y,stride);
  stop = CurrentTime();
  printf("%10.2e s %10.2f MFLOPS; daxpy y[i] += a*x[i]\n",(stop-start),(vlen*2/(stop-start))/1.0e6);

  for(i=0;i<vlen;i++)
    y[i] += a*x[i];
  stop = CurrentTime();
  printf("%10.2e s %10.2f MFLOPS; y[i] += a*x[i]\n",(stop-start),(vlen*2/(stop-start))/1.0e6);

  start = CurrentTime();
  for(i=0;i<vlen;i++)
    y[i] = c + a*x[i];
  stop = CurrentTime();
  printf("%10.2e s %10.2f MFLOPS; y[i] = c + a*x[i]\n",(stop-start),(vlen*2/(stop-start))/1.0e6);

  start = CurrentTime();
  for(i=0;i<vlen;i++)
    y[i] = c + a*x[i] + b*y[i];
  stop = CurrentTime();
  printf("%10.2e s %10.2f MFLOPS; y[i] = c + a*x[i] + b*y[i]\n",(stop-start),(vlen*4/(stop-start))/1.0e6);

  start = CurrentTime();
  for(i=0;i<vlen;i++)
    y[i] = c + a*x[i] + b*y[i]+c*x[i]*y[i];
  stop = CurrentTime();
  printf("%10.2e s %10.2f MFLOPS y[i] = c + a*x[i] + b*y[i] + c*x[i]*y[i]\n",(stop-start),(vlen*7/(stop-start))/1.0e6);

  start = CurrentTime();
  for(i=0;i<vlen;i++)
    y[i] = c + a*x[i] + b*y[i]+c*x[i]*y[i] + a*x[i]*x[i];
  stop = CurrentTime();
  printf("%10.2e s %10.2f MFLOPS; y[i] = c + a*x[i] + b*y[i] + c*x[i]*y[i] + a*x[i]*x[i]\n",(stop-start),(vlen*10/(stop-start))/1.0e6);

  start = CurrentTime();
  for(i=0;i<vlen;i++)
    y[i] = c + a*x[i] + b*y[i]+c*x[i]*y[i] + a*x[i]*x[i] + b*y[i]*y[i];
  stop = CurrentTime();
  printf("%10.2e s %10.2f MFLOPS; y[i] = c + a*x[i] + b*y[i] + c*x[i]*y[i] + a*x[i]*x[i] + b*y[i]*y[i]\n",(stop-start),(vlen*13/(stop-start))/1.0e6);

  start = CurrentTime();
  for (i=0;i<vlen;i++)
    y[i] = c + a*x[i] + b*y[i] + c*x[i]*y[i] + a*x[i]*x[i] + b*y[i]*y[i] + c*x[i]*x[i]*x[i];
  stop = CurrentTime();
  printf("%10.2e s %10.2f MFLOPS; y[i] = c + a*x[i] + b*y[i] + c*x[i]*y[i] + a*x[i]*x[i] + b*y[i]*y[i] + c*x[i]*x[i]*x[i]\n",(stop-start),(vlen*17/(stop-start))/1.0e6);

  start = CurrentTime();
  for (i=0;i<vlen;i++)
    y[i] = c + a*x[i] + b*y[i] + c*x[i]*y[i] + a*x[i]*x[i] + b*y[i]*y[i] + c*x[i]*x[i]*x[i] + a*x[i]*y[i]*y[i] + b*x[i]*z[i]*x[i];
  stop = CurrentTime();
  printf("%10.2e s %10.2f MFLOPS; y[i] = c + a*x[i] + b*y[i] + c*x[i]*y[i] + a*x[i]*x[i] + b*y[i]*y[i] + c*x[i]*x[i]*x[i] + a*x[i]*y[i]*y[i] + b*x[i]*z[i]*x[i]\n",(stop-start),(vlen*25/(stop-start))/1.0e6);

  /* 3 vectors */

  start = CurrentTime();
  for(i=0;i<vlen;i++)
    y[i] = c + a*x[i] + b*z[i];
  stop = CurrentTime();
  printf("%10.2e s %10.2f MFLOPS; y[i] = c + a*x[i] + b*z[i]\n",(stop-start),(vlen*4/(stop-start))/1.0e6);

  start = CurrentTime();
  for(i=0;i<vlen;i++)
    y[i] = c + a*x[i] + b*z[i]+c*x[i]*z[i];
  stop = CurrentTime();
  printf("%10.2e s %10.2f MFLOPS y[i] = c + a*x[i] + b*z[i] + c*x[i]*z[i]\n",(stop-start),(vlen*7/(stop-start))/1.0e6);

  start = CurrentTime();
  for(i=0;i<vlen;i++)
    y[i] = c + a*x[i] + b*z[i]+c*x[i]*z[i] + a*x[i]*x[i];
  stop = CurrentTime();
  printf("%10.2e s %10.2f MFLOPS; y[i] = c + a*x[i] + b*z[i] + c*x[i]*z[i] + a*x[i]*x[i]\n",(stop-start),(vlen*10/(stop-start))/1.0e6);

  start = CurrentTime();
  for(i=0;i<vlen;i++)
    y[i] = c + a*x[i] + b*z[i]+c*x[i]*z[i] + a*x[i]*x[i] + b*z[i]*z[i];
  stop = CurrentTime();
  printf("%10.2e s %10.2f MFLOPS; y[i] = c + a*x[i] + b*z[i] + c*x[i]*z[i] + a*x[i]*x[i] + b*z[i]*z[i]\n",(stop-start),(vlen*13/(stop-start))/1.0e6);

  start = CurrentTime();
  for (i=0;i<vlen;i++)
    y[i] = c + a*x[i] + b*z[i] + c*x[i]*z[i] + a*x[i]*x[i] + b*z[i]*z[i] + c*x[i]*x[i]*x[i];
  stop = CurrentTime();
  printf("%10.2e s %10.2f MFLOPS; y[i] = c + a*x[i] + b*z[i] + c*x[i]*z[i] + a*x[i]*x[i] + b*z[i]*z[i] + c*x[i]*x[i]*x[i]\n",(stop-start),(vlen*17/(stop-start))/1.0e6);

  start = CurrentTime();
  for (i=0;i<vlen;i++)
    y[i] = c + a*x[i] + b*z[i] + c*x[i]*z[i] + a*x[i]*x[i] + b*z[i]*z[i] + c*x[i]*x[i]*x[i] + a*x[i]*z[i]*z[i] + b*x[i]*z[i]*x[i];
  stop = CurrentTime();
  printf("%10.2e s %10.2f MFLOPS; y[i] = c + a*x[i] + b*z[i] + c*x[i]*z[i] + a*x[i]*x[i] + b*z[i]*z[i] + c*x[i]*x[i]*x[i] + a*x[i]*z[i]*z[i] + b*x[i]*z[i]*x[i]\n",(stop-start),(vlen*25/(stop-start))/1.0e6);


  /* 3 vectors */
  start = CurrentTime();
  for(i=0;i<vlen;i++)
    y[i] = c + a*u[i] + b*z[i]+c*x[i]*z[i];
  stop = CurrentTime();
  printf("%10.2e s %10.2f MFLOPS y[i] = c + a*u[i] + b*z[i] + c*x[i]*z[i]\n",(stop-start),(vlen*7/(stop-start))/1.0e6);

  start = CurrentTime();
  for(i=0;i<vlen;i++)
    y[i] = c + a*u[i] + b*z[i]+c*x[i]*z[i] + a*u[i]*x[i];
  stop = CurrentTime();
  printf("%10.2e s %10.2f MFLOPS; y[i] = c + a*u[i] + b*z[i] + c*x[i]*z[i] + a*u[i]*x[i]\n",(stop-start),(vlen*10/(stop-start))/1.0e6);

  start = CurrentTime();
  for(i=0;i<vlen;i++)
    y[i] = c + a*u[i] + b*z[i]+c*x[i]*z[i] + a*u[i]*x[i] + b*z[i]*z[i];
  stop = CurrentTime();
  printf("%10.2e s %10.2f MFLOPS; y[i] = c + a*u[i] + b*z[i] + c*x[i]*z[i] + a*u[i]*x[i] + b*z[i]*z[i]\n",(stop-start),(vlen*13/(stop-start))/1.0e6);

  start = CurrentTime();
  for (i=0;i<vlen;i++)
    y[i] = c + a*u[i] + b*z[i] + c*x[i]*z[i] + a*u[i]*x[i] + b*z[i]*z[i] + c*u[i]*x[i]*u[i];
  stop = CurrentTime();
  printf("%10.2e s %10.2f MFLOPS; y[i] = c + a*u[i] + b*z[i] + c*x[i]*z[i] + a*u[i]*x[i] + b*z[i]*z[i] + c*u[i]*x[i]*u[i]\n",(stop-start),(vlen*17/(stop-start))/1.0e6);

  start = CurrentTime();
  for (i=0;i<vlen;i++)
    y[i] = c + a*u[i] + b*z[i] + c*x[i]*z[i] + a*u[i]*x[i] + b*z[i]*z[i] + c*u[i]*x[i]*u[i] + a*x[i]*z[i]*z[i] + b*u[i]*z[i]*x[i];
  stop = CurrentTime();
  printf("%10.2e s %10.2f MFLOPS; y[i] = c + a*u[i] + b*z[i] + c*x[i]*z[i] + a*u[i]*x[i] + b*z[i]*z[i] + c*u[i]*x[i]*u[i] + a*x[i]*z[i]*z[i] + b*u[i]*z[i]*x[i]\n",(stop-start),(vlen*25/(stop-start))/1.0e6);


  start = CurrentTime();
  free(x);
  free(y);
  free(z);
  free(u);
  free(v);
  free(w);
  stop = CurrentTime();
  printf("%10.2e s %10.2e seconds per double; free; \n",(stop-start),(stop-start)/(6.0*vlen));
}
}

//****************************************************
#pragma mark -
#pragma mark * local ( static ) function implementations *
//----------------------------------------------------


//
// Returns the current time in seconds
// 
//cek took from apple
static double CurrentTime(void)
{
    static double scale = 0.0;
	
    if (0.0 == scale) {
        mach_timebase_info_data_t info;
        mach_timebase_info(&info);
        scale = info.numer / info.denom * 1e-9;
    }
	
    return mach_absolute_time() * scale;
}	// CurrentTime

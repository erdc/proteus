#include "transportCoefficients.h"
#include <math.h>
/**
   \defgroup transportCoefficients transportCoefficients

   A library of coefficients for nonlinear PDE's.
   
   @{
*/

void linearADR_ConstantCoefficientsEvaluate(const int nPoints,
                                            const int nSpace,
                                            const double M,
                                            const double *A,
                                            const double *B,
                                            const double C,
                                            const double t,
                                            const double *x,
                                            const double *u,
                                            double *m,
                                            double *dm,
                                            double *f,
                                            double *df,
                                            double *a,
                                            double *da,
                                            double *phi,
                                            double *dphi,
                                            double *r,
                                            double *dr)
{
  int k,I,J;
  const int nSpace2=nSpace*nSpace;
  for (k=0;k<nPoints;k++)
    {
      m[k]=M*u[k];
      dm[k]=M;
      
      
      for (I=0;I<nSpace;I++)
        {
          f[k*nSpace+I]=B[I]*u[k];
          df[k*nSpace+I]=B[I];
          for (J=0;J<nSpace;J++)
            {
              a[k*nSpace2+I*nSpace+J]=A[I*nSpace+J];
              da[k*nSpace2+I*nSpace+J]=0.0;
            }
        }
      
      phi[k]=u[k];
      dphi[k]=1.0;
      
      r[k]=C*u[k];
      dr[k]=C;
    }
}

void nonlinearADR_pqrstEvaluate(const int nPoints,
                                const int nSpace,
                                const double M,
                                const double* A,
                                const double* B,
                                const double C,
                                const double p_pow,
                                const double q_pow,
                                const double r_pow,
                                const double s_pow,
                                const double t_pow,
                                const double t,
                                const double *x,
                                const double *u,
                                double *m,
                                double *dm,
                                double *f,
                                double *df,
                                double *a,
                                double *da,
                                double *phi,
                                double *dphi,
                                double *r,
                                double *dr)
{
  int k,I,J;
  const int nSpace2=nSpace*nSpace;
  double uPlus,
    tmp_f, 
    tmp_df,
    tmp_a,
    tmp_da;
  const double pM1_pow=p_pow-1.0,
    qM1_pow=q_pow-1.0,
    rM1_pow=r_pow-1.0,
    sM1_pow=s_pow-1.0,
    tM1_pow=t_pow-1.0;
  for (k=0;k<nPoints;k++)
    {
      uPlus = u[k] > 0.0 ? u[k] : 0.0;
      m[k]  = p_pow > 1.0 ?       M*pow(uPlus,p_pow)   : M*u[k];
      dm[k] = p_pow > 1.0 ? p_pow*M*pow(uPlus,pM1_pow) : M;

      tmp_f  = q_pow > 1.0 ?       pow(uPlus,q_pow)   : u[k];
      tmp_df = q_pow > 1.0 ? q_pow*pow(uPlus,qM1_pow) : 1.0;

      tmp_a  = t_pow > 0.0 ?       pow(uPlus,t_pow)   : 1.0;
      tmp_da = t_pow > 0.0 ? t_pow*pow(uPlus,tM1_pow) : 0.0;

      for (I=0;I<nSpace;I++)
        {
          f[k*nSpace+I]  = B[I]*tmp_f;
          df[k*nSpace+I] = B[I]*tmp_df;
          for (J=0;J<nSpace;J++)
            {
              a[k*nSpace2+I*nSpace+J]  = A[I*nSpace + J]*tmp_a;
              da[k*nSpace2+I*nSpace+J] = A[I*nSpace + J]*tmp_da;
            }
        }
      phi[k]  = r_pow > 1.0 ?       pow(uPlus,r_pow)   : u[k];
      dphi[k] = r_pow > 1.0 ? r_pow*pow(uPlus,rM1_pow) : 1.0;
      
      r[k]  =  s_pow > 1.0 ?       C*pow(uPlus,s_pow)   : C*u[k];
      dr[k] =  s_pow > 1.0 ? s_pow*C*pow(uPlus,sM1_pow) : C;
    }
}

void nonlinearADR_pqrstDualEvaluate(const int nPoints,
                                    const int nSpace,
                                    const double M,
                                    const double* A,
                                    const double* B,
                                    const double C,
                                    const double p1,
                                    const double q1,
                                    const double r1,
                                    const double s1,
                                    const double t1,
                                    const double p2,
                                    const double q2,
                                    const double r2,
                                    const double s2,
                                    const double t2,
                                    const double t,
                                    const double *x,
                                    const double *u,
                                    double *m,
                                    double *dm,
                                    double *f,
                                    double *df,
                                    double *a,
                                    double *da,
                                    double *phi,
                                    double *dphi,
                                    double *r,
                                    double *dr)
{
  int k,I,J;
  const int nSpace2=nSpace*nSpace;
  double max_1mu_0,atmp,datmp;
  const double p2M1=p2-1.0,
    q2M1=q2-1.0,
    r2M1=r2-1.0,
    s2M1=s2-1.0,
    t2M1=t2-1.0;

  nonlinearADR_pqrstEvaluate(nPoints,
                             nSpace,
                             M,
                             A,
                             B,
                             C,
                             p1,q1,r1,s1,t1,
                             t,
                             x,
                             u,
                             m,dm,
                             f,df,
                             a,da,
                             phi,dphi,
                             r,
                             dr);

  for (k=0; k < nPoints; k++)
    {
      max_1mu_0 = 1.0-u[k] > 0.0 ? 1.0-u[k] : 0.0;

      if (p2 > 1.0)
        {
          m[k] *= pow(max_1mu_0,p2);
          dm[k] *= pow(max_1mu_0,p2M1)*p2;
        }

      if (q2 > 1.0)
        {
          for (I=0; I < nSpace; I++)
            {
              f[k*nSpace+I] *= pow(max_1mu_0,q2);
              df[k*nSpace+I] *= pow(max_1mu_0,q2M1)*q2;
            }
        }

      if (t2 > 0.0)
        {
          atmp = pow(max_1mu_0,t2);
          datmp = pow(max_1mu_0,t2M1);
          for (I=0; I < nSpace; I++)
            {
              for (J=0; J < nSpace; J++)
                {
                  a[k*nSpace2+I*nSpace+J] *= atmp;
                  da[k*nSpace2+I*nSpace+J] *= datmp*t2;
                }
            }
        }

      if (r2 > 1.0)
        {
          phi[k] *= pow(max_1mu_0,r2);
          dphi[k] *= pow(max_1mu_0,r2M1)*r2;
        }

      if (s2 > 1.0)
        {
          r[k] *= pow(max_1mu_0,s2*r[k]);
          dr[k] *= pow(max_1mu_0,s2M1)*s2;
        }
    }
}

void unitSquareRotationEvaluate(const int nPoints,
                                const int nSpace,
                                const double *x,
                                const double *u,
                                double *m,
                                double *dm,
                                double *f,
                                double *df)
{
  double vx, vy;
  int k;
  for (k=0; k < nPoints; k++)
    {
      m[k] = u[k];
      dm[k] = 1.0;

      vx = 2.0*M_PI*(x[k*3+1] - 0.5);
      vy = 2.0*M_PI*(0.5     - x[k*3]);
      f[k*nSpace] =   vx*u[k];
      f[k*nSpace+1] = vy*u[k];
      df[k*nSpace] = vx;
      df[k*nSpace+1] = vy;
    }
}

void rotatingPulseVelEvaluate(const int nPoints,
                              const int nSpace,
                              const double self_a,
                              const double *x,
                              const double *u,
                              double *m,
                              double *dm,
                              double *f,
                              double *df,
                              double *a,
                              double *da,
                              double *phi,
                              double *dphi)
{
  double vx,vy;
  int k,I,J;
  const int nSpace2 = nSpace*nSpace;

  for (k=0; k < nPoints; k++)
    {
      m[k] = phi[k] = u[k];
      dm[k] = dphi[k] = 1.0;

      vx = -4.0*(x[k*3+1] - 0.5);
      vy = 4.0*(x[k*3] - 0.5);
      f[k*nSpace] = vx*u[k];
      f[k*nSpace+1] = vy*u[k];
      df[k*nSpace] = vx;
      df[k*nSpace+1] = vy;

      for (I=0; I < nSpace; I++)
        {
          a[k*nSpace2 + I*nSpace + I] = self_a;

          for (J=0; J < nSpace; J++)
            {
              da[k*nSpace2 + I*nSpace + J] = 0.0;
            }
        }
    }
}

void disRotatingPulseVelEvaluate(const int nPoints,
                                 const int nSpace,
                                 const double self_a,
                                 const double *x,
                                 const double *u,
                                 double *m,
                                 double *dm,
                                 double *f,
                                 double *df,
                                 double *a,
                                 double *da,
                                 double *phi,
                                 double *dphi)
{
  double X,Y,vx,vy;
  int k,I,J;
  const int nSpace2 = nSpace*nSpace;

  for (k=0; k < nPoints; k++)
    {
      m[k] = phi[k] = u[k];
      dm[k] = dphi[k] = 1.0;

      X = x[k*3+1] - 0.5;
      Y = x[k*3] - 0.5;
      vx = -4.0*X;
      vy = 4.0*Y;
      f[k*nSpace] = vx*u[k];
      f[k*nSpace+1] = vy*u[k];
      df[k*nSpace] = vx;
      df[k*nSpace+1] = vy;
      if (sqrt(X*X+Y*Y) < 0.25)
        {
          f[k*nSpace] *= 0.001;
          f[k*nSpace+1] *= 0.001;
          df[k*nSpace] *= 0.001;
          df[k*nSpace+1] *= 0.001;
        }

      for (I=0; I < nSpace; I++)
        {
          a[k*nSpace2 + I*nSpace + I] = self_a;

          for (J=0; J < nSpace; J++)
            {
              da[k*nSpace2 + I*nSpace + J] = 0.0;
            }
        }
    }
}

void disVelEvaluate(const int nPoints,
                    const int nSpace,
                    const double self_a,
                    const double *x,
                    const double *u,
                    double *m,
                    double *dm,
                    double *f,
                    double *df,
                    double *a,
                    double *da,
                    double *phi,
                    double *dphi)
{
  int k,I,J;
  const int nSpace2 = nSpace*nSpace;

  for (k=0; k < nPoints; k++)
    {
      m[k] = u[k];
      phi[k] = 0.0;
      dm[k] = dphi[k] = 1.0;

      f[k*nSpace] = u[k];
      f[k*nSpace+1] = 0.0;
      df[k*nSpace] = 1.0;
      df[k*nSpace+1] = 0.0;
      if (x[k*3+1] > 0.5)
        {
          f[k*nSpace] *= 0.25;
          df[k*nSpace] *= 0.25;
        }

      for (I=0; I < nSpace; I++)
        {
          a[k*nSpace2 + I*nSpace + I] = self_a;

          for (J=0; J < nSpace; J++)
            {
              da[k*nSpace2 + I*nSpace + J] = 0.0;
            }
        }
    }
}

void burgersDiagonalVelEvaluate(const int nPoints,
                                const int nSpace,
                                const double self_a,
                                const double *u,
                                double *m,
                                double *dm,
                                double *f,
                                double *df,
                                double *a,
                                double *da,
                                double *phi,
                                double *dphi)
{
  double maxu0;
  int k,I,J;
  const int nSpace2 = nSpace*nSpace;

  for (k=0; k < nPoints; k++)
    {
      m[k] = phi[k] = u[k];
      dm[k] = dphi[k] = 1.0;
      maxu0 = u[k] > 0.0 ? u[k] : 0.0;

      for (I=0; I < nSpace; I++)
        {
          a[k*nSpace2 + I*nSpace + I] = self_a;
          f[k*nSpace+I] = pow(maxu0,2.0) * 0.5;
          df[k*nSpace+I] = maxu0;

          for (J=0; J < nSpace; J++)
            {
              da[k*nSpace2 + I*nSpace + J] = 0.0;
            }
        }
    }
}

void twophasePotentialFlowEvaluate(int nPoints,
                                   int nSpace,
                                   double *M,
                                   double *A,
                                   double *B,
                                   double *Bcon,
                                   double *C,
                                   double t,
                                   double *x,
                                   double *u,
                                   double *m,
                                   double *dm,
                                   double *f,
                                   double *df,
                                   double *a,
                                   double *da,
                                   double *phi,
                                   double *dphi,
                                   double *r,
                                   double *dr)
{
  int k,I,J;
  const int nSpace2=nSpace*nSpace;
  for (k=0;k<nPoints;k++)
    {
      m[k]=M[k]*u[k];
      dm[k]=M[k];
      for (I=0;I<nSpace;I++)
        {
          f[k*nSpace+I]=B[k*nSpace+I]*u[k] + Bcon[k*nSpace+I];
          df[k*nSpace+I]=B[k*nSpace+I];
          for (J=0;J<nSpace;J++)
            {
              a[k*nSpace2+I*nSpace+J]=A[k*nSpace2 + I*nSpace+J];
              da[k*nSpace2+I*nSpace+J]=0.0;
            }
        }
      
      phi[k]=u[k];
      dphi[k]=1.0;
      
      r[k]=C[k]*u[k];
      dr[k]=C[k];
    }
}

void twophasePotentialFlowUpdateFreeSurface(int nPoints,
                                            int nSpace,
                                            double eps,
                                            double* u_levelSet,
                                            double M1, double M2, double *M,
                                            double* A1, double* A2, double *A,
                                            double* B1, double* B2, double *B,
                                            double* Bcon1, double* Bcon2, double *Bcon,
                                            double C1, double C2, double *C)
{
  int k,I,J;
  const int nSpace2=nSpace*nSpace;
  double H,oneMinusH;
  for (k=0;k<nPoints;k++)
    {
      if (u_levelSet[k] > eps)
        {
          M[k] = M1;
          for (I=0;I<nSpace;I++)
            {
              B[k*nSpace+I]=B1[I];
              Bcon[k*nSpace+I]=Bcon1[I];
              for (J=0;J<nSpace;J++)
                {
                  A[k*nSpace2+I*nSpace+J]=A1[I*nSpace+J];
                }
            }
          C[k] = C1;
        }
      else if (u_levelSet[k] < -eps)
        {
          M[k] = M2;
          for (I=0;I<nSpace;I++)
            {
              B[k*nSpace+I]=B2[I];
              Bcon[k*nSpace+I]=Bcon2[I];
              for (J=0;J<nSpace;J++)
                {
                  A[k*nSpace2+I*nSpace+J]=A2[I*nSpace+J];
                }
            }
          C[k] = C2;
        }
      else
        {
          H = 0.5*(1.0 + u_levelSet[k]/eps + sin((M_PI*u_levelSet[k])/eps)/M_PI);
          oneMinusH=1.0-H;
          M[k] = oneMinusH*M2 + H*M1;
          for (I=0;I<nSpace;I++)
            {
              B[k*nSpace+I]=oneMinusH*B2[I] + H*B1[I];
              Bcon[k*nSpace+I]=oneMinusH*Bcon2[I] + H*Bcon1[I];
              for (J=0;J<nSpace;J++)
                {
                  A[k*nSpace2+I*nSpace+J]=oneMinusH*A2[I*nSpace+J] + H*A1[I*nSpace+J];
                }
            }
          C[k]=oneMinusH*C2+H*C1;
        }
    }
}

void twophaseLevelSetCoefficientsUpdateVelocity(int nPoints,
                                                int nSpace,
                                                double v_scale,
                                                double* vIn,
                                                double* vOut)
{
  int i,I;
  for  (i=0;i<nPoints;i++)
    for (I=0;I<nSpace;I++)
      vOut[i*nSpace+I]=v_scale*vIn[i*nSpace+I];
}

void twophaseLevelSetCoefficientsEvaluate(int nPoints,
                                          int nSpace,
                                          double* B,
                                          double  t,
                                          double* x,
                                          double* u,
                                          double* grad_u,
                                          double* m, double* dm,
                                          double* h, double* dh,
                                          double* rh)
{
  int i,I;
  for (i=0;i<nPoints;i++)
    {
      rh[i]=0.0;
      h[i]=0.0;
      m[i]=u[i];
      dm[i]=1.0;
      for (I=0;I<nSpace;I++)
        {
          h[i] += B[i*nSpace+I]*grad_u[i*nSpace+I];
          dh[i*nSpace+I]=B[i*nSpace+I];
        }
    }
}

void twophaseLevelSetCoefficientsEvaluateCI(int nPoints,
                                            int nSpace,
                                            double* B,
                                            double  t,
                                            double* x,
                                            double* u,
                                            double* m, double* dm,
                                            double* f, double* df,
                                            double* a, double* da,
                                            double* phi, double* dphi,
                                            double* r, double* dr)
{
  int i,I;
  for (i=0;i<nPoints;i++)
    {
      m[i]=u[i];
      dm[i]=1.0;
      for (I=0;I<nSpace;I++)
        {
          f[i*nSpace+I] = B[i*nSpace+I]*u[i];
          df[i*nSpace+I]=B[i*nSpace+I];
        }
    }
}

void twophaseSignedDistanceCoefficientsUpdateSignFunction(int nPoints,
                                                          double eps,
                                                          double* u_levelSet,
                                                          double* S)
{
  int k;
  double H;
  for (k=0;k<nPoints;k++)
    {
      if (u_levelSet[k] > eps)
	S[k]= 1.0;
      else if (u_levelSet[k] < -eps)
        S[k] = -1.0;
      else
        {
          H = 0.5*(1.0 + u_levelSet[k]/eps + sin((M_PI*u_levelSet[k])/eps)/M_PI);
          S[k]=2.0*H - 1.0;
        }
    }
/*   for  (k=0;k<nPoints;k++) */
/*     S[k] = u_levelSet[k]/sqrt(u_levelSet[k]*u_levelSet[k] +  eps*eps); */
/*     if(u_levelSet[k] > 0.0) */
/*       S[k] = 1.0; */
/*     else */
/*       S[k] = -1.0; */
}

void twophaseSignedDistanceCoefficientsEvaluate(int nPoints,
                                                int nSpace,
                                                double* S,
                                                double* u,
                                                double* grad_u,
                                                double* m, double* dm,
                                                double* h, double* dh,
                                                double* rh)
{
  int i,I;
  for (i=0;i<nPoints;i++)
    {
      m[i]=u[i];
      dm[i]=1.0;
      rh[i]=-1.0;
      h[i]=0.0;
      for (I=0;I<nSpace;I++)
        h[i] += grad_u[i*nSpace+I]*grad_u[i*nSpace+I];
      h[i] = sqrt(h[i]);
      for (I=0;I<nSpace;I++)
	if(h[i]>= fabs(S[i]*grad_u[i*nSpace+I])*1.0e-8)
	  dh[i*nSpace+I] = (S[i]*grad_u[i*nSpace+I])/h[i];
	else
          dh[i*nSpace+I] = (S[i]*grad_u[i*nSpace+I])/1.0e-8;
      h[i]+=rh[i];
      h[i]*=S[i];
      rh[i]*=S[i];
    }
}
/* /\** Coefficients for the mass conservative head-based Richards' equation using Mualem-Van Genuchten. */
/*  *\/ */
/* void conservativeHeadRichardsMualemVanGenuchtenEvaluate(const int nPoints, */
/*                                                         const int nSpace, */
/*                                                         const double rho, */
/*                                                         const double* gravity, */
/*                                                         const double* alpha, */
/*                                                         const double* n, */
/*                                                         const double* m, */
/*                                                         const double* thetaS, */
/*                                                         const double* thetaR, */
/*                                                         const double* thetaSR, */
/*                                                         const double* KWs, */
/*                                                         double *u, */
/*                                                         double *mass, */
/*                                                         double *dmass, */
/*                                                         double *f, */
/*                                                         double *df, */
/*                                                         double *a, */
/*                                                         double *da, */
/*                                                         double *phi, */
/*                                                         double *dphi) */
/* { */
/*   int k,I,J; */
/*   const int nSpace2=nSpace*nSpace; */
/*   register double psiC,alphaPsiC,alphaPsiC_n,alphaPsiC_nM1,onePlus_alphaPsiC_n,sBar,sBarByOnePlus_alphaPsiC_n,sqrt_sBar,sqrt_1minusSbar,thetaW,DsBar_DpC,DthetaW_DpC,vBar,uBar,krW,krN, */
/*     alphaPsiC_nM2,sBarBy_onePlus_alphaPsiC_n_2,DDsBar_DDpC,DDthetaW_DDpC,DkrW_DpC,rho2=rho*rho; */
/*   double KW[9],KN[9],DKW_DpC[9]; */
/*   for (k=0;k<nPoints;k++) */
/*     { */
/*       psiC = -u[k]; */
/*       if (psiC > 0.0) */
/*         { */
/*           alphaPsiC = alpha[k]*psiC; */
/*           alphaPsiC_n = pow(alphaPsiC,n[k]); */
/*           alphaPsiC_nM1 = alphaPsiC_n/alphaPsiC; */
/*           alphaPsiC_nM2 =   alphaPsiC_nM1/alphaPsiC;                 */
/*           onePlus_alphaPsiC_n = 1.0 + alphaPsiC_n; */
/*           sBar = pow(onePlus_alphaPsiC_n,-m[k]); */
/*           sBarByOnePlus_alphaPsiC_n = sBar/onePlus_alphaPsiC_n; */
/*           sqrt_sBar = sqrt(sBar); */
/*           sqrt_1minusSbar = sqrt(1.0 - sBar); */
/*           thetaW = thetaSR[k]*sBar + thetaR[k]; */
/*           DsBar_DpC = -alpha[k]*(n[k]-1.0)*alphaPsiC_nM1  */
/*             *sBarByOnePlus_alphaPsiC_n; */
/*           DthetaW_DpC = thetaSR[k] * DsBar_DpC;  */
/*           vBar = 1.0-alphaPsiC_nM1*sBar; */
/*           uBar = alphaPsiC_nM1*sBar; */
/*           krW = sqrt_sBar*vBar*vBar; */
/*           DkrW_DpC = (0.5/sqrt_sBar)*DsBar_DpC*vBar*vBar */
/*             - */
/*             2.0*sqrt_sBar*vBar* */
/*             (alpha[k]*(n[k]-1.0)*alphaPsiC_nM2*sBar */
/*              + alphaPsiC_nM1 * DsBar_DpC); */
/*           for (I=0;I<nSpace;I++) */
/*             for(J=0;J<nSpace;J++) */
/*               { */
/*                 KW[I*nSpace + J] = KWs[k*nSpace2 + I*nSpace + J]*krW;  */
/*                 DKW_DpC[I*nSpace + J] = KWs[k*nSpace2 + I*nSpace + J]*DkrW_DpC;  */
/*               } */
/*         } */
/*       else */
/*         { */
/*           sBar = 1.0; */
/*           thetaW = thetaS[k]; */
/*           DsBar_DpC = 0.0; */
/*           DthetaW_DpC = 0.0; */
/*           krW = 1.0; */
/*           for (I=0;I<nSpace;I++) */
/*             for(J=0;J<nSpace;J++) */
/*               { */
/*                 KW[I*nSpace + J] = KWs[k*nSpace2 + I*nSpace + J];  */
/*                 DKW_DpC[I*nSpace + J] = 0.0;  */
/*               } */
/*         } */
/*       mass[k] = rho*thetaW; */
/*       dmass[k] = -rho*DthetaW_DpC; */
/*       for (I=0;I<nSpace;I++) */
/*         { */
/*           f[k*nSpace+I] = 0.0; */
/*           df[k*nSpace+I] = 0.0; */
/*           for (J=0;J<nSpace;J++) */
/*             { */
/*               f[k*nSpace+I]  += rho2*KW[I*nSpace + J]*gravity[J]; */
/*               df[k*nSpace+I] -= rho2*DKW_DpC[I*nSpace + J]*gravity[J]; */
/*               a[k*nSpace2+I*nSpace+J]  = rho*KW[I*nSpace + J]; */
/*               da[k*nSpace2+I*nSpace+J] = -rho*DKW_DpC[I*nSpace + J]; */
/*             } */
/*         } */
/*     } */
/* } */

/** Coefficients for the mass conservative head-based Richards' equation using Mualem-Van Genuchten.
 */
void conservativeHeadRichardsMualemVanGenuchtenHomEvaluate(const int nPoints,
                                                           const int nSpace,
                                                           const double rho,
                                                           const double* gravity,
                                                           const double alpha,
                                                           const double n,
                                                           const double m,
                                                           const double thetaR,
                                                           const double thetaSR,
                                                           const double KWs,
                                                           double *u,
                                                           double *mass,
                                                           double *dmass,
                                                           double *f,
                                                           double *df,
                                                           double *a,
                                                           double *da)
{
  int k,I;
  const int nSpace2=nSpace*nSpace;
  register double psiC,
    pcBar,pcBar_n,pcBar_nM1,pcBar_nM2,
    onePlus_pcBar_n,
    sBar,sqrt_sBar,DsBar_DpsiC,
    thetaW,DthetaW_DpsiC,
    vBar,vBar2,DvBar_DpsiC,
    KW,DKW_DpsiC,
    rho2=rho*rho,
    thetaS=thetaR+thetaSR;
  for (k=0;k<nPoints;k++)
    {
      psiC = -u[k];
      if (psiC > 0.0)
        {
          pcBar = alpha*psiC;
          pcBar_nM2 = pow(pcBar,n-2);
          pcBar_nM1 = pcBar_nM2*pcBar;
          pcBar_n   = pcBar_nM1*pcBar;
          onePlus_pcBar_n = 1.0 + pcBar_n;
          
          sBar = pow(onePlus_pcBar_n,-m);
          /* using -mn = 1-n */
          DsBar_DpsiC = alpha*(1.0-n)*(sBar/onePlus_pcBar_n)*pcBar_nM1;
          
          vBar = 1.0-pcBar_nM1*sBar;
          vBar2 = vBar*vBar;
          DvBar_DpsiC = -alpha*(n-1.0)*pcBar_nM2*sBar - pcBar_nM1*DsBar_DpsiC;
          
          thetaW = thetaSR*sBar + thetaR;
          DthetaW_DpsiC = thetaSR * DsBar_DpsiC; 
          
          sqrt_sBar = sqrt(sBar);
          KW= KWs*sqrt_sBar*vBar2;
          DKW_DpsiC= KWs*
            ((0.5/sqrt_sBar)*DsBar_DpsiC*vBar2
             +
             2.0*sqrt_sBar*vBar*DvBar_DpsiC);
        }
      else
        {
          thetaW        = thetaS;
          DthetaW_DpsiC = 0.0;
          KW            = KWs; 
          DKW_DpsiC     = 0.0;
        }
      mass[k] = rho*thetaW;
      dmass[k] = -rho*DthetaW_DpsiC;
      for (I=0;I<nSpace;I++)
        {
          f[k*nSpace+I]  = rho2*KW*gravity[I];
          df[k*nSpace+I] = -rho2*DKW_DpsiC*gravity[I];
          a[k*nSpace2+I*nSpace+I]  = rho*KW;
          da[k*nSpace2+I*nSpace+I] = -rho*DKW_DpsiC;
        }
    }
}

/** Coefficients for the mass conservative saturation-based Richards' equation using Mualem-Van Genuchten.
 */
void conservativeSatRichardsMualemVanGenuchtenHomEvaluate(const int nPoints,
                                                          const int nSpace,
                                                          const double rho,
                                                          const double* gravity,
                                                          const double alpha,
                                                          const double n,
                                                          const double m,
                                                          const double thetaR,
                                                          const double thetaSR,
                                                          const double KWs,
                                                          double *u,
                                                          double *mass,
                                                          double *dmass,
                                                          double *f,
                                                          double *df,
                                                          double *a,
                                                          double *da,
                                                          double *phi,
                                                          double *dphi)
{
  int k,I;
  const int nSpace2=nSpace*nSpace;
  register double psiC,
    pcBar,pcBar_n,pcBar_nM1,pcBar_nM2,
    onePlus_pcBar_n,
    sBar,sqrt_sBar,DsBar_DpsiC,
    thetaW,DthetaW_DpsiC,
    vBar,vBar2,DvBar_DpsiC,
    KW,DKW_DpsiC,
    rho2=rho*rho,
    thetaS=thetaR+thetaSR,
    sR = thetaR/thetaS;
  for (k=0;k<nPoints;k++)
    {
      sBar = u[k];
      thetaW = thetaS*u[k];
      if (u[k] < 1.0)
        {
          /* piggy back  on head  based formulas */
          psiC = pow(pow(sBar,-1.0/m)-1.0,1.0/n)/alpha;
          
          pcBar = alpha*psiC;
          pcBar_nM2 = pow(pcBar,n-2);
          pcBar_nM1 = pcBar_nM2*pcBar;
          pcBar_n   = pcBar_nM1*pcBar;
          onePlus_pcBar_n = 1.0 + pcBar_n;
          
          /* using -mn = 1-n */
          DsBar_DpsiC = alpha*(1.0-n)*(sBar/onePlus_pcBar_n)*pcBar_nM1;
          
          vBar = 1.0-pcBar_nM1*sBar;
          vBar2 = vBar*vBar;
          DvBar_DpsiC = -alpha*(n-1.0)*pcBar_nM2*sBar - pcBar_nM1*DsBar_DpsiC;
          
          sqrt_sBar = sqrt(sBar);
          KW= KWs*sqrt_sBar*vBar2;
          DKW_DpsiC= KWs*
            ((0.5/sqrt_sBar)*DsBar_DpsiC*vBar2
             +
             2.0*sqrt_sBar*vBar*DvBar_DpsiC);
        }
      else
        {
          KW            = KWs; 
          DKW_DpsiC     = 0.0;
        }
      mass[k] = rho*thetaW;
      dmass[k] = rho*thetaS;
      phi[k] = -psiC;
      dphi[k] = -1.0/DsBar_DpsiC;
      for (I=0;I<nSpace;I++)
        {
          f[k*nSpace+I]  = rho2*KW*gravity[I];
          df[k*nSpace+I] = rho2*DKW_DpsiC*gravity[I]/DsBar_DpsiC;
          a[k*nSpace2+I*nSpace+I]  = rho*KW;
          da[k*nSpace2+I*nSpace+I] = rho*DKW_DpsiC/DsBar_DpsiC;
        }
    }
}

/* /\** Coefficients for the mass conservative Saturation equation for the fractional flow form of two-phase flow equations using Mualem-Van Genuchten. */
/*  *\/ */
/* void conservativeTwophaseSaturationMualemVanGenuchtenHomEvaluate(const int nPoints, */
/*                                                                  const int nSpace, */
/*                                                                  const double rho, */
/*                                                                  const double* gravity, */
/*                                                                  const double alpha, */
/*                                                                  const double n, */
/*                                                                  const double m, */
/*                                                                  const double thetaR, */
/*                                                                  const double thetaSR, */
/*                                                                  const double KWs, */
/*                                                                  const double viscosityRatio, */
/*                                                                  const double densityRatio, */
/*                                                                  const double* v, */
/*                                                                  double *u, */
/*                                                                  double *mass, */
/*                                                                  double *dmass, */
/*                                                                  double *f, */
/*                                                                  double *df, */
/*                                                                  double *a, */
/*                                                                  double *da) */
/* { */
/*   int k,I; */
/*   const int nSpace2=nSpace*nSpace; */
/*   register double psiC, */
/*     pcBar,pcBar_n,pcBar_nM1,pcBar_nM2, */
/*     onePlus_pcBar_n, */
/*     sBar,sqrt_sBar,DsBar_DpsiC, */
/*     thetaW,DthetaW_DpsiC, */
/*     vBar,vBar2,DvBar_DpsiC, */
/*     KW,DKW_DpsiC, */
/*     rho2=rho*rho, */
/*     thetaS=thetaR+thetaSR; */
/*   for (k=0;k<nPoints;k++) */
/*     { */
/*       sBar = fmin(fmax((u[k] - sIR)/sMir,0.0),1.0); */
/*       DsBar_Du = 1.0/sMIR; */

/*       onePlus_alphaPsiC_n = pow(sBar,1.0/-m); */
/*       alphaPsiC_n = onePlus_alphaPsiC_n - 1.0; */
/*       alphaPsiC = pow(alphaPsiC_n,1.0/n); */
      
/*       psiC = alphaPsiC/alpha; */
  
/*       alphaPsiC_nM1 = alphaPsiC_n/alphaPsiC; */
/*       sBarByOnePlus_alphaPsiC_n = sBar/onePlus_alphaPsiC_n; */
/*       DsBar_DpC = -alpha[i]*(n[i]-1.0)*alphaPsiC_nM1  */
/*         *sBarByOnePlus_alphaPsiC_n; */
      
/*       if(psiC<=0.0)  */
/*         { */
/*           DsBar_DpC = 0.0; */
/*         } */
/* } */


/* inline void VanGenuchten2p::calculateDerivatives() */
/* { */
/*   alphaPsiC_nM2 =   alphaPsiC_nM1/alphaPsiC;       */
  
/*   sBarBy_onePlus_alphaPsiC_n_2 = sBarByOnePlus_alphaPsiC_n */
/*     /onePlus_alphaPsiC_n; */
/*   DDsBar_DDpC =  alpha[i]*alpha[i]*(n[i]-1) */
/*     *((2*n[i]-1)*alphaPsiC_nM1*alphaPsiC_nM1 */
/*       *sBarBy_onePlus_alphaPsiC_n_2 */
/*       - */
/*       (n[i]-1)*alphaPsiC_nM2 */
/*       *sBarByOnePlus_alphaPsiC_n); */

/*   if (psiC <= 0.0) */
/*     { */
/*       DDsBar_DDpC = 0.0; */
/*     } */
/* } */
/*       psiC = -u[k]; */
/*       if (psiC > 0.0) */
/*         { */
/*           pcBar = alpha*psiC; */
/*           pcBar_nM2 = pow(pcBar,n-2); */
/*           pcBar_nM1 = pcBar_nM2*pcBar; */
/*           pcBar_n   = pcBar_nM1*pcBar; */
/*           onePlus_pcBar_n = 1.0 + pcBar_n; */
          
/*           sBar = pow(onePlus_pcBar_n,-m); */
/*           /\* using -mn = 1-n *\/ */
/*           DsBar_DpsiC = alpha*(1.0-n)*(sBar/onePlus_pcBar_n)*pcBar_nM1; */
          
/*           vBar = 1.0-pcBar_nM1*sBar; */
/*           vBar2 = vBar*vBar; */
/*           DvBar_DpsiC = -alpha*(n-1.0)*pcBar_nM2*sBar - pcBar_nM1*DsBar_DpsiC; */
          
/*           thetaW = thetaSR*sBar + thetaR; */
/*           DthetaW_DpsiC = thetaSR * DsBar_DpsiC;  */
          
/*           sqrt_sBar = sqrt(sBar); */
/*           krw = sqrt_sBar*vBar2; */
/*           Dkrw_DpsiC= KWs* */
/*             ((0.5/sqrt_sBar)*DsBar_DpsiC*vBar2 */
/*              + */
/*              2.0*sqrt_sBar*vBar*DvBar_DpsiC); */
/*         } */
/*       else */
/*         { */
/*           thetaW        = thetaS; */
/*           DthetaW_DpsiC = 0.0; */
/*           KW            = KWs;  */
/*           DKW_DpsiC     = 0.0; */
/*         } */
/*       mass[k] = rho*thetaW; */
/*       dmass[k] = -rho*DthetaW_DpsiC; */
/*       for (I=0;I<nSpace;I++) */
/*         { */
/*           f[k*nSpace+I]  = rho2*KW*gravity[I]; */
/*           df[k*nSpace+I] = -rho2*DKW_DpsiC*gravity[I]; */
/*           a[k*nSpace2+I*nSpace+I]  = rho*KW; */
/*           da[k*nSpace2+I*nSpace+I] = -rho*DKW_DpsiC; */
/*         } */
/*     } */
/* } */

/** @} */

void valFromDOF(int eN, int k, int j, int t,double* u,double* dof,int* l2g,double* v)
{
  u[eN*nQuadraturePoints_element*nComponents+
    k*nComponents+
    t] 
    += 
    dof[l2g[eN*nDOF_trial_element+
            j]*nComponents+
        t]
    *
    v[eN*nQuadraturePoints_element*nDOF_trial_element+
      k*nDOF_trial_element+
      j];
}

void gradFromDOF(int eN, int k, int j, int t, int I,double* grad_u,double* dof,int* l2g,double* grad_v)
{
  grad_u[eN*nQuadraturePoints_element*nComponents*nSpace+
         k*nComponents*nSpace+
         t*nSpace+
         I]
    +=
    dof[l2g[eN*nDOF_trial_element+
            j]*nComponents+
        t]
    *
    grad_v[eN*nQuadraturePoints_element*nDOF_trial_element*nSpace+
           k*nDOF_trial_element*nSpace+
           j*nSpace+
           I];
}

void valFromDOF_trace_ext(int ebNE,
                          int eN,
                          int ebN_local,
                          int k,
                          int j,
                          int t,
                          const int* exteriorElementBoundariesArray,
                          const int* elementBoundaryElementsArray,
                          const int* elementBoundaryLocalElementBoundariesArray,
                          double* psi,
                          double* vArray)
{
  vArray[ebNE*nElementBoundaryQuadraturePoints_elementBoundary*nDOF_element+
         k*nDOF_element+
         j]
    =
    psi[ebN_local*nElementBoundaryQuadraturePoints_elementBoundary*nDOF_element + 
        k*nDOF_element + 
        j];
}

void gradFromDOF_trace_ext(int ebNE,
                           int eN,
                           int k,
                           int t,
                           int I,
                           int* l2g,
                           double* dof,
                           double* grad_v,
                           double* grad_u)
{
  grad_u[ebNE*nQuadraturePoints_elementBoundary*nComponents*nSpace+
         k*nComponents*nSpace+
         t*nSpace+
         I] 
    += 
    dof[l2g[eN*nDOF_trial_element+j]*nComponents+
        t]
    *
    grad_v[ebNE*nQuadraturePoints_elementBoundary*nDOF_trial_element*nSpace+
           k*nDOF_trial_element*nSpace+
           j*nSpace+
           I];
}

void evaluateCoefficients(int K,
                          const double eps_rho,
                          const double eps_mu,
                          const double sigma,
                          const double rho_0,
                          const double nu_0,
                          const double rho_1,
                          const double nu_1,
                          const double* g,
                          const double* phi,
                          const double* n,
                          const double* kappa,
                          const double *p,
                          const double *grad_p,
                          const double *u,
                          const double *v,
                          const double *w,
                          double *mom_u_acc,
                          double *dmom_u_acc_u,
                          double *mom_v_acc,
                          double *dmom_v_acc_v,
                          double *mom_w_acc,
                          double *dmom_w_acc_w,
                          double *mass_adv,
                          double *dmass_adv_u,
                          double *dmass_adv_v,
                          double *dmass_adv_w,
                          double *mom_u_adv,
                          double *dmom_u_adv_u,
                          double *dmom_u_adv_v,
                          double *dmom_u_adv_w,
                          double *mom_v_adv,
                          double *dmom_v_adv_u,
                          double *dmom_v_adv_v,
                          double *dmom_v_adv_w,
                          double *mom_w_adv,
                          double *dmom_w_adv_u,
                          double *dmom_w_adv_v,
                          double *dmom_w_adv_w,
                          double *mom_u_diff_ten,
                          double *mom_v_diff_ten,
                          double *mom_w_diff_ten,
                          double *mom_uv_diff_ten,
                          double *mom_uw_diff_ten,
                          double *mom_vu_diff_ten,
                          double *mom_vw_diff_ten,
                          double *mom_wu_diff_ten,
                          double *mom_wv_diff_ten,
                          double *mom_u_source,
                          double *mom_v_source,
                          double *mom_w_source,
                          double *mom_u_ham,
                          double *dmom_u_ham_grad_p,
                          double *mom_v_ham,
                          double *dmom_v_ham_grad_p,
                          double *mom_w_ham,
                          double *dmom_w_ham_grad_p)
{
  double rho,nu,mu,H_rho,d_rho,H_mu,d_mu,norm_n;
  H_rho = smoothedHeaviside(eps_rho,phi[K]);
  d_rho = smoothedDirac(eps_rho,phi[K]);
  H_mu = smoothedHeaviside(eps_mu,phi[K]);
  d_mu = smoothedDirac(eps_mu,phi[K]);
  
  rho = rho_0*(1.0-H_rho)+rho_1*H_rho;
  nu  = nu_0*(1.0-H_mu)+nu_1*H_mu;
  mu  = rho_0*nu_0*(1.0-H_mu)+rho_1*nu_1*H_mu;
  
  //u momentum accumulation
  mom_u_acc[K]=u[K];
  dmom_u_acc_u[K]=1.0;
  
  //v momentum accumulation
  mom_v_acc[K]=v[K];
  dmom_v_acc_v[K]=1.0;
  
  //w momentum accumulation
  mom_w_acc[K]=w[K];
  dmom_w_acc_w[K]=1.0;
  
  
  //mass advective flux
  mass_adv[K*3+0]=u[K];
  mass_adv[K*3+1]=v[K];
  mass_adv[K*3+2]=w[K];
  
  dmass_adv_u[K*3+0]=1.0;
  dmass_adv_v[K*3+1]=1.0;
  dmass_adv_w[K*3+2]=1.0;
  
  //u momentum advective flux
  mom_u_adv[K*3+0]=u[K]*u[K];
  mom_u_adv[K*3+1]=u[K]*v[K];
  mom_u_adv[K*3+2]=u[K]*w[K];
  
  dmom_u_adv_u[K*3+0]=2.0*u[K];
  dmom_u_adv_u[K*3+1]=v[K];
  dmom_u_adv_u[K*3+2]=w[K];
  
  dmom_u_adv_v[K*3+1]=u[K];
  
  dmom_u_adv_w[K*3+2]=u[K];
  
  //v momentum advective_flux
  mom_v_adv[K*3+0]=v[K]*u[K];
  mom_v_adv[K*3+1]=v[K]*v[K];
  mom_v_adv[K*3+2]=v[K]*w[K];
  
  dmom_v_adv_u[K*3+0]=v[K];
  
  dmom_v_adv_w[K*3+2]=v[K];
  
  dmom_v_adv_v[K*3+0]=u[K];
  dmom_v_adv_v[K*3+1]=2.0*v[K];
  dmom_v_adv_v[K*3+2]=w[K];
  
  //w momentum advective_flux
  mom_w_adv[K*3+0]=w[K]*u[K];
  mom_w_adv[K*3+1]=w[K]*v[K];
  mom_w_adv[K*3+2]=w[K]*w[K];
  
  dmom_w_adv_u[K*3+0]=w[K];
  
  dmom_w_adv_v[K*3+1]=w[K];
  
  dmom_w_adv_w[K*3+0]=u[K];
  dmom_w_adv_w[K*3+1]=v[K];
  dmom_w_adv_w[K*3+2]=2.0*w[K];
  
  //u momentum diffusion tensor
  mom_u_diff_ten[K*9+0] = 2.0*nu;
  mom_u_diff_ten[K*9+4] = nu;
  mom_u_diff_ten[K*9+8] = nu;
  
  mom_uv_diff_ten[K*9+3]=nu;
  
  mom_uw_diff_ten[K*9+6]=nu;
  
  //v momentum diffusion tensor
  mom_v_diff_ten[K*9+0] = nu;
  mom_v_diff_ten[K*9+4] = 2.0*nu;
  mom_v_diff_ten[K*9+8] = nu;
  
  mom_vu_diff_ten[K*9+1]=nu;
  
  mom_vw_diff_ten[K*9+7]=nu;
  
  //w momentum diffusion tensor
  mom_w_diff_ten[K*9+0] = nu;
  mom_w_diff_ten[K*9+4] = nu;
  mom_w_diff_ten[K*9+8] = 2.0*nu;
  
  mom_wu_diff_ten[K*9+2]=nu;
  
  mom_wv_diff_ten[K*9+5]=nu;
  
  //momentum sources
  norm_n = sqrt(n[K*3+0]*n[K*3+0]+n[K*3+1]*n[K*3+1]+n[K*3+2]*n[K*3+2]);
  mom_u_source[K] = -g[0] - d_mu*sigma*kappa[K]*n[K*3+0]/(rho*(norm_n+1.0e-8));
  mom_v_source[K] = -g[1] - d_mu*sigma*kappa[K]*n[K*3+1]/(rho*(norm_n+1.0e-8));
  mom_w_source[K] = -g[2] - d_mu*sigma*kappa[K]*n[K*3+2]/(rho*(norm_n+1.0e-8));
   
  //u momentum Hamiltonian (pressure)
  mom_u_ham[K] = grad_p[K*3+0]/rho;
  dmom_u_ham_grad_p[K*3+0]=1.0/rho;
  
  //v momentum Hamiltonian (pressure)
  mom_v_ham[K] = grad_p[K*3+1]/rho;
  dmom_v_ham_grad_p[K*3+1]=1.0/rho;
  
  //w momentum Hamiltonian (pressure)
  mom_w_ham[K] = grad_p[K*3+2]/rho;
  dmom_w_ham_grad_p[K*3+2]=1.0/rho;
}
  
void backwardEuler(int eN, int k, double dt, double* m_old, double* m, double* dm, double* mt, double* dmt)
{  
  mt[eN*nQuadraturePoints_element+
     k]
    =
    (
     m[eN*nQuadraturePoints_element+
       k]
     -m_old[eN*nQuadraturePoints_element+
            k]
     )/dt;
  dmt[eN*nQuadraturePoints_element+
      k]
    = dm[eN*nQuadraturePoints_element+
         k]/dt;
}

void updateAdvection_strong(int eN,
			    int k,
			    int I,
			    double* df,
			    double* grad_u,
			    double* strong_residual)
{
  strong_residual[eN*nQuadraturePoints_element+
                  k] 
    +=
    df[eN*nQuadraturePoints_element*nSpace + 
       k*nSpace + 
       I]
    *
    grad_u[eN*nQuadraturePoints_element*nSpace + 
           k*nSpace + 
           I];
}

void updateAdvectionJacobian_strong(int eN,
                                    int k,
                                    int j,
                                    int I,
                                    double* df,
                                    double* grad_v,
                                    double* dstrong_residual)
{
  dstrong_residual[eN*nQuadraturePoints_element*nDOF_trial_element+
                   k*nDOF_trial_element + 
                   j] 
    +=
    df[eN*nQuadraturePoints_element*nSpace + 
       k*nSpace + 
       I]
    *
    grad_v[eN*nQuadraturePoints_element*nSpace*nDOF_trial_element + 
           k*nSpace*nDOF_trial_element +
           j*nSpace+
           I];
}

void updateAdvection_adjoint(int eN,
			     int k,
			     int i,
			     int I,
			     double* df,
			     double* grad_w_dV,
			     double* Lstar_w_dV)
{
  Lstar_w_dV[eN*nQuadraturePoints_element*nDOF_test_element + 
             k*nDOF_test_element + 
             i] 
    -=  
    df[eN*nQuadraturePoints_element*nSpace + 
       k*nSpace + 
       I]
    *
    grad_w_dV[eN*nQuadraturePoints_element*nDOF_test_element*nSpace + 
              k*nDOF_test_element*nSpace + 
              i*nSpace + 
              I];
}

void updateHamiltonian_strong(int eN,
			      int k,
			      int I,
			      double* dH,
			      double* grad_u,
			      double* strong_residual)
{
  strong_residual[eN*nQuadraturePoints_element+k] 
    += 
    dH[eN*nQuadraturePoints_element*nSpace+
       k*nSpace+
       I]
    *
    grad_u[eN*nQuadraturePoints_element*nSpace+
           k*nSpace+
           I];
}

void updateMass_strong(int eN,
		       int k,
		       double* mt,
		       double* strong_residual)
{
  strong_residual[eN*nQuadraturePoints_element+
		  k] 
    +=
    mt[eN*nQuadraturePoints_element+
       k]; 
}

void updateReaction_strong(int eN,
			   int k,
			   double* r,
			   double* strong_residual)
{
  strong_residual[eN*nQuadraturePoints_element+
		  k] 
    += 
    r[eN*nQuadraturePoints_element+
      k]; 
}

void calculateSubgridError_tau(int eN,
                               int k,
                               int nSpace,
                               double  hFactor,
                               double* elementDiameter,
                               double* dmt,
                               double* dm,
                               double* f,
                               double* a,
                               double* tau0,
                               double* tau1,
                               double* cfl)
{
  double h,oneByAbsdt,density,viscosity,Re_max=0.0,CFL_max=0.0,nrm_v;
  h = hFactor*elementDiameter[eN];
  density = dm[eN*nQuadraturePoints_element+
               k];
  viscosity =  a[eN*nQuadraturePoints_element*nSpace2 + 
                 k*nSpace2+
                 nSpace+1]; 
  nrm_v=0.0;
  for(I=0;I<nSpace;I++)
    nrm_v+=f[eN*nQuadraturePoints_element*nSpace+
             k*nSpace+
             I]*
      f[eN*nQuadraturePoints_element*nSpace+
        k*nSpace+
        I];
  nrm_v = sqrt(nrm_v);
  Re_max = fmax(nrm_v*h/(viscosity/density),Re_max);
  CFL_max = fmax(nrm_v/h,CFL_max);
  cfl[eN*nQuadraturePoints_element+k] = nrm_v/h;
  oneByAbsdt =  fabs(dmt[eN*nQuadraturePoints_element+k]);
  
  tau0[eN*nQuadraturePoints_element+k] = 1.0/(4.0*viscosity/(h*h) +
                                              2.0*density*nrm_v/h +
                                              oneByAbsdt);
  
  tau1[eN*nQuadraturePoints_element+k] = 4.0*viscosity +
    2.0*density*nrm_v*h+
    oneByAbsdt*h*h;
}

void calculateSubgridError_tauRes(int eN,
                                  int k,
                                  double* tau0,
                                  double* tau1,
                                  double* pdeResidualP,
                                  double* dpdeResidualP_du,
                                  double* dpdeResidualP_dv,
                                  double* dpdeResidualP_dw,
                                  double* pdeResidualU,
                                  double* dpdeResidualU_dp,
                                  double* dpdeResidualU_du,
                                  double* dpdeResidualU_dv,
                                  double* dpdeResidualU_dw,
                                  double* pdeResidualV,
                                  double* dpdeResidualV_dp,
                                  double* dpdeResidualV_du,
                                  double* dpdeResidualV_dv,
                                  double* dpdeResidualV_dw,
                                  double* pdeResidualW,
                                  double* dpdeResidualW_dp,
                                  double* dpdeResidualW_du,
                                  double* dpdeResidualW_dv,
                                  double* dpdeResidualW_dw,
                                  double* subgridErrorP,
                                  double* dsubgridErrorP_du,
                                  double* dsubgridErrorP_dv,
                                  double* dsubgridErrorP_dw,
                                  double* subgridErrorU,
                                  double* dsubgridErrorU_dp,
                                  double* dsubgridErrorU_du,
                                  double* dsubgridErrorU_dv,
                                  double* dsubgridErrorU_dw,
                                  double* subgridErrorV,
                                  double* dsubgridErrorV_dp,
                                  double* dsubgridErrorV_du,
                                  double* dsubgridErrorV_dv,
                                  double* dsubgridErrorV_dw,
                                  double* subgridErrorW,
                                  double* dsubgridErrorW_dp,
                                  double* dsubgridErrorW_du,
                                  double* dsubgridErrorW_dv,
                                  double* dsubgridErrorW_dw)
{
  /* GLS momentum */
  subgridErrorU[eN*nQuadraturePoints_element+
                k] =
    tau0[eN*nQuadraturePoints_element+k]
    *
    pdeResidualU[eN*nQuadraturePoints_element+
                 k];
  subgridErrorV[eN*nQuadraturePoints_element+
                k] = 
    tau0[eN*nQuadraturePoints_element+k]
    *
    pdeResidualV[eN*nQuadraturePoints_element+
                 k];
  subgridErrorW[eN*nQuadraturePoints_element+
                k] = 
    tau0[eN*nQuadraturePoints_element+k]
    *
    pdeResidualW[eN*nQuadraturePoints_element+
                 k];
  /* GLS pressure */
  subgridErrorP[eN*nQuadraturePoints_element+
                k] =
    tau1[eN*nQuadraturePoints_element+k]
    *pdeResidualP[eN*nQuadraturePoints_element+
                  k];
  for (j=0;j<nDOF_trial_element;j++)
    {
      /* GLS  momentum*/
      /* u */
      dsubgridErrorU_dp[eN*nQuadraturePoints_element*nDOF_trial_element+
                        k*nDOF_trial_element+
                        j] =
        tau0[eN*nQuadraturePoints_element+k]
        *
        dpdeResidualU_dp[eN*nQuadraturePoints_element*nDOF_trial_element+
                         k*nDOF_trial_element+
                         j];
      dsubgridErrorU_du[eN*nQuadraturePoints_element*nDOF_trial_element+
                        k*nDOF_trial_element+
                        j] = 
        tau0[eN*nQuadraturePoints_element+k]
        *
        dpdeResidualU_du[eN*nQuadraturePoints_element*nDOF_trial_element+
                         k*nDOF_trial_element+
                         j];
      dsubgridErrorU_dv[eN*nQuadraturePoints_element*nDOF_trial_element+
                        k*nDOF_trial_element+
                        j] = 
        tau0[eN*nQuadraturePoints_element+k]
        *
        dpdeResidualU_dv[eN*nQuadraturePoints_element*nDOF_trial_element+
                         k*nDOF_trial_element+
                         j];
      dsubgridErrorU_dw[eN*nQuadraturePoints_element*nDOF_trial_element+
                        k*nDOF_trial_element+
                        j] = 
        tau0[eN*nQuadraturePoints_element+k]
        *
        dpdeResidualU_dw[eN*nQuadraturePoints_element*nDOF_trial_element+
                         k*nDOF_trial_element+
                         j];
      /* v */
      dsubgridErrorV_dp[eN*nQuadraturePoints_element*nDOF_trial_element+
                        k*nDOF_trial_element+
                        j] = 
        tau0[eN*nQuadraturePoints_element+k]
        *
        dpdeResidualV_dp[eN*nQuadraturePoints_element*nDOF_trial_element+
                         k*nDOF_trial_element+
                         j];
      dsubgridErrorV_du[eN*nQuadraturePoints_element*nDOF_trial_element+
                        k*nDOF_trial_element+
                        j] = 
        tau0[eN*nQuadraturePoints_element+k]
        *
        dpdeResidualV_du[eN*nQuadraturePoints_element*nDOF_trial_element+
                         k*nDOF_trial_element+
                         j];
      dsubgridErrorV_dv[eN*nQuadraturePoints_element*nDOF_trial_element+
                        k*nDOF_trial_element+
                        j] = 
        tau0[eN*nQuadraturePoints_element+k]
        *
        dpdeResidualV_dv[eN*nQuadraturePoints_element*nDOF_trial_element+
                         k*nDOF_trial_element+
                         j];
      dsubgridErrorV_dw[eN*nQuadraturePoints_element*nDOF_trial_element+
                        k*nDOF_trial_element+
                        j] = 
        tau0[eN*nQuadraturePoints_element+k]
        *
        dpdeResidualV_dw[eN*nQuadraturePoints_element*nDOF_trial_element+
                         k*nDOF_trial_element+
                         j];
      /* w */
      dsubgridErrorW_dp[eN*nQuadraturePoints_element*nDOF_trial_element+
                        k*nDOF_trial_element+
                        j] = 
        tau0[eN*nQuadraturePoints_element+k]
        *
        dpdeResidualW_dp[eN*nQuadraturePoints_element*nDOF_trial_element+
                                 k*nDOF_trial_element+
                         j];
      dsubgridErrorW_du[eN*nQuadraturePoints_element*nDOF_trial_element+
                        k*nDOF_trial_element+
                        j] = 
        tau0[eN*nQuadraturePoints_element+k]
        *
        dpdeResidualW_du[eN*nQuadraturePoints_element*nDOF_trial_element+
                         k*nDOF_trial_element+
                         j];
      dsubgridErrorW_dv[eN*nQuadraturePoints_element*nDOF_trial_element+
                        k*nDOF_trial_element+
                        j] = 
        tau0[eN*nQuadraturePoints_element+k]
        *
        dpdeResidualW_dv[eN*nQuadraturePoints_element*nDOF_trial_element+
                         k*nDOF_trial_element+
                         j];
      dsubgridErrorW_dw[eN*nQuadraturePoints_element*nDOF_trial_element+
                        k*nDOF_trial_element+
                        j] = 
        tau0[eN*nQuadraturePoints_element+k]
        *
        dpdeResidualW_dw[eN*nQuadraturePoints_element*nDOF_trial_element+
                         k*nDOF_trial_element+
                         j];
      /* GLS pressure */
      dsubgridErrorP_du[eN*nQuadraturePoints_element*nDOF_trial_element+
                        k*nDOF_trial_element+
                        j]
        =
        tau1[eN*nQuadraturePoints_element+k]
        *
        dpdeResidualP_du[eN*nQuadraturePoints_element*nDOF_trial_element+
                         k*nDOF_trial_element+
                         j];
      dsubgridErrorP_dv[eN*nQuadraturePoints_element*nDOF_trial_element+
                        k*nDOF_trial_element+
                        j]
        =
        tau1[eN*nQuadraturePoints_element+k]
        *
        dpdeResidualP_dv[eN*nQuadraturePoints_element*nDOF_trial_element+
                         k*nDOF_trial_element+
                         j];
      dsubgridErrorP_dw[eN*nQuadraturePoints_element*nDOF_trial_element+
                        k*nDOF_trial_element+
                        j]
        =
        tau1[eN*nQuadraturePoints_element+k]
        *
        dpdeResidualP_dw[eN*nQuadraturePoints_element*nDOF_trial_element+
                         k*nDOF_trial_element+
                         j];
    }
}

void calculateNumericalDiffusionResGradQuad(int eN,
                                            int k,
                                            int nSpace,
                                            double shockCapturingDiffusion,
                                            double* elementDiameter,
                                            double* strong_residual,
                                            double* grad_u,
                                            double* numDiff)
{
  int I;
  double h,
    num,
    den,
    n_grad_u;
  h = elementDiameter[eN];
  n_grad_u = 0.0;
  for (I=0;I<nSpace;I++)
    {
      n_grad_u += grad_u[eN*nQuadraturePoints_element*nSpace+k*nSpace+I]*grad_u[eN*nQuadraturePoints_element*nSpace+k*nSpace+I];
    }
  num = shockCapturingDiffusion*0.5*h*fabs(strong_residual[eN*nQuadraturePoints_element+k]);
  den = sqrt(n_grad_u) + 1.0e-8;
  numDiff[eN*nQuadraturePoints_element+k] = num/den;
}


#ifndef COMPKERNEL_H
#define COMPKERNEL_H
#include <cmath>

//I separated the space mapping part of the kernel so I could partially specializethe template on NSPACE
template<const int NSPACE, const int NDOF_MESH_TRIAL_ELEMENT>
class CompKernelSpaceMapping
{
public:
  inline void calculateMapping_element(const int eN,
				       const int k,
				       double* mesh_dof,
				       int* mesh_l2g,
				       double* mesh_trial_ref,
				       double* mesh_grad_trial_ref,
				       double* jac,
				       double& jacDet,
				       double* jacInv,
				       double& x,
				       double& y,
				       double& z);
  
  inline void calculateMapping_elementBoundary(const int eN,
					       const int ebN_local,
					       const int kb,
					       const int ebN_local_kb,
					       double* mesh_dof,
					       int* mesh_l2g,
					       double* mesh_trial_trace_ref,
					       double* mesh_grad_trial_trace_ref,
					       double* boundaryJac_ref,
					       double* jac,
					       double& jacDet,
					       double* jacInv,
					       double* boundaryJac,
					       double* metricTensor,
					       double& metricTensorDetSqrt,
					       double* normal_ref,
					       double* normal,
					       double& x,
					       double& y,
					       double& z);
};

//specialization for 3D
template<const int NDOF_MESH_TRIAL_ELEMENT>
class CompKernelSpaceMapping<3,NDOF_MESH_TRIAL_ELEMENT>
{
public:
  inline void calculateMapping_element(const int eN,
				       const int k,
				       double* mesh_dof,
				       int* mesh_l2g,
				       double* mesh_trial_ref,
				       double* mesh_grad_trial_ref,
				       double* jac,
				       double& jacDet,
				       double* jacInv,
				       double& x,
				       double& y,
				       double& z)
  {
    const int X=0,Y=1,Z=2,
      XX=0,XY=1,XZ=2,
      YX=3,YY=NDOF_MESH_TRIAL_ELEMENT,YZ=5,
      ZX=6,ZY=7,ZZ=8;
    
    register double Grad_x[3],Grad_y[3],Grad_z[3],oneOverJacDet;
    
    //
    //mapping of reference element to physical element
    //
    x=0.0;y=0.0;z=0.0;
    for (int I=0;I<3;I++)
      {
	Grad_x[I]=0.0;Grad_y[I]=0.0;Grad_z[I]=0.0;
      }
    for (int j=0;j<NDOF_MESH_TRIAL_ELEMENT;j++)
      {
	int eN_j=eN*NDOF_MESH_TRIAL_ELEMENT+j;
	x += mesh_dof[mesh_l2g[eN_j]*3+0]*mesh_trial_ref[k*NDOF_MESH_TRIAL_ELEMENT+j];
	y += mesh_dof[mesh_l2g[eN_j]*3+1]*mesh_trial_ref[k*NDOF_MESH_TRIAL_ELEMENT+j];
	z += mesh_dof[mesh_l2g[eN_j]*3+2]*mesh_trial_ref[k*NDOF_MESH_TRIAL_ELEMENT+j];	      
	for (int I=0;I<3;I++)
	  {
	    Grad_x[I] += mesh_dof[mesh_l2g[eN_j]*3+0]*mesh_grad_trial_ref[k*NDOF_MESH_TRIAL_ELEMENT*3+j*3+I];
	    Grad_y[I] += mesh_dof[mesh_l2g[eN_j]*3+1]*mesh_grad_trial_ref[k*NDOF_MESH_TRIAL_ELEMENT*3+j*3+I];
	    Grad_z[I] += mesh_dof[mesh_l2g[eN_j]*3+2]*mesh_grad_trial_ref[k*NDOF_MESH_TRIAL_ELEMENT*3+j*3+I];
	  }
      }
    jac[XX] = Grad_x[X];//node[X]*grad[X];
    jac[XY] = Grad_x[Y];//node[X]*grad[Y];
    jac[XZ] = Grad_x[Z];//node[X]*grad[Z];
    jac[YX] = Grad_y[X];//node[Y]*grad[X];
    jac[YY] = Grad_y[Y];//node[Y]*grad[Y];
    jac[YZ] = Grad_y[Z];//node[Y]*grad[Z];
    jac[ZX] = Grad_z[X];//node[Z]*grad[X];
    jac[ZY] = Grad_z[Y];//node[Z]*grad[Y];
    jac[ZZ] = Grad_z[Z];//node[Z]*grad[Z];
    jacDet = 
      jac[XX]*(jac[YY]*jac[ZZ] - jac[YZ]*jac[ZY]) -
      jac[XY]*(jac[YX]*jac[ZZ] - jac[YZ]*jac[ZX]) +
      jac[XZ]*(jac[YX]*jac[ZY] - jac[YY]*jac[ZX]);
    oneOverJacDet = 1.0/jacDet;
    jacInv[XX] = oneOverJacDet*(jac[YY]*jac[ZZ] - jac[YZ]*jac[ZY]);
    jacInv[YX] = oneOverJacDet*(jac[YZ]*jac[ZX] - jac[YX]*jac[ZZ]);
    jacInv[ZX] = oneOverJacDet*(jac[YX]*jac[ZY] - jac[YY]*jac[ZX]);
    jacInv[XY] = oneOverJacDet*(jac[ZY]*jac[XZ] - jac[ZZ]*jac[XY]);
    jacInv[YY] = oneOverJacDet*(jac[ZZ]*jac[XX] - jac[ZX]*jac[XZ]);
    jacInv[ZY] = oneOverJacDet*(jac[ZX]*jac[XY] - jac[ZY]*jac[XX]);
    jacInv[XZ] = oneOverJacDet*(jac[XY]*jac[YZ] - jac[XZ]*jac[YY]);
    jacInv[YZ] = oneOverJacDet*(jac[XZ]*jac[YX] - jac[XX]*jac[YZ]);
    jacInv[ZZ] = oneOverJacDet*(jac[XX]*jac[YY] - jac[XY]*jac[YX]);
  }
  
  inline void calculateMapping_elementBoundary(const int eN,
					       const int ebN_local,
					       const int kb,
					       const int ebN_local_kb,
					       double* mesh_dof,
					       int* mesh_l2g,
					       double* mesh_trial_trace_ref,
					       double* mesh_grad_trial_trace_ref,
					       double* boundaryJac_ref,
					       double* jac,
					       double& jacDet,
					       double* jacInv,
					       double* boundaryJac,
					       double* metricTensor,
					       double& metricTensorDetSqrt,
					       double* normal_ref,
					       double* normal,
					       double& x,
					       double& y,
					       double& z)
  {
    const int X=0,Y=1,Z=2,
      XX=0,XY=1,XZ=2,
      YX=3,YY=NDOF_MESH_TRIAL_ELEMENT,YZ=5,
      ZX=6,ZY=7,ZZ=8,
      XHX=0,XHY=1,
      YHX=2,YHY=3,
      ZHX=NDOF_MESH_TRIAL_ELEMENT,ZHY=5,
      HXHX=0,HXHY=1,
      HYHX=2,HYHY=3;
    const int ebN_local_kb_nSpace = ebN_local_kb*3;
  
    register double Grad_x_ext[3],Grad_y_ext[3],Grad_z_ext[3],oneOverJacDet,norm_normal=0.0;
    // 
    //calculate mapping from the reference element to the physical element
    // 
    x=0.0;y=0.0;z=0.0;
    for (int I=0;I<3;I++)
      {
	Grad_x_ext[I] = 0.0;
	Grad_y_ext[I] = 0.0;
	Grad_z_ext[I] = 0.0;
      }
    for (int j=0;j<NDOF_MESH_TRIAL_ELEMENT;j++) 
      { 
	int eN_j = eN*NDOF_MESH_TRIAL_ELEMENT+j;
	int ebN_local_kb_j = ebN_local_kb*NDOF_MESH_TRIAL_ELEMENT+j;
	int ebN_local_kb_j_nSpace = ebN_local_kb_j*3;
	x += mesh_dof[mesh_l2g[eN_j]*3+0]*mesh_trial_trace_ref[ebN_local_kb_j]; 
	y += mesh_dof[mesh_l2g[eN_j]*3+1]*mesh_trial_trace_ref[ebN_local_kb_j]; 
	z += mesh_dof[mesh_l2g[eN_j]*3+2]*mesh_trial_trace_ref[ebN_local_kb_j]; 
	for (int I=0;I<3;I++)
	  {
	    Grad_x_ext[I] += mesh_dof[mesh_l2g[eN_j]*3+0]*mesh_grad_trial_trace_ref[ebN_local_kb_j_nSpace+I];
	    Grad_y_ext[I] += mesh_dof[mesh_l2g[eN_j]*3+1]*mesh_grad_trial_trace_ref[ebN_local_kb_j_nSpace+I]; 
	    Grad_z_ext[I] += mesh_dof[mesh_l2g[eN_j]*3+2]*mesh_grad_trial_trace_ref[ebN_local_kb_j_nSpace+I];
	  } 
      }
    //Space Mapping Jacobian
    jac[XX] = Grad_x_ext[X];//node[X]*grad[X];
    jac[XY] = Grad_x_ext[Y];//node[X]*grad[Y];
    jac[XZ] = Grad_x_ext[Z];//node[X]*grad[Z];
    jac[YX] = Grad_y_ext[X];//node[Y]*grad[X];
    jac[YY] = Grad_y_ext[Y];//node[Y]*grad[Y];
    jac[YZ] = Grad_y_ext[Z];//node[Y]*grad[Z];
    jac[ZX] = Grad_z_ext[X];//node[Z]*grad[X];
    jac[ZY] = Grad_z_ext[Y];//node[Z]*grad[Y];
    jac[ZZ] = Grad_z_ext[Z];//node[Z]*grad[Z];
    jacDet = 
      jac[XX]*(jac[YY]*jac[ZZ] - jac[YZ]*jac[ZY]) -
      jac[XY]*(jac[YX]*jac[ZZ] - jac[YZ]*jac[ZX]) +
      jac[XZ]*(jac[YX]*jac[ZY] - jac[YY]*jac[ZX]);
    oneOverJacDet = 1.0/jacDet;
    jacInv[XX] = oneOverJacDet*(jac[YY]*jac[ZZ] - jac[YZ]*jac[ZY]);
    jacInv[YX] = oneOverJacDet*(jac[YZ]*jac[ZX] - jac[YX]*jac[ZZ]);
    jacInv[ZX] = oneOverJacDet*(jac[YX]*jac[ZY] - jac[YY]*jac[ZX]);
    jacInv[XY] = oneOverJacDet*(jac[ZY]*jac[XZ] - jac[ZZ]*jac[XY]);
    jacInv[YY] = oneOverJacDet*(jac[ZZ]*jac[XX] - jac[ZX]*jac[XZ]);
    jacInv[ZY] = oneOverJacDet*(jac[ZX]*jac[XY] - jac[ZY]*jac[XX]);
    jacInv[XZ] = oneOverJacDet*(jac[XY]*jac[YZ] - jac[XZ]*jac[YY]);
    jacInv[YZ] = oneOverJacDet*(jac[XZ]*jac[YX] - jac[XX]*jac[YZ]);
    jacInv[ZZ] = oneOverJacDet*(jac[XX]*jac[YY] - jac[XY]*jac[YX]);
    //normal
    norm_normal=0.0;
    for (int I=0;I<3;I++)
      normal[I] = 0.0;
    for (int I=0;I<3;I++)
      {
	for (int J=0;J<3;J++)
	  normal[I] += jacInv[J*3+I]*normal_ref[ebN_local_kb_nSpace+J];
	norm_normal+=normal[I]*normal[I];
      }
    norm_normal = sqrt(norm_normal);
    for (int I=0;I<3;I++)
      normal[I] /= norm_normal;
    //metric tensor and determinant
    boundaryJac[XHX] = jac[XX]*boundaryJac_ref[ebN_local*6+XHX]+jac[XY]*boundaryJac_ref[ebN_local*6+YHX]+jac[XZ]*boundaryJac_ref[ebN_local*6+ZHX];
    boundaryJac[XHY] = jac[XX]*boundaryJac_ref[ebN_local*6+XHY]+jac[XY]*boundaryJac_ref[ebN_local*6+YHY]+jac[XZ]*boundaryJac_ref[ebN_local*6+ZHY];
    boundaryJac[YHX] = jac[YX]*boundaryJac_ref[ebN_local*6+XHX]+jac[YY]*boundaryJac_ref[ebN_local*6+YHX]+jac[YZ]*boundaryJac_ref[ebN_local*6+ZHX];
    boundaryJac[YHY] = jac[YX]*boundaryJac_ref[ebN_local*6+XHY]+jac[YY]*boundaryJac_ref[ebN_local*6+YHY]+jac[YZ]*boundaryJac_ref[ebN_local*6+ZHY];
    boundaryJac[ZHX] = jac[ZX]*boundaryJac_ref[ebN_local*6+XHX]+jac[ZY]*boundaryJac_ref[ebN_local*6+YHX]+jac[ZZ]*boundaryJac_ref[ebN_local*6+ZHX];
    boundaryJac[ZHY] = jac[ZX]*boundaryJac_ref[ebN_local*6+XHY]+jac[ZY]*boundaryJac_ref[ebN_local*6+YHY]+jac[ZZ]*boundaryJac_ref[ebN_local*6+ZHY];
  
    metricTensor[HXHX] = boundaryJac[XHX]*boundaryJac[XHX]+boundaryJac[YHX]*boundaryJac[YHX]+boundaryJac[ZHX]*boundaryJac[ZHX];
    metricTensor[HXHY] = boundaryJac[XHX]*boundaryJac[XHY]+boundaryJac[YHX]*boundaryJac[YHY]+boundaryJac[ZHX]*boundaryJac[ZHY];
    metricTensor[HYHX] = boundaryJac[XHY]*boundaryJac[XHX]+boundaryJac[YHY]*boundaryJac[YHX]+boundaryJac[ZHY]*boundaryJac[ZHX];
    metricTensor[HYHY] = boundaryJac[XHY]*boundaryJac[XHY]+boundaryJac[YHY]*boundaryJac[YHY]+boundaryJac[ZHY]*boundaryJac[ZHY];
  
    metricTensorDetSqrt=sqrt(metricTensor[HXHX]*metricTensor[HYHY]- metricTensor[HXHY]*metricTensor[HYHX]);
  }
};

template<const int NSPACE, const int NDOF_MESH_TRIAL_ELEMENT, const int NDOF_TRIAL_ELEMENT, const int NDOF_TEST_ELEMENT>
class CompKernel
{
public:
  CompKernelSpaceMapping<NSPACE,NDOF_MESH_TRIAL_ELEMENT> mapping;

  inline void valFromDOF(const double* dof,const int* l2g_element,const double* trial_ref,double& val)
  {
    val=0.0;
    for (int j=0;j<NDOF_TRIAL_ELEMENT;j++)
      val+=dof[l2g_element[j]]*trial_ref[j];
  }

  inline void gradFromDOF(const double* dof,const int* l2g_element,const double* grad_trial,double* grad)
  {
    for(int I=0;I<NSPACE;I++)
      grad[I] = 0.0;
    for (int j=0;j<NDOF_TRIAL_ELEMENT;j++)
      for(int I=0;I<NSPACE;I++)
	grad[I] += dof[l2g_element[j]]*grad_trial[j*NSPACE+I];
  }

  inline void gradTrialFromRef(const double* grad_trial_ref, const double* jacInv, double* grad_trial)
  {
    for (int j=0;j<NDOF_TRIAL_ELEMENT;j++)
      for(int I=0;I<NSPACE;I++)
	grad_trial[j*NSPACE+I] = 0.0;
    for (int j=0;j<NDOF_TRIAL_ELEMENT;j++)
      for(int I=0;I<NSPACE;I++)
	for(int J=0;J<NSPACE;J++)
	  grad_trial[j*NSPACE+I] += jacInv[J*NSPACE+I]*grad_trial_ref[j*NSPACE+J];
  }
  
  inline void gradTestFromRef(const double* grad_test_ref, const double* jacInv, double* grad_test)
  {
    for (int i=0;i<NDOF_TEST_ELEMENT;i++)
      for(int I=0;I<NSPACE;I++)
	grad_test[i*NSPACE+I] = 0.0;
    for (int i=0;i<NDOF_TEST_ELEMENT;i++)
      for(int I=0;I<NSPACE;I++)
	for(int J=0;J<NSPACE;J++)
	  grad_test[i*NSPACE+I] += jacInv[J*NSPACE+I]*grad_test_ref[i*NSPACE+J];
  }
  
  inline void backwardEuler(const double& dt, const double& m_old, const double& m, const double& dm, double& mt, double& dmt)
  {  
    mt =(m-m_old)/dt;
    dmt = dm/dt;
  }

  inline void bdf(const double& alpha, const double& beta, const double& m, const double& dm, double& mt, double& dmt)
  {  
    mt =alpha*m + beta;
    dmt = alpha*dm;
  }

  inline double Mass_weak(const double& mt, const double& w_dV)
  {
    return mt*w_dV;
  }

  inline double MassJacobian_weak(const double& dmt,
				  const double& v,
				  const double& w_dV)
  {
    return dmt*v*w_dV;
  }

  inline double Mass_strong(const double& mt)
  {
    return mt; 
  }

  inline double MassJacobian_strong(const double& dmt,
				    const double& v)
  {
    return dmt*v;
  }

  inline double Mass_adjoint(const double& dmt,
			     const double& w_dV)
  {
    return dmt*w_dV;
  }

  inline double Advection_weak(const double  f[NSPACE],
			       const double grad_w_dV[NSPACE])
  {
    double tmp=0.0;
    for(int I=0;I<NSPACE;I++)
      tmp -= f[I]*grad_w_dV[I];
    return tmp;
  }

  inline double AdvectionJacobian_weak(const double df[NSPACE],
				       const double& v,
				       const double grad_w_dV[NSPACE])
  {
    double tmp=0.0;
    for(int I=0;I<NSPACE;I++)
      tmp -= df[I]*v*grad_w_dV[I];
    return tmp;
  }

  inline double Advection_strong(const double df[NSPACE],
				 const double grad_u[NSPACE])
  {
    double tmp=0.0;
    for(int I=0;I<NSPACE;I++)
      tmp += df[I]*grad_u[I];
    return tmp;
  }

  inline double AdvectionJacobian_strong(const double df[NSPACE],
					 const double grad_v[NSPACE])
  {
    double tmp=0.0;
    for(int I=0;I<NSPACE;I++)
      tmp += df[I]*grad_v[I];
    return tmp;
  }

  inline double Advection_adjoint(const double df[NSPACE],
				  const double grad_w_dV[NSPACE])
  {
    double tmp=0.0;
    for(int I=0;I<NSPACE;I++)
      tmp -= df[I]*grad_w_dV[I];
    return tmp;
  }

  inline double Hamiltonian_weak(const double& H,
				 const double& w_dV)
  {
    return H*w_dV;
  }

  inline double HamiltonianJacobian_weak(const double dH[NSPACE],
					 const double grad_v[NSPACE],
					 const double& w_dV)
  {
    double tmp=0.0;
    for(int I=0;I<NSPACE;I++)
      tmp += dH[I]*grad_v[I]*w_dV;
    return tmp;
  }

  inline double Hamiltonian_strong(const double dH[NSPACE],
				   const double grad_u[NSPACE])
  {
    double tmp=0.0;
    for(int I=0;I<NSPACE;I++)
      tmp += dH[I]*grad_u[I];
    return tmp;
  }

  inline double HamiltonianJacobian_strong(const double dH[NSPACE],
					   const double grad_v[NSPACE])
  {
    double tmp=0.0;
    for(int I=0;I<NSPACE;I++)
      tmp += dH[I]*grad_v[I];
    return tmp;
  }

  inline double Hamiltonian_adjoint(const double dH[NSPACE],
				    const double grad_w_dV[NSPACE])
  {
    double tmp=0.0;
    for(int I=0;I<NSPACE;I++)
      tmp -= dH[I]*grad_w_dV[I];
    return tmp;
  }

  inline double Diffusion_weak(int* rowptr,
			       int* colind,
			       double* a,
			       const double grad_phi[NSPACE],
			       const double grad_w_dV[NSPACE])
  {
    double tmp=0.0;
    for(int I=0;I<NSPACE;I++)
      for (int m=rowptr[I];m<rowptr[I+1];m++)
	tmp += a[m]*grad_phi[colind[m]]*grad_w_dV[I];
    return tmp;
  }

  inline double DiffusionJacobian_weak(int* rowptr,
				       int* colind,
				       double* a,
				       double* da,
				       const double grad_phi[NSPACE],
				       const double grad_w_dV[NSPACE],
				       const double& dphi,
				       const double& v,
				       const double grad_v[NSPACE])
  {
    double daProduct=0.0,dphiProduct=0.0;
    for (int I=0;I<NSPACE;I++)
      for (int m=rowptr[I];m<rowptr[I+1];m++)
	{
	  daProduct += da[m]*grad_phi[colind[m]]*grad_w_dV[I];
	  dphiProduct += a[m]*grad_v[colind[m]]*grad_w_dV[I];
	}
    return daProduct*v+dphiProduct*dphi;
  }

  inline double SimpleDiffusionJacobian_weak(int* rowptr,
					     int* colind,
					     double* a,
					     const double grad_v[NSPACE],
					     const double grad_w_dV[NSPACE])
  {
    double dphiProduct=0.0;
    for (int I=0;I<NSPACE;I++)
      for (int m=rowptr[I];m<rowptr[I+1];m++)
	{
	  dphiProduct += a[m]*grad_v[colind[m]]*grad_w_dV[I];
	}
    return dphiProduct;
  }

  inline double Reaction_weak(const double& r,
			      const double& w_dV)
  {
    return r*w_dV;
  }
  
  inline double ReactionJacobian_weak(const double& dr,
				      const double& v,
				      const double& w_dV)
  {
    return dr*v*w_dV;
  }

  inline double Reaction_strong(const double& r)
  {
    return r; 
  }

  inline double ReactionJacobian_strong(const double& dr,
					const double& v)
  {
    return dr*v;
  }

  inline double Reaction_adjoint(const double& dr,
				 const double& w_dV)
  {
    return dr*w_dV;
  }

  inline void calculateNumericalDiffusion(const double& shockCapturingDiffusion,
					  const double& elementDiameter,
					  const double& strong_residual,
					  const double grad_u[NSPACE],
					  double& numDiff)
  {
    double h,
      num,
      den,
      n_grad_u;
    h = elementDiameter;
    n_grad_u = 0.0;
    for (int I=0;I<NSPACE;I++)
      n_grad_u += grad_u[I]*grad_u[I];
    num = shockCapturingDiffusion*0.5*h*fabs(strong_residual);
    den = sqrt(n_grad_u) + 1.0e-8;
    numDiff = num/den;
  }

  inline double SubgridError(const double& error,
			     const double& Lstar_w_dV)
  {
    return -error*Lstar_w_dV;
  }

  inline double SubgridErrorJacobian(const double& derror,
				     const double& Lstar_w_dV)
  {
    return -derror*Lstar_w_dV;
  }

  inline double NumericalDiffusion(const double& numDiff,
				   const double grad_u[NSPACE],
				   const double grad_w_dV[NSPACE])
  {
    double tmp=0.0;
    for (int I=0;I<NSPACE;I++)
      tmp +=  numDiff*grad_u[I]*grad_w_dV[I];
    return tmp;
  }

  inline double NumericalDiffusionJacobian(const double& numDiff,
					   const double grad_v[NSPACE],
					   const double grad_w_dV[NSPACE])
  {
    double tmp=0.0;
    for (int I=0;I<NSPACE;I++)
      tmp += numDiff*grad_v[I]*grad_w_dV[I];
    return tmp;
  }

 

  inline double ExteriorElementBoundaryFlux(const double& flux,
					    const double& w_dS)
  {
    return flux*w_dS;
  }

  inline double ExteriorNumericalAdvectiveFluxJacobian(const double& dflux_left,
						       const double& v)
  {
    return dflux_left*v;
  }

  inline void calculateMapping_element(const int eN,
				       const int k,
				       double* mesh_dof,
				       int* mesh_l2g,
				       double* mesh_trial_ref,
				       double* mesh_grad_trial_ref,
				       double* jac,
				       double& jacDet,
				       double* jacInv,
				       double& x,
				       double& y,
				       double& z)
  {
    mapping.calculateMapping_element(eN,k,mesh_dof,mesh_l2g,mesh_trial_ref,mesh_grad_trial_ref,jac,jacDet,jacInv,x,y,z);
  }

  inline 
  void calculateMapping_elementBoundary(const int eN,
					const int ebN_local,
					const int kb,
					const int ebN_local_kb,
					double* mesh_dof,
					int* mesh_l2g,
					double* mesh_trial_trace_ref,
					double* mesh_grad_trial_trace_ref,
					double* boundaryJac_ref,
					double* jac,
					double& jacDet,
					double* jacInv,
					double* boundaryJac,
					double* metricTensor,
					double& metricTensorDetSqrt,
					double* normal_ref,
					double* normal,
					double& x,
					double& y,
					double& z)
  {
    mapping.calculateMapping_elementBoundary(eN,
					     ebN_local,
					     kb,
					     ebN_local_kb,
					     mesh_dof,
					     mesh_l2g,
					     mesh_trial_trace_ref,
					     mesh_grad_trial_trace_ref,
					     boundaryJac_ref,
					     jac,
					     jacDet,
					     jacInv,
					     boundaryJac,
					     metricTensor,
					     metricTensorDetSqrt,
					     normal_ref,
					     normal,
					     x,
					     y,
					     z);
  }
};

#endif

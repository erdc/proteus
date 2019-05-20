#ifndef EQUIVALENT_POLYNOMIALS_H
#define EQUIVALENT_POLYNOMIALS_H
#include <cmath>
#include <cassert>
#include "equivalent_polynomials_coefficients.h"

namespace equivalent_polynomials
{
  template<int nSpace>
  inline void _calculate_normal(double* phys_nodes_cut, double* level_set_normal);

  template<>
  inline void _calculate_normal<1>(double* phys_nodes_cut, double* level_set_normal)
  {
    const unsigned int nSpace(1);
    level_set_normal[0] = 1.0;
  }

  template<>
  inline void _calculate_normal<2>(double* phys_nodes_cut, double* level_set_normal)
  {
    const unsigned int nSpace(2);
    double level_set_tangent[2];
    for(int I=0; I < nSpace; I++)
      {
        level_set_tangent[I] = phys_nodes_cut[3+I] - phys_nodes_cut[I];//nodes always 3D
      }
    level_set_normal[0] =   level_set_tangent[1];
    level_set_normal[1] = - level_set_tangent[0];
    double norm = std::sqrt(level_set_normal[0]*level_set_normal[0] +
                            level_set_normal[1]*level_set_normal[1]);
    for(int I=0; I < nSpace; I++)
      level_set_normal[I] /=  norm;
  }

  template<>
  inline void _calculate_normal<3>(double* phys_nodes_cut, double* level_set_normal)
  {
    const unsigned int nSpace(3);
    double level_set_tangent_a[3], level_set_tangent_b[3];
    for(int I=0; I < nSpace; I++)
      {
        level_set_tangent_a[I] = phys_nodes_cut[3+I] - phys_nodes_cut[I];
        level_set_tangent_b[I] = phys_nodes_cut[6+I] - phys_nodes_cut[I];
      }
    level_set_normal[0] =   level_set_tangent_a[1]*level_set_tangent_b[2]
      - level_set_tangent_a[2]*level_set_tangent_b[1];
    level_set_normal[1] = - level_set_tangent_a[0]*level_set_tangent_b[2]
      + level_set_tangent_a[2]*level_set_tangent_b[0];
    level_set_normal[2] =   level_set_tangent_a[0]*level_set_tangent_b[1]
      - level_set_tangent_a[1]*level_set_tangent_b[0];
    double norm = std::sqrt(level_set_normal[0]*level_set_normal[0] +
                            level_set_normal[1]*level_set_normal[1] +
                            level_set_normal[2]*level_set_normal[2]);
    for(int I=0; I < nSpace; I++)
      level_set_normal[I] /=  norm;
  }

  template<int nP>
  inline void _calculate_polynomial_1D(double* xi, double* C_H, double* C_ImH, double* C_D, double& _H, double& _ImH, double& _D)
  {
    _H   = 0.0;
    _ImH = 0.0;
    _D   = 0.0;
    unsigned int iDOF    = 0;
    register double x_pow_i=1.0;
    for(int i=0; i < nP+1; i++, iDOF++)
      {
        register double psi=x_pow_i;
        _H   += C_H[iDOF]*psi;
        _ImH += C_ImH[iDOF]*psi;
        _D   += C_D[iDOF]*psi;
        x_pow_i *= xi[0];
      }
  }

  template<int nP>
  inline void _calculate_polynomial_2D(double* xi, double* C_H, double* C_ImH, double* C_D, double& _H, double& _ImH, double& _D)
  {
    _H   = 0.0;
    _ImH = 0.0;
    _D   = 0.0;
    unsigned int iDOF    = 0;
    register double x_pow_i=1.0;
    for(int i=0; i < nP+1; i++)
      {
        register double y_pow_j=1.0;
        for(int j=0; j < nP+1-i; j++, iDOF++)
          {
            register double psi=x_pow_i*y_pow_j;
            _H   += C_H[iDOF]*psi;
            _ImH += C_ImH[iDOF]*psi;
            _D   += C_D[iDOF]*psi;
            y_pow_j *= xi[1];
          }
        x_pow_i *= xi[0];
      }
  }

  template<int nP>
  inline void _calculate_polynomial_3D(double* xi, double* C_H, double* C_ImH, double* C_D, double& _H, double& _ImH, double& _D)
  {
    _H   = 0.0;
    _ImH = 0.0;
    _D   = 0.0;
    unsigned int iDOF    = 0;
    register double x_pow_i=1.0;
    for(int i=0; i < nP+1; i++)
      {
        register double y_pow_j=1.0;
        for(int j=0; j < nP+1-i; j++)
          {
            register double x_pow_i_y_pow_j=x_pow_i*y_pow_j, z_pow_k=1.0;
            for(int k=0; k < nP+1-i-j; k++, iDOF++)
              {
                register double psi=x_pow_i_y_pow_j*z_pow_k;
                _H   += C_H[iDOF]*psi;
                _ImH += C_ImH[iDOF]*psi;
                _D   += C_D[iDOF]*psi;
                z_pow_k *= xi[2];
              }
            y_pow_j *= xi[1];
          }
        x_pow_i *= xi[0];
      }
  }

  template<int nSpace, int nP, int nQ>
  class Simplex
  {
  public:
    Simplex()
    {
      if (nSpace == 1)
        assert(nDOF == nP+1);
      else if (nSpace == 2)
        assert(nDOF == (nP+1)*(nP+2)/2);
      else if (nSpace == 3)
        assert(nDOF == (nP+1)*(nP+2)*(nP+3)/6);
      else
        assert(false);
      _set_Ainv<nSpace,nP>(Ainv);
    }
    inline void calculate(const double* phi_dof, const double* phi_nodes, const double* xi_r);
    inline void set_quad(unsigned int q)
    {
      assert(q >= 0);
      assert(q < nQ);
      D = _D[q];
      if(inside_out)
        {
          H   = _ImH[q];
          ImH = _H[q];
        }
      else
        {
          H   = _H[q];
          ImH = _ImH[q];
        }
    }
    double H, ImH, D;
    inline double* get_H(){return _H;};
    inline double* get_ImH(){return _ImH;};
    inline double* get_D(){return _D;};
  private:
    static const unsigned int nN=nSpace+1;
    unsigned int root_node, permutation[nN];
    double phi[nN], nodes[nN*3];
    double Jac[nSpace*nSpace], inv_Jac[nSpace*nSpace];
    double level_set_normal[nSpace], X_0[nSpace], phys_nodes_cut[nSpace*3];
    bool inside_out;
    static const unsigned int nDOF=((nSpace-1)/2)*(nSpace-2)*(nP+1)*(nP+2)*(nP+3)/6 + (nSpace-1)*(3-nSpace)*(nP+1)*(nP+2)/2 + (2-nSpace)*((3-nSpace)/2)*(nP+1);
    double Ainv[nDOF*nDOF];
    double C_H[nDOF], C_ImH[nDOF], C_D[nDOF];
    inline int _calculate_permutation(const double* phi_dof, const double* phi_nodes);
    inline void _calculate_cuts();
    inline void _calculate_C();
    double _H[nQ], _ImH[nQ], _D[nQ];
  };

  template<int nSpace, int nP, int nQ>
  inline void Simplex<nSpace,nP,nQ>::_calculate_C()
  {
    register double b_H[nDOF], b_ImH[nDOF], b_dH[nDOF*nSpace];
    _calculate_b<nSpace,nP>(X_0, b_H, b_ImH, b_dH);

    register double Jt_dphi_dx[nSpace];
    for (int I=0; I < nSpace; I++)
      {
        Jt_dphi_dx[I] = 0.0;
        for(int J=0; J < nSpace; J++)
          Jt_dphi_dx[I] += Jac[J*nSpace + I]*level_set_normal[J];
      }
    for (int i=0; i < nDOF; i++)
      {
        C_H[i] = 0.0;
        C_ImH[i] = 0.0;
        C_D[i] = 0.0;
        for (int j=0; j < nDOF; j++)
          {
            C_H[i]   += Ainv[i*nDOF + j]*b_H[j];
            C_ImH[i] += Ainv[i*nDOF + j]*b_ImH[j];
            for (int I=0; I < nSpace; I++)
              {
                C_D[i]   -= Ainv[i*nDOF + j]*b_dH[j*nSpace+I]/Jt_dphi_dx[I];
              }
          }
      }
  }

  const int XX(0),XY(1),XZ(2),YX(3),YY(4),YZ(5),ZX(6),ZY(7),ZZ(8);

  template<int nSpace>
  inline double det(const double* A);

  template<>
  inline double det<1>(const double* A)
  {
    return A[0];
  }

  template<>
  inline double det<2>(const double* A)
  {
    return A[0]*A[3] - A[1]*A[2];
  }

  template<>
  inline double det<3>(const double* A)
  {
    return A[XX]*(A[YY]*A[ZZ] - A[YZ]*A[ZY]) -
      A[XY]*(A[YX]*A[ZZ] - A[YZ]*A[ZX]) +
      A[XZ]*(A[YX]*A[ZY] - A[YY]*A[ZX]);
  }

  template<int nSpace>
  inline double inv(const double* A, double* Ainv);

  template<>
  inline double inv<1>(const double* A, double* Ainv)
  {
    Ainv[0] = 1.0/A[0];
  }

  template<>
  inline double inv<2>(const double* A, double* Ainv)
  {
    double oneOverDet = 1.0/det<2>(A);
    Ainv[0] = oneOverDet*A[3];
    Ainv[1] = -oneOverDet*A[1];
    Ainv[2] = -oneOverDet*A[2];
    Ainv[3] = oneOverDet*A[0];
  }

  template<>
  inline double inv<3>(const double* A, double* Ainv)
  {
    double oneOverDet = 1.0/det<3>(A);
    Ainv[XX] = oneOverDet*(A[YY]*A[ZZ] - A[YZ]*A[ZY]);
    Ainv[YX] = oneOverDet*(A[YZ]*A[ZX] - A[YX]*A[ZZ]);
    Ainv[ZX] = oneOverDet*(A[YX]*A[ZY] - A[YY]*A[ZX]);
    Ainv[XY] = oneOverDet*(A[ZY]*A[XZ] - A[ZZ]*A[XY]);
    Ainv[YY] = oneOverDet*(A[ZZ]*A[XX] - A[ZX]*A[XZ]);
    Ainv[ZY] = oneOverDet*(A[ZX]*A[XY] - A[ZY]*A[XX]);
    Ainv[XZ] = oneOverDet*(A[XY]*A[YZ] - A[XZ]*A[YY]);
    Ainv[YZ] = oneOverDet*(A[XZ]*A[YX] - A[XX]*A[YZ]);
    Ainv[ZZ] = oneOverDet*(A[XX]*A[YY] - A[XY]*A[YX]);
  }

  template<int nSpace, int nP, int nQ>
  inline int Simplex<nSpace,nP,nQ>::_calculate_permutation(const double* phi_dof, const double* phi_nodes)
  {
    int p_i, pcount=0, n_i, ncount=0, z_i, zcount=0;
    root_node=0;
    inside_out=false;
    for (int i=0; i < nN; i++)
      {
        if(phi_dof[i] > 0.0)
          {
            p_i = i;
            pcount  += 1;
          }
        else if(phi_dof[i] < 0.0)
          {
            n_i = i;
            ncount += 1;
          }
        else
          {
            z_i = i;
            zcount += 1;
          }
      }
    if(pcount == nN)
      {
        return 1;
      }
    else if(ncount == nN)
      {
        return -1;
      }
    else if(pcount == 1)
      {
        if (zcount == nN-1)//interface is on an element boundary, only integrate once
          return 1;
        else
          {
            if (nSpace > 1)
              {
                root_node = p_i;
                inside_out = true;
              }
            else
              root_node = n_i;
          }
      }
    else if(ncount == 1)
      {
        root_node = n_i;
      }
    else
      {
        assert(zcount < nN-1);
        if(pcount)
          return 1;
        else if (ncount)
          return -1;
        else
          assert(false);
      }
    for(int i=0; i < nN; i++)
      {
        permutation[i] = (root_node+i)%nN;
      }
    for(int i=0; i < nN; i++)
      {
        phi[i] = phi_dof[permutation[i]];
        for(int I=0; I < 3; I++)
          {
            nodes[i*3 + I] = phi_nodes[permutation[i]*3 + I];//nodes always 3D
          }
      }
    for(int i=0; i < nN - 1; i++)
      for(int I=0; I < nSpace; I++)
        Jac[I*nSpace+i] = nodes[(1+i)*3 + I] - nodes[I];
    double det_Jac = det<nSpace>(Jac);
    if(det_Jac < 0.0)
      {
        double tmp = permutation[nN-1];
        permutation[nN-1] = permutation[nN-2];
        permutation[nN-2] = tmp;
        for(int i=0; i < nN; i++)
          {
            phi[i] = phi_dof[permutation[i]];
            for(int I=0; I < 3; I++)
              {
                nodes[i*3 + I] = phi_nodes[permutation[i]*3 + I];//nodes always 3D
              }
          }
        for(int i=0; i < nN-1; i++)
          for(int I=0; I < nSpace; I++)
            Jac[I*nSpace+i] = nodes[(1+i)*3 + I] - nodes[I];
        det_Jac = det<nSpace>(Jac);
        assert(det_Jac > 0);
      }
    inv<nSpace>(Jac, inv_Jac);
    return 0;
  }

  template<int nSpace, int nP, int nQ>
  inline void Simplex<nSpace,nP,nQ>::_calculate_cuts()
  {
    for (int i=0; i < nN-1;i++)
      {
        if(phi[i+1]*phi[0] < 0.0)
          {
            X_0[i] = 0.5 - 0.5*(phi[i+1] + phi[0])/(phi[i+1]-phi[0]);
            assert(X_0[i] <=1.0);
            assert(X_0[i] >=0.0);
            for (int I=0; I < 3; I++)
              {
                phys_nodes_cut[i*3 + I] = (1-X_0[i])*nodes[I] + X_0[i]*nodes[(i+1)*3+I];
              }
          }
        else
          {
            assert(phi[i+1] == 0.0);
            X_0[i] = 1.0;
            for (int I=0; I < 3; I++)
              {
                phys_nodes_cut[i*3 + I] = nodes[(i+1)*3+I];
              }
          }
      }
  }

  template<int nSpace, int nP, int nQ>
  inline void Simplex<nSpace,nP,nQ>::calculate(const double* phi_dof, const double* phi_nodes, const double* xi_r)
  {
    int icase = _calculate_permutation(phi_dof, phi_nodes);//permuation, Jac,inv_Jac...
    if(icase == 1)
      {
        for (int q=0; q < nQ; q++)
          {
            _H[q] = 1.0;
            _ImH[q] = 0.0;
            _D[q] = 0.0;
          }
        return;
      }
    else if(icase == -1)
      {
        for (int q=0; q < nQ; q++)
          {
            _H[q] = 0.0;
            _ImH[q] = 1.0;
            _D[q] = 0.0;
          }
        return;
      }
    _calculate_cuts();//X_0, array of interface cuts on reference simplex
    _calculate_normal<nSpace>(phys_nodes_cut, level_set_normal);//normal to interface
    _calculate_C();//coefficients of equiv poly
    //compute the default affine map based on phi_nodes[0]
    double Jac_0[nSpace*nSpace];
    for(int i=0; i < nN - 1; i++)
      for(int I=0; I < nSpace; I++)
        Jac_0[i*nSpace+I] = phi_nodes[(1+i)*3 + I] - phi_nodes[I];
    for(int q=0; q < nQ; q++)
      {
        //Due to the permutation, the quadrature points on the reference may be rotated
        //map reference to physical simplex, then back to permuted reference
        register double x[nSpace], xi[nSpace];
        for (int I=0; I < nSpace; I++)
          {
            x[I]=0.0;
            for (int J=0; J < nSpace;J++)
              {
                x[I] += Jac_0[I*nSpace + J]*xi_r[q*3 + J];
              }
            x[I] += phi_nodes[I] - nodes[I];//translate back to reference about nodes[0]
          }
        for (int I=0; I < nSpace; I++)
          {
            xi[I] = 0.0;
            for (int J=0; J < nSpace; J++)
              xi[I] += inv_Jac[I*nSpace + J]*x[J];
          }
        if (nSpace == 1)
          _calculate_polynomial_1D<nP>(xi,C_H,C_ImH,C_D,_H[q],_ImH[q],_D[q]);
        else if (nSpace == 2)
          _calculate_polynomial_2D<nP>(xi,C_H,C_ImH,C_D,_H[q],_ImH[q],_D[q]);
        else if (nSpace == 3)
          _calculate_polynomial_3D<nP>(xi,C_H,C_ImH,C_D,_H[q],_ImH[q],_D[q]);
      }
    set_quad(0);
  }
}//equivalent_polynomials

#endif

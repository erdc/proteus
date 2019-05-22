#ifndef EQUIVALENT_POYNOMIALS_UTILS
#define EQUIVALENT_POYNOMIALS_UTILS
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
  inline void inv(const double* A, double* Ainv);

  template<>
  inline void inv<1>(const double* A, double* Ainv)
  {
    Ainv[0] = 1.0/A[0];
  }

  template<>
  inline void inv<2>(const double* A, double* Ainv)
  {
    double oneOverDet = 1.0/det<2>(A);
    Ainv[0] = oneOverDet*A[3];
    Ainv[1] = -oneOverDet*A[1];
    Ainv[2] = -oneOverDet*A[2];
    Ainv[3] = oneOverDet*A[0];
  }

  template<>
  inline void inv<3>(const double* A, double* Ainv)
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
}//equivalent_polynomials
#endif

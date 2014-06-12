#include <apf.h>
#include <apfVector.h>
#include <apfMatrix.h>

using namespace apf;
#ifndef _H_EIGEN
#define _H_EIGEN

  long eigen (Matrix3x3, Matrix3x3&, Vector3& ,int checkOrth=0);
  double trace (Matrix3x3);
  double trace2 (Matrix3x3);
  double det (Matrix3x3);
  long FindCubicRoots(const double[4], Vector3&);
  long FindQuadraticRoots(const double, const double, Vector3&);
  long NullSpace(const double *, double *, double, long);
  int checkUnitaryOthoganal(Matrix3x3 e,int &factor);
  double dotProd(Vector3 ,Vector3);
  void crossProd(Vector3,Vector3,Vector3&);
  void diffVt(Vector3,Vector3,Vector3);
  int normVt(Vector3,Vector3&);
#endif


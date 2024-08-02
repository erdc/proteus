#!/bin/env python
import math
import sympy
from sympy.integrals import intpoly
import numpy as np
from sympy.abc import I, J, K, N
from sympy import Sum
from sympy.geometry import Triangle, Point
from sympy.printing.cxxcode import cxxcode
x,y,z,x_0,y_0,z_0 = sympy.symbols('x,y,z,x_0,y_0,z_0')

f=open("equivalent_polynomials_coefficients.h",'w')
f.write("""#ifndef EQUIVALENT_POLYNOMIALS_COEFFICIENTS_H
#define EQUIVALENT_POLYNOMIALS_COEFFICIENTS_H

namespace equivalent_polynomials
{
  template<int nSpace, int nP>
  inline void _set_Ainv(double* Ainv);

  template<int nSpace, int nP>
  inline void _calculate_b(double* X_0, double* b_H, double* b_ImH, double* b_dH);

""")

for nSpace in range(1,4):
    for order in range(1,5):
        nDOF1D=order + 1
        if nSpace == 3:
            basis = [x**i*y**j*z**k for i in range(nDOF1D) for j in range(nDOF1D-i) for k in range(nDOF1D-i-j)]
            nDOF = len(basis)
            test_nDOF = sympy.factor(Sum(Sum(Sum(1,(K,0,N-I-J-1)),(J,0,N-I-1)),(I,0,N-1)).doit())
            unit_tet = [[(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)],
                        [1, 2, 3], [2, 3, 0], [3, 0, 1], [0, 1, 2]]
        elif nSpace == 2:
            basis = [x**i*y**j for i in range(nDOF1D) for j in range(nDOF1D-i)]
            nDOF = len(basis)
            test_nDOF = sympy.factor(Sum(Sum(1,(J,0,N-I-1)),(I,0,N-1)).doit())
            unit_triangle = Triangle(Point(0, 0), Point(0, 1), Point(1, 0))
        elif nSpace == 1:
            basis = [x**i for i in range(nDOF1D)]
            nDOF = len(basis)
            test_nDOF = sympy.factor(Sum(1,(I,0,N-1)).doit())
        assert(test_nDOF.evalf(subs={'N':nDOF1D}) == nDOF)
        A = np.zeros((nDOF,nDOF),'d')
        for i in range(nDOF):
            for j in range(nDOF):
                if nSpace == 3:
                    A[i,j] = intpoly.polytope_integrate(unit_tet,sympy.expand(basis[j]*basis[i])).evalf()
                elif nSpace == 2:
                    A[i,j] = intpoly.polytope_integrate(unit_triangle,sympy.expand(basis[j]*basis[i])).evalf()
                elif nSpace == 1:
                    A[i,j] = sympy.integrate(basis[j]*basis[i],(x,0,1)).evalf()
        Ainv = np.linalg.inv(A)
        b_H=[]
        b_dH_x=[]
        b_dH_y=[]
        b_dH_z=[]
        b_1mH=[]
        if nSpace == 3:
            sub_tet = [[(0, 0, 0), (x_0, 0, 0), (0, y_0, 0), (0, 0, z_0)],
                       [1, 2, 3], [2, 3, 0], [3, 0, 1], [0, 1, 2]]
            for i in range(nDOF):
                b_1mH.append(sympy.simplify(intpoly.polytope_integrate(sub_tet, basis[i])))
                b_H.append(sympy.simplify(intpoly.polytope_integrate(unit_tet, basis[i]) 
                                          - sympy.simplify(intpoly.polytope_integrate(sub_tet, basis[i]))))
                b_dH_x.append(sympy.simplify(b_H[i].diff(x_0)))
                b_dH_y.append(sympy.simplify(b_H[i].diff(y_0)))
                b_dH_z.append(sympy.simplify(b_H[i].diff(z_0)))
        elif nSpace == 2:
            sub_triangle = Triangle(Point(0, 0), Point(0, y_0), Point(x_0, 0))
            for i in range(nDOF):
                b_1mH.append(intpoly.polytope_integrate(sub_triangle, basis[i]))
                b_H.append(intpoly.polytope_integrate(unit_triangle, basis[i]) 
                           -intpoly.polytope_integrate(sub_triangle, basis[i]))
                b_dH_x.append(sympy.simplify(b_H[i].diff(x_0)))
                b_dH_y.append(sympy.simplify(b_H[i].diff(y_0)))
        elif nSpace == 1:
            for i in range(nDOF):
                b_1mH.append(sympy.integrate(basis[i],(x,0,x_0)))
                b_H.append(sympy.integrate(basis[i],(x,x_0,1)))
                b_dH_x.append(b_H[i].diff(x_0))

        f.write("""  template<>
  inline void _set_Ainv<{0:d},{1:d}>(double* Ainv)
  {{
""".format(nSpace, order))
        for i in range(nDOF):
            for j in range(nDOF):
                f.write("    Ainv[{0:d}] = {1};\n".format(i*len(b_H)+j,repr(Ainv[i,j])))
        f.write("""  }}

  template<>
  inline void _calculate_b<{0:d},{1:d}>(double* X_0, double* b_H, double* b_ImH, double* b_dH)
  {{
""".format(nSpace, order))
        if nSpace == 3:
            f.write("""    const double x_0(X_0[0]), y_0(X_0[1]), z_0(X_0[2]);
""")
        if nSpace == 2:
            f.write("""    const double x_0(X_0[0]), y_0(X_0[1]);
""")
        if nSpace == 1:
            f.write("""    const double x_0(X_0[0]);
""")
        for i in range(len(b_H)):
            f.write("    b_H[{0:d}] = {1:s};\n".format(i,cxxcode(b_H[i])))
            f.write("    b_ImH[{0:d}] = {1:s};\n".format(i,cxxcode(b_1mH[i])))
            f.write("    b_dH[{0:d}] = {1:s};\n".format(i*nSpace,cxxcode(b_dH_x[i])))
            if nSpace > 1:
                f.write("    b_dH[{0:d}] = {1:s};\n".format(i*nSpace+1,cxxcode(b_dH_y[i])))
            if nSpace > 2:
                f.write("    b_dH[{0:d}] = {1:s};\n".format(i*nSpace+2,cxxcode(b_dH_z[i])))
        f.write("""  }

""")
f.write("""}//equivalent_polynomials
#endif
""")
f.close()
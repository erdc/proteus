#!/bin/env python
import math
import sympy
from sympy.integrals import intpoly
import numpy as np
from sympy.abc import I, J, K, N
from sympy import Sum
from sympy.geometry import Triangle, Point
from sympy.matrices import Matrix
from sympy.printing.cxxcode import cxxcode
x,y,z = sympy.symbols('x,y,z')
n01_x,n01_y,n01_z = sympy.symbols('n01_x,n01_y,n01_z')
n02_x,n02_y,n02_z = sympy.symbols('n02_x,n02_y,n02_z')
n31_x,n31_y,n31_z = sympy.symbols('n31_x,n31_y,n31_z')
n32_x,n32_y,n32_z = sympy.symbols('n32_x,n32_y,n32_z')

f=open("equivalent_polynomials_coefficients_quad.h",'w')
f.write("""#ifndef EQUIVALENT_POLYNOMIALS_COEFFICIENTS_QUAD_H
#define EQUIVALENT_POLYNOMIALS_COEFFICIENTS_QUAD_H

namespace equivalent_polynomials
{
  template<int nP>
  inline void _calculate_b(double theta01, double theta02, double theta31, double theta32, 
                           double phi0, double phi1, double phi2, double phi3, 
                           double* b_H, double* b_ImH, double* b_D);

""")
use_simplify=True

def simplify(expression):
    if use_simplify:
        return sympy.simplify(expression)
    else:
        return expression
nSpace=3
maxOrder = 5
nDOF1D=maxOrder + 1
unit_tet = [[(0, 0, 0),
             (1, 0, 0),
             (0, 1, 0),
             (0, 0, 1)],
            [1, 2, 3],
            [2, 3, 0],
            [3, 0, 1],
            [0, 1, 2]]
n0=(0,0,0)
n1=(1,0,0)
n2=(0,1,0)
n3=(0,0,1)
theta01,theta02,theta31,theta32 = sympy.symbols('theta01,theta02,theta31,theta32')
phi0,phi1,phi2,phi3 = sympy.symbols('phi0,phi1,phi2,phi3')
n01 = (theta01,0,0)
n02 = (0,theta02,0)
n31 = (theta31, 0,        1-theta31)
n32 = (0,        theta32, 1-theta32)
sub_tet0 = [[n1,#0
             n01,#1
             n02,#2
             n31],#3
            [3, 2, 1], #n31, n02, n01
            [2, 3, 0], #n02, n31, n1
            [0, 3, 1], #n1, n31, n01
            [0, 1, 2]] #n1, n01, n02
sub_tet1 = [[n2,#0
             n02,#1
             n32,#2
             n31],#3
            [1, 3, 2], #n02, n31, n32
            [0, 2, 3], #n2, n32, n31
            [3, 1, 0], #n31, n02, n2
            [0, 1, 2]] #n2, n02, n32
sub_tet2 = [[n1,#0
             n2,#1
             n02,#2
             n31],#3
            [2, 3, 1], #n02, n31, n2
            [2, 0, 3], #n02, n1, n31
            [3, 0, 1], #n31, n1, n2
            [1, 0, 2]] #n2, n1, n02

b_H={}
b_1mH={}
b_D={}
for i in range(nDOF1D):
    for j in range(nDOF1D-i):
        for k in range(nDOF1D-i-j):
            basis_ijk = x**i*y**j*z**k
            print("basis function ", i, j, k)
            T0 = simplify(intpoly.polytope_integrate(sub_tet0, basis_ijk))
            T1 = simplify(intpoly.polytope_integrate(sub_tet1, basis_ijk))
            T2 = simplify(intpoly.polytope_integrate(sub_tet2, basis_ijk))
            frag = simplify(T0 + T1 + T2)
            b_H[(i,j,k)] = frag 
            b_1mH[(i,j,k)] = intpoly.polytope_integrate(unit_tet, basis_ijk) - frag
            theta01_expr = 0.5 - 0.5*(phi1 + phi0)/(phi1-phi0)
            theta02_expr = 0.5 - 0.5*(phi2 + phi0)/(phi2-phi0)
            theta31_expr = 0.5 - 0.5*(phi1 + phi3)/(phi1-phi3)
            theta32_expr = 0.5 - 0.5*(phi2 + phi3)/(phi2-phi3)
            b_D[(i,j,k)] = frag.diff(theta01,simplify=use_simplify)*(theta01_expr.diff(phi0,simplify=False) + theta01_expr.diff(phi1,simplify=False)) + \
                           frag.diff(theta02,simplify=use_simplify)*(theta02_expr.diff(phi0,simplify=False) + theta02_expr.diff(phi2,simplify=False)) + \
                           frag.diff(theta31,simplify=use_simplify)*(theta31_expr.diff(phi3,simplify=False) + theta31_expr.diff(phi1,simplify=False)) + \
                           frag.diff(theta32,simplify=use_simplify)*(theta32_expr.diff(phi3,simplify=False) + theta32_expr.diff(phi2,simplify=False))

for order in range(1,maxOrder):
    nDOF1D=order + 1
    print("Order ", order)
    f.write("""  template<>
inline void _calculate_b<{0:d}>(double theta01, double theta02, double theta31, double theta32,
                                double phi0, double phi1, double phi2, double phi3,
                                double* b_H, double* b_ImH, double* b_D)
  {{
""".format(order))
    n=0
    for i in range(nDOF1D):
        for j in range(nDOF1D-i):
            for k in range(nDOF1D-i-j):
                f.write("    b_H[{0:d}] = {1:s};\n".format(n,cxxcode(b_H[(i,j,k)])))
                f.write("    b_ImH[{0:d}] = {1:s};\n".format(n,cxxcode(b_1mH[(i,j,k)])))
                f.write("    b_D[{0:d}] = {1:s};\n".format(n,cxxcode(b_D[(i,j,k)])))
                n+=1
    f.write("""  }

""")

f.write("""}//equivalent_polynomials
#endif
""")
f.close()
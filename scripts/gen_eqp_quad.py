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
                           double* n01, double* n02, double* n31, double* n32, 
                           double* b_H, double* b_ImH, double* b_D);

""")
use_simplify=True

def simplify(expression):
    if use_simplify:
        return sympy.simplify(expression)
    else:
        return expression

for nSpace in range(3,4):
    for order in range(1,5):
        nDOF1D=order + 1
        if nSpace == 3:
            basis = [x**i*y**j*z**k for i in range(nDOF1D) for j in range(nDOF1D-i) for k in range(nDOF1D-i-j)]
            nDOF = len(basis)
            test_nDOF = sympy.factor(Sum(Sum(Sum(1,(K,0,N-I-J-1)),(J,0,N-I-1)),(I,0,N-1)).doit())
            unit_tet = [[(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)],
                        [1, 2, 3], [2, 3, 0], [3, 0, 1], [0, 1, 2]]
        assert(test_nDOF.evalf(subs={'N':nDOF1D}) == nDOF)
        b_H=[]
        b_1mH=[]
        b_D=[]
        if nSpace == 3:
            n0=(0,0,0)
            n1=(1,0,0)
            n2=(0,1,0)
            n3=(0,0,1)
            theta01,theta02,theta31,theta32 = sympy.symbols('theta01,theta02,theta31,theta32')
            n01 = (theta01,0,0)
            n02 = (0,theta02,0)
            n31 = (theta31, 0,        1-theta31)
            n32 = (0,        theta32, 1-theta32)
            sub_tet0 = [[n1,#0
                         n01,#1
                         n02,#2
                         n31],#3
                        [3, 2, 1], 
                        [2, 3, 0], 
                        [0, 3, 1], 
                        [0, 2, 1]]
            sub_tet1 = [[n2,#0
                         n02,#1
                         n32,#2
                         n31],#3
                        [1, 3, 2], 
                        [0, 2, 3], 
                        [3, 1, 0], 
                        [0, 1, 2]]
            sub_tet2 = [[n1,#0
                         n2,#1
                         n02,#2
                         n31],#3
                        [1, 3, 2], 
                        [2, 0, 3], 
                        [3, 0, 1], 
                        [1, 0, 2]]

            interface = Triangle(Point(0, 0), Point(0, 1), Point(1, 0))
            n01_x,n01_y,n01_z=sympy.symbols('n01_x,n01_y,n01_z')
            n02_x,n02_y,n02_z=sympy.symbols('n02_x,n02_y,n02_z')
            n31_x,n31_y,n31_z=sympy.symbols('n31_x,n31_y,n31_z')
            n32_x,n32_y,n32_z=sympy.symbols('n32_x,n32_y,n32_z')
            v01 = Matrix([n01_x,n01_y,n01_z])
            v02 = Matrix([n02_x,n02_y,n02_z])
            v31 = Matrix([n31_x,n31_y,n31_z])
            v32 = Matrix([n32_x,n32_y,n32_z])
            
            normalDir1 = (v31-v01).cross(v32-v01)
            J1=normalDir1.norm()
            nx1 = normalDir1[0]
            ny1 = normalDir1[1]
            nz1 = normalDir1[2]
            phi1 = (nx1*x + ny1*y + nz1*z - v01.dot(normalDir1))/J1
            X1 = v31[0]*x + v32[0]*y + v01[0]*(1-x-y)
            Y1 = v31[1]*x + v32[1]*y + v01[1]*(1-x-y)
            Z1 = v31[2]*x + v32[2]*y + v01[2]*(1-x-y)
            
            normalDir2 = (v32-v01).cross(v02-v01)
            J2=normalDir2.norm()
            nx2 = normalDir2[0]
            ny2 = normalDir2[1]
            nz2 = normalDir2[2]
            phi2 = (nx2*x + ny2*y + nz2*z - v01.dot(normalDir2))/J2
            X2 = v32[0]*x + v02[0]*y + v01[0]*(1-x-y)
            Y2 = v32[1]*x + v02[1]*y + v01[1]*(1-x-y)
            Z2 = v32[2]*x + v02[2]*y + v01[2]*(1-x-y)
            for i in range(nDOF):
                frag = simplify(simplify(intpoly.polytope_integrate(sub_tet0, basis[i])) +
                                simplify(intpoly.polytope_integrate(sub_tet1, basis[i])) +
                                simplify(intpoly.polytope_integrate(sub_tet2, basis[i])))
                b_1mH.append(frag)
                b_H.append(simplify(intpoly.polytope_integrate(unit_tet, basis[i])
                                    - frag))
                subtris = simplify(J1*intpoly.polytope_integrate(interface,sympy.expand(basis[i].subs([(x,X1),
                                                                                                       (y,Y1),
                                                                                                       (z,Z1)])))
                                   +
                                   J2*intpoly.polytope_integrate(interface,sympy.expand(basis[i].subs([(x,X2),
                                                                                                       (y,Y2),
                                                                                                       (z,Z2)]))))
                b_D.append(subtris)

        f.write("""  template<>
inline void _calculate_b<{1:d}>(double theta01, double theta02, double theta31, double theta32, 
                                double* n01, double* n02, double* n31, double* n32, 
                                double* b_H, double* b_ImH, double* b_D)
  {{
""".format(nSpace, order))
        if nSpace == 3:
            f.write("""    const double n01_x(n01[0]),n01_y(n01[1]),n01_z(n01[2]);
    const double n02_x(n02[0]),n02_y(n02[1]),n02_z(n02[2]);
    const double n31_x(n31[0]),n31_y(n31[1]),n31_z(n31[2]);
    const double n32_x(n32[0]),n32_y(n32[1]),n32_z(n32[2]);
""")
        for i in range(len(b_H)):
            f.write("    b_H[{0:d}] = {1:s};\n".format(i,cxxcode(b_H[i])))
            f.write("    b_ImH[{0:d}] = {1:s};\n".format(i,cxxcode(b_1mH[i])))
            f.write("    b_D[{0:d}] = {1:s};\n".format(i,cxxcode(b_D[i])))
        f.write("""  }

""")
f.write("""}//equivalent_polynomials
#endif
""")
f.close()

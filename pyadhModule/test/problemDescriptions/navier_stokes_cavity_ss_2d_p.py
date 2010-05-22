from pyadh import *
from pyadh.default_p import *
##\page Tests Test Problems
#\ref navier_stokes_cavity_ss_2d_p.py "Incompressible Navier-Stokes in a driven cavity at steady state"
#

##\ingroup test
#\file navier_stokes_cavity_ss_2d_p.py
#
#\brief Incompressible Navier-Stokes flow in a driven cavity at steady state
#
#The incompressible Navier-Stokes equations are described by the pyadh::TransportCoefficients::NavierStokes class.
#
#The boundary conditions are given by
#
#\f{eqnarray*}
#p(x=0.5,y=0,t) &=& 0 \\
#u(x,y,t) &=& \left\{ \begin{array}{rl}
#1.0e-5 & 0\leq x < 0.5;\quad y=0,1 \\
#-1.0e-5 & 0.5 \leq x \leq 1.0;\quad y=0,1 \\
#0&x=0,1
#\end{array} \right. \\
#v(x,y,t) &=& \left\{ \begin{array}{rl}
#-1.0e-5 & x=0,1; \quad 0\leq y \leq 0.5 \\
#1.0e-5 & x=0,1; \quad 0.5 < y \leq 1.0 \\
#0&y=0,1
#\end{array} \right. \\
#\f}
#
#\image html navier_stokes_cavity_ss_v.jpg
#\image html navier_stokes_cavity_ss_p.jpg
nd = 2

initialConditions = None

analyticalSolution = None

coefficients = NavierStokes(g=[0.0,0.0],nd=nd)

beltSpeed = 1.0e-5

def getDBC_pressure(x,flag):
    if x[1] == 0.0:
        if x[0] == 0.5:
            return lambda x,t: 0.0

def getDBC_u(x,flag):
    #set x velocity to 1 at ends and to 0 at top and bottom
    if x[0] == 0.0:
        return lambda x,t: 0.0
    elif x[0] == 1.0:
        return lambda x,t: 0.0
    elif x[1] == 0.0:
        if x[0] < 0.5:
            return lambda x,t: beltSpeed
        else:
            return lambda x,t: -beltSpeed
    elif x[1] == 1.0:
        if x[0] < 0.5:
            return lambda x,t: beltSpeed
        else:
            return lambda x,t: -beltSpeed

def getDBC_v(x,flag):
    if x[1] == 0.0:
        return lambda x,t: 0.0
    elif x[1] == 1.0:
        return lambda x,t: 0.0
    elif x[0] == 0.0:
        if x[1] > 0.5:
            return lambda x,t: beltSpeed
        else:
            return lambda x,t: -beltSpeed
    elif x[0] == 1.0:
        if x[1] > 0.5:
            return lambda x,t: beltSpeed
        else:
            return lambda x,t: -beltSpeed

dirichletConditions = {0:getDBC_pressure,
                       1:getDBC_u,
                       2:getDBC_v}

fluxBoundaryConditions = {0:'noFlow',
                          1:'noFlow',
                          2:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}


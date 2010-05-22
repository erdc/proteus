from pyadh import *
from pyadh.default_p import *
"""
A driven cavity problem for Stokes flow.
"""
##\page Tests Test Problems 
# \ref stokes_cavity_2d_p.py "Stokes flow in a driven cavity"
#

##\ingroup test
#\file stokes_cavity_2d_p.py
#\brief A driven cavity problem for Stokes flow.
#
#The model equations are defined by the src::TransportCoefficients::Stokes class.
#The initial/boundary conditions are given by
#\f{eqnarray*}
#\mathvf{v}(\mathbf{x},0) &=& 0 \\
#v_x(\mathbf{x},t) &=& \left\{ \begin{array}{lr}
#0 & (x,y) \in \{0,1\} \times (0,1) \\
#0.1 & (x,y) \in (0,1/2) \times \{0\} \\
#-0.1 & (x,y) \in (1/2,1) \times \{0\} \\
#0.1 & (x,y) \in (0,1/2) \times \{1\} \\
#-0.1 & (x,y) \in (1/2,1) \times \{1\} \\
#\end{array} \right. \\
#v_y(\mathbf{x},t) &=& \left\{ \begin{array}{lr}
#0 & (x,y) \in (0,1) \times \{0,1\} \\
#0.1 & (x,y) \in \{0\} \times (0,1/2)  \\
#-0.1 & (x,y) \in  \{0\} \times (1/2,1)  \\
#0.1 & (x,y) \in \{1\} \times (0,1/2)  \\
#-0.1 & (x,y) \in  \{1\} \times (1/2,1)  \\
#\end{array} \right. \\
#\f}
nd = 2

L=(1.0,1.0,1.0)

analyticalSolution = None

coefficients = Stokes(g=[0.0,0.0],steady=False)

#now define the Dirichlet boundary conditions

def getDBC_pressure(x):
    if x[1] == 0.0:
        if x[0] == 0.5:
            return lambda x,t: 0.0

def getDBC_u(x):
    #set x velocity to 1 at ends and to 0 at top and bottom
    if x[0] == 0.0:
        return lambda x,t: 0.0
    elif x[0] == 1.0:
        return lambda x,t: 0.0
    elif x[1] == 0.0:
        if x[0] < 0.5:
            return lambda x,t: 0.1
        else:
            return lambda x,t: -0.1
    elif x[1] == 1.0:
        if x[0] < 0.5:
            return lambda x,t: 0.1
        else:
            return lambda x,t: -0.1

def getDBC_v(x):
    if x[1] == 0.0:
        return lambda x,t: 0.0
    elif x[1] == 1.0:
        return lambda x,t: 0.0
    elif x[0] == 0.0:
        if x[1] > 0.5:
            return lambda x,t: 0.1
        else:
            return lambda x,t: -0.1
    elif x[0] == 1.0:
        if x[1] > 0.5:
            return lambda x,t: 0.1
        else:
            return lambda x,t: -0.1

dirichletConditions = {0:getDBC_pressure,
                       1:getDBC_u,
                       2:getDBC_v}

class ZeroIC:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions  = {0:ZeroIC(),1:ZeroIC(),2:ZeroIC()}

fluxBoundaryConditions = {0:'noFlow',
                          1:'noFlow',
                          2:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

T = 1.0e4

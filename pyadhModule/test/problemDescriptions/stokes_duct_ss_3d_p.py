from pyadh import *
from pyadh.default_p import *
"""
Stokes flow in a duct with square cross-section.
"""
##\page Tests Test Problems 
# \ref stokes_duct_ss_3d_p.py "Stokes flow in a square duct at steady state"
#

##\ingroup test
#\file stokes_duct_ss_3d_p.py
#
#\brief Stokes flow in a duct with square cross-section.
#
#The model equations are defined by the Stokes class. The
#boundary conditions are
#\f{eqnarray*}
#p(1,y,z,t) &=& (1-z) \rho  g \\
#v_x(0,y,z,t) &=& 0.1 \\
#v_x(x,y,z,t) &=& 0.0 \quad y,z \in \{0,1\}\\
#\frac{\partial v_x}{\partial x}(1,y,z,t) &=& 0 \\
#v_y(x,y,z,t) &=& 0.0 \mbox{ on } \partial \Omega\\
#v_z(x,y,z,t) &=& 0.0 \mbox{ on } \partial \Omega\\
#\f}
nd = 3
L=(1.0,1.0,1.0)

initialConditions = None

analyticalSolution = None

coefficients = Stokes(g=[0.0,0.0,0.0],nd=nd,steady=True)
#coefficients = NavierStokes(g=[0.0,0.0,0.0],nd=nd)

#now define the Dirichlet boundary conditions
inflow = 0.0001

def getDBC_pressure_duct(x):
    if x[0] == L[0]:
        return lambda x,t: (x[2]-L[2])*coefficients.rho*coefficients.g[2]
    else:
        pass
#     #lid driven cavity
#     if (x[0] == L[0]
#         and x[1] == L[1]
#         and x[2] == L[2]):
#         return lambda x,t: 0.0
#     else:
#         pass

def getDBC_u_duct(x):
    if (x[1] == 0.0 or
        x[1] == L[1] or
        x[2] == 0.0 or
        x[2] == L[2]):
        return lambda x,t: 0.0
    else:
        pass
#     #lid driven cavity
#     if (x[0] == 0.0 or
#         x[0] == L[0] or
#         x[1] == 0.0 or
#         x[1] == L[1] or
#         x[2] == 0.0 or
#         x[2] == L[2]):
#         return lambda x,t: 0.0
#     else:
#         pass

def getDBC_v_duct(x):
    if (x[1] == 0.0 or
        x[1] == L[1] or
        x[2] == 0.0 or
        x[2] == L[2]):
        return lambda x,t: 0.0
    else:
        pass
#     #lid driven cavity
#     if (x[0] == L[0] or
#         x[1] == 0.0 or
#         x[1] == L[1] or
#         x[2] == 0.0 or
#         x[2] == L[2]):
#         return lambda x,t: 0.0
#     elif x[0] == 0.0:
#         return lambda x,t: inflow

def getDBC_w_duct(x):
    if (x[1] == 0.0 or
        x[1] == L[1] or
        x[2] == 0.0 or
        x[2] == L[2]):
        return lambda x,t: 0.0
    else:
        pass
#     #lid driven cavity
#     if (x[0] == 0.0 or
#         x[0] == L[0] or
#         x[1] == 0.0 or
#         x[1] == L[1] or
#         x[2] == 0.0 or
#         x[2] == L[2]):
#         return lambda x,t: 0.0
#     else:
#         pass

dirichletConditions = {0:getDBC_pressure_duct,
                       1:getDBC_u_duct,
                       2:getDBC_v_duct,
                       3:getDBC_w_duct}


def getAFBC_p_duct(x):
    if (x[1] == L[1] or
        x[1] == 0.0 or
        x[2] == L[2] or
        x[2] == 0.0):
        return lambda x,t: 0.0
    elif (x[0] == 0.0):
        return lambda x,t: -inflow
    else:
        pass
#     #lid driven cavity
#     if (x[0] == 0.0 or
#         x[0] == L[0] or
#         x[1] == 0.0 or
#         x[1] == L[1] or
#         x[2] == 0.0 or
#         x[2] == L[2]):
#         return lambda x,t: 0.0

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow',
                          3:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_duct}

diffusiveFluxBoundaryConditions = {0:{}}


from pyadh import *
from pyadh.default_p import *
"""
Stokes flow around a sphere
"""
##\page Tests Test Problems 
# \ref stokes_sphere_ss_3d_p.py "Stokes flow around sphere
#

##\ingroup test
#\file stokes_sphere_ss_3d_p.py
#
#\brief Stokes flow in a sphere packing
#
#The model equations are defined by the Stokes class. The
#boundary conditions are
#\f{eqnarray*}
#p(L,y,z,t) &=& (1-z) \rho  g \\
#p(0,y,z,t) &=& 10+(1-z) \rho g \\
#v_x(x,y,z,t) &=& 0.0 \mbox{ on } y=[0,L],z=[0,L]\\
#v_y(x,y,z,t) &=& 0.0 \mbox{ on } y=[0,L],z=[0,L]\\
#v_z(x,y,z,t) &=& 0.0 \mbox{ on } y=[0,L],z=[0,L]\\
#\f}
polyfile = "sphere_in_box"
nd = 3

initialConditions = None

analyticalSolution = None

coefficients = Stokes(g=[0.0,0.0,0.0],nd=nd)
inflow = 0.1
radius = 0.001
nspheres=1
L = (0.008,0.008,0.008)
center = (L[0]*0.5,L[1]*0.5,L[2]*0.5)
#now define the Dirichlet boundary conditions
EPS = 1.0e-4
def onSphere(x):
    dist = sqrt((x[0]-center[0])**2+(x[1]-center[1])**2+(x[2]-center[2])**2)
    return dist <= radius+EPS

def getDBC_pressure_pome(x):
    if x[0] >= L[0]-EPS:
        return lambda x,t: (L[2]-x[2])*abs(coefficients.g[2])*coefficients.rho
    if x[0] <= EPS:
        return lambda x,t: 10.0+(L[2]-x[2])*abs(coefficients.g[2])*coefficients.rho
def getDBC_u_pome(x):
    if (x[1] <= EPS or
        x[1] >= L[1]-EPS or
        x[2] <= EPS or
        x[2] >= L[2]-EPS or
        onSphere(x)):
        return lambda x,t: 0.0
    
def getDBC_v_pome(x):
    if (x[1] <= EPS or
        x[1] >= L[1]-EPS or
        x[2] <= EPS or
        x[2] >= L[2]-EPS or
        onSphere(x)):
        return lambda x,t: 0.0
    
def getDBC_w_pome(x):
    if (x[1] <= EPS or
        x[1] >= L[1]-EPS or
        x[2] <= EPS or
        x[2] >= L[2]-EPS or
        onSphere(x)):
        return lambda x,t: 0.0
dirichletConditions = {0:getDBC_pressure_pome,
                       1:getDBC_u_pome,
                       2:getDBC_v_pome,
                       3:getDBC_w_pome}


def getAFBC_p_pome(x):
    if (x[1] >= L[1]-EPS or
        x[1] <= EPS or
        x[2] >= L[2]-EPS or
        x[2] <= EPS or
        onSphere(x)):
        return lambda x,t: 0.0

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow',
                          3:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_pome}

diffusiveFluxBoundaryConditions = {0:{}}


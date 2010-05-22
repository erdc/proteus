from pyadh import *
from pyadh.default_p import *
"""
Navier Stokes flow around a sphere
"""
##\page Tests Test Problems 
# \ref navier_stokes_sphere_ss_3d_p.py "Navier-Stokes flow around sphere
#

##\ingroup test
#\file navier_stokes_sphere_ss_3d_p.py
#
#\brief Navier_stokes flow in a sphere packing
#
#The model equations are defined by the NavierStokes class. The
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

#coefficients = NavierStokes(g=[0.0,0.0,0.0],nd=nd)
coefficients = TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                             sigma=0.0,
                                             rho_0=998.2,nu_0=1.004e-6,
                                             rho_1=998.2,nu_1=1.004e-6,
                                             g=[0.0,0.0,0.0],
                                             nd=3,
                                             LS_model=None,
                                             KN_model=None,
                                             epsFact_density=None,
                                             stokes=False)
inflow = 0.1
radius = 0.001
inflowPbc = False
inflowP= 1.
nspheres=1
L = (0.008,0.008,0.008)
center = (L[0]*0.5,L[1]*0.5,L[2]*0.5)
#now define the Dirichlet boundary conditions
EPS = 1.0e-6
def onSphere(x):
    dist = sqrt((x[0]-center[0])**2+(x[1]-center[1])**2+(x[2]-center[2])**2)
    return dist <= radius+EPS

T=100.0*L[0]/inflow
tnList = [0.0,0.1*L[0]/inflow,T]

def velRamp(t):
    residence_time = L[0]/inflow
    if t < residence_time:
        return 1.0-exp(-25.0*t/residence_time)
    else:
        return 1.0
def inflowU(t):
    return velRamp(t)*inflow
def getDBC_pressure_pome(x,flag):
    if x[0] >= L[0]-EPS:
        return lambda x,t: (L[2]-x[2])*abs(coefficients.g[2])*coefficients.rho
    if x[0] <= EPS and inflowPbc:
        return lambda x,t: inflowP+(L[2]-x[2])*abs(coefficients.g[2])*coefficients.rho
def getDBC_u_pome(x,flag):
    if (x[1] <= EPS or
        x[1] >= L[1]-EPS or
        x[2] <= EPS or
        x[2] >= L[2]-EPS or
        onSphere(x)):
        return lambda x,t: 0.0
    if x[0] <= EPS and not inflowPbc:
        return lambda x,t: inflowU(t)
def getDBC_v_pome(x,flag):
    if (x[1] <= EPS or
        x[1] >= L[1]-EPS or
        x[2] <= EPS or
        x[2] >= L[2]-EPS or
        onSphere(x)):
        return lambda x,t: 0.0
    if x[0] <= EPS and not inflowPbc:
        return lambda x,t: 0.0 
def getDBC_w_pome(x,flag):
    if (x[1] <= EPS or
        x[1] >= L[1]-EPS or
        x[2] <= EPS or
        x[2] >= L[2]-EPS or
        onSphere(x)):
        return lambda x,t: 0.0
    if x[0] <= EPS and not inflowPbc:
        return lambda x,t: 0.0 
dirichletConditions = {0:getDBC_pressure_pome,
                       1:getDBC_u_pome,
                       2:getDBC_v_pome,
                       3:getDBC_w_pome}


def getAFBC_p_pome(x,flag):
    if (x[1] >= L[1]-EPS or
        x[1] <= EPS or
        x[2] >= L[2]-EPS or
        x[2] <= EPS or
        onSphere(x)):
        return lambda x,t: 0.0
    if x[0] <= EPS and not inflowPbc:
        return lambda x,t: -inflowU(t) 

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow',
                          3:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_pome}

diffusiveFluxBoundaryConditions = {0:{}}

T=1.01e1

from pyadh import *
from pyadh.default_p import *
"""
Incompressible Navier-Stokes flow in a driven cavity.
"""

##\page Tests Test Problems 
# \ref navier_stokes_cavity_2d_p.py "Incompressible Navier-Stokes in a driven cavity."
#

##\ingroup test
#\file navier_stokes_cavity_2d_p.py
#
#\brief Incompressible Navier-Stokes flow in a driven cavity.
#
#The Navier-Stokes equations implemented in the pyadh::TransportCoefficients::NavierStokes class.
#
#The domain is a unit square, and the initial/boundary conditions are
#\f{eqnarray*}
#\mathbf{v}(x,y,t=0) &=& 0.0 \mbox{ on } \Omega\\
#p(x,y,t=0) &=& 0.0 \mbox{ on } \Omega\\
#p(x=0.5,y=0.0,t) &=& 0.0 \\
#u(x,y=1.0,t) &=& 10000.0 \\
#u(x,y=0.0,t) &=& u(0.0,y,t) = u(1.0,y,t) = 0.0 \\
#v(x,y,t) &=& 0.0 \mbox{ on } \partial \Omega\\
#\f}
#
#\image html navier_stokes_cavity_v.jpg
#
# <A href="https://juanita.wes.army.mil/~cekees/pyadh-doc/images/navier_stokes_cavity.avi"> AVI Animation of velocity</A>
#
nd = 3

L=(1.0,0.5,1.0)

#coefficients = NavierStokes(g=[0.0,9.8])
#speed = 1.0e3*coefficients.nu
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
Re = 100.0
speed = Re*coefficients.nu_0
theta=math.pi/4.0
v_x = speed*math.cos(theta)
v_y = speed*math.sin(theta)
print "Re = %12.5e" % (Re,)
#now define the Dirichlet boundary conditions

def getDBC_pressure(x,flag):
    if x[1] == 0.0:
        if x[0] == 0.5:
            return lambda x,t: 0.0

def getDBC_u(x,flag):
    if x[0] in [0.0,L[0]]:
        return lambda x,t: 0.0
    elif x[1] in [0.0,L[1]]:
        return lambda x,t: 0.0
    elif x[2] == 0.0:
        return lambda x,t: 0.0
    elif x[2] == L[2]:
        return lambda x,t: v_x

def getDBC_v(x,flag):
    if x[0] in [0.0,L[0]]:
        return lambda x,t: 0.0
    elif x[1] in [0.0,L[1]]:
        return lambda x,t: 0.0
    elif x[2] == 0.0:
        return lambda x,t: 0.0
    elif x[2] == L[2]:
        return lambda x,t: v_y

def getDBC_w(x,flag):
    if x[0] in [0.0,L[0]]:
        return lambda x,t: 0.0
    elif x[1] in [0.0,L[1]]:
        return lambda x,t: 0.0
    elif x[2] in [0.0,L[2]]:
        return lambda x,t: 0.0

dirichletConditions = {0:getDBC_pressure,
                       1:getDBC_u,
                       2:getDBC_v,
                       3:getDBC_w}

class ZeroIC:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions  = {0:ZeroIC(),1:ZeroIC(),2:ZeroIC(),3:ZeroIC()}

fluxBoundaryConditions = {0:'noFlow',
                          1:'noFlow',
                          2:'noFlow',
                          3:'noFlow'}

def getAFBC_p_tube(x,flag):
    if x[0] in [0.0,L[0]]:
        return lambda x,t: 0.0
    elif x[1] in [0.0,L[1]]:
        return lambda x,t: 0.0
    elif x[2] in [0.0,L[2]]:
        return lambda x,t: 0.0
def getAFBC_u_tube(x,flag):
    pass
def getAFBC_v_tube(x,flag):
    pass
def getAFBC_w_tube(x,flag):
    pass

advectiveFluxBoundaryConditions =  {0:getAFBC_p_tube,
                                    1:getAFBC_u_tube,
                                    2:getAFBC_v_tube,
                                    3:getAFBC_w_tube}

def getDFBC_tube(x,flag):
    pass

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_tube},
                                   2:{2:getDFBC_tube},
                                   3:{3:getDFBC_tube}}

T = 100.0*L[0]/speed

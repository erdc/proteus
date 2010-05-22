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
nd = 2

L=(1.0,1.0,1.0)

coefficients = TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                             sigma=0.0,
                                             rho_0=1.0,nu_0=1.0e-3,
                                             rho_1=1.0,nu_1=1.0e-3,
                                             g=[0.0,-9.8],
                                             nd=2,
                                             LS_model=None,
                                             KN_model=None,
                                             epsFact_density=None,
                                             stokes=True)
Re = 20000.0
speed = Re*coefficients.nu_0
residence_time = 1.0/speed

Re = 10000.0
speed = Re*coefficients.nu_0
residence_time = 0.0

Re = 100.0
speed = Re*coefficients.nu_0
residence_time = 0.0

#speed = 1.0
#residence_time = 0.0

#now define the Dirichlet boundary conditions

def velRamp(t):
    if t < residence_time:
        return 1.0-exp(-25.0*t/residence_time)
    else:
        return 1.0

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
        return lambda x,t: 0.0
    elif x[1] == 1.0:
        return lambda x,t: speed*velRamp(t)

def getDBC_v(x,flag):
    if x[1] == 0.0:
        return lambda x,t: 0.0
    elif x[1] == 1.0:
        return lambda x,t: 0.0
    elif x[0] == 0.0:
        return lambda x,t: 0.0
    elif x[0] == 1.0:
        return lambda x,t: 0.0

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

def getZeroBC(x,flag):
    return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getZeroBC,1:getZeroBC,2:getZeroBC}

diffusiveFluxBoundaryConditions = {0:{},1:{1:getZeroBC},2:{2:getZeroBC}}

T = 1000.0*L[0]/speed


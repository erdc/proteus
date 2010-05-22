from pyadh import *
from pyadh.default_p import *
import step2d
from pyadh import AnalyticalSolutions
"""
Incompressible Navier-Stokes flow around a square obstacle in 2D.
"""
## \page Tests Test Problems
#\ref navier_stokes_squareobstacle_2d_p.py "Incompressible Navier-Stokes flow around a square obstacle"
# 

##\ingroup test
# \file navier_stokes_squareobstacle_2d_p.py
# @{
#
# \brief Incompressible flow around a square with no-slip boundary conditions on the square.
#
# The governing equations are described by the pyadh::TransportCoefficients::NavierStokes class. The domain and initial/boundary conditions are
#
# \f{eqnarray*}
# \Omega &=& \Omega_E \backslash \Omega_O  \\
# \Omega_E &=& [0,20] \times [0,10]   \\
# \Omega_O &=& [2,3] \times [5,6]   \\
# u&=&0 \mbox{ on } \partial \Omega_O   \\
# v&=&0 \mbox{ on } \partial \Omega_O   \\
# u &=& 100 \mbox{ on } \{0\} \times [0,10]    \\
# v &=& 0 \mbox{ on } [0,20] \times \{0,10\}   \\
# \frac{\partial u}{\partial x} &=& \{0\} \mbox{ on } 20 \times [0,10]   \\
# \frac{\partial v}{\partial x} &=& \{0\} \mbox{ on } 20 \times [0,10]   \\
# \f}
#
# The flow generates vortices behind the obstacle, which break off periodically and are carried down stream (a Karman vortex street).
#
#\image html squareobstacle_p.jpg
#\image html squareobstacle_v.jpg
#
# <A href="https://adh.usace.army.mil/pyadh-images/squareobstacle.avi"> AVI Animation of velocity relative to free stream velocity</A>
#
nd = 2

#T = 9.999e0
T = 1.0
upstream_height=5.2
upstream_start =4.9
upstream_length=5
downstream_height=5.2+4.9
downstream_length=150.0
polyfile = "step2d"
boundaryTags = step2d.genPoly(fileprefix=polyfile,
                              upstream_height=upstream_height,
                              downstream_height=downstream_height,
                              upstream_length=upstream_length,
                              downstream_length=downstream_length,
                              step_fun=step2d.linear_profile,
                              n_points_on_step=2,
                              step_length=0.0)

Re=700.0
coefficients = TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                             sigma=0.0,
                                             rho_0=1.0,nu_0=1.0/Re,
                                             rho_1=1.0,nu_1=1.0/Re,
                                             g=[0.0,0.0],
                                             nd=nd,
                                             LS_model=None,
                                             KN_model=None,
                                             epsFact_density=None,
                                             stokes=False)
max_upstream_speed = 3.0*Re*coefficients.nu/(2.0*upstream_height)
uProfile = AnalyticalSolutions.PlanePoiseuilleFlow_u2(plane_theta=math.pi/2.0,
                                                      plane_phi=math.pi/2.0,
                                                      v_theta=0.0,
                                                      v_phi=math.pi/2.0,
                                                      v_norm=0.0,
                                                      mu=coefficients.rho*coefficients.nu,
                                                      grad_p=-6.0*Re*coefficients.nu**2/upstream_height**3,
                                                      L=[1.0,upstream_height,1.0])

def velRamp(t):
    return 1.0
#     if t < 25.0/(tubeEnd/inflow):
#         return 1.0-exp(-t*tubeEnd/inflow)
#     else:
#         return 1.0
def getDBC_pressure_tube(x,flag):
    if flag == boundaryTags['downstream']:
        return lambda x,t: -(downstream_height-x[1])*coefficients.rho*coefficients.g[1]
    else:
        pass

bottom = [boundaryTags['upstream_bottom'],boundaryTags['step_bottom'],boundaryTags['downstream_bottom']]

def getDBC_u_tube(x,flag):
    if flag == boundaryTags['upstream']:
        return lambda x,t: uProfile.uOfX(x-numpy.array([0.0,upstream_start,0.0]))*velRamp(t)
    elif (flag == boundaryTags['top'] or
          (flag in bottom)):
        return lambda x,t: 0.0

def getDBC_v_tube(x,flag):
    if (flag == boundaryTags['top'] or
        (flag in bottom)):
        return lambda x,t: 0.0

dirichletConditions = {0:getDBC_pressure_tube,
                       1:getDBC_u_tube,
                       2:getDBC_v_tube}


# def getPDBC(x):
#    if x[1] < 1.0e-8 or x[1] >= tubeTop - 1.0e-8:
#        print tuple(x)
#        print "link",tuple(numpy.array([x[0],0.0,0.0]))
#        return numpy.array([x[0],0.0,0.0])

# periodicDirichletConditions = {0:getPDBC,1:getPDBC,2:getPDBC}

def getAFBC_p_tube(x,flag):
    #top and bottom
    if (flag == boundaryTags['top'] or
        (flag in bottom)):
        return lambda x,t: 0.0

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_tube}

diffusiveFluxBoundaryConditions = {0:{}}

class Steady_p:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return -(downstream_height-x[1])*coefficients.rho*coefficients.g[1]

class Steady_u:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return uProfile.uOfX(x-numpy.array([0.0,upstream_start,0.0]))

class Steady_v:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:Steady_p(),
                     1:Steady_u(),
                     2:Steady_v()}

## @}

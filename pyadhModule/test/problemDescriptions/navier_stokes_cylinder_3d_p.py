from pyadh import *
from pyadh.default_p import *
import cylinder3d
import sys
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
nd = 3

rho=1.0
nu=1.0e-3
inflow_height=0.41
bottom_width = 0.41
bottom_length=2.5
cylinder_radius=0.1/2.0
cylinder_center = (0.45+0.1/2.0,
                   0.15+0.1/2.0)
#REs = raw_input("Enter Reynolds Number\n")
REs = "1.0"
name = "cylinder3d_"+REs+"_"
RE = float(REs)
Ubar = nu*RE/(2.0*cylinder_radius)
Um = 3.0*Ubar/2.0
T = 10.0*bottom_length/Um
print "REYNOLDS NUMBER = ",Ubar*2.0*cylinder_radius/nu

genMesh=True
polyfile = "cylinder3d"
circularCrossSection = True
if circularCrossSection:
    #circular
    boundaryTags = cylinder3d.genPoly(polyfile,
                                      cross_section=cylinder3d.circular_cross_section,
                                      height=inflow_height,
                                      length=bottom_length,
                                      width=bottom_width,
                                      radius=cylinder_radius,
                                      center=cylinder_center,
                                      n_points_on_obstacle=2*11-2)
else:
    #square
    boundaryTags = cylinder3d.genPoly(polyfile,
                                      cross_section=cylinder3d.circular_cross_section,
                                      height=inflow_height,
                                      length=bottom_length,
                                      width=bottom_width,
                                      radius=math.sqrt(2.0*cylinder_radius**2),
                                      center=cylinder_center,
                                      n_points_on_obstacle=4,
                                      thetaOffset=math.pi/4.0)

analyticalSolution = {}
from pyadh import RANS2P
useOpt=True#False
if useOpt:
    LevelModelType = RANS2P.OneLevelRANS2P
bcsTimeDependent=False#True

coefficients = TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                             sigma=0.0,
                                             rho_0=rho,nu_0=nu,
                                             rho_1=rho,nu_1=nu,
                                             g=[0.0,0.0,0.0],
                                             nd=nd,
                                             LS_model=None,
                                             KN_model=None,
                                             epsFact_density=1.5,#None,
                                             stokes=False)
mu = coefficients.rho_0*coefficients.nu_0
#coefficients = NavierStokes(g=[0.0,0.0],nd=nd)
#coefficients = Stokes(g=[0.0,-9.8],nd=nd)

#now define the Dirichlet boundary conditions
residence_time = bottom_length/Um
bcsTimeDependent=False
def velRamp(t):
    if t < residence_time:
        return 1.0-exp(-25.0*t/residence_time)
    else:
        return 1.0

def getDBC_pressure_tube(x,flag):
    if flag == boundaryTags['downstream']:
        return lambda x,t: -(inflow_height-x[1])*coefficients.rho*coefficients.g[1]
    else:
        pass

noSlip = [boundaryTags['obstacle'],boundaryTags['bottom'],boundaryTags['top'],boundaryTags['front'],boundaryTags['back']]
          
def getDBC_u_tube(x,flag):
    if flag == boundaryTags['upstream']:
        return lambda x,t: Um#velRamp(t)*16.0*Um*x[1]*x[2]*(bottom_width-x[1])*(inflow_height-x[2])/(inflow_height**4)
    elif flag in noSlip:
        return lambda x,t: 0.0

def getDBC_v_tube(x,flag):
    if flag == boundaryTags['upstream']:
        return lambda x,t: 0.0
    elif flag in noSlip:
        return lambda x,t: 0.0

def getDBC_w_tube(x,flag):
    if flag == boundaryTags['upstream']:
        return lambda x,t: 0.0
    elif flag in noSlip:
        return lambda x,t: 0.0

dirichletConditions = {0:getDBC_pressure_tube,
                       1:getDBC_u_tube,
                       2:getDBC_v_tube,
                       3:getDBC_w_tube}


# def getPDBC(x):
#    if x[1] < 1.0e-8 or x[1] >= tubeTop - 1.0e-8:
#        print tuple(x)
#        print "link",tuple(numpy.array([x[0],0.0,0.0]))
#        return numpy.array([x[0],0.0,0.0])

# periodicDirichletConditions = {0:getPDBC,1:getPDBC,2:getPDBC}

def getAFBC_p_tube(x,flag):
    if flag == boundaryTags['upstream']:
        return lambda x,t: -Um
    if flag in noSlip:
        return lambda x,t: 0.0

def getAFBC_tube(x,flag):
    return None

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow',
                          3:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_tube,
                                    1:getAFBC_tube,
                                    2:getAFBC_tube,
                                    3:getAFBC_tube}

def getDFBC_tube(x,flag):
    if flag == boundaryTags['downstream']:
        return lambda x,t: 0.0

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_tube},
                                   2:{2:getDFBC_tube},
                                   3:{3:getDFBC_tube}}

class Steady_p:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return -(inflow_height-x[1])*coefficients.rho*coefficients.g[1]

class Steady_u:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

class Steady_v:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0
    
class Steady_w:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:Steady_p(),
                     1:Steady_u(),
                     2:Steady_v(),
                     3:Steady_w()}


## @}

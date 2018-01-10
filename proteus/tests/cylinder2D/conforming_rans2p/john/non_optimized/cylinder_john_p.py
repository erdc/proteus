from proteus import *
from proteus.default_p import *
import obstacle2d
import sys
"""
Incompressible Navier-Stokes flow around a cylinder, 2d version

@Article{Shafer_etal_96,
  author =       {M. Shafer and S. Turek and F. Durst and E. Krause and R. Rannacher},
  title =        {Benchmark Computations of  Laminar Flow Around a Cylinder},
  journal =      {Notes on Numerical Fluid Mechanics},
  year =         1996,
  volume =    52,
  pages =     {547-566}}
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


rho=1.0
nu=1.0e-3
inflow_height=0.41 #0.41
fl_H = inflow_height
bottom_length=2.2
cylinder_radius=0.1/2.0
cylinder_center = [0.15+0.1/2.0,0.15+0.1/2.0]
RE = 120.0
name = "cylinder2d_john_"#+str(RE)+"_"
Ubar = 1.0 #nu*RE/(2.0*cylinder_radius)
#Um = 3.0*Ubar/2.0
T = 8.0 #20.0*bottom_length/Um

#polyfile = "cylinder"
 
# boundaryTags = obstacle2d.genPoly(polyfile,
#                                   cross_section=obstacle2d.circular_cross_section,
#                                   height=inflow_height,
#                                   length=bottom_length,
#                                   radius=cylinder_radius,
#                                   center = cylinder_center,
#  
#                                 n_points_on_obstacle=2*41-2)
from cylinder2dDomain import *
domain = cylinder2d(height = inflow_height,
                    length = bottom_length,
                    radius = cylinder_radius,
                    center = cylinder_center,
                    n_points_on_obstacle = 2*41-2,
                    cross_section = circular_cross_section)
boundaryTags=domain.boundaryTags
domain.writePoly("cylinder")
initialConditions = None

analyticalSolution = None

coefficients = TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                             sigma=0.0,
                                             rho_0=rho,nu_0=nu,
                                             rho_1=rho,nu_1=nu,
                                             g=[0.0,0.0],
                                             nd=2,
                                             LS_model=None,
                                             KN_model=None,
                                             epsFact_density=None,
                                             stokes=False)
#coefficients = NavierStokes(g=[0.0,0.0],nd=nd)
#coefficients = StokesP(g=[0.0,-9.8],nd=nd)
mu = coefficients.rho_0*coefficients.nu_0

#now define the Dirichlet boundary conditions

#spin up the inflow
def velRamp(t):
    #residence_time = bottom_length/Um
    return 0.41**-2*sin(pi*t/8.0)#1.0-exp(-25.0*t/residence_time)
   

def getDBC_pressure_tube(x,flag):
    if flag == boundaryTags['downstream']:
        return lambda x,t: 0.0 #-(inflow_height-x[1])*coefficients.rho*coefficients.g[1]
    else:
        pass

noSlip = [boundaryTags['obstacle'],boundaryTags['bottom'],boundaryTags['top']]
          
def getDBC_u_tube(x,flag):
    if flag == boundaryTags['upstream']:
        return lambda x,t: velRamp(t)*6.0*x[1]*(fl_H-x[1])
    elif flag in noSlip:
        return lambda x,t: 0.0

def getDBC_v_tube(x,flag):
    if flag == boundaryTags['upstream']:
        return lambda x,t: 0.0
    elif flag in noSlip:
        return lambda x,t: 0.0

dirichletConditions = {0:getDBC_pressure_tube,
                       1:getDBC_u_tube,
                       2:getDBC_v_tube}

def getAFBC_p_tube(x,flag):
    if flag in noSlip:
        return lambda x,t: 0.0

def getDFBC_downstream(x,flag):
    if flag == boundaryTags['downstream']:
        return lambda x,t: 0.0

def getNoBC(x,flag):
    pass

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_tube,
                                    1:getNoBC,
                                    2:getNoBC}

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_downstream},
                                   2:{2:getDFBC_downstream}}

class Steady_p:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0 #-(inflow_height-x[1])*coefficients.rho*coefficients.g[1]

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

initialConditions = {0:Steady_p(),
                     1:Steady_u(),
                     2:Steady_v()}

## @}

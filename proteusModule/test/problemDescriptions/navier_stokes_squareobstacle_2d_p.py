from pyadh import *
from pyadh.default_p import *
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

polyfile = "squareobstacle"

initialConditions = None

analyticalSolution = None

coefficients = NavierStokes(g=[0.0,0.0],nd=nd)
#coefficients = Stokes(g=[0.0,-9.8],nd=nd)

#now define the Dirichlet boundary conditions
inflow = 100.0
tubeTop = 11.0
tubeEnd = 20.0
squareTop = 6.0
squareBottom = 5.0
squareLeft = 2.0
squareRight = 3.0
triangleOptions="q30Dena%f" % (0.5*(0.05*tubeTop)**2,)

def velRamp(t):
    if t < 25.0/(tubeEnd/inflow):
        return 1.0-exp(-t*tubeEnd/inflow)
    else:
        return 1.0
def getDBC_pressure_tube(x):
    #at end of tube to hydrostatic
    if x[0] == tubeEnd:
        return lambda x,t: -(tubeTop-x[1])*coefficients.rho*coefficients.g[1]
    else:
        pass

def getDBC_u_tube(x):
    #inflow
    if (x[0] == 0.0):
        return lambda x,t: inflow*velRamp(t)
    #obstacle
    elif (x[1] <= squareTop and
          x[1] >= squareBottom and
          x[0] >= squareLeft and
          x[0] <= squareRight):
        return lambda x,t: 0.0
    #top and bottom
#     if (x[1] == 0.0 or
#         x[1] == tubeTop):
#         return lambda x,t: 0.0

def getDBC_v_tube(x):
    #top and bottom
    if (x[1] == 0.0 or
        x[1] == tubeTop):
        return lambda x,t: 0.0
    #obstacle
    elif (x[1] <= squareTop and
          x[1] >= squareBottom and
          x[0] >= squareLeft and
          x[0] <= squareRight):
        return lambda x,t: 0.0
#    elif (x[0] == 0.0):
#        return lambda x,t: 0.0
dirichletConditions = {0:getDBC_pressure_tube,
                       1:getDBC_u_tube,
                       2:getDBC_v_tube}


# def getPDBC(x):
#    if x[1] < 1.0e-8 or x[1] >= tubeTop - 1.0e-8:
#        print tuple(x)
#        print "link",tuple(numpy.array([x[0],0.0,0.0]))
#        return numpy.array([x[0],0.0,0.0])

# periodicDirichletConditions = {0:getPDBC,1:getPDBC,2:getPDBC}

def getAFBC_p_tube(x):
    #top and bottom
    if (x[1] == tubeTop or
        x[1] == 0.0):
        return lambda x,t: 0.0
    #inflow
    elif (x[0] == 0.0):
        return lambda x,t: -inflow*velRamp(t)
    #obstacle
    elif (x[1] <= squareTop and
          x[1] >= squareBottom and
          x[0] >= squareLeft and
          x[0] <= squareRight):
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
        return -(tubeTop-x[1])*coefficients.rho*coefficients.g[1]

class Steady_u:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0
        #return inflow

class Steady_v:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:Steady_p(),
                     1:Steady_u(),
                     2:Steady_v()}

## @}

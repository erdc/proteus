from pyadh import *
from pyadh.default_p import *
"""
Incompressible Navier-Stokes flow between parallel plates (plane Couette flow).
"""

##\page Tests Test Problems 
# \ref navier_stokes_couette_ss_2d_p.py "Incompressible Navier-Stokes between parallel plates (plane Couette flow)"
#

##\ingroup test
#\file navier_stokes_couette_ss_2d_p.py
#
#\brief Incompressible Navier-Stokes flow between parallel plates (plane Couette flow).
#
#The Navier-Stokes equations implemented in the pyadh::TransportCoefficients::NavierStokes class.
#
#The domain is a rectangle, \f$[0,L_x] \time [0,L_y]\f$, and the boundary conditions are
#\f{eqnarray*}
#p(x=0,,y=0.0,t) &=& 0.0 \\
#u(x,y=1.0,t) &=& 10000.0 \\
#u(x,y=0.0,t) &=& u(0.0,y,t) = u(1.0,y,t) = 0.0 \\
#v(x,y,t) &=& 0.0 \mbox{ on } \partial \Omega\\
#\f}
#
#\image html navier_stokes_couette_v.jpg
#
#
nd = 2

plateSeperation = 1.0 #m
plateLength = 1.0
upperPlateVelocity = 1.0 #m/s

L=(plateLength,plateSeperation,1.0)

pSol = AnalyticalSolutions.PlaneCouetteFlow_p(plateSeperation,upperPlateVelocity)
uSol = AnalyticalSolutions.PlaneCouetteFlow_u(plateSeperation,upperPlateVelocity)
vSol = AnalyticalSolutions.PlaneCouetteFlow_v(plateSeperation,upperPlateVelocity)

analyticalSolution = {0:pSol,
                      1:uSol,
                      2:vSol}

coefficients = NavierStokes(g=[0.0,0.0],nd=2)

print "Re = %12.5e" % (upperPlateVelocity/coefficients.nu)

#now define the Dirichlet boundary conditions

def getDBC_pressure(x,flag):
    if (x[0] in [0.0,plateLength] or
        x[1] in [0.0,plateSeperation]):
        return lambda x,t: pSol.uOfX(x)

def getDBC_u(x,flag):
    if (x[0] in [0.0,plateLength] or
        x[1] in [0.0,plateSeperation]):
        return lambda x,t: uSol.uOfX(x)

def getDBC_v(x,flag):
    if (x[0] in [0.0,plateLength] or
        x[1] in [0.0,plateSeperation]):
        return lambda x,t: vSol.uOfX(x)

dirichletConditions = {0:getDBC_pressure,
                       1:getDBC_u,
                       2:getDBC_v}

initialConditions  = None

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow'}

def getAFBC_p_tube(x,flag):
    #top and bottom
    if x[1] in [0.0,plateSeperation]:
        return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC_p_tube}

diffusiveFluxBoundaryConditions = {0:{}}

#T = 1.0e5

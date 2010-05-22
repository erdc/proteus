from pyadh import *
from pyadh.default_p import *
"""
Incompressible Navier-Stokes flow between parallel plates (plane Poiseuille flow).
"""

##\page Tests Test Problems 
# \ref navier_stokes_couette_ss_2d_p.py "Incompressible Navier-Stokes between parallel plates (plane Poiseuille flow)"
#

##\ingroup test
#\file navier_stokes_couette_ss_2d_p.py
#
#\brief Incompressible Navier-Stokes flow between parallel plates (plane Poiseuille flow).
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

#this are the parameters required to define the problem (and analytical solution)
plateSeperation = 1.0 #m
#assumed infinite but we need to define something
plateLength = 1.0 #m

#domain dimensions
L=(plateLength,plateSeperation,1.0)

coefficients = TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                             sigma=0.0,
                                             rho_0=998.2,nu_0=1.004e-6,
                                             rho_1=998.2,nu_1=1.004e-6,
                                             g=[0.0,0.0],
                                             nd=2,
                                             LS_model=None,
                                             KN_model=None,
                                             epsFact_density=None,
                                             stokes=True)
# coefficients = ThreephaseNavierStokes_ST_LS_SO(epsFact=1.5,
#                  sigma=0.0,
#                  rho_0=998.2,nu_0=1.004e-6,
#                  rho_1=998.2,nu_1=1.004e-6,
#                  rho_s=998.2,nu_s=1.004e-6,
#                  g=[0.0,0.0],
#                  nd=2,
#                  LS_model=None,
#                  KN_model=None,
#                  epsFact_density=1.5,
#                  stokes=True)
mu = coefficients.rho_0*coefficients.nu_0
# coefficients = NavierStokes(g=[0.0,0.0],nd=2)
#coefficients = Stokes(g=[0.0,0.0],nd=2)
#coefficients = TwophaseNavierStokes_ST_LS_SO(g=[9.8,0.0],nd=2)
#calculated dynamic viscosity
mu = coefficients.rho*coefficients.nu


plane_theta=math.pi/2.0#+math.pi/6.0
plane_phi=math.pi/2.0
v_theta=0.0#+math.pi/6.0
v_phi=math.pi/2.0
v_norm=1.0#1.0e-4
grad_p=-1.0
#grad_p=0.0
#construct analytical solution objects
#pSol = AnalyticalSolutions.PlanePoiseuilleFlow_p(plateSeperation,mu,grad_p)
#uSol = AnalyticalSolutions.PlanePoiseuilleFlow_u(plateSeperation,mu,grad_p)
#vSol = AnalyticalSolutions.PlanePoiseuilleFlow_v(plateSeperation,mu,grad_p)

pSol = AnalyticalSolutions.PlanePoiseuilleFlow_p2(plane_theta,
                                                  plane_phi,
                                                  v_theta,
                                                  v_phi,
                                                  v_norm,
                                                  mu,
                                                  grad_p,
                                                  L)
uSol = AnalyticalSolutions.PlanePoiseuilleFlow_u2(plane_theta,
                                                  plane_phi,
                                                  v_theta,
                                                  v_phi,
                                                  v_norm,
                                                  mu,
                                                  grad_p,
                                                  L)
vSol = AnalyticalSolutions.PlanePoiseuilleFlow_v2(plane_theta,
                                                  plane_phi,
                                                  v_theta,
                                                  v_phi,
                                                  v_norm,
                                                  mu,
                                                  grad_p,
                                                  L)
#in 3D we could add wSol = AnalyticalSolutions.PlanePoiseuilleFlow_v(plateSeperation,mu,grad_p) since v=w=0

#load analytical solutions into dictionary of solutions for the component
analyticalSolution = {0:pSol,
                      1:uSol,
                      2:vSol}


#print "Re = %12.5e" % (upperPlateVelocity/coefficients.nu)

#now define the Dirichlet boundary conditions
#use the analytical solutions to get the Dirichlet boundary conditions

def getDBC_pressure(x,flag):
    #if x is on the boundary
#     if (x[0] in [0.0,plateLength] or
#         x[1] in [0.0,plateSeperation]):
#         return lambda x,t: pSol.uOfX(x)
#     if (x[1] in [0.0,plateSeperation]):
#         return lambda x,t: pSol.uOfX(x)
    if (x[0] in [0.0,plateLength]):
        return lambda x,t: pSol.uOfX(x)
#     #if x is on the boundary
#     if (x[0] == 0.0):
#         return lambda x,t: 1.0
#     if (x[0] == L[0]):
#         return lambda x,t: 0.0

def getDBC_u(x,flag):
    if x[1] in [0.0,plateSeperation]:
        return lambda x,t: 0.0
#     if (x[0] in [0.0] or
#         x[1] in [0.0]):
#         return lambda x,t: uSol.uOfX(x)
#     if (x[0] in [0.0,plateLength] or
#         x[1] in [0.0,plateSeperation]):
#         return lambda x,t: uSol.uOfX(x)
#     if (x[1] in [0.0,plateSeperation]):
#         return lambda x,t: 0.0

def getDBC_v(x,flag):
    if x[1] in [0.0,plateSeperation]:
        return lambda x,t: 0.0
#     if (x[0] in [0.0] or
#         x[1] in [0.0]):
#         return lambda x,t: vSol.uOfX(x)
#     if (x[0] in [0.0,plateLength] or
#         x[1] in [0.0,plateSeperation]):
#         return lambda x,t: vSol.uOfX(x)
#     if (x[1] in [0.0,plateSeperation]):
#         return lambda x,t: 0.0

dirichletConditions = {0:getDBC_pressure,
                       1:getDBC_u,
                       2:getDBC_v}

initialConditions  = {0:pSol,
                     1:uSol,
                     2:vSol}

fluxBoundaryConditions = {0:'setFlow',
                          1:'setFlow',
                          2:'setFlow'}

def getAFBC_p_tube(x,flag):
    #pass
    #top and bottom
    if x[1] in [0.0,plateSeperation]:
        return lambda x,t: 0.0
#    if x[0] in [0.0,plateLength]:
#        return lambda x,t: 0.0

def getAFBC_u_tube(x,flag):
    #pass
    #top and bottom
#     if x[0] == 0.0:
#         return lambda x,t: 0.0
    if x[0] in [0.0,plateLength]:
        return lambda x,t: 0.0

def getAFBC_v_tube(x,flag):
    #pass
    #top and bottom
#     if x[0] == 0.0:
#         return lambda x,t: 0.0
    if x[0] in [0.0,plateLength]:
        return lambda x,t: 0.0
    
advectiveFluxBoundaryConditions =  {0:getAFBC_p_tube,
                                    1:getAFBC_u_tube,
                                    2:getAFBC_v_tube}

def getDFBC(x,flag):
    if x[0] in [0.0,plateLength]:
        return lambda x,t: 0.0
    #pass

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC,2:getDFBC},
                                   2:{2:getDFBC,1:getDFBC}}

def getPDBC_p(x,flag):
    pass
#   if x[0] in [0.0,plateLength]:
#       if x[1] > 0.0 and x[1] < plateSeperation:
#           print tuple(x)
#           print "link",tuple(numpy.array([0.0,x[1],0.0]))
#           return numpy.array([0.0,x[1],0.0])
#   if x[0] in [0.0,plateLength]:
#       print tuple(x)
#       print "link",tuple(numpy.array([x[0],0.0,0.0]))
#       return numpy.array([x[0],0.0,0.0])
def getPDBC(x,flag):
   if x[0] in [0.0,plateLength]:
       if x[1] > 0.0 and x[1] < plateSeperation:
           print tuple(x)
           print "link",tuple(numpy.array([0.0,x[1],0.0]))
           return numpy.array([0.0,x[1],0.0])

periodicDirichletConditions = {0:getPDBC_p,1:getPDBC,2:getPDBC}

#end time of dynamic simulation
#T = 1.0e5

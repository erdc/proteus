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
nd = 3

#domain dimensions
L=(1.0,1.0,1.0)

coefficients = TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                             sigma=0.0,
                                             rho_0=998.2,nu_0=1.004e-6,
                                             rho_1=998.2,nu_1=1.004e-6,
                                             g=[0.0,0.0,0.0],
                                             nd=nd,
                                             LS_model=None,
                                             KN_model=None,
                                             epsFact_density=None,
                                             stokes=True)
mu = coefficients.rho_0*coefficients.nu_0
# coefficients = NavierStokes(g=[0.0,0.0],nd=2)
#coefficients = Stokes(g=[0.0,0.0,0.0],nd=3)
#coefficients = TwophaseNavierStokes_ST_LS_SO(g=[9.8,0.0],nd=2)
#calculated dynamic viscosity
#mu = coefficients.rho*coefficients.nu

plane_theta=math.pi/6.0
plane_phi=math.pi/2.0-math.pi/3.0
v_theta=math.pi/3.0
v_phi=math.pi/2.0-math.pi/3.0
v_norm=1.0
grad_p=-1.0#-1.0e-2
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
wSol = AnalyticalSolutions.PlanePoiseuilleFlow_w2(plane_theta,
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
                      2:vSol,
                      3:wSol}


#print "Re = %12.5e" % (upperPlateVelocity/coefficients.nu)

#now define the Dirichlet boundary conditions
#use the analytical solutions to get the Dirichlet boundary conditions

def getDBC_pressure(x,flag):
    if (x[0] in [0.0,L[0]] or
        x[1] in [0.0,L[1]] or
        x[2] in [0.0,L[2]]):
        return lambda x,t: pSol.uOfX(x)

def getDBC_u(x,flag):
    if (x[0] in [0.0,L[0]] or
        x[1] in [0.0,L[1]] or
        x[2] in [0.0,L[2]]):
        return lambda x,t: uSol.uOfX(x)

def getDBC_v(x,flag):
    if (x[0] in [0.0,L[0]] or
        x[1] in [0.0,L[1]] or
        x[2] in [0.0,L[2]]):
        return lambda x,t: vSol.uOfX(x)

def getDBC_w(x,flag):
    if (x[0] in [0.0,L[0]] or
        x[1] in [0.0,L[1]] or
        x[2] in [0.0,L[2]]):
        return lambda x,t: wSol.uOfX(x)

dirichletConditions = {0:getDBC_pressure,
                       1:getDBC_u,
                       2:getDBC_v,
                       3:getDBC_w}

#initialConditions  = {0:pSol,
#                      1:uSol,
#                      2:vSol,
#                      3:wSol}

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow',
                          3:'outFlow'}

def getAFBC_p_tube(x,flag):
    pass
#     #top and bottom
#     if x[1] in [0.0,plateSeperation]:
#         return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC_p_tube}

def getDFBC(x,flag):
    pass

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC},
                                   2:{2:getDFBC},
                                   3:{3:getDFBC}}


#end time of dynamic simulation
#T = 1.0e5

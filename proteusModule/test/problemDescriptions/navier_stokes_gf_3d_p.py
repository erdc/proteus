from pyadh import *
from pyadh.default_p import *
from pyadh import RANS2P

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
nu=1.0#1.0e-3

genMesh=True

# gf = InputTranslators.GF("Container06.GF",
#                          boundingBox=True,
#                          insertHoles=True)
# polyfile="Container06_3d"#gf.polyfile

# domain = Domain.PiecewiseLinearComplexDomain(polyfile)
# domain.writeAsymptote("Container06")
# domain.writePLY("Container06")
# length_scale = max([gf.Lx,gf.Ly,gf.Lz])
# inflow_height= gf.Lz
# bottom_width = gf.Ly
# bottom_length= gf.Lx
from vesselPoly3dDomain import *
domain = vesselPoly3d(fileprefix="Container06_test",vesselprefix="Container06_3d",outerBoxFactor=[10.0,10.0,0.1],offsetFactor=[0.6,0.0,0.0])
domain.writePoly("Container06_test")
length_scale = max(domain.L)
inflow_height= domain.L[2]
bottom_width = domain.L[1]
bottom_length= domain.L[0]

#REs = raw_input("Enter Reynolds Number\n")
<<<<<<< .mine
length_scale = max([gf.Lx,gf.Ly,gf.Lz])
inflow_height= gf.Lz
bottom_width = gf.Ly
bottom_length= gf.Lx
REs="10.0"
name = "ns_gf_3d"+REs+"_"
=======
REs="1000.0"
name = "ns_gf_3d_3p"+REs+"_"
>>>>>>> .r3170
RE = float(REs)
Ubar = nu*RE/(2.0*length_scale)
Um = 9.0*Ubar/4.0
T = 1#20.0*length_scale/Um
print "REYNOLDS NUMBER = ",Ubar*2.0*1.0/nu


initialConditions = None

analyticalSolution = None

def waterLineSignedDistanceAndNormal(t,x,phi,n,elevation=domain.L[2]*0.5):
    """
    try to use a level set representation for water line
    """
    npoints = len(phi.flat)
    for k in range(npoints):
        u = 0.0; v = 0.0; w = elevation-x.flat[k*3+2]
        phi.flat[k] = w
        nx = 0.0; ny = 0.0; nz = -1.0
        n.flat[k*3+0]=nx; n.flat[k*3+1]=ny; n.flat[k*3+2] = nz

rho_a = rho*1.0e-1; nu_a = nu
# coefficients = ThreephaseNavierStokes_ST_LS_SO(epsFact=0.0,
#                                                sigma=0.0,
#                                                rho_0=rho,nu_0=nu,
#                                                rho_1=rho_a,nu_1=nu_a,
#                                                g=[0.0,0.0,0.0],
#                                                nd=nd,
#                                                LS_model=None,
#                                                KN_model=None,
#                                                epsFact_density=None,
#                                                stokes=False,
#                                                defaultFluidProfile=waterLineSignedDistanceAndNormal)
LevelModelType = RANS2P.OneLevelRANS2P
coefficients = TwophaseNavierStokes_ST_LS_SO(epsFact=3.0,
                                             sigma=0.0,
                                             rho_0 = rho,
                                             nu_0 = nu,
                                             rho_1 = rho,
                                             nu_1 = nu,
                                             g=[0.0,0.0,-9.8],
                                             nd=nd,
                                             LS_model=None,
                                             epsFact_density=3.0,
                                             stokes=False)

# coefficients = TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
#                                              sigma=0.0,
#                                              rho_0=rho,nu_0=nu,
#                                              rho_1=rho,nu_1=nu,
#                                              g=[0.0,0.0,0.0],
#                                              nd=nd,
#                                              LS_model=None,
#                                              KN_model=None,
#                                              epsFact_density=None,
#                                              stokes=False)
mu = coefficients.rho_0*coefficients.nu_0
#coefficients = NavierStokes(g=[0.0,0.0],nd=nd)
#coefficients = Stokes(g=[0.0,-9.8],nd=nd)

#now define the Dirichlet boundary conditions

def velRamp(t):
    residence_time = bottom_length/Um
    if t < residence_time and False:#mwf hack turn off for steady state
        return 1.0-exp(-25.0*t/residence_time)
    else:
        return 1.0

def getDBC_pressure_tube(x,flag):
    if flag == domain.boundaryTags['right']:
        return lambda x,t: -(inflow_height-x[2])*coefficients.rho*coefficients.g[2]
    else:
        pass

noSlip = [3,4,5,6,7]

def getDBC_u_tube(x,flag):
    if flag == 1:
        return lambda x,t: velRamp(t)*Um
    elif flag in noSlip:
        return lambda x,t: 0.0

def getDBC_v_tube(x,flag):
    if flag == 1:
        return lambda x,t: 0.0
    elif flag in noSlip:
        return lambda x,t: 0.0

def getDBC_w_tube(x,flag):
    if flag == 1:
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
    if flag == 1:
        return lambda x,t: -velRamp(t)*Um
    if flag in noSlip:
        return lambda x,t: 0.0

def getNoBC(x,flag):
    pass

def getDFBC_vel_tube(x,flag):
    if flag == 2:
        return lambda x,t: 0.0

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow',
                          3:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_tube,1:getNoBC,2:getNoBC,3:getNoBC}

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_vel_tube,2:getDFBC_vel_tube,3:getDFBC_vel_tube},
                                   2:{1:getDFBC_vel_tube,2:getDFBC_vel_tube,3:getDFBC_vel_tube},
                                   3:{1:getDFBC_vel_tube,2:getDFBC_vel_tube,3:getDFBC_vel_tube}}
diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_vel_tube},
                                   2:{2:getDFBC_vel_tube},
                                   3:{3:getDFBC_vel_tube}}

class Steady_p:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return -(inflow_height-x[2])*coefficients.rho*coefficients.g[2]

class Steady_u:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return Um

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

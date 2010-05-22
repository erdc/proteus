from pyadh import *
from pyadh.default_p import *
import cylinder2d
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

def genCylinderPoly_fib(fileprefix,
                        height=1.0,
                        length=1.0,
                        radius=0.25,
                        center=(0.5,0.5)):
    boundaries=['upstream',
                'downstream',
                'bottom',
                'top']
                #'obstacle']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    boundaryTags['obstacle'] = -100 #unreachable
    #work around the domain from (0.0,0.0) going counterclockwise
    vertices = [('upstream_bottom', (0.0,0.0),boundaryTags['bottom']),
                ('downstream_bottom',(length,0.0),boundaryTags['bottom']),
                ('downstream_top',(length,height),boundaryTags['top']),
                ('upstream_top',(0.0,height),boundaryTags['top'])]
    nv = len(vertices)
    segments = [(0,1,boundaryTags['bottom']),
                (1,2,boundaryTags['downstream']),
                (2,3,boundaryTags['top']),
                (3,0,boundaryTags['upstream'])]

    poly = open(fileprefix+'.poly','w')
    poly.write('%d %d %d %d \n' % (len(vertices),2,0,1))
    poly.write("#vertices \n")
    for vN,v in enumerate(vertices):
        poly.write('%d %12.5e %12.5e %d #%s \n' % (vN+1,v[1][0],v[1][1],v[2],v[0]))
    poly.write('%d %d \n' % (len(vertices),1))
    poly.write("#segments \n")
    for sN,s in enumerate(segments):
        poly.write('%d %d %d %d \n' % (sN+1,s[0]+1,s[1]+1,s[2]))
    poly.write("#holes \n")
    poly.write('%d\n' % (0,))

    poly.write("#regions \n")
    poly.write('%d\n' % (1,))
    poly.write('%d %12.5e %12.5e %d\n' % (1,
                                          vertices[0][1][0]+length*1.0e-8,
                                          vertices[0][1][1]+length*1.0e-8,
                                          0+1))
    poly.close()
    return boundaryTags



nd = 2


rho=1.0
nu=1.0e-3
inflow_height=0.41
bottom_length=2.2
cylinder_radius=0.1/2.0
cylinder_center = (0.15+cylinder_radius,0.15+cylinder_radius)
Um = 0.01
#Um = 0.3
Um = 1.5
#REs = raw_input("Enter Reynolds Number\n")
REs = "0.5"#"120.0"
name = "cylinder2d_fib_"+REs+"_"
RE = float(REs)
Ubar = nu*RE/(2.0*cylinder_radius)
Um = 3.0*Ubar/2.0
#T = 20.0*bottom_length/Um
T = 5.0e0*bottom_length/Um
#Ubar = 2.0*Um/3.0

def cylinderSignedDistanceAndNormal(t,x,phi,n):
    """
    signed distance from cylinder for solid boundary
    """
    npoints = len(phi.flat)
    for k in range(npoints):
        u = x.flat[k*3+0] - cylinder_center[0]; v = x.flat[k*3+1] - cylinder_center[1]
        phi.flat[k] = sqrt(u**2 + v**2) - cylinder_radius
        nx = u/(sqrt(u**2 + v**2)); ny = v/(sqrt(u**2 + v**2))
        n.flat[k*2+0]=nx; n.flat[k*2+1] = ny;


print "REYNOLDS NUMBER = ",Ubar*2.0*cylinder_radius/nu
polyfile = "cylinder_fib"

boundaryTags = genCylinderPoly_fib(polyfile,
                                  height=inflow_height,
                                  length=bottom_length,
                                  radius=cylinder_radius,
                                  center = cylinder_center)

initialConditions = None

analyticalSolution = None

coefficients = ThreephaseNavierStokes_ST_LS_SO(epsFact=1.5,#1.5,
                                               sigma=0.0,
                                               rho_0=rho,nu_0=nu,
                                               rho_1=rho,nu_1=nu,
                                               rho_s=rho,nu_s=nu,
                                               g=[0.0,0.0],
                                               nd=2,
                                               LS_model=None,
                                               KN_model=None,
                                               epsFact_density=None,
                                               defaultSolidProfile=cylinderSignedDistanceAndNormal,
                                               stokes=False,
                                               boundaryPenaltyCoef=1.e3,
                                               volumePenaltyCoef=1.e3)
mu = coefficients.rho_0*coefficients.nu_0
#coefficients = NavierStokes(g=[0.0,0.0],nd=nd)
#coefficients = StokesP(g=[0.0,-9.8],nd=nd)

#now define the Dirichlet boundary conditions

def velRamp(t):
    residence_time = bottom_length/Um
    if t < residence_time:
        return 1.0-exp(-25.0*t/residence_time)
    else:
        return 1.0

def getDBC_pressure_tube(x,flag):
    if flag == boundaryTags['downstream']:
        return lambda x,t: -(inflow_height-x[1])*coefficients.rho*coefficients.g[1]
    else:
        pass

noSlip = [boundaryTags['obstacle'],boundaryTags['bottom'],boundaryTags['top']]
          
def getDBC_u_tube(x,flag):
    if flag == boundaryTags['upstream']:
        return lambda x,t: velRamp(t)*4.0*Um*x[1]*(inflow_height-x[1])/(inflow_height**2)
    elif flag in noSlip:
        return lambda x,t: 0.0
    #mwf play with strong bc's in cylinder
    #phi,n = cylinderSignedDistanceAndNormal(x)
    #if phi <= 0.0:
    #    return lambda x,t:0.0
def getDBC_v_tube(x,flag):
    if flag == boundaryTags['upstream']:
        return lambda x,t: 0.0
    elif flag in noSlip:
        return lambda x,t: 0.0
    #mwf play with strong bc's in cylinder
    #phi,n = cylinderSignedDistanceAndNormal(x)
    #if phi <= 0.0:
    #   return lambda x,t:0.0

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

advectiveFluxBoundaryConditions =  {0:getAFBC_p_tube,1:getNoBC,2:getNoBC}

diffusiveFluxBoundaryConditions = {0:{},1:{1:getDFBC_downstream},2:{2:getDFBC_downstream}}

# advectiveFluxBoundaryConditions =  {0:getAFBC_p_tube}

# diffusiveFluxBoundaryConditions = {0:{},1:{1:getNoBC},2:{2:getNoBC}}

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

initialConditions = {0:Steady_p(),
                     1:Steady_u(),
                     2:Steady_v()}


## @}

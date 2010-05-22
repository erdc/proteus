from pyadh import *
from pyadh.default_p import *
import cylinder3d
import sys
"""
Incompressible Navier-Stokes flow over a cylinder in 3D.
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
                        width=1.0,
                        radius=0.25,
                        center=(0.5,0.5),
                        cross_section=cylinder3d.circular_cross_section,
                        thetaOffset=0.0):
    n_points_on_obstacle = 0
    boundaries=['upstream',
                'downstream',
                'bottom',
                'top',
                'front',
                'back']
                #'obstacle']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    boundaryTags['obstacle'] = -100 #unreachable
    #work around the domain from (0.0,0.0) going counterclockwise
    vertices = [('upstream_bottom', (0.0,0.0),boundaryTags['bottom']),
                ('downstream_bottom',(length,0.0),boundaryTags['bottom']),
                ('downstream_top',(length,height),boundaryTags['top']),
                ('upstream_top',(0.0,height),boundaryTags['top'])]
    nv = len(vertices)
    #now need to convert rep to 3D
    vertices3d={}
    for vN,v in enumerate(vertices):
        vertices3d[v[0]+'_front']=(vN,(v[1][0],0.0,v[1][1]),v[2])
    for vN,v in enumerate(vertices):
        vertices3d[v[0]+'_back']=(vN+len(vertices),(v[1][0],width,v[1][1]),v[2])

    poly = open(fileprefix+'.poly','w')
    poly.write('%d %d %d %d \n' % (len(vertices3d),3,0,1))
    poly.write("#vertices \n")
    a=vertices3d.values()
    a.sort()
    for v in a:
        poly.write('%d %12.5e %12.5e %12.5e %d #%s \n' % (1+v[0],v[1][0],v[1][1],v[1][2],v[2],v[-1]))
    poly.write("#facets \n")
    poly.write('%d %d \n' % (6+n_points_on_obstacle,1))
    #upstream
    poly.write('%d %d %d \n' % (1,0,boundaryTags['upstream'])) #1 polygon, 0 holes, tag
    poly.write('%d %d %d %d %d \n' % (4,
                                      1+vertices3d['upstream_bottom_front'][0],
                                      1+vertices3d['upstream_bottom_back'][0],
                                      1+vertices3d['upstream_top_back'][0],
                                      1+vertices3d['upstream_top_front'][0]))
    #downstream
    poly.write('%d %d %d \n' % (1,0,boundaryTags['downstream'])) #1 polygon, 0 holes, tag
    poly.write('%d %d %d %d %d \n' % (4,
                                      1+vertices3d['downstream_bottom_front'][0],
                                      1+vertices3d['downstream_bottom_back'][0],
                                      1+vertices3d['downstream_top_back'][0],
                                      1+vertices3d['downstream_top_front'][0]))
    #top
    poly.write('%d %d %d \n' % (1,0,boundaryTags['top'])) #1 polygon, 0 holes, tag
    poly.write('%d %d %d %d %d \n' % (4,
                                      1+vertices3d['upstream_top_front'][0],
                                      1+vertices3d['downstream_top_front'][0],
                                      1+vertices3d['downstream_top_back'][0],
                                      1+vertices3d['upstream_top_back'][0]))
    #bottom
    poly.write('%d %d %d \n' % (1,0,boundaryTags['bottom'])) #1 polygon, 0 holes, tag
    poly.write('%d %d %d %d %d \n' % (4,
                                      1+vertices3d['upstream_bottom_front'][0],
                                      1+vertices3d['downstream_bottom_front'][0],
                                      1+vertices3d['downstream_bottom_back'][0],
                                      1+vertices3d['upstream_bottom_back'][0]))
    #front
    poly.write('%d %d %d \n' % (1,0,boundaryTags['front'])) #1 polygon, 0 holes, tag
    poly.write('%d %d %d %d %d \n' % (4,
                                      1+vertices3d['upstream_bottom_front'][0],
                                      1+vertices3d['downstream_bottom_front'][0],
                                      1+vertices3d['downstream_top_front'][0],
                                      1+vertices3d['upstream_top_front'][0]))
    #back
    poly.write('%d %d %d \n' % (1,0,boundaryTags['back'])) #1 polygon, 0 holes, tag
    poly.write('%d %d %d %d %d \n' % (4,
                                      1+vertices3d['upstream_bottom_back'][0],
                                      1+vertices3d['downstream_bottom_back'][0],
                                      1+vertices3d['downstream_top_back'][0],
                                      1+vertices3d['upstream_top_back'][0]))
    #finished with facets
    poly.write("#holes \n")
    poly.write('%d\n' % (0,))
    poly.write("#regions \n")
    poly.write('%d\n' % (1,))
    #use a perturbation of the lower bottomhand corner to get a pont inside
    poly.write('%d %12.5e %12.5e %12.5e %d\n' % (1,
                                                 1.0e-8*height,
                                                 1.0e-8*height,
                                                 1.0e-8*height,
                                                 0+1))
    poly.close()
    return boundaryTags

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
    #poly.write("#holes \n")
    #poly.write('%d\n' % (0,))

    poly.write("#regions \n")
    poly.write('%d\n' % (1,))
    poly.write('%d %12.5e %12.5e %d\n' % (1,
                                          vertices[0][1][0]+length*1.0e-8,
                                          vertices[0][1][1]+length*1.0e-8,
                                          0+1))
    poly.close()
    return boundaryTags



nd = 3


rho=1.0
nu=1.0e-3
inflow_height=0.41
bottom_width = 0.41
bottom_length=2.5
cylinder_radius=0.1/2.0#0.1/2.0
cylinder_center = (0.45+cylinder_radius,
                   0.15+cylinder_radius)
                   

#REs = raw_input("Enter Reynolds Number\n")
REs = "1.0"#"120.0"
name = "cylinder3d_fib_"+REs+"_"
RE = float(REs)
Ubar = nu*RE/(2.0*cylinder_radius)
Um = 9.0*Ubar/4.0
#T = 20.0*bottom_length/Um
T = 5.0e0*bottom_length/Um
print "REYNOLDS NUMBER = ",Ubar*2.0*cylinder_radius/nu
#Ubar = 2.0*Um/3.0

def cylinderSignedDistanceAndNormal(t,x,phi,n):
    """
    signed distance from cylinder for solid boundary note cylinder oriented along y axis
    """
    npoints = len(phi.flat)
    for k in range(npoints):
        u = x.flat[k*3+0] - cylinder_center[0]; w = x.flat[k*3+2] - cylinder_center[1]; v = 0.0;
        phi.flat[k] = sqrt(u**2 + w**2) - cylinder_radius
        nx = u/(sqrt(u**2 + w**2)); ny = 0.0; nz=w/(sqrt(u**2 +  w**2))
        n.flat[k*3+0]=nx; n.flat[k*3+1] = ny; n.flat[k*3+2] = nz; 
#    u = x[0] - cylinder_center[0]; w = x[2] - cylinder_center[1]; v = 0.0;
#    phi = sqrt(u**2 + w**2) - cylinder_radius
#    nx = u/(sqrt(u**2 + w**2)); ny = 0.0; nz=w/(sqrt(u**2 +  w**2))
#    return phi,numpy.array([nx,ny,nz],'d')

polyfile = "cylinder3d_fib"

boundaryTags = genCylinderPoly_fib(polyfile,
                                   height=inflow_height,
                                   length=bottom_length,
                                   radius=cylinder_radius,
                                   width=bottom_width,
                                   center = cylinder_center)

initialConditions = None

analyticalSolution = None

coefficients = ThreephaseNavierStokes_ST_LS_SO(epsFact=0.5,#1.5,
                                               sigma=0.0,
                                               rho_0=rho,nu_0=nu,
                                               rho_1=rho,nu_1=nu,
                                               rho_s=rho,nu_s=nu,
                                               g=[0.0,0.0,0.0],
                                               nd=3,
                                               LS_model=None,
                                               KN_model=None,
                                               epsFact_density=None,
                                               defaultSolidProfile=cylinderSignedDistanceAndNormal,
                                               stokes=True,
                                               boundaryPenaltyCoef=1.e3,
                                               volumePenaltyCoef=1.e3)
mu = coefficients.rho_0*coefficients.nu_0

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

noSlip = [boundaryTags['obstacle'],boundaryTags['bottom'],boundaryTags['top'],boundaryTags['front'],boundaryTags['back']]
          
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
def getDBC_w_tube(x,flag):
    if flag == boundaryTags['upstream']:
        return lambda x,t: 0.0
    elif flag in noSlip:
        return lambda x,t: 0.0

dirichletConditions = {0:getDBC_pressure_tube,
                       1:getDBC_u_tube,
                       2:getDBC_v_tube,
                       3:getDBC_w_tube}


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
                          2:'outFlow',
                          3:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_tube,1:getNoBC,2:getNoBC,3:getNoBC}

diffusiveFluxBoundaryConditions = {0:{},1:{1:getDFBC_downstream},2:{2:getDFBC_downstream},3:{2:getDFBC_downstream}}


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

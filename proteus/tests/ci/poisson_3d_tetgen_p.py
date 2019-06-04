from __future__ import division
from builtins import range
from builtins import object
from past.utils import old_div
from proteus import *
from proteus.default_p import *
"""
Heterogeneous Poisson's equation, -div(a(x)u) = f(x), on unit domain [0,1]x[0,1]x[0,1]
"""

##\page Tests Test Problems
# \ref poisson_3d_p.py "Heterogeneous Poisson's equation, -div(a(x)u) = f(x), on unit domain [0,1]x[0,1]x[0,1]"
#

##\ingroup test
#\file poisson_3d_p.py
#
#\brief Heterogenous Poisson's equations in 3D unit domain [0,1]x[0,1]x[0,1]

#----------------------------------------------------
# Domain - mesh - quadrature
#----------------------------------------------------
#space dimension
nd = 3

hull_length = 0.5
hull_beam   = 0.5
hull_draft  = 0.5

L=(2.0*hull_length,
   2.0*hull_beam, 
   2.0*hull_draft)

x_ll = (0.0,
        0.0,
        0.0)

hull_center = (0.5*hull_length,
               0.5*hull_beam,
               0.5*hull_draft)

nLevels = 1

he = old_div(L[0],10.0)
#he = hull_draft/1.0
#he = hull_draft/6.0
genMesh=True#False
vessel = None
#vessel = 'cube'
#vessel = 'wigley'
boundaryTags = { 'bottom': 1, 'front':2, 'right':3, 'back': 4, 'left':5, 'top':6, 'obstacle':7}
if vessel is 'wigley-gmsh':
    domain = Domain.MeshTetgenDomain(fileprefix="mesh")
    domain.boundaryTags = boundaryTags
else:
    vertices=[[x_ll[0],x_ll[1],x_ll[2]],#0
              [x_ll[0]+L[0],x_ll[1],x_ll[2]],#1
              [x_ll[0]+L[0],x_ll[1]+L[1],x_ll[2]],#2
              [x_ll[0],x_ll[1]+L[1],x_ll[2]],#3
              [x_ll[0],x_ll[1],x_ll[2]+L[2]],#4
              [x_ll[0]+L[0],x_ll[1],x_ll[2]+L[2]],#5
              [x_ll[0]+L[0],x_ll[1]+L[1],x_ll[2]+L[2]],#6
              [x_ll[0],x_ll[1]+L[1],x_ll[2]+L[2]]]#7
    vertexFlags=[boundaryTags['left'],
                 boundaryTags['right'],
                 boundaryTags['right'],
                 boundaryTags['left'],
                 boundaryTags['left'],
                 boundaryTags['right'],
                 boundaryTags['right'],
                 boundaryTags['left']]
    facets=[[[0,1,2,3]],
            [[0,1,5,4]],
            [[1,2,6,5]],
            [[2,3,7,6]],
            [[3,0,4,7]],
            [[4,5,6,7]]]
    facetFlags=[boundaryTags['bottom'],
                boundaryTags['front'],
                boundaryTags['right'],
                boundaryTags['back'],
                boundaryTags['left'],
                boundaryTags['top']]
    regions=[[x_ll[0]+0.5*L[0],x_ll[1]+0.5*L[1],x_ll[2]+0.5*L[2]]]
    regionFlags=[1.0]
    holes=[]
    if vessel is 'wigley':
        from math import log
        he_hull = old_div(log(64.0*he+1.0),64.0)
        #print he,he_hull
        #he_hull = he
        n_points_length = int(ceil(old_div(hull_length,he_hull)))+1
        n_points_draft  = 2*int(ceil(old_div(hull_draft,he_hull)))+1
        #print "points",n_points_length,n_points_draft
        dx = old_div(hull_length,float(n_points_length-1))
        dz = 2.0*hull_draft/float(n_points_draft-1)
        #print "he",he,dx,dz
        #grid on right half of hull
        for i in range(n_points_length):
            for j in range(n_points_draft):
                x = i*dx - 0.5*hull_length
                z = j*dz - hull_draft
                zStar = min(0.0,z)
                y = 0.5*hull_beam*(1.0 - 4.0*(old_div(x,hull_length))**2) * (1.0 - (old_div(zStar,hull_draft))**2)
                vertices.append([x+hull_center[0],
                                 y+hull_center[1],
                                 z+hull_center[2]])
                vertexFlags.append(boundaryTags['obstacle'])
        def vN_right(i,j):
            return 8 + i*n_points_draft+j
        for i in range(n_points_length-1):
            for j in range(n_points_draft-1):
                if i < old_div(n_points_length,2):
                    facets.append([[vN_right(i,j),vN_right(i+1,j+1),vN_right(i+1,j)]])
                    facetFlags.append(boundaryTags['obstacle'])
                    facets.append([[vN_right(i,j),vN_right(i,j+1),vN_right(i+1,j+1)]])
                    facetFlags.append(boundaryTags['obstacle'])
                else:
                    facets.append([[vN_right(i,j),vN_right(i,j+1),vN_right(i+1,j)]])
                    facetFlags.append(boundaryTags['obstacle'])
                    facets.append([[vN_right(i,j+1),vN_right(i+1,j+1),vN_right(i+1,j)]])
                    facetFlags.append(boundaryTags['obstacle'])                
        #grid on left half of hull
        for i in range(1,n_points_length-1):
            for j in range(1,n_points_draft):
                x = i*dx - 0.5*hull_length
                z = j*dz - hull_draft
                zStar = min(0.0,z)
                y = 0.5*hull_beam*(1.0 - 4.0*(old_div(x,hull_length))**2) * (1.0 - (old_div(zStar,hull_draft))**2)
                vertices.append([x+hull_center[0],
                                 hull_center[1] - y,
                                 z+hull_center[2]])
                vertexFlags.append(boundaryTags['obstacle'])
        def vN_left(i,j):
            if i== 0 or j==0:
                return vN_right(i,j)
            if i == (n_points_length-1):# or j==(n_points_draft-1):
                return vN_right(i,j)
            else:
                return 8 + n_points_length*n_points_draft+(i-1)*(n_points_draft-1)+j-1
        for i in range(n_points_length-1):
            for j in range(n_points_draft-1):
                if i < old_div(n_points_length,2):
                    facets.append([[vN_left(i,j),vN_left(i+1,j+1),vN_left(i+1,j)]])
                    facetFlags.append(boundaryTags['obstacle'])
                    facets.append([[vN_left(i,j),vN_left(i,j+1),vN_left(i+1,j+1)]])
                    facetFlags.append(boundaryTags['obstacle'])
                else:
                    facets.append([[vN_left(i,j),vN_left(i,j+1),vN_left(i+1,j)]])
                    facetFlags.append(boundaryTags['obstacle'])
                    facets.append([[vN_left(i,j+1),vN_left(i+1,j+1),vN_left(i+1,j)]])
                    facetFlags.append(boundaryTags['obstacle'])
        topFacet=[]
        for i in range(n_points_length):
            topFacet.append(vN_right(i,n_points_draft-1))
        for i in range(n_points_length-2,0,-1):
            topFacet.append(vN_left(i,n_points_draft-1))
        facets.append([topFacet])
        facetFlags.append(boundaryTags['obstacle'])
        #for v in vertices: print v
        #for f in facets: print f
        holes.append(hull_center)
    if vessel is 'cube':
        nStart = len(vertices)
        vertices.append([hull_center[0] - 0.5*hull_length,
                         hull_center[1] - 0.5*hull_beam,
                         hull_center[2] - 0.5*hull_draft])
        vertexFlags.append(boundaryTags['obstacle'])
        vertices.append([hull_center[0] - 0.5*hull_length,
                         hull_center[1] + 0.5*hull_beam,
                         hull_center[2] - 0.5*hull_draft])
        vertexFlags.append(boundaryTags['obstacle'])
        vertices.append([hull_center[0] + 0.5*hull_length,
                         hull_center[1] + 0.5*hull_beam,
                         hull_center[2] - 0.5*hull_draft])
        vertexFlags.append(boundaryTags['obstacle'])
        vertices.append([hull_center[0] + 0.5*hull_length,
                         hull_center[1] - 0.5*hull_beam,
                         hull_center[2] - 0.5*hull_draft])
        vertexFlags.append(boundaryTags['obstacle'])
        vertices.append([hull_center[0] - 0.5*hull_length,
                         hull_center[1] - 0.5*hull_beam,
                         hull_center[2] + 0.5*hull_draft])
        vertexFlags.append(boundaryTags['obstacle'])
        vertices.append([hull_center[0] - 0.5*hull_length,
                         hull_center[1] + 0.5*hull_beam,
                         hull_center[2] + 0.5*hull_draft])
        vertexFlags.append(boundaryTags['obstacle'])
        vertices.append([hull_center[0] + 0.5*hull_length,
                         hull_center[1] + 0.5*hull_beam,
                         hull_center[2] + 0.5*hull_draft])
        vertexFlags.append(boundaryTags['obstacle'])
        vertices.append([hull_center[0] + 0.5*hull_length,
                         hull_center[1] - 0.5*hull_beam,
                         hull_center[2] + 0.5*hull_draft])
        vertexFlags.append(boundaryTags['obstacle'])
        facets.append([[nStart,nStart+1,nStart+2,nStart+3]])#1
        facetFlags.append(boundaryTags['obstacle'])
        facets.append([[nStart,nStart+1,nStart+5,nStart+4]])#2
        facetFlags.append(boundaryTags['obstacle'])
        facets.append([[nStart+1,nStart+2,nStart+6,nStart+5]])#3
        facetFlags.append(boundaryTags['obstacle'])
        facets.append([[nStart+2,nStart+3,nStart+7,nStart+6]])#4
        facetFlags.append(boundaryTags['obstacle'])
        facets.append([[nStart+3,nStart,nStart+4,nStart+7]])#5
        facetFlags.append(boundaryTags['obstacle'])
        facets.append([[nStart+4,nStart+5,nStart+6,nStart+7]])#6
        facetFlags.append(boundaryTags['obstacle'])
        holes.append(hull_center)
    domain = Domain.PiecewiseLinearComplexDomain(vertices=vertices,
                                                 vertexFlags=vertexFlags,
                                                 facets=facets,
                                                 facetFlags=facetFlags,
                                                 regions=regions,
                                                 regionFlags=regionFlags,
                                                 holes=holes)
    #go ahead and add a boundary tags member 
    domain.boundaryTags = boundaryTags
    from proteus import Comm
    comm = Comm.get()
    if vessel:
        polyfile="mesh_"+vessel
    else:
        polyfile="meshNoVessel"
    if comm.isMaster():
        domain.writePoly(polyfile)
    else:
        domain.polyfile=polyfile

restrictFineSolutionToAllMeshes=False
parallelPartitioningType = MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0

#steady-state so no initial conditions
initialConditions = None
#use sparse diffusion representation
sd=True
#identity tensor for defining analytical heterogeneity functions
Ident = numpy.zeros((nd,nd),'d')
Ident[0,0]=1.0; Ident[1,1] = 1.0; Ident[2,2]=1.0

#for computing exact 'Darcy' velocity
class velEx(object):
    def __init__(self,duex,aex):
        self.duex = duex
        self.aex = aex
    def uOfX(self,X):
        du = self.duex.duOfX(X)
        A  = numpy.reshape(self.aex(X),(3,3))
        return -numpy.dot(A,du)
    def uOfXT(self,X,T):
        return self.uOfX(X)


##################################################
#define coefficients a(x)=[a_{ij}] i,j=0,2, right hand side f(x)  and analytical solution u(x)
#u = x*x + y*y + z*z, a_00 = x + 5, a_11 = y + 5.0 + a_22 = z + 10.0
#f = -2*x -2*(5+x) -2*y-2*(5+y) -2*z-2*(10+z)
#

def a5(x):
    return numpy.array([[x[0] + 5.0,0.0,0.0],[0.0,x[1] + 5.0,0.0],[0.0,0.0,x[2]+10.0]],'d')
def f5(x):
    return -2.0*x[0] -2*(5.+x[0]) -2.*x[1]-2.*(5.+x[1]) -2.*x[2]-2.*(10+x[2])
#'manufactured' analytical solution
class u5Ex(object):
    def __init__(self):
        pass
    def uOfX(self,x):
        return x[0]**2+x[1]**2+x[2]**2
    def uOfXT(self,X,T):
        return self.uOfX(X)
    def duOfX(self,X):
        du = 2.0*numpy.reshape(X[0:3],(3,))
        return du
    def duOfXT(self,X,T):
        return self.duOfX(X)

#dirichlet boundary condition functions on (x=0,y,z), (x,y=0,z), (x,y=1,z), (x,y,z=0), (x,y,z=1)
def getDBC5(x,flag):
    if flag in [boundaryTags['bottom'],boundaryTags['top'],boundaryTags['front'],boundaryTags['back'],boundaryTags['left']]:
        return lambda x,t: u5Ex().uOfXT(x,t)
def getAdvFluxBC5(x,flag):
    pass
#specify flux on (x=1,y,z)
def getDiffFluxBC5(x,flag):
    if flag == boundaryTags['right']:
        n = numpy.zeros((nd,),'d'); n[0]=1.0
        return lambda x,t: numpy.dot(velEx(u5Ex(),a5).uOfXT(x,t),n)
    elif flag == 0:
        return lambda x,t: 0.0

#dirichlet boundary condition functions on (x=0,y,z), (x,y=0,z), (x,y=1,z), (x,y,z=0), (x,y,z=1)
# def getDBC5(x,flag):
#     if x[0] in [0.0] or x[1] in [0.0,1.0] or x[2] in [0.0,1.0]:
#         return lambda x,t: u5Ex().uOfXT(x,t)
# def getAdvFluxBC5(x,flag):
#     pass
# #specify flux on (x=1,y,z)
# def getDiffFluxBC5(x,flag):
#     if x[0] == 1.0:
#         n = numpy.zeros((nd,),'d'); n[0]=1.0
#         return lambda x,t: numpy.dot(velEx(u5Ex(),a5).uOfXT(x,t),n)
#     if not (x[0] in [0.0] or x[1] in [0.0,1.0] or x[2] in [0.0,1.0]):
#         return lambda x,t: 0.0

#store a,f in dictionaries since coefficients class allows for one entry per component
aOfX = {0:a5}; fOfX = {0:f5}

#one component
nc = 1
#load analytical solution, dirichlet conditions, flux boundary conditions into the expected variables
analyticalSolution = {0:u5Ex()}
analyticalSolutionVelocity = {0:velEx(analyticalSolution[0],aOfX[0])}
#
dirichletConditions = {0:getDBC5}
advectiveFluxBoundaryConditions =  {0:getAdvFluxBC5}
diffusiveFluxBoundaryConditions = {0:{0:getDiffFluxBC5}}
fluxBoundaryConditions = {0:'setFlow'} #options are 'setFlow','noFlow','mixedFlow'


#equation coefficient names
coefficients = TransportCoefficients.PoissonEquationCoefficients(aOfX,fOfX,nc,nd)
#
coefficients.variableNames=['u0']

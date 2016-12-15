
# coding: utf-8

# # LA Test Problem
#
# This notebook demonstrates how to
# - Set up a scalar, nonlinear hyperbolic PDE
# - Set up an initial-boundary value problem
# - Choose a particular set of numerics
# - Post-process the solution

# In[ ]:

#get_ipython().magic(u'matplotlib notebook')
import  math
import  numpy as  np
from proteus.iproteus import *
Profiling.logLevel=5
Profiling.verbose=True


# # Defining an equation
# The equation we want to solve is
# \begin{equation*}
# m_t + \nabla \cdot \mathbf{f} = 0
# \end{equation*}
# where $u$ is the unknown solution and the coefficients have the specific  forms
# \begin{align}
# m(u) &= u \\
# f(u) &= (sin(u), cos(u))
# \end{align}

# In[ ]:

class LA(TransportCoefficients.TC_base):
    """
    The coefficients of the linear advection-diffusion equation
    """
    def __init__(self):
        TransportCoefficients.TC_base.__init__(self,
                                               nc=1, #number of components
                                               variableNames=['u'],
                                               mass      = {0:{0:'linear'}},
                                               advection = {0:{0:'nonlinear'}})#,
                                               #diffusion = {0:{0:{0:'constant'}}},
                                               #potential = {0:{0:'u'}})

    def evaluate(self,t,c):
        c[('m',0)][:]         = c[('u',0)]
        c[('dm',0,0)][:]      = 1.0
        c[('f',0)][...,0]     = -2.0*math.pi*c['x'][...,1]*c[('u',0)]
        c[('f',0)][...,1]     =  2.0*math.pi*c['x'][...,0]*c[('u',0)]
        c[('df',0,0)][...,0]  = -2.0*math.pi*c['x'][...,1]
        c[('df',0,0)][...,1]  =  2.0*math.pi*c['x'][...,0]
        #c[('df_sge',0,0)] = c[('df',0,0)]
        #c[('da_sge',0,0,0)] = c[('da',0,0,0)]
        #c[('dphi_sge',0,0)] = c[('dm',0,0)]
        #c[('grad(phi)_sge',0)] = c[('grad(u)',0)]

# # Physics

# In[ ]:

from proteus import default_p as p
#physics
p.name = "la"
p.nd = 2; #Two dimensions
box=Domain.RectangularDomain(L=(2.0,2.0),
                             x=(-1.0,-1.0),
                             name="la");
box.writePoly("la")
p.domain=Domain.PlanarStraightLineGraphDomain(fileprefix="la")

p.T=1.0

p.coefficients=LA()

def getDBC(x,flag):
    return None
#    if flag in ['left','right','top','bottom']:
#        return lambda x,t: 0.0

def getAFBC(x,flag):
    if flag in ['left','right','top','bottom']:
        return lambda x,t: 0.0

p.dirichletConditions = {0:getDBC}
p.advectiveFluxBoundaryConditions = {0:getAFBC}

class IC:
    def __init__(self):
        self.x0=[0.3,0.0]
        self.r0=0.25
    def uOfXT(self,x,t):
        return 0.5*(1.0 - math.tanh(((x[0]-self.x0[0])**2 +
                                     (x[1]-self.x0[1])**2)/self.r0**2 - 1.0))

class IC_3B:
    def __init__(self,
                 xd=(0,0.5),
                 xc=(0,-0.5),
                 xh=(-0.5,0),
                 r0=0.3):
        self.xd=np.array(xd+(0.0,))
        self.xc=np.array(xc+(0.0,))
        self.xh=np.array(xh+(0.0,))
        self.r0=r0
    def g(self, r):
        from math import cos, pi
        return (1.+cos(pi*min(r/self.r0,1.)))/4.0
    def uOfXT(self,x, t):
        from math import fabs
        import numpy as  np
        from numpy.linalg import norm
        x=np.array(x)
        rd = norm(x-self.xd)
        rc = norm(x-self.xc)
        rh = norm(x-self.xh)
        if ( rd <= self.r0 and
             (fabs(x[0]) >= 0.05 or x[1] >= 0.7)):
            return 1.0
        elif rc <= self.r0:
            return 1. - rc/self.r0
        elif rh <= self.r0:
            return self.g(rh)
        else:
            return 0.0

p.initialConditions  = {0:IC_3B()}
#p.initialConditions  = {0:IC()}


# # Numerics

# In[ ]:

from proteus import default_n as n
from proteus.Gauges  import LineGauges
import numpy as np
import proteus as pr
import sys
timeOrder = int(sys.argv[1])
spaceOrder = int(sys.argv[2])
he = 0.02*0.5**float(sys.argv[3])
n.triangleOptions="VApq33Dena%8.8f" % (he**2/2.0,)
if timeOrder is 1:
    n.timeIntegration = pr.TimeIntegration.LinearSSPRKintegration
    n.nStagesTime = 1
    #n.timeIntegration = pr.TimeIntegration.BackwardEuler_cfl
elif timeOrder is 2:
    n.timeIntegration = pr.TimeIntegration.LinearSSPRKintegration
    #n.timeIntegration = pr.TimeIntegration.VBDF
    n.timeOrder = 2
    n.nStagesTime = 2
elif timeOrder is 3:
    n.timeIntegration = pr.TimeIntegration.LinearSSPRKintegration
    #n.timeIntegration = pr.TimeIntegration.VBDF
    n.timeOrder = 3
    n.nStagesTime = 3
n.stepController = pr.StepControl.Min_dt_RKcontroller
#n.stepController = pr.StepControl.Min_dt_controller
n.runCFL=0.33

if  spaceOrder is 1:
    n.femSpaces = {0:pr.FemTools.C0_AffineLinearOnSimplexWithNodalBasis}
elif spaceOrder is 2:
    n.femSpaces = {0:pr.FemTools.C0_AffineQuadraticOnSimplexWithNodalBasis}
n.elementQuadrature = pr.Quadrature.SimplexLobattoQuadrature(p.nd,1)
n.elementBoundaryQuadrature = pr.Quadrature.SimplexLobattoQuadrature(p.nd-1,1)
n.subgridError = pr.SubgridError.AdvectionLag_ASGS(p.coefficients,
                                                   p.nd,
                                                   lag=False)
#n.subgridError=None
n.shockCapturing = pr.ShockCapturing.ResGradQuadDelayLag_SC(p.coefficients,
                                                            p.nd,
                                                            shockCapturingFactor=0.9,
                                                            nStepsToDelay=1,
                                                            lag=False)
#n.shockCapturing=None
n.numericalFluxType = pr.NumericalFlux.Advection_DiagonalUpwind_exterior
n.tnList=[float(i)/10.0 for i in range(11)]
n.matrix = pr.LinearAlgebraTools.SparseMatrix
n.multilevelLinearSolver = pr.LinearSolvers.LU#KSP_petsc4py
n.l_atol_res = 0.0
n.linTolFac = 0.001
n.tolFac = 0.0
n.nl_atol_res = 1.0e-4
n.maxNonlinearIts =500
n.maxLineSearches =0
n.parallelPartitioningType = pr.MeshTools.MeshParallelPartitioningTypes.node
n.nLayersOfOverlapForParallel = 0
n.periodicDirichletConditions = None
#lineGauges  = LineGauges(gauges=((('u',),#fields in gauge set
#                                  (#lines for these fields
#                                      ((-2.0, -2.5, 0.0),(2.0, 1.5, 0.0)),
#                                  ),#end  lines
#                              ),#end gauge set
#                             ),#end gauges
#                         fileName="profile.csv")
#n.auxiliaryVariables=[lineGauges]


# # Operator Splitting
# Trivial since this is a scalar PDE

# In[ ]:

from proteus import default_s,default_so
so = default_so
p.name += "k{0}_p{1}_h{2}".format(timeOrder,spaceOrder,he)

so.name = p.name
so.sList=[default_s]
so.tnList = n.tnList
from proteus.Archiver import ArchiveFlags
so.archiveFlag = ArchiveFlags.EVERY_USER_STEP
so.systemStepControllerType = SplitOperator.Sequential_MinAdaptiveModelStep

# # Initialize Numerical Solution Object

# In[ ]:

ns = NumericalSolution.NS_base(so,[p],[n],so.sList,opts)


# # Calculate Solution

# In[ ]:
failed = ns.calculateSolution('la_k{0}_p{1}_h{2}'.format(timeOrder,spaceOrder,he))
assert(not failed)

print ns.modelList[0].levelModelList[-1].u[0].dof.max()
# In[ ]:
plots=False
if spaceOrder == 1 and plots:
    from matplotlib import pyplot as plt
    fig = plt.figure()
    x = ns.modelList[0].levelModelList[-1].mesh.nodeArray[:,0]
    y = ns.modelList[0].levelModelList[-1].mesh.nodeArray[:,1]
    triangles = ns.modelList[0].levelModelList[-1].mesh.elementNodesArray
    u = ns.modelList[0].levelModelList[-1].u[0].dof
    plt.axis('equal')
    plt.tricontour(x,y,triangles,u,np.linspace(math.pi/4.0,14.0*math.pi/4.0,51)[1:],colors='k')
    plt.savefig(p.name+"_Contour.png")

    # In[ ]:

    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    from matplotlib  import pyplot as plt

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    x = ns.modelList[0].levelModelList[-1].mesh.nodeArray[:,0]
    y = ns.modelList[0].levelModelList[-1].mesh.nodeArray[:,1]
    triangles = ns.modelList[0].levelModelList[-1].mesh.elementNodesArray
    u = ns.modelList[0].levelModelList[-1].u[0].dof
    ax.plot_trisurf(x,y,u,cmap=cm.jet, linewidth=0.2)
    plt.savefig(p.name+"_Surface.png")

# In[ ]:
#fig = plt.figure()
#column = np.loadtxt("profile.csv",skiprows=1,delimiter=",")
#z = np.linspace(0.0,math.sqrt(2.0)*4.0,column[-1,1:].shape[0])
#plt.plot(z,column[-1,1:])
#fig.savefig("laLine.png")

# In[ ]:

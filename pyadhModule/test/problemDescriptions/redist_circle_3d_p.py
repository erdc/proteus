from pyadh import *
from pyadh.default_p import *
from math import *

## \page Tests Test Problems 
# \ref redist_circle_3d_p.py "Hamilton Jacobi LS redistancing (circle)"
# 

##\ingroup test
# \file redist_circle_3d_p.py
# @{
#  \brief Solve the level set redistancing equation for a fixed interface. 
#
#  The redistancing function is the steady state solution of
# \f{eqnarray*}
# \phi_t + S\|\nabla \phi\| &=& S \\ 
# \Omega &=& [0,1] \times [0,1] \\
#  S &=& 2\left(H(\phi^0)-\frac{1}{2}\right) \\
# \f}
# where \f$H(u)\f$ is the Heaviside function. In the implementation, a smoothed
# version is used
# \f{eqnarray*}
# H_{\epsilon}(u) = \left\{ \begin{array}{ll}
# 0 & u < -\epsilon \\
# 1 & u > \epsilon  \\
# \frac{1}{2}(1 + \frac{u}{\epsilon} + \frac{1}{\pi}\sin(\frac{\pi u}{\epsilon})) & \mbox{otherwise}
# \end{array} \right.
# \f}
#
# In this case, the fixed interface is a circle. The initial condition is the
# squared signed distance to the circle. By default, the initial condition is
#
# \f{eqnarray*}
# \phi^0 &=& \left(x-\frac{1}{2}\right)^2 + \left(y - \frac{1}{2}\right)^2 - \frac{1}{16} 
# \f}
#
# \image html  save_redist_circle_3d_c0p1_phi0.jpg "initial condition"
# \image latex save_redist_circle_3d_c0p1_phi0.eps "initial condition"
# \image html  save_redist_circle_3d_c0p1_phi.jpg "signed distance solution"
# \image latex save_redist_circle_3d_c0p1_phi.eps "signed distance solution"
#
#
#

nd = 3

#mwf hack test unstructured mesh
#polyfile = "unstSquare"


class ConstantVelocityCircleSquared:
    def __init__(self,radius=0.1,b=[1.,1.],startX=0.25,startY=0.5,startZ=0.5):
        self.radius = radius
        self.b      = b
        self.startX = startX
        self.startY = startY
        self.startZ = startZ
    def uOfXT(self,x,t):
        centerX = self.b[0]*0.0 + self.startX
        centerY = self.b[1]*0.0 + self.startY
        centerZ = self.b[1]*0.0 + self.startZ
        return (x[0]-centerX)**2 + (x[1]-centerY)**2 + (x[2]-centerZ)**2-self.radius**2

class ConstantVelocityCircle:
    def __init__(self,radius=0.1,b=[1.,1.],startX=0.25,startY=0.5,startZ=0.5):
        self.radius = radius
        self.b      = b
        self.startX = startX
        self.startY = startY
        self.startZ = startZ
    def uOfXT(self,x,t):
        centerX = self.b[0]*0.0 + self.startX
        centerY = self.b[1]*0.0 + self.startY
        centerZ = self.b[1]*0.0 + self.startZ
        return sqrt((x[0]-centerX)**2 + (x[1]-centerY)**2) + (x[2]-centerZ)**2-self.radius

x0 = 0.50
y0 = 0.50
z0 = 0.50
r0 = 1./8.
b0 = Numeric.array([-1.,0,0.0],Numeric.Float)

analyticalSolution = {0:ConstantVelocityCircle(r0,b0,x0,y0,z0)}

def Heaviside(p,eps=1.0e-2):
    if p < -eps:
        return 0.0
    if p > eps:
        return 1.0
    return 0.5*(1.+p/eps + 1./math.pi*math.sin(math.pi*p/eps))
def Se(p,eps=1.0e-2):
    return 2.0*(Heaviside(p,eps)-0.5)

class RedistanceLevelSetOld(TransportCoefficients.TC_base):
    from pyadh.ctransportCoefficients import redistanceLevelSetCoefficientsEvaluate
    def __init__(self,eps=1.0e0):
        mass={0:{0:'linear'}}
#mwf need for outflow bcs 
#        advection={0:{0:'linear'}}
        advection={}
        diffusion={}
        potential={}
        reaction={0:{0:'constant'}}
        hamiltonian={0:{0:'nonlinear'}}
        self.phi0xt = {0:ConstantVelocityCircleSquared(r0,b0,x0)}
        self.phi_q = None
        self.phi_ebq = None
        TransportCoefficients.TC_base.__init__(self,
                                               1,
                                               mass,
                                               advection,
                                               diffusion,
                                               potential,
                                               reaction,
                                               hamiltonian)
        self.eps = eps
    def evaluate(self,t,c):
        phi = None
        if self.phi_q == None:
            #assume first c through is q for now
            self.phi_q= Numeric.zeros(c[('m',0)].shape,Numeric.Float)
            #evaulate phi_q
            for i in range(len(c[('m',0)].flat)):
                self.phi_q.flat[i]=self.phi0xt[0].uOfXT(c['x'].flat[3*i:3*(i+1)],0.)
            #end i
            phi = self.phi_q
        elif self.phi_q.shape == c[('m',0)].shape:
            phi = self.phi_q
        elif self.phi_ebq == None:
            self.phi_ebq= Numeric.zeros(c[('m',0)].shape,Numeric.Float)
            #evaluate at ebq
            for i in range(len(c[('m',0)].flat)):
                self.phi_ebq.flat[i]=self.phi0xt[0].uOfXT(c['x'].flat[3*i:3*(i+1)],0.)
            #end i
            phi = self.phi_ebq
        else:
            phi = self.phi_ebq
        #end
        self.redistanceLevelSetCoefficientsEvaluate(self.eps,
                                                    phi,
                                                    c[('u',0)],
                                                    c[('grad(u)',0)],
                                                    c[('m',0)],
                                                    c[('dm',0,0)],
                                                    c[('H',0)],
                                                    c[('dH',0,0)],
                                                    c[('r',0)])
        
        
##         if ('m',0) in c.keys():
##             c[('m',0)].flat[:] = c[('u',0)].flat[:]
##         if ('dm',0,0) in c.keys():
##             c[('dm',0,0)].flat[:] = 1.0
##         if ('f',0) in c.keys():
##              c[('f',0)].flat[:] = 0.0
##         if ('df',0,0) in c.keys():
##             c[('df',0,0)].flat[:] = 0.0
##         for i in range(len(c[('u',0)].flat)):
##             normgrad=sqrt(Numeric.dot(c[('grad(u)',0)].flat[nd*i:nd*(i+1)],
##                                        c[('grad(u)',0)].flat[nd*i:nd*(i+1)]))
             
##             Sei = Se(phi.flat[i])
##             c[('H',0)].flat[i]=Sei*normgrad
##             c[('dH',0,0)].flat[nd*i:nd*(i+1)]=Sei*c[('grad(u)',0)].flat[nd*i:nd*(i+1)]/(normgrad+1.0e-10)
##             c[('r',0)].flat[i]=-Sei
##         #end for
    #end def
                                          

coefficients = RedistanceLevelSet(epsFact=1.5,u0=ConstantVelocityCircleSquared(r0,b0,x0,y0))

coefficients.variableNames=['u']

#now define the Dirichlet boundary conditions

def getDBC(x):
    pass
    #if (x[1] == 0.0):
    #    return lambda x,t: 0.0
    #if (x[0] == 0.0 or
    #    x[0] == 1.0 or
    #    x[1] == 0.0 or
    #    x[1] == 1.0):
    #    return lambda x,t: 0.0
    
dirichletConditions = {0:getDBC}

def getZeroLSDBC(vt):
        vt.dirichletNodeSetList[0] = []
        vt.dirichletGlobalNodeSet[0]=set()
        vt.dirichletValues[0]={}
        for eN in range(vt.mesh.nElements_global):
            vt.dirichletNodeSetList[0].append(set())
            signU = 0
            j0=0
            while ((signU == 0) and
                   (j0 < vt.nDOF_trial_element[0])):
                J0 = vt.u[0].femSpace.dofMap.l2g[eN,j0]
                if vt.u[0].dof[J0] < 0.0:
                    signU = -1
                elif  vt.u[0].dof[J0] > 0.0:
                    signU = 1
                else:
                    vt.dirichletNodeSetList[0][eN].add(j0)
                    vt.dirichletValues[0][(eN,j0)]=float(vt.u[0].dof[J0])
                    vt.dirichletGlobalNodeSet[0].add(J0)
                j0 += 1
            for j in range(j0,vt.nDOF_trial_element[0]):
                J = vt.u[0].femSpace.dofMap.l2g[eN,j]
                if (((vt.u[0].dof[J] < 0.0) and
                     (signU == 1)) or
                    ((vt.u[0].dof[J] > 0.0) and
                     (signU == -1))):
                    for jj in range(vt.nDOF_trial_element[0]):
                        JJ = vt.u[0].femSpace.dofMap.l2g[eN,jj]
                        vt.dirichletNodeSetList[0][eN].add(jj)
                        vt.dirichletValues[0][(eN,jj)]=float(vt.u[0].dof[JJ])
                        vt.dirichletGlobalNodeSet[0].add(JJ)
                    break
                elif (vt.u[0].dof[J] == 0.0):
                    vt.dirichletNodeSetList[0][eN].add(j)
                    vt.dirichletValues[0][(eN,j)]=float(vt.u[0].dof[J])
                    vt.dirichletGlobalNodeSet[0].add(J)
        for eN in range(vt.mesh.nElements_global):
            for j in range(vt.nDOF_trial_element[0]):
                J = vt.u[0].femSpace.dofMap.l2g[eN,j]
                if J in vt.dirichletGlobalNodeSet[0]:
                    vt.dirichletNodeSetList[0][eN].add(j)
                    vt.dirichletValues[0][(eN,j)]=float(vt.u[0].dof[J])

#weakDirichletConditions = {0:getZeroLSDBC}
#weakDirichletConditions = {0:coefficients.setZeroLSweakDirichletBCs2}
weakDirichletConditions = {0:coefficients.setZeroLSweakDirichletBCs}
#weakDirichletConditions = None

initialConditions  = {0:ConstantVelocityCircleSquared(r0,b0,x0,y0)}

fluxBoundaryConditions = {0:'noFlow'}
#fluxBoundaryConditions = {0:'outFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

T = 0.75e1

## @}

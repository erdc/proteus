from pyadh import *
from pyadh.default_p import *

nd = 2
Ident = Numeric.zeros((nd,nd),Numeric.Float)
Ident[0,0]=1.0; Ident[1,1] = 1.0
gravity = Numeric.zeros((nd,),Numeric.Float)
gravity[1] = -1.0

rho = 1.0; 
K   = 1.0; 
zB   = 0.0; dp   = 1.0
waterTable = 0.5
psiB       = rho*(waterTable-zB)
def Heaviside(p,eps=1.0e-2):
    if p <= -eps:
        return 0.0
    if p > eps:
        return 1.0
    return 0.5*(1.+p/eps + 1./math.pi*math.sin(math.pi*p/eps))

class DarcyFlowIM(TransportCoefficients.TC_base):
     from pyadh.ctransportCoefficients import darcySharpInterfaceFlowEvaluate
     from pyadh.ctransportCoefficients import darcySharpInterfaceFlowImEvaluate
     def __init__(self,K,rho,gravity,LS_model=1):
          nc = 1
          self.K   = K
          self.rho = rho
          self.gravity = gravity
          self.K0  = 0.0
          self.rho0= 0.0
          self.LS_model = LS_model
          self.heaviEps = 1.0e-1
          mass     = {}
          advection = {}
          diffusion = {}
          potential = {}
          reaction  = {}
          hamiltonian = {}
          #Darcy mass balance
          i = 0
          mass[i]     = {i : 'linear'}
          diffusion[i] = {i : {i:'constant'}}
          potential[i] = {i : 'u'}
          advection[i] = {i : 'constant'}
          reaction[i] = {i : 'constant'}
          TransportCoefficients.TC_base.__init__(self,
                                                 nc,
                                                 mass,
                                                 advection,
                                                 diffusion,
                                                 potential,
                                                 reaction,
                                                 hamiltonian)
          self.useC = True
          self.q_uls  = None
          self.ebq_uls= None
          self.udof_uls=None
     #end init
     def attachModels(self,modelList):
          self.q_uls = modelList[self.LS_model].q[('u',0)]
          self.ebq_uls = modelList[self.LS_model].ebq[('u',0)]
          self.uls= modelList[self.LS_model].u[0]
     def initializeElementQuadrature(self,t,cq):
          #initialize so it can run without a flow model
          self.q_uls = Numeric.ones(cq[('u',0)].shape,Numeric.Float)
     def initializeElementBoundaryQuadrature(self,t,cebq,cebq_global):
          self.ebq_uls = Numeric.ones(cebq[('u',0)].shape,Numeric.Float)
    
     def evaluate(self,t,c):
          if self.q_uls.shape == c[('u',0)].shape:
               uls = self.q_uls
          else:
               uls = self.ebq_uls
          #ImEvaluate doesn't smear density
          self.darcySharpInterfaceFlowImEvaluate(self.K0,self.rho0,
                                                 self.K,self.rho,
                                                 self.heaviEps,
                                                 self.gravity,
                                                 c[('u',0)],
                                                 c[('grad(u)',0)],
                                                 uls,#level set function
                                                 c[('phi',0)],#potential
                                                 c[('a',0,0)],
                                                 c[('f',0)],
                                                 c[('r',0)],
                                                 c[('m',0)],
                                                 c[('dphi',0,0)],
                                                 c[('da',0,0,0)],
                                                 c[('df',0,0)],
                                                 c[('dr',0,0)],
                                                 c[('dm',0,0)])
          #ci = 0
          ##c[('m',ci)].flat[:] = c[('u',ci)].flat[:]
          ##c[('dm',ci,ci)].flat[:] = 1.0
          #c[('phi',ci)].flat[:] = c[('u',ci)].flat[:]
          #c[('r',ci)].flat[:] = 0.0
          #for i in range(len(c[('u',ci)].flat)):
          #     c[('a',ci,ci)].flat[nd*nd*i:nd*nd*(i+1)] = Ident.flat[:]*self.K*self.rho
          #     c[('f',ci)].flat[nd*i:nd*(i+1)]=Heaviside(uls.flat[i],eps=0.0)*gravity.flat[:]*self.K*(self.rho**2)
     #eval
#           
class HydrostaticEx:
    def __init__(self,rho=1.0,waterTable=0.5,psiB=0.0,zB=0.0):
        self.waterTable=waterTable
        self.rho   =rho
        self.psiB  =psiB
        self.zB    =zB
        self.psiI  =self.psiB + self.rho*(self.zB-self.waterTable)
        #mwf debug
        #print """psiB=%g psiI=%g\n""" % (self.psiB,self.psiI)
    def uOfX(self,X):
        if X[1] <= self.waterTable + self.zB:
            return self.psiB + self.rho*(self.zB+gravity[1]*X[1])
        else:
            return 0.0
    #end
    def uOfXT(self,X,t):
        return self.uOfX(X)
    #end
    def duOfX(self,X):
        du = Numeric.zeros((nd,),Numeric.Float)
        return du
#end HydrostaticEx
class HydrostaticVelEx:
    def __init__(self):
        pass
    def uOfX(self,X):
        du = Numeric.zeros((nd,),Numeric.Float)
        return du
    #end
    def uOfXT(self,X,t):
        return self.uOfX(X)
    #end
    def duOfX(self,X):
        du = Numeric.zeros((nd,),Numeric.Float)
        return du
#end HydrostaticEx

class PerturbedWaterTable_psi:
    def __init__(self,rho=1.0,waterTable=0.5,psiB=0.0,zB=0.0,dp=1.0):
        self.rho   =rho
        self.zB    =zB
        self.dp    = dp
        self.waterTablePert=lambda x : self.dp*x*(1.0-x) + waterTable + self.zB
        self.psiI  = lambda x : 0.0
    def uOfX(self,X):
        zP = self.waterTablePert(X[0])
        if X[1] <= zP:
            return 0.0 + self.rho*(zP+gravity[1]*X[1])
        else:
            return 0.
    #end
    def uOfXT(self,X,t):
        return self.uOfX(X)
    #end
#end PerturbedWaterTable_psi

def getDBCHorizFlow(X):
    #looks like works in relax to steady state
    if X[0] in [0.0,1.0]:
        return lambda x,t: HydrostaticEx(rho,waterTable,psiB,zB).uOfXT(x,t)
    
def getDBCBottom(X):
    #looks like works too for relax to steady state
    if X[1] in [0.0]:
        return lambda x,t: HydrostaticEx(rho,waterTable,psiB,zB).uOfXT(x,t)
#dirichlet
def getDBCTopCorner(X):
    #looks like works for relax to steady state
    if X[0] in [1.0] and X[1] in [1.0]:
        return lambda x,t: HydrostaticEx(rho,waterTable,psiB,zB).uOfXT(x,t)
#dirichlet
def getDBCHydro(X):
    #looks like works in relax to steady state
    if X[0] in [0.0,1.0]:
        return lambda x,t: HydrostaticSinglePhase_psi(rho).uOfXT(x,t)
#dirichlet

analyticalSolution = {0:HydrostaticEx(rho,waterTable,psiB,zB)}
analyticalSolutionVelocity = {0:HydrostaticVelEx()}

initialConditions  = None #{0:PerturbedWaterTable_psi(rho,waterTable,psiB,zB,dp)}
dirichletConditions= {0:getDBCBottom}#{0:getDBCHorizFlow}

#has to have same numbering and space for LS and Flow model!
def getZeroDBC(vt):
     vt.dirichletNodeSetList[0] = []
     vt.dirichletGlobalNodeSet[0]=set()
     vt.dirichletValues[0]={}
     uls = vt.coefficients.uls
     for eN in range(uls.femSpace.mesh.nElements_global):
          vt.dirichletNodeSetList[0].append(set())
          for j in range(uls.femSpace.referenceFiniteElement.localFunctionSpace.dim):
               J = uls.femSpace.dofMap.l2g[eN,j]
               if uls.dof[J] <= 0.0:
                    #mwf debug
                    #print """setting zero p eN=%d j=%d J=%d \n""" % (eN,j,J)
                    vt.dirichletNodeSetList[0][eN].add(j)
                    vt.dirichletValues[0][(eN,j)]=0.0# float(uls.dof[J])
                    vt.dirichletGlobalNodeSet[0].add(J)

weakDirichletConditions = {0:getZeroDBC}
#weakDirichletConditions = None


coefficients = DarcyFlowIM(K,rho,gravity)

coefficients.variableNames=['psi']
fluxBoundaryConditions = {0:'noFlow'}
advectiveFluxBoundaryConditions =  {}
diffusiveFluxBoundaryConditions = {0:{}}




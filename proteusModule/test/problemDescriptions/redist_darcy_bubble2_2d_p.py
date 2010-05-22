from pyadh import *
from pyadh.default_p import *
from math import *

nd = 2
class RedistanceLevelSetHack(TC_base):
    from pyadh.ctransportCoefficients import redistanceLevelSetCoefficientsEvaluate
    def __init__(self,epsFact=2.0,nModel=1,u0=None,rdModel=None):
        variableNames=['phid']
        nc=1
        mass={0:{0:'linear'}}
        hamiltonian={0:{0:'nonlinear'}}
        advection={}
        diffusion={}
        potential={}
        reaction={0:{0:'constant'}}
        TC_base.__init__(self,
                         nc,
                         mass,
                         advection,
                         diffusion,
                         potential,
                         reaction,
                         hamiltonian,
                         variableNames)
        self.nModel = nModel
        self.rdModel= rdModel
        self.epsFact=epsFact
        self.q_u0   = None
        self.ebq_u0 = None
        self.dof_u0 = None
        self.hmin   = 1.0
        self.u0 = u0
    def attachModels(self,modelList):
        self.q_u0 = modelList[self.nModel].q[('u',0)]
        self.ebq_u0 = modelList[self.nModel].ebq[('u',0)]
        self.dof_u0 = modelList[self.nModel].u[0].dof
        if self.rdModel != None:
            modelList[self.rdModel].u[0].dof = self.dof_u0 #force RD to be shallow copy?
    def initializeMesh(self,mesh):
        self.eps = self.epsFact*mesh.h
    def initializeElementQuadrature(self,t,cq):
        self.q_u0 = Numeric.zeros(cq[('u',0)].shape,Numeric.Float)
        if self.u0 != None:
            for i in range(len(cq[('u',0)].flat)):
                self.q_u0.flat[i]=self.u0.uOfXT(cq['x'].flat[3*i:3*(i+1)],0.)
        self.dof_u0 = None
    def initializeElementBoundaryQuadrature(self,t,cebq,cebq_global):
        self.ebq_u0 = Numeric.zeros(cebq[('u',0)].shape,Numeric.Float)
        if self.u0 != None:
            for i in range(len(cebq[('u',0)].flat)):
                self.ebq_u0.flat[i]=self.u0.uOfXT(cebq['x'].flat[3*i:3*(i+1)],0.)
    def getICDofs(self,cj):
        #mwf debug
        #print """In coefficients getICDofs() dof_u0= %s """ % self.dof_u0
        return self.dof_u0
    def evaluate(self,t,c):
        if c[('H',0)].shape == self.q_u0.shape:
            u0 = self.q_u0
        else:
            u0 = self.ebq_u0
        self.redistanceLevelSetCoefficientsEvaluate(self.eps,
                                                    u0,
                                                    c[('u',0)],
                                                    c[('grad(u)',0)],
                                                    c[('m',0)],
                                                    c[('dm',0,0)],
                                                    c[('H',0)],
                                                    c[('dH',0,0)],
                                                    c[('r',0)])
#end class

coefficients = RedistanceLevelSet(epsFact=1.5,nModelId=1,rdModelId=2)
#coefficients = RedistanceLevelSetHack(epsFact=1.5,nModel=1,rdModel=-1)

#now define the Dirichlet boundary conditions

def getDBC(x):
    pass
    
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
                    #level set passes through this node
                    vt.dirichletNodeSetList[0][eN].add(j0)
                    vt.dirichletValues[0][(eN,j0)]=float(vt.u[0].dof[J0])
                    vt.dirichletGlobalNodeSet[0].add(J0)
                j0 += 1
            if signU != 0:
                for j in range(j0,vt.nDOF_trial_element[0]):
                    J = vt.u[0].femSpace.dofMap.l2g[eN,j]
                    if (((vt.u[0].dof[J] < 0.0) and
                         (signU == 1)) or
                        ((vt.u[0].dof[J] > 0.0) and
                         (signU == -1))):
                        #level set passes through element interior
                        #so freeze all nodes
                        for jj in range(vt.nDOF_trial_element[0]):
                            JJ = vt.u[0].femSpace.dofMap.l2g[eN,jj]
                            vt.dirichletNodeSetList[0][eN].add(jj)
                            vt.dirichletValues[0][(eN,jj)]=vt.u[0].dof[JJ]
                            vt.dirichletGlobalNodeSet[0].add(JJ)
                        break#all nodes set so no point in continuing
                    elif (vt.u[0].dof[J] == 0.0):
                        #level set passes through this node at least
                        vt.dirichletNodeSetList[0][eN].add(j)
                        vt.dirichletValues[0][(eN,j)]=float(vt.u[0].dof[J])
                        vt.dirichletGlobalNodeSet[0].add(J)
        for eN in range(vt.mesh.nElements_global):
            for j in range(vt.nDOF_trial_element[0]):
                J = vt.u[0].femSpace.dofMap.l2g[eN,j]
                if J in vt.dirichletGlobalNodeSet[0]:
                    vt.dirichletNodeSetList[0][eN].add(j)
                    vt.dirichletValues[0][(eN,j)]=vt.u[0].dof[J]

weakDirichletConditions = {0:getZeroLSDBC}
#weakDirichletConditions = None


initialConditions  = None

fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

T = 0.75e4

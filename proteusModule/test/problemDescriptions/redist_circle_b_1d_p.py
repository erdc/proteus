from pyadh import *
from pyadh.default_p import *

nd = 1

analyticalSolutions = None
polyfile = None

class RDLevelSetCoefficients(TransportCoefficients.TC_base):
    from pyadh.ctransportCoefficients import redistanceLevelSetCoefficientsEvaluate
    def __init__(self,epsFact=2,lsModelId=0,rdModelId=1):
        variableNames=['phid']
        nc=1
        mass={0:{0:'linear'}}
        hamiltonian={0:{0:'nonlinear'}}
        #for outflow
        #advection={0:{0:'linear'}}
        advection={}
        diffusion={}
        potential={}
        reaction={0:{0:'constant'}}
        TransportCoefficients.TC_base.__init__(self,
                                               nc,
                                               mass,
                                               advection,
                                               diffusion,
                                               potential,
                                               reaction,
                                               hamiltonian,
                                               variableNames)
        self.lsModelId = lsModelId #model that has the basic level set function
        self.rdModelId = rdModelId #the redistancing model
        self.epsFact=epsFact
        self.q_u0   = None
        self.ebq_u0 = None
        self.dof_u0 = None  #level set dofs
        self.dof_uRD= None  #redistance dofs
        #mwf debug
        #print """in RD_LS init dof_u0= %s """ % self.dof_u0
    def attachModels(self,modelList):
        self.q_u0 = modelList[self.lsModelId].q[('u',0)]
        self.ebq_u0 = modelList[self.lsModelId].ebq[('u',0)]
        self.dof_u0 = modelList[self.lsModelId].u[0].dof
        #mwf debug
        #print """in RD_LS attachModels dof_u0= %s """ % self.dof_u0
        #force the RD model to have shallow copy of the LS model dofs
        modelList[self.rdModelId].u[0].dof = modelList[self.lsModelId].u[0].dof
        #or try to hold off until simulation starts (skip initial condition setting?)
        self.dof_uRD = modelList[self.rdModelId].u[0].dof
        #pick eps to be multiple of
        imin = Numeric.argmin(modelList[self.lsModelId].mesh.elementDiametersArray)
        hmin = modelList[self.lsModelId].mesh.elementDiametersArray[imin]
        #mwf debug
        print """RD_LS_1D_circle coef hmin= %g epsFact= %g eps= %g """ % (hmin,self.epsFact,
                                                                          hmin*self.epsFact)
        self.eps = self.epsFact*hmin
        #save model for doing update
    def getICDofs(self,cj):
        #mwf debug
        #print """In coefficients getICDofs() dof_u0= %s """ % self.dof_u0
        return self.dof_u0
    def evaluate(self,t,c):
        if c[('H',0)].shape == self.q_u0.shape:
            u0 = self.q_u0
        else:
            u0 = self.ebq_u0
        #mwf debug
        #print """In RDLS1D coefficients evaluate u0= %s """ % u0
        #print """In coefficients evaluate dof_uRD= %s """ % self.dof_uRD
        self.redistanceLevelSetCoefficientsEvaluate(self.eps,
                                                    u0,
                                                    c[('u',0)],
                                                    c[('grad(u)',0)],
                                                    c[('m',0)],
                                                    c[('dm',0,0)],
                                                    c[('H',0)],
                                                    c[('dH',0,0)],
                                                    c[('r',0)])
        if c.has_key(('dr',0,0)):
            #mwf
            #print """RD_LS_1D zeroing dr """ 
            c[('dr',0,0)].flat[:]=0.0                               
        #for outflow bcs
        if c.has_key(('f',0)):
            c[('f',0)].flat[:] = 0.0
        if c.has_key(('df',0,0)):
            c[('df',0,0)].flat[:] = 0.0
#end class
                                                    
coefficients = RDLevelSetCoefficients()

def getDBC(x):
    pass

def getZeroLSDBC(vt,eps=0.0):
        
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
                if vt.u[0].dof[J0] < -eps:
                    signU = -1
                elif  vt.u[0].dof[J0] > eps:
                    signU = 1
                else:
                    vt.dirichletNodeSetList[0][eN].add(j0)
                    vt.dirichletValues[0][(eN,j0)]=float(vt.u[0].dof[J0])
                    vt.dirichletGlobalNodeSet[0].add(J0)
                j0 += 1
            for j in range(j0,vt.nDOF_trial_element[0]):
                J = vt.u[0].femSpace.dofMap.l2g[eN,j]
                if (((vt.u[0].dof[J] < -eps) and
                     (signU == 1)) or
                    ((vt.u[0].dof[J] > eps) and
                     (signU == -1))):
                    for jj in range(vt.nDOF_trial_element[0]):
                        JJ = vt.u[0].femSpace.dofMap.l2g[eN,jj]
                        vt.dirichletNodeSetList[0][eN].add(jj)
                        vt.dirichletValues[0][(eN,jj)]=float(vt.u[0].dof[JJ])
                        vt.dirichletGlobalNodeSet[0].add(JJ)
                    break
                elif (abs(vt.u[0].dof[J]) <= eps):
                    vt.dirichletNodeSetList[0][eN].add(j)
                    vt.dirichletValues[0][(eN,j)]=float(vt.u[0].dof[J])
                    vt.dirichletGlobalNodeSet[0].add(J)
        for eN in range(vt.mesh.nElements_global):
            for j in range(vt.nDOF_trial_element[0]):
                J = vt.u[0].femSpace.dofMap.l2g[eN,j]
                if J in vt.dirichletGlobalNodeSet[0]:
                    vt.dirichletNodeSetList[0][eN].add(j)
                    vt.dirichletValues[0][(eN,j)]=float(vt.u[0].dof[J])

        #end for eN
        #mwf debug
        print """RDLSgetLSDBC dirichletNodeSet[0]= %s """ % vt.dirichletGlobalNodeSet[0]
#end def
weakDirichletConditions = None
#weakDirichletConditions = {0:getZeroLSDBC}


dirichletConditions = {0:getDBC}

class DummyIC:
    def __init__(self,dbc):
        self.dbc=dbc
    def uOfXT(self,x,t):
        return 0.0

initialConditions  = {0:DummyIC(getDBC)}
dummyInitialConditions = True

#fluxBoundaryConditions = {0:'outFlow'}
fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

T = 100

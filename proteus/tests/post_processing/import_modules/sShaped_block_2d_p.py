from proteus import *
from proteus.default_p import *
from . import blockDomain

"""
Heterogenous Poisson's equation in 2D on an sShaped domain
"""

#########################################################

# variables set from context
name = 'poisson_bdm1_test'

#########################################################


nd = 2
right = 1.0; top = 1.0
nblocks = [4,4]

polyfile = "blockDomain"
L,boundaryTags,regions = blockDomain.genPoly(polyfileBase=polyfile,
                                             nx=nblocks[0],
                                             ny=nblocks[1],
                                             Lx=right,Ly=top)

#left and right heads
head_left = 1.0
head_right= 0.0
initialConditions = None

Ident = numpy.zeros((nd,nd),'d')
Ident[0,0]=1.0; Ident[1,1] = 1.0


#try low conductivity block in upper right corner
hydraulicConductivities = {}
#no sources
sources = {}
def hydraulicConductivity_0(x,t):
    return numpy.array([[10.0,1.0],[1.0,5.0]])
def hydraulicConductivity_1(x,t):
    return numpy.array([[1.0e-4,0.0],[0.0,1.0e-4]])
def nosource(x,t):
    return 0.0
for j in range(nblocks[1]):
    for i in range(nblocks[0]):
        sources[regions[(i,j)][-1]] = nosource
        if j==0 and i==0:
            hydraulicConductivities[regions[(i,j)][-1]] = hydraulicConductivity_1
        elif j==1 and (i==0 or i==2 or i==3):
            hydraulicConductivities[regions[(i,j)][-1]] = hydraulicConductivity_1
        elif j==2 and (i==0 or i==2 or i==3):
            hydraulicConductivities[regions[(i,j)][-1]] = hydraulicConductivity_1
        elif j==3 and (i==2 or i==3):
            hydraulicConductivities[regions[(i,j)][-1]] = hydraulicConductivity_1
        else:
            hydraulicConductivities[regions[(i,j)][-1]] = hydraulicConductivity_0

def headBCs(x,tag):
    if tag == boundaryTags['left'] and (L[1] >= x[1] >= L[1]*0.75):
        return lambda x,t: head_left
    if tag == boundaryTags['right']:
        return lambda x,t: head_right
def noFlowBCs(x,tag):
    if tag == boundaryTags['top']:
        return lambda x,t: 0.0
    if tag == boundaryTags['bottom']:
        return lambda x,t: 0.0
    if tag == boundaryTags['left'] and x[1] < L[1]*0.75:
        return lambda x,t: 0.0

class velEx(object):
    def __init__(self,duex,aex):
        self.duex = duex
        self.aex = aex
    def uOfX(self,X):
        du = self.duex.duOfX(X)
        A  = numpy.reshape(self.aex(X),(2,2))
        return -numpy.dot(A,du)
    def uOfXT(self,X,T):
        return self.uOfX(X)



##################################################
class SinglePhaseDarcyCoefficients(TC_base):
    def __init__(self,a_types,source_types,nc=1,nd=2,
                 timeVaryingCoefficients=False):
        self.a_types = a_types
        self.source_types = source_types
        self.nd = nd
        self.timeVaryingCoefficients=timeVaryingCoefficients
        mass = {}
        advection = {}
        diffusion = {}
        potential = {}
        reaction  = {}
        hamiltonian = {}
        for i in range(nc):
            diffusion[i] = {i : {i:'constant'}}
            reaction[i]  = {i : 'constant'}
            advection[i] = {i : 'constant'} #now include for gravity type terms
            potential[i] = {i : 'u'}
        #end i
        TC_base.__init__(self,
                         nc,
                         mass,
                         advection,
                         diffusion,
                         potential,
                         reaction,
                         hamiltonian,
                         useSparseDiffusion=True,
                         sparseDiffusionTensors={0:(numpy.arange(start=0,stop=nd**2+1,step=nd,dtype='i'),
                                                  numpy.array([list(range(nd)) for row in range(nd)],dtype='i').flatten())})
    def initializeMesh(self,mesh):
        self.elementMaterialTypes = mesh.elementMaterialTypes
        self.exteriorElementBoundaryTypes =  numpy.zeros((mesh.nExteriorElementBoundaries_global),'i')
        for ebNE in range(mesh.nExteriorElementBoundaries_global):
            ebN = mesh.exteriorElementBoundariesArray[ebNE]
            eN  = mesh.elementBoundaryElementsArray[ebN,0]
            self.exteriorElementBoundaryTypes[ebNE] = self.elementMaterialTypes[eN]
        self.elementBoundaryTypes = numpy.zeros((mesh.nElementBoundaries_global,2),'i')
        self.elementBoundariesArray = mesh.elementBoundariesArray
        for ebN in range(mesh.nElementBoundaries_global):
            eN_left = mesh.elementBoundaryElementsArray[ebN,0]
            eN_right= mesh.elementBoundaryElementsArray[ebN,1]
            self.elementBoundaryTypes[ebN,0] = self.elementMaterialTypes[eN_left]
            if eN_right >= 0:
                self.elementBoundaryTypes[ebN,1] = self.elementMaterialTypes[eN_right]
            else:
                self.elementBoundaryTypes[ebN,1] = self.elementMaterialTypes[eN_left]

    def initializeElementQuadrature(self,t,cq):
        for ci in range(self.nc):
            cq[('f',ci)].flat[:] = 0.0
            for eN in range(cq['x'].shape[0]):
                material=self.elementMaterialTypes[eN]
                for k in range(cq['x'].shape[1]):
                    cq[('a',ci,ci)][eN,k,:] = self.a_types[material](cq['x'][eN,k],t).flat
                    cq[('r',ci)][eN,k]        =-self.source_types[material](cq['x'][eN,k],t)
        #ci
    def initializeElementBoundaryQuadrature(self,t,cebq,cebq_global):
        nd = self.nd
        #use harmonic average for a, arith average for f
        for ci in range(self.nc):
            if ('f',ci) in cebq_global: cebq_global[('f',ci)].flat[:] = 0.0
            if ('f',ci) in cebq: cebq[('f',ci)].flat[:] = 0.0
            for ebN in range(cebq_global['x'].shape[0]):
                material_left = self.elementBoundaryTypes[ebN,0]
                material_right= self.elementBoundaryTypes[ebN,1]
                for k in range(cebq_global['x'].shape[1]):
                    if ('r',ci) in cebq_global:
                        cebq_global[('r',ci)][eN,k] =-0.5*(self.source_types[material_left](cebq_global['x'][ebN,k],t)+
                                                           self.source_types[material_right](cebq_global['x'][ebN,k],t))
                    if ('a',ci,ci) in cebq_global:
                        for i in range(nd):
                            for j in range(nd):
                                x = cebq_global['x'][ebN,k];
                                numer = 2.0*self.a_types[material_left](x,t)[i,j]*self.a_types[material_right](x,t)[i,j]
                                denom = self.a_types[material_left](x,t)[i,j] + self.a_types[material_right](x,t)[i,j] + 1.0e-20
                                cebq_global[('a',ci,ci)][eN,k,i*nd+j] = numer/denom
            for eN in range(cebq['x'].shape[0]):
                for ebN_local in range(cebq['x'].shape[1]):
                    ebN = self.elementBoundariesArray[eN,ebN_local]
                    material_left = self.elementBoundaryTypes[ebN,0]
                    material_right= self.elementBoundaryTypes[ebN,1]
                    for k in range(cebq['x'].shape[2]):
                        x = cebq['x'][eN,ebN_local,k]
                        if ('r',ci) in cebq:
                            cebq[('r',ci)][eN,ebN_local,k] =-0.5*(self.source_types[material_left](x,t)+
                                                                  self.source_types[material_right](x,t))
                        if ('a',ci,ci) in cebq:
                            for i in range(nd):
                                for j in range(nd):
                                    numer = 2.0*self.a_types[material_left](x,t)[i,j]*self.a_types[material_right](x,t)[i,j]
                                    denom = self.a_types[material_left](x,t)[i,j] + self.a_types[material_right](x,t)[i,j] + 1.0e-20
                                    cebq[('a',ci,ci)][eN,ebN_local,k,i*nd+j] = numer/denom
                    #
                #
            #
        #
    def initializeGlobalExteriorElementBoundaryQuadrature(self,t,cebqe):
        nd = self.nd
        for ci in range(self.nc):
            if ('f',ci) in cebqe: cebqe[('f',ci)].flat[:] = 0.0
            for ebNE in range(cebqe['x'].shape[0]):
                material = self.exteriorElementBoundaryTypes[ebNE]
                for k in range(cebqe['x'].shape[1]):
                    x = cebqe['x'][ebNE,k]
                    if ('r',ci) in cebqe:
                        cebqe[('r',ci)][ebNE,k] = -self.source_types[material](x,t)
                    if ('a',ci,ci) in cebqe:
                        cebqe[('a',ci,ci)][ebNE,k,:] = self.a_types[material](x,t).flat
                #
            #
        #
    def evaluate(self,t,c):
        pass #need to put in eval for time varying coefficients
        #mwf debug
        #print "eval c[('a',ci,ci)]= %s" % c[('a',0,0)]
    #end def
########################################

analyticalSolution = None
dirichletConditions = {0:headBCs}

analyticalSolutionVelocity = None

coefficients = SinglePhaseDarcyCoefficients(hydraulicConductivities,sources,nc=1,nd=2)

fluxBoundaryConditions = {0:'setFlow'}

advectiveFluxBoundaryConditions =  {0:noFlowBCs}#{0:dummyBCs}

diffusiveFluxBoundaryConditions = {0:{0:noFlowBCs}}

from pyadh import *
from pyadh.default_p import *

name="gw_2d"

nd = 2

points_on_grain=40
points_on_boundary= 30

# size of 2D rectangular domain
block_size=[1.0, 1.0]

#creat poly file
from block2dDomain import*
domain =block2D(block_size=block_size, points_on_grain=points_on_grain,
                       points_on_boundary=points_on_boundary)
domain.writePoly("block2D")

boundaryTags=domain.boundaryFlags

Ident = numpy.zeros((nd,nd),'d')
Ident[0,0]=1.0; Ident[1,1] = 1.0

#left, right, front, and back heads
head_left = 1.0
head_right= 0.0
head_back = 1.0
head_front = 0.0
initialConditions = None

#define hydraulic conductivity
def hydraulicConductivity(x,t):
    return numpy.array([[10.0, 1.0], [1.0, 100.0]])

def headBCs(x,tag):
    if tag == boundaryTags['left']:
        return lambda x,t: head_left
    if tag == boundaryTags['right']:
        return lambda x,t: head_right
    if tag == boundaryTags['front']:
        return lambda x,t: head_front
    if tag == boundaryTags['back']:
        return lambda x,t: head_back
def noFlowBCs(x,tag):
    if tag == boundaryTags['back']:
        return lambda x,t: 0.0
    if tag == boundaryTags['front']:
        return lambda x,t: 0.0
#################################################
class SinglePhaseDarcyCoefficients(TC_base):
    def __init__(self, a_0,nc=1,nd=2,
                 timeVaryingCoefficients=False):
        self.a_0 = a_0
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
                         hamiltonian)

     
    def initializeMesh(self,mesh):
        self.elementBoundariesArray = mesh.elementBoundariesArray

    def initializeElementQuadrature(self,t,cq):
        for ci in range(self.nc):
            cq[('f',ci)].flat[:] = 0.0
            for eN in range(cq['x'].shape[0]):
                for k in range(cq['x'].shape[1]):
                    cq[('a',ci,ci)][eN,k,:] = self.a_0(cq['x'][eN,k],t).flat
                                
    def initializeElementBoundaryQuadrature(self,t,cebq,cebq_global):
        nd = self.nd
        #use harmonic average for a, arith average for f
        for ci in range(self.nc):
            if cebq_global.has_key(('f',ci)): cebq_global[('f',ci)].flat[:] = 0.0
            if cebq.has_key(('f',ci)): cebq[('f',ci)].flat[:] = 0.0
            for ebN in range(cebq_global['x'].shape[0]):
                for k in range(cebq_global['x'].shape[1]):
                    if cebq_global.has_key(('a',ci,ci)):
                        for i in range(nd):
                            for j in range(nd):
                                x = cebq_global['x'][ebN,k];
                                numer = 2.0*self.a_0(x,t)[i,j]*self.a_0(x,t)[i,j]
                                denom = self.a_0(x,t)[i,j] + self.a_0(x,t)[i,j] + 1.0e-20
                                cebq_global[('a',ci,ci)][eN,k,i*nd+j] = numer/denom
            for eN in range(cebq['x'].shape[0]):
                for ebN_local in range(cebq['x'].shape[1]):
                    ebN = self.elementBoundariesArray[eN,ebN_local]
                    for k in range(cebq['x'].shape[2]):
                        x = cebq['x'][eN,ebN_local,k]
                        if cebq.has_key(('a',ci,ci)):
                            for i in range(nd):
                                for j in range(nd):
                                    numer = 2.0*self.a_0(x,t)[i,j]*self.a_0(x,t)[i,j]
                                    denom = self.a_0(x,t)[i,j] + self.a_0(x,t)[i,j] + 1.0e-20
                                    cebq[('a',ci,ci)][eN,ebN_local,k,i*nd+j] = numer/denom

    def initializeGlobalExteriorElementBoundaryQuadrature(self,t,cebqe):
        nd = self.nd
        for ci in range(self.nc):
            if cebqe.has_key(('f',ci)): cebqe[('f',ci)].flat[:] = 0.0
            for ebNE in range(cebqe['x'].shape[0]):
                for k in range(cebqe['x'].shape[1]):
                    x = cebqe['x'][ebNE,k]
                    if cebqe.has_key(('a',ci,ci)):
                        cebqe[('a',ci,ci)][ebNE,k,:] = self.a_0(x,t).flat
                
            
        
    def evaluate(self,t,c):
        pass 

########################################

analyticalSolution = None
dirichletConditions = {0:headBCs}

analyticalSolutionVelocity = None

coefficients = SinglePhaseDarcyCoefficients(hydraulicConductivity,nc=1,nd=nd)
   

fluxBoundaryConditions = {0:'setFlow'}

advectiveFluxBoundaryConditions =  {0:noFlowBCs}

diffusiveFluxBoundaryConditions = {0:{0:noFlowBCs}}   

        
        

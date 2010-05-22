from pyadh import *
from pyadh.default_p import *
"""
Heterogenous Poisson's equation in 2D using perm field from UTEP
"""

##\page Tests Test Problems 
# \ref poisson_het_2d_p.py "Heterogeneous Poisson's equation"
#

##\ingroup test
#\file poisson_het_2d_p.py
#
#\brief Heterogeneous Poisson's equation in 2D
name="GroundwaterFlow"
L=(1,1,1) #?km
nd = 2

polyfile = "UTEPexample128x128"
permFile = "UTEPpermField128x128.dat"
wellFile = "UTEPwellRates128x128.dat"
#
#polyfile = "UTEPexample16x16"
#permFile = "UTEPpermField16x16.dat"
#wellFile = "UTEPwellRates16x16.dat"

initialConditions = None

Ident = numpy.zeros((nd,nd),'d')
Ident[0,0]=1.0; Ident[1,1] = 1.0

class velEx:
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
##################################################
#u = x + y, a_00 = sin^2(2\pi x) + 5, a_11 = cos^2(2\pi y) + 5.0
#f = 4\pi*(-cos(2\pi*x)sin(2\pi*x) + cos(2\pi*y)sin(2\pi*y))
#dirichlet bc's everywhere

Kx = 5.0e-3# / 2.0**8 #m/s
Ky = Kx#1.0e-6 #m/s
sigma_x = 2.0/L[0]#3.0/L[0]
sigma_y = 2.0/L[1]#6.0/L[1]
def a4(x):
    #return numpy.array([[Kx,0.0],[0.0,Ky]],'d')
    #return numpy.array([[Kx*100.0*(sin(2.0*pi*sigma_x*x[0])**2 + 0.01),0.0],[0.0,Ky*100.0*(cos(2.0*pi*sigma_y*x[1])**2 + 0.01)]],'d')
    return numpy.array([[Kx*(cos(2.0*pi*sigma_x*x[0]) + 1.1)**3 * (cos(2.0*pi*sigma_y*x[1]) + 1.1)**3,0.0],[0.0,Ky*(cos(2.0*pi*sigma_x*x[0]) + 1.1)**3 * (cos(2.0*pi*sigma_y*x[1]) + 1.1)**3]],'d')
def f4(x):
    return 0.0

def getDBC4(x,flag):
    if x[0]==0.0 and x[1] == 0.0:
        return lambda x,t: 10.0 #m
    if x[0]==L[0] and x[1] == L[1]:
        return lambda x,t: 0.0 #m

def getDBCnwCorner(x,flag):
    if x[0] == 0.0 and x[1] >= L[1]-.10*L[1]:
        return lambda x,t: 0.0 #?
#
        
##################################################


class PoissonEquationCoefficients(TransportCoefficients.TC_base):
    """
    Variable coefficient Poisson's equation.
    """
    def __init__(self,aOfX,fOfX,nc=1,permFileName=None,wellFileName=None):
        self.aOfX = aOfX
        self.fOfX = fOfX
        mass = {}
        advection = {}
        diffusion = {}
        potential = {}
        reaction  = {}
        hamiltonian = {}
        for i in range(nc):
            diffusion[i] = {i : {i:'constant'}}
            reaction[i]  = {i : 'constant'}
            potential[i] = {i : 'u'}
        #end i
        sdInfo  = {(0,0):(numpy.array([0,1,2],dtype='i'),
                          numpy.array([0,1],dtype='i'))}
        TransportCoefficients.TC_base.__init__(self,
                                               nc,
                                               mass,
                                               advection,
                                               diffusion,
                                               potential,
                                               reaction,
                                               hamiltonian,
                                               ['phi'],
                                               sparseDiffusionTensors=sdInfo)
        if permFileName != None:
            self.permeabilityTypes = numpy.loadtxt(permFileName,unpack=True)
        else:
            self.permeabilityTypes = numpy.array([1.0],'d')
        if wellFileName != None:
            self.wellRates = numpy.loadtxt(wellFileName,unpack=True)
        else:
            self.wellRates = numpy.array([0.0],'d')
        print "Total Pumping------------------------------------------",sum(self.wellRates)
    def initializeMesh(self,mesh):
        self.elementMaterialTypes = mesh.elementMaterialTypes
        self.mesh = mesh
    def initializeElementQuadrature(self,t,cq):
        for eN in range(cq[('phi',0)].shape[0]):
            cq[('a',0,0)][eN,:,0]= self.permeabilityTypes[self.elementMaterialTypes[eN]]
            cq[('a',0,0)][eN,:,1]= self.permeabilityTypes[self.elementMaterialTypes[eN]]
            cq[('r',0)][eN,:]      = -self.wellRates[self.elementMaterialTypes[eN]]
        #if want to visualize permeability uncomment this command and add 
        #simFlags['plotQuantities'].append("q:('visPerm',0)") to batch file
        cq[('visPerm',0)] = cq[('a',0,0)][:,:,0]
    def initializeElementBoundaryQuadrature(self,t,cebq,cebq_global):
        for eN in range(cebq[('phi',0)].shape[0]):
            #make this harmonic mean?
            cebq[('a',0,0)][eN,:,:,0]= self.permeabilityTypes[self.elementMaterialTypes[eN]]
            cebq[('a',0,0)][eN,:,:,1]= self.permeabilityTypes[self.elementMaterialTypes[eN]]
            cebq[('r',0)][eN,:,:]      = -self.wellRates[self.elementMaterialTypes[eN]]
        #
    def initializeGlobalExteriorElementBoundaryQuadrature(self,t,cebqe):
        for ebNE in range(cebqe[('phi',0)].shape[0]):
            ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
            eN = self.mesh.elementBoundaryElementsArray[ebN,0]
            cebqe[('a',0,0)][ebNE,:,0]= self.permeabilityTypes[self.elementMaterialTypes[eN]]
            cebqe[('a',0,0)][ebNE,:,1]= self.permeabilityTypes[self.elementMaterialTypes[eN]]
    def evaluate(self,t,c):
        pass


nc = 1

#ex 4
aOfX = {0:a4,1:a4}; fOfX = {0:f4,1:f4}
analyticalSolution = None

#dirichletConditions = {0:getDBC4}

dirichletConditions = {0:getDBCnwCorner}

aOfX = {0:a4}
fOfX = {0:f4}

analyticalSolutionVelocity = None

coefficients = PoissonEquationCoefficients(aOfX,fOfX,nc,permFileName=permFile,
                                           wellFileName=wellFile)
   

def zeroInflow(x,flag):
    if (x[0] == 0.0 or
        x[0] == L[0] or
        x[1] == 0.0 or
        x[1] == L[1]):
        return lambda x,t: 0.0

fluxBoundaryConditions = {0:'setFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{0:zeroInflow}}



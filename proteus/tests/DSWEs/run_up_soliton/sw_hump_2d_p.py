from proteus import *
from proteus.default_p import *
from proteus.mprans import DSW2DCV
from proteus.Domain import RectangularDomain
import numpy as np

nd=2

L=(10.0,1.0)
g = 9.81

# for bathymetry
SS = 9.0 #SS for slope starting point 
TH = 0.09  #TH for topography height

# for initial/exact solution 
x0 = 2.0 #where solitary wave starts
h1=  0.1 + .005
h2=  0.14
D = np.sqrt(g * h2)


T = 4.0
nDTout=150

domain = RectangularDomain(L=L,x=[0,0,0])
mannings=0.0

cE=1.0
LUMPED_MASS_MATRIX=1
LINEAR_FRICTION=1

bt = domain.boundaryTags
bt['front'] = bt['bottom']
bt['back'] = bt['top']
domain.writePoly("tank2d")

######################
##### BATHYMETRY #####
######################
def bathymetry_function(X):
    x=X[0] 
    bath = np.piecewise(x, [x < SS, x >= SS], [lambda x: TH, lambda x: TH + 1.0/19.85 * (x-SS)])
    #bath = TH
    return bath

##############################
##### INITIAL CONDITIONS #####
##############################
def solitary(X,t):
    xi = X[0] - D*t-x0
    z1 = 3.0*(h2-h1)
    z2 = h2 * h1**2
    z = np.sqrt(z1 / z2)
    soliton =  h1 + (h2 - h1) * 1.0/(np.cosh(xi/2.0 * z)**2)
    return soliton

class water_height_at_t0:
    def uOfXT(self,X,t):
        #return solitary(X,t)
        return max(solitary(X,t)-bathymetry_function(X),0.0)

class mom_at_t0:
    def uOfXT(self,X,t):
        h =  max(solitary(X,t)-bathymetry_function(X),0.0)
        return D*(h - h1)
       
class heta_at_t0:
    def uOfXT(self,X,t):
        h =  max(solitary(X,t)-bathymetry_function(X),0.0)
        return  h**2
       
class Zero:
    def uOfXT(self,x,t):
        return 0.0

analyticalSolution = {0:water_height_at_t0(),
                      1:Zero(),
                      2:Zero(),
                      3:Zero(),
                      4:Zero()}

initialConditions = {0:water_height_at_t0(),
                     1:mom_at_t0(),
                     2:Zero(),
                     3:heta_at_t0(),
                     4:Zero()}

###################################
##### FOR BOUNDARY CONDITIONS #####
###################################
def getDBC_h(x,flag):
    if x[0]==0: # or x[0]==L[0]:
        return lambda x,t: h1 - TH#bathymetry_function([0])
    elif x[0]==L[0]:
        return lambda x,t: 0.0

def getDBC_hu(x,flag):
    None
   # if x[0]==0: #or x[0]==L[0]:
   #    return lambda x,t: 0.

def getDBC_hv(x,flag):
    return lambda x,t: 0.0

def getDBC_heta(x,flag):
     None
 #    if x[0] == 0: #or x[0]==L[0]:
 #       return lambda x,t: (h1 - TH)**2#bathymetry_function([0]))**2
 #  elif x[0] == L[0]:
 #     return lambda x,t: 0.0

def getDBC_hw(x,flag):
    None
    #return lambda x,t: 0.0

dirichletConditions = {0:getDBC_h,
                       1:getDBC_hu,
                       2:getDBC_hv,
                       3:getDBC_heta,
                       4:getDBC_hw}

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow',
                          3:'outFlow',
                          4:'outFlow'}
advectiveFluxBoundaryConditions =  {}
diffusiveFluxBoundaryConditions = {}

#########################################
##### CREATE MODEL AND COEFFICIENTS #####
#########################################
bathymetry={0:bathymetry_function}
LevelModelType = DSW2DCV.LevelModel
coefficients = DSW2DCV.Coefficients(g=g,bathymetry=bathymetry,cE=cE,LUMPED_MASS_MATRIX=LUMPED_MASS_MATRIX,LINEAR_FRICTION=LINEAR_FRICTION,mannings=mannings)

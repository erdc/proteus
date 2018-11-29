from proteus import *
from proteus.default_p import *
from proteus.mprans import DSW2DCV
from proteus.Domain import RectangularDomain
import numpy as np

nd=2

# set up stufff #
L=(45.0,0.8)
X_coords = (-35.0,10.0)
domain = RectangularDomain(L=L,x=[-35.0,0,0])
T=20.1
nDTout=200
cE=1001.0
LUMPED_MASS_MATRIX=1
LINEAR_FRICTION=0

mannings = 0.014

bt = domain.boundaryTags
bt['front'] = bt['bottom']
bt['back'] = bt['top']
domain.writePoly("tank2d")

###############################
#  CONSTANTS NEEDED FOR SETUP #
###############################
g = 9.81
h0 = 1.0
a = 0.30  # amplitude
slope = 1.0 / 19.850
k_wavenumber = np.sqrt(3.0 * a/(4.0 * h0**3))  # wavenumber
z = np.sqrt(3.0 * a * h0) / (2.0 * h0 * np.sqrt(h0 * (1.0 + a)))
L_wave = 2.0 / k_wavenumber * np.arccosh(np.sqrt(1.0 / 0.050))  # wavelength of solitary wave
c = np.sqrt(g * (1.0 + a) * h0)
x0 = - h0/slope - L_wave/2.0  # location of the toe of the beach

def solitary_wave(x,t):
    sechSqd = (1.00/np.cosh( z*(x-x0-c*t)))**2.00
    soliton = a * h0 * sechSqd
    return soliton

######################
##### BATHYMETRY #####
######################
def bathymetry_function(X):
    x = X[0]
    return numpy.maximum(slope*x,-h0)

##############################
##### INITIAL CONDITIONS #####
##############################
class water_height_at_t0(object):
    def uOfXT(self,X,t):
        eta = solitary_wave(X[0],0)
        h = eta-bathymetry_function(X)
        hp = max(h,0.)
        return hp

class x_mom_at_t0(object):
    def uOfXT(self,X,t):
        eta = solitary_wave(X[0],0)
        h = eta-bathymetry_function(X)
        hp = max(h,0.)
        Umom = hp * c * eta / (h0 + eta)
        return Umom

class y_mom_at_t0(object):
    def uOfXT(self,X,t):
        h = water_height_at_t0().uOfXT(X,t)
        return 0.

class heta_at_t0(object):
    def uOfXT(self,X,t):
        h = water_height_at_t0().uOfXT(X,t)
        return h**2

class hw_at_t0(object):
    def uOfXT(self,X,t):
        eta = solitary_wave(X[0],0)
        h = eta-bathymetry_function(X)
        hp = max(h,0.)
        hprime = -2.0 * z * eta * np.tanh(z*(X[0]-x0-c*t))
        hw = hp * (-c * h0 * eta * hprime /(h0 + eta)**2)
        return hw

initialConditions = {0:water_height_at_t0(),
                     1:x_mom_at_t0(),
                     2:y_mom_at_t0(),
                     3:heta_at_t0(),
                     4:hw_at_t0()}

###################################
##### FOR BOUNDARY CONDITIONS #####
###################################
def water_height_DBC(X,flag):
    if X[0]==X_coords[0]:
        return lambda x,t: water_height_at_t0().uOfXT(X,0.0)
    elif X[0]==X_coords[1]:
        return lambda x,t: water_height_at_t0().uOfXT(X,0.0)

def x_mom_DBC(X,flag):
    if X[0]==X_coords[0]:
        return lambda X,t: x_mom_at_t0().uOfXT(X,0.0)

def y_mom_DBC(X,flag):
    return lambda x,t: 0.0

def heta_DBC(X,flag):
    if X[0]==X_coords[0]:
        return lambda x,t: heta_at_t0().uOfXT(X,0.0)

def hw_DBC(X,flag):
    if X[0]==X_coords[0]:
        return lambda x,t: hw_at_t0().uOfXT(X,0.0)

dirichletConditions = {0:water_height_DBC,
                       1:x_mom_DBC,
                       2:y_mom_DBC,
                       3:heta_DBC,
                       4:hw_DBC}


#########################################
##### CREATE MODEL AND COEFFICIENTS #####
#########################################
bathymetry={0:bathymetry_function}
LevelModelType = DSW2DCV.LevelModel
coefficients = DSW2DCV.Coefficients(g=g,bathymetry=bathymetry,cE=cE,LUMPED_MASS_MATRIX=LUMPED_MASS_MATRIX,LINEAR_FRICTION=LINEAR_FRICTION,mannings=mannings)

from pyadh import *
from pyadh.default_p import *

import santillanaDamBreakDomain

nd=2
genMesh=True
polyfile='sw_dam_break_santillana'
up_len= 2.0; mid_len =10.0; down_len=2.0
up_width=3.0; mid_width=1.0
up_height=2.0
def bathymetry(x):
    if x[0] < up_len:
        return up_height
    if x[0] > up_len+mid_len:
        return 0.0
    return up_height-up_height*(x[0]-up_len)/mid_len
def bathymetryGrad(x):
    if x[0] < up_len:
        return (0.0,0.0)
    if x[0] > up_len+mid_len:
        return (0.0,0.0)
    return (-up_height/mid_len,0.0)
    
domain = santillanaDamBreakDomain.twoReservoirDomain(up_width=up_width,
                                                   mid_width=mid_width,
                                                   up_len=up_len,
                                                   mid_len=mid_len,
                                                   down_len=down_len,
                                                   name=polyfile)

domain.writeAsymptote(polyfile)
domain.writePoly(polyfile)

g = 9.8
H = 1.0
U = 1.0
T=2.0
dt_init = 1.0e-4
coefficients = ShallowWater(g=g,
                            nd=nd,
                            bathymetryFunc=bathymetry,
                            bathymetryGradientFunc=bathymetryGrad)

class ContractionIC:
    def __init__(self,H0=H,HB=0.0):
        self.H0=H0
        self.HB=HB
    def uOfXT(self,x,t):
        if x[0] > up_len:
            return self.HB
        return self.H0

class ZeroIC:
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:ContractionIC(),
                     1:ZeroIC(),
                     2:ZeroIC()}

def getDBC_h(x,flag):
    return None

def getDBC_hu(x,flag):
    if flag in [domain.boundaryTags['upstream'],domain.boundaryTags['downstream'],domain.boundaryTags['wall_vert']]:
        return lambda x,t:0.0

def getDBC_hv(x,flag):
    if flag not in  [domain.boundaryTags['upstream'],domain.boundaryTags['downstream'],domain.boundaryTags['wall_vert']]:
        return lambda x,t: 0.0

dirichletConditions = {0:getDBC_h,
                       1:getDBC_hu,
                       2:getDBC_hv}

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow'}

def getAFBC_h(x,flag):
    if flag in domain.boundaryTags.values():
        return lambda x,t: 0.0
def getAFBC_hu(x,flag):
    pass
    
def getAFBC_hv(x,flag):
    pass

advectiveFluxBoundaryConditions =  {0:getAFBC_h,
                                    1:getAFBC_hu,
                                    2:getAFBC_hv}

diffusiveFluxBoundaryConditions = {0:{},1:{},2:{}}




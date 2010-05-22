from pyadh import *
from pyadh.default_p import *
import contraction2d
nd=2

genMesh=True
polyfile='sw_contraction2d'
boundaryTags = contraction2d.genPoly(polyfile,
                                     curve_fun=contraction2d.linear_profile)

g = 9.8
H = 1.0
U = 1.0
T=1.0#10.0
dt_init = 1.0e-4
coefficients = ShallowWater(g=g,
                            nd=nd)

class ContractionIC:
    def __init__(self,H0=H):
        self.H0=H0
    def uOfXT(self,x,t):
        return self.H0

class ZeroIC:
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:ContractionIC(),
                     1:ZeroIC(),
                     2:ZeroIC()}

def getDBC_h(x,flag):
#     if flag ==  boundaryTags['upstream']:
#         return lambda x,t: H
    if flag ==  boundaryTags['downstream']:
        return lambda x,t: H
#    if flag ==  boundaryTags['downstream']:
#        return lambda x,t: H

def getDBC_hu(x,flag):
#     if flag == boundaryTags['upstream']:
#         return lambda x,t: H*U
    if flag == boundaryTags['downstream']:
        return lambda x,t: H*U

def getDBC_hv(x,flag):
#     if flag == boundaryTags['upstream']:
#         return lambda x,t: 0.0
    if flag == boundaryTags['downstream']:
        return lambda x,t: 0.0

dirichletConditions = {0:getDBC_h,
                       1:getDBC_hu,
                       2:getDBC_hv}

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow'}

def getAFBC_h(x,flag):
    if flag == boundaryTags['upstream']:
        return lambda x,t: -1.0
    if not (flag == boundaryTags['upstream'] or 
            flag == boundaryTags['downstream']):
        return lambda x,t: 0.0

def getAFBC_hu(x,flag):
    pass
#    if flag == boundaryTags['upstream']:
#        return lambda x,t: -1.0
#     if not (flag == boundaryTags['upstream'] or 
#             flag == boundaryTags['downstream']):
#         return lambda x,t: 0.0
    
def getAFBC_hv(x,flag):
    pass
#    if flag == boundaryTags['upstream']:
#        return lambda x,t: 0.0
#     if not (flag == boundaryTags['upstream'] or 
#             flag == boundaryTags['downstream']):
#         return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC_h,
                                    1:getAFBC_hu,
                                    2:getAFBC_hv}

diffusiveFluxBoundaryConditions = {0:{},1:{},2:{}}




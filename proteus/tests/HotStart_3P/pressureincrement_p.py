from __future__ import absolute_import
from builtins import object
from math import *
from proteus import *
from proteus.default_p import *
from .NS_hotstart import *

#domain = ctx.domain
#nd = ctx.nd
name = "pressureincrement"

from proteus.mprans import PresInc
coefficients=PresInc.Coefficients(rho_f_min = (1.0-1.0e-8)*rho_1,
                                  rho_s_min = (1.0-1.0e-8)*rho_s,
                                  nd = nd,
                                  modelIndex=1,
                                  fluidModelIndex=0,
                                  fixNullSpace=fixNullSpace_PresInc, 
                                  INTEGRATE_BY_PARTS_DIV_U=INTEGRATE_BY_PARTS_DIV_U_PresInc)
LevelModelType = PresInc.LevelModel

def getDBC_phi(x,flag):
   None
def getAdvectiveFlux_qt(x,flag):
   if manufactured_solution==1: #u.n!=0                                                                                                                                                                                                                                            
       if (flag==1): #left boundary                                                                                                                                                                                                                                                
           return lambda x,t: -np.sin(x[0])*np.sin(x[1]+t)
       elif (flag==2): # right boundary                                                                                                                                                                                                                                            
           return lambda x,t: np.sin(x[0])*np.sin(x[1]+t)
       elif (flag==3): # bottom boundary                                                                                                                                                                                                                                          
           return lambda x,t: -np.cos(x[0])*np.cos(x[1]+t)
       else:
           return lambda x,t: np.cos(x[0])*np.cos(x[1]+t)
   else: #u.n=0                                                                                                                                                                                                                                                                    
       return lambda x,t: 0.
def getDiffusiveFlux_phi(x,flag):
   return lambda x,t: 0.
class getIBC_phi(object):
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:getIBC_phi()}
dirichletConditions = {0:getDBC_phi}
advectiveFluxBoundaryConditions = {0:getAdvectiveFlux_qt}
diffusiveFluxBoundaryConditions = {0:{0:getDiffusiveFlux_phi}}

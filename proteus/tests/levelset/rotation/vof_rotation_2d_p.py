from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from rotation2D import *
from proteus.mprans import VOF
name=soname+"_vof"

"""
The non-conservative level set description of a bubble in a two-phase flow
"""

LevelModelType = VOF.LevelModel

##\ingroup test
#\file vof_rotation_2d_p.py
#
# \todo finish vof_rotation_2d_p.py

if applyRedistancing:
    coefficients = VOF.Coefficients(LS_model=0,V_model=0,RD_model=1,ME_model=2,checkMass=checkMass,
                                   epsFact=epsFact_vof,useMetrics=useMetrics)
elif not onlyVOF:
    coefficients = VOF.Coefficients(LS_model=0,V_model=0,RD_model=None,ME_model=1,checkMass=checkMass,
                                    epsFact=epsFact_vof,useMetrics=useMetrics)
else:
    coefficients = VOF.Coefficients(RD_model=None,ME_model=0,checkMass=checkMass,
                                    epsFact=epsFact_vof,useMetrics=useMetrics)


def Heaviside(phi):
    if phi > 0:
        return 1.0
    elif phi < 0:
        return 0.0
    else:
        return 0.5

class Rotation_phi:
    def __init__(self,L=[1.0,
                         1.0],
                 center=[0.5,
                         0.5],
                 radius=0.45):#45/100 = 9/20
        self.radius = radius
        self.slotwidth = radius*3.0/9.0
        self.slotlength = radius
        self.xc = center
        self.xnw = [center[0] - 0.5*self.slotwidth,center[1] - (radius - self.slotlength)]
        self.xne = [center[0] + 0.5*self.slotwidth,center[1] - (radius - self.slotlength)]
        self.xsw = [center[0] - 0.5*self.slotwidth,center[1] - (radius)]
        self.xse = [center[0] + 0.5*self.slotwidth,center[1] - (radius)]
    def uOfXT(self,x,t):
        from math import sqrt
        dist = lambda u,v: sqrt( (u[0] - v[0])**2 + (u[1] - v[1])**2)
        phic = dist(self.xc,x) - self.radius
        phine = -dist(self.xne,x)
        phinw = -dist(self.xnw,x)
        phise = dist(self.xse,x)
        phisw = dist(self.xsw,x)
        phin = self.xnw[1] - x[1]
        phis = -(self.xsw[1] - x[1])
        phie = self.xne[0] - x[0]
        phiw = -(self.xnw[0] - x[0])
        if x[1] >= self.xnw[1]:
            if x[0] < self.xnw[0]:
                phi = max(phic,phinw)
            else:
                if x[0] < self.xne[0]:
                    phi = max(phic,phin)
                else:
                    phi = max(phic,phine)
        elif x[1] >= self.xsw[1]:
            if x[0] < self.xnw[0]:
                phi = max(phic,phiw)
            else:
                if x[0] < self.xne[0]:
                    phi = min([phin,phie,phiw])
                else:
                    phi = max(phic,phie)
        else:
            if x[0] < self.xsw[0]:
                phi = phic
            else:
                if x[0] < self.xse[0]:
                    phi = min(phisw,phise)
                else:
                    phi = phic
        return smoothedHeaviside(epsFactHeaviside*he,-phi)

class Rotation_phi_cylinder:
    def __init__(self,center=[0.5,0.75,0.5],radius=0.15):
        self.radius  = radius
        self.center  = center
    def uOfX(self,X):
        zalesak=True
        if not zalesak:
            dx = X[0]-self.center[0]; dy = X[1]-self.center[1];
            dBubble = self.radius - sqrt(dx**2 + dy**2)
            return smoothedHeaviside(epsFactHeaviside*he,dBubble)#Heaviside(dBubble)
        else:
            X = x[0]
            Y = x[1]
            # distance to center of the disk
            r = math.sqrt((X-self.xc)**2 + (Y-self.yc)**2)
            dist = self.radius - r # distance to circle
            # coordinates to slot coorners
            xslot1 = self.xc-0.025
            xslot2 = self.xc+0.025
            yslot1 = 0.75 - np.sqrt(self.radius**2-0.025**2)
            yslot2 = 0.85

            #distance to the boundary of the slot
            aux1 = np.abs(X-xslot1)
            aux2 = np.abs(X-xslot2)
            aux3 = np.abs(Y-yslot2)
            aux4 = np.abs(Y-yslot1)

            if (Y > yslot1): #above disk
                if (xslot1 < X and X <= xslot2 and y <= yslot2): #inside slot
                    return smoothedHeaviside(epsFactHeaviside*he,-np.min([aux1,aux2,aux3]))
                else: #Not inside slot
                    if X <= xslot1: #left of slot
                        if Y <= yslot2: #down top boundary of slot
                            return smoothedHeaviside(epsFactHeaviside*he,np.min([dist,aux1]))
                        else: #above top boundary of slot
                            return smoothedHeaviside(epsFactHeaviside*he,np.min([dist,np.sqrt(aux1**2+aux3**2)]))
                    elif X >= xslot2: #right of slot
                        if Y <= yslot2: #down top boundary of slot
                            return smoothedHeaviside(epsFactHeaviside*he,np.min([dist,aux2]))
                        else: #above top boundary of slot
                            return smoothedHeaviside(epsFactHeaviside*he,np.min([dist,np.sqrt(aux2**2+aux3**2)]))
                    else: #x-coordiante is within slot
                        return smoothedHeaviside(epsFactHeaviside*he,np.min([dist,aux3]))
            else: #below disk
                if X > xslot1 and X < xslot2: #x-coordinate is within slot
                    return  smoothedHeaviside(epsFactHeaviside*he,-np.min([np.sqrt(aux1**2 + aux4**2), np.sqrt(aux2**2 + aux4**2)]))
                else:
                    return smoothedHeaviside(epsFactHeaviside*he,dist)
            
    #end
    def uOfXT(self,X,t):
        return self.uOfX(X)
    #end
#end Rotation_phi

analyticalSolutions = {0:Rotation_phi(L=L,
                                      center=[0.0,
                                              0.5],
                                      radius=0.25)}

def getDBC(x,flag):
    pass

dirichletConditions = {0:getDBC}

initialConditions  = {0:analyticalSolutions[0]}

fluxBoundaryConditions = {0:'outFlow'}

#cek made no flux since v.n = 0 for this v
def getAFBC(x,flag):
   return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC}

diffusiveFluxBoundaryConditions = {0:{}}

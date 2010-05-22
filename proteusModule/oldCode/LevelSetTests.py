#!/usr/bin/env python
from math import *
from EGeometry import *
from ScalarTransport import *
"""
Classes for defining level set test problems 
"""

## \defgroup LevelSetTests LevelSetTests
#
# Classes for defining level set test problems                                                           # 
# @{   

class UnitSquareRotation(ScalarTransportCoefficients):
    from transportCoefficients import unitSquareRotationEvaluate
    def __init__(self):
        self.useC=True
    def evaluate(self,
                 t,
                 x,
                 u,
                 m,dm,
                 f,df,
                 a,da,
                 phi,dphi,
                 r,dr):
        if self.useC:
            self.unitSquareRotationEvaluate(u.shape[0],
                                            f.shape[1],
                                            t,
                                            x,
                                            u,
                                            m,dm,
                                            f,df,
                                            a,da,
                                            phi,dphi,
                                            r,dr)
        else:
            m[:] = u[:]
            dm[:] = 1.0
            #mwf debug
            #print 'in UnitSquareRot u=\n',u,' m= ',m

            for i in range(u.shape[0]):
                vx = 2*pi*(x[i][1] - 0.5)
                vy = 2*pi*(0.5     - x[i][0]) 
                #mwf debug
                vx = 1.0
                vy = 1.0
                f[i][0] = vx*u[i]
                f[i][1] = vy*u[i]
                df[i][0] = vx
                df[i][1] = vy
            #end for
    #end evaluate
#end unit square rotation
class UnitSquareRotationWithDiffusion(ScalarTransportCoefficients):
    def __init__(self,A0=0.0):
        self.A0 = A0
    def evaluate(self,
                 t,
                 x,
                 u,
                 m,dm,
                 f,df,
                 a,da,
                 phi,dphi,
                 r,dr):
        m[:]  = u
        dm[:]  = 1.0
        phi[:] = u
        dphi[:]=1.0
        a[:]   =self.A0
        da[:]  =0.0
        for i in range(u.shape[0]):
            vx = 2*pi*(x[i][1] - 0.5)
            vy = 2*pi*(0.5     - x[i][0]) 
            f[i][0] = vx*u[i]
            f[i][1] = vy*u[i]
            df[i][0] = vx
            df[i][1] = vy
        #end for
    #end evaluate
#end unit square rotation
class RotatingCone2D:
    def __init__(self,radius):
        self.radius = radius
    #end def
    def uOfXT(self,x,t):
        centerX = 0.25*sin(2*pi*t) + 0.5
        centerY = 0.25*cos(2*pi*t) + 0.5
        coneX = x[0] - centerX
        coneY = x[1] - centerY
        if sqrt(coneX**2 + coneY**2) < self.radius:
            return 0.25*(1.0 + cos(pi*coneX/self.radius))*(1.0+cos(pi*coneY/self.radius))
        else:
            return 0.0
        #end else
    #end def uOfXT
#end class def
 
def getDBC_hom(x):
    if x == 0.0:
        return lambda x,t: 0.0
    if x == 1.0:
        return lambda x,t: 0.0
#end def

def getHomogeneousDBC2D(x):
    if x[X] == 0.0:
        return lambda x,t: 0.0
    elif x[X] == 1.0:
        return lambda x,t: 0.0
    elif x[Y] == 0.0:
        return lambda x,t: 0.0
    elif x[Y] == 1.0:
        return lambda x,t: 0.0
    #end else
#end def

#cek
def getNoDBC(x):
    pass

## @}

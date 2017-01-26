from proteus.TransportCoefficients import *

class LAD(TC_base):
    """
    The coefficients of the linear advection-diffusion equation
    """
    def __init__(self,M,A,B):
        TC_base.__init__(self, 
                         nc=1, #number of components
                         variableNames=['u'],
                         mass      = {0:{0:'linear'}},
                         advection = {0:{0:'linear'}},
                         diffusion = {0:{0:{0:'constant'}}},
                         potential = {0:{0:'u'}},
                         reaction  = {0:{0:'linear'}})
        self.M=M;
        self.A=A;
        self.B=B;
    
    def evaluate(self,t,c):
        c[('m',0)][:]         = self.M*c[('u',0)]  
        c[('dm',0,0)][:]      = self.M
        c[('f',0)][...,0]     = self.B[0]*c[('u',0)]
        c[('f',0)][...,1]     = self.B[1]*c[('u',0)]
        c[('df',0,0)][...,0]  = self.B[0]
        c[('df',0,0)][...,1]  = self.B[1]
        c[('a',0,0)][...,0,0] = self.A[0][0]
        c[('a',0,0)][...,1,1] = self.A[1][1]

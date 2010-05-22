from pyadh.TransportCoefficients import *
import numpy
class LAD(TC_base):
    def __init__(self,M,A,B):
        self.M=M;self.A=A;self.B=B;
        #diagonal
        sdInfo = {(0,0):(numpy.arange(start=0,stop=3,step=1,dtype='i'),
                         numpy.arange(start=0,stop=2,step=1,dtype='i'))}
        #full
        #sdInfo = {(0,0):(numpy.array([0,2,4]),
        #                 numpy.array([0,1,0,1]))}
        
        TC_base.__init__(self, nc=1,variableNames=['u'],
                         mass={0:{0:'linear'}},
                         advection={0:{0:'linear'}},
                         diffusion={0:{0:{0:'constant'}}},
                         potential={0:{0:'u'}},
                         reaction={0:{0:'linear'}},
                         hamiltonian={},
                         sparseDiffusionTensors=sdInfo)
    def evaluate(self,t,c):
        c[('m',0)][:]=self.M*c[('u',0)];  
        c[('dm',0,0)][:]=self.M
        for i,u in enumerate(c[('u',0)].flat):
             c[('f',0)].flat[i*2+0]=self.B[0]*u
             c[('f',0)].flat[i*2+1]=self.B[1]*u
             c[('df',0,0)].flat[i*2+0]=self.B[0]
             c[('df',0,0)].flat[i*2+1]=self.B[1]
             c[('a',0,0)].flat[i*2+0]=self.A[0][0]
             c[('a',0,0)].flat[i*2+1]=self.A[1][1]
             #c[('a',0,0)].flat[i*4+0]=self.A[0][0]
             #c[('a',0,0)].flat[i*4+3]=self.A[1][1]


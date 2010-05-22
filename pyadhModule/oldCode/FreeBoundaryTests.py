from math import *
import cPickle
import copy
import os
from MeshTools import *
from FemTools import *
from QuadTools import *
from LinearSolvers import *
from NonlinearSolvers import *
import Gnuplot
from NormTools import *
from AnalyticalSolutions import *
import ScalarTransport		
import HamiltonJacobi
from LatexReport import *
import transportCoefficients
"""
A module for testing solvers for free boundary problems.
"""
## \defgroup FreeBoundaryTests FreeBoundaryTests
#
# A module for testing solvers for free boundary problems. 
#
# @{

class TwophasePotentialFlow(ScalarTransport.ScalarTransportCoefficients):
    from transportCoefficients import twophasePotentialFlowUpdateFreeSurface
    from transportCoefficients import twophasePotentialFlowEvaluate
    def __init__(self,
                 M1,M2,
                 A1,A2,
                 B1,B2,
                 Bcon1,Bcon2,
                 C1,C2):
        ScalarTransport.ScalarTransportCoefficients.__init__(self,
                                                             mass='linear',
                                                             advection='linear',
                                                             diffusion='constant',
                                                             potential='linear',
                                                             reaction='linear')
        self.M1=M1
        self.M2=M2
        self.A1=A1
        self.A2=A2
        self.B1=B1
        self.B2=B2
        self.Bcon1=Bcon1
        self.Bcon2=Bcon2
        self.C1=C1
        self.C2=C2
        self.useC=True
    def initializeFreeSurface(self,uq_levelSetList,epsilonIn):
        self.shapeDict={}
        self.MqList=[]
        self.BqList=[]
        self.BconqList=[]
        self.AqList=[]
        self.CqList=[]
        for s,uq_levelSet in enumerate(uq_levelSetList):
            self.shapeDict[uq_levelSet.shape] = s
            self.MqList.append(Numeric.zeros(uq_levelSet.shape,Numeric.Float))
            self.BqList.append(Numeric.zeros((uq_levelSet.shape[0],self.B1.shape[0]),
                                  Numeric.Float))
            self.BconqList.append(Numeric.zeros((uq_levelSet.shape[0],self.B1.shape[0]),
                                  Numeric.Float))
            self.AqList.append(Numeric.zeros((uq_levelSet.shape[0],self.A1.shape[0],self.A1.shape[1]),
                                  Numeric.Float))
            self.CqList.append(Numeric.zeros(uq_levelSet.shape,Numeric.Float))
        self.epsilon=epsilonIn
        self.updateFreeSurface(uq_levelSetList)
    def updateFreeSurface(self,uq_levelSetList):
        for uq_levelSet in uq_levelSetList:
            s = self.shapeDict[uq_levelSet.shape]
            if self.useC:
                self.twophasePotentialFlowUpdateFreeSurface(uq_levelSet.shape[0],
                                                            self.B1.shape[0],
                                                            self.epsilon,
                                                            uq_levelSet,
                                                            self.M1,self.M2,self.MqList[s],
                                                            self.A1,self.A2,self.AqList[s],
                                                            self.B1,self.B2,self.BqList[s],
                                                            self.Bcon1,self.Bcon2,self.BconqList[s],
                                                            self.C1,self.C2,self.CqList[s])
            else:
                for i in range(uq_levelSet.shape[0]):
                    if uq_levelSet[i] > self.epsilon:
                        self.MqList[s][i]=self.M1
                        self.AqList[s][i,:,:]=self.A1
                        self.BqList[s][i,:]=self.B1
                        self.BconqList[s][i,:]=self.Bcon1
                        self.CqList[s][i]=self.C1
                    elif uq_levelSet[i] < -self.epsilon:
                        self.MqList[s][i]=self.M2
                        self.AqList[s][i,:,:]=self.A2
                        self.BqList[s][i,:]=self.B2
                        self.BconqList[s][i,:]=self.Bcon2
                        self.CqList[s][i]=self.C2
                    else:
                        H = 0.5*(1.0 + uq_levelSet[i]/self.epsilon + math.sin((math.pi*uq_levelSet[i])/self.epsilon)/math.pi)
                        oneMinusH=1.0-H
                        self.MqList[s][i]=oneMinusH*self.M2 + H*self.M1
                        self.AqList[s][i,:,:]=oneMinusH*self.A2 + H*self.A1
                        self.BqList[s][i,:]=oneMinusH*self.B2 + H*self.B1
                        self.BconqList[s][i,:]=oneMinusH*self.Bcon2 + H*self.Bcon1
                        self.CqList[s][i]=oneMinusH*self.C2+H*self.C1

    def evaluate(self,
                 t,
                 x,
                 u,
                 m,dm,
                 f,df,
                 a,da,
                 phi,dphi,
                 r,dr):
        s = self.shapeDict[u.shape]
        if self.useC:
            self.twophasePotentialFlowEvaluate(x.shape[0],
                                               f.shape[1],
                                               self.MqList[s],
                                               self.AqList[s],
                                               self.BqList[s],
                                               self.BconqList[s],
                                               self.CqList[s],
                                               t,
                                               x,
                                               u,
                                               m,dm,
                                               f,df,
                                               a,da,
                                               phi,dphi,
                                               r,dr)
        else:
            for i in range(u.shape[0]):
                m[i] = self.MqList[s][i]*u[i]
                dm[i] = self.MqList[s][i]
                f[i,:]=self.BqList[s][i,:]*u[i] + self.BconqList[s][i,:]
                df[i,:]=self.BqList[s][i,:]
                a[i,:,:]=self.AqList[s][i,:,:]
                da[i,:,:]=0.0
                phi[i]=u[i]
                dphi[i]=1.0
                r[i]=self.CqList[s][i]*u[i]
                dr[i]=self.CqList[s][i]

class TwophaseLevelSetCoefficients(HamiltonJacobi.HamiltonJacobiCoefficients):
    from transportCoefficients import twophaseLevelSetCoefficientsUpdateVelocity
    from transportCoefficients import twophaseLevelSetCoefficientsEvaluate
    def __init__(self,velocityScaling):
        self.vq_scale = velocityScaling
        self.useC=True
    def initializeVelocity(self,vqList):
        self.shapeDict={}
        self.BqList=[]
        for s,vq in enumerate(vqList):
            self.shapeDict[vq.shape]=s
            self.BqList.append(Numeric.zeros(vq.shape,
                                             Numeric.Float))
        self.updateVelocity(vqList)
    def updateVelocity(self,vqList):
        for vq in vqList:
            s = self.shapeDict[vq.shape]
            if self.useC:
                self.twophaseLevelSetCoefficientsUpdateVelocity(vq.shape[0],
                                                                vq.shape[1],
                                                                self.vq_scale,
                                                                vq,
                                                                self.BqList[s])
            else:
                for i in range(vq.shape[0]):
                    for  I in range(vq.shape[1]):
                        self.BqList[s][i,I]=vq[i,I]*self.vq_scale
    def evaluate(self,
                 t,
                 x,
                 u,grad_u,
                 m,dm,
                 h,dh,
                 rh):
        s = self.shapeDict[dh.shape]
        if self.useC:
            self.twophaseLevelSetCoefficientsEvaluate(u.shape[0],
                                                      dh.shape[1],
                                                      self.BqList[s],
                                                      t,
                                                      x,
                                                      u,grad_u,
                                                      m,dm,
                                                      h,dh,
                                                      rh)
        else:
            rh[:]=0.0
            h[:]=0.0
            for i in range(u.shape[0]):
                m[i]=u[i]
                dm[i]=1.0
                for  I in range(grad_u.shape[1]):
                    h[i]+=self.BqList[s][i,I]*grad_u[i,I]
                    dh[i,I]=self.BqList[s][i,I]

class ciTwoPhaseLevelSetCoefficients(ScalarTransport.ScalarTransportCoefficients):
    from transportCoefficients import twophaseLevelSetCoefficientsUpdateVelocity
    from transportCoefficients import twophaseLevelSetCoefficientsEvaluateCI
    def __init__(self,velocityScaling):
        ScalarTransport.ScalarTransportCoefficients.__init__(self,
                                                             mass='linear',
                                                             advection='linear',
                                                             diffusion=None,
                                                             potential=None,
                                                             reaction=None)
        self.vq_scale = velocityScaling
        self.useC=True
    def initializeVelocity(self,vqList):
        self.shapeDict={}
        self.BqList=[]
        for s,vq in enumerate(vqList):
            self.shapeDict[vq.shape]=s
            self.BqList.append(Numeric.zeros(vq.shape,
                                             Numeric.Float))
        self.updateVelocity(vqList)
    def updateVelocity(self,vqList):
        for vq in vqList:
            s = self.shapeDict[vq.shape]
            if self.useC:
                self.twophaseLevelSetCoefficientsUpdateVelocity(vq.shape[0],
                                                                vq.shape[1],
                                                                self.vq_scale,
                                                                vq,
                                                                self.BqList[s])
            else:
                for i in range(vq.shape[0]):
                    for  I in range(vq.shape[1]):
                        self.BqList[s][i,I]=vq[i,I]*self.vq_scale
    def evaluate(self,
                 t,
                 x,
                 u,
                 m,dm,
                 f,df,
                 a,da,
                 phi,dphi,
                 r,dr):
        s = self.shapeDict[f.shape]
        if self.useC:
            self.twophaseLevelSetCoefficientsEvaluateCI(u.shape[0],
                                                        f.shape[1],
                                                        self.BqList[s].flat,
                                                        t,
                                                        x,
                                                        u,
                                                        m,dm,
                                                        f,df,
                                                        a,da,
                                                        phi,dphi,
                                                        r,dr)
        else:
            for i in range(u.shape[0]):
                m[i]=u[i]
                dm[i]=1.0
                for  I in range(f.shape[1]):
                    f[i,I]=self.BqList[s][i,I]*u[i]
                    df[i,I]=self.BqList[s][i,I]

class TwophaseSignedDistanceCoefficients(HamiltonJacobi.HamiltonJacobiCoefficients):
    def __init__(self):
        self.mass='linear'
        self.hamiltonian='nonlinear'
    def initializeSignFunction(self,uqList,grad_uqList,epsilonIn):
        self.SqList=[]
        self.shapeDict={}
	self.epsilon=epsilonIn
        for s,uq in enumerate(uqList):
            self.shapeDict[uq.shape]=s
            self.SqList.append(Numeric.zeros(uq.shape,
                                             Numeric.Float))
    def updateSignFunction(self,uqList,grad_uqList):
        for uq,grad_uq in zip(uqList,grad_uqList):
            s = self.shapeDict[uq.shape]
            transportCoefficients.twophaseSignedDistanceUpdateSignFunction(uq.shape[0],
                                                                           self.epsilon,
                                                                           uq,
                                                                           self.SqList[s])
#              s = self.shapeDict[uq.shape]
#              for i in range(uq.shape[0]):
#                  self.SqList[s][i] =  uq[i]/(math.sqrt(uq[i]*uq[i] + diameter*diameter))
#                 norm_grad_u_2=0.0
#                 for I in range(grad_uq.shape[1]):
#                     norm_grad_u_2+=grad_uq[i,I]*grad_uq[i,I]
#                 self.SqList[s][i] =  uq[i]/(math.sqrt(uq[i]*uq[i] + norm_grad_u_2*diameter*diameter))
#  		if  uq[i] > 0.0:
#  		    self.SqList[s][i]=1.0
#  		elif uq[i] < 0.0:
#  		    self.SqList[s][i]=-1.0
#  		else:
#  		    self.SqList[s][i]=0.0
    def evaluate(self,
                 t,
                 x,
                 u,grad_u,
                 m,dm,
                 h,dh,
                 rh):
        s = self.shapeDict[u.shape]
        transportCoefficients.twophaseSignedDistanceCoefficientsEvaluate(u.shape[0],
                                                                         dh.shape[1],
                                                                         self.SqList[s],
                                                                         u,
                                                                         grad_u,
                                                                         m,dm,
                                                                         h,dh,
                                                                         rh)
#          for i in range(grad_u.shape[0]):
#              print grad_u.shape[1]
#              print "|| grad(u) || - 1", (l2Norm(grad_u[i])-1.0)
#              print "h",h[i]
 #         rh[:]=-1.0
#          dm[:]=1.0
#  	h[:]=0.0
#          for i in range(u.shape[0]):
#              m[i]=u[i]
#              for I in range(grad_u.shape[1]):
#                  h[i] += grad_u[i,I]*grad_u[i,I]
#              h[i]=math.sqrt(h[i])
#              for  I in range(grad_u.shape[1]):
#                  if  h[i] != 0.0:
#                      dh[i,I]=(self.SqList[s][i]*grad_u[i,I])/h[i]
#                  else:
#                      dh[i,I]=(self.SqList[s][i]*grad_u[i,I])/1.0e-8
#              h[i]+=rh[i]
#              h[i]*=self.SqList[s][i]
#              rh[i]*=self.SqList[s][i]

def runTests():
    def redistanceLoop():
        for l in range(nLevels):
            signedDistance.modelList[l].u.dof[:] = freeSurface.modelList[l].u.dof
            signedDistance.uList[l][:] = freeSurface.modelList[l].u.dof
            signedDistance.modelList[l].coefficients.updateSignFunction([freeSurface.modelList[l].q['u'].flat,
                                                                         freeSurface.modelList[l].ebq['u'].flat],
                                                                        [Numeric.reshape(freeSurface.modelList[l].q['grad(u)'].flat,
                                                                                         (freeSurface.modelList[l].nQuadraturePoints_global,
                                                                                          freeSurface.modelList[l].nSpace_global)),
                                                                         Numeric.reshape(freeSurface.modelList[l].ebq['grad(u)'].flat,
                                                                                         (freeSurface.modelList[l].nElementBoundaryQuadraturePoints_global,
                                                                                          freeSurface.modelList[l].nSpace_global))])
            signedDistance.modelList[l].updateZeroLevelSetDirichletConditions()
        signedDistance.initializeTimeIntegration()
        signedDistance.chooseDT()
	signedDistance.updateTimeHistory()
        its = 0
	for l in range(nLevels):
		signedDistance.modelList[l].trialSpace.writeFunctionEnsight(signedDistance.modelList[l].u,'signedDistance'+`l`,append=True)
	redistanceTimeValues.append(signedDistance.modelList[-1].T)
        signedDistance.modelList[-1].getResidual(u = signedDistance.uList[-1],
                                                 r = signedDistance.rList[-1])
        rnorm0 = l2Norm(signedDistance.rList[-1])
        rnorm = rnorm0
        dt=None
        print "rnorm0 %12.5e" % rnorm0
        nonlinearSolverFailed=False
        while ( (rnorm > 10.0*atol) 
                and its < 50
                and nonlinearSolverFailed==False):
            its +=1
            signedDistance.chooseDT(DTSET=dt)
            print "DT %12.5e" % signedDistance.DT
            nonlinearSolverFailed = signedDistanceNonlinearSolver.solveMultilevel(uList = signedDistance.uList,
                                                                                  rList = signedDistance.rList)
            signedDistance.updateTimeHistory()
	    redistanceTimeValues.append(signedDistance.modelList[-1].T)
	    for l in range(nLevels):
		signedDistance.modelList[l].trialSpace.writeFunctionEnsight(signedDistance.modelList[l].u,'signedDistance'+`l`,append=True)
            signedDistance.modelList[-1].getResidual(u = signedDistance.uList[-1],
                                                     r = signedDistance.rList[-1])
            rnormLast = rnorm
            rnorm = l2Norm(signedDistance.rList[-1])
            print "rnorm %12.5e" % rnorm
            dt = signedDistance.DT*min(1.5,rnormLast/rnorm)
 #             for l in range(nLevels):
#                  signedDistance.modelList[l].coefficients.updateSignFunction([signedDistance.modelList[l].q['u'].flat,
#                                                                               signedDistance.modelList[l].ebq['u'].flat],
#                                                                              [Numeric.reshape(signedDistance.modelList[l].q['grad(u)'].flat,
#                                                                                               (signedDistance.modelList[l].nQuadraturePoints_global,
#                                                                                                signedDistance.modelList[l].nSpace_global)),
#                                                                               Numeric.reshape(signedDistance.modelList[l].ebq['grad(u)'].flat,
#                                                                                               (signedDistance.modelList[l].nElementBoundaryQuadraturePoints_global,
#                                                                                                signedDistance.modelList[l].nSpace_global))])
        #print signedDistanceNonlinearSolver.info()
        #if nonlinearSolverFailed == False:
        for l in range(nLevels):
            freeSurface.modelList[l].u.dof[:] = signedDistance.modelList[l].u.dof
            freeSurface.modelList[l].setFreeDOF(freeSurface.uList[l])
            freeSurface.modelList[l].updateCoefficients()
        #else:
            #print "Redistancing Failed"
    #First define the mesh independent problem definition
    #consisting of the following information
    #see examples below for definitions
    nd={}  #number of spatial dimensions
    getFlowDirichletConditions={}
    getFreeSurfaceDirichletConditions={}
    getSignedDistanceDirichletConditions={}
    flowCoefficients={}
    nciFreeSurfaceCoefficients={}
    ciFreeSurfaceCoefficients={}
    flowAnalyticalSolution={}
    freeSurfaceAnalyticalSolution={}
    flowTimeIntegration={}
    nciFreeSurfaceTimeIntegration={}
    ciFreeSurfaceTimeIntegration={}
    signedDistanceTimeIntegration={}
    getFlowInitialConditions={}
    getFreeSurfaceInitialConditions={}
    runCFL = {}
    T = {}
    flowFullNewtonFlag = {}
    freeSurfaceFullNewtonFlag = {}
    signedDistanceFullNewtonFlag = {}
    readMesh = {}
    flowNBC = {}
    flowNBCobject = {}
    #
    #porous media properties
    #
    lengthScale  = 1.0    #m
    permeability = 1.0e-7 #m^2
    viscosity_1  = 8.9e-4 #kg/(m*s)
    viscosity_2  = 7.1e-6 #kg/(m*s)
    #viscosity_2=viscosity_1
    density_1    = 997.0  #kg/m^3
    density_2    = 1.205  #kg/m^3
    #density_2=density_1
    gravity      = 9.8    #m/s^2
    porosity     = 0.5    #-
    #make non-dimensional
    K_1 = density_1*gravity*permeability/(viscosity_1*sqrt(gravity*lengthScale))
    K_2 = density_1*gravity*permeability/(viscosity_2*sqrt(gravity*lengthScale))
    dimensionless_density_1 = 1.0
    dimensionless_density_2 = density_2/density_1
    dimensionless_gravity  = 1.0
    #Define test problems
    testProblems = []
    #
    # pressure  boundary conditions on entire boundary consistant with finalWaterDepth
    # initial conditions  consistent  with initialWaterDepth "up" is given by vector dg
    #
    test = 'RisingWaterTable2D-BE-North'
    testProblems.append(test)
    flowTimeIntegration[test]=ScalarTransport.NoIntegration
    nciFreeSurfaceTimeIntegration[test]=HamiltonJacobi.BackwardEuler
    ciFreeSurfaceTimeIntegration[test]=ScalarTransport.BackwardEuler
    signedDistanceTimeIntegration[test]=HamiltonJacobi.BackwardEuler
    readMesh[test] = False
    runCFL[test]= 2.5
    nd[test]=2
    finalWaterDepth = 0.5
    initialWaterDepth = 0.25
    theta = math.pi/2.0
    north = Numeric.array([math.cos(theta),math.sin(theta)])
    def getFreeSurfaceDBC_SteadyState2D(x):
        pass
    class FlowSolution_SteadyState2D(SteadyState):
        def __init__(self,waterDepthIn,dg=north):
            self.waterDepth=waterDepthIn
            self.dg = dg
            self.top = Numeric.dot(self.dg,Numeric.array([1.0,1.0]))
        def uOfX(self,x):
            h = Numeric.dot(self.dg,x)
            if h >= self.waterDepth:
                return (self.top-h)*dimensionless_density_2*dimensionless_gravity
            else:
                return ((self.top-self.waterDepth)*dimensionless_density_2 +(self.waterDepth-h)*dimensionless_density_1)*dimensionless_gravity
    flowAnalyticalSolution[test]=risingWaterTableSol_north=FlowSolution_SteadyState2D(finalWaterDepth,north)
    class FreeSurfaceSolution_SteadyState2D(SteadyState):
        def __init__(self,waterDepthIn,dg=north):
            self.waterDepth=waterDepthIn
            self.dg=dg
        def uOfX(self,x):
            return self.waterDepth - Numeric.dot(self.dg,x)
    freeSurfaceAnalyticalSolution[test]=FreeSurfaceSolution_SteadyState2D(finalWaterDepth,north)
    def getRisingWaterTable2DDBC_north(xIn):
        if (xIn[0] == 0.0 or
            xIn[0] == 1.0 or
            xIn[1] == 0.0 or
            xIn[1] == 1.0):
            return lambda x,t: risingWaterTableSol_north.uOfX(x)
#          if (xIn[1] == 0.0 or
#              xIn[1] == 1.0):
#              return lambda x,t: risingWaterTableSol.uOfX(x)
    getFlowDirichletConditions[test]=getRisingWaterTable2DDBC_north
    getFreeSurfaceDirichletConditions[test]=getFreeSurfaceDBC_SteadyState2D
    flowCoefficients[test]=TwophasePotentialFlow(M1=0.0,M2=0.0,
                                                 B1=Numeric.array([0.0,0.0]),
                                                 B2=Numeric.array([0.0,0.0]),
                                                 Bcon1=-K_1*dimensionless_density_1*dimensionless_gravity*north,
                                                 Bcon2=-K_2*dimensionless_density_2*dimensionless_gravity*north,
                                                 A1=Numeric.array([[K_1,0.0],[0.0,K_1]]),
                                                 A2=Numeric.array([[K_2,0.0],[0.0,K_2]]),
                                                 C1=0.0,C2=0.0)
    flowCoefficients[test].mass=None
    flowCoefficients[test].advection='constant'
    flowCoefficients[test].diffusion='constant'
    flowCoefficients[test].potential='linear'
    flowCoefficients[test].reaction=None
    nciFreeSurfaceCoefficients[test]=TwophaseLevelSetCoefficients(1.0/porosity)
    nciFreeSurfaceCoefficients[test].mass='linear'
    nciFreeSurfaceCoefficients[test].hamiltonian='linear'
    ciFreeSurfaceCoefficients[test]=ciTwoPhaseLevelSetCoefficients(1.0/porosity)
    ciFreeSurfaceCoefficients[test].mass='linear'
    ciFreeSurfaceCoefficients[test].advection='linear'
    T[test]=0.5#2.0
    getFlowInitialConditions[test] = FlowSolution_SteadyState2D(initialWaterDepth,north)
    getFreeSurfaceInitialConditions[test] = FreeSurfaceSolution_SteadyState2D(initialWaterDepth,north)
    flowFullNewtonFlag[test] = False
    freeSurfaceFullNewtonFlag[test]=False
    signedDistanceFullNewtonFlag[test]=True
    flowNBC[test] = 'noFlow'
    def getSignedDistanceDBC(x):
        pass
    getSignedDistanceDirichletConditions[test]=getSignedDistanceDBC
    #
    # same as above with different dg
    #
    test = 'RisingWaterTable2D-BE-Northeast'
    testProblems.append(test)
    flowTimeIntegration[test]=ScalarTransport.NoIntegration
    nciFreeSurfaceTimeIntegration[test]=HamiltonJacobi.BackwardEuler
    ciFreeSurfaceTimeIntegration[test]=ScalarTransport.BackwardEuler
    signedDistanceTimeIntegration[test]=HamiltonJacobi.BackwardEuler
    readMesh[test] = False
    runCFL[test]= 2.5
    nd[test]=2
    finalWaterDepth = 0.5
    initialWaterDepth = 0.25
    theta = math.pi/4.0
    northeast = Numeric.array([math.cos(theta),math.sin(theta)])
    flowAnalyticalSolution[test]=risingWaterTableSol_northeast=FlowSolution_SteadyState2D(finalWaterDepth,northeast)
    freeSurfaceAnalyticalSolution[test]=FreeSurfaceSolution_SteadyState2D(finalWaterDepth,northeast)
    def getRisingWaterTable2DDBC_northeast(xIn):
        if (xIn[0] == 0.0 or
            xIn[0] == 1.0 or
            xIn[1] == 0.0 or
            xIn[1] == 1.0):
            return lambda x,t: risingWaterTableSol_northeast.uOfX(x)
    getFlowDirichletConditions[test]=getRisingWaterTable2DDBC_northeast
    getFreeSurfaceDirichletConditions[test]=getFreeSurfaceDBC_SteadyState2D
    flowCoefficients[test]=TwophasePotentialFlow(M1=0.0,M2=0.0,
                                                 B1=Numeric.array([0.0,0.0]),
                                                 B2=Numeric.array([0.0,0.0]),
                                                 Bcon1=-K_1*dimensionless_density_1*dimensionless_gravity*northeast,
                                                 Bcon2=-K_2*dimensionless_density_2*dimensionless_gravity*northeast,
                                                 A1=Numeric.array([[K_1,0.0],[0.0,K_1]]),
                                                 A2=Numeric.array([[K_2,0.0],[0.0,K_2]]),
                                                 C1=0.0,C2=0.0)
    flowCoefficients[test].mass=None
    flowCoefficients[test].advection='constant'
    flowCoefficients[test].diffusion='constant'
    flowCoefficients[test].potential='linear'
    flowCoefficients[test].reaction=None
    nciFreeSurfaceCoefficients[test]=TwophaseLevelSetCoefficients(1.0/porosity)
    nciFreeSurfaceCoefficients[test].mass='linear'
    nciFreeSurfaceCoefficients[test].hamiltonian='linear'
    ciFreeSurfaceCoefficients[test]=ciTwoPhaseLevelSetCoefficients(1.0/porosity)
    ciFreeSurfaceCoefficients[test].mass='linear'
    ciFreeSurfaceCoefficients[test].advection='linear'
    T[test]=0.5#2.0
    getFlowInitialConditions[test] = FlowSolution_SteadyState2D(initialWaterDepth,northeast)
    getFreeSurfaceInitialConditions[test] = FreeSurfaceSolution_SteadyState2D(initialWaterDepth,northeast)
    flowFullNewtonFlag[test] = False
    freeSurfaceFullNewtonFlag[test]=False
    signedDistanceFullNewtonFlag[test]=True
    flowNBC[test] = 'noFlow'
    getSignedDistanceDirichletConditions[test]=getSignedDistanceDBC
    #
    # same as above with different dg
    #
    test = 'RisingWaterTable2D-BE-East'
    testProblems.append(test)
    flowTimeIntegration[test]=ScalarTransport.NoIntegration
    nciFreeSurfaceTimeIntegration[test]=HamiltonJacobi.BackwardEuler
    ciFreeSurfaceTimeIntegration[test]=ScalarTransport.BackwardEuler
    signedDistanceTimeIntegration[test]=HamiltonJacobi.BackwardEuler
    readMesh[test] = False
    runCFL[test]= 2.5
    nd[test]=2
    finalWaterDepth = 0.5
    initialWaterDepth = 0.25
    theta = 0.0
    east = Numeric.array([math.cos(theta),math.sin(theta)])
    flowAnalyticalSolution[test]=risingWaterTableSol_east=FlowSolution_SteadyState2D(finalWaterDepth,east)
    freeSurfaceAnalyticalSolution[test]=FreeSurfaceSolution_SteadyState2D(finalWaterDepth,east)
    def getRisingWaterTable2DDBC_east(xIn):
        if (xIn[0] == 0.0 or
            xIn[0] == 1.0 or
            xIn[1] == 0.0 or
            xIn[1] == 1.0):
            return lambda x,t: risingWaterTableSol_east.uOfX(x)
    getFlowDirichletConditions[test]=getRisingWaterTable2DDBC_east
    getFreeSurfaceDirichletConditions[test]=getFreeSurfaceDBC_SteadyState2D
    flowCoefficients[test]=TwophasePotentialFlow(M1=0.0,M2=0.0,
                                                 B1=Numeric.array([0.0,0.0]),
                                                 B2=Numeric.array([0.0,0.0]),
                                                 Bcon1=-K_1*dimensionless_density_1*dimensionless_gravity*east,
                                                 Bcon2=-K_2*dimensionless_density_2*dimensionless_gravity*east,
                                                 A1=Numeric.array([[K_1,0.0],[0.0,K_1]]),
                                                 A2=Numeric.array([[K_2,0.0],[0.0,K_2]]),
                                                 C1=0.0,C2=0.0)
    flowCoefficients[test].mass=None
    flowCoefficients[test].advection='constant'
    flowCoefficients[test].diffusion='constant'
    flowCoefficients[test].potential='linear'
    flowCoefficients[test].reaction=None
    nciFreeSurfaceCoefficients[test]=TwophaseLevelSetCoefficients(1.0/porosity)
    nciFreeSurfaceCoefficients[test].mass='linear'
    nciFreeSurfaceCoefficients[test].hamiltonian='linear'
    ciFreeSurfaceCoefficients[test]=ciTwoPhaseLevelSetCoefficients(1.0/porosity)
    ciFreeSurfaceCoefficients[test].mass='linear'
    ciFreeSurfaceCoefficients[test].advection='linear'
    T[test]=0.5#2.0
    getFlowInitialConditions[test] = FlowSolution_SteadyState2D(initialWaterDepth,east)
    getFreeSurfaceInitialConditions[test] = FreeSurfaceSolution_SteadyState2D(initialWaterDepth,east)
    flowFullNewtonFlag[test] = False
    freeSurfaceFullNewtonFlag[test]=False
    signedDistanceFullNewtonFlag[test]=True
    flowNBC[test] = 'noFlow'
    getSignedDistanceDirichletConditions[test]=getSignedDistanceDBC
    #
    # same as previous test with intial and final depths swapped
    #
    test = 'RecedingWaterTable2D-BE'
    testProblems.append(test)
    flowTimeIntegration[test]=ScalarTransport.NoIntegration
    nciFreeSurfaceTimeIntegration[test]=HamiltonJacobi.BackwardEuler
    ciFreeSurfaceTimeIntegration[test]=ScalarTransport.BackwardEuler
    signedDistanceTimeIntegration[test]=HamiltonJacobi.BackwardEuler
    readMesh[test] = False
    runCFL[test]= 2.5
    nd[test]=2
    finalWaterDepth = 0.50
    initialWaterDepth = 0.75
    flowAnalyticalSolution[test]=recedingWaterTableSol=FlowSolution_SteadyState2D(finalWaterDepth)
    freeSurfaceAnalyticalSolution[test]=FreeSurfaceSolution_SteadyState2D(finalWaterDepth)
    def getRecedingWaterTable2DDBC(xIn):
        if (xIn[0] == 0.0 or
            xIn[0] == 1.0 or
            xIn[1] == 0.0 or
            xIn[1] == 1.0):
            return lambda x,t: recedingWaterTableSol.uOfX(x)
#          if (xIn[1] == 0.0 or
#              xIn[1] == 1.0):
#              return lambda x,t: recedingWaterTableSol.uOfX(x)
    getFlowDirichletConditions[test]=getRecedingWaterTable2DDBC
    getFreeSurfaceDirichletConditions[test]=getFreeSurfaceDBC_SteadyState2D
    flowCoefficients[test]=TwophasePotentialFlow(M1=0.0,M2=0.0,
                                                 B1=Numeric.array([0.0,0.0]),
                                                 B2=Numeric.array([0.0,0.0]),
                                                 Bcon1=Numeric.array([0.0,-K_1*dimensionless_density_1*dimensionless_gravity]),
                                                 Bcon2=Numeric.array([0.0,-K_2*dimensionless_density_2*dimensionless_gravity]),
                                                 A1=Numeric.array([[K_1,0.0],[0.0,K_1]]),
                                                 A2=Numeric.array([[K_2,0.0],[0.0,K_2]]),
                                                 C1=0.0,C2=0.0)
    flowCoefficients[test].mass=None
    flowCoefficients[test].advection='constant'
    flowCoefficients[test].diffusion='constant'
    flowCoefficients[test].potential='linear'
    flowCoefficients[test].reaction=None
    nciFreeSurfaceCoefficients[test]=TwophaseLevelSetCoefficients(1.0/porosity)
    nciFreeSurfaceCoefficients[test].mass='linear'
    nciFreeSurfaceCoefficients[test].hamiltonian='linear'
    ciFreeSurfaceCoefficients[test]=ciTwoPhaseLevelSetCoefficients(1.0/porosity)
    ciFreeSurfaceCoefficients[test].mass='linear'
    ciFreeSurfaceCoefficients[test].advection='linear'
    T[test]=0.5#2.0
    getFlowInitialConditions[test] = FlowSolution_SteadyState2D(initialWaterDepth)
    getFreeSurfaceInitialConditions[test] = FreeSurfaceSolution_SteadyState2D(initialWaterDepth)
    flowFullNewtonFlag[test] = False
    freeSurfaceFullNewtonFlag[test]=False
    signedDistanceFullNewtonFlag[test]=True
    flowNBC[test] = 'noFlow'
    getSignedDistanceDirichletConditions[test]=getSignedDistanceDBC
    #
    # No flow conditions on sides and bottom
    #
    test = 'PerturbedBucket2D-BE'
    testProblems.append(test)
    flowTimeIntegration[test]=ScalarTransport.NoIntegration
    nciFreeSurfaceTimeIntegration[test]=HamiltonJacobi.BackwardEuler
    ciFreeSurfaceTimeIntegration[test]=ScalarTransport.BackwardEuler
    signedDistanceTimeIntegration[test]=HamiltonJacobi.BackwardEuler
    readMesh[test] = False
    runCFL[test]= 2.5
    waterDepth=0.5
    nd[test]=2
    def getBucketBC2D(xIn):
        if xIn[1] == 1.0:
            return lambda x,t: 0.0
    flowAnalyticalSolution[test]=FlowSolution_SteadyState2D(waterDepth)
    freeSurfaceAnalyticalSolution[test]=FreeSurfaceSolution_SteadyState2D(waterDepth)
    getFlowDirichletConditions[test]=getBucketBC2D
    getFreeSurfaceDirichletConditions[test]=getFreeSurfaceDBC_SteadyState2D
    flowCoefficients[test]=TwophasePotentialFlow(M1=0.0,M2=0.0,
                                                 B1=Numeric.array([0.0,0.0]),
                                                 B2=Numeric.array([0.0,0.0]),
                                                 Bcon1=Numeric.array([0.0,-K_1*dimensionless_density_1*dimensionless_gravity]),
                                                 Bcon2=Numeric.array([0.0,-K_2*dimensionless_density_2*dimensionless_gravity]),
                                                 A1=Numeric.array([[K_1,0.0],[0.0,K_1]]),
                                                 A2=Numeric.array([[K_2,0.0],[0.0,K_2]]),
                                                 C1=0.0,C2=0.0)
    flowCoefficients[test].mass=None
    flowCoefficients[test].advection='constant'
    flowCoefficients[test].diffusion='constant'
    flowCoefficients[test].potential='linear'
    flowCoefficients[test].reaction=None
    nciFreeSurfaceCoefficients[test]=TwophaseLevelSetCoefficients(1.0/porosity)
    nciFreeSurfaceCoefficients[test].mass='linear'
    nciFreeSurfaceCoefficients[test].hamiltonian='linear'
    ciFreeSurfaceCoefficients[test]=ciTwoPhaseLevelSetCoefficients(1.0/porosity)
    ciFreeSurfaceCoefficients[test].mass='linear'
    ciFreeSurfaceCoefficients[test].advection='linear'
    T[test]=0.5#2.0
    getFlowInitialConditions[test] = FlowSolution_SteadyState2D(initialWaterDepth)
    class sinusoidalFreeSurface2D(SteadyState):
        def __init__(self,waterDepthIn,surfaceWaveNumber=1,surfaceAmplitude=0.1):
            self.waterDepth=waterDepthIn
            self.k=surfaceWaveNumber
            self.amp=surfaceAmplitude
        def uOfX(self,x):
            return self.waterDepth + self.amp*math.sin(self.k*math.pi*x[0])- x[1]
    getFreeSurfaceInitialConditions[test] = sinusoidalFreeSurface2D(waterDepth,surfaceWaveNumber=1,surfaceAmplitude=0.1)
    flowFullNewtonFlag[test] = False
    freeSurfaceFullNewtonFlag[test]=False
    signedDistanceFullNewtonFlag[test]=True
    flowNBC[test] = 'noFlow'
    getSignedDistanceDirichletConditions[test]=getSignedDistanceDBC
    #
    # No flow conditions on sides and bottom
    #
    test = 'InvertedPerturbedBucket2D-BE'
    testProblems.append(test)
    flowTimeIntegration[test]=ScalarTransport.NoIntegration
    nciFreeSurfaceTimeIntegration[test]=HamiltonJacobi.BackwardEuler
    ciFreeSurfaceTimeIntegration[test]=ScalarTransport.BackwardEuler
    signedDistanceTimeIntegration[test]=HamiltonJacobi.BackwardEuler
    readMesh[test] = False
    runCFL[test]= 2.5
    waterDepth=0.5
    nd[test]=2
    def getBucketBC2D(xIn):
        if xIn[1] == 1.0:
            return lambda x,t: 0.0
    flowAnalyticalSolution[test]=FlowSolution_SteadyState2D(waterDepth)
    freeSurfaceAnalyticalSolution[test]=FreeSurfaceSolution_SteadyState2D(waterDepth)
    getFlowDirichletConditions[test]=getBucketBC2D
    getFreeSurfaceDirichletConditions[test]=getFreeSurfaceDBC_SteadyState2D
    flowCoefficients[test]=TwophasePotentialFlow(M1=0.0,M2=0.0,
                                                 B1=Numeric.array([0.0,0.0]),
                                                 B2=Numeric.array([0.0,0.0]),
                                                 Bcon1=Numeric.array([0.0,-K_1*dimensionless_density_1*dimensionless_gravity]),
                                                 Bcon2=Numeric.array([0.0,-K_2*dimensionless_density_2*dimensionless_gravity]),
                                                 A1=Numeric.array([[K_1,0.0],[0.0,K_1]]),
                                                 A2=Numeric.array([[K_2,0.0],[0.0,K_2]]),
                                                 C1=0.0,C2=0.0)
    flowCoefficients[test].mass=None
    flowCoefficients[test].advection='constant'
    flowCoefficients[test].diffusion='constant'
    flowCoefficients[test].potential='linear'
    flowCoefficients[test].reaction=None
    nciFreeSurfaceCoefficients[test]=TwophaseLevelSetCoefficients(1.0/porosity)
    nciFreeSurfaceCoefficients[test].mass='linear'
    nciFreeSurfaceCoefficients[test].hamiltonian='linear'
    ciFreeSurfaceCoefficients[test]=ciTwoPhaseLevelSetCoefficients(1.0/porosity)
    ciFreeSurfaceCoefficients[test].mass='linear'
    ciFreeSurfaceCoefficients[test].advection='linear'
    T[test]=0.5#2.0
    getFlowInitialConditions[test] = FlowSolution_SteadyState2D(initialWaterDepth)
    class invertedSinusoidalFreeSurface2D(SteadyState):
        def __init__(self,waterDepthIn,surfaceWaveNumber=1,surfaceAmplitude=0.1):
            self.waterDepth=waterDepthIn
            self.k=surfaceWaveNumber
            self.amp=surfaceAmplitude
        def uOfX(self,x):
            return min(-(self.waterDepth + self.amp*math.sin(self.k*math.pi*x[0] + math.pi)- x[1]),self.waterDepth + 2.0*self.amp - x[1])
    getFreeSurfaceInitialConditions[test] = invertedSinusoidalFreeSurface2D(waterDepth,surfaceWaveNumber=1,surfaceAmplitude=0.1)
    flowFullNewtonFlag[test] = False
    freeSurfaceFullNewtonFlag[test]=False
    signedDistanceFullNewtonFlag[test]=True
    flowNBC[test] = 'noFlow'
    getSignedDistanceDirichletConditions[test]=getSignedDistanceDBC
    #
    # pressure  boundary conditions on entire boundary consistant with finalWaterDepth
    # initial conditions  consistent  with initialWaterDepth
    #
    test = 'RisingWaterTable2D-FE'
    testProblems.append(test)
    flowTimeIntegration[test]=ScalarTransport.NoIntegration
    nciFreeSurfaceTimeIntegration[test]=HamiltonJacobi.ForwardEuler
    ciFreeSurfaceTimeIntegration[test]=ScalarTransport.ForwardEuler
    signedDistanceTimeIntegration[test]=HamiltonJacobi.BackwardEuler
    readMesh[test] = False
    runCFL[test]= 0.1
    nd[test]=2
    finalWaterDepth = 0.5
    initialWaterDepth = 0.25
    def getFreeSurfaceDBC_SteadyState2D(x):
        pass
    freeSurfaceAnalyticalSolution[test]=FreeSurfaceSolution_SteadyState2D(finalWaterDepth)
    getFlowDirichletConditions[test]=getRisingWaterTable2DDBC_north
    getFreeSurfaceDirichletConditions[test]=getFreeSurfaceDBC_SteadyState2D
    flowCoefficients[test]=TwophasePotentialFlow(M1=0.0,M2=0.0,
                                                 B1=Numeric.array([0.0,0.0]),
                                                 B2=Numeric.array([0.0,0.0]),
                                                 Bcon1=Numeric.array([0.0,-K_1*dimensionless_density_1*dimensionless_gravity]),
                                                 Bcon2=Numeric.array([0.0,-K_2*dimensionless_density_2*dimensionless_gravity]),
                                                 A1=Numeric.array([[K_1,0.0],[0.0,K_1]]),
                                                 A2=Numeric.array([[K_2,0.0],[0.0,K_2]]),
                                                 C1=0.0,C2=0.0)
    flowCoefficients[test].mass=None
    flowCoefficients[test].advection='constant'
    flowCoefficients[test].diffusion='constant'
    flowCoefficients[test].potential='linear'
    flowCoefficients[test].reaction=None
    nciFreeSurfaceCoefficients[test]=TwophaseLevelSetCoefficients(1.0/porosity)
    nciFreeSurfaceCoefficients[test].mass='linear'
    nciFreeSurfaceCoefficients[test].hamiltonian='linear'
    ciFreeSurfaceCoefficients[test]=ciTwoPhaseLevelSetCoefficients(1.0/porosity)
    ciFreeSurfaceCoefficients[test].mass='linear'
    ciFreeSurfaceCoefficients[test].advection='linear'
    T[test]=0.5#2.0
    getFlowInitialConditions[test] = FlowSolution_SteadyState2D(initialWaterDepth)
    getFreeSurfaceInitialConditions[test] = FreeSurfaceSolution_SteadyState2D(initialWaterDepth)
    flowFullNewtonFlag[test] = False
    freeSurfaceFullNewtonFlag[test]=False
    signedDistanceFullNewtonFlag[test]=True
    flowNBC[test] = 'noFlow'
    getSignedDistanceDirichletConditions[test]=getSignedDistanceDBC
    #
    # same as previous test with intial and final depths swapped
    #
    test = 'RecedingWaterTable2D-FE'
    testProblems.append(test)
    flowTimeIntegration[test]=ScalarTransport.NoIntegration
    nciFreeSurfaceTimeIntegration[test]=HamiltonJacobi.ForwardEuler
    ciFreeSurfaceTimeIntegration[test]=ScalarTransport.ForwardEuler
    signedDistanceTimeIntegration[test]=HamiltonJacobi.BackwardEuler
    readMesh[test] = False
    runCFL[test] = 0.1
    nd[test]=2
    finalWaterDepth = 0.50
    initialWaterDepth = 0.75
    flowAnalyticalSolution[test]=recedingWaterTableSol=FlowSolution_SteadyState2D(finalWaterDepth)
    freeSurfaceAnalyticalSolution[test]=FreeSurfaceSolution_SteadyState2D(finalWaterDepth)
    getFlowDirichletConditions[test]=getRecedingWaterTable2DDBC
    getFreeSurfaceDirichletConditions[test]=getFreeSurfaceDBC_SteadyState2D
    flowCoefficients[test]=TwophasePotentialFlow(M1=0.0,M2=0.0,
                                                 B1=Numeric.array([0.0,0.0]),
                                                 B2=Numeric.array([0.0,0.0]),
                                                 Bcon1=Numeric.array([0.0,-K_1*dimensionless_density_1*dimensionless_gravity]),
                                                 Bcon2=Numeric.array([0.0,-K_2*dimensionless_density_2*dimensionless_gravity]),
                                                 A1=Numeric.array([[K_1,0.0],[0.0,K_1]]),
                                                 A2=Numeric.array([[K_2,0.0],[0.0,K_2]]),
                                                 C1=0.0,C2=0.0)
    flowCoefficients[test].mass=None
    flowCoefficients[test].advection='constant'
    flowCoefficients[test].diffusion='constant'
    flowCoefficients[test].potential='linear'
    flowCoefficients[test].reaction=None
    nciFreeSurfaceCoefficients[test]=TwophaseLevelSetCoefficients(1.0/porosity)
    nciFreeSurfaceCoefficients[test].mass='linear'
    nciFreeSurfaceCoefficients[test].hamiltonian='linear'
    ciFreeSurfaceCoefficients[test]=ciTwoPhaseLevelSetCoefficients(1.0/porosity)
    ciFreeSurfaceCoefficients[test].mass='linear'
    ciFreeSurfaceCoefficients[test].advection='linear'
    T[test]=2.0
    getFlowInitialConditions[test] = FlowSolution_SteadyState2D(initialWaterDepth)
    getFreeSurfaceInitialConditions[test] = FreeSurfaceSolution_SteadyState2D(initialWaterDepth)
    flowFullNewtonFlag[test] = False
    freeSurfaceFullNewtonFlag[test]=False
    signedDistanceFullNewtonFlag[test]=True
    flowNBC[test] = 'noFlow'
    getSignedDistanceDirichletConditions[test]=getSignedDistanceDBC
    #
    # No flow conditions on sides and bottom
    #
    test = 'PerturbedBucket2D-FE'
    testProblems.append(test)
    flowTimeIntegration[test]=ScalarTransport.NoIntegration
    nciFreeSurfaceTimeIntegration[test]=HamiltonJacobi.ForwardEuler
    ciFreeSurfaceTimeIntegration[test]=ScalarTransport.ForwardEuler
    signedDistanceTimeIntegration[test]=HamiltonJacobi.BackwardEuler
    readMesh[test] = False
    runCFL[test] = 0.1
    waterDepth=0.5
    nd[test]=2
    flowAnalyticalSolution[test]=FlowSolution_SteadyState2D(waterDepth)
    freeSurfaceAnalyticalSolution[test]=FreeSurfaceSolution_SteadyState2D(waterDepth)
    getFlowDirichletConditions[test]=getBucketBC2D
    getFreeSurfaceDirichletConditions[test]=getFreeSurfaceDBC_SteadyState2D
    flowCoefficients[test]=TwophasePotentialFlow(M1=0.0,M2=0.0,
                                                 B1=Numeric.array([0.0,0.0]),
                                                 B2=Numeric.array([0.0,0.0]),
                                                 Bcon1=Numeric.array([0.0,-K_1*dimensionless_density_1*dimensionless_gravity]),
                                                 Bcon2=Numeric.array([0.0,-K_2*dimensionless_density_2*dimensionless_gravity]),
                                                 A1=Numeric.array([[K_1,0.0],[0.0,K_1]]),
                                                 A2=Numeric.array([[K_2,0.0],[0.0,K_2]]),
                                                 C1=0.0,C2=0.0)
    flowCoefficients[test].mass=None
    flowCoefficients[test].advection='constant'
    flowCoefficients[test].diffusion='constant'
    flowCoefficients[test].potential='linear'
    flowCoefficients[test].reaction=None
    nciFreeSurfaceCoefficients[test]=TwophaseLevelSetCoefficients(1.0/porosity)
    nciFreeSurfaceCoefficients[test].mass='linear'
    nciFreeSurfaceCoefficients[test].hamiltonian='linear'
    ciFreeSurfaceCoefficients[test]=ciTwoPhaseLevelSetCoefficients(1.0/porosity)
    ciFreeSurfaceCoefficients[test].mass='linear'
    ciFreeSurfaceCoefficients[test].advection='linear'
    T[test]=2.0
    getFlowInitialConditions[test] = FlowSolution_SteadyState2D(initialWaterDepth)
    getFreeSurfaceInitialConditions[test] = sinusoidalFreeSurface2D(waterDepth,surfaceWaveNumber=1,surfaceAmplitude=0.1)
    flowFullNewtonFlag[test] = False
    freeSurfaceFullNewtonFlag[test]=False
    signedDistanceFullNewtonFlag[test]=True
    flowNBC[test] = 'noFlow'
    getSignedDistanceDirichletConditions[test]=getSignedDistanceDBC
    #
    # No flow conditions on sides and bottom
    #
    test = 'InvertedPerturbedBucket2D-FE'
    testProblems.append(test)
    flowTimeIntegration[test]=ScalarTransport.NoIntegration
    nciFreeSurfaceTimeIntegration[test]=HamiltonJacobi.ForwardEuler
    ciFreeSurfaceTimeIntegration[test]=ScalarTransport.ForwardEuler
    signedDistanceTimeIntegration[test]=HamiltonJacobi.BackwardEuler
    readMesh[test] = False
    runCFL[test] = 0.1
    waterDepth=0.5
    nd[test]=2
    def getBucketBC2D(xIn):
        if xIn[1] == 1.0:
            return lambda x,t: 0.0
    flowAnalyticalSolution[test]=FlowSolution_SteadyState2D(waterDepth)
    freeSurfaceAnalyticalSolution[test]=FreeSurfaceSolution_SteadyState2D(waterDepth)
    getFlowDirichletConditions[test]=getBucketBC2D
    getFreeSurfaceDirichletConditions[test]=getFreeSurfaceDBC_SteadyState2D
    flowCoefficients[test]=TwophasePotentialFlow(M1=0.0,M2=0.0,
                                                 B1=Numeric.array([0.0,0.0]),
                                                 B2=Numeric.array([0.0,0.0]),
                                                 Bcon1=Numeric.array([0.0,-K_1*dimensionless_density_1*dimensionless_gravity]),
                                                 Bcon2=Numeric.array([0.0,-K_2*dimensionless_density_2*dimensionless_gravity]),
                                                 A1=Numeric.array([[K_1,0.0],[0.0,K_1]]),
                                                 A2=Numeric.array([[K_2,0.0],[0.0,K_2]]),
                                                 C1=0.0,C2=0.0)
    flowCoefficients[test].mass=None
    flowCoefficients[test].advection='constant'
    flowCoefficients[test].diffusion='constant'
    flowCoefficients[test].potential='linear'
    flowCoefficients[test].reaction=None
    nciFreeSurfaceCoefficients[test]=TwophaseLevelSetCoefficients(1.0/porosity)
    nciFreeSurfaceCoefficients[test].mass='linear'
    nciFreeSurfaceCoefficients[test].hamiltonian='linear'
    ciFreeSurfaceCoefficients[test]=ciTwoPhaseLevelSetCoefficients(1.0/porosity)
    ciFreeSurfaceCoefficients[test].mass='linear'
    ciFreeSurfaceCoefficients[test].advection='linear'
    T[test]=2.0
    getFlowInitialConditions[test] = FlowSolution_SteadyState2D(initialWaterDepth)
    getFreeSurfaceInitialConditions[test] = invertedSinusoidalFreeSurface2D(waterDepth,surfaceWaveNumber=1,surfaceAmplitude=0.1)
    flowFullNewtonFlag[test] = False
    freeSurfaceFullNewtonFlag[test]=False
    signedDistanceFullNewtonFlag[test]=True
    flowNBC[test] = 'noFlow'
    getSignedDistanceDirichletConditions[test]=getSignedDistanceDBC
    #
    # simple weir problem
    #
    test = 'simp_lg_rsw'
    testProblems.append(test)
    flowTimeIntegration[test]=ScalarTransport.NoIntegration
    nciFreeSurfaceTimeIntegration[test]=HamiltonJacobi.BackwardEuler
    ciFreeSurfaceTimeIntegration[test]=ScalarTransport.BackwardEuler
    signedDistanceTimeIntegration[test]=HamiltonJacobi.BackwardEuler
    readMesh[test] = False
    runCFL[test] = 2.5
    waterDepth_rsw_inflow=0.5
    totalFlow_rsw = 1.0
    bottomDepth_rsw_outflow = 0.25
    nd[test]=2
    #weirStart=3.0/8.0
    #weirStop =5.0/8.0
    #def bottomProfile(x):
    #    return math.sin(2.0*math.pi*(x[0]-weirStart)/(weirStop-weirStart))**2 * (1.0-x[1]) + x[1]
    def getSimpRSW_flowDBC(xIn):
        #potential on left boundary of water
#        if xIn[0] == 0.0:
#            if xIn[1] <= waterDepth_rsw_inflow:
#                return lambda x,t: 0.0
#        if xIn[0] == 0.0:
#           if  xIn[1] >=waterDepth_rsw_inflow:
#               return lambda x,t: 0.0
        #potential on left boundary 
 	if xIn[0] == 0.0:
 	    return lambda x,t: 0.0
        #potential on top boundary
#        if xIn[1] == 1.0:
#            return lambda x,t: 0.0
        #potential on right boundary above original water surface
#	if xIn[0] == 1.0:
#	    if xIn[1] >= waterDepth_rsw_inflow:
#		return lambda x,t: 0.0
        #potential on right boundary
# 	if xIn[0] == 1.0:
# 	    if xIn[1] >= bottomDepth_rsw_outflow:
# 		return lambda x,t: 0.0
    def getSimpRSW_freeSurfaceDBC(xIn):
        if xIn[0] == 0.0:
            return lambda x,t: waterDepth_rsw_inflow-x[1]
    def getSimpRSW_signedDistanceDBC(xIn):
        pass
    flowNBC[test] = 'setFlow'
    class RSW_FluxBC:
        def __init__(self):
            self.totalFlow_rsw_f1 = totalFlow_rsw
            self.bottomDepth_rsw_outflow = bottomDepth_rsw_outflow
            self.waterDepth_rsw_inflow = waterDepth_rsw_inflow
            self.v1_inflow = self.totalFlow_rsw_f1/self.waterDepth_rsw_inflow
	    self.v2_inflow = self.v1_inflow
            self.totalFlow_rsw_f2 = self.v2_inflow*(1.0-self.waterDepth_rsw_inflow)
            #self.v2_inflow = 0.0
            #self.totalFlow_rsw_f2 = 0.0
        def setExteriorFlux(self,nExteriorElementBoundaries_global,
                            nElementBoundaryQuadraturePoints_elementBoundary,
                            nElementBoundaries_element,
                            nSpace_global,
                            exteriorElementBoundariesArray,
                            elementBoundaryElementsArray,
                            elementBoundaryLocalElementBoundariesArray,
                            x,
                            n,
                            dx_f,
                            dx_a,
                            advectiveFlux_dx_f,
                            diffusiveFlux_dx_a):
            #calculate the outflow  surface area of each phase
            f1_area=0.0
            f2_area=0.0
            for ebNE in range(nExteriorElementBoundaries_global):
                ebN = exteriorElementBoundariesArray[ebNE]
                eN_global   = elementBoundaryElementsArray[ebN,0]
                ebN_element = elementBoundaryLocalElementBoundariesArray[ebN,0]
                normals   = n[eN_global,ebN_element]
                for k in range(nElementBoundaryQuadraturePoints_elementBoundary):
                    if x[eN_global,ebN_element,k,0] == 1.0:
                        if x[eN_global,ebN_element,k,1] > self.bottomDepth_rsw_outflow:
                            if self.phi[eN_global,ebN_element,k] > 0.0:
                                f1_area += dx_a[ebN,k]
                            else:
                                f2_area += dx_a[ebN,k]
            if f1_area != 0.0:
                v1_outflow = self.totalFlow_rsw_f1/f1_area
            else:
                v1_outflow  = 0.0
#             if f2_area != 0.0:
#                 v2_outflow = self.totalFlow_rsw_f2/f2_area
#             else:
#                 v2_outflow = 0.0
            v2_outflow = v1_outflow
#            self.v2_inflow = v2_outflow
##             print "f1_area",f1_area
##             print "v1_outflow",v1_outflow
##             print "f2_area",f2_area
##             print "v2_outflow",v2_outflow
            for ebNE in range(nExteriorElementBoundaries_global):
                ebN = exteriorElementBoundariesArray[ebNE]
                eN_global   = elementBoundaryElementsArray[ebN,0]
                ebN_element = elementBoundaryLocalElementBoundariesArray[ebN,0]
                normals   = n[eN_global,ebN_element]
                for k in range(nElementBoundaryQuadraturePoints_elementBoundary):
                    if x[eN_global,ebN_element,k,0] == 1.0:
                        if x[eN_global,ebN_element,k,1] > self.bottomDepth_rsw_outflow:
                            if self.phi[eN_global,ebN_element,k] > 0.0:
                                diffusiveFlux_dx_a[ebN,k] = v1_outflow*dx_a[ebN,k]
                            else:
                                diffusiveFlux_dx_a[ebN,k] = v2_outflow*dx_a[ebN,k]
                        else:
                            diffusiveFlux_dx_a[ebN,k] = 0.0
                        advectiveFlux_dx_f[ebN,k] = 0.0
		    elif x[eN_global,ebN_element,k,0] == 0.0:
			if self.phi[eN_global,ebN_element,k] > 0.0:
			    diffusiveFlux_dx_a[ebN,k] = -self.v1_inflow*dx_a[ebN,k]
			else:
			    diffusiveFlux_dx_a[ebN,k] = -self.v2_inflow*dx_a[ebN,k]
                        advectiveFlux_dx_f[ebN,k] = 0.0
		    else:
                        diffusiveFlux_dx_a[ebN,k] = 0.0
                        advectiveFlux_dx_f[ebN,k] = 0.0
    flowNBCobject[test] = RSW_FluxBC
    flowAnalyticalSolution[test]=None
    freeSurfaceAnalyticalSolution[test]=None
    getFlowDirichletConditions[test]=getSimpRSW_flowDBC
    getFreeSurfaceDirichletConditions[test]=getSimpRSW_freeSurfaceDBC
    getSignedDistanceDirichletConditions[test]=getSimpRSW_signedDistanceDBC
    flowCoefficients[test]=TwophasePotentialFlow(M1=0.0,
                                                 M2=0.0,
                                                 B1=Numeric.array([0.0,0.0]),
                                                 B2=Numeric.array([0.0,0.0]),
                                                 Bcon1=Numeric.array([0.0,0.0]),
                                                 Bcon2=Numeric.array([0.0,0.0]),
                                                 A1=Numeric.array([[1.0,0.0],[0.0,1.0]]),
                                                 A2=Numeric.array([[1.0,0.0],[0.0,1.0]]),
                                                 C1=0.0,C2=0.0)
    flowCoefficients[test].mass=None
    flowCoefficients[test].advection=None
    flowCoefficients[test].diffusion='constant'
    flowCoefficients[test].potential='linear'
    flowCoefficients[test].reaction=None
    nciFreeSurfaceCoefficients[test]=TwophaseLevelSetCoefficients(1.0)
    nciFreeSurfaceCoefficients[test].mass='linear'
    nciFreeSurfaceCoefficients[test].hamiltonian='linear'
    ciFreeSurfaceCoefficients[test]=ciTwoPhaseLevelSetCoefficients(1.0)
    ciFreeSurfaceCoefficients[test].mass='linear'
    ciFreeSurfaceCoefficients[test].advection='linear'
    T[test]=2.0
    class rswFlowIC(SteadyState):
        def __init__(self):
            pass
        def uOfX(self,x):
            return 0.0
    class rswFreeSurfaceIC(SteadyState):
        def __init__(self,waterDepth):
            self.waterDepth = waterDepth
        def uOfX(self,x):
            return self.waterDepth - x[1]  
    getFlowInitialConditions[test] = rswFlowIC()
    getFreeSurfaceInitialConditions[test] = rswFreeSurfaceIC(waterDepth_rsw_inflow)
    flowFullNewtonFlag[test] = False
    freeSurfaceFullNewtonFlag[test]=False
    signedDistanceFullNewtonFlag[test]=True
    #
    # simple weir problem-reversed
    #
    test = 'simp_lg_rsw_r'
    testProblems.append(test)
    flowTimeIntegration[test]=ScalarTransport.NoIntegration
    nciFreeSurfaceTimeIntegration[test]=HamiltonJacobi.BackwardEuler
    ciFreeSurfaceTimeIntegration[test]=ScalarTransport.BackwardEuler
    signedDistanceTimeIntegration[test]=HamiltonJacobi.BackwardEuler
    readMesh[test] = False
    runCFL[test] = 2.5
    waterDepth_rsw_outflow_r=0.25
    waterDepth_rsw_inflow_r=0.5
    totalFlow_rsw_r = 1.0
    bottomDepth_rsw_inflow_r = 0.25
    nd[test]=2
    def getSimpRSW_flowDBC_r(xIn):
        if xIn[0] ==0.0:
            return lambda x,t: 0.0
    def getSimpRSW_freeSurfaceDBC_r(xIn):
        if xIn[0] == 0.0:
            return lambda x,t: min(waterDepth_rsw_inflow_r - x[1],x[1] - bottomDepth_rsw_inflow_r)
    def getSimpRSW_signedDistanceDBC_r(xIn):
        pass
    flowNBC[test] = 'setFlow'
    class RSW_FluxBC_r:
        def __init__(self):
            self.totalFlow_rsw_f1 = totalFlow_rsw_r
            self.bottomDepth_rsw_inflow = bottomDepth_rsw_inflow_r
            self.waterDepth_rsw_inflow  = waterDepth_rsw_inflow_r
            v1_inflow = self.totalFlow_rsw_f1/(self.waterDepth_rsw_inflow - self.bottomDepth_rsw_inflow)
            v2_inflow = v1_inflow
            #self.totalFlow_rsw_f2 = v2_inflow*(1.0-self.waterDepth_rsw_inflow)
            self.totalFlow_rsw_f2 = 0.0
        def setExteriorFlux(self,nExteriorElementBoundaries_global,
                            nElementBoundaryQuadraturePoints_elementBoundary,
                            nElementBoundaries_element,
                            nSpace_global,
                            exteriorElementBoundariesArray,
                            elementBoundaryElementsArray,
                            elementBoundaryLocalElementBoundariesArray,
                            x,
                            n,
                            dx_f,
                            dx_a,
                            advectiveFlux_dx_f,
                            diffusiveFlux_dx_a):
            #calculate the outflow  surface area of each phase
            f1_area=0.0
            f2_area=0.0
            for ebNE in range(nExteriorElementBoundaries_global):
                ebN = exteriorElementBoundariesArray[ebNE]
                eN_global   = elementBoundaryElementsArray[ebN,0]
                ebN_element = elementBoundaryLocalElementBoundariesArray[ebN,0]
                normals   = n[eN_global,ebN_element]
                for k in range(nElementBoundaryQuadraturePoints_elementBoundary):
                    X = x[eN_global,ebN_element,k]
                    if X[0] == 1.0:
                        if self.phi[eN_global,ebN_element,k] >= 0.0:
                            f1_area += dx_a[ebN,k]
                        else:
                            f2_area += dx_a[ebN,k]
            if f1_area > 0.0:
                v1_outflow = self.totalFlow_rsw_f1/f1_area
            else:
                v1_outflow  = 0.0
            if f2_area > 0.0:
                v2_outflow = self.totalFlow_rsw_f2/f2_area
            else:
                v2_outflow = 0.0
            for ebNE in range(nExteriorElementBoundaries_global):
                ebN = exteriorElementBoundariesArray[ebNE]
                eN_global   = elementBoundaryElementsArray[ebN,0]
                ebN_element = elementBoundaryLocalElementBoundariesArray[ebN,0]
                normals   = n[eN_global,ebN_element]
                for k in range(nElementBoundaryQuadraturePoints_elementBoundary):
                    X = x[eN_global,ebN_element,k]
                    if X[0] == 1.0:
                        if self.phi[eN_global,ebN_element,k] >= 0.0:
                            diffusiveFlux_dx_a[ebN,k] = v1_outflow*dx_a[ebN,k]
                        else:
                            diffusiveFlux_dx_a[ebN,k] = v2_outflow*dx_a[ebN,k]
                        advectiveFlux_dx_f[ebN,k] = 0.0
                    else:
                        diffusiveFlux_dx_a[ebN,k] = 0.0
                        advectiveFlux_dx_f[ebN,k] = 0.0
    flowNBCobject[test] = RSW_FluxBC_r
    flowAnalyticalSolution[test]=None
    freeSurfaceAnalyticalSolution[test]=None
    getFlowDirichletConditions[test]=getSimpRSW_flowDBC_r
    getFreeSurfaceDirichletConditions[test]=getSimpRSW_freeSurfaceDBC_r
    getSignedDistanceDirichletConditions[test]=getSimpRSW_signedDistanceDBC_r
    flowCoefficients[test]=TwophasePotentialFlow(M1=0.0,
                                                 M2=0.0,
                                                 B1=Numeric.array([0.0,0.0]),
                                                 B2=Numeric.array([0.0,0.0]),
                                                 Bcon1=Numeric.array([0.0,0.0]),
                                                 Bcon2=Numeric.array([0.0,0.0]),
                                                 A1=Numeric.array([[1.0,0.0],[0.0,1.0]]),
                                                 A2=Numeric.array([[1.0,0.0],[0.0,1.0]]),
                                                 C1=0.0,C2=0.0)
    flowCoefficients[test].mass=None
    flowCoefficients[test].advection=None
    flowCoefficients[test].diffusion='constant'
    flowCoefficients[test].potential='linear'
    flowCoefficients[test].reaction=None
    nciFreeSurfaceCoefficients[test]=TwophaseLevelSetCoefficients(1.0)
    nciFreeSurfaceCoefficients[test].mass='linear'
    nciFreeSurfaceCoefficients[test].hamiltonian='linear'
    ciFreeSurfaceCoefficients[test]=ciTwoPhaseLevelSetCoefficients(1.0)
    ciFreeSurfaceCoefficients[test].mass='linear'
    ciFreeSurfaceCoefficients[test].advection='linear'
    T[test]=2.0
    class rswFlowIC_r(SteadyState):
        def __init__(self):
            pass
        def uOfX(self,x):
            return 0.0
    class rswFreeSurfaceIC_r(SteadyState):
        def __init__(self,waterDepth,inflowBottomDepth):
            self.waterDepth = waterDepth
            self.inflowBottomDepth = inflowBottomDepth
        def uOfX(self,x):
            return min(self.waterDepth - x[1],x[1] - self.inflowBottomDepth)  
    getFlowInitialConditions[test] = rswFlowIC_r()
    getFreeSurfaceInitialConditions[test] = rswFreeSurfaceIC_r(waterDepth_rsw_inflow_r,bottomDepth_rsw_inflow_r)
    flowFullNewtonFlag[test] = False
    freeSurfaceFullNewtonFlag[test]=False
    signedDistanceFullNewtonFlag[test]=True
    #
    #
    #
    for nTest,test in enumerate(testProblems):
        print nTest,test
    for test in testProblems[:6]:
        nLevels=3
        tolFac = 0.01
        linTolFac = 0.001
        DG = False
	CI = False#True
	polyOrder = 1
        tOrder = 1
	redistance = True
        heavySideSmoothingFactor = 1.5
        signedDistanceSmoothingFactor = 1.5#1.0e-8
        if DG:
	    if polyOrder == 1:
		FemSpaceTypeCon = DG_AffineLinearOnSimplexWithNodalBasis
		FemSpaceTypeNonCon = C0_AffineLinearOnSimplexWithNodalBasis
	    elif polyOrder == 2:
		FemSpaceTypeCon = DG_AffineQuadraticOnSimplexWithNodalBasis
		FemSpaceTypeNonCon = C0_AffineQuadraticOnSimplexWithNodalBasis
	    else:
		print "polyOrder should be 1 or 2"
		return
            #
            # flow
            #
            flowConservativeFlux = None
            flowNumericalFlux = True
            flowStabilization= None
            flowShockCapturing=None
	    flowShockCapturingDiffusion=0.0
            #
            # ci free surface
            #
            ciFreeSurfaceConservativeFlux = None
            ciFreeSurfaceNumericalFlux = True
            ciFreeSurfaceStabilization= None
            ciFreeSurfaceShockCapturing=None
	    ciFreeSurfaceShockCapturingDiffusion=0.0
            #
            #neither  of the next to are  implemented as DG yet
            #
            #
            # nci free surface
            #
            nciFreeSurfaceConservativeFlux = None
            nciFreeSurfaceNumericalFlux = None
            nciFreeSurfaceStabilization= None
            nciFreeSurfaceShockCapturing=None
	    nciFreeSurfaceShockCapturingDiffusion=0.0
            #
            # signed distance
            #
            signedDistanceNumericalFlux = None
            signedDistanceStabilization= None
            signedDistanceShockCapturing=None
	    signedDistanceShockCapturingDiffusion=0.0
        else:
	    if polyOrder == 1:
		FemSpaceTypeCon = C0_AffineLinearOnSimplexWithNodalBasis
		FemSpaceTypeNonCon = C0_AffineLinearOnSimplexWithNodalBasis
	    elif polyOrder == 2:
		FemSpaceTypeCon = C0_AffineLinearOnSimplexWithNodalBasis
		FemSpaceTypeNonCon = C0_AffineQuadraticOnSimplexWithNodalBasis
	    else:
		print "polyOrder should be 1 or 2"
		return
            #
            #flow
            #
            #flowConservativeFlux = 'pwl'
            flowConservativeFlux = None
            flowNumericalFlux = None
            flowStabilization= None
            flowShockCapturing= None
            flowShockCapturingDiffusion = 0.0
            #
            # ci free surface
            #
            ciFreeSurfaceStabilization='2'
	    ciFreeSurfaceShockCapturing= '2'
            ciFreeSurfaceShockCapturingDiffusion = 0.09
	    ciFreeSurfaceConservativeFlux = None
	    ciFreeSurfaceNumericalFlux = None
            #
            # nci free surface
            #
            nciFreeSurfaceStabilization='2'
            nciFreeSurfaceShockCapturing= '2'
            nciFreeSurfaceShockCapturingDiffusion = 0.09
	    nciFreeSurfaceConservativeFlux = None
	    nciFreeSurfaceNumericalFlux = None
            #
            # signed distance
            #
            signedDistanceStabilization='2'
            signedDistanceShockCapturing= '2'
            signedDistanceShockCapturingDiffusion = 0.09
	    signedDistanceConservativeFlux = None
	    signedDistanceNumericalFlux = None
        quadratureOrder=3
        computeSpaceTimeError=False
        nn=3 #number of nodes on the coarsest mesh
        print "Starting Test "+`test`
        print "Setting up flow quadrature"
        flowQuadrature={}
        gq = SimplexGaussQuadrature(nd[test])
        #gq = SimplexLobattoQuadrature(nd[test])
        gq.setOrder(quadratureOrder)
        for integral in ScalarTransport.OneLevelScalarTransport.integralKeys:
            flowQuadrature[integral] = gq 
        if flowStabilization != None:
            flowQuadrature['stab'] = gq 
        if flowShockCapturing != None:
            flowQuadrature['numDiff'] = gq 
        flowElementBoundaryQuadrature={}
        ebgq = SimplexGaussQuadrature(nd[test]-1)
        ebgq.setOrder(quadratureOrder)
        for elementBoundaryIntegral in ScalarTransport.OneLevelScalarTransport.elementBoundaryIntegralKeys:
            flowElementBoundaryQuadrature[elementBoundaryIntegral] = ebgq 
        print "Setting up CI free surface quadrature"
        ciFreeSurfaceQuadrature={}
        gq = SimplexGaussQuadrature(nd[test])
        #gq = SimplexLobattoQuadrature(nd[test])
        gq.setOrder(quadratureOrder)
        for integral in ScalarTransport.OneLevelScalarTransport.integralKeys:
            ciFreeSurfaceQuadrature[integral] = gq 
        if ciFreeSurfaceStabilization != None:
            ciFreeSurfaceQuadrature['stab'] = gq 
        if ciFreeSurfaceShockCapturing != None:
            ciFreeSurfaceQuadrature['numDiff'] = gq 
        ciFreeSurfaceElementBoundaryQuadrature={}
        ebgq = SimplexGaussQuadrature(nd[test]-1)
        ebgq.setOrder(quadratureOrder)
        for elementBoundaryIntegral in ScalarTransport.OneLevelScalarTransport.elementBoundaryIntegralKeys:
            ciFreeSurfaceElementBoundaryQuadrature[elementBoundaryIntegral] = ebgq 
        print "Setting up NCI free surface quadrature"
        nciFreeSurfaceQuadrature={}
        gq.setOrder(quadratureOrder)
        for integral in HamiltonJacobi.OneLevelHamiltonJacobi.integralKeys:
            nciFreeSurfaceQuadrature[integral] = gq 
        nciFreeSurfaceElementBoundaryQuadrature={}
        ebgq = SimplexGaussQuadrature(nd[test]-1)
        ebgq.setOrder(quadratureOrder)
        for elementBoundaryIntegral in HamiltonJacobi.OneLevelHamiltonJacobi.elementBoundaryIntegralKeys:
            nciFreeSurfaceElementBoundaryQuadrature[elementBoundaryIntegral] = ebgq 
        print "Setting up signed distance quadrature"
        signedDistanceQuadrature={}
        gq.setOrder(quadratureOrder)
        for integral in HamiltonJacobi.OneLevelHamiltonJacobi.integralKeys:
            signedDistanceQuadrature[integral] = gq 
        signedDistanceElementBoundaryQuadrature={}
        ebgq = SimplexGaussQuadrature(nd[test]-1)
        ebgq.setOrder(quadratureOrder)
        for elementBoundaryIntegral in HamiltonJacobi.OneLevelHamiltonJacobi.elementBoundaryIntegralKeys:
            signedDistanceElementBoundaryQuadrature[elementBoundaryIntegral] = ebgq 
        #
        #define the mesh hierarchy
        #
        print "Setting up MultilevelMesh"
        mlMesh = None
        mlMeshFileName = "mlMesh%dD.%d" % (nd[test],nLevels)
        try:
            mlMeshFile = open(mlMeshFileName,'rb')
            print "reading mesh"
            mlMesh = cPickle.load(mlMeshFile)
            print "done reading mesh"
        except:
            print "generating mesh"
            mlMeshFile = open(mlMeshFileName,'wb')
            if nd[test]==1:
                mlMesh = MultilevelEdgeMesh(nn,1,1,refinementLevels=nLevels)
            elif nd[test]==2:
                mlMesh = MultilevelTriangularMesh(nn,nn,1,
                                                  refinementLevels=nLevels)
            elif nd[test]==3:
                mlMesh = MultilevelTetrahedralMesh(nn,nn,nn,
                                                   refinementLevels=nLevels)
            cPickle.dump(mlMesh,mlMeshFile,protocol=cPickle.HIGHEST_PROTOCOL)
            print "done generating mesh"        
        print "Setting up MultilevelScalarTransport"
	#matType = 'dense'
	matType = 'csr'
        tolList=[]
        linTolList=[]
        matlab = os.popen('/Applications/MATLAB71/bin/matlab -nosplash -nodesktop -nojvm','w')
#         (matlab,matlabout) = os.popen2('/Applications/MATLAB71/bin/matlab -nosplash -nodesktop -nojvm','w')
#         #(matlab,matlabout,matlaberr) = os.popen3('/Applications/MATLAB71/bin/matlab -nosplash -nodesktop -nojvm','w')
        for l in range(nLevels):
            mlMesh.meshList[l].computeGeometricInfo()
            tolList.append(tolFac*(mlMesh.meshList[l].h**2))
            linTolList.append(linTolFac*(mlMesh.meshList[l].h**2))
            matlab.write("vertices_%s = [" % l)
            for nN in range(mlMesh.meshList[l].nNodes_global-1):
                matlab.write(`mlMesh.meshList[l].nodeArray[nN,0]`+","+
                             `mlMesh.meshList[l].nodeArray[nN,1]`+","+
                             `mlMesh.meshList[l].nodeArray[nN,2]`+";")
            matlab.write(`mlMesh.meshList[l].nodeArray[nN+1,0]`+","+
                         `mlMesh.meshList[l].nodeArray[nN+1,1]`+","+
                         `mlMesh.meshList[l].nodeArray[nN+1,2]`+"];\n")
            matlab.write("faces_%s = [" % l)
            for eN in range(mlMesh.meshList[l].nElements_global-1):
                for nN in range(mlMesh.meshList[l].nNodes_element-1):
                    matlab.write(`mlMesh.meshList[l].elementNodesArray[eN,nN]+1`+",")
                matlab.write(`mlMesh.meshList[l].elementNodesArray[eN,nN+1]+1`+";")
            for nN in range(mlMesh.meshList[l].nNodes_element-1):
                matlab.write(`mlMesh.meshList[l].elementNodesArray[eN+1,nN]+1`+",")
            matlab.write(`mlMesh.meshList[l].elementNodesArray[eN+1,nN+1]+1`+"];\n")
            matlab.write("figure; title(\'Mesh, level = %d\');patch(\'Vertices\',vertices_%s,\'Faces\',faces_%s,\'FaceColor\',\'none\',\'EdgeColor\',\'red\')\n" % (l,l,l))
            matlab.flush()
        atol = min(tolList)
        lin_atol = min(linTolList)
        flow = ScalarTransport.MultilevelScalarTransport(
            nd[test],
            mlMesh,
            FemSpaceTypeCon,
            FemSpaceTypeCon,
            matType,
            getFlowDirichletConditions[test],
            flowCoefficients[test],
            flowQuadrature,
            flowElementBoundaryQuadrature,
            fluxBoundaryConditions=flowNBC[test],
            stabilization = flowStabilization,
            shockCapturing = flowShockCapturing,
            shockCapturingDiffusion = flowShockCapturingDiffusion,
            conservativeFlux = flowConservativeFlux,
            numericalFlux = flowNumericalFlux,
            TimeIntegrationClass = flowTimeIntegration[test],
	    tOrder = tOrder)
	if CI == True:
	    freeSurface = ScalarTransport.MultilevelScalarTransport(
		nd[test],
		mlMesh,
		FemSpaceTypeCon,
		FemSpaceTypeCon,
		matType,
		getFreeSurfaceDirichletConditions[test],
		ciFreeSurfaceCoefficients[test],
		ciFreeSurfaceQuadrature,
		ciFreeSurfaceElementBoundaryQuadrature,
		fluxBoundaryConditions='outFlow',
		stabilization  = ciFreeSurfaceStabilization,
		shockCapturing = ciFreeSurfaceShockCapturing,
		shockCapturingDiffusion = ciFreeSurfaceShockCapturingDiffusion,
		conservativeFlux = ciFreeSurfaceConservativeFlux,
		numericalFlux = ciFreeSurfaceNumericalFlux,
		TimeIntegrationClass = ciFreeSurfaceTimeIntegration[test],
		tOrder = tOrder)
	else:
	    freeSurface = HamiltonJacobi.MultilevelHamiltonJacobi(
		nd[test],
		mlMesh,
		FemSpaceTypeNonCon,
		FemSpaceTypeNonCon,
		matType,
		getFreeSurfaceDirichletConditions[test],
		coefficients = nciFreeSurfaceCoefficients[test],
		quadrature = nciFreeSurfaceQuadrature,
                elementBoundaryQuadrature = nciFreeSurfaceElementBoundaryQuadrature,
		stabilization = nciFreeSurfaceStabilization,
		shockCapturing = nciFreeSurfaceShockCapturing,
		shockCapturingDiffusion = nciFreeSurfaceShockCapturingDiffusion,
		TimeIntegrationClass = nciFreeSurfaceTimeIntegration[test],
                calculateElementBoundaryValues = True,
		tOrder=tOrder)
        if flowNBC[test] == 'setFlow':
            for l in range(nLevels):           
                flow.modelList[l].flowBC_setter = flowNBCobject[test]()
                flow.modelList[l].flowBC_setter.phi = freeSurface.modelList[l].ebq['u']
	signedDistance = HamiltonJacobi.MultilevelHamiltonJacobi(
            nd[test],
            mlMesh,
            FemSpaceTypeNonCon,
            FemSpaceTypeNonCon,
            matType,
            getSignedDistanceDirichletConditions[test],
            coefficients = TwophaseSignedDistanceCoefficients(),
            quadrature = signedDistanceQuadrature,
            elementBoundaryQuadrature = signedDistanceElementBoundaryQuadrature,
            stabilization = signedDistanceStabilization,
            shockCapturing = signedDistanceShockCapturing,
            shockCapturingDiffusion = signedDistanceShockCapturingDiffusion,
            TimeIntegrationClass = signedDistanceTimeIntegration[test],
            calculateElementBoundaryValues = True,
	    tOrder = tOrder)
        for l in range(nLevels):
            freeSurface.modelList[l].coefficients.initializeVelocity(
                [Numeric.reshape(flow.modelList[l].q['velocity'],
                                 (flow.modelList[l].nQuadraturePoints_global,
                                  flow.modelList[l].nSpace_global)),
                 Numeric.reshape(flow.modelList[l].ebq['velocity'],
                                 (flow.modelList[l].nElementBoundaryQuadraturePoints_global,
                                  flow.modelList[l].nSpace_global))])
            flow.modelList[l].coefficients.initializeFreeSurface(
                [freeSurface.modelList[l].q['u'].flat,
                 freeSurface.modelList[l].ebq['u'].flat],
                heavySideSmoothingFactor*freeSurface.modelList[l].mesh.h)
            signedDistance.modelList[l].coefficients.initializeSignFunction([freeSurface.modelList[l].q['u'].flat,
									     freeSurface.modelList[l].ebq['u'].flat],
                                                                            [Numeric.reshape(freeSurface.modelList[l].q['grad(u)'].flat,
                                                                                             (freeSurface.modelList[l].nQuadraturePoints_global,
                                                                                              freeSurface.modelList[l].nSpace_global)),
                                                                             Numeric.reshape(freeSurface.modelList[l].ebq['grad(u)'].flat,
                                                                                             (freeSurface.modelList[l].nElementBoundaryQuadraturePoints_global,
                                                                                              freeSurface.modelList[l].nSpace_global))],
                                                                            signedDistanceSmoothingFactor*freeSurface.modelList[l].mesh.h)
        print "Setting up Flow Solvers"
	print "Setting up Linear Solver"
	(flowLinearSolver,directSolverFlag) = multilevelLinearSolverChooser(
	    linearOperatorList = flow.jacobianList,
	    multilevelLinearSolverType = 'SparseLU',
	    #multilevelLinearSolverType = 'DenseLU',#'NI','MGM','Jacobi','GaussSeidel','StarILU'
            computeSolverRates=True,
            printSolverInfo=False,
	    levelLinearSolverType = 'SparseLU',#'DenseLU','MGM','Jacobi','GaussSeidel','StarILU'
            computeLevelSolverRates=True,
            printLevelSolverInfo=False,
	    smootherType = None,#'Jacobi','GaussSeidel','StarILU'
            computeSmootherRates=True,
            printSmootherInfo=False,
	    prolongList = flow.meshTransfers.prolongList,
	    restrictList = flow.meshTransfers.restrictList,
	    connectivityListList = [flow.modelList[l].freeNodeStarList for l in range(nLevels)],
	    relativeToleranceList = linTolList,
	    absoluteTolerance = min(linTolList),
	    solverMaxIts = 50,
	    cycles=3,
	    preSmooths=3,
	    postSmooths=3)
        print "Setting up NonlinearSolver"
	flowNonlinearSolver = multilevelNonlinearSolverChooser(
	    nonlinearOperatorList = flow.modelList,
	    jacobianList = flow.jacobianList,
	    multilevelNonlinearSolverType = 'NLNI',#'Newton','FAS','NLJacobi','NLGaussSeidel','NLJacobi'
	    #multilevelNonlinearSolverType = 'Newton',#'FAS','NLJacobi','NLGaussSeidel','NLJacobi'
            computeSolverRates=True,
            printSolverInfo=False,
	    relativeToleranceList = tolList,
	    absoluteTolerance = min(tolList),
	    levelNonlinearSolverType='Newton',#'FAS','NLJacobi','NLGaussSeidel','NLJacobi'
            computeLevelSolverRates=True,
	    printLevelSolverInfo=False,
	    smootherType = None,#'NLJacobi','NLGaussSeidel','NLStarILU'
            computeSmootherRates=True,
	    printSmootherInfo=False,
	    preSmooths=3,
	    postSmooths=3,
	    cycles=3,
	    maxSolverIts=50,
	    prolong_bcList = flow.meshTransfers.prolong_bcList,
	    restrict_bcList = flow.meshTransfers.restrict_bcList,
	    restrict_bcSumList = flow.meshTransfers.restrict_bcSumList,
	    prolongList = flow.meshTransfers.prolongList,
	    restrictList = flow.meshTransfers.restrictList,
	    restrictionRowSumList = flow.meshTransfers.restrictSumList,
	    connectionListList=[flow.modelList[l].freeNodeStarList for l in range(nLevels)],
	    linearSolverList=flowLinearSolver.solverList,
	    linearDirectSolverFlag=directSolverFlag,
	    solverFullNewtonFlag=flowFullNewtonFlag[test],
	    levelSolverFullNewtonFlag=flowFullNewtonFlag[test],
	    smootherFullNewtonFlag=flowFullNewtonFlag[test])
        print "Setting up FreeSurface Transport Solvers"
	print "Setting up Linear Solver"
	(freeSurfaceLinearSolver,directSolverFlag) = multilevelLinearSolverChooser(
	    linearOperatorList = freeSurface.jacobianList,
	    multilevelLinearSolverType = 'SparseLU',
	    #multilevelLinearSolverType = 'DenseLU',#'NI','MGM','Jacobi','GaussSeidel','StarILU'
            computeSolverRates=True,
            printSolverInfo=False,
	    levelLinearSolverType = 'SparseLU',#'DenseLU','MGM','Jacobi','GaussSeidel','StarILU'
            computeLevelSolverRates=True,
            printLevelSolverInfo=False,
	    smootherType = None,#'Jacobi','GaussSeidel','StarILU'
            computeSmootherRates=True,
            printSmootherInfo=False,
	    prolongList = freeSurface.meshTransfers.prolongList,
	    restrictList = freeSurface.meshTransfers.restrictList,
	    connectivityListList = [freeSurface.modelList[l].freeNodeStarList for l in range(nLevels)],
	    relativeToleranceList = linTolList,
	    absoluteTolerance = min(linTolList),
	    solverMaxIts = 50,
	    cycles=3,
	    preSmooths=3,
	    postSmooths=3)
        print "Setting up NonlinearSolver"
	freeSurfaceNonlinearSolver = multilevelNonlinearSolverChooser(
	    nonlinearOperatorList = freeSurface.modelList,
	    jacobianList = freeSurface.jacobianList,
	    multilevelNonlinearSolverType = 'NLNI',#'Newton','FAS','NLJacobi','NLGaussSeidel','NLJacobi'
	    #multilevelNonlinearSolverType = 'Newton',#'FAS','NLJacobi','NLGaussSeidel','NLJacobi'
            computeSolverRates=True,
            printSolverInfo=False,
	    relativeToleranceList = tolList,
	    absoluteTolerance = min(tolList),
	    levelNonlinearSolverType='Newton',#'FAS','NLJacobi','NLGaussSeidel','NLJacobi'
            computeLevelSolverRates=True,
	    printLevelSolverInfo=False,
	    smootherType = None,#'NLJacobi','NLGaussSeidel','NLStarILU'
            computeSmootherRates=True,
	    printSmootherInfo=False,
	    preSmooths=3,
	    postSmooths=3,
	    cycles=3,
	    maxSolverIts=50,
	    prolong_bcList = freeSurface.meshTransfers.prolong_bcList,
	    restrict_bcList = freeSurface.meshTransfers.restrict_bcList,
	    restrict_bcSumList = freeSurface.meshTransfers.restrict_bcSumList,
	    prolongList = freeSurface.meshTransfers.prolongList,
	    restrictList = freeSurface.meshTransfers.restrictList,
	    restrictionRowSumList = freeSurface.meshTransfers.restrictSumList,
	    connectionListList=[freeSurface.modelList[l].freeNodeStarList for l in range(nLevels)],
	    linearSolverList=freeSurfaceLinearSolver.solverList,
	    linearDirectSolverFlag=directSolverFlag,
	    solverFullNewtonFlag=freeSurfaceFullNewtonFlag[test],
	    levelSolverFullNewtonFlag=freeSurfaceFullNewtonFlag[test],
	    smootherFullNewtonFlag=freeSurfaceFullNewtonFlag[test])
        print "Setting up Signed Distance Solvers"
	print "Setting up Linear Solver"
	(signedDistanceLinearSolver,directSolverFlag) = multilevelLinearSolverChooser(
	    linearOperatorList = signedDistance.jacobianList,
	    multilevelLinearSolverType = 'SparseLU',
	    #multilevelLinearSolverType = 'DenseLU',#'NI','MGM','Jacobi','GaussSeidel','StarILU'
            computeSolverRates=True,
            printSolverInfo=False,
	    levelLinearSolverType = 'SparseLU',#'DenseLU','MGM','Jacobi','GaussSeidel','StarILU'
            computeLevelSolverRates=True,
            printLevelSolverInfo=False,
	    smootherType = None,#'Jacobi','GaussSeidel','StarILU'
            computeSmootherRates=True,
            printSmootherInfo=False,
	    prolongList = signedDistance.meshTransfers.prolongList,
	    restrictList = signedDistance.meshTransfers.restrictList,
	    connectivityListList = [signedDistance.modelList[l].freeNodeStarList for l in range(nLevels)],
	    relativeToleranceList = linTolList,
	    absoluteTolerance = min(linTolList),
	    solverMaxIts = 50,
	    cycles=3,
	    preSmooths=3,
	    postSmooths=3)
        print "Setting up NonlinearSolver"
	signedDistanceNonlinearSolver = multilevelNonlinearSolverChooser(
	    nonlinearOperatorList = signedDistance.modelList,
	    jacobianList = signedDistance.jacobianList,
	    multilevelNonlinearSolverType = 'NLNI',#'Newton','FAS','NLJacobi','NLGaussSeidel','NLJacobi'
	    #multilevelNonlinearSolverType = 'Newton',#'FAS','NLJacobi','NLGaussSeidel','NLJacobi'
            computeSolverRates=True,
            printSolverInfo=False,
	    relativeToleranceList = tolList,
	    absoluteTolerance = min(tolList),
	    levelNonlinearSolverType='Newton',#'FAS','NLJacobi','NLGaussSeidel','NLJacobi'
            computeLevelSolverRates=True,
	    printLevelSolverInfo=False,
	    smootherType = None,#'NLJacobi','NLGaussSeidel','NLStarILU'
            computeSmootherRates=True,
	    printSmootherInfo=False,
	    preSmooths=3,
	    postSmooths=3,
	    cycles=3,
	    maxSolverIts=50,
	    prolong_bcList = signedDistance.meshTransfers.prolong_bcList,
	    restrict_bcList = signedDistance.meshTransfers.restrict_bcList,
	    restrict_bcSumList = signedDistance.meshTransfers.restrict_bcSumList,
	    prolongList = signedDistance.meshTransfers.prolongList,
	    restrictList = signedDistance.meshTransfers.restrictList,
	    restrictionRowSumList = signedDistance.meshTransfers.restrictSumList,
	    connectionListList=[signedDistance.modelList[l].freeNodeStarList for l in range(nLevels)],
	    linearSolverList=signedDistanceLinearSolver.solverList,
	    linearDirectSolverFlag=directSolverFlag,
	    solverFullNewtonFlag=signedDistanceFullNewtonFlag[test],
	    levelSolverFullNewtonFlag=signedDistanceFullNewtonFlag[test],
	    smootherFullNewtonFlag=signedDistanceFullNewtonFlag[test])
        print "Running Solver"
        try:
            os.mkdir(test)
        except:
            pass
        os.chdir(test)
        tn = 0.0
	redistance_tn = 0.0
        nSteps = 0
        freeSurface.setInitialConditions(getFreeSurfaceInitialConditions[test],tn)
	signedDistance.setInitialConditions(getFreeSurfaceInitialConditions[test],tn)
        flow.setInitialConditions(getFlowInitialConditions[test],tn)
	for l in range(nLevels):
	    if polyOrder == 1:
		mlMesh.meshList[l].writeMeshEnsight('flow'+`l`,'flow'+`l`)
		mlMesh.meshList[l].writeMeshEnsight('freeSurface'+`l`,'freeSurface'+`l`)
		mlMesh.meshList[l].writeMeshEnsight('signedDistance'+`l`,'signedDistance'+`l`)
	    elif polyOrder == 2:
		flow.modelList[l].trialSpace.writeMeshEnsight('flow'+`l`,'flow'+`l`)
		freeSurface.modelList[l].trialSpace.writeMeshEnsight('freeSurface'+`l`,'freeSurface'+`l`)
		signedDistance.modelList[l].trialSpace.writeMeshEnsight('signedDistance'+`l`,'signedDistance'+`l`)
	    flow.modelList[l].u.name='p'
	    flow.modelList[l].trialSpace.writeFunctionEnsight(flow.modelList[l].u,'flow'+`l`,append=False)
	    freeSurface.modelList[l].u.name='phi'
	    freeSurface.modelList[l].trialSpace.writeFunctionEnsight(freeSurface.modelList[l].u,'freeSurface'+`l`,append=False)
	    signedDistance.modelList[l].u.name='phi'
	    signedDistance.modelList[l].trialSpace.writeFunctionEnsight(signedDistance.modelList[l].u,'signedDistance'+`l`,append=False)
	redistanceTimeValues=[0.0]
	if redistance:
            redistanceLoop()
        for l in range(nLevels):
	    flow.modelList[l].coefficients.updateFreeSurface([freeSurface.modelList[l].q['u'].flat,
							      freeSurface.modelList[l].ebq['u'].flat])
	flowNonlinearSolver.solveMultilevel(uList = flow.uList,
					    rList = flow.rList)
	for l in range(nLevels):
	    flow.modelList[l].getFlowVelocity()
	    if flowConservativeFlux == 'pwl':
		flow.modelList[l].getConservationFluxPWL()
	    elif flowConservativeFlux == 'pwc':
		flow.modelList[l].getConservationFluxPWC()
	    else:
		flow.modelList[l].getFlowVelocityElementBoundaries()
	    freeSurface.modelList[l].coefficients.updateVelocity([Numeric.reshape(flow.modelList[l].q['velocity'],
										  (flow.modelList[l].nQuadraturePoints_global,
										   flow.modelList[l].nSpace_global)),
								  Numeric.reshape(flow.modelList[l].ebq['velocity'],
										  (flow.modelList[l].nElementBoundaryQuadraturePoints_global,
										   flow.modelList[l].nSpace_global))])
	    freeSurface.modelList[l].updateCoefficients()
            matlab.write("phi_%s = [" % l)
            for nN in range(flow.modelList[l].mesh.nNodes_global-1):
                matlab.write(`flow.modelList[l].u.dof[nN]`+";")
            matlab.write(`flow.modelList[l].u.dof[nN+1]`+"] \n")
            matlab.write("figure; title(\'potential, level = %s\'); patch(\'Vertices\',vertices_%s,\'Faces\',faces_%s,\'FaceVertexCData\',phi_%s,\'FaceColor\',\'interp\')\n" % (l,l,l,l))
            matlab.flush()
            matlab.write("u_%s = [" % l)
            for nN in range(freeSurface.modelList[l].mesh.nNodes_global-1):
                matlab.write(`freeSurface.modelList[l].u.dof[nN]`+";")
            matlab.write(`freeSurface.modelList[l].u.dof[nN+1]`+"] \n")
            matlab.write("figure; title(\'free surface, level = %s\'); patch(\'Vertices\',vertices_%s,\'Faces\',faces_%s,\'FaceVertexCData\',u_%s,\'FaceColor\',\'interp\')\n" % (l,l,l,l))
            matlab.flush()
        import sys
        tstring=None
        eraseTime='\b\b\b\b\b\b\b\b\b\b\b\b'
	DTSET=None
	timeValues = [tn]
        freeSurface.modelList[-1].timeIntegration.runCFL = runCFL[test]
        signedDistance.modelList[-1].timeIntegration.runCFL = 0.9
        while (tn < T[test]):
	    freeSurface.chooseDT(DTSET)
	    if tn + freeSurface.DT > T[test]:
		freeSurface.chooseDT(DTSET=T[test]-tn)
	    freeSurface.initializeTimeIntegration()
	    tn += freeSurface.DT
#  	    if tstring != None:
#  		sys.stdout.write(eraseTime)
#  	    else:
#  		sys.stdout.write('T = %12.5e, tn = ' % T[test])
            sys.stdout.write('T = %12.5e, tn = ' % T[test])
# 	    tstring='%12.5e' % (tn,)
	    tstring='%12.5e\n' % (tn,)
	    sys.stdout.write(tstring)
	    sys.stdout.flush()
	    nSteps += 1
	    testOut = test + ('%4.4i' % nSteps)
            if CI == True:
                for l in range(nLevels):
                    freeSurface.modelList[l].setInflowFlux()
	    freeSurfaceNonlinearSolver.updateJacobian()
	    freeSurfaceNonlinearSolver.solveMultilevel(uList=freeSurface.uList,rList=freeSurface.rList)
            freeSurface.updateTimeHistory()
	    for l in range(nLevels-1,0,-1):
		freeSurface.meshTransfers.scaled_restrict_bcList[l].matvec(freeSurface.modelList[l].u.dof,
									   freeSurface.modelList[l-1].u.dof)
	    if redistance:
                redistanceLoop()
		for l in range(nLevels-1,0,-1):
		    freeSurface.meshTransfers.scaled_restrict_bcList[l].matvec(signedDistance.modelList[l].u.dof,
									       signedDistance.modelList[l-1].u.dof)
                    freeSurface.modelList[l].u.dof[:] = signedDistance.modelList[l].u.dof
                    freeSurface.modelList[l].setFreeDOF(freeSurface.uList[l])
            freeSurface.modelList[l].updateCoefficients()
	    for l in range(nLevels):
		flow.modelList[l].coefficients.updateFreeSurface([freeSurface.modelList[l].q['u'].flat,
								  freeSurface.modelList[l].ebq['u'].flat])
            flowNonlinearSolver.updateJacobian()
	    flowNonlinearSolver.solveMultilevel(uList = flow.uList,rList = flow.rList)
	    for l in range(nLevels):
		flow.modelList[l].getFlowVelocity()
		if flowConservativeFlux == 'pwl':
		    flow.modelList[l].getConservationFluxPWL()
		elif flowConservativeFlux == 'pwc':
		    flow.modelList[l].getConservationFluxPWC()
		else:
		    flow.modelList[l].getFlowVelocityElementBoundaries()
		freeSurface.modelList[l].coefficients.updateVelocity([Numeric.reshape(flow.modelList[l].q['velocity'],
										      (flow.modelList[l].nQuadraturePoints_global,
										       flow.modelList[l].nSpace_global)),
								      Numeric.reshape(flow.modelList[l].ebq['velocity'],
										      (flow.modelList[l].nElementBoundaryQuadraturePoints_global,
										       flow.modelList[l].nSpace_global))])
		freeSurface.modelList[l].updateCoefficients()
	    timeValues.append(tn)
	    for l in range(nLevels):
		flow.modelList[l].trialSpace.writeFunctionEnsight(flow.modelList[l].u,'flow'+`l`,append=True)
		freeSurface.modelList[l].trialSpace.writeFunctionEnsight(freeSurface.modelList[l].u,'freeSurface'+`l`,append=True)
        for l in range(nLevels):
            matlab.write("phi_%s = [" % l)
            for nN in range(flow.modelList[l].mesh.nNodes_global-1):
                matlab.write(`flow.modelList[l].u.dof[nN]`+";")
            matlab.write(`flow.modelList[l].u.dof[nN+1]`+"]; \n")
            matlab.write("figure; title(\'potential, level = %s\'); patch(\'Vertices\',vertices_%s,\'Faces\',faces_%s,\'FaceVertexCData\',phi_%s,\'FaceColor\',\'interp\'); colorbar; \n" % (l,l,l,l))
            matlab.flush()
            matlab.write("xebq_%s = [" % l)
            maxVel = 0.0
            for ebN in range(flow.modelList[l].mesh.nElementBoundaries_global):
                for k in range(flow.modelList[l].nElementBoundaryQuadraturePoints_elementBoundary):
                    maxVel = max(maxVel,math.fabs(math.sqrt(flow.modelList[l].ebq_global['velocity'][ebN,k,0]**2 +
                                                            flow.modelList[l].ebq_global['velocity'][ebN,k,1]**2)))
                    matlab.write(`flow.modelList[l].ebq_global['x'][ebN,k,0]`+" ")
            matlab.write(";")
            scale = (0.25*flow.modelList[l].mesh.h)/maxVel
            print "scale",scale,"maxVel",maxVel
            for ebN in range(flow.modelList[l].mesh.nElementBoundaries_global):
                for k in range(flow.modelList[l].nElementBoundaryQuadraturePoints_elementBoundary):
                    matlab.write(`flow.modelList[l].ebq_global['x'][ebN,k,0]+flow.modelList[l].ebq_global['velocity'][ebN,k,0]*scale`+" ")
            matlab.write("]; \n")
            
            matlab.write("yebq_%s = [" % l)
            for ebN in range(flow.modelList[l].mesh.nElementBoundaries_global):
                for k in range(flow.modelList[l].nElementBoundaryQuadraturePoints_elementBoundary):
                    matlab.write(`flow.modelList[l].ebq_global['x'][ebN,k,1]`+" ")
            matlab.write(";")
            for ebN in range(flow.modelList[l].mesh.nElementBoundaries_global):
                for k in range(flow.modelList[l].nElementBoundaryQuadraturePoints_elementBoundary):
                    matlab.write(`flow.modelList[l].ebq_global['x'][ebN,k,1]+flow.modelList[l].ebq_global['velocity'][ebN,k,1]*scale`+" ")
            matlab.write("]; \n")
            matlab.write("zebq_%s = [" % l)
            for ebN in range(flow.modelList[l].mesh.nElementBoundaries_global):
                for k in range(flow.modelList[l].nElementBoundaryQuadraturePoints_elementBoundary):
                    matlab.write("0.0 ")
            matlab.write(";")
            for ebN in range(flow.modelList[l].mesh.nElementBoundaries_global):
                for k in range(flow.modelList[l].nElementBoundaryQuadraturePoints_elementBoundary):
                    matlab.write("0.0 ")
            matlab.write("]; \n")
            matlab.write("hold on; plot3(xebq_%s,yebq_%s,zebq_%s,'-k'); xlim([-0.1,1.1]); ylim([-0.1,1.1]); hold off; \n" % (l,l,l) )
            matlab.flush()
            matlab.write("u_%s = [" % l)
            for nN in range(freeSurface.modelList[l].mesh.nNodes_global-1):
#                matlab.write(`freeSurface.modelList[l].u.dof[nN]`+";")
                if freeSurface.modelList[l].u.dof[nN] > 0:
                    matlab.write(`0`+";")
                else:
                    matlab.write(`1`+";")
            if freeSurface.modelList[l].u.dof[nN+1] > 0:
                matlab.write(`0`+"]; \n")
            else:
                matlab.write(`1`+"]; \n")
#            matlab.write(`freeSurface.modelList[l].u.dof[nN+1]`+"]; \n")
            matlab.write("figure; title(\'free surface, level = %s\'); patch(\'Vertices\',vertices_%s,\'Faces\',faces_%s,\'FaceVertexCData\',u_%s,\'FaceColor\',\'interp\')\n" % (l,l,l,l))
            matlab.flush()
	print flowNonlinearSolver.info()
        print flowLinearSolver.info()
	print freeSurfaceNonlinearSolver.info()
        print freeSurfaceLinearSolver.info()
        raw_input('Please press return to continue... \n')
        matlab.close()
        #print matlabout.read()
        #print matlaberr.read()
	for l in range(nLevels):
	    flow.modelList[l].trialSpace.endTimeSeriesEnsight(timeValues,'flow'+`l`,'flow'+`l`)
	    freeSurface.modelList[l].trialSpace.endTimeSeriesEnsight(timeValues,'freeSurface'+`l`,'freeSurface'+`l`)
	    signedDistance.modelList[l].trialSpace.endTimeSeriesEnsight(redistanceTimeValues,'signedDistance'+`l`,'signedDistance'+`l`)
        os.chdir('..')

## @}

if __name__ == '__main__':
    runTests()
    #import profile
    #import pstats
    #profile.run('runTests()','FreeBoundaryProf')
    #p = pstats.Stats('FreeBoundaryProf')
    #p.sort_stats('cumulative').print_stats(20)
    #p.sort_stats('time').print_stats(20)

import proteus
import proteus.mprans.RANS2P as RANS2P
from proteus.mprans.cBEAMS import *
import numpy as np
import beamFEM
from ArchiveBeams import *
from mpi4py import MPI

class Coefficients(proteus.mprans.RANS2P.Coefficients):
    def __init__(self,
                 epsFact=1.5,
                 sigma=72.8,
                 rho_0=998.2,nu_0=1.004e-6,
                 rho_1=1.205,nu_1=1.500e-5,
                 g=[0.0,0.0,-9.8],
                 nd=3,
                 ME_model=0,
                 LS_model=None,
                 VF_model=None,
                 KN_model=None,
                 Closure_0_model=None, #Turbulence closure model
                 Closure_1_model=None, #Second possible Turbulence closure model
                 epsFact_density=None,
                 stokes=False,
                 sd=True,
                 movingDomain=False,
                 useVF=0.0,
                 useRBLES=0.0,
                 useMetrics=0.0,
                 useConstant_he=False,
                 dragAlpha=0.0,
                 dragBeta =0.0,
                 setParamsFunc=None,      #uses setParamsFunc if given
                 dragAlphaTypes=None, #otherwise can use element constant values
                 dragBetaTypes=None, #otherwise can use element constant values
                 porosityTypes=None,
                 killNonlinearDrag=False,
                 waveFlag=None,
                 waveHeight=0.01,
                 waveCelerity=1.0,
                 waveFrequency=1.0,
                 waveNumber=2.0,
                 waterDepth=0.5,
                 Omega_s=[[0.45,0.55],[0.2,0.4],[0.0,1.0]],
                 epsFact_source=1.,
                 epsFact_solid=None,
                 eb_adjoint_sigma=1.0,
                 eb_penalty_constant=10.0,
                 forceStrongDirichlet=False,
                 turbulenceClosureModel=0, #0=No Model, 1=Smagorinksy, 2=Dynamic Smagorinsky, 3=K-Epsilon, 4=K-Omega
                 smagorinskyConstant=0.1,
                 barycenters=None,
                 beamLocation = [],
                 beamLength=[],
                 beamRadius=[],
                 EI = [],
                 GJ = [],
                 nBeamElements=4,
                 beam_quadOrder=3,
                 beamFilename = "Beams",
                 beam_useSparse = False,
                 beam_Cd = 1.2,
                 beam_nlTol= 1.0e-5,
                 beamRigid = True,
                 yVertical = False):

        RANS2P.Coefficients.__init__(self,
                 epsFact,
                 sigma,
                 rho_0,nu_0,
                 rho_1,nu_1,
                 g,
                 nd,
                 ME_model,
                 LS_model,
                 VF_model,
                 KN_model,
                 Closure_0_model, #Turbulence closure model
                 Closure_1_model, #Second possible Turbulence closure model
                 epsFact_density,
                 stokes,
                 sd,
                 movingDomain,
                 useVF,
                 useRBLES,
                 useMetrics,
                 useConstant_he,
                 dragAlpha,
                 dragBeta,
                 setParamsFunc,      #uses setParamsFunc if given
                 dragAlphaTypes, #otherwise can use element constant values
                 dragBetaTypes, #otherwise can use element constant values
                 porosityTypes,
                 killNonlinearDrag,
                 waveFlag,
                 waveHeight,
                 waveCelerity,
                 waveFrequency,
                 waveNumber,
                 waterDepth,
                 Omega_s,
                 epsFact_source,
                 epsFact_solid,
                 eb_adjoint_sigma,
                 eb_penalty_constant,
                 forceStrongDirichlet,
                 turbulenceClosureModel, #0=No Model, 1=Smagorinksy, 2=Dynamic Smagorinsky, 3=K-Epsilon, 4=K-Omega
                 smagorinskyConstant,
                 barycenters)
        self.beamFilename = beamFilename
        self.beamLength=beamLength
        self.beamRadius=np.array(beamRadius)
        self.beamLocation=beamLocation
        self.beamVolume=math.pi*self.beamRadius*self.beamRadius*self.beamLength
        self.nBeamElements=nBeamElements
        self.EI = EI
        self.GJ = GJ
        self.beam_quadOrder=beam_quadOrder
        self.beam_nlTol=beam_nlTol
        self.beam_useSparse=beam_useSparse
        self.beam_Cd = beam_Cd
        self.beam_nlTol = beam_nlTol
        self.beamRigid= beamRigid
        self.beamIsLocal = np.zeros((len(self.beamLocation),1), dtype=np.bool)
        self.yVertical = yVertical

    def attachModels(self,modelList):
        RANS2P.Coefficients.attachModels(self,modelList)
        self.initializeBeams()
        self.netBeamDrag = np.array([0.0])

    def initializeElementQuadrature(self,t,cq):
        RANS2P.Coefficients.initializeElementQuadrature(self,t,cq)

        for i in range(len(self.beamIsLocal)):
            n1 = np.less(cq['x'][:,0]-self.beamLocation[i][0], self.beamLength[i])
            n1 = np.logical_and(n1, np.greater(cq['x'][:,0]-self.beamLocation[i][0], -self.beamLength[i]))
            if self.yVertical:
                 n1 = np.logical_and(n1,np.less(cq['x'][:,2]-self.beamLocation[i][1], self.beamLength[i]))
                 n1 = np.logical_and(n1,np.greater(cq['x'][:,2]-self.beamLocation[i][1], -self.beamLength[i]))
                 n1 = np.logical_and(n1,np.less(cq['x'][:,1], self.beamLength[i]))
            else:
                n1 = np.logical_and(n1,np.less(cq['x'][:,1]-self.beamLocation[i][1], self.beamLength[i]))
                n1 = np.logical_and(n1,np.greater(cq['x'][:,1]-self.beamLocation[i][1], -self.beamLength[i]))
                n1 = np.logical_and(n1,np.less(cq['x'][:,2], self.beamLength[i]))
                                     
            if n1.any():
                self.beamIsLocal[i] = True
            

    def postStep(self,t,firstStep=False):
        
        self.beamDrag=np.array([0.0,0.0,0.0])
        if not self.beamRigid:
            self.updateBeams(t)
        else:
            self.updateBeamLoad()
        if self.comm.isMaster():
            print "wettedAreas"
            print self.wettedAreas[:]
            print "Forces_p"
            print self.netForces_p[:,:]
            print "Forces_v"
            print self.netForces_v[:,:]
            #self.wettedAreaHistory.write("%21.16e\n" % (self.wettedAreas[-1],))
            #self.forceHistory_p.write("%21.16e %21.16e %21.16e\n" %tuple(self.netForces_p[-1,:]))
            #self.forceHistory_p.flush()
            #self.forceHistory_v.write("%21.16e %21.16e %21.16e\n" %tuple(self.netForces_v[-1,:]))
            #self.forceHistory_v.flush()
            #self.momentHistory.write("%21.15e %21.16e %21.16e\n" % tuple(self.netMoments[-1,:]))
            #self.momentHistory.flush()
            #self.forceHistory_drag.write("%21.16e %21.16e %21.16e  %21.16e\n" %tuple([t,self.beamDrag[0], self.beamDrag[1], self.beamDrag[2]]))
            #self.forceHistory_drag.flush()
            self.vel_avg = self.vel_avg/self.rho_0
            self.velocityAverage.write("%21.16e %21.16e %21.16e  %21.16e\n" %tuple([t,self.vel_avg[0], self.vel_avg[1], self.vel_avg[2]]))
            self.velocityAverage.flush()
            self.beamDragHistory.write("%21.16e %21.16e\n" %tuple([t,self.netBeamDrag[0]]))
            self.beamDragHistory.flush()
            
            if self.beamRigid==False:
                self.deflectionHistory.write("%21.16e %21.16e %21.16e %21.16e\n" %tuple([t, self.avgHeight, self.avgDeflection, self.avgAngle]))
                self.deflectionHistory.flush()



    def updateBeamLoad(self):
        from proteus.flcbdfWrappers import globalSum

        for I in range(self.nBeams):
            if I%self.comm.size() == self.comm.rank():
                if self.nd==3:
                    self.Beam_Solver[I].updateLoads(self.q1[I,:],self.q2[I,:], self.q3[I,:])
                elif self.nd==2:
                    self.Beam_Solver[I].updateLoads(self.q3[I,:],self.q2[I,:], self.q1[I,:])
                self.Beam_Solver[I].updateQs(endLoad=[0.0,0.0,0.0],scale=1.0)
                if self.centerBox[I]==True:
                    if self.nd==3:
                        self.beamDrag[0]+= self.Beam_Solver[I].Q1[0]
                        self.beamDrag[1]+= self.Beam_Solver[I].Q2[0]
                        self.beamDrag[2]+= self.Beam_Solver[I].Q3[0]
                    elif self.nd==2:
                        self.beamDrag[0]+= self.Beam_Solver[I].Q3[0]
                        self.beamDrag[1]+= self.Beam_Solver[I].Q2[0]
                        self.beamDrag[2]+= self.Beam_Solver[I].Q1[0]
        #for i in range(3):
        #    self.beamDrag.flat[i] = globalSum(self.beamDrag.flat[i])
        cbeamDrag = np.copy(self.beamDrag)
        self.comm2.Allreduce([self.beamDrag,MPI.DOUBLE], [cbeamDrag, MPI.DOUBLE], op=MPI.SUM)
        self.beamDrag = cbeamDrag
                    
    def updateBeams(self,t):
        from proteus.flcbdfWrappers import globalSum
        self.beamDrag = np.array([0.0,0.0,0.0])
        loadSteps = 20
        xv_hold=np.copy(self.xv)
        yv_hold=np.copy(self.yv)
        zv_hold=np.copy(self.zv)
        xq_hold=np.copy(self.xq)
        yq_hold=np.copy(self.yq)
        zq_hold=np.copy(self.zq)
        for I in range(self.nBeams):
            self.xv[I].flat = 0.0 
            self.yv[I].flat = 0.0
            self.zv[I].flat = 0.0
            self.xq[I].flat = 0.0 
            self.yq[I].flat = 0.0
            self.zq[I].flat = 0.0
            if I%self.comm.size() == self.comm.rank():
                if self.nd==3:
                    self.Beam_Solver[I].updateLoads(self.q1[I,:],self.q2[I,:], self.q3[I,:])
                elif self.nd==2:                   
                    self.Beam_Solver[I].updateLoads(self.q3[I,:],self.q2[I,:], self.q1[I,:])
                for j in range(loadSteps):
                    self.Beam_Solver[I].Phi.flat[:]=0.0
                    self.Beam_Solver[I].updateQs(endLoad=[0.0,0.0,0.0],scale=float(j)/float(loadSteps))
                    counter = 0
                    go = True
                    Phi=np.copy(self.Beam_Solver[I].Phi)
                    while go==True:
                        self.Beam_Solver[I].calculateGradient_Hessian()
                        self.Beam_Solver[I].setBCs()
                        self.Beam_Solver[I].reduceOrder()
                        error=self.Beam_Solver[I].calculateResidual()
                        self.Beam_Solver[I].updateSolution()
                        go=self.Beam_Solver[I].checkConvergence()
                        counter += 1
                        if counter >= 100:
                            go = False
                            self.Beam_Solver[I].Phi=np.copy(Phi)
                if self.nd==3:
                    xv , yv, zv =self.Beam_Solver[I].updateCoords()
                    if np.min(zv) > -1.0e-10:
                        self.xv[I,:] = xv
                        self.yv[I,:] = yv
                        self.zv[I,:] = zv
                        self.xq[I,:].flat[:], self.yq[I,:].flat[:], self.zq[I,:].flat[:] = self.Beam_Solver[I].getCoords_at_Quad()
                    else:
                        self.xv[I,:] = xv_hold[I,:]
                        self.yv[I,:] = yv_hold[I,:]
                        self.zv[I,:] = zv_hold[I,:]
                        self.xq[I,:] = xq_hold[I,:]
                        self.yq[I,:] = yq_hold[I,:]
                        self.zq[I,:] = zq_hold[I,:]
                    self.beamDrag[0]+= self.Beam_Solver[I].Q1[0]
                    self.beamDrag[1]+= self.Beam_Solver[I].Q2[0]
                    self.beamDrag[2]+= self.Beam_Solver[I].Q3[0]
                elif self.nd==2:
                    self.zv[I,:], self.yv[I,:], self.xv[I,:] = self.Beam_Solver[I].updateCoords()
                    self.zq[I,:].flat[:], self.yq[I,:].flat[:], self.xq[I,:].flat[:] = self.Beam_Solver[I].getCoords_at_Quad()
                    self.xq[I,:]+=0.2
                    self.xv[I,:]+=0.2
                    self.beamDrag[0]+= self.Beam_Solver[I].Q3[0]
                    self.beamDrag[1]+= self.Beam_Solver[I].Q2[0]
                    self.beamDrag[2]+= self.Beam_Solver[I].Q1[0]
        #self.beamDrag[0] = globalSum(self.beamDrag[0])
        #self.beamDrag[1] = globalSum(self.beamDrag[1])
        #self.beamDrag[2] = globalSum(self.beamDrag[2])
        cbeamDrag = np.copy(self.beamDrag)
        self.comm2.Allreduce([self.beamDrag,MPI.DOUBLE], [cbeamDrag, MPI.DOUBLE], op=MPI.SUM)
        self.beamDrag = cbeamDrag
        #self.beamDrag = self.comm2.allreduce(self.beamDrag, op=MPI.SUM)
        #for i in range(self.nBeams*(self.nBeamElements+1)):
        #    #self.xv.flat[i] = globalSum(self.xv.flat[i])
        #    #self.yv.flat[i] = globalSum(self.yv.flat[i])
        #    #self.zv.flat[i] = globalSum(self.zv.flat[i])
        cxv = np.copy(self.xv)
        cyv = np.copy(self.yv)
        czv = np.copy(self.zv)
        self.comm2.Allreduce([self.xv, MPI.DOUBLE], [cxv, MPI.DOUBLE], op=MPI.SUM)
        self.comm2.Allreduce([self.yv, MPI.DOUBLE], [cyv, MPI.DOUBLE], op=MPI.SUM)
        self.comm2.Allreduce([self.zv, MPI.DOUBLE], [czv, MPI.DOUBLE], op=MPI.SUM)
        self.xv = cxv
        self.yv = cyv
        self.zv = czv
                
        #self.xv = self.comm2.allreduce(self.xv, op=MPI.SUM)
        #self.yv = self.comm2.allreduce(self.yv, op=MPI.SUM)
        #self.zv = self.comm2.allreduce(self.zv, op=MPI.SUM)

        #for i in range(self.nBeams*self.nBeamElements*self.beam_quadOrder):
        #    #self.xq.flat[i] = globalSum(self.xq.flat[i])
        #    #self.yq.flat[i] = globalSum(self.yq.flat[i])
        #    #self.zq.flat[i] = globalSum(self.zq.flat[i])
        #    self.xq = self.comm2.allreduce(self.xq, op=MPI.SUM)
        #    self.yq = self.comm2.allreduce(self.yq, op=MPI.SUM)
        #    self.zq = self.comm2.allreduce(self.zq, op=MPI.SUM)                
        cxq = np.copy(self.xq)
        cyq = np.copy(self.yq)
        czq = np.copy(self.zq)
        self.comm2.Allreduce([self.xq, MPI.DOUBLE], [cxq, MPI.DOUBLE], op=MPI.SUM)
        self.comm2.Allreduce([self.yq, MPI.DOUBLE], [cyq, MPI.DOUBLE], op=MPI.SUM)
        self.comm2.Allreduce([self.zq, MPI.DOUBLE], [czq, MPI.DOUBLE], op=MPI.SUM)
        self.xq = cxq
        self.yq = cyq
        self.zq = czq
        if self.nBeams > 0 and self.comm.isMaster():
            Archive_time_step(Beam_x=self.xv,
                              Beam_y=self.yv,
                              Beam_z=self.zv,
                              nBeams=self.nBeams,
                              filename=self.beamFilename,
                              t=t,
                              tStep=self.tStep)
            self.tStep += 1
            #self.avgHeight=0.0
            #self.avgDeflection=0.0
            #self.avgAngle=0.0
            #nBeamsLocal=0
            #for I in range(self.nBeams):
            #    if self.centerBox[I]==True:
            #        nBeamsLocal+=1
            #        self.avgHeight+=self.zv[I,-1]
            #        self.avgDeflection += abs(self.xv[I,-1]-self.xv[I,0])
            #        self.avgAngle+= math.degrees(math.atan((self.zv[I,-1]-self.zv[I,-2])/(self.xv[I,-1]-self.xv[I,-2])))#self.Beam_Solver[I].Phi[-3]
            self.avgHeight = np.sum(self.zv[:,-1])/float(self.nBeams)
            self.avgDeflection = np.sum(np.abs(self.xv[:,-1]-self.xv[:,0]))/float(self.nBeams)
            self.avgAngle = np.sum(np.rad2deg(np.arctan((self.zv[:,-1]-self.zv[:,-2])/(self.xv[:,-1]-self.xv[:,-2]))))/float(self.nBeams)
                

                    
    def initializeBeams(self):
        comm = Comm.get()
        self.comm=comm
        self.comm2 = MPI.COMM_WORLD

        if self.comm.isMaster():
            self.netForceHistory = open("netForceHistory.txt","w")
            self.forceHistory_drag = open("forceHistory_drag.txt","w")
            self.velocityAverage=open("velocityAverage.txt","w")
            self.dragCoefficientHistory=open("dragCoefficientHistory.txt", "w")
            self.beamDragHistory=open("dragBeamHistory.txt", "w")
            if self.beamRigid==False:
                self.deflectionHistory=open("deflectionHistory.txt", "w")

        bBox=[.25*.69, .69-.25*.69]
        self.centerBox=[]

        self.nBeams = len(self.beamLocation)
        # eulerian coords for beam at mesh vertices
        self.xv=np.zeros((self.nBeams, self.nBeamElements+1))
        self.yv=np.zeros((self.nBeams, self.nBeamElements+1))
        self.zv=np.zeros((self.nBeams, self.nBeamElements+1))

        # distributed loads for beam at quadrature points #mesh vertices
        self.q1=np.zeros((self.nBeams, self.nBeamElements,self.beam_quadOrder)) #np.zeros((self.nBeams, self.nBeamElements+1))
        self.q2=np.zeros((self.nBeams, self.nBeamElements,self.beam_quadOrder)) #np.zeros((self.nBeams, self.nBeamElements+1))
        self.q3=np.zeros((self.nBeams, self.nBeamElements,self.beam_quadOrder))#np.zeros((self.nBeams, self.nBeamElements+1))

        #element diameters
        self.Beam_h = np.zeros((self.nBeams, self.nBeamElements))
        
        #eulerian coords for beam at quadrataure points
        self.xq=np.zeros((self.nBeams, self.nBeamElements,self.beam_quadOrder))
        self.yq=np.zeros((self.nBeams, self.nBeamElements,self.beam_quadOrder))
        self.zq=np.zeros((self.nBeams, self.nBeamElements,self.beam_quadOrder))

        #integration constant for beams
        self.dV_beam = np.zeros((self.nBeams, self.nBeamElements,self.beam_quadOrder))

        self.Beam_mesh=np.zeros((self.nBeams, self.nBeamElements+1))
        
        self.Beam_Solver=[]
        boxCount=0
        for i in range(self.nBeams):
            if self.beamLocation[i][1] >= bBox[0] and self.beamLocation[i][1] <= bBox[1]:
                self.centerBox.append(True)
                boxCount+=1
            else:
                self.centerBox.append(False)
            self.Beam_Solver.append(beamFEM.FEMTools(L=self.beamLength[i],
                                          nElements = self.nBeamElements,
                                          quadOrder= self.beam_quadOrder,
                                          EI = self.EI[i],
                                          GJ = self.GJ[i],
                                          nlTol = self.beam_nlTol,
                                          useSparse = self.beam_useSparse,
                                          beamLocation=self.beamLocation[i]))
            self.Beam_mesh[i,:], self.Beam_h[i,:] = self.Beam_Solver[i].structuredMesh()
            self.Beam_Solver[i].GaussQuad()
            self.Beam_Solver[i].basisFunctions()
            self.Beam_Solver[i].initializePhi()
            self.Beam_Solver[i].initializeCoords()
            if self.nd==3:
		self.xv[i,:], self.yv[i,:], self.zv[i,:] = self.Beam_Solver[i].updateCoords()
            elif self.nd==2:
		self.zv[i,:], self.yv[i,:], self.xv[i,:] = self.Beam_Solver[i].updateCoords()
                self.xv[i,:]+= 0.2
            #self.xv[i,:], self.yv[i,:], self.zv[i,:] = self.Beam_Solver[i].updateCoords()
            for j in range(self.nBeamElements):
                for k in range(self.beam_quadOrder):
                    self.dV_beam[i,j,k] = 0.5*self.Beam_h[i,j]*self.Beam_Solver[i].w[k]
           
            if self.nd==3:
                self.xq[i,:].flat[:], self.yq[i,:].flat[:], self.zq[i,:].flat[:] = self.Beam_Solver[i].getCoords_at_Quad()
            elif self.nd==2:
                self.zq[i,:].flat[:], self.yq[i,:].flat[:], self.xq[i,:].flat[:] = self.Beam_Solver[i].getCoords_at_Quad()
                self.xq[i,:] += 0.2

            self.tStep=0
            if self.nBeams > 0 and self.comm.isMaster():
                Archive_time_step(Beam_x=self.xv,
                                  Beam_y=self.yv,
                                  Beam_z=self.zv,
                                  nBeams=self.nBeams,
                                  filename=self.beamFilename,
                                  t=0.0,
                                  tStep=self.tStep)
            self.tStep+=1
        print boxCount
   
        
    pass

class LevelModel(proteus.mprans.RANS2P.LevelModel):

    def __init__(self,
                 uDict,
                 phiDict,
                 testSpaceDict,
                 matType,
                 dofBoundaryConditionsDict,
                 dofBoundaryConditionsSetterDict,
                 coefficients,
                 elementQuadrature,
                 elementBoundaryQuadrature,
                 fluxBoundaryConditionsDict=None,
                 advectiveFluxBoundaryConditionsSetterDict=None,
                 diffusiveFluxBoundaryConditionsSetterDictDict=None,
                 stressTraceBoundaryConditionsSetterDictDict=None,
                 stabilization=None,
                 shockCapturing=None,
                 conservativeFluxDict=None,
                 numericalFluxType=None,
                 TimeIntegrationClass=None,
                 massLumping=False,
                 reactionLumping=False,
                 options=None,
                 name='RANS2P',
                 reuse_trial_and_test_quadrature=True,
                 sd = True,
                 movingDomain=False):
        RANS2P.LevelModel.__init__(self,
                                   uDict,
                                   phiDict,
                                   testSpaceDict,
                                   matType,
                                   dofBoundaryConditionsDict,
                                   dofBoundaryConditionsSetterDict,
                                   coefficients,
                                   elementQuadrature,
                                   elementBoundaryQuadrature,
                                   fluxBoundaryConditionsDict,
                                   advectiveFluxBoundaryConditionsSetterDict,
                                   diffusiveFluxBoundaryConditionsSetterDictDict,
                                   stressTraceBoundaryConditionsSetterDictDict,
                                   stabilization,
                                   shockCapturing,
                                   conservativeFluxDict,
                                   numericalFluxType,
                                   TimeIntegrationClass,
                                   massLumping,
                                   reactionLumping,
                                   options,
                                   name,
                                   reuse_trial_and_test_quadrature,
                                   sd,
                                   movingDomain)
        compKernelFlag = 0
        self.beams = cBEAMS_base(self.nSpace_global,
                                       self.nQuadraturePoints_element,
                                       self.u[0].femSpace.elementMaps.localFunctionSpace.dim,
                                       self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.dim,
                                       self.testSpace[0].referenceFiniteElement.localFunctionSpace.dim,
                                       self.nElementBoundaryQuadraturePoints_elementBoundary,
                                       compKernelFlag)
    

                 
    
    def calculateAuxiliaryQuantitiesAfterStep(self):
        RANS2P.LevelModel.calculateAuxiliaryQuantitiesAfterStep(self)
        self.beamStep()

    def beamStep(self):
        self.coefficients.vel_avg=np.array([0.0,0.0,0.0])
        self.coefficients.mom_u_source_aux.flat[:]=0.0
        self.coefficients.mom_v_source_aux.flat[:]=0.0
        self.coefficients.mom_w_source_aux.flat[:]=0.0
        self.q1 = numpy.zeros(self.coefficients.q1.shape,'d')#.flat[:]
        self.q2 = numpy.zeros(self.coefficients.q1.shape,'d')#.flat[:]
        self.q3 = numpy.zeros(self.coefficients.q1.shape,'d')#.flat[:]
        self.coefficients.netBeamDrag.flat[:]=0.0
        if self.coefficients.yVertical:
            self.beams.calculateBeams(
                self.mesh.nElements_global,
                self.coefficients.rho_0,
                self.coefficients.rho_1,
                self.coefficients.q_phi,
                self.q['x'],
                self.q[('velocity',0)],
                self.q['dV'],
                self.coefficients.mom_u_source_aux,
                self.coefficients.mom_v_source_aux,
                self.coefficients.mom_w_source_aux,
                self.coefficients.nBeams,
                self.coefficients.nBeamElements,
                self.coefficients.beam_quadOrder,
                self.coefficients.beam_Cd,
                self.coefficients.beamRadius,
                self.coefficients.xq, #.flat[:],
                self.coefficients.zq, #.flat[:],
                self.coefficients.yq, #.flat[:],
                self.coefficients.Beam_h, #.flat[:],
                self.coefficients.dV_beam, #.flat[:],
                self.q1,
                self.q3,
                self.q2,
                self.coefficients.vel_avg,
                self.coefficients.netBeamDrag,
                self.coefficients.beamIsLocal)
        else:
            self.beams.calculateBeams(
                self.mesh.nElements_global,
                self.coefficients.rho_0,
                self.coefficients.rho_1,
                self.coefficients.q_phi,
                self.q['x'],
                self.q[('velocity',0)],
                self.q['dV'],
                self.coefficients.mom_u_source_aux,
                self.coefficients.mom_v_source_aux,
                self.coefficients.mom_w_source_aux,
                self.coefficients.nBeams,
                self.coefficients.nBeamElements,
                self.coefficients.beam_quadOrder,
                self.coefficients.beam_Cd,
                self.coefficients.beamRadius,
                self.coefficients.xq, #.flat[:],
                self.coefficients.yq, #.flat[:],
                self.coefficients.zq, #.flat[:],
                self.coefficients.Beam_h, #.flat[:],
                self.coefficients.dV_beam, #.flat[:],
                self.q1,
                self.q2,
                self.q3,
                self.coefficients.vel_avg,
                self.coefficients.netBeamDrag,
                self.coefficients.beamIsLocal)
        #import pdb
        #pdb.set_trace()
 
        from proteus.flcbdfWrappers import globalSum
        for i in range(self.coefficients.nBeams):
            for j in range(self.coefficients.nBeamElements):
                for k in range(self.coefficients.beam_quadOrder):
                    self.coefficients.q1[i,j,k] = globalSum(self.q1[i,j,k])#globalSum(self.q1[i*(self.coefficients.nBeamElements*self.coefficients.beam_quadOrder)+j*self.coefficients.beam_quadOrder+k])
                    self.coefficients.q2[i,j,k] = globalSum(self.q2[i,j,k])#globalSum(self.q2[i*(self.coefficients.nBeamElements*self.coefficients.beam_quadOrder)+j*self.coefficients.beam_quadOrder+k])
                    self.coefficients.q3[i,j,k] = globalSum(self.q3[i,j,k])#globalSum(self.q3[i*(self.coefficients.nBeamElements*self.coefficients.beam_quadOrder)+j*self.coefficients.beam_quadOrder+k])
        for i in range(3):
            self.coefficients.vel_avg[i]=globalSum(self.coefficients.vel_avg[i])

        self.coefficients.vel_avg=self.coefficients.vel_avg/0.1472
        self.coefficients.netBeamDrag[0] = globalSum(self.coefficients.netBeamDrag[0])

        

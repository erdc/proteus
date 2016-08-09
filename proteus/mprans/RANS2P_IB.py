import proteus
import proteus.mprans.RANS2P as RANS2P
from proteus.mprans.cRANS2P_IB import *
import numpy as np
import beamFEM
from ArchiveBeams import *

class Coefficients(proteus.mprans.RANS2P.Coefficients):
    def __init__(self,
                 epsFact=1.5,
                 sigma=72.8,
                 rho_0=998.2,nu_0=1.004e-6,
                 rho_1=1.205,nu_1=1.500e-5,
                 g=[0.0,0.0,-9.8],
                 nd=3,
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
                 dragAlpha=0.01,
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
                 beamRigid = True):

        RANS2P.Coefficients.__init__(self,
                 epsFact=1.5,
                 sigma=72.8,
                 rho_0=998.2,nu_0=1.004e-6,
                 rho_1=1.205,nu_1=1.500e-5,
                 g=[0.0,0.0,-9.8],
                 nd=3,
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
                 dragAlpha=0.01,
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
                 barycenters=None)
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

    def attachModels(self,modelList):
        RANS2P.Coefficients.attachModels(self,modelList)
        self.initializeBeams()
        self.netBeamDrag = np.array([0.0])

    def initializeElementQuadrature(self,t,cq):
        RANS2P.Coefficients.initializeElementQuadrature(self,t,cq)
        self.q_dragBeam1 = numpy.zeros(cq[('u',1)].shape,'d')
        self.q_dragBeam2 = numpy.zeros(cq[('u',1)].shape,'d')
        self.q_dragBeam3 = numpy.zeros(cq[('u',1)].shape,'d')

    def initializeGlobalExteriorElementBoundaryQuadrature(self,t,cebqe):
        RANS2P.Coefficients.initializeGlobalExteriorElementBoundaryQuadrature(self,t,cebqe)
        self.ebqe_dragBeam1 = numpy.zeros(cebqe[('u',1)].shape,'d')
        self.ebqe_dragBeam2 = numpy.zeros(cebqe[('u',1)].shape,'d')
        self.ebqe_dragBeam3 = numpy.zeros(cebqe[('u',1)].shape,'d')

    def initializeMesh(self, mesh):
        RANS2P.Coefficients.initializeMesh(self, mesh)
        if self.comm.isMaster():
            self.netForceHistory = open("netForceHistory.txt","w")
            self.forceHistory_drag = open("forceHistory_drag.txt","w")
            self.velocityAverage=open("velocityAverage.txt","w")
            self.dragCoefficientHistory=open("dragCoefficientHistory.txt", "w")
            if self.beamRigid==False:
                self.deflectionHistory=open("deflectionHistory.txt", "w")

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
            #self.velocityAverage.write("%21.16e %21.16e %21.16e  %21.16e\n" %tuple([t,self.vel_avg[0], self.vel_avg[1], self.vel_avg[2]]))
            #self.velocityAverage.flush()
            
            if self.beamRigid==False:
                self.deflectionHistory.write("%21.16e %21.16e %21.16e %21.16e\n" %tuple([t, self.avgHeight, self.avgDeflection, self.avgAngle]))
                self.deflectionHistory.flush()

    def delta_h(self, r):
        if r <= -2.0:
            return 0.0
        elif r <= -1.0:
            return 1.0/8.0*(5.0+2.0*r-(-7.0-12.0*r-4*r*r)**0.5)
        elif r <= 0.0:
            return 1.0/8.0*(3.0+2.0*r+(1.0-4.0*r-4.0*r*r)**0.5)
        elif r <= 1.0:
            return 1.0/8.0*(3.0-2.0*r+(1.0+4.0*r-4.0*r*r)**0.5)
        elif r <= 2.0:
            return 1.0/8.0*(5.0-2.0*r-(-7.0+12.0*r-4.0*r*r)**0.5)
        else:
            return 0.0

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
        for i in range(3):
            self.beamDrag.flat[i] = globalSum(self.beamDrag.flat[i])
 
                                                                                                    
                    
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
        self.beamDrag[0] = globalSum(self.beamDrag[0])
        self.beamDrag[1] = globalSum(self.beamDrag[1])
        self.beamDrag[2] = globalSum(self.beamDrag[2])
                                     
        for i in range(self.nBeams*(self.nBeamElements+1)):
            self.xv.flat[i] = globalSum(self.xv.flat[i])
            self.yv.flat[i] = globalSum(self.yv.flat[i])
            self.zv.flat[i] = globalSum(self.zv.flat[i])
        for i in range(self.nBeams*self.nBeamElements*self.beam_quadOrder):
            self.xq.flat[i] = globalSum(self.xq.flat[i])
            self.yq.flat[i] = globalSum(self.yq.flat[i])
            self.zq.flat[i] = globalSum(self.zq.flat[i])
                
        
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
    # def __init__(self,
    #              uDict,
    #              phiDict,
    #              testSpaceDict,
    #              matType,
    #              dofBoundaryConditionsDict,
    #              dofBoundaryConditionsSetterDict,
    #              coefficients,
    #              elementQuadrature,
    #              elementBoundaryQuadrature,
    #              fluxBoundaryConditionsDict=None,
    #              advectiveFluxBoundaryConditionsSetterDict=None,
    #              diffusiveFluxBoundaryConditionsSetterDictDict=None,
    #              stressTraceBoundaryConditionsSetterDictDict=None,
    #              stabilization=None,
    #              shockCapturing=None,
    #              conservativeFluxDict=None,
    #              numericalFluxType=None,
    #              TimeIntegrationClass=None,
    #              massLumping=False,
    #              reactionLumping=False,
    #              options=None,
    #              name='RANS2P',
    #              reuse_trial_and_test_quadrature=True,
    #              sd = True,
    #              movingDomain=False):
    #      RANS2P.LevelModel.__init__(self,
    #                                       uDict,
    #                                       phiDict,
    #                                       testSpaceDict,
    #                                       matType,
    #                                       dofBoundaryConditionsDict,
    #                                       dofBoundaryConditionsSetterDict,
    #                                       coefficients,
    #                                       elementQuadrature,
    #                                       elementBoundaryQuadrature,
    #                                       fluxBoundaryConditionsDict=None,
    #                                       advectiveFluxBoundaryConditionsSetterDict=None,
    #                                       diffusiveFluxBoundaryConditionsSetterDictDict=None,
    #                                       stressTraceBoundaryConditionsSetterDictDict=None,
    #                                       stabilization=None,
    #                                       shockCapturing=None,
    #                                       conservativeFluxDict=None,
    #                                       numericalFluxType=None,
    #                                       TimeIntegrationClass=None,
    #                                       massLumping=False,
    #                                       reactionLumping=False,
    #                                       options=None,
    #                                       name='RANS2P',
    #                                       reuse_trial_and_test_quadrature=True,
    #                                       sd = True,
    #                                       movingDomain=False)
    #      self.rans2p_ib = cRANS2P_IB_base(self.nSpace_global,
    #                                    self.nQuadraturePoints_element,
    #                                    self.u[0].femSpace.elementMaps.localFunctionSpace.dim,
    #                                    self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.dim,
    #                                    self.testSpace[0].referenceFiniteElement.localFunctionSpace.dim,
    #                                    self.nElementBoundaryQuadraturePoints_elementBoundary,
    #                                    compKernelFlag)

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
        self.eb_adjoint_sigma = coefficients.eb_adjoint_sigma
        useConstant_he=coefficients.useConstant_he#this is a hack to test the effect of using a constant smoothing width
        self.postProcessing = True
        #
        #set the objects describing the method and boundary conditions
        #
        self.movingDomain=coefficients.movingDomain
        self.tLast_mesh=None
        #
        #cek todo clean up these flags in the optimized version
        self.bcsTimeDependent=options.bcsTimeDependent
        self.bcsSet=False
        self.name=name
        self.sd=sd
        self.lowmem=True
        self.timeTerm=True#allow turning off  the  time derivative
        self.testIsTrial=True
        self.phiTrialIsTrial=True            
        self.u = uDict
        self.Hess=False
        if isinstance(self.u[0].femSpace,C0_AffineQuadraticOnSimplexWithNodalBasis):
            self.Hess=True
        self.ua = {}#analytical solutions
        self.phi  = phiDict
        self.dphi={}
        self.matType = matType
        #mwf try to reuse test and trial information across components if spaces are the same
        self.reuse_test_trial_quadrature = reuse_trial_and_test_quadrature#True#False
        if self.reuse_test_trial_quadrature:
            for ci in range(1,coefficients.nc):
                assert self.u[ci].femSpace.__class__.__name__ == self.u[0].femSpace.__class__.__name__, "to reuse_test_trial_quad all femSpaces must be the same!"
        ## Simplicial Mesh
        self.mesh = self.u[0].femSpace.mesh #assume the same mesh for  all components for now
        self.testSpace = testSpaceDict
        self.dirichletConditions = dofBoundaryConditionsDict
        self.dirichletNodeSetList=None #explicit Dirichlet  conditions for now, no Dirichlet BC constraints
        self.coefficients = coefficients
        self.coefficients.initializeMesh(self.mesh)
        self.nc = self.coefficients.nc
        self.stabilization = stabilization
        self.shockCapturing = shockCapturing
        self.conservativeFlux = conservativeFluxDict #no velocity post-processing for now
        self.fluxBoundaryConditions=fluxBoundaryConditionsDict
        self.advectiveFluxBoundaryConditionsSetterDict=advectiveFluxBoundaryConditionsSetterDict
        self.diffusiveFluxBoundaryConditionsSetterDictDict = diffusiveFluxBoundaryConditionsSetterDictDict
        #determine whether  the stabilization term is nonlinear
        self.stabilizationIsNonlinear = False
        #cek come back
        if self.stabilization != None:
            for ci in range(self.nc):
        	if coefficients.mass.has_key(ci):
        	    for flag in coefficients.mass[ci].values():
        		if flag == 'nonlinear':
        		    self.stabilizationIsNonlinear=True
        	if  coefficients.advection.has_key(ci):
        	    for  flag  in coefficients.advection[ci].values():
        		if flag == 'nonlinear':
        		    self.stabilizationIsNonlinear=True
        	if  coefficients.diffusion.has_key(ci):
        	    for diffusionDict in coefficients.diffusion[ci].values():
        		for  flag  in diffusionDict.values():
        		    if flag != 'constant':
        			self.stabilizationIsNonlinear=True
        	if  coefficients.potential.has_key(ci):
        	    for flag in coefficients.potential[ci].values():
        		if  flag == 'nonlinear':
        		    self.stabilizationIsNonlinear=True
        	if coefficients.reaction.has_key(ci):
        	    for flag in coefficients.reaction[ci].values():
        		if  flag == 'nonlinear':
        		    self.stabilizationIsNonlinear=True
        	if coefficients.hamiltonian.has_key(ci):
        	    for flag in coefficients.hamiltonian[ci].values():
        		if  flag == 'nonlinear':
        		    self.stabilizationIsNonlinear=True
        #determine if we need element boundary storage
        self.elementBoundaryIntegrals = {}
        for ci  in range(self.nc):
            self.elementBoundaryIntegrals[ci] = ((self.conservativeFlux != None) or 
                                                 (numericalFluxType != None) or 
                                                 (self.fluxBoundaryConditions[ci] == 'outFlow') or
                                                 (self.fluxBoundaryConditions[ci] == 'mixedFlow') or
                                                 (self.fluxBoundaryConditions[ci] == 'setFlow'))
        #
        #calculate some dimensions
        #
        self.nSpace_global    = self.u[0].femSpace.nSpace_global #assume same space dim for all variables
        self.nDOF_trial_element     = [u_j.femSpace.max_nDOF_element for  u_j in self.u.values()]
        self.nDOF_phi_trial_element     = [phi_k.femSpace.max_nDOF_element for  phi_k in self.phi.values()]
        self.n_phi_ip_element = [phi_k.femSpace.referenceFiniteElement.interpolationConditions.nQuadraturePoints for  phi_k in self.phi.values()]
        self.nDOF_test_element     = [femSpace.max_nDOF_element for femSpace in self.testSpace.values()]
        self.nFreeDOF_global  = [dc.nFreeDOF_global for dc in self.dirichletConditions.values()]
        self.nVDOF_element    = sum(self.nDOF_trial_element)
        self.nFreeVDOF_global = sum(self.nFreeDOF_global) 
        #
        NonlinearEquation.__init__(self,self.nFreeVDOF_global)
        #
        #build the quadrature point dictionaries from the input (this
        #is just for convenience so that the input doesn't have to be
        #complete)
        #
        elementQuadratureDict={}
        elemQuadIsDict = isinstance(elementQuadrature,dict)
        if elemQuadIsDict: #set terms manually
            for I in self.coefficients.elementIntegralKeys:
                if elementQuadrature.has_key(I):
                    elementQuadratureDict[I] = elementQuadrature[I]
                else:
                    elementQuadratureDict[I] = elementQuadrature['default']
        else:
            for I in self.coefficients.elementIntegralKeys:
                elementQuadratureDict[I] = elementQuadrature
        if self.stabilization != None:
            for I in self.coefficients.elementIntegralKeys:
                if elemQuadIsDict:
                    if elementQuadrature.has_key(I):
                        elementQuadratureDict[('stab',)+I[1:]] = elementQuadrature[I]
                    else:
                        elementQuadratureDict[('stab',)+I[1:]] = elementQuadrature['default']
                else:
                    elementQuadratureDict[('stab',)+I[1:]] = elementQuadrature
        if self.shockCapturing != None:
            for ci in self.shockCapturing.components:
                if elemQuadIsDict:
                    if elementQuadrature.has_key(('numDiff',ci,ci)):
                        elementQuadratureDict[('numDiff',ci,ci)] = elementQuadrature[('numDiff',ci,ci)]
                    else:
                        elementQuadratureDict[('numDiff',ci,ci)] = elementQuadrature['default']
                else:
                    elementQuadratureDict[('numDiff',ci,ci)] = elementQuadrature
        if massLumping:
            for ci in self.coefficients.mass.keys():
                elementQuadratureDict[('m',ci)] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global,1)
            for I in self.coefficients.elementIntegralKeys:
                elementQuadratureDict[('stab',)+I[1:]] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global,1)
        if reactionLumping:
            for ci in self.coefficients.mass.keys():
                elementQuadratureDict[('r',ci)] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global,1)
            for I in self.coefficients.elementIntegralKeys:
                elementQuadratureDict[('stab',)+I[1:]] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global,1)
        elementBoundaryQuadratureDict={}
        if isinstance(elementBoundaryQuadrature,dict): #set terms manually
            for I in self.coefficients.elementBoundaryIntegralKeys:
                if elementBoundaryQuadrature.has_key(I):
                    elementBoundaryQuadratureDict[I] = elementBoundaryQuadrature[I]
                else:
                    elementBoundaryQuadratureDict[I] = elementBoundaryQuadrature['default']
        else:
            for I in self.coefficients.elementBoundaryIntegralKeys: 
                elementBoundaryQuadratureDict[I] = elementBoundaryQuadrature
        #
        # find the union of all element quadrature points and
        # build a quadrature rule for each integral that has a
        # weight at each point in the union
        (self.elementQuadraturePoints,self.elementQuadratureWeights,
         self.elementQuadratureRuleIndeces) = Quadrature.buildUnion(elementQuadratureDict)
        self.nQuadraturePoints_element = self.elementQuadraturePoints.shape[0]
        self.nQuadraturePoints_global = self.nQuadraturePoints_element*self.mesh.nElements_global
        #
        #Repeat the same thing for the element boundary quadrature
        #
        (self.elementBoundaryQuadraturePoints,
         self.elementBoundaryQuadratureWeights,
         self.elementBoundaryQuadratureRuleIndeces) = Quadrature.buildUnion(elementBoundaryQuadratureDict)
        self.nElementBoundaryQuadraturePoints_elementBoundary = self.elementBoundaryQuadraturePoints.shape[0]
        self.nElementBoundaryQuadraturePoints_global = (self.mesh.nElements_global*
                                                        self.mesh.nElementBoundaries_element*
                                                        self.nElementBoundaryQuadraturePoints_elementBoundary)
        #
        #simplified allocations for test==trial and also check if space is mixed or not
        #
        self.q={}
        self.ebq={}
        self.ebq_global={}
        self.ebqe={}
        self.phi_ip={}
        #mesh
        self.ebqe['x'] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,3),'d')
        self.ebq_global[('totalFlux',0)] = numpy.zeros((self.mesh.nElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebq_global[('velocityAverage',0)] = numpy.zeros((self.mesh.nElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),'d')
        self.q[('u',1)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('u',2)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('u',3)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('m',1)] = self.q[('u',1)]
        self.q[('m',2)] = self.q[('u',2)]
        self.q[('m',3)] = self.q[('u',3)]
        self.q[('m_last',1)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('m_last',2)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('m_last',3)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('m_tmp',1)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('m_tmp',2)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('m_tmp',3)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('mt',1)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('mt',2)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('mt',3)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('dV_u',1)] = (1.0/self.mesh.nElements_global)*numpy.ones((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('dV_u',2)] = (1.0/self.mesh.nElements_global)*numpy.ones((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('dV_u',3)] = (1.0/self.mesh.nElements_global)*numpy.ones((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('f',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nSpace_global),'d')
        self.q[('velocity',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nSpace_global),'d')
        self.q['velocity_solid'] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nSpace_global),'d')
        self.q['phi_solid'] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q['x'] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,3),'d')
        self.q[('cfl',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('numDiff',1,1)] =  numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('numDiff',2,2)] =  numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('numDiff',3,3)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.ebqe[('u',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('u',1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('u',2)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('u',3)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('advectiveFlux_bc_flag',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'i')
        self.ebqe[('advectiveFlux_bc_flag',1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'i')
        self.ebqe[('advectiveFlux_bc_flag',2)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'i')
        self.ebqe[('advectiveFlux_bc_flag',3)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'i')
        self.ebqe[('diffusiveFlux_bc_flag',1,1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'i')
        self.ebqe[('diffusiveFlux_bc_flag',2,2)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'i')
        self.ebqe[('diffusiveFlux_bc_flag',3,3)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'i')
        self.ebqe[('advectiveFlux_bc',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('advectiveFlux_bc',1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('advectiveFlux_bc',2)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('advectiveFlux_bc',3)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('diffusiveFlux_bc',1,1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe['penalty'] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('diffusiveFlux_bc',2,2)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('diffusiveFlux_bc',3,3)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('velocity',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),'d')
        self.ebqe[('velocity',1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),'d')
        self.ebqe[('velocity',2)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),'d')
        self.ebqe[('velocity',3)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),'d')
        #VRANS start, defaults to RANS 
        self.q[('r',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q['eddy_viscosity'] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        #VRANS end
        #RANS 2eq Models start
        self.q[('grad(u)',1)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nSpace_global),'d')
        self.q[('grad(u)',2)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nSpace_global),'d')
        self.q[('grad(u)',3)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nSpace_global),'d')
        #probably don't need ebqe gradients
        self.ebqe[('grad(u)',1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),'d')
        self.ebqe[('grad(u)',2)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),'d')
        self.ebqe[('grad(u)',3)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),'d')
        #RANS 2eq Models end
        self.points_elementBoundaryQuadrature= set()
        self.scalars_elementBoundaryQuadrature= set([('u',ci) for ci in range(self.nc)])
        self.vectors_elementBoundaryQuadrature= set()
        self.tensors_elementBoundaryQuadrature= set()
        #use post processing tools to get conservative fluxes, None by default
        if self.postProcessing:
            self.q[('v',0)] = numpy.zeros(
                (self.mesh.nElements_global,
                 self.nQuadraturePoints_element,
                 self.nDOF_trial_element[0]),
                'd')
            self.q['J'] = numpy.zeros(
                (self.mesh.nElements_global,
                 self.nQuadraturePoints_element,
                 self.nSpace_global,
                 self.nSpace_global),
                'd')
            self.q['det(J)'] = numpy.zeros(
                (self.mesh.nElements_global,
                 self.nQuadraturePoints_element),
                'd')
            self.q['inverse(J)'] = numpy.zeros(
                (self.mesh.nElements_global,
                 self.nQuadraturePoints_element,
                 self.nSpace_global,
                 self.nSpace_global),
                'd')
            self.ebq[('v',0)] = numpy.zeros(
                (self.mesh.nElements_global,
                 self.mesh.nElementBoundaries_element,
                 self.nElementBoundaryQuadraturePoints_elementBoundary,
                 self.nDOF_trial_element[0]),
                'd')
            self.ebq[('w',0)] = numpy.zeros(
                (self.mesh.nElements_global,
                 self.mesh.nElementBoundaries_element,
                 self.nElementBoundaryQuadraturePoints_elementBoundary,
                 self.nDOF_trial_element[0]),
                'd')
            self.ebq['x'] = numpy.zeros(
                (self.mesh.nElements_global,
                 self.mesh.nElementBoundaries_element,
                 self.nElementBoundaryQuadraturePoints_elementBoundary,
                 3),
                'd')
            self.ebq['hat(x)'] = numpy.zeros(
                (self.mesh.nElements_global,
                 self.mesh.nElementBoundaries_element,
                 self.nElementBoundaryQuadraturePoints_elementBoundary,
                 3),
                'd')
            self.ebq['inverse(J)'] = numpy.zeros(
                (self.mesh.nElements_global,
                 self.mesh.nElementBoundaries_element,
                 self.nElementBoundaryQuadraturePoints_elementBoundary,
                 self.nSpace_global,
                 self.nSpace_global),
                'd')
            self.ebq['g'] = numpy.zeros(
                (self.mesh.nElements_global,
                 self.mesh.nElementBoundaries_element,
                 self.nElementBoundaryQuadraturePoints_elementBoundary,
                 self.nSpace_global-1,
                 self.nSpace_global-1),
                'd')
            self.ebq['sqrt(det(g))'] = numpy.zeros(
                (self.mesh.nElements_global,
                 self.mesh.nElementBoundaries_element,
                 self.nElementBoundaryQuadraturePoints_elementBoundary),
                'd')
            self.ebq['n'] = numpy.zeros(
                (self.mesh.nElements_global,
                 self.mesh.nElementBoundaries_element,
                 self.nElementBoundaryQuadraturePoints_elementBoundary,
                 self.nSpace_global),
                'd')
            self.ebq[('dS_u',0)] = numpy.zeros(
                (self.mesh.nElements_global,
                 self.mesh.nElementBoundaries_element,
                 self.nElementBoundaryQuadraturePoints_elementBoundary),
                'd')
            self.ebqe['dS'] = numpy.zeros(
                (self.mesh.nExteriorElementBoundaries_global,
                 self.nElementBoundaryQuadraturePoints_elementBoundary),
                'd')
            self.ebqe[('dS_u',0)] = self.ebqe['dS']
            self.ebqe['n'] = numpy.zeros(
                (self.mesh.nExteriorElementBoundaries_global,
                 self.nElementBoundaryQuadraturePoints_elementBoundary,
                 self.nSpace_global),
                'd')
            self.ebqe['inverse(J)'] = numpy.zeros(
                (self.mesh.nExteriorElementBoundaries_global,
                 self.nElementBoundaryQuadraturePoints_elementBoundary,
                 self.nSpace_global,
                 self.nSpace_global),
                'd')
            self.ebqe['g'] = numpy.zeros(
                (self.mesh.nExteriorElementBoundaries_global,
                 self.nElementBoundaryQuadraturePoints_elementBoundary,
                 self.nSpace_global-1,
                 self.nSpace_global-1),
                'd')
            self.ebqe['sqrt(det(g))'] = numpy.zeros(
                (self.mesh.nExteriorElementBoundaries_global,
                 self.nElementBoundaryQuadraturePoints_elementBoundary),
                'd')
            self.ebq_global['n'] = numpy.zeros(
                (self.mesh.nElementBoundaries_global,
                 self.nElementBoundaryQuadraturePoints_elementBoundary,
                 self.nSpace_global),
                'd')
            self.ebq_global['x'] = numpy.zeros(
                (self.mesh.nElementBoundaries_global,
                 self.nElementBoundaryQuadraturePoints_elementBoundary,
                 3),
                'd')
        #
        #show quadrature
        #
        log("Dumping quadrature shapes for model %s" % self.name,level=9)
        log("Element quadrature array (q)", level=9)
        for (k,v) in self.q.iteritems(): log(str((k,v.shape)),level=9)
        log("Element boundary quadrature (ebq)",level=9) 
        for (k,v) in self.ebq.iteritems(): log(str((k,v.shape)),level=9)
        log("Global element boundary quadrature (ebq_global)",level=9)
        for (k,v) in self.ebq_global.iteritems(): log(str((k,v.shape)),level=9)
        log("Exterior element boundary quadrature (ebqe)",level=9)
        for (k,v) in self.ebqe.iteritems(): log(str((k,v.shape)),level=9)
        log("Interpolation points for nonlinear diffusion potential (phi_ip)",level=9)
        for (k,v) in self.phi_ip.iteritems(): log(str((k,v.shape)),level=9)
        #
        # allocate residual and Jacobian storage
        #
        #
        # allocate residual and Jacobian storage
        #
        self.elementResidual = [numpy.zeros(
                (self.mesh.nElements_global,
                 self.nDOF_test_element[ci]),
                'd')]
        self.inflowBoundaryBC = {}
        self.inflowBoundaryBC_values = {}
        self.inflowFlux = {}
        for cj in range(self.nc):
            self.inflowBoundaryBC[cj] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,),'i')
            self.inflowBoundaryBC_values[cj] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nDOF_trial_element[cj]),'d')
            self.inflowFlux[cj] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.internalNodes = set(range(self.mesh.nNodes_global))
        #identify the internal nodes this is ought to be in mesh
        ##\todo move this to mesh
        for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
            ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
            eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
            ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
            for i in range(self.mesh.nNodes_element):
                if i != ebN_element:
                    I = self.mesh.elementNodesArray[eN_global,i]
                    self.internalNodes -= set([I])
        self.nNodes_internal = len(self.internalNodes)
        self.internalNodesArray=numpy.zeros((self.nNodes_internal,),'i')
        for nI,n in enumerate(self.internalNodes):
            self.internalNodesArray[nI]=n
        #
        del self.internalNodes
        self.internalNodes = None
        log("Updating local to global mappings",2)
        self.updateLocal2Global()
        log("Building time integration object",2)
        log(memory("inflowBC, internalNodes,updateLocal2Global","OneLevelTransport"),level=4)
        #mwf for interpolating subgrid error for gradients etc
        if self.stabilization and self.stabilization.usesGradientStabilization:
            self.timeIntegration = TimeIntegrationClass(self,integrateInterpolationPoints=True)
        else:
             self.timeIntegration = TimeIntegrationClass(self)
           
        if options != None:
            self.timeIntegration.setFromOptions(options)
        log(memory("TimeIntegration","OneLevelTransport"),level=4)
        log("Calculating numerical quadrature formulas",2)
        self.calculateQuadrature()

        self.setupFieldStrides()

        comm = Comm.get()
        self.comm=comm
        if comm.size() > 1:
            assert numericalFluxType != None and numericalFluxType.useWeakDirichletConditions,"You must use a numerical flux to apply weak boundary conditions for parallel runs"

        log("initalizing numerical flux")
        log(memory("stride+offset","OneLevelTransport"),level=4)
        if numericalFluxType != None:
            if options == None or options.periodicDirichletConditions == None:
                self.numericalFlux = numericalFluxType(self,
                                                       dofBoundaryConditionsSetterDict,
                                                       advectiveFluxBoundaryConditionsSetterDict,
                                                       diffusiveFluxBoundaryConditionsSetterDictDict)
            else:
                self.numericalFlux = numericalFluxType(self,
                                                       dofBoundaryConditionsSetterDict,
                                                       advectiveFluxBoundaryConditionsSetterDict,
                                                       diffusiveFluxBoundaryConditionsSetterDictDict,
                                                       options.periodicDirichletConditions)
        else:
            self.numericalFlux = None
        #set penalty terms
        log("initializing numerical flux penalty")
        self.numericalFlux.penalty_constant = self.coefficients.eb_penalty_constant
        #cek todo move into numerical flux initialization
        if self.ebq_global.has_key('penalty'):
            for ebN in range(self.mesh.nElementBoundaries_global):
                for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                    self.ebq_global['penalty'][ebN,k] = self.numericalFlux.penalty_constant/(self.mesh.elementBoundaryDiametersArray[ebN]**self.numericalFlux.penalty_power)
        #penalty term
        #cek move  to Numerical flux initialization
        if self.ebqe.has_key('penalty'):
            for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
                ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
                for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                    self.ebqe['penalty'][ebNE,k] = self.numericalFlux.penalty_constant/self.mesh.elementBoundaryDiametersArray[ebN]**self.numericalFlux.penalty_power
        log(memory("numericalFlux","OneLevelTransport"),level=4)
        self.elementEffectiveDiametersArray  = self.mesh.elementInnerDiametersArray
        log("setting up post-processing")
        from proteus import PostProcessingTools
        self.velocityPostProcessor = PostProcessingTools.VelocityPostProcessingChooser(self)  
        log(memory("velocity postprocessor","OneLevelTransport"),level=4)
        #helper for writing out data storage
        log("initializing archiver")
        from proteus import Archiver
        self.elementQuadratureDictionaryWriter = Archiver.XdmfWriter()
        self.elementBoundaryQuadratureDictionaryWriter = Archiver.XdmfWriter()
        self.exteriorElementBoundaryQuadratureDictionaryWriter = Archiver.XdmfWriter()
        log("flux bc objects")
        for ci,fbcObject  in self.fluxBoundaryConditionsObjectsDict.iteritems():
            self.ebqe[('advectiveFlux_bc_flag',ci)] = numpy.zeros(self.ebqe[('advectiveFlux_bc',ci)].shape,'i')
            for t,g in fbcObject.advectiveFluxBoundaryConditionsDict.iteritems():
                if self.coefficients.advection.has_key(ci):
                    self.ebqe[('advectiveFlux_bc',ci)][t[0],t[1]] = g(self.ebqe[('x')][t[0],t[1]],self.timeIntegration.t)
                    self.ebqe[('advectiveFlux_bc_flag',ci)][t[0],t[1]] = 1
            for ck,diffusiveFluxBoundaryConditionsDict in fbcObject.diffusiveFluxBoundaryConditionsDictDict.iteritems():
                self.ebqe[('diffusiveFlux_bc_flag',ck,ci)] = numpy.zeros(self.ebqe[('diffusiveFlux_bc',ck,ci)].shape,'i')
                for t,g in diffusiveFluxBoundaryConditionsDict.iteritems():
                    self.ebqe[('diffusiveFlux_bc',ck,ci)][t[0],t[1]] = g(self.ebqe[('x')][t[0],t[1]],self.timeIntegration.t)
                    self.ebqe[('diffusiveFlux_bc_flag',ck,ci)][t[0],t[1]] = 1
        self.numericalFlux.setDirichletValues(self.ebqe)
        if self.movingDomain:
            self.MOVING_DOMAIN=1.0
        else:
            self.MOVING_DOMAIN=0.0
        if self.mesh.nodeVelocityArray==None:
            self.mesh.nodeVelocityArray = numpy.zeros(self.mesh.nodeArray.shape,'d')
        #cek/ido todo replace python loops in modules with optimized code if possible/necessary
        log("dirichlet conditions")
        self.forceStrongConditions=coefficients.forceStrongDirichlet
        self.dirichletConditionsForceDOF = {}
        if self.forceStrongConditions:
            for cj in range(self.nc):
                self.dirichletConditionsForceDOF[cj] = DOFBoundaryConditions(self.u[cj].femSpace,dofBoundaryConditionsSetterDict[cj],weakDirichletConditions=False)
        log("final allocations")
        compKernelFlag = 0
        if self.coefficients.useConstant_he:
            self.elementDiameter = self.mesh.elementDiametersArray.copy()
            self.elementDiameter[:] = max(self.mesh.elementDiametersArray)
        else:
            self.elementDiameter = self.mesh.elementDiametersArray
        if self.nSpace_global == 2:
            import copy
            self.u[3] = copy.deepcopy(self.u[2])
            self.timeIntegration.m_tmp[3] = self.timeIntegration.m_tmp[2].copy()
            self.timeIntegration.beta_bdf[3] = self.timeIntegration.beta_bdf[2].copy()
            self.coefficients.sdInfo[(1,3)] = (numpy.array([0,1,2],dtype='i'),
                                  numpy.array([0,1],dtype='i'))
            self.coefficients.sdInfo[(2,3)] = (numpy.array([0,1,2],dtype='i'),
                                  numpy.array([0,1],dtype='i'))
            self.coefficients.sdInfo[(3,0)] = (numpy.array([0,1,2],dtype='i'),
                                  numpy.array([0,1],dtype='i'))
            self.coefficients.sdInfo[(3,1)] = (numpy.array([0,1,2],dtype='i'),
                                  numpy.array([0,1],dtype='i'))
            self.coefficients.sdInfo[(3,2)] = (numpy.array([0,1,2],dtype='i'),
                                  numpy.array([0,1],dtype='i'))
            self.coefficients.sdInfo[(3,3)] = (numpy.array([0,1,2],dtype='i'),
                                  numpy.array([0,1],dtype='i'))
            self.offset.append(self.offset[2])
            self.stride.append(self.stride[2])
            self.numericalFlux.isDOFBoundary[3] = self.numericalFlux.isDOFBoundary[2].copy()
            self.numericalFlux.ebqe[('u',3)] = self.numericalFlux.ebqe[('u',2)].copy()
            log("calling cRANS2P2D_base ctor")
            self.rans2p = cRANS2P2D_base(self.nSpace_global,
                                         self.nQuadraturePoints_element,
                                         self.u[0].femSpace.elementMaps.localFunctionSpace.dim,
                                         self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.dim,
                                         self.testSpace[0].referenceFiniteElement.localFunctionSpace.dim,
                                         self.nElementBoundaryQuadraturePoints_elementBoundary,
                                         compKernelFlag)
        else:
            log("calling  cRANS2P_base ctor")
            self.rans2p = cRANS2P_IB_base(self.nSpace_global,
                                       self.nQuadraturePoints_element,
                                       self.u[0].femSpace.elementMaps.localFunctionSpace.dim,
                                       self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.dim,
                                       self.testSpace[0].referenceFiniteElement.localFunctionSpace.dim,
                                       self.nElementBoundaryQuadraturePoints_elementBoundary,
                                       compKernelFlag)
        
    def getResidual(self,u,r):
        """
        Calculate the element residuals and add in to the global residual
        """

        #Load the unknowns into the finite element dof
        self.timeIntegration.calculateCoefs()
        self.timeIntegration.calculateU(u)
        self.setUnknowns(self.timeIntegration.u)
        #cek todo put in logic to skip if BC's don't depend on t or u
        #hack
        if self.bcsTimeDependent or not self.bcsSet:
            self.bcsSet=True
            #Dirichlet boundary conditions
            self.numericalFlux.setDirichletValues(self.ebqe)
            #Flux boundary conditions
            for ci,fbcObject  in self.fluxBoundaryConditionsObjectsDict.iteritems():
                for t,g in fbcObject.advectiveFluxBoundaryConditionsDict.iteritems():
                    if self.coefficients.advection.has_key(ci):
                        self.ebqe[('advectiveFlux_bc',ci)][t[0],t[1]] = g(self.ebqe[('x')][t[0],t[1]],self.timeIntegration.t)
                        self.ebqe[('advectiveFlux_bc_flag',ci)][t[0],t[1]] = 1
                for ck,diffusiveFluxBoundaryConditionsDict in fbcObject.diffusiveFluxBoundaryConditionsDictDict.iteritems():
                    for t,g in diffusiveFluxBoundaryConditionsDict.iteritems():
                        self.ebqe[('diffusiveFlux_bc',ck,ci)][t[0],t[1]] = g(self.ebqe[('x')][t[0],t[1]],self.timeIntegration.t)
                        self.ebqe[('diffusiveFlux_bc_flag',ck,ci)][t[0],t[1]] = 1
        r.fill(0.0)
        self.Ct_sge = 4.0
        self.Cd_sge = 36.0
        #TODO how to request problem specific evaluations from coefficient class
        if 'evaluateForcingTerms' in dir(self.coefficients):
            self.coefficients.evaluateForcingTerms(self.timeIntegration.t,self.q,self.mesh,
                                                   self.u[0].femSpace.elementMaps.psi,self.mesh.elementNodesArray)
        self.coefficients.wettedAreas[:]  = 0.0
        self.coefficients.netForces_p[:,:]  = 0.0
        self.coefficients.netForces_v[:,:]  = 0.0
        self.coefficients.netMoments[:,:] = 0.0

        if self.forceStrongConditions:
            for cj in range(len(self.dirichletConditionsForceDOF)):
                for dofN,g in self.dirichletConditionsForceDOF[cj].DOFBoundaryConditionsDict.iteritems():
                    self.u[cj].dof[dofN] = g(self.dirichletConditionsForceDOF[cj].DOFBoundaryPointDict[dofN],self.timeIntegration.t)
        self.r = r
        #self.beamStep()
        # self.coefficients.beamDrag=np.array([0.0,0.0,0.0])
        
        # self.coefficients.updateBeamLoad()
        # print self.coefficients.beamDrag
        # import pdb
        # pdb.set_trace()

        self.rans2p.calculateResidual(#element
            self.u[0].femSpace.elementMaps.psi,
            self.u[0].femSpace.elementMaps.grad_psi,
            self.mesh.nodeArray,
            self.mesh.nodeVelocityArray,
            self.MOVING_DOMAIN,
            self.mesh.elementNodesArray,
            self.elementQuadratureWeights[('u',0)],
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            self.u[1].femSpace.psi,
            self.u[1].femSpace.grad_psi,
            self.u[1].femSpace.psi,
            self.u[1].femSpace.grad_psi,
            #element boundary
            self.u[0].femSpace.elementMaps.psi_trace,
            self.u[0].femSpace.elementMaps.grad_psi_trace,
            self.elementBoundaryQuadratureWeights[('u',0)],
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[1].femSpace.psi_trace,
            self.u[1].femSpace.grad_psi_trace,
            self.u[1].femSpace.psi_trace,
            self.u[1].femSpace.grad_psi_trace,
            self.u[0].femSpace.elementMaps.boundaryNormals,
            self.u[0].femSpace.elementMaps.boundaryJacobians,
            #physics
            self.eb_adjoint_sigma,
            self.elementDiameter,#mesh.elementDiametersArray,
            self.mesh.nodeDiametersArray,
            self.stabilization.hFactor,
            self.mesh.nElements_global,
            self.mesh.nElementBoundaries_owned,
            self.coefficients.useRBLES,
            self.coefficients.useMetrics,
            self.timeIntegration.alpha_bdf,
            self.coefficients.epsFact_density,
            self.coefficients.epsFact,
            self.coefficients.sigma,
            self.coefficients.rho_0,
            self.coefficients.nu_0,
            self.coefficients.rho_1,
            self.coefficients.nu_1,
            self.coefficients.smagorinskyConstant,
            self.coefficients.turbulenceClosureModel,
            self.Ct_sge,
            self.Cd_sge,
            self.shockCapturing.shockCapturingFactor,
            self.numericalFlux.penalty_constant,
            #VRANS start
            self.coefficients.epsFact_solid,
            self.coefficients.q_phi_solid,
            self.coefficients.q_velocity_solid,
            self.coefficients.q_porosity,
            self.coefficients.q_dragAlpha,
            self.coefficients.q_dragBeta,
            self.q[('r',0)],
            self.coefficients.q_turb_var[0],
            self.coefficients.q_turb_var[1],
            self.coefficients.q_turb_var_grad[0],
            #VRANS end
            self.u[0].femSpace.dofMap.l2g,
            self.u[1].femSpace.dofMap.l2g,
            self.u[0].dof,
            self.u[1].dof,
            self.u[2].dof,
            self.u[3].dof,
            self.coefficients.g,
            self.coefficients.useVF,
            self.coefficients.q_vf,
            self.coefficients.q_phi,
            self.coefficients.q_n,
            self.coefficients.q_kappa,            
            self.timeIntegration.m_tmp[1],
            self.timeIntegration.m_tmp[2],
            self.timeIntegration.m_tmp[3],
            self.q[('f',0)],
            self.timeIntegration.beta_bdf[1],
            self.timeIntegration.beta_bdf[2],
            self.timeIntegration.beta_bdf[3],
            self.stabilization.v_last,
            self.q[('cfl',0)],
            self.q[('numDiff',1,1)], 
            self.q[('numDiff',2,2)], 
            self.q[('numDiff',3,3)],
            self.shockCapturing.numDiff_last[1],
            self.shockCapturing.numDiff_last[2],
            self.shockCapturing.numDiff_last[3],
            self.coefficients.sdInfo[(1,1)][0],self.coefficients.sdInfo[(1,1)][1],
            self.coefficients.sdInfo[(1,2)][0],self.coefficients.sdInfo[(1,2)][1],
            self.coefficients.sdInfo[(1,3)][0],self.coefficients.sdInfo[(1,3)][1],
            self.coefficients.sdInfo[(2,2)][0],self.coefficients.sdInfo[(2,2)][1],
            self.coefficients.sdInfo[(2,1)][0],self.coefficients.sdInfo[(2,1)][1],
            self.coefficients.sdInfo[(2,3)][0],self.coefficients.sdInfo[(2,3)][1],
            self.coefficients.sdInfo[(3,3)][0],self.coefficients.sdInfo[(3,3)][1],
            self.coefficients.sdInfo[(3,1)][0],self.coefficients.sdInfo[(3,1)][1],
            self.coefficients.sdInfo[(3,2)][0],self.coefficients.sdInfo[(3,2)][1],
            self.offset[0],self.offset[1],self.offset[2],self.offset[3],
            self.stride[0],self.stride[1],self.stride[2],self.stride[3],
            r,
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray,
            self.coefficients.ebqe_vf,
            self.coefficients.bc_ebqe_vf,
            self.coefficients.ebqe_phi,
            self.coefficients.bc_ebqe_phi,
            self.coefficients.ebqe_n,
            self.coefficients.ebqe_kappa,
            #VRANS start
            self.coefficients.ebqe_porosity,
            self.coefficients.ebqe_turb_var[0],
            self.coefficients.ebqe_turb_var[1],
            #VRANS end
            self.numericalFlux.isDOFBoundary[0],
            self.numericalFlux.isDOFBoundary[1],
            self.numericalFlux.isDOFBoundary[2],
            self.numericalFlux.isDOFBoundary[3],
            self.ebqe[('advectiveFlux_bc_flag',0)],
            self.ebqe[('advectiveFlux_bc_flag',1)],
            self.ebqe[('advectiveFlux_bc_flag',2)],
            self.ebqe[('advectiveFlux_bc_flag',3)],
            self.ebqe[('diffusiveFlux_bc_flag',1,1)],
            self.ebqe[('diffusiveFlux_bc_flag',2,2)],
            self.ebqe[('diffusiveFlux_bc_flag',3,3)],
            self.numericalFlux.ebqe[('u',0)],
            self.ebqe[('advectiveFlux_bc',0)],
            self.ebqe[('advectiveFlux_bc',1)],
            self.ebqe[('advectiveFlux_bc',2)],
            self.ebqe[('advectiveFlux_bc',3)],
            self.numericalFlux.ebqe[('u',1)],
            self.ebqe[('diffusiveFlux_bc',1,1)],
            self.ebqe['penalty'],
            self.numericalFlux.ebqe[('u',2)],
            self.ebqe[('diffusiveFlux_bc',2,2)],
            self.numericalFlux.ebqe[('u',3)],
            self.ebqe[('diffusiveFlux_bc',3,3)],
            self.q['x'],
            self.q[('velocity',0)],
            self.ebqe[('velocity',0)],
            self.ebq_global[('totalFlux',0)],
            self.elementResidual[0],
            self.mesh.elementBoundaryMaterialTypes,
            self.coefficients.barycenters,
            self.coefficients.wettedAreas,
            self.coefficients.netForces_p,
            self.coefficients.netForces_v,
            self.coefficients.netMoments,
            self.coefficients.q_dragBeam1,
            self.coefficients.q_dragBeam2,
            self.coefficients.q_dragBeam3,
            self.coefficients.ebqe_dragBeam1,
            self.coefficients.ebqe_dragBeam2,
            self.coefficients.ebqe_dragBeam3)
	from proteus.flcbdfWrappers import globalSum
        for i in range(self.coefficients.netForces_p.shape[0]):
            self.coefficients.wettedAreas[i] = globalSum(self.coefficients.wettedAreas[i])
            for I in range(3):
                self.coefficients.netForces_p[i,I]  = globalSum(self.coefficients.netForces_p[i,I]) 
                self.coefficients.netForces_v[i,I]  = globalSum(self.coefficients.netForces_v[i,I]) 
                self.coefficients.netMoments[i,I] = globalSum(self.coefficients.netMoments[i,I]) 
	if self.forceStrongConditions:#
	    for cj in range(len(self.dirichletConditionsForceDOF)):#
		for dofN,g in self.dirichletConditionsForceDOF[cj].DOFBoundaryConditionsDict.iteritems():
                     r[self.offset[cj]+self.stride[cj]*dofN] = 0
        cflMax=globalMax(self.q[('cfl',0)].max())*self.timeIntegration.dt
        log("Maximum CFL = " + str(cflMax),level=2)
        if self.stabilization:
            self.stabilization.accumulateSubgridMassHistory(self.q)
        log("Global residual",level=9,data=r)
        #mwf decide if this is reasonable for keeping solver statistics
        self.nonlinear_function_evaluations += 1
    def getJacobian(self,jacobian):
	cfemIntegrals.zeroJacobian_CSR(self.nNonzerosInJacobian,
				       jacobian)
        if self.nSpace_global == 2:
            self.csrRowIndeces[(0,3)]  = self.csrRowIndeces[(0,2)]
            self.csrColumnOffsets[(0,3)] = self.csrColumnOffsets[(0,2)]
            self.csrRowIndeces[(1,3)] = self.csrRowIndeces[(0,2)]
            self.csrColumnOffsets[(1,3)] = self.csrColumnOffsets[(0,2)]
            self.csrRowIndeces[(2,3)] = self.csrRowIndeces[(0,2)]
            self.csrColumnOffsets[(2,3)] = self.csrColumnOffsets[(0,2)]
            self.csrRowIndeces[(3,0)] = self.csrRowIndeces[(2,0)]
            self.csrColumnOffsets[(3,0)] = self.csrColumnOffsets[(2,0)]
            self.csrRowIndeces[(3,1)] = self.csrRowIndeces[(2,0)]
            self.csrColumnOffsets[(3,1)] = self.csrColumnOffsets[(2,0)]
            self.csrRowIndeces[(3,2)] = self.csrRowIndeces[(2,0)]
            self.csrColumnOffsets[(3,2)] = self.csrColumnOffsets[(2,0)]
            self.csrRowIndeces[(3,3)] = self.csrRowIndeces[(2,0)]
            self.csrColumnOffsets[(3,3)] = self.csrColumnOffsets[(2,0)]
            self.csrColumnOffsets_eb[(0,3)] = self.csrColumnOffsets[(0,2)]
            self.csrColumnOffsets_eb[(1,3)] = self.csrColumnOffsets[(0,2)]
            self.csrColumnOffsets_eb[(2,3)] = self.csrColumnOffsets[(0,2)]
            self.csrColumnOffsets_eb[(3,0)] = self.csrColumnOffsets[(0,2)]
            self.csrColumnOffsets_eb[(3,1)] = self.csrColumnOffsets[(0,2)]
            self.csrColumnOffsets_eb[(3,2)] = self.csrColumnOffsets[(0,2)]
            self.csrColumnOffsets_eb[(3,3)] = self.csrColumnOffsets[(0,2)]
  
        self.rans2p.calculateJacobian(#element
            self.u[0].femSpace.elementMaps.psi,
            self.u[0].femSpace.elementMaps.grad_psi,
            self.mesh.nodeArray,
            self.mesh.nodeVelocityArray,
            self.MOVING_DOMAIN,
            self.mesh.elementNodesArray,
            self.elementQuadratureWeights[('u',0)],
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            self.u[1].femSpace.psi,
            self.u[1].femSpace.grad_psi,
            self.u[1].femSpace.psi,
            self.u[1].femSpace.grad_psi,
            #element boundary
            self.u[0].femSpace.elementMaps.psi_trace,
            self.u[0].femSpace.elementMaps.grad_psi_trace,
            self.elementBoundaryQuadratureWeights[('u',0)],
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[1].femSpace.psi_trace,
            self.u[1].femSpace.grad_psi_trace,
            self.u[1].femSpace.psi_trace,
            self.u[1].femSpace.grad_psi_trace,
            self.u[0].femSpace.elementMaps.boundaryNormals,
            self.u[0].femSpace.elementMaps.boundaryJacobians,
            self.eb_adjoint_sigma,
            self.elementDiameter,#mesh.elementDiametersArray,
            self.mesh.nodeDiametersArray,
            self.stabilization.hFactor,
            self.mesh.nElements_global,
            self.coefficients.useRBLES,
            self.coefficients.useMetrics,
            self.timeIntegration.alpha_bdf,
            self.coefficients.epsFact_density,
            self.coefficients.epsFact,
            self.coefficients.sigma,
            self.coefficients.rho_0,
            self.coefficients.nu_0,
            self.coefficients.rho_1,
            self.coefficients.nu_1,
            self.coefficients.smagorinskyConstant,
            self.coefficients.turbulenceClosureModel,
            self.Ct_sge,
            self.Cd_sge,
            self.shockCapturing.shockCapturingFactor,
            self.numericalFlux.penalty_constant,
            #VRANS start
            self.coefficients.epsFact_solid,
            self.coefficients.q_phi_solid,
            self.coefficients.q_velocity_solid,
            self.coefficients.q_porosity,
            self.coefficients.q_dragAlpha,
            self.coefficients.q_dragBeta,
            self.q[('r',0)],
            self.coefficients.q_turb_var[0],
            self.coefficients.q_turb_var[1],
            self.coefficients.q_turb_var_grad[0],
            #VRANS end
            self.u[0].femSpace.dofMap.l2g,
            self.u[1].femSpace.dofMap.l2g,
            self.u[0].dof,
            self.u[1].dof,
            self.u[2].dof,
            self.u[3].dof,
            self.coefficients.g,
            self.coefficients.useVF,
            self.coefficients.q_vf,
            self.coefficients.q_phi,
            self.coefficients.q_n,
            self.coefficients.q_kappa,
            self.timeIntegration.beta_bdf[1],
            self.timeIntegration.beta_bdf[2],
            self.timeIntegration.beta_bdf[3],
            self.stabilization.v_last,
            self.q[('cfl',0)],
            self.shockCapturing.numDiff_last[1],
            self.shockCapturing.numDiff_last[2],
            self.shockCapturing.numDiff_last[3],
            self.coefficients.sdInfo[(1,1)][0],self.coefficients.sdInfo[(1,1)][1],
            self.coefficients.sdInfo[(1,2)][0],self.coefficients.sdInfo[(1,2)][1],
            self.coefficients.sdInfo[(1,3)][0],self.coefficients.sdInfo[(1,3)][1],
            self.coefficients.sdInfo[(2,2)][0],self.coefficients.sdInfo[(2,2)][1],
            self.coefficients.sdInfo[(2,1)][0],self.coefficients.sdInfo[(2,1)][1],
            self.coefficients.sdInfo[(2,3)][0],self.coefficients.sdInfo[(2,3)][1],
            self.coefficients.sdInfo[(3,3)][0],self.coefficients.sdInfo[(3,3)][1],
            self.coefficients.sdInfo[(3,1)][0],self.coefficients.sdInfo[(3,1)][1],
            self.coefficients.sdInfo[(3,2)][0],self.coefficients.sdInfo[(3,2)][1],
            self.csrRowIndeces[(0,0)],self.csrColumnOffsets[(0,0)],
            self.csrRowIndeces[(0,1)],self.csrColumnOffsets[(0,1)],
            self.csrRowIndeces[(0,2)],self.csrColumnOffsets[(0,2)],
            self.csrRowIndeces[(0,3)],self.csrColumnOffsets[(0,3)],
            self.csrRowIndeces[(1,0)],self.csrColumnOffsets[(1,0)],
            self.csrRowIndeces[(1,1)],self.csrColumnOffsets[(1,1)],
            self.csrRowIndeces[(1,2)],self.csrColumnOffsets[(1,2)],
            self.csrRowIndeces[(1,3)],self.csrColumnOffsets[(1,3)],
            self.csrRowIndeces[(2,0)],self.csrColumnOffsets[(2,0)],
            self.csrRowIndeces[(2,1)],self.csrColumnOffsets[(2,1)],
            self.csrRowIndeces[(2,2)],self.csrColumnOffsets[(2,2)],
            self.csrRowIndeces[(2,3)],self.csrColumnOffsets[(2,3)],
            self.csrRowIndeces[(3,0)],self.csrColumnOffsets[(3,0)],
            self.csrRowIndeces[(3,1)],self.csrColumnOffsets[(3,1)],
            self.csrRowIndeces[(3,2)],self.csrColumnOffsets[(3,2)],
            self.csrRowIndeces[(3,3)],self.csrColumnOffsets[(3,3)],
            jacobian,
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray,
            self.coefficients.ebqe_vf,
            self.coefficients.bc_ebqe_vf,
            self.coefficients.ebqe_phi,
            self.coefficients.bc_ebqe_phi,
            self.coefficients.ebqe_n,
            self.coefficients.ebqe_kappa,
            #VRANS start
            self.coefficients.ebqe_porosity,
            self.coefficients.ebqe_turb_var[0],
            self.coefficients.ebqe_turb_var[1],
            #VRANS end
            self.numericalFlux.isDOFBoundary[0],
            self.numericalFlux.isDOFBoundary[1],
            self.numericalFlux.isDOFBoundary[2],
            self.numericalFlux.isDOFBoundary[3],
            self.ebqe[('advectiveFlux_bc_flag',0)],
            self.ebqe[('advectiveFlux_bc_flag',1)],
            self.ebqe[('advectiveFlux_bc_flag',2)],
            self.ebqe[('advectiveFlux_bc_flag',3)],
            self.ebqe[('diffusiveFlux_bc_flag',1,1)],
            self.ebqe[('diffusiveFlux_bc_flag',2,2)],
            self.ebqe[('diffusiveFlux_bc_flag',3,3)],
            self.numericalFlux.ebqe[('u',0)],
            self.ebqe[('advectiveFlux_bc',0)],
            self.ebqe[('advectiveFlux_bc',1)],
            self.ebqe[('advectiveFlux_bc',2)],
            self.ebqe[('advectiveFlux_bc',3)],
            self.numericalFlux.ebqe[('u',1)],
            self.ebqe[('diffusiveFlux_bc',1,1)],
            self.ebqe['penalty'],
            self.numericalFlux.ebqe[('u',2)],
            self.ebqe[('diffusiveFlux_bc',2,2)],
            self.numericalFlux.ebqe[('u',3)],
            self.ebqe[('diffusiveFlux_bc',3,3)],
            self.csrColumnOffsets_eb[(0,0)],
            self.csrColumnOffsets_eb[(0,1)],
            self.csrColumnOffsets_eb[(0,2)],
            self.csrColumnOffsets_eb[(0,3)],
            self.csrColumnOffsets_eb[(1,0)],
            self.csrColumnOffsets_eb[(1,1)],
            self.csrColumnOffsets_eb[(1,2)],
            self.csrColumnOffsets_eb[(1,3)],
            self.csrColumnOffsets_eb[(2,0)],
            self.csrColumnOffsets_eb[(2,1)],
            self.csrColumnOffsets_eb[(2,2)],
            self.csrColumnOffsets_eb[(2,3)],
            self.csrColumnOffsets_eb[(3,0)],
            self.csrColumnOffsets_eb[(3,1)],
            self.csrColumnOffsets_eb[(3,2)],
            self.csrColumnOffsets_eb[(3,3)],
            self.coefficients.q_dragBeam1,
            self.coefficients.q_dragBeam2,
            self.coefficients.q_dragBeam3,
            self.coefficients.ebqe_dragBeam1,
            self.coefficients.ebqe_dragBeam2,
            self.coefficients.ebqe_dragBeam3)

        if not self.forceStrongConditions and max(numpy.linalg.norm(self.u[1].dof,numpy.inf),numpy.linalg.norm(self.u[2].dof,numpy.inf),numpy.linalg.norm(self.u[3].dof,numpy.inf)) < 1.0e-8:
            self.pp_hasConstantNullSpace=True
        else:
            self.pp_hasConstantNullSpace=False

        #Load the Dirichlet conditions directly into residual
        if self.forceStrongConditions:
            scaling = 1.0#probably want to add some scaling to match non-dirichlet diagonals in linear system 
            for cj in range(self.nc):
                for dofN in self.dirichletConditionsForceDOF[cj].DOFBoundaryConditionsDict.keys():
                    global_dofN = self.offset[cj]+self.stride[cj]*dofN
                    for i in range(self.rowptr[global_dofN],self.rowptr[global_dofN+1]):
                        if (self.colind[i] == global_dofN):
                            #print "RBLES forcing residual cj = %s dofN= %s global_dofN= %s was self.nzval[i]= %s now =%s " % (cj,dofN,global_dofN,self.nzval[i],scaling)
                            self.nzval[i] = scaling
                        else:
                            self.nzval[i] = 0.0
                            #print "RBLES zeroing residual cj = %s dofN= %s global_dofN= %s " % (cj,dofN,global_dofN)
        log("Jacobian ",level=10,data=jacobian)
        #mwf decide if this is reasonable for solver statistics
        self.nonlinear_function_jacobian_evaluations += 1
        return jacobian

    def calculateAuxiliaryQuantitiesAfterStep(self):
        #RANS2P.LevelModel.calculateAuxiliaryQuantitiesAfterStep(self)
        self.beamStep()

    def beamStep(self):
        self.coefficients.vel_avg=np.array([0.0,0.0,0.0])
        self.coefficients.q_dragBeam1.flat[:]=0.0
        self.coefficients.q_dragBeam2.flat[:]=0.0
        self.coefficients.q_dragBeam3.flat[:]=0.0
        self.coefficients.ebqe_dragBeam1.flat[:]=0.0
        self.coefficients.ebqe_dragBeam2.flat[:]=0.0
        self.coefficients.ebqe_dragBeam3.flat[:]=0.0
        self.q1 = numpy.zeros(self.coefficients.q1.shape,'d')#.flat[:]
        self.q2 = numpy.zeros(self.coefficients.q1.shape,'d')#.flat[:]
        self.q3 = numpy.zeros(self.coefficients.q1.shape,'d')#.flat[:]
        self.coefficients.netBeamDrag.flat[:]=0.0
        self.rans2p.calculateBeams(#element
            self.u[0].femSpace.elementMaps.psi,
            self.u[0].femSpace.elementMaps.grad_psi,
            self.mesh.nodeArray,
            self.mesh.nodeVelocityArray,
            self.MOVING_DOMAIN,
            self.mesh.elementNodesArray,
            self.elementQuadratureWeights[('u',0)],
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            self.u[1].femSpace.psi,
            self.u[1].femSpace.grad_psi,
            self.u[1].femSpace.psi,
            self.u[1].femSpace.grad_psi,
            #element boundary
            self.u[0].femSpace.elementMaps.psi_trace,
            self.u[0].femSpace.elementMaps.grad_psi_trace,
            self.elementBoundaryQuadratureWeights[('u',0)],
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[1].femSpace.psi_trace,
            self.u[1].femSpace.grad_psi_trace,
            self.u[1].femSpace.psi_trace,
            self.u[1].femSpace.grad_psi_trace,
            self.u[0].femSpace.elementMaps.boundaryNormals,
            self.u[0].femSpace.elementMaps.boundaryJacobians,
            #physics
            self.eb_adjoint_sigma,
            self.elementDiameter,#mesh.elementDiametersArray,
            self.mesh.nodeDiametersArray,
            self.stabilization.hFactor,
            self.mesh.nElements_global,
            self.mesh.nElementBoundaries_owned,
            self.coefficients.useRBLES,
            self.coefficients.useMetrics,
            self.timeIntegration.alpha_bdf,
            self.coefficients.epsFact_density,
            self.coefficients.epsFact,
            self.coefficients.sigma,
            self.coefficients.rho_0,
            self.coefficients.nu_0,
            self.coefficients.rho_1,
            self.coefficients.nu_1,
            self.coefficients.smagorinskyConstant,
            self.coefficients.turbulenceClosureModel,
            self.Ct_sge,
            self.Cd_sge,
            self.shockCapturing.shockCapturingFactor,
            self.numericalFlux.penalty_constant,
            #VRANS start
            self.coefficients.epsFact_solid,
            self.coefficients.q_phi_solid,
            self.coefficients.q_velocity_solid,
            self.coefficients.q_porosity,
            self.coefficients.q_dragAlpha,
            self.coefficients.q_dragBeta,
            self.q[('r',0)],
            self.coefficients.q_turb_var[0],
            self.coefficients.q_turb_var[1],
            self.coefficients.q_turb_var_grad[0],
            #VRANS end
            self.u[0].femSpace.dofMap.l2g,
            self.u[1].femSpace.dofMap.l2g,
            self.u[0].dof,
            self.u[1].dof,
            self.u[2].dof,
            self.u[3].dof,
            self.coefficients.g,
            self.coefficients.useVF,
            self.coefficients.q_vf,
            self.coefficients.q_phi,
            self.coefficients.q_n,
            self.coefficients.q_kappa,            
            self.timeIntegration.m_tmp[1],
            self.timeIntegration.m_tmp[2],
            self.timeIntegration.m_tmp[3],
            self.q[('f',0)],
            self.timeIntegration.beta_bdf[1],
            self.timeIntegration.beta_bdf[2],
            self.timeIntegration.beta_bdf[3],
            self.stabilization.v_last,
            self.q[('cfl',0)],
            self.q[('numDiff',1,1)], 
            self.q[('numDiff',2,2)], 
            self.q[('numDiff',3,3)],
            self.shockCapturing.numDiff_last[1],
            self.shockCapturing.numDiff_last[2],
            self.shockCapturing.numDiff_last[3],
            self.coefficients.sdInfo[(1,1)][0],self.coefficients.sdInfo[(1,1)][1],
            self.coefficients.sdInfo[(1,2)][0],self.coefficients.sdInfo[(1,2)][1],
            self.coefficients.sdInfo[(1,3)][0],self.coefficients.sdInfo[(1,3)][1],
            self.coefficients.sdInfo[(2,2)][0],self.coefficients.sdInfo[(2,2)][1],
            self.coefficients.sdInfo[(2,1)][0],self.coefficients.sdInfo[(2,1)][1],
            self.coefficients.sdInfo[(2,3)][0],self.coefficients.sdInfo[(2,3)][1],
            self.coefficients.sdInfo[(3,3)][0],self.coefficients.sdInfo[(3,3)][1],
            self.coefficients.sdInfo[(3,1)][0],self.coefficients.sdInfo[(3,1)][1],
            self.coefficients.sdInfo[(3,2)][0],self.coefficients.sdInfo[(3,2)][1],
            self.offset[0],self.offset[1],self.offset[2],self.offset[3],
            self.stride[0],self.stride[1],self.stride[2],self.stride[3],
            self.r,
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray,
            self.coefficients.ebqe_vf,
            self.coefficients.bc_ebqe_vf,
            self.coefficients.ebqe_phi,
            self.coefficients.bc_ebqe_phi,
            self.coefficients.ebqe_n,
            self.coefficients.ebqe_kappa,
            #VRANS start
            self.coefficients.ebqe_porosity,
            self.coefficients.ebqe_turb_var[0],
            self.coefficients.ebqe_turb_var[1],
            #VRANS end
            self.numericalFlux.isDOFBoundary[0],
            self.numericalFlux.isDOFBoundary[1],
            self.numericalFlux.isDOFBoundary[2],
            self.numericalFlux.isDOFBoundary[3],
            self.ebqe[('advectiveFlux_bc_flag',0)],
            self.ebqe[('advectiveFlux_bc_flag',1)],
            self.ebqe[('advectiveFlux_bc_flag',2)],
            self.ebqe[('advectiveFlux_bc_flag',3)],
            self.ebqe[('diffusiveFlux_bc_flag',1,1)],
            self.ebqe[('diffusiveFlux_bc_flag',2,2)],
            self.ebqe[('diffusiveFlux_bc_flag',3,3)],
            self.numericalFlux.ebqe[('u',0)],
            self.ebqe[('advectiveFlux_bc',0)],
            self.ebqe[('advectiveFlux_bc',1)],
            self.ebqe[('advectiveFlux_bc',2)],
            self.ebqe[('advectiveFlux_bc',3)],
            self.numericalFlux.ebqe[('u',1)],
            self.ebqe[('diffusiveFlux_bc',1,1)],
            self.ebqe['penalty'],
            self.numericalFlux.ebqe[('u',2)],
            self.ebqe[('diffusiveFlux_bc',2,2)],
            self.numericalFlux.ebqe[('u',3)],
            self.ebqe[('diffusiveFlux_bc',3,3)],
            self.q['x'],
            self.q[('velocity',0)],
            self.ebqe[('velocity',0)],
            self.ebq_global[('totalFlux',0)],
            self.elementResidual[0],
            self.mesh.elementBoundaryMaterialTypes,
            self.coefficients.barycenters,
            self.coefficients.wettedAreas,
            self.coefficients.netForces_p,
            self.coefficients.netForces_v,
            self.coefficients.netMoments,
            self.coefficients.q_dragBeam1,
            self.coefficients.q_dragBeam2,
            self.coefficients.q_dragBeam3,
            self.coefficients.ebqe_dragBeam1,
            self.coefficients.ebqe_dragBeam2,
            self.coefficients.ebqe_dragBeam3,
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
            self.coefficients.netBeamDrag)
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

        

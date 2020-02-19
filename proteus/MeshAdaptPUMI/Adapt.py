from __future__ import division
from builtins import str
from builtins import range
import proteus
import sys
import numpy
from proteus.Profiling import logEvent
from proteus import MeshTools
from proteus import SimTools


def reconstructMesh(domain,mesh):

   if hasattr(domain,"PUMIMesh") and not isinstance(domain,proteus.Domain.PUMIDomain) :

     logEvent("Reconstruct based on Proteus, convert PUMI mesh to Proteus")

     nd = domain.nd
     from scipy import spatial
     meshVertexTree = spatial.cKDTree(mesh.nodeArray)
     meshVertex2Model= [0]*mesh.nNodes_owned

     assert domain.vertices, "model vertices (domain.vertices) were not specified"
     assert domain.vertexFlags, "model classification (domain.vertexFlags) needs to be specified"

     for idx,vertex in enumerate(domain.vertices):
       if(nd==2 and len(vertex) == 2): #there might be a smarter way to do this
         vertex.append(0.0) #need to make a 3D coordinate
       closestVertex = meshVertexTree.query(vertex)
       meshVertex2Model[closestVertex[1]] = 1

     isModelVert = numpy.asarray(meshVertex2Model).astype("i")

     meshBoundaryConnectivity = numpy.zeros((mesh.nExteriorElementBoundaries_global,2+nd),dtype=numpy.int32)
     for elementBdyIdx in range(len(mesh.exteriorElementBoundariesArray)):
       exteriorIdx = mesh.exteriorElementBoundariesArray[elementBdyIdx]
       meshBoundaryConnectivity[elementBdyIdx][0] = mesh.elementBoundaryMaterialTypes[exteriorIdx]
       meshBoundaryConnectivity[elementBdyIdx][1] = mesh.elementBoundaryElementsArray[exteriorIdx][0]
       meshBoundaryConnectivity[elementBdyIdx][2] = mesh.elementBoundaryNodesArray[exteriorIdx][0]
       meshBoundaryConnectivity[elementBdyIdx][3] = mesh.elementBoundaryNodesArray[exteriorIdx][1]
       if(nd==3):
         meshBoundaryConnectivity[elementBdyIdx][4] = mesh.elementBoundaryNodesArray[exteriorIdx][2]

     domain.PUMIMesh.reconstructFromProteus2(mesh.cmesh,isModelVert,meshBoundaryConnectivity)

def PUMI_reallocate(solver,mesh):
   p0 = solver.pList[0]
   n0 = solver.nList[0]
   if solver.TwoPhaseFlow:
       nLevels = p0.myTpFlowProblem.general['nLevels']
       nLayersOfOverlapForParallel = p0.myTpFlowProblem.general['nLayersOfOverlapForParallel']
       parallelPartitioningType = MeshTools.MeshParallelPartitioningTypes.element
       domain = p0.myTpFlowProblem.domain
       domain.MeshOptions.setParallelPartitioningType('element')
   else:
       nLevels = n0.nLevels
       nLayersOfOverlapForParallel = n0.nLayersOfOverlapForParallel
       parallelPartitioningType = n0.parallelPartitioningType
       domain = p0.domain

   logEvent("Generating %i-level mesh from PUMI mesh" % (nLevels,))
   if domain.nd == 3:
     mlMesh = MeshTools.MultilevelTetrahedralMesh(
         0,0,0,skipInit=True,
         nLayersOfOverlap=nLayersOfOverlapForParallel,
         parallelPartitioningType=parallelPartitioningType)
   if domain.nd == 2:
     mlMesh = MeshTools.MultilevelTriangularMesh(
         0,0,0,skipInit=True,
         nLayersOfOverlap=nLayersOfOverlapForParallel,
         parallelPartitioningType=parallelPartitioningType)
   if solver.comm.size()==1:
       mlMesh.generateFromExistingCoarseMesh(
           mesh,nLevels,
           nLayersOfOverlap=nLayersOfOverlapForParallel,
           parallelPartitioningType=parallelPartitioningType)
   else:
       mlMesh.generatePartitionedMeshFromPUMI(
           mesh,nLevels,
           nLayersOfOverlap=nLayersOfOverlapForParallel)
   solver.mlMesh_nList=[]
   for p in solver.pList:
       solver.mlMesh_nList.append(mlMesh)
   if (domain.PUMIMesh.size_field_config() == "isotropicProteus"):
       mlMesh.meshList[0].subdomainMesh.size_field = numpy.ones((mlMesh.meshList[0].subdomainMesh.nNodes_global,1),'d')*1.0e-1
   if (domain.PUMIMesh.size_field_config() == 'anisotropicProteus'):
       mlMesh.meshList[0].subdomainMesh.size_scale = numpy.ones((mlMesh.meshList[0].subdomainMesh.nNodes_global,3),'d')
       mlMesh.meshList[0].subdomainMesh.size_frame = numpy.ones((mlMesh.meshList[0].subdomainMesh.nNodes_global,9),'d')

   #may want to trigger garbage collection here
   solver.modelListOld = solver.modelList
   logEvent("Allocating models on new mesh")
   solver.allocateModels()

def PUMI_recomputeStructures(solver,modelListOld):

   ##This section is to correct any differences in the quadrature point field from the old model

   #Shock capturing lagging needs to be matched

   import copy

   ###Details for solution transfer
   #To get shock capturing lagging correct, the numDiff array needs to be computed correctly with the u^{n} solution.
   #numDiff depends on the PDE residual and can depend on the subgrid error (SGE)
   #the PDE residual depends on the alpha and beta_bdf terms which depend on m_tmp from u^{n-1} as well as VOF or LS fields.
   #getResidual() is used to populate m_tmp, numDiff.
   #The goal is therefore to populate the nodal fields with the old solution, get m_tmp properly and lagged sge properly.
   #Mimic the solver stagger with a new loop to repopulate the nodal fields with u^{n} solution. This is necessary because NS relies on the u^{n-1} field for VOF/LS

   ###This loop stores the current solution (u^n) and loads in the previous timestep solution (u^{n-1}

   for m,mOld in zip(solver.modelList, modelListOld):
       for lm, lu, lr, lmOld in zip(m.levelModelList, m.uList, m.rList, mOld.levelModelList):
           #lm.coefficients.postAdaptStep() #MCorr needs this at the moment
           lm.u_store = lm.u.copy()
           for ci in range(0,lm.coefficients.nc):
               lm.u_store[ci] = lm.u[ci].copy()
           lm.dt_store = copy.deepcopy(lm.timeIntegration.dt)
           for ci in range(0,lm.coefficients.nc):
               lm.u[ci].dof[:] = lm.u[ci].dof_last
           lm.setFreeDOF(lu)

   #All solution fields are now in state u^{n-1} and used to get m_tmp and u_sge
   for m,mOld in zip(solver.modelList, modelListOld):
       for lm, lu, lr, lmOld in zip(m.levelModelList, m.uList, m.rList, mOld.levelModelList):
           lm.getResidual(lu,lr)

           #This gets the subgrid error history correct
           if(modelListOld[0].levelModelList[0].stabilization.lag and ((modelListOld[0].levelModelList[0].stabilization.nSteps - 1) > modelListOld[0].levelModelList[0].stabilization.nStepsToDelay) ):
               solver.modelList[0].levelModelList[0].stabilization.nSteps = solver.modelList[0].levelModelList[0].stabilization.nStepsToDelay
               solver.modelList[0].levelModelList[0].stabilization.updateSubgridErrorHistory()

           #update the eddy-viscosity history
           lm.calculateAuxiliaryQuantitiesAfterStep()


   #shock capturing depends on m_tmp or m_last (if lagged). m_tmp is modified by mass-correction and is pushed into m_last during updateTimeHistory().
   #This leads to a situation where m_last comes from the mass-corrected solutions so post-step is needed to get this behavior.
   #If adapt is called after the first time-step, then skip the post-step for the old solution
   if( (abs(solver.systemStepController.t_system_last - solver.tnList[1])> 1e-12 and  abs(solver.systemStepController.t_system_last - solver.tnList[0])> 1e-12 )
     or solver.opts.hotStart):

       for idx in [3,4]:
           model = solver.modelList[idx]
           solver.preStep(model)
           solver.setWeakDirichletConditions(model)
           model.stepController.setInitialGuess(model.uList,model.rList)
           solverFailed = model.solver.solveMultilevel(uList=model.uList,
                                               rList=model.rList,
                                               par_uList=model.par_uList,
                                               par_rList=model.par_rList)
           solver.postStep(model)

   for m,mOld in zip(solver.modelList, modelListOld):
       for lm, lu, lr, lmOld in zip(m.levelModelList, m.uList, m.rList, mOld.levelModelList):
           lm.timeIntegration.postAdaptUpdate(lmOld.timeIntegration)

           if(hasattr(lm.timeIntegration,"dtLast") and lm.timeIntegration.dtLast is not None):
               lm.timeIntegration.dt = lm.timeIntegration.dtLast

   ###This loop reloads the current solution and the previous solution into proper places
   for m,mOld in zip(solver.modelList, modelListOld):
       for lm, lu, lr, lmOld in zip(m.levelModelList, m.uList, m.rList, mOld.levelModelList):
           for ci in range(0,lm.coefficients.nc):
               lm.u[ci].dof[:] = lm.u_store[ci].dof
               lm.u[ci].dof_last[:] = lm.u_store[ci].dof_last

           lm.setFreeDOF(lu)
           lm.getResidual(lu,lr)

           #This gets the subgrid error history correct
           if(modelListOld[0].levelModelList[0].stabilization.lag and modelListOld[0].levelModelList[0].stabilization.nSteps > modelListOld[0].levelModelList[0].stabilization.nStepsToDelay):
               solver.modelList[0].levelModelList[0].stabilization.nSteps = solver.modelList[0].levelModelList[0].stabilization.nStepsToDelay
               solver.modelList[0].levelModelList[0].stabilization.updateSubgridErrorHistory()
   ###

   ###need to re-distance and mass correct
   if( (abs(solver.systemStepController.t_system_last - solver.tnList[0])> 1e-12) or solver.opts.hotStart  ):
       for idx in [3,4]:
           model = solver.modelList[idx]
           solver.preStep(model)
           solver.setWeakDirichletConditions(model)
           model.stepController.setInitialGuess(model.uList,model.rList)
           solverFailed = model.solver.solveMultilevel(uList=model.uList,
                                                       rList=model.rList,
                                                       par_uList=model.par_uList,
                                                       par_rList=model.par_rList)
           solver.postStep(model)

   for m,mOld in zip(solver.modelList, modelListOld):
       for lm, lu, lr, lmOld in zip(m.levelModelList, m.uList, m.rList, mOld.levelModelList):

         lm.timeIntegration.postAdaptUpdate(lmOld.timeIntegration)
         lm.timeIntegration.dt = lm.dt_store

   ###Shock capturing update happens with the time history update
         if(lmOld.shockCapturing and lmOld.shockCapturing.nStepsToDelay is not None and lmOld.shockCapturing.nSteps > lmOld.shockCapturing.nStepsToDelay):
               lm.shockCapturing.nSteps=lm.shockCapturing.nStepsToDelay
               lm.shockCapturing.updateShockCapturingHistory()

         #update the eddy-viscosity history
         lm.calculateAuxiliaryQuantitiesAfterStep()

def PUMI2Proteus(solver,domain):
   #p0 = solver.pList[0] #This can probably be cleaned up somehow
   #n0 = solver.nList[0]
   p0 = solver.pList[0]
   n0 = solver.nList[0]

   modelListOld = solver.modelListOld
   logEvent("Attach auxiliary variables to new models")
   #(cut and pasted from init, need to cleanup)
   solver.simOutputList = []
   solver.auxiliaryVariables = {}
   solver.newAuxiliaryVariables = {}
   if solver.simFlagsList is not None:
       for p, n, simFlags, model, index in zip(
               solver.pList,
               solver.nList,
               solver.simFlagsList,
               solver.modelList,
               list(range(len(solver.pList)))):
           solver.simOutputList.append(
               SimTools.SimulationProcessor(
                   flags=simFlags,
                   nLevels=n.nLevels,
                   pFile=p,
                   nFile=n,
                   analyticalSolution=p.analyticalSolution))
           model.simTools = solver.simOutputList[-1]

           #Code to refresh attached gauges. The goal is to first purge
           #existing point gauge node associations as that may have changed
           #If there is a line gauge, then all the points must be deleted
           #and remade.
           from collections import OrderedDict
           for av in n.auxiliaryVariables:
             if hasattr(av,'adapted'):
               av.adapted=True
               for point, l_d in av.points.items():
                 if 'nearest_node' in l_d:
                   l_d.pop('nearest_node')
               if(av.isLineGauge or av.isLineIntegralGauge): #if line gauges, need to remove all points
                 av.points = OrderedDict()
               if(av.isGaugeOwner):
                 if(solver.comm.rank()==0 and not av.file.closed):
                   av.file.close()
                 for item in av.pointGaugeVecs:
                   item.destroy()
                 for item in av.pointGaugeMats:
                   item.destroy()
                 for item in av.dofsVecs:
                   item.destroy()

                 av.pointGaugeVecs = []
                 av.pointGaugeMats = []
                 av.dofsVecs = []
                 av.field_ids=[]
                 av.isGaugeOwner=False
           ##reinitialize auxiliaryVariables
           solver.auxiliaryVariables[model.name]= [av.attachModel(model,solver.ar[index]) for av in n.auxiliaryVariables]
   else:
       for p,n,s,model,index in zip(
               solver.pList,
               solver.nList,
               solver.sList,
               solver.modelList,
               list(range(len(solver.pList)))):
           solver.simOutputList.append(SimTools.SimulationProcessor(pFile=p,nFile=n))
           model.simTools = solver.simOutputList[-1]
           model.viewer = Viewers.V_base(p,n,s)
           solver.auxiliaryVariables[model.name]= [av.attachModel(model,solver.ar[index]) for av in n.auxiliaryVariables]
   for avList in list(solver.auxiliaryVariables.values()):
       for av in avList:
           av.attachAuxiliaryVariables(solver.auxiliaryVariables)

   logEvent("Transfering fields from PUMI to Proteus")
   for m in solver.modelList:
     for lm in m.levelModelList:
       coef = lm.coefficients
       if coef.vectorComponents is not None:
         vector=numpy.zeros((lm.mesh.nNodes_global,3),'d')
         domain.PUMIMesh.transferFieldToProteus(
                coef.vectorName.encode('utf-8'), vector)
         for vci in range(len(coef.vectorComponents)):
           lm.u[coef.vectorComponents[vci]].dof[:] = vector[:,vci]
         domain.PUMIMesh.transferFieldToProteus(
                coef.vectorName.encode('utf-8')+b"_old", vector)
         for vci in range(len(coef.vectorComponents)):
           lm.u[coef.vectorComponents[vci]].dof_last[:] = vector[:,vci]
         domain.PUMIMesh.transferFieldToProteus(
                coef.vectorName.encode('utf-8')+b"_old_old", vector)
         for vci in range(len(coef.vectorComponents)):
           lm.u[coef.vectorComponents[vci]].dof_last_last[:] = vector[:,vci]

         del vector
       for ci in range(coef.nc):
         if coef.vectorComponents is None or \
            ci not in coef.vectorComponents:
           scalar=numpy.zeros((lm.mesh.nNodes_global,1),'d')
           domain.PUMIMesh.transferFieldToProteus(
               coef.variableNames[ci].encode('utf-8'), scalar)
           lm.u[ci].dof[:] = scalar[:,0]
           domain.PUMIMesh.transferFieldToProteus(
               coef.variableNames[ci].encode('utf-8')+b"_old", scalar)
           lm.u[ci].dof_last[:] = scalar[:,0]
           domain.PUMIMesh.transferFieldToProteus(
               coef.variableNames[ci].encode('utf-8')+b"_old_old", scalar)
           lm.u[ci].dof_last_last[:] = scalar[:,0]

           del scalar

   logEvent("Attaching models on new mesh to each other")
   for m,ptmp,mOld in zip(solver.modelList, solver.pList, modelListOld):
       for lm, lu, lr, lmOld in zip(m.levelModelList, m.uList, m.rList,mOld.levelModelList):
           #save_dof=[]
           #for ci in range(lm.coefficients.nc):
           #    save_dof.append( lm.u[ci].dof.copy())
           #    lm.u[ci].dof_last = lm.u[ci].dof.copy()
           lm.setFreeDOF(lu)
           #for ci in range(lm.coefficients.nc):
           #    assert((save_dof[ci] == lm.u[ci].dof).all())
           lm.calculateSolutionAtQuadrature()
           lm.timeIntegration.tLast = lmOld.timeIntegration.tLast
           lm.timeIntegration.t = lmOld.timeIntegration.t
           lm.timeIntegration.dt = lmOld.timeIntegration.dt
           assert(lmOld.timeIntegration.tLast == lm.timeIntegration.tLast)
           assert(lmOld.timeIntegration.t == lm.timeIntegration.t)
           assert(lmOld.timeIntegration.dt == lm.timeIntegration.dt)
       m.stepController.dt_model = mOld.stepController.dt_model
       m.stepController.t_model = mOld.stepController.t_model
       m.stepController.t_model_last = mOld.stepController.t_model_last
       m.stepController.substeps = mOld.stepController.substeps

   #if first time-step / initial adapt & not hotstarted
   if(abs(solver.systemStepController.t_system_last - solver.tnList[0])< 1e-12 and not solver.opts.hotStart):
       for index,p,n,m,simOutput in zip(range(len(solver.modelList)),solver.pList,solver.nList,solver.modelList,solver.simOutputList):
           if p.initialConditions is not None:
               logEvent("Setting initial conditions for "+p.name)
               m.setInitialConditions(p.initialConditions,solver.tnList[0])


   #Attach models and do sample residual calculation. The results are usually irrelevant.
   #What's important right now is to re-establish the relationships between data structures.
   #The necessary values will be written in later.
   for m,ptmp,mOld in zip(solver.modelList, solver.pList, modelListOld):
       logEvent("Attaching models to model "+ptmp.name)
       m.attachModels(solver.modelList)
   logEvent("Evaluating residuals and time integration")

   for m,ptmp,mOld in zip(solver.modelList, solver.pList, modelListOld):
       for lm, lu, lr, lmOld in zip(m.levelModelList, m.uList, m.rList, mOld.levelModelList):
           lm.timeTerm=True
           lm.getResidual(lu,lr)
           lm.timeIntegration.initializeTimeHistory(resetFromDOF=True)
           lm.initializeTimeHistory()
           lm.timeIntegration.initializeSpaceHistory()
           lm.getResidual(lu,lr)
           #lm.estimate_mt() #function is empty in all models
       assert(m.stepController.dt_model == mOld.stepController.dt_model)
       assert(m.stepController.t_model == mOld.stepController.t_model)
       assert(m.stepController.t_model_last == mOld.stepController.t_model_last)
       logEvent("Initializing time history for model step controller")
       if(not solver.opts.hotStart):
         m.stepController.initializeTimeHistory()
   #p0.domain.initFlag=True #For next step to take initial conditions from solution, only used on restarts

       #m.stepController.initializeTimeHistory()
   #domain.initFlag=True #For next step to take initial conditions from solution, only used on restarts

   solver.systemStepController.modelList = solver.modelList
   solver.systemStepController.exitModelStep = {}
   solver.systemStepController.controllerList = []
   for model in solver.modelList:
       solver.systemStepController.exitModelStep[model] = False
       if model.levelModelList[-1].timeIntegration.isAdaptive:
           solver.systemStepController.controllerList.append(model)
           solver.systemStepController.maxFailures = model.stepController.maxSolverFailures

   #this sets the timeIntegration time, which might be unnecessary for restart
   if(solver.opts.hotStart):
     solver.systemStepController.stepSequence=[(solver.systemStepController.t_system,m) for m in solver.systemStepController.modelList]
   else:
     solver.systemStepController.choose_dt_system()

   #Don't do anything if this is the initial adapt
   if(abs(solver.systemStepController.t_system_last - solver.tnList[0])> 1e-12  or
     (abs(solver.systemStepController.t_system_last - solver.tnList[0]) < 1e-12 and solver.opts.hotStart)):
       PUMI_recomputeStructures(solver,modelListOld)

       #something different is needed for initial conditions
       #do nothing if archive sequence step because there will be an archive
       #if solver.archiveFlag != ArchiveFlags.EVERY_SEQUENCE_STEP:
       #  solver.tCount+=1
       #  for index,model in enumerate(solver.modelList):
       #    #import pdb; pdb.set_trace()
       #    solver.archiveSolution(
       #      model,
       #      index,
       #      #solver.systemStepController.t_system_last+1.0e-6)
       #      solver.systemStepController.t_system)

       #This logic won't account for if final step doesn't match frequency or if adapt isn't being called
       if((solver.PUMIcheckpointer.frequency>0) and ( (domain.PUMIMesh.nAdapt()!=0) and (domain.PUMIMesh.nAdapt() % solver.PUMIcheckpointer.frequency==0 ) or solver.systemStepController.t_system_last==solver.tnList[-1])):

         solver.PUMIcheckpointer.checkpoint(solver.systemStepController.t_system_last)

   #del modelListOld to free up memory
   del modelListOld
   import gc;
   gc.disable()
   gc.collect()
   solver.comm.barrier()

def PUMI_transferFields(solver):
   p0 = solver.pList[0].ct
   n0 = solver.nList[0].ct

   if solver.TwoPhaseFlow:
       domain = p0.myTpFlowProblem.domain
       rho_0 = p0.myTpFlowProblem.physical_parameters['densityA']
       nu_0 = p0.myTpFlowProblem.physical_parameters['kinematicViscosityA']
       rho_1 = p0.myTpFlowProblem.physical_parameters['densityB']
       nu_1 = p0.myTpFlowProblem.physical_parameters['kinematicViscosityB']
       g = p0.myTpFlowProblem.physical_parameters['gravity']
       epsFact_density = p0.myTpFlowProblem.clsvof_parameters['epsFactHeaviside']
   else:
       domain = p0.domain
       rho_0  = p0.rho_0
       nu_0   = p0.nu_0
       rho_1  = p0.rho_1
       nu_1   = p0.nu_1
       g      = p0.g
       epsFact_density = p0.epsFact_density
   logEvent("Copying coordinates to PUMI")
   domain.PUMIMesh.transferFieldToPUMI(b"coordinates",
       solver.modelList[0].levelModelList[0].mesh.nodeArray)

   #I want to compute the density and viscosity arrays here
   #arrays are length = number of elements and will correspond to density at center of element
   rho_transfer = numpy.zeros((solver.modelList[0].levelModelList[0].mesh.nElements_owned),'d')      
   nu_transfer = numpy.zeros((solver.modelList[0].levelModelList[0].mesh.nElements_owned),'d')      
   #get quadrature points at element centroid and evaluate at shape functions
   from proteus import Quadrature
   transferQpt = Quadrature.SimplexGaussQuadrature(p0.domain.nd,1)
   qpt_centroid = numpy.asarray([transferQpt.points[0]])
   materialSpace = solver.nList[0].femSpaces[0](solver.modelList[0].levelModelList[0].mesh.subdomainMesh,p0.domain.nd)
   materialSpace.getBasisValuesRef(qpt_centroid)
   
   #obtain the level-set or vof value at each element centroid
   #pass through heaviside function to get material property
   from proteus.ctransportCoefficients import smoothedHeaviside
   
   IEN = solver.modelList[2].levelModelList[0].u[0].femSpace.dofMap.l2g
   for (eID, dofs) in enumerate(IEN):
       phi_val = 0.0
       for idx in range(len(dofs)):
           phi_val += materialSpace.psi[0][idx]*solver.modelList[2].levelModelList[0].u[0].dof[dofs[idx]]
       #rho_transfer[eID] = phi_val
       
       #heaviside
       h_phi=0.0;
       for idx in range(len(dofs)):
           h_phi += (materialSpace.psi[0][idx])*(solver.modelList[2].levelModelList[0].mesh.nodeDiametersArray[dofs[idx]]);
       eps_rho = p0.epsFact_density*h_phi
       smoothed_phi_val = smoothedHeaviside(eps_rho,phi_val)

       rho_transfer[eID] = (1.0-smoothed_phi_val)*solver.pList[0].ct.rho_0 + smoothed_phi_val*solver.pList[0].ct.rho_1
       nu_transfer[eID] = (1.0-smoothed_phi_val)*solver.pList[0].ct.nu_0 + smoothed_phi_val*solver.pList[0].ct.nu_1

   solver.modelList[0].levelModelList[0].mesh.elementMaterial = numpy.zeros((solver.modelList[0].levelModelList[0].mesh.nElements_owned),'d') 
   solver.modelList[0].levelModelList[0].mesh.elementMaterial[:] = rho_transfer[:]

   #put the solution field as uList
   #VOF and LS needs to reset the u.dof array for proper transfer
   #but needs to be returned to the original form if not actually adapting....be careful with the following statements, unsure if this doesn't break something else
   import copy
   for m in solver.modelList:
       for lm in m.levelModelList:
           lm.u_store = lm.u.copy()
           for ci in range(0,lm.coefficients.nc):
               lm.u_store[ci] = lm.u[ci].copy()

   solver.modelList[1].levelModelList[0].setUnknowns(solver.modelList[1].uList[0])
   solver.modelList[2].levelModelList[0].setUnknowns(solver.modelList[2].uList[0])

   logEvent("Copying DOF and parameters to PUMI")
   for m in solver.modelList:
     for lm in m.levelModelList:
       coef = lm.coefficients
       if coef.vectorComponents is not None:
         vector=numpy.zeros((lm.mesh.nNodes_global,3),'d')
         for vci in range(len(coef.vectorComponents)):
           vector[:,vci] = lm.u[coef.vectorComponents[vci]].dof[:]

         domain.PUMIMesh.transferFieldToPUMI(
             coef.vectorName.encode('utf-8'), vector)
         #Transfer dof_last
         for vci in range(len(coef.vectorComponents)):
           vector[:,vci] = lm.u[coef.vectorComponents[vci]].dof_last[:]
         domain.PUMIMesh.transferFieldToPUMI(
                coef.vectorName.encode('utf-8')+b"_old", vector)
         #Transfer dof_last_last
         for vci in range(len(coef.vectorComponents)):
           vector[:,vci] = lm.u[coef.vectorComponents[vci]].dof_last_last[:]
         p0.domain.PUMIMesh.transferFieldToPUMI(
                coef.vectorName.encode('utf-8')+b"_old_old", vector)

         del vector
       for ci in range(coef.nc):
         if coef.vectorComponents is None or \
            ci not in coef.vectorComponents:
           scalar=numpy.zeros((lm.mesh.nNodes_global,1),'d')
           scalar[:,0] = lm.u[ci].dof[:]
           domain.PUMIMesh.transferFieldToPUMI(
               coef.variableNames[ci].encode('utf-8'), scalar)

           #Transfer dof_last
           scalar[:,0] = lm.u[ci].dof_last[:]
           domain.PUMIMesh.transferFieldToPUMI(
                coef.variableNames[ci].encode('utf-8')+b"_old", scalar)
           #Transfer dof_last_last
           scalar[:,0] = lm.u[ci].dof_last_last[:]
           p0.domain.PUMIMesh.transferFieldToPUMI(
                coef.variableNames[ci].encode('utf-8')+b"_old_old", scalar)

           del scalar

   scalar=numpy.zeros((lm.mesh.nNodes_global,1),'d')

   del scalar
   #Get Physical Parameters
   #Can we do this in a problem-independent  way?

   rho = numpy.array([rho_0,
                      rho_1])
   nu = numpy.array([nu_0,
                     nu_1])
   g = numpy.asarray(g)

   #This condition is to account for adapting before the simulation started
   if(hasattr(solver,"tn")):
       #deltaT = solver.tn-solver.tn_last
       #is actually the time step for next step, solver.tn and solver.tn_last refer to entries in tnList
       #deltaT = solver.systemStepController.dt_system 
       deltaT=solver.modelList[0].levelModelList[0].timeIntegration.dtLast
       T_current = solver.systemStepController.t_system_last
       deltaT_next = solver.systemStepController.dt_system 
   else:
       deltaT = 0
       deltaT_next = 0.0
       T_current = 0.0

   epsFact = epsFact_density
   #domain.PUMIMesh.transferPropertiesToPUMI(rho,nu,g,deltaT,epsFact)
   domain.PUMIMesh.transferPropertiesToPUMI(rho_transfer,nu_transfer,g,deltaT,deltaT_next,T_current,epsFact)

   del rho, nu, g, epsFact


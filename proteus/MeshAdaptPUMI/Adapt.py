from __future__ import division
from builtins import str
from builtins import range
import proteus
import sys
import numpy
from proteus import Profiling


def reconstructMesh(domain,mesh):

   if hasattr(domain,"PUMIMesh") and not isinstance(domain,proteus.Domain.PUMIDomain) :

     Profiling.logEvent("Reconstruct based on Proteus, convert PUMI mesh to Proteus")

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



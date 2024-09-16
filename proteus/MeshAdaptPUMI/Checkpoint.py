import proteus
import sys
import numpy
from proteus import Profiling

#it should probably be associated with the PUMI domain somehow
#The current implementation assumes we're using NS, VOF, LS, RD, MCorr setup with lagging and Backwards Euler.
#Future work on this module should include creating an abstract class from which variations based on the models and numerical accuracy can be created
#Have the dictionary submodels be labeled by physical model names like "twp_navier_stokes"

class Checkpointer:
    "This class is meant to handle the checkpointing process for adapted meshes. Information that's needed to be loaded into hotstart needs to be output and then read in to be handled for data reconstruction"
    def __init__(self,NSobject,frequency=10):
      self.NSobject = NSobject
      self.counter = 0
      self.frequency = frequency
    def checkpoint(self,hotStartTime):
      self.transferInfo()
      self.saveMesh()
      modelListOld=self.EncodeModel(hotStartTime)

      #pickling is apparently unsafe so we use json to try storing modelListOld
      filename = "checkpointInfo"+str(self.counter)
      f = open(filename, 'w')
      import json
      #json.dump(modelListOld.__dict__,f)
      json.dump(modelListOld,f)
      f.close()
    def transferInfo(self):
      self.NSobject.PUMI_transferFields()
    def saveMesh(self):
      fileName="checkpoint"+str(self.counter)+"_.smb"
      self.NSobject.pList[0].domain.PUMIMesh.writeMesh(fileName)

    def EncodeModel(self,hotStartTime):
      "Grab only necessary components from modelListOld so far consistent only with first-order time integrator" 
      #def __init__(self,modelListOld,hotStartTime):
      modelListOld = self.NSobject.modelListOld
      saveModel = {}
      saveModel["tCount"] = self.NSobject.tCount+1 #+1 just because of how indexing works in h5 file
      saveModel["counter"] = self.counter
      saveModel["numModels"] = len(modelListOld)
      saveModel["hotStartTime"] = hotStartTime
      saveModel["nAdapt"] = self.NSobject.pList[0].domain.PUMIMesh.nAdapt()
      saveModel["checkpoint_status"] = ""

      if(hasattr(self.NSobject,"tn") and (self.NSobject.systemStepController.t_system_last < self.NSobject.tn)):
        saveModel["checkpoint_status"] = "midsystem"
        saveModel["tCount"] = self.NSobject.tCount+2 #don't know how to justify this yet but it's what is needed
      else:
        saveModel["checkpoint_status"] = "endsystem"

      saveModel["systemStepController"]=[]
      controllerAttribute={}
      controllerAttribute["dt_system"]=self.NSobject.systemStepController.dt_system
      controllerAttribute["dt_system_fixed"]=self.NSobject.systemStepController.dt_system_fixed
      controllerAttribute["t_system_last"]=self.NSobject.systemStepController.t_system_last
      controllerAttribute["t_system"]=self.NSobject.systemStepController.t_system
      saveModel["systemStepController"].append(controllerAttribute)

      saveModel["stepController"]=[]
      saveModel["timeIntegration"]=[]
      saveModel["shockCapturing"]=[]
      saveModel["stabilization"]=[]
      for i in range(0,len(modelListOld)):
        #step controller
        subModel={}
        subModel["dt_model"]= modelListOld[i].stepController.dt_model
        subModel["t_model"] = modelListOld[i].stepController.t_model
        subModel["t_model_last"] = modelListOld[i].stepController.t_model_last
        subModel["substeps"]=modelListOld[i].stepController.substeps
        saveModel["stepController"].append(subModel)
        
        #time integration
        subModel={}
        subModel["dt"] = modelListOld[i].levelModelList[0].timeIntegration.dt
        subModel["t"] = modelListOld[i].levelModelList[0].timeIntegration.t
        if(hasattr(modelListOld[i].levelModelList[0].timeIntegration,'dtLast')):
          subModel["dtLast"] = modelListOld[i].levelModelList[0].timeIntegration.dtLast
        else:
          subModel["dtLast"] = None
        saveModel["timeIntegration"].append(subModel)

        #shock capturing
        subModel={}
        if(modelListOld[i].levelModelList[0].shockCapturing is not None):
          subModel["nSteps"]=modelListOld[i].levelModelList[0].shockCapturing.nSteps
          subModel["nStepsToDelay"]= modelListOld[i].levelModelList[0].shockCapturing.nStepsToDelay
        saveModel["shockCapturing"].append(subModel)

      #Assuming the 0th model is RANS2P
      #stabilization
      subModel={}
      subModel["nSteps"]= modelListOld[0].levelModelList[0].stabilization.nSteps
      saveModel["stabilization"].append(subModel)
      return saveModel

    def DecodeModel(self,filename):
      "create a modelListOld that can interact with the post-adapt restart capabilities" 
      f = open(filename, 'r')
      import json
      previousInfo = json.load(f)
      f.close()

      systemStepController = previousInfo["systemStepController"][0]

      self.NSobject.systemStepController.dt_system = systemStepController["dt_system"]  
      self.NSobject.systemStepController.dt_system_fixed = systemStepController["dt_system_fixed"]  
      self.NSobject.systemStepController.t_system_last = systemStepController["t_system_last"]  
      self.NSobject.systemStepController.t_system = systemStepController["t_system"]  

      numModels = previousInfo["numModels"]
      stepController=previousInfo["stepController"]
      timeIntegration=previousInfo["timeIntegration"]
      shockCapturing=previousInfo["shockCapturing"]
      stabilization=previousInfo["stabilization"]
      self.counter = previousInfo["counter"]+1
      
      for i in range(0,numModels):

        self.NSobject.modelList[i].stepController.dt_model = stepController[i]["dt_model"]
        self.NSobject.modelList[i].stepController.t_model = stepController[i]["t_model"]
        self.NSobject.modelList[i].stepController.t_model_last = stepController[i]["t_model_last"]
        self.NSobject.modelList[i].stepController.substeps = stepController[i]["substeps"]

        self.NSobject.modelList[i].levelModelList[0].timeIntegration.dt = timeIntegration[i]["dt"]
        self.NSobject.modelList[i].levelModelList[0].timeIntegration.t = timeIntegration[i]["t"]
        self.NSobject.modelList[i].levelModelList[0].timeIntegration.dtLast = timeIntegration[i]["dtLast"]


        if(self.NSobject.modelList[i].levelModelList[0].shockCapturing is not None):
          self.NSobject.modelList[i].levelModelList[0].shockCapturing.nSteps = shockCapturing[i]["nSteps"]
          self.NSobject.modelList[i].levelModelList[0].shockCapturing.nStepsToDelay = shockCapturing[i]["nStepsToDelay"]

      self.NSobject.modelList[0].levelModelList[0].stabilization.nSteps = stabilization[0]["nSteps"]

      self.NSobject.pList[0].domain.PUMIMesh.set_nAdapt(previousInfo["nAdapt"])

from __future__ import division
from builtins import str
from builtins import range
import proteus
import sys
import numpy
from proteus import Profiling

#it should probably be associated with the PUMI domain somehow
#The current implementation assumes we're using NS, VOF, LS, RD, MCorr setup with lagging and Backwards Euler.
#Future work on this module should include creating an abstract class from which variations based on the models and numerical accuracy can be created
#Also, output items as dictionaries instead of lists, which would be easier to encode and decode
class Checkpointer:
    "This class is meant to handle the checkpointing process for adapted meshes. Information that's needed to be loaded into hotstart needs to be output and then read in to be handled for data reconstruction"
    def __init__(self,NSobject,frequency=10):
      self.A = "Hello"
      self.B = "World"
      self.NSobject = NSobject
      self.counter = 0
      self.frequency = frequency
    def checkpoint(self):
      self.transferInfo()
      self.saveMesh()
      modelListOld=self.EncodeModel(self.NSobject.systemStepController.t_system_last)

      #pickling is apparently unsafe so we use json to try storing modelListOld
      filename = "checkpointInfo"+str(self.counter)
      f = open(filename, 'w')
      import json
      #json.dump(modelListOld.__dict__,f)
      json.dump(modelListOld,f)
      f.close()
      self.counter+=1
    def transferInfo(self):
      self.NSobject.PUMI_transferFields()
    def saveMesh(self):
      fileName="checkpoint"+str(self.counter)+"_.smb"
      self.NSobject.pList[0].domain.PUMIMesh.writeMesh(fileName)
    def DecodeModel(self,filename):
      "create a modelListOld that can interact with the post-adapt restart capabilities" 
      f = open(filename, 'r')
      import json
      previousInfo = json.load(f)
      f.close()

      numModels = previousInfo["numModels"]
      stepController=previousInfo["stepController"]
      timeIntegration=previousInfo["timeIntegration"]
      shockCapturing=previousInfo["shockCapturing"]
      stabilization=previousInfo["stabilization"]
      self.counter = previousInfo["counter"]+1
      
      for i in range(0,numModels):
        self.NSobject.modelList[i].stepController.dt_model = stepController[0][0]
        self.NSobject.modelList[i].stepController.t_model = stepController[0][1]
        self.NSobject.modelList[i].stepController.t_model_last = stepController[0][2]
        self.NSobject.modelList[i].stepController.substeps = stepController[0][3]

        self.NSobject.modelList[i].levelModelList[0].timeIntegration.dt = timeIntegration[0][0]
        self.NSobject.modelList[i].levelModelList[0].timeIntegration.t = timeIntegration[0][1]
        self.NSobject.modelList[i].levelModelList[0].timeIntegration.dt = timeIntegration[0][2]
        self.NSobject.modelList[i].levelModelList[0].timeIntegration.dtLast = timeIntegration[0][3]

        if(self.NSobject.modelList[i].levelModelList[0].shockCapturing is not None):
          self.NSobject.modelList[i].levelModelList[0].shockCapturing.nSteps = shockCapturing[0][0]
          self.NSobject.modelList[i].levelModelList[0].shockCapturing.nStepsToDelay = shockCapturing[0][1]

      self.NSobject.modelList[0].levelModelList[0].stabilization.nSteps = stabilization[0]

#Maybe I can form an abstract class for encoding decoding?
    def EncodeModel(self,hotStartTime):
      "Grab only necessary components from modelListOld so far consistent only with first-order time integrator" 
      #def __init__(self,modelListOld,hotStartTime):
      modelListOld = self.NSobject.modelListOld
      saveModel = {}
      saveModel["counter"] = self.counter
      saveModel["numModels"] = len(modelListOld)
      saveModel["hotStartTime"] = hotStartTime
      saveModel["stepController"]=[]
      saveModel["timeIntegration"]=[]
      saveModel["shockCapturing"]=[]
      saveModel["stabilization"]=[]
      for i in range(0,len(modelListOld)):
        saveModel["stepController"].append([modelListOld[i].stepController.dt_model,modelListOld[i].stepController.t_model,modelListOld[i].stepController.t_model_last, modelListOld[i].stepController.substeps])
        if(hasattr(modelListOld[i].levelModelList[0].timeIntegration,'dtLast')):
          saveModel["timeIntegration"].append([modelListOld[i].levelModelList[0].timeIntegration.dt,modelListOld[i].levelModelList[0].timeIntegration.t,modelListOld[i].levelModelList[0].timeIntegration.dt,modelListOld[i].levelModelList[0].timeIntegration.dtLast])
        else:
          saveModel["timeIntegration"].append([modelListOld[i].levelModelList[0].timeIntegration.dt,modelListOld[i].levelModelList[0].timeIntegration.t,modelListOld[i].levelModelList[0].timeIntegration.dt])
        if(modelListOld[i].levelModelList[0].shockCapturing is not None):
          saveModel["shockCapturing"].append([modelListOld[i].levelModelList[0].shockCapturing.nSteps, modelListOld[i].levelModelList[0].shockCapturing.nStepsToDelay])
        else:
          saveModel["shockCapturing"].append([]) #need to make numbering consistent

      #Assuming the 0th model is RANS2P
      saveModel["stabilization"].append(modelListOld[0].levelModelList[0].stabilization.nSteps)
      return saveModel

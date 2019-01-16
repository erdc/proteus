"""
Optimized  Two-Phase Reynolds Averaged Navier-Stokes
"""
from __future__ import division
from builtins import str
from builtins import range
import proteus
import sys
import numpy
from proteus import Profiling

#it should probably be associated with the PUMI domain somehow
class Checkpointer:
    "This class is meant to handle the checkpointing process for adapted meshes. Information that's needed to be loaded into hotstart needs to be output and then read in to be handled for data reconstruction"
    def __init__(self,NSobject,frequency=10):
      self.A = "Hello"
      self.B = "World"
      self.NSobject = NSobject
      self.counter = 0
      self.frequency = frequency
    def doSomething(self):
      self.transferInfo()
      self.saveMesh()
      self.counter+=1

      #simplifyModelListOld(self.NSobject.modelListOld)
      modelListOld=EncodeModel(self.NSobject.modelListOld,self.NSobject.systemStepController.t_system_last)
      #pickling is apparently unsafe so we use json to try storing modelListOld
      filename = "checkpointInfo"+str(self.counter)
      f = open(filename, 'w')
      import json
      json.dump(modelListOld.__dict__,f)
      f.close()
      #f2 = open(filename,'w')
      #f2.write("Time for hotstart is "+str(self.NSobject.systemStepController.t_system_last))
      #f2.close()
    def transferInfo(self):
      self.NSobject.PUMI_transferFields()
    def saveMesh(self):
      fileName="checkpoint"+str(self.counter)+"_.smb"
      self.NSobject.pList[0].domain.PUMIMesh.writeMesh(fileName)
    #def loadCheckPoint(self,filename):
      #I want to load the mesh

#Maybe I can form an abstract class for encoding decoding?
class EncodeModel:
      "Grab only necessary components from modelListOld so far consistent only with first-order time integrator" 
      def __init__(self,modelListOld,hotStartTime):
        self.numModels = len(modelListOld)
        self.hotStartTime = hotStartTime
        self.stepController=[]
        self.timeIntegration=[]
        self.shockCapturing=[]
        for i in range(0,self.numModels):
          self.stepController.append([modelListOld[i].stepController.dt_model,modelListOld[i].stepController.t_model,modelListOld[i].stepController.t_model_last, modelListOld[i].stepController.substeps])
          if(hasattr(modelListOld[i].levelModelList[0].timeIntegration,'dtLast')):
            self.timeIntegration.append([modelListOld[i].levelModelList[0].timeIntegration.dt,modelListOld[i].levelModelList[0].timeIntegration.t,modelListOld[i].levelModelList[0].timeIntegration.dt,modelListOld[i].levelModelList[0].timeIntegration.dtLast])
          else:
            self.timeIntegration.append([modelListOld[i].levelModelList[0].timeIntegration.dt,modelListOld[i].levelModelList[0].timeIntegration.t,modelListOld[i].levelModelList[0].timeIntegration.dt])
          if(modelListOld[i].levelModelList[0].shockCapturing is not None):
            self.shockCapturing.append([modelListOld[i].levelModelList[0].shockCapturing.nSteps, modelListOld[i].levelModelList[0].shockCapturing.nStepsToDelay])

          #self.stepController.append(modelListOld[i].stepController.__dict__)
          #self.timeIntegration.append(modelListOld[i].levelModelList[0].timeIntegration.__dict__)
          #if(modelListOld[i].levelModelList[0].shockCapturing is not None):
          #  self.shockCapturing.append(modelListOld[i].levelModelList[0].shockCapturing.__dict__)
        #import pdb; pdb.set_trace()
        #import sys; sys.exit()

class DecodeModel:
      "create a modelListOld that can interact with the post-adapt restart capabilities" 
      def __init__(self,modelListOld):
        self.numModels = modelListOld.numModels
        #self.hotStartTime = hotStartTime
        self.stepController={}
        self.timeIntegration={}
        self.shockCapturing={}
        #for i in range(0,self.numModels):
        #  self.stepController.append([modelListOld[i].stepController.dt_model,modelListOld[i].stepController.t_model,modelListOld[i].stepController.t_model_last, modelListOld[i].stepController.substeps])
        #  if(hasattr(modelListOld[i].levelModelList[0].timeIntegration,'dtLast')):
        #    self.timeIntegration.append([modelListOld[i].levelModelList[0].timeIntegration.dt,modelListOld[i].levelModelList[0].timeIntegration.t,modelListOld[i].levelModelList[0].timeIntegration.dt,modelListOld[i].levelModelList[0].timeIntegration.dtLast])
        #  else:
        #    self.timeIntegration.append([modelListOld[i].levelModelList[0].timeIntegration.dt,modelListOld[i].levelModelList[0].timeIntegration.t,modelListOld[i].levelModelList[0].timeIntegration.dt])
        #  if(modelListOld[i].levelModelList[0].shockCapturing is not None):
        #    self.shockCapturing.append([modelListOld[i].levelModelList[0].shockCapturing.nSteps, modelListOld[i].levelModelList[0].shockCapturing.nStepsToDelay])


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
      print self.counter
      print self.A
      print self.B 
      print self.NSobject
      self.transferInfo()
      self.saveMesh()
      self.counter+=1
    def transferInfo(self):
      self.NSobject.PUMI_transferFields()
    def saveMesh(self):
      fileName="checkpoint"+str(self.counter)+"_.smb"
      self.NSobject.pList[0].domain.PUMIMesh.writeMesh(fileName)


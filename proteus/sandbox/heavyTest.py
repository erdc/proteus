#!/usr/bin/env python
import numpy
import shelve

lightData = None
heavyData = None

class Sol:
    def __init__(self,nn):
        self.sol = numpy.zeros((nn,nn,nn),'d')
        self.DT = 1.0
        self.T=0.0
        self.nStep = 0
    def takeStep(self):
        self.T += self.DT
        self.nStep +=1
        self.sol+=1
    def __getstate__(self):
        heavyData['sol%i' % (self.T,)] = self.sol
        heavyData.sync()
        odict = self.__dict__.copy()
        del odict['sol']
        return odict
    def __setstate__(self,idict):
        self.__dict__.update(idict)
        self.sol=None
    def retrieveValue(self,name,nStep,heavyData):
        self.sol = heavyData[name+ ('%i' % (nStep,))]
    def reinitialize(self,heavyData):
        self.retrieveValue('sol',self.nStep,heavyData)
if __name__ == "__main__":
    lightData = shelve.open('lightData',flag='n',protocol=0)
    heavyData = shelve.open('heavyData',flag='n',protocol=0)
    s = Sol(10)
    lightData['sol%i' % (s.nStep,)] = s
    for i in range(10):
        s.takeStep()
        lightData['sol%i' % (s.nStep,)] = s
        lightData.sync()
    lightData.close()

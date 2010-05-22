from pyadh import *
from pyadh.default_p import *
from lwi_dike import *

"""
The non-conservative level set description of a flume in a two-phase flow
"""

##\ingroup test
#\file ls_flume_2d_p.py
#
# \todo finish ls_flume_2d_p.py

class VolumeAveragedVOFCoefficients(VOFCoefficients):
    from pyadh.ctransportCoefficients import VolumeAveragedVOFCoefficientsEvaluate
    from pyadh.cfemIntegrals import copyExteriorElementBoundaryValuesFromElementBoundaryValues
    def __init__(self,LS_model=-1,V_model=0,RD_model=-1,ME_model=1,EikonalSolverFlag=0,checkMass=False,epsFact=0.0,
                 setParamsFunc=None):
        VOFCoefficients.__init__(self,
                                 LS_model=LS_model,
                                 V_model=V_model,
                                 RD_model=RD_model,
                                 ME_model=ME_model,
                                 EikonalSolverFlag=EikonalSolverFlag,
                                 checkMass=checkMass,
                                 epsFact=epsFact)
        self.setParamsFunc   = setParamsFunc
        self.flowCoefficients=None
        self.q_porosity = None; self.ebq_porosity = None; self.ebqe_porosity = None
    #
    def attachModels(self,modelList):
        VOFCoefficients.attachModels(self,modelList)
        self.flowCoefficients = modelList[self.flowModelIndex].coefficients
        if hasattr(self.flowCoefficients,'q_porosity'):
            self.q_porosity = self.flowCoefficients.q_porosity
        else:
            self.q_porosity = Numeric.ones(modelList[self.modelIndex].q[('u',0)].shape,
                                           Numeric.Float)
            if self.setParamsFunc != None:
                self.setParamsFunc(modelList[self.modelIndex].q['x'],self.q_porosity)
            #
        #
        if hasattr(self.flowCoefficients,'ebq_porosity'):
            self.ebq_porosity = self.flowCoefficients.ebq_porosity
        else:
            self.ebq_porosity = Numeric.ones(modelList[self.modelIndex].ebq[('u',0)].shape,
                                             Numeric.Float)
            if self.setParamsFunc != None:
                self.setParamsFunc(modelList[self.modelIndex].ebq['x'],self.ebq_porosity)
            #
        #
        if hasattr(self.flowCoefficients,'ebqe_porosity'):
            self.ebqe_porosity = self.flowCoefficients.ebqe_porosity
        elif modelList[self.LS_modelIndex].numericalFlux != None:
            self.ebqe_porosity = Numeric.ones(modelList[self.LS_modelIndex].numericalFlux.ebqe[('u',0)].shape,
                                              Numeric.Float)
            if self.setParamsFunc != None:
                self.setParamsFunc(modelList[self.LS_modelIndex].numericalFlux.ebqe['x'],self.ebqe_porosity)
            #
        else:
            self.ebqe_porosity = Numeric.ones((modelList[self.modelIndex].mesh.nExteriorElementBoundaries_global,
                                               modelList[self.modelIndex].nElementBoundaryQuadraturePoints_elementBoundary),
                                              Numeric.Float)
            self.copyExteriorElementBoundaryValuesFromElementBoundaryValues(modelList[self.modelIndex].mesh.exteriorElementBoundariesArray,
                                                                            modelList[self.modelIndex].mesh.elementBoundaryElementsArray,
                                                                            modelList[self.modelIndex].mesh.elementBoundaryLocalElementBoundariesArray,
                                                                            self.ebq_porosity,
                                                                            self.ebqe_porosity)
        #
        
        
    #
    def evaluate(self,t,c):

        if c[('f',0)].shape == self.q_v.shape:
            v = self.q_v
            phi = self.q_phi
            porosity  = self.q_porosity
        elif c[('f',0)].shape == self.ebq_v.shape:
            v = self.ebq_v
            phi = self.ebq_phi
            porosity  = self.ebq_porosity
        elif ((self.ebqe_v != None and self.ebqe_phi != None) and c[('f',0)].shape == self.ebqe_v.shape):
            v = self.ebqe_v
            phi = self.ebqe_phi
            porosity  = self.ebqe_porosity
        else:
            v=None
            phi=None
            porosity=None
        if v != None:
            self.VolumeAveragedVOFCoefficientsEvaluate(self.eps,
                                                       v,
                                                       phi,
                                                       porosity,
                                                       c[('u',0)],
                                                       c[('m',0)],
                                                       c[('dm',0,0)],
                                                       c[('f',0)],
                                                       c[('df',0,0)])
    #
#
coefficients = VolumeAveragedVOFCoefficients(LS_model=2,V_model=1,RD_model=4,ME_model=3,
                                             setParamsFunc=setObstaclePorosity)

analyticalSolutions = None
dirichletConditions = {0:getDBC_vof}
initialConditions  = {0:Flat_H()}
fluxBoundaryConditions = {0:'outFlow'}
def getAFBC_H(x):
    pass
def getDFBC_H(x):
    pass

advectiveFluxBoundaryConditions =  {0:getAFBC_H}

diffusiveFluxBoundaryConditions = {0:{}}

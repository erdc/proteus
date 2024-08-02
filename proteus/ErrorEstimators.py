"""
Classes for a posteriori error estimation

.. inheritance-diagram:: proteus.ErrorEstimators
   :parts: 1
"""
from .Profiling import logEvent

class HierarchicalMeshEstimator(object):
    def __init__(self,mlTransport):
        self.mlTransport = mlTransport
    def calculate(self):
        import numpy
        from . import Norms
        from . import FemTools
        t=0.0
        nLevels = len(self.mlTransport.uList)
        assert nLevels > 1, "Need at least two grids for hierarchical mesh estimate"
        coefficients = self.mlTransport.modelList[-1].coefficients
        mCoarse = self.mlTransport.modelList[-2]
        uCoarse = self.mlTransport.modelList[-2].u
        mFine = self.mlTransport.modelList[-1]
        uFine = self.mlTransport.modelList[-1].u
        proj_uCoarse  = [FemTools.FiniteElementFunction(mFine.u[ci].femSpace) for ci in range(coefficients.nc)]
        proj_uCoarse_q   = [numpy.zeros(mFine.q[('u',ci)].shape,'d') for ci in range(coefficients.nc)]
        elementError = [numpy.zeros((mFine.mesh.nElements_global,),'d') for ci in range(coefficients.nc)]
        elementTagArray = numpy.zeros((mFine.mesh.nElements_global,),'i')
        elementTagArray.flat[:]=0.0
        localRefinement=False
        error = 0.0
        for ci in range(coefficients.nc):
            self.mlTransport.meshTransfers.prolong_bcListDict[ci][-1].matvec(uCoarse[ci].dof,
                                                                             proj_uCoarse[ci].dof)
            #load Dirichlet conditions in
            for dofN,g in mFine.dirichletConditions[ci].DOFBoundaryConditionsDict.items():
                proj_uCoarse[ci].dof[dofN] = g(mFine.dirichletConditions[ci].DOFBoundaryPointDict[dofN],t)
            proj_uCoarse[ci].getValues(mFine.q['v',ci],proj_uCoarse_q[ci])
            error += Norms.L2errorSFEM_local(mFine.q[('dV_u',ci)],
                                             mFine.q[('u',ci)],
                                             proj_uCoarse_q[ci],
                                             elementError[ci])
            for eN in range(mFine.mesh.nElements_global):
                if (mFine.mesh.nElements_global*elementError[ci][eN]/error) > 5.0:
                    elementTagArray[eN] = 1
                    localRefinement=True
        if not localRefinement:
            elementTagArray.flat[:]=1
        logEvent("error "+repr(error),level=3)#mwf debug turn off,"elementTagArray",elementTagArray
        return (error,elementTagArray)

"""
Modules for computing numerical solutions of differential equations
"""

try:
    import pkg_resources
    pkg_resources.declare_namespace(__name__)
except ImportError:
    import pkgutil
    __path__ = pkgutil.extend_path(__path__, __name__)

__version__ = '1.8.3'

__all__ = ["Archiver",
           "Domain",
           "InputTranslators",
           "SplitOperator",
           "StepControl",
           "NumericalSolution",
           "Comm",
           "AnalyticalSolutions",
           "TransportCoefficients",
           "TwophaseDarcyCoefficients",
           "DiagUtils",
           "EGeometry",
           "ErrorEstimators",
           "FemTools",
           "LatexReport",
           "LinearAlgebraTools",
           "LinearSolvers",
           "MeshTools",
           "NonlinearSolvers",
           "Norms",
           "NumericalFlux",
           "PostProcessingTools",
           "Profiling",
           "Quadrature",
           "RefUtils",
           "ShockCapturing",
           "SimTools",
           "SubgridError",
           "SubsurfaceTransportCoefficients",
           "TimeIntegration",
           "Transport",
           "Viewers",
           "AuxiliaryVariables",
           "deim_utils",
           "default_p",
           "default_n",
           "default_s",
           "default_so",
           "InputTranslators",
           "cfemIntegrals",
           "cshockCapturing",
           "csmoothers",
           "csubgridError",
           "ctimeIntegration",
           "ctransportCoefficients",
           "clapack",
           "superluWrappers",
           "cmeshTools",
           "cnumericalFlux",
           "cTwophaseDarcyCoefficients",
           "cADR",
           "deim_utils",
           "WaveTools",
           "Context",
           "BoundaryConditions",
           "SpatialTools",
           "defaults"]

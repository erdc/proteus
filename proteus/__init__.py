"""
Modules for computing numerical solutions of differential equations
"""

try:
    import pkg_resources
    pkg_resources.declare_namespace(__name__)
except ImportError:
    import pkgutil
    __path__ = pkgutil.extend_path(__path__, __name__)

__version__ = '1.5.1.dev0'

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
           "StupidHeap",
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
           "ADR",
           "deim_utils",
           "WaveTools",
           "Context",
           "BoundaryConditions",
           "SpatialTools",
           "defaults"]

def test(verbose=False, cleanup=True):
    """Run all proteus tests

    Parameters
    ----------
    verbose : bool
              Print verbose testing information
    cleanup : bool
              Remove the temporary directory containing output
    """
    from os import path
    from tempfile import mkdtemp
    from shutil import rmtree
    import pytest
    flags="--boxed "
    if verbose:
        flags+="-v "
    original_dir = os.get_cwd()
    tmp_dir = mkdtemp()
    os.chdir(tmp_dir)
    pytest.main(flags+path.join(path.dirname(__file__),'tests'))
    os.chdir(original_dir)
    if cleanup:
        rmtree(tmp_dir)

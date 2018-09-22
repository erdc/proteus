from __future__ import division
from past.utils import old_div
from proteus import FemTools as ft
from proteus import MeshTools as mt

class SWFlowProblem:
    """ SWFlowProblem """

    def __init__(self,
                 sw_model=0, #0: shallow water, 1: disp shallow water
                 # TIME STEPPING #
                 cfl=0.33,
                 outputStepping=None,
                 # DOMAIN AND MESH #
                 structured=False,
                 he=None,
                 nnx=None,
                 nny=None,
                 domain=None,
                 AdH_file=None,
                 # BATHYMETRY #
                 bathymetry=None,
                 triangleFlag=1,
                 # INITIAL CONDITIONS #
                 initialConditions=None,
                 # BOUNDARY CONDITIONS #
                 boundaryConditions=None,
                 reflectingBCs=False,
                 # OTHERS #
                 useSuperlu=None):
        """ Constructor for structured meshes  """
        # ***** SET OF ASSERTS ***** #
        assert sw_model in [0,1], "sw_model={0,1} for shallow water equations or dispersive shallow water equations respectively"
        assert cfl <= 1, "Choose cfl <= 1"
        assert isinstance (outputStepping,OutputStepping), "Provide an object from the OutputStepping class"
        assert type(he)==float , "Provide (float) he (characteristic mesh size)"
        if structured:
            assert type(nnx)==int and type(nny)==int, "Provide (int) nnx and nny"
        if domain is None:
            assert AdH_file is not None, "If domain is None then provide an AdH File"
        else:
            assert callable(bathymetry), "Bathymetry must be a function"
        assert triangleFlag in [0,1,2], "triangleFlag must be 1, 2 or 3"
        assert type(initialConditions)==dict, "Provide dict of initial conditions"
        assert type(boundaryConditions)==dict, "Provide dict of boundary conditions"
        self.assert_initialConditions(sw_model,initialConditions)
        self.assert_boundaryConditions(sw_model,boundaryConditions)

        # ***** SAVE PARAMETERS ***** #
        self.sw_model=sw_model
        self.cfl=cfl
        self.outputStepping=outputStepping.getOutputStepping()
        self.he=he
        self.nnx=nnx
        self.nny=nny
        self.nnz=1
        self.domain=domain
        self.AdH_file=AdH_file
        self.bathymetry=bathymetry
        self.triangleFlag=triangleFlag
        self.initialConditions=initialConditions
        self.boundaryConditions=boundaryConditions
        self.reflectingBCs=reflectingBCs
        self.useSuperlu = useSuperlu

        # ***** CREATE FINITE ELEMENT SPACES ***** #
        self.FESpace = FESpace().getFESpace()

        # ***** DEFINE PHYSICAL AND NUMERICAL PARAMETERS ***** #
        self.physical_parameters = default_physical_parameters
        self.swe_parameters = default_swe_parameters
        self.dswe_parameters = default_dswe_parameters

    def assert_initialConditions(self,sw_model,initialConditions):
        assert 'water_height' in initialConditions, 'Provide water_height in ICs'
        assert 'x_mom' in initialConditions, 'Provide y_mom in ICs'
        assert 'y_mom' in initialConditions, 'Provide x_mom in ICs'
        if sw_model==1: # dispersive SWEs
            assert 'h_times_eta' in initialConditions, 'Provide auxiliary function h*eta in ICs'
            assert 'h_times_w' in initialConditions, 'Provide auxiliary function h*w in ICs'
    #
    def assert_boundaryConditions(self,sw_model,boundaryConditions):
        # check dirichlet BCs
        assert 'water_height' in boundaryConditions, 'Provide water_height in BCs'
        assert 'x_mom' in boundaryConditions, 'Provide y_mom in BCs'
        assert 'y_mom' in boundaryConditions, 'Provide x_mom in BCs'
        if sw_model==1: # dispersive SWEs
            assert 'h_times_eta' in boundaryConditions, 'Provide auxiliary function h*eta in BCs'
            assert 'h_times_w' in boundaryConditions, 'Provide auxiliary function h*w in BCs'
            
class OutputStepping:
    """
    OutputStepping handles how often the solution is outputted.
    """
    def __init__(self,
                 final_time,
                 dt_init=0.001,
                 dt_output=None,
                 nDTout=None):
        self.final_time=final_time
        self.dt_init=dt_init
        assert not (dt_output is None and nDTout is None), "Provide dt_output or nDTout"
        self.dt_output=dt_output
        self.nDTout = nDTout
    #
    def getOutputStepping(self):
        # COMPUTE dt_init #
        dt_init = min(0.1 * self.dt_output, self.dt_init)
        if self.nDTout is None:
            self.nDTout = int(round(old_div(self.final_time, self.dt_output)))
        else:
            self.dt_output = float(self.final_time)/float(self.nDTout)
        #
        return {'final_time': self.final_time,
                'dt_init': dt_init,
                'dt_output': self.dt_output,
                'nDTout': self.nDTout}
#
class FESpace:
    """
    Create FE Spaces.
    """
    def __init__(self):
        pass

    def getFESpace(self):
        basis = ft.C0_AffineLinearOnSimplexWithNodalBasis # p1 space
        # QUADRATURE RULE #
        elementQuadrature = ft.SimplexGaussQuadrature(2, 3)
        elementBoundaryQuadrature = ft.SimplexGaussQuadrature(1, 3)
        return {'basis': basis,
                'elementQuadrature': elementQuadrature,
                'elementBoundaryQuadrature': elementBoundaryQuadrature}

# ***************************************** #
# ********** PHYSICAL PARAMETERS ********** #
# ***************************************** #
default_physical_parameters ={'gravity': 9.8,
                              'LINEAR_FRICTION': 0,
                              'mannings': 0.0}

# ****************************************** #
# ********** NUMERICAL PARAMETERS ********** #
# ****************************************** #
default_swe_parameters = {'LUMPED_MASS_MATRIX': 0,
                          'cfl': 0.33,
                          'SSPOrder': 3,
                          'cE': 1}
default_dswe_parameters = {}

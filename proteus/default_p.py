"""
The default values for problem definition modules

.. todo::

  Clean up default_p.py
"""
from MeshTools import *
from FemTools import *
from TransportCoefficients import *
from Transport import *

name = None
"""The name of the model

:type: None or string
"""

nd = 1
"""The number of spatial dimensions of the model domain"""

domain = None
"""The domain object

:type: None or proteus.Domain.D_base
"""

movingDomain=False
"""Boolean to indicate moving domain"""

polyfile = None
"""The filename of a polyfile giving the domain"""

meshfile = None
"""The filename of a mesh giving the mesh/domain"""

genMesh = True
"""Boolean to trigger mesh generation"""

L=(1.0,1.0,1.0)
"""Tuple of dimensions for simple box shaped domain"""

x0=(0.0, 0.0, 0.0)
"""Tuple of coordinates for corner of simple box shaped domain"""

analyticalSolution = {}
"""Dictionary of analytical solutions for each component

Each element should be an object of type :class:`proteus.AnalyticalSolutions.AS_base`"""

coefficients = None
"""Transport coefficients object

The object should be of type :class:`proteus.TransportCoefficients.TC_base`"""

dirichletConditions = {}
"""Dictionary of Dirichlet conditions for each component"""

periodicDirichletConditions = None
"""Dictionary of periodic boundary conditions for each component"""

fluxBoundaryConditions = {}
r"""Dictionary of flux boundary condition flags for each component ('outFlow','noFlow','setFlow','mixedFlow')"""

advectiveFluxBoundaryConditions =  {}
"""Dictionary of advective flux boundary conditions setter functions"""

diffusiveFluxBoundaryConditions = {}
"""Dictionary of diffusive flux boundary conditions setter functions"""

stressFluxBoundaryConditions = {}
"""Dictionary of stress tensor flux boundary conditions setter functions"""

initialConditions = None
"""Dictionary of initial condition function objects"""

weakDirichletConditions = None
"""Dictionary of weak Dirichlet constraint setters"""

dummyInitialConditions = False #mwf temporary hack for RD level sets

boundaryCreatesNullSpace = False
"""Boolean indicating whether the boundary conditions imply a null space exits"""

#what to do after each coupling step?
def finalizeStepDummy(c):
    pass

finalizeStep = finalizeStepDummy

T=1.0
"""End of time interval"""

sd = True
"""Use sparse representation of diffusion tensors"""

LevelModelType = OneLevelTransport

#!/usr/bin/env python
"""

Test module for generating unit_simplex elements

"""
from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
Profiling.logLevel=7
Profiling.verbose=True
import numpy.testing as npt
from nose.tools import ok_ as ok
from nose.tools import eq_ as eq
from nose.tools import set_trace

def generate_reference_triangle():
    '''
    This modeul tests the generation of reference simplex meshes.
    '''
    unit_simplex_domain_2d = Domain.unitSimplex(2)

    assert unit_simplex_domain_2d.x[0] == 0.0
    assert unit_simplex_domain_2d.x[1] == 0.0
    assert unit_simplex_domain_2d.L[0] == 1.0
    assert unit_simplex_domain_2d.L[1] == 1.0

    unit_simplex_2d = MeshTools.buildReferenceSimplex(nd=2)

    assert len(unit_simplex_2d.nodeArray) == 3
    assert any(i in unit_simplex_2d.nodeArray[0] for i in [0.,0.,0.])
    assert any(i in unit_simplex_2d.nodeArray[1] for i in [0.,1.,0.])
    assert any(i in unit_simplex_2d.nodeArray[2] for i in [1.,0.,0.])
    assert unit_simplex_2d.globalMesh == unit_simplex_2d, \
        'Reference simplex should only have a single mesh global mesh.'
    assert unit_simplex_2d.nodeOffsets_subdomain_owned[1] == 3, \
        'Mesh partitioner has likely not been called. '

def generate_reference_simplex():

    unit_simplex_domain_3d = Domain.unitSimplex(3)

    assert unit_simplex_domain_3d.x[0] == 0.0
    assert unit_simplex_domain_3d.x[1] == 0.0
    assert unit_simplex_domain_3d.x[2] == 0.0
    assert unit_simplex_domain_3d.L[0] == 1.0
    assert unit_simplex_domain_3d.L[1] == 1.0
    assert unit_simplex_domain_3d.L[2] == 1.0
    
    polyfile="reference_simplex"
    unit_simplex_domain_3d.writePoly(polyfile)

    # TEST FAILS
    assert len(mesh.nodeArray) == 4
#    assert len(mesh.edgeArray) = 


if __name__ == '__main__':
#    from proteus import Comm
#    comm = Comm.init()
#    import nose
#    nose.main()
    generate_reference_triangle()
#    generate_reference_simplex()

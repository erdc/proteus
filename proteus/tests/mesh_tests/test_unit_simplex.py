#!/usr/bin/env python
"""

Test module for generating unit_simplex elements

"""
import pytest

from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
Profiling.logLevel=7
Profiling.verbose=True

@pytest.mark.MeshTools
class TestReferenceSimplex():

    @classmethod
    def setup_class(cls):
        pass

    @classmethod
    def teardown_class(cls):
        pass

    def setup_method(self,method):
        pass

    def teardown_method(self,method):
        ''' Tear down the test problem '''
        filenames = ['reference_element.ele', 'reference_element.face',
                     'reference_element.node', 'reference_element.poly']
        for file in filenames:
            if os.path.exists(file):
                try:
                    os.remove(file)
                except OSError, e:
                    print ("Error: %s - %s." %(e.filename,e.strerror))
            else:
                pass
    
    def test_generate_reference_triangle(self):
        ''' Basic tests for the unit triangle function '''
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

    def test_generate_reference_simplex(self):
        """ Basic tests for the unit simplex function """
        
        unit_simplex_domain_3d = Domain.unitSimplex(3)

        assert len(unit_simplex_domain_3d.vertices) == 4
        assert unit_simplex_domain_3d.x[0] == 0.0
        assert unit_simplex_domain_3d.x[1] == 0.0
        assert unit_simplex_domain_3d.x[2] == 0.0
        assert unit_simplex_domain_3d.L[0] == 1.0
        assert unit_simplex_domain_3d.L[1] == 1.0
        assert unit_simplex_domain_3d.L[2] == 1.0
        
        unit_simplex_3d = MeshTools.buildReferenceSimplex(nd=3)
        
        assert len(unit_simplex_3d.nodeArray) == 4
        assert any(i in unit_simplex_3d.nodeArray[0] for i in [0.,0.,0.])
        assert any(i in unit_simplex_3d.nodeArray[1] for i in [0.,0.,1.])
        assert any(i in unit_simplex_3d.nodeArray[2] for i in [0.,1.,0.])
        assert any(i in unit_simplex_3d.nodeArray[3] for i in [1.,0.,0.])
        



if __name__ == '__main__':
    pass

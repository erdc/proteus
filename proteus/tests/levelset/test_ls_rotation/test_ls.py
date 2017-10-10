#!/usr/bin/env python
"""
Test level set solver
"""

from proteus.iproteus import *
from proteus import Comm
from proteus import Context
import tables

comm = Comm.get()
Profiling.logLevel=2
Profiling.verbose=True
import numpy as np

class TestLS():

    @classmethod
    def setup_class(cls):
        pass

    @classmethod
    def teardown_class(cls):
        pass

    def setup_method(self,method):
        """Initialize the test problem. """
        self._scriptdir = os.path.dirname(__file__)
        self.sim_names = []
        self.aux_names = []
        self.opts =Context.Options([
('useTextArchive', False,""),
('generatePartitionedMeshFromFiles', False,""),
('petscOptions', None,""),
('subdomainArchives', False,""),
('viewer', False,""),
('gatherArchive', False,""),
('dataDir', '',""),
('save_dof', False,""),
('probDir', '.',""),
('global_sync', True,""),
('viewLevels', False,""),
('logLevel', 5,""),
('hotStart', False,""),
('petscOptionsFile', None,""),
('writeVPP', False,""),
('profile', False,""),
('contextOptionsFile', None,""),
('batchFileName', '',""),
('inspect', '',""),
('setupOnly', False,""),
('ensight', False,""),
('cacheArchive', False,""),
('wait', False,""),
('viewMesh', False,""),
('memHardLimit', -1.0,""),
('logAllProcesses', False,""),
('debug', True,""),
('interactive', '', ""),],mutable=True)

    def teardown_method(self,method):
        pass

    def test_ex1(self):

        Context.contextOptionsString="T=2"

        import rotation2D
        from rotation2D import *
        from ls_rotation_2d_so import *

        self.pList = [__import__(pnList[0][0])]
        self.nList = [__import__(pnList[0][1])]
        self.sList = [default_s]
        self.so = default_so
        self.so.name = soname
        self.so.tnList = tnList

        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               self.pList,
                                               self.nList,
                                               self.sList,
                                               self.opts)
        self.sim_names.append(ns.modelList[0].name)
        ns.calculateSolution(soname)
        # COMPARE VS SAVED FILES #
        expected_path = 'comparison_files/rotation_c0p1cg_vbdf_2_level_1_ls.h5'
        expected = tables.open_file(os.path.join(self._scriptdir,expected_path))
        actual = tables.open_file(os.path.join(self._scriptdir,
                                               'rotation_c0p1cg_vbdf_2_level_1.h5'))
        assert np.allclose(expected.root.u_t2,
                           actual.root.u_t2,
                           atol=1e-10)
        expected.close()
        actual.close()

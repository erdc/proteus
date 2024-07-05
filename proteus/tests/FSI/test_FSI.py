#!/usr/bin/env python
import os
import pytest
from proteus.iproteus import opts, default_s
from proteus import Profiling, NumericalSolution, Comm
import unittest
import numpy as np
import numpy.testing as npt
from importlib import import_module
from petsc4py import PETSc
import importlib
import pytest

modulepath = os.path.dirname(os.path.abspath(__file__))
comm = Comm.init()

class TestIBM(unittest.TestCase):

    def teardown_method(self, method):
        """ Tear down function """
        FileList = ['addedmass2D.xmf',
                    'addedmass2D.h5',
                    'addedmass3D.xmf',
                    'addedmass3D.h5',
                    'record_rectangle1.csv',
                    'record_rectangle1_Aij.csv',
                    'record_cuboid1.csv',
                    'record_cuboid1_Aij.csv',
                    'mesh.ele',
                    'mesh.edge',
                    'mesh.node',
                    'mesh.neigh',
                    'mesh.face',
                    'mesh.poly',
                    'forceHistory_p.txt',
                    'forceHistory_v.txt',
                    'momentHistory.txt',
                    'wettedAreaHistory.txt',
                    'proteus.log'
                    'addedmass2D_so.log'
                    'addedmass3D_so.log'
                    ]
        for file in FileList:
            if os.path.isfile(file) and comm.isMaster():
                os.remove(file)
            else:
                pass

    @pytest.mark.skip(reason="need to redo after history revision")                         
    def test_fallingCylinderIBM_ball(self):
        from . import fallingCylinder
        from proteus import defaults
        case = fallingCylinder
        case.m.rans2p.p.coefficients.use_ball_as_particle = 1.
        case.myTpFlowProblem.initializeAll()
        so = case.myTpFlowProblem.so
        so.name = 'fallingCylinderIBM'
        pList = []
        nList = []
        for (pModule,nModule) in so.pnList:
            pList.append(pModule)
            nList.append(nModule)
        if so.sList == []:
            for i in range(len(so.pnList)):
                s = default_s
                so.sList.append(s)
        Profiling.logLevel = 7
        Profiling.verbose = True
        # PETSc solver configuration
        OptDB = PETSc.Options()
        dirloc = os.path.dirname(os.path.realpath(__file__))
        with open(os.path.join(dirloc, "petsc.options.superlu_dist")) as f:
            all = f.read().split()
            i=0
            while i < len(all):
                if i < len(all)-1:
                    if all[i+1][0]!='-':
                        print("setting ", all[i].strip(), all[i+1])
                        OptDB.setValue(all[i].strip('-'),all[i+1])
                        i=i+2
                    else:
                        print("setting ", all[i].strip(), "True")
                        OptDB.setValue(all[i].strip('-'),True)
                        i=i+1
                else:
                    print("setting ", all[i].strip(), "True")
                    OptDB.setValue(all[i].strip('-'),True)
                    i=i+1
        ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
        ns.calculateSolution('fallingCylinderIBM')
        pos = case.body.getPosition()

        npt.assert_almost_equal(pos, np.array([1.5, 1.98645, 0.]), decimal=5)
        #self.teardown_method(self)

    @pytest.mark.skip(reason="need to redo after history revision")                         
    def test_fallingCylinderIBM_sdf(self):
        from proteus import defaults
        from . import fallingCylinder
        importlib.reload(fallingCylinder)
        case = fallingCylinder
        case.m.rans2p.p.coefficients.use_ball_as_particle = 0.
        case.myTpFlowProblem.initializeAll()
        so = case.myTpFlowProblem.so
        so.name = 'fallingCylinderIBM2'
        pList = []
        nList = []
        for (pModule,nModule) in so.pnList:
            pList.append(pModule)
            nList.append(nModule)
        if so.sList == []:
            for i in range(len(so.pnList)):
                s = default_s
                so.sList.append(s)
        Profiling.logLevel = 7
        Profiling.verbose = True
        # PETSc solver configuration
        OptDB = PETSc.Options()
        dirloc = os.path.dirname(os.path.realpath(__file__))
        with open(os.path.join(dirloc, "petsc.options.superlu_dist")) as f:
            all = f.read().split()
            i=0
            while i < len(all):
                if i < len(all)-1:
                    if all[i+1][0]!='-':
                        print("setting ", all[i].strip(), all[i+1])
                        OptDB.setValue(all[i].strip('-'),all[i+1])
                        i=i+2
                    else:
                        print("setting ", all[i].strip(), "True")
                        OptDB.setValue(all[i].strip('-'),True)
                        i=i+1
                else:
                    print("setting ", all[i].strip(), "True")
                    OptDB.setValue(all[i].strip('-'),True)
                    i=i+1
        ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
        ns.calculateSolution('fallingCylinderIBM2')
        pos = case.body.getPosition()

        npt.assert_almost_equal(pos, np.array([1.5, 1.98645, 0.]), decimal=5)
        #self.teardown_method(self)

    @pytest.mark.skipif(os.sys.platform == "darwin", reason="does not run on macOS")    
    def test_floatingCylinderALE(self):
        from proteus import defaults
        from . import floatingCylinder
        case = floatingCylinder
        case.myTpFlowProblem.initializeAll()
        so = case.myTpFlowProblem.so
        so.name = 'floatingCylinderALE'
        pList = []
        nList = []
        for (pModule,nModule) in so.pnList:
            pList.append(pModule)
            nList.append(nModule)
        if so.sList == []:
            for i in range(len(so.pnList)):
                s = default_s
                so.sList.append(s)
        Profiling.logLevel = 7
        Profiling.verbose = True
        # PETSc solver configuration
        OptDB = PETSc.Options()
        dirloc = os.path.dirname(os.path.realpath(__file__))
        with open(os.path.join(dirloc, "petsc.options.superlu_dist")) as f:
            all = f.read().split()
            i=0
            while i < len(all):
                if i < len(all)-1:
                    if all[i+1][0]!='-':
                        print("setting ", all[i].strip(), all[i+1])
                        OptDB.setValue(all[i].strip('-'),all[i+1])
                        i=i+2
                    else:
                        print("setting ", all[i].strip(), "True")
                        OptDB.setValue(all[i].strip('-'),True)
                        i=i+1
                else:
                    print("setting ", all[i].strip(), "True")
                    OptDB.setValue(all[i].strip('-'),True)
                    i=i+1
        ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
        ns.calculateSolution('floatingCylinderALE')
        pos = case.body.getPosition()

        npt.assert_almost_equal(pos, np.array([0.5, 0.5074055958, 0.]), decimal=5)
        #self.teardown_method(self)

    @pytest.mark.skipif(os.sys.platform == "darwin", reason="does not run on macOS")
    def test_floatingCubeALE(self):
        from proteus import defaults
        from . import floatingCube
        case = floatingCube
        case.myTpFlowProblem.initializeAll()
        so = case.myTpFlowProblem.so
        so.name = 'floatingCubeALE'
        pList = []
        nList = []
        for (pModule,nModule) in so.pnList:
            pList.append(pModule)
            nList.append(nModule)
        if so.sList == []:
            for i in range(len(so.pnList)):
                s = default_s
                so.sList.append(s)
        Profiling.logLevel = 7
        Profiling.verbose = True
        # PETSc solver configuration
        OptDB = PETSc.Options()
        dirloc = os.path.dirname(os.path.realpath(__file__))
        with open(os.path.join(dirloc, "petsc.options.superlu_dist")) as f:
            all = f.read().split()
            i=0
            while i < len(all):
                if i < len(all)-1:
                    if all[i+1][0]!='-':
                        print("setting ", all[i].strip(), all[i+1])
                        OptDB.setValue(all[i].strip('-'),all[i+1])
                        i=i+2
                    else:
                        print("setting ", all[i].strip(), "True")
                        OptDB.setValue(all[i].strip('-'),True)
                        i=i+1
                else:
                    print("setting ", all[i].strip(), "True")
                    OptDB.setValue(all[i].strip('-'),True)
                    i=i+1
        ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
        ns.calculateSolution('floatingCubeALE')
        pos = case.body.getPosition()

        npt.assert_almost_equal(pos, np.array([0.49945, 0.50047, 0.50807]), decimal=5)
        #self.teardown_method(self)

if __name__ == "__main__":
    unittest.main()

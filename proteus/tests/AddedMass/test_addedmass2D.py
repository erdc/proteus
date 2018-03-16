#!/usr/bin/env python
import os
from proteus.iproteus import *
import unittest
import numpy.testing as npt
from importlib import import_module
from petsc4py import PETSc

class TestAddedMass2D(unittest.TestCase):

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
            if os.path.isfile(file):
                os.remove(file)
            else:
                pass

    def test_AddedMass_2D(self):
        from proteus import default_so
        reload(default_so)
        import addedmass2D_so
        pList = []
        nList = []
        for (p,n) in addedmass2D_so.pnList:
            pList.append(import_module("."+p,"proteus.tests.AddedMass"))
            nList.append(import_module("."+n,"proteus.tests.AddedMass"))
            if pList[-1].name == None:
                pList[-1].name = p
        so = addedmass2D_so
        so.name = "addedmass2D"
        if so.sList == []:
            for i in range(len(so.pnList)):
                s = default_s
                so.sList.append(s)
        Profiling.logLevel=7
        Profiling.verbose=True
        # PETSc solver configuration
        OptDB = PETSc.Options()
        dirloc = os.path.dirname(os.path.realpath(__file__))
        with open(os.path.join(dirloc, "petsc.options.superlu_dist")) as f:
            all = f.read().split()
            i=0
            while i < len(all):
                if i < len(all)-1:
                    if all[i+1][0]!='-':
                        print "setting ", all[i].strip(), all[i+1]
                        OptDB.setValue(all[i].strip('-'),all[i+1])
                        i=i+2
                    else:
                        print "setting ", all[i].strip(), "True"
                        OptDB.setValue(all[i].strip('-'),True)
                        i=i+1
                else:
                    print "setting ", all[i].strip(), "True"
                    OptDB.setValue(all[i].strip('-'),True)
                    i=i+1
        ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
        ns.calculateSolution('addedmass2D')
        Aij = ns.so.ct.body.Aij
        npt.assert_almost_equal(Aij[0,0], 500.17986551764352, decimal=5)
        npt.assert_almost_equal(Aij[1,1], 1299.5165614689747, decimal=4)
        self.teardown_method(self)

if __name__ == "__main__":
    unittest.main()

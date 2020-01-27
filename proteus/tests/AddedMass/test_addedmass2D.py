#!/usr/bin/env python
from __future__ import print_function
from __future__ import absolute_import
from builtins import range
import os
from proteus.iproteus import *
import unittest
import numpy as np
import numpy.testing as npt
from importlib import import_module
from petsc4py import PETSc

modulepath = os.path.dirname(os.path.abspath(__file__))

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
        from proteus import defaults
        from . import addedmass2D as am2D
        am2D.myTpFlowProblem.initializeAll()
        so = am2D.myTpFlowProblem.so
        so.name = 'addedMass2D'
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
        ns.calculateSolution('addedmass2D')
        Aij = am2D.body.Aij

        #np.savetxt('Aij_sol2D.csv', Aij, delimiter=',')
        Aij_sol = np.genfromtxt(os.path.join(modulepath, 'Aij_sol2D.csv'), delimiter=',')
        npt.assert_almost_equal(Aij, Aij_sol, decimal=5)
        self.teardown_method(self)

if __name__ == "__main__":
    unittest.main()

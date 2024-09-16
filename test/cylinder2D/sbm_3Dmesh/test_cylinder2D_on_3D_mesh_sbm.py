"""Tests for 2d flow around a cylinder with shifted boundary method on a 3D mesh"""
from proteus.iproteus import *
from proteus import Comm, defaults
from proteus import Context
import h5py
import importlib
import os
import numpy as np

comm = Comm.get()
Profiling.logLevel = 7
Profiling.verbose = False

modulepath = os.path.dirname(os.path.abspath(__file__))

class Test_sbm_cylinder2D_on_mesh3D(object):

    @classmethod
    def setup_class(cls):
        cls._scriptdir = os.path.dirname(os.path.abspath(__file__))
    @classmethod
    def teardown_class(cls):
        pass

    def setup_method(self, method):
        pass

    def teardown_method(self, method):
        """ Tear down function """
        FileList = ['cylinder_sbm_mesh3D_T001_P1_sbm_3Dmesh.h5', 'cylinder_sbm_mesh3D_T001_P1_sbm_3Dmesh.xmf',
                   ]
        for file in FileList:
            if os.path.isfile(file):
                os.remove(file)
            else:
                pass



    def test_ex1(self):
        self.compare_name = "T001_P1_sbm_3Dmesh"
        self.example_setting("T=0.01 spaceOrder=1 onlySaveFinalSolution=True")
        self.teardown_method(self)

    # really slow
#     def test_ex2(self):
#         self.compare_name = "T001_P2_sbm_3Dmesh"
#         self.example_setting("T=0.01 spaceOrder=2 onlySaveFinalSolution=True")


    def example_setting(self, pre_setting):
        from petsc4py import PETSc
        import sys
        OptDB = PETSc.Options()
        OptDB.clear()

        Context.contextOptionsString = pre_setting

        my_so = defaults.load_system("cylinder_so",modulepath)

        opts.profile = False
        opts.gatherArchive = True
        
        pList=[]
        nList=[]
        sList=[]
        for (pModule,nModule) in my_so.pnList:
            pList.append(defaults.load_physics(pModule,modulepath))
            sys.modules[pModule]=pList[-1]
            nList.append(defaults.load_numerics(nModule,modulepath))
            sys.modules[nModule]=nList[-1]
            if pList[-1].name == None:
                pList[-1].name = pModule

        if my_so.sList == []:
            for i in range(len(my_so.pnList)):
                s = default_s
                sList.append(s)
        else:
            sList = my_so.sList

        my_so.name += "_sbm_mesh3D_"+self.compare_name #save data with different filename
        try:
            ns = proteus.NumericalSolution.NS_base(my_so,
                                               pList,
                                               nList,
                                               sList,
                                               opts)
        except:
            assert 0, "NS setup failed"
        try:
            ns.calculateSolution(my_so.name)
        except:
            assert 0, "NS calculation failed"

        actual = h5py.File('cylinder_sbm_mesh3D_T001_P1_sbm_3Dmesh'+'.h5','r')
        expected_path = 'comparison_files/' + 'comparison_u_t2.csv'
        #write comparison file
        #np.array(actual.root.u_t2).tofile(os.path.join(self._scriptdir, expected_path),sep=",")
        np.testing.assert_almost_equal(np.fromfile(os.path.join(self._scriptdir, expected_path),sep=","),np.array(actual['u_t2']),decimal=10)

        actual.close()

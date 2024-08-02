"""Tests for 2d flow around an immersed boundary cylinder with rans3p"""
import pytest
from proteus.iproteus import *
from proteus import Comm
from proteus import Context
import h5py
import importlib


comm = Comm.get()
Profiling.logLevel = 7
Profiling.verbose = False
import numpy as np


class Test_ibm(object):

    @classmethod
    def setup_class(cls):
        cls._scriptdir = os.path.dirname(os.path.abspath(__file__))
        sys.path.insert(0,cls._scriptdir)
    @classmethod
    def teardown_class(cls):
        sys.path.remove(cls._scriptdir)
        pass

    def setup_method(self, method):
        """Initialize the test problem. """
        self.aux_names = []

    def teardown_method(self, method):
        """ Tear down function """
        FileList = ['cylinder_ibm_T1_rans3p.h5', 'cylinder_ibm_T1_rans3p.xmf'
                    ]
        for file in FileList:
            if os.path.isfile(file):
                os.remove(file)
            else:
                pass

#     def test_ex1(self):
#         self.compare_name = "T8P2"
#         self.example_setting("T=8.0 vspaceOrder=2 onlySaveFinalSolution=True")

    @pytest.mark.skip(reason="need to redo after history revision")                         
    def test_ex2(self):
        self.compare_name = "T1_rans3p"
        self.example_setting("T=0.01 onlySaveFinalSolution=True")


    def example_setting(self, pre_setting):
        from importlib import reload
        Context.contextOptionsString = pre_setting
        from . import cylinder_so as my_so
        reload(my_so)
        # defined in iproteus
        opts.profile = False
        opts.gatherArchive = True        
        pList=[]
        nList=[]
        sList=[]
        for (pModule,nModule) in my_so.pnList:
            pList.append(
                importlib.import_module("."+pModule,
                                        "proteus.tests.cylinder2D.ibm_method"))
            nList.append(
                importlib.import_module("."+nModule,
                                        "proteus.tests.cylinder2D.ibm_method"))
            if pList[-1].name == None:
                pList[-1].name = pModule
            reload(pList[-1])  # Serious error
            reload(nList[-1])

        if my_so.sList == []:
            for i in range(len(my_so.pnList)):
                s = default_s
                sList.append(s)
        else:
            sList = my_so.sList
        my_so.name += "_ibm_"+self.compare_name #save data with different filename
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(my_so,
                                               pList,
                                               nList,
                                               sList,
                                               opts)
        self.aux_names.append(ns.modelList[0].name)
        ns.calculateSolution(my_so.name)
        actual = h5py.File( my_so.name + '.h5')
        expected_path = 'comparison_files/' + 'comparison_' + self.compare_name + '_u_t2.csv'
        #write comparison file
        #np.array(actual.root.u_t2).tofile(os.path.join(self._scriptdir, expected_path),sep=",")
        np.testing.assert_almost_equal(np.fromfile(os.path.join(self._scriptdir, expected_path),sep=","),np.array(actual['u_t2']).flatten(),decimal=10)
        actual.close()

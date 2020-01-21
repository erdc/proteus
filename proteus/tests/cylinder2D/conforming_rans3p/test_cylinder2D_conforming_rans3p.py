"""Tests for 2d flow around a cylinder with a conforming mesh and rans3p"""
from builtins import range
from builtins import object
from proteus.iproteus import *
from proteus import Comm
from proteus import Context
import tables
import importlib

comm = Comm.get()
Profiling.logLevel = 7
Profiling.verbose = False
import numpy as np


class Test_rans3p(object):

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
        pass

#     def test_ex1(self):
#         self.compare_name = "T8P2"
#         self.example_setting("T=8.0 vspaceOrder=2 onlySaveFinalSolution=True")

    def test_ex2(self):
        self.compare_name = "T1P1"
        self.example_setting("T=0.1 vspaceOrder=1 onlySaveFinalSolution=True")


    def example_setting(self, pre_setting):
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
                                        "proteus.tests.cylinder2D.conforming_rans3p"))
            nList.append(
                importlib.import_module("."+nModule,
                                        "proteus.tests.cylinder2D.conforming_rans3p"))
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
        my_so.name += "_rans3p_"+self.compare_name #save data with different filename
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(my_so,
                                               pList,
                                               nList,
                                               sList,
                                               opts)
        self.aux_names.append(ns.modelList[0].name)
        ns.calculateSolution(my_so.name)
        # COMPARE VS SAVED FILES #
        actual=tables.open_file( my_so.name + '.h5')
        expected_path = 'comparison_files/' + 'comparison_' + self.compare_name + '_u_t2.csv'
        #write comparison file
        #np.array(actual.root.u_t2).tofile(os.path.join(self._scriptdir, expected_path),sep=",")
        np.testing.assert_almost_equal(np.fromfile(os.path.join(self._scriptdir, expected_path),sep=","),np.array(actual.root.u_t2).flatten(),decimal=10)

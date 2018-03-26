"""Tests for 2d flow around a cylinder with a conforming mesh and rans3p"""
from proteus.iproteus import *
from proteus import Comm
from proteus import Context
import tables
import importlib

comm = Comm.get()
Profiling.logLevel = 7
Profiling.verbose = False
import numpy as np


class Test_HotStart_rans3p():

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


    def test_hotstart_p1(self):
        self.compare_name = "T01P1_hotstart"
        self.example_setting("T=0.1 vspaceOrder=1 onlySaveFinalSolution=True",h5_filename="solution_p1")
        self.example_setting("T=0.1 vspaceOrder=1 onlySaveFinalSolution=True isHotStart=True", h5_filename="solution_p1", check_result=True, isHotstart=True,hotstart_t=0.1)

    def test_hotstart_p2(self):
        self.compare_name = "T01P2_hotstart"
        self.example_setting("T=0.1 vspaceOrder=2 onlySaveFinalSolution=True",h5_filename="solution_p2")
        self.example_setting("T=0.1 vspaceOrder=2 onlySaveFinalSolution=True isHotStart=True", h5_filename="solution_p2", check_result=True, isHotstart=True,hotstart_t=0.1)


    def example_setting(self, pre_setting, h5_filename, check_result=False, isHotstart=False, hotstart_t=0.0):
        Context.contextOptionsString = pre_setting
        from . import NS_hotstart_so as my_so
        reload(my_so)
        # defined in iproteus
        opts.profile = False
        opts.gatherArchive = True
        opts.hotStart = isHotstart
        opts.hotStartTime = hotstart_t
        
        pList=[]
        nList=[]
        sList=[]
        for (pModule,nModule) in my_so.pnList:
            pList.append(
                importlib.import_module("."+pModule,
                                        "proteus.tests.HotStart_3P"))
            nList.append(
                importlib.import_module("."+nModule,
                                        "proteus.tests.HotStart_3P"))
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
        my_so.name = h5_filename#"_hotstart_"+self.compare_name #save data with different filename
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(my_so,
                                               pList,
                                               nList,
                                               sList,
                                               opts)
        self.aux_names.append(ns.modelList[0].name)
        ns.calculateSolution(my_so.name)
        if check_result:
            # COMPARE VS SAVED FILES #
            expected_path = 'comparison_files/' + self.compare_name + '.h5'
            with tables.open_file(os.path.join(self._scriptdir, expected_path)) as expected, \
                    tables.open_file( my_so.name + '.h5') as actual:
                assert np.allclose(expected.root.u_t2,
                                   actual.root.u_t2,
                                   atol=1e-10)

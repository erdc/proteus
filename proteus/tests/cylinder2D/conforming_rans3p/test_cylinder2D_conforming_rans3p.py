from proteus.iproteus import *
from proteus import Comm
import tables
comm = Comm.get()
Profiling.logLevel = 7
Profiling.verbose = False
import numpy as np
import importlib

# from . import cylinder
# from . import cylinder_so as my_so
# from . import twp_navier_stokes_n as ns_n
# from . import twp_navier_stokes_p as ns_p
# from . import pressureincrement_n as pInc_n
# from . import pressureincrement_p as pInc_p
# from . import pressureInitial_n as pInit_n
# from . import pressureInitial_p as pInit_p
# from . import pressure_n as p_n
# from . import pressure_p as p_p

class TestIT():

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
        self.example_setting("T=1.0 vspaceOrder=1 onlySaveFinalSolution=True")


    def example_setting(self, pre_setting):
        Context.contextOptionsString = pre_setting
        cylinder = importlib.import_module('cylinder')
        my_so = importlib.import_module('cylinder_so')
        reload(cylinder)
        reload(my_so)
        # defined in iproteus
        opts.profile = False
        opts.gatherArchive = True
        
        pList=[]
        nList=[]
        sList=[]
        for (pModule,nModule) in my_so.pnList:
            pList.append(importlib.import_module(pModule))
            nList.append(importlib.import_module(nModule))
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
        expected_path = 'comparison_files/' + self.compare_name + '.h5'
        with tables.open_file(os.path.join(self._scriptdir, expected_path)) as expected, \
                tables.open_file(os.path.join(self._scriptdir + '/' + my_so.name + '.h5')) as actual:
            assert np.allclose(expected.root.u_t2,
                               actual.root.u_t2,
                               atol=1e-10)

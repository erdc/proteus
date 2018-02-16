from proteus.iproteus import *
from proteus import Comm
from proteus import Context
import tables


comm = Comm.get()
Profiling.logLevel = 7
Profiling.verbose = False
import numpy as np


class TestIT():

    @classmethod
    def setup_class(cls):
        pass

    @classmethod
    def teardown_class(cls):
        pass

    def setup_method(self, method):
        """Initialize the test problem. """
        self.aux_names = []
        self.meshdir = os.path.dirname(os.path.abspath(__file__))
        self._scriptdir = os.path.dirname(os.path.abspath(__file__))

    def teardown_method(self, method):
        return
        filenames = []
        for aux_name in self.aux_names:
            filenames.extend([aux_name + '.' + ext for ext in ['h5', 'xmf']])
        filenames.append('proteus.log')
        for f in filenames:
            if os.path.exists(f):
                try:
                    os.remove(f)
                except OSError, e:
                    print ("Error: %s - %s" % (e.filename, e.strerror))
            else:
                pass

#     def test_ex1(self):
#         self.compare_name = "T8P2"
#         self.example_setting("T=8.0 vspaceOrder=2 onlySaveFinalSolution=True")

    def test_ex2(self):
        self.compare_name = "T1_rans2p"
        self.example_setting("T=1.0 onlySaveFinalSolution=True")


    def example_setting(self, pre_setting):
        Context.contextOptionsString = pre_setting

        import cylinder2d as cylinder
        import cylinder_so as my_so
        reload(cylinder)  # Serious error
        reload(my_so)  # Serious error

        # defined in iproteus
        opts.profile = False
        opts.gatherArchive = True
        
        pList=[]
        nList=[]
        sList=[]

        for (pModule,nModule) in my_so.pnList:
            pList.append(__import__(pModule))
            if pList[-1].name == None:
                pList[-1].name = pModule
            nList.append(__import__(nModule))
            reload(pList[-1])  # Serious error
            reload(nList[-1])

        if my_so.sList == []:
            for i in range(len(my_so.pnList)):
                s = default_s
                sList.append(s)
        else:
            sList = my_so.sList
        my_so.name += "_rans2p_"+self.compare_name #save data with different filename
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

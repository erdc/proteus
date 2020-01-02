"""Tests for 2d flow around a cylinder with shifted boundary method on a 3D mesh"""
from builtins import range
from builtins import object
from proteus.iproteus import *
from proteus import Comm, defaults
from proteus import Context
import tables
import importlib


comm = Comm.get()
Profiling.logLevel = 7
Profiling.verbose = False
import numpy as np

from proteus.tests import Norms

L2_norm_cylinder_u_baseline=0.267723275157835

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
        FileList = ['',
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
        Context.contextOptionsString = pre_setting

        my_so = defaults.load_system("cylinder_so",modulepath)

        opts.profile = False
        opts.gatherArchive = True
        
        pList=[]
        nList=[]
        sList=[]
        for (pModule,nModule) in my_so.pnList:
            pList.append(defaults.load_physics(pModule,modulepath))
            nList.append(defaults.load_numerics(nModule,modulepath))
            if pList[-1].name == None:
                pList[-1].name = pModule

        if my_so.sList == []:
            for i in range(len(my_so.pnList)):
                s = default_s
                sList.append(s)
        else:
            sList = my_so.sList

        my_so.name += "_sbm_mesh3D_"+self.compare_name #save data with different filename
        ns = proteus.NumericalSolution.NS_base(my_so,
                                               pList,
                                               nList,
                                               sList,
                                               opts)
        ns.calculateSolution(my_so.name)

        actual = tables.open_file('cylinder_sbm_mesh3D_T001_P1_sbm_3Dmesh'+'.h5','r')
        L2_norm_u = Norms.get_L2_norm(actual,actual.root.u_t2)
        #print("%.15f" % L2_norm_u); import sys; sys.exit()
        np.testing.assert_almost_equal(L2_norm_u,L2_norm_cylinder_u_baseline)
        actual.close()


#!/usr/bin/env python
"""
Test level set solver
"""

from proteus.iproteus import *
from proteus import Comm
from proteus import Context
import tables


comm = Comm.get()
Profiling.logLevel = 2
Profiling.verbose = False
import numpy as np


class TestLS():

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

    def test_ex1(self):
        self.example_setting("nLevels=1 timeIntegration_ls='be'")

    def test_ex2(self):
        self.example_setting("nLevels=3 timeIntegration_ls='be'")

    def example_setting(self, pre_setting):

        Context.contextOptionsString = pre_setting

        import burgers2D
        import vof_burgers_2d_so as my_so
        reload(burgers2D)  # Serious error
        reload(my_so)  # Serious error

        # defined in iproteus
        opts.profile = False
        opts.gatherArchive = True

        simFlagsList = [{}]
        simFlagsList[0]['errorQuantities'] = ['u']
        # compute error in soln and glob. mass bal
        simFlagsList[0]['errorTypes'] = ['numericalSolution']
        # compute L2 norm in space or H0 or ...
        simFlagsList[0]['errorNorms'] = ['L1', 'L2']
        simFlagsList[0]['errorTimes'] = ['Last']  # 'All', 'Last'
        simFlagsList[0]['echo'] = True
        simFlagsList[0]['storeTimes'] = []
        simFlagsList[0]['storeQuantities'] = ['meshsize', 'errorData']
        simFlagsList[0]['dataDir'] = '.'
        simFlagsList[0]['dataFile'] = my_so.ct.datafile

        pList = [__import__(my_so.pnList[0][0])]
        nList = [__import__(my_so.pnList[0][1])]
        reload(pList[0])  # Serious error
        reload(nList[0])

        sList = [default_s]
        so = default_so
        so.name = my_so.name
        so.tnList = my_so.tnList

        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(so,
                                               pList,
                                               nList,
                                               sList,
                                               opts,
                                               simFlagsList)
        self.aux_names.append(ns.modelList[0].name)
        ns.calculateSolution(my_so.soname)
        ns.postStep(ns.modelList[0])
        # COMPARE VS SAVED FILES #
        expected_path = 'comparison_files/' + my_so.name + '.h5'
        with tables.open_file(os.path.join(self._scriptdir, expected_path)) as expected, \
                tables.open_file(os.path.join(self._scriptdir + '/' + my_so.name + '.h5')) as actual:
            assert np.allclose(expected.root.u_t2,
                               actual.root.u_t2,
                               atol=1e-10)

        import get_convergence as comp
        (error, convergence_rate) = comp.get_convergence_rate(my_so.ct.datafile)
        print my_so.ct.datafile
        print "\n error: ", error
        print "\n convergence: ", convergence_rate

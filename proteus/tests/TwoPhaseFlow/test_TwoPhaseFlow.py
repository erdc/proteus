#!/usr/bin/env python
"""
Test module for TwoPhaseFlow
"""
import os
import pytest
import tables
import numpy as np
import proteus.defaults
from proteus import Context
from proteus import default_so
from proteus.iproteus import *
import sys

Profiling.logLevel=1
Profiling.verbose=True

class TestTwoPhaseFlow(object):

    def setup_method(self,method):
        self._scriptdir = os.path.dirname(__file__)

    def runProblem(self,name,module,ns_model):
        #sys.argv.append('-f')
        #sys.argv.append(name+'.py')

        # paths
        path_models = proteus.__path__[0]+"/TwoPhaseFlow/models/"
        #path_utils = proteus.__path__[0]+"/TwoPhaseFlow/utils/"
        
        # Load so file
        #so_name = "TwoPhaseFlow_so.py"
        #so = proteus.defaults.load_system(so_name[:-3],path=path_utils)
        so = proteus.defaults.System_base()
        so.name = name
        
        # Read context 
        case = __import__(name)
        Context.setFromModule(case)
        ct = Context.get()

        # Create pnList
        if ns_model==0: #rans2p
            pnList = [("rans2p_p", "rans2p_n"),
                      ("clsvof_p", "clsvof_n")]
        else: #rans3p
            pnList = [("clsvof_p", "clsvof_n"),#0
                      ("rans3p_p", "rans3p_n"),#1
                      ("pressureincrement_p", "pressureincrement_n"),#2
                      ("pressure_p", "pressure_n"),#3
                      ("pressureInitial_p", "pressureInitial_n")]#4

        # System step controller 
        if ns_model==1: #rans3p
            PINIT_model=4
            modelSpinUpList = [PINIT_model]
            class Sequential_MinAdaptiveModelStepPS(default_so.Sequential_MinAdaptiveModelStep):
                def __init__(self,modelList,system=default_so.defaultSystem,stepExact=True):
                    default_so.Sequential_MinAdaptiveModelStep.__init__(self,
                                                                        modelList,
                                                                        system,
                                                                        stepExact)
                    self.modelList = modelList[:len(pnList)-1]
            #
            so.systemStepControllerType = Sequential_MinAdaptiveModelStepPS
        else:
            so.systemStepControllerType = default_so.Sequential_MinAdaptiveModelStep
        # create tnList
        so.tnList=[0,0.001,0.1]

        # Erase name from sys.argv
        #sys.argv = sys.argv[:-2]
        
        pList=[]
        nList=[]
        sList=[]
        
        # p and n lists
        for (pModule,nModule) in pnList:
            pList.append(proteus.defaults.load_physics(pModule,path=path_models))
            if pList[-1].name == None:
                pList[-1].name = pModule
            nList.append(proteus.defaults.load_numerics(nModule,path=path_models))
        # s list
        if sList == []:
            for i in range(len(pnList)):
                s = default_s
                sList.append(s)
        else:
            sList = sList
        #
        # Calculate solution
        ns = proteus.NumericalSolution.NS_base(so,pList,nList,sList,opts)
        ns.calculateSolution(name)

        # COMPARE VS SAVED FILES #
        expected_path = 'comparison_files/' + name + '.h5'
        expected = tables.open_file(os.path.join(self._scriptdir,expected_path))
        actual = tables.open_file(name + '.h5','r')
        assert np.allclose(expected.root.phi_t2,actual.root.phi_t2,atol=1e-8)
        assert np.allclose(expected.root.velocity_t2,actual.root.velocity_t2,atol=1e-8)
        expected.close()
        actual.close()
        
    def test_fallingBubble(self):
        from . import fallingBubble
        self.runProblem("fallingBubble",fallingBubble,ns_model=0)

    def test_marin(self):
        from . import marin
        self.runProblem("marin",marin,ns_model=0)

    def test_quiescentTank(self):
        from . import quiescentTank
        self.runProblem("quiescentTank",quiescentTank,ns_model=0)

    def test_risingBubble(self):
        from . import risingBubble
        self.runProblem("risingBubble",risingBubble,ns_model=1)

    def test_TwoDimBucklingFlow(self):
        from . import TwoDimBucklingFlow
        self.runProblem("TwoDimBucklingFlow",TwoDimBucklingFlow,ns_model=1)

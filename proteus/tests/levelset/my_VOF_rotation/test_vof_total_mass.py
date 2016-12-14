#!/usr/bin/env python
"""
Test module for VOF rotation with EV
"""
import os,sys,inspect

cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
if cmd_folder not in sys.path:
    sys.path.insert(0,cmd_folder)

cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() )) [0],"import_modules")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0,cmd_subfolder)

from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
Profiling.logLevel=2
Profiling.verbose=True
import numpy as np
import vof_rotation_2d_test_template as vf

class TestVOFrotationEV():

    @classmethod
    def setup_class(cls):
        pass

    @classmethod
    def teardown_class(cls):
        pass

    def setup_method(self,method):
        """Initialize the test problem. """
        reload(vf)
        self.sim_names = []
        self.aux_names = []


    def teardown_method(self,method):
        filenames = []
        for sim_name in self.sim_names:
            filenames.extend([sim_name+'.'+post for post in ['xmf','h5']])
        for aux_name in self.aux_names:
            filenames.extend(aux_name)

        for f in filenames:
            if os.path.exists(f):
                try:
                    os.remove(f)
                except OSError, e:
                    print ("Error: %s - %s." %(e.filename, e.strerror ))
            else:
                pass
    def test_vof_total_mass_T1m1(self):

        run_dir = os.path.dirname(os.path.abspath(__file__))
        ref_total_mass = np.loadtxt(os.path.join(run_dir,'comparison_files','total_mass_comp_0_rotation_c0p1_SSP33_level_3_vof_T1m1.txt'))
        #set the time step
        vf.p.T = 0.1
        vf.n.DT = vf.p.T/float(vf.n.nDTout)
        vf.so.DT = vf.n.DT
        vf.so.tnList = [i*vf.n.DT for i  in range(vf.n.nDTout+1)]

        ns = proteus.NumericalSolution.NS_base(vf.so,[vf.p],[vf.n],vf.so.sList,opts)
        sim_name = ns.modelList[0].name
        aux = ns.auxiliaryVariables[ns.modelList[0].name][0]
        self.sim_names.append(sim_name)
        self.aux_names.append(aux.ofile.name)
        ns.calculateSolution('test_vof_total_mass_T1')
        aux.ofile.close() #have to close manually for now, would be good to have a hook for this
        sim_total_mass = np.loadtxt('total_mass_comp_0_rotation_c0p1_SSP33_level_3_vof.txt')

        failed = np.allclose(ref_total_mass, sim_total_mass,
                             rtol=1e-05, atol=1e-07, equal_nan=True)
        #return ref_total_mass,sim_total_mass

    def test_vof_total_mass_T1(self):
        run_dir = os.path.dirname(os.path.abspath(__file__))
        ref_total_mass = np.loadtxt(os.path.join(run_dir,'comparison_files','total_mass_comp_0_rotation_c0p1_SSP33_level_3_vof.txt'))
        #set the time step
        vf.p.T = 1.0
        vf.n.DT = vf.p.T/float(vf.n.nDTout)
        vf.so.DT = vf.n.DT
        vf.so.tnList = [i*vf.n.DT for i  in range(vf.n.nDTout+1)]

        ns = proteus.NumericalSolution.NS_base(vf.so,[vf.p],[vf.n],vf.so.sList,opts)
        sim_name = ns.modelList[0].name
        aux = ns.auxiliaryVariables[ns.modelList[0].name][0]
        self.sim_names.append(sim_name)
        self.aux_names.append(aux.ofile.name)

        ns.calculateSolution('test_vof_total_mass_T1')
        aux.ofile.close() #have to close manually for now, would be good to have a hook for this
        sim_total_mass = np.loadtxt('total_mass_comp_0_rotation_c0p1_SSP33_level_3_vof.txt')

        failed = np.allclose(ref_total_mass, sim_total_mass,
                             rtol=1e-05, atol=1e-07, equal_nan=True)
        #return ref_total_mass,sim_total_mass

if __name__ == '__main__':
    pass
    # reload(vf)
    # ref,sim = TestVOFrotationEV().test_vof_total_mass_T1()
    # print ref.shape
    # print sim.shape
    # ref1 = np.loadtxt('total_mass_save_T01.txt')
    # print ref1.shape
    # np.max(ref1[:,2]-sim[:,2])

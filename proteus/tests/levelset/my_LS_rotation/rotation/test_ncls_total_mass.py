#!/usr/bin/env python
"""
Test module for NCLS rotation with EV
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
import ncls_rotation_2d_test_template as ncls


class TestNCLSrotationEV():

    @classmethod
    def setup_class(cls):
        pass

    @classmethod
    def teardown_class(cls):
        pass

    def setup_method(self,method):
        """Initialize the test problem. """
        reload(ncls)
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

    def test_ncls_T1_SSP33_ElementBasedEV(self):
        """
        This is a regression test. 
        We check the total mass, int(u^2), int(|u|), max(u) at every time step against some reference. 
        The following method is considered 
        2nd order MPP edge based EV stabilization with SSP33
        These are the flags used for NCLS.h in this benchmark 
        #define FIX_BOUNDARY_KUZMINS 1
        #define QUANTITIES_OF_INTEREST 0
        """        
        run_dir = os.path.dirname(os.path.abspath(__file__))
        
        #set the time step
        ncls.p.T = 0.1
        ncls.n.nDTout = 10
        ncls.n.DT = ncls.p.T/float(ncls.n.nDTout)
        ncls.so.DT = ncls.n.DT
        ncls.so.tnList = [i*ncls.n.DT for i  in range(ncls.n.nDTout+1)]

        #force SSP33
        ncls.timeIntegration_ncls="SSP33"
        ncls.n.timeOrder = 3
        ncls.n.nStagesTime = 3
        ncls.rot2D.soname=ncls.rot2D.soname.replace("FE","SSP33")
        ncls.p.name = ncls.p.name.replace("FE","SSP33")
        ncls.so.name = ncls.rot2D.soname

        #############################################
        # SECOND ORDER NON-MPP ELEMENT BASED METHOD #
        #############################################
        """ 
        2nd ORDER KUZMINS METHOD (MPP via FCT)
        """
        
        ncls.p.coefficients = ncls.rot2D.MyCoefficients(epsFact=ncls.rot2D.epsFactHeaviside,checkMass=ncls.rot2D.checkMass,useMetrics=ncls.rot2D.useMetrics,ME_model=0,
                                                        EDGE_VISCOSITY=1,
                                                        ENTROPY_VISCOSITY=1,
                                                        POWER_SMOOTHNESS_INDICATOR=1,
                                                        LUMPED_MASS_MATRIX=0,
                                                        FCT=1)
        
        ns = proteus.NumericalSolution.NS_base(ncls.so,[ncls.p],[ncls.n],ncls.so.sList,opts)
        sim_name = ns.modelList[0].name
        # READ REFERENCE 
        ref_total_mass = np.loadtxt(os.path.join(run_dir,'comparison_files','total_mass_comp_0_T1m1_rotation_c0p1_SSP33_EdgeBasedEV_level_3_ncls.txt'))
        # end of reading reference 

        aux = ns.auxiliaryVariables[ns.modelList[0].name][0]
        self.sim_names.append(sim_name)
        self.aux_names.append(aux.ofile.name)
        ns.calculateSolution('test_ncls_total_mass_T1')
        aux.ofile.close() #have to close manually for now, would be good to have a hook for this

        sim_total_mass = np.loadtxt('total_mass_comp_0_'+sim_name+'.txt')

        print "****************************************************************"
        print "*****... SECOND ORDER NON-MPP ELEMENT BASED METHOD ...**********"
        print "****************************************************************"
        # Regression test: check total_mass, int(u^2), int(|u|) and max(u)

        print ref_total_mass
        print sim_total_mass
        failed = np.allclose(ref_total_mass, sim_total_mass,
                             rtol=1e-05, atol=1e-07, equal_nan=True)
        print(failed)
        assert(failed)

if __name__ == '__main__':
    pass

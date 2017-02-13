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

    def test_vof_total_mass_FE_EdgeBasedEV(self):
        """
        Test total mass for Forward Euler Integration running for any final time. 
        The following (4) methods are considered: 
           1. 1st order MPP KUZMIN's method with boundary treatment 
           2. 2nd order MPP ...
           3. 2nd order non-MPP KUZMIN's with entropy viscosity, consistent mass matrix and art compression 
           4. 2nd order MPP KUZMIN's via FCT with entropy viscosity, consistent mass matrix and art comrpession
        This test check some correctness of different components of the full algorithm 
        These are the flags used for VOF.h in this benchmark 
        #define KUZMINS_METHOD 1
        #define INTEGRATE_BY_PARTS 1
        #define QUANTITIES_OF_INTEREST 0
        #define FIX_BOUNDARY_KUZMINS 1        
        """
        run_dir = os.path.dirname(os.path.abspath(__file__))

        #set the time step
        vf.p.T = 1.0
        vf.n.nDTout = 10
        vf.n.DT = vf.p.T/float(vf.n.nDTout)
        vf.so.DT = vf.n.DT
        vf.so.tnList = [i*vf.n.DT for i  in range(vf.n.nDTout+1)]

        #force F2orwardEuler
        vf.timeIntegration_vof="FE"
        vf.n.timeOrder = 1
        vf.n.nStagesTime = 1
        vf.rot2D.soname=vf.rot2D.soname.replace("SSP33","FE")
        vf.p.name = vf.p.name.replace("SSP33","FE")
        vf.so.name = vf.rot2D.soname

        #################################
        # 1. FIRST ORDER KUZMINS METHOD #
        #################################
        """ 
        1st ORDER KUZMINS METHOD 
        """
        
        vf.p.coefficients = vf.rot2D.MyCoefficients(epsFact=vf.rot2D.epsFactHeaviside,checkMass=vf.rot2D.checkMass,useMetrics=vf.rot2D.useMetrics,ME_model=0,
                                                    EDGE_VISCOSITY=1,
                                                    ENTROPY_VISCOSITY=0,
                                                    POWER_SMOOTHNESS_INDICATOR=0,
                                                    LUMPED_MASS_MATRIX=1,
                                                    FCT=0,
                                                    cK=0)

        ns = proteus.NumericalSolution.NS_base(vf.so,[vf.p],[vf.n],vf.so.sList,opts)
        sim_name = ns.modelList[0].name
        # READ REFERENCE 
        init_mass = np.loadtxt(os.path.join(run_dir,'comparison_files','init_mass.txt'))
        # end of loading reference

        aux = ns.auxiliaryVariables[ns.modelList[0].name][0]
        self.sim_names.append(sim_name)
        sim_total_mass = np.loadtxt('total_mass_comp_0_'+sim_name+'.txt')
        self.aux_names.append(aux.ofile.name)
        ns.calculateSolution('test_vof_total_mass')
        aux.ofile.close() #have to close manually for now, would be good to have a hook for this
        sim_total_mass = np.loadtxt('total_mass_comp_0_'+sim_name+'.txt')
        
        # Check initial vs final mass
        print "*************************************************"
        print "*****... FIRST ORDER KUZMINS METHOD ...**********"
        print "*************************************************"
        final_mass = sim_total_mass[:,0][-1]        
        print init_mass, final_mass
        failed = np.allclose(init_mass, final_mass, rtol=1e-05, atol=1e-07, equal_nan=True) 
        print failed
        assert(failed) 

        ########################################
        # 2. SECOND ORDER KUZMINS METHOD (MPP) #
        ########################################
        """ 
        2nd ORDER KUZMINS METHOD 
        """
        
        vf.p.coefficients = vf.rot2D.MyCoefficients(epsFact=vf.rot2D.epsFactHeaviside,checkMass=vf.rot2D.checkMass,useMetrics=vf.rot2D.useMetrics,ME_model=0,
                                                    EDGE_VISCOSITY=1,
                                                    ENTROPY_VISCOSITY=0,
                                                    POWER_SMOOTHNESS_INDICATOR=2,
                                                    LUMPED_MASS_MATRIX=1,
                                                    FCT=0,
                                                    cK=0)

        ns = proteus.NumericalSolution.NS_base(vf.so,[vf.p],[vf.n],vf.so.sList,opts)
        sim_name = ns.modelList[0].name
        
        aux = ns.auxiliaryVariables[ns.modelList[0].name][0]
        self.sim_names.append(sim_name)
        sim_total_mass = np.loadtxt('total_mass_comp_0_'+sim_name+'.txt')
        self.aux_names.append(aux.ofile.name)
        ns.calculateSolution('test_vof_total_mass')
        aux.ofile.close() #have to close manually for now, would be good to have a hook for this
        sim_total_mass = np.loadtxt('total_mass_comp_0_'+sim_name+'.txt')
        
        # Check initial vs final mass
        print "********************************************************"
        print "*****... SECOND ORDER KUZMINS METHOD (MPP) ...**********"
        print "********************************************************"
        final_mass = sim_total_mass[:,0][-1]
        print init_mass, final_mass
        failed = np.allclose(init_mass, final_mass, rtol=1e-05, atol=1e-07, equal_nan=True) 
        print failed
        assert(failed) 

        ############################################
        # 3. SECOND ORDER KUZMINS METHOD (Non-MPP) #
        ############################################
        """ 
        2nd ORDER KUZMINS METHOD (Non-MPP)
        """
        
        vf.p.coefficients = vf.rot2D.MyCoefficients(epsFact=vf.rot2D.epsFactHeaviside,checkMass=vf.rot2D.checkMass,useMetrics=vf.rot2D.useMetrics,ME_model=0,
                                                    EDGE_VISCOSITY=1,
                                                    ENTROPY_VISCOSITY=1, #NOTE: ENTROPY VISCOSITY IS ACTIVATED
                                                    POWER_SMOOTHNESS_INDICATOR=2,
                                                    LUMPED_MASS_MATRIX=0, # NOTE: CONSISTENT MASS MATRIX
                                                    FCT=0,
                                                    cK=0.25) #NOTE: ARTIFICIAL COMPRESSION 

        ns = proteus.NumericalSolution.NS_base(vf.so,[vf.p],[vf.n],vf.so.sList,opts)
        sim_name = ns.modelList[0].name

        aux = ns.auxiliaryVariables[ns.modelList[0].name][0]
        self.sim_names.append(sim_name)
        sim_total_mass = np.loadtxt('total_mass_comp_0_'+sim_name+'.txt')
        self.aux_names.append(aux.ofile.name)
        ns.calculateSolution('test_vof_total_mass')
        aux.ofile.close() #have to close manually for now, would be good to have a hook for this
        sim_total_mass = np.loadtxt('total_mass_comp_0_'+sim_name+'.txt')
        
        # Check initial vs final mass
        print "************************************************************"
        print "*****... SECOND ORDER KUZMINS METHOD (Non-MPP) ...**********"
        print "************************************************************"
        final_mass = sim_total_mass[:,0][-1]
        print init_mass, final_mass
        failed = np.allclose(init_mass, final_mass, rtol=1e-05, atol=1e-07, equal_nan=True) 
        print failed
        assert(failed) 

        ################################################
        # 4. SECOND ORDER KUZMINS METHOD (MPP via FCT) #
        ################################################
        """ 
        2nd ORDER KUZMINS METHOD (MPP via FCT)
        """
        
        vf.p.coefficients = vf.rot2D.MyCoefficients(epsFact=vf.rot2D.epsFactHeaviside,checkMass=vf.rot2D.checkMass,useMetrics=vf.rot2D.useMetrics,ME_model=0,
                                                    EDGE_VISCOSITY=1,
                                                    ENTROPY_VISCOSITY=1,
                                                    POWER_SMOOTHNESS_INDICATOR=2,
                                                    LUMPED_MASS_MATRIX=0,
                                                    FCT=1, #NOTE: FCT IS USED
                                                    cK=0.25)

        ns = proteus.NumericalSolution.NS_base(vf.so,[vf.p],[vf.n],vf.so.sList,opts)
        sim_name = ns.modelList[0].name

        aux = ns.auxiliaryVariables[ns.modelList[0].name][0]
        self.sim_names.append(sim_name)
        sim_total_mass = np.loadtxt('total_mass_comp_0_'+sim_name+'.txt')
        self.aux_names.append(aux.ofile.name)
        ns.calculateSolution('test_vof_total_mass')
        aux.ofile.close() #have to close manually for now, would be good to have a hook for this
        sim_total_mass = np.loadtxt('total_mass_comp_0_'+sim_name+'.txt')
        
        # Check initial vs final mass
        print "****************************************************************"
        print "*****... SECOND ORDER KUZMINS METHOD (MPP via FCT) ...**********"
        print "****************************************************************"
        final_mass = sim_total_mass[:,0][-1]
        print init_mass, final_mass
        failed = np.allclose(init_mass, final_mass, rtol=1e-05, atol=1e-07, equal_nan=True) 
        print failed
        assert(failed) 

    def test_vof_total_mass_SSP33_EdgeBasedEV(self):
        """
        Test total mass for SSP33 Integration running for any final time. 
        The following method is considered
           2nd order MPP KUZMIN's via FCT with entropy viscosity, consistent mass matrix and art comrpession
        This test check some correctness of different components of the full algorithm 
        These are the flags used for VOF.h in this benchmark 
        #define KUZMINS_METHOD 1
        #define INTEGRATE_BY_PARTS 1
        #define QUANTITIES_OF_INTEREST 0
        #define FIX_BOUNDARY_KUZMINS 1        
        """
        run_dir = os.path.dirname(os.path.abspath(__file__))

        #set the time step
        vf.p.T = 1.0
        vf.n.nDTout = 10
        vf.n.DT = vf.p.T/float(vf.n.nDTout)
        vf.so.DT = vf.n.DT
        vf.so.tnList = [i*vf.n.DT for i  in range(vf.n.nDTout+1)]

        #force SSP33
        vf.timeIntegration_vof="SSP33"
        vf.n.timeOrder = 3
        vf.n.nStagesTime = 3
        vf.rot2D.soname=vf.rot2D.soname.replace("FE","SSP33")
        vf.p.name = vf.p.name.replace("FE","SSP33")
        vf.so.name = vf.rot2D.soname

        #############################################
        # SECOND ORDER KUZMINS METHOD (MPP via FCT) #
        #############################################
        """ 
        2nd ORDER KUZMINS METHOD (MPP via FCT)
        """
        
        vf.p.coefficients = vf.rot2D.MyCoefficients(epsFact=vf.rot2D.epsFactHeaviside,checkMass=vf.rot2D.checkMass,useMetrics=vf.rot2D.useMetrics,ME_model=0,
                                                    EDGE_VISCOSITY=1,
                                                    ENTROPY_VISCOSITY=1,
                                                    POWER_SMOOTHNESS_INDICATOR=2,
                                                    LUMPED_MASS_MATRIX=0,
                                                    FCT=1, #NOTE: FCT IS USED
                                                    cK=0.25)

        ns = proteus.NumericalSolution.NS_base(vf.so,[vf.p],[vf.n],vf.so.sList,opts)
        sim_name = ns.modelList[0].name
        # READ REFERENCE
        init_mass = np.loadtxt(os.path.join(run_dir,'comparison_files','init_mass.txt'))
        # end of loading reference

        aux = ns.auxiliaryVariables[ns.modelList[0].name][0]
        self.sim_names.append(sim_name)
        sim_total_mass = np.loadtxt('total_mass_comp_0_'+sim_name+'.txt')
        self.aux_names.append(aux.ofile.name)
        ns.calculateSolution('test_vof_total_mass')
        aux.ofile.close() #have to close manually for now, would be good to have a hook for this
        sim_total_mass = np.loadtxt('total_mass_comp_0_'+sim_name+'.txt')
        
        # Check initial vs final mass
        print "****************************************************************"
        print "*****... SECOND ORDER KUZMINS METHOD (MPP via FCT) ...**********"
        print "****************************************************************"
        final_mass = sim_total_mass[:,0][-1]
        print init_mass, final_mass
        failed = np.allclose(init_mass, final_mass, rtol=1e-05, atol=1e-07, equal_nan=True) 
        print failed
        assert(failed) 

    def test_vof_T1_SSP33_EdgeBasedEV(self):
        """
        This is a regression test. 
        We check the total mass, int(u^2), int(|u|), max(u) at every time step against some reference. 
        The following method is considered 
           2nd order MPP KUZMIN's via FCT with entropy viscosity, consistent mass matrix and art comrpession
        These are the flags used for VOF.h in this benchmark 
        #define KUZMINS_METHOD 1
        #define INTEGRATE_BY_PARTS 1
        #define QUANTITIES_OF_INTEREST 0
        #define FIX_BOUNDARY_KUZMINS 1
        """
        run_dir = os.path.dirname(os.path.abspath(__file__))
        
        #set the time step
        vf.p.T = 1.0
        vf.n.nDTout = 10
        vf.n.DT = vf.p.T/float(vf.n.nDTout)
        vf.so.DT = vf.n.DT
        vf.so.tnList = [i*vf.n.DT for i  in range(vf.n.nDTout+1)]

        #force SSP33
        vf.timeIntegration_vof="SSP33"
        vf.n.timeOrder = 3
        vf.n.nStagesTime = 3
        vf.rot2D.soname=vf.rot2D.soname.replace("FE","SSP33")
        vf.p.name = vf.p.name.replace("FE","SSP33")
        vf.so.name = vf.rot2D.soname

        #############################################
        # SECOND ORDER KUZMINS METHOD (MPP via FCT) #
        #############################################
        """ 
        2nd ORDER KUZMINS METHOD (MPP via FCT)
        """
        
        vf.p.coefficients = vf.rot2D.MyCoefficients(epsFact=vf.rot2D.epsFactHeaviside,checkMass=vf.rot2D.checkMass,useMetrics=vf.rot2D.useMetrics,ME_model=0,
                                                    EDGE_VISCOSITY=1,
                                                    ENTROPY_VISCOSITY=1,
                                                    POWER_SMOOTHNESS_INDICATOR=2,
                                                    LUMPED_MASS_MATRIX=0,
                                                    FCT=1, #NOTE: FCT IS USED
                                                    cK=0.25)
        
        ns = proteus.NumericalSolution.NS_base(vf.so,[vf.p],[vf.n],vf.so.sList,opts)
        sim_name = ns.modelList[0].name
        # READ REFERENCE 
        ref_total_mass = np.loadtxt(os.path.join(run_dir,'comparison_files','total_mass_comp_0_T1_rotation_c0p1_SSP33_EdgeBasedEV_level_3_vof.txt'))
        #ref_total_mass = np.loadtxt(os.path.join(run_dir,'comparison_files','total_mass_comp_0_T1_'+sim_name+'.txt'))
        # end of reading reference 

        aux = ns.auxiliaryVariables[ns.modelList[0].name][0]
        self.sim_names.append(sim_name)
        self.aux_names.append(aux.ofile.name)
        ns.calculateSolution('test_vof_total_mass_T1')
        aux.ofile.close() #have to close manually for now, would be good to have a hook for this

        sim_total_mass = np.loadtxt('total_mass_comp_0_'+sim_name+'.txt')

        failed = np.allclose(ref_total_mass, sim_total_mass,
                             rtol=1e-05, atol=1e-07, equal_nan=True)
        print(failed)
        assert(failed)

    def test_vof_T1_SSP33_ElementBasedEV(self):
        """
        This is a regression test. In addition, we check correctness of the method by checking the final mass.
        We check the total mass, int(u^2), int(|u|), max(u) at every time step against some reference. 
        The following method is considered 
           2nd order Non-MPP element based EV stabilization with artificial compression 
        These are the flags used for VOF.h in this benchmark 
        #define KUZMINS_METHOD 1
        #define INTEGRATE_BY_PARTS 1
        #define QUANTITIES_OF_INTEREST 0
        #define FIX_BOUNDARY_KUZMINS 1
        """        
        run_dir = os.path.dirname(os.path.abspath(__file__))
        
        #set the time step
        vf.p.T = 1.0
        vf.n.nDTout = 10
        vf.n.DT = vf.p.T/float(vf.n.nDTout)
        vf.so.DT = vf.n.DT
        vf.so.tnList = [i*vf.n.DT for i  in range(vf.n.nDTout+1)]

        #force SSP33
        vf.timeIntegration_vof="SSP33"
        vf.n.timeOrder = 3
        vf.n.nStagesTime = 3
        vf.rot2D.soname=vf.rot2D.soname.replace("FE","SSP33")
        vf.p.name = vf.p.name.replace("FE","SSP33")
        vf.so.name = vf.rot2D.soname

        #############################################
        # SECOND ORDER NON-MPP ELEMENT BASED METHOD #
        #############################################
        """ 
        2nd ORDER KUZMINS METHOD (MPP via FCT)
        """
        
        vf.p.coefficients = vf.rot2D.MyCoefficients(epsFact=vf.rot2D.epsFactHeaviside,checkMass=vf.rot2D.checkMass,useMetrics=vf.rot2D.useMetrics,ME_model=0,
                                                    EDGE_VISCOSITY=0,
                                                    ENTROPY_VISCOSITY=1,
                                                    POWER_SMOOTHNESS_INDICATOR=2, #This is irrelevant
                                                    LUMPED_MASS_MATRIX=0,
                                                    FCT=0, #NOTE: NO FCT IS USED
                                                    cK=0.25, 
                                                    cE=1.0, #NOTE: For the element based EV these constants are important
                                                    cMax=0.1)
        
        ns = proteus.NumericalSolution.NS_base(vf.so,[vf.p],[vf.n],vf.so.sList,opts)
        sim_name = ns.modelList[0].name
        # READ REFERENCE 
        ref_total_mass = np.loadtxt(os.path.join(run_dir,'comparison_files','total_mass_comp_0_T1_rotation_c0p1_SSP33_ElementBasedEV_level_3_vof.txt'))
        # end of reading reference 

        aux = ns.auxiliaryVariables[ns.modelList[0].name][0]
        self.sim_names.append(sim_name)
        self.aux_names.append(aux.ofile.name)
        ns.calculateSolution('test_vof_total_mass_T1')
        aux.ofile.close() #have to close manually for now, would be good to have a hook for this

        sim_total_mass = np.loadtxt('total_mass_comp_0_'+sim_name+'.txt')

        print "****************************************************************"
        print "*****... SECOND ORDER NON-MPP ELEMENT BASED METHOD ...**********"
        print "****************************************************************"
        # Regression test: check total_mass, int(u^2), int(|u|) and max(u)
        failed = np.allclose(ref_total_mass, sim_total_mass,
                             rtol=1e-05, atol=1e-07, equal_nan=True)
        print(failed)
        assert(failed)

        # Correctness: check final mass vs reference initial mass
        init_mass = np.loadtxt(os.path.join(run_dir,'comparison_files','init_mass.txt'))
        final_mass = sim_total_mass[:,0][-1]
        print init_mass, final_mass
        failed = np.allclose(init_mass, final_mass, rtol=1e-05, atol=1e-07, equal_nan=True) 
        print failed
        assert(failed) 

if __name__ == '__main__':
    pass

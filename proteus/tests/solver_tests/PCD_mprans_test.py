#!/usr/bin/env python
"""

Test module for the convection-diffusion operator.

"""
import os,sys,inspect
cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
if cmd_folder not in sys.path:
    sys.path.insert(0,cmd_folder)

cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],
                                                              "import_modules")))
cmd_subfolder_0 = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],
                                                              "comparison_files")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0,cmd_subfolder)

from proteus.iproteus import *
from proteus import Comm
from proteus import LinearAlgebraTools
comm = Comm.get()
#Profiling.logLevel=7
#Profiling.verbose=True
import numpy.testing as npt
from nose.tools import ok_ as ok
from nose.tools import eq_ as eq
from nose.tools import set_trace
from petsc4py import PETSc as p4pyPETSc
import pytest

from scipy.sparse import csr_matrix
import petsc4py
import numpy as np

class TestTempNSEDrivenCavity():
    """ This class runs a small NSE test problem """
    @classmethod
    def setup_class(self):
        import twp_navier_stokes_DC_2d_p
        import twp_navier_stokes_DC_2d_n
        pList = [twp_navier_stokes_DC_2d_p]
        nList = [twp_navier_stokes_DC_2d_n]    
        so = default_so
        so.tnList = [0.,1.]
        so.name = pList[0].name
        so.sList=[default_s]
        opts.verbose=True
        opts.profile=True
        self._setPETSc()
        self.ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
    
    @classmethod
    def teardown_class(self):
        """Tear down function """
        FileList = ['mprans_test.xmf',
                    'mprans_test.h5',
                    'Cp.m',
                    'Cp',
                    'Cp.info',
                    'rdomain.edge',
                    'rdomain.ele',
                    'rdomain.neig',
                    'rdomain.node',
                    'rdomain.poly',
                    'proteus.log']
        for file in FileList:
            if os.path.isfile(file):
                os.remove(file)

    @pytest.mark.skip(reason="WIP")
    def test_01_FullRun(self):
        self.ns.calculateSolution('test_nse')
        assert(0==1)

    @classmethod
    def _setPETSc(self):
        petsc4py.PETSc.Options().setValue("ksp_type","fgmres")
        petsc4py.PETSc.Options().setValue("ksp_atol",1e-20)
        petsc4py.PETSc.Options().setValue("ksp_atol",1e-6)
        petsc4py.PETSc.Options().setValue("pc_fieldsplit_type","schur")
        petsc4py.PETSc.Options().setValue("pc_fieldsplit_schur_fact_type","full")
        petsc4py.PETSc.Options().setValue("fieldsplit_velocity_ksp_type","preonly")
        petsc4py.PETSc.Options().setValue("fieldsplit_velocity_pc_type","lu")
        petsc4py.PETSc.Options().setValue("fieldsplit_pressure_ksp_type","fgmres")
        petsc4py.PETSc.Options().setValue("fieldsplit_pressure_ksp_max_it",3)
        petsc4py.PETSc.Options().setValue("fieldsplit_pressure_ksp_atol",1e-2)
        petsc4py.PETSc.Options().setValue("fieldsplit_pressure_ksp_rtol",1e-2)

if __name__ == '__main__':
    from proteus import Comm
    comm = Comm.init()

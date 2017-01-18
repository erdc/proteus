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
Profiling.logLevel=7
Profiling.verbose=True
import numpy.testing as npt
# from nose.tools import ok_ as ok
# from nose.tools import eq_ as eq
# from nose.tools import set_trace
from petsc4py import PETSc as p4pyPETSc

from scipy.sparse import csr_matrix
import petsc4py
import numpy as np
import pytest

@pytest.mark.skip(reason="uknown regression - WIP")
def test_pcd_shell():
    '''  Tests the pcd_shell operators produce correct output. '''
    from proteus import LinearAlgebraTools
    vals_Qp = [1., 3., 1.5, -2.1, 4., 3.]
    col_idx_Qp = [0, 2, 0, 1, 1, 2]
    row_idx_Qp = [0, 2, 4, 6]
    vals_Fp = [3.2, 1.1, 6.3, 1., -5.1]
    col_idx_Fp = [0, 1, 0, 2, 0]
    row_idx_Fp = [0, 2, 4, 5]
    vals_Ap = [3., 1., 2., 7., 4., 3., 8., 1.1]
    col_idx_Ap = [0, 1, 2, 0, 1, 0, 1, 2]
    row_idx_Ap = [0, 3, 5, 8]
    size_n = len(row_idx_Qp) - 1
    petsc_matQp = p4pyPETSc.Mat().createAIJ(size = (size_n,size_n), 
                                            csr = (row_idx_Qp,col_idx_Qp,vals_Qp))
    petsc_matFp = p4pyPETSc.Mat().createAIJ(size = (size_n,size_n), 
                                            csr = (row_idx_Fp,col_idx_Fp,vals_Fp))
    petsc_matAp = p4pyPETSc.Mat().createAIJ(size = (size_n,size_n), 
                                            csr = (row_idx_Ap,col_idx_Ap,vals_Ap))
    PCD_shell = LinearAlgebraTools.PCDInv_shell(petsc_matQp,petsc_matFp,petsc_matAp)
    x_vec = np.ones(size_n)
    y_vec = np.zeros(size_n)
    x_PETSc_vec = p4pyPETSc.Vec().createWithArray(x_vec)
    y_PETSc_vec = p4pyPETSc.Vec().createWithArray(y_vec)
    A = None
    PCD_shell.apply(A,x_PETSc_vec,y_PETSc_vec)
    true_vec = np.mat('[0.52511998  -0.13008364  -0.03200969]')
    assert np.allclose(y_vec,true_vec)

def test_BAinvBt_shell():
    ''' Test the B_Ainv_Bt_shell '''
    from proteus import LinearAlgebraTools
    vals_A =    [5.5,3.6,1.8,7.1,7.9,2.1,1.0]
    col_idx_A = [0  ,1  ,2  ,1  ,2  ,0  ,2 ]
    row_idx_A = [0, 3, 5, 7]
    vals_B    =  [1.1, 6.3, 7.3, 3.6, 6.3, 2.4, 6.2, 1.2, 1.0]
    col_idx_B =  [0  , 2  , 0  , 1  , 2  , 1  , 2  , 0  , 2  ]
    row_idx_B =  [0, 2, 5, 7, 9]
    num_B_rows = len(row_idx_B)-1
    num_B_cols = len(row_idx_A)-1
    petsc_matA = p4pyPETSc.Mat().createAIJ(size = (num_B_cols,num_B_cols) ,
                                           csr = (row_idx_A, col_idx_A, vals_A))
    petsc_matB = p4pyPETSc.Mat().createAIJ(size = (num_B_rows,num_B_cols) ,
                                           csr = (row_idx_B, col_idx_B, vals_B))
    BinvABt_shell = LinearAlgebraTools.B_Ainv_Bt_shell(petsc_matA,petsc_matB)
    x_vec = np.ones(num_B_rows)
    y_vec = np.zeros(num_B_rows)
    x_PETSc_vec = p4pyPETSc.Vec().createWithArray(x_vec)
    y_PETSc_vec = p4pyPETSc.Vec().createWithArray(y_vec)
    A = None
    BinvABt_shell.mult(A,x_PETSc_vec,y_PETSc_vec)
    true_solu = np.mat('[64.60463468 60.7745635 35.16771019 15.33818394]')
    assert np.allclose(y_vec,true_solu)

def test_lsc_shell():
    ''' Test for the lsc operator shell '''
    vals_A =     [5.5,7.1,1.0]
    col_idx_A =  [0 , 1 , 2  ]
    row_idx_A =  [0, 1, 2, 3]
    vals_B    =  [1.1, 6.3, 7.3, 3.6, 6.3]
    col_idx_B =  [0  , 2  , 0  , 1  , 2  ]
    row_idx_B =  [0, 2, 5]
    vals_F  =    [3.2, 1.1, 6.3, 1., -5.1]
    col_idx_F  = [0, 1, 0, 2, 0]
    row_idx_F  = [0, 2, 4, 5]
    num_B_rows = len(row_idx_B) - 1
    num_B_cols = len(row_idx_A) - 1
    petsc_matA = p4pyPETSc.Mat().createAIJ(size = (num_B_cols,num_B_cols), 
                                            csr = (row_idx_A,col_idx_A,vals_A))
    petsc_matB = p4pyPETSc.Mat().createAIJ(size = (num_B_rows,num_B_cols), 
                                            csr = (row_idx_B,col_idx_B,vals_B))
    petsc_matF = p4pyPETSc.Mat().createAIJ(size = (num_B_cols,num_B_cols), 
                                            csr = (row_idx_F,col_idx_F,vals_F))
    LSC_shell = LinearAlgebraTools.LSCInv_shell(petsc_matA,petsc_matB,petsc_matF)
    x_vec = np.ones(num_B_rows)
    y_vec = np.zeros(num_B_rows)
    x_PETSc_vec = p4pyPETSc.Vec().createWithArray(x_vec)
    y_PETSc_vec = p4pyPETSc.Vec().createWithArray(y_vec)
    A = None
    LSC_shell.apply(A,x_PETSc_vec,y_PETSc_vec)
    true_solu = np.mat('[-0.01096996,0.00983216]')
    assert np.allclose(y_vec,true_solu)

class TestNSEDrivenCavity():
    """ This class runs a small NSE test problem """
    @classmethod
    def setup_class(self):
        import nseDrivenCavity_2d_p
        import nseDrivenCavity_2d_n
        pList = [nseDrivenCavity_2d_p]
        nList = [nseDrivenCavity_2d_n]    
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
        """ Tear down function """
        FileList = ['drivenCavityNSETrial.h5',
                    'drivenCavityNSETrial.xmf',
                    'Cp.m',
                    'Cp',
                    'Cp.info',
                    'proteus.log',
                    'proteus_default.log',
                    'rdomain.edge',
                    'rdomain.ele',
                    'rdomain.neig',
                    'rdomain.node',
                    'rdomain.poly']
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
    test_pcd_shell()
    test_BAinvBt_shell()
    test_lsc_shell()
#    test = TestNSEDrivenCavity()
#    test.setUp()
#    test.test_01_FullRun()
#    test.tearDown()

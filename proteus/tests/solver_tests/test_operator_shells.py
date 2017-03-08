#!/usr/bin/env python
"""

Test module for solver shell operator.

"""
import proteus.test_utils.TestTools
from proteus.iproteus import *
from proteus import Comm
from proteus import LinearAlgebraTools

import os
import sys
import inspect
import petsc4py
import numpy as np
import pytest

from petsc4py import PETSc as p4pyPETSc
from scipy.sparse import csr_matrix

proteus.test_utils.TestTools.addSubFolders( inspect.currentframe() )

class TestOperatorShells(proteus.test_utils.TestTools.BasicTest):

    def setup_method(self,method):
        self.petsc_options = p4pyPETSc.Options()

    def teardown_method(self,method):
        # ARB TODO - There must be a better way to do this...
        rm_lst = []
        aux_names = ['innerPCDsolver_Ap_',
                     'innerPCDsolver_Qp_',
                     'innerLSCsolver_BTinvBt_',
                     'innerLSCsolver_T_']
        for name in aux_names:
            rm_lst.extend([name+g for g in ['ksp_type','pc_type']])
        for element in rm_lst:
            self.petsc_options.delValue(element)

    def test_pcd_shell(self):
        '''  Tests the pcd_shell operators produce correct output. '''
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
        self.petsc_options.setValue('innerPCDsolver_Ap_ksp_type','preonly')
        self.petsc_options.setValue('innerPCDsolver_Qp_ksp_type','preonly')
        self.petsc_options.setValue('innerPCDsolver_Ap_pc_type','lu')
        self.petsc_options.setValue('innerPCDsolver_Qp_pc_type','lu')
        PCD_shell = LinearAlgebraTools.PCDInv_shell(petsc_matQp , petsc_matFp , petsc_matAp)
        x_vec = np.ones(size_n)
        y_vec = np.zeros(size_n)
        x_PETSc_vec = p4pyPETSc.Vec().createWithArray(x_vec)
        y_PETSc_vec = p4pyPETSc.Vec().createWithArray(y_vec)
        A = None
        PCD_shell.apply(A,x_PETSc_vec,y_PETSc_vec)
        true_vec = np.mat('[0.52511998  -0.13008364  -0.03200969]')
        assert np.allclose(y_vec,true_vec)

    def test_BAinvBt_shell(self):
        ''' Test the B_Ainv_Bt_shell '''
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

    def test_lsc_shell(self):
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
        self.petsc_options.setValue('innerLSCsolver_BTinvBt_ksp_type','preonly')
        self.petsc_options.setValue('innerLSCsolver_T_ksp_type','preonly')
        self.petsc_options.setValue('innerLSCsolver_BTinvBt_pc_type','lu')
        self.petsc_options.setValue('innerLSCsolver_T_pc_type','lu')
        LSC_shell = LinearAlgebraTools.LSCInv_shell(petsc_matA,petsc_matB,petsc_matF)
        x_vec = np.ones(num_B_rows)
        y_vec = np.zeros(num_B_rows)
        x_PETSc_vec = p4pyPETSc.Vec().createWithArray(x_vec)
        y_PETSc_vec = p4pyPETSc.Vec().createWithArray(y_vec)
        A = None
        LSC_shell.apply(A,x_PETSc_vec,y_PETSc_vec)
        true_solu = np.mat('[-0.01096996,0.00983216]')
        assert np.allclose(y_vec,true_solu)

    def test_twp_pcd_shell_1(self):
        ''' test example for a two-phase pcd shell with no temproal component '''
        Qp_rho_vals = [0.1,0.4,7.1,0.4,1.2,3.4,8]
        Qp_rho_col_idx = [0, 1, 1, 2, 0, 1, 2]
        Qp_rho_row_idx = [0,2,4,7]
        Qp_rho_num_rows = len(Qp_rho_row_idx) - 1
        Qp_rho_petsc = p4pyPETSc.Mat().createAIJ(size = (Qp_rho_num_rows,Qp_rho_num_rows),
                                                 csr = (Qp_rho_row_idx,Qp_rho_col_idx,Qp_rho_vals))

        Qp_visc_vals = [1.1, 2.2, 0.3, 4.0, 8.3, 0.3, 6., 5.]
        Qp_visc_col_idx = [0, 1, 2, 0, 1, 2, 0, 2]
        Qp_visc_row_idx = [0, 3, 6, 8]
        Qp_visc_num_rows = len(Qp_visc_row_idx) - 1
        Qp_visc_petsc = p4pyPETSc.Mat().createAIJ(size = (Qp_visc_num_rows,Qp_visc_num_rows),
                                                 csr = (Qp_visc_row_idx,Qp_visc_col_idx,Qp_visc_vals))

        Np_vals = [3.2, 1.1, 6.3, 1.0, -5.1]
        Np_col_idx = [0, 1, 1, 2, 2]
        Np_row_idx = [0, 2, 4, 5]
        Np_num_rows = len(Np_row_idx) - 1
        Np_petsc = p4pyPETSc.Mat().createAIJ(size = (Np_num_rows,Np_num_rows),
                                                 csr = (Np_row_idx,Np_col_idx,Np_vals))

        Ap_vals = [0.6, 8.2, 1.0, 0.2, 1.2, 2.5]
        Ap_col_idx = [0, 1, 2, 1, 0, 2]
        Ap_row_idx = [0, 3, 4, 6]
        Ap_num_rows = len(Ap_row_idx) - 1
        Ap_petsc = p4pyPETSc.Mat().createAIJ(size = (Ap_num_rows,Ap_num_rows),
                                                 csr = (Ap_row_idx,Ap_col_idx,Ap_vals))

        self.petsc_options.setValue('innerTPPCDsolver_Qp_visc_ksp_type','preonly')
        self.petsc_options.setValue('innerTPPCDsolver_Qp_dens_ksp_type','preonly')
        self.petsc_options.setValue('innerTPPCDsolver_Ap_rho_ksp_type','preonly')
        self.petsc_options.setValue('innerTPPCDsolver_Qp_visc_pc_type','lu')
        self.petsc_options.setValue('innerTPPCDsolver_Qp_dens_pc_type','lu')
        self.petsc_options.setValue('innerTPPCDsolver_Ap_rho_pc_type','lu')

        TPPCD_shell_1 = LinearAlgebraTools.TwoPhase_PCDInv_shell(Qp_visc_petsc,
                                                                 Qp_rho_petsc,
                                                                 Ap_petsc,
                                                                 Np_petsc)

        TPPCD_shell_2 = LinearAlgebraTools.TwoPhase_PCDInv_shell(Qp_visc_petsc,
                                                                 Qp_rho_petsc,
                                                                 Ap_petsc,
                                                                 Np_petsc,
                                                                 True,
                                                                 4)        

        x_vec = np.ones(Ap_num_rows)
        y_vec = np.zeros(Ap_num_rows)
        x_vec_petsc = p4pyPETSc.Vec().createWithArray(x_vec)
        y_vec_petsc = p4pyPETSc.Vec().createWithArray(y_vec)
        A = None
        TPPCD_shell_1.apply(A,x_vec_petsc,y_vec_petsc)
        true_sol = np.mat('[211.32623352, 1.58465852, -96.29629194]')
        assert np.allclose(y_vec,true_sol)
        TPPCD_shell_2.apply(A,x_vec_petsc,y_vec_petsc)
        true_sol = np.mat('[127.15956685,2.83465852,-55.79629194]')
        assert np.allclose(y_vec,true_sol)
        
if __name__ == '__main__':
    pass

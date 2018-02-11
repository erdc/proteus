#!/usr/bin/env python
"""

Test module for solver shell operator.

"""
import proteus.test_utils.TestTools
from proteus.iproteus import *
from proteus import Comm
from proteus import LinearAlgebraTools as LAT

import os
import sys
import inspect
import petsc4py
import numpy as np
import pytest

from petsc4py import PETSc as p4pyPETSc
from scipy.sparse import csr_matrix

proteus.test_utils.TestTools.addSubFolders( inspect.currentframe() )

@pytest.fixture()
def create_simple_saddle_point_problem(request):
    """Builds simple matrices and vectors for saddle point shell tests.

    Returns
    -------
    output_lst : lst
        This function returns an output list with three sublists.  See
        the notes below for a description of these lists.

    Notes
    -----
    The output of this function returns xxx sublists.  The first
    contains a list of matrices.  The second returns a list
    of vectors.  The last returns a list of sizes.
    """
    class Output_Storage(object):
        """Storage class for matrix objects. """
        def __init__(self,
                     petsc_matF,
                     petsc_matD,
                     petsc_matB,
                     petsc_matBt,
                     petsc_matC,
                     x_vec,
                     y_vec,
                     num_p_unkwn,
                     num_v_unkwn):
            self.petsc_matF = petsc_matF
            self.petsc_matD = petsc_matD
            self.petsc_matB = petsc_matB
            self.petsc_matBt = petsc_matBt
            self.petsc_matC = petsc_matC
            self.x_vec = x_vec
            self.y_vec = y_vec
            self.num_p_unkwn = num_p_unkwn
            self.num_v_unkwn = num_v_unkwn

    vals_F  =    [3.2, 1.1, 6.3, 1., -5.1]
    col_idx_F  = [0, 1, 0, 2, 0]
    row_idx_F  = [0, 2, 4, 5]

    vals_D =     [5.5,7.1,1.0]
    col_idx_D =  [0 , 1 , 2  ]
    row_idx_D =  [0, 1, 2, 3]

    vals_B    =  [1.1, 6.3, 7.3, 3.6, 6.3]
    col_idx_B =  [0  , 2  , 0  , 1  , 2  ]
    row_idx_B =  [0, 2, 5]

    vals_Bt   =  [1.1, 7.3, 3.6, 6.3, 6.3]
    col_idx_Bt = [0, 1, 1, 0, 1]
    row_idx_Bt = [0, 2, 3, 5]

    vals_C = [1.2, 2.1, 3.3]
    col_idx_C = [0, 1, 1]
    row_idx_C = [0, 2, 3]

    num_p_unkwn = len(row_idx_B) - 1
    num_v_unkwn = len(row_idx_F) - 1

    petsc_matF = LAT.csr_2_petsc(size = (num_v_unkwn,num_v_unkwn),
                                 csr = (row_idx_F,col_idx_F,vals_F))
    petsc_matD = LAT.csr_2_petsc(size = (num_v_unkwn,num_v_unkwn),
                                 csr = (row_idx_D,col_idx_D,vals_D))
    petsc_matB = LAT.csr_2_petsc(size = (num_p_unkwn,num_v_unkwn),
                                 csr = (row_idx_B,col_idx_B,vals_B))
    petsc_matBt = LAT.csr_2_petsc(size = (num_v_unkwn,num_p_unkwn),
                                  csr = (row_idx_Bt,col_idx_Bt,vals_Bt))
    petsc_matC = LAT.csr_2_petsc(size = (num_p_unkwn,num_p_unkwn),
                                 csr = (row_idx_C,col_idx_C,vals_C))

    x_vec = np.ones(num_p_unkwn)
    y_vec = np.zeros(num_p_unkwn)
    x_PETSc_vec = p4pyPETSc.Vec().createWithArray(x_vec)
    y_PETSc_vec = p4pyPETSc.Vec().createWithArray(y_vec)

    output_data = Output_Storage(petsc_matF,
                                 petsc_matD,
                                 petsc_matB,
                                 petsc_matBt,
                                 petsc_matC,
                                 x_PETSc_vec,
                                 y_PETSc_vec,
                                 num_p_unkwn,
                                 num_v_unkwn)

    yield output_data

@pytest.mark.LinearAlgebraTools
class TestOperatorShells(proteus.test_utils.TestTools.BasicTest):

    def setup_method(self,method):
        self.petsc_options = p4pyPETSc.Options()
        self._scriptdir = os.path.dirname(__file__)

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

    def test_lsc_shell(self, create_simple_saddle_point_problem):
        ''' Test for the lsc operator shell '''
        fixture_data = create_simple_saddle_point_problem

        self.petsc_options.setValue('innerLSCsolver_BTinvBt_ksp_type','preonly')
        self.petsc_options.setValue('innerLSCsolver_T_ksp_type','preonly')
        self.petsc_options.setValue('innerLSCsolver_BTinvBt_pc_type','lu')
        self.petsc_options.setValue('innerLSCsolver_T_pc_type','lu')

        LSC_shell = LAT.LSCInv_shell(fixture_data.petsc_matD,
                                     fixture_data.petsc_matB,
                                     fixture_data.petsc_matBt,
                                     fixture_data.petsc_matF)
        LSC_shell.apply(None,
                        fixture_data.x_vec,
                        fixture_data.y_vec)
        true_solu = np.mat('[-0.01096996,0.00983216]')
        assert np.allclose(fixture_data.y_vec.getArray(),
                           true_solu)

    def test_tppcd_shell(self):
        ''' Test for the two-phase pcd operator shell '''
        Qp_visc = LAT.petsc_load_matrix(os.path.join(self._scriptdir,'import_modules/Qp_visc'))
        Qp_dens = LAT.petsc_load_matrix(os.path.join(self._scriptdir,'import_modules/Qp_dens'))
        Ap_rho = LAT.petsc_load_matrix(os.path.join(self._scriptdir, 'import_modules/Ap_rho'))
        Np_rho = LAT.petsc_load_matrix(os.path.join(self._scriptdir, 'import_modules/Np_rho'))
        alpha = True
        delta_t = 0.001
        x_vec = LAT.petsc_load_vector(os.path.join(self._scriptdir,'import_modules/input_vec_tppcd'))

        self.petsc_options.setValue('innerTPPCDsolver_Ap_rho_ksp_type','preonly')
        self.petsc_options.setValue('innerTPPCDsolver_Ap_rho_ksp_constant_null_space', '')
        self.petsc_options.setValue('innerTPPCDsolver_Ap_rho_pc_type','hypre')
        self.petsc_options.setValue('innerTPPCDsolver_Ap_rho_pc_hypre_type','boomeramg')                
        TPPCD_shell = LAT.TwoPhase_PCDInv_shell(Qp_visc,
                                                Qp_dens,
                                                Ap_rho,
                                                Np_rho,
                                                alpha,
                                                delta_t,
                                                5)
        y_vec = x_vec.copy()
        y_vec.zeroEntries()
        A = None
        TPPCD_shell.apply(A,x_vec,y_vec)
        true_solu = LAT.petsc_load_vector(os.path.join(self._scriptdir, 'import_modules/tp_pcd_y_output'))
        assert np.allclose(y_vec.getArray(),true_solu.getArray())

    def test_SpInv_shell(self, create_simple_saddle_point_problem):
        """Test :math:`S_{p}` shell has correct behavior. """
        fixture_data = create_simple_saddle_point_problem

        self.petsc_options.setValue('innerSpsolver_ksp_type', 'preonly')
        self.petsc_options.setValue('innerSpsolver_pc_type', 'hypre')
        self.petsc_options.setValue('innerSpsolver_pc_hypre_type', 'boomeramg')

        SpInv_shell = LAT.SpInv_shell(fixture_data.petsc_matF,
                                      fixture_data.petsc_matC,
                                      fixture_data.petsc_matBt,
                                      fixture_data.petsc_matB,
                                      constNullSpace=False)
        SpInv_shell.apply(None,
                          fixture_data.x_vec,
                          fixture_data.y_vec)
        true_solu = np.mat('[1.0655362, -0.30354183]')
        assert np.allclose(fixture_data.y_vec.getArray(),
                           true_solu)

if __name__ == '__main__':
    pass

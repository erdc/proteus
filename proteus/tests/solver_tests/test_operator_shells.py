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

from petsc4py import PETSc
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
    x_PETSc_vec = PETSc.Vec().createWithArray(x_vec)
    y_PETSc_vec = PETSc.Vec().createWithArray(y_vec)

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

def setup_LSC_shell(petsc_options, fixture_data):
    petsc_options.clear()
    for k in petsc_options.getAll(): petsc_options.delValue(k)
    petsc_options.setValue('innerLSCsolver_BTinvBt_ksp_type','preonly')
    petsc_options.setValue('innerLSCsolver_T_ksp_type','preonly')
    petsc_options.setValue('innerLSCsolver_BTinvBt_pc_type','lu')
    petsc_options.setValue('innerLSCsolver_T_pc_type','lu')

    return LAT.LSCInv_shell(fixture_data.petsc_matD,
                            fixture_data.petsc_matB,
                            fixture_data.petsc_matBt,
                            fixture_data.petsc_matF)

@pytest.fixture()
def create_simple_petsc_matrix(request):
    vals_F  =    [3.2, 1.1, 6.3, 1., -5.1]
    col_idx_F  = [0, 1, 0, 1, 2]
    row_idx_F  = [0, 2, 4, 5]
    num_v_unkwn = 3
    petsc_matF = LAT.csr_2_petsc(size = (num_v_unkwn,num_v_unkwn),
                                 csr = (row_idx_F,col_idx_F,vals_F))
    yield petsc_matF

@pytest.mark.LinearAlgebraTools
class TestOperatorShells(proteus.test_utils.TestTools.BasicTest):

    def setup_method(self,method):
        self.petsc_options = PETSc.Options()
        self.petsc_options.clear()
        for k in self.petsc_options.getAll(): self.petsc_options.delValue(k)
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
        self.petsc_options.clear()
        for k in self.petsc_options.getAll(): self.petsc_options.delValue(k)

    def test_lsc_shell(self, create_simple_saddle_point_problem):
        ''' Test for the lsc operator shell '''
        fixture_data = create_simple_saddle_point_problem
        LSC_shell = setup_LSC_shell(self.petsc_options,
                                    fixture_data)
        LSC_shell.apply(None,
                        fixture_data.x_vec,
                        fixture_data.y_vec)
        true_solu = np.array([-0.01096996,0.00983216],'d')
        np.testing.assert_allclose(fixture_data.y_vec.getArray(),
                                   true_solu,
                                   rtol=1e-8, atol=1e-8)

    def test_tppcd_shell(self):
        ''' Test for the two-phase pcd operator shell '''
        Qp_visc = LAT.petsc_load_matrix(os.path.join(self._scriptdir,'import_modules/Qp_visc.bin'))
        Qp_dens = LAT.petsc_load_matrix(os.path.join(self._scriptdir,'import_modules/Qp_dens.bin'))
        Ap_rho = LAT.petsc_load_matrix(os.path.join(self._scriptdir, 'import_modules/Ap_rho.bin'))
        Np_rho = LAT.petsc_load_matrix(os.path.join(self._scriptdir, 'import_modules/Np_rho.bin'))
        alpha = True
        delta_t = 0.001
        x_vec = LAT.petsc_load_vector(os.path.join(self._scriptdir,'import_modules/input_vec_tppcd.bin'))
        self.petsc_options.clear()
        for k in self.petsc_options.getAll(): self.petsc_options.delValue(k)
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
                                                0,
                                                laplace_null_space=True)
        y_vec = x_vec.copy()
        y_vec.zeroEntries()
        A = None
        TPPCD_shell.apply(A,x_vec,y_vec)
        #np.savetxt(os.path.join(self._scriptdir,'comparison_files/tp_pcd_y_output.csv'),y_vec.getArray(),delimiter=',')
        true_solu = np.loadtxt(os.path.join(self._scriptdir,'comparison_files/tp_pcd_y_output.csv'),delimiter=',')
        np.testing.assert_allclose(y_vec.getArray(),
                                   true_solu,
                                   rtol=1e-8, atol=1e-8)

    def test_tppcd_shell_with_chebyshev_iteration(self):
        ''' Test for the lsc operator shell '''
        Qp_visc = LAT.petsc_load_matrix(os.path.join(self._scriptdir,'import_modules/Qp_visc.bin'))
        Qp_dens = LAT.petsc_load_matrix(os.path.join(self._scriptdir,'import_modules/Qp_dens.bin'))
        Ap_rho = LAT.petsc_load_matrix(os.path.join(self._scriptdir, 'import_modules/Ap_rho.bin'))
        Np_rho = LAT.petsc_load_matrix(os.path.join(self._scriptdir, 'import_modules/Np_rho.bin'))
        alpha = True
        delta_t = 0.001
        x_vec = LAT.petsc_load_vector(os.path.join(self._scriptdir,'import_modules/input_vec_tppcd.bin'))
        self.petsc_options.clear()
        for k in self.petsc_options.getAll(): self.petsc_options.delValue(k)
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
                                                5,
                                                laplace_null_space=True)
        y_vec = x_vec.copy()
        y_vec.zeroEntries()
        A = None
        TPPCD_shell.apply(A,x_vec,y_vec)
        #np.savetxt(os.path.join(self._scriptdir,'comparison_files/tp_pcd_y_output_cheb.csv'),y_vec.getArray(),delimiter=',')
        true_solu = np.loadtxt(os.path.join(self._scriptdir,'comparison_files/tp_pcd_y_output_cheb.csv'),delimiter=',')
        np.testing.assert_allclose(y_vec.getArray(),
                                       true_solu,
                                       rtol=1e-8, atol=1e-8)

    def test_SpInv_shell(self, create_simple_saddle_point_problem):
        """Test :math:`S_{p}` shell has correct behavior. """
        fixture_data = create_simple_saddle_point_problem
        self.petsc_options.clear()
        for k in self.petsc_options.getAll(): self.petsc_options.delValue(k)
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
        true_solu = np.array([1.0655362, -0.30354183],'d')
        np.testing.assert_allclose(fixture_data.y_vec.getArray(),
                                       true_solu,
                                       rtol=1e-8, atol=1e-8)

    def test_tppcd_shell_with_dirichlet_pressure(self):
        ''' Test for the lsc operator shell '''
        Qp_visc = LAT.petsc_load_matrix(os.path.join(self._scriptdir,'import_modules/Qp_visc.bin'))
        Qp_dens = LAT.petsc_load_matrix(os.path.join(self._scriptdir,'import_modules/Qp_dens.bin'))
        Ap_rho = LAT.petsc_load_matrix(os.path.join(self._scriptdir, 'import_modules/Ap_rho.bin'))
        Np_rho = LAT.petsc_load_matrix(os.path.join(self._scriptdir, 'import_modules/Np_rho.bin'))
        alpha = True
        delta_t = 0.001
        x_vec = LAT.petsc_load_vector(os.path.join(self._scriptdir,'import_modules/input_vec_tppcd.bin'))
        self.petsc_options.clear()
        for k in self.petsc_options.getAll(): self.petsc_options.delValue(k)
        self.petsc_options.setValue('innerTPPCDsolver_Ap_rho_ksp_type','preonly')
        self.petsc_options.setValue('innerTPPCDsolver_Ap_rho_ksp_constant_null_space', '')
        self.petsc_options.setValue('innerTPPCDsolver_Ap_rho_pc_type','hypre')
        self.petsc_options.setValue('innerTPPCDsolver_Ap_rho_pc_hypre_type','boomeramg')

        dirichlet_nodes = [3, 12, 15, 21, 33]
        TPPCD_shell = LAT.TwoPhase_PCDInv_shell(Qp_visc,
                                                Qp_dens,
                                                Ap_rho,
                                                Np_rho,
                                                alpha = alpha,
                                                delta_t = delta_t,
                                                num_chebyshev_its = 5,
                                                strong_dirichlet_DOF = dirichlet_nodes)
        # Test the index set is set correctly
        unknown_indices = np.arange(TPPCD_shell.getSize())
        known_indices_mask = np.ones(TPPCD_shell.getSize(),dtype=bool)
        known_indices_mask[dirichlet_nodes] = False
        unknown_is = unknown_indices[known_indices_mask]
        np.testing.assert_equal(unknown_is,
                                TPPCD_shell.unknown_dof_is.getIndices())

        y_vec = x_vec.copy()
        y_vec.zeroEntries()
        A = None
        TPPCD_shell.apply(A,x_vec,y_vec)
        assert np.array_equal(y_vec[dirichlet_nodes], [0.,0.,0.,0.,0.])
        #np.savetxt(os.path.join(self._scriptdir,'comparison_files/tppcd_y_dirichlet_dof.csv'),y_vec.getArray(),delimiter=',')
        true_solu = np.loadtxt(os.path.join(self._scriptdir,'comparison_files/tppcd_y_dirichlet_dof.csv'),delimiter=',')
        np.testing.assert_allclose(y_vec.getArray(),
                                       true_solu,
                                       rtol=1e-8, atol=1e-8)

def test_tmp_vec_creation():
    A = LAT.InvOperatorShell._create_tmp_vec(4)
    assert A.norm() == 0.0
    assert A.getSize() == 4

@pytest.mark.dev
def test_create_petsc_ksp_obj(create_simple_petsc_matrix):
    def getSize():
        return 3

    petsc_options = PETSc.Options()
    petsc_options.clear()
    for k in petsc_options.getAll(): petsc_options.delValue(k)
    petsc_options.setValue('test_F_ksp_type','preonly')
    petsc_options.setValue('test_F_pc_type','lu')
    petsc_matF = create_simple_petsc_matrix

    InvOpShell = LAT.InvOperatorShell()
    InvOpShell.const_null_space =False
    InvOpShell.options = petsc_options
    InvOpShell.strong_dirichlet_DOF = [1]
    InvOpShell.getSize = getSize
    A = InvOpShell.create_petsc_ksp_obj('test_F_',
                                        petsc_matF)
    assert A.getClassName()=='KSP'
    assert A.getOperators()[0].equal(petsc_matF)

    InvOpShell._set_dirichlet_idx_set()
    np.testing.assert_array_equal(InvOpShell.unknown_dof_is.getIndices(),
                                  np.array([0,2]))

if __name__ == '__main__':
    pass

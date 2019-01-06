#cython: boundscheck=False, wraparound=False, nonecheck=False, initializedcheck=False

# ARB - the cython directives above help improve the code speed.  it is worth noting
# however that these directives can lead to segfaults if the code is not implemented
# correctly.  perhaps they should only be enabled when not in debug mode?

import numpy as np
import cython
cimport numpy as np

# cdef extern from "proteus_superlu.h"  -- ARB: need to use this or some other non-hardcoded approach ...
cdef extern from "../linux2/include/slu_ddefs.h":
    ctypedef enum _fact_t 'fact_t':
        _DOFACT 'DFACT'
        _SamePattern_SameRowPerm "SamePattern_SameRowPerm"
    ctypedef enum _Stype_t 'Stype_t':
        _SLU_NC 'SLU_NC'
        _SLU_NCP 'SLU_NCP'
        _SLU_NR 'SLU_NR'
        _SLU_SC 'SLU_SC'
        _SLU_SCP 'SLU_SCP'
        _SLU_SR 'SLU_SR'
        _SLU_DN 'SLU_DN'
        _SLU_NR_loc 'SLU_NR_loc'
    ctypedef enum _Dtype_t 'Dtype_t':
        _SLU_S 'SLU_S'
        _SLU_D 'SLU_D'
        _SLU_C 'SLU_C'
        _SLU_Z 'SLU_Z'
    ctypedef enum _Mtype_t 'Mtype_t':
        _SLU_GE 'SLU_GE'
        _SLU_TRLU 'SLU_TRLU'
        _SLU_TRUU 'SLU_TRUU'
        _SLU_TRL 'SLU_TRL'
        _SLU_TRU 'SLU_TRU'
        _SLU_SYL 'SLU_SYL'
        _SLU_SYU 'SLU_SYU'
        _SLU_HEL 'SLU_HEL'
        _SLU_HEU 'SLU_HEU'
    ctypedef enum _trans_t "trans_t":
        _NOTRANS 'NOTRANS'
        _TRANS 'TRANS'
        _CONJ 'CONJ'
    # ctypedef struct _NRformat 'NRformat':
    #     int nnz
    #     void *nzval
    #     int *colind
    #     int *rowptr
    ctypedef struct _SuperLUStat_t 'SuperLUStat_t':
        pass
    ctypedef struct _GlobalLU_t "GlobalLU_t":
        pass
    ctypedef struct _superlu_options_t "superlu_options_t":
        pass
    ctypedef struct _SuperMatrix 'SuperMatrix':
        _Stype_t Stype
        _Dtype_t Dtype
        _Mtype_t Mtype
        int nrow
        int ncol
        void * Store
    void cdgstrs "dgstrs"(_trans_t, _SuperMatrix *, _SuperMatrix *, int *, int *, _SuperMatrix *, _SuperLUStat_t*, int *)
    void cdgstrf "dgstrf"(_superlu_options_t *, _SuperMatrix *, int, int, int *, void *, int, int *, int *, _SuperMatrix *, _SuperMatrix *, _GlobalLU_t*, _SuperLUStat_t *, int *)
    void cStatInit "StatInit"(_SuperLUStat_t *)
    void cset_default_options "set_default_options"(_superlu_options_t *)
    void cget_perm_c "get_perm_c"(int, _SuperMatrix *, int *)
    void cDestroy_CompCol_Permuted "Destroy_CompCol_Permuted"(_SuperMatrix *)
    void cDestroy_SuperNode_Matrix "Destroy_SuperNode_Matrix"(_SuperMatrix *)
    void cDestroy_CompCol_Matrix "Destroy_CompCol_Matrix"(_SuperMatrix *)
    void csp_preorder "sp_preorder"(_superlu_options_t *, _SuperMatrix *, int *, int *, _SuperMatrix *)

class SparseMatrix(object):

    def __init__(self,
                  nr,
                  nc,
                  nnz,
                  nzvals,
                  colind,
                  rowptr):
        self.nr = nr ; self.nc = nc
        self.nzvals = nzvals
        self.colind = colind
        self.rowptr = rowptr
        self._cSparseMatrix = cSparseMatrix(nr,
                                            nc,
                                            nnz,
                                            self.nzvals,
                                            self.colind,
                                            self.rowptr)

    def matvec(self, x, y):
        """
        Compute the sparse matrix-vector product y = Ax

        Arguments
        ---------
        x (input) :  numpy array
        y (output) : numpy array
        """
        SparseMatrix_matvec(self._cSparseMatrix, x, y)

    def fwrite(self, filename, base):
        """ Write the sparse matrix to a file

        Arguments
        ---------
        filename : str
            The output filename
        base : int
            ?!Possibly something to do with parallel?!
        """
        pass

    def getCSRrepresentation(self):
        """ Get the CSR representation of the sparse matrix.

        Returns
        -------
        csr_data : tuple of nparray
            (rowptr, colptr, vals)
        """
        return (self.rowptr, self.colind, self.nzvals)

    def getSubMatCSRrepresentation(self,
                                   range_start,
                                   range_end):
        """  Get the CSR representation for a submatrix.

        Arguments
        ---------
        range_start : int
        range_end : int

        Returns
        -------
        csr_data : tuple of nparray
            (rowptr, colind, nzvals)
        """
        assert range_start >= 0 ; assert range_end <= self.nr
        assert range_end > range_start
        _rows = range_end - range_start
        assert _rows <= self.nr

        nnz = self.rowptr[range_end] - self.rowptr[range_start]

        rowptr = self.rowptr[range_start : range_start + _rows + 1]
        colind = self.colind[self.rowptr[range_start] : self.rowptr[range_start] + nnz]
        nzvals = self.nzvals[self.rowptr[range_start] : self.rowptr[range_start] + nnz]

        return rowptr, colind, nzvals

cdef struct _NRformat:
    int nnz
    np.float64_t [:] nzval
    np.int32_t [:] colind
    np.int32_t [:] rowptr

cdef class cSparseMatrix(object):

    cdef int dim[2]
    cdef _NRformat A

    def __cinit__(self,
                 int nr,
                 int nc,
                 int nnz,
                 np.float64_t [:] nzval,
                 np.int32_t [:] colind,
                 np.int32_t [:] rowptr):
        self.dim[0] = nr ; self.dim[1] = nc
        self.A.nnz = nnz
        self.A.nzval = nzval
        self.A.colind = colind
        self.A.rowptr = rowptr        

cdef void SparseMatrix_matvec(cSparseMatrix sm,
                              np.float64_t [:] xp,
                              np.float64_t [:] yp):
    cdef np.float64_t tmp = 0.
    cdef int i, k

    for i in range(sm.dim[0]):
        tmp = 0.
        for k in range(sm.A.rowptr[i], sm.A.rowptr[i+1]):
            tmp += sm.A.nzval[k] * xp[sm.A.colind[k]]
        yp[i] = tmp

# cdef class SparseFactor(object):

#     cdef _superlu_options_t options

#     cdef _SuperMatrix A
#     cdef _SuperMatrix AC
#     cdef _SuperMatrix L
#     cdef _SuperMatrix U
#     cdef _SuperMatrix X

#     cdef _GlobalLU_t Glu
#     cdef _SuperLUStat_t stat

#     cdef int *perm_c
#     cdef int *perm_r
#     cdef int *etree

#     cdef unsigned int use_same_perm_c
#     cdef unsigned int use_same_sparsity

#     cdef public int dim

#     def __init__(self, dim):
#         """
#         Arguments
#         ---------
#         dim : int
#             Dimension of the sparse factor.
#         """
#         cStatInit(&self.stat)
#         cset_default_options(&self.options)
#         self._set_mat_types()
#         self.dim = dim
#         self.A.nrow = dim ; self.A.ncol = dim
#         self.AC.nrow = dim ; self.AC.ncol = dim
#         self.L.nrow = dim ; self.L.ncol = dim
#         self.U.nrow = dim ; self.U.ncol = dim
#         self.X.nrow = dim ; self.X.ncol = 1
#         self.use_same_perm_c = 0
#         self.use_same_sparsity = 0

#     cdef _set_mat_types(self):
#         self.A.Stype = _SLU_NC
#         self.A.Dtype = _SLU_D
#         self.A.Mtype = _SLU_GE

#         self.AC.Stype = _SLU_NCP
#         self.AC.Dtype = _SLU_D
#         self.AC.Mtype = _SLU_GE

#         self.L.Stype = _SLU_NC
#         self.L.Dtype = _SLU_D
#         self.L.Mtype = _SLU_TRLU

#         self.U.Stype = _SLU_NC
#         self.U.Dtype = _SLU_D
#         self.U.Mtype = _SLU_TRU

#         self.X.Stype = _SLU_DN
#         self.X.Dtype = _SLU_D
#         self.X.Mtype = _SLU_GE

# def sparseFactorPrepare(sparse_matrix,
#                         sparseFactor):
#     """ Python wrapper for superlu Sparse Factor Prepare function.

#     Arguments
#     ---------
#     sparse_matrix : petsc_mat
#         one level transport matrix object?
#     sparseFactor: superluWrappers.SparseFactor

#     """
#     superluWrappersSparseFactorPrepare(sparse_matrix,
#                                        sparseFactor)

# cdef void superluWrappersSparseFactorPrepare(sparse_matrix,
#                                              SparseFactor sparseFactor):
#     cdef int permc_spec = 3
#     cdef int n
#     cdef int relax=1
#     cdef int panel_size = 10
#     cdef int lwork = 0
#     cdef int info = 0
#     cdef void *work = NULL

#     sparseFactor.nnz = mat.A.nnz
#     sparseFactor.nzval = mat.A.nzval
#     sparseFactor.colptr = mat.A.rowptr
#     sparseFactor.rowind = mat.A.colind

#     if sparseFactor.use_same_perm_c == 0:
#         cget_perm_c(permc_spec,
#                     &sparseFactor.A,
#                     sparseFactor.perm_c)
#         sparseFactor.use_same_perm_c =1

#     if sparseFactor.use_same_sparsity == 0:
#         if sparseFactor.AC.Store != NULL:
#             cDestroy_CompCol_Permuted(&sparseFactor.AC)
#             cDestroy_SuperNode_Matrix(&sparseFactor.L)
#             cDestroy_CompCol_Matrix(&sparseFactor.U)
#         csp_preorder(&sparseFactor.options,
#                      &sparseFactor.A,
#                      sparseFactor.perm_c,
#                      sparseFactor.etree,
#                      &sparseFactor.AC)
#         sparseFactor.use_same_sparsity = 1
#     else:
#         pass
#         #     # ARB - need to lookinto this part in more detail,
#         #     # i'm not following all the typecasts?
#         # sparseFactor.options.Fact = _SamePattern_SameRowPerm
#         # n = sparseFactor.A.ncol
#         # for i in range(n):
#         #     sparseFactor.AC.Store.colberg[sparseFactor.perm_c[i]] = sparseFactor.A.Store.colptr[i]
#         #     sparseFactor.AC.Store.colend[sparseFactor.per_c[i]] = sparseFactor.colptr[i+1]
#     cdgstrf(&sparseFactor.options,
#             &sparseFactor.AC,
#             relax,
#             panel_size,
#             sparseFactor.etree,
#             work,
#             lwork,
#             sparseFactor.perm_c,
#             sparseFactor.perm_r,
#             &sparseFactor.L,
#             &sparseFactor.U,
#             &sparseFactor.Glu,
#             &sparseFactor.stat,
#             &info)

# cdef void superluWrappersSparseFactorSolve(SparseFactor sparseFactor,
#                                            x):
#     cdef _trans_t trans = _TRANS
#     cdef int info = 0

#     sparseFactor.storeX.nzval = x
#     cdgstrs(trans,
#             &sparseFactor.L,
#             &sparseFactor.U,
#             sparseFactor.perm_c,
#             sparseFactor.perm_r,
#             &sparseFactor.X,
#             &sparseFactor.stat,
#             &info)

# def sparseFactorSolve(sparseFactor,
#                       x):
#     superluWrappersSparseFactorSolve(sparseFactor,
#                                      x)

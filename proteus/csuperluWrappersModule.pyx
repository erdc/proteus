import numpy as np
cimport numpy as np

# cdef extern from "proteus_superlu.h"  -- ARB: need to use this or some other non-hardcoded approach ...
cdef extern from "../linux2/include/slu_ddefs.h":
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
    void cdgstrs "dgstrs"(trans_t, _SuperMatrix *, _SuperMatrix *, int *, int *, _SuperMatrix *, _SuperLUStat_t*, int *)
    void cdgstrf "dgstrf"(_superlu_options_t *, _SuperMatrix *, int, int, int *, void *, int, int *, int *, _SuperMatrix *, _SuperMatrix *, _GlobalLU_t*, _SuperLUStat_t *, int *)
    void cStatInit "StatInit"(_SuperLUStat_t *)
    void cset_default_options "set_default_options"(_superlu_options_t *)

class SparseMatrix(object):
    
    def __init__(self,
                 nr,
                 nc,
                 nnz,
                 nzval,
                 colind,
                 rowptr):
        self.dim = (nr, nc)

    def matvec(self, x, y):
        """
o        Compute the sparse matrix-vector product y = Ax

        Arguments
        ---------
        x :  numpy array

        Returns
        -------
        y : numpy array
        """
        pass

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

    def getCSRrepresentation():
        """ Get the CSR representation of the sparse matrix.

        Returns
        -------
        csr_data : tuple of nparray
            (rowptr, colptr, vals)
        """
        pass

    def getSubMatCSRrepresentation():
        """  Get the CSR representation for a submatrix.
        
        Arguments
        ---------
        range_start : int
        range_end : int

        Returns
        -------
        csr_data : tuple of nparray
            (rowptr, colptr, vals)
        """
        pass

cdef class SparseFactor(object):

    cdef _superlu_options_t options
    
    cdef _SuperMatrix A
    cdef _SuperMatrix AC
    cdef _SuperMatrix L
    cdef _SuperMatrix U
    cdef _SuperMatrix X

    cdef _SuperLUStat_t stat
    
    cdef unsigned int *perm_c
    cdef unsigned int *perm_r
    cdef unsigned int *etree

    cdef unsigned int use_same_perm_c
    cdef unsigned int use_same_sparsity
    
    cdef public int dim

    def __init__(self, dim):
        """
        Arguments
        ---------
        dim : int
            Dimension of the sparse factor.
        """
        cStatInit(&self.stat)
        cset_default_options(&self.options)
        self._set_mat_types()
        self.dim = dim        
        self.A.nrow = dim ; self.A.ncol = dim
        self.AC.nrow = dim ; self.AC.ncol = dim
        self.L.nrow = dim ; self.L.ncol = dim
        self.U.nrow = dim ; self.U.ncol = dim        
        self.X.nrow = dim ; self.X.ncol = 1
        self.use_same_perm_c = 0
        self.use_same_sparsity = 0

    cdef _set_mat_types(self):
        self.A.Stype = _SLU_NC
        self.A.Dtype = _SLU_D
        self.A.Mtype = _SLU_GE
    
        self.AC.Stype = _SLU_NCP
        self.AC.Dtype = _SLU_D
        self.AC.Mtype = _SLU_GE

        self.L.Stype = _SLU_NC
        self.L.Dtype = _SLU_D
        self.L.Mtype = _SLU_TRLU

        self.U.Stype = _SLU_NC
        self.U.Dtype = _SLU_D
        self.U.Mtype = _SLU_TRU

        self.X.Stype = _SLU_DN
        self.X.Dtype = _SLU_D
        self.X.Mtype = _SLU_GE

def sparseFactorPrepare(mat,
                        sparseFactor):
    """ Python wrapper for superlu Sparse Factor Prepare function.

    Arguments
    ---------
    mat : petsc_mat
        one level transport matrix object?
    sparseFactor: superluWrappers.SparseFactor

    """
    pass

cdef superluWrappersSparseFactorPrepare(mat,
                                        sparseFactor):
    sparseFactor.nnz = mat.A.nnz

def sparseFactorSolve():
    pass

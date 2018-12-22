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
    ctypedef struct SuperLUStat_t:
        pass
    ctypedef struct _SuperMatrix 'SuperMatrix':
        _Stype_t Stype
        _Dtype_t Dtype
        _Mtype_t Mtype
        int nrow
        int ncol
        void * Store
    void cdgstrs "dgstrs"(trans_t, _SuperMatrix *, _SuperMatrix *, int *, int *, _SuperMatrix *, SuperLUStat_t*, int *)

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
        Compute the sparse matrix-vector product y = Ax

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

class SparseFactor(object):

    def __init__(self, dim):
        """
        Arguments
        ---------
        dim : int
            Dimension of the sparse factor.
        """
        self.dim = dim
        self.A = _SuperMatrix()

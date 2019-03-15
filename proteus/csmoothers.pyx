import numpy as np
import cython
cimport numpy as np

class ASMFactor(object):

    def __init__(self, L):
        self.L = L
        self._cASMFactor = cASMFactor(self.L._cSparseMatrix)

cdef class cASMFactor(object):

    def __cinit__(self,
                  superluWrappers.cSparseMatrix L):
        cdef int rval = 0
        cdef SuperMatrix AS
        AS.Stype = superluWrappers._SLU_NR
        AS.Dtype = superluWrappers._SLU_D
        AS.Mtype = superluWrappers._SLU_GE
        AS.nrow = L.nr
        AS.ncol = L.nc
        AS.Store = &L.A        
        rval = casm_NR_init(&AS,
                            &self.subdomain_dim,
                            &self.l2g_L,
                            &self.subdomain_L,
                            &self.subdomain_R,
                            &self.subdomain_dX,
                            &self.subdomain_pivots)
        assert rval == 0

    def __dealloc__(self):
        casm_NR_free(self.N,
                     self.subdomain_dim,
                     self.l2g_L,
                     self.subdomain_L,
                     self.subdomain_R,
                     self.subdomain_dX,
                     self.subdomain_pivots)
            

class BASMFactor(object):

    def __init__(self, L, bs):
        self.L = L
        self.bs = bs
        self._cBASMFactor = cBASMFactor(self.L._cSparseMatrix,
                                        self.bs)

cdef class cBASMFactor(object):

    def __cinit__(self,
                  superluWrappers.cSparseMatrix L,
                  int bs):
        cdef int rval = 0
        cdef SuperMatrix AS
        AS.Stype = superluWrappers._SLU_NR
        AS.Dtype = superluWrappers._SLU_D
        AS.Mtype = superluWrappers._SLU_GE
        AS.nrow = L.nr
        AS.ncol = L.nc
        AS.Store = &L.A
        rval = cbasm_NR_init(bs,
                             &AS,
                             &self.subdomain_dim,
                             &self.l2g_L,
                             &self.subdomain_L,
                             &self.subdomain_R,
                             &self.subdomain_dX,
                             &self.subdomain_pivots,
                             &self.subdomain_col_pivots)
        assert rval == 0

    def __dealloc__(self):
        cbasm_NR_free(self.N,
                      self.subdomain_dim,
                      self.l2g_L,
                      self.subdomain_L,
                      self.subdomain_R,
                      self.subdomain_dX,
                      self.subdomain_pivots,
                      self.subdomain_col_pivots)
    
def jacobi_NR_prepare(A, w, tol, M):
    """

    Arguments
    ---------
    A : superluWrappers.SparseMatrix
    w : double
    tol : double
    M : np.array
    """
    smootherWrappersjacobi_NR_prepare(A._cSparseMatrix, w, tol, M)

cdef void smootherWrappersjacobi_NR_prepare(superluWrappers.cSparseMatrix sm,
                                            double w,
                                            double tol,
                                            DDATA M):
    cdef SuperMatrix AS
    AS.Stype = superluWrappers._SLU_NR
    AS.Dtype = superluWrappers._SLU_D
    AS.Mtype = superluWrappers._SLU_GE
    AS.nrow = sm.nr
    AS.ncol = sm.nc
    AS.Store = &sm.A
    cjacobi_NR_prepare(&AS, w, tol, &M[0])

def jacobi_NR_solve(A, M, R, node_order, dX):
    """
    
    Arguments
    ---------
    A : superluWrappers.SparseMatrix
    M : np.array double
    R : np.array double
    node_order : np.array int
    dX : np.array double
    """
    smootherWrappersjacobi_NR_solve(A._cSparseMatrix, M, R, node_order, dX)

cdef void smootherWrappersjacobi_NR_solve(superluWrappers.cSparseMatrix sm,
                                          DDATA M,
                                          DDATA R,
                                          IDATA node_order,
                                          DDATA dX):
    cdef SuperMatrix AS
    AS.Stype = superluWrappers._SLU_NR
    AS.Dtype = superluWrappers._SLU_D
    AS.Mtype = superluWrappers._SLU_GE
    AS.nrow = sm.nr
    AS.ncol = sm.nc
    AS.Store = &sm.A
    cjacobi_NR_solve(&AS, &M[0], &R[0], &node_order[0], &dX[0])

def nl_jacobi_NR_solve(A, R, node_order, w, tol, dX):
    """

    Arguments
    ---------
    A : superluWrappers.SparseMatrix
    R : np.array double
    node_order : np.array int
    w : np.float
    tol : np.float
    dX : np.array double
    """
    smootherWrappersnl_jacobi_NR_solve(A._cSparseMatrix, R, node_order, w, tol, dX)

cdef void smootherWrappersnl_jacobi_NR_solve(superluWrappers.cSparseMatrix sm,
                                             DDATA R,
                                             IDATA node_order,
                                             double w,
                                             double tol,
                                             DDATA dX):
    cdef SuperMatrix AS
    AS.Stype = superluWrappers._SLU_NR
    AS.Dtype = superluWrappers._SLU_D
    AS.Mtype = superluWrappers._SLU_GE
    AS.nrow = sm.nr
    AS.ncol = sm.nc
    AS.Store = &sm.A
    cnl_jacobi_NR_solve(&AS, &R[0], &node_order[0], w, tol, &dX[0])

def gauss_seidel_NR_preare(A, w, tol, M):
    """

    Arguments
    ---------
    A :
    w :
    tol :
    M :
    """
    smootherWrappersgauss_seidel_NR_prepare(A._cSparseMatrix, w, tol, M)

cdef void smootherWrappersgauss_seidel_NR_prepare(superluWrappers.cSparseMatrix sm,
                                                  double w,
                                                  double tol,
                                                  DDATA M):
    cdef SuperMatrix AS
    AS.Stype = superluWrappers._SLU_NR
    AS.Dtype = superluWrappers._SLU_D
    AS.Mtype = superluWrappers._SLU_GE
    AS.nrow = sm.nr
    AS.ncol = sm.nc
    AS.Store = &sm.A
    cgauss_seidel_NR_prepare(&AS, w, tol, &M[0])

def gauss_seidel_NR_solve(A, M, R, node_order, dX):
    """

    Arguments
    ---------
    A : superluWrappers.SparseMatrix
    M : np.array double
    R : np.array double
    node_order : np.array int
    dX : np.array double
    """
    smootherWrappersgauss_seidel_NR_solve(A._cSparseMatrix, M, R, node_order, dX)

cdef void smootherWrappersgauss_seidel_NR_solve(superluWrappers.cSparseMatrix sm,
                                                DDATA M,
                                                DDATA R,
                                                IDATA node_order,
                                                DDATA dX):
    cdef SuperMatrix AS
    AS.Stype = superluWrappers._SLU_NR
    AS.Dtype = superluWrappers._SLU_D
    AS.Mtype = superluWrappers._SLU_GE
    AS.nrow = sm.nr
    AS.ncol = sm.nc
    AS.Store = &sm.A
    cgauss_seidel_NR_solve(&AS, &M[0], &R[0], &node_order[0], &dX[0])

def nl_gauss_seidel_NR_solve(A, R, node_order, w, tol, dX):
    """
    
    Arguments
    ---------
    A : superluWrappers.SparseMatrix
    R : np.array double
    node_order : np.array int
    w : double
    tol : double
    dX : np.array double
    """
    smootherWrappers_nl_gauss_seidel_NR_solve(A._cSparseMatrix,
                                              R,
                                              node_order,
                                              w,
                                              tol,
                                              dX)

cdef smootherWrappers_nl_gauss_seidel_NR_solve(superluWrappers.cSparseMatrix sm,
                                               DDATA R,
                                               IDATA node_order,
                                               double w,
                                               double tol,
                                               DDATA dX):
    cdef SuperMatrix AS
    AS.Stype = superluWrappers._SLU_NR
    AS.Dtype = superluWrappers._SLU_D
    AS.Mtype = superluWrappers._SLU_GE
    AS.nrow = sm.nr
    AS.ncol = sm.nc
    AS.Store = &sm.A
    cnl_gauss_seidel_NR_solve(&AS, &R[0], &node_order[0], w, tol, &dX[0])

def asm_NR_prepare(A, asmFactor):
    """

    Arguments
    ---------
    A :
    asmFactor :
    """
    smootherWrappers_asm_NR_prepare(A._cSparseMatrix,
                                    asmFactor._cASMFactor)
    
cdef void smootherWrappers_asm_NR_prepare(superluWrappers.cSparseMatrix sm,
                                          cASMFactor asmFactor):
    cdef SuperMatrix AS
    AS.Stype = superluWrappers._SLU_NR
    AS.Dtype = superluWrappers._SLU_D
    AS.Mtype = superluWrappers._SLU_GE
    AS.nrow = sm.nr
    AS.ncol = sm.nc
    AS.Store = &sm.A
    casm_NR_prepare(&AS,
                    asmFactor.subdomain_dim,
                    asmFactor.l2g_L,
                    asmFactor.subdomain_L,
                    asmFactor.subdomain_pivots)

def asm_NR_solve(A, w, asmFactor, node_order, R, dX):
    """
    
    Arguments
    ---------
    A : superluWrappers.sparseMatrix
    w : double
    asmFactor : csmoothers.asmFactor
    node_order : np.array int
    R : np.array double
    dX : np.array double
    """
    smootherWrappers_asm_NR_solve(A._cSparseMatrix,
                                  w,
                                  asmFactor._cASMFactor,
                                  node_order,
                                  R,
                                  dX)

cdef void smootherWrappers_asm_NR_solve(superluWrappers.cSparseMatrix sm,
                                        double w,
                                        cASMFactor asmFactor,
                                        IDATA node_order,
                                        DDATA R,
                                        DDATA dX):
    cdef SuperMatrix AS
    AS.Stype = superluWrappers._SLU_NR
    AS.Dtype = superluWrappers._SLU_D
    AS.Mtype = superluWrappers._SLU_GE
    AS.nrow = sm.nr
    AS.ncol = sm.nc
    AS.Store = &sm.A
    casm_NR_solve(&AS,
                  w,
                  asmFactor.subdomain_L,
                  asmFactor.subdomain_dim,
                  asmFactor.l2g_L,
                  &R[0],
                  asmFactor.subdomain_R,
                  &node_order[0],
                  asmFactor.subdomain_dX,
                  &dX[0],
                  asmFactor.subdomain_pivots)

def basm_NR_prepare(A, basmFactor):
    """

    Arguments
    ---------
    A :
    basmFactor :
    """
    smootherWrappers_basm_NR_prepare(A._cSparseMatrix,
                                     basmFactor._cBASMFactor)

cdef void smootherWrappers_basm_NR_prepare(superluWrappers.cSparseMatrix sm,
                                           cBASMFactor basmFactor):
    cdef SuperMatrix AS
    AS.Dtype = superluWrappers._SLU_D
    AS.Mtype = superluWrappers._SLU_GE
    AS.nrow = sm.nr
    AS.ncol = sm.nc
    AS.Store = &sm.A
    cbasm_NR_prepare(basmFactor.bs,
                     basmFactor.N,
                     &AS,
                     basmFactor.subdomain_dim,
                     basmFactor.l2g_L,
                     basmFactor.subdomain_L,
                     basmFactor.subdomain_pivots,
                     basmFactor.subdomain_col_pivots)
    

def basm_NR_solve(A, w, basmFactor, node_order, R, dX):
    """
    
    Arguments
    ---------
    A :
    w :
    basmFactor :
    node_order :
    R :
    dX :
    """
    smootherWrappers_basm_NR_solve(A._cSparseMatrix,
                                   w,
                                   basmFactor._cBASMFactor,
                                   node_order,
                                   R,
                                   dX)

cdef void smootherWrappers_basm_NR_solve(superluWrappers.cSparseMatrix sm,
                                         double w,
                                         cBASMFactor basmFactor,
                                         IDATA node_order,
                                         DDATA R,
                                         DDATA dX):
    cdef SuperMatrix AS
    AS.Stype = superluWrappers._SLU_NR
    AS.Dtype = superluWrappers._SLU_D
    AS.Mtype = superluWrappers._SLU_GE
    AS.nrow = sm.nr
    AS.ncol = sm.nc
    AS.Store = &sm.A
    cbasm_NR_solve(basmFactor.bs,
                   basmFactor.N,
                   &AS,
                   w,
                   basmFactor.subdomain_L,
                   basmFactor.subdomain_dim,
                   basmFactor.l2g_L,
                   &R[0],
                   basmFactor.subdomain_R,
                   &node_order[0],
                   basmFactor.subdomain_dX,
                   &dX[0],
                   basmFactor.subdomain_pivots,
                   basmFactor.subdomain_col_pivots)

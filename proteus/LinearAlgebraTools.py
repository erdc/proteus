"""
Tools for n-dimensional linear algebra

Vectors are just numpy arrays, as are dense matrices. Sparse matrices
are CSR matrices. Parallel vector and matrix are built on top of those
representations using PETSc.

.. inheritance-diagram:: proteus.LinearAlgebraTools
   :parts: 1
"""
import numpy
import math
import sys
import superluWrappers
from .superluWrappers import *
from .Profiling import logEvent
from petsc4py import PETSc as p4pyPETSc
from . import flcbdfWrappers

def _petsc_view(obj, filename):
    """Saves object to disk using a PETSc binary viewer.
    """
    viewer = p4pyPETSc.Viewer().createBinary(filename, 'w')
    viewer(obj)
    viewer2 = p4pyPETSc.Viewer().createASCII(filename+".m", 'w')
    viewer2.pushFormat(1)
    viewer2(obj)
    viewer2.popFormat()

def petsc_load_matrix(filename):
    """ This function loads a PETSc matrix from a binary format.
    (Eg. what is saved using the petsc_view function).

    Parameters
    ----------
    filename : str
        This is the name of the binary with the file stored.

    Returns
    -------
    matrix : petsc4py matrix
        The matrix that is stored in the binary file.
    """
    try:
        viewer = p4pyPETSc.Viewer().createBinary(filename,'r')
        output = p4pyPETSc.Mat().load(viewer)
    except:
        print("Either you've entered an invalid file name or your object is not a matrix (try petsc_load_vector).")
    return output

def petsc_load_vector(filename):
    """ This function loads a PETSc vector from a binary format.
    (Eg. what is saved using the petsc_view function).

    Parameters
    ----------
    filename : str
        This is the name of the binary with the file stored.

    Returns
    -------
    matrix : petsc4py vector
        The matrix that is stored in the binary file.
    """
    try:
        viewer = p4pyPETSc.Viewer().createBinary(filename,'r')
        output = p4pyPETSc.Vec().load(viewer)
    except:
        print("Either you've entered an invalid file name or your object is not a vector (try petsc_load_matrix).")
    return output
    
def _pythonCSR_2_dense(rowptr,colptr,data,nr,nc,output=False):
    """ Takes python CSR datatypes and makes a dense matrix """
    dense_matrix = numpy.zeros(shape = (nr,nc), dtype='float')
    for idx in range(len(rowptr)-1):
        row_vals = data[rowptr[idx]:rowptr[idx+1]]
        for val_idx,j in enumerate(colptr[rowptr[idx]:rowptr[idx+1]]):
            dense_matrix[idx][j] = row_vals[val_idx]
    if output is not False:
        numpy.save(output,dense_matrix)
    return dense_matrix

def superlu_sparse_2_dense(sparse_matrix,output=False):
    """ Converts a sparse superluWrapper into a dense matrix.

    Parameters
    ----------
    sparse_matrix : 
    output : str
        Out file name to store the matrix.

    Returns
    -------
    dense_matrix : numpy array
        A numpy array storing the dense matrix.

    Notes
    -----
    This function should not be used for large matrices.
    """
    rowptr = sparse_matrix.getCSRrepresentation()[0]
    colptr = sparse_matrix.getCSRrepresentation()[1]
    data   = sparse_matrix.getCSRrepresentation()[2]
    nr     = sparse_matrix.shape[0]
    nc     = sparse_matrix.shape[1]
    return _pythonCSR_2_dense(rowptr,colptr,data,nr,nc,output)

def petsc4py_sparse_2_dense(sparse_matrix,output=False):
    """ Converts a PETSc4Py matrix to a dense numpyarray.

    Parameters
    ----------
    sparse_matrix : PETSc4py matrix
    output : str
        Output file name to store the matrix.

    Returns
    -------
    dense_matrix : numpy array
        A numpy array with the dense matrix.
    
    Notes
    -----
    This function is very inefficient for large matrices.
    """
    rowptr = sparse_matrix.getValuesCSR()[0]
    colptr = sparse_matrix.getValuesCSR()[1]
    data   = sparse_matrix.getValuesCSR()[2]
    nr     = sparse_matrix.getSize()[0]
    nc     = sparse_matrix.getSize()[1]
    return _pythonCSR_2_dense(rowptr,colptr,data,nr,nc,output)

def superlu_2_petsc4py(sparse_superlu):
    """ Copy a sparse superlu matrix to a sparse petsc4py matrix

    Parameters
    ----------
    sparse_matrix : :class:`proteus.superluWrappers.SparseMatrix`

    Returns
    -------
    sparse_matrix : PETSc4py matrix
    """
    rowptr, colind, nzval = sparse_superlu.getCSRrepresentation()
    A_rowptr = rowptr.copy()
    A_colind = colind.copy()
    A_nzval  = nzval.copy()
    nr       = sparse_superlu.shape[0]
    nc       = sparse_superlu.shape[1]
    A_petsc4py = p4pyPETSc.Mat().createAIJWithArrays((nr,nc),
                                                     (A_rowptr,
                                                      A_colind,
                                                      A_nzval))
    return A_petsc4py

class ParVec:
    """
    A parallel vector built on top of daetk's wrappers for petsc
    """
    def __init__(self,array,blockSize,n,N,nghosts=None,subdomain2global=None,blockVecType="simple"):#"block"
        import flcbdfWrappers
        self.dim_proc=n*blockSize
        if nghosts is None:
            if blockVecType=="simple":
                self.cparVec=flcbdfWrappers.ParVec(blockSize,n,N,-1,None,array,0)
            else:
                self.cparVec=flcbdfWrappers.ParVec(blockSize,n,N,-1,None,array,1)
        else:
            assert nghosts >= 0, "The number of ghostnodes must be non-negative"
            assert subdomain2global.shape[0] == (n+nghosts), ("The subdomain2global map is the wrong length n=%i,nghosts=%i,shape=%i \n" % (n,n+nghosts,subdomain2global.shape[0]))
            assert len(array.flat) == (n+nghosts)*blockSize, ("%i  != (%i+%i)*%i \n"%(len(array.flat),  n,nghosts,blockSize))
            if blockVecType=="simple":
                self.cparVec=flcbdfWrappers.ParVec(blockSize,n,N,nghosts,subdomain2global,array,0)
            else:
                self.cparVec=flcbdfWrappers.ParVec(blockSize,n,N,nghosts,subdomain2global,array,1)
        self.nghosts = nghosts
    def scatter_forward_insert(self):
       self.cparVec.scatter_forward_insert()
    def scatter_reverse_add(self):
       self.cparVec.scatter_reverse_add()


class ParVec_petsc4py(p4pyPETSc.Vec):
    """
    Parallel vector using petsc4py's wrappers for PETSc
    WIP -- This function builds the local to global mapping for the PETSc parallel vectors.  At this
    point it only works when the variables can be interwoven (eg. stablized elements where velocity and
    pressure come from the same space).  We would like to extend this functionality to include finite
    element spaces that cannot be interwoven such as Taylor Hood.
    """
    def __init__(self,array=None,bs=None,n=None,N=None,nghosts=None,subdomain2global=None,blockVecType="simple",ghosts=None,
                                                 proteus2petsc_subdomain=None,
                                                 petsc2proteus_subdomain=None):
        p4pyPETSc.Vec.__init__(self)
        if array is None:
            return#when duplicating for petsc usage
        self.proteus2petsc_subdomain=proteus2petsc_subdomain
        self.petsc2proteus_subdomain=petsc2proteus_subdomain
        blockSize = max(1,bs)
        self.dim_proc = n*blockSize
        self.nghosts = nghosts
        self.blockVecType = blockVecType
        assert self.blockVecType == "simple", "petsc4py wrappers require self.blockVecType=simple"
        self.proteus_array = array
        if nghosts is None:
            if blockVecType == "simple":
                self.createWithArray(array,size=(blockSize*n,blockSize*N),bsize=1)
            else:
                self.createWithArray(array,size=(blockSize*n,blockSize*N),bsize=blockSize)
            self.subdomain2global=subdomain2global
            self.petsc_l2g = None
            self.setUp()
        else:
            assert nghosts >= 0, "The number of ghostnodes must be non-negative"
            assert subdomain2global.shape[0] == (n+nghosts), ("The subdomain2global map is the wrong length n=%i,nghosts=%i,shape=%i \n" % (n,n+nghosts,subdomain2global.shape[0]))
            assert len(array.flat) == (n+nghosts)*blockSize, "%i  != (%i+%i)*%i \n" % (len(array.flat),  n,nghosts,blockSize)
            if blockVecType == "simple":
                if ghosts is None:
                    ghosts = numpy.zeros((blockSize*nghosts),'i')
                    for j in range(blockSize):
                        ghosts[j::blockSize]=subdomain2global[n:]*blockSize+j
                self.createGhostWithArray(ghosts,array,size=(blockSize*n,blockSize*N),bsize=1)
                if blockSize > 1: #have to build in block dofs
                    subdomain2globalTotal = numpy.zeros((blockSize*subdomain2global.shape[0],),'i')
                    for j in range(blockSize):
                        subdomain2globalTotal[j::blockSize]=subdomain2global*blockSize+j
                    self.subdomain2global=subdomain2globalTotal
                else:
                    self.subdomain2global=subdomain2global
            else:
                #TODO need to debug
                ghosts = subdomain2global[n:]
                self.createGhostWithArray(ghosts,array,size=(blockSize*n,blockSize*N),bsize=blockSize)
                self.subdomain2global = subdomain2global
            self.setUp()
            #self.petsc_l2g = p4pyPETSc.LGMap()
            #self.petsc_l2g.create(self.subdomain2global)
            #self.setLGMap(self.petsc_l2g)
        self.setFromOptions()
    def scatter_forward_insert(self):
        if self.proteus2petsc_subdomain is not None:
            self.proteus_array[:] = self.proteus_array[self.petsc2proteus_subdomain]
        self.ghostUpdateBegin(p4pyPETSc.InsertMode.INSERT,p4pyPETSc.ScatterMode.FORWARD)
        self.ghostUpdateEnd(p4pyPETSc.InsertMode.INSERT,p4pyPETSc.ScatterMode.FORWARD)
        if self.proteus2petsc_subdomain is not None:
            self.proteus_array[:] = self.proteus_array[self.proteus2petsc_subdomain]
    def scatter_reverse_add(self):
        if self.proteus2petsc_subdomain is not None:
            self.proteus_array[:] = self.proteus_array[self.petsc2proteus_subdomain]
        self.ghostUpdateBegin(p4pyPETSc.InsertMode.ADD_VALUES,p4pyPETSc.ScatterMode.REVERSE)
        self.ghostUpdateEnd(p4pyPETSc.InsertMode.ADD_VALUES,p4pyPETSc.ScatterMode.REVERSE)
        if self.proteus2petsc_subdomain is not None:
            self.proteus_array[:] = self.proteus_array[self.proteus2petsc_subdomain]

    def save(self, filename):
        """Saves to disk using a PETSc binary viewer.
        """
        _petsc_view(self, filename)

class ParMat_petsc4py(p4pyPETSc.Mat):
    """
    Parallel matrix based on petsc4py's wrappers for PETSc.
    """
    def __init__(self,ghosted_csr_mat=None,par_bs=None,par_n=None,par_N=None,par_nghost=None,subdomain2global=None,blockVecType="simple",pde=None, par_nc=None, par_Nc=None, proteus_jacobian=None, nzval_proteus2petsc=None):
        p4pyPETSc.Mat.__init__(self)
        if ghosted_csr_mat is None:
            return#when duplicating for petsc usage
        self.pde = pde
        if par_nc is None:
            par_nc = par_n
        if par_Nc is None:
            par_Nc = par_N
        self.proteus_jacobian=proteus_jacobian
        self.nzval_proteus2petsc = nzval_proteus2petsc
        self.ghosted_csr_mat=ghosted_csr_mat
        self.blockVecType = blockVecType
        assert self.blockVecType == "simple", "petsc4py wrappers require self.blockVecType=simple"
        self.create(p4pyPETSc.COMM_WORLD)
        self.blockSize = max(1,par_bs)
        if self.blockSize > 1 and blockVecType != "simple":
            ## \todo fix block aij in ParMat_petsc4py
            self.setType('baij')
            self.setSizes([[self.blockSize*par_n,self.blockSize*par_N],[self.blockSize*par_nc,self.blockSize*par_Nc]],bsize=self.blockSize)
            self.setBlockSize(self.blockSize)
            self.subdomain2global = subdomain2global #no need to include extra block dofs?
        else:
            self.setType('aij')
            self.setSizes([[par_n*self.blockSize,par_N*self.blockSize],[par_nc*self.blockSize,par_Nc*self.blockSize]],bsize=1)
            if self.blockSize > 1: #have to build in block dofs
                subdomain2globalTotal = numpy.zeros((self.blockSize*subdomain2global.shape[0],),'i')
                for j in range(self.blockSize):
                    subdomain2globalTotal[j::self.blockSize]=subdomain2global*self.blockSize+j
                self.subdomain2global=subdomain2globalTotal
            else:
                self.subdomain2global=subdomain2global
        from proteus import Comm
        comm = Comm.get()
        logEvent("ParMat_petsc4py comm.rank= %s blockSize = %s par_n= %s par_N=%s par_nghost=%s par_jacobian.getSizes()= %s "
                 % (comm.rank(),self.blockSize,par_n,par_N,par_nghost,self.getSizes()))
        self.csr_rep = ghosted_csr_mat.getCSRrepresentation()
        if self.proteus_jacobian is not None:
            self.proteus_csr_rep = self.proteus_jacobian.getCSRrepresentation()
        if self.blockSize > 1:
            blockOwned = self.blockSize*par_n
            self.csr_rep_local = ghosted_csr_mat.getSubMatCSRrepresentation(0,blockOwned)
        else:
            self.csr_rep_local = ghosted_csr_mat.getSubMatCSRrepresentation(0,par_n)
        self.petsc_l2g = p4pyPETSc.LGMap()
        self.petsc_l2g.create(self.subdomain2global)
        self.setUp()
        self.setLGMap(self.petsc_l2g)
        #
        self.colind_global = self.petsc_l2g.apply(self.csr_rep_local[1]) #prealloc needs global indices
        self.setPreallocationCSR([self.csr_rep_local[0],self.colind_global,self.csr_rep_local[2]])
        self.setFromOptions()

    def save(self, filename):
        """Saves to disk using a PETSc binary viewer.
        """
        _petsc_view(self, filename)

def Vec(n):
    """
    Build a vector of length n (using numpy)

    For example::
    
      >>> Vec(3)
      array([ 0.,  0.,  0.])

    """
    return numpy.zeros((n,),'d')


def Mat(m,n):
    """
    Build an m x n matrix (using numpy)

    For example::

      >>> Mat(2,3)
      array([[ 0.,  0.,  0.],
            [ 0.,  0.,  0.]])

    """
    return numpy.zeros((m,n),'d')


def SparseMatFromDict(nr,nc,aDict):
    """
    Build a nr x nc sparse matrix from a dictionary representation
    """
    import superluWrappers
    indeces = aDict.keys()
    indeces.sort()
    nnz     = len(indeces)
    nzval   = numpy.zeros((nnz,),'d')
    rowptr  = numpy.zeros((nr+1,),'i')
    colind  = numpy.zeros((nnz,),'i')
    i=0
    k=0
    rowptr[i]=0
    for ij in indeces:
        nzval[k] = aDict[ij]
        colind[k] = ij[1]
        if ij[0] > i:
            i += 1
            rowptr[i]=k
        k+=1
    rowptr[i+1] = k
    return (SparseMat(nr,nc,nnz,nzval,colind,rowptr),nzval)


def SparseMat(nr,nc,nnz,nzval,colind,rowptr):
    """ Build a nr x nc sparse matrix from the CSR data structures

    Parameters
    ----------
    nr : int
        The number of rows.
    nc : int
        The number of columns.
    nnz : int
        The number of non-zero matrix entries.
    nzval : numpy array
        Array with non-zero matrix entries.
    colind : numpy array of 32bit integers
        CSR column array.
    rowptr : numpy array of 32bit integers
        CSR row pointer.

    Returns
    -------
    sparse_matrix : :class:`proteus.superluWrappers.SparseMatrix`
        superlu sparse matrix in CSR format.

    Note
    ----
    For the superluWrapper, both the colind and rowptr should use
    32-bit integer data types.
    """
    if (colind.dtype != 'int32' or rowptr.dtype != 'int32'):
        print('ERROR - colind and rowptr must be "int32" numpy arrays for ' \
              'superluWrappers')
        sys.exit(1)
    return superluWrappers.SparseMatrix(nr,nc,nnz,nzval,colind,rowptr)

class SparseMatShell:
    """ Build a parallel matrix shell from CSR data structures.

    Parameters
    ----------
    ghosted_csr_mat: :class: `proteus.superluWrappers.SparseMatrix`
    """
    def __init__(self,ghosted_csr_mat):
        self.ghosted_csr_mat=ghosted_csr_mat
        self.par_b = None
        self.xGhosted = None
        self.yGhosted = None
    def create(self, A):
        pass
    def mult(self, A, x, y):
        assert self.par_b is not None, "The parallel RHS vector par_b must be " \
                            "initialized before using the mult function"
        logEvent("Using SparseMatShell in LinearSolver matrix multiply")
        if self.xGhosted is None:
            self.xGhosted = self.par_b.duplicate()
            self.yGhosted = self.par_b.duplicate()
        self.xGhosted.setArray(x.getArray())
        self.xGhosted.ghostUpdateBegin(p4pyPETSc.InsertMode.INSERT,p4pyPETSc.ScatterMode.FORWARD)
        self.xGhosted.ghostUpdateEnd(p4pyPETSc.InsertMode.INSERT,p4pyPETSc.ScatterMode.FORWARD)
        self.yGhosted.zeroEntries()
        with self.xGhosted.localForm() as xlf, self.yGhosted.localForm() as ylf:
            self.ghosted_csr_mat.matvec(xlf.getArray(),ylf.getArray())
        y.setArray(self.yGhosted.getArray())

class OperatorShell:
    """ A base class for operator shells """
    def __init__(self):
        pass
    def create(self,A):
        pass

class ProductOperatorShell(OperatorShell):
    """ A base class for shell operators that apply multiplcation. 
    
    Operators derived from this class should have working multiplication
    functions.
    """
    def __init__(self):
        pass
    def mult(self, A, x, y):
        raise NotImplementedError('You need to define a multiply' \
                                  'function for your shell')

class InvOperatorShell(OperatorShell):
    """ A base class for inverse operator shells 
    
    Operators derived from this class should have working apply
    functions.
    """
    def __init__(self):
        pass
    def apply(self, A, x, y):
        raise NotImplementedError('You need to define an apply' \
                                  'function for your shell')

    def _create_tmp_vec(self,size):
        """ Creates an empty vector of given size. 
        
        Arguments
        ---------
        size : int
            Size of the temporary vector.

        Returns
        -------
        vec : PETSc vector
        """
        tmp = p4pyPETSc.Vec().create()
        tmp.setType('mpi')
        tmp.setSizes(size)
        return tmp

    def _create_constant_nullspace(self):
        """Initialize a constant null space. """
        self.const_null_space = p4pyPETSc.NullSpace().create(comm=p4pyPETSc.COMM_WORLD,
                                                             vectors = (),
                                                             constant = True)
    
    def _converged_trueRes(self,ksp,its,rnorm):
        """ Function handle to feed to ksp's setConvergenceTest  """
        ksp.buildResidual(self.r_work)
        truenorm = self.r_work.norm()
        if its == 0:
            self.rnorm0 = truenorm
            # ARB - Leaving these log events in for future debugging purposes.
            # logEvent("NumericalAnalytics KSP_LSC_LaplaceResidual: %12.5e" %(truenorm) )
            # logEvent("NumericalAnalytics KSP_LSC_LaplaceResidual(relative): %12.5e" %(truenorm / self.rnorm0) )
            # logEvent("        KSP it %i norm(r) = %e  norm(r)/|b| = %e ; atol=%e rtol=%e " % (its,
            #                                                                                   truenorm,
            #                                                                                   (truenorm/ self.rnorm0),
            #                                                                                   ksp.atol,
            #                                                                                   ksp.rtol))
            return False
        else:
            # ARB - Leaving these log events in for future debugging purposes.
            # logEvent("NumericalAnalytics KSP_LSC_LaplaceResidual: %12.5e" %(truenorm) )
            # logEvent("NumericalAnalytics KSP_LSC_LaplaceResidual(relative): %12.5e" %(truenorm / self.rnorm0) )
            # logEvent("        KSP it %i norm(r) = %e  norm(r)/|b| = %e ; atol=%e rtol=%e " % (its,
            #                                                                                   truenorm,
            #                                                                                   (truenorm/ self.rnorm0),
            #                                                                                   ksp.atol,
            #                                                                                   ksp.rtol))
            if truenorm < self.rnorm0*ksp.rtol:
                return p4pyPETSc.KSP.ConvergedReason.CONVERGED_RTOL
            if truenorm < ksp.atol:
                return p4pyPETSc.KSP.ConvergedReason.CONVERGED_ATOL
        return False

class MatrixShell(ProductOperatorShell):
    """ A shell class for a matrix. """
    def __init__(self,A):
        """
        Specifies a basic matrix shell.

        Parameters
        ----------
        A : matrix
            A petsc4py matrix object
        """
        self.A = A
    def mult(self,A,x,y):
        """
        Multiply the matrix and x.

        Parameters
        ----------
        A : matrix
            Dummy place holder for PETSc compatibility
        x : vector

        Returns
        -------
        y : vector
        """
        self.A.mult(x,y)

class MatrixInvShell(InvOperatorShell):
    """ A PETSc shell class for a inverse operator. """
    def __init__(self, A):
        """ Initializes operators and solvers for inverse operator.

        Parameters
        ----------
        A : PETSc matrix
            This is the matrix object used to construct the inverse.
        """
        self.A = A
        self.ksp = p4pyPETSc.KSP().create()
        self.ksp.setOperators(self.A,self.A)
        self.ksp.setType('preonly')
        self.ksp.pc.setType('lu')
        self.ksp.setUp()
    def apply(self,A,x,y):
        """ Apply the inverse pressure mass matrix.

        Parameters
        ----------
        A : matrix
            Dummy place holder for PETSc compatibility
        x : vector

        Returns
        -------
        y : vector
        """
        self.ksp.solve(x,y)
        
def l2Norm(x):
    """
    Compute the parallel :math:`l_2` norm
    """
    return math.sqrt(flcbdfWrappers.globalSum(numpy.dot(x,x)))


def l1Norm(x):
    """
    Compute the parallel :math:`l_1` norm
    
    The :math:`l_1` norm of a vector :math:`\mathbf{x} \in
    \mathbb{R}^n` is
    
    .. math:: 
    
       \| \mathbf{x} \|_{1} = \sum_{i=0} |x_i|
    
    If Python is running in parallel, then the sum is over all
    dimensions on all processors so that the input must not contain
    "ghost" entries.
    
    This implemtation works for a distributed array with no ghost
    components (each component must be on a single processor).
    
    :param x: numpy array of length n
    :return: float
    """
    return flcbdfWrappers.globalSum(numpy.sum(numpy.abs(x)))


def lInfNorm(x):
    """
    Compute the parallel :math:`l_{\infty}` norm

    The :math:`l_{\infty}` norm of a vector :math:`\mathbf{x} \in
    \mathbb{R}^n` is

    .. math::

       \|x\|_{\infty} = \max_i |x_i|
       
    This implemtation works for a distributed array with no ghost
    components (each component must be on a single processor).
    
    :param x: numpy array of length n
    :return: float
    """
    return flcbdfWrappers.globalMax(numpy.linalg.norm(x,numpy.inf))


def wDot(x,y,h):
    """
    Compute the parallel weighted dot product of vectors x and y using
    weight vector h.
    
    The weighted dot product is defined for a weight vector
    :math:`\mathbf{h}` as

    .. math:: 

       (\mathbf{x},\mathbf{y})_h = \sum_{i} h_{i} x_{i} y_{i}
    
    All weight vector components should be positive.

    :param x,y,h: numpy arrays for vectors and weight 
    :return: the weighted dot product
    """
    return flcbdfWrappers.globalSum(numpy.sum(x*y*h))

def wl2Norm(x,h):
    """
    Compute the parallel weighted l_2 norm with weight h
    """
    return math.sqrt(flcbdfWrappers.globalSum(wDot(x,x,h)))


def wl1Norm(x,h):
    """
    Compute the parallel weighted l_1 norm with weight h
    """
    return flcbdfWrappers.globalSum(numpy.sum(numpy.abs(h*x)))


def wlInfNorm(x,h):
    """
    Compute the parallel weighted l_{\infty} norm with weight h
    """
    return flcbdfWrappers.globalMax(numpy.linalg.norm(h*x,numpy.inf))

def energyDot(x,y,A):
    """
    Compute the "energy" dot product x^t A y (not parallel)
    """
    return numpy.dot(numpy.dot(x,A),y)

def energyNorm(x,A):
    """
    Compute the "energy" norm x^t A x (not parallel)
    """
    return math.sqrt(energyDot(x,x,A))

def l2NormAvg(x):
    """
    Compute the arithmetic averaged l_2 norm (root mean squared norm)
    """
    scale = 1.0/flcbdfWrappers.globalSum(len(x.flat))
    return math.sqrt(scale*flcbdfWrappers.globalSum(numpy.dot(x,x)))


rmsNorm = l2NormAvg


def l2Norm_local(x):
    """
    Compute the l_2 norm for just local (processor) system  (not parallel)
    """
    return math.sqrt(numpy.dot(x,x))


class WeightedNorm:
    """
    Compute the weighted norm for time step control (not currently parallel)
    """
    def __init__(self,shape,atol,rtol):
        self.shape = shape
        self.dim = sum(self.shape)
        self.atol= atol
        self.rtol= rtol
        self.weight = numpy.ones(shape,'d')
        self.tmp    = numpy.ones(shape,'d')
    def setWeight(self,y):
        self.weight[:] = numpy.absolute(y)
        self.weight   *= self.rtol
        self.weight   += self.atol
    def norm(self,y,type):
        self.tmp[:] = y
        self.tmp /= self.weight
        value = numpy.linalg.norm(self.tmp.flat,type)
        return value/self.dim


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    

# def test_MGV():
#     n=2**8 + 1
#     h =1.0/(n-1.0)
#     freq=10
#     u = numpy.random.uniform(0,1,(n))
#     u[0]=0.0
#     u[n-1]=0.0
#     x = numpy.arange(0,1.0+h,h)
#     AList=[]
#     N=n
#     pList=[]
#     rList=[]
#     resList=[]
#     while N >= 3:
#         resList.append(Vec(N-2))
#         A = dict()#SparseMat(N-2,N-2,3*(N-2),sym=True)
#         H = 1.0/(N-1.0)
#         #beginAssembly(A)
#         for i in range(N-2):
#             A[(i,i)] = 2.0/H**2
#             if i > 0:
#                 A[(i,i-1)] = -1.0/H**2
#             if i < N-3:
#                 A[(i,i+1)] = -1.0/H**2
#         #endAssembly(A)
#         AList.append(SparseMatFromDict(N-2,N-2,A)[0])
#         cN = (N - 1)/2 + 1
#         r = dict()#SparseMat(cN-2,N-2,3*(N-2))
#         p = dict()#SparseMat(N-2,cN-2,3*(N-2))
#         for i in range(cN-2):
#             r[(i,2*i)]   = 1.0/4.0
#             r[(i,2*i+1)] = 2.0/4.0
#             r[(i,2*i+2)] = 1.0/4.0
#             p[(2*i,i)] = 1.0/2.0
#             p[(2*i+1,i)]= 2.0/2.0
#             p[(2*i+2,i)]= 1.0/2.0
#         #r.to_csr()
#         print cN-2,N-2,r.keys()
#         if cN-2 > 0:
#             rList.append(SparseMatFromDict(cN-2,N-2,r)[0])
#         else:
#             rList.append(None)
#         #p.to_csr()
#         pList.append(SparseMatFromDict(N-2,cN-2,p)[0])
#         N = cN
#     class Jacobi:
#         def __init__(self,A):
#             self.A=A
#             self.n=A.shape[0]
#             self.M=Vec(self.n)
#             for i in range(self.n):
#                 self.M[i]=1.0/A[i,i]
#             self.res=Vec(self.n)
#             self.dx=Vec(self.n)
#         def apply(self,w,jits,b,x):
#             self.A.matvec(x,self.res)
#             self.res-=b
#             for it in range(jits):
#                 self.dx[:] = self.M*self.res
#                 self.dx*=w
#                 x -= self.dx
#                 self.A.matvec(x,self.res)
#                 self.res -= b
#     jacobiList=[]
#     for A in AList:
#         jacobiList.append(Jacobi(A))
#     jits = 3
#     w = 2.0/3.0
#     class MGV:
#         def __init__(self,smootherList,AList,pList,rList,resList):
#             self.AList = AList
#             self.pList = pList
#             self.rList = rList
#             self.resList = resList
#             self.xList=[]
#             self.vList=[]
#             self.bList=[]
#             self.gpList=[]
#             for res in resList:
#                 self.xList.append(Vec(len(res)))
#                 self.vList.append(Vec(len(res)))
#                 self.bList.append(Vec(len(res)))
#             self.smootherList = smootherList

#         def apply(self,w,nsPre,nsPost,level,b,x):
#             logEvent("Level = "+`level`)
#             if level == len(self.AList)-1:
#                 self.smootherList[level].apply(1.0,1,b,x)
#             else:
#                 #smooth
#                 self.smootherList[level].apply(w,nsPre,b,x)
#                 #restrict the defect
#                 self.rList[level].matvec(self.smootherList[level].res,self.bList[level+1])
#                 #V-cycle on the error equation
#                 self.xList[level+1][:]=0.0
#                 self.apply(w,nsPre,nsPost,level+1,self.bList[level+1],self.xList[level+1])
#                 #prolong
#                 self.pList[level].matvec(self.xList[level+1],self.vList[level])
#                 #correct
#                 x-=self.vList[level]
#                 #smooth
#                 self.smootherList[level].apply(w,nsPost,b,x)
#                 self.resList[level][:]=self.smootherList[level].res
#     mgv = MGV(jacobiList,AList,pList,rList,resList)
#     rnorm=1.0
#     mgits = 0
#     while rnorm > 1.0e-10 and mgits < 20:
#         mgits +=1
#         mgv.apply(w,jits,jits,0,f[1:n-1],u[1:n-1])
#         rnorm = l2Norm(resList[0])

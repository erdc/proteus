"""
Tools for n-dimensional linear algebra

Vectors are just numpy arrays, as are dense matrices. Sparse matrices
are CSR matrices. Parallel vector and matrix are built on top of those
representations using PETSc.

\todo LinearAlgebraTools: make better use of numpy.linalg and petsc4py to provide the needed functionality and improve test suite
"""
import numpy
import math
from .superluWrappers import *
from . import flcbdfWrappers
from Profiling import logEvent
#PETSc import, forces comm init if not already done
from petsc4py import PETSc as p4pyPETSc
from . import Comm
Comm.set_isInitialized()
#end PETSc import


class ParVec:
    """
    A parallel vector built on top of daetk's wrappers for petsc
    """
    def __init__(self,array,blockSize,n,N,nghosts=None,subdomain2global=None,blockVecType="simple"):#"block"
        import flcbdfWrappers
        self.dim_proc=n*blockSize
        if nghosts==None:
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
    """
    def __init__(self,array=None,bs=None,n=None,N=None,nghosts=None,subdomain2global=None,blockVecType="simple"):
        if array == None:
            return#cek hack, don't know why init gets called by PETSc.Vec duplicate function
        p4pyPETSc.Vec.__init__(self)
        blockSize = max(1,bs)
        self.dim_proc = n*blockSize
        self.nghosts = nghosts
        self.blockVecType = blockVecType
        assert self.blockVecType == "simple", "petsc4py wrappers require self.blockVecType=simple"
        self.proteus_array = array
        if nghosts == None:
            if blockVecType == "simple":
                self.createWithArray(array,size=(blockSize*n,blockSize*N),bsize=1)
            else:
                self.createWithArray(array,size=(blockSize*n,blockSize*N),bsize=blockSize)
            self.subdomain2global=subdomain2global
            self.petsc_l2g = None
        else:
            assert nghosts >= 0, "The number of ghostnodes must be non-negative"
            assert subdomain2global.shape[0] == (n+nghosts), ("The subdomain2global map is the wrong length n=%i,nghosts=%i,shape=%i \n" % (n,n+nghosts,subdomain2global.shape[0]))
            assert len(array.flat) == (n+nghosts)*blockSize
            if blockVecType == "simple":
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
                self.petsc_l2g = p4pyPETSc.LGMap()
                self.petsc_l2g.create(self.subdomain2global)
                self.setLGMap(self.petsc_l2g)

            else:
                #TODO need to debug
                ghosts = subdomain2global[n:]
                self.createGhostWithArray(ghosts,array,size=(blockSize*n,blockSize*N),bsize=blockSize)
                self.subdomain2global = subdomain2global
                self.petsc_l2g = p4pyPETSc.LGMap()
                self.petsc_l2g.create(self.subdomain2global)
                self.setLGMap(self.petsc_l2g)
        self.setFromOptions()
    def scatter_forward_insert(self):
        self.ghostUpdateBegin(p4pyPETSc.InsertMode.INSERT,p4pyPETSc.ScatterMode.FORWARD)
        self.ghostUpdateEnd(p4pyPETSc.InsertMode.INSERT,p4pyPETSc.ScatterMode.FORWARD)
    def scatter_reverse_add(self):
        self.ghostUpdateBegin(p4pyPETSc.InsertMode.ADD_VALUES,p4pyPETSc.ScatterMode.REVERSE)
        self.ghostUpdateEnd(p4pyPETSc.InsertMode.ADD_VALUES,p4pyPETSc.ScatterMode.REVERSE)


class ParMat_petsc4py(p4pyPETSc.Mat):
    """
    Parallel matrix based on petsc4py's wrappers for PETSc.
    """
    def __init__(self,ghosted_csr_mat,par_bs,par_n,par_N,par_nghost,subdomain2global,blockVecType="simple"):
        p4pyPETSc.Mat.__init__(self)
        self.ghosted_csr_mat=ghosted_csr_mat
        self.blockVecType = blockVecType
        assert self.blockVecType == "simple", "petsc4py wrappers require self.blockVecType=simple"
        self.create(p4pyPETSc.COMM_WORLD)
        blockSize = max(1,par_bs)
        if blockSize >= 1 and blockVecType != "simple":
            ## \todo fix block aij in ParMat_petsc4py
            self.setType('baij')
            self.setSizes([[blockSize*par_n,blockSize*par_N],[blockSize*par_n,blockSize*par_N]],bsize=blockSize)
            self.setBlockSize(blockSize)
            self.subdomain2global = subdomain2global #no need to include extra block dofs?
        else:
            self.setType('aij')
            self.setSizes([[par_n*blockSize,par_N*blockSize],[par_n*blockSize,par_N*blockSize]],bsize=1)
            if blockSize > 1: #have to build in block dofs
                subdomain2globalTotal = numpy.zeros((blockSize*subdomain2global.shape[0],),'i')
                for j in range(blockSize):
                    subdomain2globalTotal[j::blockSize]=subdomain2global*blockSize+j
                self.subdomain2global=subdomain2globalTotal
            else:
                self.subdomain2global=subdomain2global
        import Comm
        comm = Comm.get()
        logEvent("ParMat_petsc4py comm.rank= %s blockSize = %s par_n= %s par_N=%s par_nghost=%s par_jacobian.getSizes()= %s "
                 % (comm.rank(),blockSize,par_n,par_N,par_nghost,self.getSizes()))
        self.csr_rep = ghosted_csr_mat.getCSRrepresentation()
        blockOwned = blockSize*par_n
        self.csr_rep_owned = ghosted_csr_mat.getSubMatCSRrepresentation(0,blockOwned)
        self.petsc_l2g = p4pyPETSc.LGMap()
        self.petsc_l2g.create(self.subdomain2global)
        self.colind_global = self.petsc_l2g.apply(self.csr_rep_owned[1]) #prealloc needs global indices
        self.setPreallocationCSR([self.csr_rep_owned[0],self.colind_global,self.csr_rep_owned[2]])
        self.setUp()
        self.setLGMap(self.petsc_l2g)
        self.setFromOptions()


def Vec(n):
    """
    Build a vector of length n (using numpy)

    For example:
    >>> Vec(3)
    array([ 0.,  0.,  0.])
    """
    return numpy.zeros((n,),'d')


def Mat(m,n):
    """
    Build an m x n matrix (using numpy)

    For example:
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
    """
    Build a nr x nc sparse matrix from the CSR data structures
    """
    import superluWrappers
    return superluWrappers.SparseMatrix(nr,nc,nnz,nzval,colind,rowptr)


class SparseMatShell:
    """
    Build a parallel matrix shell using the subdomain CSR data structures (must have overlapping subdomains)
    """
    def __init__(self,ghosted_csr_mat):
        self.ghosted_csr_mat=ghosted_csr_mat
        self.par_b = None
        self.xGhosted = None
        self.yGhosted = None
    def create(self, A):
        pass
    def mult(self, A, x, y):
        logEvent("Using SparseMatShell in LinearSolver matrix multiply")
        if self.xGhosted == None:
            self.xGhosted = self.par_b.duplicate()
            self.yGhosted = self.par_b.duplicate()
        self.xGhosted.setArray(x.getArray())
        self.xGhosted.ghostUpdateBegin(p4pyPETSc.InsertMode.INSERT,p4pyPETSc.ScatterMode.FORWARD)
        self.xGhosted.ghostUpdateEnd(p4pyPETSc.InsertMode.INSERT,p4pyPETSc.ScatterMode.FORWARD)
        self.yGhosted.zeroEntries()
        with self.xGhosted.localForm() as xlf, self.yGhosted.localForm() as ylf:
            self.ghosted_csr_mat.matvec(xlf,ylf)
        y.setArray(self.yGhosted.getArray())


def l2Norm(x):
    """
    Compute the parallel l_2 norm
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
    return scale*math.sqrt(flcbdfWrappers.globalSum(numpy.dot(x,x)))


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
    

def test_MGV():
    '''Non-working function (fix imports below)'''
    #    from LinearAlgebraTools import *
    #    import Gnuplot
    #    from Gnuplot import *
    #    from math import *
    #    from RandomArray import *
    gf = Gnuplot.Gnuplot()
    gf("set terminal x11")
    ginit = Gnuplot.Gnuplot()
    ginit("set terminal x11")
    gsol = Gnuplot.Gnuplot()
    gsol("set terminal x11")
    gres = Gnuplot.Gnuplot()
    gres("set terminal x11")
    n=2**8 + 1
    h =1.0/(n-1.0)
    freq=10
    u = uniform(0,1,(n))
    u[0]=0.0
    u[n-1]=0.0
    x = numpy.arange(0,1.0+h,h)
    f = Gnuplot.Gnuplot()
    f("set terminal x11")
    f = (freq*2*pi)**2*numpy.sin(freq*2*pi*x)
    AList=[]
    N=n
    pList=[]
    rList=[]
    resList=[]
    while N >= 3:
        resList.append(Vec(N-2))
        A = SparseMat(N-2,N-2,3*(N-2),sym=True)
        H = 1.0/(N-1.0)
        beginAssembly(A)
        for i in range(N-2):
            A[i,i] = 2.0/H**2
            if i > 0:
                A[i,i-1] = -1.0/H**2
            if i < N-3:
                A[i,i+1] = -1.0/H**2
        endAssembly(A)
        AList.append(A)
        cN = (N - 1)/2 + 1
        r = SparseMat(cN-2,N-2,3*(N-2))
        p = SparseMat(N-2,cN-2,3*(N-2))
        for i in range(cN-2):
            r[i,2*i]   = 1.0/4.0
            r[i,2*i+1] = 2.0/4.0
            r[i,2*i+2] = 1.0/4.0
            p[2*i,i] = 1.0/2.0
            p[2*i+1,i]= 2.0/2.0
            p[2*i+2,i]= 1.0/2.0
        r.to_csr()
        rList.append(r)
        p.to_csr()
        pList.append(p)
        N = cN
    class Jacobi:
        def __init__(self,A):
            self.A=A
            self.n=A.shape[0]
            self.M=Vec(self.n)
            for i in range(self.n):
                self.M[i]=1.0/A[i,i]
            self.res=Vec(self.n)
            self.dx=Vec(self.n)
        def apply(self,w,jits,b,x):
            self.A.matvec(x,self.res)
            self.res-=b
            for it in range(jits):
                self.dx[:] = self.M*self.res
                self.dx*=w
                x -= self.dx
                self.A.matvec(x,self.res)
                self.res -= b
    jacobiList=[]
    for A in AList:
        jacobiList.append(Jacobi(A))
    jits = 3
    w = 2.0/3.0
    class MGV:
        def __init__(self,smootherList,AList,pList,rList,resList):
            self.AList = AList
            self.pList = pList
            self.rList = rList
            self.resList = resList
            self.xList=[]
            self.vList=[]
            self.bList=[]
            self.gpList=[]
            for res in resList:
                self.xList.append(Vec(len(res)))
                self.vList.append(Vec(len(res)))
                self.bList.append(Vec(len(res)))
                self.gpList.append(Gnuplot.Gnuplot(debug=1)("set terminal x11"))
            self.smootherList = smootherList

        def apply(self,w,nsPre,nsPost,level,b,x):
            logEvent("Level = "+`level`)
            if level == len(self.AList)-1:
                self.smootherList[level].apply(1.0,1,b,x)
            else:
                #smooth
                self.smootherList[level].apply(w,nsPre,b,x)
                #restrict the defect
                self.rList[level].matvec(self.smootherList[level].res,self.bList[level+1])
                #V-cycle on the error equation
                self.xList[level+1][:]=0.0
                self.apply(w,nsPre,nsPost,level+1,self.bList[level+1],self.xList[level+1])
                #self.gpList[level].plot(Gnuplot.Data(self.smootherList[level+1].res,title='residual'))
                #prolong
                self.pList[level].matvec(self.xList[level+1],self.vList[level])
                #correct
                x-=self.vList[level]
                #smooth
                self.smootherList[level].apply(w,nsPost,b,x)
                self.resList[level][:]=self.smootherList[level].res
    mgv = MGV(jacobiList,AList,pList,rList,resList)
    rnorm=1.0
    mgits = 0
    while rnorm > 1.0e-10 and mgits < 20:
        mgits +=1
        mgv.apply(w,jits,jits,0,f[1:n-1],u[1:n-1])
        rnorm = l2Norm(resList[0])
    gsol.plot(Gnuplot.Data(x,u,title='numerical solution'),
              Gnuplot.Data(x,numpy.sin(freq*2*pi*x),title='exact solution'))
    #gres.plot(Gnuplot.Data(x[1:n-1],mgv.smootherList[0].res,title='final residual'))
    raw_input('Please press return to continue... \n')

"""
Tools for n-dimensional linear algebra

Vectors will just be numpy arrays
as will dense matrices. 
"""
import numpy
from superluWrappers import *
import flcbdfWrappers
from Profiling import logEvent

class ParVec:
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
			assert len(array.flat) == (n+nghosts)*blockSize
			if blockVecType=="simple":
				self.cparVec=flcbdfWrappers.ParVec(blockSize,n,N,nghosts,subdomain2global,array,0)
			else:
				self.cparVec=flcbdfWrappers.ParVec(blockSize,n,N,nghosts,subdomain2global,array,1)
		self.nghosts = nghosts
        def scatter_forward_insert(self):
                self.cparVec.scatter_forward_insert()
        def scatter_reverse_add(self):
                self.cparVec.scatter_reverse_add()

class ParVec2:
        def __init__(self,array,blockSize,n,N,subdomain2global):
                import flcbdfWrappers
		self.arrayAll = array
		self.arrayPetsc = numpy.zeros((blockSize*n,),'d')
		self.subdomain2global = subdomain2global
		self.blockSize=blockSize
                self.dim_proc=n*blockSize
		self.cparVec=flcbdfWrappers.ParVec(blockSize,n,N,-1,None,self.arrayPetsc,0)
        def scatter_forward_insert(self):
                self.cparVec.scatter_forward_insertAll(self.subdomain2global,self.arrayAll)
        def scatter_reverse_add(self):
                self.cparVec.scatter_reverse_addAll(self.subdomain2global,self.arrayAll)

from petsc4py import PETSc
class ParVec_petsc4py(PETSc.Vec):
	"""
	wrapper for petsc4py's interface to PETSc Vec
	TODO:

	"""
	def __init__(self,array=None,bs=None,n=None,N=None,nghosts=None,subdomain2global=None,blockVecType="simple"):
		if array == None:
			return#cek hack, don't know why init gets called by PETSc.Vec duplicate function
		PETSc.Vec.__init__(self)
		blockSize = max(1,bs)
		self.dim_proc = n*blockSize
		self.nghosts = nghosts
		self.blockVecType = blockVecType
		assert self.blockVecType == "simple", "petsc4py wrappers require self.blockVecType=simple"
		self.pyadh_array = array
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
				#keep for debugging for a little while
				#ghosts2 = numpy.zeros((blockSize*nghosts),'i')
				#for i in range(nghosts):
				#	for j in range(blockSize):
				#		ghosts2[i*blockSize+j]=subdomain2global[n+i]*blockSize+j
				#if len(ghosts2) > 0 or len(ghosts) > 0:
				#	assert numpy.max(numpy.absolute(ghosts2-ghosts)) < 1.0e-8
				
				self.createGhostWithArray(ghosts,array,size=(blockSize*n,blockSize*N),bsize=1)
				if blockSize > 1: #have to build in block dofs
					subdomain2globalTotal = numpy.zeros((blockSize*subdomain2global.shape[0],),'i')
					for j in range(blockSize):
						subdomain2globalTotal[j::blockSize]=subdomain2global*blockSize+j
					#keep for debugging for a minute
					#subdomain2globalTotal2= numpy.zeros((blockSize*subdomain2global.shape[0],),'i')
					#for i in range(subdomain2global.shape[0]):
					#	for j in range(blockSize):
					#		subdomain2globalTotal2[i*blockSize+j]=subdomain2global[i]*blockSize+j
					#if len(subdomain2globalTotal) > 0 or len(subdomain2globalTotal2) > 0:
					#	assert numpy.max(numpy.absolute(subdomain2globalTotal2-subdomain2globalTotal)) < 1.0e-8
					self.subdomain2global=subdomain2globalTotal
				else:
					self.subdomain2global=subdomain2global
				self.petsc_l2g = PETSc.LGMap()
				self.petsc_l2g.create(self.subdomain2global)
				self.setLGMap(self.petsc_l2g)

			else:
				#TODO need to debug
				ghosts = subdomain2global[n:]
				self.createGhostWithArray(ghosts,array,size=(blockSize*n,blockSize*N),bsize=blockSize)
				self.subdomain2global = subdomain2global
				self.petsc_l2g = PETSc.LGMap()
				self.petsc_l2g.create(self.subdomain2global)
				self.setLGMap(self.petsc_l2g)
		#
		self.setFromOptions()
	#
        def scatter_forward_insert(self):
		#mwf debug
		#import Comm
		#comm = Comm.get()
		#print "rank= %s ParVec3 into scatter_forward view " % (comm.rank()) 
		#self.view()
		addValues = False; scatterReverse=False
                self.ghostUpdateBegin(addValues,scatterReverse)
		self.ghostUpdateEnd(addValues,scatterReverse)
		#print "rank= %s ParVec3 outof scatter_forward " % (comm.rank())
		#self.view()
		#comm.barrier()
        def scatter_reverse_add(self):
		addValues = True; scatterReverse=True
                self.ghostUpdateBegin(addValues,scatterReverse)
		self.ghostUpdateEnd(addValues,scatterReverse)

			
class ParMat_petsc4py(PETSc.Mat):
	"""
	wrapper for petsc4py interface to PETSc Mat
	TODO:
	
	"""
	def __init__(self,ghosted_csr_mat,par_bs,par_n,par_N,par_nghost,subdomain2global,blockVecType="simple"):
		PETSc.Mat.__init__(self)
		self.ghosted_csr_mat=ghosted_csr_mat
		self.blockVecType = blockVecType
		assert self.blockVecType == "simple", "petsc4py wrappers require self.blockVecType=simple"
		self.create(PETSc.COMM_WORLD)
		blockSize = max(1,par_bs)
		if blockSize >= 1 and blockVecType != "simple":
			#need to debug
			self.setType('baij')
			self.setSizes([[blockSize*par_n,blockSize*par_N],[blockSize*par_n,blockSize*par_N]],bsize=blockSize)
			self.setBlockSize(blockSize)
			self.subdomain2global = subdomain2global #no need to include extra block dofs?
		else:
			self.setType('aij')
			self.setSizes([[par_n*blockSize,par_N*blockSize],[par_n*blockSize,par_N*blockSize]],bsize=blockSize)
			if blockSize > 1: #have to build in block dofs
				subdomain2globalTotal = numpy.zeros((blockSize*subdomain2global.shape[0],),'i')
				for j in range(blockSize):
					subdomain2globalTotal[j::blockSize]=subdomain2global*blockSize+j
				#keep for debugging for a minute
				#subdomain2globalTotal2= numpy.zeros((blockSize*subdomain2global.shape[0],),'i')
				#for i in range(subdomain2global.shape[0]):
				#	for j in range(blockSize):
				#		subdomain2globalTotal2[i*blockSize+j]=subdomain2global[i]*blockSize+j
				#if len(subdomain2globalTotal) > 0 or len(subdomain2globalTotal2) > 0:
				#	assert numpy.max(numpy.absolute(subdomain2globalTotal2-subdomain2globalTotal)) < 1.0e-8
				self.subdomain2global=subdomain2globalTotal
			else:
				self.subdomain2global=subdomain2global
		#
		import Comm
		comm = Comm.get()
		print "ParMat_petsc4py comm.rank= %s blockSize = %s par_n= %s par_N=%s par_nghost=%s par_jacobian.getSizes()= %s " % (comm.rank(),blockSize,
																      par_n,par_N,par_nghost,
																      self.getSizes())
		
		#includes overlap
		self.csr_rep = ghosted_csr_mat.getCSRrepresentation()
		blockOwned = blockSize*par_n
		self.csr_rep_owned = ghosted_csr_mat.getSubMatCSRrepresentation(0,blockOwned)
		self.petsc_l2g = PETSc.LGMap()
		self.petsc_l2g.create(self.subdomain2global)
		self.setLGMap(self.petsc_l2g)
		self.setFromOptions()
	
def Vec(n):
        return numpy.zeros((n,),'d')

def Mat(n,m):
        return numpy.zeros((n,m),'d')

def SparseMatFromDict(nr,nc,aDict):
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
##         print 'nr',nr
##         print 'nc',nc
##         print 'nnz',nnz
##         print 'rowptr',rowptr
##         print 'colind',colind
##         print 'nzval',nzval
        return (SparseMat(nr,nc,nnz,nzval,colind,rowptr),nzval)
def SparseMat(nr,nc,nnz,nzval,colind,rowptr):
        import superluWrappers
        return superluWrappers.SparseMatrix(nr,nc,nnz,nzval,colind,rowptr)

class SparseMatShell:
	from petsc4py import PETSc
	def __init__(self,ghosted_csr_mat):
		self.ghosted_csr_mat=ghosted_csr_mat
		self.par_b = None
		self.xGhosted = None
		self.yGhosted = None
	#cek define matshell
	def create(self, A):
		pass
	def mult(self, A, x, y):
		logEvent("Using SparseMatShell in LinearSolver matrix multiply")
		if self.xGhosted == None:
			self.xGhosted = self.par_b.duplicate()
			self.yGhosted = self.par_b.duplicate()
		self.xGhosted.setArray(x.getArray())
		self.xGhosted.ghostUpdateBegin(PETSc.InsertMode.INSERT,PETSc.ScatterMode.FORWARD)
		self.xGhosted.ghostUpdateEnd(PETSc.InsertMode.INSERT,PETSc.ScatterMode.FORWARD)
		self.yGhosted.zeroEntries()
		self.ghosted_csr_mat.matvec(self.xGhosted.getLocalForm().getArray(),self.yGhosted.getLocalForm().getArray())
		y.setArray(self.yGhosted.getArray())

## def SparseMat_old(n,m,storage):
##         return spmatrix.ll_mat(n,m,storage)

def l2Norm(x):
        return numpy.sqrt(flcbdfWrappers.globalSum(numpy.dot(x,x)))

def l1Norm(x):
        return flcbdfWrappers.globalSum(numpy.sum(numpy.abs(x)))

def lInfNorm(x):
        return flcbdfWrappers.globalMax(numpy.linalg.norm(x,numpy.inf))

def wDot(x,y,h):
        return flcbdfWrappers.globalSum(numpy.sum(x*y*h))

def wl2Norm(x,h):
        return numpy.sqrt(flcbdfWrappers.globalSum(wDot(x,x,h)))

def wl1Norm(x,h):
        return numpy.sum(flcbdfWrappers.globalSum(numpy.abs(h*x)))

def wlInfNorm(x,h):
        return flcbdfWrappers.globalMax(numpy.max(h*numpy.abs(x)))

def energyDot(x,y,A):
        flcbdfWrappers.globalSum(numpy.dot(A*x,y))

def energyNorm(x,A):
        numpy.sqrt(energyDot(x,x,A))

def l2NormAvg(x):
	#mwf what if scale by number of unknowns?
	scale = 1.0/flcbdfWrappers.globalSum(len(x.flat))
        return scale*numpy.sqrt(flcbdfWrappers.globalSum(numpy.dot(x,x)))
def l2Norm_local(x):
	"""
	l2 norm for just local (processor) system 
	"""
	return numpy.sqrt(numpy.dot(x,x))
class WeightedNorm:
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
## @}

if __name__ == '__main__':
        from LinearAlgebraTools import *
        import Gnuplot
        from Gnuplot import *
        from math import *
        from RandomArray import *
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
  

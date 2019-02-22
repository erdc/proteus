# A type of -*- python -*- file
import numpy as np
cimport numpy as np

cdef class FLCBDF_integrator:
    def choose_dt(self, double t, double tout):
        return self.flcbdf.chooseDT(t,tout)

    def set_dt(self, double DT):
        self.flcbdf.setDT(DT)
        return DT

    def set_order(self, int k):
        self.flcbdf.useFixedOrder(k)
        return k

    def initialize_dt(self,
                       double t0,
                       double tOut,
                       np.ndarray y,
                       np.ndarray yPrime):
        cdef Vec *yVec = new Vec(REF, <double*>(y.data), self.flcbdf.yn.ldim_)
        cdef Vec *yPrimeVec = new Vec(REF, <double*>(yPrime.data),self.flcbdf.yn.ldim_)
        DT = self.flcbdf.chooseInitialStepSize(t0,tOut,yVec[0],yPrimeVec[0])
        del yVec
        del yPrimeVec
        return DT

    def setInitialGuess(self, np.ndarray y):
        yn = np.asarray(<double[:self.flcbdf.yn.ldim_]> self.flcbdf.yn.p_)
        y[:] = yn
        
    def lastStepErrorOk(self, np.ndarray y):
        cdef Vec *yVec = new Vec(REF,<double*>(y.data),self.flcbdf.yn.ldim_)
        Ok =  not self.flcbdf.errorForStepTooLarge(yVec[0])
        del yVec
        return int(Ok)

    def calculate_yprime(self,
                          np.ndarray y,
                          np.ndarray Dy,
                          np.ndarray yprime,
                          np.ndarray Dyprime):
        if self.yVec is NULL:
            self.yVec = new Vec(REF, <double*>(y.data), self.flcbdf.yn.ldim_)
            self.DyVec = new Vec(REF, <double*>(Dy.data), self.flcbdf.yn.ldim_)
            self.yprimeVec = new Vec(REF, <double*>(yprime.data), self.flcbdf.yn.ldim_)
            self.DyprimeVec = new Vec(REF, <double*>(Dyprime.data), self.flcbdf.yn.ldim_)
        self.flcbdf.calculate_yprime(self.yVec[0], self.DyVec[0], self.yprimeVec[0], self.DyprimeVec[0])

    def stepTaken(self, np.ndarray y):
        cdef Vec *yVec = new Vec(REF,<double*>(y.data),self.flcbdf.yn.ldim_)
        self.flcbdf.estimateError(yVec[0])
        del yVec
        
    def retryStep_errorFailure(self):
        cdef double h
        h = self.flcbdf.retryStep_errorFailure()
        return h

    def retryStep_solverFailure(self):
        cdef double h
        h = self.flcbdf.retryStep_solverFailure()
        return h

    def initializeTimeHistory(self,
                               np.ndarray y,
                               np.ndarray yPrime):
        cdef Vec *yVec = new Vec(REF,<double*>(y.data),self.flcbdf.yn.ldim_)
        cdef Vec *yPrimeVec = new Vec(REF,<double*>(yPrime.data),self.flcbdf.yn.ldim_)
        self.flcbdf.initializeTimeHistory(yVec[0],yPrimeVec[0])
        del yVec
        del yPrimeVec
        
    def setTolerances(self,
                      double atol,
                      double rtol,
                      np.ndarray dV):
        cdef Vec *dV_Vec  = new Vec(REF,<double*>(dV.data),self.flcbdf.yn.ldim_)
        self.wNormp.setTolerances(atol,rtol,dV_Vec[0])
        del dV_Vec
        
    def getCurrentAlpha(self):
        cdef double alpha
        alpha = self.flcbdf.getCurrentAlpha()
        return alpha

    def __cinit__(self,
                   np.ndarray y,
                   object yName):
        initialized=True
        cdef int ignore1
        cdef char** ignore2
        self.petscSys = new Sys(ignore1, ignore2, NULL, NULL)
        cdef int dim=y.size
        self.sizeVec = new Vec(REF,<double*>(y.data),dim)#this is so we have a parallel vector in the registry
        self.sizeVec.setExample()
        self.wNormp = new WeightedL2Norm(dim)
        if not isinstance(yName, bytes):
            yName = yName.encode()
        dataFilename = "data_{0:d}_{1:s}.txt".format(dim,yName)
        self.data = new FullDataFile(0.0, dataFilename)
        self.flcbdf = new FLCBDF_lite(dim, self.wNormp[0], self.data[0])

    def __dealloc__(self):
        del self.wNormp
        del self.petscSys
        del self.data
        del self.flcbdf
        del self.sizeVec

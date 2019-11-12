# A type of -*- python -*- file
import numpy as np
cimport numpy as np
cimport equivalent_polynomials as eqp
#Note: this Simplex class is for testing equivalent polynomials in python
#It uses a simplistic approach to dealing with the Simplex template
#It is not intended to take full advantage of the C++ implementation

cdef extern from *:
    ctypedef int nSpace1T "1"
    ctypedef int nSpace2T "2"
    ctypedef int nSpace3T "3"
    ctypedef int nP1T "1"
    ctypedef int nP2T "2"
    ctypedef int nP3T "3"
    ctypedef int nQT "50"
#nQ=50 will provide enough space for testing most quadrature rules
#but will be slow


cdef class Simplex:
    cdef np.ndarray xiBuffer
    cdef int nSpace
    cdef int nP
    cdef int nQ
    cdef int q
    cdef bool inside_out
    cdef np.ndarray _H
    cdef np.ndarray _ImH
    cdef np.ndarray _D
    #instatiate template classes for 1,2,3D and P1-P3
    cdef eqp.cSimplex[nSpace1T,nP1T,nQT] s11
    cdef eqp.cSimplex[nSpace1T,nP2T,nQT] s12
    cdef eqp.cSimplex[nSpace1T,nP3T,nQT] s13
    cdef eqp.cSimplex[nSpace2T,nP1T,nQT] s21
    cdef eqp.cSimplex[nSpace2T,nP2T,nQT] s22
    cdef eqp.cSimplex[nSpace2T,nP3T,nQT] s23
    cdef eqp.cSimplex[nSpace3T,nP1T,nQT] s31
    cdef eqp.cSimplex[nSpace3T,nP2T,nQT] s32
    cdef eqp.cSimplex[nSpace3T,nP3T,nQT] s33
    def __cinit__(self, nSpace, nP, nQ):
        self.xiBuffer=np.zeros((50,3),'d')
        self.nSpace = nSpace
        self.nP = nP
        self.nQ=nQ
        self.q=0
    def calculate(self, np.ndarray phi_dof, np.ndarray phi_nodes, np.ndarray xi):
        self.xiBuffer[:xi.shape[0]]=xi
        if (self.nSpace,self.nP) == (1,1):
            icase = self.s11.calculate(<double*>(phi_dof.data), <double*>(phi_nodes.data), <double*>(self.xiBuffer.data))
            self._H = np.asarray(<double[:self.nQ]>self.s11.get_H())
            self._ImH = np.asarray(<double[:self.nQ]>self.s11.get_ImH())
            self._D = np.asarray(<double[:self.nQ]>self.s11.get_D())
            self.inside_out = self.s11.inside_out
        elif (self.nSpace,self.nP) == (1,2):
            icase = self.s12.calculate(<double*>(phi_dof.data), <double*>(phi_nodes.data), <double*>(self.xiBuffer.data))
            self._H = np.asarray(<double[:self.nQ]>self.s12.get_H())
            self._ImH = np.asarray(<double[:self.nQ]>self.s12.get_ImH())
            self._D = np.asarray(<double[:self.nQ]>self.s12.get_D())
            self.inside_out = self.s12.inside_out
        elif (self.nSpace,self.nP) == (1,3):
            icase = self.s13.calculate(<double*>(phi_dof.data), <double*>(phi_nodes.data), <double*>(self.xiBuffer.data))
            self._H = np.asarray(<double[:self.nQ]>self.s13.get_H())
            self._ImH = np.asarray(<double[:self.nQ]>self.s13.get_ImH())
            self._D = np.asarray(<double[:self.nQ]>self.s13.get_D())
            self.inside_out = self.s13.inside_out
        elif (self.nSpace,self.nP) == (2,1):
            icase = self.s21.calculate(<double*>(phi_dof.data), <double*>(phi_nodes.data), <double*>(self.xiBuffer.data))
            self._H = np.asarray(<double[:self.nQ]>self.s21.get_H())
            self._ImH = np.asarray(<double[:self.nQ]>self.s21.get_ImH())
            self._D = np.asarray(<double[:self.nQ]>self.s21.get_D())
            self.inside_out = self.s21.inside_out
        elif (self.nSpace,self.nP) == (2,2):
            icase = self.s22.calculate(<double*>(phi_dof.data), <double*>(phi_nodes.data), <double*>(self.xiBuffer.data))
            self._H = np.asarray(<double[:self.nQ]>self.s22.get_H())
            self._ImH = np.asarray(<double[:self.nQ]>self.s22.get_ImH())
            self._D = np.asarray(<double[:self.nQ]>self.s22.get_D())
            self.inside_out = self.s22.inside_out
        elif (self.nSpace,self.nP) == (2,3):
            icase = self.s23.calculate(<double*>(phi_dof.data), <double*>(phi_nodes.data), <double*>(self.xiBuffer.data))
            self._H = np.asarray(<double[:self.nQ]>self.s23.get_H())
            self._ImH = np.asarray(<double[:self.nQ]>self.s23.get_ImH())
            self._D = np.asarray(<double[:self.nQ]>self.s23.get_D())
            self.inside_out = self.s23.inside_out
        if (self.nSpace,self.nP) == (3,1):
            icase = self.s31.calculate(<double*>(phi_dof.data), <double*>(phi_nodes.data), <double*>(self.xiBuffer.data))
            self._H = np.asarray(<double[:self.nQ]>self.s31.get_H())
            self._ImH = np.asarray(<double[:self.nQ]>self.s31.get_ImH())
            self._D = np.asarray(<double[:self.nQ]>self.s31.get_D())
            self.inside_out = self.s31.inside_out
        elif (self.nSpace,self.nP) == (3,2):
            icase = self.s32.calculate(<double*>(phi_dof.data), <double*>(phi_nodes.data), <double*>(self.xiBuffer.data))
            self._H = np.asarray(<double[:self.nQ]>self.s32.get_H())
            self._ImH = np.asarray(<double[:self.nQ]>self.s32.get_ImH())
            self._D = np.asarray(<double[:self.nQ]>self.s32.get_D())
            self.inside_out = self.s32.inside_out
        elif (self.nSpace,self.nP) == (3,3):
            icase = self.s33.calculate(<double*>(phi_dof.data), <double*>(phi_nodes.data), <double*>(self.xiBuffer.data))
            self._H = np.asarray(<double[:self.nQ]>self.s33.get_H())
            self._ImH = np.asarray(<double[:self.nQ]>self.s33.get_ImH())
            self._D = np.asarray(<double[:self.nQ]>self.s33.get_D())
            self.inside_out = self.s33.inside_out
    def set_quad(self, int q):
        self.q=q
    @property
    def H(self):
        if self.inside_out:
            return self._ImH[self.q]
        else:
            return self._H[self.q]
    @property
    def ImH(self):
        if self.inside_out:
            return self._H[self.q]
        else:
            return self._ImH[self.q]
    @property
    def D(self):
        return self._D[self.q]

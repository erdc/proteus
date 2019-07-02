# A type of -*- python -*- file
import numpy as np
cimport numpy as np
cimport equivalent_polynomials as eqp
#Note: this Simplex class is for testing equivalent polynomials in python
#It uses a simplistic approach to dealing with the Simplex template
#It is not intended to take full advantage of the C++ implementation

cdef extern from *:
    ctypedef int nSpace1 "1"
    ctypedef int nSpace2 "2"
    ctypedef int nSpace3 "3"
    ctypedef int nP1 "1"
    ctypedef int nP2 "2"
    ctypedef int nP3 "3"
    ctypedef int nQ "50"
#nQ=50 will provide enough space for testing most quadrature rules
#but will be slow

#instatiate template classes for 1,2,3D and P1-P3
cdef eqp.cSimplex[nSpace1,nP1,nQ] s11
cdef eqp.cSimplex[nSpace1,nP2,nQ] s12
cdef eqp.cSimplex[nSpace1,nP3,nQ] s13
cdef eqp.cSimplex[nSpace2,nP1,nQ] s21
cdef eqp.cSimplex[nSpace2,nP2,nQ] s22
cdef eqp.cSimplex[nSpace2,nP3,nQ] s23
cdef eqp.cSimplex[nSpace3,nP1,nQ] s31
cdef eqp.cSimplex[nSpace3,nP2,nQ] s32
cdef eqp.cSimplex[nSpace3,nP3,nQ] s33
cdef double eps_dummy=0;
cdef double phi_dummy=0;

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
    def __cinit__(self, nSpace, nP, nQ):
        self.xiBuffer=np.zeros((50,3),'d')
        self.nSpace = nSpace
        self.nP = nP
        self.nQ=nQ
        self.q=0
    def calculate(self, np.ndarray phi_dof, np.ndarray phi_nodes, np.ndarray xi):
        self.xiBuffer[:xi.shape[0]]=xi
        if (self.nSpace,self.nP) == (1,1):
            s11.calculate(<double*>(phi_dof.data), <double*>(phi_nodes.data), <double*>(self.xiBuffer.data))
            self._H = np.asarray(<double[:self.nQ]>s11.get_H(eps_dummy, phi_dummy))
            self._ImH = np.asarray(<double[:self.nQ]>s11.get_ImH(eps_dummy, phi_dummy))
            self._D = np.asarray(<double[:self.nQ]>s11.get_D(eps_dummy, phi_dummy))
            self.inside_out = s11.inside_out
        elif (self.nSpace,self.nP) == (1,2):
            s12.calculate(<double*>(phi_dof.data), <double*>(phi_nodes.data), <double*>(self.xiBuffer.data))
            self._H = np.asarray(<double[:self.nQ]>s12.get_H(eps_dummy, phi_dummy))
            self._ImH = np.asarray(<double[:self.nQ]>s12.get_ImH(eps_dummy, phi_dummy))
            self._D = np.asarray(<double[:self.nQ]>s12.get_D(eps_dummy, phi_dummy))
            self.inside_out = s12.inside_out
        elif (self.nSpace,self.nP) == (1,3):
            s13.calculate(<double*>(phi_dof.data), <double*>(phi_nodes.data), <double*>(self.xiBuffer.data))
            self._H = np.asarray(<double[:self.nQ]>s13.get_H(eps_dummy, phi_dummy))
            self._ImH = np.asarray(<double[:self.nQ]>s13.get_ImH(eps_dummy, phi_dummy))
            self._D = np.asarray(<double[:self.nQ]>s13.get_D(eps_dummy, phi_dummy))
            self.inside_out = s13.inside_out
        elif (self.nSpace,self.nP) == (2,1):
            s21.calculate(<double*>(phi_dof.data), <double*>(phi_nodes.data), <double*>(self.xiBuffer.data))
            self._H = np.asarray(<double[:self.nQ]>s21.get_H(eps_dummy, phi_dummy))
            self._ImH = np.asarray(<double[:self.nQ]>s21.get_ImH(eps_dummy, phi_dummy))
            self._D = np.asarray(<double[:self.nQ]>s21.get_D(eps_dummy, phi_dummy))
            self.inside_out = s21.inside_out
        elif (self.nSpace,self.nP) == (2,2):
            s22.calculate(<double*>(phi_dof.data), <double*>(phi_nodes.data), <double*>(self.xiBuffer.data))
            self._H = np.asarray(<double[:self.nQ]>s22.get_H(eps_dummy, phi_dummy))
            self._ImH = np.asarray(<double[:self.nQ]>s22.get_ImH(eps_dummy, phi_dummy))
            self._D = np.asarray(<double[:self.nQ]>s22.get_D(eps_dummy, phi_dummy))
            self.inside_out = s22.inside_out
        elif (self.nSpace,self.nP) == (2,3):
            s23.calculate(<double*>(phi_dof.data), <double*>(phi_nodes.data), <double*>(self.xiBuffer.data))
            self._H = np.asarray(<double[:self.nQ]>s23.get_H(eps_dummy, phi_dummy))
            self._ImH = np.asarray(<double[:self.nQ]>s23.get_ImH(eps_dummy, phi_dummy))
            self._D = np.asarray(<double[:self.nQ]>s23.get_D(eps_dummy, phi_dummy))
            self.inside_out = s23.inside_out
        if (self.nSpace,self.nP) == (3,1):
            s31.calculate(<double*>(phi_dof.data), <double*>(phi_nodes.data), <double*>(self.xiBuffer.data))
            self._H = np.asarray(<double[:self.nQ]>s31.get_H(eps_dummy, phi_dummy))
            self._ImH = np.asarray(<double[:self.nQ]>s31.get_ImH(eps_dummy, phi_dummy))
            self._D = np.asarray(<double[:self.nQ]>s31.get_D(eps_dummy, phi_dummy))
            self.inside_out = s31.inside_out
        elif (self.nSpace,self.nP) == (3,2):
            s32.calculate(<double*>(phi_dof.data), <double*>(phi_nodes.data), <double*>(self.xiBuffer.data))
            self._H = np.asarray(<double[:self.nQ]>s32.get_H(eps_dummy, phi_dummy))
            self._ImH = np.asarray(<double[:self.nQ]>s32.get_ImH(eps_dummy, phi_dummy))
            self._D = np.asarray(<double[:self.nQ]>s32.get_D(eps_dummy, phi_dummy))
            self.inside_out = s32.inside_out
        elif (self.nSpace,self.nP) == (3,3):
            s33.calculate(<double*>(phi_dof.data), <double*>(phi_nodes.data), <double*>(self.xiBuffer.data))
            self._H = np.asarray(<double[:self.nQ]>s33.get_H(eps_dummy, phi_dummy))
            self._ImH = np.asarray(<double[:self.nQ]>s33.get_ImH(eps_dummy, phi_dummy))
            self._D = np.asarray(<double[:self.nQ]>s33.get_D(eps_dummy, phi_dummy))
            self.inside_out = s33.inside_out
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

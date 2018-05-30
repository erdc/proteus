cimport numpy as np
import numpy as np
from libcpp cimport bool


cdef class cCoefficients:
    cdef public:
        object pyCoefficients
        double C

    def __cinit__(self):
        self.C = 1.

    def attachPyCoefficients(self,
                             object pyCoefficients):
        self.pyCoefficients = pyCoefficients
        print("$$$$$$$$", pyCoefficients)
        print("$$$$$$$$$3", self.pyCoefficients)

    def preStep(self):
        pc = self.pyCoefficients
        print("$$$$$$$$$2", pc, self.pyCoefficients)
        self.cppPreStep(q_uOfX=pc.uOfXTatQuadrature,
                        q_J=pc.model.q['abs(det(J))'],
                        q_weights=pc.model.elementQuadratureWeights[('u', 0)],
                        areas_=pc.areas,
                        q_rci=pc.model.q[('r', 0)],
                        q_fci=pc.model.q[('f', 0)],
                        t=pc.t)

    cdef cppPreStep(self,
                    double[:,:] q_uOfX,
                    double[:,:] q_J,
                    double[:] q_weights,
                    double[:] areas_,
                    double[:, :] q_rci,
                    double[:, :, :] q_fci,
                    double t):
        cdef double integral_1_over_f = 0.
        cdef int N_eN = q_J.shape[0]
        cdef int nE = 0
        cdef double area = 0
        for e in xrange(len(q_J)):
            area = 0
            for k in xrange(len(q_weights)):
                area += q_J[e, k]*q_weights[k]
                integral_1_over_f += q_J[e, k]*q_weights[k]/q_uOfX[e, k]
            areas_[e] = area
            nE += 1
        self.C = integral_1_over_f/nE  # update scaling coefficient

        for eN in range(len(q_rci)):
            for k in range(len(q_rci[eN])):
                for kk in range(len(q_fci[eN, k])):
                    q_fci[eN, k, kk] = 0.0
                q_rci[eN, k] = -(1./(q_uOfX[eN, k]*self.C)-1./areas_[eN])

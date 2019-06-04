# A type of -*- python -*- file
cdef extern from "pskRelations.h":
    cdef cppclass SimplePSK:
        pass
    cdef cppclass PskSpline:
        pass
    cdef cppclass VGM:
        double psic
        double krw
        double krn
        double dpsic
        double dkrw
        double dkrn
        double Se
        double dSe_dpsic
        VGM()
        VGM(const VGM)
        VGM(double* rwork, double* iwork)
        VGM(double* rwork)
        void setTolerances(const double* rwork_tol)
        void calc(const double& Sw)
        void calc_from_psic(const double& psic)

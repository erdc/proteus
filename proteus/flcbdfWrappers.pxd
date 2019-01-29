# A type of -*- python -*- file
cimport numpy as np
from libcpp cimport bool

cdef extern from "DaetkPetscSys.h" namespace "Daetk::Petsc::Sys":
    cdef bool initialized

cdef extern from "DaetkPetscSys.h" namespace "Daetk::Petsc":
    cdef cppclass Sys:
        Sys()
        Sys(int& argc, char **argv, char* help, char* file)

cdef extern from "DaetkPetscVec.h" namespace "Daetk::Petsc::Vec":
    cdef enum CopyType :  NOT_REF, REF

cdef extern from "DaetkPetscVec.h" namespace "Daetk::Petsc":
    cdef cppclass Vec:
        enum CopyType :  NOT_REF, REF
        Vec()
        Vec(CopyType t, double* a, int localDim)
        void setExample()
        int dim()
        int ldim_
        double* p_

cdef extern from "FullDataFile.h" namespace "Daetk":
    cdef cppclass FullDataFile:
        FullDataFile(const double t0, const char* filename)

cdef extern from "WeightedRMSNorm.h" namespace "Daetk":
    cdef cppclass VectorNorm:
        pass
    cdef cppclass WeightedL2Norm(VectorNorm):
        WeightedL2Norm()
        WeightedL2Norm(int dim)
        void setTolerances(const double& atolR,const double& rtolR, const Vec& dV);

cdef extern from "FLCBDF_lite.h" namespace "Daetk":
    cdef cppclass FLCBDF_lite:
        FLCBDF_lite()
        FLCBDF_lite(int dim, WeightedL2Norm&, FullDataFile&)
        void useFixedStep(const double step)
        void useFixedOrder(const int order)
        double chooseDT(const double& , const double& )
        double setDT(const double& hin)
        double chooseInitialStepSize(const double& t, 
                                   const double& tout,
                                   const Vec& y,
                                   const Vec& yPrime)
        bool initializationPhaseStep()
        bool step()
        double estimateError(const Vec& y)
        bool errorForStepTooLarge(const Vec&)
        bool checkError(const Vec& y)
        bool solverFailure()
        bool calculate_yprime(const Vec& y,
                              const Vec& Dy,
                              Vec& yprime,
                              Vec& Dyprime)
        bool initializeTimeHistory(const Vec& yIn, const Vec& yPrimeIn)
        double retryStep_solverFailure()
        double retryStep_errorFailure()
        double getCurrentAlpha() const
        Vec yn,ynprime,yCminusyP,tempvec

cdef class FLCBDF_integrator:
    cdef Sys *petscSys
    cdef WeightedL2Norm *wNormp
    cdef FullDataFile *data
    cdef FLCBDF_lite *flcbdf
    cdef Vec *yVec, *DyVec, *yprimeVec, *DyprimeVec, *sizeVec

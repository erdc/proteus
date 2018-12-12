# A type of -*- python -*- file
cdef extern from "flcbdfWrappersModule.h":

cdef class FLCBDF_integrator:
    cdef Daetk::Petsc::Vec* sizeVec,*yVec,*DyVec,*yprimeVec,*DyprimeVec
    cdef Daetk::Petsc::Sys* petscSys
    cdef Daetk::WeightedL2Norm* wNorm
    cdef Daetk::FullDataFile* data
    cdef Daetk::FLCBDF_lite* flcbdf

cdef class DaetkPetscSys:
    cdef Daetk::Petsc::Sys* petscSys

#get rid of ParVec
#cdef class ParVec:
#    cdef double* array
#    cdef Vec v
#    cdef PyArrayObject* numpy_array

#get rid of ParMat

#get rid of TrueResidualTestCtx
#get rid of CKSP

cdef class FLCBDF_integrator:

cdef class DaetkPetscSys:

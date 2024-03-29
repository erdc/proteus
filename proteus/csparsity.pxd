# A type of -*- python -*- file
# cython: language_level=3
cdef extern from "sparsity.h" namespace "proteus":
    cdef cppclass SparsityInfo:
        int nrows
        int nnz
        int *rowptr
        int *colind
        double *nzval
        SparsityInfo()
        void findNonzeros(int nElements_global,
                          int nDOF_test_element,
                          int nDOF_trial_element,
                          int* nFreeDOF_test,
                          int* freeGlobal_test,
                          int* nFreeDOF_trial,
                          int* freeGlobal_trial,
                          int offset_test,
                          int stride_test,
                          int offset_trial,
                          int stride_trial,
                          int hasNumericalFlux,
                          int hasDiffusionInMixedForm,
                          int needNumericalFluxJacobian,
                          int nElementBoundaries_element,
                          int* elementNeighborsArray,
                          int nInteriorElementBoundaries_global,
                          int* interiorElementBoundariesArray,
                          int* elementBoundaryElementsArray,
                          int* elementBoundaryLocalElementBoundariesArray,
                          int hasFluxBoundaryConditions,
                          int nExteriorElementBoundaries_global,
                          int* exteriorElementBoundariesArray,
                          int hasOutflowBoundary,
                          int needOutflowJacobian);
        void getOffsets_CSR(int nElements_global,
                            int nDOF_test_element,
                            int nDOF_trial_element,
                            int* nFreeDOF_test,
                            int* freeGlobal_test,
                            int* nFreeDOF_trial,
                            int* freeGlobal_trial,
                            int offset_test,
                            int stride_test,
                            int offset_trial,
                            int stride_trial,
                            int hasNumericalFlux,
                            int hasDiffusionInMixedForm,
                            int needNumericalFluxJacobian,
                            int nElementBoundaries_element,
                            int* elementNeighborsArray,
                            int nInteriorElementBoundaries_global,
                            int* interiorElementBoundariesArray,
                            int* elementBoundaryElementsArray,
                            int* elementBoundaryLocalElementBoundariesArray,
                            int hasFluxBoundaryConditions,
                            int nExteriorElementBoundaries_global,
                            int* exteriorElementBoundariesArray,
                            int hasOutflowBoundary,
                            int needOutflowJacobian,
                            int* rowptr,
                            int* csrRowIndeces,
                            int* csrColumnOffsets,
                            int* csrColumnOffsets_eNebN,
                            int* csrColumnOffsets_eb,
                            int* csrColumnOffsets_eb_eNebN)
        void getCSR()

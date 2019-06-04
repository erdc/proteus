#ifndef SPARSITYINFO_H
#define SPARSITYINFO_H
#include <map>
#include <set>
#include <cassert>
#include <algorithm>

namespace proteus
{
  class SparsityInfo
  {
  public:
    SparsityInfo();
    ~SparsityInfo();
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
                        int* csrColumnOffsets_eb_eNebN);
    void getCSR();
    int nrows;
    int nnz;
    int *rowptr, *colind;
    double* nzval;
  private:
    std::map<int,std::set<int> > columnIndecesMap;
    std::map<int,std::map<int,int> > columnOffsetsMap;
  };
}
#endif

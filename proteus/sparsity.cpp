#include "sparsity.h"

namespace proteus
{
  void SparsityInfo::findNonzeros(int nElements_global,
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
                                  int needOutflowJacobian)
  {      
    //elements
    for(int eN=0;eN<nElements_global;eN++)
      {
        for(int ii=0;ii<nFreeDOF_test[eN];ii++)
          {
            int I = offset_test + stride_test*freeGlobal_test[eN*nDOF_test_element+ii];
            for(int jj=0;jj<nFreeDOF_trial[eN];jj++)
              {
                int J = offset_trial + stride_trial*freeGlobal_trial[eN*nDOF_trial_element+jj];
                columnIndecesMap[I].insert(J);
              }
          }
        //get element neighbor DOF for mixed form diffusion
        if (hasNumericalFlux &&
            hasDiffusionInMixedForm &&
            needNumericalFluxJacobian)
          {
            for(int ebN=0;ebN<nElementBoundaries_element;ebN++)
              {
                int eN_ebN = elementNeighborsArray[eN*nElementBoundaries_element + ebN];
                if (eN_ebN >= 0)
                  {
                    for (int ii=0;ii<nFreeDOF_test[eN];ii++)
                      {
                        int I = offset_test + stride_test*freeGlobal_test[eN*nDOF_test_element+ii];
                        for (int jj=0;jj<nFreeDOF_trial[eN_ebN];jj++)
                          {
                            int J = offset_trial + stride_trial*freeGlobal_trial[eN_ebN*nDOF_trial_element+jj];
                            columnIndecesMap[I].insert(J);
                          }
                      }
                  }
              }
          }
      }
    //interior element boundaries
    if (hasNumericalFlux &&
        needNumericalFluxJacobian)
      {
        for(int ebNI=0;ebNI<nInteriorElementBoundaries_global;ebNI++)
          {
            int ebN = interiorElementBoundariesArray[ebNI],
              left_eN_global   = elementBoundaryElementsArray[ebN*2 + 0],
              right_eN_global  = elementBoundaryElementsArray[ebN*2 + 1];
            //left_ebN_element  = elementBoundaryLocalElementBoundariesArray[ebN*2 + 0],
            //right_ebN_element = elementBoundaryLocalElementBoundariesArray[ebN*2 + 1];
            for(int ii=0;ii<nFreeDOF_test[left_eN_global];ii++)
              {
                int left_I = offset_test+stride_test*freeGlobal_test[left_eN_global*nDOF_test_element + ii];
                for (int jj=0;jj<nFreeDOF_trial[left_eN_global];jj++)
                  {
                    int left_J = offset_trial+stride_trial*freeGlobal_trial[left_eN_global*nDOF_trial_element + jj];
                    columnIndecesMap[left_I].insert(left_J);
                  }
                for (int jj=0;jj<nFreeDOF_trial[right_eN_global];jj++)
                  {
                    int right_J = offset_trial+stride_trial*freeGlobal_trial[right_eN_global*nDOF_trial_element + jj];
                    columnIndecesMap[left_I].insert(right_J);
                  }
              }
            for(int ii=0;ii<nFreeDOF_test[right_eN_global];ii++)
              {
                int right_I = offset_test+stride_test*freeGlobal_test[right_eN_global*nDOF_test_element+ii];
                for(int jj=0;jj<nFreeDOF_trial[left_eN_global];jj++)
                  {
                    int left_J = offset_trial+stride_trial*freeGlobal_trial[left_eN_global*nDOF_trial_element+jj];
                    columnIndecesMap[right_I].insert(left_J);
                  }
                for(int jj=0;jj<nFreeDOF_trial[right_eN_global];jj++)
                  {
                    int right_J = offset_trial+stride_trial*freeGlobal_trial[right_eN_global*nDOF_trial_element+jj];
                    columnIndecesMap[right_I].insert(right_J);
                  }
              }
            if(hasDiffusionInMixedForm)
              {
                for(int ebN_eN=0;ebN_eN<nElementBoundaries_element;ebN_eN++)
                  {
                    int left_eN_ebN = elementNeighborsArray[left_eN_global*nElementBoundaries_element+ebN_eN],
                      right_eN_ebN = elementNeighborsArray[right_eN_global*nElementBoundaries_element+ebN_eN];
                    for(int ii=0;ii<nFreeDOF_test[left_eN_global];ii++)
                      {
                        int left_I = offset_test+stride_test*freeGlobal_test[left_eN_global*nDOF_test_element+ii];
                        if(left_eN_ebN >= 0)
                          {
                            for (int jj=0;jj<nFreeDOF_trial[left_eN_ebN];jj++)
                              {
                                int left_J = offset_trial+stride_trial*freeGlobal_trial[left_eN_ebN*nDOF_trial_element+jj];
                                columnIndecesMap[left_I].insert(left_J);
                              }
                          }
                        if(right_eN_ebN >= 0)
                          {
                            for(int jj=0;jj<nFreeDOF_trial[right_eN_ebN];jj++)
                              {
                                int right_J = offset_trial+stride_trial*freeGlobal_trial[right_eN_ebN*nDOF_trial_element+jj];
                                columnIndecesMap[left_I].insert(right_J);
                              }
                          }
                      }
                    for(int ii=0;ii<nFreeDOF_test[right_eN_global];ii++)
                      {
                        int right_I = offset_test+stride_test*freeGlobal_test[right_eN_global*nDOF_test_element+ii];
                        if(left_eN_ebN >= 0)
                          {
                            for(int jj=0;jj<nFreeDOF_trial[left_eN_ebN];jj++)
                              {
                                int left_J = offset_trial+stride_trial*freeGlobal_trial[left_eN_ebN*nDOF_trial_element+jj];
                                columnIndecesMap[right_I].insert(left_J);
                              }
                          }
                        if(right_eN_ebN >= 0)
                          {
                            for(int jj=0;jj<nFreeDOF_trial[right_eN_ebN];jj++)
                              {
                                int right_J = offset_trial+stride_trial*freeGlobal_trial[right_eN_ebN*nDOF_trial_element+jj];
                                columnIndecesMap[right_I].insert(right_J);
                              }
                          }
                      }
                  }
              }
          }
      }
    //exterior element boundaries
    if ((hasNumericalFlux &&
         needNumericalFluxJacobian) ||
        (hasOutflowBoundary &&
         needOutflowJacobian))
      {
        for(int ebNE=0;ebNE<nExteriorElementBoundaries_global;ebNE++)
          {
            int ebN = exteriorElementBoundariesArray[ebNE],
              eN_global   = elementBoundaryElementsArray[ebN*2+0];
            //ebN_element  = elementBoundaryLocalElementBoundariesArray[ebN*2+0];
            for (int ii=0;ii<nFreeDOF_test[eN_global];ii++)
              {
                int I = offset_test+stride_test*freeGlobal_test[eN_global*nDOF_test_element+ii];
                for (int jj=0;jj<nFreeDOF_trial[eN_global];jj++)
                  {
                    int J = offset_trial+stride_trial*freeGlobal_trial[eN_global*nDOF_trial_element+jj];
                    columnIndecesMap[I].insert(J);
                  }
              }
            if(hasNumericalFlux &&
               hasDiffusionInMixedForm)
              {
                for(int ebN_eN=0;ebN_eN < nElementBoundaries_element;ebN_eN++)
                  {
                    int eN_ebN = elementNeighborsArray[eN_global*nElementBoundaries_element+ebN_eN];
                    for (int ii=0;ii<nFreeDOF_test[eN_global];ii++)
                      {
                        int I = offset_test + stride_test*freeGlobal_test[eN_global*nDOF_test_element+ii];
                        if(eN_ebN >= 0)
                          {
                            for(int jj=0;jj<nFreeDOF_trial[eN_ebN];jj++)
                              {
                                int J = offset_trial+stride_trial*freeGlobal_trial[eN_ebN*nDOF_trial_element+jj];
                                columnIndecesMap[I].insert(J);
                              }
                          }
                      }
                  }
              }
          }
      }
    //debug
    //   for (std::map<int, std::set<int> >::iterator mit = columnIndecesMap.begin();mit != columnIndecesMap.end();mit++)
    //     {
    //       for(std::set<int>::iterator sit = mit->second.begin();sit!=mit->second.end();sit++)
    //         std::cout<<*sit<<'\t';
    //       std::cout<<std::endl;
    //     }
  }

  void SparsityInfo::getOffsets_CSR(int nElements_global,
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
  {
    //elements
    for(int eN=0;eN<nElements_global;eN++)
      {
        for(int ii=0;ii<nFreeDOF_test[eN];ii++)
          {
            int I = offset_test + stride_test*freeGlobal_test[eN*nDOF_test_element+ii];
            csrRowIndeces[eN*nDOF_test_element+ii] = rowptr[I];
            for(int jj=0;jj<nFreeDOF_trial[eN];jj++)
              {
                int J = offset_trial + stride_trial*freeGlobal_trial[eN*nDOF_trial_element+jj];
                csrColumnOffsets[eN*nDOF_test_element*nDOF_trial_element+
                                 ii*nDOF_trial_element+
                                 jj] = columnOffsetsMap[I][J];
              }
          }
        //get element neighbor DOF for mixed form diffusion
        if (hasNumericalFlux &&
            hasDiffusionInMixedForm &&
            needNumericalFluxJacobian)
          {
            for(int ebN=0;ebN<nElementBoundaries_element;ebN++)
              {
                int eN_ebN = elementNeighborsArray[eN*nElementBoundaries_element + ebN];
                if (eN_ebN >= 0)
                  {
                    for (int ii=0;ii<nFreeDOF_test[eN];ii++)
                      {
                        int I = offset_test + stride_test*freeGlobal_test[eN*nDOF_test_element+ii];
                        for (int jj=0;jj<nFreeDOF_trial[eN_ebN];jj++)
                          {
                            int J = offset_trial + stride_trial*freeGlobal_trial[eN_ebN*nDOF_trial_element+jj];
                            csrColumnOffsets_eNebN[eN*nElementBoundaries_element*nDOF_test_element*nDOF_trial_element+
                                                   ebN*nDOF_test_element*nDOF_trial_element+
                                                   ii*nDOF_test_element+
                                                   jj] = columnOffsetsMap[I][J];
                          }
                      }
                  }
              }
          }
      }
    //interior element boundaries
    if (hasNumericalFlux &&
        needNumericalFluxJacobian)
      {
        for(int ebNI=0;ebNI<nInteriorElementBoundaries_global;ebNI++)
          {
            int ebN = interiorElementBoundariesArray[ebNI],
              left_eN_global   = elementBoundaryElementsArray[ebN*2 + 0],
              right_eN_global  = elementBoundaryElementsArray[ebN*2 + 1];
            //left_ebN_element  = elementBoundaryLocalElementBoundariesArray[ebN*2 + 0],
            //right_ebN_element = elementBoundaryLocalElementBoundariesArray[ebN*2 + 1];
            for(int ii=0;ii<nFreeDOF_test[left_eN_global];ii++)
              {
                int left_I = offset_test+stride_test*freeGlobal_test[left_eN_global*nDOF_test_element + ii];
                for (int jj=0;jj<nFreeDOF_trial[left_eN_global];jj++)
                  {
                    int left_J = offset_trial+stride_trial*freeGlobal_trial[left_eN_global*nDOF_trial_element + jj];
                    csrColumnOffsets_eb[ebN*2*2*nDOF_test_element*nDOF_trial_element+
                                        0*2*nDOF_test_element*nDOF_trial_element+
                                        0*nDOF_test_element*nDOF_trial_element+
                                        ii*nDOF_trial_element+
                                        jj] = columnOffsetsMap[left_I][left_J];
                  }
                for (int jj=0;jj<nFreeDOF_trial[right_eN_global];jj++)
                  {
                    int right_J = offset_trial+stride_trial*freeGlobal_trial[right_eN_global*nDOF_trial_element + jj];
                    csrColumnOffsets_eb[ebN*2*2*nDOF_test_element*nDOF_trial_element+
                                        0*2*nDOF_test_element*nDOF_trial_element+
                                        1*nDOF_test_element*nDOF_trial_element+
                                        ii*nDOF_trial_element+
                                        jj] = columnOffsetsMap[left_I][right_J];
                  }
              }
            for(int ii=0;ii<nFreeDOF_test[right_eN_global];ii++)
              {
                int right_I = offset_test+stride_test*freeGlobal_test[right_eN_global*nDOF_test_element+ii];
                for(int jj=0;jj<nFreeDOF_trial[left_eN_global];jj++)
                  {
                    int left_J = offset_trial+stride_trial*freeGlobal_trial[left_eN_global*nDOF_trial_element+jj];
                    csrColumnOffsets_eb[ebN*2*2*nDOF_test_element*nDOF_trial_element+
                                        1*2*nDOF_test_element*nDOF_trial_element+
                                        0*nDOF_test_element*nDOF_trial_element+
                                        ii*nDOF_trial_element+
                                        jj] = columnOffsetsMap[right_I][left_J];
                  }
                for(int jj=0;jj<nFreeDOF_trial[right_eN_global];jj++)
                  {
                    int right_J = offset_trial+stride_trial*freeGlobal_trial[right_eN_global*nDOF_trial_element+jj];
                    csrColumnOffsets_eb[ebN*2*2*nDOF_test_element*nDOF_trial_element+
                                        1*2*nDOF_test_element*nDOF_trial_element+
                                        1*nDOF_test_element*nDOF_trial_element+
                                        ii*nDOF_trial_element+
                                        jj] = columnOffsetsMap[right_I][right_J];
                  }
              }
            if(hasDiffusionInMixedForm)
              {
                for(int ebN_eN=0;ebN_eN<nElementBoundaries_element;ebN_eN++)
                  {
                    int left_eN_ebN = elementNeighborsArray[left_eN_global*nElementBoundaries_element+ebN_eN],
                      right_eN_ebN = elementNeighborsArray[right_eN_global*nElementBoundaries_element+ebN_eN];
                    for(int ii=0;ii<nFreeDOF_test[left_eN_global];ii++)
                      {
                        int left_I = offset_test+stride_test*freeGlobal_test[left_eN_global*nDOF_test_element+ii];
                        if(left_eN_ebN >= 0)
                          {
                            for (int jj=0;jj<nFreeDOF_trial[left_eN_ebN];jj++)
                              {
                                int left_J = offset_trial+stride_trial*freeGlobal_trial[left_eN_ebN*nDOF_trial_element+jj];
                                csrColumnOffsets_eb_eNebN[ebN*2*2*nElementBoundaries_element*nDOF_test_element*nDOF_trial_element+
                                                          0*2*nElementBoundaries_element*nDOF_test_element*nDOF_trial_element+
                                                          0*nElementBoundaries_element*nDOF_test_element*nDOF_trial_element+
                                                          ebN_eN*nDOF_test_element*nDOF_trial_element+
                                                          ii*nDOF_trial_element+
                                                          jj] = columnOffsetsMap[left_I][left_J];
                              }
                          }
                        if(right_eN_ebN >= 0)
                          {
                            for(int jj=0;jj<nFreeDOF_trial[right_eN_ebN];jj++)
                              {
                                int right_J = offset_trial+stride_trial*freeGlobal_trial[right_eN_ebN*nDOF_trial_element+jj];
                                csrColumnOffsets_eb_eNebN[ebN*2*2*nElementBoundaries_element*nDOF_test_element*nDOF_trial_element+
                                                          0*2*nElementBoundaries_element*nDOF_test_element*nDOF_trial_element+
                                                          1*nElementBoundaries_element*nDOF_test_element*nDOF_trial_element+
                                                          ebN_eN*nDOF_test_element*nDOF_trial_element+
                                                          ii*nDOF_trial_element+
                                                          jj] = columnOffsetsMap[left_I][right_J];
                              }
                          }
                      }
                    for(int ii=0;ii<nFreeDOF_test[right_eN_global];ii++)
                      {
                        int right_I = offset_test+stride_test*freeGlobal_test[right_eN_global*nDOF_test_element+ii];
                        if(left_eN_ebN >= 0)
                          {
                            for(int jj=0;jj<nFreeDOF_trial[left_eN_ebN];jj++)
                              {
                                int left_J = offset_trial+stride_trial*freeGlobal_trial[left_eN_ebN*nDOF_trial_element+jj];
                                csrColumnOffsets_eb_eNebN[ebN*2*2*nElementBoundaries_element*nDOF_test_element*nDOF_trial_element+
                                                          1*2*nElementBoundaries_element*nDOF_test_element*nDOF_trial_element+
                                                          0*nElementBoundaries_element*nDOF_test_element*nDOF_trial_element+
                                                          ebN_eN*nDOF_test_element*nDOF_trial_element+
                                                          ii*nDOF_trial_element+
                                                          jj] = columnOffsetsMap[right_I][left_J];
                              }
                          }
                        if(right_eN_ebN >= 0)
                          {
                            for(int jj=0;jj<nFreeDOF_trial[right_eN_ebN];jj++)
                              {
                                int right_J = offset_trial+stride_trial*freeGlobal_trial[right_eN_ebN*nDOF_trial_element+jj];
                                csrColumnOffsets_eb_eNebN[ebN*2*2*nElementBoundaries_element*nDOF_test_element*nDOF_trial_element+
                                                          1*2*nElementBoundaries_element*nDOF_test_element*nDOF_trial_element+
                                                          1*nElementBoundaries_element*nDOF_test_element*nDOF_trial_element+
                                                          ebN_eN*nDOF_test_element*nDOF_trial_element+
                                                          ii*nDOF_trial_element+
                                                          jj] = columnOffsetsMap[right_I][right_J];
                              }
                          }
                      }
                  }
              }
          }
      }
    //exterior element boundaries
    if ((hasNumericalFlux &&
         needNumericalFluxJacobian) ||
        (hasOutflowBoundary &&
         needOutflowJacobian))
      {
        for(int ebNE=0;ebNE<nExteriorElementBoundaries_global;ebNE++)
          {
            int ebN = exteriorElementBoundariesArray[ebNE],
              eN_global   = elementBoundaryElementsArray[ebN*2+0];
            //ebN_element  = elementBoundaryLocalElementBoundariesArray[ebN*2+0];
            for (int ii=0;ii<nFreeDOF_test[eN_global];ii++)
              {
                int I = offset_test+stride_test*freeGlobal_test[eN_global*nDOF_test_element+ii];
                for (int jj=0;jj<nFreeDOF_trial[eN_global];jj++)
                  {
                    int J = offset_trial+stride_trial*freeGlobal_trial[eN_global*nDOF_trial_element+jj];
                    csrColumnOffsets_eb[ebN*2*2*nDOF_test_element*nDOF_trial_element+
                                        0*2*nDOF_test_element*nDOF_trial_element+
                                        0*nDOF_test_element*nDOF_trial_element+
                                        ii*nDOF_trial_element+
                                        jj] = columnOffsetsMap[I][J];
                  }
              }
            if(hasNumericalFlux &&
               hasDiffusionInMixedForm)
              {
                for(int ebN_eN=0;ebN_eN < nElementBoundaries_element;ebN_eN++)
                  {
                    int eN_ebN = elementNeighborsArray[eN_global*nElementBoundaries_element+ebN_eN];
                    for (int ii=0;ii<nFreeDOF_test[eN_global];ii++)
                      {
                        int I = offset_test + stride_test*freeGlobal_test[eN_global*nDOF_test_element+ii];
                        if(eN_ebN >= 0)
                          {
                            for(int jj=0;jj<nFreeDOF_trial[eN_ebN];jj++)
                              {
                                int J = offset_trial+stride_trial*freeGlobal_trial[eN_ebN*nDOF_trial_element+jj];
                                csrColumnOffsets_eb_eNebN[ebN*2*2*nElementBoundaries_element*nDOF_test_element*nDOF_trial_element+
                                                          0*2*nElementBoundaries_element*nDOF_test_element*nDOF_trial_element+
                                                          0*nElementBoundaries_element*nDOF_test_element*nDOF_trial_element+
                                                          ebN_eN*nDOF_test_element*nDOF_trial_element+
                                                          ii*nDOF_trial_element+
                                                          jj] = columnOffsetsMap[I][J];
                              }
                          }
                      }
                  }
              }
          }
      }
  }
  
  void SparsityInfo::getCSR()
  {
    //debug
    //   for (std::map<int, std::set<int> >::iterator mit = self->columnIndecesMap.begin();mit != self->columnIndecesMap.end();mit++)
    //     {
    //       for(std::set<int>::iterator sit = mit->second.begin();sit!=mit->second.end();sit++)
    //         std:://cout<<*sit<<'\t';
    //       std::cout<<std::endl;
    //     }
    nrows = int(columnIndecesMap.size())+1;
    rowptr = new int[nrows];
    rowptr[0] = 0;
    for(int I=1;I< columnIndecesMap.size()+1;I++)
      rowptr[I]=rowptr[I-1] + int(columnIndecesMap[I-1].size());
    nnz = rowptr[columnIndecesMap.size()];
    colind = new int[nnz];
    nzval = new double[nnz];
    int max_nonzeros=0;
    for(int I=0;I< columnIndecesMap.size();I++)
      {
        int offset=0;
        for(std::set<int>::iterator sit=columnIndecesMap[I].begin();sit != columnIndecesMap[I].end();sit++)
          {
            columnOffsetsMap[I][*sit] = offset;
            assert(rowptr[I]+offset < nnz);
            colind[rowptr[I]+offset]=*sit;
            offset++;
          }
        std::sort(&colind[rowptr[I]],&colind[rowptr[I+1]]);
        max_nonzeros = std::max(max_nonzeros,rowptr[I+1] - rowptr[I]);
      }
    //std::cout<<"Proteus: Maximum nonzeros in any row is "<<max_nonzeros<<std::endl;
  }
  
  SparsityInfo::SparsityInfo():
    rowptr(NULL),
    colind(NULL),
    nzval(NULL)
  {}
  
  SparsityInfo::~SparsityInfo()
  {
    if (rowptr != NULL)
      delete [] rowptr;
    if (colind != NULL)
      delete [] colind;
    if (nzval != NULL)
      delete [] nzval;
  }
}//proteus

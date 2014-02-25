#include "MeshAdaptPUMI.h"
#include "mMeshIO.h"

using namespace std;
int MeshAdaptPUMIDrvr::TransferSolutionToPUMI(double* inArray, int nVar, int nN) {
   
   int size;
   PUMI_Part_GetNumEnt(PUMI_Part, PUMI_VERTEX, PUMI_ALLTOPO, &size);
   assert(size==nN);
   
   PUMI_Mesh_CreateTag(PUMI_MeshInstance, "Solution", SCUtil_DBL, nVar, SolutionTag);

   pMeshEnt meshEnt;
   pPartEntIter EntIt;

   PUMI_PartEntIter_Init (PUMI_Part, PUMI_VERTEX, PUMI_ALLTOPO, EntIt);
   int isEnd = 0;
   while (!isEnd) {
     PUMI_PartEntIter_GetNext(EntIt, meshEnt);

     int vtxID = PUMI_MeshEnt_ID(meshEnt);
     double* sol  = new double[nVar];
     for(int j=0; j<nVar; j++) {
       sol[j]=inArray[j*nN+vtxID];
     }
     
     PUMI_MeshEnt_SetDblArrTag(PUMI_MeshInstance, meshEnt, SolutionTag, sol, nVar);

     delete [] sol;
     PUMI_PartEntIter_IsEnd(EntIt, &isEnd);
   }
   PUMI_PartEntIter_Del(EntIt);

//   exportMeshToVTK(PUMI_MeshInstance, "pumi");   
}

int MeshAdaptPUMIDrvr::TransferSolutionToProteus(double* outArray, int nVar, int nN) {
   
   pMeshEnt meshEnt;
   pPartEntIter EntIt;

   PUMI_PartEntIter_Init (PUMI_Part, PUMI_VERTEX, PUMI_ALLTOPO, EntIt);
   int isEnd = 0;
   outArray = new double[nVar*nN];
   while (!isEnd) {
     PUMI_PartEntIter_GetNext(EntIt, meshEnt);

     int vtxID = PUMI_MeshEnt_ID(meshEnt);
     
     double *sol = new double [nVar];
     int ncount;
     PUMI_MeshEnt_GetDblArrTag(PUMI_MeshInstance, meshEnt, SolutionTag, &sol, &ncount);
     for(int j=0; j<ncount; j++) {
       outArray[j*nN+vtxID] = sol[j]; 
     }
     PUMI_PartEntIter_IsEnd(EntIt, &isEnd);
   }

}

#include "MeshAdaptPUMI.h"
#include "mMeshIO.h"

using namespace std;
int MeshAdaptPUMIDrvr::TransferSolutionToPUMI(double* inArray, int nVar, int nN) {
   
   int size;
   PUMI_Part_GetNumEnt(PUMI_Part, PUMI_VERTEX, PUMI_ALLTOPO, &size);
   assert(size==nN);
   
   numVar = nVar;
   PUMI_Mesh_CreateTag(PUMI_MeshInstance, "Solution", SCUtil_DBL, nVar, SolutionTag);
   PUMI_Mesh_SetAutoTagMigrOn(PUMI_MeshInstance, SolutionTag, PUMI_ALLTYPE);

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

   printf("Converted Proteus solution to PUMI tags\n");

//   exportMeshToVTK(PUMI_MeshInstance, "pumi");   
}

int MeshAdaptPUMIDrvr::TransferSolutionToProteus(double* outArray, int nVar, int nN) {
   
   pMeshEnt meshEnt;
   pPartEntIter EntIt;

   PUMI_PartEntIter_Init (PUMI_Part, PUMI_VERTEX, PUMI_ALLTOPO, EntIt);
   int isEnd = 0;
//   outArray = new double[nVar*nN];
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
   PUMI_Mesh_DelTag (PUMI_MeshInstance, SolutionTag, 1);

   printf("Converted PUMI tags to Proteus solution\n");
}

int MeshAdaptPUMIDrvr::DeleteMeshEntIDs() {

  int partDim;
  PUMI_Part_GetDim(PUMI_Part, &partDim);
  for (int type=partDim; type>-1; type--) {
    int isEnd = 0;
    int eN = 0;
    pPartEntIter EntIt;
    pMeshEnt meshEnt;
    PUMI_PartEntIter_Init (PUMI_Part, type, PUMI_ALLTOPO, EntIt);
    while (!isEnd) {
      PUMI_PartEntIter_GetNext(EntIt, meshEnt);

      PUMI_MeshEnt_DelID(meshEnt);
      PUMI_PartEntIter_IsEnd(EntIt, &isEnd);
    } //region loop
    PUMI_PartEntIter_Del(EntIt);

  }
}

#include "MeshAdaptPUMI.h"
#include "mMeshIO.h"

int MeshAdaptPUMIDrvr::CalculateSizeField() {

  PUMI_Mesh_CreateTag(PUMI_MeshInstance, "NodeMeshSize", SCUtil_DBL, 1, SFTag);
  PUMI_Mesh_SetAutoTagMigrOn(PUMI_MeshInstance, SFTag, PUMI_ALLTYPE);

  printf("Calculating size field\n");
  
  double max_size=0.06;
  double min_size=0.01; 
  pPartEntIter EntIt;
  pMeshEnt meshEnt;
  PUMI_PartEntIter_Init (PUMI_Part, PUMI_VERTEX, PUMI_ALLTOPO, EntIt);
  int isEnd = 0;
  while (!isEnd) {
     PUMI_PartEntIter_GetNext(EntIt, meshEnt);

//     int owned;
//     PUMI_MeshEnt_IsOwned(meshEnt, PUMI_Part, &owned);
//     if(owned) { 
       double* sol = new double [numVar];
       int ncount;
       PUMI_MeshEnt_GetDblArrTag(PUMI_MeshInstance, meshEnt, SolutionTag, &sol, &ncount);    
       assert(ncount==numVar);

       double* size = new double[1];
       if(sol[4]>0.0 && sol[4]<1.0) {
         size[0] = 0.01;
       }  
       else {
         size[0] = 0.04;
       }
       PUMI_MeshEnt_SetDblArrTag(PUMI_MeshInstance, meshEnt, SFTag, size, 1);

       delete [] size;
       delete [] sol;
//     }
     PUMI_PartEntIter_IsEnd(EntIt, &isEnd);
  }
  PUMI_PartEntIter_Del(EntIt);
  
//  CommuSizeField();

//  PUMI_Mesh_WriteToFile(PUMI_MeshInstance, "Dambreak_debug.smb", 1);  
  exportMeshToVTK(PUMI_MeshInstance, "pumi.vtk"); 
  return 0;
}

int MeshAdaptPUMIDrvr::CommuSizeField() {
  pPartEntIter EntIt;
  pMeshEnt meshEnt;
  PUMI_PartEntIter_InitPartBdry(PUMI_Part, -1, PUMI_VERTEX, PUMI_ALLTOPO, EntIt);
  int isEnd = 0;
  int eN = 0;

  PCU_Comm_Start (PCU_GLOBAL_METHOD);
//  std::pair<pMeshEnt, int>* msg_send = (std::pair<pMeshEnt, int>)malloc(sizeof(std::pair<pMeshEnt, int>));
  std::pair<pMeshEnt, double>* msg_send = new std::pair<pMeshEnt, double>;
  size_t msg_size = sizeof(std::pair<pMeshEnt, double>);

  while (!isEnd) {
    PUMI_PartEntIter_GetNext(EntIt, meshEnt);

    int owned;
    PUMI_MeshEnt_IsOwned(meshEnt, PUMI_Part, &owned);
    if(owned) {
      double* ownSize = new double[1];
      int count;
      PUMI_MeshEnt_GetDblArrTag(PUMI_MeshInstance, meshEnt, SFTag, &ownSize, &count);
      msg_send->second =  ownSize[0];

      std::vector<std::pair<int, pMeshEnt> > VecRemCpy;
      PUMI_MeshEnt_GetAllRmt(meshEnt, VecRemCpy);
      for (int iRC = 0; iRC < VecRemCpy.size(); ++iRC)
      {
        msg_send->first = VecRemCpy[iRC].second;
        PCU_Comm_Write (VecRemCpy[iRC].first, (void*)msg_send, msg_size);
      }
    }
    PUMI_PartEntIter_IsEnd(EntIt, &isEnd);
    eN++;
  }
  PUMI_PartEntIter_Del(EntIt);
  MPI_Barrier(MPI_COMM_WORLD);

  delete msg_send;
  PCU_Comm_Send ();
  size_t recv_size;
  int pid_from;
  void* msg_recv;
  while (PCU_Comm_Read (&pid_from,&msg_recv,&recv_size)) {
     std::pair<pMeshEnt, double>& msg_pair = *(static_cast<std::pair<pMeshEnt, double> *>(msg_recv));
     pMeshEnt copyEnt = msg_pair.first;
     PUMI_MeshEnt_SetDblArrTag(PUMI_MeshInstance, copyEnt, SFTag, &msg_pair.second, 1);
  }
  return 0;
}

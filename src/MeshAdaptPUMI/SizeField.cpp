#include "MeshAdaptPUMI.h"
#include "mMeshIO.h"
#include "apf.h"
#include "apfVector.h"
#include "apfPUMI.h"
#include "apfSPR.h"
#include "apfMesh.h"

struct commData {
  int numNeigh;
  double Size;
};

pTag NumNeighTag;
pTag NeighSFTag;

int MeshAdaptPUMIDrvr::CalculateSizeField(pMAdapt MA_Drvr) {

  PUMI_Mesh_CreateTag(PUMI_MeshInstance, "NodeMeshSize", SCUtil_DBL, 1, SFTag);
  PUMI_Mesh_SetAutoTagMigrOn(PUMI_MeshInstance, SFTag, PUMI_ALLTYPE);

  printf("Calculating size field\n");

  double epsilon=0.03;
  
  pPartEntIter EntIt;
  pMeshEnt meshEnt;
  PUMI_PartEntIter_Init (PUMI_Part, PUMI_VERTEX, PUMI_ALLTOPO, EntIt);
  int isEnd = 0;
  while (!isEnd) {
     PUMI_PartEntIter_GetNext(EntIt, meshEnt);

     double* sol = new double [numVar];
     int ncount;
     PUMI_MeshEnt_GetDblArrTag(PUMI_MeshInstance, meshEnt, SolutionTag, &sol, &ncount);    
     assert(ncount==numVar);
     
     double* size = new double[1];
     double phi = sqrt(sol[5]*sol[5]);
     if(fabs(phi)<epsilon) {
       size[0] = hmin;
     } 
     else if(phi<3*epsilon) {
       size[0]=(hmin+hmax)/2;
     } 
     else {
       size[0] = hmax;
     }
     PUMI_MeshEnt_SetDblArrTag(PUMI_MeshInstance, meshEnt, SFTag, size, 1);
     delete [] size;
     delete [] sol;
     PUMI_PartEntIter_IsEnd(EntIt, &isEnd);
  }
  PUMI_PartEntIter_Del(EntIt);
  
  for(int i=0; i<1; i++)
     SmoothSizeField();

  //set size field to the meshadaptor
  PUMI_PartEntIter_Init (PUMI_Part, PUMI_VERTEX, PUMI_ALLTOPO, EntIt);
  isEnd=0;
  int sizeCounter=0;
  double dVtxSize;
  while (!isEnd)
  {
    PUMI_PartEntIter_GetNext(EntIt, meshEnt);
    if(SCUtil_SUCCESS == PUMI_MeshEnt_GetDblTag (PUMI_MeshInstance, meshEnt, SFTag, &dVtxSize)) {
      MA_SetIsoVtxSize(MA_Drvr, (pVertex)meshEnt, dVtxSize);   // sets size field from tag data
      sizeCounter++;
    }
    PUMI_PartEntIter_IsEnd(EntIt, &isEnd);
  }
  PUMI_PartEntIter_Del(EntIt);
  std::cerr<<" - set size field for "<<sizeCounter<<" vertices\n";
  

//  PUMI_Mesh_WriteToFile(PUMI_MeshInstance, "Dambreak_debug.smb", 1);  
  exportMeshToVTK(PUMI_MeshInstance, "pumi.vtk"); 
  return 0;
}

int MeshAdaptPUMIDrvr::CalculateAnisoSizeField(pMAdapt MA_Drvr, apf::Field* f) {
  
  apf::Mesh* apf_mesh = apf::createMesh(PUMI_Part);
  apf::Field* gradphi = recoverGradientByVolume(f);
  apf::Field* metric = createLagrangeField(apf_mesh,"sizeMetric",apf::MATRIX,1);

  pPartEntIter EntIt;
  pMeshEnt meshEnt;
  int isEnd = 0;
  PUMI_PartEntIter_Init(PUMI_Part,PUMI_VERTEX,PUMI_ALLTOPO,EntIt);
  while (!isEnd) {
    PUMI_PartEntIter_GetNext(EntIt, meshEnt);
    
    double sizeMetric[3][3];
    apf::MeshEntity* e = apf::castEntity(reinterpret_cast<mEntity*>(meshEnt));

    apf::Vector3 gphi;
    double vof;
    apf::Matrix3x3 mtx;
    getVector(gradphi, e, 0, gphi);
    vof = getScalar(voff, e, 0);

    double factor=hmin-hmax;

    double abs_phi = gphi[0]*gphi[0]+gphi[1]*gphi[1]+gphi[2]*gphi[2];
    apf::Vector3 normal_gphi = gphi.normalize();

    for(int i=0;i<3;++i) {
      for(int j=0;j<3;++j) {
        sizeMetric[i][j]=0.0;
        mtx[i][j]=0.0;
      }
    }     
//    sizeMetric[0][0]=hmax; sizeMetric[1][1]=hmax; sizeMetric[2][2]=hmax;
//    mtx[0][0]=hmax; mtx[1][1]=hmax; mtx[2][2]=hmax;
//
    apf::Vector3 y_axis=apf::Vector3(0.0,1.0,0.0); 
    apf::Vector3 x_axis = apf::Vector3(1.0,0.0,0.0);
    apf::Vector3 dir2=cross(normal_gphi,y_axis);
    if(dir2.getLength()<1e-3)
       dir2=cross(normal_gphi,x_axis);
    apf::Vector3 dir3=cross(normal_gphi,dir2);
    
    if(vof>0.0 && vof<1.0) {
//      for(int i=0;i<3;++i) {
        for(int j=0;j<3;++j) {
//          sizeMetric[i][j]=sizeMetric[i][j]+factor*gphi[i]*gphi[j]/abs_phi;
//          mtx[i][j]=mtx[i][j]+factor*gphi[i]*gphi[j]/abs_phi;
          sizeMetric[0][j]=hmin*normal_gphi[j];
          sizeMetric[1][j]=hmax*dir2[j];
          sizeMetric[2][j]=hmax*dir3[j];
        }
//      }
    } else {
      sizeMetric[0][0]=hmax; sizeMetric[1][1]=hmax; sizeMetric[2][2]=hmax;
      mtx[0][0]=hmax; mtx[1][1]=hmax; mtx[2][2]=hmax;
    } 
    MA_SetAnisoVtxSize(MA_Drvr, pVertex(meshEnt), sizeMetric);
//    apf::setMatrix(metric, e, 0, mtx);
 
    PUMI_PartEntIter_IsEnd(EntIt, &isEnd);
  }
  PUMI_PartEntIter_Del(EntIt);
  
  apf:: writeVtkFiles("sizeField", apf_mesh);
  return 0;
} 

int MeshAdaptPUMIDrvr::SmoothSizeField() {

  PUMI_Mesh_CreateTag(PUMI_MeshInstance, "NumberOfNeighbors", SCUtil_INT, 1, NumNeighTag);
  PUMI_Mesh_CreateTag(PUMI_MeshInstance, "NeighborSF", SCUtil_DBL, 1, NeighSFTag);

  pPartEntIter EntIt;
  pMeshEnt meshEnt;
  int isEnd = 0;

  PUMI_PartEntIter_Init(PUMI_Part,PUMI_VERTEX,PUMI_ALLTOPO,EntIt);
  while (!isEnd) {
    PUMI_PartEntIter_GetNext(EntIt, meshEnt);
    int ownedSelf;
    PUMI_MeshEnt_IsOwned(meshEnt, PUMI_Part, &ownedSelf);
    double* Size = new double[1];
    int ncount;
    
    int numNeighVert=0;
    if(ownedSelf) {
      numNeighVert++;
      PUMI_MeshEnt_GetDblArrTag(PUMI_MeshInstance, meshEnt, SFTag, &Size, &ncount);
    } else {
      Size[0]=0.0;
    }

    std::vector<pMeshEnt> vecEdge;
    PUMI_MeshEnt_GetAdj(meshEnt, PUMI_EDGE, 0, vecEdge);
    int iNumEdge = vecEdge.size();
    for (int iEdge = 0; iEdge < iNumEdge; ++iEdge) {
       std::vector<pMeshEnt> vecVtx;
       PUMI_MeshEnt_GetAdj(vecEdge[iEdge], PUMI_VERTEX, 0, vecVtx);
       int iNumVtx = vecVtx.size();
       assert(iNumVtx==2);
       for (int iVtx = 0; iVtx < iNumVtx; ++iVtx) {
          if(vecVtx[iVtx]!=meshEnt) {
            int ownedNeigh;
            PUMI_MeshEnt_IsOwned(vecVtx[iVtx], PUMI_Part, &ownedNeigh);
            if(ownedSelf || ownedNeigh){
              double* SizeNeigh = new double[1];
              PUMI_MeshEnt_GetDblArrTag(PUMI_MeshInstance, vecVtx[iVtx], SFTag, &SizeNeigh, &ncount);
              Size[0]=Size[0]+SizeNeigh[0];
              numNeighVert++;
              delete [] SizeNeigh;
            }
          }
       }
    }
    PUMI_MeshEnt_SetDblArrTag(PUMI_MeshInstance,meshEnt,NeighSFTag,Size,1);
    PUMI_MeshEnt_SetIntTag(PUMI_MeshInstance,meshEnt,NumNeighTag,numNeighVert);
    delete [] Size;
    PUMI_PartEntIter_IsEnd(EntIt, &isEnd);
  }
  PUMI_PartEntIter_Del(EntIt);

  //part boundary vertices
  PUMI_PartEntIter_InitPartBdry(PUMI_Part, -1, PUMI_VERTEX, PUMI_ALLTOPO, EntIt);
  PCU_Comm_Start (PCU_GLOBAL_METHOD);

  std::pair<pMeshEnt, commData>* msg_send = new std::pair<pMeshEnt, commData>;
  size_t msg_size = sizeof(std::pair<pMeshEnt, commData>);
  isEnd=0;
  while (!isEnd) {
    PUMI_PartEntIter_GetNext(EntIt, meshEnt);

    int numNeighVert=0;
    int owned;
    PUMI_MeshEnt_IsOwned(meshEnt, PUMI_Part, &owned);

    double* ownSize = new double[1];
    int count;
    PUMI_MeshEnt_GetDblArrTag(PUMI_MeshInstance, meshEnt, NeighSFTag, &ownSize, &count);
    PUMI_MeshEnt_GetIntTag(PUMI_MeshInstance, meshEnt, NumNeighTag, &numNeighVert);

    msg_send->second.Size =  ownSize[0];
    msg_send->second.numNeigh = numNeighVert;
    std::vector<std::pair<int, pMeshEnt> > VecRemCpy;
    PUMI_MeshEnt_GetAllRmt(meshEnt, VecRemCpy);
    for (int iRC = 0; iRC < VecRemCpy.size(); ++iRC)
    {
       msg_send->first = VecRemCpy[iRC].second;
       PCU_Comm_Write (VecRemCpy[iRC].first, (void*)msg_send, msg_size);
    }
    delete [] ownSize;
    PUMI_PartEntIter_IsEnd(EntIt, &isEnd);
  }
  PUMI_PartEntIter_Del(EntIt);
  MPI_Barrier(MPI_COMM_WORLD);

  delete msg_send;
  PCU_Comm_Send ();

  size_t recv_size;
  int pid_from;
  void* msg_recv;
  while (PCU_Comm_Read (&pid_from,&msg_recv,&recv_size)) {
     std::pair<pMeshEnt, commData>& msg_pair = *(static_cast<std::pair<pMeshEnt, commData> *>(msg_recv));
     pMeshEnt copyEnt = msg_pair.first;

     double* ownSize= new double[1];
     int numNeigh;
     int count;
     PUMI_MeshEnt_GetDblArrTag(PUMI_MeshInstance, copyEnt, NeighSFTag, &ownSize, &count);
     PUMI_MeshEnt_GetIntTag(PUMI_MeshInstance, copyEnt, NumNeighTag, &numNeigh);
//     printf("Size: %lf numNeigh: %d\n",msg_pair.second.Size, msg_pair.second.numNeigh);
     ownSize[0]=ownSize[0]+msg_pair.second.Size;
     numNeigh=numNeigh+msg_pair.second.numNeigh;

//     PUMI_MeshEnt_DelTag(copyEnt, NeighSFTag);
//     PUMI_MeshEnt_DelTag(copyEnt, NumNeighTag);
     PUMI_MeshEnt_SetDblArrTag(PUMI_MeshInstance, copyEnt, NeighSFTag, ownSize, 1);
     PUMI_MeshEnt_SetIntTag(PUMI_MeshInstance, copyEnt, NumNeighTag, numNeigh);
     delete [] ownSize;
  }

  PUMI_PartEntIter_Init(PUMI_Part,PUMI_VERTEX,PUMI_ALLTOPO,EntIt);
  isEnd=0;
  while (!isEnd) {
    PUMI_PartEntIter_GetNext(EntIt, meshEnt);
    double* ownSize= new double[1];
    double* orgSize= new double[1];
    int numNeigh;
    int count;
    PUMI_MeshEnt_GetDblArrTag(PUMI_MeshInstance, meshEnt, NeighSFTag, &ownSize, &count);
    PUMI_MeshEnt_GetIntTag(PUMI_MeshInstance, meshEnt, NumNeighTag, &numNeigh);    
    PUMI_MeshEnt_GetDblArrTag(PUMI_MeshInstance, meshEnt, SFTag, &orgSize, &count);
    
    if(numNeigh==0) printf("Error, no neighbors\n");
    double* avgSize = new double[1];
    avgSize[0]=ownSize[0]/numNeigh;

//    PUMI_MeshEnt_DelTag(meshEnt, SFTag);
    PUMI_MeshEnt_SetDblArrTag(PUMI_MeshInstance, meshEnt, SFTag, avgSize, 1);

    delete [] ownSize;
    delete [] orgSize;
    delete [] avgSize;

    PUMI_PartEntIter_IsEnd(EntIt, &isEnd);
  }
  PUMI_PartEntIter_Del(EntIt);
  PUMI_Mesh_DelTag (PUMI_MeshInstance, NeighSFTag, 1);
  PUMI_Mesh_DelTag (PUMI_MeshInstance, NumNeighTag, 1);

  return 0;
}

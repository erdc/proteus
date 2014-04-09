#include "MeshAdaptPUMI.h"
#include "mMeshIO.h"
#include "apf.h"
#include "apfVector.h"
#include "apfPUMI.h"
#include "apfSPR.h"
#include "apfMesh.h"
#include "Eigen.h"

struct commData {
  int numNeigh;
  double Size[20];
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
/*     
     if(sol[4]>0.0 && sol[4]<1.0) {
        size[0]=hmin;
     } else {
        size[0]=hmax;
     }     
*/     
     PUMI_MeshEnt_SetDblArrTag(PUMI_MeshInstance, meshEnt, SFTag, size, 1);
     delete [] size;
     delete [] sol;
     PUMI_PartEntIter_IsEnd(EntIt, &isEnd);
  }
  PUMI_PartEntIter_Del(EntIt);
  
  for(int i=0; i<3; i++)
     SmoothField(SFTag, 1);

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
   
  PUMI_Mesh_CreateTag(PUMI_MeshInstance, "NodeMeshDir", SCUtil_DBL, 9, SFDirTag);
  PUMI_Mesh_CreateTag(PUMI_MeshInstance, "NodeMeshSize", SCUtil_DBL, 3, SFTag);
  PUMI_Mesh_SetAutoTagMigrOn(PUMI_MeshInstance, SFDirTag, PUMI_ALLTYPE);
  PUMI_Mesh_SetAutoTagMigrOn(PUMI_MeshInstance, SFTag, PUMI_ALLTYPE);

  double epsilon=4*hmin+exp(-nAdapt/4)*hmin*5;
  printf("epsilon value for level set: %lf\n", epsilon);
  apf::Mesh* apf_mesh = apf::createMesh(PUMI_Part);
  apf::Field* gradphi = recoverGradientByVolume(f);
  apf::Field* grad2phi = recoverGradientByVolume(gradphi);
  apf::Field* metric = createLagrangeField(apf_mesh,"sizeMetric",apf::MATRIX,1);
  apf::Field* sizes = createLagrangeField(apf_mesh,"sizes",apf::VECTOR,1);
  apf::Field* gradnorml = createLagrangeField(apf_mesh,"gradnorm",apf::VECTOR,1);
  apf::Field* hess = createLagrangeField(apf_mesh,"hess",apf::MATRIX,1);
  apf::Field* cur = createLagrangeField(apf_mesh,"curve",apf::SCALAR,1);
  apf::Field* direction = createLagrangeField(apf_mesh,"dir",apf::MATRIX,1);

  apf::MeshIterator* it = apf_mesh->begin(0);
  apf::MeshEntity* v;
  double hminsq = hmin*hmin;
  while ((v = apf_mesh->iterate(it)))
  { 
    apf::Vector3 gphi;
    apf::Matrix3x3 g2phi, g2phit;
    getMatrix(grad2phi, v, 0, g2phi);
    g2phit = transpose(g2phi);
    apf::Matrix3x3 gphiprod; 
    gphiprod = tensorProduct(gphi, gphi);
    apf::Matrix3x3 h;
    h=(g2phi+g2phit)/2;
//    h=h/0.001+gphiprod/hminsq;
    setMatrix(hess, v, 0, h);
    getVector(gradphi, v, 0, gphi);
    apf::Vector3 normal_gphi = gphi.normalize();
    setVector(gradnorml, v, 0, normal_gphi);
  }
  apf_mesh->end(it);

  apf::Field* g2norml = recoverGradientByVolume(gradnorml);
  
  Matrix3x3 dir;
  Vector3 size;
  pPartEntIter EntIt;
  pMeshEnt meshEnt;
  int isEnd = 0;
  PUMI_PartEntIter_Init(PUMI_Part,PUMI_VERTEX,PUMI_ALLTOPO,EntIt);
  while (!isEnd) {
    PUMI_PartEntIter_GetNext(EntIt, meshEnt);
    
    double sizeMetric[3][3];
    apf::MeshEntity* e = apf::castEntity(reinterpret_cast<mEntity*>(meshEnt));

    apf::Vector3 gphi, hphi1, hphi2, hphi3;
    double phi,vof;
    apf::Matrix3x3 mtx, he;
    getMatrix(hess, e, 0, he);
    eigen(he, dir, size, 1);
    setMatrix(direction, e, 0, dir);
    setVector(sizes, e, 0, size);
    for(int i=0; i<3; ++i) {
      size[i]=fabs(size[i]);
//      printf("sizes: %lf\n", size[i]);
//      size[i]=sqrt(0.001/size[i]);
    }

//    getMatrix(g2norml, e, 0, he);
    getVector(gradphi, e, 0, gphi);
    vof = getScalar(voff, e, 0);
    phi = getScalar(phif, e, 0);
    apf::Vector3 normal_gphi = gphi.normalize();
    getMatrix(g2norml, e, 0, he);
    double curve = he[0][0]+he[1][1]+he[2][2];

    double dot = dotProd(dir[0], normal_gphi);
    dot=acos(dot/(dir[0].getLength()*normal_gphi.getLength()));
//    printf("angle between gradient and 1st hessian eigenvector %lf\n", dot*180/3.142);

    double* sfdir = new double[9]; int k=0;
    double *sf = new double[3];
    for(int i=0;i<3;++i) {
      for(int j=0;j<3;++j) {
//        sizeMetric[i][j]=dir[i][j]*size[i];
        sizeMetric[i][j]=0.0;
        sfdir[k]=0.0;
        k++;
      }
    }
    setMatrix(metric, e, 0, sizeMetric);
    apf::Vector3 y_axis=apf::Vector3(0.0,1.0,0.0); 
    apf::Vector3 x_axis = apf::Vector3(1.0,0.0,0.0);
    apf::Vector3 z_axis = apf::Vector3(0.0,0.0,1.0);
    apf::Vector3 dir2=cross(normal_gphi,y_axis);
    apf::Vector3 normal_dir2=dir2.normalize();
    if(normal_dir2.getLength()<1e-3) {
       printf("cross with y axis is zero!\n");
       dir2=cross(normal_gphi,x_axis);
    }
    normal_dir2=dir2.normalize();
    apf::Vector3 dir3=cross(normal_gphi,dir2);
    apf::Vector3 normal_dir3=dir3.normalize();

///*   
//    for(int i=0;i<3;++i) {
      for(int j=0;j<3;++j) {
       sfdir[j]=normal_gphi[j];
       sfdir[j+3]=normal_dir2[j];
       sfdir[j+6]=normal_dir3[j];
//         sfdir[i*3+j]=dir[i][j];
      }
//    }
    setVector(sizes, e, 0, size);
///*    
//      printf("sizes: %lf %lf\n", size[1], size[2]);    
//    if(vof<1.0 && vof>0.0) {
    if(sqrt(phi*phi)<3*epsilon) {
      sf[0]=hmin;
      sf[1]=sqrt(0.0004/size[1]);
      sf[2]=sqrt(0.0004/size[2]);
//        sf[0]=hmin; sf[1]=sqrt(0.001/(size[1]+size[2])); sf[2]=hmax/4;
    } else {
      sf[0]=sf[1]=sf[2]=hmax;
//      sizeMetric[0][0]=hmax; sizeMetric[1][1]=hmax; sizeMetric[2][2]=hmax;
    }
//*/   
    for(int i=0;i<3;++i) {
      if(sf[i]<hmin) 
        sf[i]=hmin;
      if(sf[i]>hmax)
        sf[i]=hmax;
    }
    PUMI_MeshEnt_SetDblArrTag(PUMI_MeshInstance, meshEnt, SFDirTag, sfdir, 9);
    PUMI_MeshEnt_SetDblArrTag(PUMI_MeshInstance, meshEnt, SFTag, sf, 3);
    delete [] sf;
    delete [] sfdir;
    PUMI_PartEntIter_IsEnd(EntIt, &isEnd);
  }
  PUMI_PartEntIter_Del(EntIt);

  for(int i=0; i<2; ++i)
     SmoothField(SFTag, 3);

  for(int i=0; i<1; ++i)
     SmoothField(SFDirTag, 9);

  isEnd=0;
  PUMI_PartEntIter_Init(PUMI_Part,PUMI_VERTEX,PUMI_ALLTOPO,EntIt);
  while (!isEnd) {
    PUMI_PartEntIter_GetNext(EntIt, meshEnt);
    double* sf = new double[3]; int ncount;
    double* sfdir = new double[9];
    double sizeMetric[3][3];
    PUMI_MeshEnt_GetDblArrTag(PUMI_MeshInstance, meshEnt, SFTag, &sf, &ncount);    
    PUMI_MeshEnt_GetDblArrTag(PUMI_MeshInstance, meshEnt, SFDirTag, &sfdir, &ncount);    
    int k=0;
    for(int i=0; i<3; ++i) {
      for(int j=0; j<3; ++j) {
        sizeMetric[i][j]=sfdir[k]*sf[i];
        k++;
      }
    }
    MA_SetAnisoVtxSize(MA_Drvr, pVertex(meshEnt), sizeMetric);
    PUMI_PartEntIter_IsEnd(EntIt, &isEnd);
    delete [] sf;
    delete [] sfdir;
  }
  PUMI_PartEntIter_Del(EntIt);
  exportMeshToVTK(PUMI_MeshInstance, "pumi.vtk"); 
  apf:: writeVtkFiles("sizeField", apf_mesh);
  PUMI_Mesh_DelTag (PUMI_MeshInstance, SFDirTag, 1);
  PUMI_Mesh_DelTag (PUMI_MeshInstance, SFTag, 1);
  return 0;
} 

int MeshAdaptPUMIDrvr::SmoothField(pTag tag, int num) {

  PUMI_Mesh_CreateTag(PUMI_MeshInstance, "NumberOfNeighbors", SCUtil_INT, 1, NumNeighTag);
  PUMI_Mesh_CreateTag(PUMI_MeshInstance, "NeighborSF", SCUtil_DBL, num, NeighSFTag);

  pPartEntIter EntIt;
  pMeshEnt meshEnt;
  int isEnd = 0;
  PUMI_PartEntIter_Init(PUMI_Part,PUMI_VERTEX,PUMI_ALLTOPO,EntIt);
  while (!isEnd) {
    PUMI_PartEntIter_GetNext(EntIt, meshEnt);
    int ownedSelf;
    PUMI_MeshEnt_IsOwned(meshEnt, PUMI_Part, &ownedSelf);
    double* Size = new double[num];
    int ncount;
    
    int numNeighVert=0;
    if(ownedSelf) {
      numNeighVert++;
      PUMI_MeshEnt_GetDblArrTag(PUMI_MeshInstance, meshEnt, tag, &Size, &ncount);
    } else {
      for(int i=0; i<num; ++i)
        Size[i]=0.0;
    }

    std::vector<pMeshEnt> vecEdge;
    PUMI_MeshEnt_GetAdj(meshEnt, PUMI_EDGE, 0, vecEdge);
    int iNumEdge = vecEdge.size();
    for (int iEdge = 0; iEdge < iNumEdge; ++iEdge) {
       std::vector<pMeshEnt> vecVtx;
       PUMI_MeshEnt_GetAdj(vecEdge[iEdge], PUMI_VERTEX, 0, vecVtx);
       int iNumVtx = vecVtx.size();
       for (int iVtx = 0; iVtx < iNumVtx; ++iVtx) {
          if(vecVtx[iVtx]!=meshEnt) {
            int ownedNeigh;
            PUMI_MeshEnt_IsOwned(vecVtx[iVtx], PUMI_Part, &ownedNeigh);
            if(ownedSelf || ownedNeigh){
              double* SizeNeigh = new double[num];
              PUMI_MeshEnt_GetDblArrTag(PUMI_MeshInstance, vecVtx[iVtx], tag, &SizeNeigh, &ncount);
              for(int i=0; i<num; ++i)
                 Size[i]=Size[i]+SizeNeigh[i];
              numNeighVert++;
              delete [] SizeNeigh;
            }
          }
       }
    }
    PUMI_MeshEnt_SetDblArrTag(PUMI_MeshInstance,meshEnt,NeighSFTag,Size,num);
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

    double* ownSize = new double[num];
    int count;
    PUMI_MeshEnt_GetDblArrTag(PUMI_MeshInstance, meshEnt, NeighSFTag, &ownSize, &count);
    PUMI_MeshEnt_GetIntTag(PUMI_MeshInstance, meshEnt, NumNeighTag, &numNeighVert);
    
    for(int i=0; i<num; ++i)
      msg_send->second.Size[i] =  ownSize[i];
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

     double* ownSize= new double[num];
     int numNeigh;
     int count;
     PUMI_MeshEnt_GetDblArrTag(PUMI_MeshInstance, copyEnt, NeighSFTag, &ownSize, &count);
     PUMI_MeshEnt_GetIntTag(PUMI_MeshInstance, copyEnt, NumNeighTag, &numNeigh);
//     printf("Size: %lf numNeigh: %d\n",msg_pair.second.Size, msg_pair.second.numNeigh);
     for(int i=0; i<num; ++i)  
       ownSize[i]=ownSize[i]+msg_pair.second.Size[i];
     numNeigh=numNeigh+msg_pair.second.numNeigh;

//     PUMI_MeshEnt_DelTag(copyEnt, NeighSFTag);
//     PUMI_MeshEnt_DelTag(copyEnt, NumNeighTag);
     PUMI_MeshEnt_SetDblArrTag(PUMI_MeshInstance, copyEnt, NeighSFTag, ownSize, num);
     PUMI_MeshEnt_SetIntTag(PUMI_MeshInstance, copyEnt, NumNeighTag, numNeigh);
     delete [] ownSize;
  }

  PUMI_PartEntIter_Init(PUMI_Part,PUMI_VERTEX,PUMI_ALLTOPO,EntIt);
  isEnd=0;
  while (!isEnd) {
    PUMI_PartEntIter_GetNext(EntIt, meshEnt);
    double* ownSize= new double[num];
    double* orgSize= new double[num];
    int numNeigh;
    int count;
    PUMI_MeshEnt_GetDblArrTag(PUMI_MeshInstance, meshEnt, NeighSFTag, &ownSize, &count);
    PUMI_MeshEnt_GetIntTag(PUMI_MeshInstance, meshEnt, NumNeighTag, &numNeigh);    
    PUMI_MeshEnt_GetDblArrTag(PUMI_MeshInstance, meshEnt, tag, &orgSize, &count);
    
    if(numNeigh==0) printf("Error, no neighbors\n");
    double* avgSize = new double[num];
    for(int i=0; i<num; ++i)  
      avgSize[i]=ownSize[i]/numNeigh;

//    PUMI_MeshEnt_DelTag(meshEnt, SFTag);
    PUMI_MeshEnt_SetDblArrTag(PUMI_MeshInstance, meshEnt, tag, avgSize, num);
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

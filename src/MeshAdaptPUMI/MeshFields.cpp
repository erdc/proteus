#include "MeshAdaptPUMI.h"
#include "mMeshIO.h"
#include "apfSPR.h"
#include "apfMesh.h"
#include "apfPUMI.h"

using namespace std;
using namespace apf;

//API To transfer dof arrays from Proteus to PUMI, called from Proteus
//input Array comes from Proteus along with number of variables and the length of the array
//It is a numpy array so the memory is layed out nicely
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
   return 0;
}

//API to transfer solution tags to Proteus dof arrays, called from Proteus
//This is called after adaptation, before starting next solve step to use
//the interpolated solution as the initial conditions
//outArray is dof array in Proteus
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
     delete [] sol;
   }
   PUMI_PartEntIter_Del(EntIt);
   PUMI_Mesh_SetAutoTagMigrOff(PUMI_MeshInstance, SolutionTag, PUMI_VERTEX);
   PUMI_Mesh_DelTag (PUMI_MeshInstance, SolutionTag, 1);
   printf("Converted PUMI tags to Proteus solution\n");
   return 0;
}

//Cleaning up Ent IDs
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
  return 0;
}

//APF support, will play more role in the future, transfer the solution tags to apf fields
//we should be using packed fields here
int MeshAdaptPUMIDrvr::getFieldFromTag(apf::Mesh* mesh, pMeshMdl mesh_pumi, const char* tag_name)
{
  pTag tag;
  int rc = PUMI_Mesh_FindTag(mesh_pumi,tag_name,tag);
  assert(rc == PUMI_SUCCESS);

  presf = apf::createLagrangeField(mesh,"pres",apf::SCALAR,1);
  velf = apf::createLagrangeField(mesh,"vel",apf::VECTOR,1);
  voff = apf::createLagrangeField(mesh,"vof",apf::SCALAR,1);
  phif = apf::createLagrangeField(mesh,"phi",apf::SCALAR,1);
  phidf = apf::createLagrangeField(mesh,"phid",apf::SCALAR,1);
  phiCorrf = apf::createLagrangeField(mesh,"phiCorr",apf::SCALAR,1);

  pPart part;
  PUMI_Mesh_GetPart(mesh_pumi,0,part);
  for (mPart::iterall it = part->beginall(0);
      it != part->endall(0); ++it)
  {
    pMeshEnt vertex = *it;
    double* v = new double[numVar];
    int really_unnecessary;
    rc = PUMI_MeshEnt_GetDblArrTag(mesh_pumi,vertex,tag,
        &v,&really_unnecessary);
    assert(rc == PUMI_SUCCESS);
    setScalar(presf,castEntity(vertex),0,v[0]);
    setScalar(voff,castEntity(vertex),0,v[4]);
    setScalar(phif,castEntity(vertex),0,v[5]);
//    PUMI_MeshEnt_DelTag(vertex,tag);
    delete [] v;
  }
  return 0;
}

//Set solution tags from the apf field, should use packed fields here as well
int MeshAdaptPUMIDrvr::getTagFromField(apf::Mesh* mesh, pMeshMdl mesh_pumi, const char* tag_name)
{
  pTag tag;
  PUMI_Mesh_CreateTag(mesh_pumi,tag_name,PUMI_DBL,1,tag);
  pPart part;
  PUMI_Mesh_GetPart(mesh_pumi,0,part);
  for (mPart::iterall it = part->beginall(0);
      it != part->endall(0); ++it)
  {
    pMeshEnt vertex = *it;
    double pres = getScalar(presf,castEntity(vertex),0);
//    Vector3 vel = getVector(velf,castEntity(vertex),0);
    double vof = getScalar(voff,castEntity(vertex),0);
    double phi = getScalar(phif,castEntity(vertex),0);
//    double phid = getScalar(phidf,castEntity(vertex),0);
//    double phiCorr = getScalar(phiCorrf,castEntity(vertex),0);
    double* v = new double[1];
    v[0]=pres;
/*    
    v[1]=vel[0]; v[2]=vel[1]; v[3]=vel[2];
    v[4]=vof; v[5]=phi; v[6]=phid; v[7]=phiCorr; 
*/    
    PUMI_MeshEnt_SetDblArrTag(mesh_pumi,vertex,tag,v,1);
    delete [] v;
  }
  return 0;
}


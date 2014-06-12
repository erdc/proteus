#include "MeshAdaptPUMI.h"

// implementation of the callback fctn
//extern int CBcounter; /* Do we need this? */
int isEdgeBasedTransfer = 1;
//for measuring procsize

double dist(double xyz1[3],double xyz2[3])
{
  double x1,x2,x3;
  x1=xyz1[0]-xyz2[0];
  x2=xyz1[1]-xyz2[1];
  x3=xyz1[2]-xyz2[2];
  return x1*x1+x2*x2+x3*x3;
}

int inverseMapE (pEdge edge, double* qpt, double* pt )
{
    double xyz[3];
    double node1[3], node2[3];
    double dist1, dist2, distance;

    pVertex nd1, nd2;
    xyz[0] = qpt[0];
    xyz[1] = qpt[1];
    xyz[2] = qpt[2];

    nd1 = E_vertex(edge, 0);
    nd2 = E_vertex(edge, 1);
    V_coord(nd1, node1);
    V_coord(nd2, node2);

    dist1 = sqrt(dist(node1, xyz));
    dist2 = sqrt(dist(node2, xyz));

    distance = dist1+dist2;
   pt[0] = dist2/distance;
   pt[1] = dist1/distance;
   return 1;
}

// arguments are :
// mesh mod. type, look at file MeshSimAdapt.h
// data containing mesh entities created, deleted or reshaped
// userData used to set a pointer to be passed in callback function

void TransferTopSCOREC(pPList oldEnts, pPList newRgn, void *userData, modType mtype, pEntity ent)
{ 
//    std::cerr<<"TransferTopSCOREC() called\n";
    if(isEdgeBasedTransfer && mtype==F_REFINE)
        return;
    if(!isEdgeBasedTransfer && (mtype==E_REFINE || mtype==F_REFINE))
        return;

    pPList fresh;
    pPList old;
    pEntity entity;

    int oldEntsSize = PList_size(oldEnts);
    int freshEntsSize = PList_size(newRgn);

    // old - holds old regions (to be deleted)
    old = PList_new();
    // freash - holds new regions (created)
    fresh = PList_new();

    // VerticestoHandle - hold vertices to be attached solution or deleted
    pPList VerticestoHandle = PList_new();
    PList_append(VerticestoHandle,(void*)ent);

    for(int i=0; i<oldEntsSize; i++){
        entity = (pEntity)PList_item(oldEnts,i);
        PList_append(old,entity);
    }

    for(int i=0; i<freshEntsSize; i++){
        entity = (pEntity)PList_item(newRgn,i);
        PList_append(fresh,entity);
    }
    
 //   if(isEdgeBasedTransfer) //fresh is not used in ESPLIT and E_REFINE
      ((MeshAdaptPUMIDrvr*)userData)->TransferBottomE(old, fresh, VerticestoHandle, mtype);
 //    else
 //      phastaTransferBottom(old, fresh, VerticestoHandle, mtype);    

    PList_delete(fresh);
    PList_delete(old);
    PList_delete(VerticestoHandle);
}

int MeshAdaptPUMIDrvr::TransferBottomE(pPList parent, pPList fresh, pPList VtxstoHandle, modType mtype) 
{
    pEdge edge;
    pVertex vtx;
    double xietazeta[2], xyz[3];

    double field_data_buffer_array[numVar];
    double* field_data = field_data_buffer_array;

    int dont_care;

    if(mtype==E_REFINE || mtype==ESPLIT){
        for(int i=0; i<PList_size(VtxstoHandle);i++) { //loop over the new vertices
            vtx = (pVertex)PList_item(VtxstoHandle, i);

            if(SCUtil_SUCCESS != PUMI_MeshEnt_GetDblArrTag (PUMI_MeshInstance, vtx, SolutionTag, &field_data, &dont_care)) {
                 V_coord(vtx, xyz);
                 edge = (pEdge)PList_item(parent, i);
                 inverseMapE(edge, xyz, xietazeta);
                 InterpolateSolutionE(edge, xietazeta, numVar, SolutionTag, field_data);
                 PUMI_MeshEnt_SetDblArrTag (PUMI_MeshInstance, vtx, SolutionTag, field_data, numVar);
            }

        }
    }

//  take care of the centroid point created by R_REFINE
    else if(mtype==R_REFINE){
        vtx=(pVertex)PList_item(VtxstoHandle, 0);

//      if vertex has field data...          
        if(vtx && (SCUtil_SUCCESS != PUMI_MeshEnt_GetDblArrTag (PUMI_MeshInstance, vtx, SolutionTag, &field_data, &dont_care))) {
            pRegion rgn = (pRegion)PList_item(parent,0);
            pPList vertices = R_vertices(rgn, 1);
            int numVtx = PList_size(vertices);

            double field_buffer_array[numVar];
            double* field = field_buffer_array;
            for(int i=0; i<numVar;i++) { field[i] = 0.0; }
            for(int j=0; j<numVtx; j++) {
                pVertex oldvtx = (pVertex)PList_item(vertices,j);
                if(SCUtil_SUCCESS != PUMI_MeshEnt_GetDblArrTag (PUMI_MeshInstance, oldvtx, SolutionTag, &field_data, &dont_care)) {
                   printf("Error in callback function, no solution attached to an old vertex, velocity\n");
                   exit(-1);
                }
                for(int i=0; i<numVar; i++)
                {
                    field[i] += field_data[i];
                }
            }
            for(int i=0; i<numVar; i++)
            {
                field[i] /= numVtx;
            }
            PUMI_MeshEnt_SetDblArrTag (PUMI_MeshInstance, vtx, SolutionTag, field, numVar);
            PList_delete(vertices);
       }
    }
    return 0;
}
 
int MeshAdaptPUMIDrvr::InterpolateSolutionE( pEdge edge, double xi[2], int field_size, pTag pTagTag, double* result) 
{
     double* vcc1 = new double[field_size];
     double* vcc2 = new double[field_size];
     for(int i=0; i<field_size; i++) {result[i]=0.0;}

     pVertex nd1 = E_vertex(edge,0);
     pVertex nd2 = E_vertex(edge,1);

     int iNumVals;
     if (SCUtil_SUCCESS != PUMI_MeshEnt_GetDblArrTag(PUMI_MeshInstance, nd1, pTagTag, &vcc1, &iNumVals))
     {
       printf("Error in InterpolateSolution: No solution attached to an old vertex \n");
       int isBdry;
       PUMI_MeshEnt_IsPartBdry(nd1, &isBdry);
       printf("is On Bdry: %d\n", isBdry); 
       V_info(nd1);
//       exit(-1);
     }

     if (SCUtil_SUCCESS != PUMI_MeshEnt_GetDblArrTag(PUMI_MeshInstance, nd2, pTagTag, &vcc2, &iNumVals))
     {
       printf("Error in InterpolateSolution: No solution attached to an old vertex \n");
       int isBdry;
       PUMI_MeshEnt_IsPartBdry(nd2, &isBdry);
       printf("is On Bdry: %d\n", isBdry); 
       V_info(nd2);
//       exit(-1);
     }

     //interpolation
     for(int i=0;i<field_size;i++)
     {
       result[i] = vcc1[i]*xi[0]+vcc2[i]*xi[1];
     }

     delete [] vcc1;
     delete [] vcc2;
     return 0;
}



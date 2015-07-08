#include "MeshAdaptPUMI.h"
#include <PCU.h>
#include <petscksp.h>

#include <apf.h>
#include <apfMesh.h>
#include <apfShape.h>
#include <apfDynamicMatrix.h>

#include <iostream>

enum {
  //PHI_IDX = 5
  VOF_IDX = 4
};

//used to attach error estimates to nodes
static void averageToEntity(apf::Field* ef, apf::Field* vf, apf::MeshEntity* ent) //taken from Dan's superconvergent patch recovery code
{
  apf::Mesh* m = apf::getMesh(ef);
  apf::Adjacent elements;
  m->getAdjacent(ent, m->getDimension(), elements);
  double s=0;
  for (std::size_t i=0; i < elements.getSize(); ++i)
    s += apf::getScalar(ef, elements[i], 0);
  s /= elements.getSize();
  apf::setScalar(vf, ent, 0, s);
  return;
}

static apf::Field* extractVOF(apf::Field* solution)
{
  apf::Mesh* m = apf::getMesh(solution);
  apf::Field* voff = apf::createLagrangeField(m,"proteus_vof",apf::SCALAR,1);
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* v;
  apf::NewArray<double> tmp(apf::countComponents(solution));
  while ((v = m->iterate(it))) {
    apf::getComponents(solution, v, 0, &tmp[0]);
    double vof = tmp[VOF_IDX];
    apf::setScalar(voff, v, 0, vof);
  }
  m->end(it);
  return voff;
}

void getViscosity(double &nu_0, double &nu_1)
{
  //***** Viscosity hard-coded hook *****//
  nu_0 = 1.004e-6; //water 
  nu_1 = 1.5e-5; //air
  return;
}


void MeshAdaptPUMIDrvr::get_local_error() 
//This function aims to compute error at each element via ERM.
//First get the mesh and impose a 2nd order field
//Then get the desired quadrature points
{
  getViscosity(nu_0,nu_1);
  //std::cout<<"Viscosity? "<<nu_0<<" "<<nu_1<<std::endl;

  //***** Get Phi First *****//
  apf::Field* voff = extractVOF(solution);
  apf::Field* visc = apf::createLagrangeField(m,"viscosity",apf::SCALAR,1);
  //*****               *****//
  
  //***** Compute the viscosity field *****//
  apf::Mesh*m = apf::getMesh(solution); 
  apf::MeshEntity* ent;
  apf::MeshIterator* iter = m->begin(0);
  double vof_val, visc_val;
  int nsd = m->getDimension();
  while(ent = m->iterate(iter)){ //loop through all elements
    vof_val=apf::getScalar(voff,ent,0);
    visc_val = nu_0*(1-vof_val)+nu_1*vof_val;
    apf::setScalar(visc, ent, 0,visc_val);
  }
  m->end(iter);
  //apf::writeVtkFiles("viscosity", m);

  //Start computing element quantities
  int numqpt; //number of quadrature points
  int nshl,nshl_sol; //number of local shape functions
  int elem_type; //what type of topology
  double weight; //value container for the weight at each qpt
  double Jdet;
  int approx_order = 1; //if using Lagrange shape functions. Serendipity is automatic
  int int_order = approx_order*2; //determines number of integration points
  apf::FieldShape* err_shape = apf::getLagrange(approx_order);
  apf::FieldShape* sol_shape = apf::getShape(solution);
  apf::EntityShape* elem_shape;
  apf::Vector3 qpt; //container for quadrature points
  apf::MeshElement* element;
  apf::Element* visc_elem, *sol_elem;

  apf::Matrix3x3 J; //actual Jacobian matrix
  apf::Matrix3x3 invJ; //inverse of Jacobian
  apf::NewArray <double> shpval, shpval_sol; //array to store shape function values at quadrature points
  apf::NewArray <apf::Vector3> shgval, shgval_sol; //array to store shape function values at quadrature points

  apf::DynamicMatrix invJ_copy;
  apf::NewArray <apf::DynamicVector> shdrv,shdrv_sol;
  apf::NewArray <apf::DynamicVector> shgval_copy, shgval_sol_copy;
  
  iter = m->begin(nsd); //loop over elements
int testcount = 0;
  while(ent = m->iterate(iter)){ //loop through all elements
    element = apf::createMeshElement(m,ent);
    numqpt=apf::countIntPoints(element,int_order); //generally p*p maximum for shape functions
    elem_type = m->getType(ent);
    nshl=apf::countElementNodes(err_shape,elem_type);
    nshl_sol=apf::countElementNodes(sol_shape,elem_type);
    if(testcount==0 && comm_rank==0)
      std::cout<<"nshls "<<nshl<<" "<<nshl_sol<<" numqpt "<<numqpt<<std::endl;
    shpval.allocate(nshl);  shgval.allocate(nshl); shgval_copy.allocate(nshl); shdrv.allocate(nshl);
    shpval_sol.allocate(nshl_sol);  shgval_sol.allocate(nshl_sol); shgval_sol_copy.allocate(nshl_sol); shdrv_sol.allocate(nshl_sol);

    //LHS Matrix Initialization
    int ndofs = nshl*nsd;
    int ndofs_sol = nshl_sol*nsd;
    Mat K; //matrix size depends on nshl, which may vary from element to element
    MatCreate(PETSC_COMM_SELF,&K);
    MatSetType(K,MATSBAIJ); //should be symmetric sparse
    MatSetSizes(K,ndofs,ndofs,ndofs,ndofs);
    MatSeqSBAIJSetPreallocation(K,ndofs,ndofs,NULL); //the local matrices are going to be dense
    
    //loop through all qpts
    for(int k=0;k<numqpt;k++){
      apf::getIntPoint(element,int_order,k,qpt); //get a quadrature point and store in qpt
      apf::getJacobian(element,qpt,J); //evaluate the Jacobian at the quadrature point
      if(J[2][2] == 0){ J[2][2]=1;} //can't compute an inverse with a 0 in a diagonal

      J = apf::transpose(J); //Is PUMI still defined in this way?
      invJ = invert(J);
      Jdet=apf::getJacobianDeterminant(J,nsd); 
      weight = apf::getIntWeight(element,int_order,k);
      invJ_copy = apf::fromMatrix(invJ);

      //first get the shape function values for error shape functions
      elem_shape = err_shape->getEntityShape(elem_type);
      elem_shape->getValues(NULL,NULL,qpt,shpval);
      elem_shape->getLocalGradients(NULL,NULL,qpt,shgval); 

      for(int i =0;i<nshl;i++){ //get the true derivative
        shgval_copy[i] = apf::fromVector(shgval[i]);
        apf::multiply(shgval_copy[i],invJ_copy,shdrv[i]); 
      }

       
      //obtain viscosity value
      visc_elem = apf::createElement(visc,element); //at vof currently
      visc_val = apf::getScalar(visc_elem,qpt);

      /* test visocosity values at qpts
      apf::Element* vof_elem; 
      vof_elem = apf::createElement(voff,element); //at vof currently
      apf::Vector3 xyz;
      apf::mapLocalToGlobal(element,qpt,xyz);

      if(comm_rank==0 && k ==0 ) std::cout<<"VOF "<<apf::getScalar(vof_elem,qpt)<< " Viscosity at qpt " << visc_val<<" at "<<xyz<<std::endl;
      */
      //Left-Hand Side
      PetscScalar term1[nshl][nshl], term2[nshl][nshl];

      //Calculate LHS Diagonal Block Term
      for(int s=0; s<nshl;s++){
        for(int t=0; t<nshl;t++){
          double temp=0;
          for(int j=0;j<nsd;j++){
            temp+=shdrv[s][j]*shdrv[t][j];
          }
          term1[s][t] = temp*weight;
        }
      } 
      int idx[nshl]; //indices for PETSc Mat insertion
      for(int i=0; i<nsd;i++){
        for(int j=0;j<nshl;j++){
          idx[j] = i*nshl+j;
        }
        MatSetValues(K,nshl,idx,nshl,idx,term1[0],ADD_VALUES);
      }

      int idxr[nshl],idxc[nshl]; //indices for PETSc rows and columns
      for(int i = 0; i< nsd;i++){
        for(int j=0; j< nsd;j++){
          for(int s=0;s<nshl;s++){
            for(int t=0;t<nshl;t++){
              term2[s][t] = shdrv[s][j]*shdrv[t][i]*weight;
            }
          }
          for(int count=0;count<nshl;count++){
            idxr[count] = i*nshl+count;
            idxc[count] = j*nshl+count;
          }
          MatSetValues(K,nshl,idxr,nshl,idxc,term2[0],ADD_VALUES);
        }
      } //end 2nd term loop

      //solution shape functions
      elem_shape = sol_shape->getEntityShape(elem_type);
      elem_shape->getValues(NULL,NULL,qpt,shpval_sol);
      elem_shape->getLocalGradients(NULL,NULL,qpt,shgval_sol); 

      for(int i =0;i<nshl_sol;i++){ //get the true derivative
        shgval_sol_copy[i] = apf::fromVector(shgval_sol[i]);
        apf::multiply(shgval_sol_copy[i],invJ_copy,shdrv_sol[i]); 
      }

      //Get RHS
      
      

    } //end loop over quadrature
   
    //to complete integration, scale by the determinant of the Jacobian

testcount++;
    MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);
    MatScale(K,Jdet); //must be done after assembly
    if(comm_rank==0 && testcount==1){ 
      std::cout<<"Jdet ?"<<Jdet<<std::endl;
      MatView(K,PETSC_VIEWER_STDOUT_SELF);
      PetscBool flg;
      MatIsSymmetric(K,0.0,&flg);
      std::cout<<"Is Symmetric "<<flg<<std::endl;
    }
    MatDestroy(&K); //destroy the matrix
  }
  m->end(iter);
  apf::destroyField(voff);
  apf::destroyField(visc);
  printf("It cleared the function.\n");
}



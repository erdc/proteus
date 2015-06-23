#include "MeshAdaptPUMI.h"
#include <PCU.h>
#include <petscksp.h>

#include <apf.h>
#include <apfMesh.h>
#include <apfShape.h>
#include <apfDynamicMatrix.h>

void averageToEntity(apf::Field* ef, apf::Field* vf, apf::MeshEntity* ent);

void MeshAdaptPUMIDrvr::get_local_error() 
//This function aims to compute error at each element via ERM.
//First get the mesh and impose a 2nd order field
//Then get the desired quadrature points
{
  int numqpt; //number of quadrature points
  int nshl; //number of local shape functions
  int elem_type; //what type of topology
  double weight; //value container for the weight at each qpt
  double Jdet;
  int approx_order = 2; //if using Lagrange shape functions. Serendipity is automatic
  apf::Mesh* m = apf::getMesh(solution); 
  apf::FieldShape* err_shape = apf::getLagrange(approx_order);
  apf::EntityShape* elem_shape;
  apf::Vector3 qpt; //container for quadrature points
  apf::MeshElement* element;

  apf::Matrix3x3 J; //actual Jacobian matrix
  apf::Matrix3x3 invJ; //inverse of Jacobian
  apf::NewArray <double> shpval; //array to store shape function values at quadrature points
  apf::NewArray <apf::Vector3> shgval; //array to store shape function values at quadrature points

  apf::DynamicMatrix invJ_copy;
  apf::NewArray <apf::DynamicVector> shdrv;
  apf::NewArray <apf::DynamicVector> shgval_copy;
  

  apf::MeshEntity* ent;
  apf::MeshIterator* iter = m->begin(m->getDimension());
  //Vec F;

  while(ent = m->iterate(iter)){ //loop through all elements
    element = apf::createMeshElement(m,ent);
    numqpt=apf::countIntPoints(element,approx_order);
    elem_type = m->getType(ent);
    nshl=apf::countElementNodes(err_shape,elem_type);

    shpval.allocate(nshl);
    shgval.allocate(nshl);
    shgval_copy.allocate(nshl); //need in Dynamic Vector form
    shdrv.allocate(nshl);
    //loop through all qpts
    for(int k=0;k<numqpt;k++){
      apf::getIntPoint(element,approx_order,k,qpt); //get a quadrature point and store in qpt
      apf::getJacobian(element,qpt,J); //evaluate the Jacobian at the quadrature point
      if(J[2][2] == 0){ J[2][2]=1;} //can't compute an inverse with a 0 in a diagonal

      J = apf::transpose(J); //Is PUMI still defined in this way?
      invJ = invert(J);
      Jdet=apf::getJacobianDeterminant(J,m->getDimension()); 
      weight = apf::getIntWeight(element,approx_order,k);
      elem_shape = err_shape->getEntityShape(elem_type);
      elem_shape->getValues(NULL,NULL,qpt,shpval);
      elem_shape->getLocalGradients(NULL,NULL,qpt,shgval); 
      invJ_copy = apf::fromMatrix(invJ);
/*
      std::cout<<"Jdet "<<Jdet<<std::endl;
      std::cout<<"invJ "<<invJ_copy<<std::endl;
      std::cout<<"What is qpt and weight? "<<qpt<<" "<<weight<<std::endl;
*/
      for(int i =0;i<nshl;i++){ //get the true derivative
        shgval_copy[i] = apf::fromVector(shgval[i]);
        //std::cout<<"What is shgval? "<<shgval_copy[i]<<std::endl;
        apf::multiply(shgval_copy[i],invJ_copy,shdrv[i]); 
        //std::cout<<"What is shdrv? "<<shdrv[i]<<std::endl;
      }

/*
    VecCreate(PETSC_COMM_SELF,&F); 
    VecSetSizes(F,5,5);
    VecSetFromOptions(F); 
    VecSetUp(F);
    for(int i=0;i<5;i++){
      VecSetValue(F,i,2*i+1,INSERT_VALUES);
    }
    VecView(F,PETSC_VIEWER_STDOUT_SELF);
    VecDestroy(&F);
*/
    }
  }
  m->end(iter);
  printf("It cleared the function.\n");
}

void averageToEntity(apf::Field* ef, apf::Field* vf, apf::MeshEntity* ent) //taken from Dan's superconvergent patch recovery code
{
  apf::Mesh* m = apf::getMesh(ef);
  apf::Adjacent elements;
  m->getAdjacent(ent, m->getDimension(), elements);
  double s=0;
  for (std::size_t i=0; i < elements.getSize(); ++i)
    s += apf::getScalar(ef, elements[i], 0);
  s /= elements.getSize();
  apf::setScalar(vf, ent, 0, s);
}


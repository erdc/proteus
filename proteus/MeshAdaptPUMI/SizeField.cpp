#include "MeshAdaptPUMI.h"
#include <apf.h>
#include <apfVector.h>
#include <apfMesh.h>
#include <apfShape.h>
#include <apfDynamicVector.h>
#include <apfCavityOp.h>
#include <string>
#include <iostream>
#include <sstream>
#include <PCU.h>
#include <samElementCount.h>
#include <queue>
#include <algorithm> 
#include <cmath>
#include <Eigen/Eigenvalues>

#define RATIO 8.0

static void SmoothField(apf::Field *f);
void gradeAnisoMesh(apf::Mesh* m,double gradingFactor);
void gradeAspectRatio(apf::Mesh* m, int idx, double gradingFactor);
//double getMetricLength(apf::Matrix3x3 metric, apf::Vector3 directionVector);
//void getGeometricMetricLength(apf::Matrix3x3 metric1, apf::Matrix3x3 metric2, apf::Vector3 directionVector,double &metricLength1, double &metricLength2);
void getGeometricMetricLength(apf::Matrix3x3 metric1, apf::Matrix3x3 metric2, apf::Vector3 directionVector,double &metricLength);

apf::Matrix3x3 getMetricIntersection(apf::Matrix3x3 metric1, apf::Matrix3x3 metric2);

/* Based on the distance from the interface epsilon can be controlled to determine
   thickness of refinement near the interface */
static double isotropicFormula(double phi, double dphi, double verr, double hmin, double hmax, double phi_s = 0, double epsFact = 0)
{
  double size;
  double dphi_size_factor;
  double v_size_factor;
  if (fabs(phi_s) < (epsFact*7.5) * hmin)
    return hmin;
  else
    return hmax;
}

static void setSizeField(apf::Mesh2 *m,apf::MeshEntity *vertex,double h,apf::MeshTag *marker,apf::Field* sizeField)
{
  int isMarked=0;
  if(m->hasTag(vertex,marker))
    isMarked=1;
  double h_new;
  if(isMarked)
    h_new = std::min(h,apf::getScalar(sizeField,vertex,0));
  else
  {
    h_new = h;
    int newMark = 1;
    m->setIntTag(vertex,marker,&newMark);
  }
  apf::setScalar(sizeField,vertex,0,h_new);

  //Parallel Communication with owning copy
  if(!m->isOwned(vertex))
  {
    apf::Copies remotes;
    m->getRemotes(vertex,remotes);
    int owningPart=m->getOwner(vertex);
    PCU_COMM_PACK(owningPart, remotes[owningPart]);
    PCU_COMM_PACK(owningPart, h_new);
  }
}


int MeshAdaptPUMIDrvr::calculateSizeField()
{
  //freeField(size_iso);
  //size_iso = apf::createLagrangeField(m, "proteus_size", apf::SCALAR, 1);
  apf::Field* interfaceBand = apf::createLagrangeField(m, "interfaceBand", apf::SCALAR, 1);
  apf::Field *phif = m->findField("phi");
  assert(phif);

  //Implementation of banded interface, edge intersection algorithm
  apf::MeshTag* vertexMarker = m->createIntTag("vertexMarker",1);
  apf::MeshIterator *it = m->begin(1);
  apf::MeshEntity *edge;

  double safetyFactor = 2.0; //need to make this user defined
  double L_band = N_interface_band*hPhi*safetyFactor;

  PCU_Comm_Begin();
  while ((edge = m->iterate(it)))
  {
    apf::Adjacent edge_adjVerts;
    m->getAdjacent(edge,0,edge_adjVerts);
    apf::MeshEntity *vertex1 = edge_adjVerts[0];
    apf::MeshEntity *vertex2 = edge_adjVerts[1];
    double phi1 = apf::getScalar(phif,vertex1,0);
    double phi2 = apf::getScalar(phif,vertex2,0);
    int caseNumber = 1;
    if(std::fabs(phi1)>L_band)
      caseNumber++;
    if(std::fabs(phi2)>L_band)
      caseNumber++;

    if(caseNumber==1 || caseNumber == 2)
    {
      setSizeField(m,vertex1,hPhi,vertexMarker,interfaceBand);
      setSizeField(m,vertex2,hPhi,vertexMarker,interfaceBand);
    }
    else
    {
      if (phi1*phi2 <0)
      {
        setSizeField(m,vertex1,hPhi,vertexMarker,interfaceBand);
        setSizeField(m,vertex2,hPhi,vertexMarker,interfaceBand);
      }
      else
      {
        setSizeField(m,vertex1,hmax,vertexMarker,interfaceBand);
        setSizeField(m,vertex2,hmax,vertexMarker,interfaceBand);
      }
    }

  }//end while

  PCU_Comm_Send();

  //Take minimum between received value and current value
  apf::MeshEntity *ent;
  double h_received;
  while(PCU_Comm_Receive())
  {
    //Note: the only receiving entities should be owning copies
    PCU_COMM_UNPACK(ent);
    PCU_COMM_UNPACK(h_received);
    //take minimum of received values
    double h_current = apf::getScalar(interfaceBand,ent,0);
    double h_final = std::min(h_current,h_received);
    apf::setScalar(interfaceBand,ent,0,h_final);
  }

  //Synchronization has all remote copies track the owning copy value
  apf::synchronize(interfaceBand);
  m->end(it);

  //Grade the Mesh
  //gradeMesh();
  
  m->destroyTag(vertexMarker);

  //add to queue
  sizeFieldList.push(interfaceBand);
  return 0;
}

void MeshAdaptPUMIDrvr::isotropicIntersect()
{
  freeField(size_iso);
  size_iso = apf::createFieldOn(m, "proteus_size", apf::SCALAR);

  apf::MeshEntity *vert;
  apf::MeshIterator *it = m->begin(0);
  
  apf::Field *field = sizeFieldList.front();
  apf::copyData(size_iso,field);
  sizeFieldList.pop();
  apf::destroyField(field);
  while(!sizeFieldList.empty())
  {
    field = sizeFieldList.front();
    while(vert = m->iterate(it))
    {
      double value1 = apf::getScalar(size_iso,vert,0);
      double value2 = apf::getScalar(field,vert,0);
      double minValue = std::min(value1,value2);
      apf::setScalar(size_iso,vert,0,minValue);
    } 
    sizeFieldList.pop();
    apf::destroyField(field);
  }
  gradeMesh();
}

//taken from Dan's superconvergent patch recovery code
void MeshAdaptPUMIDrvr::averageToEntity(apf::Field *ef, apf::Field *vf,
                                        apf::MeshEntity *ent)
//Function used to convert a region/element field into a vertex field via averaging.
//For each vertex, the average value of the adjacent elements is taken
//Input:
//  ef is the element field
//  ent is the target vertex
//Output:
//  vf is the vertex field
{
  apf::Mesh *m = apf::getMesh(ef);
  apf::Adjacent elements;
  m->getAdjacent(ent, m->getDimension(), elements);
  double s = 0;
  for (std::size_t i = 0; i < elements.getSize(); ++i)
    s += apf::getScalar(ef, elements[i], 0);
  s /= elements.getSize();
  apf::setScalar(vf, ent, 0, s);
  return;
}

void minToEntity(apf::Field* ef, apf::Field* vf,
    apf::MeshEntity* ent)
{
  apf::Mesh* m = apf::getMesh(ef);
  apf::Adjacent elements;
  m->getAdjacent(ent, m->getDimension(), elements);
  double s=0;
  for (std::size_t i=0; i < elements.getSize(); ++i){
    if(i==0)
      s = apf::getScalar(ef, elements[i], 0);
    else if(apf::getScalar(ef,elements[i],0) < s)
      s= apf::getScalar(ef,elements[i],0);
  }
  apf::setScalar(vf, ent, 0, s);
  return;
}

void MeshAdaptPUMIDrvr::volumeAverageToEntity(apf::Field *ef, apf::Field *vf,
                                              apf::MeshEntity *ent)
//Serves the same purpose as averageToEntity but considers a volume-weighted average
{
  apf::Mesh *m = apf::getMesh(ef);
  apf::Adjacent elements;
  apf::MeshElement *testElement;
  m->getAdjacent(ent, m->getDimension(), elements);
  double s = 0;
  double VolumeTotal = 0;
  for (std::size_t i = 0; i < elements.getSize(); ++i)
  {
    testElement = apf::createMeshElement(m, elements[i]);
    s += apf::getScalar(ef, elements[i], 0)*apf::measure(testElement);
    VolumeTotal += apf::measure(testElement);
    if (comm_rank == 0)
    {
      std::cout << "What is s " << s << " Volume? " << apf::measure(testElement) << " scale? " << apf::getScalar(ef, elements[i], 0) << std::endl;
    }
    apf::destroyMeshElement(testElement);
  }
  s /= VolumeTotal;
  if (comm_rank == 0)
  {
    std::cout << "What is s final? " << s << std::endl;
  }
  apf::setScalar(vf, ent, 0, s);
  return;
}

void errorAverageToEntity(apf::Field *ef, apf::Field *vf, apf::Field* err, apf::MeshEntity *ent)
//Serves the same purpose as averageToEntity but considers a error-weighted average
{
  apf::Mesh *m = apf::getMesh(ef);
  apf::Adjacent elements;
  m->getAdjacent(ent, m->getDimension(), elements);
  double s = 0;
  double errorTotal = 0;
  for (std::size_t i = 0; i < elements.getSize(); ++i)
  {
    s += apf::getScalar(ef, elements[i], 0)*apf::getScalar(err,elements[i],0);
    errorTotal += apf::getScalar(err,elements[i],0);
  }
  s /= errorTotal;
/*
  if (comm_rank == 0)
  {
    std::cout << "What is s final? " << s << std::endl;
  }
*/
  apf::setScalar(vf, ent, 0, s);
  return;
}


static apf::Field *extractSpeed(apf::Field *velocity)
//Function used to convert the velocity field into a speed field
{
  apf::Mesh *m = apf::getMesh(velocity);
  apf::Field *speedF = apf::createLagrangeField(m, "proteus_speed", apf::SCALAR, 1);
  apf::MeshIterator *it = m->begin(0);
  apf::MeshEntity *v;
  apf::Vector3 vel_vect;
  while ((v = m->iterate(it)))
  {
    apf::getVector(velocity, v, 0, vel_vect);
    double speed = vel_vect.getLength();
    apf::setScalar(speedF, v, 0, speed);
  }
  m->end(it);
  return speedF;
}

static apf::Matrix3x3 hessianFormula(apf::Matrix3x3 const &g2phi)
//Function used to output a symmetric hessian
{
  apf::Matrix3x3 g2phit = apf::transpose(g2phi);
  return (g2phi + g2phit) / 2;
}

static apf::Field *computeHessianField(apf::Field *grad2phi)
//Function that writes a hessian field
{
  apf::Mesh *m = apf::getMesh(grad2phi);
  apf::Field *hessf = createLagrangeField(m, "proteus_hess", apf::MATRIX, 1);
  apf::MeshIterator *it = m->begin(0);
  apf::MeshEntity *v;
  while ((v = m->iterate(it)))
  {
    apf::Matrix3x3 g2phi;
    apf::getMatrix(grad2phi, v, 0, g2phi);
    apf::Matrix3x3 hess = hessianFormula(g2phi);
    apf::setMatrix(hessf, v, 0, hess);
  }
  m->end(it);
  return hessf;
}

static apf::Matrix3x3 computeGradPhiMetric(apf::Field *gradphi, apf::Field *grad2phi, apf::Field *size_iso, double eps_u,apf::MeshEntity* v,double hPhi)
//Function used to generate a field of metric tensors at vertices. It is meant to be an instantiation or member of computeMetricField .
//Currently needs development.
{
    apf::Mesh *m = apf::getMesh(grad2phi);
    apf::Matrix3x3 g2phi;
    apf::getMatrix(grad2phi, v, 0, g2phi);
    apf::Vector3 gphi;
    apf::getVector(gradphi, v, 0, gphi);
    apf::Matrix3x3 gphigphit(gphi[0] * gphi[0], gphi[0] * gphi[1], gphi[0] * gphi[2],
                             gphi[0] * gphi[1], gphi[1] * gphi[1], gphi[1] * gphi[2],
                             gphi[0] * gphi[2], gphi[1] * gphi[2], gphi[2] * gphi[2]);

    apf::Matrix3x3 hess = hessianFormula(g2phi);
    //apf::Matrix3x3 metric = gphigphit/(apf::getScalar(size_iso,v,0)*apf::getScalar(size_iso,v,0))+ hess/eps_u; 
    //apf::Matrix3x3 metric = gphigphit/(apf::getScalar(size_iso,v,0)*apf::getScalar(size_iso,v,0));
    //apf::Matrix3x3 metric = gphigphit/(hPhi*hPhi);
    //apf::Matrix3x3 metric = hess;
    apf::Matrix3x3 metric = gphigphit/(hPhi*hPhi)+ hess/eps_u;
    
/*
    //std::cout<<"metric is \n"<<metric<<"\n gphi "<< gphi<<" metric length "<<getMetricLength(metric,gphi)<<std::endl; 
    double metricLength =0.0;
    //double metricLength2 =0.0;
    apf::Matrix3x3 identity(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);
    //getGeometricMetricLength(metric, identity, gphi, metricLength1, metricLength2);
    getGeometricMetricLength(identity, identity, gphi, metricLength);
    getMetricIntersection(metric,identity);
    std::abort();
*/
  return metric;
}

static apf::Field *computeMetricField(apf::Field *gradphi, apf::Field *grad2phi, apf::Field *size_iso, double eps_u,double hPhi)
//Function used to generate a field of metric tensors at vertices. It is meant to be an umbrella function that can compute Hessians too.
//Currently needs development.
{
  apf::Mesh *m = apf::getMesh(grad2phi);
  apf::Field *metricf = createLagrangeField(m, "proteus_metric", apf::MATRIX, 1);
  apf::MeshIterator *it = m->begin(0);
  apf::MeshEntity *v;
  apf::Field* phif = m->findField("phi");
  while ((v = m->iterate(it)))
  {
    apf::Matrix3x3 metric;

    //metric = computeGradPhiMetric(gradphi,grad2phi, size_iso, eps_u, v);
    double epsilon = RATIO * hPhi*2.0;
    if( (fabs(apf::getScalar(phif,v,0)) < epsilon ))
    {
      metric = computeGradPhiMetric(gradphi,grad2phi, size_iso, eps_u, v,hPhi);
    }
    else
    {
      for(int i=0; i<3; i++)
      {
        for(int j=0;j<3; j++)
        {
          if(i==j)
            metric[i][j] = 1.0;
          else
            metric[i][j] = 0.0;
        }
      }
    }

    if(m->getDimension()==2){
      metric[0][2]=0.0;
      metric[1][2]=0.0;
      metric[2][2]=1.0;
      metric[2][0]=0.0;
      metric[2][1]=0.0;
    }


/*
    apf::Matrix3x3 g2phi;
    apf::getMatrix(grad2phi, v, 0, g2phi);

    apf::Matrix3x3 hess = hessianFormula(g2phi);
    apf::Matrix3x3 metric = hess;
*/
    apf::setMatrix(metricf, v, 0, metric);
  }
  m->end(it);
  return metricf;
}


// Gaussian, Mean and principal curvatures based on Hessian and gradient of phi.
static void curveFormula(apf::Matrix3x3 const &h, apf::Vector3 const &g,
                         apf::Vector3 &curve)
{
  double a = (h[1][1] + h[2][2]) * g[0] * g[0] + (h[0][0] + h[2][2]) * g[1] * g[1] + (h[0][0] + h[1][1]) * g[2] * g[2];

  double b = g[0] * g[1] * h[0][1] + g[0] * g[2] * h[0][2] + g[1] * g[2] * h[1][2];

  double Km = (a - 2 * b) / pow(g * g, 1.5);

  double c = g[0] * g[0] * (h[1][1] * h[2][2] - h[1][2] * h[1][2]) + g[1] * g[1] * (h[0][0] * h[2][2] - h[0][2] * h[0][2]) + g[2] * g[2] * (h[0][0] * h[1][1] - h[0][1] * h[0][1]);

  double d = g[0] * g[1] * (h[0][2] * h[1][2] - h[0][1] * h[2][2]) + g[1] * g[2] * (h[0][1] * h[0][2] - h[1][2] * h[0][0]) + g[0] * g[2] * (h[0][1] * h[1][2] - h[0][2] * h[1][1]);

  double Kg = (c + 2 * d) / pow(g * g, 2);

  curve[0] = Km;                      //Mean curvature= (k_1+k_2)/2, Not used in normal direction
  curve[1] = Km + sqrt(Km * Km - Kg); //k_1, First principal curvature  (maximum curvature)
  curve[2] = Km - sqrt(Km * Km - Kg); //k_2, Second principal curvature (minimum curvature)
}

static apf::Field *getCurves(apf::Field *hessians, apf::Field *gradphi)
{
  apf::Mesh *m = apf::getMesh(hessians);
  apf::Field *curves;
  curves = apf::createLagrangeField(m, "proteus_curves", apf::VECTOR, 1);
  apf::MeshIterator *it = m->begin(0);
  apf::MeshEntity *v;
  while ((v = m->iterate(it)))
  {
    apf::Matrix3x3 hessian;
    apf::getMatrix(hessians, v, 0, hessian);
    apf::Vector3 gphi;
    apf::getVector(gradphi, v, 0, gphi);
    apf::Vector3 curve;
    curveFormula(hessian, gphi, curve);
    apf::setVector(curves, v, 0, curve);
  }
  m->end(it);
  return curves;
}

static void clamp(double &v, double min, double max)
//Function used to restrict the mesh size to user specified minimums and maximums
{
  v = std::min(max, v);
  v = std::max(min, v);
}

static void clampField(apf::Field *field, double min, double max)
//Function that loops through a field and clamps the values
{
  apf::Mesh *m = apf::getMesh(field);
  int numcomps = apf::countComponents(field);
  double components[numcomps];
  apf::MeshEntity *v;
  apf::MeshIterator *it = m->begin(0);
  while ((v = m->iterate(it)))
  {
    //double tempValue = apf::getScalar(field,v,0);
    apf::getComponents(field, v, 0, &components[0]);
    for (int i = 0; i < numcomps; i++)
      clamp(components[i], min, max);
    //apf::setScalar(field,v,0,tempValue);
    apf::setComponents(field, v, 0, &components[0]);
  }
  m->end(it);
}

static void scaleFormula(double phi, double hmin, double hmax,
                         int adapt_step,
                         apf::Vector3 const &curves,
                         apf::Vector3 &scale)
//Function used to set the size scale vector for the interface size field configuration
{
  double epsilon = 7.0 * hmin;
  if (fabs(phi) < epsilon)
  {
    scale[0] = hmin;
    scale[1] = sqrt(0.002 / fabs(curves[1]));
    scale[2] = sqrt(0.002 / fabs(curves[2]));
  }
  else if (fabs(phi) < 4 * epsilon)
  {
    scale[0] = 2 * hmin;
    scale[1] = 2 * sqrt(0.002 / fabs(curves[1]));
    scale[2] = 2 * sqrt(0.002 / fabs(curves[2]));
  }
  else
  {
    scale = apf::Vector3(1, 1, 1) * hmax;
  }

  //for (int i = 0; i < 3; ++i)
  //  clamp(scale[i], hmin, hmax);
}

static void scaleFormulaERM(double phi, double hmin, double hmax, double h_dest,
                            apf::Vector3 const &curves,
                            double lambda[3], double eps_u, apf::Vector3 &scale,int nsd,double maxAspect)
//Function used to set the size scale vector for the anisotropic ERM size field configuration
//Inputs:
// phi is is the distance to the interface
// hmin is the minimum mesh size
// hmax is the maximum mesh size
// h_dest is the computed new mesh size
// curves is a vector that denotes the curvature of the interface
// lambda is the ordered set of eigenvalues from an eigendecomposition of the metric tensor
// eps_u is a tolerance for the distance away from the interface
//Output:
// scale is the mesh size in each direction for a vertex
{
/*
  double epsilon = 7.0 * hmin;
*/

  //This ensures that the eigenvalues are not too small
  double lambdamin = lambda[0]/maxAspect;
  if (lambda[1] < 1e-10)
  {
    lambda[1] = lambdamin;
    lambda[2] = lambdamin;
  }
  if (lambda[2] < 1e-10)
  {
    lambda[2] = lambdamin;
  }

/*
  if(nsd == 2){
    scale[0] = h_dest;  
    scale[1] = sqrt(lambda[0]/lambda[1])*scale[0];
    scale[2] = 1.0;
  }
  else{
    scale[0] = h_dest;  
    scale[1] = sqrt(lambda[0]/lambda[1])*scale[0];
    scale[2] = sqrt(lambda[0]/lambda[2])*scale[0];
  }
*/

    double epsilon = RATIO * hmin*2.0;
    if (fabs(phi) < epsilon)
    {
      scale[0] = hmin;  
      scale[1] = sqrt(lambda[0]/lambda[1])*scale[0];
      scale[2] = sqrt(lambda[0]/lambda[2])*scale[0];
    }
    else{
      scale[0] = h_dest;
      scale[1] = h_dest;
      scale[2] = h_dest;
    }


  if(scale[1]/scale[0] > maxAspect)
    scale[1] = maxAspect*scale[0];
  if(scale[2]/scale[0] > maxAspect)
    scale[2] = maxAspect*scale[0];
  if(scale[1]/scale[0] > maxAspect || scale[2]/scale[0] > maxAspect){
    std::cout<<"Scales reached maximum aspect ratio\n";
  }

}

static void scaleFormulaPhiVMS(double phi, double hmin, double hmax, double h_dest,
                            apf::Vector3 const &curves,
                            double lambda[3], double eps_u, apf::Vector3 &scale,int nsd,double maxAspect)
//Function used to set the size scale vector for the anisotropic ERM size field configuration
//Inputs:
// phi is is the distance to the interface
// hmin is the minimum mesh size
// hmax is the maximum mesh size
// h_dest is the computed new mesh size
// curves is a vector that denotes the curvature of the interface
// lambda is the ordered set of eigenvalues from an eigendecomposition of the metric tensor
// eps_u is a tolerance for the distance away from the interface
//Output:
// scale is the mesh size in each direction for a vertex
{

  //This ensures that the eigenvalues are not too small
  double lambdamin = lambda[0]/maxAspect;
  if (lambda[1] < 1e-10)
  {
    lambda[1] = lambdamin;
    lambda[2] = lambdamin;
  }
  if (lambda[2] < 1e-10)
  {
    lambda[2] = lambdamin;
  }

/*
  if(nsd == 2){
    scale[0] = hmin;
    scale[1] = sqrt(0.002 / fabs(curves[1]));
    scale[2] = sqrt(0.002 / fabs(curves[2]));
  }
  else{

    double epsilon = 7.0 * hmin;
    if (fabs(phi) < epsilon)
    {
      scale[0] = hmin;
      scale[1] = sqrt(0.002 / fabs(curves[1]));
      scale[2] = sqrt(0.002 / fabs(curves[2]));
    }
    else{
      scale[0] = h_dest;
      scale[1] = h_dest;
      scale[2] = h_dest;
    }
  }
*/
    double epsilon = RATIO * hmin*2.0;
    if (fabs(phi) < epsilon)
    {
      scale[0] = hmin;
      scale[1] = sqrt(0.002 / fabs(curves[1]));
      scale[2] = sqrt(0.002 / fabs(curves[2]));
    }
    else{
      scale[0] = h_dest;
      scale[1] = h_dest;
      scale[2] = h_dest;
    }


  if(scale[1]/scale[0] > maxAspect)
    scale[1] = maxAspect*scale[0];
  if(scale[2]/scale[0] > maxAspect)
    scale[2] = maxAspect*scale[0];
  if(scale[1]/scale[0] > maxAspect || scale[2]/scale[0] > maxAspect){
    std::cout<<"Scales reached maximum aspect ratio\n";
  }

}



static apf::Field *getSizeScales(apf::Field *phif, apf::Field *curves,
                                 double hmin, double hmax, int adapt_step)
//Function that gets size scales in the interface size field configuration
{
  apf::Mesh *m = apf::getMesh(phif);
  apf::Field *scales;
  scales = apf::createLagrangeField(m, "proteus_size_scale", apf::VECTOR, 1);
  apf::MeshIterator *it = m->begin(0);
  apf::MeshEntity *v;
  while ((v = m->iterate(it)))
  {
    double phi = apf::getScalar(phif, v, 0);
    apf::Vector3 curve;
    apf::getVector(curves, v, 0, curve);
    apf::Vector3 scale;

    scaleFormula(phi, hmin, hmax, adapt_step, curve, scale);
    apf::setVector(scales, v, 0, scale);
  }
  m->end(it);
  return scales;
}

struct SortingStruct
{
  apf::Vector3 v;
  double wm;
  bool operator<(const SortingStruct &other) const
  {
    return wm < other.wm;
  }
};

static apf::Field *getSizeFrames(apf::Field *hessians, apf::Field *gradphi)
//Function that gets size frames, or the metric tensor, in the interface size field configuration
{
  apf::Mesh *m = apf::getMesh(gradphi);
  apf::Field *frames;
  frames = apf::createLagrangeField(m, "proteus_size_frame", apf::MATRIX, 1);
  apf::MeshIterator *it = m->begin(0);
  apf::MeshEntity *v;
  while ((v = m->iterate(it)))
  {
    apf::Vector3 gphi;
    apf::getVector(gradphi, v, 0, gphi);
    apf::Vector3 dir;
    if (gphi.getLength() > 1e-16)
      dir = gphi.normalize();
    else
      dir = apf::Vector3(1, 0, 0);
    apf::Matrix3x3 hessian;
    apf::getMatrix(hessians, v, 0, hessian);
    apf::Vector3 eigenVectors[3];
    double eigenValues[3];
    apf::eigen(hessian, eigenVectors, eigenValues);
    SortingStruct ssa[3];
    for (int i = 0; i < 3; ++i)
    {
      ssa[i].v = eigenVectors[i];
      ssa[i].wm = std::fabs(eigenValues[i]);
    }
    std::sort(ssa, ssa + 3);
    assert(ssa[2].wm >= ssa[1].wm);
    assert(ssa[1].wm >= ssa[0].wm);
    double firstEigenvalue = ssa[2].wm;
    apf::Matrix3x3 frame;
    frame[0] = dir;
    if (firstEigenvalue > 1e-16)
    {
      apf::Vector3 firstEigenvector = ssa[2].v;
      frame[1] = apf::reject(firstEigenvector, dir);
      frame[2] = apf::cross(frame[0], frame[1]);
      if (frame[2].getLength() < 1e-16)
        frame = apf::getFrame(dir);
    }
    else
      frame = apf::getFrame(dir);
    for (int i = 0; i < 3; ++i)
      frame[i] = frame[i].normalize();
    frame = apf::transpose(frame);
    apf::setMatrix(frames, v, 0, frame);
  }
  m->end(it);
  return frames;
}

static apf::Field *getERMSizeFrames(apf::Field *hessians, apf::Field *gradphi, apf::Field *frame_comps[3])
//Function that gets size frames, or the metric tensor, in the interface size field configuration
{
  apf::Mesh *m = apf::getMesh(gradphi);
  apf::Field *frames;
  frames = m->findField("proteus_size_frame");
  apf::MeshIterator *it = m->begin(0);
  apf::MeshEntity *v;
  while ((v = m->iterate(it)))
  {
    apf::Matrix3x3 frame(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    apf::Vector3 gphi;
    apf::getVector(gradphi, v, 0, gphi);
/*
    apf::Vector3 dir;
    if (gphi.getLength() > 1e-16)
      dir = gphi.normalize();
    else
      dir = apf::Vector3(1, 0, 0);
*/

    //get eigen values and eigenvectors from hessian
    apf::Matrix3x3 hessian;
    apf::getMatrix(hessians, v, 0, hessian);
    apf::Vector3 eigenVectors[3];
    double eigenValues[3];
    apf::eigen(hessian, eigenVectors, eigenValues);
    SortingStruct ssa[3];
    for (int i = 0; i < 3; ++i)
    {
      ssa[i].v = eigenVectors[i];
      ssa[i].wm = std::fabs(eigenValues[i]);
    }

    //sort eigenvalues and eigenvectors
    std::sort(ssa, ssa + 3);
    assert(ssa[2].wm >= ssa[1].wm);
    assert(ssa[1].wm >= ssa[0].wm);
    double firstEigenvalue = ssa[2].wm;
    assert(firstEigenvalue > 1e-12);
    //frame[0] = dir;
    //if (firstEigenvalue > 1e-16)
    //{
      frame[0] = ssa[2].v;
      frame[1] = ssa[1].v;
      frame[2] = ssa[0].v;
    //}
    //else
    //  frame = apf::getFrame(dir);
    
/*
    apf::Vector3 test(1.0,0.0,0.0);
    frame = apf::getFrame(test);
    apf::setMatrix(frames,v,0,frame);
*/

    //normalize eigenvectors
    for (int i = 0; i < 3; ++i)
      frame[i] = frame[i].normalize();
    frame = apf::transpose(frame);
    apf::setMatrix(frames, v, 0, frame);
    apf::setVector(frame_comps[0], v, 0, frame[0]);
    apf::setVector(frame_comps[1], v, 0, frame[1]);
    apf::setVector(frame_comps[2], v, 0, frame[2]);
  }
  m->end(it);
  return frames;
}

int MeshAdaptPUMIDrvr::calculateAnisoSizeField()
//High level function that obtains the size scales and the size frames for anistropic interface-based adapt
{
  apf::Field *phif = m->findField("phi");
  assert(phif);
  apf::Field *gradphi = apf::recoverGradientByVolume(phif);
  apf::Field *grad2phi = apf::recoverGradientByVolume(gradphi);
  apf::Field *hess = computeHessianField(grad2phi);
  apf::destroyField(grad2phi);
  apf::Field *curves = getCurves(hess, gradphi);
  freeField(size_scale);

  size_scale = getSizeScales(phif, curves, hmin, hmax, nAdapt);
  apf::destroyField(curves);
  freeField(size_frame);
  size_frame = getSizeFrames(hess, gradphi);
  apf::destroyField(hess);

  apf::destroyField(gradphi);
/*
  for (int i = 0; i < 2; ++i)
    SmoothField(size_scale);
*/

  return 0;
}

struct Smoother : public apf::CavityOp
{
  Smoother(apf::Field *f) : apf::CavityOp(apf::getMesh(f))
  {
    field = f;
    int nc = apf::countComponents(f);
    newField = apf::createPackedField(mesh, "proteus_smooth_new", nc);
    sum.setSize(nc);
    value.setSize(nc);
    nApplied = 0;
  }
  ~Smoother()
  {
    copyData(field, newField);
    apf::destroyField(newField);
  }
  virtual Outcome setEntity(apf::MeshEntity *e)
  {
    if (apf::hasEntity(newField, e))
      return SKIP;
    if (!this->requestLocality(&e, 1))
      return REQUEST;
    vertex = e;
    return OK;
  }
  virtual void apply()
  {
    /* the smoothing function used here is the average of the
   vertex value and neighboring vertex values, with the
   center vertex weighted equally as the neighbors */
    apf::Up edges;
    mesh->getUp(vertex, edges);
    apf::getComponents(field, vertex, 0, &sum[0]);
    for (int i = 0; i < edges.n; ++i)
    {
      apf::MeshEntity *ov = apf::getEdgeVertOppositeVert(
          mesh, edges.e[i], vertex);
      apf::getComponents(field, ov, 0, &value[0]);
      sum += value;
    }
    sum /= edges.n + 1;
    apf::setComponents(newField, vertex, 0, &sum[0]);
    ++nApplied;
  }
  apf::Field *field;
  apf::Field *newField;
  apf::MeshEntity *vertex;
  apf::DynamicVector sum;
  apf::DynamicVector value;
  apf::MeshTag *emptyTag;
  int nApplied;
};

static void SmoothField(apf::Field *f)
{
  Smoother op(f);
  op.applyToDimension(0);
}

void getTargetError(apf::Mesh* m, apf::Field* errField, double &target_error,double totalError){
  //Implemented for 3D and for serial case only so far
  //Need to communicate target error in parallel
  if(PCU_Comm_Self()>0) 
    std::cout<<"WARNING/ERROR:Parallel implementation is not completed yet\n";
  if(m->getDimension()==2){
    target_error = totalError/sqrt(m->count(m->getDimension()));
    if(PCU_Comm_Self()==0)
      std::cout<<"The estimated target error is "<<target_error<<std::endl;
    return;
  }
  apf::Field* interfaceField = m->findField("vof");
  apf::Field* targetField = apf::createField(m,"targetError",apf::SCALAR,apf::getVoronoiShape(m->getDimension(),1));
  apf::MeshEntity* ent;
  apf::MeshIterator* it = m->begin(m->getDimension());
  apf::MeshElement* element;
  apf::Element* vofElem;
  std::vector <double> errVect;
  while( (ent = m->iterate(it))){
    element = apf::createMeshElement(m, ent);
    vofElem = apf::createElement(interfaceField,element);
    double vofVal = apf::getScalar(vofElem,apf::Vector3(1./3.,1./3.,1./3.));

    if(vofVal < 0.9 && vofVal > 0.1){ //at the interface
      double errorValue = apf::getScalar(errField,ent,0);
      errVect.push_back(errorValue);    
      apf::setScalar(targetField,ent,0,errorValue);
    }
    else{
      apf::setScalar(targetField,ent,0,0.0);
    }
  }
  m->end(it);
  if(PCU_Comm_Self()==0)
    std::cout<<"Past creation of vector\n";
  if(errVect.size()==0){
    target_error = totalError/sqrt(m->count(m->getDimension()));
  }
  else{
    std::ofstream myfile;
    myfile.open("interfaceErrors.txt", std::ios::app );
    for(int i=0;i<errVect.size();i++){
      myfile << errVect[i]<<std::endl;
    }
    myfile.close();
    std::sort(errVect.begin(),errVect.end());
    int vectorSize = errVect.size();
    if(vectorSize %2 ==0){
      int idx1 = vectorSize/2-1;
      target_error = (errVect[idx1]+errVect[idx1+1])/2; //get average
    }
    else
      target_error = errVect[(vectorSize-1)/2];
  }
  if(PCU_Comm_Self()==0)
    std::cout<<"The estimated target error is "<<target_error<<std::endl;
  //std::abort();
}

int MeshAdaptPUMIDrvr::getERMSizeField(double err_total)
//High level function that obtains the size scales and the size frames for ERM-based adapt and uses the computed total error
{

  freeField(size_frame);
  freeField(size_scale);
  freeField(size_iso);

  //Initialize fields and needed types/variables
  apf::Field* errField;
  //apf::Mesh* m;
  if(size_field_config=="ERM")
    errField = m->findField("ErrorRegion");
  else if(size_field_config=="VMS" || size_field_config=="combined")
    errField = m->findField("VMSH1");
  assert(errField); 
  //apf::Mesh *m = apf::getMesh(vmsErrH1);
  //apf::getMesh(errField);
  apf::MeshIterator *it;
  apf::MeshEntity *v;
  apf::MeshElement *element;
  apf::MeshEntity *reg;
  //size_iso = apf::createLagrangeField(m, "proteus_size", apf::SCALAR, 1);
  apf::Field *errorSize = apf::createLagrangeField(m, "errorSize", apf::SCALAR, 1);

  if (adapt_type_config == "anisotropic"){
    size_scale = apf::createLagrangeField(m, "proteus_size_scale", apf::VECTOR, 1);
    size_frame = apf::createLagrangeField(m, "proteus_size_frame", apf::MATRIX, 1);
  }
  apf::Field *errorSize_reg = apf::createField(m, "iso_size", apf::SCALAR, apf::getConstant(nsd));
  apf::Field *clipped_vtx = apf::createLagrangeField(m, "iso_clipped", apf::SCALAR, 1);
  
  //Get total number of elements
  int numel = 0;
  int nsd = m->getDimension();
  numel = m->count(nsd);
  PCU_Add_Ints(&numel, 1);

  //if target error is not specified, choose one based on interface or based on equidistribution assumption
  if(target_error==0){
    if(m->findField("vof")!=NULL)
      getTargetError(m,errField,target_error,err_total);
    else
      target_error = err_total/sqrt(m->count(nsd));
  }
   
  // Get domain volume
  // should only need to be computed once unless geometry is complex
  if (domainVolume == 0)
  {
    double volTotal = 0.0;
    it = m->begin(nsd);
    while (reg = m->iterate(it))
    {
      volTotal += apf::measure(m, reg);
    }
    PCU_Add_Doubles(&volTotal, 1);
    domainVolume = volTotal;
    assert(domainVolume > 0);
  }
  //compute the new size field over elements
  it = m->begin(nsd);
  double err_curr = 0.0;
  double errRho_curr = 0.0;
  double errRho_target = target_error / sqrt(domainVolume);
  apf::Vector3 err_vect;
  while (reg = m->iterate(it))
  {
    double h_old;
    double h_new;
    element = apf::createMeshElement(m, reg);

    if (m->getDimension() == 2)
      h_old = sqrt(apf::measure(element) * 4 / sqrt(3));
    else
      //h_old = pow(apf::measure(element) * 6 * sqrt(2), 1.0 / 3.0); //edge of a regular tet
      h_old = apf::computeShortestHeightInTet(m,reg);
    //err_curr = apf::getScalar(vmsErrH1, reg, 0);
    err_curr = apf::getScalar(errField, reg, 0);
    //err_curr = err_vect[0];
    //errRho_curr = apf::getScalar(errRho_reg, reg, 0);
    //h_new = h_old*errRho_target/errRho_curr;
    //h_new = h_old*sqrt(apf::measure(element))/sqrt(domainVolume)*target_error/err_curr;
    //
    //error-to-size relationship should be different between anisotropic and isotropic cases
    //consider moving this to where size frames are computed to get aspect ratio info
    
/* temporarily turn off until this is tested
    if (adapt_type_config == "anisotropic")
      if(target_error/err_curr <= 1)
        h_new = h_old * pow((target_error / err_curr),2.0/(2.0*(1.0)+1.0)); //refinement
      else
        h_new = h_old * pow((target_error / err_curr),2.0/(2.0*(1.0)+3.0)); //coarsening
    else //isotropic
*/
      h_new = h_old * pow((target_error / err_curr),2.0/(2.0*(1.0)+nsd));

    apf::setScalar(errorSize_reg, reg, 0, h_new);
    apf::destroyMeshElement(element);
  }
  m->end(it);

  //Transfer size field from elements to vertices through averaging
  it = m->begin(0);
  while ((v = m->iterate(it)))
  {
    //averageToEntity(errorSize_reg, errorSize, v);
    //volumeAverageToEntity(errorSize_reg, errorSize, v);
    errorAverageToEntity(errorSize_reg, errorSize,errField, v);
    //minToEntity(errorSize_reg, errorSize, v);
  }
  m->end(it);


  //Get the anisotropic size frame
  if (adapt_type_config == "anisotropic")
  {
    if(comm_rank==0)
      std::cout<<"Entering anisotropic loop to compute size scales and frames\n";
    double eps_u = 0.002; //distance from the interface
    apf::Field *phif = m->findField("phi");
    apf::Field *gradphi = apf::recoverGradientByVolume(phif);
    apf::Field *grad2phi = apf::recoverGradientByVolume(gradphi);

    apf::Field *speedF = extractSpeed(m->findField("velocity"));
    apf::Field *gradSpeed = apf::recoverGradientByVolume(speedF);
    apf::Field *grad2Speed = apf::recoverGradientByVolume(gradSpeed);

    apf::Field *hess = computeHessianField(grad2phi);
    apf::Field *curves = getCurves(hess, gradphi);
    //apf::Field* metricf = computeGradPhiMetricField(gradphi,grad2phi,errorSize,eps_u);
    apf::Field *metricf = computeMetricField(gradphi, grad2phi, errorSize, eps_u,hPhi);
    //apf::Field *metricf = computeMetricField(gradSpeed, grad2Speed, errorSize, eps_u);
    apf::Field *frame_comps[3] = {apf::createLagrangeField(m, "frame_0", apf::VECTOR, 1), apf::createLagrangeField(m, "frame_1", apf::VECTOR, 1), apf::createLagrangeField(m, "frame_2", apf::VECTOR, 1)};
    //getERMSizeFrames(metricf, gradSpeed, frame_comps);

    //Set the size scale for vertices
    it = m->begin(0);
    apf::Vector3 scale;
    while ((v = m->iterate(it)))
    {
      double tempScale = apf::getScalar(errorSize, v, 0);
      if (tempScale < hmin)
        apf::setScalar(clipped_vtx, v, 0, -1);
      else if (tempScale > hmax)
        apf::setScalar(clipped_vtx, v, 0, 1);
      else
        apf::setScalar(clipped_vtx, v, 0, 0);
      clamp(tempScale, hmin, hmax);
      apf::setScalar(errorSize,v,0,tempScale);
    }
    it = m->begin(0);
    while( (v = m->iterate(it)) ){
      double phi = apf::getScalar(phif, v, 0);
      apf::Vector3 curve;
      apf::getVector(curves, v, 0, curve);

      //metricf is the hessian
      apf::Matrix3x3 metric;
      apf::getMatrix(metricf, v, 0, metric);

      apf::Vector3 eigenVectors[3];
      double eigenValues[3];
      apf::eigen(metric, eigenVectors, eigenValues);
      // Sort the eigenvalues and corresponding vectors
      // Larger eigenvalues means a need for a finer mesh
      SortingStruct ssa[3];
      for (int i = 0; i < 3; ++i)
      {
        ssa[i].v = eigenVectors[i];
        ssa[i].wm = std::fabs(eigenValues[i]);
      }
      std::sort(ssa, ssa + 3);

      assert(ssa[2].wm >= ssa[1].wm);
      assert(ssa[1].wm >= ssa[0].wm);

      double lambda[3] = {ssa[2].wm, ssa[1].wm, ssa[0].wm};

      //anisotropic scales based on eigenvalues of metric
      //scaleFormulaERM(phi, hPhi, hmax, apf::getScalar(errorSize, v, 0), curve, lambda, eps_u, scale,nsd,maxAspect);
      scaleFormulaPhiVMS(phi, hPhi, hmax, apf::getScalar(errorSize, v, 0), curve, lambda, eps_u, scale,nsd,maxAspect);

      apf::setVector(size_scale, v, 0, scale);
      //get frames

      apf::Matrix3x3 frame(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);

      //get eigen values and eigenvectors from hessian
      double firstEigenvalue = ssa[2].wm;
      assert(firstEigenvalue > 1e-12);
      frame[0] = ssa[2].v;
      frame[1] = ssa[1].v;
      frame[2] = ssa[0].v;
    
      //normalize eigenvectors
      for (int i = 0; i < 3; ++i)
        frame[i] = frame[i].normalize();
      frame = apf::transpose(frame);
      apf::setMatrix(size_frame, v, 0, frame);
    }
    m->end(it);

    //Do simple size and aspect ratio grading
    /*
    gradeAnisoMesh(m,gradingFactor);
    if(comm_rank==0)
      std::cout<<"Finished grading size 0\n";
    gradeAspectRatio(m,1,gradingFactor);
    if(comm_rank==0)
      std::cout<<"Finished grading size 1\n";
    gradeAspectRatio(m,2,gradingFactor);
    if(comm_rank==0)
      std::cout<<"Finished grading size 2\n";
    */

    apf::synchronize(size_scale);

    apf::destroyField(gradphi);
    apf::destroyField(grad2phi);
    apf::destroyField(curves);
    apf::destroyField(hess);

/*
    if(logging_config=="on"){
      char namebuffer[20];
      sprintf(namebuffer,"pumi_preadapt_aniso_%i",nAdapt);
      apf::writeVtkFiles(namebuffer, m);
    }
*/

    apf::destroyField(metricf);
    apf::destroyField(frame_comps[0]);
    apf::destroyField(frame_comps[1]);
    apf::destroyField(frame_comps[2]);
    apf::destroyField(speedF);
    apf::destroyField(gradSpeed);
    apf::destroyField(grad2Speed);
    apf::destroyField(errorSize); //should be moved over to an intersecting function when available
  }
  else
  {
    it = m->begin(0);
    while ((v = m->iterate(it)))
    {
      double tempScale = apf::getScalar(errorSize, v, 0);
      if (tempScale < hmin)
        apf::setScalar(clipped_vtx, v, 0, -1);
      else if (tempScale > hmax)
        apf::setScalar(clipped_vtx, v, 0, 1);
      else
        apf::setScalar(clipped_vtx, v, 0, 0);
      clamp(tempScale, hmin, hmax);
      apf::setScalar(errorSize, v, 0, tempScale);
    }
    //gradeMesh();
    apf::synchronize(errorSize);
    m->end(it);
    if (target_element_count != 0)
    {
      sam::scaleIsoSizeField(errorSize, target_element_count);
      clampField(errorSize, hmin, hmax);
      //gradeMesh();
      //SmoothField(errorSize);
    }
    sizeFieldList.push(errorSize);
  } 

  //Destroy locally required fields
  apf::destroyField(errorSize_reg);
  apf::destroyField(clipped_vtx);
  if(comm_rank==0)
    std::cout<<"Finished Size Field\n";
  return 0;
}

int MeshAdaptPUMIDrvr::testIsotropicSizeField()
//Function that tests MeshAdapt by generating an isotropic sizefield based on hmin
{ 
    freeField(size_iso);
    size_iso = apf::createLagrangeField(m, "proteus_size",apf::SCALAR,1);
    apf::MeshIterator* it = m->begin(0);
    apf::MeshEntity* v;
    while(v = m->iterate(it)){
      double phi = hmin;
      clamp(phi,hmin,hmax);
      apf::setScalar(size_iso,v,0,phi);
    }
    return 0;
}

int gradeSizeModify(apf::Mesh* m, double gradingFactor, 
    double size[2], apf::Adjacent edgAdjVert, 
    apf::Adjacent vertAdjEdg,
    std::queue<apf::MeshEntity*> &markedEdges,
    apf::MeshTag* isMarked,
    int fieldType,
    int vecPos, //which idx of sizeVec to modify
    int idxFlag)
{
    //Determine a switching scheme depending on which vertex needs a modification
    int idx1,idx2;
    if(idxFlag == 0){
      idx1=0;
      idx2=1;
    } 
    else{
      idx1=1;
      idx2 = 0;
    } 
    
    int marker[3] = {0,1,0}; 
    double marginVal = 0.01;
    int needsParallel=0;

    if(fieldType == apf::SCALAR){
      apf::Field* size_iso = m->findField("proteus_size");

      if(size[idx1]>(gradingFactor*size[idx2])*(1+marginVal))
      {
        if(m->isOwned(edgAdjVert[idx1]))
        {
          size[idx1] = gradingFactor*size[idx2];
          apf::setScalar(size_iso,edgAdjVert[idx1],0,size[idx1]);
          m->getAdjacent(edgAdjVert[idx1], 1, vertAdjEdg);
          for (std::size_t i=0; i<vertAdjEdg.getSize();++i){
            m->getIntTag(vertAdjEdg[i],isMarked,&marker[2]);
            //if edge is not already marked
            if(!marker[2]){
              m->setIntTag(vertAdjEdg[i],isMarked,&marker[1]);
              markedEdges.push(vertAdjEdg[i]);
            }
          }
        } //end isOwned
        else
        { //Pack information to owning processor
          needsParallel=1;
          apf::Copies remotes;
          m->getRemotes(edgAdjVert[idx1],remotes);
          double newSize = gradingFactor*size[idx2];
          int owningPart=m->getOwner(edgAdjVert[idx1]);
          PCU_COMM_PACK(owningPart, remotes[owningPart]);
          PCU_COMM_PACK(owningPart,newSize);
        }
      }

    }//end if apf::SCALAR
    else{
      apf::Field* size_scale = m->findField("proteus_size_scale");
      apf::Vector3 sizeVec;
      if(size[idx1]>(gradingFactor*size[idx2])*(1+marginVal)){
        size[idx1] = gradingFactor*size[idx2];
        apf::getVector(size_scale,edgAdjVert[idx1],0,sizeVec);
        if(vecPos > 0){
          sizeVec[vecPos] = size[idx1]*sizeVec[0]; //realize the new aspect ratio
        }
        else{
          sizeVec[0] = size[idx1];
        }
        apf::setVector(size_scale,edgAdjVert[idx1],0,sizeVec);
        m->getAdjacent(edgAdjVert[idx1], 1, vertAdjEdg);
        for (std::size_t i=0; i<vertAdjEdg.getSize();++i){
          m->getIntTag(vertAdjEdg[i],isMarked,&marker[2]);
          //if edge is not already marked
          if(!marker[2]){
            m->setIntTag(vertAdjEdg[i],isMarked,&marker[1]);
            markedEdges.push(vertAdjEdg[i]);
          }
        }
      }
    }
  return needsParallel;
}

void markEdgesInitial(apf::Mesh* m, std::queue<apf::MeshEntity*> &markedEdges,double gradingFactor)
{
  //marker structure for 0) not marked 1) marked 2)storage
  int marker[3] = {0,1,0}; 

  double size[2];
  apf::MeshTag* isMarked = m->findTag("isMarked");
  apf::Field* size_iso = m->findField("proteus_size");
  apf::Adjacent edgAdjVert;
  apf::MeshEntity* edge;
  apf::MeshIterator* it = m->begin(1);
  while((edge=m->iterate(it))){
    m->getAdjacent(edge, 0, edgAdjVert);
    for (std::size_t i=0; i < edgAdjVert.getSize(); ++i){
      size[i]=apf::getScalar(size_iso,edgAdjVert[i],0);
    }
    if( (size[0] > gradingFactor*size[1]) || (size[1] > gradingFactor*size[0]) ){
      //add edge to a queue 
      markedEdges.push(edge);
      //tag edge to indicate that it is part of queue 
      m->setIntTag(edge,isMarked,&marker[1]); 
    }
    else{
      m->setIntTag(edge,isMarked,&marker[0]); 
    }
  }
  m->end(it); 
}

int serialGradation(apf::Mesh* m, std::queue<apf::MeshEntity*> &markedEdges,double gradingFactor)
{
  double size[2];
  //marker structure for 0) not marked 1) marked 2)storage
  int marker[3] = {0,1,0}; 
  apf::MeshTag* isMarked = m->findTag("isMarked");
  apf::Field* size_iso = m->findField("proteus_size");
  apf::Adjacent edgAdjVert;
  apf::Adjacent vertAdjEdg;
  apf::MeshEntity* edge;
  apf::MeshIterator* it = m->begin(1);
  int needsParallel=0;

  //perform serial gradation while packing necessary info for parallel
  while(!markedEdges.empty()){ 
    edge = markedEdges.front();
    m->getAdjacent(edge, 0, edgAdjVert);
    for (std::size_t i=0; i < edgAdjVert.getSize(); ++i){
      size[i] = apf::getScalar(size_iso,edgAdjVert[i],0);
    }

    needsParallel+=gradeSizeModify(m, gradingFactor, size, edgAdjVert, 
      vertAdjEdg, markedEdges, isMarked, apf::SCALAR,0, 0);
    needsParallel+=gradeSizeModify(m, gradingFactor, size, edgAdjVert, 
      vertAdjEdg, markedEdges, isMarked, apf::SCALAR,0, 1);

    m->setIntTag(edge,isMarked,&marker[0]);
    markedEdges.pop();
  }
  return needsParallel;
}

int MeshAdaptPUMIDrvr::gradeMesh()
//Function to grade isotropic mesh through comparison of edge vertex size ratios
//For simplicity, we do not bother with accounting for entities across partitions
{
  //
  if(comm_rank==0)
    std::cout<<"Starting grading\n";
  apf::MeshEntity* edge;
  apf::Adjacent edgAdjVert;
  apf::Adjacent vertAdjEdg;
  //double gradingFactor = 1.5;
  double size[2];
  std::queue<apf::MeshEntity*> markedEdges;
  apf::MeshTag* isMarked = m->createIntTag("isMarked",1);

  //marker structure for 0) not marked 1) marked 2)storage
  int marker[3] = {0,1,0}; 

  apf::MeshIterator* it;
  markEdgesInitial(m,markedEdges,gradingFactor);

  int needsParallel=1;
  int nCount=1;
  while(needsParallel)
  {
    PCU_Comm_Begin();
    needsParallel = serialGradation(m,markedEdges,gradingFactor);

    PCU_Add_Ints(&needsParallel,1);
    if(comm_rank==0)
      std::cerr<<"Sending size info for gradation"<<std::endl;
    PCU_Comm_Send(); 

    apf::MeshEntity* ent;
    double receivedSize;
    double currentSize;
    double newSize;

    //Need a container to get all entitites that need to be updated on remotes
    std::queue<apf::MeshEntity*> updateRemoteVertices;

    apf::Copies remotes;
    //owning copies are receiving
    while(PCU_Comm_Receive())
    {
      PCU_COMM_UNPACK(ent);
      PCU_COMM_UNPACK(receivedSize);

      if(!m->isOwned(ent)){
        std::cout<<"THERE WAS AN ERROR"<<std::endl;
        std::exit(1);
      }

      currentSize = apf::getScalar(size_iso,ent,0);
      newSize = std::min(receivedSize,currentSize);
      apf::setScalar(size_iso,ent,0,newSize);
      
      //add adjacent edges into Q
      m->getAdjacent(ent, 1, vertAdjEdg);
      for (std::size_t i=0; i<vertAdjEdg.getSize();++i)
      {
        edge = vertAdjEdg[i];
        m->getIntTag(vertAdjEdg[i],isMarked,&marker[2]);
        if(!marker[2])
        {
          markedEdges.push(edge);
          //tag edge to indicate that it is part of queue 
          m->setIntTag(edge,isMarked,&marker[1]);
        }
      }
      updateRemoteVertices.push(ent);
    }

    PCU_Comm_Begin();

    while(!updateRemoteVertices.empty())
    { 
      ent = updateRemoteVertices.front();
      //get remote copies and send updated mesh sizes
      m->getRemotes(ent,remotes);
      currentSize = apf::getScalar(size_iso,ent,0);
      for(apf::Copies::iterator iter=remotes.begin(); iter!=remotes.end();++iter)
      {
        PCU_COMM_PACK(iter->first, iter->second);
      }
      updateRemoteVertices.pop();
    }

    PCU_Comm_Send();
    //while remote copies are receiving
    while(PCU_Comm_Receive())
    {
      //unpack
      PCU_COMM_UNPACK(ent);
      //PCU_COMM_UNPACK(receivedSize);
      assert(!m->isOwned(ent));

      if(m->isOwned(ent)){
        std::cout<<"Problem occurred\n";
        std::exit(1);
      }
/*
      currentSize = apf::getScalar(size_iso,ent,0);
      newSize = std::min(receivedSize,currentSize);
      apf::setScalar(size_iso,ent,0,newSize);
*/

      //add adjacent edges into Q
      m->getAdjacent(ent, 1, vertAdjEdg);
      for (std::size_t i=0; i<vertAdjEdg.getSize();++i)
      {
        edge = vertAdjEdg[i];
        m->getIntTag(vertAdjEdg[i],isMarked,&marker[2]);
        if(!marker[2])
        {
          markedEdges.push(edge);
          //tag edge to indicate that it is part of queue 
          m->setIntTag(edge,isMarked,&marker[1]);
        }
      }
    }
    apf::synchronize(size_iso);

  } //end outer while

  //Cleanup of edge marker field
  it = m->begin(1);
  while((edge=m->iterate(it))){
    m->removeTag(edge,isMarked);
  }
  m->end(it); 
  m->destroyTag(isMarked);

  //apf::synchronize(size_iso);
  if(comm_rank==0)
    std::cout<<"Completed grading\n";
  return needsParallel;
}

void gradeAnisoMesh(apf::Mesh* m)
//Function to grade anisotropic mesh through comparison of edge vertex aspect ratios and minimum sizes
//For simplicity, we do not bother with accounting for entities across partitions
{
  //
  //if(comm_rank==0)
  //  std::cout<<"Starting anisotropic grading\n";
  apf::MeshIterator* it = m->begin(1);
  apf::MeshEntity* edge;
  apf::Adjacent edgAdjVert;
  apf::Adjacent vertAdjEdg;
  double gradingFactor = 1.3;
  double size[2];
  apf::Vector3 sizeVec;
  std::queue<apf::MeshEntity*> markedEdges;
  apf::MeshTag* isMarked = m->createIntTag("isMarked",1);
  apf::Field* size_scale = m->findField("proteus_size_scale");

  //marker structure for 0) not marked 1) marked 2)storage
  int marker[3] = {0,1,0}; 

  while((edge=m->iterate(it))){
    m->getAdjacent(edge, 0, edgAdjVert);
    for (std::size_t i=0; i < edgAdjVert.getSize(); ++i){
      apf::getVector(size_scale,edgAdjVert[i],0,sizeVec);
      size[i]=sizeVec[0];
    }
    if( (size[0] > gradingFactor*size[1]) || (size[1] > gradingFactor*size[0]) ){
      //add edge to a queue 
      markedEdges.push(edge);
      //tag edge to indicate that it is part of queue 
      m->setIntTag(edge,isMarked,&marker[1]); 

    }
    else{
      m->setIntTag(edge,isMarked,&marker[0]); 
    }
  }
  m->end(it); 
  while(!markedEdges.empty()){
    edge = markedEdges.front();
    m->getAdjacent(edge, 0, edgAdjVert);
    for (std::size_t i=0; i < edgAdjVert.getSize(); ++i){
      apf::getVector(size_scale,edgAdjVert[i],0,sizeVec);
      size[i]=sizeVec[0];
    }
    gradeSizeModify(m, gradingFactor, size, edgAdjVert, 
      vertAdjEdg, markedEdges, isMarked, apf::VECTOR,0, 0);
    gradeSizeModify(m, gradingFactor, size, edgAdjVert, 
      vertAdjEdg, markedEdges, isMarked, apf::VECTOR,0, 1);

/*
    if(size[0]>gradingFactor*size[1]){
      size[0] = gradingFactor*size[1];
      apf::getVector(size_scale,edgAdjVert[0],0,sizeVec);
      sizeVec[0] = size[0];
      apf::setVector(size_scale,edgAdjVert[0],0,sizeVec);
      m->getAdjacent(edgAdjVert[0], 1, vertAdjEdg);
      for (std::size_t i=0; i<vertAdjEdg.getSize();++i){
        m->getIntTag(vertAdjEdg[i],isMarked,&marker[2]);
        //if edge is not already marked
        if(!marker[2]){
          m->setIntTag(vertAdjEdg[i],isMarked,&marker[1]);
          markedEdges.push(vertAdjEdg[i]);
        }
      }
    }
    if(size[1]>gradingFactor*size[0]){
      size[1] = gradingFactor*size[0];
      apf::getVector(size_scale,edgAdjVert[1],0,sizeVec);
      sizeVec[0] = size[1];
      apf::setVector(size_scale,edgAdjVert[1],0,sizeVec);
      m->getAdjacent(edgAdjVert[1], 1, vertAdjEdg);
      for (std::size_t i=0; i<vertAdjEdg.getSize();++i){
        m->getIntTag(vertAdjEdg[i],isMarked,&marker[2]);
        //if edge is not already marked
        if(!marker[2]){
          m->setIntTag(vertAdjEdg[i],isMarked,&marker[1]);
          markedEdges.push(vertAdjEdg[i]);
        }
      }
    }
*/
    m->setIntTag(edge,isMarked,&marker[0]);
    markedEdges.pop();
  }
  it = m->begin(1);
  while((edge=m->iterate(it))){
    m->removeTag(edge,isMarked);
  }
  m->end(it); 
  m->destroyTag(isMarked);
  apf::synchronize(size_scale);
  //if(comm_rank==0)
  //  std::cout<<"Completed minimum size grading\n";
}

void gradeAspectRatio(apf::Mesh* m,int idx)
//Function to grade anisotropic mesh through comparison of edge vertex aspect ratios and minimum sizes
//For simplicity, we do not bother with accounting for entities across partitions
{
  std::cout<<"Entered function\n"; 
  apf::MeshIterator* it = m->begin(1);
  apf::MeshEntity* edge;
  apf::Adjacent edgAdjVert;
  apf::Adjacent vertAdjEdg;
  double gradingFactor = 1.3;
  double size[2];
  apf::Vector3 sizeVec;
  std::queue<apf::MeshEntity*> markedEdges;
  apf::MeshTag* isMarked = m->createIntTag("isMarked",1);
  apf::Field* size_scale = m->findField("proteus_size_scale");

  //marker structure for 0) not marked 1) marked 2)storage
  int marker[3] = {0,1,0}; 

  while((edge=m->iterate(it))){
    m->getAdjacent(edge, 0, edgAdjVert);
    for (std::size_t i=0; i < edgAdjVert.getSize(); ++i){
      apf::getVector(size_scale,edgAdjVert[i],0,sizeVec);
      size[i]=sizeVec[idx]/sizeVec[0];
    }
    if( (size[0] > gradingFactor*size[1]) || (size[1] > gradingFactor*size[0]) ){
      //add edge to a queue 
      markedEdges.push(edge);
      //tag edge to indicate that it is part of queue 
      m->setIntTag(edge,isMarked,&marker[1]); 

    }
    else{
      m->setIntTag(edge,isMarked,&marker[0]); 
    }
  }
  m->end(it); 

  std::cout<<"Got queue of size "<<markedEdges.size()<<std::endl; 
  while(!markedEdges.empty()){
    edge = markedEdges.front();
    m->getAdjacent(edge, 0, edgAdjVert);
    for (std::size_t i=0; i < edgAdjVert.getSize(); ++i){
      apf::getVector(size_scale,edgAdjVert[i],0,sizeVec);
      size[i]=sizeVec[idx]/sizeVec[0];
    }
    gradeSizeModify(m, gradingFactor, size, edgAdjVert, 
      vertAdjEdg, markedEdges, isMarked, apf::VECTOR, idx, 0);
    gradeSizeModify(m, gradingFactor, size, edgAdjVert, 
      vertAdjEdg, markedEdges, isMarked, apf::VECTOR, idx, 1);

    m->setIntTag(edge,isMarked,&marker[0]);
    markedEdges.pop();
  }
  it = m->begin(1);
  while((edge=m->iterate(it))){
    m->removeTag(edge,isMarked);
  }
  m->end(it); 
  m->destroyTag(isMarked);
  apf::synchronize(size_scale);
}

void gradeAnisoMesh(apf::Mesh* m,double gradingFactor)
//Function to grade anisotropic mesh through comparison of edge vertex aspect ratios and minimum sizes
//For simplicity, we do not bother with accounting for entities across partitions
{
  //
  //if(comm_rank==0)
  //  std::cout<<"Starting anisotropic grading\n";
  apf::MeshIterator* it = m->begin(1);
  apf::MeshEntity* edge;
  apf::Adjacent edgAdjVert;
  apf::Adjacent vertAdjEdg;
  //double gradingFactor = 1.3;
  double size[2];
  apf::Vector3 sizeVec;
  std::queue<apf::MeshEntity*> markedEdges;
  apf::MeshTag* isMarked = m->createIntTag("isMarked",1);
  apf::Field* size_scale = m->findField("proteus_size_scale");

  //marker structure for 0) not marked 1) marked 2)storage
  int marker[3] = {0,1,0}; 

  while((edge=m->iterate(it))){
    m->getAdjacent(edge, 0, edgAdjVert);
    for (std::size_t i=0; i < edgAdjVert.getSize(); ++i){
      apf::getVector(size_scale,edgAdjVert[i],0,sizeVec);
      size[i]=sizeVec[0];
    }
    if( (size[0] > gradingFactor*size[1]) || (size[1] > gradingFactor*size[0]) ){
      //add edge to a queue 
      markedEdges.push(edge);
      //tag edge to indicate that it is part of queue 
      m->setIntTag(edge,isMarked,&marker[1]); 

    }
    else{
      m->setIntTag(edge,isMarked,&marker[0]); 
    }
  }
  m->end(it); 
  while(!markedEdges.empty()){
    edge = markedEdges.front();
    m->getAdjacent(edge, 0, edgAdjVert);
    for (std::size_t i=0; i < edgAdjVert.getSize(); ++i){
      apf::getVector(size_scale,edgAdjVert[i],0,sizeVec);
      size[i]=sizeVec[0];
    }
    gradeSizeModify(m, gradingFactor, size, edgAdjVert, 
      vertAdjEdg, markedEdges, isMarked, apf::VECTOR,0, 0);
    gradeSizeModify(m, gradingFactor, size, edgAdjVert, 
      vertAdjEdg, markedEdges, isMarked, apf::VECTOR,0, 1);

/*
    if(size[0]>gradingFactor*size[1]){
      size[0] = gradingFactor*size[1];
      apf::getVector(size_scale,edgAdjVert[0],0,sizeVec);
      sizeVec[0] = size[0];
      apf::setVector(size_scale,edgAdjVert[0],0,sizeVec);
      m->getAdjacent(edgAdjVert[0], 1, vertAdjEdg);
      for (std::size_t i=0; i<vertAdjEdg.getSize();++i){
        m->getIntTag(vertAdjEdg[i],isMarked,&marker[2]);
        //if edge is not already marked
        if(!marker[2]){
          m->setIntTag(vertAdjEdg[i],isMarked,&marker[1]);
          markedEdges.push(vertAdjEdg[i]);
        }
      }
    }
    if(size[1]>gradingFactor*size[0]){
      size[1] = gradingFactor*size[0];
      apf::getVector(size_scale,edgAdjVert[1],0,sizeVec);
      sizeVec[0] = size[1];
      apf::setVector(size_scale,edgAdjVert[1],0,sizeVec);
      m->getAdjacent(edgAdjVert[1], 1, vertAdjEdg);
      for (std::size_t i=0; i<vertAdjEdg.getSize();++i){
        m->getIntTag(vertAdjEdg[i],isMarked,&marker[2]);
        //if edge is not already marked
        if(!marker[2]){
          m->setIntTag(vertAdjEdg[i],isMarked,&marker[1]);
          markedEdges.push(vertAdjEdg[i]);
        }
      }
    }
*/
    m->setIntTag(edge,isMarked,&marker[0]);
    markedEdges.pop();
  }
  it = m->begin(1);
  while((edge=m->iterate(it))){
    m->removeTag(edge,isMarked);
  }
  m->end(it); 
  m->destroyTag(isMarked);
  apf::synchronize(size_scale);
  //if(comm_rank==0)
  //  std::cout<<"Completed minimum size grading\n";
}

void gradeAspectRatio(apf::Mesh* m,int idx,double gradingFactor)
//Function to grade anisotropic mesh through comparison of edge vertex aspect ratios and minimum sizes
//For simplicity, we do not bother with accounting for entities across partitions
{
  if(PCU_Comm_Self()==0)
    std::cout<<"Entered function\n"; 
  apf::MeshIterator* it = m->begin(1);
  apf::MeshEntity* edge;
  apf::Adjacent edgAdjVert;
  apf::Adjacent vertAdjEdg;
  //double gradingFactor = 1.3;
  double size[2];
  apf::Vector3 sizeVec;
  std::queue<apf::MeshEntity*> markedEdges;
  apf::MeshTag* isMarked = m->createIntTag("isMarked",1);
  apf::Field* size_scale = m->findField("proteus_size_scale");

  //marker structure for 0) not marked 1) marked 2)storage
  int marker[3] = {0,1,0}; 

  while((edge=m->iterate(it))){
    m->getAdjacent(edge, 0, edgAdjVert);
    for (std::size_t i=0; i < edgAdjVert.getSize(); ++i){
      apf::getVector(size_scale,edgAdjVert[i],0,sizeVec);
      size[i]=sizeVec[idx]/sizeVec[0];
    }
    if( (size[0] > gradingFactor*size[1]) || (size[1] > gradingFactor*size[0]) ){
      //add edge to a queue 
      markedEdges.push(edge);
      //tag edge to indicate that it is part of queue 
      m->setIntTag(edge,isMarked,&marker[1]); 

    }
    else{
      m->setIntTag(edge,isMarked,&marker[0]); 
    }
  }
  m->end(it); 

  if(PCU_Comm_Self()==0)
    std::cout<<"Got queue of size "<<markedEdges.size()<<std::endl; 
  while(!markedEdges.empty()){
    edge = markedEdges.front();
    m->getAdjacent(edge, 0, edgAdjVert);
    for (std::size_t i=0; i < edgAdjVert.getSize(); ++i){
      apf::getVector(size_scale,edgAdjVert[i],0,sizeVec);
      size[i]=sizeVec[idx]/sizeVec[0];
    }
    gradeSizeModify(m, gradingFactor, size, edgAdjVert, 
      vertAdjEdg, markedEdges, isMarked, apf::VECTOR, idx, 0);
    gradeSizeModify(m, gradingFactor, size, edgAdjVert, 
      vertAdjEdg, markedEdges, isMarked, apf::VECTOR, idx, 1);

    m->setIntTag(edge,isMarked,&marker[0]);
    markedEdges.pop();
  }
  it = m->begin(1);
  while((edge=m->iterate(it))){
    m->removeTag(edge,isMarked);
  }
  m->end(it); 
  m->destroyTag(isMarked);
  apf::synchronize(size_scale);
}



double getMetricLength(apf::Matrix3x3 metric, apf::Vector3 directionVector)
{
  double length =0.0;
  length = directionVector*(metric*directionVector);
  return sqrt(length);
}

//void getGeometricMetricLength(apf::Matrix3x3 metric1, apf::Matrix3x3 metric2, apf::Vector3 directionVector,double &metricLength1, double &metricLength2)
double getGeometricMetricLength(apf::Matrix3x3 metric1, apf::Matrix3x3 metric2, apf::Vector3 directionVector)
//Compute geometric interpolated metric length of edge w/ respect to both end points
{
  double edgeLength = directionVector.getLength();
  double l1 = edgeLength/getMetricLength(metric1,directionVector);
  double l2 = edgeLength/getMetricLength(metric2,directionVector);
  double a = l1/l2;
  //double a2 = l2/l1;

  double metricLength;
  if(fabs(a-1.0)<1e-13)
    metricLength = l1;
  else
    metricLength = l1*(a-1)/(a*std::log(a));
  //metricLength2 = l2*(a2-1)/(a2*std::log(a2));

/*
  std::cout<<"L1 "<<l1<<std::endl;
  std::cout<<"L2 "<<l2<<std::endl;
  std::cout<<"a "<<a<<std::endl;
  std::cout<<"a2 "<<a2<<std::endl;
  std::cout<<"metricLength 1 "<<metricLength1<<std::endl;
  std::cout<<"metricLength 2 "<<metricLength2<<std::endl;
*/
  return metricLength;
}

apf::Matrix3x3 getMetricIntersection(apf::Matrix3x3 metric1, apf::Matrix3x3 metric2)
{

  //if metric1 = metric 2 then the intersection is itself
 
  apf::Matrix3x3 diff = metric1-metric2;
  if(sqrt(apf::getInnerProduct(diff,diff)) < 1e-12)
    return metric1;

  apf::Matrix3x3 N = apf::invert(metric1)*metric2;
  //std::cout<<"N "<<N<<std::endl;
  apf::Vector3 eigenVectors[3];
  double eigenValues[3];
  Eigen::MatrixXd A(3,3);
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      A(i,j) = N[i][j];
    }
  }
  Eigen::EigenSolver<Eigen::MatrixXd> es(A);
/*
  std::cout << "Here is a random 6x6 matrix, A:" << std::endl << A << std::endl << std::endl;
  std::cout << "The eigenvalues of A are:" << std::endl << es.eigenvalues() << std::endl;
  std::cout << "The matrix of eigenvectors, V, is:" << std::endl << es.eigenvectors() << std::endl << std::endl;
*/

  for(int i=0;i<3;i++)
  {
    eigenValues[i] = es.eigenvalues()[i].real();
//    std::cout<<"eigenvalues "<<eigenValues[i]<<std::endl;
    for(int j=0;j<3;j++)
    {
      eigenVectors[i][j] = es.eigenvectors().col(i)[j].real();
    }
  }
//  std::cout<<"eigenVectors "<<eigenVectors[0]<<" "<<eigenVectors[1]<<" "<<eigenVectors[2]<<std::endl;

  //apf::eigen(N, eigenVectors, eigenValues);
//  std::cout<<"Past the eigen call\n";
 
  apf::Matrix3x3 R; //this will be used to store the eigenvectors
  apf::Matrix3x3 R_t; //transpose

  //This loads eigenvectors as rows, but it should be columns, so we transpose
  for(int i=0;i<3;i++)
    R_t[i] = eigenVectors[i];

  R = apf::transpose(R_t); 

  apf::Matrix3x3 finalMat;
  apf::Matrix3x3 lambdaMat = R_t*metric1*R;
  apf::Matrix3x3 muMat = R_t*metric2*R;
  for(int i=0;i<3;i++) 
  {
    for(int j=0;j<3;j++)
    {
      if(i == j)
        finalMat[i][j] = std::max(std::fabs(lambdaMat[i][j]),std::fabs(muMat[i][j]));
      else
        finalMat[i][j] = 0.0;
    }
  }

  apf::Matrix3x3 intersect =apf::transpose(apf::invert(R))*finalMat*apf::invert(R);
/*
  std::cout<<"N "<<N<<std::endl;
  std::cout<<"R "<<R<<std::endl;
  std::cout<<"eigenVectors "<<eigenVectors[0]<<std::endl;
  std::cout<<"metric 1 "<<metric1<<std::endl;
  std::cout<<"lambdaMat "<<lambdaMat<<std::endl;
  std::cout<<"muMat "<<muMat<<std::endl;
  std::cout<<"finalMat "<<finalMat<<std::endl;
  std::cout<<"intersect "<<apf::transpose(apf::invert(R))*finalMat*apf::invert(R)<<std::endl;
*/
  return intersect;
}

int metricIntersection_core(apf::Mesh2* m, apf::MeshEntity* edge, apf::Matrix3x3 intersect[],double gradingFactor)
{
//  1. look at the metric on each vertex
//  2. compute the scaled metric based on metric length for each vertex
//  3. compute the intersections
//  4. compare the intersections with the original metric via Frobenius norm ; More generally, this is where we determine if an ellipsoid is in another ellipsoid

  apf::Matrix3x3 metric[2];
  apf::Matrix3x3 scaledMetric[2];
  apf::Field* metricField = m->findField("metric");
  apf::Adjacent edgAdjVert;

  //1.
  m->getAdjacent(edge, 0, edgAdjVert);
  for (std::size_t i=0; i < edgAdjVert.getSize(); ++i)
  {
    apf::getMatrix(metricField,edgAdjVert[i],0,metric[i]);
  }

  //2.

  //get the points to compute a direction for the edge
  apf::Vector3 pts[2];
  for (std::size_t i=0; i < edgAdjVert.getSize(); ++i)
  {
    m->getPoint(edgAdjVert[i],0,pts[i]);
  }
  apf::Vector3 directionVector = pts[1]-pts[0];
  double metricLength = getGeometricMetricLength(metric[0], metric[1], directionVector);
  double scalingFactor = 1.0+metricLength*std::log(gradingFactor);
  scalingFactor = scalingFactor*scalingFactor;
  scalingFactor = 1.0/scalingFactor;

  for (std::size_t i=0; i < edgAdjVert.getSize(); ++i)
  {
    scaledMetric[i]=metric[i]*scalingFactor;
  }

  //3.
  for (std::size_t i=0; i < edgAdjVert.getSize(); ++i)
  {
    if(i==0)
      intersect[0] = getMetricIntersection(metric[0],scaledMetric[1]);
    else if (i==1)
      intersect[1] = getMetricIntersection(metric[1],scaledMetric[0]);
  }

  //4.
  
  int needToIntersect = 0;
  for (std::size_t i=0; i < edgAdjVert.getSize(); ++i)
  {
    apf::Matrix3x3 diff =metric[i]-intersect[i];
    //std::cout<<"metric \n"<<metric[i]<<std::endl;
    //std::cout<<"intersect \n"<<intersect[i]<<std::endl;
    double diff_norm = sqrt(apf::getInnerProduct(diff,diff));
    if(diff_norm > 1e-10)
      needToIntersect++;
  }
  //std::abort();

  return needToIntersect;
}


void markEdgesInitialMetric(apf::Mesh2* m, std::queue<apf::MeshEntity*> &markedEdges,double gradingFactor)
//This function looks at each edge and determines whether or not applying simultaneous intersection will be necessary
//For each edge

//This will actually waste some resources making the comparisons, but this is simpler to understand
{
  //marker structure for 0) not marked 1) marked 2)storage
  int marker[2] = {0,1}; 

  apf::Matrix3x3 intersect[2];
  apf::MeshTag* isMarked = m->findTag("isMarked");
  apf::Adjacent edgAdjVert;
  apf::MeshEntity* edge;
  apf::MeshIterator* it = m->begin(1);
  while((edge=m->iterate(it)))
  {
    int needToIntersect = metricIntersection_core(m, edge, intersect,gradingFactor);
    if(needToIntersect>0) 
    {
      //add edge to a queue 
      markedEdges.push(edge);
      //tag edge to indicate that it is part of queue 
      m->setIntTag(edge,isMarked,&marker[1]); 
    }
    else
    {
      m->setIntTag(edge,isMarked,&marker[0]); 
    }
  } //end while
  m->end(it); 
}



int serialMetricGradation(apf::Mesh2* m, std::queue<apf::MeshEntity*> &markedEdges,double gradingFactor)
{
  //marker structure for 0) not marked 1) marked 2)storage
  int marker[2] = {0,1}; 
  int markerTag; //used to check marker

  apf::MeshTag* isMarked = m->findTag("isMarked");
  apf::Field* metricField = m->findField("metric");

  apf::Adjacent edgAdjVert;
  apf::Adjacent vertAdjEdg;
  apf::MeshEntity* edge;
  int needsParallel=0;
  apf::Matrix3x3 intersect[2];

  //perform serial gradation while packing necessary info for parallel
  std::cout<<"Enter queue loop\n";
  while(!markedEdges.empty())
  { 
    edge = markedEdges.front();
    m->getAdjacent(edge, 0, edgAdjVert);
    
    int needToIntersect = metricIntersection_core(m, edge, intersect,gradingFactor);
    //if need to intersect, set intersect as the metric
    if(needToIntersect>0)
    {

      for (std::size_t i=0; i < edgAdjVert.getSize(); ++i)
      {
        apf::setMatrix(metricField,edgAdjVert[i],0,intersect[i]);
      }

/*
      //Need to look at add adjacent edges to the queue
      for (std::size_t i=0; i < edgAdjVert.getSize(); ++i)
      {
        m->getAdjacent(edgAdjVert[i], 1, vertAdjEdg);
        //Loop over adjacent edges to vertex
        std::cout<<i<<" edgAdjVert size "<<edgAdjVert.getSize()<<" vertAdjEdg size "<<vertAdjEdg.getSize()<<std::endl;
        for (std::size_t j=0; i < vertAdjEdg.getSize(); ++j)
        {
          std::cout<<"Is this the problem?\n";
          std::cout<<"type is "<<m->getType(vertAdjEdg[j])<<" type of base "<<m->getType(edgAdjVert[i])<<std::endl;
          apf::Vector3 pts; 
          m->getPoint(edgAdjVert[i],0,pts);
          std::cout<<"points "<<pts<<"  "<<std::endl;

          if(m->getType(vertAdjEdg[j])==1)
          {
            m->getIntTag(vertAdjEdg[j],isMarked,&markerTag);
            std::cout<<"no?\n";
            if(!markerTag)
            {
              //push edge to queue and tag edge to indicate that it is part of queue 
              markedEdges.push(vertAdjEdg[j]);
              m->setIntTag(vertAdjEdg[j],isMarked,&marker[1]); 
            }
          }
        }
      } //end for - add adjacent edges to queue
*/
      apf::Adjacent edgAdjEdg;
      apf::getBridgeAdjacent(m,edge,0,1,edgAdjEdg);
      for (std::size_t i=0; i < edgAdjEdg.getSize(); ++i)
      {
        m->getIntTag(edgAdjEdg[i],isMarked,&markerTag);
        if(!markerTag)
        {
        //push edge to queue and tag edge to indicate that it is part of queue 
          markedEdges.push(edgAdjEdg[i]);
          m->setIntTag(edgAdjEdg[i],isMarked,&marker[1]); 
        }

      } //end for adjacent edges
      
    } //end if
  
    //Will need a section to determine if parallel communication is necessary
    m->setIntTag(edge,isMarked,&marker[0]); 
    markedEdges.pop();
    std::cout <<" queue size is "<<markedEdges.size()<<std::endl;
  }
  return needsParallel;
}



void MeshAdaptPUMIDrvr::gradeMetric()
//Based on Alauzet, Frédéric. "Size gradation control of anisotropic meshes." Finite Elements in Analysis and Design 46.1-2 (2010): 181-202.
{
  if(PCU_Comm_Self()==0)
    std::cout<<"Entered metric gradation function\n"; 

  std::queue<apf::MeshEntity*> markedEdges;
  apf::MeshTag* isMarked = m->createIntTag("isMarked",1);

  markEdgesInitialMetric(m,markedEdges,gradingFactor);
  std::cout<<"Past markEdgesInitialMetric\n";

  int needsParallel=1;
  int nCount=1;


  //initialize queue
  //mark edges that need intersection
  //while queue is non-empty
  //  Compute intersection
  //  Check via Frobenius Norm if metric need intersection on vertex 0
  //  If yes, replace metric with intersection and set flag to 1
  //  Check vertex 1
  //  

  while(needsParallel)
  {
    needsParallel = serialMetricGradation(m,markedEdges,gradingFactor);
  }
  std::cout<<"Past serialMetricGradation\n";
  //Cleanup of edge marker field
  apf::MeshEntity* edge;
  apf::MeshIterator* it;

  it = m->begin(1);
  while((edge=m->iterate(it)))
  {
    m->removeTag(edge,isMarked);
  }
  m->end(it); 
  m->destroyTag(isMarked);

}

void MeshAdaptPUMIDrvr::anisotropicIntersect()
{
  //create metric fields 
  
  apf::Field* metricField = apf::createLagrangeField(m, "metric", apf::MATRIX, 1);
  apf::Field* size_scales = m->findField("proteus_size_scale");
  apf::Field* size_frames = m->findField("proteus_size_frame");

  apf::Matrix3x3 metric;
  apf::Matrix3x3 frame;
  apf::Vector3 scale;

  apf::MeshEntity* v;
  apf::MeshIterator* it = m->begin(0);
  while( (v = m->iterate(it)) )
  {
    apf::getVector(size_scales,v,0,scale);
    apf::getMatrix(size_frames,v,0,frame);
    apf::Matrix3x3 scaleMatrix(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);
    scaleMatrix[0][0] = scale[0];
    scaleMatrix[1][1] = scale[1];
    scaleMatrix[2][2] = scale[2];
    
    metric = apf::transpose(frame)*scaleMatrix*frame;
    //symmetrize metric
    //metric = (metric+apf::transpose(metric))/2.0;
    apf::setMatrix(metricField,v,0,metric);
  }
  m->end(it);
  
  //intersect them - don't need this yet as we only have one size field description
  //
  char namebuffer[50];
  sprintf(namebuffer,"pumi_preGrade_%i",nAdapt);
  apf::writeVtkFiles(namebuffer, m);
 
  
  //grade the metric field
  gradeMetric();

  sprintf(namebuffer,"pumi_postGrade_%i",nAdapt);
  apf::writeVtkFiles(namebuffer, m);
  //std::abort();

  //deconstruct metric fields back into frame and size

  it = m->begin(0);
  while( (v = m->iterate(it)) )
  {
    apf::getMatrix(metricField,v,0,metric);
  
    //eigendecomposition
    //sort and set
    apf::Vector3 eigenVectors[3];
    double eigenValues[3];
    apf::eigen(metric, eigenVectors, eigenValues);
    SortingStruct ssa[3];
    for (int i = 0; i < 3; ++i)
    {
      ssa[i].v = eigenVectors[i];
      ssa[i].wm = std::fabs(eigenValues[i]);
    }
    std::sort(ssa, ssa + 3);

    for (int i = 0; i < 3; ++i)
    {
      scale[i] = ssa[2-i].wm;
      frame[i] = ssa[2-i].v;
      frame[i] = frame[i].normalize();
    }

    frame = apf::transpose(frame);


    apf::setVector(size_scales,v,0,scale);
    apf::setMatrix(size_frames,v,0,frame);
    

  }
  m->end(it);

  apf::destroyField(metricField);
}




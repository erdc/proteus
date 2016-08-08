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

static void SmoothField(apf::Field* f);

/* Based on the distance from the interface epsilon can be controlled to determine
   thickness of refinement near the interface */
static double isotropicFormula(double phi, double hmin, double hmax)
{
  static double const epsilon = 0.02;
  phi = fabs(phi);
  double size;
  if (fabs(phi) < epsilon)
    size = hmin;
  else if (phi < 3 * epsilon)
    size = (hmin + hmax) / 2;
  else
    size = hmax;
  return size;
}

int MeshAdaptPUMIDrvr::calculateSizeField()
{
  freeField(size_iso);
  size_iso = apf::createLagrangeField(m, "proteus_size",apf::SCALAR,1);
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* v;
  apf::Field* phif = m->findField("phi");
  assert(phif);
  while ((v = m->iterate(it))) {
    double phi = apf::getScalar(phif, v, 0);
    double size = isotropicFormula(phi, hmin, hmax);
    apf::setScalar(size_iso, v, 0, size);
  }
  m->end(it);
  for(int i=0; i < 3; i++)
    SmoothField(size_iso);
  return 0;
}


//taken from Dan's superconvergent patch recovery code
void MeshAdaptPUMIDrvr::averageToEntity(apf::Field* ef, apf::Field* vf,
    apf::MeshEntity* ent)
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

void MeshAdaptPUMIDrvr::volumeAverageToEntity(apf::Field* ef, apf::Field* vf,
    apf::MeshEntity* ent)
{
  apf::Mesh* m = apf::getMesh(ef);
  apf::Adjacent elements;
  apf::MeshElement* testElement;
  m->getAdjacent(ent, m->getDimension(), elements);
  double s=0;
  double invVolumeTotal=0;
  for (std::size_t i=0; i < elements.getSize(); ++i){
      testElement = apf::createMeshElement(m,elements[i]);
      s+= apf::getScalar(ef,elements[i],0)/apf::measure(testElement);
      invVolumeTotal += 1.0/apf::measure(testElement);
if(comm_rank==0){
  std::cout<<"What is s "<<s<<" Volume? "<<apf::measure(testElement)<<" scale? "<<apf::getScalar(ef,elements[i],0)<<std::endl;
}
      apf::destroyMeshElement(testElement);
  }
  s /= invVolumeTotal;
if(comm_rank==0){
  std::cout<<"What is s final? "<<s<<std::endl;
}
  apf::setScalar(vf, ent, 0, s);
  return;
}


static apf::Field* extractSpeed(apf::Field* velocity)
{
  apf::Mesh* m = apf::getMesh(velocity);
  apf::Field* speedF = apf::createLagrangeField(m,"proteus_speed",apf::SCALAR,1);
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* v;
  apf::Vector3 vel_vect;
  while ((v = m->iterate(it))) {
    apf::getVector(velocity, v, 0, vel_vect);
    double speed = vel_vect.getLength();
    apf::setScalar(speedF, v, 0, speed);
  }
  m->end(it);
  return speedF;
}

static apf::Matrix3x3 hessianFormula(apf::Matrix3x3 const& g2phi)
{
  apf::Matrix3x3 g2phit = apf::transpose(g2phi);
  //return (g2phi + g2phit) / 2;
  return g2phi;
}

static apf::Field* computeHessianField(apf::Field* grad2phi)
{
  apf::Mesh* m = apf::getMesh(grad2phi);
  apf::Field* hessf = createLagrangeField(m,"proteus_hess",apf::MATRIX,1);
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* v;
  while ((v = m->iterate(it))) {
    apf::Matrix3x3 g2phi;
    apf::getMatrix(grad2phi, v, 0, g2phi);
    apf::Matrix3x3 hess = hessianFormula(g2phi);
    apf::setMatrix(hessf, v, 0, hess);
  }
  m->end(it);
  return hessf;
}

static apf::Field* computeMetricField(apf::Field* gradphi, apf::Field*grad2phi,apf::Field* size_iso,double eps_u){
  apf::Mesh* m = apf::getMesh(grad2phi);
  apf::Field* metricf = createLagrangeField(m,"proteus_metric",apf::MATRIX,1);
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* v;
  while ((v = m->iterate(it))) {
    apf::Matrix3x3 g2phi;
    apf::getMatrix(grad2phi, v, 0, g2phi);
    apf::Vector3 gphi;
    apf::getVector(gradphi,v,0,gphi);
    apf::Matrix3x3 gphigphit(gphi[0]*gphi[0], gphi[0]*gphi[1], gphi[0]*gphi[2],
                             gphi[0]*gphi[1], gphi[1]*gphi[1], gphi[1]*gphi[2],
                             gphi[0]*gphi[2], gphi[1]*gphi[2], gphi[2]*gphi[2]); 
    apf::Matrix3x3 hess = hessianFormula(g2phi);
    apf::Matrix3x3 metric = hess;
    //apf::Matrix3x3 metric = gphigphit/(apf::getScalar(size_iso,v,0)*apf::getScalar(size_iso,v,0))+ hess/eps_u;
    //apf::Matrix3x3 metric(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);
    apf::setMatrix(metricf, v, 0, metric);
  }
  m->end(it);
  return metricf;
}

// Gaussian, Mean and principal curvatures based on Hessian and gradient of phi. 
static void curveFormula(apf::Matrix3x3 const& h, apf::Vector3 const& g,
    apf::Vector3& curve)
{
  double a =   (h[1][1] + h[2][2]) * g[0] * g[0]
             + (h[0][0] + h[2][2]) * g[1] * g[1]
             + (h[0][0] + h[1][1]) * g[2] * g[2];

  double b =   g[0] * g[1] * h[0][1]
             + g[0] * g[2] * h[0][2]
             + g[1] * g[2] * h[1][2];

  double Km = (a - 2 * b) / pow(g * g, 1.5);

  double c =   g[0] * g[0] * (h[1][1] * h[2][2] - h[1][2] * h[1][2])
             + g[1] * g[1] * (h[0][0] * h[2][2] - h[0][2] * h[0][2])
             + g[2] * g[2] * (h[0][0] * h[1][1] - h[0][1] * h[0][1]);

  double d =   g[0] * g[1] * (h[0][2] * h[1][2] - h[0][1] * h[2][2])
             + g[1] * g[2] * (h[0][1] * h[0][2] - h[1][2] * h[0][0])
             + g[0] * g[2] * (h[0][1] * h[1][2] - h[0][2] * h[1][1]);

  double Kg = (c + 2 * d) / pow(g * g, 2);

  curve[0] = Km;                  //Mean curvature= (k_1+k_2)/2, Not used in normal direction
  curve[1] = Km+ sqrt(Km*Km- Kg); //k_1, First principal curvature  (maximum curvature)
  curve[2] = Km- sqrt(Km*Km- Kg); //k_2, Second principal curvature (minimum curvature)
}

static apf::Field* getCurves(apf::Field* hessians, apf::Field* gradphi)
{
  apf::Mesh* m = apf::getMesh(hessians);
  apf::Field* curves;
  curves = apf::createLagrangeField(m, "proteus_curves", apf::VECTOR, 1);
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* v;
  while ((v = m->iterate(it))) {
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

static void clamp(double& v, double min, double max)
{
  v = std::min(max, v);
  v = std::max(min, v);
}

static void scaleFormula(double phi, double hmin, double hmax,
    int adapt_step,
    apf::Vector3 const& curves,
    apf::Vector3& scale)
{
  double epsilon = 7.0* hmin; 
  if (fabs(phi) < epsilon) {
     scale[0] = hmin;
     scale[1] = sqrt(0.002/ fabs(curves[1]));
     scale[2] = sqrt(0.002/ fabs(curves[2]));
  }else if(fabs(phi) < 4 * epsilon){
     scale[0] = 2 * hmin;
     scale[1] = 2 * sqrt(0.002/ fabs(curves[1])); 
     scale[2] = 2 * sqrt(0.002/ fabs(curves[2])); 
  }else{
     scale = apf::Vector3(1,1,1) * hmax; 
 }	  


  for (int i = 0; i < 3; ++i)
    clamp(scale[i], hmin, hmax);
}

static void scaleFormulaERM(double phi, double hmin, double hmax, double h_dest, 
                            apf::Vector3 const& curves,
                            double lambda[3], double eps_u, apf::Vector3& scale,std::string adapt_type)
{

  if(adapt_type=="isotropic"){
    scale = apf::Vector3(1,1,1) * h_dest;
    for(int i=0;i<3;i++)
      clamp(scale[i], hmin, hmax);
  }
  else if(adapt_type=="anisotropic") { 
    double epsilon = 7.0* hmin; 
    double lambdamin = 1.0/(hmin*hmin);
    if(lambda[1] < 1e-10){lambda[1]=lambdamin; lambda[2]=lambdamin;}
    if(lambda[2] < 1e-10){lambda[2]=lambdamin;}
///* useful
    scale[0] = h_dest*pow((lambda[1]*lambda[2])/(lambda[0]*lambda[0]),1.0/6.0);
    scale[1] = sqrt(lambda[0]/lambda[1])*scale[0];
    scale[2] = sqrt(lambda[0]/lambda[2])*scale[0];
//*/
/*
    if(fabs(phi)<epsilon){
      scale[0] = h_dest;
      scale[1] = sqrt(eps_u/fabs(curves[2]));
      scale[2] = sqrt(eps_u/fabs(curves[1]));
    }
    else
      scale = apf::Vector3(1,1,1) * h_dest; 
*/
  }
  else{
    std::cerr << "unknown adapt type config " << adapt_type << '\n';
    abort();
  }
} 


static apf::Field* getSizeScales(apf::Field* phif, apf::Field* curves,
    double hmin, double hmax, int adapt_step)
{
  apf::Mesh* m = apf::getMesh(phif);
  apf::Field* scales;
  scales = apf::createLagrangeField(m, "proteus_size_scale", apf::VECTOR, 1);
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* v;
  while ((v = m->iterate(it))) {
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
  bool operator<(const SortingStruct& other) const
  {
    return wm < other.wm;
  }
};

static apf::Field* getSizeFrames(apf::Field* hessians, apf::Field* gradphi)
{
  apf::Mesh* m = apf::getMesh(gradphi);
  apf::Field* frames;
  frames = apf::createLagrangeField(m, "proteus_size_frame", apf::MATRIX, 1);
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* v;
  while ((v = m->iterate(it))) {
    apf::Vector3 gphi;
    apf::getVector(gradphi, v, 0, gphi);
    apf::Vector3 dir;
    if (gphi.getLength() > 1e-16)
      dir = gphi.normalize();
    else
      dir = apf::Vector3(1,0,0);
    apf::Matrix3x3 hessian;
    apf::getMatrix(hessians, v, 0, hessian);
    apf::Vector3 eigenVectors[3];
    double eigenValues[3];
    apf::eigen(hessian, eigenVectors, eigenValues);
    SortingStruct ssa[3];
    for (int i = 0; i < 3; ++i) {
      ssa[i].v = eigenVectors[i];
      ssa[i].wm = std::fabs(eigenValues[i]);
    }
    std::sort(ssa, ssa + 3);
    assert(ssa[2].wm >= ssa[1].wm);
    assert(ssa[1].wm >= ssa[0].wm);
    double firstEigenvalue = ssa[2].wm;
    apf::Matrix3x3 frame;
    frame[0] = dir;
    if (firstEigenvalue > 1e-16) {
      apf::Vector3 firstEigenvector = ssa[2].v;
      frame[1] = apf::reject(firstEigenvector, dir);
      frame[2] = apf::cross(frame[0], frame[1]);
      if (frame[2].getLength() < 1e-16)
        frame = apf::getFrame(dir);
    } else
      frame = apf::getFrame(dir);
    for (int i = 0; i < 3; ++i)
      frame[i] = frame[i].normalize();
    frame = apf::transpose(frame);
    apf::setMatrix(frames, v, 0, frame);
  }
  m->end(it);
  return frames;
}

static apf::Field* getERMSizeFrames(apf::Field* hessians, apf::Field* gradphi,apf::Field* frame_comps[3],std::string adapt_type)
{
  apf::Mesh* m = apf::getMesh(gradphi);
  apf::Field* frames;
  frames = apf::createLagrangeField(m, "proteus_size_frame", apf::MATRIX, 1);
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* v;
  while ((v = m->iterate(it))) {
    apf::Matrix3x3 frame(1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0);
    if(adapt_type=="isotropic")
    {
    }
    else
    {
    apf::Vector3 gphi;
    apf::getVector(gradphi, v, 0, gphi);
    apf::Vector3 dir;
    if (gphi.getLength() > 1e-16)
      dir = gphi.normalize();
    else
      dir = apf::Vector3(1,0,0);
    apf::Matrix3x3 hessian;
    apf::getMatrix(hessians, v, 0, hessian);
    apf::Vector3 eigenVectors[3];
    double eigenValues[3];
    apf::eigen(hessian, eigenVectors, eigenValues);
    SortingStruct ssa[3];
    for (int i = 0; i < 3; ++i) {
      ssa[i].v = eigenVectors[i];
      ssa[i].wm = std::fabs(eigenValues[i]);
    }
    std::sort(ssa, ssa + 3);
    assert(ssa[2].wm >= ssa[1].wm);
    assert(ssa[1].wm >= ssa[0].wm);
    double firstEigenvalue = ssa[2].wm;
    frame[0] = dir;
    if (firstEigenvalue > 1e-16) {
      frame[0] = ssa[2].v;
      frame[1] = ssa[1].v;
      frame[2] = ssa[0].v;
    } else
      frame = apf::getFrame(dir);
    for (int i = 0; i < 3; ++i)
      frame[i] = frame[i].normalize();
    frame = apf::transpose(frame);
    }
    apf::setMatrix(frames, v, 0, frame);
    apf::setVector(frame_comps[0], v, 0, frame[0]);
    apf::setVector(frame_comps[1], v, 0, frame[1]);
    apf::setVector(frame_comps[2], v, 0, frame[2]);
  }
  m->end(it);
  return frames;
}

int MeshAdaptPUMIDrvr::calculateAnisoSizeField()
{
  apf::Field* phif = m->findField("phi");
  assert(phif);
  apf::Field* gradphi = apf::recoverGradientByVolume(phif);
  apf::Field* grad2phi = apf::recoverGradientByVolume(gradphi);
  apf::Field* hess = computeHessianField(grad2phi);
  apf::destroyField(grad2phi);
  apf::Field* curves = getCurves(hess, gradphi);
  freeField(size_scale);
  
  size_scale = getSizeScales(phif, curves, hmin, hmax, nAdapt);
  apf::destroyField(curves);
  freeField(size_frame);
  size_frame = getSizeFrames(hess, gradphi);
  apf::destroyField(hess);

  apf::destroyField(gradphi);
  for (int i = 0; i < 2; ++i)
    SmoothField(size_scale);
 
  return 0;
}

struct Smoother : public apf::CavityOp
{
  Smoother(apf::Field* f):
    apf::CavityOp(apf::getMesh(f))
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
  virtual Outcome setEntity(apf::MeshEntity* e)
  {
    if (apf::hasEntity(newField, e))
      return SKIP;
    if ( ! this->requestLocality(&e, 1))
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
    for (int i = 0; i < edges.n; ++i) {
      apf::MeshEntity* ov = apf::getEdgeVertOppositeVert(
          mesh, edges.e[i], vertex);
      apf::getComponents(field, ov, 0, &value[0]);
      sum += value;
    }
    sum /= edges.n + 1;
    apf::setComponents(newField, vertex, 0, &sum[0]);
    ++nApplied;
  }
  apf::Field* field;
  apf::Field* newField;
  apf::MeshEntity* vertex;
  apf::DynamicVector sum;
  apf::DynamicVector value;
  apf::MeshTag* emptyTag;
  int nApplied;
};

static void SmoothField(apf::Field* f)
{
  Smoother op(f);
  op.applyToDimension(0);
}

int MeshAdaptPUMIDrvr::getERMSizeField(double err_total)
{
  double eps_u = 0.002; //distance from the interface
  double tolerance = 0.1;
  double alpha;
  if(target_error!=0){
    alpha = target_error/err_total;
//    if(alpha>1)
//      alpha=1.0;
  }
  else 
    alpha = 0.125;//tolerance/rel_err_total; //refinement constant

  freeField(size_frame);
  freeField(size_scale);
  freeField(size_iso);

//set the desired size field over regions
  apf::Mesh* m = apf::getMesh(err_reg);
  size_scale = apf::createLagrangeField(m, "proteus_size_scale", apf::VECTOR, 1);
  apf::MeshIterator* it;
  apf::MeshEntity* v;
  int numel = 0;
  int nsd = m->getDimension();
  it = m->begin(nsd);
  apf::Field* size_iso_reg = apf::createField(m, "iso_size",apf::SCALAR,apf::getConstant(nsd));
  size_iso = apf::createLagrangeField(m, "proteus_size",apf::SCALAR,1);

  //Get total number of elements
  numel = m->count(nsd);
  PCU_Add_Ints(&numel, 1);

  double err_dest = alpha*err_total/sqrt(numel);
if(comm_rank==0) std::cout<<"refinement ratio "<<alpha<<" error destination "<<err_dest<<" numel "<<numel<<std::endl;
  double err_curr = 0.0;
  apf::Vector3 err_vect;
  //compute the new size field
  apf::MeshElement* element;
  apf::MeshEntity* reg;
  it = m->begin(nsd); 
  while(reg=m->iterate(it)){
    double h_old;
    double h_new;
    element = apf::createMeshElement(m,reg);
    //h_old = pow(apf::measure(element),1.0/3.0);
    h_old = pow(apf::measure(element)*6*sqrt(2),1.0/3.0); //edge of a regular tet
    apf::getVector(err_reg,reg,0,err_vect);
    err_curr = err_vect[0];
    if(alpha==1)
      h_new = h_old;
    else
      h_new = h_old*(err_dest/err_curr);
    apf::setScalar(size_iso_reg,reg,0,h_new);
  }
  apf::destroyMeshElement(element);
  m->end(it);
/*
  it = m->begin(0);
  apf::Adjacent regions;
  apf::Adjacent edges;
  while((v=m->iterate(it))){
    double h_old, h_new;
    m->getAdjacent(v,nsd,regions);
    m->getAdjacent(v,1,edges);
    err_curr= 0;
    for(int i=0;i<regions.getSize();i++){
      apf::getVector(err_reg,regions[i],0,err_vect);
      err_curr += err_vect[0];
    }
    err_curr /= regions.getSize();

    h_old = 0;
    for(int i=0; i<edges.getSize();i++){
      element = apf::createMeshElement(m,edges[i]);
      h_old+=apf::measure(element); 
      apf::destroyMeshElement(element);
    }
    h_old /= edges.getSize();
    h_new = h_old*sqrt(err_dest/err_curr);
    apf::setScalar(size_iso,v,0,h_new);
  }
  m->end(it);
*/
//*
  it = m->begin(0);
  while((v=m->iterate(it))){
    //minToEntity(size_iso_reg, size_iso, v);
    //volumeAverageToEntity(size_iso_reg, size_iso, v);
    averageToEntity(size_iso_reg, size_iso, v);
  }
  m->end(it); 
//*/

//Get the anisotropic size frame
  apf::Field* phif = m->findField("phi");
  apf::Field* gradphi = apf::recoverGradientByVolume(phif);
  apf::Field* grad2phi = apf::recoverGradientByVolume(gradphi);
  apf::Field* speedF = extractSpeed(m->findField("velocity"));
  apf::Field* gradSpeed = apf::recoverGradientByVolume(speedF);
  apf::Field* grad2Speed = apf::recoverGradientByVolume(gradSpeed);
  apf::Field* hess = computeHessianField(grad2phi);
  apf::Field* curves = getCurves(hess, gradphi);
  //apf::Field* metricf = computeMetricField(gradphi,grad2phi,size_iso,eps_u);
  apf::Field* metricf = computeMetricField(gradSpeed,grad2Speed,size_iso,eps_u);

  apf::Field* frame_comps[3] = {apf::createLagrangeField(m, "frame_0", apf::VECTOR, 1),apf::createLagrangeField(m, "frame_1", apf::VECTOR, 1),apf::createLagrangeField(m, "frame_2", apf::VECTOR, 1)};
  //size_frame = getERMSizeFrames(hess, gradphi,frame_comps);
  //size_frame = getERMSizeFrames(metricf, gradphi,frame_comps);
  size_frame = getERMSizeFrames(grad2Speed,gradSpeed,frame_comps,adapt_type_config);
//

//Clipped Field

  apf::Field* clipped_vtx = apf::createLagrangeField(m, "iso_clipped",apf::SCALAR,1);

//Set the size scale for vertices
  it = m->begin(0);
  apf::Vector3 scale;
  while ((v = m->iterate(it))) {
    double vtx_vol=0;
    double phi = apf::getScalar(phif, v, 0);
    apf::Vector3 curve;
    apf::getVector(curves, v, 0, curve);

    //apf::Matrix3x3 hessian;
    //apf::getMatrix(hess, v, 0, hessian);
    apf::Matrix3x3 metric;
    apf::getMatrix(metricf, v, 0, metric);
    apf::Vector3 eigenVectors[3];
    double eigenValues[3];
    //apf::eigen(hessian, eigenVectors, eigenValues);
    apf::eigen(metric, eigenVectors, eigenValues);
    SortingStruct ssa[3];
    for (int i = 0; i < 3; ++i) {
      ssa[i].v = eigenVectors[i];
      ssa[i].wm = std::fabs(eigenValues[i]);
    }
    std::sort(ssa, ssa + 3);
/*
    assert(ssa[2].wm >= ssa[1].wm);
    assert(ssa[1].wm >= ssa[0].wm);
*/
    double lambda[3] = {ssa[2].wm, ssa[1].wm, ssa[0].wm};

    if(apf::getScalar(size_iso,v,0) < hmin)
      apf::setScalar(clipped_vtx,v,0,-1);
    else if(apf::getScalar(size_iso,v,0) > hmax)
      apf::setScalar(clipped_vtx,v,0,1);
    else
      apf::setScalar(clipped_vtx,v,0,0);

    //if(comm_rank==0)std::cout<<"Original Lambdas "<<lambda[0]<<" "<<lambda[1]<<" "<<lambda[2]<<std::endl;
    scaleFormulaERM(phi,hmin,hmax,apf::getScalar(size_iso,v,0),curve,lambda,eps_u,scale,adapt_type_config);
    //if(comm_rank==0) std::cout<<"Scales "<<scale[0]<<" "<<scale[1]<<" "<<scale[2]<<"Lambdas "<<lambda[0]<<" "<<lambda[1]<<" "<<lambda[2]<<std::endl;
    apf::setVector(size_scale,v,0,scale);
  }
  m->end(it);
/*
  SmoothField(size_scale);
  SmoothField(size_scale);
  SmoothField(size_scale);
*/

  if(logging_config=="on"){
    char namebuffer[20];
    sprintf(namebuffer,"pumi_preadapt_%i",nAdapt);
    apf::writeVtkFiles(namebuffer, m);
  }
  apf::destroyField(size_iso_reg); //will throw error if not destroyed
  apf::destroyField(clipped_vtx);
  apf::destroyField(grad2phi);
  //apf::destroyField(phif);
  apf::destroyField(curves);
  apf::destroyField(hess);
  apf::destroyField(metricf);
  apf::destroyField(gradphi);
  apf::destroyField(frame_comps[0]); apf::destroyField(frame_comps[1]); apf::destroyField(frame_comps[2]);
  apf::destroyField(speedF);
  apf::destroyField(gradSpeed);
  apf::destroyField(grad2Speed);

  freeField(size_iso); //no longer necessary
  return 0;
}

int MeshAdaptPUMIDrvr::testIsotropicSizeField()
{
    size_scale = apf::createLagrangeField(m, "proteus_size",apf::VECTOR,1);
    size_frame = apf::createLagrangeField(m, "proteus_size_frame", apf::MATRIX, 1);
    apf::MeshIterator* it = m->begin(0);
    apf::MeshEntity* v;
    //apf::Field* phif = m->findField("phi");
    while(v = m->iterate(it)){
      double phi = hmin;//apf::getScalar(phif, v, 0);
      clamp(phi,hmin,hmax);
      apf::Vector3 scale = (apf::Vector3(1.0,1.0,1.0))*phi;
      apf::setVector(size_scale, v, 0, scale);
      apf::Matrix3x3 frame(1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0);
      apf::setMatrix(size_frame, v, 0, frame);
    }
    for(int i=0; i < 3; i++)
      SmoothField(size_scale);
    char namebuffer[20];
    sprintf(namebuffer,"pumi_adapt_%i",nAdapt);
    apf::writeVtkFiles(namebuffer, m);
}




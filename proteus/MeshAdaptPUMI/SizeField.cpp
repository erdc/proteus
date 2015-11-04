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

int MeshAdaptPUMIDrvr::getERMSizeField(double err_total)
{
  double alpha = 0.6; //refinement constant

  free(size_frame);
  free(size_scale);
  free(size_iso);

  apf::Mesh* m = apf::getMesh(err_reg);
  size_frame = apf::createLagrangeField(m, "proteus_size_frame", apf::MATRIX, 1);
  size_scale = apf::createLagrangeField(m, "proteus_size_scale", apf::VECTOR, 1);
  apf::MeshIterator* it;
  int numel = 0;
  int nsd = m->getDimension();
  it = m->begin(nsd);
  apf::Field* size_iso_reg = apf::createField(m, "iso_size",apf::SCALAR,apf::getConstant(nsd));
  size_iso = apf::createLagrangeField(m, "proteus_size",apf::SCALAR,1);

  apf::MeshEntity* reg;
  while( reg=m->iterate(it)){
    numel++; 
  }
  m->end(it);
  it = m->begin(nsd); 
  double err_dest = alpha*err_total/sqrt(numel);
  double err_curr = 0.0;
  //compute the new size field
  double h_old = 0.0;
  double h_new;
  apf::MeshElement* element;
///*
  while(reg=m->iterate(it)){
    element = apf::createMeshElement(m,reg);
    h_old = pow(apf::measure(element),1.0/3.0);
    err_curr = apf::getScalar(err_reg,reg,0);
std::cout<<"Err curr "<<err_curr<<" Err_dest "<<err_dest<<std::endl;
    h_new = h_old*pow(err_dest/err_curr,0.5);
    if(h_new>hmax) h_new = hmax;
    if(h_new<hmin) h_new = hmin;
    std::cout<<"Old length "<<h_old<<" New length "<<h_new<<std::endl;
    apf::setScalar(size_iso_reg,reg,0,h_new);
  }
  apf::destroyMeshElement(element);
  m->end(it);
//*/
  apf::MeshEntity* v;
  it = m->begin(0);
  apf::Vector3 scale(0.1,0.1,0.1);
  apf::Matrix3x3 identity(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0); //isotropic size frame
  while ((v = m->iterate(it))) {
    averageToEntity(size_iso_reg, size_iso, v);
    scale = apf::Vector3(1.0,1.0,1.0)* apf::getScalar(size_iso,v,0); //isotropic
    //scale = apf::Vector3(1.0,1.0,1.0)*(hmin+hmax)/2.0; //for debugging purposes only
    apf::setVector(size_scale,v,0,scale);
    apf::setMatrix(size_frame,v,0,identity);
  }
  m->end(it);
  //for (int i = 0; i < 2; ++i)
  SmoothField(size_scale);
  apf::destroyField(size_iso_reg); //will throw error if not destroyed
  apf::destroyField(err_reg);
  freeField(size_iso); //no longer necessary
  return 0;
}

static apf::Matrix3x3 hessianFormula(apf::Matrix3x3 const& g2phi)
{
  apf::Matrix3x3 g2phit = apf::transpose(g2phi);
  return (g2phi + g2phit) / 2;
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

  double Km = 0.5* (a - 2 * b) / pow(g * g, 1.5);

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
  apf::destroyField(phif);
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

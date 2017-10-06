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

static void SmoothField(apf::Field *f);

/* Based on the distance from the interface epsilon can be controlled to determine
   thickness of refinement near the interface */
static double isotropicFormula(double phi, double dphi, double verr, double hmin, double hmax, double phi_s = 0)
{
  double size;
  double dphi_size_factor;
  double v_size_factor;
  //This is just a hack for now. This disable the refinement over phi and does it over phi_s
  // if (phi_s != 0.0)
  // {
  if (fabs(phi_s) < 5.0 * hmin)
    return hmin;
  else
    return hmax;
  // }
  // else
  // {
  //   if (fabs(phi) < 5.0 * hmin)
  //   {
  //     dphi_size_factor = fmax(hmin / 10.0, fmin(1.0, pow(((hmin / 1000.0) / fabs(dphi + 1.0e-8)), 1.0 / 2.0)));
  //     size = hmin * dphi_size_factor;
  //   }
  //   else
  //     size = hmax;

  //   size = fmax(hmin / 100.0, fmin(size, 0.001 / (verr + 1.0e-8)));

  //   return size;
  // }
}

int MeshAdaptPUMIDrvr::calculateSizeField()
{
  freeField(size_iso);
  size_iso = apf::createLagrangeField(m, "proteus_size", apf::SCALAR, 1);
  apf::MeshIterator *it = m->begin(0);
  apf::MeshEntity *v;
  apf::Field *phif = m->findField("phi");
  assert(phif);
  ////////////////////////////////////////
  apf::Field *phisError = m->findField("phi_s");
  assert(phisError);
  /////////////////////////////////////////
  apf::Field *phiCorr = m->findField("phiCorr");
  assert(phiCorr);
  apf::Field *velocityError = m->findField("velocityError");
  assert(phiCorr);
  while ((v = m->iterate(it)))
  {
    double phi = apf::getScalar(phif, v, 0);
    double phi_s = apf::getScalar(phisError, v, 0);
    // double dphi = apf::getScalar(phiCorr, v, 0);
    // double verr = apf::getScalar(velocityError, v, 0);
    double size = isotropicFormula(0.0, 0.0, 0.0, hmin, hmax, phi_s);
    apf::setScalar(size_iso, v, 0, size);
  }
  m->end(it);
  /*
    If you just smooth then hmax will just diffuse into the hmin band
    and you won't really get a band around phi=0 with uniform diameter
    hmin. Instead, reset to hmin after each smooth within the band in
    order to ensure the band uses hmin. Iterate on that process until
    changes in the smoothed size are less than 50% of hmin.
   */
/*
  double err_h_max=hmax;
  //int its=0;
  //while (err_h_max > 0.5*hmin && its < 200)
  for(int its=0;its<200;its++)
    {
      its++;
      SmoothField(size_iso);
      err_h_max=0.0;
      it = m->begin(0);
      while ((v = m->iterate(it))) {
	double phi = apf::getScalar(phif, v, 0);
	double dphi = apf::getScalar(phiCorr, v, 0);
	double verr = apf::getScalar(velocityError, v, 0);
	double size_current = apf::getScalar(size_iso, v, 0);
	double size = fmin(size_current,isotropicFormula(phi, dphi, verr, hmin, hmax));
	err_h_max = fmax(err_h_max,fabs(size_current-size));
	apf::setScalar(size_iso, v, 0, size);
      }
      m->end(it);
    }
  PCU_Barrier();    
*/
  gradeMesh();
  return 0;
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

/*
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
*/

void MeshAdaptPUMIDrvr::volumeAverageToEntity(apf::Field *ef, apf::Field *vf,
                                              apf::MeshEntity *ent)
//Serves the same purpose as averageToEntity but considers a volume-weighted average
{
  apf::Mesh *m = apf::getMesh(ef);
  apf::Adjacent elements;
  apf::MeshElement *testElement;
  m->getAdjacent(ent, m->getDimension(), elements);
  double s = 0;
  double invVolumeTotal = 0;
  for (std::size_t i = 0; i < elements.getSize(); ++i)
  {
    testElement = apf::createMeshElement(m, elements[i]);
    s += apf::getScalar(ef, elements[i], 0) / apf::measure(testElement);
    invVolumeTotal += 1.0 / apf::measure(testElement);
    if (comm_rank == 0)
    {
      std::cout << "What is s " << s << " Volume? " << apf::measure(testElement) << " scale? " << apf::getScalar(ef, elements[i], 0) << std::endl;
    }
    apf::destroyMeshElement(testElement);
  }
  s /= invVolumeTotal;
  if (comm_rank == 0)
  {
    std::cout << "What is s final? " << s << std::endl;
  }
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

static apf::Field *computeMetricField(apf::Field *gradphi, apf::Field *grad2phi, apf::Field *size_iso, double eps_u)
//Function used to generate a field of metric tensors at vertices. It is meant to be an umbrella function that can compute Hessians too.
//Currently needs development.
{
  apf::Mesh *m = apf::getMesh(grad2phi);
  apf::Field *metricf = createLagrangeField(m, "proteus_metric", apf::MATRIX, 1);
  apf::MeshIterator *it = m->begin(0);
  apf::MeshEntity *v;
  while ((v = m->iterate(it)))
  {
    apf::Matrix3x3 g2phi;
    apf::getMatrix(grad2phi, v, 0, g2phi);
    apf::Vector3 gphi;
    apf::getVector(gradphi, v, 0, gphi);
    apf::Matrix3x3 gphigphit(gphi[0] * gphi[0], gphi[0] * gphi[1], gphi[0] * gphi[2],
                             gphi[0] * gphi[1], gphi[1] * gphi[1], gphi[1] * gphi[2],
                             gphi[0] * gphi[2], gphi[1] * gphi[2], gphi[2] * gphi[2]);
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

  for (int i = 0; i < 3; ++i)
    clamp(scale[i], hmin, hmax);
}

static void scaleFormulaERM(double phi, double hmin, double hmax, double h_dest,
                            apf::Vector3 const &curves,
                            double lambda[3], double eps_u, apf::Vector3 &scale)
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
  double epsilon = 7.0 * hmin;
  double lambdamin = 1.0 / (hmin * hmin);
  if (lambda[1] < 1e-10)
  {
    lambda[1] = lambdamin;
    lambda[2] = lambdamin;
  }
  if (lambda[2] < 1e-10)
  {
    lambda[2] = lambdamin;
  }
  ///* useful
  scale[0] = h_dest * pow((lambda[1] * lambda[2]) / (lambda[0] * lambda[0]), 1.0 / 6.0);
  scale[1] = sqrt(lambda[0] / lambda[1]) * scale[0];
  scale[2] = sqrt(lambda[0] / lambda[2]) * scale[0];
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
  frames = apf::createLagrangeField(m, "proteus_size_frame", apf::MATRIX, 1);
  apf::MeshIterator *it = m->begin(0);
  apf::MeshEntity *v;
  while ((v = m->iterate(it)))
  {
    apf::Matrix3x3 frame(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
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
    frame[0] = dir;
    if (firstEigenvalue > 1e-16)
    {
      frame[0] = ssa[2].v;
      frame[1] = ssa[1].v;
      frame[2] = ssa[0].v;
    }
    else
      frame = apf::getFrame(dir);
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
  for (int i = 0; i < 2; ++i)
    SmoothField(size_scale);

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

int MeshAdaptPUMIDrvr::getERMSizeField(double err_total)
//High level function that obtains the size scales and the size frames for ERM-based adapt and uses the computed total error
{
  freeField(size_frame);
  freeField(size_scale);
  freeField(size_iso);

  //Initialize fields and needed types/variables
  apf::Mesh *m = apf::getMesh(err_reg);
  apf::MeshIterator *it;
  apf::MeshEntity *v;
  apf::MeshElement *element;
  apf::MeshEntity *reg;
  size_iso = apf::createLagrangeField(m, "proteus_size", apf::SCALAR, 1);
  apf::Field *size_iso_reg = apf::createField(m, "iso_size", apf::SCALAR, apf::getConstant(nsd));
  apf::Field *clipped_vtx = apf::createLagrangeField(m, "iso_clipped", apf::SCALAR, 1);

  //Get total number of elements
  int numel = 0;
  int nsd = m->getDimension();
  numel = m->count(nsd);
  PCU_Add_Ints(&numel, 1);

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
      h_old = pow(apf::measure(element) * 6 * sqrt(2), 1.0 / 3.0); //edge of a regular tet
    apf::getVector(err_reg, reg, 0, err_vect);
    err_curr = err_vect[0];
    errRho_curr = apf::getScalar(errRho_reg, reg, 0);
    //h_new = h_old*errRho_target/errRho_curr;
    //h_new = h_old*sqrt(apf::measure(element))/sqrt(domainVolume)*target_error/err_curr;
    if (target_error == 0)
      target_error = err_total / sqrt(numel);
    h_new = h_old * (target_error / err_curr);
    apf::setScalar(size_iso_reg, reg, 0, h_new);
    apf::destroyMeshElement(element);
  }
  m->end(it);

  //Transfer size field from elements to vertices through averaging
  it = m->begin(0);
  while ((v = m->iterate(it)))
  {
    averageToEntity(size_iso_reg, size_iso, v);
  }
  m->end(it);

  //Get the anisotropic size frame
  if (adapt_type_config == "anisotropic")
  {
    size_scale = apf::createLagrangeField(m, "proteus_size_scale", apf::VECTOR, 1);
    double eps_u = 0.002; //distance from the interface
    apf::Field *phif = m->findField("phi");
    apf::Field *gradphi = apf::recoverGradientByVolume(phif);
    apf::Field *grad2phi = apf::recoverGradientByVolume(gradphi);
    apf::Field *speedF = extractSpeed(m->findField("velocity"));
    apf::Field *gradSpeed = apf::recoverGradientByVolume(speedF);
    apf::Field *grad2Speed = apf::recoverGradientByVolume(gradSpeed);
    apf::Field *hess = computeHessianField(grad2phi);
    apf::Field *curves = getCurves(hess, gradphi);
    //apf::Field* metricf = computeMetricField(gradphi,grad2phi,size_iso,eps_u);
    apf::Field *metricf = computeMetricField(gradSpeed, grad2Speed, size_iso, eps_u);
    apf::Field *frame_comps[3] = {apf::createLagrangeField(m, "frame_0", apf::VECTOR, 1), apf::createLagrangeField(m, "frame_1", apf::VECTOR, 1), apf::createLagrangeField(m, "frame_2", apf::VECTOR, 1)};
    size_frame = getERMSizeFrames(metricf, gradSpeed, frame_comps);

    //Set the size scale for vertices
    it = m->begin(0);
    apf::Vector3 scale;
    while ((v = m->iterate(it)))
    {
      double vtx_vol = 0;
      double phi = apf::getScalar(phif, v, 0);
      apf::Vector3 curve;
      apf::getVector(curves, v, 0, curve);
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

      if (apf::getScalar(size_iso, v, 0) < hmin)
        apf::setScalar(clipped_vtx, v, 0, -1);
      else if (apf::getScalar(size_iso, v, 0) > hmax)
        apf::setScalar(clipped_vtx, v, 0, 1);
      else
        apf::setScalar(clipped_vtx, v, 0, 0);

      scaleFormulaERM(phi, hmin, hmax, apf::getScalar(size_iso, v, 0), curve, lambda, eps_u, scale);
      apf::setVector(size_scale, v, 0, scale);
    }
    m->end(it);
    apf::synchronize(size_scale);
    apf::destroyField(grad2phi);
    apf::destroyField(curves);
    apf::destroyField(hess);
    apf::destroyField(metricf);
    apf::destroyField(gradphi);
    apf::destroyField(frame_comps[0]);
    apf::destroyField(frame_comps[1]);
    apf::destroyField(frame_comps[2]);
    apf::destroyField(speedF);
    apf::destroyField(gradSpeed);
    apf::destroyField(grad2Speed);
  }
  else
  {
    it = m->begin(0);
    while ((v = m->iterate(it)))
    {
      double tempScale = apf::getScalar(size_iso, v, 0);
      if (tempScale < hmin)
        apf::setScalar(clipped_vtx, v, 0, -1);
      else if (tempScale > hmax)
        apf::setScalar(clipped_vtx, v, 0, 1);
      else
        apf::setScalar(clipped_vtx, v, 0, 0);
      clamp(tempScale, hmin, hmax);
      apf::setScalar(size_iso, v, 0, tempScale);
    }
    apf::synchronize(size_iso);
    m->end(it);
    if (target_element_count != 0)
    {
      sam::scaleIsoSizeField(size_iso, target_element_count);
      clampField(size_iso, hmin, hmax);
      SmoothField(size_iso);
    }
  }

  //Destroy locally required fields
  apf::destroyField(size_iso_reg);
  apf::destroyField(clipped_vtx);

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
}

int MeshAdaptPUMIDrvr::gradeMesh()
//Function to grade isotropic mesh through comparison of edge vertex size ratios
//For simplicity, we do not bother with accounting for entities across partitions
{
  //
  if(comm_rank==0)
    std::cout<<"Starting grading\n";
  apf::MeshIterator* it = m->begin(1);
  apf::MeshEntity* edge;
  apf::Adjacent edgAdjVert;
  apf::Adjacent vertAdjEdg;
  double gradingFactor = 1.5;
  double size[2];
  std::queue<apf::MeshEntity*> markedEdges;
  apf::MeshTag* isMarked = m->createIntTag("isMarked",1);

  //marker structure for 0) not marked 1) marked 2)storage
  int marker[3] = {0,1,0}; 

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
  while(!markedEdges.empty()){
    edge = markedEdges.front();
    m->getAdjacent(edge, 0, edgAdjVert);
    for (std::size_t i=0; i < edgAdjVert.getSize(); ++i){
      size[i] = apf::getScalar(size_iso,edgAdjVert[i],0);
    }
    if(size[0]>gradingFactor*size[1]){
      size[0] = gradingFactor*size[1];
      apf::setScalar(size_iso,edgAdjVert[0],0,size[0]);
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
      apf::setScalar(size_iso,edgAdjVert[1],0,size[1]);
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
    m->setIntTag(edge,isMarked,&marker[0]);
    markedEdges.pop();
  }

  it = m->begin(1);
  while((edge=m->iterate(it))){
    m->removeTag(edge,isMarked);
  }
  m->end(it); 
  m->destroyTag(isMarked);
  apf::synchronize(size_iso);
  if(comm_rank==0)
    std::cout<<"Completed grading\n";
}


#include "MeshAdaptPUMI.h"
#include "apf.h"
#include "apfVector.h"
#include "apfPUMI.h"
#include "apfSPR.h"
#include "apfMesh.h"

enum {
  PHI_IDX = 5;
};

static void SmoothField(apf::Field* f);

/* Based on the distance from the interface
   epsilon can be controlled to determine
   thickness of refinement near the interface */
static double isotropicFormula(double* solution, double hmin, double hmax)
{
  static double const epsilon = 0.02;
  double phi = sqrt(solution[PHI_IDX] * solution[PHI_IDX]);
  double size;
  if (fabs(phi) < epsilon)
    size = hmin;
  else if (phi < 3 * epsilon)
    size = (hmin + hmax) / 2;
  else
    size = hmax;
  return size;
}

int MeshAdaptPUMIDrvr::CalculateSizeField()
{
  size_iso = apf::createLagrangeField(mesh, "proteus_size",apf::SCALAR,1);
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* v;
  apf::NewArray<double> sol(apf::countComponents(solution));
  while ((v = m->iterate(it))) {
    apf::getComponents(solution, v, 0, &sol[0]);
    double size = isotropicFormula(&sol[0], hmin, hmax);
    apf::setScalar(size_iso, v, 0, size);
  }
  m->end(it);
  for(int i=0; i < 3; i++)
    SmoothField(size_iso);
  apf::writeVtkFIles("pumi_size", mesh);
  return 0;
}

static apf::Field* extractPhi(apf::Field* solution)
{
  apf::Mesh* m = apf::getMesh(solution);
  apf::Field* phif = apf::createLagrangeField(mesh,"proteus_phi",apf::SCALAR,1);
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* v;
  apf::NewArray<double> tmp(apf::countComponents(solution));
  while ((v = m->iterate(it))) {
    apf::getComponents(solution, v, 0, &tmp[0]);
    double phi = tmp[PHI_IDX];
    apf::setScalar(phif, v, 0, phi);
  }
  m->end(it);
  return phif;
}

static apf::Matrix3x3 hessianFormula(apf::Matrix3x3 const& g2phi)
{
  apf::Matrix3x3 g2phit = apf::transpose(g2phi);
  return (g2phi + g2phit) / 2;
}

static apf::Field* computeHessianField(apf::Field* grad2phi)
{
  apf::Mesh* m = apf::getMesh(grad2phi);
  apf::Field* hessf = createLagrangeField(apf_mesh,"proteus_hess",apf::MATRIX,1);
  apf::MeshIterator* it = apf_mesh->begin(0);
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

/* curves based on hessian and gradient of phi.
   vector and matrix operators should be able
   to simplify this from its current component-based form */
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

  curve[0] = Km;
  curve[1] = Km + sqrt(Km * Km - Kg);
  curve[2] = Km - sqrt(Km * Km - Kg);
}

static apf::Field* getCurves(apf::Field* hessians, apf::Field* gradphi)
{
  apf::Field* curves;
  curves = apf::createLagrangeField(mesh, "proteus_curves", apf::VECTOR, 1);
  apf::Mesh* m = apf::getMesh(hessians);
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* v;
  while ((e = m->iterate(it))) {
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
  double epsilon = 2 * hmin + exp( -adapt_step / 4) * hmin * 2.5;
  if (fabs(phi) < 3 * epsilon) {
    scale[0] = hmin;
    scale[1] = sqrt(.0004 / fabs(curves[1]));
    scale[2] = sqrt(.0004 / fabs(curves[2]));
  } else {
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
  while ((e = m->iterate(it))) {
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

static apf::Field* getSizeFrames(apf::Field* gradphi)
{
  apf::Mesh* m = apf::getMesh(gradphi);
  apf::Field* frames;
  frames = apf::createLagrangeField(m, "proteus_size_frame", apf::MATRIX, 1);
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* v;
  while ((v = m->iterate(it))) {
    apf::Vector3 gphi;
    apf::getVector(gradphi, v, 0, gphi);
    apf::Matrix3x3 frame;
/* the anisotropic frame is just aligned
   to have grad(phi) along the first axis */
    frame = apf::transpose(apf::getFrame(gphi));
    apf::setMatrix(frames, v, 0, frame);
  }
  m->end(it);
  return scales;
}

/*
 Function to calculate anisotropic size field
 The idea is to set hmin in level set gradient's direction, the other 2 sizes are
 set based on the local curvatures (calculated from Hessian)
 Currently directions are chosen perpendicular to the gradient, since for the dambreak
 y axis is the thickness, it is a natural option as one of the directions, but will need to
 have a more general formulation
 */
int MeshAdaptPUMIDrvr::CalculateAnisoSizeField(pMAdapt MA_Drvr, apf::Field* f)
{
  apf::Field* phif = extractPhi(solution);
  apf::Field* gradphi = apf::recoverGradientByVolume(phif);
  apf::Field* grad2phi = apf::recoverGradientByVolume(gradphi);
  apf::Field* hess = computeHessianField(grad2phi);
  apf::Field* curves = getCurves(hess, gradphi);
  size_scale = getSizeScales(phif, curves, hmin, hmax, nAdapt);
  size_frame = getSizeFrames(gradphi);
  for (int i = 0; i < 2; ++i)
    SmoothField(size_scale);
}

static void getSelfAndNeighbors(apf::mesh* m, apf::meshEntity* v, apf::Up& vs)
{
  apf::Up es;
  m->getUp(v, es);
  for (int i = 0; i < es.n; ++i)
    vs.e[i] = apf::getEdgeVertOppositeVert(m, es[i], v);
  vs.e[vs.n] = v;
  ++vs.n;
}

/* the smoothing function used here is the average of the
   vertex value and neighboring vertex values, with the
   center vertex weighted equally as the neighbors */
static void SmoothField(apf::Field* f)
{
  apf::Mesh* m = apf::getMesh(f);
  int nc = apf::countComponents(f);
  apf::Field* sumf = apf::createPackedField(m, "proteus_sum", nc);
  apf::Field* numf = apf::createLagrangeField(m,"proteus_num",apf::SCALAR,1);
  apf::DynamicVector sum(nc);
  apf::DynamicVector val(nc);
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* v;
  while ((v = m->iterate(it))) {
    apf::Up vs;
    getSelfAndNeighbors(m, v, vs);
    for (int i = 0; i < nc; ++i)
      sum[i] = 0;
    int num = 0;
    for (int i = 0; i < vs.n; ++i)
      if (m->isOwned(vs.e[i])) {
        apf::getComponents(f, v, 0, &val[0]);
        sum += val;
        ++num;
      }
    apf::setComponents(sumf, v, 0, &sum[0]);
    apf::setScalar(numf, v, 0, num);
  }
  m->end(it);
  apf::accumulate(sumf);
  apf::accumulate(numf);
  apf::MeshIterator* it = m->begin(0);
  while ((v = m->iterate(it))) {
    apf::getComponents(sumf, v, 0, &sum[0]);
    num = apf::getScalar(numf, v, 0);
    sum /= num;
    apf::setComponents(f, v, 0, &sum[0]);
  }
  m->end(it);
  apf::destroyField(sumf);
  apf::destroyField(numf);
}

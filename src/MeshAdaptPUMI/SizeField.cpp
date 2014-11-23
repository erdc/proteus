#include "MeshAdaptPUMI.h"
#include <apf.h>
#include <apfVector.h>
#include <apfMesh.h>
#include <apfDynamicVector.h>
#include <string>
#include <iostream>
#include <sstream>

enum {
  PHI_IDX = 5
};

static void SmoothField(apf::Field* f);

/* Based on the distance from the interface epsilon can be controlled to determine
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
  freeField(size_iso);
  size_iso = apf::createLagrangeField(m, "proteus_size",apf::SCALAR,1);
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
  apf::writeVtkFiles("pumi_size", m);
  return 0;
}

static apf::Field* extractPhi(apf::Field* solution)
{
  apf::Mesh* m = apf::getMesh(solution);
  apf::Field* phif = apf::createLagrangeField(m,"proteus_phi",apf::SCALAR,1);
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
  int i;
  double m;
  };

static apf::Field* getSizeFrames(apf::Field* hessians, apf::Field* gradphi)     
{
  apf::Mesh* m = apf::getMesh(gradphi);
  apf::Field* frames;
  frames = apf::createLagrangeField(m, "proteus_size_frame", apf::MATRIX, 1);
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* v;
  while ((v = m->iterate(it))) {

//  Finding the eigenvalues and eigenvectors of the Hessian of the distance field.
//  Largest eigenvalue and the corresponding eigenvector, represent the principal curvature and principal
//  direction; respectively.	
    apf::Matrix3x3 hessian;                                                     
    apf::getMatrix(hessians, v, 0, hessian);                                    
    apf::Vector3 eigenVector[3];         
    double eigenValue[3];                                                       
	int n= apf::eigen(hessian, eigenVector, eigenValue);                        
 
//  Sorting EigenValues by decreasing absolute value 
    SortingStruct s[3] =
    {{0,fabs(eigenValue[0])},{1,fabs(eigenValue[1])},{2,fabs(eigenValue[2])}};
    if (s[2].m > s[1].m)
       std::swap(s[1],s[2]);
    if (s[1].m > s[0].m)
       std::swap(s[0],s[1]);
    if (s[2].m > s[1].m)
       std::swap(s[1],s[2]);
    int largest = s[0].i;
    int medium  = s[1].i;
    int smallest= s[2].i;

// 1st Frame Size (normal):         First dferivative of the phi and not the Hessian (to be more accurate)	
// 2nd Frame Size (1st tangential): Hessian's eigenvector corresponding to Hessian's larget eigenvalue
// 3rd Frame Size (2nd tangenital): Cross product of the normal vector and the 1st tangential vector	
	apf::Vector3 gphi;
    apf::getVector(gradphi, v, 0, gphi);
	apf::Matrix3x3 frameHess;
    frameHess[0] = gphi;                            
	frameHess[1] = eigenVector[largest];            
    frameHess[2] = cross(frameHess[0],frameHess[1]);

	apf::Matrix3x3 frame;
    for (int i = 0; i < 3; ++i)       
       frame[i] = frameHess[i].normalize(); 
    frame = apf::transpose(frame);
    apf::setMatrix(frames, v, 0, frame);
  }
  m->end(it);
  return frames;
}

int MeshAdaptPUMIDrvr::CalculateAnisoSizeField()
{
  apf::Field* phif = extractPhi(solution);
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
  size_frame = getSizeFrames(hess,gradphi); 

  apf::destroyField(hess);                    
  apf::destroyField(gradphi);
  for (int i = 0; i < 2; ++i)
    SmoothField(size_scale);
 
  //apf::writeVtkFiles("pumi_size", m);
  std::stringstream s1;
  s1 << "SizeField_t" << nAdapt;
  std::string str1 = s1.str();
  apf::writeVtkFiles(str1.c_str(), m);

  std::stringstream s2;
  s2 << "SizeField_t" << nAdapt<<".smb";
  std::string str2 = s2.str();
  m->writeNative(str2.c_str());

  return 0;
}

static void getSelfAndNeighbors(apf::Mesh* m, apf::MeshEntity* v, apf::Up& vs)
{
  apf::Up es;
  m->getUp(v, es);
  vs.n = es.n;
  for (int i = 0; i < es.n; ++i)
    vs.e[i] = apf::getEdgeVertOppositeVert(m, es.e[i], v);
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
        apf::getComponents(f, vs.e[i], 0, &val[0]);
        sum += val;
        ++num;
      }
    apf::setComponents(sumf, v, 0, &sum[0]);
    apf::setScalar(numf, v, 0, num);
  }
  m->end(it);
  apf::accumulate(sumf);
  apf::accumulate(numf);
  it = m->begin(0);
  while ((v = m->iterate(it))) {
    apf::getComponents(sumf, v, 0, &sum[0]);
    double num = apf::getScalar(numf, v, 0);
    sum /= num;
    apf::setComponents(f, v, 0, &sum[0]);
  }
  m->end(it);
  apf::destroyField(sumf);
  apf::destroyField(numf);
}

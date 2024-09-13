#include <iostream>
#include <algorithm>
#include <vector>
#include <valarray>
#include <map>
#include "FMMandFSW.h"

double pos(double a)
{return std::max<double>(a,0.0);}
double neg(double a)
{return std::min<double>(a,0.0);}

int main()
{
  using namespace std;

  int nx = 11;
  double Lx = 1.0; double dx = Lx/(nx-1);
  Mesh mesh;
  std::cout<<"creating mesh nx= "<<nx<<" Lx= "<<Lx<<std::endl;
  initializeMesh(mesh);
  
  edgeMeshElements(nx,mesh);
  regularEdgeMeshNodes(nx,Lx,mesh);
  constructElementBoundaryElementsArray_edge(mesh);

  std::cout<<"creating FMMEikonalSolver1d ... "<<std::endl;
  FMMEikonalSolver1d fmm(&mesh);

  //try redistancing from 1d circle
  double xc = 0.5; double radius = 0.1;
  std::valarray<double> phi0(nx),T(nx);
  for (int nN = 0; nN < mesh.nNodes_global; nN++)
    {
      double x = mesh.nodeArray[nN*3+0]; double d2 = (x-xc)*(x-xc) - radius*radius;
      phi0[nN] = d2;
    }
  
  std::cout<<"phi0 :"<<std::endl;
  for (int nN = 0; nN < mesh.nNodes_global; nN++)
    {
      std::cout<<"x= "<<mesh.nodeArray[nN*3+0]<<" phi["<<nN<<"] = "<<phi0[nN]<<std::endl;
    }
  double zeroTol=1.0e-4; double trialTol = 1.0e-1; int verbose = 10;
  
//   std::cout<<"initializingKnownPointsUsingMagnitude  ..."<<std::endl;
//   fmm.initializeKnownPointsUsingMagnitude(&phi0[0],&T[0],zeroTol,verbose);
//   std::cout<<"Known= "<<std::endl;
//   std::copy(fmm.Known.begin(),fmm.Known.end(),std::ostream_iterator<int>(std::cout," , "));
//   std::cout<<std::endl;
//   std::cout<<"T= "<<std::endl;
//   std::copy(&T[0],&T[0]+nx,std::ostream_iterator<double>(std::cout," , "));
//   std::cout<<std::endl;

//   std::cout<<"initializingKnownPointsUsingFront  ..."<<std::endl;
//   fmm.initializeKnownPointsUsingFrontIntersection(&phi0[0],&T[0],zeroTol,verbose);
//   std::cout<<"Known= "<<std::endl;
//   std::copy(fmm.Known.begin(),fmm.Known.end(),std::ostream_iterator<int>(std::cout," , "));
//   std::cout<<std::endl;
//   std::cout<<"T= "<<std::endl;
//   std::copy(&T[0],&T[0]+nx,std::ostream_iterator<double>(std::cout," , "));
//   std::cout<<std::endl;

  //typedef double (*PF)(double);
  //PF pf = pos;
  std::valarray<double> phi0p = phi0.apply(pos);
  //pf    = neg; 
  std::valarray<double> phi0m = phi0.apply(neg);
  std::valarray<double> Tp(T);
  std::valarray<double> Tm(T);
  std::valarray<double> nodalSpeeds(1.0,nx);

  fmm.solve(&phi0p[0],&nodalSpeeds[0],&Tp[0],zeroTol,trialTol,verbose);
  std::cout<<"phi0p= "<<std::endl;
  std::copy(&phi0p[0],&phi0p[0]+nx,std::ostream_iterator<double>(std::cout," , "));
  std::cout<<std::endl;
  std::cout<<"Tp= "<<std::endl;
  std::copy(&Tp[0],&Tp[0]+nx,std::ostream_iterator<double>(std::cout," , "));
  std::cout<<std::endl;

  fmm.solve(&phi0m[0],&nodalSpeeds[0],&Tm[0],zeroTol,trialTol,verbose);
  std::cout<<"phi0m= "<<std::endl;
  std::copy(&phi0m[0],&phi0m[0]+nx,std::ostream_iterator<double>(std::cout," , "));
  std::cout<<std::endl;
  std::cout<<"Tm= "<<std::endl;
  std::copy(&Tm[0],&Tm[0]+nx,std::ostream_iterator<double>(std::cout," , "));
  std::cout<<std::endl;

  T = Tp - Tm;
  std::cout<<"T= "<<std::endl;
  std::copy(&T[0],&T[0]+nx,std::ostream_iterator<double>(std::cout," , "));
  std::cout<<std::endl;

  std::cout<<"cleaning up"<<std::endl;
  deleteMesh(mesh);
  return 0;
}

#ifndef FMM_AND_FSW_H
#define FMM_AND_FSW_H

#include <vector>
#include <valarray>
#include <cmath>
#include "mesh.h"


class UnstructuredLocalUpwindSolvers
{
  /***********************************************************************
     holding place for "upwind" Eikonal solvers used for local updates in
     FMM and FSW solvers. 1d is just basic first order approach. 2d and 3d
     are from Qian Zhan Siam paper 

     Need to think about using this to eliminate use of separate FMM 1d,2d,3d solvers
   ***********************************************************************/
 public:
  UnstructuredLocalUpwindSolvers();
  virtual ~UnstructuredLocalUpwindSolvers();

  double solve1d(double* x_N, double* x_N0, double T_N0, 
		 double speed, int verbose);
  double solve2d(const double* x_C, const double* x_A, const double * x_B,
		 double T_A,  double T_B, 
		 double speed, int verbose);
  double solve2din3d(const double* x_C, const double* x_A, const double * x_B,
		     double T_A,  double T_B, 
		     double speed, int verbose);
  double solve3d(const double * x_D,const double* x_1, const double* x_2, 
		 const double * x_3,
		 double T_1,  double T_2,  double T_3, 
		 double speed, int verbose,
		 double areaTol=1.0e-8);
protected:
  void cross3d(double x0, double x1, double x2, double y0, double y1, double y2,
	       double& xy0, double& xy1, double& xy2)
  {
    xy0 = x1*y2 - x2*y1; xy1 = x2*y0 - x0*y2; xy2 = x0*y1 - x1*y0;
  }
  double parArea3d(double x0, double x1, double x2, double y0, double y1, double y2)
  {
    double xy0,xy1,xy2,norm;
    cross3d(x0,x1,x2,y0,y1,y2,xy0,xy1,xy2);
    norm = sqrt(xy0*xy0 + xy1*xy1 + xy2*xy2);
    return norm;
  }
  //could have common place for constants and enums ...
  const double UNINITIALIZED;
};


class FMMEikonalSolverBase
{
public:
  enum STATUS {FAR = -1, TRIAL = 0, KNOWN= 1};
  enum INIT_TYPE {MAGNITUDE,FRONT_AND_MAGNITUDE};

  FMMEikonalSolverBase(Mesh* meshIn=0,int nSpace=1,
		       INIT_TYPE initIn = MAGNITUDE,
		       bool forcePositiveInitialValuesIn=true);
  virtual ~FMMEikonalSolverBase();
  //initialize T, set Status and Known data structures
  //leave additional layer of indirection in case want more sophisticated
  //initialization later
  virtual bool initializeKnownPoints(const double * phi0, double * T, 
				     double zeroTol = 1.0e-4, int verbose = 0);

  //initialize T, set Status and Known data structures
  //This version sets point to known iff its magnitude is less than zeroTol
  virtual bool initializeKnownPointsUsingMagnitude(const double * phi0, double * T, 
						   double zeroTol = 1.0e-4, int verbose = 0);
  //initialize T, set Status and Known data structures
  //This version sets point to known if they are <= 0.0 or the zero level set
  //intersects their element
  virtual bool initializeKnownPointsUsingFrontIntersection(const double * phi0, double * T, 
							   double zeroTol = 1.0e-4, int verbose = 0);
  virtual bool initializeKnownPointsAlaWeakDirBCs(const double* phi0, double * T,
						  double zeroTol = 0.0, int verbose=0);
  //basic FMM algorithm
  virtual bool solve(const double* phi0, const double * nodalSpeeds, double * T,
		     double zeroTol = 1.0e-4, double trialTol = 1.0e-1,
		     int initTypeFlag = -1,// -1 -->ignore, 0 --> magn. 1--> frontInt 
		     int verbose = 0);

  virtual bool localUpdate(int N, double & T_N, const double * T,
			   double speed = 1.0, int verbose = 0) = 0;
  //keep public for now
  Mesh* mesh;
  std::vector<STATUS> Status;
  std::vector<int> Known;
  const int nSpace;
  //
  UnstructuredLocalUpwindSolvers localSolvers;
  //had some problems with static const initialization
  const double UNINITIALIZED;
  INIT_TYPE initFlag;
  bool forcePositiveInitialValues;
};

class FMMEikonalSolver1d : public FMMEikonalSolverBase
{
 public:
  
  FMMEikonalSolver1d(Mesh* meshIn=0,
		     INIT_TYPE initIn = MAGNITUDE);
  virtual ~FMMEikonalSolver1d();

  //go through node star and try to update Travel time at node N using known points
  //takes min over node star
  virtual bool localUpdate(int N, double & T_N, const double * T,
			   double speed = 1.0, int verbose = 0);
  //-------------------- data members --------------------
};

class FMMEikonalSolver2d : public FMMEikonalSolverBase
{
 public:
  
  FMMEikonalSolver2d(Mesh* meshIn=0,
		     INIT_TYPE initIn = MAGNITUDE);
  virtual ~FMMEikonalSolver2d();
  //go through node star and try to update Travel time at node N using known points
  //takes min over node star
  virtual bool localUpdate(int N, double & T_N, const double * T,
			   double speed = 1.0, int verbose = 0);
  //-------------------- data members --------------------
};


class FMMEikonalSolver3d : public FMMEikonalSolverBase
{
 public:
  
  FMMEikonalSolver3d(Mesh* meshIn=0,
		     INIT_TYPE initIn = FRONT_AND_MAGNITUDE);
  virtual ~FMMEikonalSolver3d();
  //go through node star and try to update Travel time at node N using known points
  //takes min over node star
  virtual bool localUpdate(int N, double & T_N, const double * T,
			   double speed = 1.0, int verbose = 0);

  //-------------------- data members --------------------
 protected:

};


class FSWEikonalSolverBase: public FMMEikonalSolverBase
{
public:
  //construct basic FSW solver using refPoints if given
  FSWEikonalSolverBase(Mesh* meshIn=0, int nSpace=1,
		       double atolIn = 1.0e-8, double rtolIn=1.0e-8,
		       int maxItsIn = 1000,
		       INIT_TYPE initIn = MAGNITUDE,
		       int nRefPointsIn = 0,
		       const double* refPointsIn = 0);
  virtual ~FSWEikonalSolverBase();

  //initialize T, set Status and Known data structures
  //This version sets point to known iff its magnitude is less than zeroTol
  virtual bool initializeKnownPointsUsingMagnitude(const double * phi0, double * T, 
						   double zeroTol = 1.0e-4, int verbose = 0);
  //initialize T, set Status and Known data structures
  //This version sets point to known if they are <= 0.0 or the zero level set
  //intersects their element
  virtual bool initializeKnownPointsUsingFrontIntersection(const double * phi0, double * T, 
							   double zeroTol = 1.0e-4, int verbose = 0);
  virtual bool initializeKnownPointsAlaWeakDirBCs(const double* phi0, double * T,
						  double zeroTol = 0.0, int verbose=0);
  //build reference point orderings, if none given in constructor assume Cartesian
  //grid
  virtual bool buildOrderings();
  //basic FSW algorithm
  virtual bool solve(const double* phi0, const double * nodalSpeeds, double * T,
		     double zeroTol = 1.0e-4, double trialTol = 1.0e-1,
		     int initTypeFlag = -1,// -1 -->ignore, 0 --> magn. 1 --> frontInt 
		     int verbose = 0);

  //points always 3d, 1d array nRefPoints*3
  int nRefPoints;
  std::valarray<double> refPoints;
  //nRefPoint vectors of length  nNodes_global giving nodes in increasing distance
  //from each reference point
  std::vector<std::vector<int> > Order;
  double iterAtol,iterRtol;
  int maxIts;
  std::valarray<double> T0;//keep around for convergence test
};


class FSWEikonalSolver1d : public FSWEikonalSolverBase
{
 public:
  
  FSWEikonalSolver1d(Mesh* meshIn=0, 
		     double atolIn = 1.0e-8, double rtolIn=1.0e-8,
		     int maxItsIn = 1000,
		     INIT_TYPE initIn = MAGNITUDE,
		     int nRefPointsIn = 0,
		     const double* refPointsIn = 0);
  virtual ~FSWEikonalSolver1d();

  //go through node star and try to update Travel time at node N using known points
  //takes min over node star
  virtual bool localUpdate(int N, double & T_N, const double * T,
			   double speed = 1.0, int verbose = 0);
  //-------------------- data members --------------------
};

class FSWEikonalSolver2d : public FSWEikonalSolverBase
{
 public:
  
  FSWEikonalSolver2d(Mesh* meshIn=0, 
		     double atolIn = 1.0e-8, double rtolIn=1.0e-8,
		     int maxItsIn = 1000,
		     INIT_TYPE initIn = MAGNITUDE,
		     int nRefPointsIn = 0,
		     const double* refPointsIn = 0);
  virtual ~FSWEikonalSolver2d();

  //go through node star and try to update Travel time at node N using known points
  //takes min over node star
  virtual bool localUpdate(int N, double & T_N, const double * T,
			   double speed = 1.0, int verbose = 0);
  //-------------------- data members --------------------
};

class FSWEikonalSolver3d : public FSWEikonalSolverBase
{
 public:
  
  FSWEikonalSolver3d(Mesh* meshIn=0, 
		     double atolIn = 1.0e-8, double rtolIn=1.0e-8,
		     int maxItsIn = 1000,
		     INIT_TYPE initIn = MAGNITUDE,
		     int nRefPointsIn = 0,
		     const double* refPointsIn = 0);
  virtual ~FSWEikonalSolver3d();

  //go through node star and try to update Travel time at node N using known points
  //takes min over node star
  virtual bool localUpdate(int N, double & T_N, const double * T,
			   double speed = 1.0, int verbose = 0);
  //-------------------- data members --------------------
};


//----------------------------------------------------------------------
// crude techniques for adjusting or reconstructing values near
// front, would be nice to put in narrow band procedure as well
//----------------------------------------------------------------------
namespace NearFrontInitialization
{
  bool localPWLreconstruction(int nSpace, Mesh* mesh, const double* phi0, double * phi0R, 
			      double zeroTol, int verbose = 0);

  //not implemented needs to be corrected
  bool localDGPWLreconstruction(int nSpace, Mesh* mesh, const int* l2g,
				const double* phi0, double * phi0R,
				int * reconTag,
				double zeroTol, int verbose = 0);
  bool copyOverOutsideBand(int N, double tol,
			   const double * in,
			   double * out);
}


#endif

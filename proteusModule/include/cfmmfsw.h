#ifndef FMM_FSW_H
#define FMM_FSW_H
#include "Python.h"

#include "mesh.h"
extern "C" 
{
//could just have one wrapper
class FMMEikonalSolverBase;
class FSWEikonalSolverBase;

typedef struct
{
  PyObject_HEAD
  FMMEikonalSolverBase* pFMMsolver;
  int nSpace;
} FMMEikonalSolver_wrapper;

extern 
int FMMEikonalSolver_wrapper_init(Mesh* cmesh,
				  int nSpace,
				  FMMEikonalSolver_wrapper * wrap); 

extern 
int FMMEikonalSolver_wrapper_del(FMMEikonalSolver_wrapper * wrap); 

extern
int FMMEikonalSolver_wrapper_copy(FMMEikonalSolver_wrapper * from, 
				  FMMEikonalSolver_wrapper * to); 

extern
int FMMEikonalSolver_wrapper_solve(FMMEikonalSolver_wrapper * self, 
				   const double * phi0,
				   const double * nodalSpeeds,
				   double * T,
				   double zeroTol,
				   double trialTol,
				   int initFlag,
				   int verbose); 

typedef struct
{
  PyObject_HEAD
  FSWEikonalSolverBase* pFSWsolver;
  int nSpace;
} FSWEikonalSolver_wrapper;

extern 
int FSWEikonalSolver_wrapper_init(Mesh* cmesh,
				  int nSpace,
				  double atol, double rtol,
				  int maxIts,
				  int initFlag,
				  int nRefPoints,
				  const double* refPoints,
				  FSWEikonalSolver_wrapper * wrap); 

extern 
int FSWEikonalSolver_wrapper_del(FSWEikonalSolver_wrapper * wrap); 

extern
int FSWEikonalSolver_wrapper_copy(FSWEikonalSolver_wrapper * from, 
				  FSWEikonalSolver_wrapper * to); 

extern
int FSWEikonalSolver_wrapper_solve(FSWEikonalSolver_wrapper * self, 
				   const double * phi0,
				   const double * nodalSpeeds,
				   double * T,
				   double zeroTol,
				   double trialTol,
				   int initFlag,
				   int verbose); 

extern 
int localPWL_wrapper(int nSpace, Mesh* cmesh, 
		     const double * phi0,
		     double * phi0R,
		     double zeroTol,
		     int verbose);
}
#endif

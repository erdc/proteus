#include "cfmmfsw.h"
#include "FMMandFSW.h"
#include "mesh.h"
#include <iostream>
#include <cassert>

extern "C"
{

int FMMEikonalSolver_wrapper_init(Mesh* cmesh,
				  int nSpace,
				  FMMEikonalSolver_wrapper * wrap)
{
  assert(wrap);
  wrap->pFMMsolver = NULL;
  if (nSpace == 1)
    {
      wrap->pFMMsolver = new FMMEikonalSolver1d(cmesh);
      wrap->nSpace = nSpace;
    }
  else if (nSpace == 2)
    {
      wrap->pFMMsolver = new FMMEikonalSolver2d(cmesh);
      wrap->nSpace = nSpace;
    }
  else if (nSpace == 3)
    {
      wrap->pFMMsolver = new FMMEikonalSolver3d(cmesh);
      wrap->nSpace = nSpace;
    }
  else
    {
      std::cout<<"nSpace= "<<nSpace<<" not implemented" <<std::endl;
      std::exit(1);
      return true;
    }
  return false;
}

int FMMEikonalSolver_wrapper_del(FMMEikonalSolver_wrapper * wrap)
{
  assert(wrap);
  if (wrap->pFMMsolver)
    delete wrap->pFMMsolver;
  return false;
}

int FMMEikonalSolver_wrapper_copy(FMMEikonalSolver_wrapper * from, 
				  FMMEikonalSolver_wrapper * to)
{
  assert(from); assert(to);
  if (to->pFMMsolver)
    delete to->pFMMsolver;
  if (from->nSpace == 1)
    {
      to->pFMMsolver = new FMMEikonalSolver1d(from->pFMMsolver->mesh);
      to->nSpace = from->nSpace;
    }
  else if (from->nSpace == 2)
    {
      to->pFMMsolver = new FMMEikonalSolver2d(from->pFMMsolver->mesh);
      to->nSpace = from->nSpace;
    }
  else if (from->nSpace == 3)
    {
      to->pFMMsolver = new FMMEikonalSolver3d(from->pFMMsolver->mesh);
      to->nSpace = from->nSpace;
    }
   else
    {
      std::cout<<"nSpace= "<<from->nSpace<<" not implemented" <<std::endl;
      std::exit(1);
      return true;
    }
 
  return false;
}

int FMMEikonalSolver_wrapper_solve(FMMEikonalSolver_wrapper * self, 
				   const double * phi0,
				   const double * nodalSpeeds,
				   double * T,
				   double zeroTol,
				   double trialTol,
				   int initFlag,
				   int verbose)
{
  int failed;
  assert(self);
  assert(self->pFMMsolver);
  //mwf hack
  //std::cout<<"FMMEikonalSolver_wrapper calling pFMMsolver->solve verbose= "<<verbose<<std::endl;
  failed = self->pFMMsolver->solve(phi0,
				   nodalSpeeds,
				   T,
				   zeroTol,
				   trialTol,
				   initFlag,
				   verbose);
  return failed;
  
}

//----------------------------------------------------------------------
int FSWEikonalSolver_wrapper_init(Mesh* cmesh,
				  int nSpace,
				  double atol, double rtol,
				  int maxIts,
				  int initFlag,
				  int nRefPoints,
				  const double* refPoints,
				  FSWEikonalSolver_wrapper * wrap)
{
  assert(wrap);
  wrap->pFSWsolver = NULL;
  FMMEikonalSolverBase::INIT_TYPE initType = FMMEikonalSolverBase::FRONT_AND_MAGNITUDE;
  if (initFlag == 0)
    initType = FMMEikonalSolverBase::MAGNITUDE;
  if (nSpace == 1)
    {
      wrap->pFSWsolver = new FSWEikonalSolver1d(cmesh,atol,rtol,maxIts,initType,
						nRefPoints,refPoints);
      wrap->nSpace = nSpace;
    }
  else if (nSpace == 2)
    {
      wrap->pFSWsolver = new FSWEikonalSolver2d(cmesh,atol,rtol,maxIts,initType,
						nRefPoints,refPoints);
      wrap->nSpace = nSpace;
    }
  else if (nSpace == 3)
    {
      wrap->pFSWsolver = new FSWEikonalSolver3d(cmesh,atol,rtol,maxIts,initType,
						nRefPoints,refPoints);
      wrap->nSpace = nSpace;
    }
  else
    {
      std::cout<<"nSpace= "<<nSpace<<" not implemented" <<std::endl;
      std::exit(1);
      return true;
    }
  return false;
}

int FSWEikonalSolver_wrapper_del(FSWEikonalSolver_wrapper * wrap)
{
  assert(wrap);
  if (wrap->pFSWsolver)
    delete wrap->pFSWsolver;
  return false;
}

int FSWEikonalSolver_wrapper_copy(FSWEikonalSolver_wrapper * from, 
				  FSWEikonalSolver_wrapper * to)
{
  assert(from); assert(to);
  if (to->pFSWsolver)
    delete to->pFSWsolver;
  if (from->nSpace == 1)
    {
      to->pFSWsolver = new FSWEikonalSolver1d(from->pFSWsolver->mesh,
					      from->pFSWsolver->iterAtol,
					      from->pFSWsolver->iterRtol,
					      from->pFSWsolver->maxIts,
					      from->pFSWsolver->initFlag,
					      from->pFSWsolver->nRefPoints,
					      &(from->pFSWsolver->refPoints[0]));
      to->nSpace = from->nSpace;
    }
  else if (from->nSpace == 2)
    {
      to->pFSWsolver = new FSWEikonalSolver2d(from->pFSWsolver->mesh,
					      from->pFSWsolver->iterAtol,
					      from->pFSWsolver->iterRtol,
					      from->pFSWsolver->maxIts,
					      from->pFSWsolver->initFlag,
					      from->pFSWsolver->nRefPoints,
					      &(from->pFSWsolver->refPoints[0]));
      to->nSpace = from->nSpace;
    }
  else if (from->nSpace == 3)
    {
      to->pFSWsolver = new FSWEikonalSolver3d(from->pFSWsolver->mesh,
					      from->pFSWsolver->iterAtol,
					      from->pFSWsolver->iterRtol,
					      from->pFSWsolver->maxIts,
					      from->pFSWsolver->initFlag,
					      from->pFSWsolver->nRefPoints,
					      &(from->pFSWsolver->refPoints[0]));
      to->nSpace = from->nSpace;
    }
  else
    {
      std::cout<<"nSpace= "<<from->nSpace<<" not implemented" <<std::endl;
      std::exit(1);
      return true;
    }
 
  return false;
}

int FSWEikonalSolver_wrapper_solve(FSWEikonalSolver_wrapper * self, 
				   const double * phi0,
				   const double * nodalSpeeds,
				   double * T,
				   double zeroTol,
				   double trialTol,
				   int initFlag,
				   int verbose)
{
  int failed;
  assert(self);
  assert(self->pFSWsolver);

  failed = self->pFSWsolver->solve(phi0,
				   nodalSpeeds,
				   T,
				   zeroTol,
				   trialTol,
				   initFlag,
				   verbose);
  return failed;
  
}

//----------------------------------------------------------------------
int localPWL_wrapper(int nSpace, Mesh* cmesh,
		     const double * phi0,
		     double * phi0R,
		     double zeroTol,
		     int verbose)
{
  int failed;
  assert(cmesh);
  failed = NearFrontInitialization::localPWLreconstruction(nSpace,cmesh,
							   phi0,phi0R,
							   zeroTol,verbose);
  return failed;
}






}//extern C

#include "tracking.h"

#include <iostream>
#include <cassert>
#include <cmath>
#include <list>

/**********************************************************************

   TODO:
     General
        setup tolerances to scale based on element diameter
     1D
        add linear in time, space velocity field 
            RK2
            RK2 element by element?
            exact integration following Tom's paper?
 **********************************************************************/
extern "C" void trackPointsConstVelocity1d(int nElements_global,               //mesh representation
					   int nNodes_global,
					   int nNodes_element,
					   int nElementBoundaries_element,
					   const double * nodeArray,
					   const int * elementNodesArray,
					   const int * elementNeighborsArray, //local boundary id is associated with node across from boundary 
					   double cVelocity,                  //characteristic speed (velocity) representation
					   double dir,                        //direction in time
					   const double *x_depart_times,            //in  -- time of departure for each point
                                           double *x_arrive_times,            //in  -- desired stopping time
					                                      //out -- stopping time for each point 
					   int nPointsToTrack,                //number of points to track
					   double zeroTolForTracking,         //ignore point if |u| < eps or |v| < eps
					   const double* x_in,                //points for tracking (always 3d)
					   int *x_element,                    //in -- element where point i is located at tIn
					                                      //out -- element where point i is located at tOut
					   double * x_out,                    //stopping location for point i at tOut
					   int * flag)                        //in: > -2  -- track point
                                                                              //in: -3    -- skip point
                                                                              //out: -1   -- point in interior at tOut
					                                      //out: -2   -- point exited domain somewhere in (tIn,tOut)
                                                                              //out: -3   -- did not track (e.g.,  or u = 0)
{
  using namespace std;
  const double LARGE = 1.2345e12;
  double v = cVelocity; //tracking with this velocity 
  const double zeroVelocityTol = zeroTolForTracking;
  const double zeroIncrementTol     = 1.0e-3*zeroTolForTracking;
  for (int k=0; k < nPointsToTrack; k++)
    {
      //1d only 
      double x = x_in[k*3+0];
      double t = x_depart_times[k];
      double tTarget = x_arrive_times[k];
      assert((tTarget-t)*dir >= 0.0);
      //current element and its nodes
      int eN = x_element[k];
      int nN_0 = elementNodesArray[eN*nNodes_element + 0], nN_1 = elementNodesArray[eN*nNodes_element + 1]; 
      double x_0 = nodeArray[nN_0*3+0], x_1 = nodeArray[nN_1*3+0];
      bool trackPoint = flag[k] >= -1 && fabs(v) > zeroVelocityTol;
      if (!trackPoint)
	{
	  x_element[k] = eN;
	  
	  x_arrive_times[k] = tTarget;
	  x_out[k*3+0] = x_in[k*3+0]; x_out[k*3+1] = x_in[k*3+1]; x_out[k*3+2] = x_in[k*3+2];
	}
      else
	{
	  bool done = false;
	  int exitedDomain = 0;
	  while (!done)
	    {
	      double dt_target = tTarget-t;
	      double dt = dt_target; //by default take full step
	      done = true;
	      int eN_exit = eN;
	      //exit times to boundaries
	      double tExit_0 = (x_0 - x)/v + t;
	      double tExit_1 = (x_1 - x)/v + t;
	      double dt_0 = LARGE*dir, dt_1 = LARGE*dir;
	      double dt_exit = LARGE*dir; int exitBoundary = -1;
	      if ((tExit_0-t)*dir > zeroIncrementTol) //possibly went out 0 with nonzero step
		{dt_0 = tExit_0-t; exitBoundary = 0;}
	      if ((tExit_1-t)*dir > zeroIncrementTol) //possibly went out 1 with nonzero step
		{dt_1 = tExit_1-t; exitBoundary = 1;}
	      //catch possibility that integration point started at boundary of element
	      //and left out of that location
	      if (exitBoundary == -1 && fabs(x_0-x) <= zeroIncrementTol)
		{
		  dt_0 = 0.0; exitBoundary = 0; dt_exit=dt_0;
		  dt = dt_exit;
		  int neigIndex = (exitBoundary + nElementBoundaries_element-1) % nElementBoundaries_element;
		  eN_exit = elementNeighborsArray[eN*nElementBoundaries_element + neigIndex];
		  if (eN_exit < 0)
		    exitedDomain = 1;
		  done = exitedDomain; //only done if left the physical domain
		  
		}
	      else if (exitBoundary == -1 && fabs(x_1-x) <= zeroIncrementTol)
		{
		  dt_1 = 0.0; exitBoundary = 1; dt_exit=dt_1; 
		  dt = dt_exit;
		  int neigIndex = (exitBoundary + nElementBoundaries_element-1) % nElementBoundaries_element;
		  eN_exit = elementNeighborsArray[eN*nElementBoundaries_element + neigIndex];
		  if (eN_exit < 0)
		    exitedDomain = 1;
		  done = exitedDomain; //only done if left the physical domain
		}
	      else
		{
		  if (fabs(dt_0) < fabs(dt_exit))
		    {
		      dt_exit = dt_0; exitBoundary = 0;
		    }
		  if (fabs(dt_1) < fabs(dt_exit))
		    {
		      dt_exit = dt_1; exitBoundary = 1;
		    }
		  if (dt_exit*dir > zeroIncrementTol && fabs(dt_exit) < fabs(dt) - zeroIncrementTol) //exited current element before reaching target time
		    {
		      assert(exitBoundary >= 0.0);
		      dt = dt_exit;
		      int neigIndex = (exitBoundary + nElementBoundaries_element-1) % nElementBoundaries_element;
		      eN_exit = elementNeighborsArray[eN*nElementBoundaries_element + neigIndex];
		      if (eN_exit < 0)
			exitedDomain = 1;
		      done = exitedDomain; //only done if left the physical domain
		    }
		}
	      //now update point and time
	      x += dt*v;
	      t += dt;
	      if (!done)
		{
		  eN = eN_exit;
		  nN_0 = elementNodesArray[eN*nNodes_element + 0]; nN_1 = elementNodesArray[eN*nNodes_element + 1]; 
		  x_0 = nodeArray[nN_0*3+0]; x_1 = nodeArray[nN_1*3+0];
		}
	    }//done

	  x_arrive_times[k] = t;
	  x_out[k*3+0] = x; x_out[k*3+1] = x_in[k*3+1]; x_out[k*3+2] = x_in[k*3+2];
	  flag[k] = -1;
	  if (exitedDomain)
	    flag[k] = -2;
	  x_element[k] = eN;
	}//need to track point
    }//points to track
}//simple element to element tracking

extern "C" void trackPointsC0P1Velocity1d(int nElements_global,               //mesh representation
					   int nNodes_global,
					   int nNodes_element,
					   int nElementBoundaries_element,
					   const double * nodeArray,
					   const int * elementNodesArray,
					   const int * elementNeighborsArray, //local boundary id is associated with node across from boundary 
					   const int * cvelocity_l2g,   //characteristic speed (velocity) representation
					   const double * cvelocity_dof,
					   double dir,                        //direction in time
					   const double *x_depart_times,            //in  -- time of departure for each point
                                           double *x_arrive_times,            //in  -- desired stopping time
					                                      //out -- stopping time for each point 
					   int nPointsToTrack,                //number of points to track
					   double zeroTolForTracking,         //ignore point if |u| < eps or |v| < eps
					   const double* x_in,                //points for tracking (always 3d)
					   int *x_element,                    //in -- element where point i is located at tIn
					                                      //out -- element where point i is located at tOut
					   double * x_out,                    //stopping location for point i at tOut
					   int * flag)                        //in: > -2  -- track point
                                                                              //in: -3    -- skip point
                                                                              //out: -1   -- point in interior at tOut
					                                      //out: -2   -- point exited domain somewhere in (tIn,tOut)
                                                                              //out: -3   -- did not track (e.g.,  or u = 0)

{
  using namespace std;
  using namespace Tracking;
  const double LARGE = 1.2345e12;
  const int nDOF_element_velocity = 2;
  const double zeroVelocityTol = zeroTolForTracking;
  const double zeroIncrementTol     = 1.0e-3*zeroTolForTracking;
  for (int k=0; k < nPointsToTrack; k++)
    {
      //1d only 
      double x = x_in[k*3+0];
      double t = x_depart_times[k];
      double tTarget = x_arrive_times[k];
      assert((tTarget-t)*dir >= 0.0);
      //current element and its nodes
      int eN = x_element[k];
      int nN_0 = elementNodesArray[eN*nNodes_element + 0], nN_1 = elementNodesArray[eN*nNodes_element + 1]; 
      double x_0 = nodeArray[nN_0*3+0], x_1 = nodeArray[nN_1*3+0];
      //use linear interpolation in physical space to get velocity at x
      double w[nDOF_element_velocity] = {0.0,0.0};
      double grad_w[nDOF_element_velocity] = {0.0,0.0};
      double x_eval[3] = {0.0,0.0,0.0};
      evaluateTestFunctionsOnElement_AffineLinear_1D(&x_in[k*3],
						     &nodeArray[nN_0*3],
						     &nodeArray[nN_1*3],
						     w);
      evaluateTestFunctionDerivativesOnElement_AffineLinear_1D(&x_in[k*3],
							       &nodeArray[nN_0*3],
							       &nodeArray[nN_1*3],
							       grad_w);
      double v = 0.0, vx = 0.0;
      for (int i = 0; i < nDOF_element_velocity; i++)
	{
	  v += w[i]*cvelocity_dof[cvelocity_l2g[eN*nDOF_element_velocity + i]];
	  vx+= grad_w[i]*cvelocity_dof[cvelocity_l2g[eN*nDOF_element_velocity + i]];
	}
      bool trackPoint = flag[k] >= -1 && fabs(v) > zeroVelocityTol;
      if (!trackPoint)
	{
	  x_element[k] = eN;

	  x_arrive_times[k] = tTarget;
	  x_out[k*3+0] = x_in[k*3+0]; x_out[k*3+1] = x_in[k*3+1]; x_out[k*3+2] = x_in[k*3+2];
	}
      else
	{
	  bool done = false;
	  int exitedDomain = 0;
	  while (!done)
	    {
	      double dt_target = tTarget-t;
	      double dt = dt_target; //by default take full step
	      done = true;
	      int eN_exit = eN;
	      //exit times to boundaries assuming |vx| = 0
	      double tExit_0 = (x_0 - x)/v + t;
	      double tExit_1 = (x_1 - x)/v + t;
	      double dt_0 = LARGE*dir, dt_1 = LARGE*dir;
	      double dt_exit = LARGE*dir; int exitBoundary = -1;
	      //account for linear variation
	      if (fabs(vx) > zeroVelocityTol)
		{
		  //can it go out x_0 or x_1 
		  double dx_0 = x_0-x;
		  if (dx_0*vx/v <= -1.0)
		    tExit_0 = t -dir*LARGE; //no way to hit x_0
		  else
		    tExit_0 = t + log(dx_0*vx/v + 1.)/vx;
		  double dx_1 = x_1-x;
		  if (dx_1*vx/v <= -1.0)
		    tExit_1 = t -dir*LARGE; //no way to hit x_1
		  else
		    tExit_1 = t + log(dx_1*vx/v + 1.)/vx;
		}

	      if ((tExit_0-t)*dir > zeroIncrementTol) //possibly went out 0 with nonzero step
		{dt_0 = tExit_0-t; exitBoundary = 0;}
	      if ((tExit_1-t)*dir > zeroIncrementTol) //possibly went out 1 with nonzero step
		{dt_1 = tExit_1-t; exitBoundary = 1;}
	      //catch possibility that integration point started at boundary of element
	      //and left out of that location
	      if (exitBoundary == -1 && fabs(x_0-x) <= zeroIncrementTol)
		{
		  dt_0 = 0.0; exitBoundary = 0; dt_exit=dt_0;
		  dt = dt_exit;
		  int neigIndex = (exitBoundary + nElementBoundaries_element-1) % nElementBoundaries_element;
		  eN_exit = elementNeighborsArray[eN*nElementBoundaries_element + neigIndex];
		  if (eN_exit < 0)
		    exitedDomain = 1;
		  done = exitedDomain; //only done if left the physical domain
		  
		}
	      else if (exitBoundary == -1 && fabs(x_1-x) <= zeroIncrementTol)
		{
		  dt_1 = 0.0; exitBoundary = 1; dt_exit=dt_1; 
		  dt = dt_exit;
		  int neigIndex = (exitBoundary + nElementBoundaries_element-1) % nElementBoundaries_element;
		  eN_exit = elementNeighborsArray[eN*nElementBoundaries_element + neigIndex];
		  if (eN_exit < 0)
		    exitedDomain = 1;
		  done = exitedDomain; //only done if left the physical domain
		}
	      else
		{
		  if (fabs(dt_0) < fabs(dt_exit))
		    {
		      dt_exit = dt_0; exitBoundary = 0;
		    }
		  if (fabs(dt_1) < fabs(dt_exit))
		    {
		      dt_exit = dt_1; exitBoundary = 1;
		    }
		  if (dt_exit*dir > zeroVelocityTol && fabs(dt_exit) < fabs(dt) -zeroIncrementTol) //exited current element before reaching target time
		    {
		      assert(exitBoundary >= 0.0);
		      dt = dt_exit;
		      int neigIndex = (exitBoundary + nElementBoundaries_element-1) % nElementBoundaries_element;
		      eN_exit = elementNeighborsArray[eN*nElementBoundaries_element + neigIndex];
		      if (eN_exit < 0)
			exitedDomain = 1;
		      done = exitedDomain; //only done if left the physical domain
		    }
		}
	      //now update point and time
	      t += dt;
	      if (fabs(vx) < zeroVelocityTol)
		{
		  x += dt*v;
		}
	      else
		{
		  x += v*(exp(vx*dt)-1.0)/vx;
		}
	      if (!done)
		{
		  eN = eN_exit;
		  nN_0 = elementNodesArray[eN*nNodes_element + 0]; nN_1 = elementNodesArray[eN*nNodes_element + 1]; 
		  x_0 = nodeArray[nN_0*3+0]; x_1 = nodeArray[nN_1*3+0];
		  x_eval[0] = x;
		  evaluateTestFunctionsOnElement_AffineLinear_1D(x_eval,
								 &nodeArray[nN_0*3],
								 &nodeArray[nN_1*3],
								 w);
		  evaluateTestFunctionDerivativesOnElement_AffineLinear_1D(x_eval,
									   &nodeArray[nN_0*3],
									   &nodeArray[nN_1*3],
									   grad_w);
		  v = 0.0; vx = 0.0;
		  for (int i = 0; i < nDOF_element_velocity; i++)
		    {
		      v += w[i]*cvelocity_dof[cvelocity_l2g[eN*nDOF_element_velocity + i]];
		      vx+= grad_w[i]*cvelocity_dof[cvelocity_l2g[eN*nDOF_element_velocity + i]];
		    }
		  
		}
	    }//done

	  x_arrive_times[k] = t;
	  x_out[k*3+0] = x; x_out[k*3+1] = x_in[k*3+1]; x_out[k*3+2] = x_in[k*3+2];
	  flag[k] = -1; 
	  if (exitedDomain) 
	    flag[k] = -2;
	  x_element[k] = eN;
	}//need to track point
    }//points to track
}//C0P1 1d

extern "C" void trackPointsRT0Velocity2d(int debugLevel,
					 int localVelocityRepresentationFlag, //2 -- 'flux' dofs, 1 -- \vec a + b \vec x
					 int nElements_global,               //mesh representation
					 int nNodes_global,
					 int nNodes_element,
					 int nElementBoundaries_element,
					 const double * nodeArray,
					 const int * elementNodesArray,
					 const int * elementNeighborsArray, //local boundary id is associated with node across from boundary 
					 const int * elementBoundariesArray, //element --> global element boundary id
					 const double * elementBoundaryBarycentersArray,  //\vec{\bar{x}}_f [nElementBoundaries_global,3]
					 const double * elementLocalBoundaryOuterNormalsArray, //\vec n_{f}      [nElements_global,nElementBoundaries_element,nd]
					 const int * cvelocity_l2g,
					 const double * cvelocity_dof,   //characteristic speed (velocity) representation
					                                //v = \sum_{i=0}^{n_d} V^i \vec N_i
                                                                        //\vec N_i = \frac{1}{d|\Omega_e|}(\vec x - \vec x_{n,i}) i=0,..,n_d
					 double dir,                    //direction in time
					 const double *x_depart_times,  //in  -- time of departure for each point
					 double *x_arrive_times,        //in  -- desired stopping time
					                                //out -- stopping time for each point 
					 int nPointsToTrack,            //number of points to track
					 double zeroTolForTracking,     //ignore point if |u| < eps or |v| < eps
					 const double* x_in,            //points for tracking (always 3d)
					 int *x_element,                //in -- element where point i is located at tIn
					                                //out -- element where point i is located at tOut
					 double * x_out,                //stopping location for point i at tOut
					 int * flag)                    //in: > -2  -- track point
                                                                        //in: -3    -- skip point
                                                                        //out: -1   -- point in interior at tOut
					                                //out: -2   -- point exited domain somewhere in (tIn,tOut)
                                                                        //out: -3   -- did not track (e.g.,  or u = 0)
{
  using namespace std;
  using namespace Tracking;
  const double LARGE = 1.2345e12;
  const int nDOF_element_velocity = 3;
  const int nSpace = 2;
  const double zeroVelocityTol = zeroTolForTracking;
  const double zeroIncrementTol     = 1.0e-3*zeroTolForTracking;
  const double insideElementTol     = 1.0e-5;
  for (int k=0; k < nPointsToTrack; k++)
    {
      if (flag[k] >= -1)
	{
	  register double x[3] = {0.0,0.0,0.0};
	  x[0] = x_in[k*3+0]; x[1] = x_in[k*3+1]; x[2] = x_in[k*3+2]; 

	  double t = x_depart_times[k];
	  double tTarget = x_arrive_times[k];
	  assert((tTarget-t)*dir >= 0.0);
	  //current element
	  int eN = x_element[k];

	  //
	  //perform sanity check that actually in element
	  /***************************************************
         recall in element check for straight sided polygon
           is that point lines within intersection of all
           negative half planes of the faces
  
         face f: (\vec x - \vec \bar{x}_f) \cdot \vec n_f = 0        
	  **************************************************/
	  for (int ebN = 0; ebN < nElementBoundaries_element; ebN++)
	    {
	      const int ebN_global = elementBoundariesArray[eN*nElementBoundaries_element+ebN];
	      double dxf = 0.0;
	      for (int I = 0; I < nSpace; I++)
		{
		  //note switching order of x-x_f in plane definition to match convention for \Delta x in 1d and Haselbacher JCP paper
		  dxf += (elementBoundaryBarycentersArray[ebN_global*3 + I] - x[I])*elementLocalBoundaryOuterNormalsArray[eN*nElementBoundaries_element*nSpace + ebN*nSpace + I];
		}
	      if (dxf < -zeroIncrementTol) //outside element
		{
		  std::cout<<"PROBLEM trackPointsRT0Velocity2d initial element check k= "<<k<<" eN= "<<eN<<" x= ["<<x[0]<<","<<x[1]<<","<<x[2]
			   <<"] ebN= "<<ebN<<" ebN_global= "<<ebN_global<<" n= ["
			   <<elementLocalBoundaryOuterNormalsArray[eN*nElementBoundaries_element*nSpace + ebN*nSpace + 0]
			   <<","
			   <<elementLocalBoundaryOuterNormalsArray[eN*nElementBoundaries_element*nSpace + ebN*nSpace + 1]
			   <<"] dxf= "<<dxf<<std::endl;
		  //move point to boundary
		  if (dxf > -insideElementTol)
		    {
		      for (int I = 0; I < nSpace; I++)
			x[I] += dxf*elementLocalBoundaryOuterNormalsArray[eN*nElementBoundaries_element*nSpace + ebN*nSpace + I];
		      //double check
		      dxf = 0.0;
		      for (int I = 0; I < nSpace; I++)
			{
			  //note switching order of x-x_f in plane definition to match convention for \Delta x in 1d and Haselbacher JCP paper
			  dxf += (elementBoundaryBarycentersArray[ebN_global*3 + I] - x[I])*elementLocalBoundaryOuterNormalsArray[eN*nElementBoundaries_element*nSpace + ebN*nSpace + I];
			}
		    }//move to element
		  assert(dxf >= -zeroIncrementTol);
		}
	    }
	  int nN_0 = elementNodesArray[eN*nNodes_element + 0], nN_1 = elementNodesArray[eN*nNodes_element + 1], nN_2 = elementNodesArray[eN*nNodes_element + 2];
	  //evaluate velocity at current location
	  double w[nDOF_element_velocity][2];
	  double grad_w[nDOF_element_velocity][2][2];
	  if (localVelocityRepresentationFlag == 1)
	    {
	      evaluateTestFunctionsOnElement_RT0_simplex_2D_axbRep(&x_in[k*3],
								   &nodeArray[nN_0*3],
								   &nodeArray[nN_1*3],
								   &nodeArray[nN_2*3],
								   w);
	      evaluateTestFunctionDerivativesOnElement_RT0_simplex_2D_axbRep(&x_in[k*3],
									     &nodeArray[nN_0*3],
									     &nodeArray[nN_1*3],
									     &nodeArray[nN_2*3],
									     grad_w);

	    }
	  else
	    {
	      evaluateTestFunctionsOnElement_RT0_simplex_2D(&x_in[k*3],
							    &nodeArray[nN_0*3],
							    &nodeArray[nN_1*3],
							    &nodeArray[nN_2*3],
							    w);
	      evaluateTestFunctionDerivativesOnElement_RT0_simplex_2D(&x_in[k*3],
								      &nodeArray[nN_0*3],
								      &nodeArray[nN_1*3],
								      &nodeArray[nN_2*3],
								      grad_w);
	    }


	  register double v[2],grad_v[2][2],
	    vx = 0.0, normv = 0.0;

	  for (int I = 0; I < nSpace; I++)
	    {
	      v[I] = 0.0; 
	      for (int i = 0; i < nDOF_element_velocity; i++)
		{
		  v[I] += w[i][I]*cvelocity_dof[cvelocity_l2g[eN*nDOF_element_velocity+i]];
		}
	      normv += v[I]*v[I];
	      for (int J = 0; J < nSpace; J++)
		{
		  grad_v[I][J] = 0.0;
		  for (int i = 0; i < nDOF_element_velocity; i++)
		    {
		      grad_v[I][J] += grad_w[i][I][J]*cvelocity_dof[cvelocity_l2g[eN*nDOF_element_velocity+i]];
		    }
		}
	    }
	  normv = sqrt(normv); 
	  assert(fabs(grad_v[0][0]-grad_v[1][1]) < zeroVelocityTol);
	  vx = grad_v[0][0]; 
	  //mwf debug
	  if (debugLevel > 0)
	    std::cout<<"trackPointsRT0Velocity2d k= "<<k<<" eN= "<<eN<<" x= ["<<x[0]<<","<<x[1]<<"], initial v= ["<<v[0]<<","<<v[1]<<"] grad_v[0][0]= "<<grad_v[0][0]<<std::endl;
	  //TODO refine this test since flag[k] >= 0 should be true always at this point
	  bool trackPoint = flag[k] >= -1 && fabs(normv) > zeroVelocityTol;
	  if (!trackPoint)
	    {
	      x_element[k] = eN;

	      x_arrive_times[k] = tTarget;
	      x_out[k*3+0] = x_in[k*3+0]; x_out[k*3+1] = x_in[k*3+1]; x_out[k*3+2] = x_in[k*3+2];
	    }
	  else
	    {
	      bool done = false;
	      int exitedDomain = 0;
	      while (!done)
		{
		  double dt_target = tTarget-t;
		  double dt = dt_target; //by default take full step
		  //mwf debug
		  if (debugLevel > 0)
		    std::cout<<"trackPointsRT0Velocity2d eN= "<<eN<<" k= "<<k<<" tTarget= "<<tTarget<<" t= "<<t<<" dt_target= "<<dt_target<<std::endl;
		  done = true;
		  int eN_exit = eN;
		  /***************************************************
	        in 2d, put check to make sure going into element
	          before calculating exit times because could have exit time along correct time axis
	          with a face but with an intersection that's outside the element if the point
	          never actually tracked into the element

                point leaving element only if 

                  |(\vec \bar{x}_f - \vec x_0)\cdot n_f| < \epsilon and
                     \vec v_0 \cdot \vec n_f > 0 
                  for some face f

                TODO: determine what can happen if v.n_f = 0, should have
                      intersection with another face that fits usual paradigm
                TODO: get tolerances straight on v.n test
		  ***************************************************/
		  bool enteredElement = true;
		  //do not reevaluate velocity if point "exits" an element without entering
		  // that is the point is on the boundary of the element (e.g., a node) but v.n > 0
                  bool evaluateVelocity = true; 
		  for (int ebN = 0; ebN < nElementBoundaries_element; ebN++)
		    {
		      const int ebN_global = elementBoundariesArray[eN*nElementBoundaries_element+ebN];
		      double dxf = 0.0, vdotn = 0.0;
		      for (int I = 0; I < nSpace; I++)
			{
			  dxf  += (elementBoundaryBarycentersArray[ebN_global*3 + I] - x[I])*elementLocalBoundaryOuterNormalsArray[eN*nElementBoundaries_element*nSpace + ebN*nSpace + I];
			  vdotn+= v[I]*elementLocalBoundaryOuterNormalsArray[eN*nElementBoundaries_element*nSpace + ebN*nSpace + I];
			}
		      if (dxf < -zeroIncrementTol) //outside element
			{
			  std::cout<<"PROBLEM trackPointsRT0Velocity2d done loop element check eN= "<<eN<<" x= ["<<x[0]<<","<<x[1]<<","<<x[2]
				   <<"] ebN= "<<ebN<<" ebN_global= "<<ebN_global<<" n= ["
				   <<elementLocalBoundaryOuterNormalsArray[eN*nElementBoundaries_element*nSpace + ebN*nSpace + 0]
				   <<","
				   <<elementLocalBoundaryOuterNormalsArray[eN*nElementBoundaries_element*nSpace + ebN*nSpace + 1]
				   <<"] dxf= "<<dxf<<std::endl;
			  std::cout<<"global Nodes on element : "
				   <<nN_0<<" --> ["<<nodeArray[nN_0*3+0]<<","<<nodeArray[nN_0*3+1]<<","<<nodeArray[nN_0*3+2]<<"] "
				   <<nN_1<<" --> ["<<nodeArray[nN_1*3+0]<<","<<nodeArray[nN_1*3+1]<<","<<nodeArray[nN_1*3+2]<<"] "
				   <<nN_2<<" --> ["<<nodeArray[nN_2*3+0]<<","<<nodeArray[nN_2*3+1]<<","<<nodeArray[nN_2*3+2]<<"] "
				   <<std::endl;
			  //move point to boundary
			  if (dxf > -insideElementTol)
			    {
			      for (int I = 0; I < nSpace; I++)
				x[I] += dxf*elementLocalBoundaryOuterNormalsArray[eN*nElementBoundaries_element*nSpace + ebN*nSpace + I];
			      //double check
			      dxf = 0.0;
			      for (int I = 0; I < nSpace; I++)
				{
				  //note switching order of x-x_f in plane definition to match convention for \Delta x in 1d and Haselbacher JCP paper
				  dxf += (elementBoundaryBarycentersArray[ebN_global*3 + I] - x[I])*elementLocalBoundaryOuterNormalsArray[eN*nElementBoundaries_element*nSpace + ebN*nSpace + I];
				}
			    }//move to element
			  assert(dxf >= -zeroIncrementTol);
			}
		      else if (fabs(dxf) < zeroIncrementTol && vdotn*dir > zeroIncrementTol) 
			{
			  enteredElement  = false;
			  evaluateVelocity= false; 
			  eN_exit = elementNeighborsArray[eN*nElementBoundaries_element + ebN];
			  dt = 0.0;
			  if (eN_exit < 0)
			    exitedDomain = 1;
			  done = exitedDomain; //only done if left the physical domain
			  //mwf debug
			  if (debugLevel > 1)
			    std::cout<<"trackPointsRT0Velocity2d on boundary did not enter domain eN= "<<eN<<" ebN= "<<ebN<<" dxf= "<<dxf
				     <<" vdotn= "<<vdotn<<" eN_exit= "<<eN_exit<<std::endl;
			}
		    }
		  if (enteredElement)
		    {
		      register double t_bnd[3],dt_bnd[3],dxf[3],vdotn[3];
		      //exit times to boundaries assuming |vx| = 0
		      for (int ebN = 0; ebN < nElementBoundaries_element; ebN++)
			{
			  const int ebN_global = elementBoundariesArray[eN*nElementBoundaries_element + ebN];
			  dxf[ebN] = 0.0, vdotn[ebN] = 0.0;
			  for (int I = 0; I < nSpace; I++)
			    {
			      dxf[ebN]  += (elementBoundaryBarycentersArray[ebN_global*3 + I] - x[I])*elementLocalBoundaryOuterNormalsArray[eN*nElementBoundaries_element*nSpace + ebN*nSpace + I];
			      vdotn[ebN]+= v[I]*elementLocalBoundaryOuterNormalsArray[eN*nElementBoundaries_element*nSpace + ebN*nSpace + I];
			    }
			  //by default assume no way to hit boundary
			  //set to invalid t exit
			  t_bnd[ebN] = t - LARGE*dir;
			  if (fabs(vdotn[ebN]) > zeroVelocityTol)
			    {
			      t_bnd[ebN] = dxf[ebN]/vdotn[ebN] + t;
			    }
			  //mwf debug
			  if (debugLevel > 2)
			    std::cout<<"eN= "<<eN<<" ebN= "<<ebN<<" ebN_global= "<<ebN_global<<" x= ["<<x[0]<<","<<x[1]<<"] dxf= "<<dxf[ebN]<<" vdotn= "<<vdotn[ebN]<<" t= "<<t<<" t_bnd= "<<t_bnd[ebN]<<std::endl;
			  //set dt to boundary to a big value so won't get chosen in min step
			  dt_bnd[ebN] = LARGE*dir;
			}
		      //account for linear variation
		      if (fabs(vx) > zeroVelocityTol)
			{
			  for (int ebN = 0; ebN < nElementBoundaries_element; ebN++)
			    {
			      if (fabs(vdotn[ebN]) > zeroVelocityTol)
				{
				  const double alpha = dxf[ebN]/vdotn[ebN];
				  if (alpha*vx <= -1.0) //no way to hit boundary, set to invalid time
				    t_bnd[ebN] = t - LARGE*dir;
				  else
				    t_bnd[ebN] = t + log(alpha*vx + 1.0)/vx;
				}
			    }
			}
		      //final exit time and boundary 
		      double dt_exit = LARGE*dir; int exitBoundary = -1;
		      //determine valid exit time steps after calculating exit times
		      for (int ebN = 0; ebN < nElementBoundaries_element; ebN++)
			{
			  if ((t_bnd[ebN]-t)*dir > zeroIncrementTol) //possibly went out ebN with nonzero step
			    {
			      dt_bnd[ebN] = t_bnd[ebN]-t;
			      //see if this valid exit time is the smallest one
			      if (fabs(dt_bnd[ebN]) < fabs(dt_exit))
				{
				  dt_exit = dt_bnd[ebN];
				  exitBoundary = ebN;
				}
			    }
			}
		      //mwf debug
		      if (debugLevel > 0)
			std::cout<<"after calculating exit times eN = "<<eN<<" dt_exit= "<<dt_exit<<" exitBoundary= "<<exitBoundary<<" dt_bnd= ["<<dt_bnd[0]<<","<<dt_bnd[1]<<","<<dt_bnd[2]<<"] "<<std::endl;

		      if (dt_exit*dir > zeroIncrementTol && fabs(dt_exit) < fabs(dt) - zeroIncrementTol) //exited current element before reaching target time
			{
			  assert(exitBoundary >= 0.0 && exitBoundary < nElementBoundaries_element);
			  dt = dt_exit;
			  int neigIndex = exitBoundary;
			  eN_exit = elementNeighborsArray[eN*nElementBoundaries_element + neigIndex];
			  if (eN_exit < 0)
			    exitedDomain = 1;
			  done = exitedDomain; //only done if left the physical domain
			  //mwf debug
			  if (debugLevel > 0)
			    {
			      std::cout<<"eN= "<<eN<<" x= ["<<x[0]<<","<<x[1]<<"] leaving through boundary "<<exitBoundary<<" globalBoundary= "<<elementBoundariesArray[eN*nElementBoundaries_element+exitBoundary]<<" to neighbor "<<eN_exit<<" v= ["<<v[0]<<","<<v[1]<<"]= vx= "<<vx
				       <<" dt= "<<dt<<" evaluateVelocity= "<<evaluateVelocity<<std::endl;
			      std::cout<<"eN elementBoundaries= ["<<elementBoundariesArray[eN*nElementBoundaries_element+0]
				       <<","<<elementBoundariesArray[eN*nElementBoundaries_element+1]
				       <<","<<elementBoundariesArray[eN*nElementBoundaries_element+2]<<"] "
				       <<"eN neighbors= ["<<elementNeighborsArray[eN*nElementBoundaries_element+0]
				       <<","<<elementNeighborsArray[eN*nElementBoundaries_element+1]
				       <<","<<elementNeighborsArray[eN*nElementBoundaries_element+2]<<"] "
				       <<std::endl;
			      if (eN_exit >= 0)
				{
				  std::cout<<"eN_exit elementBoundaries= ["<<elementBoundariesArray[eN_exit*nElementBoundaries_element+0]
					   <<","<<elementBoundariesArray[eN_exit*nElementBoundaries_element+1]
					   <<","<<elementBoundariesArray[eN_exit*nElementBoundaries_element+2]<<"] "
					   <<"eN_exit neighbors= ["<<elementNeighborsArray[eN_exit*nElementBoundaries_element+0]
					   <<","<<elementNeighborsArray[eN_exit*nElementBoundaries_element+1]
					   <<","<<elementNeighborsArray[eN_exit*nElementBoundaries_element+2]<<"] "
					   <<std::endl;
				}
			    }
			}
		    }//end entered element
		  //now update point and time
		  t += dt;
		  if (fabs(vx) < zeroVelocityTol)
		    {
		      for (int I = 0; I < nSpace; I++)
			x[I] += dt*v[I];
		    }
		  else
		    {
		      //mwf debug
		      if (debugLevel > 1)
			std::cout<<"updating x with linear formula x= ["<<x[0]<<","<<x[1]<<"] v= ["<<v[0]<<","<<v[1]<<"]= vx= "<<vx<<" dt= "<<dt;
		      for (int I = 0; I < nSpace; I++)
			x[I] += v[I]*(exp(vx*dt)-1.0)/vx;
		      //mwf debug
		      if (debugLevel > 1)
			std::cout<<" --> ["<<x[0]<<","<<x[1]<<"]"<<std::endl;
		    }
		  if (!done)
		    {
		      eN = eN_exit;

		      nN_0 = elementNodesArray[eN*nNodes_element + 0]; nN_1 = elementNodesArray[eN*nNodes_element + 1]; nN_2 = elementNodesArray[eN*nNodes_element + 2];
		      if (localVelocityRepresentationFlag == 1)
			{
			  evaluateTestFunctionsOnElement_RT0_simplex_2D_axbRep(x,
									       &nodeArray[nN_0*3],
									       &nodeArray[nN_1*3],
									       &nodeArray[nN_2*3],
									       w);
			  evaluateTestFunctionDerivativesOnElement_RT0_simplex_2D_axbRep(x,
											 &nodeArray[nN_0*3],
											 &nodeArray[nN_1*3],
											 &nodeArray[nN_2*3],
											 grad_w);

			}
		      else
			{
			  evaluateTestFunctionsOnElement_RT0_simplex_2D(x,
									&nodeArray[nN_0*3],
									&nodeArray[nN_1*3],
									&nodeArray[nN_2*3],
									w);
			  evaluateTestFunctionDerivativesOnElement_RT0_simplex_2D(x,
										  &nodeArray[nN_0*3],
										  &nodeArray[nN_1*3],
										  &nodeArray[nN_2*3],
										  grad_w);
			}

		      //evaluate velocity at current location
		      if (evaluateVelocity)
			{
			  normv = 0.0;
			  for (int I = 0; I < nSpace; I++)
			    {
			      v[I] = 0.0; 
			      for (int i = 0; i < nDOF_element_velocity; i++)
				{
				  v[I] += w[i][I]*cvelocity_dof[cvelocity_l2g[eN*nDOF_element_velocity+i]];
				}
			      normv += v[I]*v[I];
			      for (int J = 0; J < nSpace; J++)
				{
				  grad_v[I][J] = 0.0;
				  for (int i = 0; i < nDOF_element_velocity; i++)
				    {
				      //mwf debug
				      if (debugLevel > 2)
					std::cout<<"done loop evaluting velocity eN= "<<eN<<" x= ["<<x[0]<<","<<x[1]<<"], t= "<<t<<" grad_w["<<i<<"]["<<I<<"]["<<J<<"]= "<<grad_w[i][I][J]
						 <<" cvelocity_dof= "<<cvelocity_dof[cvelocity_l2g[eN*nDOF_element_velocity+i]]<<std::endl;
				      grad_v[I][J] += grad_w[i][I][J]*cvelocity_dof[cvelocity_l2g[eN*nDOF_element_velocity+i]];
				    }
				}
			    }
			  normv = sqrt(normv); 
			  //take advantage of RT0 velocity properties to simplify
			  assert(fabs(grad_v[0][0]-grad_v[1][1]) < zeroVelocityTol);
			  vx = grad_v[0][0]; 
			}//need to evaluate velocity
		      //mwf debug
		      if (debugLevel > 0 )
			std::cout<<"trackPointsRT0Velocity2d k= "<<k<<" eN= "<<eN<<" x= ["<<x[0]<<","<<x[1]<<"], t= "<<t<<" v= ["<<v[0]<<","<<v[1]<<"] grad_v[0][0]= "<<grad_v[0][0]<<std::endl;
		    }
		}//done

	      x_arrive_times[k] = t;
	      x_out[k*3+0] = x[0]; x_out[k*3+1] = x[1]; x_out[k*3+2] = x[2];
	      flag[k] = -1;
	      if (exitedDomain)
		flag[k] = -2;
	      x_element[k] = eN;
	    }//need to track point
	}//need to track point
    }//points to track

 }//simple element to element tracking

extern "C" void trackPointsRT0Velocity2dWithTrajectories(int debugLevel,
					 int localVelocityRepresentationFlag, //2 -- 'flux' dofs, 1 -- \vec a + b \vec x
					 int nElements_global,               //mesh representation
					 int nNodes_global,
					 int nNodes_element,
					 int nElementBoundaries_element,
					 const double * nodeArray,
					 const int * elementNodesArray,
					 const int * elementNeighborsArray, //local boundary id is associated with node across from boundary 
					 const int * elementBoundariesArray, //element --> global element boundary id
					 const double * elementBoundaryBarycentersArray,  //\vec{\bar{x}}_f [nElementBoundaries_global,3]
					 const double * elementLocalBoundaryOuterNormalsArray, //\vec n_{f}      [nElements_global,nElementBoundaries_element,nd]
					 const int * cvelocity_l2g,
					 const double * cvelocity_dof,   //characteristic speed (velocity) representation
					                                //v = \sum_{i=0}^{n_d} V^i \vec N_i
                                                                        //\vec N_i = \frac{1}{d|\Omega_e|}(\vec x - \vec x_{n,i}) i=0,..,n_d
					 double dir,                    //direction in time
					 const double *x_depart_times,  //in  -- time of departure for each point
					 double *x_arrive_times,        //in  -- desired stopping time
					                                //out -- stopping time for each point 
					 int nPointsToTrack,            //number of points to track
					 double zeroTolForTracking,     //ignore point if |u| < eps or |v| < eps
					 const double* x_in,            //points for tracking (always 3d)
					 int *x_element,                //in -- element where point i is located at tIn
					                                //out -- element where point i is located at tOut
					 double * x_out,                //stopping location for point i at tOut
					 int * flag,                    //in: > -2  -- track point
                                                                        //in: -3    -- skip point
                                                                        //out: -1   -- point in interior at tOut
					                                //out: -2   -- point exited domain somewhere in (tIn,tOut)
                                                                        //out: -3   -- did not track (e.g.,  or u = 0)
					 int & n_traj,                  //out: number of points in path
					 int & n_tracked,               //out: number of points actually tracked
					 int*& offsets_traj, 		 //out: starting locations for particle trajectories
					 double *& x_traj,               //out: x for path taken by point
					 double *& t_traj,               //out: t for path taken by point
					 int *& e_traj)                  //out: element id for path taken by point
 
{
  using namespace std;
  using namespace Tracking;
  const double LARGE = 1.2345e12;
  const int nDOF_element_velocity = 3;
  const int nSpace = 2;
  const double zeroVelocityTol = zeroTolForTracking;
  const double zeroIncrementTol     = 1.0e-3*zeroTolForTracking;
  const double insideElementTol     = 1.0e-5;
  std::list<double> x_list,t_list;
  std::list<int> e_list,offsets_list;
  n_traj = 0; n_tracked = 0;
  
  for (int k=0; k < nPointsToTrack; k++)
    {
      if (flag[k] >= -1)
	{
	  register double x[3] = {0.0,0.0,0.0};
	  x[0] = x_in[k*3+0]; x[1] = x_in[k*3+1]; x[2] = x_in[k*3+2]; 

	  double t = x_depart_times[k];
	  double tTarget = x_arrive_times[k];
	  assert((tTarget-t)*dir >= 0.0);
	  //current element
	  int eN = x_element[k];

	  //update trajectory
	  offsets_list.push_back(n_traj);
	  n_tracked++;
	  n_traj++;
	  t_list.push_back(t);
	  for (int I=0; I < 3; I++)
	    x_list.push_back(x[I]);
	  e_list.push_back(eN);
	  //
	  //perform sanity check that actually in element
	  /***************************************************
         recall in element check for straight sided polygon
           is that point lines within intersection of all
           negative half planes of the faces
  
         face f: (\vec x - \vec \bar{x}_f) \cdot \vec n_f = 0        
	  **************************************************/
	  for (int ebN = 0; ebN < nElementBoundaries_element; ebN++)
	    {
	      const int ebN_global = elementBoundariesArray[eN*nElementBoundaries_element+ebN];
	      double dxf = 0.0;
	      for (int I = 0; I < nSpace; I++)
		{
		  //note switching order of x-x_f in plane definition to match convention for \Delta x in 1d and Haselbacher JCP paper
		  dxf += (elementBoundaryBarycentersArray[ebN_global*3 + I] - x[I])*elementLocalBoundaryOuterNormalsArray[eN*nElementBoundaries_element*nSpace + ebN*nSpace + I];
		}
	      if (dxf < -zeroIncrementTol) //outside element
		{
		  std::cout<<"PROBLEM trackPointsRT0Velocity2d initial element check k= "<<k<<" eN= "<<eN<<" x= ["<<x[0]<<","<<x[1]<<","<<x[2]
			   <<"] ebN= "<<ebN<<" ebN_global= "<<ebN_global<<" n= ["
			   <<elementLocalBoundaryOuterNormalsArray[eN*nElementBoundaries_element*nSpace + ebN*nSpace + 0]
			   <<","
			   <<elementLocalBoundaryOuterNormalsArray[eN*nElementBoundaries_element*nSpace + ebN*nSpace + 1]
			   <<"] dxf= "<<dxf<<std::endl;
		  //move point to boundary
		  if (dxf > -insideElementTol)
		    {
		      for (int I = 0; I < nSpace; I++)
			x[I] += dxf*elementLocalBoundaryOuterNormalsArray[eN*nElementBoundaries_element*nSpace + ebN*nSpace + I];
		      //double check
		      dxf = 0.0;
		      for (int I = 0; I < nSpace; I++)
			{
			  //note switching order of x-x_f in plane definition to match convention for \Delta x in 1d and Haselbacher JCP paper
			  dxf += (elementBoundaryBarycentersArray[ebN_global*3 + I] - x[I])*elementLocalBoundaryOuterNormalsArray[eN*nElementBoundaries_element*nSpace + ebN*nSpace + I];
			}
		    }//move to element
		  assert(dxf >= -zeroIncrementTol);
		}
	    }
	  int nN_0 = elementNodesArray[eN*nNodes_element + 0], nN_1 = elementNodesArray[eN*nNodes_element + 1], nN_2 = elementNodesArray[eN*nNodes_element + 2];
	  //evaluate velocity at current location
	  double w[nDOF_element_velocity][2];
	  double grad_w[nDOF_element_velocity][2][2];
	  if (localVelocityRepresentationFlag == 1)
	    {
	      evaluateTestFunctionsOnElement_RT0_simplex_2D_axbRep(&x_in[k*3],
								   &nodeArray[nN_0*3],
								   &nodeArray[nN_1*3],
								   &nodeArray[nN_2*3],
								   w);
	      evaluateTestFunctionDerivativesOnElement_RT0_simplex_2D_axbRep(&x_in[k*3],
									     &nodeArray[nN_0*3],
									     &nodeArray[nN_1*3],
									     &nodeArray[nN_2*3],
									     grad_w);

	    }
	  else
	    {
	      evaluateTestFunctionsOnElement_RT0_simplex_2D(&x_in[k*3],
							    &nodeArray[nN_0*3],
							    &nodeArray[nN_1*3],
							    &nodeArray[nN_2*3],
							    w);
	      evaluateTestFunctionDerivativesOnElement_RT0_simplex_2D(&x_in[k*3],
								      &nodeArray[nN_0*3],
								      &nodeArray[nN_1*3],
								      &nodeArray[nN_2*3],
								      grad_w);
	    }


	  register double v[2],grad_v[2][2],
	    vx = 0.0, normv = 0.0;

	  for (int I = 0; I < nSpace; I++)
	    {
	      v[I] = 0.0; 
	      for (int i = 0; i < nDOF_element_velocity; i++)
		{
		  v[I] += w[i][I]*cvelocity_dof[cvelocity_l2g[eN*nDOF_element_velocity+i]];
		}
	      normv += v[I]*v[I];
	      for (int J = 0; J < nSpace; J++)
		{
		  grad_v[I][J] = 0.0;
		  for (int i = 0; i < nDOF_element_velocity; i++)
		    {
		      grad_v[I][J] += grad_w[i][I][J]*cvelocity_dof[cvelocity_l2g[eN*nDOF_element_velocity+i]];
		    }
		}
	    }
	  normv = sqrt(normv); 
	  assert(fabs(grad_v[0][0]-grad_v[1][1]) < zeroVelocityTol);
	  vx = grad_v[0][0]; 
	  //mwf debug
	  if (debugLevel > 0)
	    std::cout<<"trackPointsRT0Velocity2d k= "<<k<<" eN= "<<eN<<" x= ["<<x[0]<<","<<x[1]<<"], initial v= ["<<v[0]<<","<<v[1]<<"] grad_v[0][0]= "<<grad_v[0][0]<<std::endl;
	  //TODO refine this test since flag[k] >= 0 should be true always at this point
	  bool trackPoint = flag[k] >= -1 && fabs(normv) > zeroVelocityTol;
	  if (!trackPoint)
	    {
	      x_element[k] = eN;

	      x_arrive_times[k] = tTarget;
	      x_out[k*3+0] = x_in[k*3+0]; x_out[k*3+1] = x_in[k*3+1]; x_out[k*3+2] = x_in[k*3+2];
	    }
	  else
	    {
	      bool done = false;
	      int exitedDomain = 0;
	      while (!done)
		{
		  double dt_target = tTarget-t;
		  double dt = dt_target; //by default take full step
		  //mwf debug
		  if (debugLevel > 0)
		    std::cout<<"trackPointsRT0Velocity2d eN= "<<eN<<" k= "<<k<<" tTarget= "<<tTarget<<" t= "<<t<<" dt_target= "<<dt_target<<std::endl;
		  done = true;
		  int eN_exit = eN;
		  /***************************************************
	        in 2d, put check to make sure going into element
	          before calculating exit times because could have exit time along correct time axis
	          with a face but with an intersection that's outside the element if the point
	          never actually tracked into the element

                point leaving element only if 

                  |(\vec \bar{x}_f - \vec x_0)\cdot n_f| < \epsilon and
                     \vec v_0 \cdot \vec n_f > 0 
                  for some face f

                TODO: determine what can happen if v.n_f = 0, should have
                      intersection with another face that fits usual paradigm
                TODO: get tolerances straight on v.n test
		  ***************************************************/
		  bool enteredElement = true;
		  //do not reevaluate velocity if point "exits" an element without entering
		  // that is the point is on the boundary of the element (e.g., a node) but v.n > 0
                  bool evaluateVelocity = true; 
		  for (int ebN = 0; ebN < nElementBoundaries_element; ebN++)
		    {
		      const int ebN_global = elementBoundariesArray[eN*nElementBoundaries_element+ebN];
		      double dxf = 0.0, vdotn = 0.0;
		      for (int I = 0; I < nSpace; I++)
			{
			  dxf  += (elementBoundaryBarycentersArray[ebN_global*3 + I] - x[I])*elementLocalBoundaryOuterNormalsArray[eN*nElementBoundaries_element*nSpace + ebN*nSpace + I];
			  vdotn+= v[I]*elementLocalBoundaryOuterNormalsArray[eN*nElementBoundaries_element*nSpace + ebN*nSpace + I];
			}
		      if (dxf < -zeroIncrementTol) //outside element
			{
			  std::cout<<"PROBLEM trackPointsRT0Velocity2d done loop element check eN= "<<eN<<" x= ["<<x[0]<<","<<x[1]<<","<<x[2]
				   <<"] ebN= "<<ebN<<" ebN_global= "<<ebN_global<<" n= ["
				   <<elementLocalBoundaryOuterNormalsArray[eN*nElementBoundaries_element*nSpace + ebN*nSpace + 0]
				   <<","
				   <<elementLocalBoundaryOuterNormalsArray[eN*nElementBoundaries_element*nSpace + ebN*nSpace + 1]
				   <<"] dxf= "<<dxf<<std::endl;
			  std::cout<<"global Nodes on element : "
				   <<nN_0<<" --> ["<<nodeArray[nN_0*3+0]<<","<<nodeArray[nN_0*3+1]<<","<<nodeArray[nN_0*3+2]<<"] "
				   <<nN_1<<" --> ["<<nodeArray[nN_1*3+0]<<","<<nodeArray[nN_1*3+1]<<","<<nodeArray[nN_1*3+2]<<"] "
				   <<nN_2<<" --> ["<<nodeArray[nN_2*3+0]<<","<<nodeArray[nN_2*3+1]<<","<<nodeArray[nN_2*3+2]<<"] "
				   <<std::endl;
			  //move point to boundary
			  if (dxf > -insideElementTol)
			    {
			      for (int I = 0; I < nSpace; I++)
				x[I] += dxf*elementLocalBoundaryOuterNormalsArray[eN*nElementBoundaries_element*nSpace + ebN*nSpace + I];
			      //double check
			      dxf = 0.0;
			      for (int I = 0; I < nSpace; I++)
				{
				  //note switching order of x-x_f in plane definition to match convention for \Delta x in 1d and Haselbacher JCP paper
				  dxf += (elementBoundaryBarycentersArray[ebN_global*3 + I] - x[I])*elementLocalBoundaryOuterNormalsArray[eN*nElementBoundaries_element*nSpace + ebN*nSpace + I];
				}
			    }//move to element
			  assert(dxf >= -zeroIncrementTol);
			}
		      else if (fabs(dxf) < zeroIncrementTol && vdotn*dir > zeroIncrementTol) 
			{
			  enteredElement  = false;
			  evaluateVelocity= false; 
			  eN_exit = elementNeighborsArray[eN*nElementBoundaries_element + ebN];
			  dt = 0.0;
			  if (eN_exit < 0)
			    exitedDomain = 1;
			  done = exitedDomain; //only done if left the physical domain
			  //mwf debug
			  if (debugLevel > 1)
			    std::cout<<"trackPointsRT0Velocity2d on boundary did not enter domain eN= "<<eN<<" ebN= "<<ebN<<" dxf= "<<dxf
				     <<" vdotn= "<<vdotn<<" eN_exit= "<<eN_exit<<std::endl;
			}
		    }
		  if (enteredElement)
		    {
		      register double t_bnd[3],dt_bnd[3],dxf[3],vdotn[3];
		      //exit times to boundaries assuming |vx| = 0
		      for (int ebN = 0; ebN < nElementBoundaries_element; ebN++)
			{
			  const int ebN_global = elementBoundariesArray[eN*nElementBoundaries_element + ebN];
			  dxf[ebN] = 0.0, vdotn[ebN] = 0.0;
			  for (int I = 0; I < nSpace; I++)
			    {
			      dxf[ebN]  += (elementBoundaryBarycentersArray[ebN_global*3 + I] - x[I])*elementLocalBoundaryOuterNormalsArray[eN*nElementBoundaries_element*nSpace + ebN*nSpace + I];
			      vdotn[ebN]+= v[I]*elementLocalBoundaryOuterNormalsArray[eN*nElementBoundaries_element*nSpace + ebN*nSpace + I];
			    }
			  //by default assume no way to hit boundary
			  //set to invalid t exit
			  t_bnd[ebN] = t - LARGE*dir;
			  if (fabs(vdotn[ebN]) > zeroVelocityTol)
			    {
			      t_bnd[ebN] = dxf[ebN]/vdotn[ebN] + t;
			    }
			  //mwf debug
			  if (debugLevel > 2)
			    std::cout<<"eN= "<<eN<<" ebN= "<<ebN<<" ebN_global= "<<ebN_global<<" x= ["<<x[0]<<","<<x[1]<<"] dxf= "<<dxf[ebN]<<" vdotn= "<<vdotn[ebN]<<" t= "<<t<<" t_bnd= "<<t_bnd[ebN]<<std::endl;
			  //set dt to boundary to a big value so won't get chosen in min step
			  dt_bnd[ebN] = LARGE*dir;
			}
		      //account for linear variation
		      if (fabs(vx) > zeroVelocityTol)
			{
			  for (int ebN = 0; ebN < nElementBoundaries_element; ebN++)
			    {
			      if (fabs(vdotn[ebN]) > zeroVelocityTol)
				{
				  const double alpha = dxf[ebN]/vdotn[ebN];
				  if (alpha*vx <= -1.0) //no way to hit boundary, set to invalid time
				    t_bnd[ebN] = t - LARGE*dir;
				  else
				    t_bnd[ebN] = t + log(alpha*vx + 1.0)/vx;
				}
			    }
			}
		      //final exit time and boundary 
		      double dt_exit = LARGE*dir; int exitBoundary = -1;
		      //determine valid exit time steps after calculating exit times
		      for (int ebN = 0; ebN < nElementBoundaries_element; ebN++)
			{
			  if ((t_bnd[ebN]-t)*dir > zeroIncrementTol) //possibly went out ebN with nonzero step
			    {
			      dt_bnd[ebN] = t_bnd[ebN]-t;
			      //see if this valid exit time is the smallest one
			      if (fabs(dt_bnd[ebN]) < fabs(dt_exit))
				{
				  dt_exit = dt_bnd[ebN];
				  exitBoundary = ebN;
				}
			    }
			}
		      //mwf debug
		      if (debugLevel > 0)
			std::cout<<"after calculating exit times eN = "<<eN<<" dt_exit= "<<dt_exit<<" exitBoundary= "<<exitBoundary<<" dt_bnd= ["<<dt_bnd[0]<<","<<dt_bnd[1]<<","<<dt_bnd[2]<<"] "<<std::endl;

		      if (dt_exit*dir > zeroIncrementTol && fabs(dt_exit) < fabs(dt) - zeroIncrementTol) //exited current element before reaching target time
			{
			  assert(exitBoundary >= 0.0 && exitBoundary < nElementBoundaries_element);
			  dt = dt_exit;
			  int neigIndex = exitBoundary;
			  eN_exit = elementNeighborsArray[eN*nElementBoundaries_element + neigIndex];
			  if (eN_exit < 0)
			    exitedDomain = 1;
			  done = exitedDomain; //only done if left the physical domain
			  //mwf debug
			  if (debugLevel > 0)
			    {
			      std::cout<<"eN= "<<eN<<" x= ["<<x[0]<<","<<x[1]<<"] leaving through boundary "<<exitBoundary<<" globalBoundary= "<<elementBoundariesArray[eN*nElementBoundaries_element+exitBoundary]<<" to neighbor "<<eN_exit<<" v= ["<<v[0]<<","<<v[1]<<"]= vx= "<<vx
				       <<" dt= "<<dt<<" evaluateVelocity= "<<evaluateVelocity<<std::endl;
			      std::cout<<"eN elementBoundaries= ["<<elementBoundariesArray[eN*nElementBoundaries_element+0]
				       <<","<<elementBoundariesArray[eN*nElementBoundaries_element+1]
				       <<","<<elementBoundariesArray[eN*nElementBoundaries_element+2]<<"] "
				       <<"eN neighbors= ["<<elementNeighborsArray[eN*nElementBoundaries_element+0]
				       <<","<<elementNeighborsArray[eN*nElementBoundaries_element+1]
				       <<","<<elementNeighborsArray[eN*nElementBoundaries_element+2]<<"] "
				       <<std::endl;
			      if (eN_exit >= 0)
				{
				  std::cout<<"eN_exit elementBoundaries= ["<<elementBoundariesArray[eN_exit*nElementBoundaries_element+0]
					   <<","<<elementBoundariesArray[eN_exit*nElementBoundaries_element+1]
					   <<","<<elementBoundariesArray[eN_exit*nElementBoundaries_element+2]<<"] "
					   <<"eN_exit neighbors= ["<<elementNeighborsArray[eN_exit*nElementBoundaries_element+0]
					   <<","<<elementNeighborsArray[eN_exit*nElementBoundaries_element+1]
					   <<","<<elementNeighborsArray[eN_exit*nElementBoundaries_element+2]<<"] "
					   <<std::endl;
				}
			    }
			}
		    }//end entered element
		  //now update point and time
		  t += dt;
		  if (fabs(vx) < zeroVelocityTol)
		    {
		      for (int I = 0; I < nSpace; I++)
			x[I] += dt*v[I];
		    }
		  else
		    {
		      //mwf debug
		      if (debugLevel > 1)
			std::cout<<"updating x with linear formula x= ["<<x[0]<<","<<x[1]<<"] v= ["<<v[0]<<","<<v[1]<<"]= vx= "<<vx<<" dt= "<<dt;
		      for (int I = 0; I < nSpace; I++)
			x[I] += v[I]*(exp(vx*dt)-1.0)/vx;
		      //mwf debug
		      if (debugLevel > 1)
			std::cout<<" --> ["<<x[0]<<","<<x[1]<<"]"<<std::endl;
		    }
		  //update trajectory
		  n_traj++;
		  t_list.push_back(t);
		  for (int I=0; I < 3; I++)
		    x_list.push_back(x[I]);
		  e_list.push_back(eN);
		  
		  if (!done)
		    {
		      eN = eN_exit;

		      nN_0 = elementNodesArray[eN*nNodes_element + 0]; nN_1 = elementNodesArray[eN*nNodes_element + 1]; nN_2 = elementNodesArray[eN*nNodes_element + 2];
		      if (localVelocityRepresentationFlag == 1)
			{
			  evaluateTestFunctionsOnElement_RT0_simplex_2D_axbRep(x,
									       &nodeArray[nN_0*3],
									       &nodeArray[nN_1*3],
									       &nodeArray[nN_2*3],
									       w);
			  evaluateTestFunctionDerivativesOnElement_RT0_simplex_2D_axbRep(x,
											 &nodeArray[nN_0*3],
											 &nodeArray[nN_1*3],
											 &nodeArray[nN_2*3],
											 grad_w);

			}
		      else
			{
			  evaluateTestFunctionsOnElement_RT0_simplex_2D(x,
									&nodeArray[nN_0*3],
									&nodeArray[nN_1*3],
									&nodeArray[nN_2*3],
									w);
			  evaluateTestFunctionDerivativesOnElement_RT0_simplex_2D(x,
										  &nodeArray[nN_0*3],
										  &nodeArray[nN_1*3],
										  &nodeArray[nN_2*3],
										  grad_w);
			}

		      //evaluate velocity at current location
		      if (evaluateVelocity)
			{
			  normv = 0.0;
			  for (int I = 0; I < nSpace; I++)
			    {
			      v[I] = 0.0; 
			      for (int i = 0; i < nDOF_element_velocity; i++)
				{
				  v[I] += w[i][I]*cvelocity_dof[cvelocity_l2g[eN*nDOF_element_velocity+i]];
				}
			      normv += v[I]*v[I];
			      for (int J = 0; J < nSpace; J++)
				{
				  grad_v[I][J] = 0.0;
				  for (int i = 0; i < nDOF_element_velocity; i++)
				    {
				      //mwf debug
				      if (debugLevel > 2)
					std::cout<<"done loop evaluting velocity eN= "<<eN<<" x= ["<<x[0]<<","<<x[1]<<"], t= "<<t<<" grad_w["<<i<<"]["<<I<<"]["<<J<<"]= "<<grad_w[i][I][J]
						 <<" cvelocity_dof= "<<cvelocity_dof[cvelocity_l2g[eN*nDOF_element_velocity+i]]<<std::endl;
				      grad_v[I][J] += grad_w[i][I][J]*cvelocity_dof[cvelocity_l2g[eN*nDOF_element_velocity+i]];
				    }
				}
			    }
			  normv = sqrt(normv); 
			  //take advantage of RT0 velocity properties to simplify
			  assert(fabs(grad_v[0][0]-grad_v[1][1]) < zeroVelocityTol);
			  vx = grad_v[0][0]; 
			}//need to evaluate velocity
		      //mwf debug
		      if (debugLevel > 0 )
			std::cout<<"trackPointsRT0Velocity2d k= "<<k<<" eN= "<<eN<<" x= ["<<x[0]<<","<<x[1]<<"], t= "<<t<<" v= ["<<v[0]<<","<<v[1]<<"] grad_v[0][0]= "<<grad_v[0][0]<<std::endl;
		    }
		}//done

	      x_arrive_times[k] = t;
	      x_out[k*3+0] = x[0]; x_out[k*3+1] = x[1]; x_out[k*3+2] = x[2];
	      flag[k] = -1;
	      if (exitedDomain)
		flag[k] = -2;
	      x_element[k] = eN;
	    }//need to track point
	}//need to track point
    }//points to track

  //return arrays with the trajectory information
  //get one past the end
  offsets_list.push_back(n_traj);
  if (x_traj)
    delete [] x_traj;
  if (t_traj)
    delete [] t_traj;
  if (e_traj)
    delete [] e_traj;
  if (offsets_traj)
    delete [] offsets_traj;
  //mwf debug
  std::cout<<"starting trajectory copy n_traj= "<<n_traj<<" t_list.size() = "<<t_list.size()
	   <<" x_list.size()= "<<x_list.size()<<" offsets_list.size()= "<<offsets_list.size()<<std::endl;
  assert(unsigned(n_traj) == t_list.size());
  assert(unsigned(3*n_traj) == x_list.size());
  assert(unsigned(n_tracked+1) == offsets_list.size());
  e_traj = new int[n_traj];
  t_traj = new double[n_traj];
  x_traj = new double[3*n_traj];
  offsets_traj = new int[n_tracked+1];
  list<int>::const_iterator e_iter = e_list.begin();
  for (int i=0 ; e_iter != e_list.end(); i++,e_iter++)
      e_traj[i] = *e_iter;
  list<double>::const_iterator t_iter = t_list.begin();
  for (int i=0;  t_iter != t_list.end(); i++,t_iter++)
    t_traj[i] = *t_iter;
  list<double>::const_iterator x_iter = x_list.begin();
  for (int i=0; x_iter != x_list.end(); i++,x_iter++)
    x_traj[i] = *x_iter;
  list<int>::const_iterator o_iter = offsets_list.begin();
  for (int i=0 ; o_iter != offsets_list.end(); i++,o_iter++)
      offsets_traj[i] = *o_iter;
  std::cout<<"done with trajectory copy"<<std::endl;
}//simple element to element tracking

/***********************************************************************
   utility routines 
 ***********************************************************************/
 /***********************************************************************
   tag external boundary nodes
 ***********************************************************************/
extern "C" void setNodeOnBoundaryArray(int nExteriorElementBoundaries_global,
				       int nNodes_elementBoundary,
				       const int * exteriorElementBoundariesArray,
				       const int * elementBoundaryNodesArray,
				       int * nodeOnBoundaryArray)
{
  //assumes nodeOnBoundaryArray initialized correctly
  for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++)
    {
      int ebN = exteriorElementBoundariesArray[ebNE];
      for (int nN_local = 0; nN_local < nNodes_elementBoundary; nN_local++)
	{
	  int nN = elementBoundaryNodesArray[ebN*nNodes_elementBoundary+nN_local];
	  nodeOnBoundaryArray[nN] = 1;
	}
    }
}
/***********************************************************************
   for generating outer normals directly from mesh
 ***********************************************************************/
//assume constant jacobian over element so can just reuse element quadrature info
void getOuterNormals_affineSimplex_1d(int nElements_global,
				      int nElementBoundaries_element,
				      int nQuadraturePoints_element,
				      const double* boundaryNormals,
				      const double* jacobianInverseArray,
				      double* unitNormalArray)
{
  int eN,ebN,k;
  const int X=0;
  double *n;
  for(eN=0;eN<nElements_global;eN++)
    for(ebN=0;ebN<nElementBoundaries_element;ebN++)
      {
	k = 0;//just need 1 quadrature point since affine
	n = unitNormalArray+
	  eN*nElementBoundaries_element*1+
	  ebN;
	n[X] = boundaryNormals[ebN];
      }
}
void getOuterNormals_affineSimplex_2d(int nElements_global,
				      int nElementBoundaries_element,
				      int nQuadraturePoints_element,
				      const double* boundaryNormals,
				      const double* jacobianInverseArray,
				      double* unitNormalArray)
{
  int eN,ebN,k;
  const int X=0,Y=1,
    XX=0,XY=1,
    YX=2,YY=3;
  const double *jacInv=NULL, *bn=NULL;
  double *n;
  register double oneOverNbn=0.0;
  for(eN=0;eN<nElements_global;eN++)
    for(ebN=0;ebN<nElementBoundaries_element;ebN++)
      {
	k = 0;//just need 1 quadrature point since affine
	jacInv = jacobianInverseArray +
	  eN*nQuadraturePoints_element*4+
	  k*4;
	n = unitNormalArray+
	  eN*nElementBoundaries_element*2+
	  ebN*2;
	bn = boundaryNormals + ebN*2;
	n[X] = jacInv[XX]*bn[X]+jacInv[YX]*bn[Y];
	n[Y] = jacInv[XY]*bn[X]+jacInv[YY]*bn[Y];
	oneOverNbn= 1.0/sqrt(n[X]*n[X]+n[Y]*n[Y]);
	n[X] *= oneOverNbn;
	n[Y] *= oneOverNbn;
      }
}
void getOuterNormals_affineSimplex_3d(int nElements_global,
				      int nElementBoundaries_element,
				      int nQuadraturePoints_element,
				      const double* boundaryNormals,
				      const double* jacobianInverseArray,
				      double* unitNormalArray)
{
  int eN,ebN,k;
  const int X=0,Y=1,Z=2,
    XX=0,XY=1,XZ=2,
    YX=3,YY=4,YZ=5,
    ZX=6,ZY=7,ZZ=8;

  const double *jacInv=NULL, *bn=NULL;
  double *n;
  register double oneOverNbn=0.0;
  for(eN=0;eN<nElements_global;eN++)
    for(ebN=0;ebN<nElementBoundaries_element;ebN++)
      {
	k = 0;//just need 1 quadrature point since affine
	jacInv = jacobianInverseArray +
	  eN*nQuadraturePoints_element*9+
	  k*9;
	n = unitNormalArray+
	  eN*nElementBoundaries_element*3+
	  ebN*3;
	bn = boundaryNormals + ebN*3;
	n[X] = (jacInv[XX]*bn[X]+jacInv[YX]*bn[Y]+jacInv[ZX]*bn[Z]);
	n[Y] = (jacInv[XY]*bn[X]+jacInv[YY]*bn[Y]+jacInv[ZY]*bn[Z]);
	n[Z] = (jacInv[XZ]*bn[X]+jacInv[YZ]*bn[Y]+jacInv[ZZ]*bn[Z]);
          
	oneOverNbn = 1.0/sqrt(n[X]*n[X]+n[Y]*n[Y]+n[Z]*n[Z]);
        
	n[X] *= oneOverNbn;
	n[Y] *= oneOverNbn;
	n[Z] *= oneOverNbn;
 	
      }
}
//brute force search through elements to find which ones contain points
void findElementLocationsBruteForce(int nSpace,
				    int nElements_global,
				    int nElementBoundaries_element,
				    int nPoints,
				    const int* elementBoundariesArray,
				    const double* elementBoundaryBarycentersArray,
				    const double* elementBoundaryOuterNormalsArray,
				    const double* points,
				    int* element_locations)
{
  double tol = 1.0e-8;
  for (int i=0; i < nPoints; i++)
    {
      int eN = 0; bool foundElement = false;
      while (eN < nElements_global && !foundElement)
	{
	  bool outside = false;
	  for (int ebN = 0; ebN < nElementBoundaries_element; ebN++)
	    {
	      const int ebN_global = elementBoundariesArray[eN*nElementBoundaries_element+ebN];
	      double dxf =0.0;
	      for (int I = 0; I < nSpace; I++)
		{
		  dxf += (points[i*3 + I]-elementBoundaryBarycentersArray[ebN_global*3 + I])*elementBoundaryOuterNormalsArray[eN*nElementBoundaries_element*nSpace + ebN*nSpace + I];
		}
	      outside = outside || (dxf > tol);
	    }
	  if (!outside) foundElement = true;
	  eN++;
	}
      assert(foundElement);
      element_locations[i] = eN-1;
    }
}
			  

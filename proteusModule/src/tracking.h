#ifndef TRACKING_H
#define TRACKING_H
#include <cmath>
#include <cassert>
namespace Tracking
{
inline
void evaluateTestFunctionsOnElement_AffineLinear_1D(const double x[3], const double x0[3], const double x1[3], double w[2])
{
  double lambda_0 = (x[0]-x1[0])/(x0[0]-x1[0]);
  w[0] = lambda_0; w[1] = 1.0-lambda_0;
}
inline
void evaluateTestFunctionDerivativesOnElement_AffineLinear_1D(const double x[3], const double x0[3], const double x1[3], double grad_w[2])
{
  double dlambda_0 = 1.0/(x0[0]-x1[0]);
  grad_w[0] = dlambda_0; grad_w[1] = -dlambda_0;
}

inline
void evaluateTestFunctionsOnElement_RT0_simplex_2D(const double x[3], const double x0[3], const double x1[3], const double x2[3],
						   double w[3][2])
{
   /***********************************************************************
    Compute local projection of velocity to RT_0, where the local basis
    representation is

      \vec N_i = \frac{1}{d|\Omega_e|}(\vec x - x_{n,i}), i=0,...,d

    where x_{n,i} is the vertex across from face i, |\Omega| is the volume of the element,
     and d is the space dimension.

    The degrees of freedom are 
      V^f = \int_{\gamma_i}\vec v\dot n_{f}\ds
  ***************************************************/

 
  /***************************************************
   area 
       J    = |x_1-x_0  x_2-x_0|
              |y_1-y_0  y_2-y_0|

       detJ = (x_1-x_0)(y_2-y_0) - (y_1-y_0)(x_2-x_0)

  ***************************************************/
  using namespace std;
  const double dx10 = x1[0]-x0[0]; const double dx20 = x2[0]-x0[0];
  const double dy10 = x1[1]-x0[1]; const double dy20 = x2[1]-x0[1];
  const double det  = dx10*dy20-dy10*dx20; const double area = fabs(0.5*det);
  assert(area > 0.0);
  const double mgrad = 0.5/area;
  for (int I=0; I < 2; I++)
    {
      w[0][I] = mgrad*(x[I]-x0[I]);
      w[1][I] = mgrad*(x[I]-x1[I]);
      w[2][I] = mgrad*(x[I]-x2[I]);
    }
}

inline
void evaluateTestFunctionDerivativesOnElement_RT0_simplex_2D(const double x[3], const double x0[3], const double x1[3], const double x2[3],
							     double grad_w[3][2][2])
{
   /***********************************************************************
    Compute local projection of velocity to RT_0, where the local basis
    representation is

      \vec N_i = \frac{1}{d|\Omega_e|}(\vec x - x_{n,i}), i=0,...,d

    where x_{n,i} is the vertex across from face i, |\Omega| is the volume of the element,
     and d is the space dimension.

    The degrees of freedom are 
      V^f = \int_{\gamma_i}\vec v\dot n_{f}\ds
  ***************************************************/

 
  /***************************************************
   area 
       J    = |x_1-x_0  x_2-x_0|
              |y_1-y_0  y_2-y_0|

       detJ = (x_1-x_0)(y_2-y_0) - (y_1-y_0)(x_2-x_0)

  ***************************************************/
  using namespace std;
  const double dx10 = x1[0]-x0[0]; const double dx20 = x2[0]-x0[0];
  const double dy10 = x1[1]-x0[1]; const double dy20 = x2[1]-x0[1];
  const double det  = dx10*dy20-dy10*dx20; const double area = fabs(0.5*det);
  assert(area > 0.0);
  const double mgrad = 0.5/area;
  for (int i=0; i < 3; i++)
    {
      grad_w[i][0][0] = mgrad; grad_w[i][0][1] = 0.0; grad_w[i][1][0] = 0.0; grad_w[i][1][1] = mgrad;
    }
}

inline
void evaluateTestFunctionsOnElement_RT0_simplex_2D_axbRep(const double x[3], const double x0[3], const double x1[3], const double x2[3],
							  double w[3][2])
{
   /***********************************************************************
    Compute local projection of velocity to RT_0, where the local basis
    representation is

    \vec v_e = \vec a_e + b_e \vec x

    for simplicity, the basis numbering convention is
 
    \vec N_i = \vec e_i i=0,...,nd-1    e_0 = [1,0], e_1=[0,1]
    \vec N_{nd} =  \vec x
  ***************************************************/

  w[0][0] = 1.0; w[0][1] = 0.0;
  w[1][0] = 0.0; w[1][1] = 1.0;
  w[2][0] = x[0];w[2][1] = x[1];

}

inline
void evaluateTestFunctionDerivativesOnElement_RT0_simplex_2D_axbRep(const double x[3], const double x0[3], const double x1[3], const double x2[3],
								    double grad_w[3][2][2])
{
   /***********************************************************************
    Compute local projection of velocity to RT_0, where the local basis
    representation is

    \vec v_e = \vec a_e + b_e \vec x

    for simplicity, the basis numbering convention is
 
    \vec N_i = \vec e_i i=0,...,nd-1    e_0 = [1,0], e_1=[0,1]
    \vec N_{nd} =  \vec x
  ***************************************************/
  for (int i=0; i < 3; i++)
    for (int I = 0; I < 2; I++)
      for (int J = 0; J < 2; J++)
	grad_w[i][I][J] = 0.0;
  for (int I = 0; I < 2; I++)
    grad_w[2][I][I] = 1.0;
	  
}
}//Tracking namespace 


extern "C"
{
  void trackPointsConstVelocity1d(int nElements_global,               //mesh representation
				  int nNodes_global,
				  int nNodes_element,
				  int nElementBoundaries_element,
				  const double * nodeArray,
				  const int * elementNodesArray,
				  const int * elementNeighborsArray, //local boundary id is associated with node across from boundary 
				  double cVelocity,                  //characteristic speed (velocity) representation
				  double dir,                        //direction in time
				  const double *x_depart_times,      //in  -- time of departure for each point
				  double *x_arrive_times,            //in  -- desired stopping time for each point
				                                     //out -- stopping time for each point 
				  int nPointsToTrack,                //number of points to track
				  double zeroTolForTracking,         //ignore point if |u| < eps or |v| < eps
				  const double* x_in,                //points for tracking (always 3d)
				  int *x_element,                    //in -- element where point i is located at tIn
				                                     //out -- element where point i is located at tOut
				  double * x_out,                    //stopping location for point i at tOut
				  int * flag);                       //in: > -2  -- track point
                                                                     //in: -3    -- skip point
                                                                     //out: -1   -- point in interior at tOut
					                             //out: -2   -- point exited domain somewhere in (tIn,tOut)
                                                                     //out: -3   -- did not track (e.g.,  or u = 0)



  void trackPointsC0P1Velocity1d(int nElements_global,               //mesh representation
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
				 int * flag);                      //in: > -2  -- track point
                                                                   //in: -3    -- skip point
                                                                   //out: -1   -- point in interior at tOut
 				                                   //out: -2   -- point exited domain somewhere in (tIn,tOut)
                                                                   //out: -3   -- did not track (e.g.,  or u = 0)


void trackPointsRT0Velocity2d(int debugLevel,
			      int localVelocityRepresentationFlag, //2 -- dofs are fluxes through faces, 1 -- \vec a + b \vec x
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
			                                      //v = \vec a_e + b_e \vec x_e locally, 
			                                      //local dofs [\vec a_e,b_e]
			      double dir,                     //direction in time
			      const double *x_depart_times,   //in  -- time of departure for each point
			      double *x_arrive_times,         //in  -- desired stopping time
			                                      //out -- stopping time for each point 
			      int nPointsToTrack,             //number of points to track
			      double zeroTolForTracking,      //ignore point if |u| < eps or |v| < eps
			      const double* x_in,             //points for tracking (always 3d)
			      int *x_element,                 //in -- element where point i is located at tIn
			                                      //out -- element where point i is located at tOut
			      double * x_out,                 //stopping location for point i at tOut
			      int * flag);                      //in: > -2  -- track point
                                                                //in: -3    -- skip point
                                                                //out: -1   -- point in interior at tOut
 				                                //out: -2   -- point exited domain somewhere in (tIn,tOut)
                                                                //out: -3   -- did not track (e.g.,  or u = 0)

void trackPointsRT0Velocity2dWithTrajectories(int debugLevel,
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
   				         int *& e_traj);
/***********************************************************************
   utility routines
   
 ***********************************************************************/
void setNodeOnBoundaryArray(int nExteriorElementBoundaries_global,
			    int nNodes_elementBoundary,
			    const int * exteriorElementBoundariesArray,
			    const int * elementBoundaryNodesArray,
			    int * nodeOnBoundaryArray);
/***********************************************************************
  for generating outer normals directly from mesh
***********************************************************************/
//assume constant jacobian over element so can just reuse element quadrature info
void getOuterNormals_affineSimplex_1d(int nElements_global,
				      int nElementBoundaries_element,
				      int nQuadraturePoints_element,
				      const double* boundaryNormals,
				      const double* jacobianInverseArray,
				      double* unitNormalArray);
void getOuterNormals_affineSimplex_2d(int nElements_global,
				      int nElementBoundaries_element,
				      int nQuadraturePoints_element,
				      const double* boundaryNormals,
				      const double* jacobianInverseArray,
				      double* unitNormalArray);
void getOuterNormals_affineSimplex_3d(int nElements_global,
				      int nElementBoundaries_element,
				      int nQuadraturePoints_element,
				      const double* boundaryNormals,
				      const double* jacobianInverseArray,
				      double* unitNormalArray);

}//extern "C"
#endif

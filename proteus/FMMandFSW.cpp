#include "FMMandFSW.h"
#include "mesh.h"
#include "stupidheap.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cassert>

//----------------------------------------------------------------------
//local solvers
//----------------------------------------------------------------------

UnstructuredLocalUpwindSolvers::UnstructuredLocalUpwindSolvers():
  UNINITIALIZED(1.0e24)
{
}
UnstructuredLocalUpwindSolvers::~UnstructuredLocalUpwindSolvers()
{
}

double UnstructuredLocalUpwindSolvers::solve1d(double* x_N, double* x_N0, double T_N0, 
					       double speed, int verbose)
{
  /***********************************************************************

    repeat simple update formula in 1d but assuming unstructured mesh
    data structures and allow calling code to upwind

    Assumes Status[N_0] == 1, or othewise acceptable 

  **********************************************************************/
  const double dx = std::fabs(x_N[0]-x_N0[0]); const double dxInv2 = 1.0/(dx*dx);
  double N_a= dxInv2, N_b= T_N0*dxInv2, N_c=dxInv2*T_N0*T_N0;
  N_b *= -2.0; N_c -= 1.0/speed;
  const double discr = N_b*N_b - 4.0*N_a*N_c;
  assert(discr >= 0.0);
  //take larger root
  const double T_N = (-N_b + std::sqrt(discr))/(2.0*N_a);
  return T_N;
}


double UnstructuredLocalUpwindSolvers::solve2d(const double * x_C, 
					       const double * x_A,
					       const double * x_B,
					       double T_A,
					       double T_B,
					       double speed,
					       int verbose)
{
  /***********************************************************************
   2d local solver from Qian Zhang etal SIAM paper version 2
   ***********************************************************************/
  using namespace std;
  register double AC0,AC1,BC0,BC1,AB0,AB1,a,b,c,
    gamma,beta,alpha,theta,f_C,h,H,T_C;
  const double pi = M_PI;//otherwise 3.141592653589793;
  //edges
  AC0 = x_C[0]-x_A[0]; AC1 = x_C[1]-x_A[1]; 
  BC0 = x_C[0]-x_B[0]; BC1 = x_C[1]-x_B[1];
  AB0 = x_B[0]-x_A[0]; AB1 = x_B[1]-x_A[1];
  //lengths
  a = sqrt(BC0*BC0 + BC1*BC1); //across from node A
  b = sqrt(AC0*AC0 + AC1*AC1); //across from node B
  c = sqrt(AB0*AB0 + AB1*AB1); //across from node C
  //angles
  gamma = acos((AC0*BC0+AC1*BC1)/(a*b)); //angle at node C rev sign on both
  beta  = acos((AC0*AB0+AC1*AB1)/(b*c)); //angle at node A
  alpha = acos((-BC0*AB0-BC1*AB1)/(a*c));//angle at node B

  f_C = 1.0/speed;
  T_C = std::min<double>(T_A+b*f_C,T_B+a*f_C);
  if (fabs(T_B-T_A) <= c*f_C)
    {
      theta = asin((T_B-T_A)/(c*f_C)); //check abs in asin is there a typo in paper?
      if ((std::max<double>(0.,alpha-0.5*pi) <= theta && theta <= 0.5*pi - beta) ||
	  (alpha-0.5*pi) <= theta && theta <= std::min<double>(0.,0.5*pi - beta))
	{
	  h      = a*sin(alpha-theta); H = b*sin(beta+theta);
	  T_C = 0.5*(h*f_C+T_B) + 0.5*(H*f_C+T_A); 
	}
    }
  return T_C;
  
}



double UnstructuredLocalUpwindSolvers::solve2din3d(const double * x_C, 
						   const double * x_A,
						   const double * x_B,
						   double T_A,
						   double T_B,
						   double speed,
						   int verbose)
{
  /***********************************************************************
   2d local solver from Qian Zhang etal SIAM paper version 2
   only differs from FMMEikonalSolver2d version in that it assumes 
   embedded in 3d
   ***********************************************************************/
  using namespace std;
  register double AC0,AC1,AC2,BC0,BC1,BC2,AB0,AB1,AB2,a,b,c,
    gamma,beta,alpha,theta,f_C,h,H,T_C;
  const double pi = M_PI;//otherwise 3.141592653589793;
  //edges
  AC0 = x_C[0]-x_A[0]; AC1 = x_C[1]-x_A[1]; AC2 = x_C[2]-x_A[2]; 
  BC0 = x_C[0]-x_B[0]; BC1 = x_C[1]-x_B[1]; BC2 = x_C[2]-x_B[2];
  AB0 = x_B[0]-x_A[0]; AB1 = x_B[1]-x_A[1]; AB2 = x_B[2]-x_A[2];
  //lengths
  a = sqrt(BC0*BC0 + BC1*BC1 + BC2*BC2); //across from node A
  b = sqrt(AC0*AC0 + AC1*AC1 + AC2*AC2); //across from node B
  c = sqrt(AB0*AB0 + AB1*AB1 + AB2*AB2); //across from node C
  //angles
  gamma = acos((AC0*BC0+AC1*BC1+AC2*BC2)/(a*b)); //angle at node C rev sign on both
  beta  = acos((AC0*AB0+AC1*AB1+AC2*AB2)/(b*c)); //angle at node A
  alpha = acos((-BC0*AB0-BC1*AB1-BC2*AB2)/(a*c));//angle at node B

  f_C = 1.0/speed;
  T_C = std::min<double>(T_A+b*f_C,T_B+a*f_C);
  if (fabs(T_B-T_A) <= c*f_C)
    {
      theta = asin((T_B-T_A)/(c*f_C)); //check abs in asin is there a typo in paper?
      if ((std::max<double>(0.,alpha-0.5*pi) <= theta && theta <= 0.5*pi - beta) ||
	  (alpha-0.5*pi) <= theta && theta <= std::min<double>(0.,0.5*pi - beta))
	{
	  h      = a*sin(alpha-theta); H = b*sin(beta+theta);
	  T_C = 0.5*(h*f_C+T_B) + 0.5*(H*f_C+T_A); 
	}
    }
  return T_C;
  
}


double UnstructuredLocalUpwindSolvers::solve3d(const double * x_D,const double* x_1, 
					       const double* x_2, const double * x_3,
					       double T_1,  double T_2,  double T_3, 
					       double speed, int verbose,
					       double areaTol)
{
  /***********************************************************************
   3d local solver from Qian Zhang etal SIAM paper 
     
    Here N_D is the node for which we're solving. 
    We calculate N_A= argmin(T[N_1],T[N_2],T[N_3])
      and label the other two vertices as N_B and N_C

   ***********************************************************************/
  using namespace std;
  const double * x_A = x_1; const double * x_B = x_2; const double * x_C= x_3;
  register double T_A=T_1,T_B=T_2,T_C=T_3,f_D=1.0/speed,T_D;
  register double AB0,AB1,AB2,AC0,AC1,AC2,AD0,AD1,AD2,lAB,lAC;

  double ns[2][3]; 
  if (T_A > T_B)
    {
      swap(T_A,T_B);
      swap(x_A,x_B);
    }
  if (T_A > T_C)
    {
      swap(T_A,T_C);
      swap(x_A,x_C);
    }

  AB0 = x_B[0]-x_A[0]; AB1 = x_B[1]-x_A[1]; AB2 = x_B[2]-x_A[2];
  AC0 = x_C[0]-x_A[0]; AC1 = x_C[1]-x_A[1]; AC2 = x_C[2]-x_A[2];
  AD0 = x_D[0]-x_A[0]; AD1 = x_D[1]-x_A[1]; AD2 = x_D[2]-x_A[2];
  lAB = sqrt(AB0*AB0 + AB1*AB1 + AB2*AB2);
  lAC = sqrt(AC0*AC0 + AC1*AC1 + AC2*AC2);
  assert(lAB > 0.0); 
  assert(lAC > 0.0);

  T_D = UNINITIALIZED;
  bool T_DunInit=true;
  //branch 1
  if (T_B-T_A <= lAB*f_D && T_C-T_A <= lAC*f_D)
    {
      /*************************************************************
        solve quadratic system
         | AB  | | n_0 |   | (T_B-T_A)/f_D |  
         |-----| |     |   |               |
         | AC  | | n_1 | = | (T_C-T_A)/f_D |
         |-----| |     |   |               |
         | n   | | n_2 |   |       1       |
        
        
        write as
        
         [N] {n_01} = ({b} -  n_2 {a}) and
        
         {n_01}.{n_01} + n^2 = 1
        
         n_01 = [n_0, n_1]^T; b = [(T_B-T_A)/f_D, (T_C-T_A)/f_D]^T; a = [AB_2, AC_2]^T
        
               | AB_0 AB_1 |
         N=    | AC_0 AC_1 |
        
        solve for n_01 in terms of n_2, 
        
         {n_01} = [N]^(-1)({b} - n_2 {a})
        
        and get quadratic equation for n_2
        
         alpha n_2^2 + beta n_2 + gamma = 0
        
         alpha = C_00*a_0*a_0 + C_11*a_1*a_1 + 2*C_01*a_0*a_1 + 1
        
         beta  = -2(C_00*a_0*b_0 + C_11*a_1*b_1 + C_01*(a_1*b_0 + a_0*b_1))
        
         gamma = C_00 b_00*b_00 + C_11*b_11*b_11 + 2*C_01*b_0*b_1 - 1
        
        Here 
         
          C =  (NN^T)^(-1), call NN^T = B get B_{ij} = \sum^{1}_{k=0} N_{ik}N_{jk}
              |            |
          C = | B_11 -B_01 |           1
              |            | -------------------
              |-B_10  B_00 | (B_00B_11 - B_01B_10)

      *************************************************************/
      int ids[3] = {0,1,2};
      int permflag= 0;
      double a0  = AB2, a1 = AC2, b0 = (T_B-T_A)/f_D, b1 = (T_C-T_A)/f_D;
      double N00 = AB0, N01= AB1, N10= AC0, N11 = AC1;
      double detN = N00*N11-N01*N10;
      
      if (fabs(detN) < 1.0e-12)
	{
	  ids[0]=1; ids[1]=2; ids[2]=0;
	  permflag=1;
          //a0 = AB[ids[2]]; a1= AC[ids[2]]; b0=(T[N_B]-T[N_A])/f_D; b1=(T[N_C]-T[N_A])/f_D
	  //N00= AB[ids[0]]; N01= AB[ids[1]]; N10= AC[ids[0]]; N11= AC[ids[1]];
	  a0  = AB0;  a1 = AC0; 
	  N00 = AB1;  N01= AB2; N10 = AC1; N11 = AC2;
	  detN = N00*N11-N01*N10;
	}
      if (fabs(detN) < 1.0e-12)
	{
	  ids[0]=2; ids[1]=0; ids[2]=1;
	  permflag=2;
          //a0 = AB[ids[2]]; a1= AC[ids[2]]; b0=(T[N_B]-T[N_A])/f_D; b1=(T[N_C]-T[N_A])/f_D
	  //N00= AB[ids[0]]; N01= AB[ids[1]]; N10= AC[ids[0]]; N11= AC[ids[1]];
	  a0 = AB1; a1 = AC1; 
	  N00= AB2; N01= AB0; N10=AC2; N11=AC0;
	  detN = N00*N11 - N01*N10;
	}
      assert(fabs(detN) > 1.0e-12);
      const double Ni00=N11/detN, Ni01= -N01/detN, Ni10 = -N10/detN, Ni11 = N00/detN;
      const double B00 = N00*N00 + N01*N01, B01 = N00*N10 + N01*N11, B10 = N10*N00 + N11*N01,
	B11 = N10*N10 + N11*N11;
      const double detB = B00*B11-B01*B10;
      assert(fabs(detB) > 0.0);
      
      const double C00 = B11/detB, C01 = -B01/detB, C10 = -B10/detB, C11 = B00/detB;
      //C should be symmetric
      assert(fabs(C01-C10) < 1.0e-8);
      const double alpha = C00*a0*a0 + C11*a1*a1 + 2.0*C01*a0*a1 + 1.0;
      const double beta  =-2.0*(C00*a0*b0 + C11*a1*b1 + C01*(a1*b0 + a0*b1));
      const double gamma = C00*b0*b0 + C11*b1*b1 + 2.0*C01*b0*b1 - 1.0;

      const double disc  = beta*beta - 4.0*alpha*gamma;

      bool nsFound = false;
      if (disc >= 0.0)
	{
	  const double ns0ids2 = (-beta + sqrt(disc))/(2.0*alpha);
	  ns[0][ids[2]] = ns0ids2;
	  ns[0][ids[0]] = Ni00*(b0-a0*ns0ids2)+Ni01*(b1-a1*ns0ids2);
	  ns[0][ids[1]] = Ni10*(b0-a0*ns0ids2)+Ni11*(b1-a1*ns0ids2);
	  //
	  const double ns1ids2 = (-beta - sqrt(disc))/(2.*alpha);
	  ns[1][ids[2]] = ns1ids2;
	  ns[1][ids[0]] = Ni00*(b0-a0*ns1ids2)+Ni01*(b1-a1*ns1ids2);
	  ns[1][ids[1]] = Ni10*(b0-a0*ns1ids2)+Ni11*(b1-a1*ns1ids2);
	  nsFound = true;
	}
      if (nsFound == true)
	{
	  for (int i = 0; i < 2; i++)
	    {
	      /***************************************************
                solve for intersection point for ray comming out of D (back) 
                in direction n with
                plane ABC, call this E
                call x_E = \vec x_D - mu \vec n_{ABC}, mu > 0,
                should have (x_E - x_A) . \vec n_{ABC} = 0 

	       ***************************************************/
	      double x_E0 = 0.0; double x_E1=0.0; double x_E2=0.0;
	      bool Efound   = false;
	      double ABCn0(0.),ABCn1(0.),ABCn2(0.);

	      cross3d(AC0,AC1,AC2,
		      AB0,AB1,AB2,
		      ABCn0,ABCn1,ABCn2);
	      
	      const double nDotABCn = ns[i][0]*ABCn0 + ns[i][1]*ABCn1 + ns[i][2]*ABCn2;
	      if (fabs(nDotABCn) > 0.0)
		{
		  const double mu = (AD0*ABCn0 + AD1*ABCn1 + AD2*ABCn2)/nDotABCn;
		  if (mu > 0.)
		    {
		      x_E0 = x_D[0]-mu*ns[i][0]; x_E1 = x_D[1]-mu*ns[i][1]; x_E2 = x_D[2]-mu*ns[i][2];
		      Efound = true;
		    }//mu > 0
		}//normal not zero
	      if (Efound)
		{
		  const double areaEAB = 0.5*(parArea3d(x_E0-x_A[0],x_E1-x_A[1],x_E2-x_A[2],
							AB0,AB1,AB2));
		  const double areaEAC = 0.5*(parArea3d(x_E0-x_A[0],x_E1-x_A[1],x_E2-x_A[2],
							AC0,AC1,AC2));
		  const double areaEBC = 0.5*(parArea3d(x_E0-x_B[0],x_E1-x_B[1],x_E2-x_B[2],
							x_C[0]-x_B[0],x_C[1]-x_B[1],x_C[2]-x_B[2]));
		  const double areaABC = 0.5*(parArea3d(AC0,AC1,AC2,
							AB0,AB1,AB2));
		  if (fabs(areaEAB+areaEAC+areaEBC-areaABC) < areaTol)
		    {
		      const double T_Db1 = T_A + f_D*fabs(AD0*ns[i][0]+AD1*ns[i][1]+AD2*ns[i][2]);
		      if (T_DunInit)
			{T_D = T_Db1; T_DunInit = false;}
		      else
			T_D  = std::min<double>(T_D,T_Db1);
		    }
		  else
		    {
		      //solve 2d update on tetrahedron surfaces
		      const double T_ABD = solve2din3d(x_D,x_A,x_B,T_A,T_B,speed,verbose);
		      const double T_ACD = solve2din3d(x_D,x_A,x_C,T_A,T_C,speed,verbose);
		      const double T_BCD = solve2din3d(x_D,x_B,x_C,T_B,T_C,speed,verbose);
		      T_D = min<double>(min<double>(T_ABD,T_ACD),T_BCD); 
		      T_DunInit = false;
		    }
		}//end intersection found
	      else
		{
		  const double T_ABD = solve2din3d(x_D,x_A,x_B,T_A,T_B,speed,verbose);
		  const double T_ACD = solve2din3d(x_D,x_A,x_C,T_A,T_C,speed,verbose);
		  const double T_BCD = solve2din3d(x_D,x_B,x_C,T_B,T_C,speed,verbose);
		  T_D = min<double>(min<double>(T_ABD,T_ACD),T_BCD); 
		  T_DunInit = false;
		}
	    }//2 wave front normals 
	}//found a wave front normal
      else
	{
	  //solve on tetrahedron surfaces
	  const double T_ABD = solve2din3d(x_D,x_A,x_B,T_A,T_B,speed,verbose);
	  const double T_ACD = solve2din3d(x_D,x_A,x_C,T_A,T_C,speed,verbose);
	  const double T_BCD = solve2din3d(x_D,x_B,x_C,T_B,T_C,speed,verbose);
	  T_D = min<double>(min<double>(T_ABD,T_ACD),T_BCD); 
	  T_DunInit = false;
	}
    }//branch 1
  else
    {
      //solve on tetrahedron surfaces
      const double T_ABD = solve2din3d(x_D,x_A,x_B,T_A,T_B,speed,verbose);
      const double T_ACD = solve2din3d(x_D,x_A,x_C,T_A,T_C,speed,verbose);
      const double T_BCD = solve2din3d(x_D,x_B,x_C,T_B,T_C,speed,verbose);
      T_D = min<double>(min<double>(T_ABD,T_ACD),T_BCD); 
      T_DunInit = false;
    }
  assert(!T_DunInit);
  return T_D;
}

//----------------------------------------------------------------------

FMMEikonalSolverBase::FMMEikonalSolverBase(Mesh* meshIn,int nSpaceIn,
					   FMMEikonalSolverBase::INIT_TYPE initIn,
					   bool forcePositiveInitialValuesIn): 
  mesh(meshIn),Status(),Known(),nSpace(nSpaceIn),localSolvers(),
  UNINITIALIZED(1.0e24),initFlag(initIn),
  forcePositiveInitialValues(forcePositiveInitialValuesIn)
{
  if (mesh)
    {
      Status.resize(mesh->nNodes_global,FAR);
      Known.reserve(mesh->nNodes_global);
    }  
}

FMMEikonalSolverBase::~FMMEikonalSolverBase()
{
}
bool FMMEikonalSolverBase::initializeKnownPoints(const double* phi0, double * T,
						 double zeroTol, int verbose)
{
  if (initFlag == FRONT_AND_MAGNITUDE)
    return initializeKnownPointsAlaWeakDirBCs(phi0,T,zeroTol,verbose);
  else
    return initializeKnownPointsUsingMagnitude(phi0,T,zeroTol,verbose);
}
bool FMMEikonalSolverBase::initializeKnownPointsUsingMagnitude(const double* phi0, double * T,
							       double zeroTol, int verbose)
{
  using namespace std;
  //mwf debug
  if (verbose > 5)
    {
      std::cout<<"initializeMag zeroTol= "<<zeroTol<<" nNodes_global= "<<mesh->nNodes_global
	       <<" Status.size()= "<<Status.size()<<std::endl;
    }
  if (forcePositiveInitialValues)
    {
      for (int I = 0; I < mesh->nNodes_global; I++)
	T[I] = fabs(phi0[I]);
    }
  else
    copy(phi0,phi0+mesh->nNodes_global,T);
  fill(Status.begin(),Status.end(),FAR);
  Known.clear();

  for (int I = 0; I < mesh->nNodes_global; I++)
    {
      if (fabs(phi0[I]) <= zeroTol)
	{
	  Status[I] = KNOWN;
	  Known.push_back(I);
	  if (verbose > 5)
	    std::cout<<"FMM init magnitude known phi["<<I<<"]= "<<phi0[I]<<std::endl;
	}
      else
	{
	  Status[I]= FAR;
	}
    }
  bool failed = Known.size() < unsigned(nSpace);  //for now need at least nSpace known points to do local solves
  return failed;
}

bool FMMEikonalSolverBase::initializeKnownPointsUsingFrontIntersection(const double* phi0, double * T,
								       double zeroTol, int verbose)
{
  using namespace std;

  if (forcePositiveInitialValues)
    {
      for (int I = 0; I < mesh->nNodes_global; I++)
	T[I] = fabs(phi0[I]);
    }
  else
    copy(phi0,phi0+mesh->nNodes_global,T);
  fill(Status.begin(),Status.end(),FAR);
  Known.clear();

  for (int I = 0; I < mesh->nNodes_global; I++)
    {
      if (fabs(phi0[I]) <= zeroTol)
	{
	  Status[I] = KNOWN;
	  Known.push_back(I);
	  if (verbose > 5)
	    std::cout<<"FMM init frontMagnitude known phi["<<I<<"]= "<<phi0[I]<<" from mag"<<std::endl;
	}
      else
	{
	  bool frontFound = false;
	  int startI = mesh->nodeElementOffsets[I];
	  int nElements_nodeI = mesh->nodeElementOffsets[I+1]-startI; 
	  for (int i = 0; i < nElements_nodeI; i++)
	    {
	      const int eN = mesh->nodeElementsArray[startI+i];
	      for (int nN_local = 0; nN_local < mesh->nNodes_element; nN_local++)
		{
		  const int nN = mesh->elementNodesArray[eN*mesh->nNodes_element+nN_local];
		  if (nN != I)
		    {
		      frontFound = frontFound || (phi0[I]*phi0[nN] <= 0.0);
		      if (verbose > 10)
			std::cout<<"initFrontMag I= "<<I<<" startI= "<<startI<<" nElements_nodeI= "
				 <<nElements_nodeI<<" eN= "<<eN<<" nN_local= "<<nN_local<<" nN= "<<nN<<" frontFound = "
				 <<frontFound<<std::endl;
		    }
		}
	    }
	  if (frontFound)
	    {
	      Status[I] = KNOWN;
	      Known.push_back(I);
	      if (verbose >5)
		{
		  std::cout<<"initFrontMag I= "<<I<<" frontFound= "<<frontFound<<" phi0["<<I<<"]= "
			   <<phi0[I]<<std::endl;
		}
	    }
	
	}
    }
  bool failed = Known.size() < unsigned(nSpace); //for now need at least nSpace known points to do local solves
  return failed;

}

bool FMMEikonalSolverBase::initializeKnownPointsAlaWeakDirBCs(const double* phi0, double * T,
							      double zeroTol, int verbose)
{
  using namespace std;

  if (forcePositiveInitialValues)
    {
      for (int I = 0; I < mesh->nNodes_global; I++)
	T[I] = fabs(phi0[I]);
    }
  else
    copy(phi0,phi0+mesh->nNodes_global,T);
  fill(Status.begin(),Status.end(),FAR);
  Known.clear();


  for (int eN = 0; eN < mesh->nElements_global; eN++)
    {
      int    signU = 0;
      int j0       = 0.0;
      double eps   = 0.0;
      while (signU == 0 && j0 < mesh->nNodes_element)
	{
	  int J0 = mesh->elementNodesArray[eN*mesh->nNodes_element+j0];//same as l2g
	  if (phi0[J0] < -eps)
	    {
	      signU = -1;
	    }
	  else if (phi0[J0] > eps)
	    {
	      signU = 1;
	    }
	  else
	    {
	      if (Status[J0] != KNOWN)
		Known.push_back(J0);
	      Status[J0] = KNOWN;
	    }
	  j0 += 1;
	}//while j0
      for (int j = j0; j < mesh->nNodes_element; j++)
	{
	  int J = mesh->elementNodesArray[eN*mesh->nNodes_element+j];
	  if ((phi0[J] < -eps && signU == 1) || (phi0[J] > eps && signU == -1))
	    {
	      for (int jj = 0; jj < mesh->nNodes_element; jj++)
		{
		  int JJ = mesh->elementNodesArray[eN*mesh->nNodes_element+jj];
		  if (Status[JJ] != KNOWN)
		    Known.push_back(JJ);
		  Status[JJ] = KNOWN;
		}
	      break;
	    }
	  else if (fabs(phi0[J]) < eps)
	    {
	      if (Status[J] != KNOWN)
		Known.push_back(J);
	      Status[J] = KNOWN;
	    }
	}//j
    }//eN

  bool failed = Known.size() < unsigned(nSpace); //for now need at least nSpace known points to do local solves
  return failed;

}

bool FMMEikonalSolverBase::solve(const double* phi0, 
				 const double * nodalSpeeds,
				 double * T, double zeroTol, 
				 double trialTol,
				 int initFlagIn,// -1 -->ignore, 0 --> magn. 1--> frontInt 
				 int verbose)
{
  /***********************************************************************
        Test first order fast marching method algorithm for eikonal equation
        
        \|\grad T \| = 1, \phi(\vec x) = 0, x \in \Gamma
        
        assuming \phi_0 describes initial location of interface Gamma and
        has reasonable values (absolute values) for T close to Gamma. Here
        T can be interpreted as the travel time from Gamma.
        
        Right now assumes global node numbers <--> global dofs but this can be
        fixed easily
        Input
        

        phi0: dof array from P1 C0 FiniteElementFunction holding initial condition
        
        T   : dof array from P1 C0 FiniteElementFunction for solution
        

        Output
        T(\vec x_n)    : travel time from initial front to node (\vec x_n) 

        Internal data structures

        Status : status of nodal point (dictionary)
             -1 --> Far
              0 --> Trial
              1 --> Known
        Trial : nodal points adjacent to front tuples (index,val) stored in heap

   ***********************************************************************/
  using namespace std;
  bool failed = false;
  int nWastedRevisit = 0;
  //take opportunity to reset some flags?
  if (initFlagIn == 0)//MAGNITUDE ONLY
    { initFlag = MAGNITUDE;}
  else if (initFlagIn == 1)
    { initFlag = FRONT_AND_MAGNITUDE;}
  //otherwise don't change value constructed with
 
  //check sizes
  if (Status.size() != unsigned(mesh->nNodes_global))
    {
      //c0 p1 required
      if (verbose > 2)
	cout<<"FMM solve resizing Status from "<<Status.size()<<" to "<<unsigned(mesh->nNodes_global)
	    <<endl;
      Status.resize(mesh->nNodes_global,FAR);
    }
  StupidHeap Trial;
  //mwf hack
  //std::cout<<"In FMMandFSW baseSolve calling init"<<std::endl;
  failed = initializeKnownPoints(phi0,T,zeroTol,verbose);
  if (failed)
    {
      cout<<"initialize failed quitting"<<endl;
      return failed;
    }
  //mwf hack
  //std::cout<<"In FMMandFSW baseSolve back from calling init"<<std::endl;
  assert(Known.size() > 0);
  //loop through known points and calculate travel time for their unKNOWN neighbors
  for (vector<int>::iterator pI = Known.begin(); pI != Known.end(); pI++)
    {
      const int I = *pI;
      const int startI = mesh->nodeElementOffsets[I];
      int nElements_starI = mesh->nodeElementOffsets[I+1]-startI;
      for (int eN_star = 0; eN_star < nElements_starI; eN_star++)
	{
	  const int eN = mesh->nodeElementsArray[startI+eN_star];
	  for (int nN_local = 0; nN_local < mesh->nNodes_element; nN_local++)
	    {
	      const int N_0 = mesh->elementNodesArray[eN*mesh->nNodes_element+nN_local];
	      if (N_0 != I && Status[N_0] != KNOWN)
		{
		  double T_N = UNINITIALIZED;
		  bool updateFailed = localUpdate(N_0,T_N,T,nodalSpeeds[N_0],verbose);
		  if (!updateFailed)
		    {
		      T[N_0] = T_N;
		      if (Status[N_0] == TRIAL)
			Trial.updateNode(N_0,T[N_0]);
		      else
			Trial.insert(N_0,T[N_0]);
		      Status[N_0] = TRIAL;
		    }
		}//neighbor that's UNKNOWN
	    }//nodes on element 
	}//elements in node star
    }//Known nodes

  while (!Trial.isEmpty())
    {
      StupidHeap::EntryType minNode = Trial.pop();
      const int I = minNode.first; const double T_I = minNode.second;
      if (verbose > 1)
	{
	  cout<<"*** FMM Start Trial loop min Trial point = ("<<I<<","<<T_I<<") Status= "
	      <<Status[I]<<endl;
	}
      if (Status[I] != KNOWN)
	{
	  Status[I] = KNOWN; Known.push_back(I);
	  T[I] = T_I;
	}
      else
	{
	  nWastedRevisit++;
	}
      //now update travel times of I's unKNOWN neighbors
      const int startI = mesh->nodeElementOffsets[I];
      const int nElements_starI = mesh->nodeElementOffsets[I+1]-startI;
      for (int eN_star = 0; eN_star < nElements_starI; eN_star++)
	{
	  const int eN = mesh->nodeElementsArray[startI+eN_star];
	  for (int nN_local = 0; nN_local < mesh->nNodes_element; nN_local++)
	    {
	      const int N_0 = mesh->elementNodesArray[eN*mesh->nNodes_element+nN_local];
	      if (N_0 != I && Status[N_0] != KNOWN)
		{
		  double T_N;
		  bool updateFailed = localUpdate(N_0,T_N,T,nodalSpeeds[N_0],verbose);
		  if (!updateFailed)
		    {
		      T[N_0] = T_N;
		      if (Status[N_0] == TRIAL)
			Trial.updateNode(N_0,T[N_0]);
		      else
			Trial.insert(N_0,T[N_0]);
		      Status[N_0] = TRIAL;
		    }
		}//neighbor that's unKNOWN
	    }//local nodes on neighboring element
	}//elements in node star
    }//Trial empty
  //sanity check
  if (verbose > 0)
    {
      if (Known.size() != unsigned(mesh->nNodes_global))
	{
	  cout<<"Problem at end of FMM 2d sovle Known.size() = "<<Known.size()
	      <<" nNodes_global= "<<mesh->nNodes_global<<endl;
	}
      for (unsigned int I = 0; I < Status.size(); I++)
	{
	  if (Status[I] != KNOWN)
	    {
	      cout<<"Problem TRIAL or FAR point found in Status["<<I<<"]= "<<Status[I]<<endl;
	    }
	}
      cout<<"FMM nWastedRevisit= "<<nWastedRevisit<<endl;
    }
  return failed;
}


//----------------------------------------------------------------------
//1d
//----------------------------------------------------------------------
FMMEikonalSolver1d::FMMEikonalSolver1d(Mesh* meshIn,
				       FMMEikonalSolverBase::INIT_TYPE initIn): 
  FMMEikonalSolverBase(meshIn,1,initIn)
{
}

FMMEikonalSolver1d::~FMMEikonalSolver1d()
{
}


bool FMMEikonalSolver1d::localUpdate(int N, double& T_N, const double * T, 
				     double speed, int verbose)
{
  /***********************************************************************
    Loop through all element neighbors of N and compute Travel time using
     KNOWN nodal values on that element
    Takes min of all possible updates

    returns True if no update possible
   ***********************************************************************/
  const int startN         = mesh->nodeElementOffsets[N]; 
  const int nElements_node = mesh->nodeElementOffsets[N+1]-startN;
  T_N = UNINITIALIZED;
  bool failure = true;
  for (int eN_star = 0; eN_star < nElements_node; eN_star++)
    {
      const int eN = mesh->nodeElementsArray[startN+eN_star];
      //need nSpace Known neighbors to update
      for (int nN_local = 0; nN_local < mesh->nNodes_element; nN_local++)
	{
	  const int N_0 = mesh->elementNodesArray[eN*mesh->nNodes_element+nN_local];
	  if (N_0 != N && Status[N_0] == KNOWN) //KNOWN neighbor found
	    {
	      const double T_eN_N = localSolvers.solve1d(mesh->nodeArray + N*3,mesh->nodeArray + N_0*3,T[N_0],
							 speed,verbose);
	      if (T_N == UNINITIALIZED)
		T_N = T_eN_N;
	      else 
		T_N = std::min(T_N,T_eN_N);
	      failure = false;
	    }
	}//node neighbors
    }//element neighbors
  return failure;
}

//----------------------------------------------------------------------
//2d
//----------------------------------------------------------------------
FMMEikonalSolver2d::FMMEikonalSolver2d(Mesh* meshIn,
				       FMMEikonalSolverBase::INIT_TYPE initIn): 
  FMMEikonalSolverBase(meshIn,2,initIn)
{
}

FMMEikonalSolver2d::~FMMEikonalSolver2d()
{
}


bool FMMEikonalSolver2d::localUpdate(int N, double& T_N, const double * T, 
				     double speed, int verbose)
{
  /***********************************************************************
    Loop through all element neighbors of N and compute Travel time using
     KNOWN nodal values on that element
    Takes min of all possible updates

    returns True if no update possible
   ***********************************************************************/
  const int startN         = mesh->nodeElementOffsets[N]; 
  const int nElements_node = mesh->nodeElementOffsets[N+1]-startN;
  int knownNeigs[2];
  T_N = UNINITIALIZED;
  bool failure = true;
  for (int eN_star = 0; eN_star < nElements_node; eN_star++)
    {
      const int eN = mesh->nodeElementsArray[startN+eN_star];
      int nneig = 0;
      knownNeigs[0] = -1; knownNeigs[1] = -1;//assume not enough known neigs
      //need nSpace Known neighbors to update
      for (int nN_local = 0; nN_local < mesh->nNodes_element; nN_local++)
	{
	  const int N_0 = mesh->elementNodesArray[eN*mesh->nNodes_element+nN_local];
	  if (N_0 != N && Status[N_0] == KNOWN) //KNOWN neighbor found
	    {
	      knownNeigs[nneig] = N_0;
	      nneig++;
	    }
	}//node neighbors
      if (nneig == 2)
	{
	  const int N_A = knownNeigs[0]; const int N_B = knownNeigs[1];
	  assert(N_A >= 0); assert(N_B >= 0);
	  const double T_eN_N = localSolvers.solve2d(mesh->nodeArray + N*3,
						     mesh->nodeArray + N_A*3,
						     mesh->nodeArray + N_B*3,
						     T[N_A],T[N_B],
						     speed,verbose);
	  if (T_N == UNINITIALIZED)
	    T_N = T_eN_N;
	  else 
	    T_N = std::min(T_N,T_eN_N);
	  failure = false;
	}//2 known neighbors found
    }//element neighbors
  return failure;
}
//----------------------------------------------------------------------
//3d
//----------------------------------------------------------------------
FMMEikonalSolver3d::FMMEikonalSolver3d(Mesh* meshIn,
				       FMMEikonalSolverBase::INIT_TYPE initIn): 
  FMMEikonalSolverBase(meshIn,3,initIn)
{
}

FMMEikonalSolver3d::~FMMEikonalSolver3d()
{
}


bool FMMEikonalSolver3d::localUpdate(int N, double& T_N, const double * T, 
				     double speed, int verbose)
{
  /***********************************************************************
    Loop through all element neighbors of N and compute Travel time using
     KNOWN nodal values on that element
    Takes min of all possible updates

    returns True if no update possible
   ***********************************************************************/
  const int startN         = mesh->nodeElementOffsets[N]; 
  const int nElements_node = mesh->nodeElementOffsets[N+1]-startN;
  int knownNeigs[3];
  T_N = UNINITIALIZED;
  bool failure = true;
  for (int eN_star = 0; eN_star < nElements_node; eN_star++)
    {
      const int eN = mesh->nodeElementsArray[startN+eN_star];
      int nneig = 0;
      knownNeigs[0] = -1; knownNeigs[1] = -1; knownNeigs[2]=-1;//assume not enough known neigs
      //need nSpace Known neighbors to update
      for (int nN_local = 0; nN_local < mesh->nNodes_element; nN_local++)
	{
	  const int N_0 = mesh->elementNodesArray[eN*mesh->nNodes_element+nN_local];
	  if (N_0 != N && Status[N_0] == KNOWN) //KNOWN neighbor found
	    {
	      knownNeigs[nneig] = N_0;
	      nneig++;
	    }
	}//node neighbors
      if (nneig == 3)
	{
	  const int N_A = knownNeigs[0]; const int N_B = knownNeigs[1]; const int N_C = knownNeigs[2];
	  assert(N_A >= 0); assert(N_B >= 0); assert(N_C >= 0);
	  const double T_eN_N = localSolvers.solve3d(mesh->nodeArray + N*3,
						     mesh->nodeArray + N_A*3,
						     mesh->nodeArray + N_B*3,
						     mesh->nodeArray + N_C*3,
						     T[N_A],T[N_B],T[N_C],
						     speed,verbose);
	  if (T_N == UNINITIALIZED)
	    T_N = T_eN_N;
	  else 
	    T_N = std::min(T_N,T_eN_N);
	  failure = false;
	}//2 known neighbors found
    }//element neighbors
  return failure;
}

//----------------------------------------------------------------------
//FSW routines
//----------------------------------------------------------------------

FSWEikonalSolverBase::FSWEikonalSolverBase(Mesh* meshIn,int nSpaceIn,
					   double atolIn, double rtolIn, int maxItsIn,
					   FMMEikonalSolverBase::INIT_TYPE initIn,
					   int nRefPointsIn,
					   const double* refPointsIn):
  FMMEikonalSolverBase(meshIn,nSpaceIn,initIn,true),//always force positive for now
  nRefPoints(nRefPointsIn),refPoints(refPointsIn,3*nRefPointsIn),Order(),
  iterAtol(atolIn),iterRtol(rtolIn),maxIts(maxItsIn),T0()
{
  bool failed(false);
  if (mesh)
    T0.resize(mesh->nNodes_global);
  failed = buildOrderings();
}

FSWEikonalSolverBase::~FSWEikonalSolverBase()
{
}

bool FSWEikonalSolverBase::buildOrderings()
{
  bool failed = false;
  if (!mesh)
    { failed = true; return failed; }
  //first set default points if none given
  assert(refPoints.size() == unsigned(3*nRefPoints));
  
  if (nRefPoints == 0)
    {
      if (nSpace == 1)
	{
	  //just brute force min/max rather than get cute with stl
	  double minX = mesh->nodeArray[0]; double maxX = minX;
	  int nN_minX(0),nN_maxX(0);
	  for (int nN = 0; nN < mesh->nNodes_global; nN++)
	    {
	      if (mesh->nodeArray[nN*3+0] < minX)
		{ minX = mesh->nodeArray[nN*3+0]; nN_minX = nN;}
	      if (mesh->nodeArray[nN*3+0] > maxX)
		{ maxX = mesh->nodeArray[nN*3+0]; nN_maxX = nN;}
	    }
	  nRefPoints = 2;
	  refPoints.resize(3*nRefPoints);
	  refPoints[3*0+0]= mesh->nodeArray[nN_minX*3+0];
	  refPoints[3*0+1]= mesh->nodeArray[nN_minX*3+1];
	  refPoints[3*0+2]= mesh->nodeArray[nN_minX*3+2];
	  //
	  refPoints[3*1+0]= mesh->nodeArray[nN_maxX*3+0];
	  refPoints[3*1+1]= mesh->nodeArray[nN_maxX*3+1];
	  refPoints[3*1+2]= mesh->nodeArray[nN_maxX*3+2];
	}//1d
      else if (nSpace == 2)
	{
	  //just brute force min/max rather than get cute with stl
	  double minX = mesh->nodeArray[0]; double maxX = minX;
	  double minY = mesh->nodeArray[1]; double maxY = minY;
	  int nN_minX(0),nN_maxX(0),nN_minY(0),nN_maxY(0);
	  for (int nN = 0; nN < mesh->nNodes_global; nN++)
	    {
	      if (mesh->nodeArray[nN*3+0] < minX)
		{ minX = mesh->nodeArray[nN*3+0]; nN_minX = nN;}
	      if (mesh->nodeArray[nN*3+0] > maxX)
		{ maxX = mesh->nodeArray[nN*3+0]; nN_maxX = nN;}
	      if (mesh->nodeArray[nN*3+1] < minY)
		{ minY = mesh->nodeArray[nN*3+1]; nN_minY = nN;}
	      if (mesh->nodeArray[nN*3+1] > maxY)
		{ maxY = mesh->nodeArray[nN*3+1]; nN_maxY = nN;}
	    }
	  nRefPoints = 4;
	  //(x_min,y_min), (x_max,y_min), (x_min,y_max), (x_max,y_max)
	  refPoints.resize(3*nRefPoints);
	  refPoints[3*0+0]= mesh->nodeArray[nN_minX*3+0];
	  refPoints[3*0+1]= mesh->nodeArray[nN_minY*3+1];
	  refPoints[3*0+2]= 0.0;//better option?
	  //
	  refPoints[3*1+0]= mesh->nodeArray[nN_maxX*3+0];
	  refPoints[3*1+1]= mesh->nodeArray[nN_minY*3+1];
	  refPoints[3*1+2]= 0.0;
	  //
	  refPoints[3*2+0]= mesh->nodeArray[nN_minX*3+0];
	  refPoints[3*2+1]= mesh->nodeArray[nN_maxY*3+1];
	  refPoints[3*2+2]= 0.0;//better option?
	  //
	  refPoints[3*3+0]= mesh->nodeArray[nN_maxX*3+0];
	  refPoints[3*3+1]= mesh->nodeArray[nN_maxY*3+1];
	  refPoints[3*3+2]= 0.0;
	}//2d
      else
	{
	  //just brute force min/max rather than get cute with stl
	  double minX = mesh->nodeArray[0]; double maxX = minX;
	  double minY = mesh->nodeArray[1]; double maxY = minY;
	  double minZ = mesh->nodeArray[2]; double maxZ = minZ;
	  int nN_minX(0),nN_maxX(0),nN_minY(0),nN_maxY(0),
	    nN_minZ(0),nN_maxZ(0);
	  for (int nN = 0; nN < mesh->nNodes_global; nN++)
	    {
	      if (mesh->nodeArray[nN*3+0] < minX)
		{ minX = mesh->nodeArray[nN*3+0]; nN_minX = nN;}
	      if (mesh->nodeArray[nN*3+0] > maxX)
		{ maxX = mesh->nodeArray[nN*3+0]; nN_maxX = nN;}
	      if (mesh->nodeArray[nN*3+1] < minY)
		{ minY = mesh->nodeArray[nN*3+1]; nN_minY = nN;}
	      if (mesh->nodeArray[nN*3+1] > maxY)
		{ maxY = mesh->nodeArray[nN*3+1]; nN_maxY = nN;}
	      if (mesh->nodeArray[nN*3+2] < minZ)
		{ minZ = mesh->nodeArray[nN*3+2]; nN_minZ = nN;}
	      if (mesh->nodeArray[nN*3+2] > maxZ)
		{ maxZ = mesh->nodeArray[nN*3+2]; nN_maxZ = nN;}
	    }
	  nRefPoints = 8;
	  //(x_min,y_min,z_min), (x_max,y_min,z_min), 
	  //(x_min,y_max,z_min), (x_max,y_max,z_min)
	  //(x_min,y_min,z_max), (x_max,y_min,z_max), 
	  //(x_min,y_max,z_max), (x_max,y_max,z_max)

	  refPoints.resize(3*nRefPoints);
	  refPoints[3*0+0]= mesh->nodeArray[nN_minX*3+0];
	  refPoints[3*0+1]= mesh->nodeArray[nN_minY*3+1];
	  refPoints[3*0+2]= mesh->nodeArray[nN_minZ*3+2];
	  //
	  refPoints[3*1+0]= mesh->nodeArray[nN_maxX*3+0];
	  refPoints[3*1+1]= mesh->nodeArray[nN_minY*3+1];
	  refPoints[3*1+2]= mesh->nodeArray[nN_minZ*3+2];
	  //
	  refPoints[3*2+0]= mesh->nodeArray[nN_minX*3+0];
	  refPoints[3*2+1]= mesh->nodeArray[nN_maxY*3+1];
	  refPoints[3*2+2]= mesh->nodeArray[nN_minZ*3+2];
	  //
	  refPoints[3*3+0]= mesh->nodeArray[nN_maxX*3+0];
	  refPoints[3*3+1]= mesh->nodeArray[nN_maxY*3+1];
	  refPoints[3*3+2]= mesh->nodeArray[nN_minZ*3+2];
	  //
	  refPoints[3*0+0]= mesh->nodeArray[nN_minX*3+0];
	  refPoints[3*0+1]= mesh->nodeArray[nN_minY*3+1];
	  refPoints[3*0+2]= mesh->nodeArray[nN_maxZ*3+2];
	  //
	  refPoints[3*1+0]= mesh->nodeArray[nN_maxX*3+0];
	  refPoints[3*1+1]= mesh->nodeArray[nN_minY*3+1];
	  refPoints[3*1+2]= mesh->nodeArray[nN_maxZ*3+2];
	  //
	  refPoints[3*2+0]= mesh->nodeArray[nN_minX*3+0];
	  refPoints[3*2+1]= mesh->nodeArray[nN_maxY*3+1];
	  refPoints[3*2+2]= mesh->nodeArray[nN_maxZ*3+2];
	  //
	  refPoints[3*3+0]= mesh->nodeArray[nN_maxX*3+0];
	  refPoints[3*3+1]= mesh->nodeArray[nN_maxY*3+1];
	  refPoints[3*3+2]= mesh->nodeArray[nN_maxZ*3+2];
	  
	}//3d
    }//built own reference points
  //now build orderings using l_1 distance
  Order.resize(nRefPoints);
  for (int i = 0; i < nRefPoints; i++)
    Order[i].reserve(mesh->nNodes_global);
  //sort according to min distance from point
  StupidHeap theap;
  for (int i = 0; i < nRefPoints; i++)
    {
      for (int nN = 0; nN < mesh->nNodes_global; nN++)
	{
	  double dlp1 = 0.0;
	  for (int id=0; id < nSpace; id++)
	    dlp1 += fabs(mesh->nodeArray[nN*3+id]-refPoints[i*3+id]);
	  theap.insert(nN,dlp1);
	}//nN
      while (! theap.isEmpty())
	{
	  StupidHeap::EntryType next = theap.pop();
	  Order[i].push_back(next.first);
	}//heap empty
      assert(Order[i].size() == unsigned(mesh->nNodes_global));
    }//i
  failed = false;
  return failed;
}

//-------------------- init routines --------------------
//only real difference is setting T values to UNINITIALIZED instead of phi0
//and always forces positive values
bool FSWEikonalSolverBase::initializeKnownPointsUsingMagnitude(const double* phi0, double * T,
							       double zeroTol, int verbose)
{
  using namespace std;
  //mwf debug
  if (verbose > 5)
    {
      std::cout<<"initializeMag zeroTol= "<<zeroTol<<" nNodes_global= "<<mesh->nNodes_global
	       <<" Status.size()= "<<Status.size()<<std::endl;
    }
  fill(T,T+mesh->nNodes_global,UNINITIALIZED);
  fill(Status.begin(),Status.end(),TRIAL);
  Known.clear();

  for (int I = 0; I < mesh->nNodes_global; I++)
    {
      if (fabs(phi0[I]) <= zeroTol)
	{
	  Status[I] = KNOWN;
	  Known.push_back(I);
	  T[I] = std::fabs(phi0[I]);
	  if (verbose > 5)
	    std::cout<<"FMM init magnitude known phi["<<I<<"]= "<<phi0[I]<<std::endl;
	}
    }
  bool failed = Known.size() < unsigned(nSpace);  //for now need at least nSpace known points to do local solves
  return failed;
}

bool FSWEikonalSolverBase::initializeKnownPointsUsingFrontIntersection(const double* phi0, double * T,
								       double zeroTol, int verbose)
{
  using namespace std;

  fill(T,T+mesh->nNodes_global,UNINITIALIZED);
  fill(Status.begin(),Status.end(),TRIAL);
  Known.clear();

  for (int I = 0; I < mesh->nNodes_global; I++)
    {
      if (fabs(phi0[I]) <= zeroTol)
	{
	  Status[I] = KNOWN;
	  Known.push_back(I);
	  T[I] = fabs(phi0[I]);
	  if (verbose > 5)
	    std::cout<<"FMM init frontMagnitude known phi["<<I<<"]= "<<phi0[I]<<" from mag"<<std::endl;
	}
      else
	{
	  bool frontFound = false;
	  int startI = mesh->nodeElementOffsets[I];
	  int nElements_nodeI = mesh->nodeElementOffsets[I+1]-startI; 
	  for (int i = 0; i < nElements_nodeI; i++)
	    {
	      const int eN = mesh->nodeElementsArray[startI+i];
	      for (int nN_local = 0; nN_local < mesh->nNodes_element; nN_local++)
		{
		  const int nN = mesh->elementNodesArray[eN*mesh->nNodes_element+nN_local];
		  if (nN != I)
		    {
		      frontFound = frontFound || (phi0[I]*phi0[nN] <= 0.0);
		      if (verbose > 10)
			std::cout<<"initFrontMag I= "<<I<<" startI= "<<startI<<" nElements_nodeI= "
				 <<nElements_nodeI<<" eN= "<<eN<<" nN_local= "<<nN_local<<" nN= "<<nN<<" frontFound = "
				 <<frontFound<<std::endl;
		    }
		}
	    }
	  if (frontFound)
	    {
	      Status[I] = KNOWN;
	      Known.push_back(I);
	      T[I] = fabs(phi0[I]);
	      if (verbose >5)
		{
		  std::cout<<"initFrontMag I= "<<I<<" frontFound= "<<frontFound<<" phi0["<<I<<"]= "
			   <<phi0[I]<<std::endl;
		}
	    }
	
	}
    }
  bool failed = Known.size() < unsigned(nSpace); //for now need at least nSpace known points to do local solves
  return failed;

}
bool FSWEikonalSolverBase::initializeKnownPointsAlaWeakDirBCs(const double* phi0, double * T,
							      double zeroTol, int verbose)
{
  using namespace std;

  fill(T,T+mesh->nNodes_global,UNINITIALIZED);
  fill(Status.begin(),Status.end(),TRIAL);
  Known.clear();


  for (int eN = 0; eN < mesh->nElements_global; eN++)
    {
      int    signU = 0;
      int j0       = 0.0;
      double eps   = 0.0;
      while (signU == 0 && j0 < mesh->nNodes_element)
	{
	  int J0 = mesh->elementNodesArray[eN*mesh->nNodes_element+j0];//same as l2g
	  if (phi0[J0] < -eps)
	    {
	      signU = -1;
	    }
	  else if (phi0[J0] > eps)
	    {
	      signU = 1;
	    }
	  else
	    {
	      if (Status[J0] != KNOWN)
		Known.push_back(J0);
	      Status[J0] = KNOWN;
	      T[J0] = fabs(phi0[J0]);
	    }
	  j0 += 1;
	}//while j0
      for (int j = j0; j < mesh->nNodes_element; j++)
	{
	  int J = mesh->elementNodesArray[eN*mesh->nNodes_element+j];
	  if ((phi0[J] < -eps && signU == 1) || (phi0[J] > eps && signU == -1))
	    {
	      for (int jj = 0; jj < mesh->nNodes_element; jj++)
		{
		  int JJ = mesh->elementNodesArray[eN*mesh->nNodes_element+jj];
		  if (Status[JJ] != KNOWN)
		    Known.push_back(JJ);
		  Status[JJ] = KNOWN;
		  T[JJ] = fabs(phi0[JJ]);
		}
	      break;
	    }
	  else if (fabs(phi0[J]) < eps)
	    {
	      if (Status[J] != KNOWN)
		Known.push_back(J);
	      Status[J] = KNOWN;
	      T[J] = fabs(phi0[J]);
	    }
	}//j
    }//eN

  bool failed = Known.size() < unsigned(nSpace); //for now need at least nSpace known points to do local solves
  return failed;

}

bool FSWEikonalSolverBase::solve(const double* phi0, const double * nodalSpeeds, double * T,
				 double zeroTol, double trialTol,
				 int initTypeFlag,// -1 -->ignore, 0 --> magn. 1--> frontInt 
				 int verbose)
{
  using namespace std;
  bool failed = false;
  assert(mesh); assert(phi0); assert(T); assert(nodalSpeeds);


  //take opportunity to reset some flags?
  if (initTypeFlag == 0)//MAGNITUDE ONLY
    { initFlag = MAGNITUDE;}
  else if (initTypeFlag == 1)
    { initFlag = FRONT_AND_MAGNITUDE;}
  //otherwise don't change value constructed with
 
  //check sizes
  if (Status.size() != unsigned(mesh->nNodes_global))
    {
      //c0 p1 required
      if (verbose > 2)
	cout<<"FSW solve resizing Status from "<<Status.size()<<" to "<<unsigned(mesh->nNodes_global)
	    <<endl;
      Status.resize(mesh->nNodes_global,FAR);
    }
  if (T0.size() != unsigned(mesh->nNodes_global))
    {
      if (verbose > 0)
	cout<<"FSW WARNING having to resize T0 from "<<T0.size()<<" to "<<mesh->nNodes_global<<endl;
      T0.resize(mesh->nNodes_global);
    }

  failed = initializeKnownPoints(phi0,T,zeroTol,verbose);
  if (failed)
    {
      cout<<"initialize failed quitting"<<endl;
      return failed;
    }
  assert(Known.size() > 0);
 
 
  double phi0norm = 0.0;
  //l2 for now
  for (int I = 0; I < mesh->nNodes_global; I++)
    phi0norm += phi0[I]*phi0[I];
  phi0norm /= mesh->nNodes_global;
  phi0norm  = sqrt(phi0norm);

  int nSweeps = 0, nIts = 0, maxSweeps = nRefPoints*maxIts;
  double err = 1.0e28;
  bool converged = false;
  while (!converged && nSweeps < maxSweeps)
    {
      if (verbose > 0)
	cout<<"FSW starting sweep #"<<nSweeps<<" err= "<<err<<endl;
      copy(T,T+mesh->nNodes_global,&T0[0]);
      int iorder = nSweeps % nRefPoints;
      for (vector<int>::iterator pN = Order[iorder].begin();
	   pN != Order[iorder].end(); pN++) //up
	{
	  if (Status[*pN] != KNOWN) //leave known points alone
	    {
	      double T_N;
	      bool T_N_failed = localUpdate(*pN,T_N,T,nodalSpeeds[*pN],
					    verbose);
	      if (!T_N_failed)
		T[*pN] = T_N;
	    }
	}//up
      for (vector<int>::reverse_iterator pN = Order[iorder].rbegin(); 
	   pN != Order[iorder].rend(); pN++) //down
	{

	  if (Status[*pN] != KNOWN) //leave known points alone
	    {
	      double T_N;
	      bool T_N_failed = localUpdate(*pN,T_N,T,nodalSpeeds[*pN],
					    verbose);
	      if (!T_N_failed)
		T[*pN] = T_N;
	    }
	}//down
      nSweeps++;
      double errSum = 0.0;
      //l2
      for (int nN = 0; nN < mesh->nNodes_global; nN++)
	errSum += (T0[nN]-T[nN])*(T0[nN]-T[nN]);
      errSum /= mesh->nNodes_global;
      err = sqrt(errSum);
      converged = err <= (iterAtol + iterRtol*phi0norm);
    }//sweeps
  nIts = nSweeps/nRefPoints;
  failed = failed || !converged;
  if (verbose > 0)
    cout<<"FSW leaving nIts= "<<nIts<<" err= "<<err<<" converged= "
	<<converged<<" failed= "<<failed<<endl;
  return failed;
}

//----------------------------------------------------------------------
//1d
//----------------------------------------------------------------------
FSWEikonalSolver1d::FSWEikonalSolver1d(Mesh* meshIn, 
				       double atolIn, double rtolIn,
				       int maxItsIn,
				       INIT_TYPE initIn,
				       int nRefPointsIn,
				       const double* refPointsIn): 
  FSWEikonalSolverBase(meshIn,1,atolIn,rtolIn,maxItsIn,initIn,
		       nRefPointsIn,refPointsIn)
{
}

FSWEikonalSolver1d::~FSWEikonalSolver1d()
{
}


bool FSWEikonalSolver1d::localUpdate(int N, double& T_N, const double * T, 
				     double speed, int verbose)
{
  /***********************************************************************
    Loop through all element neighbors of N and compute Travel time using
     nodal values on that element if they are less than T[N] since
     have a monotonic update formula
    Takes min of all possible updates

    returns True if no update possible
   ***********************************************************************/
  const int startN         = mesh->nodeElementOffsets[N]; 
  const int nElements_node = mesh->nodeElementOffsets[N+1]-startN;
  T_N = UNINITIALIZED;
  bool failure = true;
  for (int eN_star = 0; eN_star < nElements_node; eN_star++)
    {
      const int eN = mesh->nodeElementsArray[startN+eN_star];
      //need nSpace Known neighbors to update
      for (int nN_local = 0; nN_local < mesh->nNodes_element; nN_local++)
	{
	  const int N_0 = mesh->elementNodesArray[eN*mesh->nNodes_element+nN_local];
	  if (N_0 != N && T[N_0] < T[N]) //earlier neighbor found
	    {
	      const double T_eN_N = localSolvers.solve1d(mesh->nodeArray + N*3,mesh->nodeArray + N_0*3,T[N_0],
							 speed,verbose);
	      if (T_N == UNINITIALIZED)
		T_N = T_eN_N;
	      else 
		T_N = std::min(T_N,T_eN_N);
	      failure = false;
	    }
	}//node neighbors
    }//element neighbors
  return failure;
}

//----------------------------------------------------------------------
//2d
//----------------------------------------------------------------------
FSWEikonalSolver2d::FSWEikonalSolver2d(Mesh* meshIn, 
				       double atolIn, double rtolIn,
				       int maxItsIn,
				       INIT_TYPE initIn,
				       int nRefPointsIn,
				       const double* refPointsIn): 
  FSWEikonalSolverBase(meshIn,2,atolIn,rtolIn,maxItsIn,initIn,
		       nRefPointsIn,refPointsIn)
{
}

FSWEikonalSolver2d::~FSWEikonalSolver2d()
{
}

bool FSWEikonalSolver2d::localUpdate(int N, double& T_N, const double * T, 
				     double speed, int verbose)
{
  /***********************************************************************
    Loop through all element neighbors of N and compute Travel time using
     earlier nodal values on that element
    Takes min of all possible updates

    returns True if no update possible
   ***********************************************************************/
  const int startN         = mesh->nodeElementOffsets[N]; 
  const int nElements_node = mesh->nodeElementOffsets[N+1]-startN;
  int earlierNeigs[2];
  T_N = UNINITIALIZED;
  bool failure = true;
  for (int eN_star = 0; eN_star < nElements_node; eN_star++)
    {
      const int eN = mesh->nodeElementsArray[startN+eN_star];
      int nneig = 0;
      earlierNeigs[0] = -1; earlierNeigs[1] = -1;//assume not enough known neigs
      //need nSpace Known neighbors to update
      for (int nN_local = 0; nN_local < mesh->nNodes_element; nN_local++)
	{
	  const int N_0 = mesh->elementNodesArray[eN*mesh->nNodes_element+nN_local];
	  if (N_0 != N && T[N_0] < T[N]) //earlier neighbor found
	    {
	      earlierNeigs[nneig] = N_0;
	      nneig++;
	    }
	}//node neighbors
      if (nneig == 2)
	{
	  const int N_A = earlierNeigs[0]; const int N_B = earlierNeigs[1];
	  assert(N_A >= 0); assert(N_B >= 0);
	  const double T_eN_N = localSolvers.solve2d(mesh->nodeArray + N*3,
						     mesh->nodeArray + N_A*3,
						     mesh->nodeArray + N_B*3,
						     T[N_A],T[N_B],
						     speed,verbose);
	  if (T_N == UNINITIALIZED)
	    T_N = T_eN_N;
	  else 
	    T_N = std::min(T_N,T_eN_N);
	  failure = false;
	}//2 earlier neighbors found
    }//element neighbors
  return failure;
}

//----------------------------------------------------------------------
//3d
//----------------------------------------------------------------------
FSWEikonalSolver3d::FSWEikonalSolver3d(Mesh* meshIn, 
				       double atolIn, double rtolIn,
				       int maxItsIn,
				       INIT_TYPE initIn,
				       int nRefPointsIn,
				       const double* refPointsIn): 
  FSWEikonalSolverBase(meshIn,3,atolIn,rtolIn,maxItsIn,initIn,
		       nRefPointsIn,refPointsIn)
{
}

FSWEikonalSolver3d::~FSWEikonalSolver3d()
{
}

bool FSWEikonalSolver3d::localUpdate(int N, double& T_N, const double * T, 
				     double speed, int verbose)
{
  /***********************************************************************
    Loop through all element neighbors of N and compute Travel time using
     earlier nodal values on that element
    Takes min of all possible updates

    returns True if no update possible
   ***********************************************************************/
  const int startN         = mesh->nodeElementOffsets[N]; 
  const int nElements_node = mesh->nodeElementOffsets[N+1]-startN;
  int earlierNeigs[3];
  T_N = UNINITIALIZED;
  bool failure = true;
  for (int eN_star = 0; eN_star < nElements_node; eN_star++)
    {
      const int eN = mesh->nodeElementsArray[startN+eN_star];
      int nneig = 0;
      earlierNeigs[0] = -1; earlierNeigs[1] = -1; earlierNeigs[2]=-1;//assume not enough known neigs
      //need nSpace Known neighbors to update
      for (int nN_local = 0; nN_local < mesh->nNodes_element; nN_local++)
	{
	  const int N_0 = mesh->elementNodesArray[eN*mesh->nNodes_element+nN_local];
	  if (N_0 != N && T[N_0] < T[N]) //earlier neighbor found
	    {
	      earlierNeigs[nneig] = N_0;
	      nneig++;
	    }
	}//node neighbors
      if (nneig == 3)
	{
	  const int N_A = earlierNeigs[0]; const int N_B = earlierNeigs[1]; const int N_C = earlierNeigs[2];
	  assert(N_A >= 0); assert(N_B >= 0); assert(N_C >= 0);
	  const double T_eN_N = localSolvers.solve3d(mesh->nodeArray + N*3,
						     mesh->nodeArray + N_A*3,
						     mesh->nodeArray + N_B*3,
						     mesh->nodeArray + N_C*3,
						     T[N_A],T[N_B],T[N_C],
						     speed,verbose);
	  if (T_N == UNINITIALIZED)
	    T_N = T_eN_N;
	  else 
	    T_N = std::min(T_N,T_eN_N);
	  failure = false;
	}//2 known neighbors found
    }//element neighbors
  return failure;
}


//----------------------------------------------------------------------
//reconstruction routines
//----------------------------------------------------------------------
namespace NearFrontInitialization
{
bool localPWLreconstruction(int nSpace, Mesh* mesh, const double* phi0, double * phi0R, 
			    double zeroTol, int verbose)
{
  using namespace std;
  assert(mesh); assert(phi0); assert(phi0R); 
  bool failed = false;
  double x0nN[3][3], xIx0[3],dX10[3],xIx01[3];
  //by default just copy over
  copy(phi0,phi0+mesh->nNodes_global,phi0R);

  if (nSpace == 2)
    {
      for (int I = 0; I < mesh->nNodes_global; I++)
	{
	  if (fabs(phi0[I]) > zeroTol)
	    {
	      const double signI = (phi0[I] >= 0.0) ? 1.0 : -1.0;
	      bool reconValFound = false;
	      double reconVal = -12345.0; //nonsense value
	      //loop through neighboring elements see if find a zero intersection

	      const int startI = mesh->nodeElementOffsets[I];
	      const int nElements_node = mesh->nodeElementOffsets[I+1]-startI;
	      for (int eN_star = 0; eN_star < nElements_node; eN_star++)
		{
		  const int eN = mesh->nodeElementsArray[startI + eN_star];
		  int nint = 0;
		  for (int nN_local = 0; nN_local < mesh->nNodes_element; nN_local++)
		    {
		      const int nN = mesh->elementNodesArray[eN*mesh->nNodes_element+nN_local];
		      if (nN != I && phi0[nN]*phi0[I] <= 0.0)
			{
			  //at least intersect through 1 edge
			  const double s0 = -phi0[I]/(phi0[nN]-phi0[I]);
			  x0nN[nint][0] = mesh->nodeArray[I*3+0]*(1.0-s0)+mesh->nodeArray[nN*3+0]*s0;
			  x0nN[nint][1] = mesh->nodeArray[I*3+1]*(1.0-s0)+mesh->nodeArray[nN*3+1]*s0;
			  x0nN[nint][2] = mesh->nodeArray[I*3+2]*(1.0-s0)+mesh->nodeArray[nN*3+2]*s0;

			  xIx0[0] = mesh->nodeArray[I*3+0]-x0nN[nint][0];
			  xIx0[1] = mesh->nodeArray[I*3+1]-x0nN[nint][1];
			  xIx0[2] = mesh->nodeArray[I*3+2]-x0nN[nint][2];
			  nint++;
			  const double phiInN = signI*sqrt(xIx0[0]*xIx0[0]+xIx0[1]*xIx0[1]+xIx0[2]*xIx0[2]);
			  if (!reconValFound)
			    reconVal = phiInN;
			  else if (fabs(reconVal) > fabs(phiInN))
			    reconVal = phiInN;
			  reconValFound = true;
			    
			}//sign difference
		    }//local nodes
		  //check if intersect through 2 edges to get minimum distance
		  if (nint == 2)
		    {
		      dX10[0] = x0nN[1][0]-x0nN[0][0]; dX10[1] = x0nN[1][1]-x0nN[0][1]; 
		      dX10[2] = x0nN[1][2]-x0nN[0][2];
		      const double ndX10 = dX10[0]*dX10[0] + dX10[1]*dX10[1] + dX10[2]*dX10[2];
		      const double s01 = (dX10[0]*(mesh->nodeArray[I*3+0]-x0nN[0][0])+
					  dX10[1]*(mesh->nodeArray[I*3+1]-x0nN[0][1])+
					  dX10[2]*(mesh->nodeArray[I*3+2]-x0nN[0][2]))/ndX10;
		      if (0.0 < s01 && s01 < 1.0)
			{
			  //intersection in element
			  xIx01[0] = mesh->nodeArray[I*3+0]-(1.0-s01)*x0nN[0][0] - s01*x0nN[1][0];
			  xIx01[1] = mesh->nodeArray[I*3+1]-(1.0-s01)*x0nN[0][1] - s01*x0nN[1][1];
			  xIx01[2] = mesh->nodeArray[I*3+2]-(1.0-s01)*x0nN[0][2] - s01*x0nN[1][2];
			  const double phiIJ01 = signI*sqrt(xIx01[0]*xIx01[0] + xIx01[1]*xIx01[1] +
							    xIx01[2]*xIx01[2]);
			  assert(reconValFound);
			  if (fabs(phiIJ01) < reconVal)
			    reconVal = phiIJ01;
			}//intersection in interor
		    }//two intersections found on this element
		}//element neighbors
	      if (reconValFound)
		phi0R[I] = reconVal;
	      if (verbose > 3 && reconValFound)
		{
		  cout<<"c localPWL I= "<<I<<" reconValFound= "<<reconValFound
		      <<" phi0["<<I<<"]= "<<phi0[I]<<" phi0R["<<I<<"]= "<<phi0R[I]<<endl;
		}
	    } // > zeroTol
	}//global nodes
    }//2d
  else
    {
      for (int I = 0; I < mesh->nNodes_global; I++)
	{
	  if (fabs(phi0[I]) > zeroTol)
	    {
	      const double signI = (phi0[I] >= 0.0) ? 1.0 : -1.0;
	      bool reconValFound = false;
	      double reconVal = -12345.0; //nonsense value
	      //loop through neighboring elements see if find a zero intersection

	      const int startI = mesh->nodeElementOffsets[I];
	      const int nElements_node = mesh->nodeElementOffsets[I+1]-startI;
	      for (int eN_star = 0; eN_star < nElements_node; eN_star++)
		{
		  const int eN = mesh->nodeElementsArray[startI + eN_star];
		  int nint = 0;
		  for (int nN_local = 0; nN_local < mesh->nNodes_element; nN_local++)
		    {
		      const int nN = mesh->elementNodesArray[eN*mesh->nNodes_element+nN_local];
		      if (nN != I && phi0[nN]*phi0[I] <= 0.0)
			{
			  //at least intersect through 1 edge
			  const double s0 = -phi0[I]/(phi0[nN]-phi0[I]);
			  x0nN[nint][0] = mesh->nodeArray[I*3+0]*(1.0-s0)+mesh->nodeArray[nN*3+0]*s0;
			  x0nN[nint][1] = mesh->nodeArray[I*3+1]*(1.0-s0)+mesh->nodeArray[nN*3+1]*s0;
			  x0nN[nint][2] = mesh->nodeArray[I*3+2]*(1.0-s0)+mesh->nodeArray[nN*3+2]*s0;

			  xIx0[0] = mesh->nodeArray[I*3+0]-x0nN[nint][0];
			  xIx0[1] = mesh->nodeArray[I*3+1]-x0nN[nint][1];
			  xIx0[2] = mesh->nodeArray[I*3+2]-x0nN[nint][2];
			  nint++;
			  const double phiInN = signI*sqrt(xIx0[0]*xIx0[0]+xIx0[1]*xIx0[1]+xIx0[2]*xIx0[2]);
			  if (!reconValFound)
			    reconVal = phiInN;
			  else if (fabs(reconVal) > fabs(phiInN))
			    reconVal = phiInN;
			  reconValFound = true;
			    
			}//sign difference
		    }//local nodes
		}//element neighbors
	      if (reconValFound)
		phi0R[I] = reconVal;
	      if (verbose > 3 && reconValFound)
		{
		  cout<<"c localPWL I= "<<I<<" reconValFound= "<<reconValFound
		      <<" phi0["<<I<<"]= "<<phi0[I]<<" phi0R["<<I<<"]= "<<phi0R[I]<<endl;
		}
	    } // > zeroTol
	}//global nodes
    }
  return failed;
}

bool localDGPWLreconstruction(int nSpace, Mesh* mesh, 
			      int nDOF_element, const int* l2g,
			      const double* phi0, double * phi0R,
			      int* reconTag,
			      double zeroTol, int verbose)
{
  /***********************************************************************
    try piecewise linear reconstruction on elements allowing for discontinuities

    ok if nDOF_element > nNodes_element as long as first nSpace+1 dofs are 
    associated with vertex values


  ***********************************************************************/
  using namespace std;
  assert(mesh); assert(phi0); assert(phi0R); 
  bool failed = false;
  //double x0nN[3][3], xIx0[3],dX10[3],xIx01[3];
  //by default just copy over
  copy(phi0,phi0+mesh->nNodes_global,phi0R);
  failed = true;
  return failed;
}

bool copyOverOutsideBand(int N, double tol,
			 const double * in,
			 double * out)
{
  bool failed = false;
  for (int i = 0; i < N; i++)
    {
      if (fabs(out[i]) >= tol)
	out[i] = in[i];
    }
  return failed;
}
}//namespace NearFrontInitialization

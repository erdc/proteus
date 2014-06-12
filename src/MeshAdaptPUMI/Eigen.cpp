#include "Eigen.h"
#include <stdio.h>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <apf.h>
#include <apfMatrix.h>
#include <apfVector.h>

using namespace std;
using namespace apf;  
#define ABS(x) ((x) < 0 ? -(x) : (x))
struct SortStruct
{
  int i;
  double m;
};

int normVt(Vector3 v1,Vector3& nv)
{
  double norm ;

  norm = v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2] ;
  norm = 1./sqrt(norm) ;
  nv[0] = v1[0]*norm ;
  nv[1] = v1[1]*norm ;
  nv[2] = v1[2]*norm ;

  return(1) ;
}
void diffVt(double *a,double *b,double *v)
{
  v[0] = a[0] - b[0] ;
  v[1] = a[1] - b[1] ;
  v[2] = a[2] - b[2] ;
}
void crossProd(Vector3 v1, Vector3 v2, Vector3& cp)
{
  cp[0] = v1[1]*v2[2] - v1[2]*v2[1] ;
  cp[1] = v1[2]*v2[0] - v1[0]*v2[2] ;
  cp[2] = v1[0]*v2[1] - v1[1]*v2[0] ;
}

double dotProd(Vector3 v1, Vector3 v2)
{
  return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2] ;
}

struct greater_abs
{
  bool operator () (const double &a, const double &b)
  {
    return fabs(a) > fabs(b);
  }
};

long eigen (Matrix3x3 pos, Matrix3x3& e, Vector3& v, int checkOrthogonality)
{  
  // characteristic polynomial of T : 
  // solve x^3 + (I[2]/I[3])*x^2 + (I[1]/I[3])*x + (I[0]/I[3])*a3 = 0
  // I1 : first invariant , trace(T)
  // I2 : second invariant , 1/2 (I1^2 -trace(T^2))
  // I3 : third invariant , det T

  double I[4];
  I[3] = 1.0;
  I[2] = - trace(pos);
  I[1] = 0.5 * (I[2]*I[2] - trace2(pos));
  I[0] = - det(pos);

//   printf (" %lf x^3 +  %lf x^2 + %lf x + %lf = 0\n",
//   I[3],I[2],I[1],I[0]);

  // solve x^3 + (I[2]/I[3])*x^2 + (I[1]/I[3])*x + (I[0]/I[3])*a3 = 0
  // solve x^3 + a1 x^2 + a2 x + a3 = 0
  long nbEigen = FindCubicRoots (I,v);

  if (fabs(v[2]) > fabs(v[1]))
    std::swap(v[1],v[2]);
  if (fabs(v[1]) > fabs(v[0]))
    std::swap(v[0],v[1]);
  if (fabs(v[2]) > fabs(v[1]))
    std::swap(v[1],v[2]);  
  
//  std::sort(v,v+3, greater_abs() );
//   printf ("nbEigen = %d %12.5E %12.5E %12.5E\n",nbEigen,v[0],v[1],v[2]);
    
  double result[12];
  int nb_vec=0;

  while(1)
    {
      double a[9] = {pos[0][0]-v[nb_vec],pos[0][1],pos[0][2],
		     pos[1][0],pos[1][1]-v[nb_vec],pos[1][2],
		     pos[2][0],pos[2][1],pos[2][2]-v[nb_vec]};
	
      // eps smaller gives better eigenvals (orig 1.0e-3) 
      double eps = 1.0e-3;
      int nb = 0;
      while (1)
	{
	  nb = NullSpace (a,result,eps,3);
	  if (nb != 0)break;
	  eps *= 2.0;
	}
      int kk=0;
      for (int i=nb_vec;i<nb+nb_vec;i++)
	{
	  e[i][0] = result[0+kk*3];
	  e[i][1] = result[1+kk*3];
	  e[i][2] = result[2+kk*3];
          Vector3 norme;
          normVt(e[i],norme);
          e[i][0]=norme[0]; e[i][1]=norme[1];e[i][2]=norme[2];
	  //printf("%d: %f (%f, %f, %f)\n",i,v[nb_vec],e[i][0],e[i][1],e[i][2]);
	  kk++;
	  if (i == 2 && checkOrthogonality) {
	    int factor;
	    if( !checkUnitaryOthoganal(e,factor) )
	      {
		  
		printf (" %lf x^3 +  %lf x^2 + %lf x + %lf = 0\n",I[3],I[2],I[1],I[0]);
		printf ("nbEigen = %d %12.5E %12.5E %12.5E\n",nbEigen,v[0],v[1],v[2]);
		for(int jj=0; jj<3; jj++ )
		  printf("%d: %f (%f, %f, %f)\n",jj,v[jj],e[jj][0],e[jj][1],e[jj][2]);
		printf("nb=%d nb_vec=%d nbEigen=%d\n",nb,nb_vec,nbEigen);
		printf("WARNING: not orthoganal (eigen)\n\n");
	      }

	    // changing the orientation of thrid vector
	    // such that it follows right hand rule
	     if(factor==-1) {
	     for(int icomp=0;icomp<3;icomp++) {
	       e[3][icomp]=factor*e[3][icomp];
	     }
	    // // cout<<"Changing orientation for third eigen-vector"<<endl;
	     }

	    return nbEigen;
	  }// if (i == 2 && checkOrthog
	}//for (int i=nb_v
      nb_vec += nb;
      if (nb_vec == 3)
	return nbEigen;
      if( nb_vec > 3 )
	return nbEigen;
//	throw;
      if (nb > 3)
	return nbEigen;
//	throw;
    }//while(1)
}
  
int checkUnitaryOthoganal(Matrix3x3 e, int &factor)
{
  int i;
  double dot;  
  Vector3 n;
  double tol=1e-14;
  double cosalpha, alpha;
    
  for( i=0; i<3; i++ ) {
    dot=dotProd(e[i],e[i]);
    if( dot < tol ) 
      { printf("the %d vector in zero length\n",i); return 0; }
    if( ABS(dot - 1.) > tol )
      { printf("the %d vector not unitary. lenthSq=%f\n",i,dot); return 0; }
  }
  dot=dotProd(e[0],e[1]);
  cosalpha=dot/sqrt(dotProd(e[0],e[0])*dotProd(e[1],e[1]));
  alpha = 57.295718*acos(cosalpha);
  if( alpha > 95 || alpha<85 ) {
    printf("first two base vectors not orthognal.  %f\n",alpha);
    return 0;
  }
  crossProd(e[0],e[1],n);
  dot=dotProd(e[2],n);

  if(dot<0.)
    factor=-1;

  cosalpha=dot/sqrt(dotProd(e[2],e[2])*dotProd(n,n));
  alpha = 57.295718*acos(cosalpha);
  if( alpha < 175 && alpha>5 ) {
    printf("third base vector not orthognal to first two.  %f\n", alpha);
    return 0;
  }
  return 1;
}
  
double trace (Matrix3x3 pos)
{
  return pos[0][0] + pos[1][1] + pos[2][2];
}

double trace2 (Matrix3x3 pos)
{
  double a00 =  pos[0][0] * pos[0][0] + 
    pos[1][0] * pos[0][1] + 
    pos[2][0] * pos[0][2]; 
  double a11 =  pos[1][0] * pos[0][1] + 
    pos[1][1] * pos[1][1] + 
    pos[1][2] * pos[2][1]; 
  double a22 =  pos[2][0] * pos[0][2] + 
    pos[2][1] * pos[1][2] + 
    pos[2][2] * pos[2][2];

  return a00 + a11 + a22;
}

double det (Matrix3x3 pos)
{
  return pos[0][0] * (pos[1][1] * pos[2][2] - pos[1][2] * pos[2][1]) -
    pos[0][1] * (pos[1][0] * pos[2][2] - pos[1][2] * pos[2][0]) +
    pos[0][2] * (pos[1][0] * pos[2][1] - pos[1][1] * pos[2][0]);
}

// solve x^2 + b x + c = 0
// x[2] is always set to be zero
long FindQuadraticRoots(const double b, const double c, Vector3& x)
{
  //    printf("Quadratic roots\n");
  x[2]=0.0;
  double delt=b*b-4.*c;
  if( delt >=0 ) {
    delt=sqrt(delt);
    x[0]=(-b+delt)/2.0;
    x[1]=(-b-delt)/2.0;
    return 3;
  }
    
  printf("Imaginary roots, impossible, delt=%f\n",delt);
  return 1;
}
  
// solve x^3 + a1 x^2 + a2 x + a3 = 0
long FindCubicRoots(const double coeff[4], Vector3& x)
{
  double a1 = coeff[2] / coeff[3];
  double a2 = coeff[1] / coeff[3];
  double a3 = coeff[0] / coeff[3];
    
  if( ABS(a3)<1.0e-8 ) 
    return FindQuadraticRoots(a1,a2,x);
    
  double Q = (a1 * a1 - 3 * a2) / 9.;
  double R = (2. * a1 * a1 * a1 - 9. * a1 * a2 + 27. * a3) / 54.;
  double Qcubed = Q * Q * Q;
  double d = Qcubed - R * R;
    
  //    printf ("d = %22.15e Q = %12.5E R = %12.5E Qcubed %12.5E\n",d,Q,R,Qcubed);

  /// three roots, 2 equal 
  if(Qcubed == 0.0 || fabs ( Qcubed - R * R ) < 1.e-8 * (fabs ( Qcubed) + fabs( R * R)) )
    {
      double theta;
      if (Qcubed <= 0.0)theta = acos(1.0);
      else if (R / sqrt(Qcubed) > 1.0)theta = acos(1.0); 
      else if (R / sqrt(Qcubed) < -1.0)theta = acos(-1.0); 
      else theta = acos(R / sqrt(Qcubed));
      double sqrtQ = sqrt(Q);
      //      printf("sqrtQ = %12.5E teta=%12.5E a1=%12.5E\n",sqrt(Q),theta,a1);
      x[0] = -2 * sqrtQ * cos( theta           / 3) - a1 / 3;
      x[1] = -2 * sqrtQ * cos((theta + 2 * M_PI) / 3) - a1 / 3;
      x[2] = -2 * sqrtQ * cos((theta + 4 * M_PI) / 3) - a1 / 3;
      return (3);
    }

  /* Three real roots */
  if (d >= 0.0) {
    double theta = acos(R / sqrt(Qcubed));
    double sqrtQ = sqrt(Q);
    x[0] = -2 * sqrtQ * cos( theta           / 3) - a1 / 3;
    x[1] = -2 * sqrtQ * cos((theta + 2 * M_PI) / 3) - a1 / 3;
    x[2] = -2 * sqrtQ * cos((theta + 4 * M_PI) / 3) - a1 / 3;
    return (3);
  }
    
  /* One real root */
  else {
    printf("IMPOSSIBLE !!!\n");

    double e = pow(sqrt(-d) + fabs(R), 1. / 3.);
    if (R > 0)
      e = -e;
    x[0] = (e + Q / e) - a1 / 3.;
    return (1);
  }
}
  
#define MAXN 32
#define R(i,j)  result[n*(i)+(j)]
  
long NullSpace(const double *a, double *result, double eps, long n)
{
  int r[MAXN], c[MAXN];
  register long i, j, k;
  int jj, kk, t;
  double max, temp;
  int ec;
    
  for (i = 0; i < n; i++)
    r[i] = c[i] = -1;                 /* Reset row and column pivot indices */
    
  // copy the input matrix if not in place
  if (result != a) 
    for (i = 0; i < n*n; i++)  
      result[i] = a[i];
  // rest of algorithm is in place wrt result[]
    
  for (i = 0; i < n; i++) {
    /* Find the biggest element in the remaining submatrix
     * for the next full pivot.
     */
    max = 0.0;
    for (k = 0; k < n; k++) {
      if (r[k] < 0) {
	for (j = 0; j < n; j++) {
	  if ((c[j] < 0) && ((temp = fabs(R(k, j))) > max)) {
	    kk = k;
	    jj = j;
	    max = temp;
	  }
	}
      }
    }
    if (max < eps)
      break;          /* Consider this and all subsequent pivots to be zero */

    c[jj] = kk;                                       /* The row */
    r[kk] = jj;                                       /* and column of the next pivot */
      
    temp = 1.0 / R(kk, jj);
    R(kk, jj) = 1.0;
    for (j = 0; j < n; j++)           /* Should this be for j != jj ? */
      R(kk, j) *= temp;               /* Row equilibration */
      
    for (k = 0; k < n; k++) { /* Row elimination */
      if (k == kk)
	continue;                     /* Don't do a thing to the pivot row */
      temp = R(k, jj);
      R(k, jj) = 0.0;
      for (j = 0; j < n; j++) {
	R(k, j) -= temp * R(kk, j);   /* Subtract row kk from row k */
	if (fabs(R(k, j)) < eps)
	  R(k, j) = 0.0;      /* Flush to zero if too small */
      }
    }
  }
    
  /* Sort into a truncated triangular matrix */
  for (j = 0; j < n; j++) {           /* For all columns... */
    while ((c[j] >= 0) && (j != c[j])) {
      for (k = 0; k < n; k++) {
	if (r[k] < 0) {
	  /* Aha! a null column vector */
	  temp = R(k, j);     /* Get it on top */
	  R(k, j) = R(k, c[j]);
	  R(k, c[j]) = temp;
	}
      }
      t = c[j];                /* Twiddle until pivots are on the diagonal */
      c[j] = c[t];
      c[t] = t;
    }
  }
    
  /* Copy the null space vectors into the top of the A matrix */
  ec = 0;
  for (k = 0; k < n; k++) {
    if (r[k] < 0) {
      R(k, k) = 1.0;                  /* Set the pivot equal to 1 */
      if (ec != k) {
	for (j = 0; j < n; j++) {
	  R(ec, j) = R(k, j);
	}
      }
      ec++;
    }
  }
  /* The first  ec  rows of the matrix  a  are the vectors which are
   * orthogonal to the columns of the matrix  a.
   */
  return (ec);
}

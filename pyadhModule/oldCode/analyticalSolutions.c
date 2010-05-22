int poissonsEquationsExp3D(double K, int nPoints, double t, double* x, double u)
{
  int i=0;
  for (i=0;i<nPoints;i++)
    {
      X = x[i*3+0];      
      Y = x[i*3+1];      
      Z = x[i*3+2];
      K*X*(1.0-X)*Y*(1.0-Y)*Z*(1.0-Z)*exp(X*X + Y*Y + Z*Z);
    }
  return 0;
}

int poissonsEquationExp3D_r(double K, int nPoints, double t, double* x, double u)
{
/*             x = X[0] */
/*             y = X[1] */
/*             z = X[2] */
/*             return  self.K_*(y*(1.0-y)*z*(1.0-z)*(4.0*(1.0-x)*x**3 - 4.0*x**2 + 6.0*x*(1.0-x) - 2.0)+ */
/*                              x*(1.0-x)*z*(1.0-z)*(4.0*(1.0-y)*y**3 - 4.0*y**2 + 6.0*y*(1.0-y) - 2.0)+ */
/*                              x*(1.0-x)*y*(1.0-y)*(4.0*(1.0-z)*z**3 - 4.0*z**2 + 6.0*z*(1.0-z) - 2.0))*exp(x**2 + y**2 + z**2)  */

}

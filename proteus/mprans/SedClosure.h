namespace proteus
{
class cppHsuSedStress
{
public:
  cppHsuSedStress(double parameterIn):
    parameter(parameterIn)
  {}
  inline double  M_sf_x(double porosity, double pf)
  {
    double M_sf_x = parameter*(porosity*porosity + pf);
    return M_sf_x; 
  }
  double parameter;
};
}

#ifndef DOMAIN_H
#define DOMAIN_H

#include <iostream>
#include <vector>
#include <utility>

using namespace std;

class D_base
{
 public:
  D_base(int dimIn = 3, string nameIn = string("Domain"), string unitsIn = string("m"));
  virtual void writeAsymptote(const char* filename)=0;
  virtual void writePoly(const char* filename)=0;
  //virtual void writeXdmf()=0;
  int dim;
  string name;
  string units;
  vector<double> x;//minx, miny, minz
  vector<double> L;//bounding box when x is origin
 };

class RectangularDomain: public D_base
{ 
 public:
  RectangularDomain(int dimIn = 3, string nameIn = string("Domain"), string unitsIn = string("m"));
  RectangularDomain(double x0 ,double Lx, string nameIn = string("Domain"), string unitsIn = string("m"));
  RectangularDomain(double x0, double y0, double Lx, double Ly, string nameIn = string("Domain"), string unitsIn = string("m"));
  RectangularDomain(double x0, double y0, double z0, double Lx, double Ly, double Lz, string nameIn = string("Domain"), string unitsIn = string("m"));
  void writePoly(const char* filename);
  void writeAsymptote(const char* filename);
  //void writeXdmf();
};

class PlanarStraightLineGraphDomain: public D_base
{
 public:
  PlanarStraightLineGraphDomain(int dimIn = 2, string nameIn = string("Domain"), string unitsIn = string("m"));
  PlanarStraightLineGraphDomain(const vector<vector<double> >& verticesIn,const vector< vector<int> >& segmentsIn, const vector< vector<double> >& holesIn, const vector<double>& regionsIn, const vector<int>& vertexFlagsIn, const vector<int>& segmentFlagsIn, string nameIn = string("Domain"),string unitsIn = string("m"));
  virtual ~PlanarStraightLineGraphDomain();

  void getSegmentPartition(vector<int>& s);
  void getBoundingBox(vector<double>& v);
  void readPoly(const char* filename);
  void writePoly(const char* filename);
  void writeAsymptote(const char* filename);
  vector<vector<double> > vertices;
  vector< vector<int> > segments;
  vector< vector<double> > holes;
  vector<double> regions;
  vector<int> nVertices;
  vector<int> nVertexAttributes;
  vector<int> vertexAttributes;
  vector<int> vertexFlags;
  vector<int> segmentFlags;
  vector<int> regionFlags;
  vector<double> v;
  vector<int> vN;
  vector<int> sN;
  vector<int> s;
  vector<int> sFlags;
  vector<int> hN;
  vector<double> h;
  vector<int> rN;
  vector<double> r;
  vector<string> line;
  int xMax, xMin;
  int yMax, yMin;
  int sMax, sMin;
  int hasVertexFlags;
  int hasSegmentFlags;
}; 

/* class TriangulatedSurfaceDomain:public D_base */
/* { */
/*  public: */
/*   TriangulatedSurfaceDomain(); */
/*   virtual ~TriangulatedSurfaceDomain(); */
/*   void readSTL(const char* filename); */
/*   void writePoly(const char* filename); */
/*   void writeAsymptote(const char* filename); */
/*   vector<double[3]> vertices; */
/*   //triangles are vectors analagous to segments in 2D  */
/*   vector<int> triangles; */
/*   vector<double[3]> nodes; */
/* }; */
/* class Mesh2DMDomain: public D_base
   {
    public:
    Mesh2DMDomain();
    virtual ~Mesh2DMDomain();
   };
   class Mesh3DMDomain: D_base;
   {
     public:
     Mesh3DMDomain();
     virtual ~Mesh3DMDomain();
   };
*/
class PiecewiseLinearComplexDomain: public D_base
{
  public:
  PiecewiseLinearComplexDomain(int dimIn = 3, string nameIn = string("Domain"), string unitsIn = string("m"));
  PiecewiseLinearComplexDomain(const vector<vector < vector <double> > >& verticesIn,const vector< vector < vector < vector <int> > > >& facetsIn, const vector<vector<double> >& holesIn, const vector<double>& regionsIn, const vector<int>& vertexFlagsIn,const vector <int>& facetFlagsIn, string nameIn = string("Domain"),string unitsIn = string("m"));
  virtual ~PiecewiseLinearComplexDomain();
  void getBoundingBox(vector <double>& v);
  void readPoly(const char* filename);
  void writePoly(const char* filename);
  void writeAsymptote(const char* filename);
  vector < vector < vector <double> > > vertices;
  vector < vector < vector <int> > > facets;
  vector<double> regions;
  vector<double> holes;
  vector<int> facetHoles;
  vector<int> vertexFlags;
  vector<int> facetFlags;
  int hasVertexFlags;
  int hasFacetFlags;
  vector<int> vN;
  vector<double> v;
  vector<int> f;
  vector<int> fN;

};
#endif

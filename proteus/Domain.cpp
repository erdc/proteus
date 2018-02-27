
#include "Domain.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstring>
#include <cassert>

using namespace std;


D_base::D_base(int dimIn, string nameIn, string unitsIn):
  dim(dimIn),
  name(nameIn),
  units(unitsIn)
{}

RectangularDomain::RectangularDomain(int dimIn, string nameIn, string unitsIn):
  D_base(dimIn,nameIn,unitsIn)
{}

RectangularDomain::RectangularDomain(double x0, double Lx, string nameIn, string unitsIn):
  D_base(1,nameIn,unitsIn)
{
  x.push_back(x0);
  L.push_back(Lx);
}

RectangularDomain::RectangularDomain(double x0, double y0, double Lx ,double Ly, string nameIn, string unitsIn):
  D_base(2,nameIn,unitsIn)
{
  x.push_back(x0);
  x.push_back(y0);
  L.push_back(Lx);
  L.push_back(Ly);
}
  
RectangularDomain::RectangularDomain(double x0, double y0, double z0, double Lx, double Ly, double Lz, string nameIn, string unitsIn):
  D_base(3,nameIn,unitsIn)
{
  x.push_back(x0);
  x.push_back(y0);
  x.push_back(z0);
  L.push_back(Lx);
  L.push_back(Ly);
  L.push_back(Lz);
}

void RectangularDomain::writePoly(const char* filename)
{
  ofstream outfile(filename);
  if(dim == 2)
    {
      
      outfile << "#vertices" << endl;
      outfile << "4 2 0 0 " << endl;
      outfile << "1" << " " << x[0] << " " << x[1] << " 1" << endl;
      outfile << "2" << " " << x[0]+L[0] << " " << x[1] << " " << "2" << endl;
      outfile << "3" << " " << x[0]+L[0] << " " << x[1]+L[1] << " " << "3" << endl;
      outfile << "4" << " " << x[0] << " " << x[1]+L[1] << " " << "4" << endl;
      outfile << "#segments" << endl;
      outfile << "4 1" << endl;
      outfile << "1 1 2 1" << endl;
      outfile << "2 2 3 2" << endl;
      outfile << "3 3 4 3" << endl;
      outfile << "4 4 1 4" << endl;
      outfile << "#holes" << endl;
      outfile << "0" << endl;
      outfile << "#regions" << endl;
      outfile << "1" << endl;
      outfile << "1" << " " << (x[0]+L[0])/2 << " " << (x[1]+L[1])/2 << " " << "1" << endl;
      outfile.close();
     }
  else if (dim == 3)
    {
      
      outfile << "#vertices" << endl;
      outfile << "8 3 0 0" << endl;
      outfile << "1" << " " << x[0] << " " << x[1] << " " << x[2] << " 1" << endl;
      outfile << "2" << " " << x[0] << " " << x[1]+L[1] << " " << x[2] << " 2" << endl;
      outfile << "3" << " " << x[0]+L[0] << " " << x[1]+L[1] << " " << x[2] << " 3" << endl;
      outfile << "4" << " " << x[0]+L[0] << " " << x[1] << " " << x[2] << " 4" << endl;
      outfile << "5" << " " << x[0] << " " << x[1] << " " << x[2]+L[2] << " 5" << endl;
      outfile << "6" << " " << x[0] << " " << x[1]+L[1] << " " << x[2]+L[2] << " 6" << endl;
      outfile << "7" << " " << x[0]+L[0] << " " << x[1]+L[1] << " " << x[2]+L[2] << " 7" << endl;
      outfile << "8" << " " << x[0]+L[0] << " " << x[1] << " " << x[2]+L[2] << " 8" << endl;
      outfile << "#facets" << endl;
      outfile << "6 1" << endl;
      outfile << "1 0 1" << endl;
      outfile << "4 1 2 3 4\t\t#bottom" << endl;
      outfile << "1 0 2" << endl;
      outfile << "4 1 5 8 4\t\t#front" << endl;
      outfile << "1 0 3" << endl;
      outfile << "4 1 5 6 2\t\t#left" << endl;
      outfile << "1 0 4" << endl;
      outfile << "4 5 6 7 8\t\t#top" << endl;
      outfile << "1 0 5" << endl;
      outfile << "4 4 8 7 3\t\t#right" << endl;
      outfile << "1 0 6" << endl;
      outfile << "4 3 7 6 2\t\t#back" << endl;
      outfile << "#holes" << endl;
      outfile << "0" << endl;
      outfile << "#regions" << endl;
      outfile << "1" << endl;
      outfile << "1" << " " << (x[0]+L[0])/2.0 << " " << (x[1]+L[1])/2.0 << " " << (x[2]+L[2])/2.9 << " 1" << endl; 
      outfile.close();
    }
}


void RectangularDomain::writeAsymptote(const char* filename)
{
  
    ofstream outfile(filename);
  

    if (dim == 1)
      {//interval1D
	outfile << "unitsize(4.0 inches/" << L[0] << ");" << endl;
	outfile << "size(6 inches);" << endl;
	outfile << "real L=" << L[0] << ";" << endl;
	outfile << "real offset=.0125L;" << endl;
	outfile << "real x=" << x[0] << ";" << endl;
	outfile << "string str=\"$" << L[0] << "\\mbox{" << units << "}$\";" << endl;
	outfile << "import graph;" << endl;
	outfile << "import palette;" << endl;
	outfile << "pen[] allPens=Wheel();" << endl;
	outfile << "pen[] myPens = new pen[3];" << endl;
	outfile << "for(int i=0;i< 3;++i)" << endl;
	outfile << "{" << endl;
	outfile << "   int iPen = round(i*allPens.length/3);" << endl;
	outfile << "   myPens[i] = allPens[iPen];" << endl;
	outfile << " }" << endl;
	outfile << "draw((x,0)--(x+L,0)^^(x,offset)--(x+L,offset),myPens[0]);" << endl;
	outfile << "draw((x,0)--(x,offset),myPens[1]);" << endl;
	outfile << "draw((x+L,0)--(x+L,offset),myPens[2]);" << endl;
	outfile << "draw(str,(x,-offset)--(x+L,-offset),black,Bars,Arrows,PenMargins)" << ";" <<  endl;
      }else if (dim == 2)
      {

	//Rectangle 2D
	outfile << "import math;" << endl;
	outfile << "import graph;" << endl;
	outfile << "unitsize(4.0 inches /" << L[0] << ");" << endl;
	outfile << "size(6 inches);" << endl;	
	outfile << "real Lx = " << L[0] << ";" << endl;
	outfile << "real Ly = " << L[1] << ";" << endl;
	outfile << "real offset = 0.0125Lx;" << endl;
	outfile << "real x = " << x[0] << ";" << endl;
	outfile << "real y = " << x[1] << ";" << endl;
	outfile << "string strx=\"$" << L[0]<< "\\mbox{" << units << "}$\";" << endl;
	outfile << "string stry=\"$" << L[1]<< "\\mbox{" << units << "}$\";" << endl;
	outfile << "import graph;" << endl;
	outfile << "import palette;" << endl;
	outfile << "pen[] allPens=Wheel();" << endl;
	outfile << "pen[] myPens = new pen[4];" << endl;
	outfile << "for(int i=0;i< 4;++i)" << endl;
	outfile << "{" << endl;
	outfile << "int iPen = round(i*allPens.length/4);" << endl;
	outfile << " myPens[i] = allPens[iPen];" << endl;
	outfile << "}" << endl;
	outfile << "draw((x,y)--(x+Lx,y),myPens[0]);" << endl;
	outfile << "draw((x,y+Ly)--(x+Lx,y+Ly),myPens[1]);" << endl;
	outfile << "draw((x,y)--(x,y+Ly),myPens[2]);" << endl;
	outfile << "draw((x+Lx,y)--(x+Lx,y+Ly),myPens[3]);" << endl;
	outfile << "draw(strx,(x,y-offset)--(x+Lx,y-offset),S,black,Bars,Arrows,PenMargins);" << endl;
	outfile << "draw(stry,(x-offset,y)--(x-offset,y+Ly),W,black,Bars,Arrows,PenMargins);" << endl;
      }else 
      {
	outfile << "import three;" << endl;
	outfile << "currentlight=adobe;" << endl;
	outfile << "currentprojection=perspective(-2,-2,1,up=Z,target=O,showtarget=true,autoadjust=true,center=false);" << endl;
	outfile << "unitsize(4.0 inches/" << L[0] <<");" << endl;
	outfile << "size(6 inches);" << endl;
	outfile << "real Lx=" << L[0] << ";" << endl;
	outfile << "real Ly=" << L[1] << ";" << endl;
	outfile << "real Lz=" << L[2] << ";" << endl;
	outfile << "real offset=.0125Lx;" << endl;
	outfile << "real x=" << x[0] << ";" << endl;
	outfile << "real y=" << x[1] << ";" << endl;
	outfile << "real z=" << x[2]<< ";" << endl;
	outfile << "triple p0=(" << x[0] << "," << x[1] << ","<< x[2] << ");" << endl;
	outfile << "triple p1=(" << x[0] << "," << x[1]+L[1] << ","<< x[2] << ");" << endl;
	outfile << "triple p2=(" << x[0]+L[0] << "," << x[1]+L[1] << "," << x[2] << ");" << endl;
	outfile << "triple p3=(" << x[0]+L[0] << "," << x[1] << "," << x[2] << ");" << endl;
	outfile << "triple p4=(" << x[0] << "," << x[1] << ","<< x[2]+L[2] << ");" << endl;
	outfile << "triple p5=(" << x[0] << "," << x[1]+L[1] << "," << x[2]+ L[2] << ");" << endl;
	outfile << "triple p6=(x+Lx,y+Ly,z+Lz);" << endl;
	outfile << "triple p7=(x+Lx,y,z+Lz);" << endl;
	outfile << "string strx=\"$" << L[0] << "\\mbox{" << units << "}$\";"<< endl;
	outfile << "string stry=\"$" << L[1] << "\\mbox{" << units << "}$\";" << endl;
	outfile << "string strz=\"$" << L[2] << "\\mbox{" << units << "}$\";" << endl;
	outfile << "draw(surface(p0--p1--p2--p3--cycle),red);" << endl;
	outfile << "draw(surface(p4--p5--p6--p7--cycle),blue);" << endl;
	outfile << "draw(surface(p0--p1--p5--p4--cycle),green);" << endl;
	outfile << "draw(surface(p1--p2--p6--p5--cycle),orange);" << endl;
	outfile << "draw(surface(p2--p3--p7--p6--cycle),purple);" << endl;
	outfile << "draw(surface(p3--p0--p4--p7--cycle),yellow);" << endl;
	outfile << "draw(strx,(x,y-offset,z-offset)--(x+Lx,y-offset,z-offset),-(Y+Z),black,Bars3(Z),Arrows3);" << endl;
	outfile << "draw(stry,(x-offset,y,z-offset)--(x-offset,y+Ly,z-offset),-(X+Z),black,Bars3(Z),Arrows3);" << endl;
	outfile << "draw(strz,(x-offset,y-offset,z)--(x-offset,y-offset,z+Lz),-(X+Z),black,Bars3(Z),Arrows3);" << endl;
	outfile << "shipout();" << endl;
      }
    outfile.close();
}


PlanarStraightLineGraphDomain::PlanarStraightLineGraphDomain(int dimIn, string nameIn, string unitsIn):
  D_base(dimIn,nameIn,unitsIn)
{
}

PlanarStraightLineGraphDomain::PlanarStraightLineGraphDomain(const vector<vector<double> >& verticesIn, const vector< vector<int> >& segmentsIn, const vector< vector<double> >& holesIn, const vector<double>& regionsIn, const vector<int>& vertexFlagsIn, const vector<int>& segmentFlagsIn, string nameIn,string unitsIn):
  vertices(verticesIn),
  segments(segmentsIn),
  holes(holesIn),
  regions(regionsIn),
  vertexFlags(vertexFlagsIn),
  segmentFlags(segmentFlagsIn),
  D_base(2,nameIn,unitsIn)
{    
}

PlanarStraightLineGraphDomain::~PlanarStraightLineGraphDomain()
{}
void PlanarStraightLineGraphDomain::writePoly(const char* filename)
{
  ofstream outfile(filename);

  if(!filename)
    {
      if(vertexFlags.size() > 0)
	hasVertexFlags = 1;
      else
	hasVertexFlags = 0;
      if(segmentFlags.size() > 0)
	hasSegmentFlags = 1;
      else
	hasSegmentFlags = 0;
      outfile << vertices.size() << "\t" << "2\t" << "0\t" << hasVertexFlags << endl;
      //write the vertices
      for(int i = 0; i < vertices.size(); i++)
	for(int j =0; j < 3; j++)
	  {
	    vN[i*3 + j] = vertices[i][j];
	    outfile << vN.size()+1 << "\t" <<  v[0] << "\t" <<  v[1] << endl;
	  }
      //import pdb; pdb.set_trace()
     //  if(vertexFlags.size() > 0)
// 	outfile << vertexFlags[vN] << endl;
//       else
// 	outfile << endl;
//       //write the segments
//       outfile << segments.size() << "\t" << hasSegmentFlags << endl;
//       for(int i = 0; i < segments.size();i++)
// 	for(int j = 0; j < 3; j++)
// 	  { 
// 	    sN[i*3+j] = segments[i][j];
// 	    outfile << sN.size()+1 <<"\t" << s[0]+1 << "\t" << s[1]+1 << endl;
// 	  }
//       if(segmentFlags.size() > 0)
// 	outfile << segmentFlags[sN] << endl;
//       else
// 	outfile << endl;
//       if(holes.size() > 0)
// 	{
// 	  outfile << holes.size() << endl;
// 	  for(int i = 0; i < holes.size(); i++)
// 	    for(int j = 0; j < 3; j++)
// 	      { 
// 		hN[i*3 + j] = holes[i][j];
// 		outfile << hN.size()+1 << "\t" << h[0] << "\t" << h[1] << endl;
// 	      }
// 	}
//       else
// 	outfile << "0" << endl;
//       if(regions.size() > 0)
// 	{
// 	  outfile << regions.size() << endl;
// 	  for(int i = 0; i < regions.size(); i++)
// 	    {
// 	      rN[i*3] = regions[i];
// 	      outfile << rN.size()+1 << "\t" << r[0] << "\t" << r[1] << endl;
// 	    }
// 	  if(regionFlags.size() > 0)
// 	    outfile << regionFlags[rN] << endl;
// 	  else
// 	    outfile << endl;
// 	}
//       else
// 	outfile << "0" << endl;
//       outfile.close();
    }
  else
    cout << "File already exists, not writing polyfile: " << filename << endl;
}

void PlanarStraightLineGraphDomain::writeAsymptote(const char* filename)
{
  ofstream outfile(filename);
 
  
  outfile << "unitsize(4.0 inches/ " << L[0] << ");" << endl;
  outfile << "size(6 inches);" << endl;
  outfile << "real Lx=" << L[0] << ";" << endl;
  outfile << "real Ly=" << L[1]<< ";" << endl;
  outfile << "real offset=0.0125Lx;" << endl;
  outfile << "real x=" << x[0] << ";" << endl;
  outfile << "real y=" << x[1] << ";" << endl;
  outfile << "string strx=\"$" << L[0] << "\\mbox{" << units << "}$\";" << endl;
  outfile << "string stry=\"$" << L[1] << "\\mbox{" << units << "}$\";" << endl;
  outfile << "draw(strx,(x,y-offset)--(x+Lx,y-offset),S,black,Bars,Arrows,PenMargins);" << endl;
  outfile << "draw(stry,(x-offset,y)--(x-offset,y+Ly),W,black,Bars,Arrows,PenMargins);" << endl;
  outfile << "import graph;" << endl;
  outfile << "import palette;" << endl;
  outfile << "pen[] allPens=wheel();" << endl;
  outfile << "pen[] myPens = new pen[" << segmentFlags.size() + 1 << "];" << endl;
  outfile << "for(int i= 0; i<" << segmentFlags.size() +1 << ";++i)" << endl;
  outfile << "{" << endl;
  outfile << "int iPen = round(i*allPens.length//" << segmentFlags.size() + 1 << ");" << endl;
  outfile << "myPens[i] = allPens[iPen];" << endl;
  outfile << "}" << endl;
  outfile << "shipout();" << endl;
  for(int i = 0; i<segments.size();i++)
    for(int j = 0; j < 3; j++)
      s[i*3 + j] = segments[i][j];
  for(int i = 0; i < segmentFlags.size();i++)
    sFlags[i*3] = segmentFlags[i]; 
  // outfile << "draw((" << vertices[s[0]] << "," << vertices[s[1]] << ")--(" <<

//   //for loop in python code
 
//   //#now loop over segments
//   //    for s,sFlag in zip(self.segments,self.segmentFlags):
//   //      fileString+="draw((%f,%f)--(%f,%f),myPens[%d]);\n" % tuple(self.vertices[s[0]]+self.vertices[s[1]]+[sFlag])

  outfile.close();
}

void PlanarStraightLineGraphDomain::getSegmentPartition(vector<int>& s)
{
 //  sMin = segmentFlags[0];
//   sMax = segmentFlags[0];
//   for(int i = 0; i < segmentFlags.size(); i++)
//     {
//       s[i] = segmentFlags[i];
     
//     }
//  if(s > sMax)
// 	sMax = s;
//       if(s < sMin)
// 	sMin = s;
 

}

void PlanarStraightLineGraphDomain::getBoundingBox(vector<double>& v)
{
  xMax = vertices[0][0];
  xMin = vertices[0][0];
  yMax = vertices[0][1];
  yMin = vertices[0][1];
  for(int i =0; i < vertices.size(); i++) 
    for(int j = 0; j <3; j++)
      { 
	v[i*3 + j] = vertices[i][j];
	if (v[0] > xMax)
	  xMax = v[0];
	if (v[1] > yMax)
	  yMax = v[1];
	if (v[0] < xMin)
	  xMin = v[0];
	if (v[1] < yMin)
	  yMin = v[1];
      }
  x.push_back(xMin);
  x.push_back(yMin);
  L.push_back(xMax-xMin);
  L.push_back(yMax-yMin);
}


void PlanarStraightLineGraphDomain::readPoly(const char* filename)
{
//   ifstream infile;
//   infile.open(filename);

//   while(!infile)
//     {
//       infile.close();
//       infile.clear();

//       cout << "Error: Invalid filename. Please re-enter: ";
//       cin >> filename;

//       infile.open(filename);
//     }
//   while(firstLine.length() == 0 || firstLine[0][0] == "#")
//     {
//       infile >> nVertices >> dim >> nVertexAttributes >> hasVertexFlags;
//     }
//   infile.close();
}




// //Constructors for TriangulatedSurfaceDomain


// TriangulatedSurfaceDomain::TriangulatedSurfaceDomain()
// {
//   //vertices = 0.0;
//   triangles = 0.0;
//   // nodes = 0.0;
//   name = "DefaultTriangulatedSurfaceDomain";
//   units = "m";
//   dim = 3;
// }


// TriangulatedSurfaceDomain::~TriangulatedSurfaceDomain()
 // {}

// void TriangulatedSurfaceDomain::writeAsymptote(const char* filename)
// {
//   ofstream outfile(filename);

// }





//PiecewiseLinearComplexDomain
PiecewiseLinearComplexDomain::PiecewiseLinearComplexDomain(int dimIn, string nameIn, string unitsIn):
  D_base(dimIn,nameIn,unitsIn)
{}

PiecewiseLinearComplexDomain::PiecewiseLinearComplexDomain(const vector<vector < vector <double> > >& verticesIn,const vector< vector < vector < vector <int> > > >& facetsIn, const vector<vector<double> >& holesIn, const vector<double>& regionsIn, const vector<int>& vertexFlagsIn,const vector <int>& facetFlagsIn, string nameIn = string("Domain"),string unitsIn = string("m")):
  vertices(verticesIn),
  facets(facetsIn),
  holes(holesIn),
  regions(regionsIn),
  vertexFlags(vertexFlagsIn),
  facetFlags(facetFlagsIn),
  D_base(3,nameIn,unitsIn)
{}
PiecewiseLinearComplexDomain::~PiecewiseLinearComplexDomain()
{}

void PiecewiseLinearComplexDomain::writeAsymptote(const char* filename)
{
  ofstream outfile(filename);


  outfile << "import three;" << endl;
  outfile << "currentlight=adobe;" << endl;
  outfile << "currentprojection=perspective(-2,-2,1,up=Z,target=(%(xc)f,%(yc)f,%(zc)f),showtarget=true,autoadjust=true,center=false);" << endl;
  outfile << "unitsize(4.0 inches/%(Lx)f);" << endl;
  outfile << "size(6 inches);" << endl;
  outfile << "real Lx=%(Lx)f;" << endl;
  outfile << "real Ly=%(Ly)f;" << endl;
  outfile << "real Lz=%(Lz)f;" << endl;
  outfile << "real offset=.0125Lx;" << endl;
  outfile << "real x=%(x)f;" << endl;
  outfile << "real y=%(y)f;" << endl;
  outfile << "real z=%(z)f;" << endl;
  outfile << "triple p0=(x,y,z);" << endl;
  outfile << "triple p1=(x,y+Ly,z);" << endl;
  outfile << "triple p2=(x+Lx,y+Ly,z);" << endl;
  outfile << "triple p3=(x+Lx,y,z);" << endl;
  outfile << "triple p4=(x,y,z+Lz);" << endl;
  outfile << "triple p5=(x,y+Ly,z+Lz);" << endl;
  outfile << "triple p6=(x+Lx,y+Ly,z+Lz);" << endl;
  outfile << "triple p7=(x+Lx,y,z+Lz);" << endl;
  outfile << "string strx=\"$%(Lx)2.2f\\mbox{%(units)s}$\";" << endl;
  outfile << "string stry=\"$%(Ly)2.2f\\mbox{%(units)s}$\";" << endl;
  outfile << "string strz=\"$%(Lz)2.2f\\mbox{%(units)s}$\";" << endl;
  outfile << "draw(strx,(x,y-offset,z-offset)--(x+Lx,y-offset,z-offset),-(Y+Z),black,Bars3(Z),Arrows3);" << endl;
  outfile << "draw(stry,(x-offset,y,z-offset)--(x-offset,y+Ly,z-offset),-(X+Z),black,Bars3(Z),Arrows3);" << endl;
  outfile << "draw(strz,(x-offset,y-offset,z)--(x-offset,y-offset,z+Lz),-(X+Y),black,Bars3(X),Arrows3);" << endl;

//   for (int i = 0; i < facets.size(); i++)
//     for(int j = 0; j < vertexList.size(); j++)
      
//       outfile << "draw(surface((" << vertices[vertexList[0]] << "," <<  vertices[vertexList[1]] << "," <<  vertices[vertexList[2]] << "))";
//   for(int k = 0; k < vN; k++)
//     //outfile << "(--(" << vertices[vN] << 
//     //Python Loop found in Domain.py  
    /*
      {'Lx':self.L[0],'Ly':self.L[1],'Lz':self.L[2],'x':self.x[0],'y':self.x[1],'z':self.x[2],'units':self.units,
      'xc':(self.x[0]+0.5*self.L[0]),
      'yc':(self.x[1]+0.5*self.L[1]),
      'zc':(self.x[2]+0.5*self.L[2])}
      f.write(fileString)
      for facet in self.facets:
      for vertexList in facet:
      f.write("draw(surface((%f,%f,%f)" % tuple(self.vertices[vertexList[0]]))
      for vN in vertexList[1:]:
      f.write("--(%f,%f,%f)" % tuple(self.vertices[vN]))
      f.write("--cycle));\n")
    
  
      outfile.close();

    */   
}
// void PiecewiseLinearComplex::getBoundingBox(vector<double> & v)
// {}		  
void PiecewiseLinearComplexDomain::writePoly(const char* filename)
{
  if(!filename)
    {
      ofstream outfile;
      
      if(vertexFlags.size() > 1)
	hasVertexFlags = 1;
      else 
	hasVertexFlags = 0;
      if(facetFlags >1)
	hasFacetFlags.size() = 1;
      else
	hasFacetFlags = 0;
      outfile << vertices.size() << " 3" << " 0 " << hasVertexFlags << endl;
      //write the vertices
      for(int i = 0; i < vertices.size(); i++)
	for(int j =0; j < 3; j++)
	  {
	    vN[i*3 + j] = vertices[i][j];
	    outfile << vN.size()+1 << "\t" <<  v[0] << "\t" <<  v[1] << endl;
	  }
      if(vertexFlags.size() > 1);
      outfile << vertexFlags[vN] << endl;
      else
	outfile << endl;
      //write facets
      outfile << facets.size() << " " << hasFacetFlags;
      for(int i = 0; i < facets.size(); i++)
	for(int j = 0; j < 3; j++)
	  {
	    fN[i*3+j] = facets[i][j];
	  }
      if(facetHoles.size() > 1)
	nFacetHoles = facetHoles[fN].size();
      else
	nFacetHoles = 0;
      if(hasFacetFlags = 1)
	outfile << f.size() << " " << nFacetHoles << " " << facetFlags[fN] << endl;
      else
	outfile << fN.size()+1 << " " << nFacetHoles << endl;
      for(int i = 0; i < f.size(); i++)
	{
	  f = segmentList[i];
	  outfile << segmentList.size()+1 << endl;
	}
	
	for(int j = 0; j < segmentList.size(); j++)
	  {
	    segmentList = vN[j];
	    outfile << vN+1 << " ";
	  }
	outfile << endl;
	if(facetHoles > 0)
	  {
	    for(int i = 0; i < facetHoles[fN].size(); i++)
	      {
		h = facetHoles[i];
		outfile << facetHoles[i] << " " << endl;
	      }
	    if(holes > 0)
	      {
		outfile << holes.size() << endl;
		for(int i = 0; i < holes.size(); i++)
		  for(int j = 0; j < 3; j++)
		    {
		      hN[i*3+j] = holes[i][j]; 
		      outfile << hN+1 << " " << h[0] << " " << h[1] << endl;
		      if(regionFlags > 0)
			outfile << regionFlags[rN] << endl;
		      else
			outfile << endl;
		    
		    }
		
	      }
	    else
	      outfile << "0" << endl;
	    
	  }
	else
	  outfile.close();
    }
  else
    cout <<"File already exits, not writing polyfile: " << filename << endl;
}

// void PiecewiseLinearComplex::readPoly(const char* filename)
// {

//   /*
//     self.polyfile = fileprefix+".poly"
//     self.name=fileprefix
//     f = open(self.polyfile,'r')
//     firstLine=f.readline().split()
//     while len(firstLine) == 0 or firstLine[0][0] == '#':
//     firstLine = f.readline()
//     nVertices = int(firstLine[0])
//     self.dim = int(firstLine[1])
//     assert(self.dim == 3)
//     nVertexAttributes = int(firstLine[2])
//     self.vertexAttributes = [[] for j in range(nVertexAttributes)]
//     hasVertexFlag = bool(firstLine[3])
//     self.vertices=[]
//     self.base=None
//     if hasVertexFlag:
//     self.vertexFlags=[]
//     for i in range(nVertices):
//     line = f.readline().split()
//     while len(line) == 0 or line[0][0] == '#':
//     line = f.readline().split()
//     if self.base is None:
//     self.base = int(line[0])
//     self.vertices.append([float(line[1]),float(line[2]),float(line[3])])
//     for j in range(nVertexAttributes):
//     self.attributes[j].append(float(line[4+j]))
//     if hasVertexFlag:
//     self.vertexFlags.append(line[4+nVertexAttributes+1])
//     facetLine = f.readline().split()
//     while len(facetLine) == 0 or facetLine[0][0] == '#':
//     facetLine = f.readline().split()
//     nFacets = int(facetLine[0])
//     hasFacetFlag = bool(facetLine[1])
//     self.facets=[]
//     self.facetHoles=[]
//     if hasFacetFlag:
//     self.facetFlags=[]
//     for i in range(nFacets):
//     line = f.readline().split()
//     while len(line) == 0 or line[0][0] == '#':
//     line = f.readline().split()
//     self.facets.append([])
//     self.facetHoles.append([])
//     nPolygons = int(line[0])
//     nHoles = int(line[1])
//     if hasFacetFlag:
//     self.facetFlags.append(int(line[2]))
//     for j in range(nPolygons):
//     line = f.readline().split()
//     while len(line) == 0 or line[0][0] == '#':
//     line = f.readline().split()
//     nSegments = int(line[0])
//     self.facets[-1].append([int(line[1+k])-self.base for k in range(nSegments)])
//     for j in range(nHoles):
//     line = f.readline().split()
//     while len(line) == 0 or line[0][0] == '#':
//     line = f.readline().split()
//     self.facetHoles[-1].append([float(line[0]),
//     float(line[1]),
//     float(line[2])])
//     holeLine = f.readline().split()
//     while len(holeLine) == 0 or holeLine[0][0] == '#':
//     holeLine = f.readline().split()
//     nHoles = int(holeLine[0])
//     self.holes=[]
//     for i in range(nHoles):
//     line =  f.readline().split()
//      while len(line) == 0 or line[0][0] == '#':
//     line = f.readline().split()
//     self.holes.append([float(line[0]),float(line[1])])
//     regionLine = f.readline().split()
//     while len(regionLine) == 0 or regionLine[0][0] == '#':
//     regionLine = f.readline().split()
//     nRegions = int(regionLine[0])
//     self.regions=[]
//     self.regionAttributes=[]
//     self.areaConstraings=[]
//     for i in range(nRegions):
//     line =  f.readline().split()
//     while len(line) == 0 or line[0][0] == '#':
//     line = f.readline().split()
//     self.regions.append([float(line[1]),float(line[2]),float(line[3])])
//     if len(line) > 4:
//     self.regionAttributes.append(int(line[4]))
//     if len(line) > 5:
//     self.areaConstraint.append(float(line[5]))                
//     self.getBoundingBox()
//   */
//   }





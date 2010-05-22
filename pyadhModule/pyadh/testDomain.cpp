#include "Domain.h"
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main()
{

  string polyfile;
  string filename;
  string filenameIn;
  string units;
  string name;
  double x0, y0, z0, Lx, Ly, Lz ;//interval variables
  vector<double> x;
  vector<double> y;
  vector<double> z;
  vector<int> a;
  vector<int> b;
  vector<int> c;
  vector<int> d;
  vector<vector<double> > vertices2;
  vector<vector<vector<double> > > vertices3;
  vector< vector<int> > segments;
  vector< vector<double> > holes;
  vector<double> regions;
  vector<int> vertexFlags;
  vector<int> segmentFlags;
  vector< vector < vector < vector<int> > > >facets;
  int dim;
  int answer;
  
  ifstream inFile;


  cout << "Enter the domain, 1 for Rectangular Domain, 2 for Planar Straight Line Graph, 3 for Triangulated Surface Domain, and 4 for Piecewise Linear Complex Domain: " << endl;
  cin >> answer;


    if (answer == 1)
      {
	//Asking for the Dimensions
	cout << "Enter 1 for 1D, 2 for 2D, and 3 for 3D asymptote rectangle:"<<endl;
	cin >> dim;


	/******* if else below has validation*******
	//Validation for Dimensions
	while (dim != 1 && dim != 2 && dim != 3)
	  {
	    cout << "Error:  Domain object must have dimension 1, 2 or 3. Please re-enter: " <<  endl;
	    cin >> dim;
	  }
  
	*/
	if(dim == 1)
	  {
	    //Asking for the points
	    cout << "Enter the x0 coordinate:" << endl;
	    cin >> x0;
	    cout << "Enter the Lx coordinate:" << endl;
	    cin >> Lx;
	    
	    RectangularDomain rect(x0, Lx, name, units);
	

	    //Entering filename
	    cout<<"Enter filename"<<endl;
	    cin >> filenameIn;
	    filename = filenameIn + ".asy";

	    cout << "Writing to "<<filename<< "..." << endl;

	    rect.writeAsymptote(filename.c_str()); //function call
    
	  }
	else if(dim == 2)
	  {
	    //Asking for the points
	    cout << "Enter the x0:" << endl;
	    cin >> x0;
	    cout << "Enter the Lx:" << endl;
	    cin >> Lx;
	    cout << "Enter the y0:" << endl;
	    cin >> y0;
	    cout << "Enter the Ly:" << endl;
	    cin >> Ly;
	    
	
	    RectangularDomain rect2(x0, y0, Lx, Ly, name, units);
 
	    //Entering filename
	    cout<<"Enter filename"<<endl;
	    cin >> filenameIn;
	    filename = filenameIn + ".asy";
	    cout << "Writing to "<< filename << "..." << endl;
	    
	    polyfile = filenameIn + ".poly";
	    
	    cout << "Writing to "<< polyfile << "..." << endl;
	    
	    rect2.writeAsymptote(filename.c_str()); //function call
	    rect2.writePoly(polyfile.c_str());
	  }
	else 
	  {
	    //Asking for the points
	    cout << "Enter the x0:" << endl;
	    cin >> x0;
	    cout << "Enter the Lx:" << endl;
	    cin >> Lx;
	    cout << "Enter the y0:" << endl;
	    cin >> y0;
	    cout << "Enter the Ly:" << endl;
	    cin >> Ly;
	    cout << "Enter the z0:" << endl;
	    cin >> z0;
	    cout << "Enter the Lz:" << endl;
	    cin >> Lz;

	RectangularDomain rect3(x0, y0, z0, Lx, Ly, Lz, name, units);
   
	    //Entering filename
	    cout<<"Enter filename"<<endl;
	    cin >> filenameIn;
	    filename = filenameIn + ".asy";
	    cout << "Writing to "<<filename<< "..." <<endl;
	    polyfile = filenameIn + ".poly";
	    cout << "Writing to" << polyfile << "..."<< endl;
	    rect3.writeAsymptote(filename.c_str());//function call
	    rect3.writePoly(polyfile.c_str());
	  }
      }
   
      
   
    else if (answer == 2)
      {
	x.push_back(0.0);
	x.push_back(0.0);
	x.push_back(1.0);
	y.push_back(0.0);
	y.push_back(1.0);
	y.push_back(0.0);

	vertices2.push_back(x);
	vertices2.push_back(y);

	a.push_back(0);
	a.push_back(1);
	a.push_back(2);
	b.push_back(1);
	b.push_back(2);
	b.push_back(0);
	
	segments.push_back(a);
	segments.push_back(b);

	PlanarStraightLineGraphDomain PSLG(vertices2,segments,holes,regions,vertexFlags,segmentFlags, name, units);

	cout << "Enter filename: ";
	cin >> filenameIn;
	filename = filenameIn + ".asy";
	polyfile = filenameIn + ".poly";
	PSLG.writeAsymptote(filename.c_str());
	PSLG.writePoly(polyfile.c_str());

      }   
//     else if(answer == 3)
//       { 
  
// 	// TriangulatedSurfaceDomain TSD( units, name);
      else if (answer == 4)
      {
	x.push_back(0.0);
	x.push_back(0.0);
	x.push_back(1.0);
	x.push_back(1.0);
	x.push_back(0.0);
	x.push_back(0.0);
	x.push_back(1.0);
	x.push_back(1.0);
	y.push_back(0.0); 
	y.push_back(1.0);
	y.push_back(1.0);
	y.push_back(0.0);
	y.push_back(0.0);
	y.push_back(1.0);
	y.push_back(1.0);
	y.push_back(0.0);
	z.push_back(0.0);
	z.push_back(0.0);	
	z.push_back(0.0);
	z.push_back(0.0);
	z.push_back(1.0);
	z.push_back(1.0);	
	z.push_back(1.0);
	z.push_back(1.0);

	vertices3.push_back(x);
	vertices3.push_back(y);
	vertices3.push_back(z);

	a.push_back(0);
	a.push_back(0);
	a.push_back(0);
	a.push_back(4);
	a.push_back(3);
	a.push_back(2);
	b.push_back(1);
	b.push_back(4);
	b.push_back(4);
	b.push_back(5);
	b.push_back(7);
	b.push_back(6);
	c.push_back(2);
	c.push_back(7);
	c.push_back(5);
	c.push_back(6);
	c.push_back(6);
	c.push_back(5);
	d.push_back(3);
	d.push_back(3);
	d.push_back(1);
	d.push_back(7);
	d.push_back(2);
	d.push_back(1);
      

	facets.push_back(a);
	facets.push_back(b);
	facets.push_back(c);
	facets.push_back(d);
	
	PiecewiseLinearComplexDomain PLC(vertices3, facets, name, units);

	cout << "Enter filename:" << endl;
	cin >> filenameIn;
	filename = filenameIn + ".asy";
	polyfile = filenameIn + ".poly";

	PLC.writeAsymptote(filename.c_str());
	PLC.writePoly(polyfile.c_str());
	}
 else
  {
    cout << "The domain entered is not supported" << endl;
  }
    
  return 0;
}

#include "mesh.h"

int main()
{
  //todo write mesh destructor so this test doesn't leak memory
  using namespace std;
  int nx,ny,nz;
  double start,stop,diff;
  string outfilename,infilename;
  cout<<"Testing simple mesh generation"<<endl;
  cout<<"Enter nx,ny, and nz"<<endl;
  cin>>nx>>ny>>nz;
  cout<<"Print mesh? (n/filename/cout)"<<endl;
  cin>>outfilename;
  Mesh mesh;
  initializeMesh(mesh);
  start = CurrentTime();
  if (ny == 1 && nz == 1)
    {
      cout<<"Generating Edge Mesh Elements"<<endl;
      edgeMeshElements(nx,mesh);
      cout<<"Generating Edge Mesh Nodes"<<endl;
      regularEdgeMeshNodes(nx,1.0,mesh);
    }
  else if (nz == 1)
    {
      cout<<"Generating Triangular Mesh Elements"<<endl;
      regularRectangularToTriangularMeshElements(nx,ny,mesh);
      cout<<"Generating Triangular Mesh Nodes"<<endl;
      regularRectangularToTriangularMeshNodes(nx,ny,1.0,1.0,mesh);
    }
  else
    {
      cout<<"Generating Tetrahedral Mesh Elements"<<endl;
      regularHexahedralToTetrahedralMeshElements(nx,ny,nz,mesh);
      cout<<"Generating Tetrahedral Mesh Nodes"<<endl;
      regularHexahedralToTetrahedralMeshNodes(nx,ny,nz,1.0,1.0,1.0,mesh);
    }
  stop = CurrentTime();
  cout<<"Elapsed time for mesh generation = "<<(stop-start)<<"s"<<endl;
  cout<<"nElements_global = "<<mesh.nElements_global<<endl;
  cout<<"nNodes_global = "<<mesh.nNodes_global<<endl;

  
  cout<<"Writing mesh"<<endl;
  start = CurrentTime();
  if(outfilename != "n")
    {
      if(outfilename != "cout")
        {
          ofstream outfile(outfilename.c_str());
          writeElements(outfile,mesh);
          writeNodes(outfile,mesh);
        }
      else
        {
          writeElements(cout,mesh);
          writeNodes(cout,mesh);
        }
    }
  stop = CurrentTime();
  cout<<"Elapsed time for writing mesh = "<<(stop-start)<<"s"<<endl;
  delete [] mesh.elementNodesArray;

  cout<<"Testing readElements"<<endl;
  cout<<"Enter mesh filename"<<endl;
  cin>>infilename;
  cout<<"Reading mesh from "<<infilename<<endl;
  start = CurrentTime();
  ifstream infile(infilename.c_str());
  readElements(infile,mesh);
  stop = CurrentTime();
  cout<<"Elapsed time for  reading mesh = "<<(stop-start)<<"s"<<endl;


  cout<<"Testing global refinement"<<endl;
  cout<<"Enter nLevels"<<endl;
  int nLevels;
  cin>>nLevels;
  start = CurrentTime();
  MultilevelMesh multilevelMesh;
  if (ny == 1 && nz==1)
    globallyRefineEdgeMesh(nLevels,mesh,multilevelMesh);
  else if (nz == 1)
    globallyRefineTriangularMesh(nLevels,mesh,multilevelMesh);
  else
    globallyRefineTetrahedralMesh(nLevels,mesh,multilevelMesh);
  stop = CurrentTime();
  cout<<"Elapsed time for global refinements = "<<(stop-start)<<"s"<<endl;
  mesh = multilevelMesh.meshArray[nLevels-1];
  constructElementBoundaryElementsArray_tetrahedron(mesh);
  computeGeometricInfo_tetrahedron(mesh);
  Mesh edgeMesh;
  initializeMesh(edgeMesh);
  cout<<"nElements_global = "<<mesh.nElements_global<<endl;
  cout<<"nNodes_global = "<<mesh.nNodes_global<<endl;
  
  cout<<"Constructing element boundary element arrays"<<endl;
  if (ny == 1 && nz == 1)
    {
      constructElementBoundaryElementsArray_edge(mesh);
      edgeMesh = mesh;
    }
  else if (nz == 1)
    {
      constructElementBoundaryElementsArray_triangle(mesh);
      //set up the edge mesh, need to make sure I've set everything
      edgeMesh.nNodes_global = mesh.nNodes_global;
      edgeMesh.nodeArray = mesh.nodeArray;
      edgeMesh.nElements_global = mesh.nElementBoundaries_global;
      edgeMesh.elementNodesArray = mesh.elementBoundaryNodesArray;
      edgeMesh.nNodes_element = mesh.nNodes_elementBoundary;
    }
  else
    {
      constructElementBoundaryElementsArray_tetrahedron(mesh);
      computeGeometricInfo_tetrahedron(mesh);
      //set up the face and edge meshes
      Mesh faceMesh;
      initializeMesh(faceMesh);
      faceMesh.nNodes_global = mesh.nNodes_global;
      faceMesh.nodeArray = mesh.nodeArray;
      faceMesh.nElements_global = mesh.nElementBoundaries_global;
      faceMesh.elementNodesArray = mesh.elementBoundaryNodesArray;
      faceMesh.nNodes_element = mesh.nNodes_elementBoundary;
      constructElementBoundaryElementsArray_triangle(faceMesh);
      edgeMesh.nNodes_global = faceMesh.nNodes_global;
      edgeMesh.nodeArray = faceMesh.nodeArray;
      edgeMesh.nElements_global = faceMesh.nElementBoundaries_global;
      edgeMesh.elementNodesArray = faceMesh.elementBoundaryNodesArray;
      edgeMesh.nNodes_element = faceMesh.nNodes_elementBoundary;
    }
  cout<<"Print mesh? (n/filename/cout)"<<endl;
  cin>>outfilename;
  start = CurrentTime();
  if(outfilename != "n")
    {
      if(outfilename != "cout")
        {
          ofstream outfile(outfilename.c_str());
          writeElements(outfile,mesh);
          writeNodes(outfile,mesh);
        }
      else
        {
          writeElements(cout,mesh);
          writeNodes(cout,mesh);
        }
    }
  stop = CurrentTime();
  cout<<"Elapsed time for writing mesh = "<<(stop-start)<<"s"<<endl;

  ofstream datFile("edgeMesh.dat");
  for(int eN=0;eN<edgeMesh.nElements_global;eN++)
    {
      datFile<<scientific<<setprecision(8)<<setw(16)<<edgeMesh.nodeArray[edgeMesh.elementNodesArray[eN*2+0]*3+0]
             <<scientific<<setprecision(8)<<setw(16)<<edgeMesh.nodeArray[edgeMesh.elementNodesArray[eN*2+0]*3+1]
             <<scientific<<setprecision(8)<<setw(16)<<edgeMesh.nodeArray[edgeMesh.elementNodesArray[eN*2+0]*3+2]<<endl
             <<scientific<<setprecision(8)<<setw(16)<<edgeMesh.nodeArray[edgeMesh.elementNodesArray[eN*2+1]*3+0]
             <<scientific<<setprecision(8)<<setw(16)<<edgeMesh.nodeArray[edgeMesh.elementNodesArray[eN*2+1]*3+1]
             <<scientific<<setprecision(8)<<setw(16)<<edgeMesh.nodeArray[edgeMesh.elementNodesArray[eN*2+1]*3+2]<<endl<<endl<<endl;
    }
  //deleteMesh(edgeMesh); mesh.nodeArray = NULL;//the edgeMesh shares the nodes
  for (int i=0;i<multilevelMesh.nLevels;i++)
    deleteMesh(multilevelMesh.meshArray[i]);
  return 0;
}

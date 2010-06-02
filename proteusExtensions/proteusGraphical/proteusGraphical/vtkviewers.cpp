#include "vtkviewers.h"
#include "vtkUnstructuredGrid.h"
#include "vtkCellType.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "stdio.h"

vtkDoubleArray* vtkPrepareScalarValueArray(int nNodes,
					   const double* scalarArray)
  
{
  double* scalarArray_tmp = const_cast<double*>(scalarArray);
  assert(scalarArray_tmp);
  
  const int doNotDelete = 1;
  vtkDoubleArray *scalars = vtkDoubleArray::New();
  scalars->SetNumberOfComponents(1);
  scalars->SetNumberOfValues(nNodes);
  scalars->SetArray(scalarArray_tmp,nNodes,doNotDelete);
  return scalars;
}

vtkDoubleArray* vtkPrepareVectorValueArray(int nNodes, 
					   int nTuples,
					   const double* vectorArray)
  
{
  double* vectorArray_tmp = const_cast<double *>(vectorArray);
  assert(vectorArray_tmp);
  
  const int doNotDelete = 1;
  vtkDoubleArray* vectors = vtkDoubleArray::New();
  vectors->SetNumberOfComponents(nTuples);
  vectors->SetNumberOfTuples(nNodes);
  vectors->SetArray(vectorArray_tmp,nNodes,doNotDelete);
  return vectors;
}

vtkPoints* vtkPrepareVTKPoints(int nNodes, const double *nodeArray)
{
  using namespace std;
  vtkPoints* points = vtkPoints::New();
  points->SetDataTypeToDouble();
  
  double* nodeArray_tmp = const_cast<double *>(nodeArray);
  assert(nodeArray_tmp);
  const int doNotDelete = 1;
  vtkDoubleArray* pcoords = vtkDoubleArray::New();
  pcoords->SetNumberOfComponents(3);
  pcoords->SetNumberOfValues(nNodes);
  pcoords->SetArray(nodeArray_tmp,nNodes*3,doNotDelete);
  points->SetData(pcoords);
  return points;
}

vtkUnstructuredGrid* vtkUnstructuredGridFromMesh(int nNodes,
						 int nElements,
						 int nSimplex,
						 const double* nodeArray,
						 const int* elementNodesArray)
  
{
  using namespace std;
  vtkPoints * points = vtkPoints::New();
  points->SetDataTypeToDouble();
  
#ifdef SET_VTK_POINTS_DIRECTLY
  points->SetNumberOfPoints(nNodes);
  for (int nN = 0; nN < nNodes; nN++)
    {
      points->InsertPoint(nN,nodeArray + nN*3);
    }
#else
  double * nodeArray_tmp = const_cast<double *>(nodeArray);
  assert(nodeArray_tmp);
  const int doNotDelete = 1;
  vtkDoubleArray * pcoords = vtkDoubleArray::New();
  pcoords->SetNumberOfComponents(3);
  pcoords->SetNumberOfValues(nNodes);
  pcoords->SetArray(nodeArray_tmp,nNodes*3,doNotDelete);
  points->SetData(pcoords);
#endif
  
  // Create the grid
  vtkUnstructuredGrid * vtkMesh = vtkUnstructuredGrid::New();
  vtkMesh->Allocate(nElements);
  vtkIdType tmp[4] = {0,0,0,0}; //can't convert from int?
  
  // Set the cell type
  int cellType = VTK_LINE;
  if (nSimplex == 3)
    {
      cellType = VTK_TRIANGLE;
    }
  else if (nSimplex == 4)
    {
      cellType = VTK_TETRA;
    }
  
  // Create the cells
  for (int eN = 0; eN < nElements; eN++)
    {
      for (int nN_local = 0; nN_local < nSimplex; nN_local++)
	{
	  tmp[nN_local] = elementNodesArray[eN*nSimplex+nN_local];
	}
      vtkMesh->InsertNextCell(cellType,nSimplex,tmp);
    }
  vtkMesh->SetPoints(points);
  points->Delete();
  
  return vtkMesh;
}

void BuildQuadraticTriangleNodeList(int numNodes, 
				    int eN, 
				    const int*l2g, 
				    const int* edgeNodesArray, 
				    vtkIdType* NodeList)
{
  int base = 0;
  int n0, n1, edgeN;
  
  //Build the node list
  for  (int i=0; i<3; i++)
    {
      NodeList[i] = l2g[eN*6+i];
    }
  //mwf looks like vtk numbers edges 
  //[node,node] --> local edge id
  //[0,1] --> 0, [1,2] --> 1 [2,0] --> 2
  //which is what's used in proteus for quadratic triangles
  for (int i=0; i<3; i++)
    {
      int iproteus = i;
      NodeList[3+i] = base + l2g[eN*6+3+iproteus];
    }
  //mwf orig doesn't work for parallel
//   //Brute force way to make sure the remaining nodes are in the right order to describe the element
//   for (int i=0; i<3; i++)
//     {
//       edgeN = l2g[eN*6+3+i] - numNodes;
//       n0 = edgeNodesArray[edgeN*2+0];
//       n1 = edgeNodesArray[edgeN*2+1];
      
//       if (n0 == l2g[eN*6+0])
// 	{
// 	  if (n1 == l2g[eN*6+1])
// 	    {
// 	      NodeList[3+0]=base+l2g[eN*6+i+3];
// 	    }
// 	  if (n1 == l2g[eN*6+2])
// 	    {
// 	      NodeList[3+2]=base+l2g[eN*6+i+3];
// 	    }
// 	}
//       else if (n1 == l2g[eN*6+0])
// 	{
// 	  if (n0 == l2g[eN*6+1])
// 	    {
// 	      NodeList[3+0]=base+l2g[eN*6+i+3];
// 	    }
// 	  if (n0 == l2g[eN*6+2])
// 	    {
// 	      NodeList[3+2]=base+l2g[eN*6+i+3];
// 	    }
// 	}
//       else if (n0 == l2g[eN*6+1])
// 	{
// 	  if (n1 == l2g[eN*6+2])
// 	    {
// 	      NodeList[3+1]=base+l2g[eN*6+i+3];
// 	    }
// 	}
//       else if (n1 == l2g[eN*6+1])
// 	{
// 	  if (n0 == l2g[eN*6+2])
// 	    {
// 	      NodeList[3+1]=base+l2g[eN*6+i+3];
// 	    }
// 	}
//       else
// 	{
// 	  printf("Element Sorting Fell Through on element %i side %i\n", eN, i);
// 	}
//     }
}

vtkUnstructuredGrid* vtkUnstructuredGridFromQuadraticTriangleMesh(int nNodes,
								  int nNodesTotal,
								  int nElements,
								  const double* nodeArray,
								  const int* l2g,
								  const int* edgeNodesArray)
{
  //Build the points list
  vtkPoints *points = vtkPrepareVTKPoints(nNodesTotal, nodeArray);
  
  //Build the grid
  int cellType = VTK_QUADRATIC_TRIANGLE;
  
  vtkUnstructuredGrid * vtkMesh = vtkUnstructuredGrid::New();
  printf("Number of elements: %i\n", nElements);
  vtkMesh->Allocate(nElements);
  vtkIdType tmp[6] = {0,0,0,0,0,0}; //can't convert from int?
  
  for (int i=0; i<nElements; i++)
    {
      BuildQuadraticTriangleNodeList(nNodes, i, l2g, edgeNodesArray, tmp);
      vtkMesh->InsertNextCell(cellType, 6, tmp);
    }
  vtkMesh->SetPoints(points);
  points->Delete();
  
  return vtkMesh;
}


void BuildQuadraticTetNodeList(int numNodes,
			       int eN,
			       const int *l2g,
			       const int* edgeNodesArray,
			       vtkIdType* NodeList)
{
  int base = 0;
  int n0, n1, edgeN;
  
  //Build the node list
  for  (int i=0; i<4; i++)
    {
      NodeList[i] = l2g[eN*10+i];
    }
  //mwf looks like vtk is expecting
  //[node,node] --> local edge id to be
  //[0,1] --> 0, [0,2] --> 2, [0,3] --> 3, [1,2] --> 1, [1,3] --> 4, [2,3]-->5
  //proteus stores 3d edge dof as
  //|(n0,n1),(n1,n2),(n2,n3)|(n0,n2),(n1,n3)|(n0,n3)|
  //[0,1] --> 0, [1,2] --> 1, [2,3] --> 2, [0,2] --> 3, [1,3] --> 4, [0,3] --> 5
  NodeList[4+0] = l2g[eN*10+4+0]; 
  NodeList[4+1] = l2g[eN*10+4+1];
  NodeList[4+2] = l2g[eN*10+4+3];
  NodeList[4+3] = l2g[eN*10+4+5];
  NodeList[4+4] = l2g[eN*10+4+4];
  NodeList[4+5] = l2g[eN*10+4+2];
  //mwf doesn't work in parallel any more
//   //brute force way to make sure the remaining nodes are in the right order to describe the element
//   for (int i=0; i<6; i++)
//     {
//       edgeN = l2g[eN*10+4+i] - numNodes;
//       n0 = edgeNodesArray[edgeN*2+0];
//       n1 = edgeNodesArray[edgeN*2+1];
      
//       if (n0 == l2g[eN*10+0])
// 	{
// 	  if (n1 == l2g[eN*10+1])
// 	    {
// 	      NodeList[4+0]=base+l2g[eN*10+i+4];
// 	    }
// 	  if (n1 == l2g[eN*10+2])
// 	    {
// 	      NodeList[4+2]=base+l2g[eN*10+i+4];
// 	    }
// 	  if (n1 == l2g[eN*10+3])
// 	    {
// 	      NodeList[4+3]=base+l2g[eN*10+i+4];
// 	    }
// 	}
//       else if (n1 == l2g[eN*10+0])
// 	{
// 	  if (n0 == l2g[eN*10+1])
// 	    {
// 	      NodeList[4+0]=base+l2g[eN*10+i+4];
// 	    }
// 	  if (n0 == l2g[eN*10+2])
// 	    {
// 	      NodeList[4+2]=base+l2g[eN*10+i+4];
// 	    }
// 	  if (n0 == l2g[eN*10+3])
// 	    {
// 	      NodeList[4+3]=base+l2g[eN*10+i+4];
// 	    }
// 	}
//       else if (n0 == l2g[eN*10+1])
// 	{
// 	  if (n1 == l2g[eN*10+2])
// 	    {
// 	      NodeList[4+1]=base+l2g[eN*10+i+4];
// 	    }
// 	  if (n1 == l2g[eN*10+3])
// 	    {
// 	      NodeList[4+4]=base+l2g[eN*10+i+4];
// 	    }
// 	}
//       else if (n1 == l2g[eN*10+1])
// 	{
// 	  if (n0 == l2g[eN*10+2])
// 	    {
// 	      NodeList[4+1]=base+l2g[eN*10+i+4];
// 	    }
// 	  if (n0 == l2g[eN*10+3])
// 	    {
// 	      NodeList[4+4]=base+l2g[eN*10+i+4];
// 	    }
// 	}
//       else if (n0 == l2g[eN*10+2])
// 	{
// 	  if (n1 == l2g[eN*10+3])
// 	    {
// 	      NodeList[4+5]=base+l2g[eN*10+i+4];
// 	    }
// 	}
//       else if (n1 == l2g[eN*10+2])
// 	{
// 	  if (n0 == l2g[eN*10+3])
// 	    {
// 	      NodeList[4+5]=base+l2g[eN*10+i+4];
// 	    }
// 	}
//       else
// 	{
// 	  printf("Element Sorting Fell Through on element %i side %i\n", eN, i);
// 	}
//     }
}

vtkUnstructuredGrid* vtkUnstructuredGridFromQuadraticTetMesh(int nNodes,
							     int nNodesTotal,
							     int nElements,
							     const double * nodeArray,
							     const int * l2g,
							     const int * edgeNodesArray)
{
  //Build the points list
  vtkPoints *points = vtkPrepareVTKPoints(nNodesTotal, nodeArray);
  
  //Build the grid
  int cellType = VTK_QUADRATIC_TETRA;
  
  vtkUnstructuredGrid * vtkMesh = vtkUnstructuredGrid::New();
  printf("Number of elements: %i\n", nElements);
  vtkMesh->Allocate(nElements);
  vtkIdType tmp[10] = {0,0,0,0,0,0,0,0,0,0}; //can't convert from int?
  
  for (int i=0; i<nElements; i++)
    {
      BuildQuadraticTetNodeList(nNodes, i, l2g, edgeNodesArray, tmp);
      vtkMesh->InsertNextCell(cellType, 10, tmp);
    }
  vtkMesh->SetPoints(points);
  points->Delete();
  
  return vtkMesh;
}	

vtkPolyData* vtkPolyDataBoundaryMesh(int nNodes,
				     int nElementBoundaries,
				     int nSimplexBoundary,
				     const double * nodeArray,
				     const int * elementBoundaryNodesArray,
				     const int * elementBoundaryMaterialTypes)
{
  using namespace std;
  vtkPoints * points = vtkPoints::New();
  points->SetDataTypeToDouble();
  const int doNotDelete = 1;
#ifdef SET_VTK_POINTS_DIRECTLY
  points->SetNumberOfPoints(nNodes);
  for (int nN = 0; nN < nNodes; nN++)
    {
      points->InsertPoint(nN,nodeArray + nN*3); 
    }
#else
  double * nodeArray_tmp = const_cast<double *>(nodeArray);
  assert(nodeArray_tmp);
  vtkDoubleArray * pcoords = vtkDoubleArray::New();
  pcoords->SetNumberOfComponents(3);
  pcoords->SetNumberOfValues(nNodes);
  pcoords->SetArray(nodeArray_tmp,nNodes*3,doNotDelete);
  points->SetData(pcoords);
#endif
  vtkPolyData * vtkBoundaryMesh = vtkPolyData::New();
  vtkBoundaryMesh->Allocate(nElementBoundaries);
  vtkBoundaryMesh->SetPoints(points);
  vtkIdType tmp[3] = {0,0,0}; //can't convert from int?
  assert(nSimplexBoundary >= 2);
  int cellType = VTK_LINE;
  if (nSimplexBoundary == 3)
    {
      cellType = VTK_TRIANGLE;
    }
  for (int ebN = 0; ebN < nElementBoundaries; ebN++)
    {
      for (int nN_local = 0; nN_local < nSimplexBoundary; nN_local++)
	tmp[nN_local] =elementBoundaryNodesArray[ebN*nSimplexBoundary+nN_local];
      vtkBoundaryMesh->InsertNextCell(cellType,nSimplexBoundary,tmp);
    }
  if (elementBoundaryMaterialTypes != 0)
    {
      int * elementBoundaryMaterialTypes_tmp = const_cast<int*>(elementBoundaryMaterialTypes);
      assert(elementBoundaryMaterialTypes_tmp);
      vtkIntArray * pflags = vtkIntArray::New();
      pflags->SetNumberOfComponents(1);
      pflags->SetNumberOfValues(nElementBoundaries);
      pflags->SetArray(elementBoundaryMaterialTypes_tmp,nElementBoundaries,doNotDelete);
      vtkBoundaryMesh->GetCellData()->SetScalars(pflags);
    }
  points->Delete();
  
  return vtkBoundaryMesh;
}

bool meshElementAndNodeArraysFromVTKUnstructuredGrid(vtkUnstructuredGrid* vtkMesh,
						     int& nNodes,
						     int& nElements,
						     int& nSimplex,
						     double*& nodeArray,
						     int*& elementNodesArray,
						     int*& elementMaterialTypes,
						     int*& nodeMaterialTypes,
						     const bool& readMaterialIds,
						     const int& defaultElementMaterialType,
						     const int& defaultNodeMaterialType)
{
  bool failed = false;
  assert(vtkMesh);
  assert(vtkMesh->IsHomogeneous());
  
  vtkIdType vtkCellId = vtkMesh->GetCellType(0);
  nSimplex = -1;
  if (vtkCellId == VTK_LINE)
    nSimplex = 2;
  else if (vtkCellId == VTK_TRIANGLE)
    nSimplex = 3;
  else if (vtkCellId == VTK_TETRA)
    nSimplex = 4;
  assert(nSimplex > 1);
  
  nElements = vtkMesh->GetNumberOfCells();
  nNodes    = vtkMesh->GetNumberOfPoints();
  //mwf debug
  std::cout<<"vtkViewers getMeshFromVTKUnst vtkCellId = "<<vtkCellId
	   <<" nSimplex = "<<nSimplex<<std::endl
	   <<"nElements = "<<nElements<<" nNodes= "<<nNodes<<std::endl; 
  if (nodeArray)
    delete [] nodeArray;
  nodeArray = new double[nNodes*3];
  if (elementNodesArray)
    delete [] elementNodesArray;
  elementNodesArray = new int[nElements*nSimplex];
  
  for (int nN = 0; nN < nNodes; nN++)
    {
      vtkMesh->GetPoint(nN,nodeArray + nN*3);
    }
  vtkIdType nPerCell(0), *pcell(0);
  for (int eN= 0; eN< nElements; eN++)
    {
      vtkMesh->GetCellPoints(eN,nPerCell,pcell);
      for (int nN_local = 0; nN_local < nSimplex; nN_local++)
	elementNodesArray[eN*nSimplex + nN_local] = pcell[nN_local];
    }
  
  if (elementMaterialTypes)
    delete [] elementMaterialTypes;
  elementMaterialTypes = new int[nElements];
  if (nodeMaterialTypes)
    delete [] nodeMaterialTypes;
  nodeMaterialTypes = new int[nNodes];
  
  //skip reading material types temporarily
  memset(elementMaterialTypes,defaultElementMaterialType,nElements*sizeof(int));
  memset(nodeMaterialTypes,defaultNodeMaterialType,nNodes*sizeof(int));
  
  return failed;
}


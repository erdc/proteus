#ifndef PROTEUS_VTK_VIEWERS_H
#define PROTEUS_VTK_VIEWERS_H
#include <cassert>
class  vtkUnstructuredGrid;
class  vtkDoubleArray;
class  vtkPoints;

vtkDoubleArray* vtkPrepareScalarValueArray(int nNodes,
					   const double * scalarArray);

vtkDoubleArray* vtkPrepareVectorValueArray(int nNodes, int nTouples,
										   const double * vectorArray);

vtkPoints* vtkPrepareVTKPoints(int nNodes, const double *nodeArray);

vtkUnstructuredGrid* vtkUnstructuredGridFromMesh(int nNodes,
						 int nElements,
						 int nSimplex,
						 const double * nodeArray,
						 const int * elementNodesArray,
						 const int* elementMaterialTypes=0,
						 const int* nodeMaterialTypes=0);

vtkUnstructuredGrid* vtkUnstructuredGridFromQuadraticTriangleMesh(int nNodes,
																  int nNodesTotal,
																  int nElements,
																  const double * nodeArray,
																  const int * l2g,
																  const int * edgeNodesArray);

vtkUnstructuredGrid* vtkUnstructuredGridFromQuadraticTetMesh(int nNodes,
															 int nNodesTotal,
															 int nElements,
															 const double * nodeArray,
															 const int * l2g,
															 const int * edgeNodesArray);

class vtkPolyData;

vtkPolyData* vtkPolyDataBoundaryMesh(int nNodes,
				     int nElementBoundaries,
				     int nBoundarySimplex,
				     const double * nodeArray,
				     const int * elementBoundaryNodesArray,
				     const int * elementBoundaryMaterialTypes);

bool meshElementAndNodeArraysFromVTKUnstructuredGrid(vtkUnstructuredGrid* vtkMesh,
						     int& nNodes,
						     int& nElements,
						     int& nSimplex,
						     double *& nodeArray,
						     int *& elementNodesArray,
						     int *& elementMaterialTypes,
						     int *& nodeMaterialTypes,
						     const bool& readMaterialIds = false,
						     const int& defaultElementMaterialType=0,
						     const int& defaultNodeMaterialType=0);

bool classifyElementMaterialPropertiesFromVTKUnstructuredGridSolid(vtkUnstructuredGrid* vtkSolid,
								   //properties of the mesh classifying
								   int nElements,
								   int newMaterialId,
								   const double* elementBarycentersArray,
								   const double* elementDiametersArray,
								   int * elementMaterialTypes,
								   double tol_hFactor,
								   int verbose =0);
bool classifyElementMaterialPropertiesFromVTKUnstructuredGridNeighborhood(vtkUnstructuredGrid* vtkSolid,
									  //properties of the mesh classifying
									  int nElements,
									  int newMaterialId,
									  const double* elementBarycentersArray,
									  const double* elementDiametersArray,
									  int * elementMaterialTypes,
									  int useAbsoluteTolerance,
									  double tolerance,
									  int verbose);
#endif

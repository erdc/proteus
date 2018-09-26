#include "DumpMesh.h"
#include <PCU.h>

/* see flcbdfWrappersModule.cpp:
   flcbdfWrappersConvertPUMIPartitionToPython
   and cmeshToolsModule.cpp:
   cmeshToolsBuildPythonMeshInterface */

static void print_int_2d(int* a, int h, int w, FILE* f)
{
  for (int i = 0; i < h; ++i) {
    for (int j = 0; j < w; ++j)
      fprintf(f, "%d ", a[i * w + j]);
    fprintf(f, "\n");
  }
}

static void print_double_2d(double* a, int h, int w, FILE* f)
{
  for (int i = 0; i < h; ++i) {
    for (int j = 0; j < w; ++j)
      fprintf(f, "%f ", a[i * w + j]);
    fprintf(f, "\n");
  }
}

static void dump_proteus_mesh_header(Mesh* m, FILE* f)
{
  fprintf(f, "nElements_global = %d\n", m->nElements_global);
  fprintf(f, "nNodes_global = %d\n", m->nNodes_global);
  fprintf(f, "nNodes_element = %d\n", m->nNodes_element);
  fprintf(f, "nNodes_elementBoundary = %d\n", m->nNodes_elementBoundary);
  fprintf(f, "nElementBoundaries_element = %d\n", m->nElementBoundaries_element);
  fprintf(f, "nElementBoundaries_global = %d\n", m->nElementBoundaries_global);
  fprintf(f, "nInteriorElementBoundaries_global = %d\n",
      m->nInteriorElementBoundaries_global);
  fprintf(f, "nExteriorElementBoundaries_global = %d\n",
      m->nExteriorElementBoundaries_global);
  fprintf(f, "max_nElements_node = %d\n", m->max_nElements_node);
  fprintf(f, "nEdges_global = %d\n", m->nEdges_global);
  fprintf(f, "max_nNodeNeighbors_node = %d\n", m->max_nNodeNeighbors_node);
}

static void dump_proteus_subdomain(Mesh* m, FILE* f)
{
  dump_proteus_mesh_header(m, f);
  fprintf(f, "elementNodesArray:\n");
  print_int_2d(m->elementNodesArray, m->nElements_global,
               m->nNodes_element, f);
  fprintf(f, "nodeElementOffsets:\n");
  if (!m->nodeElementOffsets)
    fprintf(f, " NULL\n");
  else {
    for (int i = 0; i <= m->nNodes_global; ++i)
      fprintf(f, "%d\n", m->nodeElementOffsets[i]);
    int total = m->nodeElementOffsets[m->nNodes_global];
    assert(total == m->nElements_global * m->nNodes_element);
    fprintf(f, "nodeElementsArray:\n");
    for (int i = 0; i < total; ++i)
      fprintf(f, "%d\n", m->nodeElementsArray[i]);
  }
  fprintf(f, "elementNeighborsArray:\n");
  print_int_2d(m->elementNeighborsArray, m->nElements_global,
               m->nElementBoundaries_element, f);
  fprintf(f, "elementBoundariesArray:\n");
  print_int_2d(m->elementBoundariesArray, m->nElements_global,
               m->nElementBoundaries_element, f);
  fprintf(f, "elementBoundaryNodesArray:\n");
  print_int_2d(m->elementBoundaryNodesArray, m->nElementBoundaries_global,
               m->nNodes_elementBoundary, f);
  fprintf(f, "elementBoundaryElementsArray:\n");
  print_int_2d(m->elementBoundaryElementsArray, m->nElementBoundaries_global,
               2, f);
  fprintf(f, "elementBoundaryLocalElementBoundariesArray:\n");
  print_int_2d(m->elementBoundaryLocalElementBoundariesArray,
               m->nElementBoundaries_global, 2, f);
  fprintf(f, "interiorElementBoundariesArray:\n");
  for (int i = 0; i < m->nInteriorElementBoundaries_global; ++i)
    fprintf(f, "%d\n", m->interiorElementBoundariesArray[i]);
  fprintf(f, "exteriorElementBoundariesArray:\n");
  for (int i = 0; i < m->nExteriorElementBoundaries_global; ++i)
    fprintf(f, "%d\n", m->exteriorElementBoundariesArray[i]);
  fprintf(f, "edgeNodesArray:\n");
  print_int_2d(m->edgeNodesArray, m->nEdges_global, 2, f);
  fprintf(f, "nodeStarOffsets:\n");
  if (!m->nodeStarOffsets)
    fprintf(f, " NULL\n");
  else {
    for (int i = 0; i <= m->nNodes_global; ++i)
      fprintf(f, "%d\n", m->nodeStarOffsets[i]);
    int total = m->nodeStarOffsets[m->nNodes_global];
    assert(total == m->nEdges_global * 2);
    fprintf(f, "nodeStarArray:\n");
    for (int i = 0; i < total; ++i)
      fprintf(f, "%d\n", m->nodeStarArray[i]);
  }
  fprintf(f, "elementMaterialTypes:\n");
  for (int i = 0; i < m->nElements_global; ++i)
    fprintf(f, "%d\n", m->elementMaterialTypes[i]);
  fprintf(f, "elementBoundaryMaterialTypes:\n");
  for (int i = 0; i < m->nElementBoundaries_global; ++i)
    fprintf(f, "%d\n", m->elementBoundaryMaterialTypes[i]);
  fprintf(f, "nodeMaterialTypes:\n");
  for (int i = 0; i < m->nNodes_global; ++i)
    fprintf(f, "%d\n", m->nodeMaterialTypes[i]);
  fprintf(f, "nodeArray:\n");
  print_double_2d(m->nodeArray, m->nNodes_global, 3, f);
  fprintf(f, "elementDiametersArray:\n");
  for (int i = 0; i < m->nElements_global; ++i)
    fprintf(f, "%f\n", m->elementDiametersArray[i]);
  fprintf(f, "elementInnerDiametersArray:\n");
  for (int i = 0; i < m->nElements_global; ++i)
    fprintf(f, "%f\n", m->elementInnerDiametersArray[i]);
  fprintf(f, "elementBoundaryDiametersArray:\n");
  for (int i = 0; i < m->nElementBoundaries_global; ++i)
    fprintf(f, "%f\n", m->elementBoundaryDiametersArray[i]);
  fprintf(f, "elementBarycentersArray:\n");
  print_double_2d(m->elementBarycentersArray, m->nElements_global, 3, f);
  print_double_2d(m->elementBoundaryBarycentersArray,
      m->nElementBoundaries_global, 3, f);
  fprintf(f, "nodeDiametersArray:\n");
  for (int i = 0; i < m->nNodes_global; ++i)
    fprintf(f, "%f\n", m->nodeDiametersArray[i]);
  fprintf(f, "nodeSupportArray:\n");
  for (int i = 0; i < m->nNodes_global; ++i)
    fprintf(f, "%f\n", m->nodeSupportArray[i]);
  fprintf(f, "h = %f\n", m->h);
  fprintf(f, "hMin = %f\n", m->hMin);
  fprintf(f, "sigmaMax = %f\n", m->sigmaMax);
  fprintf(f, "volume = %f\n", m->volume);
}

static void dump_proteus_parallel_arrays(Mesh* m, FILE* f)
{
  int size = PCU_Comm_Peers();
  fprintf(f, "elementOffsets_subdomain_owned:\n");
  for (int i = 0; i <= size; ++i)
    fprintf(f, "%d\n", m->elementOffsets_subdomain_owned[i]);
  fprintf(f, "elementNumbering_subdomain2global:\n");
  for (int i = 0; i < m->subdomainp->nElements_global; ++i)
    fprintf(f, "%d\n", m->elementNumbering_subdomain2global[i]);
  fprintf(f, "nodeOffsets_subdomain_owned:\n");
  for (int i = 0; i <= size; ++i)
    fprintf(f, "%d\n", m->nodeOffsets_subdomain_owned[i]);
  fprintf(f, "nodeNumbering_subdomain2global:\n");
  for (int i = 0; i < m->subdomainp->nNodes_global; ++i)
    fprintf(f, "%d\n", m->nodeNumbering_subdomain2global[i]);
  fprintf(f, "elementBoundaryOffsets_subdomain_owned:\n");
  for (int i = 0; i <= size; ++i)
    fprintf(f, "%d\n", m->elementBoundaryOffsets_subdomain_owned[i]);
  fprintf(f, "elementBoundaryNumbering_subdomain2global:\n");
  for (int i = 0; i < m->subdomainp->nElementBoundaries_global; ++i)
    fprintf(f, "%d\n", m->elementBoundaryNumbering_subdomain2global[i]);
  fprintf(f, "edgeOffsets_subdomain_owned:\n");
  for (int i = 0; i <= size; ++i)
    fprintf(f, "%d\n", m->edgeOffsets_subdomain_owned[i]);
  fprintf(f, "edgeNumbering_subdomain2global:\n");
  for (int i = 0; i < m->subdomainp->nEdges_global; ++i)
    fprintf(f, "%d\n", m->edgeNumbering_subdomain2global[i]);

}

void dump_proteus_mesh(Mesh* m, FILE* f)
{
  fprintf(stderr,"%d dumping mesh file\n", PCU_Comm_Self());
  if (PCU_Comm_Peers() > 1) {
    fprintf(f, "PARALLEL HEADER\n");
    dump_proteus_mesh_header(m, f);
    fprintf(f, "PARALLEL ARRAYS\n");
    dump_proteus_parallel_arrays(m, f);
    fprintf(f, "SUBDOMAIN\n");
    dump_proteus_subdomain(m->subdomainp, f);
  } else {
    dump_proteus_subdomain(m, f);
  }
}

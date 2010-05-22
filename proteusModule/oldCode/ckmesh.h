struct Mesh
{
  int nNodes_global,
    nElements_global,
    nNodes_element,
    nNodes_elementBoundary,
    nElementBoundaries_element,
    nElementBoundaries_global,
	nInteriorElementBoundaries_global,
    nExteriorElementBoundaries_global,
    max_nElements_node;

  double *nodeArray;

  int *nodes_element,
    *elements_node,
    *elementBoundaries,
    *nodes_elementBoundary,
    *elements_elementBoundary,
    *interiorElementBoundaries,
    *exteriorElementBoundaries;

/* we will want to add geometric  information to this */
};

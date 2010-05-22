typedef struct Mesh
{
int	nElements_global,
	nNodes_global,
	nNodes_element,
	nElementBoundaries_element,
	nElementBoundaries_global,
	nInteriorElementBoundaries_global,
	nExteriorElementBoundaries_global,
	max_nElements_node;

int	*elementNodesArray,
	*nodeElementsArray,
	*elementBoundaryElementsArray,
	*interiorElementBoundariesArray,
	*exteriorElementBoundariesArray;

double	*nodeArray;

/* we will want to add geometric  information to this */
};


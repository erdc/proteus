#! /usr/bin/env python
"""
A test script for MeshTools and cmeshTools
"""
import Comm
import MeshTools
import cmeshTools



if __name__ == '__main__':
    import optparse
    import sys
    if sys.version_info[1] >= 5:
        import cProfile as profiler
    else:
        import profile as profiler
    #
    import pstats
    usage = "usage: %prog [options] pFile.py [nFile.py]"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-P", "--petsc-options",
                      help="Options for  PETSc",
                      action="store",
                      type="string",
                      dest="petscOptions",
                      default=None)
    (opts,args) = parser.parse_args()
    if opts.petscOptions is not None:
        sys.argv = sys.argv[:-1]+opts.petscOptions.split()
        print(sys.argv)
    Comm.argv = sys.argv
    comm = Comm.init()
    #cmesh = cmeshTools.CMesh()
#     multilevelMesh = MeshTools.MultilevelEdgeMesh(3,1,1,
#                                                   1.0,1.0,1.0,
#                                                   5)
#     mesh = multilevelMesh.meshList[-1]
#     print 'done'

#     #mesh.buildFromC(cmesh)
#     print ("nElements_global,                              ",mesh.nElements_global,'\n',
#            "nNodes_global,                                 ",mesh.nNodes_global,'\n',
#            "nNodes_element,                                ",mesh.nNodes_element,'\n',
#            "nNodes_elementBoundary,                        ",mesh.nNodes_elementBoundary,'\n',
#            "nElementBoundaries_element,                    ",mesh.nElementBoundaries_element,'\n',
#            "nElementBoundaries_global,               ",mesh.nElementBoundaries_global,'\n',
#            "nInteriorElementBoundaries_global,             ",mesh.nInteriorElementBoundaries_global,'\n',
#            "nExteriorElementBoundaries_global,             ",mesh.nExteriorElementBoundaries_global,'\n',
#            "max_nElements_node,                            ",mesh.max_nElements_node,'\n',
#            "elementNodesArray,                             ",mesh.elementNodesArray,'\n',
#            "nodeElementsArray,                             ",mesh.nodeElementsArray,'\n',
#            "elementNeighborsArray,                         ",mesh.elementNeighborsArray,'\n',
#            "elementBoundaryNodesArray,                     ",mesh.elementBoundaryNodesArray,'\n',
#            "elementBoundaryElementsArray,                  ",mesh.elementBoundaryElementsArray,'\n',
#            "elementBoundaryLocalElementBoundariesArray,    ",mesh.elementBoundaryLocalElementBoundariesArray,'\n',
#            "interiorElementBoundariesArray,                ",mesh.interiorElementBoundariesArray,'\n',
#            "exteriorElementBoundariesArray,                ",mesh.exteriorElementBoundariesArray,'\n',
#            "elementMaterialTypes,                          ",mesh.elementMaterialTypes,'\n',
#            "nodeArray                                      ",mesh.nodeArray)

#    mesh = MeshTools.TriangularMesh()
#    mesh.generateTriangularMeshFromRectangularGrid(3,3,1.0,1.0)
#    print mesh.elementNodesArray
#    print mesh.nodeArray
#    print mesh.edgeNodesArray
#    mesh.writeEdgesGnuplot2('subdomain_meshEdges'+`comm.rank()`)
#    mesh.viewMeshGnuplotPipePar(['subdomain_meshEdges'+`comm.rank()`])
#     #mesh.buildFromC(cmesh)
#     print ("nElements_global,                              ",mesh.nElements_global,'\n',
#            "nNodes_global,                                 ",mesh.nNodes_global,'\n',
#            "nNodes_element,                                ",mesh.nNodes_element,'\n',
#            "nNodes_elementBoundary,                        ",mesh.nNodes_elementBoundary,'\n',
#            "nElementBoundaries_element,                    ",mesh.nElementBoundaries_element,'\n',
#            "nElementBoundaries_global,               ",mesh.nElementBoundaries_global,'\n',
#            "nInteriorElementBoundaries_global,             ",mesh.nInteriorElementBoundaries_global,'\n',
#            "nExteriorElementBoundaries_global,             ",mesh.nExteriorElementBoundaries_global,'\n',
#            "max_nElements_node,                            ",mesh.max_nElements_node,'\n',
#            "elementNodesArray,                             ",mesh.elementNodesArray,'\n',
#            "nodeElementsArray,                             ",mesh.nodeElementsArray,'\n',
#            "elementNeighborsArray,                         ",mesh.elementNeighborsArray,'\n',
#            "elementBoundaryNodesArray,                     ",mesh.elementBoundaryNodesArray,'\n',
#            "elementBoundaryElementsArray,                  ",mesh.elementBoundaryElementsArray,'\n',
#            "elementBoundaryLocalElementBoundariesArray,    ",mesh.elementBoundaryLocalElementBoundariesArray,'\n',
#            "interiorElementBoundariesArray,                ",mesh.interiorElementBoundariesArray,'\n',
#            "exteriorElementBoundariesArray,                ",mesh.exteriorElementBoundariesArray,'\n',
#            "elementMaterialTypes,                          ",mesh.elementMaterialTypes,'\n',
#            "nodeArray                                      ",mesh.nodeArray)
#     mesh = MeshTools.TetrahedralMesh()
#     mesh = MeshTools.TetrahedralMesh()
#     mesh.generateTetrahedralMeshFromRectangularGrid(3,3,3,1.0,1.0,1.0)
#     #mesh.buildFromC(cmesh)
#     print ("nElements_global,                              ",mesh.nElements_global,'\n',
#            "nNodes_global,                                 ",mesh.nNodes_global,'\n',
#            "nNodes_element,                                ",mesh.nNodes_element,'\n',
#            "nNodes_elementBoundary,                        ",mesh.nNodes_elementBoundary,'\n',
#            "nElementBoundaries_element,                    ",mesh.nElementBoundaries_element,'\n',
#            "nElementBoundaries_global,               ",mesh.nElementBoundaries_global,'\n',
#            "nInteriorElementBoundaries_global,             ",mesh.nInteriorElementBoundaries_global,'\n',
#            "nExteriorElementBoundaries_global,             ",mesh.nExteriorElementBoundaries_global,'\n',
#            "max_nElements_node,                            ",mesh.max_nElements_node,'\n',
#            "elementNodesArray,                             ",mesh.elementNodesArray,'\n',
#            "nodeElementsArray,                             ",mesh.nodeElementsArray,'\n',
#            "elementNeighborsArray,                         ",mesh.elementNeighborsArray,'\n',
#            "elementBoundaryNodesArray,                     ",mesh.elementBoundaryNodesArray,'\n',
#            "elementBoundaryElementsArray,                  ",mesh.elementBoundaryElementsArray,'\n',
#            "elementBoundaryLocalElementBoundariesArray,    ",mesh.elementBoundaryLocalElementBoundariesArray,'\n',
#            "interiorElementBoundariesArray,                ",mesh.interiorElementBoundariesArray,'\n',
#            "exteriorElementBoundariesArray,                ",mesh.exteriorElementBoundariesArray,'\n',
#            "elementMaterialTypes,                          ",mesh.elementMaterialTypes,'\n',
#            "nodeArray                                      ",mesh.nodeArray)
#     grid3d = MeshTools.RectangularGrid(3,3,3,1.0,1.0,1.0)
#     mesh = MeshTools.TetrahedralMesh()
#     mesh.rectangularToTetrahedral(grid3d)
#     print ("nElements_global,                              ",mesh.nElements_global,'\n',
#            "nNodes_global,                                 ",mesh.nNodes_global,'\n',
#            "nNodes_element,                                ",mesh.nNodes_element,'\n',
#            "nNodes_elementBoundary,                        ",mesh.nNodes_elementBoundary,'\n',
#            "nElementBoundaries_element,                    ",mesh.nElementBoundaries_element,'\n',
#            "nElementBoundaries_global,               ",mesh.nElementBoundaries_global,'\n',
#            "nInteriorElementBoundaries_global,             ",mesh.nInteriorElementBoundaries_global,'\n',
#            "nExteriorElementBoundaries_global,             ",mesh.nExteriorElementBoundaries_global,'\n',
#            "max_nElements_node,                            ",mesh.max_nElements_node,'\n',
#            "elementNodesArray,                             ",mesh.elementNodesArray,'\n',
#            "nodeElementsArray,                             ",mesh.nodeElementsArray,'\n',
#            "elementNeighborsArray,                         ",mesh.elementNeighborsArray,'\n',
#            "elementBoundaryNodesArray,                     ",mesh.elementBoundaryNodesArray,'\n',
#            "elementBoundaryElementsArray,                  ",mesh.elementBoundaryElementsArray,'\n',
#            "elementBoundaryLocalElementBoundariesArray,    ",mesh.elementBoundaryLocalElementBoundariesArray,'\n',
#            "interiorElementBoundariesArray,                ",mesh.interiorElementBoundariesArray,'\n',
#            "exteriorElementBoundariesArray,                ",mesh.exteriorElementBoundariesArray,'\n',
#            "elementMaterialTypes,                          ",mesh.elementMaterialTypes,'\n',
#            "nodeArray                                      ",mesh.nodeArray)
    print("Testing multilevel mesh")
    multilevelMesh = MeshTools.MultilevelTetrahedralMesh(10,10,10,1.0,1.0,1.0,1)
    print(multilevelMesh.meshList[-1].nodeArray)
    print(multilevelMesh.meshList[-1].edgeNodesArray)
#     for i in range(multilevelMesh.nLevels):
#         multilevelMesh.meshList[i].partitionMesh()
    print(multilevelMesh.meshList[-1].nodeArray)
    print(multilevelMesh.meshList[-1].edgeNodesArray)
    multilevelMesh.meshList[-1].subdomainMesh.writeEdgesGnuplot2('subdomain_meshEdges'+repr(comm.rank()))
    comm.barrier()
    if comm.isMaster():
        filenames = ['subdomain_meshEdges'+repr(i) for i in range(comm.size())]
        print(filenames)
        multilevelMesh.meshList[-1].subdomainMesh.viewMeshGnuplotPipePar(filenames)

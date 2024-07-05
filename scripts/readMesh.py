#! /usr/bin/env python
import proteus
from proteus.MeshTools import *

## \ingroup scripts
#
# \file readMesh.py
#
# \brief A script for reading a .3dm file for  viewing or exctracint the boundary mesh.

def readMesh():
    from optparse import OptionParser
    usage  = "usage: %readMesh  [options] meshFilename"
    parser = OptionParser(usage=usage)
    parser.add_option("-b", "--boundaryMesh",
                      help="print boundary mesh to 3dm file",
                      action="store_true",
                      dest="boundaryMesh",
                      default=False)
    parser.add_option("-e", "--ensight",
                      help="print meshes to ensight format",
                      action="store_true",
                      dest="ensight",
                      default=False)
    parser.add_option("-m", "--matlab",
                      help="print edges to files for plotting in matlab",
                      action="store_true",
                      dest="matlab",
                      default=False)
    parser.add_option("-M", "--view-matlab",
                      help="plot edges in matlab",
                      action="store_true",
                      dest="viewMatlab",
                      default=False)
    parser.add_option("-g", "--gnuplot",
                      help="print edges to files for plotting in gnuplot",
                      action="store_true",
                      dest="gnuplot",
                      default=False)
    parser.add_option("-G", "--view-gnuplot",
                      help="plot edges in gnuplot",
                      action="store_true",
                      dest="viewGnuplot",
                      default=False)
    (opts, args) = parser.parse_args()
    print(args)
    if len(args)==1:
        meshFilename = args[0]
    else:
        print(usage)
        return 1
    #mesh=TetrahedralMesh()
    #mesh.readMeshADH(filename)
    mesh = None
    meshIn = open(meshFilename+'.3dm','r')
    firstLine = meshIn.readline()
    firstWords = firstLine.split()
    meshIn.close()
    if firstWords[0] == 'MESH3D':
        mesh = Mesh3DM(meshFilename)
    elif firstWords[0] == 'MESH2D':
        mesh = Mesh2DM(meshFilename)
    else:
        print(firstWords[0])
    #print mesh.meshInfo()
    if opts.boundaryMesh:
        mesh.buildTriangleArrays()
        mesh.writeBoundaryMeshADH(meshFilename)
    if opts.ensight:
        if opts.boundaryMesh:
            mesh.writeBoundaryMeshEnsight(meshFilename)
        mesh.writeMeshEnsight(meshFilename)
    if opts.matlab or opts.viewMatlab:
        mesh.writeEdgesMatlab('meshEdges')
    if opts.viewMatlab:
        mesh.viewMeshMatlab('meshEdges')
    if opts.gnuplot or opts.viewGnuplot:
        mesh.writeEdgesGnuplot('meshEdges')
    if opts.viewGnuplot:
        mesh.viewMeshGnuplot('meshEdges')
    return 0
if __name__ == '__main__':
    import profile
    import pstats
    import gc
    gc.enable()
#     gc.set_debug(gc.DEBUG_STATS)
#      profile.run('readMesh()','readMeshProf')
#      p = pstats.Stats('readMeshProf')
#      p.sort_stats('cumulative').print_stats(20)
#      p.sort_stats('time').print_stats(20)
    readMesh()
    input('Please press return to continue... \n')

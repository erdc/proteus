#! /usr/bin/env python
from numpy import *
import proteus
from proteus.MeshTools import *
from proteus.FemTools import *
from string import *

## \ingroup scripts
#
# \file adh2ensight.py
#
# \brief Convert adh mesh and data files to ensight format.


def adh2ensight():
    from optparse import OptionParser
    usage = "usage: %adh2ensight [options] meshFilename ensightFilenameRoot dataFilename [dataFilename2 ...]"
    parser = OptionParser(usage=usage)
    parser.add_option("-b", "--boundaryMesh",
                      help="print boundary mesh to 3dm file",
                      action="store_true",
                      dest="boundaryMesh",
                      default=False)
    (opts, args) = parser.parse_args()
    if len(args) >= 3:
        meshFilename = args[0]
        ensightFilename = args[1]
        funFilenameList=[]
        for filename in args[2:]:
            funFilenameList.append(filename)
    else:
        print(usage)
        exit
    caseOut=open(ensightFilename+'.case','w')
    caseOut.write('FORMAT\n'+'type: ensight gold\n')
    caseOut.close()
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
    mesh.writeMeshEnsight(ensightFilename,ensightFilename)
    if opts.boundaryMesh:
        mesh.buildTriangleArrays()
        mesh.writeBoundaryMeshADH(ensightFilename)
        mesh.writeBoundaryMeshEnsight(ensightFilename)
    tList=[]
    ftypeList=[]
    descriptionList=[]
    for filename in funFilenameList:
        funIn = open(filename,'r')
        f = None
        fType= None
        for i in range(7):
            line = funIn.readline()
            print(("Line %i" % i ) + " = " + line.strip())
            words = line.split()
            if i==0 and words[0] != 'DATASET':
                print("%s is not an ADH data file" % filename)
            if words[0] == 'BEGSCL':
                ftype = 'scalar'
                f = numpy.zeros((mesh.nNodes_global,),numpy.float_)
            elif words[0] == 'BEGVEC':
                ftype = 'vector'
                f = numpy.zeros((mesh.nNodes_global,3),numpy.float_)
            if words[0] == 'ND':
                nNodesInDataFile = int(words[1])
                assert(nNodesInDataFile == mesh.nNodes_global)
            if words[0] == 'NAME':
                description = ''
                for w in words[1:]: description += w
                description = description.strip('\"')
                nCopies=1
                ensightDescription=description
                while descriptionList.count(ensightDescription):
                    nCopies+=1
                    ensightDescription=description+"-%i" % nCopies
                descriptionList.append(ensightDescription)
        ftypeList.append(ftype)
        tList=[]
        while line.strip() != 'ENDDS':
            print("Reading "+line.strip())
            words = line.split()
            prematureEOF=False
            if ftype == 'scalar':
                for i in range(mesh.nNodes_global):
                    line = funIn.readline()
                    try:
                        f[i] = float(line)
                    except:
                        prematureEOF=True
                        break
                if prematureEOF:
                    break
                scalarOut=open(filename+'.scl%2.2i' % len(tList),'w')
                scalarOut.write(ensightDescription+'\n')
                scalarOut.write('part\n'+'%10i\n' % 1)
                scalarOut.write('coordinates\n')
                for fi in f:
                    scalarOut.write('%12.5E\n' % fi)
                scalarOut.close()
                if opts.boundaryMesh:
                    scalarOut=open(filename+'Boundary.scl%2.2i' % len(tList),'w')
                    scalarOut.write(ensightDescription+'\n')
                    scalarOut.write('part\n'+'%10i\n' % 1)
                    scalarOut.write('coordinates\n')
                    for nN in range(mesh.nExteriorNodes_global):
                        fi = f[mesh.exteriorNodeArray[nN]]
                        scalarOut.write('%12.5E\n' % fi)
                    scalarOut.close()
            elif ftype =='vector':
                for i in range(mesh.nNodes_global):
                    line = funIn.readline()
                    v = line.split()
                    try:
                        f[i,0] = float(v[0])
                        f[i,1] = float(v[1])
                        f[i,2] = float(v[2])
                    except:
                        prematureEOF=True
                if prematureEOF:
                    break
                vectorOut=open(filename+'.vec%2.2i' % len(tList),'w')
                vectorOut.write(ensightDescription+'\n')
                vectorOut.write('part\n'+'%10i\n' % 1)
                vectorOut.write('coordinates\n')
                for fi in f:
                    vectorOut.write('%12.5E\n' % fi[0])
                for fi in f:
                    vectorOut.write('%12.5E\n' % fi[1])
                for fi in f:
                    vectorOut.write('%12.5E\n' % fi[2])
                vectorOut.close()
                if opts.boundaryMesh:
                    vectorOut=open(filename+'Boundary.vec%2.2i' % len(tList),'w')
                    vectorOut.write(ensightDescription+'\n')
                    vectorOut.write('part\n'+'%10i\n' % 1)
                    vectorOut.write('coordinates\n')
                    for nN in range(mesh.nExteriorNodes_global):
                        fi = f[mesh.exteriorNodeArray[nN]]
                        vectorOut.write('%12.5E\n' % fi[0])
                    for nN in range(mesh.nExteriorNodes_global):
                        fi = f[mesh.exteriorNodeArray[nN]]
                        vectorOut.write('%12.5E\n' % fi[1])
                    for nN in range(mesh.nExteriorNodes_global):
                        fi = f[mesh.exteriorNodeArray[nN]]
                        vectorOut.write('%12.5E\n' % fi[2])
                    vectorOut.close()
            line = funIn.readline()
            tList.append(float(words[-1]))
        funIn.close()
        print("Read %i timesteps of %s" % (len(tList),ensightDescription))
        #check time steps
        for i in range(1,len(tList)-1):
            if tList[i-1] == tList[i] or tList[i] == tList[i+1]:
                print("the time series contains a time step below the output precision, adding a perturbation to preven an error  in ensight")
                tList[i] = 0.5*(tList[i-1] + tList[i+1])
    caseOut=open(ensightFilename+'.case','a')
    caseOut.write('VARIABLE\n')
    if opts.boundaryMesh:
        caseBoundaryOut=open(ensightFilename+'Boundary.case','a')
        caseBoundaryOut.write('VARIABLE\n')
    ts=1
    for filename,ftype,description in zip(funFilenameList,ftypeList,descriptionList):
        if ftype == 'scalar':
            caseOut.write('scalar per node: '+repr(ts)+' '+description+' '+filename+'.scl**\n')
            if opts.boundaryMesh:
                caseBoundaryOut.write('scalar per node: '+repr(ts)+
                                      ' '+description+' '+filename+'Boundary.scl**\n')
        elif ftype == 'vector':
            caseOut.write('vector per node: '+repr(ts)+' '+description+' '+filename+'.vec**\n')
            if opts.boundaryMesh:
                caseBoundaryOut.write('vector per node: '+repr(ts)+
                                      ' '+description+' '+filename+'Boundary.vec**\n')
    lines = ('TIME\n'+'time set: '+repr(ts)+' '+description+'\n'+
            'number of steps: '+ repr(len(tList))+'\n'+
            'filename start number: 0\n'+
            'filename increment: 1\n'+
            'time values:')
    caseOut.write(lines)
    if opts.boundaryMesh:
        caseBoundaryOut.write(lines)
    for tn in tList:
        caseOut.write(' %12.5E' % tn)
        caseOut.write('\n')
        if opts.boundaryMesh:
            caseBoundaryOut.write(' %12.5E' % tn)
            caseBoundaryOut.write('\n')
    caseOut.close()
    if opts.boundaryMesh:
        caseBoundaryOut.close()

if __name__ == '__main__':
##     import sys
##     import profile
##     profile.run('readFun(sys.argv[1],sys.argv[2],sys.argv[3])','readFunProf')
##     import pstats
##     p = pstats.Stats('readFunProf')
##     p.sort_stats('cumulative').print_stats(20)
##     p.sort_stats('time').print_stats(20)
    adh2ensight()

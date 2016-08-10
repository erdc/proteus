#!/usr/bin/env python
"""
A collection of functions for manipulating triangleWrappers interface to
Shewchuk's Triangle package. There are a couple of different groups of
functions.

One is for writing out triangleWrappers data arrays in Triangle's fileformat

.. inheritance-diagram:: proteus.TriangleUtils
   :parts: 1
"""
#import standard Python modules
import sys,os
import numpy
import triangleWrappers
import TriangleFileUtils


########################################################################
#some global constants and useful things?
########################################################################
#where is showme on the system?
PROTEUS_PACKAGES = os.getenv('PROTEUS_PACKAGES',os.getenv('HOME')+'/src/proteus-packages')
showmeCmdBase = PROTEUS_PACKAGES+'/triangle/bin/showme'

########################################################################
#this set of functions is for writing out triangulations
########################################################################
def writeOutTriangulation(tri,filebase="mesh",nbase=0,verbose=0):
    """
    collect necessary steps to write out files for triangle data structures

    """
    failed = False
    if tri == None:
        failed = True
        return failed
    if not nbase in [0,1]:
        print 'must have vertex base numbering = 0 or 1, got nbase= ',nbase
        failed = True
        return failed
    #end

    nodes = triangleWrappers.getPoints(tri)
    patts = triangleWrappers.getPointAttributes(tri)
    pmarks= triangleWrappers.getPointMarkers(tri)
    if verbose > 0:
        print 'writing node file'
    printNodeFile(nodes,filebase,patts,pmarks,nbase=nbase)

    elems = triangleWrappers.getTriangles(tri)
    eatts = triangleWrappers.getTriangleAttributes(tri)
    if verbose > 0:
        print 'writing elements file'
    printElemFile(elems,filebase,eatts,nbase=nbase)

    segms  = triangleWrappers.getSegments(tri)
    segmmks= triangleWrappers.getSegmentMarkers(tri)
    holes  = triangleWrappers.getHoles(tri)
    regions= triangleWrappers.getRegions(tri)
    if verbose > 0:
        print 'writing poly file'
    printPolyFile(nodes,segms,filebase,patts,pmarks,
                  segmmks,holes,regions,nbase=nbase)

    edges  = triangleWrappers.getEdges(tri)
    edgemks= triangleWrappers.getEdgeMarkers(tri)
    if verbose > 0:
        print 'writing edges file'
    printEdgeFile(edges,filebase,edgemks,nbase=nbase)

    neigs  = triangleWrappers.getNeighbors(tri)
    if verbose > 0:
        print 'writing neighbors file'
    printNeighborFile(neigs,filebase,nbase=nbase)

    return failed

#end writeNodesFile
#print out data structures in format triangle is expecting
def printNodeFile(nodes,filebase='mesh',attrib=None,markers=None,nbase=0):
    """
    print out nodes in triangle format, nbase is the numbering base

    file format that triangle is expecting

    First line:  <# of vertices> <dimension (must be 2)>  <# of attributes> <# of boundary markers (0 or 1)>

    Remaining lines: <vertex #> <x> <y>  [attributes] [boundary marker]
    """
    failed = False
    if nodes == None: #ok to do nothing of nodes is empy
        return failed
    #end empty areas check
    if not nbase in [0,1]:
        print 'must have vertex base numbering = 0 or 1, got nbase= ',nbase
        failed = True
        return failed
    #end
    filename = filebase+'.node'
    nout = open(filename,'w')

    format = """
#file format for nodes
#<# of vertices> <dim (must be 2)>  <# of attr> <# of boundary markers (0 or 1)>
#Remaining lines: <vertex #> <x> <y>  [attributes] [boundary marker]
"""
    nout.write(format)
    nvert   = nodes.shape[0]
    sdim    = 2
    nattrib = 0
    if attrib != None:
        nattrib = attrib.shape[1]
    #end if
    nmarkers= 0
    if markers != None:
        nmarkers = 1
    #end if
    line = """%d %d %d %d \n""" % (nvert,sdim,nattrib,nmarkers)
    nout.write(line)
    for i in range(nvert):
        line = """%d %g %g """ % (nbase+i,nodes[i,0],nodes[i,1])
        for j in range(nattrib):
            line +=""" %g """ % attrib[i,j]
        #end j
        if nmarkers > 0:
            line +=""" %d """ % markers[i]
        #end if
        line += '\n'
        nout.write(line)
    #end i
    nout.close()

    return failed
#end printNodeFile
#print out data structures in format triangle is expecting
def printElemFile(elems,filebase='mesh',attrib=None,nbase=0):
    """
    print out triangles/elements in the mesh

    file format that triangle is expecting
    First line: <# of triangles> <nodes per triangle>  <# of attributes>
    Remaining lines: <triangle #> <node> <node>  <node> ... [attributes]
    """
    failed = False
    if elems == None: #ok to do nothing of nodes is empy
        return failed
    #end empty areas check
    if not nbase in [0,1]:
        print 'must have vertex base numbering = 0 or 1, got nbase= ',nbase
        failed = True
        return failed
    #end
    filename = filebase+'.ele'
    fout = open(filename,'w')

    format = """
#file format for elements
#First line: <# of triangles> <nodes per triangle>  <# of attributes>
#Remaining lines: <triangle #> <node> <node>  <node> ... [attributes]

"""
    fout.write(format)
    nelem   = elems.shape[0]
    nodesPerTri = elems.shape[1] #should be 3 or 6
    nattrib = 0
    if attrib != None:
        nattrib = attrib.shape[1]
    #end if
    line = """%d %d %d \n""" % (nelem,nodesPerTri,nattrib)
    fout.write(line)
    for i in range(nelem):
        line = """%d %d %d %d """ % (nbase+i,elems[i,0],elems[i,1],elems[i,2])
        if nodesPerTri > 3:
            line += """%d %d %d """ % (elems[i,3],elems[i,4],elems[i,5])
        for j in range(nattrib):
            line +=""" %g """ % attrib[i,j]
        #end j
        line += '\n'
        fout.write(line)
    #end i
    fout.close()

    return failed
#end printElemFile
#print out data structures in format triangle is expecting
def printPolyFile(nodes,segments,filebase='mesh',
                  pattrib=None,pmarkers=None,
                  smarkers=None,
                  holes=None,regions=None,
                  nbase=0):
    """
    print out mesh info in triangle's poly format, nbase is the numbering base

    First line: <# of vertices> <dimension (must be 2)>  <# of attributes> <# of boundary markers (0 or 1)>
    Following lines: <vertex #> <x> <y> [attributes] [boundary marker]
    One line: <# of segments> <# of boundary markers (0 or 1)>
    Following lines: <segment #> <endpoint> <endpoint> [boundary marker]
    One line: <# of holes>
    Following lines: <hole #> <x> <y>
    Optional line: <# of regional attributes and/or area constraints>
    Optional following lines: <region #> <x> <y> <attribute> <maximum area>

    """
    failed = False
    if nodes == None: #ok to do nothing of nodes is empy
        return failed
    if segments == None:
        if nodes == None:
            return failed
        else:
            failed = True
            return failed
        #end if
    #end if
    #end empty areas check
    if not nbase in [0,1]:
        print 'must have vertex base numbering = 0 or 1, got nbase= ',nbase
        failed = True
        return failed
    #end
    filename = filebase+'.poly'
    fout = open(filename,'w')

    format = """
#First line: <# of vertices> <dim (2)>  <# of pt. attrs> <# of pt. bound markers (0 or 1)>
#Following lines: <vertex #> <x> <y> [attributes] [boundary marker]
#One line: <# of segments> <# of boundary markers (0 or 1)>
#Following lines: <segment #> <endpoint> <endpoint> [boundary marker]
#One line: <# of holes>
#Following lines: <hole #> <x> <y>
#Optional line: <# of regional attributes and/or area constraints>
#Optional following lines: <region #> <x> <y> <attribute> <maximum area>

"""
    fout.write(format)
    nvert   = nodes.shape[0]
    sdim    = 2
    npattrib= 0
    if pattrib != None:
        npattrib = pattrib.shape[1]
    #end if
    npmarkers= 0
    if pmarkers != None:
        npmarkers = 1
    #end if
    #write points out
    line = """%d %d %d %d \n""" % (nvert,sdim,npattrib,npmarkers)
    fout.write(line)
    for i in range(nvert):
        line = """%d %g %g """ % (nbase+i,nodes[i,0],nodes[i,1])
        for j in range(npattrib):
            line +=""" %g """ % pattrib[i,j]
        #end j
        if npmarkers > 0:
            line +=""" %d """ % pmarkers[i]
        #end if
        line += '\n'
        fout.write(line)
    #end i
    #write segments out
    nsegments = segments.shape[0]
    nsmarkers = 0
    if smarkers != None:
        nsmarkers = 1
    line = """%d %d \n""" % (nsegments,nsmarkers)
    fout.write(line)
    for i in range(nsegments):
        line = """%d %d %d""" % (nbase+i,segments[i,0],segments[i,1])
        if nsmarkers > 0:
            line += """ %d """ % smarkers[i]
        #end if
        line += '\n'
        fout.write(line)
    #end i segment loop

    nholes = 0
    if holes != None:
        nholes = holes.shape[0]
    #end if
    line ="""%d\n""" % nholes
    fout.write(line)
    for i in range(nholes):
        line = """%d %g %g """ % (nbase+i,holes[i,0],holes[i,1])
        line += '\n'
        fout.write(line)
    #end i for holes
    nregions = 0
    if regions != None:
        nregions = regions.shape[0]
    #end if
    line = "#no regions specified"

    if nregions > 0:
        line = """%d\n""" % nregions
    #end if
    fout.write(line)
    for i in range(nregions):
        line = """%d %g %g %g %g """ % (nbase+i,regions[i,0],regions[i,1],
                                        regions[i,2],regions[i,3])
        line += '\n'
        fout.write(line)
    #end i loop for regions
    fout.close()

    return failed
#end printPolyFile
#print out data structures in format triangle is expecting
def printAreaFile(elemAreas,filebase='mesh',nbase=0):
    """
    print out area constraints

    file format that triangle is expecting
    First line: <# of triangles>
    Remaining lines: <triangle #> <maximum area>

    note, negative area means that the area is unconstrained
    """
    failed = False
    if elemAreas == None: #ok to do nothing if elemAreas is empty
        return failed
    #end elemAreas empty
    if not nbase in [0,1]:
        print 'must have vertex base numbering = 0 or 1, got nbase= ',nbase
        failed = True
        return failed
    #end
    filename = filebase+'.area'
    fout = open(filename,'w')

    format = """
#file format for areas
#First line: <# of triangles>
#Remaining lines: <triangle #> <maximum area>

"""
    fout.write(format)
    nelem   = elemAreass.shape[0]
    line = """%d\n""" % nelem
    fout.write(line)
    for i in range(nelem):
        line = """%d %d %g """ % (nbase+i,elemAreas[i])
        line += '\n'
        fout.write(line)
    #end i
    fout.close()

    return failed
#end printAreaFile

def printEdgeFile(edges,filebase='mesh',markers=None,nbase=0):
    """
    print out edges in the mesh

    file format that triangle is expecting
    First line: <# of edges> <# of boundary markers (0 or 1)>
    Following lines: <edge #> <endpoint> <endpoint> [boundary marker]
    """
    failed = False
    if edges == None: #ok to do nothing if elemAreas is empty
        return failed
    #end empty arrays check

    if not nbase in [0,1]:
        print 'must have vertex base numbering = 0 or 1, got nbase= ',nbase
        failed = True
        return failed
    #end
    filename = filebase+'.edge'
    fout = open(filename,'w')

    format = """
#file format for edges
#First line: <# of edges> <# of boundary markers (0 or 1)>
#Following lines: <edge #> <endpoint> <endpoint> [boundary marker]
"""
    fout.write(format)
    nedges   = edges.shape[0]
    nmarkers = 0
    if markers != None:
        nmarkers = 1
    #end if
    line = """%d %d \n""" % (nedges,nmarkers)
    fout.write(line)
    for i in range(nedges):
        line = """%d %d %d """ % (nbase+i,edges[i,0],edges[i,1])
        if nmarkers > 0:
            line += """ %d """ % markers[i]
        #end j
        line += '\n'
        fout.write(line)
    #end i
    fout.close()

    return failed
#end printElemFile

def printNeighborFile(neigs,filebase='mesh',nbase=0):
    """
    print out elment neighbors  in the mesh

    file format that triangle is expecting
    First line: <# of triangles>  <# of neighbors per triangle (always 3)>
    Following lines: <triangle #> <neighbor> <neighbor> <neighbor>
    """
    failed = False
    if neigs == None: #ok to do nothing if elemAreas is empty
        return failed
    #end
    if not nbase in [0,1]:
        print 'must have vertex base numbering = 0 or 1, got nbase= ',nbase
        failed = True
        return failed
    #end
    filename = filebase+'.neig'
    fout = open(filename,'w')

    format = """
#file format for neighbors
#First line: <# of triangles>  <# of neighbors per triangle (always 3)>
#Following lines: <triangle #> <neighbor> <neighbor> <neighbor>
"""
    fout.write(format)
    nelems   = neigs.shape[0]
    nneigsPerElem = 3
    line = """%d %d \n""" % (nelems,nneigsPerElem)
    fout.write(line)
    for i in range(nelems):
        line = """%d %d %d %d""" % (nbase+i,neigs[i,0],neigs[i,1],neigs[i,2])
        line += '\n'
        fout.write(line)
    #end i
    fout.close()

    return failed
#end printNeigFile


########################################################################
#this set of functions is for reading a triangulation in
########################################################################

class TriangleInputFileReader:
    """
    collect format strings, specifiers and routines in a class so that
    I can keep track of them better
    """
    def __init__(self,verbose=0):
        """
        just specifies a collection of fixed format strings for the various
        input file formats allowed
        """
        #use some utilities for doing this stuff
        fUtils = TriangleFileUtils
        #allowed input file types
        self.inputFileTypes = ['node','ele','poly']
        self.fileExtensions = ['.'+type for type in self.inputFileTypes]

        #data types to be read in
        self.inputDataTypes = ['node','triangle','segment','hole','region']
        #keep track of what I'm supposed to be reading
        self.recordFormatDoc = {}
        for type in self.inputDataTypes:
            self.recordFormatDoc[type] = """format doc string not set"""
        #end for

        #now go through and specify the details of different formats
        self.recordInfo     = {}
        for type in self.inputDataTypes: #set defaults
            self.recordInfo[type] = {}
            #number of possible records per entry
            self.recordInfo[type]['nRecords'] = 0
            #underlying data type for each record
            self.recordInfo[type]['types']    = []
            #default values for each type to say it is uninitialized
            self.recordInfo[type]['defaultValues'] = []
            #how to convert a data value of this type to array entry
            self.recordInfo[type]['conversions']   = []
        #end type initialization

        ###go ahead and set the necessary information for each type
        #nodes
        self.recordFormatDoc['node'] = """
#<vertex #> <x> <y>  [attributes] [boundary marker]
        """
        self.recordInfo['node']['nRecords'] = 3
        self.recordInfo['node']['types'] = ['d','d','i']
        self.recordInfo['node']['defaultValues'] = [-12345.0,-12345.0,-12345]
        self.recordInfo['node']['conversions'].append(lambda x : float(x))
        self.recordInfo['node']['conversions'].append(lambda x : float(x))
        self.recordInfo['node']['conversions'].append(lambda x : int(x))

        #elems (triangles)
        self.recordFormatDoc['triangle'] = """
#<triangle #> <node> <node>  <node> ... [attributes]
        """
        self.recordInfo['triangle']['nRecords'] = 2
        self.recordInfo['triangle']['types']    = ['i','d']
        self.recordInfo['triangle']['defaultValues'] = [-12345,-12345.0]
        self.recordInfo['triangle']['conversions'].append(lambda x : int(x))
        self.recordInfo['triangle']['conversions'].append(lambda x : float(x))

        #segments
        self.recordFormatDoc['segment'] = """
#<segment #> <endpoint> <endpoint> [boundary marker]
        """
        self.recordInfo['segment']['nRecords'] = 2
        self.recordInfo['segment']['types']    = ['i','i']
        self.recordInfo['segment']['defaultValues'] = [-12345,-12345]
        self.recordInfo['segment']['conversions'].append(lambda x : int(x))
        self.recordInfo['segment']['conversions'].append(lambda x : int(x))

        #holes
        self.recordFormatDoc['hole'] = """
#<hole #> <x> <y>
        """
        self.recordInfo['hole']['nRecords'] = 1
        self.recordInfo['hole']['types']    = ['d']
        self.recordInfo['hole']['defaultValues'] = [-12345.0]
        self.recordInfo['hole']['conversions'].append(lambda x : float(x))

        #regions
        self.recordFormatDoc['region'] = """
#<region #> <x> <y> <attribute> <maximum area>
        """
        self.recordInfo['region']['nRecords'] = 1
        self.recordInfo['region']['types']    = ['d']
        self.recordInfo['region']['defaultValues'] = [-12345.0]
        self.recordInfo['region']['conversions'].append(lambda x : float(x))

        ## end record info block


        ##now set the line format for each data type that one might read in
        self.initialLineFormat = {}
        for type in self.inputDataTypes:
            dinfo = self.recordInfo[type]
            self.initialLineFormat[type]=fUtils.generateDefaultFormatPattern(dinfo,
                                                                             verbose)
        #end for
        #have to set the line formats for segments,holes, and regions manually
        #actual format is <# of segments> <# of boundary markers (0 or 1)>
        self.initialLineFormat['segment'] = r'^\s*(\d+)\s+(\d+)\s*$'
        #hole section just contains the number of holes
        self.initialLineFormat['hole'] = r'^\s*(\d+)\s*$'
        #region section just contains the number of holes
        self.initialLineFormat['region'] = r'^\s*(\d+)\s*$'

        #because of format discrepancies, keep a different record info object
        #for processing initial line in file
        self.recordInfoInit = {}
        for type in self.inputDataTypes:
            #get a deep copy
            self.recordInfoInit[type] = self.recordInfo[type].copy()
        #end
        #fix segment, hole, region to match expected format
        self.recordInfoInit['segment']['nRecords'] = 1
        self.recordInfoInit['segment']['types'] = ['i']
        self.recordInfoInit['hole']['nRecords'] = 0
        self.recordInfoInit['hole']['types']    = []
        self.recordInfoInit['region']['nRecords'] = 0
        self.recordInfoInit['region']['types']    = []

    #end init

    def readNodes(self,filebase,commchar='#',nbase=0,verbose=0):
        """
        read nodes specified in Triangle format
        see recordFormatDoc for format info

        returns dataInfo holding basic info about data read and array of data read

        data[0] = nodes
        data[1] = node attributes
        data[1] = node markers
        """
        type = 'node'
        filename = filebase+'.node'
        fout = open(filename,'r')
        lines = fout.readlines()
        fout.close()

        dataInfo = TriangleFileUtils.findInitialFormatLine(lines,
                                                           self.recordInfoInit[type],
                                                           self.initialLineFormat[type],
                                                           commchar=commchar,
                                                           nbase=nbase,verbose=verbose)

        ndstart   = dataInfo['dataStartLineNumber']
        data,ilast= TriangleFileUtils.readSimpleFormattedDataEntries(lines[ndstart:],
                                                                     self.recordInfo[type],
                                                                     dataInfo,
                                                                     commchar=commchar,
                                                                     nbase=nbase,
                                                                     verbose=verbose)

        output = {}
        output['nodes']          = data[0]
        output['nodeAttributes'] = data[1]
        output['nodeMarkers']    = data[2]
        return dataInfo,output
    #end readNodes
    def readTriangles(self,filebase,commchar='#',nbase=0,verbose=0):
        """
        read triangles specified in Triangle format
        see recordFormatDoc for format info

        returns dataInfo holding basic info about data read and array of data read

        data[0] = triangles
        data[1] = triangle attributes

        """
        type = 'triangle'
        filename = filebase+'.ele'
        fout = open(filename,'r')
        lines = fout.readlines()
        fout.close()

        dataInfo = TriangleFileUtils.findInitialFormatLine(lines,
                                                           self.recordInfoInit[type],
                                                           self.initialLineFormat[type],
                                                           commchar=commchar,
                                                           nbase=nbase,verbose=verbose)

        ndstart   = dataInfo['dataStartLineNumber']
        data,ilast= TriangleFileUtils.readSimpleFormattedDataEntries(lines[ndstart:],
                                                                     self.recordInfo[type],
                                                                     dataInfo,
                                                                     commchar=commchar,
                                                                     nbase=nbase,
                                                                     verbose=verbose)

        output = {}
        output['triangles']          = data[0]
        output['triangleAttributes'] = data[1]

        return dataInfo,output

    #end readTriangles
    def readPoly(self,filebase,commchar='#',nbase=0,verbose=0):
        """
        read nodes, segments, holes, and regions specified in Triangle format
        see recordFormatDoc for format info

        returns dataInfo holding basic info about data read and array of data read

        this function is more tedious because of variations in input formats

        """
        filename = filebase+'.poly'
        fout = open(filename,'r')
        lines = fout.readlines()
        fout.close()

        ##start with nodes
        type = 'node'
        nodeInfo = TriangleFileUtils.findInitialFormatLine(lines,
                                                           self.recordInfoInit[type],
                                                           self.initialLineFormat[type],
                                                           commchar=commchar,
                                                           nbase=nbase,verbose=verbose)

        nodestart   = nodeInfo['dataStartLineNumber']
        #mwf debug print "recordInfo ",self.recordInfo[type]
        nodeData,nlast=TriangleFileUtils.readSimpleFormattedDataEntries(lines[nodestart:],
                                                                        self.recordInfo[type],
                                                                        nodeInfo,
                                                                        commchar=commchar,
                                                                        nbase=nbase,
                                                                        verbose=verbose)
        ##mwf debug print "node type after read simple",nodeData[2].dtype
        ##now segments
        type = 'segment'
        #start reading after node section
        slines   = lines[nlast:]
        segInfo0 = TriangleFileUtils.findInitialFormatLine(slines,
                                                           self.recordInfoInit[type],
                                                           self.initialLineFormat[type],
                                                           commchar=commchar,
                                                           nbase=nbase,verbose=verbose)

        #now add segment size into the expected record sizes per entry
        segInfo = segInfo0.copy()
        #segment size is 2 (beginning node and ending node)
        #segInfo0['recordSizes'] should be 0 or 1
        assert(segInfo0['recordSizes'][0] == 0 or
               segInfo0['recordSizes'][0] == 1)
        segSize   = 2
        markerSize= int(segInfo0['recordSizes'][0])
        segInfo['recordSizes'] = [segSize,markerSize]
        #first value per line is the entry number
        #record I is in entries recLocPerLine[I]:recLocPerLine[I+1]
        segInfo['recordLocationsPerLine'] = numpy.array([0,1,1+segSize,
                                                           1+segSize+markerSize],
                                                          'i')
        segdstart   = segInfo['dataStartLineNumber']
        segData,slast=TriangleFileUtils.readSimpleFormattedDataEntries(slines[segdstart:],
                                                                       self.recordInfo[type],
                                                                       segInfo,
                                                                       commchar=commchar,
                                                                       nbase=nbase,
                                                                       verbose=verbose)

        ##now holes
        type = 'hole'
        #start reading where segment section ended (cound from beginning of segment
        hlines= slines[slast+segdstart:]

        holeInfo0 = TriangleFileUtils.findInitialFormatLine(hlines,
                                                            self.recordInfoInit[type],
                                                            self.initialLineFormat[type],
                                                            commchar=commchar,
                                                            nbase=nbase,verbose=verbose)

        #now add hole size into the expected record sizes per entry
        holeInfo = holeInfo0.copy()
        #hole size is 2 (x and y)
        sdim = 2
        holeInfo['recordSizes'] = [sdim] #no optional sizing, just may
                                         #or may not be data included
        holeInfo['recordLocationsPerLine'] = numpy.array([0,1,1+sdim],'i')

        holedstart   = holeInfo['dataStartLineNumber']
        holeData,hlast=TriangleFileUtils.readSimpleFormattedDataEntries(hlines[holedstart:],
                                                                        self.recordInfo[type],
                                                                        holeInfo,
                                                                        commchar=commchar,
                                                                        nbase=nbase,
                                                                        verbose=verbose)

        ##last read regional constraints
        type = 'region'
        #start where holes stopped
        rlines = hlines[hlast+holedstart:]
        regionInfo0 = TriangleFileUtils.findInitialFormatLine(rlines,
                                                           self.recordInfoInit[type],
                                                           self.initialLineFormat[type],
                                                           commchar=commchar,
                                                           nbase=nbase,verbose=verbose)


        #now add hole size into the expected record sizes per entry
        regionInfo = regionInfo0.copy()
        #region size is 4 (x,y, attribute, area constraint)
        #actually, Triangle allows user to omit one of the last two values!
        rsize = 4
        regionInfo['recordSizes'] = [rsize] #no optional sizing, just may
                                         #or may not be data included
        regionInfo['recordLocationsPerLine'] = numpy.array([0,1,1+rsize],'i')

        regdstart   = regionInfo['dataStartLineNumber']
#         regionData,rlast=TriangleFileUtils.readSimpleFormattedDataEntries(rlines[regdstart:],
#                                                                           self.recordInfo[type],
#                                                                           regionInfo,
#                                                                           commchar=commchar,
#                                                                           nbase=nbase,
#                                                                           verbose=verbose)
        regionData,rlast=\
           TriangleFileUtils.readSimpleFormattedDataEntriesLastOptional(rlines[regdstart:],
                                                                        self.recordInfo[type],
                                                                        regionInfo,
                                                                        commchar=commchar,
                                                                        nbase=nbase,
                                                                        verbose=verbose)

        outputInfo = {}
        outputData = {}
        #info
        outputInfo['node']   = nodeInfo
        outputInfo['segment']= segInfo
        outputInfo['hole']   = holeInfo
        outputInfo['region'] = regionInfo
        #data
        outputData['node']   = {}
        if not len(nodeData) == 3:
            for val in ['nodes','nodeAttributes','nodeMarkers']:
                outputData['node'][val] = None
        else: #end wrong size
            outputData['node']['nodes']          = nodeData[0]
            outputData['node']['nodeAttributes'] = nodeData[1]
            outputData['node']['nodeMarkers']    = nodeData[2]
        #end if
        outputData['segment']= {}
        if not len(segData) == 2:
            for val in ['segments','segmentMarkers']:
                outputData['segment'][val] = None
        else: #end wrong size
            outputData['segment']['segments']      = segData[0]
            outputData['segment']['segmentMarkers']= segData[1]
        #end
        outputData['hole']   = {}
        outputData['hole']['holes'] = None
        if not holeData == None:
            outputData['hole']['holes']  = holeData[0]

        outputData['region'] = {}
        outputData['region']['regions'] = None
        if not regionData == None:
            outputData['region']['regions'] = regionData[0]

        return outputInfo,outputData
    #end readPoly
#end input file reader

"""
Classes for archiving numerical solution data
"""
import Profiling
import Comm
import numpy
import os
from xml.etree.ElementTree import *

log = Profiling.logEvent

def indentXML(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        for e in elem:
            indentXML(e, level+1)
            if not e.tail or not e.tail.strip():
                e.tail = i + "  "
        if not e.tail or not e.tail.strip():
            e.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i

class ArchiveFlags:
    EVERY_MODEL_STEP     = 0
    EVERY_USER_STEP      = 1
    EVERY_SEQUENCE_STEP  = 2
    UNDEFINED            =-1
#
class AR_base:
    def __init__(self,dataDir,filename,useTextArchive=False,gatherAtClose=True,hotStart=False,readOnly=False):
        import os.path
        comm=Comm.get()
        self.comm=comm
        self.dataDir=dataDir
        self.filename=filename
        if hotStart:
            self.filename+="hot"
        self.comm=comm
        self.rank = comm.rank()
        self.size = comm.size()
        import datetime
        #filename += datetime.datetime.now().isoformat()
        try:
            import tables
            self.hasTables=True
        except:
            self.hasTables=False
        self.xmlHeader = "<?xml version=\"1.0\" ?>\n<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n"
        if hotStart:
            self.ar_old =  AR_base(dataDir,filename,useTextArchive,gatherAtClose=False,hotStart=False,readOnly=True)
            self.xmlFile=open(os.path.join(self.dataDir,filename+"hot"+str(self.rank)+".xmf"),"w")
            self.tree=ElementTree(Element("Xdmf",
                                          {"Version":"2.0",
                                           "xmlns:xi":"http://www.w3.org/2001/XInclude"}))
            if self.hasTables and not useTextArchive:
                self.hdfFilename=filename+"hot"+str(self.rank)+".h5"
                self.hdfFile=tables.openFile(os.path.join(self.dataDir,self.hdfFilename),
                                             mode = "w",
                                             title = filename+" Data")
                self.dataItemFormat="HDF"
            else:
                self.textDataDir=filename+"hot_Data"
                if not os.path.exists(self.textDataDir):
                    try:
                        os.mkdir(self.textDataDir)
                    except:
                        self.textDataDir=""
                self.hdfFile=None
                self.dataItemFormat="XML"
        elif readOnly:
            self.xmlFile=open(os.path.join(self.dataDir,filename+str(self.rank)+".xmf"),"r")
            self.tree=ElementTree(file=self.xmlFile)
            if self.hasTables and not useTextArchive:
                self.hdfFilename=filename+str(self.rank)+".h5"
                self.hdfFile=tables.openFile(os.path.join(self.dataDir,self.hdfFilename),
                                             mode = "r",
                                             title = filename+" Data")
                self.dataItemFormat="HDF"
            else:
                self.textDataDir=filename+"_Data"
                assert(os.path.exists(self.textDataDir))
                self.hdfFile=None
                self.dataItemFormat="XML"            
        else:
            self.xmlFile=open(os.path.join(self.dataDir,filename+str(self.rank)+".xmf"),"w")
            self.tree=ElementTree(Element("Xdmf",
                                          {"Version":"2.0",
                                           "xmlns:xi":"http://www.w3.org/2001/XInclude"}))
            if self.hasTables and not useTextArchive:
                self.hdfFilename=filename+str(self.rank)+".h5"
                self.hdfFile=tables.openFile(os.path.join(self.dataDir,self.hdfFilename),
                                             mode = "w",
                                             title = filename+" Data")
                self.dataItemFormat="HDF"
            else:
                self.textDataDir=filename+"_Data"
                if not os.path.exists(self.textDataDir):
                    try:
                        os.mkdir(self.textDataDir)
                    except:
                        self.textDataDir=""
                self.hdfFile=None
                self.dataItemFormat="XML"
        #
        self.gatherAtClose = gatherAtClose
    def clear_xml(self):
        self.xmlFile.seek(0)
        self.xmlFile.truncate()
    def close(self):
        log("Closing Archive")
        self.clear_xml()
        self.xmlFile.write(self.xmlHeader)
        indentXML(self.tree.getroot())
        self.tree.write(self.xmlFile)
        self.xmlFile.close()
        if self.hdfFile != None:
            self.hdfFile.close()
        try:
            if self.gatherAtClose:
                self.allGather()
        except:
            pass
    def allGather(self):
        if self.rank==0:
            #replace the bottom level grid with a spatial collection
            XDMF_all=self.tree.getroot()
            Domain_all=XDMF_all[-1]
            for TemporalGridCollection in Domain_all:
                Grids = TemporalGridCollection[:]
                del TemporalGridCollection[:]
                for Grid in Grids:
                    SpatialCollection=SubElement(TemporalGridCollection,"Grid",{"GridType":"Collection",
                                                                                "CollectionType":"Spatial"})
                    SpatialCollection.append(Grid[0])#append Time in Spatial Collection
                    del Grid[0]#delete Time in grid
                    SpatialCollection.append(Grid) #append Grid without Time
            for i in range(1,self.size):
                xmlFile=open(os.path.join(self.dataDir,self.filename+str(i)+".xmf"),"r")
                tree = ElementTree(file=xmlFile)
                XDMF=tree.getroot()
                Domain=XDMF[-1]
                for TemporalGridCollection,TemporalGridCollection_all in zip(Domain,Domain_all):
                    SpatialGridCollections = TemporalGridCollection_all.findall('Grid')
                    for Grid,Grid_all in zip(TemporalGridCollection,SpatialGridCollections):
                        del Grid[0]#Time
                        Grid_all.append(Grid)
                xmlFile.close()
            f = open(os.path.join(self.dataDir,self.filename+"all_"+str(self.size)+".xmf"),"w")
            indentXML(self.tree.getroot())
            self.tree.write(f)
            f.close()
    def sync(self):
        log("Syncing Archive",level=3)
        self.clear_xml()
        self.xmlFile.write(self.xmlHeader)
        indentXML(self.tree.getroot())
        self.tree.write(self.xmlFile)
        self.xmlFile.flush()
        if self.hdfFile != None:
            self.hdfFile.flush()

XdmfArchive=AR_base

########################################################################
#for writing out various quantities in Xdmf format
#import xml.etree.ElementTree as ElementTree
#from xml.etree.ElementTree import SubElement
import numpy
class XdmfWriter:
    """
    collect functionality for writing data to Xdmf format

    Writer is supposed to keep track of grid collection (temporal collection)

    as well as current grid under grid collection where data belong,
    since data are associated with a grid of specific type
    (e.g., P1 Lagrange, P2 Lagrange, elementQuadrature dictionary, ...) 
    """
    def __init__(self,shareSingleGrid=True,arGridCollection=None,arGrid=None,arTime=None):
        self.arGridCollection = arGridCollection #collection of "grids" (at least one per time level)
        self.arGrid           = arGrid #grid in collection that data should be associated with 
        self.arTime           = arTime #time level for data
        self.shareSingleGrid  = shareSingleGrid

    def setGridCollectionAndGridElements(self,init,ar,arGrid,t,spaceSuffix):
        """
        attempt at boiler plate code to grab current arGridCollection and grid for
        a given type and time level t
        returns gridName to use in writing mesh
        """
        if init:
            #ar should have domain now and mesh should have gridCollection
            #but want own grid collection to start
            self.arGridCollection = SubElement(ar.domain,"Grid",{"Name":"Mesh"+spaceSuffix,
                                                                 "GridType":"Collection",
                                                                 "CollectionType":"Temporal"})
        elif self.arGridCollection == None:#try to get existing grid collection
            for child in ar.domain:
                if child.tag == "Grid" and child.attrib["Name"] == "Mesh"+spaceSuffix:
                    self.arGridCollection = child
        assert self.arGridCollection != None
        
        #see if dgp1 grid exists with current time?
        
        if self.shareSingleGrid:
            gridName = "Grid"+spaceSuffix
            if arGrid != None:
                self.arGrid = arGrid
                gt = arGrid.find("Time")
                self.arTime = gt
            #brute force search through child grids of arGridCollection
            #grids = self.arGridCollection.findall("Grid")
            #for g in grids:
            #    for gt in g:
            #        if gt.tag == "Time" and gt.attrib["Value"] == str(t):
            #            self.arTime = gt
            #            self.arGrid = g
            #            break
            #end brute force search
            
        else:
            gridName = "Grid"+spaceSuffix+name
        return gridName

    def writeMeshXdmf_elementQuadrature(self,ar,mesh,spaceDim,x,t=0.0,
                                       init=False,meshChanged=False,arGrid=None,tCount=0):
        return self.writeMeshXdmf_quadraturePoints(ar,mesh,spaceDim,x,t=t,quadratureType="q",
                                                   init=init,meshChanged=meshChanged,arGrid=arGrid,tCount=tCount)
    def writeMeshXdmf_elementBoundaryQuadrature(self,ar,mesh,spaceDim,x,t=0.0,
                                                init=False,meshChanged=False,arGrid=None,tCount=0):
        return self.writeMeshXdmf_quadraturePoints(ar,mesh,spaceDim,x,t=t,quadratureType="ebq_global",
                                                   init=init,meshChanged=meshChanged,arGrid=arGrid,tCount=tCount)
    def writeMeshXdmf_exteriorElementBoundaryQuadrature(self,ar,mesh,spaceDim,x,t=0.0,
                                                        init=False,meshChanged=False,arGrid=None,tCount=0):
        return self.writeMeshXdmf_quadraturePoints(ar,mesh,spaceDim,x,t=t,quadratureType="ebqe",
                                                   init=init,meshChanged=meshChanged,arGrid=arGrid,tCount=tCount)
    
    def writeScalarXdmf_elementQuadrature(self,ar,u,name,tCount=0,init=True):
        qualifiedName = name.replace(' ','_')+"_elementQuadrature"
        return self.writeScalarXdmf_quadrature(ar,u,qualifiedName,tCount=tCount,init=init)
    def writeVectorXdmf_elementQuadrature(self,ar,u,name,tCount=0,init=True):
        qualifiedName = name.replace(' ','_')+"_elementQuadrature"
        return self.writeVectorXdmf_quadrature(ar,u,qualifiedName,tCount=tCount,init=init)
    def writeTensorXdmf_elementQuadrature(self,ar,u,name,tCount=0,init=True):
        qualifiedName = name.replace(' ','_')+"_elementQuadrature"
        return self.writeTensorXdmf_quadrature(ar,u,qualifiedName,tCount=tCount,init=init)
    #
    def writeScalarXdmf_elementBoundaryQuadrature(self,ar,u,name,tCount=0,init=True):
        qualifiedName = name.replace(' ','_')+"_elementBoundaryQuadrature"
        return self.writeScalarXdmf_quadrature(ar,u,qualifiedName,tCount=tCount,init=init)
    def writeVectorXdmf_elementBoundaryQuadrature(self,ar,u,name,tCount=0,init=True):
        qualifiedName = name.replace(' ','_')+"_elementBoundaryQuadrature"
        return self.writeVectorXdmf_quadrature(ar,u,qualifiedName,tCount=tCount,init=init)
    def writeTensorXdmf_elementBoundaryQuadrature(self,ar,u,name,tCount=0,init=True):
        qualifiedName = name.replace(' ','_')+"_elementBoundaryQuadrature"
        return self.writeTensorXdmf_quadrature(ar,u,qualifiedName,tCount=tCount,init=init)
    #
    def writeScalarXdmf_exteriorElementBoundaryQuadrature(self,ar,u,name,tCount=0,init=True):
        qualifiedName = name.replace(' ','_')+"_exteriorElementBoundaryQuadrature"
        return self.writeScalarXdmf_quadrature(ar,u,qualifiedName,tCount=tCount,init=init)
    def writeVectorXdmf_exteriorElementBoundaryQuadrature(self,ar,u,name,tCount=0,init=True):
        qualifiedName = name.replace(' ','_')+"_exteriorElementBoundaryQuadrature"
        return self.writeVectorXdmf_quadrature(ar,u,qualifiedName,tCount=tCount,init=init)
    def writeTensorXdmf_exteriorElementBoundaryQuadrature(self,ar,u,name,tCount=0,init=True):
        qualifiedName = name.replace(' ','_')+"_exteriorElementBoundaryQuadrature"
        return self.writeTensorXdmf_quadrature(ar,u,qualifiedName,tCount=tCount,init=init)
    

    def writeMeshXdmf_quadraturePoints(self,ar,mesh,spaceDim,x,t=0.0,quadratureType="q",
                                       init=False,meshChanged=False,arGrid=None,tCount=0):
        """
        write out quadrature point mesh for quadrature points that are unique to
        either elements or elementBoundaries

        quadrature type
           q  --> element
           ebq_global --> element boundary
           ebqe       --> exterior element boundary
        """
        if quadratureType == "q":
            spaceSuffix = "_elementQuadrature"
        elif quadratureType == "ebq_global":
            spaceSuffix = "_elementBoundaryQuadrature"
        elif quadratureType == "ebqe":
            spaceSuffix = "_exteriorElementBoundaryQuadrature"
        else:
            raise RuntimeError("quadratureType = %s not recognized" % quadratureType)
        
        #write out basic geometry if not already done?
        mesh.writeMeshXdmf(ar,"Spatial_Domain",t,init,meshChanged,tCount=tCount)
        #now try to write out a mesh that is a collection of points per element
        
        gridName = self.setGridCollectionAndGridElements(init,ar,arGrid,t,spaceSuffix)

        assert len(x.shape) == 3 #make sure have right type of dictionary
        if self.arGrid == None or self.arTime.get('Value') != str(t):
            Xdmf_ElementTopology = "Polyvertex"
            Xdmf_NumberOfElements= x.shape[0]
            Xdmf_NodesPerElement = x.shape[1]
            Xdmf_NodesGlobal     = Xdmf_NumberOfElements*Xdmf_NodesPerElement
            
            self.arGrid = SubElement(self.arGridCollection,"Grid",{"Name":gridName,"GridType":"Uniform"})
            self.arTime = SubElement(self.arGrid,"Time",{"Value":str(t)})
            topology    = SubElement(self.arGrid,"Topology",
                                     {"Type":Xdmf_ElementTopology,
                                      "NumberOfElements":str(Xdmf_NumberOfElements),
                                      "NodesPerElement":str(Xdmf_NodesPerElement)})
            elements    = SubElement(topology,"DataItem",
                                     {"Format":ar.dataItemFormat,
                                      "DataType":"Int",
                                      "Dimensions":"%i %i" % (Xdmf_NumberOfElements,Xdmf_NodesPerElement)})
            geometry    = SubElement(self.arGrid,"Geometry",{"Type":"XYZ"})
            nodes       = SubElement(geometry,"DataItem",
                                     {"Format":ar.dataItemFormat,
                                      "DataType":"Float",
                                      "Precision":"8",
                                      "Dimensions":"%i %i" % (Xdmf_NodesGlobal,3)})
            if ar.hdfFile != None:
                elements.text = ar.hdfFilename+":/elements"+spaceSuffix+`tCount`
                nodes.text    = ar.hdfFilename+":/nodes"+spaceSuffix+`tCount`
                if init or meshChanged:
                    #this will fail if elements_dgp1 already exists
                    #q_l2g = numpy.zeros((Xdmf_NumberOfElements,Xdmf_NodesPerElement),'i')
                    #brute force to start
                    #for eN in range(Xdmf_NumberOfElements):
                    #    for nN in range(Xdmf_NodesPerElement):
                    #        q_l2g[eN,nN] = eN*Xdmf_NodesPerElement + nN
                    #
                    q_l2g = numpy.arange(Xdmf_NumberOfElements*Xdmf_NodesPerElement,dtype='i').reshape((Xdmf_NumberOfElements,Xdmf_NodesPerElement))
                    ar.hdfFile.createArray("/",'elements'+spaceSuffix+`tCount`,q_l2g)
                    ar.hdfFile.createArray("/",'nodes'+spaceSuffix+`tCount`,x.flat[:])
            else:
                SubElement(elements,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/elements"+spaceSuffix+`tCount`+".txt"})
                SubElement(nodes,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/nodes"+spaceSuffix+`tCount`+".txt"})
                if init or meshChanged:
                    #this will fail if elements_dgp1 already exists
                    #q_l2g = numpy.zeros((Xdmf_NumberOfElements,Xdmf_NodesPerElement),'i')
                    #brute force to start
                    #for eN in range(Xdmf_NumberOfElements):
                    #    for nN in range(Xdmf_NodesPerElement):
                    #        q_l2g[eN,nN] = eN*Xdmf_NodesPerElement + nN
                    #
                    q_l2g = numpy.arange(Xdmf_NumberOfElements*Xdmf_NodesPerElement,dtype='i').reshape((Xdmf_NumberOfElements,Xdmf_NodesPerElement))
                    numpy.savetxt(ar.textDataDir+"/elements"+spaceSuffix+`tCount`+".txt",q_l2g,fmt='%d')
                    numpy.savetxt(ar.textDataDir+"/nodes"+spaceSuffix+`tCount`+".txt",x.flat[:])
                    
                #
            #hdfile
        #need to write a grid
        return self.arGrid
    #def
    def writeScalarXdmf_quadrature(self,ar,u,name,tCount=0,init=True):
        assert len(u.shape) == 2
        Xdmf_NodesGlobal = u.shape[0]*u.shape[1]
        attribute = SubElement(self.arGrid,"Attribute",{"Name":name,
                                                        "AttributeType":"Scalar",
                                                        "Center":"Node"})
        values    = SubElement(attribute,"DataItem",
                               {"Format":ar.dataItemFormat,
                                "DataType":"Float",
                                "Precision":"8",
                                "Dimensions":"%i" % (Xdmf_NodesGlobal,)})
        if ar.hdfFile != None:
            values.text = ar.hdfFilename+":/"+name+str(tCount)
            ar.hdfFile.createArray("/",name+str(tCount),u.flat[:])
        else:
            numpy.savetxt(ar.textDataDir+"/"+name+str(tCount)+".txt",u.flat[:])
            SubElement(values,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/"+name+str(tCount)+".txt"}) 
        
    
    def writeVectorXdmf_quadrature(self,ar,u,name,tCount=0,init=True):
        assert len(u.shape) == 3
        Xdmf_NodesGlobal = u.shape[0]*u.shape[1]
        Xdmf_NumberOfComponents = u.shape[2]
        Xdmf_StorageDim = 3
        attribute = SubElement(self.arGrid,"Attribute",{"Name":name,
                                                        "AttributeType":"Vector",
                                                        "Center":"Node"})
        values    = SubElement(attribute,"DataItem",
                               {"Format":ar.dataItemFormat,
                                "DataType":"Float",
                                "Precision":"8",
                                "Dimensions":"%i %i" % (Xdmf_NodesGlobal,Xdmf_StorageDim)})#force 3d vector since points 3d
        tmp = numpy.zeros((Xdmf_NodesGlobal,Xdmf_StorageDim),'d')
        tmp[:,:Xdmf_NumberOfComponents]=numpy.reshape(u.flat,(Xdmf_NodesGlobal,Xdmf_NumberOfComponents))
        
        if ar.hdfFile != None:
            values.text = ar.hdfFilename+":/"+name+str(tCount)
            ar.hdfFile.createArray("/",name+str(tCount),tmp)
        else:
            numpy.savetxt(ar.textDataDir+"/"+name+str(tCount)+".txt",tmp)
            SubElement(values,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/"+name+str(tCount)+".txt"}) 
        
    
    def writeTensorXdmf_quadrature(self,ar,u,name,tCount=0,init=True):
        """
        TODO make faster tmp creation
        """
        assert len(u.shape) == 4
        Xdmf_NodesGlobal = u.shape[0]*u.shape[1]
        Xdmf_NumberOfComponents = u.shape[2]*u.shape[3] #Xdmf requires 9 though
        attribute = SubElement(self.arGrid,"Attribute",{"Name":name,
                                                        "AttributeType":"Tensor",
                                                        "Center":"Node"})
        values    = SubElement(attribute,"DataItem",
                               {"Format":ar.dataItemFormat,
                                "DataType":"Float",
                                "Precision":"8",
                                "Dimensions":"%i %i" % (Xdmf_NodesGlobal,9)})#force 3d vector since points 3d
        #mwf brute force
        tmp = numpy.zeros((Xdmf_NodesGlobal,9),'d')
        for k in range(Xdmf_NodesGlobal):
            for i in range(u.shape[2]):
                for j in range(u.shape[3]):
                    tmp.flat[k*9 + i*3 + j]=u.flat[k*Xdmf_NumberOfComponents + i*u.shape[2] + j]
        
        if ar.hdfFile != None:
            values.text = ar.hdfFilename+":/"+name+str(tCount)
            ar.hdfFile.createArray("/",name+str(tCount),tmp)
        else:
            numpy.savetxt(ar.textDataDir+"/"+name+str(tCount)+".txt",tmp)
            SubElement(values,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/"+name+str(tCount)+".txt"}) 
        
    
    def writeMeshXdmf_DGP1Lagrange(self,ar,name,mesh,spaceDim,dofMap,CGDOFMap,t=0.0,
                                   init=False,meshChanged=False,arGrid=None,tCount=0):
        """
        TODO Not tested yet
        """
        #assert False, "Not tested"
        #write out basic geometry if not already done?
        mesh.writeMeshXdmf(ar,"Spatial_Domain",t,init,meshChanged,tCount=tCount)
        #now try to write out a mesh that matches DG P1 layout

        spaceSuffix = "_dgp1_Lagrange"
        gridName = self.setGridCollectionAndGridElements(init,ar,arGrid,t,spaceSuffix)

        if self.arGrid == None or self.arTime.get('Value') != str(t):
            if spaceDim == 1:
                Xdmf_ElementTopology = "Polyline"
            elif spaceDim == 2:
                Xdmf_ElementTopology = "Triangle"
            elif spaceDim == 3:
                Xdmf_ElementTopology = "Tetrahedron"
            self.arGrid = SubElement(self.arGridCollection,"Grid",{"Name":gridName,"GridType":"Uniform"})
            self.arTime = SubElement(self.arGrid,"Time",{"Value":str(t)})
            topology    = SubElement(self.arGrid,"Topology",
                                     {"Type":Xdmf_ElementTopology,
                                      "NumberOfElements":str(mesh.nElements_global)})
            elements    = SubElement(topology,"DataItem",
                                     {"Format":ar.dataItemFormat,
                                      "DataType":"Int",
                                      "Dimensions":"%i %i" % (mesh.nElements_global,mesh.nNodes_element)})
            geometry    = SubElement(self.arGrid,"Geometry",{"Type":"XYZ"})
            nodes       = SubElement(geometry,"DataItem",
                                     {"Format":ar.dataItemFormat,
                                      "DataType":"Float",
                                      "Precision":"8",
                                      "Dimensions":"%i %i" % (dofMap.nDOF,3)})
            if ar.hdfFile != None:
                elements.text = ar.hdfFilename+":/elements"+spaceSuffix+`tCount`
                nodes.text    = ar.hdfFilename+":/nodes"+spaceSuffix+`tCount`
                if init or meshChanged:
                    #this will fail if elements_dgp1 already exists
                    ar.hdfFile.createArray("/",'elements'+spaceSuffix+`tCount`,dofMap.l2g)
                    #bad
                    dgnodes = numpy.zeros((dofMap.nDOF,3),'d')
                    for eN in range(mesh.nElements_global):
                        for nN in range(mesh.nNodes_element):
                            dgnodes[dofMap.l2g[eN,nN],:]=mesh.nodeArray[mesh.elementNodesArray[eN,nN]]
                    #make more pythonic loop
                    ar.hdfFile.createArray("/",'nodes'+spaceSuffix+`tCount`,dgnodes)
            else:
                SubElement(elements,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/elements"+spaceSuffix+`tCount`+".txt"})
                SubElement(nodes,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/nodes"+spaceSuffix+`tCount`+".txt"})
                if init or meshChanged:
                    numpy.savetxt(ar.textDataDir+"/elements"+spaceSuffix+`tCount`+".txt",dofMap.l2g,fmt='%d')
                    #bad
                    dgnodes = numpy.zeros((dofMap.nDOF,3),'d')
                    for eN in range(mesh.nElements_global):
                        for nN in range(mesh.nNodes_element):
                            dgnodes[dofMap.l2g[eN,nN],:]=mesh.nodeArray[mesh.elementNodesArray[eN,nN]]
                    #make more pythonic loop
                    numpy.savetxt(ar.textDataDir+"/nodes"+spaceSuffix+`tCount`+".txt",dgnodes)
                    
                #
            #hdfile
        #need to write a grid
        return self.arGrid
    #def
    def writeMeshXdmf_DGP2Lagrange(self,ar,name,mesh,spaceDim,dofMap,CGDOFMap,t=0.0,
                                   init=False,meshChanged=False,arGrid=None,tCount=0):
        #write out basic geometry if not already done?
        mesh.writeMeshXdmf(ar,"Spatial_Domain",t,init,meshChanged,tCount=tCount)
        #now try to write out a mesh that matches DG P1 layout
        #first duplicate geometry points etc then try to save space
        #mwf debug
        #import pdb
        #pdb.set_trace()
        spaceSuffix = "_dgp2_Lagrange"
        gridName = self.setGridCollectionAndGridElements(init,ar,arGrid,t,spaceSuffix)

        if self.arGrid == None or self.arTime.get('Value') != str(t):
            if spaceDim == 1:
                Xdmf_ElementTopology = "Edge_3"
            elif spaceDim == 2:
                Xdmf_ElementTopology = "Tri_6"
            elif spaceDim == 3:
                Xdmf_ElementTopology = "Tet_10"
            self.arGrid = SubElement(self.arGridCollection,"Grid",{"Name":gridName,"GridType":"Uniform"})
            self.arTime = SubElement(self.arGrid,"Time",{"Value":str(t)})
            topology    = SubElement(self.arGrid,"Topology",
                                     {"Type":Xdmf_ElementTopology,
                                      "NumberOfElements":str(mesh.nElements_global)})
            elements    = SubElement(topology,"DataItem",
                                     {"Format":ar.dataItemFormat,
                                      "DataType":"Int",
                                      "Dimensions":"%i %i" % dofMap.l2g.shape})
            geometry    = SubElement(self.arGrid,"Geometry",{"Type":"XYZ"})
            #try to use fancy functions later
            nodes       = SubElement(geometry,"DataItem",
                                     {"Format":ar.dataItemFormat,
                                      "DataType":"Float",
                                      "Precision":"8",
                                      "Dimensions":"%i %i" % (dofMap.nDOF,3)})
            if ar.hdfFile != None:
                elements.text = ar.hdfFilename+":/elements"+spaceSuffix+`tCount`
                nodes.text    = ar.hdfFilename+":/nodes"+spaceSuffix+`tCount`
                if init or meshChanged:
                    #this will fail if elements_dgp1 already exists
                    ar.hdfFile.createArray("/",'elements'+spaceSuffix+`tCount`,dofMap.l2g)
                    #bad
                    dgnodes = numpy.zeros((dofMap.nDOF,3),'d')
                    for eN in range(mesh.nElements_global):
                        for i in range(dofMap.l2g.shape[1]):
                            #for nN in range(mesh.nNodes_element):
                            #dgnodes[dofMap.l2g[eN,nN],:]=mesh.nodeArray[mesh.elementNodesArray[eN,nN]]
                            #now changed lagrange nodes to hold all nodes
                            dgnodes[dofMap.l2g[eN,i],:]= CGDOFMap.lagrangeNodesArray[CGDOFMap.l2g[eN,i],:]
                        #next loop over extra dofs on element and write out
                        #for nN in range(mesh.nNodes_element,dofMap.l2g.shape[1]):
                        #    dgnodes[dofMap.l2g[eN,nN],:]= CGDOFMap.lagrangeNodesArray[CGDOFMap.l2g[eN,nN]-mesh.nNodes_global,:]
                            
                    #make more pythonic loop
                    ar.hdfFile.createArray("/",'nodes'+spaceSuffix+`tCount`,dgnodes)
            else:
                SubElement(elements,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/elements"+spaceSuffix+`tCount`+".txt"})
                SubElement(nodes,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/nodes"+spaceSuffix+`tCount`+".txt"})
                if init or meshChanged:
                    numpy.savetxt(ar.textDataDir+"/elements"+spaceSuffix+`tCount`+".txt",dofMap.l2g,fmt='%d')
                    #bad
                    dgnodes = numpy.zeros((dofMap.nDOF,3),'d')
                    for eN in range(mesh.nElements_global):
                        for i in range(dofMap.l2g.shape[1]):
                            #for nN in range(mesh.nNodes_element):
                            #dgnodes[dofMap.l2g[eN,nN],:]=mesh.nodeArray[mesh.elementNodesArray[eN,nN]]
                            dgnodes[dofMap.l2g[eN,i],:]= CGDOFMap.lagrangeNodesArray[CGDOFMap.l2g[eN,i],:]
                        #next loop over extra dofs on element and write out
                        #for nN in range(mesh.nNodes_element,dofMap.l2g.shape[1]):
                        #    dgnodes[dofMap.l2g[eN,nN],:]= CGDOFMap.lagrangeNodesArray[CGDOFMap.l2g[eN,nN]-mesh.nNodes_global,:]
                    
                    #make more pythonic loop
                    numpy.savetxt(ar.textDataDir+"/nodes"+spaceSuffix+`tCount`+".txt",dgnodes)
                    
                #
            #hdfile
        #need to write a grid
        return self.arGrid

    def writeMeshXdmf_C0P2Lagrange(self,ar,name,mesh,spaceDim,dofMap,t=0.0,
                                   init=False,meshChanged=False,arGrid=None,tCount=0):
        """
        TODO: test new lagrangeNodes convention for 2d,3d, and concatNow=False
        """
        #write out basic geometry if not already done?
        mesh.writeMeshXdmf(ar,"Spatial_Domain",t,init,meshChanged,tCount=tCount)
        spaceSuffix = "_c0p2_Lagrange"
        gridName = self.setGridCollectionAndGridElements(init,ar,arGrid,t,spaceSuffix)
        if self.arGrid == None or self.arTime.get('Value') != str(t):
            if spaceDim == 1:
                Xdmf_ElementTopology = "Edge_3"
            elif spaceDim == 2:
                Xdmf_ElementTopology = "Tri_6"
            elif spaceDim == 3:
                Xdmf_ElementTopology = "Tet_10"
            self.arGrid = SubElement(self.arGridCollection,"Grid",{"Name":gridName,"GridType":"Uniform"})
            self.arTime = SubElement(self.arGrid,"Time",{"Value":str(t)})
            lagrangeNodesArray = dofMap.lagrangeNodesArray
            topology = SubElement(self.arGrid,"Topology",
                                  {"Type":Xdmf_ElementTopology,
                                   "NumberOfElements":str(mesh.nElements_global)})
            elements = SubElement(topology,"DataItem",
                                  {"Format":ar.dataItemFormat,
                                   "DataType":"Int",
                                   "Dimensions":"%i %i" % dofMap.l2g.shape})
            geometry = SubElement(self.arGrid,"Geometry",{"Type":"XYZ"})
            concatNow = True
            if concatNow:
                allNodes    = SubElement(geometry,"DataItem",
                                         {"Format":ar.dataItemFormat,
                                          "DataType":"Float",
                                          "Precision":"8",
                                          "Dimensions":"%i %i" % (dofMap.nDOF,3)})
            else:
                allNodes    = SubElement(geometry,"DataItem",
                                         {"Function":"JOIN( $0 ; $1 )",
                                          "DataType":"Float",
                                          "Precision":"8",
                                          "Dimensions":"%i %i" % (dofMap.nDOF,3)})
                #nodes    = SubElement(allNodes,"DataItem",
                #                      {"Format":ar.dataItemFormat,
                #                       "DataType":"Float",
                #                       "Precision":"8",
                #                       "Dimensions":"%i %i" % mesh.elementNodesArray.shape})
                lagrangeNodes    = SubElement(allNodes,"DataItem",
                                              {"Format":ar.dataItemFormat,
                                               "DataType":"Float",
                                               "Precision":"8",
                                               "Dimensions":"%i %i" % lagrangeNodesArray.shape})
            if ar.hdfFile != None:
                elements.text = ar.hdfFilename+":/elements"+spaceSuffix+`tCount`
                if concatNow:
                    allNodes.text = ar.hdfFilename+":/nodes"+spaceSuffix+`tCount`
                    import copy
                    if spaceDim == 3:#
                        #mwf 10/19/09 orig no longer works because of changes in orderings for
                        #parallel
                        elements=copy.deepcopy(dofMap.l2g)
                        #proteus stores 3d dof as
                        #|n0,n1,n2,n3|(n0,n1),(n1,n2),(n2,n3)|(n0,n2),(n1,n3)|(n0,n3)|
                        #looks like xdmf wants them as 
                        #|n0,n1,n2,n3|(n0,n1),(n1,n2),(n0,n2) (n0,n3),(n1,n3) (n2,n3)|
                        for eN in range(mesh.nElements_global):
                            elements[eN,4+2] = dofMap.l2g[eN,4+3]
                            elements[eN,4+3] = dofMap.l2g[eN,4+5]
                            elements[eN,4+5] = dofMap.l2g[eN,4+2]
                    else:
                        elements=dofMap.l2g
                    if init or meshChanged:
                        ar.hdfFile.createArray("/",'elements'+spaceSuffix+`tCount`,elements)
                        #print "element nodes",elements
                        #mwf orig ar.hdfFile.createArray("/",'nodes'+spaceSuffix+`tCount`,numpy.concatenate((mesh.nodeArray,lagrangeNodesArray)))
                        #ar.hdfFile.createArray("/",'nodes'+spaceSuffix+`tCount`,numpy.concatenate((mesh.nodeArray.flat[:3*mesh.nNodes_owned],
                        #                                                                           lagrangeNodesArray.flat[:3*mesh.nElements_owned],
                        #                                                                           mesh.nodeArray.flat[3*mesh.nNodes_owned:3*mesh.nNodes_global],
                        #                                                                           lagrangeNodesArray.flat[3*mesh.nElements_owned:3*mesh.nElements_global])))
                        ar.hdfFile.createArray("/",'nodes'+spaceSuffix+`tCount`,lagrangeNodesArray)
                        #mwf debug 
                        #print "nodes ",numpy.concatenate((mesh.nodeArray,lagrangeNodesArray))
                        #print "nodes ",numpy.concatenate((mesh.nodeArray.flat[:3*mesh.nNodes_owned],
                        #                                  lagrangeNodesArray.flat[:3*mesh.nElements_owned],
                        #                                  mesh.nodeArray.flat[3*mesh.nNodes_owned:3*mesh.nNodes_global],
                        #                                  lagrangeNodesArray.flat[3*mesh.nElements_owned:3*mesh.nElements_global]))
                else:
                    nodes.text = ar.hdfFilename+":/nodes"+spaceSuffix+`tCount`
                    lagrangeNodes.text = ar.hdfFilename+":/lagrangeNodes"+spaceSuffix+`tCount`
                    if init or meshChanged:
                        ar.hdfFile.createArray("/",'elements'+spaceSuffix+`tCount`,dofMap.l2g)
                        #ar.hdfFile.createArray("/",'nodes'+spaceSuffix+`tCount`,mesh.nodeArray)
                        ar.hdfFile.createArray("/",'lagrangeNodes'+spaceSuffix+`tCount`,lagrangeNodesArray)

            else:
                SubElement(elements,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/elements"+spaceSuffix+`tCount`+".txt"})
                if concatNow:
                    SubElement(allNodes,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/nodes"+spaceSuffix+`tCount`+".txt"})
                    if init or meshChanged:
                        numpy.savetxt(ar.textDataDir+"/elements"+spaceSuffix+`tCount`+".txt",dofMap.l2g,fmt='%d')
                        numpy.savetxt(ar.textDataDir+"/nodes"+spaceSuffix+`tCount`+".txt",lagrangeNodesArray)
                else:
                    #SubElement(nodes,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/nodes"+spaceSuffix+`tCount`+".txt"})
                    SubElement(lagrangeNodes,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/lagrangeNodes"+spaceSuffix+`tCount`+".txt"})
                    if init or meshChanged:
                        numpy.savetxt(ar.textDataDir+"/elements"+spaceSuffix+`tCount`+".txt",dofMap.l2g,fmt='%d')
                        #numpy.savetxt(ar.textDataDir+"/nodes"+spaceSuffix+`tCount`+".txt",mesh.nodeArray)
                        numpy.savetxt(ar.textDataDir+"/lagrangeNodes"+spaceSuffix+`tCount`+".txt",lagrangeNodesArray)
                    #

                #
            #
        #need to write a grid
        return self.arGrid

    def writeMeshXdmf_C0Q2Lagrange(self,ar,name,mesh,spaceDim,dofMap,t=0.0,init=False,meshChanged=False,arGrid=None,tCount=0):
        """
        TODO: test new lagrangeNodes convention for 2d,3d, and concatNow=False
        """
        #write out basic geometry if not already done?
        #mesh.writeMeshXdmf(ar,"Spatial_Domain",t,init,meshChanged,tCount=tCount)
        spaceSuffix = "_c0q2_Lagrange"
        gridName = self.setGridCollectionAndGridElements(init,ar,arGrid,t,spaceSuffix)
        if self.arGrid == None or self.arTime.get('Value') != str(t):
            if spaceDim == 1:
                 print "No writeMeshXdmf_C0Q2Lagrange for 1D" 
                 return 0
            elif spaceDim == 2:
                 print "No writeMeshXdmf_C0Q2Lagrange for 2D" 
                 return 0
            elif spaceDim == 3:
                Xdmf_ElementTopology = "Hexahedron"

                e2s=[[0,8,20,11,12,21,26,24], [8,1,9,20,21,13,22,26], [11,20,10,3,24,26,23,15], [20,9,2,10,26,22,14,23],
                     [12,21,26,24,4,16,25,19], [21,13,22,26,16,5,17,25], [24,26,23,15,19,25,18,7], [26,22,14,23,25,17,6,18] ]

                l2g = numpy.zeros((8*mesh.nElements_global,8),'i')
                for eN in range(mesh.nElements_global):
                   dofs=dofMap.l2g[eN,:]
                   for i in range(8): #loop over subelements
                      for j in range(8): # loop over nodes of subelement
                         l2g[8*eN+i,j] = dofs[e2s[i][j]]

            self.arGrid = SubElement(self.arGridCollection,"Grid",{"Name":gridName,"GridType":"Uniform"})
            self.arTime = SubElement(self.arGrid,"Time",{"Value":str(t)})

            lagrangeNodesArray = dofMap.lagrangeNodesArray      

            topology = SubElement(self.arGrid,"Topology",
                                  {"Type":Xdmf_ElementTopology,
                                   "NumberOfElements":str(l2g.shape[0])})
            elements = SubElement(topology,"DataItem",
                                  {"Format":ar.dataItemFormat,
                                   "DataType":"Int",
                                   "Dimensions":"%i %i" % l2g.shape})
            geometry = SubElement(self.arGrid,"Geometry",{"Type":"XYZ"})

            concatNow = True
            if concatNow:

                allNodes = SubElement(geometry,"DataItem",
                                         {"Format":ar.dataItemFormat,
                                          "DataType":"Float",
                                          "Precision":"8",
                                          "Dimensions":"%i %i" % lagrangeNodesArray.shape})
            else:
                allNodes    = SubElement(geometry,"DataItem",
                                         {"Function":"JOIN( $0 ; $1 )",
                                          "DataType":"Float",
                                          "Precision":"8",
                                          "Dimensions":"%i %i" % lagrangeNodesArray.shape})
                #nodes    = SubElement(allNodes,"DataItem",
                #                      {"Format":ar.dataItemFormat,
                #                       "DataType":"Float",
                #                       "Precision":"8",
                #                       "Dimensions":"%i %i" % mesh.elementNodesArray.shape})
                lagrangeNodes    = SubElement(allNodes,"DataItem",
                                              {"Format":ar.dataItemFormat,
                                               "DataType":"Float",
                                               "Precision":"8",
                                               "Dimensions":"%i %i" % lagrangeNodesArray.shape}) 

            if ar.hdfFile != None:
                elements.text = ar.hdfFilename+":/elements"+spaceSuffix+`tCount`
                if concatNow:
                    allNodes.text = ar.hdfFilename+":/nodes"+spaceSuffix+`tCount`
                    import copy
                    if init or meshChanged:
                        ar.hdfFile.createArray("/",'elements'+spaceSuffix+`tCount`,l2g)
                        #print "element nodes",elements
                        #mwf orig ar.hdfFile.createArray("/",'nodes'+spaceSuffix+`tCount`,numpy.concatenate((mesh.nodeArray,lagrangeNodesArray)))
                        #ar.hdfFile.createArray("/",'nodes'+spaceSuffix+`tCount`,numpy.concatenate((mesh.nodeArray.flat[:3*mesh.nNodes_owned],
                        #                                                                           lagrangeNodesArray.flat[:3*mesh.nElements_owned],
                        #                                                                           mesh.nodeArray.flat[3*mesh.nNodes_owned:3*mesh.nNodes_global],
                        #                                                                           lagrangeNodesArray.flat[3*mesh.nElements_owned:3*mesh.nElements_global])))
                        ar.hdfFile.createArray("/",'nodes'+spaceSuffix+`tCount`,lagrangeNodesArray)
                        #mwf debug 
                        #print "nodes ",numpy.concatenate((mesh.nodeArray,lagrangeNodesArray))
                        #print "nodes ",numpy.concatenate((mesh.nodeArray.flat[:3*mesh.nNodes_owned],
                        #                                  lagrangeNodesArray.flat[:3*mesh.nElements_owned],
                        #                                  mesh.nodeArray.flat[3*mesh.nNodes_owned:3*mesh.nNodes_global],
                        #                                  lagrangeNodesArray.flat[3*mesh.nElements_owned:3*mesh.nElements_global]))
                else:
                    nodes.text = ar.hdfFilename+":/nodes"+spaceSuffix+`tCount`
                    lagrangeNodes.text = ar.hdfFilename+":/lagrangeNodes"+spaceSuffix+`tCount`
                    if init or meshChanged:
                        ar.hdfFile.createArray("/",'elements'+spaceSuffix+`tCount`,l2g)
                        #ar.hdfFile.createArray("/",'nodes'+spaceSuffix+`tCount`,mesh.nodeArray)
                        ar.hdfFile.createArray("/",'lagrangeNodes'+spaceSuffix+`tCount`,lagrangeNodesArray)

            else:
                SubElement(elements,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/elements"+spaceSuffix+`tCount`+".txt"})
                if concatNow:
                    SubElement(allNodes,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/nodes"+spaceSuffix+`tCount`+".txt"})
                    if init or meshChanged:
                        numpy.savetxt(ar.textDataDir+"/elements"+spaceSuffix+`tCount`+".txt",l2g,fmt='%d')
                        numpy.savetxt(ar.textDataDir+"/nodes"+spaceSuffix+`tCount`+".txt",lagrangeNodesArray)
                else:
                    #SubElement(nodes,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/nodes"+spaceSuffix+`tCount`+".txt"})
                    SubElement(lagrangeNodes,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/lagrangeNodes"+spaceSuffix+`tCount`+".txt"})
                    if init or meshChanged:
                        numpy.savetxt(ar.textDataDir+"/elements"+spaceSuffix+`tCount`+".txt",l2g,fmt='%d')
                        #numpy.savetxt(ar.textDataDir+"/nodes"+spaceSuffix+`tCount`+".txt",mesh.nodeArray)
                        numpy.savetxt(ar.textDataDir+"/lagrangeNodes"+spaceSuffix+`tCount`+".txt",lagrangeNodesArray)
                    #

                #
            #
        #need to write a grid
        return self.arGrid


    
    def writeMeshXdmf_CrouzeixRaviartP1(self,ar,mesh,spaceDim,dofMap,t=0.0,
                                        init=False,meshChanged=False,arGrid=None,tCount=0):
        """
        Write out nonconforming P1 approximation
        Write out as a (discontinuous) Lagrange P1 function to make visualization easier
        and dof's as face centered data on original grid
        """
        #write out basic geometry if not already done?
        mesh.writeMeshXdmf(ar,"Spatial_Domain",t,init,meshChanged,tCount=tCount)
        #now try to write out a mesh that matches DG P1 layout

        spaceSuffix = "_ncp1_CrouzeixRaviart"
        gridName = self.setGridCollectionAndGridElements(init,ar,arGrid,t,spaceSuffix)

        if self.arGrid == None or self.arTime.get('Value') != str(t):
            if spaceDim == 1:
                Xdmf_ElementTopology = "Polyline"
            elif spaceDim == 2:
                Xdmf_ElementTopology = "Triangle"
            elif spaceDim == 3:
                Xdmf_ElementTopology = "Tetrahedron"
            Xdmf_NodesPerElement = spaceDim+1
            Xdmf_NumberOfElements= mesh.nElements_global
            self.arGrid = SubElement(self.arGridCollection,"Grid",{"Name":gridName,"GridType":"Uniform"})
            self.arTime = SubElement(self.arGrid,"Time",{"Value":str(t)})
            topology    = SubElement(self.arGrid,"Topology",
                                     {"Type":Xdmf_ElementTopology,
                                      "NumberOfElements":str(Xdmf_NumberOfElements)})
            elements    = SubElement(topology,"DataItem",
                                     {"Format":ar.dataItemFormat,
                                      "DataType":"Int",
                                      "Dimensions":"%i %i" % (Xdmf_NumberOfElements,Xdmf_NodesPerElement)})
            geometry    = SubElement(self.arGrid,"Geometry",{"Type":"XYZ"})
            nodes       = SubElement(geometry,"DataItem",
                                     {"Format":ar.dataItemFormat,
                                      "DataType":"Float",
                                      "Dimensions":"%i %i" % (Xdmf_NumberOfElements*Xdmf_NodesPerElement,3)})
            if ar.hdfFile != None:
                elements.text = ar.hdfFilename+":/elements"+spaceSuffix+`tCount`
                nodes.text    = ar.hdfFilename+":/nodes"+spaceSuffix+`tCount`
                if init or meshChanged:
                    #simple dg l2g mapping
                    dg_l2g = numpy.arange(Xdmf_NumberOfElements*Xdmf_NodesPerElement,dtype='i').reshape((Xdmf_NumberOfElements,Xdmf_NodesPerElement))
                    ar.hdfFile.createArray("/",'elements'+spaceSuffix+`tCount`,dg_l2g)

                    dgnodes = numpy.reshape(mesh.nodeArray[mesh.elementNodesArray],(Xdmf_NumberOfElements*Xdmf_NodesPerElement,3))
                    ar.hdfFile.createArray("/",'nodes'+spaceSuffix+`tCount`,dgnodes)
            else:
                SubElement(elements,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/elements"+spaceSuffix+`tCount`+".txt"})
                SubElement(nodes,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/nodes"+spaceSuffix+`tCount`+".txt"})
                if init or meshChanged:
                    dg_l2g = numpy.arange(Xdmf_NumberOfElements*Xdmf_NodesPerElement,dtype='i').reshape((Xdmf_NumberOfElements,Xdmf_NodesPerElement))
                    numpy.savetxt(ar.textDataDir+"/elements"+spaceSuffix+`tCount`+".txt",dg_l2g,fmt='%d')

                    dgnodes = numpy.reshape(mesh.nodeArray[mesh.elementNodesArray],(Xdmf_NumberOfElements*Xdmf_NodesPerElement,3))
                    numpy.savetxt(ar.textDataDir+"/nodes"+spaceSuffix+`tCount`+".txt",dgnodes)
                    
                #
            #hdfile
        #need to write a grid
        return self.arGrid
    #def
    def writeFunctionXdmf_DGP1Lagrange(self,ar,u,tCount=0,init=True):
        attribute = SubElement(self.arGrid,"Attribute",{"Name":u.name,
                                                        "AttributeType":"Scalar",
                                                        "Center":"Node"})
        values    = SubElement(attribute,"DataItem",
                               {"Format":ar.dataItemFormat,
                                "DataType":"Float",
                                "Precision":"8",
                                "Dimensions":"%i" % (u.nDOF_global,)})
        if ar.hdfFile != None:
            values.text = ar.hdfFilename+":/"+u.name+str(tCount)
            ar.hdfFile.createArray("/",u.name+str(tCount),u.dof)
        else:
            numpy.savetxt(ar.textDataDir+"/"+u.name+str(tCount)+".txt",u.dof)
            SubElement(values,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/"+u.name+str(tCount)+".txt"}) 
        
    def writeFunctionXdmf_DGP2Lagrange(self,ar,u,tCount=0,init=True):
        attribute = SubElement(self.arGrid,"Attribute",{"Name":u.name,
                                                        "AttributeType":"Scalar",
                                                        "Center":"Node"})
        values    = SubElement(attribute,"DataItem",
                               {"Format":ar.dataItemFormat,
                                "DataType":"Float",
                                "Precision":"8",
                                "Dimensions":"%i" % (u.nDOF_global,)})
        if ar.hdfFile != None:
            values.text = ar.hdfFilename+":/"+u.name+str(tCount)
            ar.hdfFile.createArray("/",u.name+str(tCount),u.dof)
        else:
            numpy.savetxt(ar.textDataDir+"/"+u.name+str(tCount)+".txt",u.dof)
            SubElement(values,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/"+u.name+str(tCount)+".txt"}) 
    def writeFunctionXdmf_CrouzeixRaviartP1(self,ar,u,tCount=0,init=True):
        Xdmf_NumberOfElements = u.femSpace.mesh.nElements_global
        Xdmf_NodesPerElement  = u.femSpace.mesh.nNodes_element
        
        name = u.name.replace(' ','_')
        #if writing as dgp1
        attribute = SubElement(self.arGrid,"Attribute",{"Name":name,
                                                        "AttributeType":"Scalar",
                                                        "Center":"Node"})
        values    = SubElement(attribute,"DataItem",
                               {"Format":ar.dataItemFormat,
                                "DataType":"Float",
                                "Dimensions":"%i" % (Xdmf_NumberOfElements*Xdmf_NodesPerElement,)})
        u_tmp = numpy.zeros((Xdmf_NumberOfElements*Xdmf_NodesPerElement,),'d')
        if u.femSpace.nSpace_global == 1:
            for eN in range(Xdmf_NumberOfElements):
                #dof associated with face id, so opposite usual C0P1 ordering here
                u_tmp[eN*Xdmf_NodesPerElement + 0] = u.dof[u.femSpace.dofMap.l2g[eN,1]]
                u_tmp[eN*Xdmf_NodesPerElement + 1] = u.dof[u.femSpace.dofMap.l2g[eN,0]]
        elif u.femSpace.nSpace_global == 2:
            for eN in range(Xdmf_NumberOfElements):
                #assume vertex associated with face across from it
                u_tmp[eN*Xdmf_NodesPerElement + 0] = u.dof[u.femSpace.dofMap.l2g[eN,1]]
                u_tmp[eN*Xdmf_NodesPerElement + 0]+= u.dof[u.femSpace.dofMap.l2g[eN,2]]
                u_tmp[eN*Xdmf_NodesPerElement + 0]-= u.dof[u.femSpace.dofMap.l2g[eN,0]]
                
                u_tmp[eN*Xdmf_NodesPerElement + 1] = u.dof[u.femSpace.dofMap.l2g[eN,0]]
                u_tmp[eN*Xdmf_NodesPerElement + 1]+= u.dof[u.femSpace.dofMap.l2g[eN,2]]
                u_tmp[eN*Xdmf_NodesPerElement + 1]-= u.dof[u.femSpace.dofMap.l2g[eN,1]]

                u_tmp[eN*Xdmf_NodesPerElement + 2] = u.dof[u.femSpace.dofMap.l2g[eN,0]]
                u_tmp[eN*Xdmf_NodesPerElement + 2]+= u.dof[u.femSpace.dofMap.l2g[eN,1]]
                u_tmp[eN*Xdmf_NodesPerElement + 2]-= u.dof[u.femSpace.dofMap.l2g[eN,2]]
        else:
            for eN in range(Xdmf_NumberOfElements):
                for i in range(Xdmf_NodesPerElement):
                    u_tmp[eN*Xdmf_NodesPerElement + i] = u.dof[u.femSpace.dofMap.l2g[eN,i]]*(1.0-float(u.femSpace.nSpace_global)) + \
                                                         sum([u.dof[u.femSpace.dofMap.l2g[eN,j]] for j in range(Xdmf_NodesPerElement) if j != i])
            
                
        if ar.hdfFile != None:
            values.text = ar.hdfFilename+":/"+name+str(tCount)
            ar.hdfFile.createArray("/",name+str(tCount),u_tmp)
        else:
            numpy.savetxt(ar.textDataDir+"/"+name+str(tCount)+".txt",u_tmp)
            SubElement(values,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/"+name+str(tCount)+".txt"}) 
        
    
        #if writing as face centered (true nc p1)
        #could be under self.arGrid or mesh.arGrid
        grid = u.femSpace.mesh.arGrid #self.arGrid
        attribute_dof = SubElement(grid,"Attribute",{"Name":name+"_dof",
                                                     "AttributeType":"Scalar",
                                                     "Center":"Face"})
        values_dof    = SubElement(attribute_dof,"DataItem",
                               {"Format":ar.dataItemFormat,
                                "DataType":"Float",
                                "Dimensions":"%i" % (u.nDOF_global,)})
        if ar.hdfFile != None:
            values_dof.text = ar.hdfFilename+":/"+name+"_dof"+str(tCount)
            ar.hdfFile.createArray("/",name+"_dof"+str(tCount),u.dof)
        else:
            numpy.savetxt(ar.textDataDir+"/"+name+"_dof"+str(tCount)+".txt",u.dof)
            SubElement(values,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/"+name+"_dof"+str(tCount)+".txt"}) 
        
    def writeVectorFunctionXdmf_nodal(self,ar,uList,components,vectorName,spaceSuffix,tCount=0,init=True):
        concatNow=True
        nDOF_global = uList[components[0]].nDOF_global
        if concatNow:
            attribute = SubElement(self.arGrid,"Attribute",{"Name":vectorName,
                                                            "AttributeType":"Vector",
                                                            "Center":"Node"})
            values    = SubElement(attribute,"DataItem",
                                   {"Format":ar.dataItemFormat,
                                    "DataType":"Float",
                                    "Dimensions":"%i %i" % (nDOF_global,3)})
            u_dof = uList[components[0]].dof
            if len(components) < 2:
                v_dof = numpy.zeros(u_dof.shape,dtype='d')
            else:
                v_dof = uList[components[1]].dof
            if len(components) < 3:
                w_dof = numpy.zeros(u_dof.shape,dtype='d')
            else:
                w_dof = uList[components[2]].dof
            velocity = numpy.column_stack((u_dof,v_dof,w_dof))
            if ar.hdfFile != None:
                values.text = ar.hdfFilename+":/"+vectorName+str(tCount)
                ar.hdfFile.createArray("/",vectorName+str(tCount),velocity)
        else:            
            attribute = SubElement(self.arGrid,"Attribute",{"Name":vectorName,
                                                            "AttributeType":"Vector",
                                                            "Center":"Node"})
            if len(components) == 2:
                values    = SubElement(attribute,"DataItem",
                                       {"ItemType":"Function",
                                        "Function":"JOIN($0 , $1 , (0.0 * $1 ))",
                                        "Dimensions":"%i %i" % (nDOF_global,3)})
            elif len(components) == 3:
                values    = SubElement(attribute,"DataItem",
                                       {"ItemType":"Function",
                                        "Function":"JOIN($0 , $1 , $2)",
                                        "Dimensions":"%i %i" % (nDOF_global,3)})
            for ci in components:
                ReferenceString="""/Xdmf/Domain/Grid[@Name=\"Mesh%s\"]/Grid[%i]/Attribute[%i]/DataItem""" % (spaceSuffix,tCount+1,ci+1)
                component = SubElement(values,"DataItem",{"Reference":ReferenceString})
    def writeVectorFunctionXdmf_CrouzeixRaviartP1(self,ar,uList,components,spaceSuffix,tCount=0,init=True):
        return self.writeVectorFunctionXdmf_nodal(ar,uList,components,"_ncp1_CrouzeixRaviart",
                                                  tCount=tCount,init=init)
    
    def writeMeshXdmf_MonomialDGPK(self,ar,mesh,spaceDim,interpolationPoints,t=0.0,
                                   init=False,meshChanged=False,arGrid=None,tCount=0):
        """
        as a first cut, archive DG monomial spaces using same approach as for element quadrature
        arrays using x = interpolation points (which are just Gaussian quadrature points) as values
        """
        
        #write out basic geometry if not already done?
        mesh.writeMeshXdmf(ar,"Spatial_Domain",t,init,meshChanged,tCount=tCount)
        #now try to write out a mesh that is a collection of points per element
        spaceSuffix = "_dgpk_Monomial"
        
        gridName = self.setGridCollectionAndGridElements(init,ar,arGrid,t,spaceSuffix)

        assert len(interpolationPoints.shape) == 3 #make sure have right type of dictionary
        if self.arGrid == None or self.arTime.get('Value') != str(t):
            Xdmf_ElementTopology = "Polyvertex"
            Xdmf_NumberOfElements= interpolationPoints.shape[0]
            Xdmf_NodesPerElement = interpolationPoints.shape[1]
            Xdmf_NodesGlobal     = Xdmf_NumberOfElements*Xdmf_NodesPerElement
            
            self.arGrid = SubElement(self.arGridCollection,"Grid",{"Name":gridName,"GridType":"Uniform"})
            self.arTime = SubElement(self.arGrid,"Time",{"Value":str(t)})
            topology    = SubElement(self.arGrid,"Topology",
                                     {"Type":Xdmf_ElementTopology,
                                      "NumberOfElements":str(Xdmf_NumberOfElements),
                                      "NodesPerElement":str(Xdmf_NodesPerElement)})
            elements    = SubElement(topology,"DataItem",
                                     {"Format":ar.dataItemFormat,
                                      "DataType":"Int",
                                      "Dimensions":"%i %i" % (Xdmf_NumberOfElements,Xdmf_NodesPerElement)})
            geometry    = SubElement(self.arGrid,"Geometry",{"Type":"XYZ"})
            nodes       = SubElement(geometry,"DataItem",
                                     {"Format":ar.dataItemFormat,
                                      "DataType":"Float",
                                      "Precision":"8",
                                      "Dimensions":"%i %i" % (Xdmf_NodesGlobal,3)})
            if ar.hdfFile != None:
                elements.text = ar.hdfFilename+":/elements"+spaceSuffix+`tCount`
                nodes.text    = ar.hdfFilename+":/nodes"+spaceSuffix+`tCount`
                if init or meshChanged:
                    q_l2g = numpy.arange(Xdmf_NumberOfElements*Xdmf_NodesPerElement,dtype='i').reshape((Xdmf_NumberOfElements,Xdmf_NodesPerElement))
                    ar.hdfFile.createArray("/",'elements'+spaceSuffix+`tCount`,q_l2g)
                    ar.hdfFile.createArray("/",'nodes'+spaceSuffix+`tCount`,interpolationPoints.flat[:])
            else:
                SubElement(elements,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/elements"+spaceSuffix+`tCount`+".txt"})
                SubElement(nodes,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/nodes"+spaceSuffix+`tCount`+".txt"})
                if init or meshChanged:
                    q_l2g = numpy.arange(Xdmf_NumberOfElements*Xdmf_NodesPerElement,dtype='i').reshape((Xdmf_NumberOfElements,Xdmf_NodesPerElement))
                    numpy.savetxt(ar.textDataDir+"/elements"+spaceSuffix+`tCount`+".txt",q_l2g,fmt='%d')
                    numpy.savetxt(ar.textDataDir+"/nodes"+spaceSuffix+`tCount`+".txt",interpolationPoints.flat[:])
                    
                #
            #hdfile
        #need to write a grid
        return self.arGrid
    #def
    def writeFunctionXdmf_MonomialDGPK(self,ar,interpolationValues,name,tCount=0,init=True):
        """
        Different than usual FemFunction Write routines since saves values at interpolation points
        need to check on way to save dofs as well
        """
        assert len(interpolationValues.shape) == 2
        Xdmf_NodesGlobal = interpolationValues.shape[0]*interpolationValues.shape[1]
        attribute = SubElement(self.arGrid,"Attribute",{"Name":name,
                                                        "AttributeType":"Scalar",
                                                        "Center":"Node"})
        values    = SubElement(attribute,"DataItem",
                               {"Format":ar.dataItemFormat,
                                "DataType":"Float",
                                "Precision":"8",
                                "Dimensions":"%i" % (Xdmf_NodesGlobal,)})
        if ar.hdfFile != None:
            values.text = ar.hdfFilename+":/"+name+str(tCount)
            ar.hdfFile.createArray("/",name+str(tCount),interpolationValues.flat[:])
        else:
            numpy.savetxt(ar.textDataDir+"/"+name+str(tCount)+".txt",interpolationValues.flat[:])
            SubElement(values,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/"+name+str(tCount)+".txt"}) 
        
    def writeVectorFunctionXdmf_MonomialDGPK(self,ar,interpolationValues,name,tCount=0,init=True):
        """
        Different than usual FemFunction Write routines since saves values at interpolation points
        need to check on way to save dofs as well
        """
        assert len(interpolationValues.shape) == 3
        Xdmf_NodesGlobal = interpolationValues.shape[0]*interpolationValues.shape[1]
        Xdmf_NumberOfComponents = interpolationValues.shape[2]
        attribute = SubElement(self.arGrid,"Attribute",{"Name":name,
                                                        "AttributeType":"Vector",
                                                        "Center":"Node"})
        values    = SubElement(attribute,"DataItem",
                               {"Format":ar.dataItemFormat,
                                "DataType":"Float",
                                "Precision":"8",
                                "Dimensions":"%i %i" % (Xdmf_NodesGlobal,3)})#force 3d vector since points 3d
        #mwf brute force
        tmp = numpy.zeros((Xdmf_NodesGlobal,3),'d')
        tmp[:,:Xdmf_NumberOfComponents]=numpy.reshape(interpolationValues.flat,(Xdmf_NodesGlobal,Xdmf_NumberOfComponents))
        
        if ar.hdfFile != None:
            values.text = ar.hdfFilename+":/"+name+str(tCount)
            ar.hdfFile.createArray("/",name+str(tCount),tmp)
        else:
            numpy.savetxt(ar.textDataDir+"/"+name+str(tCount)+".txt",tmp)
            SubElement(values,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/"+name+str(tCount)+".txt"}) 
        
    
    def writeMeshXdmf_DGP0(self,ar,mesh,spaceDim,
                           t=0.0,init=False,meshChanged=False,arGrid=None,tCount=0):
        """
        as a first cut, archive piecewise constant space
        """
        
        #write out basic geometry if not already done?
        meshSpaceSuffix = "Spatial_Domain"
        mesh.writeMeshXdmf(ar,meshSpaceSuffix,t,init,meshChanged,tCount=tCount)
        #now try to write out a mesh that a constant per element
        spaceSuffix = "_dgp0"
        gridName = self.setGridCollectionAndGridElements(init,ar,arGrid,t,spaceSuffix)

        if self.arGrid == None or self.arTime.get('Value') != str(t):
            if spaceDim == 1:
                Xdmf_ElementTopology = "Polyline"
            elif spaceDim == 2:
                Xdmf_ElementTopology = "Triangle"
            elif spaceDim == 3:
                Xdmf_ElementTopology = "Tetrahedron"
            self.arGrid = SubElement(self.arGridCollection,"Grid",{"Name":gridName,"GridType":"Uniform"})
            self.arTime = SubElement(self.arGrid,"Time",{"Value":str(t)})
            topology    = SubElement(self.arGrid,"Topology",
                                     {"Type":Xdmf_ElementTopology,
                                      "NumberOfElements":str(mesh.nElements_global)})
            elements    = SubElement(topology,"DataItem",
                                     {"Format":ar.dataItemFormat,
                                      "DataType":"Int",
                                      "Dimensions":"%i %i" % (mesh.nElements_global,mesh.nNodes_element)})
            geometry    = SubElement(self.arGrid,"Geometry",{"Type":"XYZ"})
            nodes       = SubElement(geometry,"DataItem",
                                     {"Format":ar.dataItemFormat,
                                      "DataType":"Float",
                                      "Precision":"8",
                                      "Dimensions":"%i %i" % (mesh.nNodes_global,3)})
            #just reuse spatial mesh entries
            if ar.hdfFile != None:
                elements.text = ar.hdfFilename+":/elements"+meshSpaceSuffix+`tCount`
                nodes.text    = ar.hdfFilename+":/nodes"+meshSpaceSuffix+`tCount`
            else:
                SubElement(elements,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/elements"+meshSpaceSuffix+".txt"})
                SubElement(nodes,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/nodes"+meshSpaceSuffix+".txt"})
            #hdfile
        #need to write a grid
        return self.arGrid
    #def
    def writeFunctionXdmf_DGP0(self,ar,u,tCount=0,init=True):
        name = u.name.replace(' ','_')
        attribute = SubElement(self.arGrid,"Attribute",{"Name":name,
                                                        "AttributeType":"Scalar",
                                                        "Center":"Cell"})
        values    = SubElement(attribute,"DataItem",
                               {"Format":ar.dataItemFormat,
                                "DataType":"Float",
                                "Precision":"8",
                                "Dimensions":"%i" % (u.femSpace.elementMaps.mesh.nElements_global,)})
        if ar.hdfFile != None:
            values.text = ar.hdfFilename+":/"+name+str(tCount)
            ar.hdfFile.createArray("/",name+str(tCount),u.dof)
        else:
            numpy.savetxt(ar.textDataDir+"/"+name+str(tCount)+".txt",u.dof)
            SubElement(values,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/"+name+str(tCount)+".txt"}) 
 
    def writeVectorFunctionXdmf_DGP0(self,ar,uList,components,vectorName,tCount=0,init=True):
        concatNow=True
        if concatNow:
            attribute = SubElement(self.arGrid,"Attribute",{"Name":vectorName,
                                                            "AttributeType":"Vector",
                                                            "Center":"Cell"})
            values    = SubElement(attribute,"DataItem",
                                   {"Format":ar.dataItemFormat,
                                    "DataType":"Float",
                                    "Precision":"8",
                                    "Dimensions":"%i %i" % (u.femSpace.elementMaps.mesh.nElements_global,3)})
            u_dof = uList[components[0]].dof
            if len(components) < 2:
                v_dof = numpy.zeros(u_dof.shape,dtype='d')
            else:
                v_dof = uList[components[1]].dof
            if len(components) < 3:
                w_dof = numpy.zeros(u_dof.shape,dtype='d')
            else:
                w_dof = uList[components[2]].dof
            velocity = numpy.column_stack((u_dof,v_dof,w_dof))
            if ar.hdfFile != None:
                values.text = ar.hdfFilename+":/"+vectorName+str(tCount)
                ar.hdfFile.createArray("/",vectorName+str(tCount),velocity)
        else:            
            attribute = SubElement(uList[components[0]].femSpace.elementMaps.mesh.arGrid,"Attribute",{"Name":vectorName,
                                                                 "AttributeType":"Vector",
                                                                 "Center":"Cell"})
            if len(components) == 2:
                values    = SubElement(attribute,"DataItem",
                                       {"ItemType":"Function",
                                        "Function":"JOIN($0 , $1 , (0.0 * $1 ))",
                                        "Dimensions":"%i %i" % (u.femSpace.elementMaps.mesh.nElements_global,3)})
            elif len(components) == 3:
                values    = SubElement(attribute,"DataItem",
                                       {"ItemType":"Function",
                                        "Function":"JOIN($0 , $1 , $2)",
                                        "Dimensions":"%i %i" % (u.femSpace.elementMaps.mesh.nElements_global,3)})
            for ci in components:
                ReferenceString="/Xdmf/Domain/Grid/Grid[%i]/Attribute[%i]/DataItem" % (tCount+1,ci+1)
                component = SubElement(values,"DataItem",{"Reference":ReferenceString})
        
    #
    def writeMeshXdmf_P1Bubble(self,ar,mesh,spaceDim,dofMap,t=0.0,
                               init=False,meshChanged=False,arGrid=None,tCount=0):
        """
        represent P1 bubble space using just vertices for now
        """
        mesh.writeMeshXdmf(ar,"Spatial_Domain",t,init,meshChanged,tCount=tCount)
        spaceSuffix = "_c0p1_Bubble%s" % tCount
        gridName = self.setGridCollectionAndGridElements(init,ar,arGrid,t,spaceSuffix)

        if self.arGrid == None or self.arTime.get('Value') != str(t):
            if spaceDim == 1:
                Xdmf_ElementTopology = "Polyline"
            elif spaceDim == 2:
                Xdmf_ElementTopology = "Triangle"
            elif spaceDim == 3:
                Xdmf_ElementTopology = "Tetrahedron"
            Xdmf_NodesPerElement = spaceDim+1 #just handle vertices for now
            Xdmf_NumberOfElements= mesh.nElements_global
            Xdmf_NumberOfNodes   = mesh.nNodes_global
            self.arGrid = SubElement(self.arGridCollection,"Grid",{"Name":gridName,"GridType":"Uniform"})
            self.arTime = SubElement(self.arGrid,"Time",{"Value":str(t)})
            topology    = SubElement(self.arGrid,"Topology",
                                     {"Type":Xdmf_ElementTopology,
                                      "NumberOfElements":str(Xdmf_NumberOfElements)})
            elements    = SubElement(topology,"DataItem",
                                     {"Format":ar.dataItemFormat,
                                      "DataType":"Int",
                                      "Dimensions":"%i %i" % (Xdmf_NumberOfElements,Xdmf_NodesPerElement)})
            geometry    = SubElement(self.arGrid,"Geometry",{"Type":"XYZ"})
            nodes       = SubElement(geometry,"DataItem",
                                     {"Format":ar.dataItemFormat,
                                      "DataType":"Float",
                                      "Precision":"8",
                                      "Dimensions":"%i %i" % (Xdmf_NumberOfNodes,3)})
            if ar.hdfFile != None:
                elements.text = ar.hdfFilename+":/elements"+spaceSuffix+`tCount`
                nodes.text    = ar.hdfFilename+":/nodes"+spaceSuffix+`tCount`
                if init or meshChanged:
                    #c0p1 mapping for now
                    ar.hdfFile.createArray("/",'elements'+spaceSuffix+`tCount`,mesh.elementNodesArray)
                    ar.hdfFile.createArray("/",'nodes'+spaceSuffix+`tCount`,mesh.nodeArray)
            else:
                SubElement(elements,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/elements"+spaceSuffix+`tCount`+".txt"})
                SubElement(nodes,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/nodes"+spaceSuffix+`tCount`+".txt"})
                if init or meshChanged:
                    numpy.savetxt(ar.textDataDir+"/elements"+spaceSuffix+`tCount`+".txt",mesh.elementNodesArray,fmt='%d')
                    numpy.savetxt(ar.textDataDir+"/nodes"+spaceSuffix+`tCount`+".txt",mesh.nodeArray)
                    
                #
            #hdfile
        #need to write a grid
        return self.arGrid
    #def
    def writeFunctionXdmf_P1Bubble(self,ar,u,tCount=0,init=True):
        #just write out nodal part right now
        Xdmf_NumberOfElements = u.femSpace.mesh.nElements_global
        Xdmf_NumberOfNodes    = u.femSpace.mesh.nNodes_global
        name = u.name.replace(' ','_')
        
        #if writing as dgp1
        attribute = SubElement(self.arGrid,"Attribute",{"Name":name,
                                                        "AttributeType":"Scalar",
                                                        "Center":"Node"})
        values    = SubElement(attribute,"DataItem",
                               {"Format":ar.dataItemFormat,
                                "DataType":"Float",
                                "Precision":"8",
                                "Dimensions":"%i" % (Xdmf_NumberOfNodes,)})
        if ar.hdfFile != None:
            values.text = ar.hdfFilename+":/"+name+str(tCount)
            ar.hdfFile.createArray("/",name+str(tCount),u.dof[0:Xdmf_NumberOfNodes])
        else:
            numpy.savetxt(ar.textDataDir+"/"+name+str(tCount)+".txt",u.dof[0:Xdmf_NumberOfNodes])
            SubElement(values,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/"+name+str(tCount)+".txt"}) 

    def writeFunctionXdmf_C0P2Lagrange(self,ar,u,tCount=0,init=True):
        attribute = SubElement(self.arGrid,"Attribute",{"Name":u.name,
                                                 "AttributeType":"Scalar",
                                                 "Center":"Node"})
        values    = SubElement(attribute,"DataItem",
                               {"Format":ar.dataItemFormat,
                                "DataType":"Float",
                                "Precision":"8",
                                "Dimensions":"%i" % (u.nDOF_global,)})
        if ar.hdfFile != None:
            values.text = ar.hdfFilename+":/"+u.name+str(tCount)
            ar.hdfFile.createArray("/",u.name+str(tCount),u.dof)
        else:
            numpy.savetxt(ar.textDataDir+"/"+u.name+str(tCount)+".txt",u.dof)
            SubElement(values,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/"+u.name+str(tCount)+".txt"}) 
        
    def writeVectorFunctionXdmf_P1Bubble(self,ar,uList,components,vectorName,spaceSuffix,tCount=0,init=True):
        concatNow=True
        nDOF_global = uList[components[0]].femSpace.mesh.nNodes_global
        if concatNow:
            attribute = SubElement(self.arGrid,"Attribute",{"Name":vectorName,
                                                            "AttributeType":"Vector",
                                                            "Center":"Node"})
            values    = SubElement(attribute,"DataItem",
                                   {"Format":ar.dataItemFormat,
                                    "DataType":"Float",
                                    "Dimensions":"%i %i" % (nDOF_global,3)})
            u_dof = uList[components[0]].dof[0:nDOF_global]
            if len(components) < 2:
                v_dof = numpy.zeros((nDOF_global,),dtype='d')
            else:
                v_dof = uList[components[1]].dof[0:nDOF_global]
            if len(components) < 3:
                w_dof = numpy.zeros((nDOF_global,),dtype='d')
            else:
                w_dof = uList[components[2]].dof[0:nDOF_global]
            velocity = numpy.column_stack((u_dof,v_dof,w_dof))
            if ar.hdfFile != None:
                values.text = ar.hdfFilename+":/"+vectorName+str(tCount)
                ar.hdfFile.createArray("/",vectorName+str(tCount),velocity)
        else:            
            attribute = SubElement(self.arGrid,"Attribute",{"Name":vectorName,
                                                            "AttributeType":"Vector",
                                                            "Center":"Node"})
            if len(components) == 2:
                values    = SubElement(attribute,"DataItem",
                                       {"ItemType":"Function",
                                        "Function":"JOIN($0 , $1 , (0.0 * $1 ))",
                                        "Dimensions":"%i %i" % (nDOF_global,3)})
            elif len(components) == 3:
                values    = SubElement(attribute,"DataItem",
                                       {"ItemType":"Function",
                                        "Function":"JOIN($0 , $1 , $2)",
                                        "Dimensions":"%i %i" % (nDOF_global,3)})
            for ci in components:
                ReferenceString="""/Xdmf/Domain/Grid[@Name=\"Mesh%s\"]/Grid[%i]/Attribute[%i]/DataItem""" % (spaceSuffix,tCount+1,ci+1)
                component = SubElement(values,"DataItem",{"Reference":ReferenceString})
 
    def writeMeshXdmf_particles(self,ar,mesh,spaceDim,x,t=0.0,
                                init=False,meshChanged=False,arGrid=None,tCount=0,
                                spaceSuffix = "_particles"):
        """
        write out arbitrary set of points on a mesh

        """
        #spaceSuffix = "_particles"
        #write out basic geometry if not already done?
        mesh.writeMeshXdmf(ar,"Spatial_Domain",t,init,meshChanged,tCount=tCount)
        #now try to write out a mesh that is a collection of points per element
        
        gridName = self.setGridCollectionAndGridElements(init,ar,arGrid,t,spaceSuffix)
        nPoints = numpy.cumprod(x.shape)[-2]
        if self.arGrid == None or self.arTime.get('Value') != str(t):
            Xdmf_ElementTopology = "Polyvertex"
            Xdmf_NumberOfElements= nPoints
            Xdmf_NodesPerElement = 1
            Xdmf_NodesGlobal     = nPoints
            
            self.arGrid = SubElement(self.arGridCollection,"Grid",{"Name":gridName,"GridType":"Uniform"})
            self.arTime = SubElement(self.arGrid,"Time",{"Value":str(t)})
            topology    = SubElement(self.arGrid,"Topology",
                                     {"Type":Xdmf_ElementTopology,
                                      "NumberOfElements":str(Xdmf_NumberOfElements),
                                      "NodesPerElement":str(Xdmf_NodesPerElement)})
            elements    = SubElement(topology,"DataItem",
                                     {"Format":ar.dataItemFormat,
                                      "DataType":"Int",
                                      "Dimensions":"%i %i" % (Xdmf_NumberOfElements,Xdmf_NodesPerElement)})
            geometry    = SubElement(self.arGrid,"Geometry",{"Type":"XYZ"})
            nodes       = SubElement(geometry,"DataItem",
                                     {"Format":ar.dataItemFormat,
                                      "DataType":"Float",
                                      "Precision":"8",
                                      "Dimensions":"%i %i" % (Xdmf_NodesGlobal,3)})
            if ar.hdfFile != None:
                elements.text = ar.hdfFilename+":/elements"+spaceSuffix+`tCount`
                nodes.text    = ar.hdfFilename+":/nodes"+spaceSuffix+`tCount`
                if init or meshChanged:
                    #this will fail if elements_dgp1 already exists
                    #q_l2g = numpy.zeros((Xdmf_NumberOfElements,Xdmf_NodesPerElement),'i')
                    #brute force to start
                    #for eN in range(Xdmf_NumberOfElements):
                    #    for nN in range(Xdmf_NodesPerElement):
                    #        q_l2g[eN,nN] = eN*Xdmf_NodesPerElement + nN
                    #
                    q_l2g = numpy.arange(Xdmf_NumberOfElements*Xdmf_NodesPerElement,dtype='i').reshape((Xdmf_NumberOfElements,Xdmf_NodesPerElement))
                    ar.hdfFile.createArray("/",'elements'+spaceSuffix+`tCount`,q_l2g)
                    ar.hdfFile.createArray("/",'nodes'+spaceSuffix+`tCount`,x.flat[:])
            else:
                SubElement(elements,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/elements"+spaceSuffix+`tCount`+".txt"})
                SubElement(nodes,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/nodes"+spaceSuffix+`tCount`+".txt"})
                if init or meshChanged:
                    #this will fail if elements_dgp1 already exists
                    #q_l2g = numpy.zeros((Xdmf_NumberOfElements,Xdmf_NodesPerElement),'i')
                    #brute force to start
                    #for eN in range(Xdmf_NumberOfElements):
                    #    for nN in range(Xdmf_NodesPerElement):
                    #        q_l2g[eN,nN] = eN*Xdmf_NodesPerElement + nN
                    #
                    q_l2g = numpy.arange(Xdmf_NumberOfElements*Xdmf_NodesPerElement,dtype='i').reshape((Xdmf_NumberOfElements,Xdmf_NodesPerElement))
                    numpy.savetxt(ar.textDataDir+"/elements"+spaceSuffix+`tCount`+".txt",q_l2g,fmt='%d')
                    numpy.savetxt(ar.textDataDir+"/nodes"+spaceSuffix+`tCount`+".txt",x.flat[:])
                    
                #
            #hdfile
        #need to write a grid
        return self.arGrid
    #def
    def writeScalarXdmf_particles(self,ar,u,name,tCount=0,init=True):
        nPoints = numpy.cumprod(u.shape)[-1]
        
        Xdmf_NodesGlobal = nPoints
        attribute = SubElement(self.arGrid,"Attribute",{"Name":name,
                                                        "AttributeType":"Scalar",
                                                        "Center":"Node"})
        values    = SubElement(attribute,"DataItem",
                               {"Format":ar.dataItemFormat,
                                "DataType":"Float",
                                "Precision":"8",
                                "Dimensions":"%i" % (Xdmf_NodesGlobal,)})
        if ar.hdfFile != None:
            values.text = ar.hdfFilename+":/"+name+str(tCount)
            ar.hdfFile.createArray("/",name+str(tCount),u.flat[:])
        else:
            numpy.savetxt(ar.textDataDir+"/"+name+str(tCount)+".txt",u.flat[:])
            SubElement(values,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/"+name+str(tCount)+".txt"}) 
        
    def writeVectorXdmf_particles(self,ar,u,name,tCount=0,init=True):
        nPoints = numpy.cumprod(u.shape)[-2]

        Xdmf_NodesGlobal = nPoints
        Xdmf_NumberOfComponents = u.shape[-1]
        attribute = SubElement(self.arGrid,"Attribute",{"Name":name,
                                                        "AttributeType":"Vector",
                                                        "Center":"Node"})
        values    = SubElement(attribute,"DataItem",
                               {"Format":ar.dataItemFormat,
                                "DataType":"Float",
                                "Precision":"8",
                                "Dimensions":"%i %i" % (Xdmf_NodesGlobal,3)})#force 3d vector since points 3d
        #mwf brute force
        tmp = numpy.zeros((Xdmf_NodesGlobal,3),'d')
        tmp[:,:Xdmf_NumberOfComponents]=numpy.reshape(u.flat,(Xdmf_NodesGlobal,Xdmf_NumberOfComponents))
        
        if ar.hdfFile != None:
            values.text = ar.hdfFilename+":/"+name+str(tCount)
            ar.hdfFile.createArray("/",name+str(tCount),tmp)
        else:
            numpy.savetxt(ar.textDataDir+"/"+name+str(tCount)+".txt",tmp)
            SubElement(values,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/"+name+str(tCount)+".txt"}) 
        
    def writeMeshXdmf_LowestOrderMixed(self,ar,mesh,spaceDim,t=0.0,init=False,meshChanged=False,arGrid=None,tCount=0,
                                       spaceSuffix = "_RT0"):
        #write out basic geometry if not already done?
        mesh.writeMeshXdmf(ar,"Spatial_Domain",t,init,meshChanged,tCount=tCount)
        #now try to write out a mesh that matches RT0 velocity as dgp1 lagrange
        gridName = self.setGridCollectionAndGridElements(init,ar,arGrid,t,spaceSuffix)

        if self.arGrid == None or self.arTime.get('Value') != str(t):
            if spaceDim == 1:
                Xdmf_ElementTopology = "Polyline"
            elif spaceDim == 2:
                Xdmf_ElementTopology = "Triangle"
            elif spaceDim == 3:
                Xdmf_ElementTopology = "Tetrahedron"
            Xdmf_NodesPerElement = spaceDim+1
            Xdmf_NumberOfElements= mesh.nElements_global
            self.arGrid = SubElement(self.arGridCollection,"Grid",{"Name":gridName,"GridType":"Uniform"})
            self.arTime = SubElement(self.arGrid,"Time",{"Value":str(t)})
            topology    = SubElement(self.arGrid,"Topology",
                                     {"Type":Xdmf_ElementTopology,
                                      "NumberOfElements":str(Xdmf_NumberOfElements)})
            elements    = SubElement(topology,"DataItem",
                                     {"Format":ar.dataItemFormat,
                                      "DataType":"Int",
                                      "Dimensions":"%i %i" % (Xdmf_NumberOfElements,Xdmf_NodesPerElement)})
            geometry    = SubElement(self.arGrid,"Geometry",{"Type":"XYZ"})
            nodes       = SubElement(geometry,"DataItem",
                                     {"Format":ar.dataItemFormat,
                                      "DataType":"Float",
                                      "Dimensions":"%i %i" % (Xdmf_NumberOfElements*Xdmf_NodesPerElement,3)})
            if ar.hdfFile != None:
                elements.text = ar.hdfFilename+":/elements"+spaceSuffix+`tCount`
                nodes.text    = ar.hdfFilename+":/nodes"+spaceSuffix+`tCount`
                if init or meshChanged:
                    #simple dg l2g mapping
                    dg_l2g = numpy.arange(Xdmf_NumberOfElements*Xdmf_NodesPerElement,dtype='i').reshape((Xdmf_NumberOfElements,Xdmf_NodesPerElement))
                    ar.hdfFile.createArray("/",'elements'+spaceSuffix+`tCount`,dg_l2g)

                    dgnodes = numpy.reshape(mesh.nodeArray[mesh.elementNodesArray],(Xdmf_NumberOfElements*Xdmf_NodesPerElement,3))
                    ar.hdfFile.createArray("/",'nodes'+spaceSuffix+`tCount`,dgnodes)
            else:
                SubElement(elements,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/elements"+spaceSuffix+`tCount`+".txt"})
                SubElement(nodes,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/nodes"+spaceSuffix+`tCount`+".txt"})
                if init or meshChanged:
                    dg_l2g = numpy.arange(Xdmf_NumberOfElements*Xdmf_NodesPerElement,dtype='i').reshape((Xdmf_NumberOfElements,Xdmf_NodesPerElement))
                    numpy.savetxt(ar.textDataDir+"/elements"+spaceSuffix+`tCount`+".txt",dg_l2g,fmt='%d')

                    dgnodes = numpy.reshape(mesh.nodeArray[mesh.elementNodesArray],(Xdmf_NumberOfElements*Xdmf_NodesPerElement,3))
                    numpy.savetxt(ar.textDataDir+"/nodes"+spaceSuffix+`tCount`+".txt",dgnodes)
                    
                #
            #hdfile
        #need to write a grid
        return self.arGrid
    #def
    def writeVectorFunctionXdmf_LowestOrderMixed(self,ar,u,tCount=0,init=True,spaceSuffix="_RT0",name="velocity"):
        Xdmf_NodesGlobal = u.shape[0]*u.shape[1]
        Xdmf_NumberOfComponents = u.shape[2]
        Xdmf_StorageDim = 3
        
        #if writing as dgp1
        attribute = SubElement(self.arGrid,"Attribute",{"Name":name,
                                                        "AttributeType":"Vector",
                                                        "Center":"Node"})
        values    = SubElement(attribute,"DataItem",
                               {"Format":ar.dataItemFormat,
                                "DataType":"Float",
                                "Precision":"8",
                                "Dimensions":"%i %i" % (Xdmf_NodesGlobal,Xdmf_StorageDim)})#force 3d vector since points 3d
        tmp = numpy.zeros((Xdmf_NodesGlobal,Xdmf_StorageDim),'d')
        tmp[:,:Xdmf_NumberOfComponents]=numpy.reshape(u.flat,(Xdmf_NodesGlobal,Xdmf_NumberOfComponents))

        if ar.hdfFile != None:
            values.text = ar.hdfFilename+":/"+name+str(tCount)
            ar.hdfFile.createArray("/",name+str(tCount),tmp)
        else:
            numpy.savetxt(ar.textDataDir+"/"+name+str(tCount)+".txt",tmp)
            SubElement(values,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/"+name+str(tCount)+".txt"}) 
        #
        

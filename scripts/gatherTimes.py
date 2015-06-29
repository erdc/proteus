#!/usr/bin/env python

import numpy
import os
from xml.etree.ElementTree import *
#rather than import proteus
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

def gatherTimes(size, filename,dataDir='.',addname="_times"):
    """
    in case archiving failed to collect results from all times
    """
    import h5py
    xmlFile = open(filename+"_all"+`size`+".xmf","r")
    h5File = h5py.File(filename+".h5","r")
    tree = ElementTree(file=xmlFile)
    xmlFile.close()
    XDMF = tree.getroot()
    Domain = XDMF[0]
    for TemporalGridCollection in Domain:
        SpatialCollection = TemporalGridCollection[-1]
        Grids = SpatialCollection[:]
        tCount = int(Grids[0].attrib['Name'])
        del TemporalGridCollection[:]
        for i in range(tCount):
            dataset_name = TemporalGridCollection.attrib['Name']+"_"+`i`
            dataset_name = dataset_name.replace(" ","_")
            grid_array = h5File["/"+dataset_name]
            SpatialCollection=SubElement(TemporalGridCollection,"Grid",{"GridType":"Collection",
                                                                        "CollectionType":"Spatial"})
            time = SubElement(SpatialCollection,"Time",{"Value":grid_array.attrs['Time'],"Name":str(i)})
            for j in range(size):
                Grid = fromstring(grid_array[j])
                SpatialCollection.append(Grid)
    xmlFile = open(filename+"_all"+`size`+"_complete.xmf","w")
    indentXML(tree.getroot())
    tree.write(xmlFile)
    xmlFile.close()

if __name__ == '__main__':
    from optparse import OptionParser
    usage = ""
    parser = OptionParser(usage=usage)
    parser.add_option("-s","--size",
                      help="number of processors for run",
                      action="store",
                      type="int",
                      dest="size",
                      default=1)

    parser.add_option("-f","--filebase",
                      help="base name for storage files",
                      action="store",
                      type="string",
                      dest="filebase",
                      default="simulation")

    (opts,args) = parser.parse_args()

    gatherTimes(opts.size, opts.filebase)

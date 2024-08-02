#!/usr/bin/env python
import numpy
import os
from xml.etree.ElementTree import *
from proteus import Comm
comm = Comm.init()
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

def clearh5(filename,dataDir='.',addname="_clean", tCount=None, global_sync=True):
    """
    in case archiving failed to collect results from all times
    """
    import h5py
    xmlFile = open(filename+".xmf","r")
    h5File = h5py.File(filename+".h5","r")
    newh5File = h5py.File(filename+"_clean"+".h5", 'w')

    tree = ElementTree(file=xmlFile)
    xmlFile.close()
    XDMF = tree.getroot()
    Domain = XDMF[0]
    for TemporalGridCollection in Domain:
        for gridChild in TemporalGridCollection.findall("Grid"):
            counter=int(gridChild.find('Time').attrib['Name'])
            if(counter >= tCount):
              print(counter)
              Domain[0].remove(gridChild)

    import re
    for fieldNames in h5File.keys():
      result = int(re.search(r'\d+', fieldNames).group())
      if(result<tCount):
        newh5File.copy(h5File["/"+fieldNames],"/"+fieldNames)
        
    xmlFile = open(filename+addname+".xmf","w")
    indentXML(tree.getroot())
    tree.write(xmlFile)
    xmlFile.close()
    h5File.close()
    newh5File.close()

if __name__ == '__main__':
    from optparse import OptionParser
    usage = ""
    parser = OptionParser(usage=usage)
    parser.add_option("-f","--filebase",
                      help="base name for storage files",
                      action="store",
                      type="string",
                      dest="filebase",
                      default="simulation")

    parser.add_option("-t","--tCount",
                      help="number of time steps",
                      action="store",
                      type="int",
                      dest="tCount",
                      default="-1")
    parser.add_option("-n","--no-global-sync",
                      help="assume hdf5 has partitioned mesh",
                      action="store_true",
                      dest="not_global_sync",
                      default=False)

    (opts,args) = parser.parse_args()

    clearh5(opts.filebase,tCount = opts.tCount, global_sync = not opts.not_global_sync)

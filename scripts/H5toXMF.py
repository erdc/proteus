

#import numpy
#import os
#from xml.etree.ElementTree import *
import tables
#from Xdmf import *

def H5toXMF(basename,size,start,finaltime,stride):

# Open XMF files

    for step in range(start,finaltime+1,stride):
        XMFfile = open(basename+"."+str(step)+".xmf","w")
        XMFfile.write(r"""<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">
<Domain>"""+"\n")

        XMFfile.write(r'      <Grid GridType="Collection" CollectionType="Spatial">'+"\n")


        for proc in range(0,size):

            filename="solution.p"+str(proc)+"."+str(step)+".h5"
            print filename
            f1 = tables.open_file(filename)
            XMFfile.write (r'<Grid GridType="Uniform">'+"\n")
            XMFfile.write(r'      <Time Value="'+str(step)+'" />'+"\n")
            for tmp in f1.root:
                if tmp.name ==  "elements":
                    XMFfile.write (r'<Topology NumberOfElements="' +str(len(tmp[:]))+ '" Type="Tetrahedron">'+"\n")
                    XMFfile.write (r' <DataItem DataType="Int" Dimensions="' +str(len(tmp[:]))+ ' 4" Format="HDF">' + filename + ':/elements</DataItem>'+"\n")
                    XMFfile.write (r'</Topology>'+"\n")

                if tmp.name ==  "nodes":
                    XMFfile.write (r'<Geometry Type="XYZ">'+"\n")
                    XMFfile.write (r'  <DataItem DataType="Float" Dimensions="' +str(len(tmp[:]))+ ' 3" Format="HDF" Precision="8">' + filename + ':/nodes</DataItem>'+"\n")
                    XMFfile.write (r'</Geometry>'+"\n")

                if tmp.name ==  "u":
                    XMFfile.write (r'<Attribute AttributeType="Scalar" Center="Node" Name="u">'+"\n")
                    XMFfile.write (r'  <DataItem DataType="Float" Dimensions="' +str(len(tmp[:]))+ '" Format="HDF" Precision="8">' + filename + ':/u</DataItem>'+"\n")
                    XMFfile.write (r'</Attribute>'+"\n")

                if tmp.name ==  "v":
                    XMFfile.write (r'<Attribute AttributeType="Scalar" Center="Node" Name="v">'+"\n")
                    XMFfile.write (r'  <DataItem DataType="Float" Dimensions="' +str(len(tmp[:]))+ '" Format="HDF" Precision="8">' + filename + ':/v</DataItem>'+"\n")
                    XMFfile.write (r'</Attribute>'+"\n")

                if tmp.name ==  "w":
                    XMFfile.write (r'<Attribute AttributeType="Scalar" Center="Node" Name="w">'+"\n")
                    XMFfile.write (r'  <DataItem DataType="Float" Dimensions="' +str(len(tmp[:]))+ '" Format="HDF" Precision="8">' + filename + ':/w</DataItem>'+"\n")
                    XMFfile.write (r'</Attribute>'+"\n")

                if tmp.name ==  "p":
                    XMFfile.write (r'<Attribute AttributeType="Scalar" Center="Node" Name="p">'+"\n")
                    XMFfile.write (r'  <DataItem DataType="Float" Dimensions="' +str(len(tmp[:]))+ '" Format="HDF" Precision="8">' + filename + ':/p</DataItem>'+"\n")
                    XMFfile.write (r'</Attribute>'+"\n")

                if tmp.name ==  "phid":
                    XMFfile.write (r'<Attribute AttributeType="Scalar" Center="Node" Name="phid">'+"\n")
                    XMFfile.write (r'  <DataItem DataType="Float" Dimensions="' +str(len(tmp[:]))+ '" Format="HDF" Precision="8">' + filename + ':/phid</DataItem>'+"\n")
                    XMFfile.write (r'</Attribute>'+"\n")

            f1.close()
            XMFfile.write('      </Grid>'+"\n")

        XMFfile.write('    </Grid>'+"\n")
        XMFfile.write('   </Domain>'+"\n")
        XMFfile.write(' </Xdmf>'+"\n")
        XMFfile.close()


if __name__ == '__main__':
    from optparse import OptionParser
    usage = ""
    parser = OptionParser(usage=usage)

    parser.add_option("-n","--size",
                      help="number of processors for run",
                      action="store",
                      type="int",
                      dest="size",
                      default=1)

    parser.add_option("-s","--stride",
                      help="stride for solution output",
                      action="store",
                      type="int",
                      dest="stride",
                      default=0)

    parser.add_option("-t","--finaltime",
                      help="finaltime",
                      action="store",
                      type="int",
                      dest="finaltime",
                      default=1000)

    parser.add_option("-f","--filebase_flow",
                      help="base name for storage files",
                      action="store",
                      type="string",
                      dest="filebase",
                      default="solution")


    (opts,args) = parser.parse_args()

    start = 0
    if opts.stride == 0 :
        start = opts.finaltime
        opts.stride = 1

    H5toXMF(opts.filebase,opts.size,start,opts.finaltime,opts.stride)

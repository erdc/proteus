

#import numpy
#import os
#from xml.etree.ElementTree import *
import tables
#from Xdmf import *

def splitH5(basename1,basename2,size,start,finaltime,stride):

# Open XMF files
     xmfFiles = {}
     for step in range(start,finaltime+1,stride):
                XMFfile = open("solution."+str(step)+".xmf","w")
                XMFfile.write(r"""<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">
  <Domain>"""+"\n")
                XMFfile.write(r'      <Grid GridType="Collection"   CollectionType="Temporal">'+"\n")
		xmfFiles[step] = XMFfile
		


     for proc in range(0,size):

# Open h5 files
     	     filename=basename1+str(proc)+".h5"
	     print filename
     	     f1 = tables.openFile(filename)  
     	     
     	     filename=basename2+str(proc)+".h5"
	     print filename
     	     f2 = tables.openFile(filename)  
     	     
	     hdfFiles = {}
     	     for step in range(start,finaltime+1,stride):	     
     	    	filename="solution.p"+str(proc)+"."+str(step)+".h5"
     	     	hdfFiles[step]=  tables.openFile(filename,
				    mode = "w",
				    title = filename+" Data")

# Loop over entries and put in appropriate file
     	     for tmp in f1.root:
	         for step in range(start,finaltime+1,stride):
     		     if tmp.name ==  "elementsSpatial_Domain"+str(step):
     			     print tmp.name
     			     hdfFiles[step].createArray("/","elements",tmp[:])
			     
			     XMFfile=xmfFiles[step]
                  	     XMFfile.write (r'<Topology NumberOfElements="' +str(len(tmp[:]))+ '" Type="Tetrahedron">'+"\n")
         		     XMFfile.write (r' <DataItem DataType="Int" Dimensions="' +str(len(tmp[:]))+ ' 4" Format="HDF">' + filename + ':/elements</DataItem>'+"\n")
        		     XMFfile.write (r'</Topology>'+"\n")  
			        			     
		     if tmp.name ==  "nodesSpatial_Domain"+str(step):
     			     print tmp.name
     			     hdfFiles[step].createArray("/","nodes",tmp[:]) 
			     
			     XMFfile=xmfFiles[step]			     
                   	     XMFfile.write (r'<Geometry Type="XYZ">'+"\n")
                  	     XMFfile.write (r'  <DataItem DataType="Float" Dimensions="' +str(len(tmp[:]))+ ' 3" Format="HDF" Precision="8">' + filename + ':/nodes</DataItem>'+"\n")
                             XMFfile.write (r'</Geometry>'+"\n")
 
     		     if tmp.name ==  "u"+str(step):
     			     print tmp.name
     			     hdfFiles[step].createArray("/","u",tmp[:])
			     
			     XMFfile=xmfFiles[step]     					
                  	     XMFfile.write (r'<Attribute AttributeType="Scalar" Center="Node" Name="u">'+"\n")
                  	     XMFfile.write (r'  <DataItem DataType="Float" Dimensions="' +str(len(tmp[:]))+ '" Format="HDF" Precision="8">' + filename + ':/u</DataItem>'+"\n")
                  	     XMFfile.write (r'</Attribute>'+"\n")

     		     if tmp.name ==  "v"+str(step):
     			     print tmp.name
     			     hdfFiles[step].createArray("/","v",tmp[:])
			     
			     XMFfile=xmfFiles[step]			     
                  	     XMFfile.write (r'<Attribute AttributeType="Scalar" Center="Node" Name="v">'+"\n")
                  	     XMFfile.write (r'  <DataItem DataType="Float" Dimensions="' +str(len(tmp[:]))+ '" Format="HDF" Precision="8">' + filename + ':/v</DataItem>'+"\n")
                  	     XMFfile.write (r'</Attribute>'+"\n")
			     
     		     if tmp.name ==  "w"+str(step):
     			     print tmp.name
     			     hdfFiles[step].createArray("/","w",tmp[:])
			     
			     XMFfile=xmfFiles[step]			     
                  	     XMFfile.write (r'<Attribute AttributeType="Scalar" Center="Node" Name="w">'+"\n")
                  	     XMFfile.write (r'  <DataItem DataType="Float" Dimensions="' +str(len(tmp[:]))+ '" Format="HDF" Precision="8">' + filename + ':/w</DataItem>'+"\n")
                  	     XMFfile.write (r'</Attribute>'+"\n")
     			     
     		     if tmp.name ==  "p"+str(step):
     			     print tmp.name
     			     hdfFiles[step].createArray("/","p",tmp[:])
			     
			     XMFfile=xmfFiles[step]			     	
                  	     XMFfile.write (r'<Attribute AttributeType="Scalar" Center="Node" Name="p">'+"\n")
                  	     XMFfile.write (r'  <DataItem DataType="Float" Dimensions="' +str(len(tmp[:]))+ '" Format="HDF" Precision="8">' + filename + ':/p</DataItem>'+"\n")
                  	     XMFfile.write (r'</Attribute>'+"\n")
					
					

             f1.close()			     
			     
       	     filename=basename2+str(proc)+".h5"
	     print filename
     	     f2 = tables.openFile(filename)     			     

     	     for tmp in f2.root:
	        for step in range(start,finaltime+1,stride):
     		     if tmp.name ==  "phid"+str(step):
     			     print tmp.name
     			     hdfFiles[step].createArray("/","phid",tmp[:])
			    
			     XMFfile=xmfFiles[step]
                  	     XMFfile.write (r'<Attribute AttributeType="Scalar" Center="Node" Name="phid">'+"\n")
                             XMFfile.write (r'  <DataItem DataType="Float" Dimensions="' +str(len(tmp[:]))+ '" Format="HDF" Precision="8">' + filename + ':/phid</DataItem>'+"\n")
                  	     XMFfile.write (r'</Attribute>'+"\n")
			     
	     f2.close()
	     
	     for step in range(start,finaltime+1,stride):
	     	xmfFiles[step]	.write('      </Grid>'+"\n")
			     	
     for step in range(start,finaltime+1,stride):
		XMFfile=xmfFiles[step]	    
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
                      dest="filebase1",
                      default="twp_navier_stokes_wigley_3d_p")
		      
    parser.add_option("-p","--filebase_phi",
                      help="base name for storage files",
                      action="store",
                      type="string",
                      dest="filebase2",
                      default="redist_wigley_3d_p")
    
    (opts,args) = parser.parse_args()
    
    start = 0 
    if opts.stride == 0 :
   	start = opts.finaltime
	opts.stride = 1
    
    splitH5(opts.filebase1,opts.filebase2,opts.size,start,opts.finaltime,opts.stride)

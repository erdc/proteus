

#import numpy
#import os
#from xml.etree.ElementTree import *
import tables
#from Xdmf import *

def splitH5all(basename1,basename2,size,start,finaltime,stride):

     for proc in range(0,size):
     	splitH5single(basename1,basename2,proc,start,finaltime,stride)

     H5toXMF("solution",size,start,finaltime,stride)
	

def splitH5single(basename1,basename2,proc,start,finaltime,stride):



# Loop over entries and put in appropriate file
     	     filename=basename1+str(proc)+".h5"
	     print filename
     	     f1 = tables.openFile(filename)  
       	     filename=basename2+str(proc)+".h5"
	     print filename
     	     f2 = tables.openFile(filename) 
	     
     	     for step in range(start,finaltime+1,stride):

     	    	filename="sol.p"+str(proc)+"."+str(step)+".h5"
     	     	hdfFile=  tables.openFile(filename,
				    mode = "w",
				    title = filename+" Data")
				    
	     
     		name =  "elementsSpatial_Domain"+str(step)
     		hdfFile.createArray("/","elements",f1.getNode("/",name)[:])
		     
		name =  "nodesSpatial_Domain"+str(step)
     		hdfFile.createArray("/","nodes",f1.getNode("/",name)[:]) 
		    
     		name =  "u"+str(step)
     		hdfFile.createArray("/","u",f1.getNode("/",name)[:])

     		name =  "v"+str(step)
     		hdfFile.createArray("/","v",f1.getNode("/",name)[:])
		     
     		name =  "w"+str(step)
     		hdfFile.createArray("/","w",f1.getNode("/",name)[:])
		     
     		name =  "p"+str(step)
     		hdfFile.createArray("/","p",f1.getNode("/",name)[:])
		
     		name =  "phid"+str(step)
     		hdfFile.createArray("/","phid",f2.getNode("/",name)[:])					

                hdfFile.close()

             f1.close()			    			     
             f2.close()
		     	    
def H5toXMF(basename,size,start,finaltime,stride):

# Open XMF files

     for step in range(start,finaltime+1,stride):
                filename = basename+"."+str(step)+".h5",
     	     	hdfFile=  tables.openFile(filename,
				    mode = "w",
				    title = filename+" Data")     
     
                XMFfile = open(basename+"."+str(step)+".xmf","w")
                XMFfile.write(r"""<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">
  <Domain>"""+"\n")
	
		XMFfile.write(r'      <Grid GridType="Collection" CollectionType="Spatial">'+"\n")	


                for proc in range(0,size):

 
     	        	f1 = tables.openFile("sol.p"+str(proc)+"."+str(step)+".h5")  
			
		        XMFfile.write (r'<Grid GridType="Uniform">'+"\n")
			XMFfile.write(r'      <Time Value="'+str(step)+'" />'+"\n")

     			data=f1.getNode("/","elements")[:]
     			hdfFile.createArray("/","elements_p"+str(proc),data)

                  	XMFfile.write (r'<Topology NumberOfElements="' +str(len(data))+ '" Type="Tetrahedron">'+"\n")
         		XMFfile.write (r' <DataItem DataType="Int" Dimensions="' +str(len(data))+ ' 4" Format="HDF">' + filename + ':/element_p'+str(proc)'</DataItem>'+"\n")
        		XMFfile.write (r'</Topology>'+"\n")  

     			data=f1.getNode("/","nodes")[:]
     			hdfFile.createArray("/","nodes_p"+str(proc),data)
					     
                   	XMFfile.write (r'<Geometry Type="XYZ">'+"\n")
                  	XMFfile.write (r'  <DataItem DataType="Float" Dimensions="' +str(len(data))+ ' 3" Format="HDF" Precision="8">' + filename + ':/nodes_p'+str(proc)'</DataItem>'+"\n")
                        XMFfile.write (r'</Geometry>'+"\n")

     			data=f1.getNode("/","u")[:]
     			hdfFile.createArray("/","u_p"+str(proc),data)
			  			
                  	XMFfile.write (r'<Attribute AttributeType="Scalar" Center="Node" Name="u">'+"\n")
                  	XMFfile.write (r'  <DataItem DataType="Float" Dimensions="' +str(len(data))+ '" Format="HDF" Precision="8">' + filename + ':/u_p'+str(proc)'</DataItem>'+"\n")
                  	XMFfile.write (r'</Attribute>'+"\n")

     			data=f1.getNode("/","v")[:]
     			hdfFile.createArray("/","v_p"+str(proc),data)
			
     			XMFfile.write (r'<Attribute AttributeType="Scalar" Center="Node" Name="v">'+"\n")
                  	XMFfile.write (r'  <DataItem DataType="Float" Dimensions="' +str(len(data))+ '" Format="HDF" Precision="8">' + filename + ':/v_'+str(proc)'</DataItem>'+"\n")
                  	XMFfile.write (r'</Attribute>'+"\n")

     			data=f1.getNode("/","w")[:]
     			hdfFile.createArray("/","w_p"+str(proc),data)
						
                  	XMFfile.write (r'<Attribute AttributeType="Scalar" Center="Node" Name="w">'+"\n")
                  	XMFfile.write (r'  <DataItem DataType="Float" Dimensions="' +str(len(data))+ '" Format="HDF" Precision="8">' + filename + ':/w_p'+str(proc)'</DataItem>'+"\n")
                  	XMFfile.write (r'</Attribute>'+"\n")

     			data=f1.getNode("/","p")[:]
     			hdfFile.createArray("/","p_p"+str(proc),data)
						
                  	XMFfile.write (r'<Attribute AttributeType="Scalar" Center="Node" Name="p">'+"\n")
                  	XMFfile.write (r'  <DataItem DataType="Float" Dimensions="' +str(len(data))+ '" Format="HDF" Precision="8">' + filename + ':/p_p'+str(proc)'</DataItem>'+"\n")
                  	XMFfile.write (r'</Attribute>'+"\n")

     			data=f1.getNode("/","phid")[:]
     			hdfFile.createArray("/","phid_p"+str(proc),data)
			
                  	XMFfile.write (r'<Attribute AttributeType="Scalar" Center="Node" Name="phid">'+"\n")
                        XMFfile.write (r'  <DataItem DataType="Float" Dimensions="' +str(len(data))+ '" Format="HDF" Precision="8">' + filename + ':/phid_p'+str(proc)'</DataItem>'+"\n")
                  	XMFfile.write (r'</Attribute>'+"\n")
			     
	     		f1.close()
                        XMFfile.write('      </Grid>'+"\n")

        	XMFfile.write('    </Grid>'+"\n")
                XMFfile.write('   </Domain>'+"\n")
                XMFfile.write(' </Xdmf>'+"\n")		
		XMFfile.close()
		hdfFile.close()


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
    if (opts.size >0) :
    	splitH5all(opts.filebase1,opts.filebase2,opts.size,start,opts.finaltime,opts.stride)
    else :
        splitH5single(opts.filebase1,opts.filebase2,-opts.size,start,opts.finaltime,opts.stride)	

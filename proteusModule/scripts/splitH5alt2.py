

#import numpy
#import os
#from xml.etree.ElementTree import *
import tables
#from Xdmf import *

def splitH5all(basename1,basename2,size,start,finaltime,stride):

     for proc in range(0,size):
     	splitH5single(basename1,basename2,proc,start,finaltime,stride)



def splitH5single(basename1,basename2,proc,start,finaltime,stride):

# Open output h5 files    	     
	     hdfFiles = {}
     	     for step in range(start,finaltime+1,stride):	     
     	    	filename="solution.p"+str(proc)+"."+str(step)+".h5"
     	     	hdfFiles[step]=  tables.openFile(filename,
				    mode = "w",
				    title = filename+" Data")

# Loop over entries and put in appropriate file
     	     filename=basename1+str(proc)+".h5"
	     print filename
     	     f1 = tables.openFile(filename)  

     	     for step in range(start,finaltime+1,stride):
     		     name =  "elementsSpatial_Domain"+str(step)
     		     hdfFiles[step].createArray("/","elements",f1.getNode("/",name)[:])
			     
		     name =  "nodesSpatial_Domain"+str(step
     		     hdfFiles[step].createArray("/","nodes",f1.getNode("/",name)[:]) 
			    
     		     name =  "u"+str(step)
     		     hdfFiles[step].createArray("/","u",f1.getNode("/",name)[:])

     		     name =  "v"+str(step)
     		     hdfFiles[step].createArray("/","v",f1.getNode("/",name)[:])
			     
     		     name =  "w"+str(step)
     		     hdfFiles[step].createArray("/","w",f1.getNode("/",name)[:])
			     
     		     name =  "p"+str(step)
     		     hdfFiles[step].createArray("/","p",f1.getNode("/",name)[:])
					

             f1.close()			     
	     
# Loop over entries and put in appropriate file			     
       	     filename=basename2+str(proc)+".h5"
	     print filename
     	     f2 = tables.openFile(filename)     			     

     	     for step in range(start,finaltime+1,stride):
     		name =  "phid"+str(step):
     		hdfFiles[step].createArray("/","phid",f2.getNode("/",name),[:])
			     
             f2.close()
		     
# Close solution files
	     for step in range(start,finaltime+1,stride):	     
     	    	hdfFiles[step].close()
	    

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

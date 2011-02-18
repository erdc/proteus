

#import numpy
#import os
#from xml.etree.ElementTree import *
import tables
#from Xdmf import *

def splitH5all(basename1,basename2,size,start,finaltime,stride):

     for proc in range(0,size):
     	splitH5single(basename1,basename2,proc,start,finaltime,stride)



def splitH5single(basename1,basename2,proc,start,finaltime,stride):


# Loop over entries and put in appropriate file
     	     filename=basename1+str(proc)+".h5"
	     print filename
     	     f1 = tables.openFile(filename)  
       	     filename=basename2+str(proc)+".h5"
	     print filename
     	     f2 = tables.openFile(filename) 
	     
     	     for step in range(start,finaltime+1,stride):

     	    	filename="solution.p"+str(proc)+"."+str(step)+".h5"
     	     	hdfFile=  tables.openFile(filename,
				    mode = "w",
				    title = filename+" Data")
				    
	     
     		name =  "elementsSpatial_Domain"+str(step)
     		hdfFile.createArray("/","elements",f1.getNode("/",name)[:])
		     
		name =  "nodesSpatial_Domain"+str(step
     		hdfFile.createArray("/","nodes",f1.getNode("/",name)[:]) 
		    
     		name =  "u"+str(step)
     		hdfFile.createArray("/","u",f1.getNode("/",name)[:])

     		name =  "v"+str(step)
     		hdfFile.createArray("/","v",f1.getNode("/",name)[:])
		     
     		name =  "w"+str(step)
     		hdfFile.createArray("/","w",f1.getNode("/",name)[:])
		     
     		name =  "p"+str(step)
     		hdfFile.createArray("/","p",f1.getNode("/",name)[:])
		
     		name =  "phid"+str(step):
     		hdfFile.createArray("/","phid",f2.getNode("/",name),[:])					

                hdfFile.close()

             f1.close()			    			     
             f2.close()
		     	    

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

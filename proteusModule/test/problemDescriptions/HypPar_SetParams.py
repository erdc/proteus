#John Chrispell, Spring 08
from math import *
L=(2.0,1.0,1.0)

T = 1.0
nd = 1

alpha = 1.0
delta =  100*((2.0/128)**2.0)

q=[1.0]   #[-1.0]


plotAlpha = False

class SimpleTest:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        if((x[0] <= 1.0) and (x[0] >= 0.0)): 
	    return -1.0*x[0]*x[0] + 2.0*x[0] + x[0]*x[0]*t - 2.0*x[0]*t
        if((x[0] <= 1.5) and (x[0] > 1.0)):
	    return -1.0*x[0]*x[0] + 2.0*x[0] + x[0]*x[0]*t - 2.0*x[0]*t
        if((x[0] <= 2.0) and (x[0] > 1.5)):
	    return -1.0*x[0]*x[0] + 2.0*x[0] + x[0]*x[0]*t - 2.0*x[0]*t
 
    def duOfXT(self,x,t):
        if((x[0] <= 1.0) and (x[0] >= 0.0)): 
	    return -2.0*x[0] + 2.0 + 2.0*x[0]*t - 2.0*t
        if((x[0] <= 1.5) and (x[0] > 1.0)):
	    return -2.0*x[0] + 2.0 + 2.0*x[0]*t - 2.0*t
        if((x[0] <= 2.0) and (x[0] > 1.5)):
	    return -2.0*x[0] + 2.0 + 2.0*x[0]*t - 2.0*t
    
    def rOfXT(self,x,t):
        #print "x = ", x
        if((x[0] <= 1.0) and (x[0] >= 0.0)): 
	    q = 1.0
	    alpha = 1.0
	    #print "Calling rOfXT int one"
	    return x[0]*x[0] - 2.0*x[0] + q*(-2.0*x[0] + 2.0 + 2.0*x[0]*t - 2.0*t) - alpha*(-2.0 + 2.0*t)
        
	if((x[0] <= 1.5) and (x[0] > 1.0)):
	    q = 1.0
	    alpha = 1.0
	    #print "Calling rOfXT int two"
	    return x[0]*x[0] - 2.0*x[0] + q*(-2.0*x[0] + 2.0 + 2.0*x[0]*t - 2.0*t) - alpha*(-2.0 + 2.0*t)
	
	if((x[0] <= 2.0) and (x[0] > 1.5)):
	    q = 1.0
	    alpha = 1.0 
	    #print "Calling rOfXT int three"
	    return x[0]*x[0] - 2.0*x[0] + q*(-2.0*x[0] + 2.0 + 2.0*x[0]*t - 2.0*t) - alpha*(-2.0 + 2.0*t)
    def drOfXT(self,x,t):
        if((x[0] <= 1.0) and (x[0] >= 0.0)):
	    q = 1.0
	    alpha = 1.0 
	    return 2.0*x[0] - 2.0 + q*(-2.0 + 2.0*t) 
        
	if((x[0] <= 1.5) and (x[0] > 1.0)):
	    q=1.0
	    alpha = 1.0
	    return 2.0*x[0] - 2.0 + q*(-2.0 + 2.0*t) 
	    
        if((x[0] <= 2.0) and (x[0] > 1.5)):
	    q = 1.0
	    alpha = 1.0
	    return 2.0*x[0] - 2.0 + q*(-2.0 + 2.0*t) 
	    
	    
class Ern_and_Proft_Test1:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        if((x[0] <= 1.0) and (x[0] >= 0.0)): 
	    return -sin(x[0] - 1) + exp(x[0] - 1.0) + t*((x[0] -1.0)**2.0) 
        if((x[0] <= 1.5) and (x[0] > 1.0)):
	    return ((x[0]-1.0)**2.0)*(1.0+t) + 1.0
        if((x[0] <= 2.0) and (x[0] > 1.5)):
	    return 1.0 + x[0]**3.0 + t*(1.0/4.0)*(x[0]**2.0) - 5.0*x[0]*x[0] + 8.0*x[0] - x[0]*t - (33.0/8.0) + (15.0/16.0)*t
 
    def duOfXT(self,x,t):
        if((x[0] <= 1.0) and (x[0] >= 0.0)): 
	    return -cos(x[0] - 1) + exp(x[0] - 1.0) + 2.0*t*(x[0] -1.0)
        if((x[0] <= 1.5) and (x[0] > 1.0)):
	    return 2*x[0] - 2.0 + 2.0*t*x[0] - 2.0*t
        if((x[0] <= 2.0) and (x[0] > 1.5)):
	    return 3.0*x[0]*x[0] + (1.0/2.0)*x[0]*t - 10.0*x[0] + 8.0 - t
    
    def rOfXT(self,x,t):
        if((x[0] <= 1.0) and (x[0] >= 0.0)): 
	    q = 1.0
	    alpha = 1.0
	    return (x[0]-1.0)**2.0 + q*( -1.0*cos(x[0]-1.0) + exp(x[0] - 1.0) + t*(x[0] - 1.0)*2.0) - alpha*( sin(x[0]-1.0) + exp(x[0] - 1.0) + 2*t)
        
	if((x[0] <= 1.5) and (x[0] > 1.0)):
	    q = 1.0
	    alpha = 0.0
	    return ((x[0]-1.0)**2.0) + q*(2.0*x[0] - 2.0 + 2.0*t*x[0] - 2.0*t) - alpha*(2.0 + 2.0*t)
	
	if((x[0] <= 2.0) and (x[0] > 1.5)):
	    q = 1.0
	    alpha = 1.0 
	    return ((1.0/4.0)*x[0]*x[0] - x[0] + (15.0/16.0)) + q*(3.0*x[0]*x[0] + (1.0/2.0)*x[0]*t -10.0*x[0] + 8.0 - t) - alpha*(6.0*x[0] + (1.0/2.0)*t -10.0)
 
    def drOfXT(self,x,t):
        if((x[0] <= 1.0) and (x[0] >= 0.0)):
	    q = 1.0
	    alpha = 1.0 
	    return 2.0*(x[0]-1) + q*(sin(x[0]-1.0) + exp(x[0]-1.0) + 2.0*t) - alpha*(cos(x[0]-1.0) + exp(x[0]-1.0))
        
	if((x[0] <= 1.5) and (x[0] > 1.0)):
	    q=1.0
	    alpha = 0.0
	    return 2.0*(x[0]-1.0) + 2.0*q*(1.0+t)
	    
        if((x[0] <= 2.0) and (x[0] > 1.5)):
	    q = 1.0
	    alpha = 1.0
	    return 0.5*x[0] - 1.0 + q*(6.0*x[0] + (1.0/2.0)*t - 10.0) - alpha*6.0
	    
rFuncClass = Ern_and_Proft_Test1
#rFuncClass = SimpleTest
#rFuncClass = None

def setParams(x_in,alpha_in):

    alpha_in.flat[:] = alpha
    
    value_lhs_a = 1.0
    value_in_ab = 0.0
    value_rhs_b = 1.0
    
    domain_split_a = 1.0
    domain_split_b = 1.5

    #have to do some mesh dependent stuff here
    if len(x_in.shape) == 3: #on element quadrature
        for eN in range(x_in.shape[0]):
            for k in range(x_in.shape[1]):
                if ((x_in[eN,k,0] <= domain_split_a) and (x_in[eN,k,0] >= 0.0)):
                    alpha_in[eN,k] = value_lhs_a
		if ((x_in[eN,k,0] > domain_split_a) and (x_in[eN,k,0] <= domain_split_b)):
		    alpha_in[eN,k] = value_in_ab
		if ((x_in[eN,k,0] > domain_split_b) and (x_in[eN,k,0] <= L[0])):
		    alpha_in[eN,k] = value_rhs_b
		    
        if plotAlpha:
            from pyadh import Viewers #need to check that it's gnuplot
            dgridx=32; dgridy=32; dgridp=16;
            for eN in range(x_in.shape[0]):
                for k in range(x_in.shape[1]):
                    Viewers.datFile.write("%12.5e %12.5e %12.5e \n" % (x_in[eN,k,0],x_in[eN,k,1],alpha_in[eN,k]))
            Viewers.datFile.write("\n \n")
            cmd = "set dgrid3d %d,%d,%g; set contour base; set term x11 %i; splot \'%s\' index %i with lines title \"%s\" \n" % (dgridx,dgridy,dgridp,
                                                                                                                                 Viewers.windowNumber,
                                                                                                                                 Viewers.datFilename,
                                                                                                                                 Viewers.plotNumber,
                                                                                                                                 'alpha')
            Viewers.cmdFile.write(cmd)
            Viewers.viewerPipe.write(cmd)
            Viewers.plotNumber+=1
            Viewers.windowNumber += 1
            
            raw_input('press return to continue')
            Viewers.windowNumber -= 1
    elif len(x_in.shape) == 4: #on element boundary quadrature
        for eN in range(x_in.shape[0]):
            for ebN in range(x_in.shape[1]):
                for k in range(x_in.shape[2]):
                    if ((x_in[eN,k,0] <= domain_split_a) and (x_in[eN,k,0] >= 0.0)):
                        alpha_in[eN,ebN,k] = value_lhs_a
		    if ((x_in[eN,k,0] > domain_split_a) and (x_in[eN,k,0] <= domain_split_b)): 
		        alpha_in[eN,ebN,k] = value_in_ab
                    if ((x_in[eN,k,0] > domain_split_b) and (x_in[eN,k,0] <= L[0])):
		        alpha_in[eN,ebN,k] = value_rhs_b
    #print "SetParams alpha_in", alpha_in 


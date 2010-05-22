from pyadh import *
from pyadh.default_p import *

nd=2

rhow = 997.0        # density wetting      (Water 997 kg/m^3)
rhon = 1.205        # density nonwetting   (Air 1.205 kg/m^3)
muw  = 1.002e-3     # viscosity nonwetting (Water := 1.002e-3 kg/m s)
mun  = 1.81e-5      # viscosity wetting    (Air := 1.81e-5 kg/ m s)
g = [0.0,-gravity,0.0]


polyfile = "forsyth2"
top = 6.5 #m
right = 8 #m
L=(right,top,1.0)

analyticalSolution = None

# # # # # # # # media types # # # # # # # # # # 
#     Type 1: k_i     = 9.33e-12 m^2 
#             K_s     = 7.89     m/d  (k_i*998.2*7.3209e10/86.4)
#             \theta_s= 0.3680
#             \theta_r= 0.3680*0.2771
#             \alpha  = 3.334    1/m
#             n       = 1.982
#     Type 2: k_i     = 5.55e-12 m^2
#             K_s     = 4.69     m/d
#             \theta_s= 0.3510
#             \theta_r= 0.3510*0.2806
#             \alpha  = 3.63     1/m
#             n       = 1.632
#     Type 3: k_i     = 4.898e-12 m^2
#             K_s     = 4.143
#             \theta_s= 0.3250
#             \theta_r= 0.3250*0.2643
#             \alpha  = 3.455     1/m
#             n       = 5.0
#     Type 4: k_i     = 4.898e-11 m^2
#             K_s     = 41.43     m/d
#             \theta_s= 0.3250
#             \theta_r= 0.3250*0.2643
#             \alpha  = 3.455     1/m
#             n       = 5.0
     

nMediaTypes  = 4
alphaVGtypes = Numeric.zeros((nMediaTypes,),Numeric.Float)
nVGtypes     = Numeric.zeros((nMediaTypes,),Numeric.Float)
thetaStypes  = Numeric.zeros((nMediaTypes,),Numeric.Float)
thetaRtypes  = Numeric.zeros((nMediaTypes,),Numeric.Float)
KsTypes      = Numeric.zeros((nMediaTypes,),Numeric.Float)

#Type ids are base 1

#Type 1
it = 0;
alphaVGtypes[it] = 3.334 ; nVGtypes[it]    = 1.982;
thetaStypes[it]  = 0.3680; thetaRtypes[it] = 0.3680*0.2771;
KsTypes[it]      = 7.89;

#Type 2
it = 1;
alphaVGtypes[it] = 3.63 ;  nVGtypes[it]    = 1.632;
thetaStypes[it]  = 0.3510; thetaRtypes[it] = 0.3510*0.2806;
KsTypes[it]      = 4.69;

#Type 3
it = 2;
alphaVGtypes[it] = 3.455;  nVGtypes[it]    = 5.0;
thetaStypes[it]  = 0.3250; thetaRtypes[it] = 0.3250*0.2643;
KsTypes[it]      = 4.143;

#Type 4
it = 3;
alphaVGtypes[it] = 3.455;  nVGtypes[it]    = 5.0;
thetaStypes[it]  = 0.3250; thetaRtypes[it] = 0.3250*0.2643;
KsTypes[it]      = 41.43;


#set dimensions for types
viscosity     = 8.9e-4  #kg/(m*s)
density       = 998.2   #kg/m^3
gravity       = 9.8     #m/s^2
m_per_s_by_m_per_d = 1.1574074e-5
lengthScale   = 1.0     #m
timeScale     = 1.0     #d #1.0/sqrt(g*lengthScale)

# Set some initial parameters that may be overwritten

model = 'VGM'       # Choose from {simp,BCB,BCM,VGB,VGM}
q = 0.0
a = 0.0
Ksw = 0.0
mvg_alpha = 5.47    # Sand 5.47 
mvg_n = 4.264       # Van Genuchten model param (see above)
bc_lambda=mvg_n-1.0 # Lambda for Brooks-Corey (BC) and m Van Genuchten Models
bc_pd = 1.0/mvg_alpha
mvg_m = 1.0-1.0/mvg_n
omega = 0.3         # volumetric fraction of porous medium that can transmit fluid (0.01 - 0.9)

T= 30./timeScale


dimensionless_conductivity  = Numeric.zeros(KsTypes.shape,Numeric.Float)
dimensionless_alpha         = Numeric.zeros(alphaVGtypes.shape,Numeric.Float)
dimensionless_density  = 1.0
dimensionless_gravity  = Numeric.array([0.0,
                                        -1.0,
                                        0.0])
pondingPressure=0.1
initialPressure=-880.0e3 / density / gravity #convert kPa to m
rechargeRate =  -2.0e-2#- means into domain m/d
rechargeXboundary = 2.25# L[0]/3.0  where recharge stops
waterTableHeight = L[1]/3.0


for i in range(len(KsTypes.flat)):
    permeability = (KsTypes[i]*m_per_s_by_m_per_d)*viscosity/(gravity*density)  #m^2
    dimensionless_conductivity[i] = (timeScale*density*gravity*permeability/(viscosity*lengthScale))/m_per_s_by_m_per_d
    print 'Ks',dimensionless_conductivity
    dimensionless_alpha[i] = alphaVGtypes[i]*lengthScale
    


#set up heterogeneity blocks, zone ids are base 1
nHetBlocks = 7
heterogeneityBlocks = {}

#zone 1 is [0,8]x[6,6.5], type 1
zone = 0
heterogeneityBlocks[zone] = {};
heterogeneityBlocks[zone]['Xm']  = [0.0,6.0,0.0]; heterogeneityBlocks[zone]['Xp'] = [right,top,1.0]
heterogeneityBlocks[zone]['type']= 0
#zone 2 is [0,8]x[5.5,6.0], type 2
zone = 1
heterogeneityBlocks[zone] = {};
heterogeneityBlocks[zone]['Xm']  = [0.0,5.5,0.0]; heterogeneityBlocks[zone]['Xp'] = [right,6.0,1.0]
heterogeneityBlocks[zone]['type']= 1
#zone 3 is [0,8]x[5.0,5.5], type 3
zone = 2
heterogeneityBlocks[zone] = {};
heterogeneityBlocks[zone]['Xm']  = [0.0,5.0,0.0]; heterogeneityBlocks[zone]['Xp'] = [right,5.5,1.0]
heterogeneityBlocks[zone]['type']= 2
#zone 4 is [0,1]x[4.0,5.0], type 3
zone = 3
heterogeneityBlocks[zone] = {};
heterogeneityBlocks[zone]['Xm']  = [0.0,4.0,0.0]; heterogeneityBlocks[zone]['Xp'] = [1.0,5.0,1.0]
heterogeneityBlocks[zone]['type']= 2
#zone 5 is [1,3]x[4.0,5.0], type 4
zone = 4
heterogeneityBlocks[zone] = {};
heterogeneityBlocks[zone]['Xm']  = [1.0,4.0,0.0]; heterogeneityBlocks[zone]['Xp'] = [3.0,5.0,1.0]
heterogeneityBlocks[zone]['type']= 3
#zone 6 is [3,8]x[4.0,5.0], type 3
zone = 5
heterogeneityBlocks[zone] = {};
heterogeneityBlocks[zone]['Xm']  = [3.0,4.0,0.0]; heterogeneityBlocks[zone]['Xp'] = [right,5.0,1.0]
heterogeneityBlocks[zone]['type']= 2
#zone 7 is [0,8]x[0.0,4.0], type 3
zone = 6
heterogeneityBlocks[zone] = {};
heterogeneityBlocks[zone]['Xm']  = [0.0,0.0,0.0]; heterogeneityBlocks[zone]['Xp'] = [right,4.0,1.0]
heterogeneityBlocks[zone]['type']= 2


plotHet = False

bndEps = 1.0e-8

def setParams(x_in,bc_lambda_in,bc_pd_in,vgm_m_in,vgm_alpha_in,Ks_in,thetaR_in,thetaS_in):
    #brute force
    
    
    #vgm_n_in.flat[:]     = nVGtypes[0]
    vgm_m_in.flat[:]     = 1.0 - (1.0/nVGtypes[0])
    vgm_alpha_in.flat[:] = alphaVGtypes[0]
    bc_pd_in.flat[:]     = 1.0/alphaVGtypes[0]
    bc_lambda_in.flat[:] = nVGtypes[0] - 1.0
    Ks_in.flat[:]        = dimensionless_conductivity[0]
    thetaR_in.flat[:]    = thetaRtypes[0]
    thetaS_in.flat[:]    = thetaStypes[0]
    
    
    #have to do some mesh dependent stuff here
    if len(x_in.shape) == 3: #on element quadrature
        for eN in range(x_in.shape[0]):
            foundAblock=False
            for nB in range(nHetBlocks):
                inBlock = True
                for k in range(x_in.shape[1]):
                    if not (x_in[eN,k,0] >= heterogeneityBlocks[nB]['Xm'][0]-bndEps and
                            x_in[eN,k,0] <= heterogeneityBlocks[nB]['Xp'][0]+bndEps and
                            x_in[eN,k,1] >= heterogeneityBlocks[nB]['Xm'][1]-bndEps and
                            x_in[eN,k,1] <= heterogeneityBlocks[nB]['Xp'][1]+bndEps):
                        inBlock = False
                        break
                if inBlock:
                    foundAblock=True
                    #mwf debug
                    #if heterogeneityBlocks[nB]['type'] == 0:
                    #    print """re_forsyth2 set eN=%d nB=%d type=%d x=%s """ % (eN,nB,heterogeneityBlocks[nB]['type'],
                    #                                                             x_in[eN,:,:])
                    for k in range(x_in.shape[1]):
                        type = heterogeneityBlocks[nB]['type']
			
			
                        #vgm_n_in[eN,k]     = nVGtypes[type]
			vgm_m_in[eN,k]     = 1.0 - (1.0/nVGtypes[type])
                        vgm_alpha_in[eN,k] = alphaVGtypes[type]
			bc_pd_in[eN,k]     = 1.0/alphaVGtypes[type]
    			bc_lambda_in[eN,k] = nVGtypes[type] - 1.0
			
			
                        Ks_in[eN,k]        = dimensionless_conductivity[type]
                        thetaR_in[eN,k]    = thetaRtypes[type]
                        thetaS_in[eN,k]    = thetaStypes[type]
			
			
                        #mwf debug
                        if type == 0:
                            assert x_in[eN,k,1] >= 6.0, "type 0 y >= 6.0"
                        if vgm_alpha_in[eN,k] == 3.334:
                            assert type == 0, "alpha 3.334 type = %d " % type
                        if x_in[eN,k,0] < 1.0 and x_in[eN,k,1] < 5.5:
                            assert type == 2, "x= %s type= %d " % (x_in[eN,k,:],type)
                        if (x_in[eN,k,0] > 1.0 and x_in[eN,k,0] < 3.0
                            and x_in[eN,k,1] > 4.0 and x_in[eN,k,1] < 5.0):
                            assert type == 3, "x= %s type= %d " % (x_in[eN,k,:],type)
                #in block
            #nB
            if not foundAblock:
                print """didn't find block eN=%d """ % eN
                for k in range(x_in.shape[1]):
                    print """x[%d,%d]=%s blockTests= """ % (eN,k,x_in[eN,k,0:2])
                    for nB in range(nHetBlocks):
                        inblock= (x_in[eN,k,0] >= heterogeneityBlocks[nB]['Xm'][0] and
                                x_in[eN,k,0] <= heterogeneityBlocks[nB]['Xp'][0] and
                                x_in[eN,k,1] >= heterogeneityBlocks[nB]['Xm'][1] and
                                x_in[eN,k,1] <= heterogeneityBlocks[nB]['Xp'][1])
                        print """BLOCK %d x >= %g x <=%g y >= %g y <= %g = %s """ % (nB,
                                                                                     heterogeneityBlocks[nB]['Xm'][0],
                                                                                     heterogeneityBlocks[nB]['Xp'][0],
                                                                                     heterogeneityBlocks[nB]['Xm'][1],
                                                                                     heterogeneityBlocks[nB]['Xp'][1],
                                                                                     inblock)
            assert foundAblock == True, "foundAblock eN=%d x=%s " % (eN,x_in[eN,:,0:2])
        if plotHet:
            from pyadh import Viewers #need to check that it's gnuplot
            dgridx=32; dgridy=32; dgridp=16;
            for eN in range(x_in.shape[0]):
                for k in range(x_in.shape[1]):
                    Viewers.datFile.write("%12.5e %12.5e %12.5e \n" % (x_in[eN,k,0],x_in[eN,k,1],Ks_in[eN,k]))
            Viewers.datFile.write("\n \n")
            cmd = "set dgrid3d %d,%d,%g; set contour base; set term x11 %i; splot \'%s\' index %i with lines title \"%s\" \n" % (dgridx,dgridy,dgridp,
                                                                                                                                 Viewers.windowNumber,
                                                                                                                                 Viewers.datFilename,
                                                                                                                                 Viewers.plotNumber,
                                                                                                                                 'Ks')
            Viewers.cmdFile.write(cmd)
            Viewers.viewerPipe.write(cmd)
            Viewers.plotNumber+=1
            Viewers.windowNumber += 1
            for eN in range(x_in.shape[0]):
                for k in range(x_in.shape[1]):
                    Viewers.datFile.write("%12.5e %12.5e %12.5e \n" % (x_in[eN,k,0],x_in[eN,k,1],vgm_m_in[eN,k]))
            Viewers.datFile.write("\n \n")
            cmd = "set dgrid3d %d,%d,%g; set contour base; set term x11 %i; splot \'%s\' index %i with lines title \"%s\" \n" % (dgridx,dgridy,dgridp,
                                                                                                                                 Viewers.windowNumber,
                                                                                                                                 Viewers.datFilename,
                                                                                                                                 Viewers.plotNumber,
                                                                                                                                 'vgm-m')
            Viewers.cmdFile.write(cmd)
            Viewers.viewerPipe.write(cmd)
            Viewers.plotNumber+=1
            Viewers.windowNumber += 1
            for eN in range(x_in.shape[0]):
                for k in range(x_in.shape[1]):
                    Viewers.datFile.write("%12.5e %12.5e %12.5e \n" % (x_in[eN,k,0],x_in[eN,k,1],vgm_alpha_in[eN,k]))
            Viewers.datFile.write("\n \n")
            cmd = "set dgrid3d %d,%d,%g; set contour base; set term x11 %i; splot \'%s\' index %i with lines title \"%s\" \n" % (dgridx,dgridy,dgridp,
                                                                                                                                 Viewers.windowNumber,
                                                                                                                                 Viewers.datFilename,
                                                                                                                                 Viewers.plotNumber,
                                                                                                                                 'vgm-alpha')
            Viewers.cmdFile.write(cmd)
            Viewers.viewerPipe.write(cmd)
            Viewers.plotNumber+=1
            Viewers.windowNumber += 1
            raw_input('press return to continue')
            Viewers.windowNumber -= 3
    #end element quad
    elif len(x_in.shape) == 4: #on element boundary quadrature
        for eN in range(x_in.shape[0]):
            foundAblock=False
            for nB in range(nHetBlocks):
                inBlock = True
                for ebN in range(x_in.shape[1]):
                    for k in range(x_in.shape[2]):
                        if not (x_in[eN,ebN,k,0] >= heterogeneityBlocks[nB]['Xm'][0]-bndEps and
                                x_in[eN,ebN,k,0] <= heterogeneityBlocks[nB]['Xp'][0]+bndEps and
                                x_in[eN,ebN,k,1] >= heterogeneityBlocks[nB]['Xm'][1]-bndEps and 
                                x_in[eN,ebN,k,1] <= heterogeneityBlocks[nB]['Xp'][1]+bndEps):
                            inBlock = False
                            break
                    if inBlock == False:
                        break
                if inBlock:
                    foundAblock=True
                    for ebN in range(x_in.shape[1]):
                        for k in range(x_in.shape[2]):
                            type = heterogeneityBlocks[nB]['type']
                            #vgm_n_in[eN,ebN,k]     = nVGtypes[type]
                            vgm_alpha_in[eN,ebN,k] = alphaVGtypes[type]			
			    vgm_m_in[eN,ebN,k]     = 1.0 - (1.0/nVGtypes[type])
			    bc_pd_in[eN,ebN,k]     = 1.0/alphaVGtypes[type]
    			    bc_lambda_in[eN,ebN,k] = nVGtypes[type] - 1.0
                            Ks_in[eN,ebN,k]        = dimensionless_conductivity[type]
                            thetaR_in[eN,ebN,k]    = thetaRtypes[type]
                            thetaS_in[eN,ebN,k]    = thetaStypes[type]


            if not foundAblock:
                print """didn't find block eN=%d """ % eN
                for ebN in range(x_in.shape[1]):
                    for k in range(x_in.shape[1]):
                        print """x[%d,%d,%d]=%s blockTests= """ % (eN,ebN,k,x_in[eN,ebN,k,0:2])
                        for nB in range(nHetBlocks):
                            inblock= (x_in[eN,ebN,k,0] >= heterogeneityBlocks[nB]['Xm'][0]-bndEps and
                                      x_in[eN,ebN,k,0] <= heterogeneityBlocks[nB]['Xp'][0]+bndEps and
                                      x_in[eN,ebN,k,1] >= heterogeneityBlocks[nB]['Xm'][1]-bndEps and
                                      x_in[eN,ebN,k,1] <= heterogeneityBlocks[nB]['Xp'][1]+bndEps)
                            print """BLOCK %d x >= %g x <=%g y >= %g y <= %g = %s """ % (nB,
                                                                                         heterogeneityBlocks[nB]['Xm'][0],
                                                                                         heterogeneityBlocks[nB]['Xp'][0],
                                                                                         heterogeneityBlocks[nB]['Xm'][1],
                                                                                         heterogeneityBlocks[nB]['Xp'][1],
                                                                                         inblock)
            assert foundAblock == True, "foundAblock eN=%d x=%s " % (eN,x_in[eN,:,:,0:2])

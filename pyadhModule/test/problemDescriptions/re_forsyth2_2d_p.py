from pyadh import *
from pyadh.default_p import *

nd = 2

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
alphaVGtypes = numpy.zeros((nMediaTypes,),'d')
nVGtypes     = numpy.zeros((nMediaTypes,),'d')
thetaStypes  = numpy.zeros((nMediaTypes,),'d')
thetaRtypes  = numpy.zeros((nMediaTypes,),'d')
thetaSRtypes = numpy.zeros((nMediaTypes,),'d')
KsTypes      = numpy.zeros((nMediaTypes,),'d')

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

thetaSRtypes = thetaStypes-thetaRtypes

#set dimensions for types
viscosity     = 8.9e-4  #kg/(m*s)
density       = 998.2   #kg/m^3
gravity       = 9.8     #m/s^2
beta          = density*gravity*4.524e-10
m_per_s_by_m_per_d = 1.1574074e-5
lengthScale   = 1.0     #m
timeScale     = 1.0     #d #1.0/sqrt(g*lengthScale)

dimensionless_conductivity  = numpy.zeros(KsTypes.shape,'d')
dimensionless_alpha         = numpy.zeros(alphaVGtypes.shape,'d')
dimensionless_density  = 1.0
dimensionless_gravity  = numpy.array([0.0,
                                        -1.0,
                                        0.0])
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

def setParams(x_in,vgm_n_in,vgm_alpha_in,Ks_in,thetaR_in,thetaSR_in):
    #brute force
    vgm_n_in.flat[:]     = nVGtypes[0]
    vgm_alpha_in.flat[:] = alphaVGtypes[0]
    Ks_in.flat[:] = dimensionless_conductivity[0]
    thetaR_in.flat[:] = thetaRtypes[0]
    thetaSR_in.flat[:] = thetaStypes[0]-thetaRtypes[0]
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
                        vgm_n_in[eN,k]     = nVGtypes[type]
                        vgm_alpha_in[eN,k] = alphaVGtypes[type]
                        Ks_in[eN,k]        = dimensionless_conductivity[type]
                        thetaR_in[eN,k]    = thetaRtypes[type]
                        thetaSR_in[eN,k]   = thetaStypes[type]-thetaRtypes[type]
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
                    Viewers.datFile.write("%12.5e %12.5e %12.5e \n" % (x_in[eN,k,0],x_in[eN,k,1],vgm_n_in[eN,k]))
            Viewers.datFile.write("\n \n")
            cmd = "set dgrid3d %d,%d,%g; set contour base; set term x11 %i; splot \'%s\' index %i with lines title \"%s\" \n" % (dgridx,dgridy,dgridp,
                                                                                                                                 Viewers.windowNumber,
                                                                                                                                 Viewers.datFilename,
                                                                                                                                 Viewers.plotNumber,
                                                                                                                                 'vgm-n')
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
                            vgm_n_in[eN,ebN,k]     = nVGtypes[type]
                            vgm_alpha_in[eN,ebN,k] = alphaVGtypes[type]
                            Ks_in[eN,ebN,k]        = dimensionless_conductivity[type]
                            thetaR_in[eN,ebN,k]    = thetaRtypes[type]
                            thetaSR_in[eN,ebN,k]   = thetaStypes[type]-thetaRtypes[type]


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

#mwf orig
#coefficients = ConservativeHeadRichardsMualemVanGenuchtenBlockHet(hydraulicConductivity=dimensionless_conductivity,
#                                                                  gravity=dimensionless_gravity,
#                                                                  density=dimensionless_density,
#                                                                  setParamsFunc=setParams)
upwindFlag = 1

# coefficients = ConservativeHeadRichardsMualemVanGenuchtenBlockHetV2withUpwind(nd,
#                                                                               KsTypes,
#                                                                               nVGtypes,
#                                                                               alphaVGtypes,
#                                                                               thetaRtypes,
#                                                                               thetaSRtypes,
#                                                                               gravity*dimensionless_gravity,
#                                                                               dimensionless_density,#need to fix compressibility
#                                                                               beta,
#                                                                               upwindFlag)
coefficients = ConservativeHeadRichardsMualemVanGenuchtenBlockHetV2withUpwind(nd,
                                                                              KsTypes,
                                                                              nVGtypes,
                                                                              alphaVGtypes,
                                                                              thetaRtypes,
                                                                              thetaSRtypes,
                                                                              gravity*dimensionless_gravity,
                                                                              dimensionless_density,#need to fix compressibility
                                                                              beta,
                                                                              upwindFlag,
                                                                              sd=1)
    
pondingPressure=0.1
initialPressure=-880.0e3 / density / gravity #convert kPa to m
rechargeRate =  -2.0e-2#- means into domain m/d
rechargeXboundary = 2.25# L[0]/3.0  where recharge stops
waterTableHeight = L[1]/3.0

eps = 1.0e-8
def getDBC_2D_Richards_HydroRight(x,flag):
    if abs(x[0]-L[0]) < eps:
        return lambda x,t: (x[1]-waterTableHeight)*dimensionless_gravity[1]*dimensionless_density
def getDBC_2D_Richards_Box(x,flag):
    pass
def getDBC_2D_Richards_top(x,flag):
    if abs(x[1]- L[1]) < eps:
        if x[0] < rechargeXboundary:
            #mwf debug
            #print """setting recharge rate at x=%s """ % x
            return lambda x,t: pondingPressure

#stupid V2 velocity postprocessing in testStuff requires explicit
#inclusion of no flux boundaries
def getRecharge_2D_UpperLeft(x,flag):
    if abs(x[1]- L[1]) < eps:
        if x[0] < rechargeXboundary:
            #mwf debug
            #print """setting recharge rate at x=%s """ % x
            return lambda x,t: rechargeRate
        else:
            return lambda x,t: 0.0
    #include no flow bcs explicitly for testing post processing
    if abs(x[0]-0.0) < eps or abs(x[1]-0.0) < eps:
        return lambda x,t:0.0

def getDummyFlux(x,flag):
    pass

#dirichletConditions = {0:getDBC_2D_Richards_HydroRight}
#dirichletConditions = {0:getDBC_2D_Richards_Box}
dirichletConditions = {0:getDBC_2D_Richards_top}

class ConstIC_2D_Richards:
    def uOfXT(self,x,t):
        if abs(x[1]-L[1]) < eps and x[0] < rechargeXboundary:
            return pondingPressure
        return initialPressure
class LinearIC_2D_Richards:
    def uOfXT(self,x,t):
        return (pondingPressure-0.0)/L[1]*(x[1]-0.0) + 0.0
class HydroIC_2D_Richards:
    def uOfXT(self,x,t):
        return (x[1]-waterTableHeight)*dimensionless_gravity[1]*dimensionless_density

initialConditions  = {0:ConstIC_2D_Richards()}
#initialConditions  = {0:HydroIC_2D_Richards()}

fluxBoundaryConditions = {0:'setFlow'}#{0:'noFlow'}

#advectiveFluxBoundaryConditions =  {0:getRecharge_2D_UpperLeft}#{}
advectiveFluxBoundaryConditions =  {0:getDummyFlux}#{}

diffusiveFluxBoundaryConditions = {0:{0:getDummyFlux}}#{0:{0:getRecharge_2D_UpperLeft}} #{0:{}}

T= 30./timeScale

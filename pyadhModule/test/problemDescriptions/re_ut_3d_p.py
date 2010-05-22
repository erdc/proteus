from pyadh import *
from pyadh.default_p import *

nd = 3

L=(15.0,15.0,15.0) #ft

analyticalSolution = None

nMediaTypes   = 2
pdBCtypes     = numpy.zeros((nMediaTypes,),'d')
lambdaBCtypes = numpy.zeros((nMediaTypes,),'d')
thetaStypes   = numpy.zeros((nMediaTypes,),'d')
thetaRtypes   = numpy.zeros((nMediaTypes,),'d')
KsTypes       = numpy.zeros((nMediaTypes,),'d')

it = 0;
thetaStypes[it]  = 0.30;
thetaRtypes[it]  = 0.03;
KsTypes[it]      = 10.0; #ft/day
pdBCtypes[it]    = 1.5;#ft
lambdaBCtypes[it]= 2.0;


it = 1;
thetaStypes[it]   = thetaStypes[0];
thetaRtypes[it]   = thetaRtypes[0];
KsTypes[it]       = 1.0; #ft/day
#KsTypes[it]       = KsTypes[0]; #ft/day
pdBCtypes[it]     = pdBCtypes[0]*math.sqrt(KsTypes[0]/KsTypes[it]);
lambdaBCtypes[it] = lambdaBCtypes[0];

#set up heterogeneity blocks, zone ids are base 1
nHetBlocks = 2
heterogeneityBlocks = {}

#7.5 < z <= 15
zone = 0
heterogeneityBlocks[zone] = {};
heterogeneityBlocks[zone]['Xm']  = [0.0,0.0,7.5]; heterogeneityBlocks[zone]['Xp'] = [L[0],L[1],L[2]]
heterogeneityBlocks[zone]['type']= 0
#0.0 <= z <=7.5
zone = 1
heterogeneityBlocks[zone] = {};
heterogeneityBlocks[zone]['Xm']  = [0.0,0.0,0.0]; heterogeneityBlocks[zone]['Xp'] = [L[0],L[1],7.5]
heterogeneityBlocks[zone]['type']= 1

plotHet = False

bndEps = 1.0e-8

def setParams(x_in,bc_lambda_in,bc_pd_in,Ks_in,thetaR_in,thetaSR_in):
    #brute force
    bc_lambda_in.flat[:]     = lambdaBCtypes[0]
    bc_pd_in.flat[:]         = pdBCtypes[0]
    Ks_in.flat[:]            = KsTypes[0]
    thetaR_in.flat[:]        = thetaRtypes[0]
    thetaSR_in.flat[:]       = thetaStypes[0]-thetaRtypes[0]
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
                    for k in range(x_in.shape[1]):
                        type = heterogeneityBlocks[nB]['type']
                        bc_lambda_in[eN,k] = lambdaBCtypes[type]
                        bc_pd_in[eN,k]     = pdBCtypes[type]
                        Ks_in[eN,k]        = KsTypes[type]
                        thetaR_in[eN,k]    = thetaRtypes[type]
                        thetaSR_in[eN,k]   = thetaStypes[type]-thetaRtypes[type]
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
                    Viewers.datFile.write("%12.5e %12.5e %12.5e \n" % (x_in[eN,k,0],x_in[eN,k,1],bc_lambda_in[eN,k]))
            Viewers.datFile.write("\n \n")
            cmd = "set dgrid3d %d,%d,%g; set contour base; set term x11 %i; splot \'%s\' index %i with lines title \"%s\" \n" % (dgridx,dgridy,dgridp,
                                                                                                                                 Viewers.windowNumber,
                                                                                                                                 Viewers.datFilename,
                                                                                                                                 Viewers.plotNumber,
                                                                                                                                 'bcb-lambda')
            Viewers.cmdFile.write(cmd)
            Viewers.viewerPipe.write(cmd)
            Viewers.plotNumber+=1
            Viewers.windowNumber += 1
            for eN in range(x_in.shape[0]):
                for k in range(x_in.shape[1]):
                    Viewers.datFile.write("%12.5e %12.5e %12.5e \n" % (x_in[eN,k,0],x_in[eN,k,1],bc_pd_in[eN,k]))
            Viewers.datFile.write("\n \n")
            cmd = "set dgrid3d %d,%d,%g; set contour base; set term x11 %i; splot \'%s\' index %i with lines title \"%s\" \n" % (dgridx,dgridy,dgridp,
                                                                                                                                 Viewers.windowNumber,
                                                                                                                                 Viewers.datFilename,
                                                                                                                                 Viewers.plotNumber,
                                                                                                                                 'bcb-pd')
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
                            bc_lambda_in[eN,ebN,k]     = lambdaBCtypes[type]
                            bc_pd_in[eN,ebN,k] = pdBCtypes[type]
                            Ks_in[eN,ebN,k]        = KsTypes[type]
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

coefficients = ConservativeHeadRichardsBrooksCoreyBurdine(hydraulicConductivity=KsTypes[0],
                                                          gravity=numpy.array([0.0,-1.0,0.0]),
                                                          density=1.0,
                                                          thetaS=thetaStypes[0],
                                                          thetaR=thetaRtypes[0],
                                                          lambdab=lambdaBCtypes[0],
                                                          pd=pdBCtypes[0])
coefficients = ConservativeHeadRichardsBrooksCoreyBurdineHet(hydraulicConductivity=1.0,
                                                            gravity=numpy.array([0.0,0.0,-1.0]),
                                                            density=1.0,
                                                            setParamsFunc=setParams)
pondingHead = 1.0 #ft
initialHead =-5.0 #ft
bottomHead  =-5.0 #ft


eps = 1.0e-8
def getDBC_3D_Richards(x):
    if x[2] >= L[2] - eps:
        return lambda x,t: pondingHead
    elif x[2] <= eps:
        return lambda x,t: bottomHead

dirichletConditions = {0:getDBC_3D_Richards}

class ConstIC_3D_Richards:
    def uOfXT(self,x,t):
        return initialHead

initialConditions  = {0:ConstIC_3D_Richards()}

fluxBoundaryConditions = {0:'noFlow'}

T= 1.0 # day

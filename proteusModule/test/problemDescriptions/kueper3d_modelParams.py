from pyadh import *
import numpy
import KueperFrindGenPoly3D

model = 'BCB'#'BCB'#'BCB' #'VGM'

# Model Parameters
rhow = 997.0        # density wetting      (Water 997 kg/m^3)
rhon = 1.205        # density nonwetting   (Air 1.205 kg/m^3)
muw  = 1.002e-3     # viscosity nonwetting (Water := 1.002e-3 kg/m s)
mun  = 1.81e-5      # viscosity wetting    (Air := 1.81e-5 kg/ m s)

Temp = 288.15       # Temperature
R    = 8.314        # Ideal Gas Constant   (N/mol k)
M    = 2.896e-2     # Molar Mass           (Air := 2.896e-2)
p_o  = 1.013e5      # (N/m^2)

gmag  = 9.8         # magnitude of gravity             (Earth 9.8 m/s^2) 
m_per_d_to_m_per_s = 1.1574074e-5



#make a "type for each block"
top   = 0.5
right = 0.7
depth = 0.7/10.0
L=(right,depth,top)
ny    = 2
source_x = (right/3.0,right*2.0/3.0)
source_y = (0.0,depth)
source_type = 12 #block id for source location
bottom_type = 0  #bottom block id

nBlocks = 20
blockLeft  = numpy.array([0.00, 0.00, 0.05, 0.10, 0.20, 0.50, 0.60, 0.20, 0.45, 0.65, 0.25, 0.10, 0.05, 
                          0.10, 0.50, 0.10, 0.20, 0.35, 0.10, 0.35])
blockRight = numpy.array([0.70, 0.05, 0.10, 0.20, 0.45, 0.65, 0.65, 0.50, 0.50, 0.70, 0.35, 0.60, 0.65, 
                          0.20, 0.60, 0.25, 0.50, 0.60, 0.60, 0.60])
blockFront = numpy.array([0.00, 0.05, 0.05, 0.05, 0.10, 0.05, 0.15, 0.05, 0.10, 0.05, 0.20, 0.30, 0.40, 
                          0.35, 0.35, 0.20, 0.35, 0.20, 0.15, 0.25])
blockBack  = numpy.array([0.05, 0.50, 0.40, 0.15, 0.15, 0.15, 0.40, 0.10, 0.15, 0.50, 0.30, 0.35, 0.50, 
                          0.40, 0.40, 0.30, 0.40, 0.25, 0.20, 0.30])


thetaS_types= numpy.array([0.41, 0.4, 0.41, 0.41, 0.41, 0.41, 0.41, 0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 
                           0.40, 0.40, 0.39, 0.39, 0.39, 0.39, 0.41])
thetaR_types= numpy.array([0.07749, 0.0312, 0.07749, 0.07749, 0.07749, 0.07749, 0.07749, 0.0312, 0.0312, 
                           0.0312, 0.0312, 0.0312, 0.0312, 0.0312, 0.0312, 0.03822, 0.03822, 0.03822, 
                           0.02691, 0.07749])
omega_types = thetaS_types[:]
sw_max_types= numpy.ones((nBlocks,),'d')
sw_min_types= thetaR_types/thetaS_types
bc_lambda_types = numpy.array([3.30, 3.86, 3.30, 3.30, 3.30, 3.30, 3.30, 3.86, 3.86, 3.86, 3.86, 3.86, 
                               3.86, 3.86, 3.86, 2.49, 2.49, 2.49, 3.51, 3.30])
bc_pd_types     = numpy.array([.3310, 0.0377, 0.3310, 0.3310, 0.3310, 0.3310, 0.3310, 0.0377, 0.0377, 
                                0.0377, 0.0377, 0.0377, 0.0377, 0.0377, 0.0377, 0.1350, 0.1350, 0.1350, 0.0443, 0.3310])
#m^2
perm_types      = numpy.array([8.19e-12, 5.04e-10, 8.19e-12, 8.19e-12, 8.19e-12, 8.19e-12, 8.19e-12, 5.04e-10, 5.04e-10, 
                               5.04e-10, 5.04e-10, 5.04e-10, 5.04e-10, 5.04e-10, 5.04e-10, 5.26e-11, 5.26e-11, 5.26e-11, 
                               2.05e-10, 8.19e-12])

#m/s
Ksw_types       = perm_types*rhow*gmag/muw

#default values for VGM
mvg_n_types     = bc_lambda_types + 1.0
mvg_alpha_types = 1.0/bc_pd_types
mvg_m_types     = 1.0 - 1.0/mvg_n_types

def seBCB(psic,pdBC,lamBC):
    if psic <= pdBC: return 1.0
    tmp1 = pow(pdBC/psic,lamBC)
    return min(max(tmp1,0.0),1.0)
def pcBCB(se,pdBC,lamBC):
    if se >= 1.0: return 0.0
    se_cutOff = max(se,1.0e-4)
    tmp1 = pow(se_cutOff,-1.0/lamBC)
    return pdBC*tmp1

#domain 
nd = 3
polyfile = "kueper3d"
boundaryTags = KueperFrindGenPoly3D.genPoly(polyfile,
                                            depth=depth,
                                            ny=ny,
                                            source_x=source_x,
                                            source_y=source_y)

def findAblock(x):
    bndEps = 1.0e-6
    foundAblock=False
    
    for nB in range(len(blockLeft)):
        inBlock = False
        if  (x[0] >= blockLeft[nB]-bndEps and
             x[0] <= blockRight[nB]+bndEps and
             x[2] >= blockFront[nB]-bndEps and
             x[2] <= blockBack[nB]+bndEps):
            inBlock = True
            break
    if inBlock:
        foundAblock = True
        return nB
    #mwf debug
    if not foundAblock:
        print "didnt find a block for x=%s" % x
        for nB in range(len(blockLeft)):
            print "nB=%d blockLeft= %s, blockRight=%s blockFront=%s blockBack=%s" % (nB,blockLeft[nB],blockRight[nB],
                                                                                     blockFront[nB],blockBack[nB])
    assert foundAblock == True, "foundAblock=%s x=%s " % (foundAblock,x)
    return None

g = [0.0,0.0,-gmag]
#boundary conditions


Se_top = 0.8
Sw_top = Se_top*(sw_max_types[source_type]-sw_min_types[source_type]) + sw_min_types[source_type]
psi_top= 0.0

waterTable = -0.1
psiTable   = 0.0
#psi_w at z=0
psi_bottom = psiTable + g[2]/gmag*(0.0-waterTable)#waterTable - 0.0
psinTable  = 0.0
psin_bottom= psinTable + rhon/rhow*g[2]/gmag*(0.0-waterTable)#waterTable - 0.0
pc_bottom  = psin_bottom-psi_bottom

Se_bottom =  pcBCB(pc_bottom,bc_pd_types[bottom_type],bc_lambda_types[bottom_type])#0.05#1.0e-3                  # effective saturation bottom
Sw_bottom= Se_bottom*(sw_max_types[bottom_type]-sw_min_types[bottom_type]) + sw_min_types[bottom_type]


eps = 1.0e-6

Sw_pad=1.0e-3

def sw_min(x,t):
    block=findAblock(x)
    return sw_min_types[block]+Sw_pad

def air0(x,t):
    block = findAblock(x)
    return - pcBCB(sw_min_types[block]+Sw_pad,bc_pd_types[block],bc_lambda_types[block])

def getDBC_sw(x,tag):
    if tag==boundaryTags['source']:
        return lambda x,t: Sw_top
    if tag in [boundaryTags['left'],boundaryTags['right']]:
        return sw_min
#     if (x[2] >= L[2]-eps and
#         source_y[0] <= x[1] and x[1] <= source_y[1] and 
#         source_x[0] <= x[0] and x[0] <= source_x[1] and tag==boundaryTags['source']):
#         return lambda x,t: Sw_top
#     if x[1] <= eps and tag==boundaryTags['bottom']:
#         return lambda x,t: Sw_bottom

def getDBC_psiw(x,tag):
    if tag==boundaryTags['source']:
        return lambda x,t: psi_top
    if tag in [boundaryTags['left'],boundaryTags['right']]:
        return air0
#     if (x[2] >= L[2]-eps and
#         source_y[0] <= x[1] and x[1] <= source_y[1] and 
#         source_x[0] <= x[0] and x[0] <= source_x[1] and tag==boundaryTags['source']):
#         return lambda x,t: psi_top
#     if x[1] <= eps and tag==boundaryTags['bottom']:
#         return lambda x,t: psi_bottom

class sw_IC:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        block = findAblock(x)
        sw = sw_min_types[block]+Sw_pad
        return sw
#         if (x[2] >= L[2]-eps and
#             source_y[0] <= x[1] and x[1] <= source_y[1] and 
#             source_x[0] <= x[0] and x[0] <= source_x[1]):
#             return Sw_top
#         #need to incorporate lens
#         psic = pc_bottom + (rhon/rhow - 1.0)*g[2]/gmag*x[2]
#         block = findAblock(x)
#         se = seBCB(psic,bc_pd_types[block],bc_lambda_types[block])
#         sw = se*(sw_max_types[block]-sw_min_types[block]) + sw_min_types[block]
#         return sw
class psiw_IC:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        block = findAblock(x)
        return - pcBCB(sw_min_types[block]+Sw_pad,bc_pd_types[block],bc_lambda_types[block])
        #return psiTable + g[2]/gmag*(x[2]-waterTable)



lengthScale   = 1.0     #m
timeScale     = 1.1574074e-05     #d #1.0/sqrt(g*lengthScale)



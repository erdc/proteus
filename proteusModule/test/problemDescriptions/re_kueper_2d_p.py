from pyadh import *
from pyadh.default_p import *
##\page Tests Test Problems
# \ref re_kueper_2d_p.py "Richards' equation, Brooks-Corey-Burdine constitutive relations, heterogeneous tank experiment"
#

##\ingroup test
#\file re_kueper_2d_p.py
#\brief Richards' equation with Brooks-Corey-Burdine Coefficients, hterogeneous tank experiment
#
#The equation formulation and coefficients are implemented in the ConservativeHeadRichardsBrooksCoreyBurdineHet
#class. The initial/boundary conditions are
#\f{eqnarray*}
#\psi(x,0) = -x \rho g \\
#\psi(0,t) = 0 \\
#\psi(10,t)= 0.1
#\f}
#

nd = 2

polyfile = "kueper"
top = 0.5
right = 0.7
analyticalSolution = None

viscosity     = 8.9e-4  #kg/(m*s)
density       = 998.2   #kg/m^3
gravity       = 9.8     #m/s^2
m_per_s_by_m_per_d = 1.1574074e-5
permeability  = (5.04*m_per_s_by_m_per_d)*viscosity/(gravity*density)  #m^2
print 'perm',permeability
thetaS        = 0.301   #-
thetaR        = 0.093   #-
mvg_alpha     = 5.47    #1/m
mvg_n         = 4.264
mvg_m         = 1.0 - 1.0/mvg_n
lengthScale   = 1.0     #m
timeScale     = 1.1574074e-05     #d #1.0/sqrt(g*lengthScale)
bc_lambda     = mvg_n-1.0
bc_pd         = 1.0/mvg_alpha
print 'pd',bc_pd
#make non-dimensional
dimensionless_conductivity  = (timeScale*density*gravity*permeability/(viscosity*lengthScale))/m_per_s_by_m_per_d
print 'Ks',dimensionless_conductivity
dimensionless_density  = 1.0
dimensionless_gravity  = numpy.array([0.0,
                                      -1.0,
                                      0.0])
dimensionless_pd    = bc_pd/lengthScale
blockLeft  = [0.00, 0.00, 0.05, 0.10, 0.20, 0.50, 0.60, 0.20, 0.45, 0.65, 0.25, 0.10, 0.05, 0.10, 0.50, 0.10, 0.20, 0.35, 0.10, 0.35]
blockRight = [0.70, 0.05, 0.10, 0.20, 0.45, 0.65, 0.65, 0.50, 0.50, 0.70, 0.35, 0.60, 0.65, 0.20, 0.60, 0.25, 0.50, 0.60, 0.60, 0.60]
blockFront = [0.00, 0.05, 0.05, 0.05, 0.10, 0.05, 0.15, 0.05, 0.10, 0.05, 0.20, 0.30, 0.40, 0.35, 0.35, 0.20, 0.35, 0.20, 0.15, 0.25]
blockBack  = [0.05, 0.50, 0.40, 0.15, 0.15, 0.15, 0.40, 0.10, 0.15, 0.50, 0.30, 0.35, 0.50, 0.40, 0.40, 0.30, 0.40, 0.25, 0.20, 0.30]
thetaS_block= [0.41, 0.4, 0.41, 0.41, 0.41, 0.41, 0.41, 0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 0.39, 0.39, 0.39, 0.39, 0.41]
thetaR_block= [0.07749, 0.0312, 0.07749, 0.07749, 0.07749, 0.07749, 0.07749, 0.0312, 0.0312, 0.0312, 0.0312, 0.0312, 0.0312, 0.0312, 0.0312, 0.03822, 0.03822, 0.03822, 0.02691, 0.07749]
bc_lambda_block= [3.30, 3.86, 3.30, 3.30, 3.30, 3.30, 3.30, 3.86, 3.86, 3.86, 3.86, 3.86, 3.86, 3.86, 3.86, 2.49, 2.49, 2.49, 3.51, 3.30]
bc_pd_block= [.3310, 0.0377, 0.3310, 0.3310, 0.3310, 0.3310, 0.3310, 0.0377, 0.0377, 0.0377, 0.0377, 0.0377, 0.0377, 0.0377, 0.0377, 0.1350, 0.1350, 0.1350, 0.0443, 0.3310]
perm_block= [8.19e-12, 5.04e-10, 8.19e-12, 8.19e-12, 8.19e-12, 8.19e-12, 8.19e-12, 5.04e-10, 5.04e-10, 5.04e-10, 5.04e-10, 5.04e-10, 5.04e-10, 5.04e-10, 5.04e-10, 5.26e-11, 5.26e-11, 5.26e-11, 2.05e-10, 8.19e-12] 
dimensionless_conductivity_block = [(timeScale*density*gravity*perm/(viscosity*lengthScale))/m_per_s_by_m_per_d for perm in perm_block]
print dimensionless_conductivity_block

plotHet = False #True
bndEps = 1.0e-8

def setParams(x_in,bcb_lambda_in,bcb_pd_in,Ks_in,thetaR_in,thetaS_in):
    bcb_lambda_in.flat[:] = bc_lambda
    bcb_pd_in.flat[:] = bc_pd
    Ks_in.flat[:] = dimensionless_conductivity
    thetaR_in.flat[:] = thetaR
    thetaS_in.flat[:] = thetaS
    #have to do some mesh dependent stuff here
    if len(x_in.shape) == 3: #on element quadrature
        for eN in range(x_in.shape[0]):
            foundAblock=False
            for nB in range(len(blockLeft)):
                inBlock = True
                for k in range(x_in.shape[1]):
                    if not (x_in[eN,k,0] >= blockLeft[nB]-bndEps and
                            x_in[eN,k,0] <= blockRight[nB]+bndEps and
                            x_in[eN,k,1] >= blockFront[nB]-bndEps and
                            x_in[eN,k,1] <= blockBack[nB]+bndEps):
                        inBlock = False
                        break
                if inBlock:
                    foundAblock = True
                    for k in range(x_in.shape[1]):
                        bcb_lambda_in[eN,k] = bc_lambda_block[nB]
                        bcb_pd_in[eN,k] = bc_pd_block[nB]
                        Ks_in[eN,k] = dimensionless_conductivity_block[nB]
                        thetaR_in[eN,k] = thetaR_block[nB]
                        thetaS_in[eN,k] = thetaS_block[nB]
        if not foundAblock:
            print """didn't find block eN=%d """ % eN
            for k in range(x_in.shape[1]):
                print """x[%d,%d]=%s blockTests= """ % (eN,k,x_in[eN,k,0:2])
                for nB in range(len(blockLeft)):
                    inblock = (x_in[eN,k,0] >= blockLeft[nB]-bndEps and
                               x_in[eN,k,0] <= blockRight[nB]+bndEps and
                               x_in[eN,k,1] >= blockFront[nB]-bndEps and
                               x_in[eN,k,1] <= blockBack[nB]+bndEps)
                    print """BLOCK %d x >= %g x <=%g y >= %g y <= %g = %s """ % (nB,
                                                                                 blockLeft[nB],
                                                                                 blockRight[nB],
                                                                                 blockFront[nB],
                                                                                 blockBack[nB],
                                                                                 inblock)
        assert foundAblock == True, "foundAblock eN=%d x=%s " % (eN,x_in[eN,:,0:2])
        from pyadh import Viewers 
        if plotHet and Viewers.viewerType == 'gnuplot':
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
                    Viewers.datFile.write("%12.5e %12.5e %12.5e \n" % (x_in[eN,k,0],x_in[eN,k,1],bcb_lambda_in[eN,k]))
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
                    Viewers.datFile.write("%12.5e %12.5e %12.5e \n" % (x_in[eN,k,0],x_in[eN,k,1],bcb_pd_in[eN,k]))
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
    elif len(x_in.shape) == 4: #on element boundary quadrature
        for eN in range(x_in.shape[0]):
            foundAblock = False
            for nB in range(len(blockLeft)):
                inBlock = True
                for ebN in range(x_in.shape[1]):
                    for k in range(x_in.shape[2]):
                        if not (x_in[eN,ebN,k,0] >= blockLeft[nB]-bndEps and
                                x_in[eN,ebN,k,0] <= blockRight[nB]+bndEps and
                                x_in[eN,ebN,k,1] >= blockFront[nB]-bndEps and
                                x_in[eN,ebN,k,1] <= blockBack[nB]+bndEps):
                            inBlock = False
                            break
                    if inBlock == False:
                        break
                if inBlock:
                    foundAblock = True
                    for ebN in range(x_in.shape[1]):
                        for k in range(x_in.shape[2]):
                            bcb_lambda_in[eN,ebN,k] = bc_lambda_block[nB]
                            bcb_pd_in[eN,ebN,k] = bc_pd_block[nB]
                            Ks_in[eN,ebN,k] = dimensionless_conductivity_block[nB]
                            thetaR_in[eN,ebN,k] = thetaR_block[nB]
                            thetaS_in[eN,ebN,k] = thetaS_block[nB]
                #in block
            #end blocks
            if not foundAblock:
                print """didn't find block eN=%d """ % eN
                for ebN in range(x_in.shape[1]):
                    for k in range(x_in.shape[1]):
                        print """x[%d,%d,%d]=%s blockTests= """ % (eN,ebN,k,x_in[eN,ebN,k,0:2])
                        for nB in range(len(blockLeft)):
                            inblock = (x_in[eN,ebN,k,0] >= blockLeft[nB]-bndEps and
                                       x_in[eN,ebN,k,0] <= blockRight[nB]+bndEps and
                                       x_in[eN,ebN,k,1] >= blockFront[nB]-bndEps and
                                       x_in[eN,ebN,k,1] <= blockBack[nB]+bndEps)
                            print """BLOCK %d x >= %g x <=%g y >= %g y <= %g = %s """ % (nB,
                                                                                         blockLeft[nB],
                                                                                         blockRight[nB],
                                                                                         blockFront[nB],
                                                                                         blockBack[nB],
                                                                                         inblock)
                #ebn
            #not found a block
            assert foundAblock == True, "foundAblock eN=%d x=%s " % (eN,x_in[eN,:,:,0:2])
        #eN
    #boundary quad


coefficients = ConservativeHeadRichardsBrooksCoreyBurdineHet(hydraulicConductivity=dimensionless_conductivity,
                                                             gravity=dimensionless_gravity,
                                                             density=dimensionless_density,
                                                             setParamsFunc=setParams)
# coefficients = ConservativeHeadRichardsBrooksCoreyBurdine(hydraulicConductivity=dimensionless_conductivity,
#                                                           gravity=dimensionless_gravity,
#                                                           density=dimensionless_density,
#                                                           thetaS=thetaS,
#                                                           thetaR=thetaR,
#                                                           lambdab = bc_lambda,
#                                                           pd = dimensionless_pd)
pondingPressure= 0.0 #-0.01 #0.1

def getDBC_2D_Richards_KueperFrind_Shock(x,flag):
    if x[1] == top:
        if (x[0] >= right/3.0 and
            x[0] <= 2.0*right/3.0):
            return lambda x,t: pondingPressure
    return None
#         else:
#             return None
#             #return lambda x,t: x[1]*dimensionless_gravity[1]*dimensionless_density
#     elif x[0] == 0.0 or x[0] == right or x[1] == 0.0:
#         return lambda x,t: x[1]*dimensionless_gravity[1]*dimensionless_density

dirichletConditions = {0:getDBC_2D_Richards_KueperFrind_Shock}

class ShockIC_2D_Richards:
    def uOfXT(self,x,t):
        if x[1] == top:
            if (x[0] >= right/3.0 and
                x[0] <= 2.0*right/3.0):
                return pondingPressure
            else:
                return x[1]*dimensionless_gravity[1]*dimensionless_density
        else:
            return x[1]*dimensionless_gravity[1]*dimensionless_density

initialConditions  = {0:ShockIC_2D_Richards()}

fluxBoundaryConditions = {0:'noFlow'}

def getAFBC_2D_Richards_KueperFrind_Shock(x,flag):
    if x[1] == top:
        if (x[0] >= right/3.0 and
            x[0] <= 2.0*right/3.0):
            return None
    return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC_2D_Richards_KueperFrind_Shock}

def getDFBC_2D_Richards_KueperFrind_Shock(x,flag):
    if x[1] == top:
        if (x[0] >= right/3.0 and
            x[0] <= 2.0*right/3.0):
            return None
    return lambda x,t: 0.0

diffusiveFluxBoundaryConditions = {0:{0:getDFBC_2D_Richards_KueperFrind_Shock}}

T = 0.5#1.0e-2/timeScale

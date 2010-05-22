#mwf try to unify through twp_darcy_modelParams from impes_modelParams_kueper_2d import *
from twp_darcy_modelParams import *

blockLeft  = [0.00, 0.00, 0.05, 0.10, 0.20, 0.50, 0.60, 0.20, 0.45, 0.65, 0.25, 0.10, 0.05, 0.10, 0.50, 0.10, 0.20, 0.35, 0.10, 0.35]
blockRight = [0.70, 0.05, 0.10, 0.20, 0.45, 0.65, 0.65, 0.50, 0.50, 0.70, 0.35, 0.60, 0.65, 0.20, 0.60, 0.25, 0.50, 0.60, 0.60, 0.60]
blockFront = [0.00, 0.05, 0.05, 0.05, 0.10, 0.05, 0.15, 0.05, 0.10, 0.05, 0.20, 0.30, 0.40, 0.35, 0.35, 0.20, 0.35, 0.20, 0.15, 0.25]
blockBack  = [0.05, 0.50, 0.40, 0.15, 0.15, 0.15, 0.40, 0.10, 0.15, 0.50, 0.30, 0.35, 0.50, 0.40, 0.40, 0.30, 0.40, 0.25, 0.20, 0.30]

thetaS_block    = [0.41, 0.4, 0.41, 0.41, 0.41, 0.41, 0.41, 0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 0.39, 0.39, 0.39, 0.39, 0.41]
thetaR_block    = [0.07749, 0.0312, 0.07749, 0.07749, 0.07749, 0.07749, 0.07749, 0.0312, 0.0312, 0.0312, 0.0312, 0.0312, 0.0312, 0.0312, 0.0312, 0.03822, 0.03822, 0.03822, 0.02691, 0.07749]
bc_lambda_block = [3.30, 3.86, 3.30, 3.30, 3.30, 3.30, 3.30, 3.86, 3.86, 3.86, 3.86, 3.86, 3.86, 3.86, 3.86, 2.49, 2.49, 2.49, 3.51, 3.30]
bc_pd_block     = [.3310, 0.0377, 0.3310, 0.3310, 0.3310, 0.3310, 0.3310, 0.0377, 0.0377, 0.0377, 0.0377, 0.0377, 0.0377, 0.0377, 0.0377, 0.1350, 0.1350, 0.1350, 0.0443, 0.3310]
perm_block      = [8.19e-12, 5.04e-10, 8.19e-12, 8.19e-12, 8.19e-12, 8.19e-12, 8.19e-12, 5.04e-10, 5.04e-10, 5.04e-10, 5.04e-10, 5.04e-10, 5.04e-10, 5.04e-10, 5.04e-10, 5.26e-11, 5.26e-11, 5.26e-11, 2.05e-10, 8.19e-12] 
Ksw_block = [rhow*gravity*perm/muw for perm in perm_block]

sw_max_block = [1.0 for i in range(len(thetaS_block))]
sw_min_block = [thetaR_block[i]/thetaS_block[i] for i in range(len(thetaS_block))]

print "Kueper_Frind Ksw= ",Ksw
print "Kueper_Frind Ksw_block= ",Ksw_block

plotHet = False
bndEps = 1.0e-8

def setParams_orig(x_in,bcb_lambda_in,bcb_pd_in,mvg_m_in,mvg_alpha_in,Ksw_in,thetaR_in,thetaS_in):
    bcb_lambda_in.flat[:] = bc_lambda
    bcb_pd_in.flat[:] = bc_pd
    mvg_m_in.flat[:] = mvg_m
    mvg_alpha_in.flat[:] = mvg_alpha
    Ksw_in.flat[:] = Ksw
    thetaR_in.flat[:] = thetaR
    thetaS_in.flat[:] = thetaS
    #cek debug
    #return
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
			mvg_m_in[eN,k] = 1.0 - (1.0/(bc_lambda_block[nB] + 1.0))
			mvg_alpha_in[eN,k] = 1.0/bc_pd_block[nB]
                        Ksw_in[eN,k] = Ksw_block[nB]
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
        if plotHet:
            from pyadh import Viewers #need to check that it's gnuplot
            dgridx=32; dgridy=32; dgridp=16;
            for eN in range(x_in.shape[0]):
                for k in range(x_in.shape[1]):
                    Viewers.datFile.write("%12.5e %12.5e %12.5e \n" % (x_in[eN,k,0],x_in[eN,k,1],Ksw_in[eN,k]))
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
	    for eN in range(x_in.shape[0]):
                for k in range(x_in.shape[1]):
                    Viewers.datFile.write("%12.5e %12.5e %12.5e \n" % (x_in[eN,k,0],x_in[eN,k,1],mvg_m_in[eN,k]))
            Viewers.datFile.write("\n \n")
            cmd = "set dgrid3d %d,%d,%g; set contour base; set term x11 %i; splot \'%s\' index %i with lines title \"%s\" \n" % (dgridx,dgridy,dgridp,
                                                                                                                                 Viewers.windowNumber,
                                                                                                                                 Viewers.datFilename,
                                                                                                                                 Viewers.plotNumber,
                                                                                                                                 'mvg-m')
            Viewers.cmdFile.write(cmd)
            Viewers.viewerPipe.write(cmd)
            Viewers.plotNumber+=1
            Viewers.windowNumber += 1
	    for eN in range(x_in.shape[0]):
                for k in range(x_in.shape[1]):
                    Viewers.datFile.write("%12.5e %12.5e %12.5e \n" % (x_in[eN,k,0],x_in[eN,k,1],mvg_alpha_in[eN,k]))
            Viewers.datFile.write("\n \n")
            cmd = "set dgrid3d %d,%d,%g; set contour base; set term x11 %i; splot \'%s\' index %i with lines title \"%s\" \n" % (dgridx,dgridy,dgridp,
                                                                                                                                 Viewers.windowNumber,
                                                                                                                                 Viewers.datFilename,
                                                                                                                                 Viewers.plotNumber,
                                                                                                                                 'mvg-alpha')
            Viewers.cmdFile.write(cmd)
            Viewers.viewerPipe.write(cmd)
            Viewers.plotNumber+=1
            Viewers.windowNumber += 1
            raw_input('press return to continue')
            Viewers.windowNumber -= 5
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
			    mvg_m_in[eN,ebN,k] = 1.0 - (1.0/(bc_lambda_block[nB] + 1.0))
			    mvg_alpha_in[eN,ebN,k] = 1.0/bc_pd_block[nB]
                            Ksw_in[eN,ebN,k] = Ksw_block[nB]
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

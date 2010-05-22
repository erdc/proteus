L=(1.0,1.0,1.0)
nd=2

model = 'VGM'       # Choose from {simp,BCB,BCM,VGB,VGM}
PorousMedia = 'sand' 


# Model Parameters
rhow = 1305.0        # density wetting      (Water 997 kg/m^3)
rhon = 998.2         # density nonwetting   (Air 1.205 kg/m^3)
muw  = 1.32e-3     # viscosity nonwetting (Water := 1.002e-3 kg/m s)
mun  = 1.005e-3      # viscosity wetting    (Air := 1.81e-5 kg/ m s)
g    = [0.0,9.8]    # gravity              (Earth 9.8 m/s^2) 

if (PorousMedia == 'sand'):
    sw_min    = 0.308
    sw_max    = 1.0
    omega     = 0.301
    Ks        = 5.040
    mvg_alpha = 5.470
    mvg_n     = 4.264
    bc_pd     = 0.120
    bc_lambda = 3.264
if (PorousMedia == 'loam'):
    sw_min    = 0.181
    sw_max    = 1.0
    omega     = 0.430
    Ks        = 0.250
    mvg_alpha = 3.600
    mvg_n     = 1.560
    bc_pd     = 0.050
    bc_lambda = 0.560
if (PorousMedia == 'clay'):
    sw_min    = 0.232
    sw_max    = 1.0
    omega     = 0.410
    Ks        = 0.062
    mvg_alpha = 1.900
    mvg_n     = 1.310
    bc_pd     = 0.044
    bc_lambda = 0.310

Ksw = Ks*1.1574074e-5   # Hydrolic Conductivity  related to medium permeability  (3.5e-10 to 3.5e-1 for Water)
mvg_m = 1.0-1.0/mvg_n
Q_m_per_d = 0.1
q    = Q_m_per_d*1.1574074e-5
q = 0.0

k = 10.0  
a = 0.025
FudgeFactor = 0.001

T = 5000.0

plotHet = False


def setParams(x_in,bcb_lambda_in,bcb_pd_in,mvg_m_in,mvg_alpha_in,Ksw_in,thetaR_in,thetaS_in):
    bcb_lambda_in.flat[:] = bc_lambda
    bcb_pd_in.flat[:] = bc_pd
    mvg_m_in.flat[:] = mvg_m
    mvg_alpha_in.flat[:] = mvg_alpha
    Ksw_in.flat[:] = Ksw
    thetaR_in.flat[:] = 0.0 #thetaR
    thetaS_in.flat[:] = omega #thetaS     
    
    bc_lambda_2 = bc_lambda # 3.30
    bc_pd_2     = bc_pd     # 0.3310
    mvg_m_2     = mvg_m     # 1.0 - (1.0/(bc_lambda_2 + 1.0))
    mvg_alpha_2 = mvg_alpha # 1.0/bc_pd_2
    #Kprm        = 8.19e-12
    Ksw_2       = Ksw       # 997.0*9.8*Kprm/1.002e-3 
    thetaR_2    = 0
    omega_2     = omega
    
    
    if len(x_in.shape) == 3: #on element quadrature
    	for eN in range(x_in.shape[0]):
            for k in range(x_in.shape[1]):
                if (x_in[eN,k,1] >= 0.5):
	            bcb_lambda_in[eN,k] = bc_lambda_2
                    bcb_pd_in[eN,k]     = bc_pd_2
                    mvg_m_in[eN,k]      = mvg_m_2
                    mvg_alpha_in[eN,k]  = mvg_alpha_2
                    Ksw_in[eN,k]        = Ksw_2
                    thetaR_in[eN,k]     = 0.0 #thetaR
                    thetaS_in[eN,k]     = omega_2
		    
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
            for ebN in range(x_in.shape[1]):
                for k in range(x_in.shape[2]):
                    if (x_in[eN,ebN,k,1] >= 0.5):	            
		        bcb_lambda_in[eN,k] = bc_lambda_2
                        bcb_pd_in[eN,k]     = bc_pd_2
                        mvg_m_in[eN,k]      = mvg_m_2
                        mvg_alpha_in[eN,k]  = mvg_alpha_2 
                        Ksw_in[eN,k]        = Ksw_2 
                        thetaR_in[eN,k]     = 0.0             #thetaR
                        thetaS_in[eN,k]     = omega_2 
			
 

# This is a driver for the Oscher Functions
from TestFunctions import *
from SearchAlgorithm import *
import Gnuplot
import numpy
import pyadh
import sys
import os
from pyadh.TransportCoefficients import *
#from Edit_TransportCoefficients import *
# Plotting and Data options
showStepPlots   = 0 # Set to 0 for off and 1 for on 
HardCopyPlot    = 0 # Set to 0 for off and 1 for on 
SaveData        = 1 # Set to 0 for off and 1 for on 
PlotCoefficents = 0
SaveCoefDatFile = 1
PlotTitle  = "" # if Title is set plot is saved. 

# Set the model parameters to be given to the impes model for two phase flow. 

model = 'VGM'        # No Choice Here
PorousMedia = 'sand' # Choose from {sand,loam,clay}

# Set up a loop to allow for running a time loop
SingleOsherRun = 'no' # 'yes' or 'no'

T_start = 10.0
T_stop  = 500.0
dt = 10.0

T = T_start
while T <= T_stop:
    # Model Parameters
    m_per_s_by_m_per_d = 1.1574074e-5
    rhow = 997.0        # density wetting      (Water 997 kg/m^3)
    muw  = 1.002e-3     # viscosity nonwetting (Water := 1.002e-3 kg/m s)


    Temp = 288.15       # Temperature
    R    = 8.314        # Ideal Gas Constant   (N/mol k)
    M    = 2.896e-2     # Molar Mass           (Air := 2.896e-2)
    p_o  = 1.013e5      # (N/m^2)

    g    = Numeric.array([9.8,0.0,0.0])     #m/s^2
    gu   = Numeric.array([1.0,0.0,0.0])    
    FudgeFactor = 0.0000

    if (PorousMedia == 'sand'):
        sw_min    = 0.308
        sw_max    = 1.0
        omega     = 0.301
        Ks        = 5.040
        mvg_alpha = 5.470
        mvg_n     = 4.264
        bc_pd     = 0.120
        bc_lambda = 3.264
	thetaS    = omega
        thetaR    = sw_min*omega   
    if (PorousMedia == 'loam'):
        sw_min    = 0.181
        sw_max    = 1.0
        omega     = 0.430
        Ks        = 0.250
        mvg_alpha = 3.600
        mvg_n     = 1.560
        bc_pd     = 0.050
        bc_lambda = 0.560	
	thetaS    = omega
        thetaR    = sw_min*omega   
    if (PorousMedia == 'clay'):
        sw_min    = 0.232
        sw_max    = 1.0
        omega     = 0.410
        Ks        = 0.062
        mvg_alpha = 1.900
        mvg_n     = 1.310
        bc_pd     = 0.044
        bc_lambda = 0.310	
	thetaS    = omega
        thetaR    = sw_min*omega   

    mvg_m     = 1.0-1.0/mvg_n
    Ksw = Ks*m_per_s_by_m_per_d  # Hydrolic Conductivity  related to medium permeability  (3.5e-10 to 3.5e-1 for Water)
     
    LHS_x   = -0.5    # Left hand side of problem domain
    RHS_x   = 0.5     # Right hand side of problem domain
    xShift  = 0.5     # Shift the xdomain over for plot matching
    Guess_s = 0.2     # Generic initial guess for fminbound
    Tol     = 1e-6    # Tolerance to pass to fminbound
    nsub    = 100     # Way of setting the refinement. 

    xLeft=LHS_x
    dx = (RHS_x - LHS_x)/nsub

    phase = 'potential' # 'saturation' or 'potential'
    phase = 'saturation'

    # Define the inital shock.
    FudgeFactor = 1e-5 
    LHS_s = 1.0-FudgeFactor
    RHS_s = 0.0+FudgeFactor

    coefficients = ConservativeSatRichardsMualemVanGenuchten(hydraulicConductivity=Ksw,
                                                          gravity=gu,#g,
                                                          density=1.0,#rhow,
                                                          thetaS=thetaS,
                                                          thetaR=thetaR,
                                                          alpha= mvg_alpha,
                                                          n = mvg_n,
                                                          m = mvg_m)
			    
    #coefficients.useC = False	    
    coefficients.useC = True  

    # Plot the different Coefficients
    # Choose from {u,m,dm,f,df,a,da,phi,dphi,r,dr,H,dH} 
    
    if(PlotCoefficents == 1):
        xaxis = 'u'
        yaxis = 'f'
	Title = qplotName
        PlotTitle = Title + "svf.eps"
        (xvec,yvec) = coefficients.PlotCoeff(6*nsub,T,xaxis,yaxis,PlotTitle)
	if(SaveCoefDatFile == 1): 
	    File = open(Title+'.dat' , 'w')
	    for i in range(len(xvec)): 
	        LineOfData = str(xvec[i]) + " " + str(yvec[i]) + "\n"
	        File.write(LineOfData)

        #yaxis = 'df'
        #xaxis = 'u'
        #PlotTitle = CodeType + "uvsdf.eps"
        #coefficients.PlotCoeff(nsub,T,xaxis,yaxis,PlotTitle)
        #xaxis = 'u'
        #yaxis = 'phi'
        #lotTitle = CodeType + "uvsphi.eps"
        #coefficients.PlotCoeff(nsub,T,xaxis,yaxis,PlotTitle)
        #xaxis = 'u'
        #yaxis = 'dphi'
        #PlotTitle = CodeType + "uvsdphi.eps"
        #coefficients.PlotCoeff(nsub,T,xaxis,yaxis,PlotTitle)



    """ Set model with f and F or with coefficients function """ 
    #f = BuckleyLeverett(0.5)
    #F = OsherFunc(LHS_s,RHS_s,fFunc=f,t=T,x=xLeft)
    F = OsherFuncCoef(LHS_s,RHS_s,fFunc=coefficients,t=T,x=xLeft)
    
    # comment jcc: Best I can tell x at the moment dosen't do anythin in the 
    # ConservativeSatRichardsMualemVanGenuchten transport coefficents c function. 
    # I set it this way for now to get around changing the function interface. 
    F.c[('x')] = 42.0
    
    solver = fminbound(FuncToMinimize=F,Tol=Tol)

    xVec=[]
    sVec=[]
    for i in range(nsub+1):
	    x = xLeft+i*dx
	    F.xi = x/T
            sVec2=[]
            fVec2=[]
            for i in range(nsub+1):
	        s2    = LHS_s+i*(RHS_s - LHS_s)/float(nsub)
	        Fval2 = F.getResidual(s2)
	        sVec2.append(s2)
	        fVec2.append(Fval2)

	    (s,Fval) = solver.solve(Guess_x=Guess_s)
	    #print "s = ",s," Fval = ",Fval
	    Guess_s = 0.5*(LHS_s+RHS_s)
	    #Guess_s = 0.25*LHS_s + 0.75*RHS_s
	    #Guess_s = Fval
	    xVec.append(x+xShift)
	    sVec.append(s)
   
	    if(showStepPlots==1):
                g = Gnuplot.Gnuplot(debug=1)
	        g.title('Plot of Function')     # (optional)
                g('set data style linespoints') # give gnuplot an arbitrary command
                d = Gnuplot.Data(Numeric.array(sVec2),Numeric.array(fVec2))
                g.xlabel('s')
                g.ylabel('f(s)')
	        dmin = Gnuplot.Data(s,Fval)
                g.plot(d,dmin)
	        raw_input('Please press return to continue\n')
	    
    g = Gnuplot.Gnuplot(debug=1)
    PlotTitle = 'Osher Solution for Model: ' + model + ' at time = ' + str(T)
    g.title(PlotTitle) # (optional)
    g('set data style lines') # give gnuplot an arbitrary command
    #xrangeString = 'set xrange ['+str(LHS_x+xShift)+':'+str(RHS_x+xShift)+']'
    #g(xrangeString)
    #g('set yrange [0:1]')
    d = Gnuplot.Data(Numeric.array(xVec),Numeric.array(sVec))
    g.xlabel('x')
    g.ylabel('s')
    g.plot(d)
    if(SingleOsherRun == 'yes'):
        raw_input('Please press return to continue...\n')

    # Plots and data storage for post processing. 
    DataName = 'Osher'+model
    if(SaveData == 1):
        OutFile = open(DataName+'.dat' , 'a')
	if (T > T_start):
	    OutFile.write("\n\n")
        for i in range(nsub+1):
            OutFile.write(str(xVec[i])+' '+str(sVec[i])+'\n')
    if(HardCopyPlot == 1):
        PlotName = DataName+'.eps'
        g.hardcopy(PlotName)
    
    T = T + dt
    print "T = ", T

OutFile.close()

if( SingleOsherRun == 'no'): 
    DatafileFullName = DataName + ".dat" 
    MakeAnimation = 'python makeOsherAnimation.py ' + DatafileFullName + " " + str(int((T_stop-T_start)/dt)) + " " + str(T_start) + " " + str(T_stop) 
    os.system(MakeAnimation)
    os.system('rm *MakeAnimations.eps')
    Cleanup = 'rm ' + DatafileFullName
    os.system(Cleanup) 


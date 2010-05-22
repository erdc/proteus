# This is a driver for the Oscher Functions
from TestFunctions import *
from SearchAlgorithm import *
import Gnuplot
import numpy
import pyadh
import sys
import os
from pyadh.TransportCoefficients import *

# Plotting and Data options
# Set to 0 for off and 1 for on
showStepPlots   = 0 # Plots the flux and shows the minimum value found across the domain.
HardCopyPlot    = 0 
SaveData        = 1 

PlotCoefficents = 0
SaveCoefDatFile = 1
PlotTitle  = ""        # if Title is set plot is saved. 

# Set the model parameters to be given to the impes model for two phase flow. 

model = 'VGM'        # Choose from {simp,BCB,BCM,VGB,VGM}
PorousMedia = 'sand' # Choose from {sand,loam,clay}

# Set up a loop to allow for running a time loop
SingleOsherRun = 'yes' # 'yes' or 'no'

T_start = 500.0
T_stop  = 500.0
dt = 10.0

T = T_start
while T <= T_stop:
    # Model Parameters
    rhow = 997.0        # density wetting      (Water 997 kg/m^3)
    rhon = 1.205        # density nonwetting   (Air 1.205 kg/m^3)
    muw  = 1.002e-3     # viscosity nonwetting (Water := 1.002e-3 kg/m s)
    mun  = 1.81e-5      # viscosity wetting    (Air := 1.81e-5 kg/ m s)

    Temp = 288.15       # Temperature
    R    = 8.314        # Ideal Gas Constant   (N/mol k)
    M    = 2.896e-2     # Molar Mass           (Air := 2.896e-2)
    p_o  = 1.013e5      # (N/m^2)

    g    = [9.8]        # gravity              (Earth 9.8 m/s^2) 
    Q_m_per_d = 0.0#0.4 
    q    = Q_m_per_d*1.1574074e-5
    q    = 0.0001
         # adjust the q_t (Total velocity) value
    a    = 1.0e-6

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

    mvg_m     = 1.0-1.0/mvg_n
    Ksw = Ks*1.1574074e-5   # Hydrolic Conductivity  related to medium permeability  (3.5e-10 to 3.5e-1 for Water)
    
    qplotName = model+PorousMedia+"q_t_1.0Ksw"
    q = 0.5*Ksw
    
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

    # Set up the modeling problem using the Transport Coefficients model. 
    coefficients = TwoPhaseFlow(q=q,
                            a=a,
			    Ksw=Ksw,
			    rhon=rhon,
			    rhow=rhow,
			    g=g,
			    mvg_alpha=mvg_alpha,
			    bc_lambda=bc_lambda,
			    bc_pd = bc_pd,
			    mvg_n = mvg_n,
			    mvg_m = mvg_m,
			    omega=omega,
			    mun=mun,
			    muw=muw,
			    sw_max=sw_max,
			    sw_min=sw_min,
			    M=M,
			    R=R,
			    Temp=Temp,
			    p_o=p_o,
			    model=model,
			    phase=phase)

			    
    #coefficients.useC = False	    
    coefficients.useC = True	    

    # Plot the different Coefficients
    # Choose from {u,m,dm,f,df,a,da,phi,dphi,r,dr,H,dH}. Note: There
    # is some code in the bottom of the "TwoPhaseFlow" coefficients class that
    # could easily be moved to other coefficents classes to make this work. 
    
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
    xrangeString = 'set xrange ['+str(LHS_x+xShift)+':'+str(RHS_x+xShift)+']'
    g(xrangeString)
    g('set yrange [0:1]')
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

OutFile.close()

if( SingleOsherRun == 'no'): 
    DatafileFullName = DataName + ".dat" 
    MakeAnimation = 'python makeOsherAnimation.py ' + DatafileFullName + " " + str(int((T_stop-T_start)/dt)) + " " + str(T_start) + " " + str(T_stop) 
    os.system(MakeAnimation)
    os.system('rm *MakeAnimations.eps')
    Cleanup = 'rm ' + DatafileFullName
    os.system(Cleanup) 


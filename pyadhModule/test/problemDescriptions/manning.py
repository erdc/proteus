from interactivePyADH import *
import pyadh
import math
import navier_stokes_manning_slope_2d_p
import navier_stokes_manning_slope_2d_n

pList[0] = navier_stokes_manning_slope_2d_p
nList[0]= navier_stokes_manning_slope_2d_n

pList[0].rho=1.0
pList[0].nu=1.0



slopes=[]
velocities=[]
Size=[0.7,1.0]
pList[0].Size=Size
heights=[]
numheights=3
for j in range(1,numheights+1):
    pList[0].Size[1]=0.05+float(j)/100.0

    for i in range(1,11):
        print 'RUN ' + str((j-1)*11+i) + ' of '+ str(numheights*11)
        A=0.0
        pList[0].slope=-float(i)*10.0**-4
        while A==0.0:
            #pList[0].angle=float(i)*math.pi/800.0
            #pList[0].slope= -math.tan(pList[0].angle)
            ###pList[0].slope=-float(i)/200.0
            from plogram2dDomain_tmp import *
            pList[0].domain=domain = plogram2D(Lx=pList[0].Size[0],Ly=pList[0].Size[1],slope=pList[0].slope)
            pList[0].domain.writePoly("plogram2D2D")
            ns = NumericalSolution.NS_base(so,pList,nList,sList,opts)
            runName= 'runName001'
            ns.calculateSolution(runName)
            A=nList[0].auxiliaryVariables[0].levelVlist[0][0]
            if A ==0.0:
                pList[0].slope=pList[0].slope+1.0e-3
            del ns



        slopes.append(pList[0].slope)
        velocities.append(math.sqrt(nList[0].auxiliaryVariables[0].levelVlist[0][0]**2 + nList[0].auxiliaryVariables[0].levelVlist[0][1]**2))
        heights.append(pList[0].Size[1])

        

        f=open('manning_results.py', 'w')
        f.write('numheights = ' + str(numheights)+ '\r' +
                'Size = ' +str(Size) + '\r' +
                'slopes = ' + str(slopes) + '\r' +
                'velocities = ' + str(velocities)+ '\r'+
                'heights = ' + str(heights) )
    f.close

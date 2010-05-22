from interactivePyADH import *
import pyadh
import navier_stokes_pome_2d_homog_p
import navier_stokes_pome_2d_homog_c0p1c0p1_n
import numpy

dir=True

pList[0] = navier_stokes_pome_2d_homog_p
nList[0]= navier_stokes_pome_2d_homog_c0p1c0p1_n

pList[0].dir=dir


ReynoldsNumbers=[]
PressureGradients=[]
inflows=[]

for i in range(101):
    pList[0].Re=0.1*float(i)
    pList[0].inflow=pList[0].Re*pList[0].nu*(1.0-pList[0].plantArea)/(2.0*pList[0].rad)
    ns = NumericalSolution.NS_base(so,pList,nList,sList,opts)
    runName= 'runName' + str(i)
    ns.calculateSolution(runName)
    
    ReynoldsNumbers.append(pList[0].Re)
    inflows.append(pList[0].inflow)
    print pList[0].Re
    #PressureGradients.append
    if dir:
         PressureGradients.append(nList[0].auxiliaryVariables[0].levelPlist[0][4])
    else:
        PressureGradients.append(nList[0].auxiliaryVariables[0].levelPlist[0][1])

f=open('forch_results.txt', 'w')
f.write('Reynolds = ' + str(ReynoldsNumbers) + '\r' +
        'inflows  =' + str(inflows) + '\r' +
        'Pressure  = ' + str(PressureGradients))
f.close

g=open('forch_results.m', 'w')
g.write('Reynolds = ' + str(ReynoldsNumbers) + ';'+ '\r' +
        'inflows  =' + str(inflows) + ';' + '\r' +
        'Pressure  = ' + str(PressureGradients)+ ';')
g.close

Reyn=numpy.zeros((len(inflows),1),'d')
inf=Reyn
Pres=Reyn

for i in range(len(ReynoldsNumbers)):
    Reyn[i]=ReynoldsNumbers[i]
    inf[i]=inflows[i]
    Pres[i]=PressureGradients[i]

h=open('forch_results.py','w')
## h.write('Pres = ' + str(Pres) + '\r'+
##         'Reyn = ' + str(Reyn) + '\r' +
##         'inf = ' + str(inf))
h.write('Reynolds = ' + str(ReynoldsNumbers) + ';'+ '\r' +
        'inflows  =' + str(inflows) + ';' + '\r' +
        'Pressure  = ' + str(PressureGradients)+ ';')
h.close

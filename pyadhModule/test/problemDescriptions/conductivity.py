

from interactivePyADH import *
import math
import numpy
import numpy.linalg
import pyadh
import stokes_pome_2d_homog_p
import stokes_pome_2d_homog_c0p1c0p1_n

pList[0] = stokes_pome_2d_homog_p
nList[0] = stokes_pome_2d_homog_c0p1c0p1_n

conductivityar = numpy.zeros((2,2), 'd')
conductivitylst = [[0.0, 0.0],[0.0,0.0]]

pList[0].dir = False
pList[0].nu=1.0
pList[0].rho=1.0
pList[0].pInflow=1000.0

ns = NumericalSolution.NS_base(so,pList,nList,sList,opts)
runName= 'run001'
ns.calculateSolution(runName)


conductivityar[0][0]=nList[0].auxiliaryVariables[0].levelVlist[0][0]
conductivityar[1][0]=nList[0].auxiliaryVariables[0].levelVlist[0][1]
conductivitylst[0][0]=nList[0].auxiliaryVariables[0].levelVlist[0][0]
conductivitylst[1][0]=nList[0].auxiliaryVariables[0].levelVlist[0][1]



pList[0].dir = True

ns = NumericalSolution.NS_base(so,pList,nList,sList,opts)
runName= 'run002'
ns.calculateSolution(runName)

conductivityar[1][1]=nList[0].auxiliaryVariables[0].levelVlist[0][1]
conductivityar[0][1]=nList[0].auxiliaryVariables[0].levelVlist[0][0]

conductivitylst[1][1]=nList[0].auxiliaryVariables[0].levelVlist[0][1]
conductivitylst[0][1]=nList[0].auxiliaryVariables[0].levelVlist[0][0]

conductivityar=0.5*(conductivityar+numpy.core.fromnumeric.transpose(conductivityar))/pList[0].pInflow

conductivitylst[0][0]=conductivityar[0][0]
conductivitylst[1][0]=conductivityar[1][0]
conductivitylst[0][1]=conductivityar[0][1]
conductivitylst[1][1]=conductivityar[1][1]
invconductivity= numpy.linalg.inv(conductivityar)
f=open('conductivity.txt', 'w')
f.write('Conductivity = ' + str(conductivityar) + ';')##  #'\r' + '\r'
##         #'InvConductivity = ' + str(invconductivity))
f.close

g=open('conductivity.m', 'w')
g.write('Conductivity = ' + str(conductivityar) + ';')##  #'\r' + '\r'
##         #'InvConductivity = ' + str(invconductivity))
g.close

h=open('conductivitycalc.py', 'w')
h.write('conductivityresult = ' + str(conductivitylst))
h.close

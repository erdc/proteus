import numpy
import numpy.linalg
import scipy
import scipy.linalg
import pylab
import conductivitycalc
import forch_results
import matplotlib
import matplotlib.pyplot

L=len(forch_results.Reynolds)

Mhat=998.2
That=1.0e6

R=numpy.zeros((L,1),'d')
P=numpy.zeros((L,1),'d')#((len(L),1),'d')
I=numpy.zeros((L,1),'d')#((len(forch_results.ReynoldsNumbers),1),'d')
H=numpy.zeros((L,1),'d')#((len(forch_results.ReynoldsNumbers),1),'d')
parab=numpy.zeros((L,1),'d')#((len(forch_results.ReynoldsNumbers),1),'d')
infsq=numpy.zeros((L,1),'d')#((len(forch_results.ReynoldsNumbers),1),'d')
iinfsq=numpy.zeros((L,2),'d')#((len(forch_results.ReynoldsNumbers),2),'d')

K=numpy.zeros((2,2), 'd')

K[0][0]=conductivitycalc.conductivityresult[0][0]
K[1][0]=conductivitycalc.conductivityresult[1][0]
K[0][1]=conductivitycalc.conductivityresult[0][1]
K[1][1]=conductivitycalc.conductivityresult[1][1]

Kinv = numpy.linalg.inv(K)

for i in range(L):
    R[i]=forch_results.Reynolds[i]
    I[i]=forch_results.inflows[i]/That
    P[i]=forch_results.Pressure[i]*Mhat/(That*That)
    H[i]=I[i]*Kinv[0][0]*Mhat/That
    parab[i]= P[i]-H[i]
    infsq[i]=I[i]*I[i]
    iinfsq[i][0]=I[i]
    iinfsq[i][1]=I[i]*I[i]


Beta=scipy.linalg.lstsq(infsq,parab)
B=Beta[0][0][0]/Mhat

KBeta=scipy.linalg.lstsq(iinfsq,P)

Fit1=numpy.zeros((L,1),'d')
Fit2=numpy.zeros((L,1),'d')

for i in range(L):
    Fit1[i]=Kinv[0][0]*Mhat/That*I[i]+Beta[0][0]*I[i]*I[i]
    Fit2[i]=KBeta[0][0][0]*I[i] +KBeta[0][1][0]*I[i]*I[i]



matplotlib.pyplot.figure(1)
matplotlib.pyplot.plot(R,P,'r', label='Navier-Stokes', linewidth=1)
matplotlib.pyplot.plot(R,H, 'b', label= 'Homogenization', linewidth=1)
matplotlib.pyplot.xlabel('Reynolds Number')
matplotlib.pyplot.ylabel('Pressure Gradient')
matplotlib.pyplot.legend(('Navier-Stokes', 'Homogenization'), loc='lower right')
matplotlib.pyplot.savefig('forch.eps')

str1='Quadr. Fit (Beta = ' +str(B) + ')'

matplotlib.pyplot.figure(2)
matplotlib.pyplot.plot(R,P,label='Navier-Stokes', linewidth=1)
matplotlib.pyplot.plot(R,Fit1,label=str1, linewidth=1)
matplotlib.pyplot.xlabel('Reynolds Number')
matplotlib.pyplot.ylabel('Pressure Gradient')
#matplotlib.pyplot.plot(R,Fit2)
matplotlib.pyplot.legend(('Navier-Stokes', str1), loc='lower right')
matplotlib.pyplot.savefig('forch2.eps')

import numpy
import numpy.linalg
import math
import scipy
import scipy.linalg
import pylab
import matplotlib
import matplotlib.pyplot
import manning_results

L=len(manning_results.slopes)




Slopes=numpy.zeros((L,1),'d')
Sizes=numpy.zeros((L,1),'d')
Velocities=numpy.zeros((L,1),'d')
Cf=numpy.zeros((L,1),'d')
A=numpy.zeros((L,3),'d')
A2=numpy.zeros((L,2),'d')
b=numpy.zeros((L,1),'d')
b2=numpy.zeros((L,1),'d')
b3=numpy.zeros((L,1),'d')
numheights=manning_results.numheights
heights=numpy.zeros((L,1),'d')


for i in range(L):
    Slopes[i]=-manning_results.slopes[i]
    Sizes[i]=manning_results.heights[i]
    Velocities[i]=manning_results.velocities[i]
    Cf[i]=1.0
    heights[i]=manning_results.heights[i]
    b[i]=math.log(Velocities[i])
    b2[i]=math.log(Velocities[i])-math.log(Slopes[i])
    b3[i]=math.log(Velocities[i])-0.5*math.log(Slopes[i])
    if Slopes[i]==0.0:
        A[i][0]==math.log(1.0e-20)
    else:
        A[i][0]=math.log(Slopes[i])
    A[i][1]=math.log(heights[i])
    A[i][2]=1.0
    A2[i][0]=math.log(heights[i])
    A2[i][1]=1.0
fit=scipy.linalg.lstsq(A,b)
fit2=scipy.linalg.lstsq(A2,b2)
fit3=scipy.linalg.lstsq(A2,b3)

L2=L/numheights
matplotlib.pyplot.figure(1)
for i in range(numheights):
    label='plot' +str(heights[i*L2])
    matplotlib.pyplot.plot(Slopes[i*L2:(i+1)*L2],Velocities[i*L2:(i+1)*L2], 'r',label=label, linewidth=1)
    matplotlib.pyplot.xlabel('-grad(H)')
    matplotlib.pyplot.ylabel('V')
#matplotlib.pyplot.legend('plot1', loc='lower right')
matplotlib.pyplot.savefig('manning.eps')

matplotlib.pyplot.figure(2)
for i in range(numheights):
    label='plot' +str(heights[i*L2])
    matplotlib.pyplot.plot(numpy.sqrt(Slopes[i*L2:(i+1)*L2])*math.sqrt(heights[i*L2]),Velocities[i*L2:(i+1)*L2], 'r',label=label, linewidth=1)
    matplotlib.pyplot.xlabel('sqrt(RS)')
    matplotlib.pyplot.ylabel('V')
#matplotlib.pyplot.legend('plot1', loc='lower right')
matplotlib.pyplot.savefig('manning2.eps')

matplotlib.pyplot.figure(3)
for i in range(numheights):
    label='plot' +str(heights[i*L2])
    matplotlib.pyplot.plot(numpy.sqrt(Slopes[i*L2:(i+1)*L2])*math.pow(heights[i*L2],2.0/3.0),Velocities[i*L2:(i+1)*L2], 'r',label=label, linewidth=1)
    matplotlib.pyplot.xlabel('R^2/3 * S^1/2')
    matplotlib.pyplot.ylabel('V')
#matplotlib.pyplot.legend('plot1', loc='lower right')
matplotlib.pyplot.savefig('manning3.eps')

matplotlib.pyplot.figure(4)
for i in range(numheights):
    label='plot' +str(heights[i*L2])
    matplotlib.pyplot.plot(Slopes[i*L2:(i+1)*L2],Velocities[i*L2:(i+1)*L2], 'r',label=label, linewidth=1)
    label='plot' +str(heights[i*L2]*50)
    matplotlib.pyplot.plot(Slopes[i*L2:(i+1)*L2],Slopes[i*L2:(i+1)*L2]*math.pow(heights[i*L2],fit2[0][0][0])*math.exp(fit2[0][1][0]), 'b',label=label, linewidth=1)
    matplotlib.pyplot.xlabel('-grad(H)')
    matplotlib.pyplot.ylabel('V')
    a= 'a =' +str(math.exp(fit2[0][1][0]))+'*R^' + str(fit2[0][0][0])
    matplotlib.pyplot.title(a)
#matplotlib.pyplot.legend('plot1', loc='lower right')
matplotlib.pyplot.savefig('manning4.eps')

matplotlib.pyplot.figure(5)
for i in range(numheights):
    label='plot' +str(heights[i*L2])
    matplotlib.pyplot.plot(Slopes[i*L2:(i+1)*L2],Velocities[i*L2:(i+1)*L2], 'r',label=label, linewidth=1)
    label='plot' +str(heights[i*L2]*50)
    matplotlib.pyplot.plot(Slopes[i*L2:(i+1)*L2],numpy.sqrt(Slopes[i*L2:(i+1)*L2])*math.pow(heights[i*L2],fit3[0][0][0])*math.exp(fit3[0][1][0]), 'b',label=label, linewidth=1)
    matplotlib.pyplot.xlabel('-grad(H)')
    matplotlib.pyplot.ylabel('V')
#matplotlib.pyplot.legend('plot1', loc='lower right')
matplotlib.pyplot.savefig('manning5.eps')

matplotlib.pyplot.figure(6)
for i in range(numheights):
    label='plot' +str(heights[i*L2])
    matplotlib.pyplot.plot(Slopes[i*L2:(i+1)*L2],abs(Velocities[i*L2:(i+1)*L2]-Slopes[i*L2:(i+1)*L2]*math.pow(heights[i*L2],fit2[0][0][0])*math.exp(fit2[0][1][0]))/Velocities[i*L2:(i+1)*L2], 'r',label=label, linewidth=1)
##     label='plot' +str(heights[i*L2]*50)
##     matplotlib.pyplot.plot(Slopes[i*L2:(i+1)*L2],Slopes[i*L2:(i+1)*L2]*math.pow(heights[i*L2],fit2[0][0][0])*math.exp(fit2[0][1][0]), 'b',label=label, linewidth=1)
    matplotlib.pyplot.xlabel('-grad(H)')
    matplotlib.pyplot.ylabel('relative error')
    a= 'a =' +str(math.exp(fit2[0][1][0]))+'*R^' + str(fit2[0][0][0])
    matplotlib.pyplot.title(a)
#matplotlib.pyplot.legend('plot1', loc='lower right')
matplotlib.pyplot.savefig('manning6.eps')

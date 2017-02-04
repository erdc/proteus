#from scipy import *
from pylab import *
from numpy import *
import collections as cll

filename='pointGauge_pressure.txt'
headerline = 0


# Read header info 

fid = open(filename,'r')
prdata=loadtxt(fid,skiprows=headerline)
fid.seek(0)
D = fid.readlines()
header = D[:headerline]
n=len(D)

a = array(zeros((int(n),5)))
step = zeros(int(n),float)
for i in range (int(n)):
 a[i,:]=D[i].split()
 step[i]=i

a = array(a,dtype=float32)

y = zeros(int(n),float)

y[:] = a[:,2]






#j=str(2)

fig1 = figure(1)
fig1.add_subplot(1,1,1)
plot(step[:],y[:],"k",lw=2)

##xlim(80,90)
##ylim(-0.05,0.05)

grid()
#title('Surface elevation at $x=$' + str(probex[int(j)])+'m')
title('Point gauge pressure evolution')
xlabel('Step')
ylabel('Pressure')

show()
savefig('plottrial.pdf')
savefig('plottrial.png',dpi=100)

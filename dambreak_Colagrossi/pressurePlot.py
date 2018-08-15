from numpy import *
from scipy import *
from pylab import *
import collections as cll
import csv

# Put relative path below
filename='pressureGauge.csv'

# Reading file
with open (filename, 'rb') as csvfile:
    data=csv.reader(csvfile, delimiter=",")
    a=[]
    time=[]
    probes=[]
    nRows=0

    for row in data:
# Time steps
        if nRows!=0:
            time.append(float(row[0]))

# Probes location ######################
        if nRows==0:                   #
            for i in row:              #
                if i!= '      time':   #
                    i=float(i[14:24])  #
                    probes.append(i)   #
########################################

        row2=[]
        for j in row:
            if j!= '      time' and nRows>0.:
                j=float(j)
                row2.append(j)
        a.append(row2)
        nRows+=1

#####################################################################################

# Choose which probes to plot  
    x = 1 
    pressure2=np.zeros(nRows-1, dtype=float)
    pressure2A=np.zeros(nRows-1, dtype=float)
    timeA=np.zeros(nRows-1, dtype=float)
    for k in range(1,nRows):
        pressure2[k-1]=(float(a[k][x]))    
        timeA[k-1]=time[k-1]*(9.81/0.6)**0.5
    pressure2A=pressure2/(998.2*9.81*0.6)
# Plot pressure in time
    import matplotlib.pyplot as plt
    plt.plot(timeA,pressure2A)
    plt.xlabel('time step [adim]')    
    plt.ylabel('pressure [adim]')
    plt.suptitle('Nondimensional pressure against time')
    plt.ylim((0,1.4))
    plt.xlim((0,12))
    plt.grid(True)
    plt.show()
    savefig('pressureAdim_in_time.png')

#####################################################################################

# Print an output file to validate the results
    maxPressureCal = max(pressure2A)
    maxPressureRef = 0.876481416000
    err = 100*abs(maxPressureRef-maxPressureCal)/maxPressureRef
    val = open('validation.txt', 'w')
    val.write('Only for gauges taken at (x,y)=(3.22,0.12)'+'\n')
    val.write('Maximum nondimensionalized pressure:'+'\n')
    val.write('Reference'+'\t'+'Simulation'+'\t'+'Error'+'\n')
    val.write(str(maxPressureRef)+'\t'+str(maxPressureCal)+'\t'+str(err))
    val.close()

#####################################################################################

# Print an output file
    info = open('probes.txt','w')
    string1=str(probes)
    string2=string1.replace('[',' ')
    string3=string2.replace(']',' ')   
    string4=string3.replace(',','\t') 
    info.write(str('x')+'\t'+string4+'\n')
    for j in range(1,nRows):
        string5=str(a[j])
        string6=string5.replace('[',' ')
        string7=string6.replace(']',' ')   
        string8=string7.replace(',','\t')   
        info.write(string8+'\n')
    info.close()





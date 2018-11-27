import numpy as np
import collections as cll
import csv
import os
import matplotlib.pyplot as plt
from proteus import WaveTools as wt
from AnalysisTools import signalFilter,zeroCrossing,reflStat

#####################################################################################
folders = ['A1','A2','A3','A4','A5','A6']
band = [2.,3.,4.,8.,8.,8.]
period = [2.,3.,4.,2.,3.,4.,]
height = [0.05,0.05,0.05,0.15,0.15,0.15]
wavelength = [3.72,5.84,7.93,3.9,6.21,8.59]
Ycoeff = [(0.04154750,0.00557540,0.00066333,0.00008224,0.00001091,0.00000154,0.00000023,0.00000007),
          (0.02547490,0.00658928,0.00138630,0.00027829,0.00005622,0.00001167,0.00000261,0.00000109),
          (0.01740417,0.00688280,0.00219367,0.00064631,0.00018577,0.00005368,0.00001682,0.00000909),
          (0.10563897,0.03899903,0.01306615,0.00457401,0.00172175,0.00070315,0.00033483,0.00024142),
          (0.05603136,0.03165234,0.01539034,0.00719950,0.00342378,0.00173118,0.00101564,0.00081894),
    (0.03386280,0.02371217,0.01455201,0.00842901,0.00485473,0.00292776,0.00199970,0.00171312)]
Bcoeff = [(0.05392302,0.00359353,0.00020733,0.00000737,-0.00000016,-0.00000005,0.00000000,0.00000000),
          (0.03973519,0.00527646,0.00071397,0.0000908,0.00001009,0.00000085,0.00000002,-0.00000001),
          (0.03117069,0.00629140,0.00134427,0.00028562,0.00005866,0.00001132,0.00000201,0.00000029),
          (0.13540388,0.02480804,0.00426381,0.00055395,0.00002809,-0.00000926,-0.00000291,-0.00000030),
          (0.08669061,0.02485208,0.00774216,0.00233218,0.00064168,0.00014869,0.00002212,-0.00000001),
          (0.06029399,0.02131785,0.00866053,0.00357527,0.00145114,0.00057285,0.00023282,0.00006823)]

ifo=-1
## Reading probes into the file
for folder in folders:
    ifo+=1
    os.chdir(folder)
    file_vof = 'column_gauges.csv'
    wave = wt.MonochromaticWaves(period=period[ifo], waveHeight=height[ifo], mwl=0.4,
                                 depth=0.4,g=np.array([0., -9.81, 0.]),
                                 waveDir=np.array([1.,0.,0.]), wavelength=wavelength[ifo],
                                 waveType="Fenton",Ycoeff=np.array(Ycoeff[ifo]),
                                 Bcoeff=np.array(Bcoeff[ifo]),Nf=len(Ycoeff[ifo]),fast=True)

    def readProbeFile(filename):
        with open (filename, 'rb') as csvfile:
            data=np.loadtxt(csvfile, delimiter=",",skiprows=1)
            time=data[:,0]
            data = data[:,1:]
            csvfile.seek(0)
            header = csvfile.readline()
            header = header.replace("time","")
            header = header.replace("[","")
            header = header.replace("]","")
            header = header.replace(","," ")
            header = header.split()
            probeType = []
            probex = []
            probey = []
            probez = []        
            for ii in range(0,len(header),4):
                probeType.append(header[ii])
                probex.append(float(header[ii+1]))
                probey.append(float(header[ii+2]))
                probez.append(float(header[ii+3]))
            probeCoord = zip(np.array(probex),np.array(probey),np.array(probez))
            datalist = [probeType,probeCoord,time,data]
            return datalist

    data_vof = readProbeFile(file_vof)
    
    #####################################################################################

    # Exctracting probes
    time = data_vof[2]
    vof = data_vof[3]
    eta_num = []
    tank_dim = [20.62,0.7]
    waterLevel = 0.4
    i_mid = 0#len(vof[0])/2-1
    for i in range(0, len(vof)):
        eta_num.append(tank_dim[1]-vof[i][i_mid]-waterLevel)
    eta_num = np.array(eta_num)
    fp = 1./period[ifo]
    minf = 0.75*fp
    maxf = 1.1*band[ifo]*fp
    time_int = np.linspace(time[0],time[-1],len(time))
    eta_num = np.interp(time_int,time,eta_num)
    eta_num = signalFilter(time,eta_num,minf, maxf, maxf, minf)
    eta_num = np.interp(time,time_int,eta_num)

    # Theoretical eta
    x = np.array(data_vof[1][2*i_mid])
    eta_th = []
    for i in range(0,len(time)):
        eta_th.append(wave.eta(x,time[i]))

    #####################################################################################

    # Plotting the probes
    plt.figure(ifo+1)
    plt.plot(time, eta_num, 'b', label='numerical')
    plt.plot(time, eta_th, 'r--', label='theoretical')
    plt.legend(loc='upper right')
    plt.xlabel('time [sec]')
    plt.ylabel('eta [m]')
    plt.xlim((76.,100.))
    plt.ylim((-0.2,0.2))
    plt.suptitle('Surface elevation against time in the middle of the tank.')
    plt.grid()
    plt.savefig('eta_NLW.png')
    
    #####################################################################################

    # Validation of the result
    S = 0.
    c = 0.
    istart = np.where(time>=60.)[0][0]
    iend = np.where(time>=120.)[0][0]
    for i in range(istart,iend):
        c = c + 1.
        S = S + (eta_th[i]-eta_num[i])**2
    err = np.sqrt(S/c)
    err = 100*err/(height[ifo]+0.4)
    val = open('validation_eta_NLW.txt', 'w')
    val.write('Eta in the middle of the tank.'+'\n')
    val.write('Gauges taken between 6s and 18s'+'\n')
    val.write('Average error (%) between the theoretical function and the simulation:'+'\n')
    val.write(str(err))
    val.close()

    #####################################################################################

    # Reflection
    dataW = readProbeFile('column_gauges.csv')
    time = dataW[2]
    L = wavelength[ifo]
    Nwaves = 10.
    T = period[ifo]
    Tend = time[-1]
    Tstart = Tend-Nwaves*T
    i_mid = len(dataW[3][0])/2-1
    time_int = np.linspace(time[0],Tend,len(time))
    data1 = np.zeros((len(time),len(dataW[3][0])),"d")
    dx_array = 0.25
    Narray = int(round(L/6./dx_array))
    data = np.zeros((len(data1),3))
    zc = []
    minf = 0.8/period[ifo]
    maxf = 1.2/period[ifo]
    for ii in range(0,3):
        data1[:,i_mid+ii*Narray] = np.interp(time_int,time,dataW[3][:,i_mid+ii*Narray])
        data[:,ii] = signalFilter(time,data1[:,i_mid+ii*Narray],minf, maxf, 1.1*maxf, 0.9*minf)
        zc.append(zeroCrossing(time,data[:,ii]))
    H1 = zc[0][1]
    H2 = zc[1][1]
    H3 = zc[2][1]
    HH = reflStat(H1,H2,H3,Narray*dx_array,L)[0]
    RR = reflStat(H1,H2,H3,Narray*dx_array,L)[2]
    print "RR = ", RR
    os.chdir("../")

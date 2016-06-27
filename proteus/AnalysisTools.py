import numpy as np
import WaveTools as WT
import math
       

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
        header =  header.replace(","," ")
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

def signalFilter(time,data,minfreq,maxfreq,costapCut = False):
    dt = (time[-1]-time[0])/(len(time)-1)
    doInterp = False
    data1 = np.zeros(data.shape,)
    

    for i in range(1,len(time)):
        dt_temp = time[i]-time[i-1]
        if dt_temp!=dt:
            doInterp = True
    if(doInterp):
#        print "Interpolating series"
        time_lin = np.linspace(time[0],time[-1],len(time))
        try:
            for ii in range(nprobes):
                data1[:,ii] = np.interp(time_lin,time,data[:,ii])
            time = time_lin
            data = data1
            nprobes = len(data[0,:])

        except:
            data1 = np.interp(time_lin,time,data[:])
            time = time_lin
            data = data1
            nprobes = -1
    nfft = len(time)
    dt = (time[-1]-time[0])/(len(time)-1)
    freq = np.fft.fftfreq(nfft,dt)   
    i1 = np.where(freq > maxfreq)[0]
    i3 = np.where(freq < -maxfreq)[0]
    i2a = []
    i2b = []
    for jj in range(1,len(freq)):
        if(freq[jj] < minfreq and freq[jj]>0):
            i2a.append(jj)
        if(freq[jj] > -minfreq and freq[jj]<0):
            i2b.append(jj)
    band1 = min(i1)-max(i2a) 
    band2 = min(i2b) - max(i3) 

    
    del data1
    data1 = np.zeros(data.shape)
    fft_x = None
    for ii in range(max(nprobes,1)):
        if(nprobes == -1):
            fft_x = np.fft.fft(data[:],nfft)
            fft_x[i1] = 0.
            fft_x[i2a] = 0.
            fft_x[i2b] = 0.
            fft_x[i3] = 0.
            if(costapCut):
                fft_x[max(i2a):min(i1)] *= WT.costap(band1,0.1)
                fft_x[max(i3):min(i2b)] *= WT.costap(band2,0.1)
            data1[:] = np.fft.ifft(fft_x)
            
        else:
            fft_x = np.fft.fft(data[:,ii],nfft)
            fft_x[i1,ii] = 0.
            fft_x[i2a,ii] = 0. 
            fft_x[i2b,ii] = 0. 
            fft_x[i3,ii] = 0. 
            if(costapCut):
                fft_x[max(i2a):min(i1),ii] *= WT.costap(band1,0.1)
                fft_x[max(i3):min(i2b),ii] *= WT.costap(band2,0.1)
            data1[:,ii] = np.fft.ifft(fft_x)


    return data1


def zeroCrossing(time,data,Tstart=0.,Tend=1e300,mode="mean",up=True,filt=True,minfreq=0.,maxfreq=1e300,costapCut=True):
    irange = np.where(time[np.where(time<Tend)[0]]>Tstart)
    time = time[irange]
    data = data[irange]
#    print time
    if(filt):
        data_filt = signalFilter(time,data,minfreq,maxfreq,costapCut)    
    trend = np.mean(data_filt)
    data_filt = data_filt - trend

    data_temp=np.zeros(data_filt.shape,)
    data_temp[0:-1] = data_filt[1:]
    zc = data_temp*data_filt
    zcPoints = np.where(zc<0)[0]
    if(up):
        if(data_filt[0]<0):
            zcPoints = zcPoints[::2]
        if(data_filt[0]>0):
            zcPoints = zcPoints[1::2]
    else:
        if(data_filt[0]<0):
            zcPoints = zcPoints[1::2]
        if(data_filt[0]>0):
            zcPoints = zcPoints[::2]

    zCH = []
    period=[]
    for zi in range(1,len(zcPoints)):
        i1 = zcPoints[zi-1]
        i2 = zcPoints[zi]
        zCH.append(max(data[i1:i2])-min(data[i1:i2]))
        period.append(time[i2]-time[i1])
    zCH = np.array(zCH)
    period = np.array(period)

    args = zCH.argsort()
    height = zCH[args]
    period = period[args]
#    print height, period

    if mode == "mean":
        height = np.mean(height)
        period = np.mean(period)
    elif mode >= 1:
        ii = len(height) - mode
        height = np.mean(height[ii:])
        period = np.mean(period[ii:])
    elif  mode < 1:
        ii = int(round(len(height*(1. - mode))))
        height = np.mean(height[ii:])
        period = np.mean(period[ii:])
    else:
        print "mode must be either an integer > 1, a float < 1 or a 'mean'"
    return [period,height]
                      

def pressureToHeight(data,Z,depth,wavelength,rho,g):
    k = 2*math.pi/wavelength
    Kp = rho*g*cosh(k*(depth+Z))/cosh(k*depth)
    return data/Kp


def ReflStat(H1,H2,H3,dx,wavelenght):
    D = 2*math.pi*dx/wavelegth
    Amp =np.array([H1/2.,H2/2.,H3/2.])
    A1 = Amp[j]*Amp[j]
    A2 = Amp[j+1]*Amp[j+1]
    A3 = Amp[j+2]*Amp[j+2]
    Lamda = (A1 + A3 - 2.*A2*cos(2*D))/(4.*sin(D)*sin(D))
    Gamma = 0.5*sqrt(
        ((2*A2-A1-A3)/(2.*sin(D)*sin(D)))**2+((A1-A3)/sin(2*D))**2)
    
    Hi = sqrt(Lamda + Gamma) + sqrt(Lamda - Gamma)
    Hr = sqrt(Lamda + Gamma) - sqrt(Lamda - Gamma)
    Rf = Hr/(Hi+1e-15)
    return [Hi,Hr,Rf]

#    i3 = np.where(freq[np.where(freq<0)[0]]
#    i4 = np.where(freq[np.where(freq<0)[0]]

    

#    WT.costap(len(fft_x),0.1)

        



#    return [time_lin,data1]
        


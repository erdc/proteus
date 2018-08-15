import numpy as np

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
    
def costapValue(f,tail):
    return 0.5*(1.-np.cos(np.pi*float(tail-f)/float(tail)))
    
def signalFilter(time,data,minfreq,maxfreq, cutoffMax, cutoffMin):
    dt = (time[-1]-time[0])/(len(time)-1)
    doInterp = False
    nfft = len(time)
    dt = (time[-1]-time[0])/(len(time)-1)
    freq = np.fft.fftfreq(nfft,dt)
    fft_x = np.fft.fft(data[:],nfft)
    ii = -1
    tailMax = cutoffMax - maxfreq
    tailMin = minfreq - cutoffMin

    if(tailMax < 0 or tailMin < 0):
        print "cutoffMax is less than maxfreq or cutoffMin larger than minfreq, this should not be the case"
    
    for ff in freq:
        ii+=1
        if ff > maxfreq:
            if ff - maxfreq < tailMax:
               f = ff -maxfreq
               fft_x[ii] = costapValue(f,tailMax)*fft_x[ii] 
            else:
                fft_x[ii] = 0.
        if ff < -maxfreq:       
            if -ff - maxfreq < tailMax:
               f = -ff -maxfreq
               fft_x[ii] = costapValue(f,tailMax)*fft_x[ii] 
            else:
               fft_x[ii] = 0.

        if (ff < minfreq and ff > -minfreq and ff!=0.):
            fd = abs(ff)
            if minfreq - fd < tailMin:
               f =  minfreq - fd
               fft_x[ii] = costapValue(f,tailMin)*fft_x[ii] 
            else:
               fft_x[ii] = 0.
    data1 = np.zeros(data.shape)    
    data1[:] = np.fft.ifft(fft_x)
    return data1

def zeroCrossing(time,data,up=True):
    trend = np.mean(data)
    data = data - trend
    data_temp=np.zeros(data.shape,)
    data_temp[0:-1] = data[1:]
    zc = data_temp*data
    zcPoints = np.where(zc<0)[0]
    if(up):
        if(data[0]<0):
            zcPoints = zcPoints[::2]
        if(data[0]>0):
            zcPoints = zcPoints[1::2]
    else:
        if(data[0]<0):
            zcPoints = zcPoints[1::2]
        if(data[0]>0):
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
    height = np.sort(zCH)
    #ii = len(height) - float(len(height))/float(mode)
    height = np.mean(height)
    period = np.mean(period)
    return [period, height]

def reflStat(H1,H2,H3,dx,wavelength):
    D = 2*np.pi*dx/wavelength
    Amp =np.array([H1/2.,H2/2.,H3/2.])
    A1 = Amp[0]*Amp[0]
    A2 = Amp[1]*Amp[1]
    A3 = Amp[2]*Amp[2]
    Lamda = (A1 + A3 - 2.*A2*np.cos(2*D))/(4.*np.sin(D)*np.sin(D))
    Gamma = 0.5*np.sqrt(
        ((2*A2-A1-A3)/(2.*np.sin(D)*np.sin(D)))**2+((A1-A3)/np.sin(2*D))**2)
    
    Hi = np.sqrt(Lamda + Gamma) + np.sqrt(Lamda - Gamma)
    Hr = np.sqrt(Lamda + Gamma) - np.sqrt(Lamda - Gamma)
    Rf = Hr/(Hi+1e-15)
    return [Hi,Hr,Rf]

#!/usr/bin/env python
from numpy import *


def  interpolate(master,slave1,slave2):
    if  master[0]<-999:
        return array([slave1[0],abs(slave1[1]),slave1[2]])

    x = master[0]

    alpha = (x-slave1[0])/(slave2[0]-slave1[0])

    y = 0.5*(abs(master[1]) + abs((1.0-alpha)*slave1[1]  + alpha*slave2[1]) )
    z = 0.5*(master[2] + (1.0-alpha)*slave1[2] + alpha*slave2[2] )

    return array([x,y,z])


def  mergeWaterline(size, step):



        # Loop over procs to load waterline
    for proc in range(size):
        data = load('waterline.' + str(proc)  + '.' + str(step)  + '.npy')

        if proc == 0:
            waterline = data
        else:
            waterline = concatenate((waterline,data))

    # Sort waterline
    waterline.view('d8,d8,d8').sort(order=['f0'], axis=0)






    # Write to file in ASCII mode
    wlfile  = open("waterline." + str(step) +".dat",'w')

    pos = array([-9999,-9999,-9999])
    neg = array([-9999,-9999,-9999])
    for wl in waterline:
        if wl[1]<0:
            cur = interpolate(pos,wl,neg)
            neg = wl
        else:
            cur = interpolate(neg,wl,pos)
            pos = wl

        wlfile.write('%12.5E %12.5E %12.5E\n'% (cur[0],cur[1],cur[2]))
    wlfile.close()


if __name__ == '__main__':
    from optparse import OptionParser

    usage = ""
    parser = OptionParser(usage=usage)

    parser.add_option("-n","--size",
                      help="number of processors for run",
                      action="store",
                      type="int",
                      dest="size",
                      default=1)

    parser.add_option("-s","--stride",
                      help="stride for solution output",
                      action="store",
                      type="int",
                      dest="stride",
                      default=0)

    parser.add_option("-t","--time",
                      help="finaltime",
                      action="store",
                      type="int",
                      dest="finaltime",
                      default=1000)


    (opts,args) = parser.parse_args()

    start = 0
    if opts.stride == 0:
        mergeWaterline(opts.size,opts.finaltime)
    elif  opts.stride > 0:
        for step in range(0,opts.finaltime+1,opts.stride):
            mergeWaterline(opts.size,step)

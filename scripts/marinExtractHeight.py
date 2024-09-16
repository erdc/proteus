#! /usr/bin/env pvpython
from paraview import servermanager
from paraview.simple import *
from optparse import OptionParser

usage = ""
parser = OptionParser(usage=usage)
parser.add_option("-f","--filebase",
                    help="filename",
                    action="store",
                    type="string",
                    dest="filename",
                    default="simulation")

parser.add_option("-r","--resolution",
                    help="filename",
                    action="store",
                    type="int",
                    dest="resolution",
                    default=100)

parser.add_option("-a","--accuracy",
                    help="filename",
                    action="store",
                    type="float",
                    dest="accuracy",
                    default=0.01)

(opts,args) = parser.parse_args()

if not paraview.servermanager.ActiveConnection:
    connection = paraview.servermanager.Connect()

reader = servermanager.sources.XDMFReader(FileNames=opts.filename)
reader.UpdatePipeline()
timesteps = reader.TimestepValues


xloc = [2.724, 2.228,1.732, 0.582]


lines=[]
for x in xloc:
    lines.append(Line(Point1=[x,0.5,0.0],Point2=[x,0.5,1.0],Resolution=opts.resolution))

probes=[]
for line in lines:
    probes.append(ProbeLocation(ProbeType=line,Input=reader))

point=PointSource(Center=[0.0,0.5,0.0],NumberOfPoints=1)
pprobe=ProbeLocation(ProbeType=point,Input=reader)

outfile = open("height.txt",'w')
for time in timesteps:
    print("Time =" + str(time))
    outfile.write(str(time))

    phi_old = 99.99
    height  = 0.0
    p =-1
    for probe in probes:
        p = p+1
        probe.UpdatePipeline (time)

        fp = servermanager.Fetch(probe)
        pdata= fp.GetPointData()
        for i in  range(opts.resolution+1):
            phi = pdata.GetArray("phid").GetTuple1(i)

            if (phi > 0.0) and (phi_old < 0.0):
                height = old_div((float(i-1) + (old_div(phi_old,(phi_old-phi)))),float(opts.resolution))
            phi_old=phi

        if (height > 0.0):
            phi = 111.0
            while abs(phi) > opts.accuracy:
                point.Center=[xloc[p],0.5,height]
                fp2 = servermanager.Fetch(pprobe)
                pdata2= fp2.GetPointData()
                phi  = pdata2.GetArray("phid").GetTuple1(0)
                height = height - phi
                print(time,height, phi)
        else:
            height = 0.0
        outfile.write("  " + str(height))
    outfile.write("\n")
outfile.close()

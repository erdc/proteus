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

(opts,args) = parser.parse_args()

if not paraview.servermanager.ActiveConnection:
    connection = paraview.servermanager.Connect()

reader = servermanager.sources.XDMFReader(FileNames=opts.filename)
reader.UpdatePipeline()
timesteps = reader.TimestepValues

points=[]
points.append(PointSource(Center=[2.3950,0.4745,0.020],NumberOfPoints=1))
points.append(PointSource(Center=[2.3950,0.4745,0.100],NumberOfPoints=1))
points.append(PointSource(Center=[2.4195,0.5255,0.161],NumberOfPoints=1))
points.append(PointSource(Center=[2.4995,0.5255,0.161],NumberOfPoints=1))

probes=[]
for point in points:
    probes.append(ProbeLocation(ProbeType=point,Input=reader))

outfile = open("pressure.txt",'w')
for time in timesteps:
    outfile.write(str(time))
    for probe in probes:
        probe.UpdatePipeline (time)

        fp = servermanager.Fetch(probe)
        pdata= fp.GetPointData()
        pressure = pdata.GetArray("p").GetTuple1(0)

        outfile.write("  " + str(pressure))
    outfile.write("\n")
outfile.close()

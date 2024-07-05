import sys
import os
import Gnuplot
import numpy as numpy
"""
Script to creat an animation of Osher solution dat file
"""
print("Making animations")

datFile = sys.argv[1]
totalNumberOfIndex = int(sys.argv[2])
T_start = float(sys.argv[3])
T_stop  = float(sys.argv[4])

FudgeFactor = 1e-5
i = 0

while (i <= totalNumberOfIndex):

    time  = T_start + i*( (T_stop-T_start)/totalNumberOfIndex )

    g = Gnuplot.Gnuplot(debug=1)
    g('set style data lines')
    g('set key bottom left')
    g('set yrange [-0.0:1.0]')
    g.xlabel('x')
    g.ylabel('S_e')
    PlotTitle = "S_e vs x at time: " + str(time)
    g.title(PlotTitle)
    g('set term postscript eps enhanced color solid')
    Fraction = float(time)/float(T_stop) + FudgeFactor
    outFile = 'set output \"' + str('Se') + 'plotFraction' + str(Fraction) + 'MakeAnimations.eps\"'
    g(outFile)
    Plot = 'plot \'' + str(datFile) + '\' index ' + str(i) + ' title \"\"'
    g(Plot)

    i += 1

i = 0
while(i < 1000000):
    # Counting delay
    i = i + 1

Convert = 'convert -loop 1 *MakeAnimations.eps OsherAnimation.gif'
print(Convert)
os.system(Convert)

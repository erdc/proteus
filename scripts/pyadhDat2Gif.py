## Automatically adapted for numpy.oldnumeric Apr 14, 2008 by -c

#! /usr/bin/env python

##
# \addtogroup scripts
# Script to creat an animation of the different parts in a gnuplot dat file

# Postprocess proteus runs:
# To use: [machine]$ python MakeAnimations.py runfile.cmd runfile.dat yrange
# example: python MakeAnimations.py TestAnimation.cmd TestAnimation.dat [-0.1:2.5]
#
# \file proteusDat2Gif.py
# @{
#   \ingroup scripts
#   \brief Convert gnuplot plotting commands and data to gif files
#
#

#John Chrispell, Summer 07
#modified by cek

import sys
import os
import Gnuplot
import numpy as numpy
import optparse
"""
Script to creat an animation of the different parts in a gnuplot dat file

Postprocess proteus runs:
To use: [machine]$ python MakeAnimations.py runfile.cmd runfile.dat yrange
example: python MakeAnimations.py TestAnimation.cmd TestAnimation.dat [-0.1:2.5]
"""

usage = "usage: %prog [options] runFile"
parser = optparse.OptionParser(usage=usage)
parser.add_option("-y", "--yrange",
                  help="Set the y axis limits",
                  action="store",
                  type="string",
                  dest="yrange",
                  default="[0.0:1.0]")
parser.add_option("-f", "--frames",
                  help="Store each frame as an eps file",
                  action="store_true",
                  dest="frames",
                  default=False)
(opts,args) = parser.parse_args()

if len(args) < 1:
    raise RuntimeError("No input file specified")

print("Making animations")
inFile = args[0]+'.cmd'
datFile = args[0]+'.dat'
YRange = 'set yrange ' + opts.yrange
Variables = set()

FudgeFactor = 1e-5
f = open(inFile,'r')
lines = f.readlines()

LineTotal = len(lines)
WordsLastLine = lines[LineTotal-1].split()
FinalTime = float(WordsLastLine[-1].strip('\'\"'))

i=0
WordsFirstLine = lines[i].split()
while ( str(WordsFirstLine[-1].strip('\'\"')) == 'Condition'):
    words = lines[i].split()
    Window = words[3].strip('\;')
    Variable = words[11].strip('\'\"\_0\:\;')
    g = Gnuplot.Gnuplot(debug=1)
    g('set style data linespoints')
    g('set key bottom left')
    g(YRange)
    g.xlabel('x')
    g.ylabel(Variable)
    PlotTitle = "Initial Condition"
    g.title(PlotTitle)
    g('set term postscript eps enhanced color solid')
    Fraction = 0.00000
    outFile = 'set output \"' + str(Variable) + 'plotFraction0.00000MakeAnimations.eps\"'
    g(outFile)
    Plot = 'plot \'' + str(datFile) + '\' index ' + str(i) + ' title \"\"'
    g(Plot)
    i += 1
    WordsFirstLine = lines[i].split()

WindowLastSeen = 0
AllWindows = 0
while (i < LineTotal -1):
    words = lines[i].split()
    Window = words[3].strip('\;')
    Variable = words[11].strip('\'\"\_0\:\;')
    time  = float(words[-1].strip('\'\"'))

    if (AllWindows == 0):
        print(Window)
        print(int(Window))
        print(WindowLastSeen)
        if(int(Window) >= WindowLastSeen):
            Variables.add(Variable)
        else:
            AllWindows = WindowLastSeen
        WindowLastSeen = int(Window)

    g = Gnuplot.Gnuplot(debug=1)
    g('set style data linespoints')
    g('set key bottom left')
    g(YRange)
    g.xlabel('x')
    g.ylabel(Variable)
    PlotTitle = "Time: " + str(time)
    g.title(PlotTitle)
    g('set term postscript eps enhanced color solid')
    Fraction = float(time)/float(FinalTime) + FudgeFactor
    outFile = 'set output \"' + str(Variable) + 'plotFraction' + str(Fraction) + 'MakeAnimations.eps\"'
    g(outFile)
    Plot = 'plot \'' + str(datFile) + '\' index ' + str(i) + ' title \"\"'
    g(Plot)

    i += 1

k = 0
AnimationNames = list(Variables)
print("Names ", AnimationNames)
while(k <= AllWindows):
    Convert = 'convert -loop 1 ' + str(AnimationNames[k])  + 'plot*MakeAnimations.eps ' + str(AnimationNames[k]) + str(k) + '.gif'
    os.system(Convert)
    k += 1

if not opts.frames:
    os.system('rm *MakeAnimations.eps')

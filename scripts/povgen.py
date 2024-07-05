#!/usr/bin/env python
"""
A script for generating povray frames of zero level set
"""
import argparse
import tables
import numpy as np
from proteus import Comm, Domain, Isosurface

parser = argparse.ArgumentParser()
parser.add_argument("prefix",
                    help="The prefix of the h5 files")
parser.add_argument("-L", type=float, default=[1.0, 1.0, 1.0], nargs='+',
                    help="extents of bounding box")
parser.add_argument("-x", type=float, default=[0.0, 0.0, 0.0], nargs='+',
                    help="lower left front corner")
parser.add_argument("-s", "--steps", type=int, default=0,
                    help="number of time steps to process")
args = parser.parse_args()
domain = Domain.RectangularDomain(L=args.L, x=args.x)
comm = Comm.init()
h5 = tables.open_file(args.prefix + ".h5", "r")
isosurface = Isosurface.Isosurface((('phi_t', (0.0,)),),
                                   domain,
                                   writeBoundary=False)

def runSteps():
    for i in range(args.steps):
        isosurface.attachHDF5(h5, i)
        isosurface.calculate(checkTime=False)

import cProfile
cProfile.run('runSteps()','isostats')
import pstats
p = pstats.Stats('isostats')
p.strip_dirs().sort_stats('time').print_stats()
p.strip_dirs().sort_stats('cumulative').print_stats()
h5.close()

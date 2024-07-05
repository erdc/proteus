#!/usr/bin/env python
import glob
cfiles = glob.glob('*.c')+ glob.glob('*.cpp')+glob.glob('*.h')
for f in cfiles:
    fh = open(f,'r')
    fstr = fh.read()
    fh.close()
    fstr = fstr.replace("PyArray_LONG",
                        "PyArray_INT")
    fh = open(f,'w')
    fh.write(fstr)
    fh.close()
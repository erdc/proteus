#!/usr/bin/env python
import glob
files = glob.glob('*.py')
print files
for f in files:
        fh = open(f,'r')
        fstr = fh.read()
        fh.close()
        fstr = fstr.replace("numpy.float_",
                            "\'d\'")
        fh = open(f,'w')
        fh.write(fstr)
        fh.close()

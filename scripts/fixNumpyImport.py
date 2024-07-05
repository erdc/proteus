#!/usr/bin/env python
import glob
files = glob.glob('*.py')
print(files)
for f in files:
    fh = open(f,'r')
    fstr = fh.read()
    fh.close()
    fstr = fstr.replace("import numpy as numpy",
                        "import numpy")
    fh = open(f,'w')
    fstr = fh.write(fstr)
    fh.close()

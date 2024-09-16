#!/usr/bin/env python
import glob
files = glob.glob('*.py') + glob.glob('*.cpp') + glob.glob('*.c') + glob.glob('*.h')  + glob.glob('*.pxd')
print(files)
for f in files:
    if f != "replacePYADH.py":
        fh = open(f,'r')
        fstr = fh.read()
        fh.close()
        fstr = fstr.replace("pyadh",
                            "proteus")
        fh = open(f,'w')
        fh.write(fstr)
        fh.close()

#!/usr/bin/env python
import glob
files = glob.glob('*.py')
print(files)
for f in files:
    fh = open(f,'r')
    fstr = fh.read()
    fh.close()
    fstr = fstr.replace(r",with='linespoints'",
                        "")
    fstr = fstr.replace(r"with='linespoints',",
                        "")
    fstr = fstr.replace(r",with='lines'",
                        "")
    fstr = fstr.replace(r"with='lines',",
                        "")
    fstr = fstr.replace(r",with='points'",
                        "")
    fstr = fstr.replace(r"with='points'",
                        "")
    fstr = fstr.replace(r",with='vectors'",
                        "")
    fstr = fstr.replace(r"with='vectors'",
                        "")
    fh = open(f,'w')
    fh.write(fstr)
    fh.close()

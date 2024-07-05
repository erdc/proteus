#!/usr/bin/env python
import glob
files = glob.glob('*.c') + glob.glob('*.cpp') + glob.glob('*.h')
print(files)
for f in files:
    fh = open(f,'r')
    fstr = fh.read()
    fh.close()
    fstr = fstr.replace("#include \"numpy/oldnumeric.h\"",
                        "#include \"numpy/arrayobject.h\"")
    fh = open(f,'w')
    fstr = fh.write(fstr)
    fh.close()

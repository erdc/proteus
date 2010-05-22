#!/usr/bin/env python
import sys
lib_dirs=[]
libs=[]
include_dirs=[]
for s in sys.argv:
    if '-L' in s:
        lib_dirs.append(s[2:])
    if '-l' in s:
        libs.append(s[2:])
    if '-I' in s:
        include_dirs.append(s[2:])
f = open('config.py','r')
fstr=f.read()
f.close()
fstr+="PYADH_PETSC_LIB_DIRS = " + `lib_dirs`+"\n"
fstr+="PYADH_PETSC_LIBS = "+`libs`+"\n"
fstr+="PYADH_PETSC_INCLUDE_DIRS = "+`include_dirs`+"\n"
f = open('config.py','w')
f.write(fstr)
f.close()

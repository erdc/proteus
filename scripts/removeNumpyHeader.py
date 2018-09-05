#!/usr/bin/env python
from __future__ import print_function
import glob
files = glob.glob('*.py')
print(files)
for f in files:
    fh = open(f,'r')
    a = fh.readlines()
    if '## Automatically adapted for numpy' in a[0]:
        fh.close()
        fh = open(f,'w')
        for l in a[2:]:
            fh.write(l)
    fh.close()

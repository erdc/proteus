#!/usr/bin/env python
import glob
#old_word = 'http://juanita.wes.army.mil/~cekees/pyadh-doc/images'
#new_word = 'https://adh.usace.army.mil/pyadh-images'
old_word = "testStuff.SSPRKNewton"
new_word = "SSPRKNewton"
files = glob.glob('*_n.py')
for fn in files:
    f = open(fn,'r')
    contents = f.read()
    f.close()
    new_contents = contents.replace(old_word,new_word)
    #print new_contents
    f = open(fn,'w')
    contents = f.write(new_contents)
    f.close()

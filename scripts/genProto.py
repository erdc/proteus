#!/usr/bin/env python
import os
import optparse

parser = optparse.OptionParser()
(opts,args) = parser.parse_args()
fileRoot = args[0][:-2]
protoFile = fileRoot+'.proto'
hFile = fileRoot+'.h'
os.system('gcc %s -aux-info %s' % (args[0],protoFile))
protof = open(protoFile,'r')
hf = open(hFile,'w')
hf.write('#ifndef '+fileRoot.upper()+'_H\n')
hf.write('#define '+fileRoot.upper()+'_H\n')
line = protof.readline()
while line:
    newline = ''
    if line.find(args[0]) != -1:
        start=line.find('/*')
        end=line.find('*/')
        while start != -1:
            start = line.find('/*',end+2)
            newline += line[end+2:start]
            end = line.find('*/',start+2)
        hf.write(newline+'\n')
    line = protof.readline()
hf.write('#endif\n')
hf.close()
hf = open(hFile,'r')
contents = hf.read()
hf.close()
hf = open(hFile,'w')
new_contents = contents.replace('const const','const')
hf.write(new_contents)
hf.close()
os.system('indent %s' % (fileRoot+'.h',))
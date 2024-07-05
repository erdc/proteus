from glob import *
from os import *
files = glob('*2D*_n.py')
for f in files:
    nf = f.lower()
    words = nf.split('_')
    nf=words[0]
    for sb in words[1:-2]:
        if sb != '2d':
            nf += '_'+sb
    nf += '_2d'
    nf += '_'+words[-2]+'_'+words[-1]
    #print nf
    #dstart=nf.find('2d')
    #nstart=nf.find('n.py')
    #nf = nf[:dstart]+nf[dstart+3:nstart]+'2d'+nf[nstart:]
    #print nf
    system('svn mv '+f+' '+nf)
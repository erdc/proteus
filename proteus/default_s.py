"""
The default values for simulation modules

.. todo::

   Clean up, straighten out sim flags/sim tools situation
"""
#
#Viewers
#
viewerType = 'vtk' #None, gnuplot

viewerPause = False

viewTimes = 'All'

#viewQuantities = ['u',"q:('f',0)","q:('f',1)","q:('f',2)","q:('u',0)","q:('u',1)","q:('u',2)","q:('m',1)","q:('m',2)"]
#viewQuantities = ['u',"q:('u',0)","q:('u',1)","q:('u',2)","q:('m',1)","q:('m',2)","q:('H',1)","q:('H',2)"]
#viewQuantities = ['u',"q:('m',0)"]
viewQuantities = ['u']

viewComponents = 'All'
#
#Archivers
#
archiverType = 'xdmf' #matlab, ensight, gnuplot

#
#Logging
#
logAllProcesses=False
flushBuffer=True
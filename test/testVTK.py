from mpi4py import MPI
from vtk import *

print "vtkcrm"
vtkcrm = vtkCompositeRenderManager()
#print "vtkcomm"
#vtkcomm = vtkMPICommunicator()
#print "vtkcomm_world"
#comm_world = vtkcomm.GetWorldCommunicator()
#print comm_world
print "vtkcontr"
#help(vtkMPIController)
vtkcontr = vtkMPIController()
vtkcontr.Initialize()#0,'',1)
print "SetCommunicator"
vtkcontr.SetCommunicator(comm_world)
print "SetController"
vtkcrm.SetController(vtkcontr)
print "done"

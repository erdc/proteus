from pyadh.default_so import *

pnList = [("navier_stokes_cylinder_2d_p","navier_stokes_cylinder_2d_c0p1c0p1_n"),
          ("movingMesh_cylinder_2d_p","movingMesh_cylinder_2d_c0p1_n")]
          

systemStepControllerType = Sequential_MinAdaptiveModelStep
systemStepControllerType = Sequential_FixedStep
systemStepControllerType = Sequential_MinFLCBDFModelStep


name="cylinder_mm"
needEBQ_GLOBAL  = True
needEBQ = True
nDTout = 100
DT = 3.0/float(nDTout)
tnList = [i*DT for i  in range(nDTout+1)]

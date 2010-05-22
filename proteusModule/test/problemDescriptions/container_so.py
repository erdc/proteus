"""
The split operator module for air/water flow around a moving rigid cylinder
"""
from pyadh.default_so import *
import container

pnList = [("twp_navier_stokes_container_3d_p" , #0
           "twp_navier_stokes_container_3d_n"),
          ("ls_container_3d_p" , #1
           "ls_container_3d_n"),
          ("vof_container_3d_p" , #2
           "vof_container_3d_n"),
          ("redist_container_3d_p" ,#3 
           "redist_container_3d_n"),
          ("ls_consrv_container_3d_p" ,#4 
           "ls_consrv_container_3d_n")]
pnList = [("twp_navier_stokes_container_3d_p" , #0
            "twp_navier_stokes_container_3d_n"),
           ("ls_container_3d_p" , #1
            "ls_container_3d_n")]
#,
#          ("movingMesh_cylinder_3d_p",#5
#           "movingMesh_cylinder_3d_c0p1_n")]

name = "container"


systemStepControllerType = Sequential_MinAdaptiveModelStep
#systemStepControllerType = ISO_fixed_MinAdaptiveModelStep

needEBQ_GLOBAL = False
needEBQ = False
useOneArchive = False
archiveFlag = ArchiveFlags.EVERY_USER_STEP
if container.nDTout != None:
    tnList = [0.0,container.dt_init]+[container.dt_init+i*(container.T-container.dt_init)/(float(container.nDTout)-1) for i in range(1,container.nDTout)]
else:
    tnList = [0.0,container.dt_init,container.T]

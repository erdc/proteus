"""
The split operator module for air/water flow around a moving rigid cylinder
"""
from pyadh.default_so import *
import threep_cylinder_md_3d

pnList = [("twp_navier_stokes_cylinder_md_3d_p" , #0
               "twp_navier_stokes_cylinder_md_3d_n"),
          ("ls_threep_cylinder_md_3d_p" , #1
           "ls_threep_cylinder_md_3d_n"),
          ("vof_threep_cylinder_md_3d_p" , #2
           "vof_threep_cylinder_md_3d_n"),
          ("redist_threep_cylinder_md_3d_p" ,#3 
           "redist_threep_cylinder_md_3d_n"),
          ("ls_consrv_threep_cylinder_md_3d_p" ,#4 
           "ls_consrv_threep_cylinder_md_3d_n")]
#,
#          ("movingMesh_cylinder_3d_p",#5
#           "movingMesh_cylinder_3d_c0p1_n")]

name = "threep_cylinder_md_3d"


systemStepControllerType = Sequential_MinAdaptiveModelStep

needEBQ_GLOBAL = False
needEBQ = False
useOneArchive = False
archiveFlag = ArchiveFlags.EVERY_USER_STEP
if threep_cylinder_md_3d.nDTout != None:
    tnList = [0.0,threep_cylinder_md_3d.dt_init]+[threep_cylinder_md_3d.dt_init+i*(threep_cylinder_md_3d.T-threep_cylinder_md_3d.dt_init)/(float(threep_cylinder_md_3d.nDTout)-1) for i in range(1,threep_cylinder_md_3d.nDTout)]
else:
    tnList = [0.0,threep_cylinder_md_3d.dt_init,threep_cylinder_md_3d.T]

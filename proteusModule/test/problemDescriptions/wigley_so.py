"""
The split operator module for air/water flow around a moving rigid cylinder
"""
from pyadh.default_so import *
import wigley

pnList = [("twp_navier_stokes_wigley_3d_p" , #0
           "twp_navier_stokes_wigley_3d_n"),
          ("ls_wigley_3d_p" , #1
           "ls_wigley_3d_n"),
          ("vof_wigley_3d_p" , #2
           "vof_wigley_3d_n"),
          ("redist_wigley_3d_p" ,#3 
           "redist_wigley_3d_n"),
          ("ls_consrv_wigley_3d_p" ,#4 
           "ls_consrv_wigley_3d_n")]
#pnList = [("twp_navier_stokes_wigley_3d_p" , #0
#            "twp_navier_stokes_wigley_3d_n"),
#           ("ls_wigley_3d_p" , #1
#            "ls_wigley_3d_n")]
#,
#          ("movingMesh_cylinder_3d_p",#5
#           "movingMesh_cylinder_3d_c0p1_n")]

name = "wigley"


systemStepControllerType = Sequential_MinAdaptiveModelStep
#systemStepControllerType = ISO_fixed_MinAdaptiveModelStep

needEBQ_GLOBAL = False
needEBQ = False
useOneArchive = False
archiveFlag = ArchiveFlags.EVERY_USER_STEP
if wigley.nDTout != None:
    tnList = [0.0,wigley.dt_init]+[wigley.dt_init+i*(wigley.T-wigley.dt_init)/(float(wigley.nDTout)-1) for i in range(1,wigley.nDTout)]
else:
    tnList = [0.0,wigley.dt_init,wigley.T]

"""
The split operator module for air/water flow around a moving rigid cylinder
"""
from pyadh.default_so import *
import threep_cylinder_md_2d

if threep_cylinder_md_2d.applyCorrection:
    pnList = [("curvature_threep_cylinder_md_2d_p" ,#0 
               "curvature_threep_cylinder_md_2d_n"),
              ("twp_navier_stokes_cylinder_md_2d_p" , #1
               "twp_navier_stokes_cylinder_md_2d_n"),
#              ("navier_stokes_cylinder_2d_p" , 
#               "navier_stokes_cylinder_2d_n"),
              ("ls_threep_cylinder_md_2d_p" , #2
               "ls_threep_cylinder_md_2d_n"),
              ("vof_threep_cylinder_md_2d_p" , #3
               "vof_threep_cylinder_md_2d_n"),
              ("redist_threep_cylinder_md_2d_p" ,#4 
               "redist_threep_cylinder_md_2d_n"),
              ("ls_consrv_threep_cylinder_md_2d_p" ,#5 
               "ls_consrv_threep_cylinder_md_2d_n"),
              ("movingMesh_cylinder_2d_p",#6
               "movingMesh_cylinder_2d_c0p1_n")]
elif threep_cylinder_md_2d.applyRedistancing:
    pnList = [("curvature_threep_cylinder_md_2d_p" , 
               "curvature_threep_cylinder_md_2d_n"),
              ("twp_navier_stokes_cylinder_md_2d_p" , 
               "twp_navier_stokes_cylinder_md_2d_n"),
              ("ls_threep_cylinder_md_2d_p" , 
               "ls_threep_cylinder_md_2d_n"),
              ("redist_threep_cylinder_md_2d_p" , 
               "redist_threep_cylinder_md_2d_n"),
              ("movingMesh_cylinder_2d_p",
               "movingMesh_cylinder_2d_c0p1_n")]
else:
    pnList = [("curvature_threep_cylinder_md_2d_p" , 
               "curvature_threep_cylinder_md_2d_n"),
              ("twp_navier_stokes_cylinder_md_2d_p" , 
               "twp_navier_stokes_cylinder_md_2d_n"),
              ("ls_threep_cylinder_md_2d_p" , 
               "ls_threep_cylinder_md_2d_n"),
              ("movingMesh_cylinder_2d_p",
               "movingMesh_cylinder_2d_c0p1_n")]
    
#pnList = [("twp_navier_stokes_cylinder_md_2d_p" , #1
#               "twp_navier_stokes_cylinder_md_2d_n"),
#          ("ls_threep_cylinder_md_2d_p" , 
#           "ls_threep_cylinder_md_2d_n")
#          ]

name = "threep_cylinder_md_2d"


systemStepControllerType = Sequential_MinAdaptiveModelStep

needEBQ_GLOBAL = True
needEBQ = True
useOneArchive = False
archiveFlag = ArchiveFlags.EVERY_USER_STEP
if threep_cylinder_md_2d.nDTout != None:
    tnList = [0.0,threep_cylinder_md_2d.dt_init]+[threep_cylinder_md_2d.dt_init+i*(threep_cylinder_md_2d.T-threep_cylinder_md_2d.dt_init)/(float(threep_cylinder_md_2d.nDTout)-1) for i in range(1,threep_cylinder_md_2d.nDTout)]
else:
    tnList = [0.0,threep_cylinder_md_2d.dt_init,threep_cylinder_md_2d.T]

import os

PROTEUS = os.getenv('PROTEUS',os.getenv('HOME')+'/src/proteus')
PROTEUS_PREFIX = os.getenv('PROTEUS_PREFIX',PROTEUS+os.getenv('PROTEUS_ARCH'))
PROTEUS_GRAPHICAL = os.getenv('PROTEUS_GRAPHICAL',os.getenv('PROTEUS')+'/proteusExtensions/proteusGraphical')
PROTEUS_GRAPHICAL_EXTRA_COMPILE_ARGS= ['-Wall']
PROTEUS_GRAPHICAL_EXTRA_LINK_ARGS=['-lblas','-lstdc++']

PROTEUS_GRAPHICAL_VTK_DEFINES = [('USE_VTK_5_9_OR_LATER',1)]
PROTEUS_GRAPHICAL_VTK_INCLUDE_DIR = os.getenv('PROTEUS_PREFIX')+'/include/paraview-4.0'
PROTEUS_GRAPHICAL_VTK_LIB_DIR = os.getenv('PROTEUS_PREFIX')+'/lib/paraview-4.0'
PROTEUS_GRAPHICAL_VTK_LIBS = ['vtkCommonCore-pv4.0','vtkCommonDataModel-pv4.0','vtkWrappingPython27Core-pv4.0']

PROTEUS_TRIANGLE_INCLUDE_DIR = PROTEUS_PREFIX+'/include'
PROTEUS_TRIANGLE_H = r'"triangle.h"'
PROTEUS_TRIANGLE_LIB_DIR = PROTEUS_PREFIX+'/lib'
PROTEUS_TRIANGLE_LIB ='tri'

from distutils.core import setup, Extension
import numpy
try:
    from config import *
except:
    print "Using defaultConfig.py"
    from defaultConfig import *

setup(name='proteusGraphical',
      version='0.0.4',
      description='Python tools for visualization',
      author='Chris Kees',
      author_email='cekees@gmail.com',
      url='http://proteus.usace.army.mil',
      packages = ['proteusGraphical'],
      package_dir={'proteusGraphical':'src'},
      ext_package='proteusGraphical',
      ext_modules=[Extension('cvtkviewers',
                             ['src/cvtkviewersModule.cpp','src/vtkviewers.cpp',PROTEUS+'/proteusModule/src/mesh.cpp',PROTEUS+'/proteusModule/src/meshio.cpp'],
                             define_macros=[('PROTEUS_TRIANGLE_H',PROTEUS_TRIANGLE_H)]+PROTEUS_GRAPHICAL_VTK_DEFINES,
                             include_dirs=[numpy.get_include(),PROTEUS_GRAPHICAL,PROTEUS+'/proteusModule/include',
                                           PROTEUS_GRAPHICAL_VTK_INCLUDE_DIR
                                           ]+
                                           [PROTEUS_TRIANGLE_INCLUDE_DIR],
                             library_dirs=[PROTEUS_GRAPHICAL_VTK_LIB_DIR]+[PROTEUS_TRIANGLE_LIB_DIR],
                             libraries=PROTEUS_GRAPHICAL_VTK_LIBS+[PROTEUS_TRIANGLE_LIB],
                             extra_compile_args=PROTEUS_GRAPHICAL_EXTRA_COMPILE_ARGS,
                             extra_link_args=PROTEUS_GRAPHICAL_EXTRA_LINK_ARGS)],
      requires=['numpy']
      )

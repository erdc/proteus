from distutils.core import setup, Extension
import numpy
try:
    from config import *
except:
    raise RuntimeError("Missing or invalid config.py file. See proteusConfig for examples")

from distutils import sysconfig
import sys
sys.path.append(os.getenv('PROTEUS')+'/externalPackages/petsc4py')
from conf.petscconf import Extension as PetscExtension
from conf.petscconf import build_ext as petsc_build_ext
from conf.petscconf import config, build, build_src
from conf.petscconf import test, sdist
from distutils import sysconfig
cv = sysconfig.get_config_vars()
cv["OPT"] = cv["OPT"].replace("-DNDEBUG","-DDEBUG")
cv["OPT"] = cv["OPT"].replace("-O3","-g")
cv["CFLAGS"] = cv["CFLAGS"].replace("-DNDEBUG","-DDEBUG")
cv["CFLAGS"] = cv["CFLAGS"].replace("-O3","-g")

class my_build_ext(petsc_build_ext):
    def build_configuration(self, arch_list):
        from distutils.util import split_quoted, execute
        #
        template, variables = self.get_config_data(arch_list)
        config_data = template % variables
        #
        build_lib   = self.build_lib
        dist_name   = self.distribution.get_name()
        config_file = os.path.join(build_lib, dist_name,'petsc.cfg')#cek hack
        #
        def write_file(filename, data):
            fh = open(filename, 'w')
            try: fh.write(config_data)
            finally: fh.close()
        execute(write_file, (config_file, config_data),
                msg='writing %s' % config_file, 
                verbose=self.verbose, dry_run=self.dry_run)
    def _copy_ext(self, ext):
        from copy import deepcopy
        extclass = ext.__class__
        fullname = self.get_ext_fullname(ext.name)
        modpath = str.split(fullname, '.')
        pkgpath = os.path.join('', *modpath[0:-1])
        name = modpath[-1]
        sources = list(ext.sources)
        newext = extclass(name, sources)
        newext.__dict__.update(deepcopy(ext.__dict__))
        newext.name = name
        pkgpath=''#cek hack
        return pkgpath, newext


if 'PROTEUS_PETSC_EXTRA_LINK_ARGS' in dir():
    PROTEUS_EXTRA_LINK_ARGS = PROTEUS_EXTRA_LINK_ARGS + PROTEUS_PETSC_EXTRA_COMPILE_ARGS
print "PROTEUS_EXTRA_LINK_ARGS",PROTEUS_EXTRA_LINK_ARGS
if 'PROTEUS_PETSC_EXTRA_COMPILE_ARGS' in dir():
    PROTEUS_EXTRA_COMPILE_ARGS = PROTEUS_EXTRA_COMPILE_ARGS + PROTEUS_PETSC_EXTRA_COMPILE_ARGS
print "PROTEUS_EXTRA_COMPILE_ARGS",PROTEUS_EXTRA_COMPILE_ARGS
setup(name='proteus',
      ext_package='proteus',
      package_dir={'proteus':'src'},
      package_data = {'proteus' : ['petsc.cfg'],},
      cmdclass     = {'config'     : config,
                      'build'      : build,
                      'build_src'  : build_src,
                      'build_ext'  : my_build_ext,
                      'test'       : test,
                      'sdist'      : sdist,
                      },
      ext_modules=[PetscExtension('flcbdfWrappers',
                                  ['src/flcbdfWrappersModule.cpp','src/mesh.cpp','src/meshio.cpp'],
                                  define_macros=[('PROTEUS_TRIANGLE_H',PROTEUS_TRIANGLE_H),
                                                 ('PROTEUS_SUPERLU_H',PROTEUS_SUPERLU_H),
                                                 ('CMRVEC_BOUNDS_CHECK',1),
                                                 ('MV_VECTOR_BOUNDS_CHECK',1),
                                                 ('PETSCVEC_BOUNDS_CHECK',1),
                                                 ('F77_POST_UNDERSCORE',1),
                                                 ('USE_BLAS',1)],
                                  include_dirs=['include',
                                                numpy.get_include(),
                                                PROTEUS_SUPERLU_INCLUDE_DIR,
                                                PROTEUS_TRIANGLE_INCLUDE_DIR] + \
                                      PROTEUS_DAETK_INCLUDE_DIR + \
                                      PROTEUS_PETSC_INCLUDE_DIRS + \
                                      [PROTEUS_MPI_INCLUDE_DIR],
                                  library_dirs=[PROTEUS_DAETK_LIB_DIR]+PROTEUS_PETSC_LIB_DIRS+[PROTEUS_MPI_LIB_DIR],
                                  libraries=['m',PROTEUS_DAETK_LIB]+PROTEUS_PETSC_LIBS+PROTEUS_MPI_LIBS,
                                  extra_link_args=PROTEUS_EXTRA_LINK_ARGS,
                                  extra_compile_args=PROTEUS_EXTRA_COMPILE_ARGS)],
      requires=['numpy']
      )

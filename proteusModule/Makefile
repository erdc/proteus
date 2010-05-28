include ${PROTEUS_PREFIX}/conf/variables

all: build_ext install

build_ext:
	${PROTEUS_PYTHON} setup.py build_ext -i
	${PROTEUS_PYTHON} setuppyx.py build_ext -i
	${PROTEUS_PYTHON} setupf.py build_ext -i

install:
	${PROTEUS_PYTHON} setuppyx.py install
	${PROTEUS_PYTHON} setupf.py install
	${PROTEUS_PYTHON} setup.py install

petsc-config:
	${PROTEUS_PYTHON} ./setPetscLib.py \'${PETSC_LIB} ${PETSC_INCLUDE}\'

docs:
	doxygen proteus-doc.conf

clean:
	${PROTEUS_PYTHON} setup.py clean

cleaner: clean
	/bin/rm -rf build proteus/*.pyc proteus/*.so proteus/*.a

mpkg:
	bdist_mpkg --readme=README --license=COPYING --welcome=INTRODUCTION

source:
	${PROTEUS_PYTHON} setup.py sdist

binary:
	${PROTEUS_PYTHON} setup.py bdist_dumb

rpm:
	${PROTEUS_PYTHON} setup.py bdist_rpm

wininst:
	${PROTEUS_PYTHON} setup.py bdist_wininst

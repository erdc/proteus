all: install

#grab environment variables from shell if they are not set

PROTEUS ?= $(shell pwd)
PROTEUS_ARCH ?= $(shell uname -s)
PROTEUS_PREFIX ?= "${PROTEUS}/${PROTEUS_ARCH}"
PROTEUS_PYTHON = ${PROTEUS_PREFIX}/bin/python

install: ${PROTEUS_PREFIX} config.py
	PROTEUS_PREFIX=${PROTEUS_PREFIX} ${PROTEUS_PYTHON} setuppyx.py install
	@echo "************************"
	@echo "done installing cython extension modules"
	@echo "************************"
	PROTEUS_PREFIX=${PROTEUS_PREFIX} ${PROTEUS_PYTHON} setupf.py install
	@echo "************************"
	@echo "done installing f2py extension modules"
	@echo "************************"
	PROTEUS_PREFIX=${PROTEUS_PREFIX} ${PROTEUS_PYTHON} setuppetsc.py build --petsc-dir=${PROTEUS_PREFIX} --petsc-arch='' install
	@echo "************************"
	@echo "done insalling petsc-based extension modules"
	@echo "************************"
	PROTEUS_PREFIX=${PROTEUS_PREFIX} ${PROTEUS_PYTHON} setup.py install
	@echo "************************"
	@echo "done installing standard extension modules"
	@echo "************************"

config.py:
	cp proteusConfig/config.py.${PROTEUS_ARCH} config.py
	@echo "config.py is now proteusConfig/config.py.${PROTEUS_ARCH}"

docs:
	doxygen proteus-doc.conf

clean:
	-PROTEUS_PREFIX=${PROTEUS_PREFIX} ${PROTEUS_PYTHON} setuppyx.py clean
	-PROTEUS_PREFIX=${PROTEUS_PREFIX} ${PROTEUS_PYTHON} setupf.py clean
	-PROTEUS_PREFIX=${PROTEUS_PREFIX} ${PROTEUS_PYTHON} setuppetsc.py clean
	-PROTEUS_PREFIX=${PROTEUS_PREFIX} ${PROTEUS_PYTHON} setup.py clean

distclean: clean
	-rm -f config.py configure.done stack.done
	-rm -rf ${PROTEUS_PREFIX}
	-rm -rf build src/*.pyc src/*.so src/*.a

${PROTEUS_PREFIX}: stack hashdist
	cp stack/examples/proteus.${PROTEUS_ARCH}.yaml stack/default.yaml
	cd stack && ${PROTEUS}/hashdist/bin/hit develop -k error -f ${PROTEUS_PREFIX}
	@echo "Stack complete, test with: make check"
	@echo "or: make parallel_check"
	@echo "Please ensure that the following is prepended  to your path"
	@echo "${PROTEUS_PREFIX}/bin"

hashdist:
	@echo "No hashdist found.  Cloning hashdist from GitHub"
	git clone https://github.com/hashdist/hashdist.git

stack:
	@echo "No stack found.  Cloning stack from GitHub"
	git clone -b proteus https://github.com/hashdist/hashstack2.git stack

check:
	@echo "************************"
	@echo "SANITY ENVIRONMENT CHECK"
	@echo PROTEUS: ${PROTEUS}
	@echo PROTEUS_ARCH: ${PROTEUS_ARCH}
	@echo PROTEUS_PREFIX: ${PROTEUS_PREFIX}
	@echo "************************"
	@echo "Hello world Check!"
	${PROTEUS_PYTHON} -c "print 'hello world'"
	@echo "************************"
	@echo "Proteus Partition Test"
	${PROTEUS_PYTHON} test/test_meshParitionFromTetgenFiles.py
	@echo "************************"

parallel_check:
	@echo "************************"
	@echo "Parallel Proteus Partition Test"
	mpirun -np 4 ${PROTEUS_PYTHON} test/test_meshParitionFromTetgenFiles.py
	@echo "************************"

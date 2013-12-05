all: install

#grab environment variables from shell if they are not set

PROTEUS ?= $(shell pwd)
PROTEUS_ARCH ?= $(shell python -c "import sys; print sys.platform")
PROTEUS_PREFIX ?= ${PROTEUS}/${PROTEUS_ARCH}
PROTEUS_PYTHON ?= ${PROTEUS_PREFIX}/bin/python

ifeq ($(PROTEUS_ARCH), darwin)
PLATFORM_ENV = MACOSX_DEPLOYMENT_TARGET=$(shell sw_vers -productVersion)
endif

PROTEUS_ENV ?= PATH="${PROTEUS_PREFIX}/bin:${PATH}" \
	PYTHONPATH=${PROTEUS_PREFIX}/lib/python2.7/site-packages \
	PROTEUS_PREFIX=${PROTEUS_PREFIX} \
	${PLATFORM_ENV}

install: ${PROTEUS_PREFIX} config.py
	${PROTEUS_ENV} ${PROTEUS_PYTHON} setuppyx.py install
	@echo "************************"
	@echo "done installing cython extension modules"
	@echo "************************"
	${PROTEUS_ENV} ${PROTEUS_PYTHON} setupf.py install
	@echo "************************"
	@echo "done installing f2py extension modules"
	@echo "************************"
	${PROTEUS_ENV} ${PROTEUS_PYTHON} setuppetsc.py build --petsc-dir=${PROTEUS_PREFIX} --petsc-arch='' install
	@echo "************************"
	@echo "done installing petsc-based extension modules"
	@echo "************************"
	${PROTEUS_ENV} ${PROTEUS_PYTHON} setup.py install
	@echo "************************"
	@echo "done installing standard extension modules"
	@echo "************************"

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
	git clone https://github.com/hashdist/hashstack2.git stack

config.py:
	@echo "No config.py file found.  Running ./configure"
	./configure

check:
	@echo "************************"
	@echo "SANITY ENVIRONMENT CHECK"
	@echo PROTEUS: ${PROTEUS}
	@echo PROTEUS_ARCH: ${PROTEUS_ARCH}
	@echo PROTEUS_PREFIX: ${PROTEUS_PREFIX}
	@echo PROTEUS_ENV: ${PROTEUS_ENV}

	@echo "************************"
	@echo "Hello world Check!"
	${PROTEUS_ENV} python -c "print 'hello world'"
	@echo "************************"
	@echo "Proteus Partition Test"
	${PROTEUS_ENV} python test/test_meshParitionFromTetgenFiles.py
	@echo "************************"

	@echo "************************"
	@echo "Parallel Proteus Partition Test"
	${PROTEUS_ENV} mpirun -np 4 python test/test_meshParitionFromTetgenFiles.py
	@echo "************************"

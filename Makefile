.PHONY: all check clean distclean doc install profile proteus update

all: install

#We use environment variables from the invoking process if they have been set,
#otherwise we try our best to determine them automatically.

SHELL=/usr/bin/env bash

PROTEUS ?= $(shell pwd)
VER_CMD = git log -1 --pretty="%H"
# shell hack for now to automatically detect Garnet front-end nodes
PROTEUS_ARCH ?= $(shell [[ $$(hostname) = garnet* ]] && echo "garnet.gnu" || python -c "import sys; print sys.platform")
PROTEUS_PREFIX ?= ${PROTEUS}/${PROTEUS_ARCH}
PROTEUS_PYTHON ?= ${PROTEUS_PREFIX}/bin/python
PROTEUS_VERSION := $(shell ${VER_CMD})

ifeq ($(PROTEUS_ARCH), darwin)
PLATFORM_ENV = MACOSX_DEPLOYMENT_TARGET=$(shell sw_vers -productVersion | sed "s/\(10.[0-9]\).*/\1/")
endif

ifeq ($(PROTEUS_ARCH), cygwin)
BOOTSTRAP = cygwin_bootstrap.done
endif

# The choice for default Fortran compiler needs to be overridden on the Garnet system
ifeq ($(PROTEUS_ARCH), garnet.gnu)
FC=ftn 
F77=ftn 
F90=ftn
endif 


PROTEUS_ENV ?= PATH="${PROTEUS_PREFIX}/bin:${PATH}" \
	PYTHONPATH=${PROTEUS_PREFIX}/lib/python2.7/site-packages \
	PROTEUS_PREFIX=${PROTEUS_PREFIX} \
	PROTEUS=${PROTEUS} \
	LD_LIBRARY_PATH="${PROTEUS_PREFIX}/lib:${LD_LIBRARY_PATH}" \
	${PLATFORM_ENV}

clean:
	-PROTEUS_PREFIX=${PROTEUS_PREFIX} ${PROTEUS_PYTHON} setuppyx.py clean
	-PROTEUS_PREFIX=${PROTEUS_PREFIX} ${PROTEUS_PYTHON} setupf.py clean
	-PROTEUS_PREFIX=${PROTEUS_PREFIX} ${PROTEUS_PYTHON} setuppetsc.py clean
	-PROTEUS_PREFIX=${PROTEUS_PREFIX} ${PROTEUS_PYTHON} setup.py clean

distclean: clean
	-rm -f config.py configure.done stack.done
	-rm -rf ${PROTEUS_PREFIX}
	-rm -rf build src/*.pyc src/*.so src/*.a

update:
	@echo "Manually trying to update all repositories to most recent stable configuration"
	git fetch origin; git checkout -q origin/stable/cygwin
	@echo "Proteus repository updated to latest stable version"
	cd stack; git fetch origin; git checkout -q origin/stable/cygwin
	@echo "Stack repository updated to latest stable version"
	cd hashdist; git fetch origin; git checkout -q origin/stable/cygwin
	@echo "Hashdist repository updated to latest stable version"
	@echo "+======================================================================================================+"
	@echo "Warning, the HEAD has been detached in all repositories"
	@echo "Type: git checkout -b branch_name to save changes" 
	@echo "+======================================================================================================+"

hashdist: 
	@echo "No hashdist found.  Cloning stable hashdist from GitHub"
	git clone -b stable/cygwin https://github.com/hashdist/hashdist.git

stack: 
	@echo "No stack found.  Cloning stable stack from GitHub"
	git clone -b stable/cygwin https://github.com/hashdist/hashstack.git stack

cygwin_bootstrap.done: stack/scripts/setup_cygstack.py stack/scripts/cygstack.txt
	python stack/scripts/setup_cygstack.py stack/scripts/cygstack.txt
	touch cygwin_bootstrap.done

profile: ${PROTEUS_PREFIX}/artifact.json

# A hashstack profile will be rebuilt if Make detects any files in the stack 
# directory newer than the profile artifact file.
${PROTEUS_PREFIX}/artifact.json: stack hashdist $(shell find stack -type f) ${BOOTSTRAP}
	@echo "************************"
	@echo "Building dependencies..."
	@echo "************************"

	@echo "Please include this information in all bug reports."
	@echo "+======================================================================================================+"
	@echo "PROTEUS          : ${PROTEUS}"
	@echo "PROTEUS_ARCH     : ${PROTEUS_ARCH}"
	@echo "PROTEUS_PREFIX   : ${PROTEUS_PREFIX}"
	@echo "PROTEUS_VERSION  : ${PROTEUS_VERSION}"
	@echo "HASHDIST_VERSION : $$(cd hashdist; ${VER_CMD})"
	@echo "HASHSTACK_VERSION: $$(cd stack; ${VER_CMD})"
	@echo "+======================================================================================================+"
	@echo ""

	cp stack/examples/proteus.${PROTEUS_ARCH}.yaml stack/default.yaml
	cd stack && ${PROTEUS}/hashdist/bin/hit develop -f -k error default.yaml ${PROTEUS_PREFIX}
        # workaround hack on Cygwin for hashdist launcher to work correctly
	-cp ${PROTEUS}/${PROTEUS_ARCH}/bin/python2.7.exe.link ${PROTEUS}/${PROTEUS_ARCH}/bin/python2.7.link

	@echo "************************"
	@echo "Dependency build complete"
	@echo "************************"

proteus: ${PROTEUS_PREFIX}/bin/proteus

${PROTEUS_PREFIX}/bin/proteus ${PROTEUS_PREFIX}/bin/proteus_env.sh: profile
	@echo "************************"
	@echo "Installing proteus scripts..."
	@echo "************************"

	echo "#!/usr/bin/env bash" > ${PROTEUS_PREFIX}/bin/proteus
	echo '${PROTEUS_ENV} python "$${@:1}"' >> ${PROTEUS_PREFIX}/bin/proteus
	chmod a+x ${PROTEUS_PREFIX}/bin/proteus

	echo "#!/usr/bin/env sh" > ${PROTEUS_PREFIX}/bin/proteus_env.sh
	echo '${PROTEUS_ENV}' >> ${PROTEUS_PREFIX}/bin/proteus_env.sh
	chmod a+x ${PROTEUS_PREFIX}/bin/proteus_env.sh

	@echo "************************"
	@echo "Proteus script successfully installed"
	@echo "************************"


#config.py file should be newer than proteusConfig/config.py.$PROTEUS_ARCH
config.py: proteusConfig/config.py.${PROTEUS_ARCH}
	@echo "************************"
	@echo "Configuring..."
	@echo "************************"
	@echo "Copying proteusConfig/config.py.$PROTEUS_ARCH to ./config.py"
	@cp proteusConfig/config.py.${PROTEUS_ARCH} config.py
	@echo "************************"
	@echo "Configure complete"
	@echo "************************"


# Proteus install should be triggered by an out-of-date hashstack profile, source tree, or modified setup files.
install: profile config.py $(shell find src -type f) $(wildcard *.py) proteus
	@echo "************************"
	@echo "Installing..."
	@echo "************************"
	@echo "Please include this information in all bug reports."
	@echo "+======================================================================================================+"
	@echo "PROTEUS          : ${PROTEUS}"
	@echo "PROTEUS_ARCH     : ${PROTEUS_ARCH}"
	@echo "PROTEUS_PREFIX   : ${PROTEUS_PREFIX}"
	@echo "PROTEUS_VERSION  : ${PROTEUS_VERSION}"
	@echo "HASHDIST_VERSION : $$(cd hashdist; ${VER_CMD})"
	@echo "HASHSTACK_VERSION: $$(cd stack; ${VER_CMD})"
	@echo "+======================================================================================================+"
	@echo ""
	${PROTEUS_ENV} ${PROTEUS_PYTHON} setuppyx.py install --prefix=${PROTEUS_PREFIX}
	@echo "************************"
	@echo "done installing cython extension modules"
	@echo "************************"
	${PROTEUS_ENV} ${PROTEUS_PYTHON} setupf.py install --prefix=${PROTEUS_PREFIX}
	@echo "************************"
	@echo "done installing f2py extension modules"
	@echo "************************"
	${PROTEUS_ENV} ${PROTEUS_PYTHON} setuppetsc.py build --petsc-dir=${PROTEUS_PREFIX} --petsc-arch='' install --prefix=${PROTEUS_PREFIX}
	@echo "************************"
	@echo "done installing petsc-based extension modules"
	@echo "************************"
	${PROTEUS_ENV} ${PROTEUS_PYTHON} setup.py install --prefix=${PROTEUS_PREFIX}
	@echo "************************"
	@echo "done installing standard extension modules"
	@echo "************************"
	@echo "Installation complete"
	@echo "************************"
	@echo ""
	@echo "Proteus was built using the following configuration:"
	@echo "Please include this information in all bug reports."
	@echo "+======================================================================================================+"
	@echo "PROTEUS          : ${PROTEUS}"
	@echo "PROTEUS_ARCH     : ${PROTEUS_ARCH}"
	@echo "PROTEUS_PREFIX   : ${PROTEUS_PREFIX}"
	@echo "PROTEUS_VERSION  : ${PROTEUS_VERSION}"
	@echo "HASHDIST_VERSION : $$(cd hashdist; ${VER_CMD})"
	@echo "HASHSTACK_VERSION: $$(cd stack; ${VER_CMD})"
	@echo "+======================================================================================================+"
	@echo "PROTEUS_VERSION  : ${PROTEUS_VERSION}" > ${PROTEUS_PREFIX}/proteus.version
	@echo "HASHDIST_VERSION : $$(cd hashdist; ${VER_CMD})" > ${PROTEUS_PREFIX}/hashdist.version
	@echo "HASHSTACK_VERSION: $$(cd stack; ${VER_CMD})" > ${PROTEUS_PREFIX}/hashstack.version
	@echo ""
	@echo "You should now verify that the install succeeded by running:"
	@echo ""
	@echo "make check"
	@echo ""

check: install
	@echo "************************"
	@echo "Sanity environment check"
	@echo PROTEUS: ${PROTEUS}
	@echo PROTEUS_ARCH: ${PROTEUS_ARCH}
	@echo PROTEUS_PREFIX: ${PROTEUS_PREFIX}
	@echo PROTEUS_ENV: ${PROTEUS_ENV}

	@echo "************************"
	@echo "Hello world Check!"
	${PROTEUS_PREFIX}/bin/proteus -c "print 'hello world'"
	@echo "************************"
	@echo "Proteus Partition Test"
	${PROTEUS_PREFIX}/bin/proteus test/test_meshParitionFromTetgenFiles.py
	@echo "************************"

	@echo "************************"
	@echo "Parallel Proteus Partition Test"
	source ${PROTEUS_PREFIX}/bin/proteus_env.sh; mpirun -np 4 ${PROTEUS_PYTHON} test/test_meshParitionFromTetgenFiles.py
	@echo "************************"

doc: install
	cd doc && ${PROTEUS_ENV} make html

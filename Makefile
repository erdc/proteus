.PHONY: all check clean distclean doc install profile proteus update FORCE

all: install

#We use environment variables from the invoking process if they have been set,
#otherwise we try our best to determine them automatically.

SHELL=/usr/bin/env bash

PROTEUS ?= $(shell python -c "from __future__ import print_function; import os; print(os.path.realpath(os.getcwd()))")
VER_CMD = git log -1 --pretty="%H"
PROTEUS_INSTALL_CMD = python setup.py install -O2
PROTEUS_DEVELOP_CMD = pip install -v -e .
# shell hack for now to automatically detect Garnet front-end nodes
PROTEUS_ARCH ?= $(shell [[ $$(hostname) = garnet* ]] && echo "garnet.gnu" || python -c "import sys; print sys.platform")
PROTEUS_PREFIX ?= ${PROTEUS}/${PROTEUS_ARCH}
PROTEUS_PYTHON ?= ${PROTEUS_PREFIX}/bin/python
PROTEUS_VERSION := $(shell ${VER_CMD})
HASHDIST_DEFAULT_VERSION := $(shell cat .hashdist_default)
HASHSTACK_DEFAULT_VERSION := $(shell cat .hashstack_default)
HASHDIST_VERSION := $(shell cd hashdist; ${VER_CMD})
HASHSTACK_VERSION := $(shell cd stack; ${VER_CMD})

define show_info
	@echo "Please include this information in all bug reports."
	@echo "+======================================================================================================+"
	@echo "PROTEUS          : ${PROTEUS}"
	@echo "PROTEUS_ARCH     : ${PROTEUS_ARCH}"
	@echo "PROTEUS_PREFIX   : ${PROTEUS_PREFIX}"
	@echo "PROTEUS_VERSION  : ${PROTEUS_VERSION}"
	@echo "HASHDIST_VERSION : ${HASHDIST_VERSION}"
	@echo "HASHSTACK_VERSION: ${HASHSTACK_VERSION}"
	@echo "+======================================================================================================+"
	@echo ""
endef

define howto
	@echo "Proteus and its dependencies have been installed!"
	@echo ""
	@echo "You may want to add the installation directory to your PATH:"
	@echo "${PROTEUS_PREFIX}/bin"
	@echo ""
	@echo "On bash or zsh:"
	@echo 'export PATH=${PROTEUS_PREFIX}/bin:$${PATH}'
	@echo ""
	@echo "You can also invoke the Python interpreter directly:"
	@echo "${PROTEUS_PREFIX}/bin/python"
	@echo ""
	@echo "You should now verify that the install succeeded by running:"
	@echo "make check"
	@echo ""
endef

ifeq ($(PROTEUS_ARCH), darwin)
PLATFORM_ENV = MACOSX_DEPLOYMENT_TARGET=$(shell sw_vers -productVersion | sed -E "s/([0-9]+\.[0-9]+).*/\1/")
endif

ifeq ($(PROTEUS_ARCH), cygwin)
BOOTSTRAP = cygwin_bootstrap.done
endif

ifdef MATLAB
MATLAB_SETUP = matlab_setup.done
endif

# The choice for default Fortran compiler needs to be overridden on the Garnet system
ifeq ($(PROTEUS_ARCH), garnet.gnu)
FC=ftn 
F77=ftn 
F90=ftn
endif 

ifdef VERBOSE
HIT_FLAGS += -v
endif

ifdef DEBUG
PROTEUS_ARCH := ${PROTEUS_ARCH}-debug
endif

PROTEUS_ENV ?= PATH="${PROTEUS_PREFIX}/bin:${PATH}"

clean:
	-PROTEUS_PREFIX=${PROTEUS_PREFIX} ${PROTEUS_PYTHON} setup.py clean

distclean: clean
	-rm -f stack.done
	-rm -rf ${PROTEUS_PREFIX}
	-rm -rf build src/*.pyc proteus/*.so proteus/*.a

update:
	@echo "Manually trying to update all repositories"
	git fetch origin; git checkout -q origin/master
	@echo "Proteus repository updated to latest versions"
	cd stack; git fetch origin; git checkout -q origin/master
	@echo "Stack repository updated to latest versions"
	cd hashdist; git fetch origin; git checkout -q origin/master
	@echo "HashDist repository updated to latest versions"
	@echo "+======================================================================================================+"
	@echo "Warning, the HEAD has been detached in all repositories"
	@echo "Type: git checkout -b branch_name to save changes" 
	@echo "+======================================================================================================+"

hashdist: 
	@echo "No hashdist found.  Cloning hashdist from GitHub"
	git clone https://github.com/hashdist/hashdist.git 
	cd hashdist && git checkout ${HASHDIST_DEFAULT_VERSION}
stack: 
	@echo "No stack found.  Cloning stack from GitHub"
	git clone https://github.com/hashdist/hashstack.git stack
	cd  stack && git checkout ${HASHSTACK_DEFAULT_VERSION}

cygwin_bootstrap.done: stack/scripts/setup_cygstack.py stack/scripts/cygstack.txt
	python stack/scripts/setup_cygstack.py stack/scripts/cygstack.txt
	touch cygwin_bootstrap.done

matlab_setup.done: stack stack/default.yaml hashdist
	@echo "User requests MATLAB install"
	@echo "MATLAB environment variable set to ${MATLAB}"
	@python setupmatlab.py stack/default.yaml ${MATLAB}; if [ $$? -ne 0 ] ; then \
	echo "+======================================================================================================+"; \
	echo "Couldn't find matlab on PATH."; \
	echo "Try"; \
	echo "    MATLAB=/path/to/matlab make"; \
	echo "+======================================================================================================+"; \
	false; fi
	touch matlab_setup.done

profile: ${PROTEUS_PREFIX}/artifact.json

stack/default.yaml: stack stack/examples/proteus.${PROTEUS_ARCH}.yaml
	# workaround since mac doesn't support '-b' and '-i' breaks travis
	-cp ${PWD}/stack/default.yaml ${PWD}/stack/default.yaml.bak
	ln -sf ${PWD}/stack/examples/proteus.${PROTEUS_ARCH}.yaml ${PWD}/stack/default.yaml


# A hashstack profile will be rebuilt if Make detects any files in the stack 
# directory newer than the profile artifact file.
${PROTEUS_PREFIX}/artifact.json: stack/default.yaml stack hashdist $(shell find stack -type f) ${BOOTSTRAP} ${MATLAB_SETUP}
	@echo "************************"
	@echo "Building dependencies..."
	@echo "************************"

	$(call show_info)

	cd stack && ${PROTEUS}/hashdist/bin/hit develop ${HIT_FLAGS} -v -f -k error default.yaml ${PROTEUS_PREFIX}

	@echo "************************"
	@echo "Dependency build complete"
	@echo "************************"

versions: ${PROTEUS_PREFIX}/versions.txt
	@echo "************************"
	@echo "Installing hashdist/hashstack versions..."
	@echo "************************"

	echo ${HASHDIST_VERSION} > ${PROTEUS_PREFIX}/hashdist_version.txt
	echo ${HASHSTACK_VERSION} > ${PROTEUS_PREFIX}/hashstack_version.txt

# this always runs
${PROTEUS_PREFIX}/versions.txt: ${PROTEUS_PREFIX}/artifact.json FORCE

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

# Proteus install should be triggered by an out-of-date hashstack profile, source tree, or modified setup files.
install: profile $(shell find proteus -type f) $(wildcard *.py) proteus
	@echo "************************"
	@echo "Installing..."
	@echo "************************"
	$(call show_info)
	${PROTEUS_ENV} ${PROTEUS_INSTALL_CMD}
	@echo "************************"
	@echo "done installing standard extension modules"
	@echo "************************"
	@echo "Installation complete"
	@echo "************************"
	@echo ""
	@echo "************************"
	@echo "Installing proteus version information..."
	@echo "************************"
	@echo ${PROTEUS_VERSION} > ${PROTEUS_PREFIX}/proteus_version.txt

	@echo "Proteus was built using the following configuration:"
	$(call show_info)
	$(call howto)

develop: proteus profile 
	@echo "************************"
	@echo "Installing development version"
	@echo "************************"
	$(call show_info)
	${PROTEUS_ENV} CFLAGS="-Wall -Wstrict-prototypes -DDEBUG" ${PROTEUS_DEVELOP_CMD}
	@echo "************************"
	@echo "Development installation complete"
	@echo "************************"
	@echo ""
	@echo "************************"
	@echo "Installing proteus version information..."
	@echo "************************"
	@echo "${PWD}" > ${PROTEUS_PREFIX}/proteus_version.txt
	@echo "Proteus was built using the following configuration:"
	$(call show_info)
	$(call howto)

check:
	@echo "************************"
	@echo "Sanity environment check"
	@echo PROTEUS: ${PROTEUS}
	@echo PROTEUS_ARCH: ${PROTEUS_ARCH}
	@echo PROTEUS_PREFIX: ${PROTEUS_PREFIX}
	@echo PROTEUS_ENV: ${PROTEUS_ENV}

	@echo "************************"
	@echo "Hello world Check!"
	${PROTEUS_PREFIX}/bin/python -c "print 'hello world'"
	@echo "************************"
	@echo "Proteus Partition Test"
	source ${PROTEUS_PREFIX}/bin/proteus_env.sh; ${PROTEUS_PREFIX}/bin/python proteus/tests/ci/test_meshPartitionFromTetgenFiles.py
	@echo "************************"

	@echo "************************"
	@echo "Parallel Proteus Partition Test"
	source ${PROTEUS_PREFIX}/bin/proteus_env.sh; mpirun -np 4 ${PROTEUS_PYTHON} proteus/tests/ci/test_meshPartitionFromTetgenFiles.py
	@echo "************************"

doc:
	@echo "************************************"
	@echo "Generating documentation with Sphinx"
	@echo "Be sure to first run"
	@echo "make develop"
	@echo "or"
	@echo "make install"
	@echo "************************************"
	cd doc && ${PROTEUS_ENV} PROTEUS=${PWD} make html
	@echo "**********************************"
	@echo "Trying to open the html at"
	@echo "../proteus-website/index.html"
	@echo "**********************************"
	-sensible-browser ../proteus-website/index.html &

test:
	@echo "************************************"
	@echo "Running travis-ci test suite"
	py.test -v proteus/tests --ignore proteus/tests/POD
	@echo "Tests complete "
	@echo "************************************"

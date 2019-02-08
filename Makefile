.PHONY: all check clean distclean doc install FORCE

all: develop

#We use environment variables from the invoking process if they have been set,
#otherwise we try our best to determine them automatically.

SHELL=/usr/bin/env bash
N=1
PROTEUS ?= $(shell python -c "from __future__ import print_function; import os; print(os.path.realpath(os.getcwd()))")
VER_CMD = git log -1 --pretty="%H"
PROTEUS_BUILD_CMD = python setup.py build_ext
PROTEUS_INSTALL_CMD = python setup.py install
PROTEUS_DEVELOP_BUILD_CMD = python setup.py build_ext -i
PROTEUS_DEVELOP_CMD = pip --disable-pip-version-check install -v -e .
#
ifeq (${N}, 1)
PROTEUS_BUILD_CMD = python -c "print('Letting install handle build_ext')"
PROTEUS_DEVELOP_BUILD_CMD = python -c "print('Letting install handle build_ext')"
endif

# automatically detect hpcmp machines
PROTEUS_ARCH ?= $(shell [[ $$(hostname) = topaz* ]] && echo "topaz" || python -c "import sys; print(sys.platform)")
PROTEUS_ARCH ?= $(shell [[ $$(hostname) = onyx* ]] && echo "onyx" || python -c "import sys; print(sys.platform)")
PROTEUS_ARCH ?= $(shell [[ $$(hostname) = copper* ]] && echo "copper" || python -c "import sys; print(sys.platform)")
PROTEUS_ARCH ?= $(shell [[ $$(hostname) = excalibur* ]] && echo "excalibur" || python -c "import sys; print(sys.platform)")
PROTEUS_ARCH ?= $(shell [[ $$(hostname) = centennial* ]] && echo "centennial" || python -c "import sys; print(sys.platform)")
PROTEUS_ARCH ?= $(shell [[ $$(hostname) = thunder* ]] && echo "thunder" || python -c "import sys; print(sys.platform)")
PROTEUS_ARCH ?= $(shell [[ $$(hostname) = gordon* ]] && echo "gordon" || python -c "import sys; print(sys.platform)")
PROTEUS_ARCH ?= $(shell [[ $$(hostname) = conrad* ]] && echo "conrad" || python -c "import sys; print(sys.platform)")
PROTEUS_PREFIX ?= ${PROTEUS}/${PROTEUS_ARCH}
PROTEUS_PYTHON ?= ${PROTEUS_PREFIX}/bin/python
PROTEUS_VERSION := $(shell ${VER_CMD})
TEST_MARKER="' '"

define show_info
	@echo "Please include this information in all bug reports."
	@echo "+======================================================================================================+"
	@echo "PROTEUS        : ${PROTEUS}"
	@echo "PROTEUS_ARCH   : ${PROTEUS_ARCH}"
	@echo "PROTEUS_PREFIX : ${PROTEUS_PREFIX}"
	@echo "PROTEUS_VERSION: ${PROTEUS_VERSION}"
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
	@echo "make test"
	@echo ""
endef

ifeq ($(PROTEUS_ARCH), darwin)
PLATFORM_ENV = MACOSX_DEPLOYMENT_TARGET=$(shell sw_vers -productVersion | sed -E "s/([0-9]+\.[0-9]+).*/\1/")
endif

ifeq ($(PROTEUS_ARCH), cygwin)
BOOTSTRAP = cygwin_bootstrap.done
endif

FC ?= gfortran
F77 ?= gfortran
F90 ?= gfortran

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
	-rm -rf build proteus/*.pyc proteus/*.so proteus/*.a proteus/MeshAdaptPUMI/*.so
	-rm -rf build proteus/mprans/*.pyc proteus/mprans/*.so proteus/mprans/*.a
	-rm -rf build proteus/mbd/*.pyc proteus/mbd/*.so proteus/mbd/*.a

stack/hit/bin/hit:
	@echo "Updating stack submodule"
	git submodule update --init stack
	@echo "Updating stack/hit submodule"
	cd stack && git submodule update --init
	@echo "Adding source cache if not done already"
	-./stack/hit/bin/hit init-home
	-./stack/hit/bin/hit remote add http://levant.hrwallingford.com/hashdist_src --objects="source"

stack:
	@echo "Updating stack submodule"
	git submodule update --init stack

air-water-vv:
	@echo "Updating air-water-vv submodule"
	git submodule update --init air-water-vv

bld_cache: stack/hit/bin/hit
	@echo "Trying to add build cache for your arch"
	HASHSTACK_BLD = $(shell lsb_release -ir | python -c "import sys; rel=dict((k.split(':')[0].split()[0],k.split(':')[1].strip().replace('.','_').lower()) for k in sys.stdin.readlines()); print('{Distributor}_{Release}'.format(**rel))")
	./stack/hit/bin/hit remote add http://levant.hrwallingford.com/hashdist_${HASHSTACK_BLD} --objects="build"

cygwin_bootstrap.done: stack/scripts/setup_cygstack.py stack/scripts/cygstack.txt
	python stack/scripts/setup_cygstack.py stack/scripts/cygstack.txt
	touch cygwin_bootstrap.done

stack/default.yaml: stack/hit/bin/hit
	@echo "Linking stack/default.yaml for this arch"
	-ln -s ${PWD}/stack/examples/proteus.${PROTEUS_ARCH}.yaml ${PWD}/stack/default.yaml

# A hashstack profile will be rebuilt if Make detects any files in the stack 
# directory newer than the profile artifact file.
${PROTEUS_PREFIX}/artifact.json: stack/default.yaml $(shell find stack -type f) ${BOOTSTRAP}
	@echo "************************"
	@echo "Building dependencies..."
	@echo "************************"

	$(call show_info)

	cd stack && ${PROTEUS}/stack/hit/bin/hit develop -j ${N} ${HIT_FLAGS} -v -f -k error default.yaml ${PROTEUS_PREFIX}

	@echo "************************"
	@echo "Dependency build complete"
	@echo "************************"

${PROTEUS_PREFIX}/bin/proteus_env.sh: ${PROTEUS_PREFIX}/artifact.json
	@echo "************************"
	@echo "Installing proteus_env.sh"
	@echo "************************"
	echo "#!/usr/bin/env sh" > ${PROTEUS_PREFIX}/bin/proteus_env.sh
	echo '${PROTEUS_ENV}' >> ${PROTEUS_PREFIX}/bin/proteus_env.sh
	chmod a+x ${PROTEUS_PREFIX}/bin/proteus_env.sh

# Proteus install should be triggered by an out-of-date hashstack profile, source tree, or modified setup files.
install: ${PROTEUS_PREFIX}/bin/proteus_env.sh stack/default.yaml ${PROTEUS_PREFIX}/artifact.json
	@echo "************************"
	@echo "Installing..."
	@echo "************************"
	$(call show_info)
	${PROTEUS_ENV} ${PROTEUS_BUILD_CMD}
	${PROTEUS_ENV} ${PROTEUS_INSTALL_CMD}
	@echo "************************"
	@echo "done installing standard extension modules"
	@echo "************************"
	@echo "installing scripts"
	cd scripts && ${PROTEUS_ENV} PROTEUS_PREFIX=${PROTEUS_PREFIX} make
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

develop: ${PROTEUS_PREFIX}/bin/proteus_env.sh stack/default.yaml ${PROTEUS_PREFIX}/artifact.json
	-ln -sf ${PROTEUS}/${PROTEUS_ARCH}/lib64/* ${PROTEUS}/${PROTEUS_ARCH}/lib
	-ln -sf ${PROTEUS}/${PROTEUS_ARCH}/lib64/cmake/* ${PROTEUS}/${PROTEUS_ARCH}/lib/cmake
	@echo "************************"
	@echo "Installing development version"
	@echo "************************"
	$(call show_info)
	${PROTEUS_ENV} ${PROTEUS_DEVELOP_BUILD_CMD}
	${PROTEUS_ENV} ${PROTEUS_DEVELOP_CMD}
	@echo "************************"
	@echo "installing scripts"
	cd scripts && ${PROTEUS_ENV} PROTEUS_PREFIX=${PROTEUS_PREFIX} make
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
	${PROTEUS_PREFIX}/bin/python -c "print('hello world')"
	@echo "************************"
	@echo "Proteus Partition Test"
	source ${PROTEUS_PREFIX}/bin/proteus_env.sh; ${PROTEUS_PREFIX}/bin/python proteus/tests/ci/test_meshPartitionFromTetgenFiles.py
	@echo "************************"

	@echo "************************"
	@echo "Parallel Proteus Partition Test"
	source ${PROTEUS_PREFIX}/bin/proteus_env.sh; mpirun -np 4 ${PROTEUS_PYTHON} proteus/tests/ci/test_meshPartitionFromTetgenFiles.py
	@echo "************************"

check_simmetrix:
	@echo "SCOREC-Serial Simmetrix Error-Estimator and Adapt Test"
	${PROTEUS_ENV} ${PROTEUS_PYTHON} proteus/MeshAdaptPUMI/test/test_MeshAdaptPUMI/test_errorAndSerialAdapt/errorCheck.py
	@echo "************************"

	@echo "SCOREC-Parallel Simmetrix Error Estimator and Adapt Test"
	${PROTEUS_ENV} mpirun -np 2 ${PROTEUS_PYTHON} proteus/MeshAdaptPUMI/test/test_MeshAdaptPUMI/test_parallelAdapt/parallelAdaptCheck.py
	@echo "************************"

	@echo "SCOREC Simmetrix Isotropic Uniform Adapt Test"
	${PROTEUS_ENV} ${PROTEUS_PYTHON} proteus/MeshAdaptPUMI/test/test_MeshAdaptPUMI/test_isotropicAdapt/isotropicCheck.py
	@echo "************************"


#doc: install
doc:
	@echo "************************************"
	@echo "Generating documentation with Sphinx"
	@echo "Be sure to first run"
	@echo "make develop"
	@echo "or"
	@echo "make install"
	@echo "************************************"
	-pip install sphinx sphinx-bootstrap-theme
	cd doc && ${PROTEUS_ENV} PROTEUS=${PWD} make html
	@echo "**********************************"
	@echo "Trying to open the html at"
	@echo "../proteus-website/index.html"
	@echo "**********************************"
	-sensible-browser ../proteus-website/index.html &

test: air-water-vv check
	@echo "**************************************************"
	@echo "Running git-lfs to get regression test data files."
	-git lfs fetch
	-git lfs checkout
	@echo "If git-lfs failed to download data, then some tests will fail, and"
	@echo "you should install git-lfs or try 'make lfs', passing all tests is needed"
	@echo "**************************************************************************"
	@echo "Running basic test suite"
	-source ${PROTEUS_PREFIX}/bin/proteus_env.sh; MPLBACKEND=Agg py.test -n ${N} --dist=loadfile --forked -v proteus/tests -m ${TEST_MARKER} --ignore proteus/tests/POD --cov=proteus
	@echo "Basic tests complete "
	@echo "************************************"
	@echo "Running air-water-vv test set 1"
	-source ${PROTEUS_PREFIX}/bin/proteus_env.sh; MPLBACKEND=Agg py.test -n ${N} --dist=loadfile --forked -v air-water-vv/Tests/1st_set -m ${TEST_MARKER}
	@echo "************************************"
	@echo "Running air-water-vv test set 2"
	-source ${PROTEUS_PREFIX}/bin/proteus_env.sh; MPLBACKEND=Agg py.test -n ${N} --dist=loadfile --forked -v air-water-vv/Tests/2nd_set -m ${TEST_MARKER}

jupyter:
	@echo "************************************"
	@echo "Enabling jupyter notebook/lab/widgets"
	source ${PROTEUS_PREFIX}/bin/proteus_env.sh
	pip install configparser==3.5.0 ipyparallel==6.1.1 ipython==5.5.0 terminado==0.8.1 jupyter==1.0.0 jupyterlab==0.31.12 ipywidgets==7.1.2 ipyleaflet==0.7.1 jupyter_dashboards==0.7.0 pythreejs==0.4.1 rise==5.2.0 cesiumpy==0.3.3 bqplot==0.10.5 hide_code==0.5.2 ipympl==0.1.0 sympy==1.1.1 transforms3d==0.3.1 ipymesh
	ipcluster nbextension enable --user
	jupyter serverextension enable --py jupyterlab --sys-prefix
	jupyter nbextension enable --py --sys-prefix widgetsnbextension
	jupyter nbextension enable --py --sys-prefix bqplot
	jupyter nbextension enable --py --sys-prefix pythreejs
	jupyter nbextension enable --py --sys-prefix ipympl
	jupyter nbextension enable --py --sys-prefix ipymesh
	jupyter nbextension enable --py --sys-prefix ipyleaflet
	jupyter nbextension install --py --sys-prefix hide_code
	jupyter nbextension enable --py --sys-prefix hide_code
	jupyter nbextension install --py --sys-prefix rise
	jupyter nbextension enable --py --sys-prefix rise
	jupyter dashboards quick-setup --sys-prefix
	jupyter nbextension install --sys-prefix --py ipyparallel
	jupyter nbextension enable --sys-prefix --py ipyparallel
	jupyter serverextension enable --sys-prefix --py ipyparallel
	ipython profile create mpi --parallel
	printf "\nc.NotebookApp.server_extensions.append('ipyparallel.nbextension')\n" >> ${HOME}/.jupyter/jupyter_notebook_config.py
	printf "\nc.IPClusterEngines.engine_launcher_class = 'MPI'\n" >> ${HOME}/.ipython/profile_mpi/ipcluster_config.py
	printf "c.LocalControllerLauncher.controller_cmd = ['python2', '-m', 'ipyparallel.controller']\n" >> ${HOME}/.ipython/profile_mpi/ipcluster_config.py
	printf "c.LocalEngineSetLauncher.engine_cmd = ['python2', '-m', 'ipyparallel.engine']\n" >> ${HOME}/.ipython/profile_mpi/ipcluster_config.py
	printf "c.MPIEngineSetLauncher.engine_cmd = ['python2', '-m', 'ipyparallel.engine']\n" >> ${HOME}/.ipython/profile_mpi/ipcluster_config.py
	jupyter labextension install @jupyter-widgets/jupyterlab-manager@0.34.0

lfs:
	pip install pyliblzma
	wget https://github.com/git-lfs/git-lfs/releases/download/v1.5.5/git-lfs-linux-amd64-1.5.5.tar.gz
	tar xzvf git-lfs-linux-amd64-1.5.5.tar.gz
	cd git-lfs-1.5.5 && PREFIX=${HOME} ./install.sh
	export PATH=${HOME}/bin:${PATH}

proteus_pkg:
	cp stack/default.yaml stack/proteus_stack.yaml
	echo "  proteus:" >> stack/proteus_stack.yaml
	sed -i '/sources:/c\#sources:' stack/pkgs/proteus.yaml
	sed -i '/- key:/c\# -key:' stack/pkgs/proteus.yaml
	sed -i '/  url:/c\#  url:' stack/pkgs/proteus.yaml
	./stack/hit/bin/hit fetch https://github.com/erdc/proteus/archive/${PROTEUS_VERSION}.zip >> stack/pkgs/proteus.yaml
	cd stack && ${PROTEUS}/stack/hit/bin/hit build -j ${N} -v proteus_stack.yaml

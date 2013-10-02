all: install

PROTEUS ?= $(shell pwd)
PROTEUS_ARCH ?= cygwin
PROTEUS_PREFIX ?= ${PROTEUS}/${PROTEUS_ARCH}
PROTEUS_ENV ?= PATH="${PROTEUS_PREFIX}/bin:${PATH}" PYTHONPATH=${PROTEUS_PREFIX}/lib/python2.7/site-packages
PROTEUS_PYTHON ?= ${PROTEUS_ENV} /usr/bin/python

install:
	${PROTEUS_PYTHON} setuppyx.py install --prefix=${PROTEUS_PREFIX}
	${PROTEUS_PYTHON} setupf.py install --prefix=${PROTEUS_PREFIX}
	${PROTEUS_PYTHON} setuppetsc.py build --petsc-dir=${PROTEUS_PREFIX} --petsc-arch='' install --prefix=${PROTEUS_PREFIX}
	${PROTEUS_PYTHON} setup.py install --prefix=${PROTEUS_PREFIX}

build_ext:
	${PROTEUS_PYTHON} setup.py build_ext
	${PROTEUS_PYTHON} setuppyx.py build_ext 
	${PROTEUS_PYTHON} setupf.py build_ext 
	${PROTEUS_PYTHON} setuppetsc.py build_ext --petsc-dir=${PROTEUS_PREFIX} --petsc-arch=''

docs:
	doxygen proteus-doc.conf

clean:
	${PROTEUS_PYTHON} setuppyx.py clean
	${PROTEUS_PYTHON} setupf.py clean
	${PROTEUS_PYTHON} setuppetsc.py clean
	${PROTEUS_PYTHON} setup.py clean

cleaner: clean
	rm -rf ${PROTEUS_PREFIX}
	rm -rf build src/*.pyc src/*.so src/*.a

newConfig:
	cd proteusConfig && cp config.py.${PROTEUS_ARCH_OLD} config.py.${PROTEUS_ARCH}

stack: hashstack bootstrap.done
	cp hashstack/config.yml.${PROTEUS_ARCH} hashstack/config.yml 
	cp hashstack/packages.yml.${PROTEUS_ARCH} hashstack/packages.yml
	cd hashstack && ./update -v && ./update --copy ${PROTEUS_PREFIX}
	@echo "Stack complete, test with: make check"
	@echo "or: make parallel_check"

hashstack: 
	@echo "No hashstack found.  Cloning hashstack from GitHub"
	git clone -b proteus-cygwin-hacks-pre1.0 https://github.com/hashdist/hashstack.git

bootstrap: bootstrap.done

bootstrap.done: hashstack/setup_cygstack.py hashstack/cygstack.txt
	python hashstack/setup_cygstack.py hashstack/cygstack.txt
	touch bootstrap.done

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

parallel_check:
	@echo "************************"
	@echo "Parallel Proteus Partition Test"
	${PROTEUS_ENV} mpirun -np 4 python test/test_meshParitionFromTetgenFiles.py
	@echo "************************"

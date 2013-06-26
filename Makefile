all: install

install:
        ${PROTEUS_PYTHON} setuppyx.py install
	${PROTEUS_PYTHON} setupf.py install
        ${PROTEUS_PYTHON} setuppetsc.py build --petsc-dir=${PROTEUS_PREFIX} --petsc-arch='' install
        ${PROTEUS_PYTHON} setup.py install

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

stack:
	echo "You must run git submodule update --init to build hashstack"
	cp hashstack/config.yml.${PROTEUS_ARCH} hashstack/config.yml 
	cp hashstack/packages.yml.${PROTEUS_ARCH} hashstack/package.yml
	cd hashstack && ./update -v && ./update --copy ${PROTEUS_PREFIX}

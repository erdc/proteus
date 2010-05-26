all: install_externalPackages install_proteusModule

clean: clean_externalPackages clean_proteusModule

install_externalPackages:
	cd externalPackages && make all

install_proteusModule:
	cd proteusModule && cp proteusConfig/config.py.${PROTEUS_ARCH} config.py && python setup.py install

clean_externalPackages:
	cd externalPackages && make clean

clean_proteusModule:
	cd proteusModule && make cleaner


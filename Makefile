all: install_externalPackages install_pyadhModule

clean: clean_externalPackages clean_pyadhModule

install_externalPackages:
	cd externalPackages && make all

install_pyadhModule:
	cd pyadhModule && cp pyadhConfig/config.py.${PYADH_ARCH} config.py && python setup.py install

clean_externalPackages:
	cd externalPackages && make clean

clean_pyadhModule:
	cd pyadhModule && make cleaner


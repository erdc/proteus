all: install_externalPackages install_proteusModule

clean: clean_externalPackages clean_proteusModule

install_externalPackages:
	cd externalPackages && make all

install_proteusModule:
	cd proteusModule && cp proteusConfig/config.py.${PROTEUS_ARCH} config.py && make petsc-config install

clean_externalPackages:
	cd externalPackages && make -k distclean

clean_proteusModule:
	cd proteusModule && make cleaner

spkg:
	rm -rf proteus-0.9.0
	mkdir proteus-0.9.0
	cp spkg-install Makefile SPKG.txt proteus-0.9.0
	svn export proteusModule proteus-0.9.0/proteusModule
	svn export externalPackages proteus-0.9.0/externalPackages
	tar cjf proteus-0.9.0.spkg proteus-0.9.0
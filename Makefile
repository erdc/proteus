all: install_externalPackages install_proteusModule

clean: clean_externalPackages clean_proteusModule
    
cleaner:
	make clean
	rm -rf ${PROTEUS_PREFIX}    

install_externalPackages:
	cd externalPackages && make all

install_proteusModule:
	cd proteusModule && cp proteusConfig/config.py.${PROTEUS_ARCH} config.py && make petsc-config install

clean_externalPackages:
	cd externalPackages && make -k distclean

clean_proteusModule:
	cd proteusModule && make cleaner

newConfig:
	cd externalPackages && make newConfig
	cd proteusModule && make newConfig

spkg:
	rm -rf proteus
	mkdir proteus
	cp spkg-install Makefile SPKG.txt proteus
	svn export proteusModule proteus/proteusModule
	svn export externalPackages proteus/externalPackages
	tar cjf proteus.spkg proteus

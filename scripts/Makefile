default: $(PROTEUS_PREFIX)/bin/gmsh2tetgen

marinHexMesh: marinHexMesh.f 
	gfortran -ffree-form marinHexMesh.f -o marinHexMesh

gmsh2tetgen:gmsh2tetgen.f 
	gfortran -ffree-form gmsh2tetgen.f -o gmsh2tetgen

$(PROTEUS_PREFIX)/bin/gmsh2tetgen: gmsh2tetgen
	@echo ${PROTEUS_PREFIX}
	@echo "Installing gmsh2tegen in "
	@echo $(PROTEUS_PREFIX)/bin
	cp gmsh2tetgen $(PROTEUS_PREFIX)/bin

gmsh2tetgen_withtags:gmsh2tetgen_withtags.f 	
	gfortran -ffree-form gmsh2tetgen_withtags.f -o gmsh2tetgen_withtags

gmsh2tetgen_mprans:gmsh2tetgen_mprans.f 	
	gfortran -ffree-form gmsh2tetgen_mprans.f -o gmsh2tetgen_mprans

gmsh2tetgen_mprans_withtags:gmsh2tetgen_mprans_withtags.f 	
	gfortran -ffree-form gmsh2tetgen_mprans_withtags.f -o gmsh2tetgen_mprans_withtags

clean:
	rm -f marinHexMesh gmsh2tetgen

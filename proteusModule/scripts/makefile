include ${PROTEUS_PREFIX}/conf/variables


default: marinHexMesh gmsh2tetgen

marinHexMesh:
	$(FC) -ffree-form marinHexMesh.f -o marinHexMesh

gmsh2tetgen:	
	$(FC) -ffree-form gmsh2tetgen.f -o gmsh2tetgen
			
clean:
	rm -f marinHexMesh gmsh2tetgen

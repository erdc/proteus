include ${PROTEUS_PREFIX}/conf/variables


default: marinHexMesh gmsh2tetgen

marinHexMesh: marinHexMesh.f 
	$(FC) -ffree-form marinHexMesh.f -o marinHexMesh

gmsh2tetgen:gmsh2tetgen.f 	
	$(FC) -ffree-form gmsh2tetgen.f -o gmsh2tetgen
			
clean:
	rm -f marinHexMesh gmsh2tetgen

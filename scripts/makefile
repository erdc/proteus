default: marinHexMesh gmsh2tetgen

marinHexMesh: marinHexMesh.f 
	$(FC) -ffree-form marinHexMesh.f -o marinHexMesh

gmsh2tetgen:gmsh2tetgen.f 	
	$(FC) -ffree-form gmsh2tetgen.f -o gmsh2tetgen

gmsh2tetgen_withtags:gmsh2tetgen_withtags.f 	
	$(FC) -ffree-form gmsh2tetgen_withtags.f -o gmsh2tetgen_withtags

gmsh2tetgen_mprans:gmsh2tetgen_mprans.f 	
	$(FC) -ffree-form gmsh2tetgen_mprans.f -o gmsh2tetgen_mprans

gmsh2tetgen_mprans_withtags:gmsh2tetgen_mprans_withtags.f 	
	$(FC) -ffree-form gmsh2tetgen_mprans_withtags.f -o gmsh2tetgen_mprans_withtags

clean:
	rm -f marinHexMesh gmsh2tetgen

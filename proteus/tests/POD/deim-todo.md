#TODO list for adding DEIM

##DEIM interpolation test using fine grid evaluation

1. Split residual into time dependent and spatial term
   similar to current spatial and mass jacobians ..math::`F=F_s + F_t` 

	- test by just calling getResidual, getSpatialResidual, and getMassResidual separately and assering getResidual == getSpatialResidual + getMassResidual *DONE*
	
2. Add archiving of spatial residual, ..math::`F_s`, to Xdmf

	- add smoke test and check output in paraview? *DONE*
	
3. Add script for generating SVD for snapshots ..math::`\mathbf{F}_s = \mathbf{V}_f \mathbf{S}_f \mathbf{W}_f` and saves results to file 

	- just test that SVD algorithm gives correct factorization *DONE*
	
4. Implement DEIM algorithm,
	- just test that applying to full basis returns a full set of indices? *DONE*
	
5  Implement script that reads in snapshot info from file, calls DEIM, computes the projection matrix
	..math:: `\mathbf{P}_F=\mathbf{V}_f(\mathbf{P}^T\mathbf{V}_f)^{-1}`

	- make sure the error goes to essentially zero if use all of the basis vectors *DONE*

	- evaluate by visualizing full evaluation compared to DEIM and computing errr on 'fine' grid


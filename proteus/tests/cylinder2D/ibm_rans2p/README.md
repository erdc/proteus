# Flow around cylinder benchmark 2D-2

This test simulates the flow (density=1, viscosity=0.001) around the
cylinder (center=(0.2, 0.2), radius=0.05) in the rectangular domain
[0,2.2]x[0,0.41]. Please refer to [1,2] for mote details.


## Preprocessing

The file cylinder.py is called first, where `he`, depending on
`Refinement`, is used to control the mesh generation. The final time
of the problem is T=4. The boundary flags are ['left', 'right', 'top',
'bottom']. For boundary condition, the 'left' side is of inflow type,
the 'right' side is of outflow type, and other sides are enforced
no-slip boundary condition. The no-slip boundary condition over the
cylinder is implemented using the domain integral with the help of
delta function with support on the surface of the cylinder in the same
spirit as Nitsche method.

## Running

This problem can be run using the following command
```bash
    parun -v -l 5 cylinder_so.py -C "Refinement=5"  -D ibm_p1
```


## Postprocessing

Several physical parameters can be used to evaluate the computation
[1]. One is the drag/lift coefficient. It is computed in
*plot_lift_drag_force_from_particle.ipynb* by using the file
*particle_forceHistory.txt*.



## Reference 

1. Turek, Schaefer; Benchmark computations of laminar flow around
cylinder; in Flow Simulation with High-Performance Computers II, Notes
on Numerical Fluid Mechanics 52, 547-566, Vieweg 1996

2. John; Higher order Finite element methods and multigrid solvers in
a benchmark problem for the 3D Navier-Stokes equations;
Int. J. Numer. Meth. Fluids 2002; 40: 775-798 DOI:10.1002/d.377

# Flow around cylinder benchmark 2D-2

This test simulates the flow (density=1, viscoity=0.001) around the
cylinder (center=(0.2, 0.2), radius=0.05) in the rectangle domain
[0,2.2]x[0,0.41]. Please refer to [1,2] for mote details. It seems
none of finite element references cited in [1] can be found.

The whole algorithm is a projection scheme, where the momentum
equation is solved using RANS3P.py, the pressure increment is computed
in PreInc.py, and the pressure at the next step is updated in
Pressure.py.

## Preprocessing

The file cylinder.py is called first, where *he*, the maximum size of
edge, is passed into *symmetric2D* to generate the mesh. To increase
the accuracy of the cylinder, there are at least 40 points on the
boundary of the cylinder.  The final time of the problem is *T=8*. The
boundary flags are ['left', 'right', 'top', 'bottom', 'obstacle'].
For boundary condition, the 'left' side is of inflow type, the 'right'
side is of outflow type, and other sides are enforced non-slip
boundary condition.

## Running

This problem can be run using the following command
```bash
    parun -v -l 5 cylinder_so.py -C "he=0.01"
```


## Postprocessing

Several physics parameters can be used to evaluate the computation
[1]. One is the drag/lift coefficient. The drag/lift coefficient is
sensitive and gives us the direction to find the bug relating the
penalty coefficient in Nitsches' term. It is computed in
*plot-lift-drag-force-from-pv.ipynb* by using the file
*forceHistory_p.txt* and *forceHistory_p.txt*.

## Reference 

1. Turek, Schaefer; Benchmark computations of laminar flow around
cylinder; in Flow Simulation with High-Performance Computers II, Notes
on Numerical Fluid Mechanics 52, 547-566, Vieweg 1996

2. John; Higher order Finite element methods and multigrid solvers in
a benchmark problem for the 3D Navier-Stokes equations;
Int. J. Numer. Meth. Fluids 2002; 40: 775-798 DOI:10.1002/d.377
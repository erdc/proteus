# Problem
Use sbm to solve poisson equation $-\Delta u = 2$ on the disk domain with center at origin and radius $R=0.7$.

The true solution is $u=0.5(R^2-r^2)$.

# How to run it

Checkout the branch `yy/sbm_poisson` and run it with the command 

```
parun -v -l 5 poisson_so.py -b L2_batch.py
```

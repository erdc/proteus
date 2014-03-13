
Simple transient groundwater problem
====================================
    
    
    S_s\pd{h}{t} + \deld \vec q = s(x,t) \\

          \vec q =-\frac{K(x,t)}{\hat{\mu}(c(x,t))} \grad h


Here K is the Saturated Conductivity, K = k\varrho_0 |g|/\mu_0 and  S_s is the specific storage

Note that we assume S_h are independent of h and t and approximate
the accumulation term in the 'conservative form' 


    \pd{m}{t} = \pd{S_h h}{t}


since this is the form that proteus uses by default
 
To run with C_0 P^1 elements on triangles

::
parun sp_gw_p.py sp_gw_c0p1_n.py -v -l 2
::

To run with non-conforming P^1 on triangles
::
parun sp_gw_p.py sp_gw_ncp1_n.py -v -l 2
::

To run in parallel, set the parallel flag to True in sp_gw_p.py and type
::
mpirun -np 2 parun sp_gw_p.py sp_gw_c0p1_n.py -v -l 2 -O petsc.options
gatherArchives.py -s 2 -f TransientSinglePhaseFlow
::


Simple transient groundwater problem
====================================
    
We are trying to solve the problem .. math::    
    S_s\frac{\partial h}{\partial t} + \nabla \cdot \vec q = s(x,t) 

    \vec q =-\frac{K(x,t)}{\hat{\mu}(c(x,t))} \nabla h


Here K is the Saturated Conductivity, :math:`K = k\varrho_0 |g|/\mu_0` and  :math:`S_s` is the specific storage

Note that we assume S_h are independent of h and t and approximate
the accumulation term in the 'conservative form' math::


    \frac{\partial m}{\partial t} = \frac{\partial S_h h}{\partial t}


since this is the form that proteus uses by default
 
To run with :math:`C_0 P^1` elements on triangles::

   parun sp_gw_p.py sp_gw_c0p1_n.py -v -l 2

To run with non-conforming P^1 on triangles::

   parun sp_gw_p.py sp_gw_ncp1_n.py -v -l 2

To run in parallel, set the parallel flag to True in sp_gw_p.py and type::

   mpirun -np 2 parun sp_gw_p.py sp_gw_c0p1_n.py -v -l 2 -O petsc.options
   gatherArchives.py -s 2 -f TransientSinglePhaseFlow


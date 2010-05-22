from pyadh import *

import step2d

#problem definitions
nd = 2
#type of model (default k-e)
useSmagorinsky=False#True
#T = 9.999e0
T = 1.0e2#5.0e1#1.0e1

#shorter domain
upstream_height=0.5
downstream_height=1.
upstream_length  =1.0
downstream_length = 5.0

#for comparison with Griebel etal 
#upstream_height=0.75 #0.5
#downstream_height=1.5#1.
#upstream_length  =7.5#1.0
#downstream_length = 20.0#5.0
polyfile = "step2d"

#now try using inflow, Re = inflow*(downstream_height - upstream_height)/nu ?
ReInflow = 3025.0#20000.0#1000.0#250.#100.0
inflow   = 1.0
nu = inflow*(downstream_height - upstream_height)/ReInflow
rho= 998.2
g = [0.0,0.0]
#k-epsilon fudge factor
c_mu  = 0.09
#Smagorinsky fudge factor
smagorinskyConstant=0.1
#defaults otherwise for k-epsilon
print "trying ReInflow = %s inflow= %s nu=%s rho=%s " % (ReInflow,inflow,nu,rho)

boundaryTags = step2d.genPoly(fileprefix=polyfile,
                              upstream_height=upstream_height,
                              downstream_height=downstream_height,
                              upstream_length=upstream_length,
                              downstream_length=downstream_length,
                              step_fun=step2d.linear_profile,
                              n_points_on_step=2,
                              step_length=0.0)

#grid options
nn=3
triangleOptions = "Aq30Dena%f" % (0.15**2 / 6.0)
nLevels=1

space_quad_order = 3

#time integration
useBackwardEuler = False
timeOrder = 2

#solvers
atol_ns = 1.0e-4
atol_k  = 1.0e-4
atol_eps= 1.0e-4
linearSolverConvergenceTest = 'r-true'

#stabilization
lag_ns_subgridError = True
lag_k_subgridError  = True
lag_epsilon_subgridError = True

ns_shockCapturingFactor=0.99
k_shockCapturingFactor =0.99
epsilon_shockCapturingFactor=0.99


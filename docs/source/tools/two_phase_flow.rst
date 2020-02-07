.. _two_phase_flow:

TwoPhaseFlow
************


..
   Initial Conditions
   -----------------

   Boundary Conditions
   -------------------

   TwoPhaseFlow problem
   --------------------

   .. code-block:: python

       import proteus.TwoPhaseFlow as tpf

       # create problem instance
       myTpFlowProblem = TpFlow.TwoPhaseFlowProblem(
           ns_model=None,
           ls_model=None,
           nd=domain.nd,
           cfl=opts.cfl,
           outputStepping=outputStepping,
           structured=False,
           he=he,
           nnx=None,
           nny=None,
           nnz=None,
           domain=my_domain,
           initialConditions=initialConditions,
           boundaryConditions=None, # set with SpatialTools,
           useSuperlu=opts.useSuperlu,
       )

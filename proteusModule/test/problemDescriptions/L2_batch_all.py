simFlags['errorQuantities']=['u']
simFlags['errorTypes']= ['numericalSolution'] #compute error in soln and glob. mass bal
simFlags['errorNorms']= ['L2'] #compute L2 norm in space
simFlags['errorTimes']= ['Last']
simFlags['echo']=True
simFlags['dataFile']       = simFlags['simulationName']+'_results'
simFlags['dataDir']        = os.getcwd()+'/results'
simFlags['storeQuantities']= ['simulationData','errorData'] #include errorData for mass bal
simFlags['storeTimes']     = ['Last']
#simFlags['plotQuantities'].append('u_exact')
#simFlags['plotQuantities'].append('velocity')
nList[0].femSpaces = {0:FemTools.C0_AffineLinearOnSimplexWithNodalBasis,
                      1:FemTools.C0_AffineLinearOnSimplexWithNodalBasis,
                      2:FemTools.C0_AffineLinearOnSimplexWithNodalBasis}
nList[0].elementQuadrature = Quadrature.SimplexGaussQuadrature(pList[0].nd,2)
nList[0].elementBoundaryQuadrature = Quadrature.SimplexGaussQuadrature(pList[0].nd-1,2)
nList[0].nn=11
pList[0].name = "p1n11tx1"
nList[0].rtol_u[1] = 1.0e-1
nList[0].rtol_u[2] = 1.0e-1
nList[0].atol_u[1] = 1.0e-1
nList[0].atol_u[2] = 1.0e-1
nList[0].atol = 1.0e-6
nList[0].subgridError.nSteps=0
pList[0].coefficients = TransportCoefficients.TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                                      sigma=0.0,
                                                      rho_0=1.0,nu_0=1.0/pList[0].Re,
                                                      rho_1=1.0,nu_1=1.0/pList[0].Re,
                                                      g=[0.0,0.0],
                                                      nd=2,
                                                      LS_model=None,
                                                      KN_model=None,
                                                      epsFact_density=None,
                                                      stokes=False)
nList[0].subgridError = SubgridError.NavierStokesASGS_velocity_pressure(pList[0].coefficients,nList[0].nd,lag=True,delayLagSteps=5)
start
pList[0].name = "p1n11tx2"
nList[0].rtol_u[1] = 1.0e-2
nList[0].rtol_u[2] = 1.0e-2
nList[0].atol_u[1] = 1.0e-2
nList[0].atol_u[2] = 1.0e-2
nList[0].atol = 1.0e-6
nList[0].subgridError.nSteps=0
pList[0].coefficients = TransportCoefficients.TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                                      sigma=0.0,
                                                      rho_0=1.0,nu_0=1.0/pList[0].Re,
                                                      rho_1=1.0,nu_1=1.0/pList[0].Re,
                                                      g=[0.0,0.0],
                                                      nd=2,
                                                      LS_model=None,
                                                      KN_model=None,
                                                      epsFact_density=None,
                                                      stokes=False)
nList[0].subgridError = SubgridError.NavierStokesASGS_velocity_pressure(pList[0].coefficients,nList[0].nd,lag=True,delayLagSteps=5)
start
pList[0].name = "p1n11tx3"
nList[0].rtol_u[1] = 1.0e-3
nList[0].rtol_u[2] = 1.0e-3
nList[0].atol_u[1] = 1.0e-3
nList[0].atol_u[2] = 1.0e-3
nList[0].atol = 1.0e-6
nList[0].subgridError.nSteps=0
pList[0].coefficients = TransportCoefficients.TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                                      sigma=0.0,
                                                      rho_0=1.0,nu_0=1.0/pList[0].Re,
                                                      rho_1=1.0,nu_1=1.0/pList[0].Re,
                                                      g=[0.0,0.0],
                                                      nd=2,
                                                      LS_model=None,
                                                      KN_model=None,
                                                      epsFact_density=None,
                                                      stokes=False)
nList[0].subgridError = SubgridError.NavierStokesASGS_velocity_pressure(pList[0].coefficients,nList[0].nd,lag=True,delayLagSteps=5)
start
pList[0].name = "p1n11tx4"
nList[0].rtol_u[1] = 1.0e-4
nList[0].rtol_u[2] = 1.0e-4
nList[0].atol_u[1] = 1.0e-4
nList[0].atol_u[2] = 1.0e-4
nList[0].atol = 1.0e-6
nList[0].subgridError.nSteps=0
pList[0].coefficients = TransportCoefficients.TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                                      sigma=0.0,
                                                      rho_0=1.0,nu_0=1.0/pList[0].Re,
                                                      rho_1=1.0,nu_1=1.0/pList[0].Re,
                                                      g=[0.0,0.0],
                                                      nd=2,
                                                      LS_model=None,
                                                      KN_model=None,
                                                      epsFact_density=None,
                                                      stokes=False)
nList[0].subgridError = SubgridError.NavierStokesASGS_velocity_pressure(pList[0].coefficients,nList[0].nd,lag=True,delayLagSteps=5)
start
pList[0].name = "p1n21tx1"
nList[0].nn=21
nList[0].rtol_u[1] = 1.0e-1
nList[0].rtol_u[2] = 1.0e-1
nList[0].atol_u[1] = 1.0e-1
nList[0].atol_u[2] = 1.0e-1
nList[0].atol = 1.0e-6
nList[0].subgridError.nSteps=0
pList[0].coefficients = TransportCoefficients.TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                                      sigma=0.0,
                                                      rho_0=1.0,nu_0=1.0/pList[0].Re,
                                                      rho_1=1.0,nu_1=1.0/pList[0].Re,
                                                      g=[0.0,0.0],
                                                      nd=2,
                                                      LS_model=None,
                                                      KN_model=None,
                                                      epsFact_density=None,
                                                      stokes=False)
nList[0].subgridError = SubgridError.NavierStokesASGS_velocity_pressure(pList[0].coefficients,nList[0].nd,lag=True,delayLagSteps=5)
start
pList[0].name = "p1n21tx2"
nList[0].rtol_u[1] = 1.0e-2
nList[0].rtol_u[2] = 1.0e-2
nList[0].atol_u[1] = 1.0e-2
nList[0].atol_u[2] = 1.0e-2
nList[0].atol = 1.0e-6
nList[0].subgridError.nSteps=0
pList[0].coefficients = TransportCoefficients.TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                                      sigma=0.0,
                                                      rho_0=1.0,nu_0=1.0/pList[0].Re,
                                                      rho_1=1.0,nu_1=1.0/pList[0].Re,
                                                      g=[0.0,0.0],
                                                      nd=2,
                                                      LS_model=None,
                                                      KN_model=None,
                                                      epsFact_density=None,
                                                      stokes=False)
nList[0].subgridError = SubgridError.NavierStokesASGS_velocity_pressure(pList[0].coefficients,nList[0].nd,lag=True,delayLagSteps=5)
start
pList[0].name = "p1n21tx3"
nList[0].rtol_u[1] = 1.0e-3
nList[0].rtol_u[2] = 1.0e-3
nList[0].atol_u[1] = 1.0e-3
nList[0].atol_u[2] = 1.0e-3
nList[0].atol = 1.0e-6
nList[0].subgridError.nSteps=0
pList[0].coefficients = TransportCoefficients.TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                                      sigma=0.0,
                                                      rho_0=1.0,nu_0=1.0/pList[0].Re,
                                                      rho_1=1.0,nu_1=1.0/pList[0].Re,
                                                      g=[0.0,0.0],
                                                      nd=2,
                                                      LS_model=None,
                                                      KN_model=None,
                                                      epsFact_density=None,
                                                      stokes=False)
nList[0].subgridError = SubgridError.NavierStokesASGS_velocity_pressure(pList[0].coefficients,nList[0].nd,lag=True,delayLagSteps=5)
start
pList[0].name = "p1n21tx4"
nList[0].rtol_u[1] = 1.0e-4
nList[0].rtol_u[2] = 1.0e-4
nList[0].atol_u[1] = 1.0e-4
nList[0].atol_u[2] = 1.0e-4
nList[0].atol = 1.0e-6
nList[0].subgridError.nSteps=0
pList[0].coefficients = TransportCoefficients.TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                                      sigma=0.0,
                                                      rho_0=1.0,nu_0=1.0/pList[0].Re,
                                                      rho_1=1.0,nu_1=1.0/pList[0].Re,
                                                      g=[0.0,0.0],
                                                      nd=2,
                                                      LS_model=None,
                                                      KN_model=None,
                                                      epsFact_density=None,
                                                      stokes=False)
nList[0].subgridError = SubgridError.NavierStokesASGS_velocity_pressure(pList[0].coefficients,nList[0].nd,lag=True,delayLagSteps=5)
start
pList[0].name = "p1n41tx1"
nList[0].nn=41
nList[0].rtol_u[1] = 1.0e-1
nList[0].rtol_u[2] = 1.0e-1
nList[0].atol_u[1] = 1.0e-1
nList[0].atol_u[2] = 1.0e-1
nList[0].atol = 1.0e-6
nList[0].subgridError.nSteps=0
pList[0].coefficients = TransportCoefficients.TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                                      sigma=0.0,
                                                      rho_0=1.0,nu_0=1.0/pList[0].Re,
                                                      rho_1=1.0,nu_1=1.0/pList[0].Re,
                                                      g=[0.0,0.0],
                                                      nd=2,
                                                      LS_model=None,
                                                      KN_model=None,
                                                      epsFact_density=None,
                                                      stokes=False)
nList[0].subgridError = SubgridError.NavierStokesASGS_velocity_pressure(pList[0].coefficients,nList[0].nd,lag=True,delayLagSteps=5)
start
pList[0].name = "p1n41tx2"
nList[0].rtol_u[1] = 1.0e-2
nList[0].rtol_u[2] = 1.0e-2
nList[0].atol_u[1] = 1.0e-2
nList[0].atol_u[2] = 1.0e-2
nList[0].atol = 1.0e-6
nList[0].subgridError.nSteps=0
pList[0].coefficients = TransportCoefficients.TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                                      sigma=0.0,
                                                      rho_0=1.0,nu_0=1.0/pList[0].Re,
                                                      rho_1=1.0,nu_1=1.0/pList[0].Re,
                                                      g=[0.0,0.0],
                                                      nd=2,
                                                      LS_model=None,
                                                      KN_model=None,
                                                      epsFact_density=None,
                                                      stokes=False)
nList[0].subgridError = SubgridError.NavierStokesASGS_velocity_pressure(pList[0].coefficients,nList[0].nd,lag=True,delayLagSteps=5)
start
pList[0].name = "p1n41tx3"
nList[0].rtol_u[1] = 1.0e-3
nList[0].rtol_u[2] = 1.0e-3
nList[0].atol_u[1] = 1.0e-3
nList[0].atol_u[2] = 1.0e-3
nList[0].atol = 1.0e-6
nList[0].subgridError.nSteps=0
pList[0].coefficients = TransportCoefficients.TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                                      sigma=0.0,
                                                      rho_0=1.0,nu_0=1.0/pList[0].Re,
                                                      rho_1=1.0,nu_1=1.0/pList[0].Re,
                                                      g=[0.0,0.0],
                                                      nd=2,
                                                      LS_model=None,
                                                      KN_model=None,
                                                      epsFact_density=None,
                                                      stokes=False)
nList[0].subgridError = SubgridError.NavierStokesASGS_velocity_pressure(pList[0].coefficients,nList[0].nd,lag=True,delayLagSteps=5)
start
pList[0].name = "p1n41tx4"
nList[0].rtol_u[1] = 1.0e-4
nList[0].rtol_u[2] = 1.0e-4
nList[0].atol_u[1] = 1.0e-4
nList[0].atol_u[2] = 1.0e-4
nList[0].atol = 1.0e-6
nList[0].subgridError.nSteps=0
pList[0].coefficients = TransportCoefficients.TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                                      sigma=0.0,
                                                      rho_0=1.0,nu_0=1.0/pList[0].Re,
                                                      rho_1=1.0,nu_1=1.0/pList[0].Re,
                                                      g=[0.0,0.0],
                                                      nd=2,
                                                      LS_model=None,
                                                      KN_model=None,
                                                      epsFact_density=None,
                                                      stokes=False)
nList[0].subgridError = SubgridError.NavierStokesASGS_velocity_pressure(pList[0].coefficients,nList[0].nd,lag=True,delayLagSteps=5)
start
pList[0].name = "p1n81tx1"
nList[0].nn=81
nList[0].rtol_u[1] = 1.0e-1
nList[0].rtol_u[2] = 1.0e-1
nList[0].atol_u[1] = 1.0e-1
nList[0].atol_u[2] = 1.0e-1
nList[0].atol = 1.0e-6
nList[0].subgridError.nSteps=0
pList[0].coefficients = TransportCoefficients.TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                                      sigma=0.0,
                                                      rho_0=1.0,nu_0=1.0/pList[0].Re,
                                                      rho_1=1.0,nu_1=1.0/pList[0].Re,
                                                      g=[0.0,0.0],
                                                      nd=2,
                                                      LS_model=None,
                                                      KN_model=None,
                                                      epsFact_density=None,
                                                      stokes=False)
nList[0].subgridError = SubgridError.NavierStokesASGS_velocity_pressure(pList[0].coefficients,nList[0].nd,lag=True,delayLagSteps=5)
start
pList[0].name = "p1n81tx2"
nList[0].rtol_u[1] = 1.0e-2
nList[0].rtol_u[2] = 1.0e-2
nList[0].atol_u[1] = 1.0e-2
nList[0].atol_u[2] = 1.0e-2
nList[0].atol = 1.0e-6
nList[0].subgridError.nSteps=0
pList[0].coefficients = TransportCoefficients.TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                                      sigma=0.0,
                                                      rho_0=1.0,nu_0=1.0/pList[0].Re,
                                                      rho_1=1.0,nu_1=1.0/pList[0].Re,
                                                      g=[0.0,0.0],
                                                      nd=2,
                                                      LS_model=None,
                                                      KN_model=None,
                                                      epsFact_density=None,
                                                      stokes=False)
nList[0].subgridError = SubgridError.NavierStokesASGS_velocity_pressure(pList[0].coefficients,nList[0].nd,lag=True,delayLagSteps=5)
start
pList[0].name = "p1n81tx3"
nList[0].rtol_u[1] = 1.0e-3
nList[0].rtol_u[2] = 1.0e-3
nList[0].atol_u[1] = 1.0e-3
nList[0].atol_u[2] = 1.0e-3
nList[0].atol = 1.0e-6
nList[0].subgridError.nSteps=0
pList[0].coefficients = TransportCoefficients.TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                                      sigma=0.0,
                                                      rho_0=1.0,nu_0=1.0/pList[0].Re,
                                                      rho_1=1.0,nu_1=1.0/pList[0].Re,
                                                      g=[0.0,0.0],
                                                      nd=2,
                                                      LS_model=None,
                                                      KN_model=None,
                                                      epsFact_density=None,
                                                      stokes=False)
nList[0].subgridError = SubgridError.NavierStokesASGS_velocity_pressure(pList[0].coefficients,nList[0].nd,lag=True,delayLagSteps=5)
start
pList[0].name = "p1n81tx4"
nList[0].rtol_u[1] = 1.0e-4
nList[0].rtol_u[2] = 1.0e-4
nList[0].atol_u[1] = 1.0e-4
nList[0].atol_u[2] = 1.0e-4
nList[0].atol = 1.0e-6
nList[0].subgridError.nSteps=0
pList[0].coefficients = TransportCoefficients.TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                                      sigma=0.0,
                                                      rho_0=1.0,nu_0=1.0/pList[0].Re,
                                                      rho_1=1.0,nu_1=1.0/pList[0].Re,
                                                      g=[0.0,0.0],
                                                      nd=2,
                                                      LS_model=None,
                                                      KN_model=None,
                                                      epsFact_density=None,
                                                      stokes=False)
nList[0].subgridError = SubgridError.NavierStokesASGS_velocity_pressure(pList[0].coefficients,nList[0].nd,lag=True,delayLagSteps=5)
start
nList[0].femSpaces = {0:FemTools.C0_AffineQuadraticOnSimplexWithNodalBasis,
                      1:FemTools.C0_AffineQuadraticOnSimplexWithNodalBasis,
                      2:FemTools.C0_AffineQuadraticOnSimplexWithNodalBasis}
nList[0].elementQuadrature = Quadrature.SimplexGaussQuadrature(pList[0].nd,5)
nList[0].elementBoundaryQuadrature = Quadrature.SimplexGaussQuadrature(pList[0].nd-1,5)
nList[0].nn=11
pList[0].name = "p2n11tx1"
nList[0].rtol_u[1] = 1.0e-1
nList[0].rtol_u[2] = 1.0e-1
nList[0].atol_u[1] = 1.0e-1
nList[0].atol_u[2] = 1.0e-1
nList[0].atol = 1.0e-6
nList[0].subgridError.nSteps=0
pList[0].coefficients = TransportCoefficients.TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                                      sigma=0.0,
                                                      rho_0=1.0,nu_0=1.0/pList[0].Re,
                                                      rho_1=1.0,nu_1=1.0/pList[0].Re,
                                                      g=[0.0,0.0],
                                                      nd=2,
                                                      LS_model=None,
                                                      KN_model=None,
                                                      epsFact_density=None,
                                                      stokes=False)
nList[0].subgridError = SubgridError.NavierStokesASGS_velocity_pressure(pList[0].coefficients,nList[0].nd,lag=True,delayLagSteps=5,hFactor=0.5)
start
pList[0].name = "p2n11tx2"
nList[0].rtol_u[1] = 1.0e-2
nList[0].rtol_u[2] = 1.0e-2
nList[0].atol_u[1] = 1.0e-2
nList[0].atol_u[2] = 1.0e-2
nList[0].atol = 1.0e-6
nList[0].subgridError.nSteps=0
pList[0].coefficients = TransportCoefficients.TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                                      sigma=0.0,
                                                      rho_0=1.0,nu_0=1.0/pList[0].Re,
                                                      rho_1=1.0,nu_1=1.0/pList[0].Re,
                                                      g=[0.0,0.0],
                                                      nd=2,
                                                      LS_model=None,
                                                      KN_model=None,
                                                      epsFact_density=None,
                                                      stokes=False)
nList[0].subgridError = SubgridError.NavierStokesASGS_velocity_pressure(pList[0].coefficients,nList[0].nd,lag=True,delayLagSteps=5,hFactor=0.5)
start
pList[0].name = "p2n11tx3"
nList[0].rtol_u[1] = 1.0e-3
nList[0].rtol_u[2] = 1.0e-3
nList[0].atol_u[1] = 1.0e-3
nList[0].atol_u[2] = 1.0e-3
nList[0].atol = 1.0e-6
nList[0].subgridError.nSteps=0
pList[0].coefficients = TransportCoefficients.TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                                      sigma=0.0,
                                                      rho_0=1.0,nu_0=1.0/pList[0].Re,
                                                      rho_1=1.0,nu_1=1.0/pList[0].Re,
                                                      g=[0.0,0.0],
                                                      nd=2,
                                                      LS_model=None,
                                                      KN_model=None,
                                                      epsFact_density=None,
                                                      stokes=False)
nList[0].subgridError = SubgridError.NavierStokesASGS_velocity_pressure(pList[0].coefficients,nList[0].nd,lag=True,delayLagSteps=5,hFactor=0.5)
start
pList[0].name = "p2n11tx4"
nList[0].rtol_u[1] = 1.0e-4
nList[0].rtol_u[2] = 1.0e-4
nList[0].atol_u[1] = 1.0e-4
nList[0].atol_u[2] = 1.0e-4
nList[0].atol = 1.0e-6
nList[0].subgridError.nSteps=0
pList[0].coefficients = TransportCoefficients.TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                                      sigma=0.0,
                                                      rho_0=1.0,nu_0=1.0/pList[0].Re,
                                                      rho_1=1.0,nu_1=1.0/pList[0].Re,
                                                      g=[0.0,0.0],
                                                      nd=2,
                                                      LS_model=None,
                                                      KN_model=None,
                                                      epsFact_density=None,
                                                      stokes=False)
nList[0].subgridError = SubgridError.NavierStokesASGS_velocity_pressure(pList[0].coefficients,nList[0].nd,lag=True,delayLagSteps=5,hFactor=0.5)
start
pList[0].name = "p2n21tx1"
nList[0].nn=21
nList[0].rtol_u[1] = 1.0e-1
nList[0].rtol_u[2] = 1.0e-1
nList[0].atol_u[1] = 1.0e-1
nList[0].atol_u[2] = 1.0e-1
nList[0].atol = 1.0e-6
nList[0].subgridError.nSteps=0
pList[0].coefficients = TransportCoefficients.TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                                      sigma=0.0,
                                                      rho_0=1.0,nu_0=1.0/pList[0].Re,
                                                      rho_1=1.0,nu_1=1.0/pList[0].Re,
                                                      g=[0.0,0.0],
                                                      nd=2,
                                                      LS_model=None,
                                                      KN_model=None,
                                                      epsFact_density=None,
                                                      stokes=False)
nList[0].subgridError = SubgridError.NavierStokesASGS_velocity_pressure(pList[0].coefficients,nList[0].nd,lag=True,delayLagSteps=5,hFactor=0.5)
start
pList[0].name = "p2n21tx2"
nList[0].rtol_u[1] = 1.0e-2
nList[0].rtol_u[2] = 1.0e-2
nList[0].atol_u[1] = 1.0e-2
nList[0].atol_u[2] = 1.0e-2
nList[0].atol = 1.0e-6
nList[0].subgridError.nSteps=0
pList[0].coefficients = TransportCoefficients.TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                                      sigma=0.0,
                                                      rho_0=1.0,nu_0=1.0/pList[0].Re,
                                                      rho_1=1.0,nu_1=1.0/pList[0].Re,
                                                      g=[0.0,0.0],
                                                      nd=2,
                                                      LS_model=None,
                                                      KN_model=None,
                                                      epsFact_density=None,
                                                      stokes=False)
nList[0].subgridError = SubgridError.NavierStokesASGS_velocity_pressure(pList[0].coefficients,nList[0].nd,lag=True,delayLagSteps=5,hFactor=0.5)
start
pList[0].name = "p2n21tx3"
nList[0].rtol_u[1] = 1.0e-3
nList[0].rtol_u[2] = 1.0e-3
nList[0].atol_u[1] = 1.0e-3
nList[0].atol_u[2] = 1.0e-3
nList[0].atol = 1.0e-6
nList[0].subgridError.nSteps=0
pList[0].coefficients = TransportCoefficients.TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                                      sigma=0.0,
                                                      rho_0=1.0,nu_0=1.0/pList[0].Re,
                                                      rho_1=1.0,nu_1=1.0/pList[0].Re,
                                                      g=[0.0,0.0],
                                                      nd=2,
                                                      LS_model=None,
                                                      KN_model=None,
                                                      epsFact_density=None,
                                                      stokes=False)
nList[0].subgridError = SubgridError.NavierStokesASGS_velocity_pressure(pList[0].coefficients,nList[0].nd,lag=True,delayLagSteps=5,hFactor=0.5)
start
pList[0].name = "p2n21tx4"
nList[0].rtol_u[1] = 1.0e-4
nList[0].rtol_u[2] = 1.0e-4
nList[0].atol_u[1] = 1.0e-4
nList[0].atol_u[2] = 1.0e-4
nList[0].atol = 1.0e-6
nList[0].subgridError.nSteps=0
pList[0].coefficients = TransportCoefficients.TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                                      sigma=0.0,
                                                      rho_0=1.0,nu_0=1.0/pList[0].Re,
                                                      rho_1=1.0,nu_1=1.0/pList[0].Re,
                                                      g=[0.0,0.0],
                                                      nd=2,
                                                      LS_model=None,
                                                      KN_model=None,
                                                      epsFact_density=None,
                                                      stokes=False)
nList[0].subgridError = SubgridError.NavierStokesASGS_velocity_pressure(pList[0].coefficients,nList[0].nd,lag=True,delayLagSteps=5,hFactor=0.5)
start
pList[0].name = "p2n41tx1"
nList[0].nn=41
nList[0].rtol_u[1] = 1.0e-1
nList[0].rtol_u[2] = 1.0e-1
nList[0].atol_u[1] = 1.0e-1
nList[0].atol_u[2] = 1.0e-1
nList[0].atol = 1.0e-6
nList[0].subgridError.nSteps=0
pList[0].coefficients = TransportCoefficients.TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                                      sigma=0.0,
                                                      rho_0=1.0,nu_0=1.0/pList[0].Re,
                                                      rho_1=1.0,nu_1=1.0/pList[0].Re,
                                                      g=[0.0,0.0],
                                                      nd=2,
                                                      LS_model=None,
                                                      KN_model=None,
                                                      epsFact_density=None,
                                                      stokes=False)
nList[0].subgridError = SubgridError.NavierStokesASGS_velocity_pressure(pList[0].coefficients,nList[0].nd,lag=True,delayLagSteps=5,hFactor=0.5)
start
pList[0].name = "p2n41tx2"
nList[0].rtol_u[1] = 1.0e-2
nList[0].rtol_u[2] = 1.0e-2
nList[0].atol_u[1] = 1.0e-2
nList[0].atol_u[2] = 1.0e-2
nList[0].atol = 1.0e-6
nList[0].subgridError.nSteps=0
pList[0].coefficients = TransportCoefficients.TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                                      sigma=0.0,
                                                      rho_0=1.0,nu_0=1.0/pList[0].Re,
                                                      rho_1=1.0,nu_1=1.0/pList[0].Re,
                                                      g=[0.0,0.0],
                                                      nd=2,
                                                      LS_model=None,
                                                      KN_model=None,
                                                      epsFact_density=None,
                                                      stokes=False)
nList[0].subgridError = SubgridError.NavierStokesASGS_velocity_pressure(pList[0].coefficients,nList[0].nd,lag=True,delayLagSteps=5,hFactor=0.5)
start
pList[0].name = "p2n41tx3"
nList[0].rtol_u[1] = 1.0e-3
nList[0].rtol_u[2] = 1.0e-3
nList[0].atol_u[1] = 1.0e-3
nList[0].atol_u[2] = 1.0e-3
nList[0].atol = 1.0e-6
nList[0].subgridError.nSteps=0
pList[0].coefficients = TransportCoefficients.TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                                      sigma=0.0,
                                                      rho_0=1.0,nu_0=1.0/pList[0].Re,
                                                      rho_1=1.0,nu_1=1.0/pList[0].Re,
                                                      g=[0.0,0.0],
                                                      nd=2,
                                                      LS_model=None,
                                                      KN_model=None,
                                                      epsFact_density=None,
                                                      stokes=False)
nList[0].subgridError = SubgridError.NavierStokesASGS_velocity_pressure(pList[0].coefficients,nList[0].nd,lag=True,delayLagSteps=5,hFactor=0.5)
start
pList[0].name = "p2n41tx4"
nList[0].rtol_u[1] = 1.0e-4
nList[0].rtol_u[2] = 1.0e-4
nList[0].atol_u[1] = 1.0e-4
nList[0].atol_u[2] = 1.0e-4
nList[0].atol = 1.0e-6
nList[0].subgridError.nSteps=0
pList[0].coefficients = TransportCoefficients.TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                                      sigma=0.0,
                                                      rho_0=1.0,nu_0=1.0/pList[0].Re,
                                                      rho_1=1.0,nu_1=1.0/pList[0].Re,
                                                      g=[0.0,0.0],
                                                      nd=2,
                                                      LS_model=None,
                                                      KN_model=None,
                                                      epsFact_density=None,
                                                      stokes=False)
nList[0].subgridError = SubgridError.NavierStokesASGS_velocity_pressure(pList[0].coefficients,nList[0].nd,lag=True,delayLagSteps=5,hFactor=0.5)
start
pList[0].name = "p2n81tx1"
nList[0].nn=81
nList[0].rtol_u[1] = 1.0e-1
nList[0].rtol_u[2] = 1.0e-1
nList[0].atol_u[1] = 1.0e-1
nList[0].atol_u[2] = 1.0e-1
nList[0].atol = 1.0e-6
nList[0].subgridError.nSteps=0
pList[0].coefficients = TransportCoefficients.TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                                      sigma=0.0,
                                                      rho_0=1.0,nu_0=1.0/pList[0].Re,
                                                      rho_1=1.0,nu_1=1.0/pList[0].Re,
                                                      g=[0.0,0.0],
                                                      nd=2,
                                                      LS_model=None,
                                                      KN_model=None,
                                                      epsFact_density=None,
                                                      stokes=False)
nList[0].subgridError = SubgridError.NavierStokesASGS_velocity_pressure(pList[0].coefficients,nList[0].nd,lag=True,delayLagSteps=5,hFactor=0.5)
start
pList[0].name = "p2n81tx2"
nList[0].rtol_u[1] = 1.0e-2
nList[0].rtol_u[2] = 1.0e-2
nList[0].atol_u[1] = 1.0e-2
nList[0].atol_u[2] = 1.0e-2
nList[0].atol = 1.0e-6
nList[0].subgridError.nSteps=0
pList[0].coefficients = TransportCoefficients.TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                                      sigma=0.0,
                                                      rho_0=1.0,nu_0=1.0/pList[0].Re,
                                                      rho_1=1.0,nu_1=1.0/pList[0].Re,
                                                      g=[0.0,0.0],
                                                      nd=2,
                                                      LS_model=None,
                                                      KN_model=None,
                                                      epsFact_density=None,
                                                      stokes=False)
nList[0].subgridError = SubgridError.NavierStokesASGS_velocity_pressure(pList[0].coefficients,nList[0].nd,lag=True,delayLagSteps=5,hFactor=0.5)
start
pList[0].name = "p2n81tx3"
nList[0].rtol_u[1] = 1.0e-3
nList[0].rtol_u[2] = 1.0e-3
nList[0].atol_u[1] = 1.0e-3
nList[0].atol_u[2] = 1.0e-3
nList[0].atol = 1.0e-6
nList[0].subgridError.nSteps=0
pList[0].coefficients = TransportCoefficients.TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                                      sigma=0.0,
                                                      rho_0=1.0,nu_0=1.0/pList[0].Re,
                                                      rho_1=1.0,nu_1=1.0/pList[0].Re,
                                                      g=[0.0,0.0],
                                                      nd=2,
                                                      LS_model=None,
                                                      KN_model=None,
                                                      epsFact_density=None,
                                                      stokes=False)
nList[0].subgridError = SubgridError.NavierStokesASGS_velocity_pressure(pList[0].coefficients,nList[0].nd,lag=True,delayLagSteps=5,hFactor=0.5)
start
pList[0].name = "p2n81tx4"
nList[0].rtol_u[1] = 1.0e-4
nList[0].rtol_u[2] = 1.0e-4
nList[0].atol_u[1] = 1.0e-4
nList[0].atol_u[2] = 1.0e-4
nList[0].atol = 1.0e-6
nList[0].subgridError.nSteps=0
pList[0].coefficients = TransportCoefficients.TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                                      sigma=0.0,
                                                      rho_0=1.0,nu_0=1.0/pList[0].Re,
                                                      rho_1=1.0,nu_1=1.0/pList[0].Re,
                                                      g=[0.0,0.0],
                                                      nd=2,
                                                      LS_model=None,
                                                      KN_model=None,
                                                      epsFact_density=None,
                                                      stokes=False)
nList[0].subgridError = SubgridError.NavierStokesASGS_velocity_pressure(pList[0].coefficients,nList[0].nd,lag=True,delayLagSteps=5,hFactor=0.5)
start
quit

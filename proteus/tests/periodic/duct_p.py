from proteus import *
from proteus.default_p import *
from proteus.mprans import RANS2P
"""
Navier-Stokes flow in a periodic duct with square cross-section
"""
nd = 2

L=(4.0, 1.0, 1.0)
if nd == 2:
    L=(4.0, 1.0)
domain = Domain.RectangularDomain(L)
boundaryTags = domain.boundaryTags
unstructured = False
if unstructured:
    domain.writePoly("duct")
    if nd == 3:
        domain = Domain.PiecewiseLinearComplexDomain()
    elif nd == 2:
        domain = Domain.PlanarStraightLineGraphDomain()
    domain.readPoly("duct")

periodic = False
weakBC=True
coordBC=False
if periodic:
    gravity = [1.0e-1, 0., 0.]
else:
    gravity = [0., 0., 0.]

LevelModelType = RANS2P.LevelModel

coefficients = RANS2P.Coefficients(epsFact=0.0,
                                   sigma=0.0,
                                   rho_0=998.2,nu_0=1.004e-6,
                                   rho_1=998.2,nu_1=1.004e-6,
                                   g=gravity,
                                   nd=nd,
                                   ME_model=0,
                                   VF_model=None,
                                   LS_model=None,
                                   Closure_0_model=None,
                                   Closure_1_model=None,
                                   epsFact_density=0.0,
                                   stokes=False,
                                   useVF=0.0,
                                   useRBLES=0.0,
                                   useMetrics=1.0,
                                   eb_adjoint_sigma=1.0,
                                   eb_penalty_constant=100.0,
                                   forceStrongDirichlet=not weakBC,
                                   turbulenceClosureModel=0,
                                   NONCONSERVATIVE_FORM=1.0)
#boundaryCreatesNullSpace = True

T = 100.0
nsave=100
dt_init = 1.0e-3
DT = (T-dt_init)/float(nsave-1)
tnList = [0.0,dt_init]+[dt_init+i*DT for i in range(nsave)]

eps=1.0e-8
if periodic:
    def getPDBC(x,flag):
        if (x[0] < eps or x[0] > L[0] - eps) and (x[1] < eps or x[1] > L[1] - eps) and (x[2] < eps or x[2] > L[2] - eps):#x,y,z corner
            return numpy.array([0.0,0.0,0.0])
        elif (x[0] < eps or x[0] > L[0] - eps) and (x[2] < eps or x[2] > L[2] - eps):#x-z edge
            return numpy.array([0.0,round(x[1],5),0.0])
        elif (x[0] < eps or x[0] > L[0] - eps) and (x[1] < eps or x[1] > L[1] - eps):#x-y edge
            return numpy.array([0.0,0.0,round(x[2],5)])
        elif (x[1] < eps) or (x[1] > L[1]-eps):#on front or back
            return numpy.array([round(x[0],5),0.0,round(x[2],5)])
        elif (x[0] < eps) or (x[0] > L[0]-eps):#on inflow or outflow (left/right)
            return numpy.array([0.0,round(x[1],5),round(x[2],5)])

    periodicDirichletConditions = {0:getPDBC,
                                   1:getPDBC,
                                   2:getPDBC,
                                   3:getPDBC}
else:
    Re = 10000
    if nd == 3:
        inflow_v = Re*coefficients.nu/L[2]
    else:
        inflow_v = Re*coefficients.nu/L[1]

if periodic:
    def getDBC_pressure_duct(x,flag):
        if not periodic:
            if x[0] > L[0] - eps:
                return lambda x,t: 0.0

    def getDBC_u_duct(x,flag):
        if not periodic:
            if x[0] < eps:
                return lambda x,t: inflow_v
            elif x[0] > L[0] - eps:
                return lambda x,t: 0.0
        if x[2] <  eps or x[2] > L[2] - eps:#top and bottom: no slip
            return lambda x,t: 0.0

    def getDBC_v_duct(x,flag):
        if not periodic:
            if x[0] > L[0] - eps:
                return lambda x,t: 0.0
            elif x[0] < eps:
                return lambda x,t: 0.0
        if x[2] <  eps or x[2] > L[2] - eps:#top and bottom: no slip
            return lambda x,t: 0.0

    def getDBC_w_duct(x,flag):
        if not periodic:
            if x[0] > L[0] - eps:
                return lambda x,t: 0.0
            elif x[0] < eps:
                return lambda x,t: 0.0
        if x[2] <  eps or x[2] > L[2] - eps:#top and bottom: no slip
            return lambda x,t: 0.0

    dirichletConditions = {0:getDBC_pressure_duct,
                           1:getDBC_u_duct,
                           2:getDBC_v_duct,
                           3:getDBC_w_duct}


    def getAFBC_p_duct(x,flag):
        if not periodic:
            if x[0] < eps:
                return lambda x,t: -inflow_v
            elif x[0] > L[0] - eps:
                return None
        if (x[2] < eps) or (x[2] > L[2] - eps):#top and bottom: no flow
            return lambda x,t: 0.0
        else:
            return lambda x,t: 0.0

    def getAFBC_u_duct(x,flag):
        if not periodic:
            if x[0] < eps:
                return None
            elif x[0] > L[0] - eps:
                return None
        if (x[2] < eps) or (x[2] > L[2] - eps):#top and bottom: no flow
            return lambda x,t: 0.0
        else:
            return lambda x,t: 0.0

    def getAFBC_v_duct(x,flag):
        if not periodic:
            if x[0] < eps:
                return None
            elif x[0] > L[0] - eps:
                return None
        if (x[2] < eps) or (x[2] > L[2] - eps):#top and bottom: no flow
            return lambda x,t: 0.0
        else:
            return lambda x,t: 0.0

    def getAFBC_w_duct(x,flag):
        if not periodic:
            if x[0] < eps:
                return None
            elif x[0] > L[0] - eps:
                return None
        if (x[2] < eps) or (x[2] > L[2] - eps):#top and bottom: no flow
            return lambda x,t: 0.0
        else:
            return lambda x,t: 0.0

    advectiveFluxBoundaryConditions =  {0:getAFBC_p_duct,
                                        1:getAFBC_u_duct,
                                        2:getAFBC_v_duct,
                                        3:getAFBC_w_duct}

    def getDFBC_duct(x,flag):
        if not periodic:
            if x[0] < eps:
                return None
            elif x[0] > L[0] - eps:
                return lambda x,t: 0.0
        if (x[2] < eps) or (x[2] > L[2] - eps):#top and bottom: no slip
            return None
        else:
            return lambda x,t: 0.0

    diffusiveFluxBoundaryConditions = {0:{},
                                       1:{1:getDFBC_duct},
                                       2:{2:getDFBC_duct},
                                       3:{3:getDFBC_duct}}

else:
    if  coordBC:
        if nd == 3:
            def onLeft(x):
                return x[0] < eps and x[2] > eps and x[2] < L[2] - eps
            def onRight(x):
                return x[0] > L[0] - eps and x[2] > eps and x[2] < L[2] - eps
            def onFront(x):
                return x[1] < eps and x[2] > eps and x[2] < L[2] - eps and x[0] > eps and x[0] < L[0] - eps
            def onBack(x):
                return x[1] > L[1] - eps and x[2] > eps and x[2] < L[2] - eps and x[0] > eps and x[0] < L[0] - eps
            def onBottom(x):
                return x[2] < eps
            def onTop(x):
                return x[2] > L[2] - eps
        elif nd == 2:
            def onLeft(x):
                return x[0] < eps and x[1] > eps and x[1] < L[1] - eps
            def onRight(x):
                return x[0] > L[0] - eps and x[1] > eps and x[1] < L[1] - eps
            def onBottom(x):
                return x[1] < eps
            def onTop(x):
                return x[1] > L[1] - eps

        def getDBC_pressure_duct(x,flag):
            if onRight(x):
                return lambda x,t: 0.0
            
        def getDBC_u_duct(x,flag):
            if onLeft(x):
                return lambda x,t: inflow_v
            if weakBC and onRight(x):
                return lambda x,t: 0.0
            if onTop(x) or onBottom(x):
                return lambda x,t: 0.0
            
        def getDBC_v_duct(x,flag):
            if onLeft(x) or onRight(x) or onTop(x) or onBottom(x):
                return lambda x,t: 0.0

        def getDBC_w_duct(x,flag):
            if onLeft(x) or onRight(x) or onTop(x) or onBottom(x):
                return lambda x,t: 0.0

        dirichletConditions = {0:getDBC_pressure_duct,
                               1:getDBC_u_duct,
                               2:getDBC_v_duct}
        if nd == 3:
            dirichletConditions[3] =  getDBC_w_duct


        def getAFBC_p_duct(x,flag):
            if onLeft(x):
                return lambda x,t: -inflow_v
            if onTop(x) or onBottom(x):
                return lambda x,t: 0.0
            if nd == 3:
                if onFront(x) or onBack(x):
                    return lambda x,t: 0.0

        def getAFBC_u_duct(x,flag):
            if nd == 3:
                if onFront(x) or onBack(x):
                    return lambda x,t: 0.0

        def getAFBC_v_duct(x,flag):
            if nd == 3:
                if onFront(x) or onBack(x):
                    return lambda x,t: 0.0

        advectiveFluxBoundaryConditions =  {0:getAFBC_p_duct,
                                            1:getAFBC_u_duct,
                                            2:getAFBC_v_duct}
        if nd == 3:
            def getAFBC_w_duct(x,flag):
                if onFront(x) or onBack(x):
                    return lambda x,t: 0.0
            advectiveFluxBoundaryConditions[3] = getAFBC_w_duct

        def getDFBC_duct(x,flag):
            if onRight(x):
                return lambda x,t: 0.0
            if nd == 3:
                if onFront(x) or onBack(x):
                    return lambda x,t: 0.0

        diffusiveFluxBoundaryConditions = {0:{},
                                           1:{1:getDFBC_duct},
                                           2:{2:getDFBC_duct}}
        if nd == 3:
            diffusiveFluxBoundaryConditions[3] = {3:getDFBC_duct}
    else:
        def getDBC_pressure_duct(x,flag):
            if flag == boundaryTags['right']:
                return lambda x,t: 0.0

        def getDBC_u_duct(x,flag):
            if flag == boundaryTags['left']:
                return lambda x,t: inflow_v
            if weakBC and flag == boundaryTags['right']:
                return lambda x,t: 0.0
            if flag in [boundaryTags['top'], boundaryTags['bottom']]:
                return lambda x,t: 0.0
        def getDBC_v_duct(x,flag):
            if flag in [boundaryTags['left'],
                        boundaryTags['right'],
                        boundaryTags['top'],
                        boundaryTags['bottom']]:
                return lambda x,t: 0.0

        def getDBC_w_duct(x,flag):
            if flag in [boundaryTags['left'],
                        boundaryTags['right'],
                        boundaryTags['top'],
                        boundaryTags['bottom']]:
                return lambda x,t: 0.0

        dirichletConditions = {0:getDBC_pressure_duct,
                               1:getDBC_u_duct,
                               2:getDBC_v_duct}
        if nd == 3:
            dirichletConditions[3] =  getDBC_w_duct


        def getAFBC_p_duct(x,flag):
            if flag == boundaryTags['left']:
                return lambda x,t: -inflow_v
            if flag in [boundaryTags['top'],
                        boundaryTags['bottom']]:
                return lambda x,t: 0.0
            if nd == 3 and flag in [boundaryTags['front'],
                                    boundaryTags['back'],
                                    0]:
                return lambda x,t: 0.0
            elif nd == 2 and flag == 0:
                return lambda x,t: 0.0

        def getAFBC_u_duct(x,flag):
            if nd == 3 and flag in [boundaryTags['front'],
                                    boundaryTags['back'],
                                    0]:
                return lambda x,t: 0.0
            elif nd == 2 and flag == 0:
                return lambda x,t: 0.0

        def getAFBC_v_duct(x,flag):
            if nd == 3 and flag in [boundaryTags['front'],
                                    boundaryTags['back'],
                                    0]:
                return lambda x,t: 0.0
            elif nd == 2 and flag == 0:
                return lambda x,t: 0.0
        def getAFBC_w_duct(x,flag):
            if flag in [boundaryTags['front'],
                        boundaryTags['back'],
                        0]:
                return lambda x,t: 0.0

        advectiveFluxBoundaryConditions =  {0:getAFBC_p_duct,
                                            1:getAFBC_u_duct,
                                            2:getAFBC_v_duct}
        if nd == 3:
            advectiveFluxBoundaryConditions[3] = getAFBC_w_duct

        def getDFBC_duct(x,flag):
            if flag == boundaryTags['right']:#outflow
                return lambda x,t: 0.0
            if nd == 3 and flag in [boundaryTags['front'],
                                    boundaryTags['back'],
                                    0]:
                return lambda x,t: 0.0
            elif nd == 2 and flag == 0:
                return lambda x,t: 0.0
        diffusiveFluxBoundaryConditions = {0:{},
                                           1:{1:getDFBC_duct},
                                           2:{2:getDFBC_duct}}
        if nd == 3:
            diffusiveFluxBoundaryConditions[3] = {3:getDFBC_duct}

import math
import numpy
from matplotlib import pyplot as plt
from proteus import (MeshTools,
                     AnalyticalSolutions,
                     canalyticalSolutions)
from proteus.Profiling import logEvent
nPoints_x = 34
nPoints_y = 34
nPoints_z = 34
dx = 1.0/(nPoints_x-1.0)
dy = 1.0/(nPoints_y-1.0)
dz = 1.0/(nPoints_z-1.0)
x = numpy.zeros((nPoints_x,nPoints_y,nPoints_z,3),'d')
for i in range(nPoints_x):
    for j in range(nPoints_y):
        for k in range(nPoints_z):
            x[i,j,k,0] = i*dx
            x[i,j,k,1] = j*dy
            x[i,j,k,2] = k*dz

def test_PlaneCouetteFlow():
    iwork = numpy.zeros((1,),'i')
    rwork = numpy.zeros((5,),'d')
    t=0.0
    rwork[0]=0.1
    rwork[1]=dy*(nPoints_y-1.0)
    rwork[2]=0.0
    rwork[3]=0.0
    rwork[4]=0.0
    ux = numpy.zeros(x.shape[:-1],'d')
    uy = numpy.zeros(x.shape[:-1],'d')
    canalyticalSolutions.PlaneCouetteFlow_u(iwork,rwork,t,x,ux)
    slice=nPoints_z/2
    fig = plt.figure()
    plt.contourf(x[:,:,slice,0],
                 x[:,:,slice,1],
                 ux[:,:,slice])
    plt.title("PlaneCouette")
    plt.savefig("PlaneCouette.png")
    fig = plt.figure()
    plt.quiver(x[:,:,slice,0],
               x[:,:,slice,1],
               ux[:,:,slice],
               uy[:,:,slice])
    plt.title("PlaneCouette Velocity")
    plt.savefig("PlaneCouetteQuiver.png")

def test_PlaneBase():
    from proteus.AnalyticalSolutions import PlaneBase
    for i in range(5):
        for j in range(5):
            for k in range(5):
                for l in range(5):
                    p = PlaneBase(plane_theta=i*math.pi/4.0,
                                  plane_phi=j*math.pi/4.0-math.pi/2.0,
                                  v_theta=k*math.pi/4.0,
                                  v_phi=l*math.pi/4.0-math.pi/2.0,
                                  v_norm=1.0,
                                  mu=1.0,
                                  grad_p=1.0)

def test_NonlinearDAE
    grid = MeshTools.RectangularGrid(500,1,1,Lx=1.0)
    X=numpy.array([n.p for n in grid.nodeList])
    Xx=numpy.array([n.p[0] for n in grid.nodeList])
    logEvent("Testing Solutions in 1D")
    logEvent("NonlinearDAE")
    from proteus.AnalyticalSolutions import NonlinearDAE
    sol = NonlinearDAE(5.0,1.0)
    y1=[sol.uOfT(x) for x in Xx]
    sol = NonlinearDAE(5.0,2.0)
    y2=[sol.uOfT(x) for x in Xx]
    sol = NonlinearDAE(5.0,3.0)
    y3=[sol.uOfT(x) for x in Xx]
    sol = NonlinearDAE(5.0,0.75)
    y075=[sol.uOfT(x) for x in Xx]
    sol = NonlinearDAE(5.0,0.25)
    y025=[sol.uOfT(x) for x in Xx]
    sol = NonlinearDAE(5.0,-0.25)
    yn025=[sol.uOfT(x) for x in Xx]
    sol = NonlinearDAE(5.0,-1.25)
    yn125=[sol.uOfT(x) for x in Xx]
    fig = plt.figure()
    plt.plot(Xx,y1,label='Solution,q=1')
    plt.plot(Xx,y2,label='Solution,q=2')
    plt.plot(Xx,y3,label='Solution,q=3')
    plt.plot(Xx,y075,label='Solution,q=0.75')
    plt.plot(Xx,y025,label='Solution,q=0.25')
    plt.plot(Xx,yn025,label='Solution,q=-0.25')
    plt.plot(Xx,yn125,label='Solution,q=-1.25')
    plt.legend()
    plt.savefig("NonlinearDAE.png")

def test_LinearAD_SteadyState():
    logEvent("AD_SteadyState")
    sol=AnalyticalSolutions.LinearAD_SteadyState()
    fig = plt.figure()
    plt.plot(Xx,[sol.uOfX(x) for x in X],label='Solution,q=1,r=1')
    plt.legend()
    plt.savefig("LinearAD_SteadyState.png")

def test_NonlinearAD_SteadyState():
    logEvent("NonlinearAD_SteadyState")
    sol=AnalyticalSolutions.NonlinearAD_SteadyState(q=2,r=1)
    fig = plt.figure()
    plt.plot(Xx,
             [sol.uOfX(x) for x in X],
             label='Solution,q=2,r=1')
    sol=AnalyticalSolutions.NonlinearAD_SteadyState(q=1,r=2)
    plt.plot(Xx,
             [sol.uOfX(x) for x in X],
             label='Solution,q=1,r=2')
    plt.legend()
    plt.savefig("NonlinearAD_SteadyStage.png")

def test_LinearADR_Sine():
    logEvent("LinearADR_Sine")
    sol=AnalyticalSolutions.LinearADR_Sine()
    fig = plt.figure()
    plt.title("LinearADR_Sine")
    plt.plot(Xx,
             [sol.uOfX(x) for x in X],
             label='Solution')
    plt.plot(Xx,
             [sol.duOfX(x)[0] for x in X],
             label='Gradient')
    plt.plot(Xx,
             [sol.advectiveFluxOfX(x)[0] for x in X],
             label='Advective Flux')
    plt.plot(Xx,
             [sol.diffusiveFluxOfX(x)[0] for x in X],
             label='Diffusive Flux')
    plt.plot(Xx,
             [sol.totalFluxOfX(x)[0] for x in X],
             label='Total Flux')
    plt.plot(Xx,
             [sol.rOfUX(sol.uOfX(x),x) for x in X],
             label='reaction')
    plt.legend()
    plt.savefig("LinearADR_Sine.png")

def test_LinearAD_DiracIC();
    logEvent("LinearAD_DiracIC")
    sol=AnalyticalSolutions.LinearAD_DiracIC()
    plt.plot(Xx,
             [sol.uOfXT(x,T=0.25) for x in X],
             label='Solution,t=0.25')
    plt.plot(Xx,
             [sol.uOfXT(x,T=0.75) for x in X],
             label='Solution,t=0.75')
    plt.plot(Xx,
             [sol.duOfXT(x,T=0.25)[0] for x in X],
             label='Gradient,t=0.25')
    plt.plot(Xx,
             [sol.duOfXT(x,T=0.75)[0] for x in X],
             label='Gradient,t=0.75')
    plt.plot(Xx,
             [sol.advectiveFluxOfXT(x,T=0.25)[0] for x in X],
             label='Advective Flux,t=0.25')
    plt.plot(Xx,
             [sol.advectiveFluxOfXT(x,T=0.75)[0] for x in X],
             label='Advective Flux,t=0.75')
    plt.plot(Xx,
             [sol.diffusiveFluxOfXT(x,T=0.25)[0] for x in X],
             label='Diffusive Flux,t=0.25')
    plt.plot(Xx,
             [sol.diffusiveFluxOfXT(x,T=0.75)[0] for x in X],
             label='Diffusive Flux,t=0.75')
    plt.plot(Xx,
             [sol.totalFluxOfXT(x,T=0.25)[0] for x in X],
             label='Total Flux,t=0.25')
    plt.plot(Xx,
             [sol.totalFluxOfXT(x,T=0.75)[0] for x in X],
             label='Total Flux,t=0.75')
    plt.legend()
    plt.savefig("LinearAD_DiracIC.png")

def test_LinearADR_Decay_DiracIC():
    logEvent("LinearADR_Decay_DiracIC")
    sol=AnalyticalSolutions.LinearADR_Decay_DiracIC()
    plt.plot(Xx,[sol.uOfXT(x,T=0.25) for x in X],
             label='Solution,t=0.25')
    plt.plot(Xx,[sol.uOfXT(x,T=0.75) for x in X],
             label='Solution,t=0.75')
    plt.plot(Xx,
             [sol.duOfXT(x,T=0.25)[0] for x in X],
             label='Gradient,t=0.25'),
    plt.plot(Xx,
             [sol.duOfXT(x,T=0.75)[0] for x in X],
             label='Gradient,t=0.75')
    plt.plot(Xx,
             [sol.advectiveFluxOfXT(x,T=0.25)[0] for x in X],
             label='Advective Flux,t=0.25')
    plt.plot(Xx,
             [sol.advectiveFluxOfXT(x,T=0.75)[0] for x in X],
             label='Advective Flux,t=0.75')
    plt.plot(Xx,
             [sol.diffusiveFluxOfXT(x,T=0.25)[0] for x in X],
             label='Diffusive Flux,t=0.25')
    plt.plot(Xx,
             [sol.diffusiveFluxOfXT(x,T=0.75)[0] for x in X],
             label='Diffusive Flux,t=0.75')
    plt.plot(Xx,
             [sol.totalFluxOfXT(x,T=0.25)[0] for x in X],
             label='Total Flux,t=0.25')
    plt.plot(Xx,
             [sol.totalFluxOfXT(x,T=0.75)[0] for x in X],
             label='Total Flux,t=0.75')
    plt.plot(Xx,
             [sol.rOfUXT(sol.uOfXT(x,T=0.25),x,T=0.25) for x in X],
             label='Reaction,T=0.25')
    plt.plot(Xx,
             [sol.rOfUXT(sol.uOfXT(x,T=0.75),x,T=0.75) for x in X],
             label='Reaction,t=0.75')
    plt.legend()
    plt.savefix("LinearADR_Decay_DiracIC.png")

def test_NonlinearADR_Decay_DiracIC():
    logEvent("NonlinearADR_Decay_DiracIC")
    sol=AnalyticalSolutions.NonlinearADR_Decay_DiracIC()
    plt.plot(Xx,
             [sol.uOfXT(x,T=0.25) for x in X],
             label='Solution,t=0.25'),
    plt.plot(Xx,
             [sol.uOfXT(x,T=0.75) for x in X],
             label='Solution,t=0.75'))
    plt.plot(Xx,
             [sol.duOfXT(x,T=0.25)[0] for x in X],
             label='Gradient,t=0.25'),
    plt.plot(Xx,
             [sol.duOfXT(x,T=0.75)[0] for x in X],
             label='Gradient,t=0.75'))
    plt.plot(Xx,
             [sol.advectiveFluxOfXT(x,T=0.25)[0] for x in X],
             label='Advective Flux,t=0.25')
    plt.plot(Xx,
             [sol.advectiveFluxOfXT(x,T=0.75)[0] for x in X],
             label='Advective Flux,t=0.75')
    plt.plot(Xx,
             [sol.diffusiveFluxOfXT(x,T=0.25)[0] for x in X],
             label='Diffusive Flux,t=0.25'),
    plt.plot(Xx,
             [sol.diffusiveFluxOfXT(x,T=0.75)[0] for x in X],
             label='Diffusive Flux,t=0.75')
    plt.plot(Xx,
             [sol.totalFluxOfXT(x,T=0.25)[0] for x in X],
             label='Total Flux,t=0.25')
    plt.plot(Xx,
             [sol.totalFluxOfXT(x,T=0.75)[0] for x in X],
             label='Total Flux,t=0.75')
    plt.plot(Xx,
             [sol.rOfUXT(sol.uOfXT(x,T=0.25),x,T=0.25) for x in X],
             label='Reaction,t=0.25')
    plt.plot(Xx,
             [sol.rOfUXT(sol.uOfXT(x,T=0.75),x,T=0.75) for x in X],
             label='Reaction,t=0.75')
    plt.legend("NonlinearADR_Decay_DiracIC.png")
    plt.savefig()

def test_NonlinearDAE_f():
    iwork = numpy.zeros((1,),'i')
    rwork = numpy.zeros((2,),'d')
    rwork[0]=5.0
    rwork[1]=2.0
    t=1.0
    u = numpy.zeros((nPoints_x,nPoints_y,nPoints_z,3),'d')
    # x is the initial u.
    u = x.copy()
    canalyticalSolutions.NonlinearDAE_f(iwork,rwork,t,x,u)
    slice = nPoints_z/2
    plt.contourf(x[:,:,0],
                 x[:,:,1],
                 u[:,:,slice,0])
    plt.savefig("NonlinearDAE_f.png")

def test_poissonsEquationExp1D():
    iwork = numpy.zeros((1,),'i')
    rwork = numpy.zeros((1,),'d')
    rwork[0]=5
    nPoints = 101
    t=0.0
    x = numpy.zeros((nPoints,3),'d')
    for i in range(nPoints):
        x[i,0] = i*(1.0/(nPoints-1.0))
    u = numpy.zeros(x.shape[0],'d')
    canalyticalSolutions.poissonsEquationExp1D(iwork,rwork,t,x,u)
    plt.plot(x[:,0],
             u,
             label='poissonsEquationExp1D')
    plt.savefig("poissonsEquationsExp1D.png")

def test_poissonsEquationsExp2D():
    iwork = numpy.zeros((1,),'i')
    rwork = numpy.zeros((1,),'d')
    rwork[0]=5
    nPoints_x = 51
    nPoints_y = 51
    t=0.0
    x = numpy.zeros((nPoints_x,nPoints_y,3),'d')
    for i in range(nPoints_x):
        for j in range(nPoints_y):
            x[i,j,0] = i*(1.0/(nPoints_x-1.0))
            x[i,j,1] = j*(1.0/(nPoints_y-1.0))
    u = numpy.zeros(x.shape[:-1],'d')
    canalyticalSolutions.poissonsEquationExp2D(iwork,rwork,t,x,u)
    plt.contourf(x[:,:,0],
                 x[:,:,1],
                 u,
                 label='poissonsEquationExp2D')
    plt.savefig("poissonsEquationExp2D.png")

def test_poissonsEquationExp3D():
    iwork = numpy.zeros((1,),'i')
    rwork = numpy.zeros((1,),'d')
    rwork[0]=5
    t=0.0
    u = numpy.zeros(x.shape[:-1],'d')
    canalyticalSolutions.poissonsEquationExp3D(iwork,rwork,t,x,u)
    slice=nPoints_x/2
    plt.contourf(x[slice,:,:,0],
                 x[slice,:,:,1],
                 u[slice,:,:])
    plt.title('poissonsEquationExp3D-x/2')
    plt.savefig('poissonsEquationExp3D_x/2.png')
    slice=nPoints_y/2
    plt.contourf(x[:,slice,:,0],
                 x[:,slice,:,1],
                 u[:,slice,:])
    plt.title('poissonsEquationExp3D-y/2')
    plt.savefig('poissonsEquationExp3D-y/2.png')
    slice=nPoints_z/2
    plt.contourf(x[:,:,slice,0],
                 x[:,:,slice,1],
                 u[:,:,slice])
    plt.title('poissonsEquationExp3D-z/2')
    plt.savefig('poissonsEquationExp3D-z/2.png')

def test_diffusionSin1D():
    iwork = numpy.zeros((1,),'i')
    iwork[0]=1
    rwork = numpy.zeros((1,),'d')
    nPoints = 101
    t=0.0
    x = numpy.zeros((nPoints,3),'d')
    for i in range(nPoints):
        x[i,0] = i*(1.0/(nPoints-1.0))
    u = numpy.zeros(x.shape[0],'d')
    canalyticalSolutions.diffusionSin1D(iwork,rwork,t,x,u)
    plt.plot(x[:,0],
             u)
    plt.title('diffusionSin1D')
    plt.savefig('diffusionSin1D')

def test_diffusionSin2D():
    iwork = numpy.zeros((2,),'i')
    iwork[0]=1
    iwork[1]=1
    rwork = numpy.zeros((1,),'d')
    nPoints_x = 51
    nPoints_y = 51
    t=0.0
    x = numpy.zeros((nPoints_x,nPoints_y,3),'d')
    for i in range(nPoints_x):
        for j in range(nPoints_y):
            x[i,j,0] = i*(1.0/(nPoints_x-1.0))
            x[i,j,1] = j*(1.0/(nPoints_y-1.0))
    u = numpy.zeros(x.shape[:-1],'d')
    canalyticalSolutions.diffusionSin2D(iwork,rwork,t,x,u)
    plt.contourf(x[:,:,0],
                 x[:,:,1],
                 u)
    plt.title('diffusionSin2D')
    plt.savefig('diffusionSin2D.png')

def test_diffusionSin3D():
    iwork = numpy.zeros((3,),'i')
    iwork[0]=1
    iwork[1]=1
    iwork[2]=1
    rwork = numpy.zeros((1,),'d')
    t=0.0
    u = numpy.zeros(x.shape[:-1],'d')
    canalyticalSolutions.diffusionSin3D(iwork,rwork,t,x,u)

def test_STflowSphere_P():
    iwork = numpy.zeros((1,),'i')
    rwork = numpy.zeros((8,),'d')
    rwork[0]=0.1
    rwork[1]=0.0
    rwork[2]=0.0
    nPoints_x = 31
    nPoints_y = 31
    nPoints_z = 31
    t=0.0
    dx = 1.0/(nPoints_x-1.0)
    dy = 1.0/(nPoints_y-1.0)
    dz = 1.0/(nPoints_z-1.0)
    rwork[3]=0.1
    rwork[4]=0.5
    rwork[5]=0.5
    rwork[6]=0.5
    rwork[7]=1.003e-3
    x = numpy.zeros((nPoints_x,nPoints_y,nPoints_z,3),'d')
    for i in range(nPoints_x):
        for j in range(nPoints_y):
            for k in range(nPoints_z):
                x[i,j,k,0] = i*dx
                x[i,j,k,1] = j*dy
                x[i,j,k,2] = k*dz
    p = numpy.zeros(x.shape[:-1],'d')
    canalyticalSolutions.STflowSphere_P(iwork,rwork,t,x,p)
    slice=nPoints_x/2
    plt.contourf(x[slice,:,:,1],
                 x[slice,:,:,2],
                 p[slice,:,:])
    plt.title('STflowSphere_P: yz plane'))
    slice=nPoints_y/2
    plt.contourf(x[:,slice,:,0],
                 x[:,slice,:,2],
                 p[:,slice,:])
    plt.title('STflowSphere_P: xz plane')
    slice=nPoints_z/2
    plt.contourf(x[:,:,slice,0],
                 x[:,:,slice,1],
                 p[:,:,slice])
     plt.title('STflowSphere_P: xy plane')
     plt.savefig("STflowSphere_P.png")
     ux = numpy.zeros(x.shape[:-1],'d')
     canalyticalSolutions.STflowSphere_Vx(iwork,rwork,t,x,ux)
     uy = numpy.zeros(x.shape[:-1],'d')
     canalyticalSolutions.STflowSphere_Vy(iwork,rwork,t,x,uy)
     uz = numpy.zeros(x.shape[:-1],'d')
     canalyticalSolutions.STflowSphere_Vz(iwork,rwork,t,x,uz)
     slice=nPoints_z/2
     plt.quiver(x[:,:,slice,0],
                x[:,:,slice,1],
                ux[:,:,slice],
                uy[:,:,slice])
     plot.title('STflowSphere : xy plane')
     plot.savefig('STflowSphere_xy.png')
     slice=nPoints_y/2
     plt.quiver(x[:,slice,:,0],
                x[:,slice,:,2],
                ux[:,slice,:],
                uz[:,slice,:])
     plt.title('STflowSphere : xz plane')
     plt.savefig('STflowSphere_xz.png')
     slice=nPoints_x/2
     plt.quiver(x[slice,:,:,1],
                x[slice,:,:,2],
                uy[slice,:,:],
                uz[slice,:,:])
     plt.title('STflowSphere : yz plane')
     plt.savefig('STflowSphere_yz.png')

def test_PoiseuilleFlow():
     iwork = numpy.zeros((1,),'i')
     rwork = numpy.zeros((7,),'d')
     t=0.0
     rwork[0]=dy*(nPoints_y-1.0)
     rwork[1]=0.001
     rwork[2]=-1.0
     # rwork[3]=83.3
     rwork[4]=0.0
     rwork[5]=0.0
     rwork[6]=0.0
     ux = numpy.zeros(x.shape[:-1],'d')
     uy = numpy.zeros(x.shape[:-1],'d')
     canalyticalSolutions.PoiseuilleFlow(iwork,rwork,t,x,ux)
     slice=nPoints_z/2
     plt.contourf(x[:,:,slice,0],
                  x[:,:,slice,1],
                  u[:,:,slice])
     plt.title("Poiseuille Flow")
     plt.savefig("PoiseuilleFlow.png")
     plt.quiver(x[:,:,slice,0],
                x[:,:,slice,1],
                ux[:,:,slice],
                uy[:,:,slice])
     plt.title("Poiseuille Flow Velocity")
     plt.savefig("PoiseuilleFlowQuiver.png")

def test_PoiseuillePipeFlow():
    iwork = numpy.zeros((1,),'i')
    rwork = numpy.zeros((7,),'d')
    t=0.0
    dx = 1.0/(nPoints_x-1.0)
    dy = 1.0/(nPoints_y-1.0)
    dz = 1.0/(nPoints_z-1.0)
    rwork[0]=dy*(nPoints_y-1.0)/2.0
    rwork[1]=0.001
    rwork[2]=-1.0
    # rwork[3]=24.5
    rwork[4]=0.5
    rwork[5]=0.5
    rwork[6]=0.5
    ux = numpy.zeros(x.shape[:-1],'d')
    uy = numpy.zeros(x.shape[:-1],'d')
    canalyticalSolutions.PoiseuillePipeFlow(iwork,rwork,t,x,ux)
    slice=nPoints_z/2
    plt.contourf(x[:,:,slice,0],
                 x[:,:,slice,1],
                 ux[:,:,slice])
    plt.title('Poiseuille Pipe Flow')
    plt.savefig("PoiseuillePipeFlow.png")
    fig = plt.figure()
    plt.quiver(x[:,:,slice,0],
               x[:,:,slice,1],
               ux[:,:,slice],
               uy[:,:,slice])
    plt.savefig("PoiseuillePipeFlowQuiver.png")

def test_PoiseuillePipeFlow_P():
     iwork = numpy.zeros((1,),'i')
     rwork = numpy.zeros((8,),'d')
     t=0.0
     dx = 1.0/(nPoints_x-1.0)
     dy = 1.0/(nPoints_y-1.0)
     dz = 1.0/(nPoints_z-1.0)
     rwork[0]=dy*(nPoints_y-1.0)/2.0
     rwork[1]=0.001
     rwork[2]=-1.0
     # rwork[3]=24.5
     rwork[4]=1.0
     rwork[5]=0.5
     rwork[6]=0.5
     rwork[7]=0.5
     ux = numpy.zeros(x.shape[:-1],'d')
     uy = numpy.zeros(x.shape[:-1],'d')
     canalyticalSolutions.PoiseuillePipeFlow_P(iwork,rwork,t,x,ux)
     slice=nPoints_z/2
     plt.contourf(x[:,:,slice,0],
                  x[:,:,slice,1],
                  ux[:,:,slice])
     plt.title('Poiseuille Pipe Flow: Pressure'))
     plt.savefig("PoiseuillePipeFlow_P.png")

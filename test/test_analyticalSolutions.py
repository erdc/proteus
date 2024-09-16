import math
import numpy as np
from proteus import (MeshTools,
                     AnalyticalSolutions,
                     canalyticalSolutions)
from proteus.Profiling import logEvent

save_plots=False
try:
    import matplotlib
    matplotlib.use('AGG')
    from matplotlib import pyplot as plt
    save_plots = True
except:
    save_plots = False
nPoints_x = 35
nPoints_y = 35
nPoints_z = 35
midpoint_x_index = (nPoints_x + 1) // 2
midpoint_y_index = (nPoints_y + 1) // 2
midpoint_z_index = (nPoints_z + 1) // 2
assert (nPoints_x + 1) % 2 == 0
assert (nPoints_y + 1) % 2 == 0
assert (nPoints_z + 1) % 2 == 0
dx = 1.0 / (nPoints_x-1.0)
dy = 1.0 / (nPoints_y-1.0)
dz = 1.0 / (nPoints_z-1.0)
x = np.zeros((nPoints_x, nPoints_y, nPoints_z, 3), 'd')
for i in range(nPoints_x):
    for j in range(nPoints_y):
        for k in range(nPoints_z):
            x[i, j, k, 0] = i*dx
            x[i, j, k, 1] = j*dy
            x[i, j, k, 2] = k*dz

grid = MeshTools.RectangularGrid(500, 1, 1, Lx=1.0)
X = np.array([n.p for n in grid.nodeList])
Xx = np.array([n.p[0] for n in grid.nodeList])


def test_PlaneCouetteFlow():
    iwork = np.zeros((1,), 'i')
    rwork = np.zeros((5,), 'd')
    t = 0.0
    rwork[0] = 0.1
    rwork[1] = dy*(nPoints_y-1.0)
    rwork[2] = 0.0
    rwork[3] = 0.0
    rwork[4] = 0.0
    ux = np.zeros(x.shape[:-1], 'd')
    uy = np.zeros(x.shape[:-1], 'd')
    canalyticalSolutions.PlaneCouetteFlow_u(iwork, rwork, t, x, ux)
    assert ux[0, 0, 0] == 0.0
    assert ux[midpoint_x_index, midpoint_y_index, midpoint_z_index] == 0.05294117647058824
    assert ux[-1, -1, -1] == 0.1
    slice = nPoints_z//2
    if save_plots:
        fig = plt.figure()
        plt.contourf(x[:, :, slice, 0],
                     x[:, :, slice, 1],
                     ux[:, :, slice])
        plt.title("PlaneCouette")
        plt.savefig("PlaneCouette.png")
        fig = plt.figure()
        plt.quiver(x[:, :, slice, 0],
                   x[:, :, slice, 1],
                   ux[:, :, slice],
                   uy[:, :, slice])
        plt.title("PlaneCouette Velocity")
        plt.savefig("PlaneCouetteQuiver.png")


def test_NonlinearDAE():
    logEvent("Testing Solutions in 1D")
    logEvent("NonlinearDAE")
    from proteus.AnalyticalSolutions import NonlinearDAE
    sol = NonlinearDAE(5.0, 1.0)
    y1 = [sol.uOfT(x) for x in Xx]
    assert y1[midpoint_x_index] == 0.8349689658715417
    sol = NonlinearDAE(5.0, 2.0)
    y2 = [sol.uOfT(x) for x in Xx]
    assert y2[midpoint_x_index] == 0.8471986417657046
    sol = NonlinearDAE(5.0, 3.0)
    y3 = [sol.uOfT(x) for x in Xx]
    assert y3[midpoint_x_index] == 0.8572655778618113
    sol = NonlinearDAE(5.0, 0.75)
    y075 = [sol.uOfT(x) for x in Xx]
    assert y075[midpoint_x_index] == 0.8314754625643243
    sol = NonlinearDAE(5.0, 0.25)
    y025 = [sol.uOfT(x) for x in Xx]
    assert y025[midpoint_x_index] == 0.8238351907493889
    sol = NonlinearDAE(5.0, -0.25)
    yn025 = [sol.uOfT(x) for x in Xx]
    assert yn025[midpoint_x_index] == 0.8151530579734992
    sol = NonlinearDAE(5.0, -1.25)
    yn125 = [sol.uOfT(x) for x in Xx]
    assert yn125[midpoint_x_index] == 0.7934541671666694
    if save_plots:
        fig = plt.figure()
        plt.plot(Xx, y1, label='Solution,q=1')
        plt.plot(Xx, y2, label='Solution,q=2')
        plt.plot(Xx, y3, label='Solution,q=3')
        plt.plot(Xx, y075, label='Solution,q=0.75')
        plt.plot(Xx, y025, label='Solution,q=0.25')
        plt.plot(Xx, yn025, label='Solution,q=-0.25')
        plt.plot(Xx, yn125, label='Solution,q=-1.25')
        plt.legend()
        plt.title("NonlinearDAE")
        plt.savefig("NonlinearDAE.png")


def test_LinearAD_SteadyState():
    logEvent("AD_SteadyState")
    sol = AnalyticalSolutions.LinearAD_SteadyState()
    y = [sol.uOfX(x) for x in X]
    assert y[midpoint_x_index] == 0.9882908500750013
    if save_plots:
        fig = plt.figure()
        plt.plot(Xx, y, label='Solution,q=1,r=1')
        plt.legend()
        plt.title("LinearAD_SteadyStage")
        plt.savefig("LinearAD_SteadyState.png")


def test_NonlinearAD_SteadyState():
    logEvent("NonlinearAD_SteadyState")
    sol = AnalyticalSolutions.NonlinearAD_SteadyState(q=2, r=1)
    y = [sol.uOfX(x) for x in X]
    assert y[midpoint_x_index] == 0.9948469998751492
    if save_plots:
        fig = plt.figure()
        plt.plot(Xx,
                 y,
                 label='Solution,q=2,r=1')
    sol = AnalyticalSolutions.NonlinearAD_SteadyState(q=1, r=2)
    y = [sol.uOfX(x) for x in X]
    assert y[midpoint_x_index] == 0.990588835252627
    if save_plots:
        fig = plt.figure()
        plt.plot(Xx,
                 [sol.uOfX(x) for x in X],
                 label='Solution,q=1,r=2')
        plt.legend()
        plt.title("NonlinearAD_SteadyState")
        plt.savefig("NonlinearAD_SteadyState.png")


def test_LinearADR_Sine():
    logEvent("LinearADR_Sine")
    sol = AnalyticalSolutions.LinearADR_Sine()
    u = [sol.uOfX(x) for x in X]
    assert u[midpoint_x_index] == 0.22471248696198207
    Du = [sol.duOfX(x)[0] for x in X]
    assert Du[midpoint_x_index] == 6.122493544431436
    aflux = [sol.advectiveFluxOfX(x)[0] for x in X]
    assert aflux[midpoint_x_index] == 0.22471248696198207
    dflux = [sol.diffusiveFluxOfX(x)[0] for x in X]
    assert dflux[midpoint_x_index] == -0.061224935444314364
    tflux = [sol.totalFluxOfX(x)[0] for x in X]
    assert tflux[midpoint_x_index] == 0.1634875515176677
    r = [sol.rOfUX(sol.uOfX(x), x) for x in X]
    assert r[midpoint_x_index] == -6.211206478443425
    if save_plots:
        fig = plt.figure()
        plt.title("LinearADR_Sine")
        plt.plot(Xx,
                 u,
                 label='Solution')
        plt.plot(Xx,
                 Du,
                 label='Gradient')
        plt.plot(Xx,
                 aflux,
                 label='Advective Flux')
        plt.plot(Xx,
                 dflux,
                 label='Diffusive Flux')
        plt.plot(Xx,
                 tflux,
                 label='Total Flux')
        plt.plot(Xx,
                 r,
                 label='reaction')
        plt.legend()
        plt.savefig("LinearADR_Sine.png")


def test_LinearAD_DiracIC():
    logEvent("LinearAD_DiracIC")
    sol = AnalyticalSolutions.LinearAD_DiracIC()
    up25 = [sol.uOfXT(x, T=0.25) for x in X]
    assert up25[midpoint_x_index] == 0.005805917122996999
    up75 = [sol.uOfXT(x, T=0.75) for x in X]
    assert up75[midpoint_x_index] == 1.362394143014481e-08
    Dup25 = [sol.duOfXT(x, T=0.25)[0] for x in X]
    assert Dup25[midpoint_x_index] == -0.24840948011219627
    Dup75 = [sol.duOfXT(x, T=0.75)[0] for x in X]
    assert Dup75[midpoint_x_index] == -6.484340861040867e-07
    afluxp25 = [sol.advectiveFluxOfXT(x, T=0.25)[0] for x in X]
    assert afluxp25[midpoint_x_index] == 0.005805917122996999
    afluxp75 = [sol.advectiveFluxOfXT(x, T=0.75)[0] for x in X]
    assert afluxp75[midpoint_x_index] == 1.362394143014481e-08
    dfluxp25 = [sol.diffusiveFluxOfXT(x, T=0.25)[0] for x in X]
    assert dfluxp25[midpoint_x_index] == 0.0024840948011219627
    dfluxp75 = [sol.diffusiveFluxOfXT(x, T=0.75)[0] for x in X]
    assert dfluxp75[midpoint_x_index] == 6.484340861040867e-09
    tfluxp25 = [sol.totalFluxOfXT(x, T=0.25)[0] for x in X]
    assert tfluxp25[midpoint_x_index] == 0.008290011924118962
    tfluxp75 = [sol.totalFluxOfXT(x, T=0.75)[0] for x in X]
    assert tfluxp75[midpoint_x_index] == 2.0108282291185675e-08
    if save_plots:
        fig = plt.figure()
        plt.plot(Xx,
                 up25,
                 label='Solution,t=0.25')
        plt.plot(Xx,
                 up75,
                 label='Solution,t=0.75')
        plt.plot(Xx,
                 Dup25,
                 label='Gradient,t=0.25')
        plt.plot(Xx,
                 Dup75,
                 label='Gradient,t=0.75')
        plt.plot(Xx,
                 afluxp25,
                 label='Advective Flux,t=0.25')
        plt.plot(Xx,
                 afluxp75,
                 label='Advective Flux,t=0.75')
        plt.plot(Xx,
                 dfluxp25,
                 label='Diffusive Flux,t=0.25')
        plt.plot(Xx,
                 dfluxp75,
                 label='Diffusive Flux,t=0.75')
        plt.plot(Xx,
                 tfluxp25,
                 label='Total Flux,t=0.25')
        plt.plot(Xx,
                 tfluxp75,
                 label='Total Flux,t=0.75')
        plt.legend()
        plt.title("LinearAD_DiracIC.png")
        plt.savefig("LinearAD_DiracIC.png")


def test_LinearADR_Decay_DiracIC():
    logEvent("LinearADR_Decay_DiracIC")
    sol = AnalyticalSolutions.LinearADR_Decay_DiracIC()
    up25 = [sol.uOfXT(x, T=0.25) for x in X]
    assert up25[midpoint_x_index] == 0.0045216528018377404
    up75 = [sol.uOfXT(x, T=0.75) for x in X]
    assert up75[midpoint_x_index] == 6.435494248102993e-09
    Dup25 = [sol.duOfXT(x, T=0.25)[0] for x in X]
    assert Dup25[midpoint_x_index] == -0.19346149763373902
    Dup75 = [sol.duOfXT(x, T=0.75)[0] for x in X]
    assert Dup75[midpoint_x_index] == -3.062985739327577e-07
    afluxp25 = [sol.advectiveFluxOfXT(x, T=0.25)[0] for x in X]
    assert afluxp25[midpoint_x_index] == 0.0045216528018377404
    afluxp75 = [sol.advectiveFluxOfXT(x, T=0.75)[0] for x in X]
    assert afluxp75[midpoint_x_index] == 6.435494248102993e-09
    dfluxp25 = [sol.diffusiveFluxOfXT(x, T=0.25)[0] for x in X]
    assert dfluxp25[midpoint_x_index] == 0.0019346149763373901
    dfluxp75 = [sol.diffusiveFluxOfXT(x, T=0.75)[0] for x in X]
    assert dfluxp75[midpoint_x_index] == 3.062985739327577e-09
    tfluxp25 = [sol.totalFluxOfXT(x, T=0.25)[0] for x in X]
    assert tfluxp25[midpoint_x_index] == 0.00645626777817513
    tfluxp75 = [sol.totalFluxOfXT(x, T=0.75)[0] for x in X]
    assert tfluxp75[midpoint_x_index] == 9.498479987430571e-09
    if save_plots:
        fig = plt.figure()
        plt.plot(Xx,
                 up25,
                 label='Solution,t=0.25')
        plt.plot(Xx,
                 up75,
                 label='Solution,t=0.75')
        plt.plot(Xx,
                 Dup25,
                 label='Gradient,t=0.25')
        plt.plot(Xx,
                 Dup75,
                 label='Gradient,t=0.75')
        plt.plot(Xx,
                 afluxp25,
                 label='Advective Flux,t=0.25')
        plt.plot(Xx,
                 afluxp75,
                 label='Advective Flux,t=0.75')
        plt.plot(Xx,
                 dfluxp25,
                 label='Diffusive Flux,t=0.25')
        plt.plot(Xx,
                 dfluxp75,
                 label='Diffusive Flux,t=0.75')
        plt.plot(Xx,
                 tfluxp25,
                 label='Total Flux,t=0.25')
        plt.plot(Xx,
                 tfluxp75,
                 label='Total Flux,t=0.75')
        plt.legend()
        plt.title("LinearADR_Decay_DiracIC.png")
        plt.savefig("LinearADR_Decay_DiracIC.png")


def test_NonlinearADR_Decay_DiracIC():
    logEvent("NonlinearADR_Decay_DiracIC")
    sol = AnalyticalSolutions.NonlinearADR_Decay_DiracIC()
    up25 = [sol.uOfXT(x, T=0.25) for x in X]
    assert up25[midpoint_x_index] == 0.0058003017280384575
    up75 = [sol.uOfXT(x, T=0.75) for x in X]
    assert up75[midpoint_x_index] == 1.3623941337338921e-08
    Dup25 = [sol.duOfXT(x, T=0.25)[0] for x in X]
    assert Dup25[midpoint_x_index] == -0.248169222231705570
    Dup75 = [sol.duOfXT(x, T=0.75)[0] for x in X]
    assert Dup75[midpoint_x_index] == -6.484340816869727e-07
    afluxp25 = [sol.advectiveFluxOfXT(x, T=0.25)[0] for x in X]
    assert afluxp25[midpoint_x_index] == 0.0058003017280384575
    afluxp75 = [sol.advectiveFluxOfXT(x, T=0.75)[0] for x in X]
    assert afluxp75[midpoint_x_index] == 1.3623941337338921e-08
    dfluxp25 = [sol.diffusiveFluxOfXT(x, T=0.25)[0] for x in X]
    assert dfluxp25[midpoint_x_index] == 0.0024816922223170556
    dfluxp75 = [sol.diffusiveFluxOfXT(x, T=0.75)[0] for x in X]
    assert dfluxp75[midpoint_x_index] == 6.484340816869727e-09
    tfluxp25 = [sol.totalFluxOfXT(x, T=0.25)[0] for x in X]
    assert tfluxp25[midpoint_x_index] == 0.008281993950355513
    tfluxp75 = [sol.totalFluxOfXT(x, T=0.75)[0] for x in X]
    assert tfluxp75[midpoint_x_index] == 2.010828215420865e-08
    if save_plots:
        fig = plt.figure()
        plt.plot(Xx,
                 up25,
                 label='Solution,t=0.25')
        plt.plot(Xx,
                 up75,
                 label='Solution,t=0.75')
        plt.plot(Xx,
                 Dup25,
                 label='Gradient,t=0.25')
        plt.plot(Xx,
                 Dup75,
                 label='Gradient,t=0.75')
        plt.plot(Xx,
                 afluxp25,
                 label='Advective Flux,t=0.25')
        plt.plot(Xx,
                 afluxp75,
                 label='Advective Flux,t=0.75')
        plt.plot(Xx,
                 dfluxp25,
                 label='Diffusive Flux,t=0.25')
        plt.plot(Xx,
                 dfluxp75,
                 label='Diffusive Flux,t=0.75')
        plt.plot(Xx,
                 tfluxp25,
                 label='Total Flux,t=0.25')
        plt.plot(Xx,
                 tfluxp75,
                 label='Total Flux,t=0.75')
        plt.legend()
        plt.title("NonlinearADR_Decay_DiracIC.png")
        plt.savefig("NonlinearADR_Decay_DiracIC.png")


def test_NonlinearDAE_f():
    iwork = np.zeros((1,), 'i')
    rwork = np.zeros((2,), 'd')
    rwork[0] = 5.0
    rwork[1] = 2.0
    t = 1.0
    f = x.copy()
    canalyticalSolutions.NonlinearDAE_f(iwork, rwork, t, x, f)
    assert f[midpoint_x_index, midpoint_y_index, midpoint_z_index][0] == -1.4013840830449829
    if save_plots:
        fig = plt.figure()
        slice = nPoints_z//2
        plt.contourf(x[:, :, slice, 0],
                     x[:, :, slice, 1],
                     f[:, :, slice, 0])
        plt.title("NonlinearDAE_f.png")
        plt.savefig("NonlinearDAE_f.png")


def test_poissonsEquationExp1D():
    iwork = np.zeros((1,), 'i')
    rwork = np.zeros((1,), 'd')
    rwork[0] = 5
    nPoints = 101
    t = 0.0
    x = np.zeros((nPoints, 3), 'd')
    for i in range(nPoints):
        x[i, 0] = i*(1.0/(nPoints-1.0))
    u = np.zeros(x.shape[0], 'd')
    canalyticalSolutions.poissonsEquationExp1D(iwork, rwork, t, x, u)
    assert u[(nPoints_x+1)//2] == 0.7623027790507058
    if save_plots:
        fig = plt.figure()
        plt.plot(x[:, 0],
                 u,
                 label='poissonsEquationExp1D')
        plt.title("poissonsEquationsExp1D.png")
        plt.savefig("poissonsEquationsExp1D.png")


def test_poissonsEquationsExp2D():
    iwork = np.zeros((1,), 'i')
    rwork = np.zeros((1,), 'd')
    rwork[0] = 5
    nPoints_x = 51
    nPoints_y = 51
    t = 0.0
    x = np.zeros((nPoints_x, nPoints_y, 3), 'd')
    for i in range(nPoints_x):
        for j in range(nPoints_y):
            x[i, j, 0] = i*(1.0/(nPoints_x-1.0))
            x[i, j, 1] = j*(1.0/(nPoints_y-1.0))
    u = np.zeros(x.shape[:-1], 'd')
    canalyticalSolutions.poissonsEquationExp2D(iwork, rwork, t, x, u)
    assert u[(nPoints_x+1)//2, (nPoints_y+1)//2] == 0.5349653114820003
    if save_plots:
        fig = plt.figure()
        plt.contourf(x[:, :, 0],
                     x[:, :, 1],
                     u)
        plt.title('poissonsEquationExp2D')
        plt.savefig("poissonsEquationExp2D.png")


def test_poissonsEquationExp3D():
    iwork = np.zeros((1,), 'i')
    rwork = np.zeros((1,), 'd')
    rwork[0] = 5
    t = 0.0
    u = np.zeros(x.shape[:-1], 'd')
    canalyticalSolutions.poissonsEquationExp3D(iwork, rwork, t, x, u)
    assert u[midpoint_x_index, midpoint_y_index, midpoint_z_index] == 0.1792429116046883
    if save_plots:
        fig = plt.figure()
        slice = (nPoints_x+1)//2
        plt.contourf(x[slice, :, :, 1],
                     x[slice, :, :, 2],
                     u[slice, :, :])
        plt.title('poissonsEquationExp3D-x/2')
        plt.savefig('poissonsEquationExp3D_xby2.png')
        fig = plt.figure()
        slice = (nPoints_y+1)//2
        plt.contourf(x[:, slice, :, 0],
                     x[:, slice, :, 2],
                     u[:, slice, :])
        plt.title('poissonsEquationExp3D-y/2')
        plt.savefig('poissonsEquationExp3D_yby2.png')
        fig = plt.figure()
        slice = (nPoints_z+1)//2
        plt.contourf(x[:, :, slice, 0],
                     x[:, :, slice, 1],
                     u[:, :, slice])
        plt.title('poissonsEquationExp3D-z/2')
        plt.savefig('poissonsEquationExp3D_zby2.png')


def test_diffusionSin1D():
    iwork = np.zeros((1,), 'i')
    iwork[0] = 1
    rwork = np.zeros((1,), 'd')
    nPoints = 101
    t = 0.0
    x = np.zeros((nPoints, 3), 'd')
    for i in range(nPoints):
        x[i, 0] = i*(1.0/(nPoints-1.0))
    u = np.zeros(x.shape[0], 'd')
    canalyticalSolutions.diffusionSin1D(iwork, rwork, t, x, u)
    # assert u[(nPoints+1)//2] == 2.4788703515727932
    np.testing.assert_approx_equal(u[(nPoints+1)//2], 2.4788703515727932,
    significant=7)
    if save_plots:
        fig = plt.figure()
        plt.plot(x[:, 0],
                 u)
        plt.title('diffusionSin1D')
        plt.savefig('diffusionSin1D')


def test_diffusionSin2D():
    iwork = np.zeros((2,), 'i')
    iwork[0] = 1
    iwork[1] = 1
    rwork = np.zeros((1,), 'd')
    nPoints_x = 51
    nPoints_y = 51
    t = 0.0
    x = np.zeros((nPoints_x, nPoints_y, 3), 'd')
    for i in range(nPoints_x):
        for j in range(nPoints_y):
            x[i, j, 0] = i*(1.0/(nPoints_x-1.0))
            x[i, j, 1] = j*(1.0/(nPoints_y-1.0))
    u = np.zeros(x.shape[:-1], 'd')
    canalyticalSolutions.diffusionSin2D(iwork, rwork, t, x, u)
    # assert u[(nPoints_x+1)//2, (nPoints_y+1)//2] == 9.895915468712143
    np.testing.assert_approx_equal(u[(nPoints_x+1)//2, (nPoints_y+1)//2], 9.895915468712143,significant=7)
    if save_plots:
        fig = plt.figure()
        plt.contourf(x[:, :, 0],
                     x[:, :, 1],
                     u)
        plt.title('diffusionSin2D')
        plt.savefig('diffusionSin2D.png')


def test_diffusionSin3D():
    iwork = np.zeros((3,), 'i')
    iwork[0] = 1
    iwork[1] = 1
    iwork[2] = 1
    rwork = np.zeros((1,), 'd')
    t = 0.0
    u = np.zeros(x.shape[:-1], 'd')
    canalyticalSolutions.diffusionSin3D(iwork, rwork, t, x, u)
    assert u[midpoint_x_index, midpoint_y_index, midpoint_z_index] == -4147.294904126027

def test_STflowSphere_P():
    iwork = np.zeros((1,), 'i')
    rwork = np.zeros((8,), 'd')
    rwork[0] = 0.1
    rwork[1] = 0.0
    rwork[2] = 0.0
    nPoints_x = 31
    nPoints_y = 31
    nPoints_z = 31
    t = 0.0
    dx = 1.0/(nPoints_x-1.0)
    dy = 1.0/(nPoints_y-1.0)
    dz = 1.0/(nPoints_z-1.0)
    rwork[3] = 0.1
    rwork[4] = 0.5
    rwork[5] = 0.5
    rwork[6] = 0.5
    rwork[7] = 1.003e-3
    x = np.zeros((nPoints_x, nPoints_y, nPoints_z, 3), 'd')
    for i in range(nPoints_x):
        for j in range(nPoints_y):
            for k in range(nPoints_z):
                x[i, j, k, 0] = i*dx
                x[i, j, k, 1] = j*dy
                x[i, j, k, 2] = k*dz
    p = np.zeros(x.shape[:-1], 'd')
    canalyticalSolutions.STflowSphere_P(iwork, rwork, t, x, p)
    slice = nPoints_x//2
    if save_plots:
        fig = plt.figure()
        plt.contourf(x[slice, :, :, 1],
                     x[slice, :, :, 2],
                     p[slice, :, :])
        plt.title('STflowSphere_P: yz plane')
        plt.savefig("STflowSphere_P_yz.png")
        fig = plt.figure()
        slice = nPoints_y//2
        plt.contourf(x[:, slice, :, 0],
                     x[:, slice, :, 2],
                     p[:, slice, :])
        plt.title('STflowSphere_P: xz plane')
        plt.savefig("STflowSphere_P_xz.png")
        fig = plt.figure()
        slice = nPoints_z//2
        plt.contourf(x[:, :, slice, 0],
                     x[:, :, slice, 1],
                     p[:, :, slice])
        plt.title('STflowSphere_P: xy plane')
        plt.savefig("STflowSphere_P_xy.png")
    ux = np.zeros(x.shape[:-1], 'd')
    canalyticalSolutions.STflowSphere_Vx(iwork, rwork, t, x, ux)
    # assert ux[midpoint_x_index, midpoint_y_index, midpoint_z_index] == 0.0422649730810374
    np.testing.assert_approx_equal(ux[midpoint_x_index, midpoint_y_index, midpoint_z_index],0.0422649730810374,significant=7)
    uy = np.zeros(x.shape[:-1], 'd')
    canalyticalSolutions.STflowSphere_Vy(iwork, rwork, t, x, uy)
    # assert uy[midpoint_x_index, midpoint_y_index, midpoint_z_index] == -0.009622504486493759
    np.testing.assert_approx_equal(uy[midpoint_x_index, midpoint_y_index, midpoint_z_index],-0.009622504486493759,significant=7)
    uz = np.zeros(x.shape[:-1], 'd')
    canalyticalSolutions.STflowSphere_Vz(iwork, rwork, t, x, uz)
    # assert uz[midpoint_x_index, midpoint_y_index, midpoint_z_index] == -0.009622504486493759
    np.testing.assert_approx_equal(uz[midpoint_x_index, midpoint_y_index, midpoint_z_index],-0.009622504486493759,significant=7)
    if save_plots:
        fig = plt.figure()
        slice = (nPoints_z+1)//2
        plt.quiver(x[:, :, slice, 0],
                   x[:, :, slice, 1],
                   ux[:, :, slice],
        uy[:, :, slice])
        plt.title('STflowSphere : xy plane')
        plt.savefig('STflowSphere_xy.png')
        fig = plt.figure()
        slice = (nPoints_y+1)//2
        plt.quiver(x[:, slice, :, 0],
                   x[:, slice, :, 2],
                   ux[:, slice, :],
                   uz[:, slice, :])
        plt.title('STflowSphere : xz plane')
        plt.savefig('STflowSphere_xz.png')
        fig = plt.figure()
        slice = (nPoints_x+1)//2
        plt.quiver(x[slice, :, :, 1],
                   x[slice, :, :, 2],
                   uy[slice, :, :],
                   uz[slice, :, :])
        plt.title('STflowSphere : yz plane')
        plt.savefig('STflowSphere_yz.png')


def test_PlanePoiseuilleFlow_u():
    iwork = np.zeros((1,), 'i')
    rwork = np.zeros((7,), 'd')
    t = 0.0
    rwork[0] = dy*(nPoints_y-1.0)
    rwork[1] = 0.001
    rwork[2] = -1.0
    # rwork[3]=83.3
    rwork[4] = 0.0
    rwork[5] = 0.0
    rwork[6] = 0.0
    ux = np.zeros(x.shape[:-1], 'd')
    uy = np.zeros(x.shape[:-1], 'd')
    canalyticalSolutions.PlanePoiseuilleFlow_u(iwork, rwork, t, x, ux)
    assert ux[midpoint_x_index, midpoint_y_index, midpoint_z_index] == 124.5674740484429
    slice = nPoints_z//2
    if save_plots:
        fig = plt.figure()
        plt.contourf(x[:, :, slice, 0],
                     x[:, :, slice, 1],
                     ux[:, :, slice])
        plt.title("Plane Poiseuille Flow")
        plt.savefig("PlanePoiseuilleFlow.png")
        fig = plt.figure()
        plt.quiver(x[:, :, slice, 0],
                   x[:, :, slice, 1],
                   ux[:, :, slice],
                   uy[:, :, slice])
        plt.title("Plane Poiseuille Flow Velocity")
        plt.savefig("PlanePoiseuilleFlowQuiver.png")


def test_PoiseuillePipeFlow():
    iwork = np.zeros((1,), 'i')
    rwork = np.zeros((7,), 'd')
    t = 0.0
    dx = 1.0/(nPoints_x-1.0)
    dy = 1.0/(nPoints_y-1.0)
    dz = 1.0/(nPoints_z-1.0)
    rwork[0] = dy*(nPoints_y-1.0)/2.0
    rwork[1] = 0.001
    rwork[2] = -1.0
    # rwork[3]=24.5
    rwork[4] = 0.5
    rwork[5] = 0.5
    rwork[6] = 0.5
    ux = np.zeros(x.shape[:-1], 'd')
    uy = np.zeros(x.shape[:-1], 'd')
    canalyticalSolutions.PoiseuillePipeFlow(iwork, rwork, t, x, ux)
    assert ux[midpoint_x_index, midpoint_y_index, midpoint_z_index] == 62.28373702422145
    slice = nPoints_z//2
    if save_plots:
        fig = plt.figure()
        plt.contourf(x[:, :, slice, 0],
                     x[:, :, slice, 1],
                     ux[:, :, slice])
        plt.title('Poiseuille Pipe Flow')
        plt.savefig("PoiseuillePipeFlow.png")
        fig = plt.figure()
        plt.quiver(x[:, :, slice, 0],
                   x[:, :, slice, 1],
                   ux[:, :, slice],
                   uy[:, :, slice])
        plt.title("PoiseuillePipeFlow")
        plt.savefig("PoiseuillePipeFlowQuiver.png")


def test_PoiseuillePipeFlow_P():
    iwork = np.zeros((1,), 'i')
    rwork = np.zeros((8,), 'd')
    t = 0.0
    dx = 1.0/(nPoints_x-1.0)
    dy = 1.0/(nPoints_y-1.0)
    dz = 1.0/(nPoints_z-1.0)
    rwork[0] = dy*(nPoints_y-1.0)/2.0
    rwork[1] = 0.001
    rwork[2] = -1.0
    # rwork[3]=24.5
    rwork[4] = 1.0
    rwork[5] = 0.5
    rwork[6] = 0.5
    rwork[7] = 0.5
    ux = np.zeros(x.shape[:-1], 'd')
    canalyticalSolutions.PoiseuillePipeFlow_P(iwork, rwork, t, x, ux)
    assert ux[midpoint_x_index, midpoint_y_index, midpoint_z_index] == 0.9705882352941176
    slice = nPoints_z//2
    if save_plots:
        fig = plt.figure()
        plt.contourf(x[:, :, slice, 0],
                     x[:, :, slice, 1],
                     ux[:, :, slice])
        plt.title('Poiseuille Pipe Flow: Pressure')
        plt.savefig("PoiseuillePipeFlow_P.png")
import numpy as np
from pytest import approx
from proteus import equivalent_polynomials as eqp

def test_1D():
    from proteus.Quadrature import GaussEdge
    polyOrders = [1,2,3]
    quadOrderMax = 9
    elements = [[0.,1.],[-2.,2.],[9.,10.]]
    phiList  = [[-1.,-1.],[1.,1.],[-1.,1.],[0.,1.],[1.,0.]]
    for nP in polyOrders:
        for qO in range(nP,quadOrderMax):
            quad = GaussEdge(order=qO)
            gf = eqp.Simplex(nSpace=1, nP=nP, nQ=len(quad.points))
            for phi in phiList:
                for e in elements:
                    dV = abs(e[1]-e[0])
                    if phi[0]*phi[1] < 0.0:
                        theta = -phi[0]/(phi[1] - phi[0])
                        x_0 = (1-theta)*e[0] + theta*e[1]
                        int_H_exact = abs(e[1] - x_0)
                        int_ImH_exact = abs(x_0 - e[0])
                        int_D_exact = 1.0
                        if phi[0] > 0.0:
                            tmp=int_H_exact
                            int_H_exact = int_ImH_exact
                            int_ImH_exact = int_H_exact
                    elif phi[0] < 0.0:
                        int_H_exact = 0.0
                        int_ImH_exact = dV
                        int_D_exact = 0.0
                    elif phi[0] > 0.0:
                        int_H_exact = dV
                        int_ImH_exact = 0.0
                        int_D_exact = 0.0
                    elif phi[0] == 0.0:
                        if phi[1] > 0.0:
                            int_H_exact = dV
                            int_ImH_exact = 0.0
                            int_D_exact = 0.0
                        if phi[1] < 0.0:
                            int_H_exact = 0.0
                            int_ImH_exact = dV
                            int_D_exact = 1.0
                    elif phi[1] == 0.0:
                        if phi[0] > 0.0:
                            int_H_exact = dV
                            int_ImH_exact = 0.0
                            int_D_exact = 0.0
                        if phi[0] < 0.0:
                            int_H_exact = 0.0
                            int_ImH_exact = dV
                            int_D_exact = 1.0
                            
                    gf.calculate(np.array(phi),
                                 np.array([[e[0],0.,0.],
                                           [e[1],0.,0.]]).transpose(),
                                 np.array(quad.points))
                    int_H=0.0
                    int_ImH=0.0
                    int_D=0.0
                    for k in range(len(quad.points)):
                        gf.set_quad(k)
                        int_H += quad.weights[k]*gf.H*dV
                        int_ImH += quad.weights[k]*gf.ImH*dV
                        int_D += quad.weights[k]*gf.D*dV
                    assert(int_H == approx(int_H_exact,1e-16))
                    assert(int_ImH == approx(int_ImH_exact,1e-15))
                    assert(int_D == approx(int_D_exact,1e-15))

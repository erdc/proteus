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
                    dV = e[1]-e[0]
                    assert(dV > 0)
                    if phi[0]*phi[1] < 0.0:
                        theta = -phi[0]/(phi[1] - phi[0])
                        x_0 = (1-theta)*e[0] + theta*e[1]
                        int_H_exact = abs(e[1] - x_0)
                        int_ImH_exact = abs(x_0 - e[0])
                        int_D_exact = 1.0
                        if phi[0] > 0.0:
                            tmp=int_H_exact
                            int_H_exact = int_ImH_exact
                            int_ImH_exact = tmp
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
                                           [e[1],0.,0.]]),
                                 np.array(quad.points))
                    int_H=0.0
                    int_ImH=0.0
                    int_D=0.0
                    for k in range(len(quad.points)):
                        gf.set_quad(k)
                        int_H += quad.weights[k]*gf.H*dV
                        int_ImH += quad.weights[k]*gf.ImH*dV
                        int_D += quad.weights[k]*gf.D*dV
                    assert(int_H == approx(int_H_exact,1e-15))
                    assert(int_ImH == approx(int_ImH_exact,1e-15))
                    assert(int_D == approx(int_D_exact,1e-15))

def test_2D():
    from proteus.Quadrature import GaussTriangle
    polyOrders = [1,2,3]
    quadOrderMax = 6
    elements = [
        np.array([[0.,0.],[0.,1.],[1.,0.]])
    ]
    phiList  = [[-1.,-1.,-1.],
                [1.,1.,1.],
                [-1.,1.,1.],
                [0.,1.,1.],
                [1.,0.,0.]]
    for nP in polyOrders:
        for qO in range(nP,quadOrderMax):
            quad = GaussTriangle(order=qO)
            gf = eqp.Simplex(nSpace=2, nP=nP, nQ=len(quad.points))
            for phi in phiList:
                for e in elements:
                    b_0 = e[2,:] - e[0,:]
                    b_1 = e[1,:] - e[0,:]
                    Jac = np.array([b_0,b_1]).transpose()
                    dV = abs(np.linalg.det(Jac))
                    area = dV/2.0
                    if phi[0]*phi[1] < 0.0 and phi[0]*phi[2] < 0.0:
                        theta0 = -phi[0]/(phi[1] - phi[0])
                        theta1 = -phi[0]/(phi[2] - phi[0])
                        x_0 = (1-theta0)*e[0,:] + theta0*e[1,:]
                        x_1 = (1-theta1)*e[0,:] + theta1*e[2,:]
                        b_0 = x_1 - e[0,:]
                        b_1 = x_0 - e[0,:]
                        Jac_0 = np.array([b_0,b_1]).transpose()
                        int_ImH_exact = np.linalg.det(Jac_0)/2.0
                        int_H_exact = area - int_ImH_exact
                        int_D_exact = np.linalg.norm(x_1 - x_0)
                        if phi[0] > 0.0:
                            tmp=int_H_exact
                            int_H_exact = int_ImH_exact
                            int_ImH_exact = tmp
                    elif phi[0] > 0.0 and phi[1] == 0.0 and phi[2] == 0.0:
                        int_H_exact = area
                        int_ImH_exact = 0.0
                        int_D_exact = np.linalg.norm(e[1,:] - e[2,:])
                    elif phi[0] < 0.0:
                        int_H_exact = 0.0
                        int_ImH_exact = area
                        int_D_exact = 0.0
                    elif phi[0] > 0.0:
                        int_H_exact = area
                        int_ImH_exact = 0.0
                        int_D_exact = 0.0
                    elif phi[0] == 0.0:
                        if phi[1] > 0.0:
                            int_H_exact = area
                            int_ImH_exact = 0.0
                            int_D_exact = 0.0
                        if phi[1] < 0.0:
                            int_H_exact = 0.0
                            int_ImH_exact = area
                            int_D_exact = 1.0
                    elif phi[1] == 0.0:
                        if phi[0] > 0.0:
                            int_H_exact = area
                            int_ImH_exact = 0.0
                            int_D_exact = 0.0
                        if phi[0] < 0.0:
                            int_H_exact = 0.0
                            int_ImH_exact = area
                            int_D_exact = 1.0
                            
                    gf.calculate(np.array(phi),
                                 np.array([[e[0,0],e[0,1],0.],
                                           [e[1,0],e[1,1],0.],
                                           [e[2,0],e[2,1],0.]]),
                                 np.array(quad.points))
                    int_H=0.0
                    int_ImH=0.0
                    int_D=0.0
                    for k in range(len(quad.points)):
                        gf.set_quad(k)
                        int_H += quad.weights[k]*gf.H*dV
                        int_ImH += quad.weights[k]*gf.ImH*dV
                        int_D += quad.weights[k]*gf.D*dV
                    assert(int_H == approx(int_H_exact,1e-15))
                    assert(int_ImH == approx(int_ImH_exact,1e-15))
                    assert(int_D == approx(int_D_exact,1e-15))

def test_3D():
    from proteus.Quadrature import GaussTetrahedron
    polyOrders = [1,2,3]
    quadOrderMax = 8
    elements = [
        np.array([[0.,0.,0.],
                  [1.,0.,0.],
                  [0.,1.,0.],
                  [0.,0.,1.]]),
       np.array([[659.82702468, 583.95327591, 217.72111178],
                 [312.90104741, 109.12205319,  82.38957441],
                 [371.40050898, 837.09283397, 117.87194133],
                 [666.41881358, 107.49166421, 419.82032267]]),
       np.array([[1.,1.,0.5],
                 [0.,0.,1.25],
                 [0.,1.,0.5],
                 [0.,0.,0.5,]])
    ]
    phiList=  [
        [-1.,-1.,-1.,-1.],
        [1.,1.,1.,1.],
        [-1.,1.,1.,1.],
        [1.,-1.,-1.,-1.],
        [-1.,1.,-1.,-1.],
        [1.,-1.,1.,1.],
        [-1.,1.,-1.,-1.],
        [1.,1.,-1.,1.],
        [-1.,-1.,1.,-1.],
        [1.,1.,1.,-1.],
        [-1.,-1.,-1.,1.],
        [-1.,1.,1.,-1.],
        [1.,-1.,-1.,1.],
        [-1.,1.,1.,-1.],
        [-2.,3.,4.,5.],
        [3.,-2.,-4.,-5.],
        [-2.,3.,-4.,-5.],
        [3.,-2.,4.,5.],
        [-2.,3.,-4.,-5.],
        [3.,4.,-2.,5.],
        [-2.,-3.,4.,-5.],
        [3.,4.,5.,-2.],
        [-2.,-3.,-4.,5.],
        [0.,1.,1.,1.],
        [9.52267983, -348.3143481, -10.31124012, -2.73674249],
        [-9.52267983, 348.3143481, -10.31124012, -2.73674249],
        [-9.52267983, -348.3143481, 10.31124012, -2.73674249],
        [-9.52267983, -348.3143481, -10.31124012, 2.73674249],
        [0.5, -0.25, 0.5, 0.5]]
    for nP in polyOrders:
        for qO in range(nP,quadOrderMax):
            quad = GaussTetrahedron(order=qO)
            print("Simplex in ",3,nP,len(quad.points))
            gf = eqp.Simplex(nSpace=3, nP=nP, nQ=len(quad.points))
            for phi in phiList:
                for e in elements:
                    b_0 = e[3,:] - e[0,:]
                    b_1 = e[2,:] - e[0,:]
                    b_2 = e[1,:] - e[0,:]
                    Jac = np.array([b_0, b_1, b_2]).transpose()
                    dV = abs(np.linalg.det(Jac))
                    volume = dV/6.0
                    if phi[0] < 0.0 and phi[1] > 0.0 and phi[2] > 0.0 and phi[3] < 0.0:
                        print("====quad cut===== ",phi)
                        print("qO ",qO)
                        print("nP ",nP)
                        theta01 = 0.5 - 0.5*(phi[1]+phi[0])/(phi[1] - phi[0])
                        theta02 = 0.5 - 0.5*(phi[2]+phi[0])/(phi[2] - phi[0])
                        theta31 = 0.5 - 0.5*(phi[1]+phi[3])/(phi[1] - phi[3])
                        theta32 = 0.5 - 0.5*(phi[2]+phi[3])/(phi[2] - phi[3])
                        x_01 = (1-theta01)*e[0,:] + theta01*e[1,:]
                        x_02 = (1-theta02)*e[0,:] + theta02*e[2,:]
                        x_31 = (1-theta31)*e[3,:] + theta31*e[1,:]
                        x_32 = (1-theta32)*e[3,:] + theta32*e[2,:]
                        b_0 = x_01 - e[1,:]
                        b_1 = x_02 - e[1,:]
                        b_2 = x_31 - e[1,:]
                        Jac_0 = np.array([b_0,b_1,b_2]).transpose()
                        int_ImH_exact = volume
                        int_H_exact = 0.0
                        int_D_exact = 0.0
                        int_H_exact += abs(np.linalg.det(Jac_0))/6.0
                        int_ImH_exact -= abs(np.linalg.det(Jac_0))/6.0
                        int_D_exact += 0.5*np.linalg.norm(np.cross(x_02 - x_01,
                                                                   x_32 - x_01))
                        b_0 = x_02 - e[2,:]
                        b_1 = x_32 - e[2,:]
                        b_2 = x_31 - e[2,:]
                        Jac_0 = np.array([b_0,b_1,b_2]).transpose()
                        int_H_exact += abs(np.linalg.det(Jac_0))/6.0
                        int_ImH_exact -= abs(np.linalg.det(Jac_0))/6.0
                        int_D_exact += 0.5*np.linalg.norm(np.cross(x_32 - x_01,
                                                                   x_31 - x_01))
                        b_0 = e[1,:] - x_02
                        b_1 = e[2,:] - x_02
                        b_2 = x_31   - x_02
                        Jac_0 = np.array([b_0,b_1,b_2]).transpose()
                        int_H_exact += abs(np.linalg.det(Jac_0))/6.0
                        int_ImH_exact -= abs(np.linalg.det(Jac_0))/6.0
                    elif phi[0] > 0.0 and phi[1] < 0.0 and phi[2] < 0.0 and phi[3] > 0.0:
                        print("====quad cut=====  inside out",phi)
                        print("qO ",qO)
                        print("nP ",nP)
                        theta01 = 0.5 - 0.5*(phi[1]+phi[0])/(phi[1] - phi[0])
                        theta02 = 0.5 - 0.5*(phi[2]+phi[0])/(phi[2] - phi[0])
                        theta31 = 0.5 - 0.5*(phi[1]+phi[3])/(phi[1] - phi[3])
                        theta32 = 0.5 - 0.5*(phi[2]+phi[3])/(phi[2] - phi[3])
                        x_01 = (1-theta01)*e[0,:] + theta01*e[1,:]
                        x_02 = (1-theta02)*e[0,:] + theta02*e[2,:]
                        x_31 = (1-theta31)*e[3,:] + theta31*e[1,:]
                        x_32 = (1-theta32)*e[3,:] + theta32*e[2,:]
                        b_0 = x_01 - e[1,:]
                        b_1 = x_02 - e[1,:]
                        b_2 = x_31 - e[1,:]
                        Jac_0 = np.array([b_0,b_1,b_2]).transpose()
                        int_ImH_exact = 0.0
                        int_H_exact = volume
                        int_D_exact = 0.0
                        int_H_exact -= abs(np.linalg.det(Jac_0))/6.0
                        int_ImH_exact += abs(np.linalg.det(Jac_0))/6.0
                        int_D_exact += 0.5*np.linalg.norm(np.cross(x_02 - x_01,
                                                                   x_32 - x_01))
                        b_0 = x_02 - e[2,:]
                        b_1 = x_32 - e[2,:]
                        b_2 = x_31 - e[2,:]
                        Jac_0 = np.array([b_0,b_1,b_2]).transpose()
                        int_H_exact -= abs(np.linalg.det(Jac_0))/6.0
                        int_ImH_exact += abs(np.linalg.det(Jac_0))/6.0
                        int_D_exact += 0.5*np.linalg.norm(np.cross(x_32 - x_01,
                                                                   x_31 - x_01))
                        b_0 = e[1,:] - x_02
                        b_1 = e[2,:] - x_02
                        b_2 = x_31   - x_02
                        Jac_0 = np.array([b_0,b_1,b_2]).transpose()
                        int_H_exact -= abs(np.linalg.det(Jac_0))/6.0
                        int_ImH_exact += abs(np.linalg.det(Jac_0))/6.0
                    elif phi[0]*phi[1] < 0.0 and phi[0]*phi[2] < 0.0 and phi[0]*phi[3] < 0.0:
                        theta0 = -phi[0]/(phi[1] - phi[0])
                        theta1 = -phi[0]/(phi[2] - phi[0])
                        theta2 = -phi[0]/(phi[3] - phi[0])
                        x_0 = (1-theta0)*e[0,:] + theta0*e[1,:]
                        x_1 = (1-theta1)*e[0,:] + theta1*e[2,:]
                        x_2 = (1-theta2)*e[0,:] + theta2*e[3,:]
                        b_0 = x_2 - e[0,:]
                        b_1 = x_1 - e[0,:]
                        b_2 = x_0 - e[0,:]
                        Jac_0 = np.array([b_0,b_1,b_2]).transpose()
                        int_ImH_exact = abs(np.linalg.det(Jac_0))/6.0
                        int_H_exact = volume - int_ImH_exact
                        int_D_exact = 0.5*np.linalg.norm(np.cross(x_2 - x_0, x_1 - x_0))
                        if phi[0] > 0.0:
                            tmp=int_H_exact
                            int_H_exact = int_ImH_exact
                            int_ImH_exact = tmp
                    elif phi[0]*phi[1] < 0.0 and phi[1]*phi[2] < 0.0 and phi[1]*phi[3] < 0.0:
                        theta0 = -phi[1]/(phi[0] - phi[1])
                        theta1 = -phi[1]/(phi[2] - phi[1])
                        theta2 = -phi[1]/(phi[3] - phi[1])
                        x_0 = (1-theta0)*e[1,:] + theta0*e[0,:]
                        x_1 = (1-theta1)*e[1,:] + theta1*e[2,:]
                        x_2 = (1-theta2)*e[1,:] + theta2*e[3,:]
                        b_0 = x_2 - e[1,:]
                        b_1 = x_1 - e[1,:]
                        b_2 = x_0 - e[1,:]
                        Jac_0 = np.array([b_0,b_1,b_2]).transpose()
                        int_ImH_exact = abs(np.linalg.det(Jac_0))/6.0
                        int_H_exact = volume - int_ImH_exact
                        int_D_exact = 0.5*np.linalg.norm(np.cross(x_2 - x_0, x_1 - x_0))
                        if phi[1] > 0.0:
                            tmp=int_H_exact
                            int_H_exact = int_ImH_exact
                            int_ImH_exact = tmp
                    elif phi[0]*phi[2] < 0.0 and phi[1]*phi[2] < 0.0 and phi[3]*phi[2] < 0.0:
                        theta0 = -phi[2]/(phi[0] - phi[2])
                        theta1 = -phi[2]/(phi[1] - phi[2])
                        theta2 = -phi[2]/(phi[3] - phi[2])
                        x_0 = (1-theta0)*e[2,:] + theta0*e[0,:]
                        x_1 = (1-theta1)*e[2,:] + theta1*e[1,:]
                        x_2 = (1-theta2)*e[2,:] + theta2*e[3,:]
                        b_0 = x_2 - e[2,:]
                        b_1 = x_1 - e[2,:]
                        b_2 = x_0 - e[2,:]
                        Jac_0 = np.array([b_0,b_1,b_2]).transpose()
                        int_ImH_exact = abs(np.linalg.det(Jac_0))/6.0
                        int_H_exact = volume - int_ImH_exact
                        int_D_exact = 0.5*np.linalg.norm(np.cross(x_2 - x_0, x_1 - x_0))
                        if phi[2] > 0.0:
                            tmp=int_H_exact
                            int_H_exact = int_ImH_exact
                            int_ImH_exact = tmp
                    elif phi[0]*phi[3] < 0.0 and phi[1]*phi[3] < 0.0 and phi[2]*phi[3] < 0.0:
                        theta0 = -phi[3]/(phi[0] - phi[3])
                        theta1 = -phi[3]/(phi[1] - phi[3])
                        theta2 = -phi[3]/(phi[2] - phi[3])
                        x_0 = (1-theta0)*e[3,:] + theta0*e[0,:]
                        x_1 = (1-theta1)*e[3,:] + theta1*e[1,:]
                        x_2 = (1-theta2)*e[3,:] + theta2*e[2,:]
                        b_0 = x_2 - e[3,:]
                        b_1 = x_1 - e[3,:]
                        b_2 = x_0 - e[3,:]
                        Jac_0 = np.array([b_0,b_1,b_2]).transpose()
                        int_ImH_exact = abs(np.linalg.det(Jac_0))/6.0
                        int_H_exact = volume - int_ImH_exact
                        int_D_exact = 0.5*np.linalg.norm(np.cross(x_2 - x_0, x_1 - x_0))
                        if phi[3] > 0.0:
                            tmp=int_H_exact
                            int_H_exact = int_ImH_exact
                            int_ImH_exact = tmp
                    elif phi[0] < 0.0:
                        int_H_exact = 0.0
                        int_ImH_exact = volume
                        int_D_exact = 0.0
                    elif phi[0] > 0.0:
                        int_H_exact = volume
                        int_ImH_exact = 0.0
                        int_D_exact = 0.0
                    elif phi[0] == 0.0:
                        if phi[1] > 0.0:
                            int_H_exact = volume
                            int_ImH_exact = 0.0
                            int_D_exact = 0.0
                        if phi[1] < 0.0:
                            int_H_exact = 0.0
                            int_ImH_exact = volume
                            int_D_exact = 1.0
                    elif phi[1] == 0.0:
                        if phi[0] > 0.0:
                            int_H_exact = volume
                            int_ImH_exact = 0.0
                            int_D_exact = 0.0
                        if phi[0] < 0.0:
                            int_H_exact = 0.0
                            int_ImH_exact = volume
                            int_D_exact = 1.0
                    phi_in = np.array(phi)
                    print("phi in ", phi_in) 
                    nodesIn = e
                    print("nodes in ",nodesIn)
                    print("quad in ",np.array(quad.points))
                    gf.calculate(phi_in,
                                 nodesIn,
                                 np.array(quad.points))
                    int_H=0.0
                    int_ImH=0.0
                    int_D=0.0
                    for k in range(len(quad.points)):
                        gf.set_quad(k)
                        print("k ", k, gf.H, gf.ImH, gf.D)
                        int_H += quad.weights[k]*gf.H*dV
                        int_ImH += quad.weights[k]*gf.ImH*dV
                        int_D += quad.weights[k]*gf.D*dV
                    print("dV ", dV)
                    print(int_H, int_H_exact)
                    print(int_ImH, int_ImH_exact)
                    print(int_D, int_D_exact)
                    assert(int_H == approx(int_H_exact,1e-12,1e-12))
                    assert(int_ImH == approx(int_ImH_exact,1e-12,1e-12))
                    assert(int_D == approx(int_D_exact,1e-11,1e-11))
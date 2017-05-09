from proteus import Comm, Profiling
import numpy as np
import numpy.testing as npt
import unittest
import random
from math import cos,sin,cosh,sinh,pi,tanh,log,atan2,asin,acos,radians
import sys,os
import logging
import pytest

import cython

comm = Comm.init()
Profiling.procID = comm.rank()
def getpath():
    path = sys.path[0]+'/'
    return path

def remove_files(filenames):
    ''' delete files in filenames list '''
    for f in filenames:
        if os.path.exists(f):
            try:
                os.remove(f)
            except OSError, e:
                print ("Error: %s - %s" %(e.filename,e.strerror))
        else:
            pass

Profiling.logEvent("Testing Bodydynamics")

class TestAuxFunctions(unittest.TestCase):
    def testForwardEler(self):
        from proteus.mprans.BodyDynamics import forward_euler
        p0 = 5.0
        v0 = 2.0
        a  = 0.5
        dt = 0.01
        v = v0 + a*dt
        p = p0 + v*dt
        p1, v1 = forward_euler(p0, v0, a, dt)
        self.assertTrue(p==p1)
        self.assertTrue(v==v1)

    def testRungeKutta(self):
        from proteus.mprans.BodyDynamics import runge_kutta
        u0 = 5.0
        v0 = 2.0
        dt = 0.01
        substeps = 20
        F = 150.0
        K = 100000.0
        C = 1000.0
        m = 50.0
        a0 = (F - C*v0 - K*u0) / m
        velCheck = True
        # Loop through substeps
        for i in range(substeps):
            # 1st step
            u1 = u0
            v1 = v0
            a1 = (F - C*v1 - K*u1) / m
            # 2nd step
            u2 = u0 + v1*dt/2.
            v2 = v0 + a1*dt/2.
            a2 = (F - C*v2 - K*u2) / m
            # 3rd step
            u3 = u0 + v2*dt/2.
            v3 = v0 + a2*dt/2.
            a3 = (F - C*v3 - K*u3) / m
            # 4th step
            u4 = u0 + v3*dt
            v4 = v0 + a3*dt
            a4 = (F - C*v4 - K*u4) / m
            # Final step
            u = u0 + dt/6.*(v1+2*v2+2*v3+v4)
            v = v0 + dt/6.*(a1+2*a2+2*a3+a4)
            a = (F - C*v - K*u) / m
            # Check for friction module, when velocity is reversed
            if velCheck:
                if v0*v < 0.0:
                    break
        uc, vc, ac = runge_kutta(u0, v0, a0, dt, substeps, F, K, C, m, velCheck)
        self.assertTrue(u==uc)
        self.assertTrue(v==vc)
        self.assertTrue(a==ac)

    def testEulerAngles(self):
        from proteus.mprans.BodyDynamics import getEulerAngles
        coords = np.array([ [1.,0.,0.],
                            [0.,1.,0.],
                            [0.,0.,1.] ])     
        anglex = radians(5.0)
        angley = radians(10.)
        anglez = radians(15.)
        # rotation around x-axis 
        rotx = np.array([  [  1.0,   0.0,          0.0 ],
                           [  0.0,   cos(anglex),  sin(anglex) ],
                           [  0.0,  -sin(anglex),  cos(anglex) ] ])       
        # rotation around y-axis 
        roty = np.array([  [  cos(angley),  0.0,  -sin(angley) ],
                           [  0.0,          1.0,   0.0 ],
                           [  sin(angley),  0.0,   cos(angley) ] ])
        # rotation around z-axis 
        rotz = np.array([  [  cos(anglez),  sin(anglez),  0.0 ],
                           [ -sin(anglez),  cos(anglez),  0.0 ],
                           [  0.0,         0.0,           1.0 ] ])
        # total 3D rotation
        rot3D = np.dot(np.dot(rotx,roty),rotz)
        coordinates = np.dot( rot3D, coords) 
        ax = atan2(coordinates[1,2], coordinates[2,2])
        ay = asin(-coordinates[0,2])
        az = atan2(coordinates[0,1], coordinates[0,0])
        angles = [ax, ay, az]
        rotation = getEulerAngles(coordinates)
        self.assertTrue(angles==rotation)

class TestRigidBody(unittest.TestCase):

    def testCalculateInit(self):
        from proteus import Domain  
        from proteus.mprans import SpatialTools as st
        from proteus.mprans import BodyDynamics as bd   
        domain = Domain.PlanarStraightLineGraphDomain()
        pos = np.array([1.0, 1.0, 0.0])
        rot = np.array([ [1.,0.,0.], [0.,1.,0.], [0.,0.,1.] ])  
        dim=(0.3, 0.385)
        caisson = st.Rectangle(domain, dim=dim, coords=[pos[0], pos[1]])
        caisson2D = bd.RigidBody(shape=caisson)
        st.assembleDomain(domain)      
        c2d = caisson2D
        nd = 2
        angle = bd.getEulerAngles(rot)

        c2d.calculate_init()

        npt.assert_equal(c2d.position, pos)
        npt.assert_equal(c2d.last_position, pos)
        npt.assert_equal(c2d.rotation, rot)
        npt.assert_equal(c2d.last_rotation, rot)
        npt.assert_equal(c2d.rotation_euler, angle)
        npt.assert_equal(c2d.last_rotation_euler, angle)

    def testStoreLastValues(self):
        from proteus import Domain  
        from proteus.mprans import SpatialTools as st
        from proteus.mprans import BodyDynamics as bd   
        domain = Domain.PlanarStraightLineGraphDomain()
        pos = np.array([1.0, 1.0, 0.0])
        rot = np.array([ [1.,0.,0.], [0.,1.,0.], [0.,0.,1.] ])  
        dim=(0.3, 0.385)
        caisson = st.Rectangle(domain, dim=dim, coords=[pos[0], pos[1]])
        caisson2D = bd.RigidBody(shape=caisson)
        st.assembleDomain(domain)      
        c2d = caisson2D
        nd = 2
        anglez = radians(5.0)
        rotz = np.array([  [  cos(anglez),  sin(anglez),  0.0 ],
                           [ -sin(anglez),  cos(anglez),  0.0 ],
                           [  0.0,         0.0,           1.0 ] ])

        c2d.calculate_init()

        # Imposing values on rigid body's parameters
        c2d.position = np.array([2.0, 3.0, 0.0])
        c2d.velocity = np.array([0.5, -0.3, 0.0])
        c2d.acceleration = np.array([0.1, 0.2, 0.0])
        c2d.rotation = rotz
        c2d.rotation_euler =  bd.getEulerAngles(rotz)
        c2d.ang_disp = np.array([0.0, 0.0, 0.001])
        c2d.ang_vel = np.array([0.0, 0.0, -0.003])
        c2d.ang_acc = np.array([0.0, 0.0, 0.005])
        c2d.F = np.array([100.0, -300.0, 0.0])
        c2d.M = np.array([0.0, 0.0, 50.0])
        c2d.pivot = c2d.barycenter
        c2d.ux = 0.5
        c2d.uy = 0.05
        c2d.uz = 0.0

        c2d._store_last_values()
        
        npt.assert_equal(c2d.position, c2d.last_position)
        npt.assert_equal(c2d.velocity, c2d.last_velocity)
        npt.assert_equal(c2d.acceleration, c2d.last_acceleration)
        npt.assert_equal(c2d.rotation, c2d.last_rotation)
        npt.assert_equal(c2d.rotation_euler, c2d.last_rotation_euler)
        npt.assert_equal(c2d.ang_disp, c2d.last_ang_disp)
        npt.assert_equal(c2d.ang_vel, c2d.last_ang_vel)
        npt.assert_equal(c2d.ang_acc, c2d.last_ang_acc)
        npt.assert_equal(c2d.F, c2d.last_F)
        npt.assert_equal(c2d.M, c2d.last_M)
        npt.assert_equal(c2d.pivot, c2d.last_pivot)
        npt.assert_equal(c2d.ux, c2d.last_ux)
        npt.assert_equal(c2d.uy, c2d.last_uy)
        npt.assert_equal(c2d.uz, c2d.last_uz)

    def testGetInertia(self):
        from proteus import Domain  
        from proteus import SpatialTools as st    
        from proteus.mprans import SpatialTools as mst
        from proteus.mprans import BodyDynamics as bd   

        # 2d case
        domain2D = Domain.PlanarStraightLineGraphDomain()
        pos = np.array([1.0, 1.0, 0.0])
        dim = (0.3, 0.385)
        caisson = mst.Rectangle(domain=domain2D, dim=dim, coords=[pos[0], pos[1]])
        caisson2D = bd.RigidBody(shape=caisson)
        mst.assembleDomain(domain2D)      
        c2d = caisson2D
        c2d.nd = 2
        Ix = ((dim[0]**2)/12.)
        Iy = ((dim[1]**2)/12.)
        It = Ix + Iy 
        mass = 50.0
        I1 = It*mass
        c2d.mass = mass
        c2d.It = It
        I2 = c2d.getInertia()
        npt.assert_equal(I1, I2)

        # 3d case
        pos, dim, caisson = [], [], []
        domain3D = Domain.PiecewiseLinearComplexDomain()
        pos = np.array([0.0, 0.0, 0.0])
        dim = (0.3, 0.385, 0.4)
        angle = radians(10.)
        Ix = ((dim[0]**2)/12.)
        Iy = ((dim[1]**2)/12.)
        Iz = ((dim[2]**2)/12.)
        It = np.array([ [Iz + Iy, 0.0, 0.0],
                        [0.0, Ix + Iz, 0.0],
                        [0.0, 0.0, Ix + Iy] ])
        mass = 50.0
        # rotational axis and versors
        ax, ay, az = axis = np.array([1., 1., 1.])
        mod = np.sqrt((ax**2)+(ay**2)+(az**2))
        axis_ver = (axis/mod)
        # shape
        caisson = mst.Cuboid(domain=domain3D, dim=dim, coords=[pos[0], pos[1], pos[2]])
        caisson.rotate(rot=angle, axis=axis)
        coords_system = caisson.coords_system
        caisson3D = bd.RigidBody(shape=caisson)
        mst.assembleDomain(domain3D)      
        # rotational operator for inertia calculation on an arbitrary axis
        rotx, roty, rotz = np.dot(axis_ver, np.linalg.inv(coords_system))
        rot = np.array([ [(rotx**2), rotx*roty, rotx*rotz],
                         [roty*rotx, (roty**2), roty*rotz],
                         [rotz*rotx, rotz*roty, (rotz**2)] ])
        #inertia = np.einsum('ij,ij->', mass*It, rot)
        inertia = np.sum(mass*It*rot)
        # testing from BodyDynamics
        c3d = caisson3D
        c3d.nd = 3
        c3d.mass = mass
        c3d.It = It
        c3d.coords_system = coords_system
        I = c3d.getInertia(vec=axis)
        npt.assert_equal(I, inertia)


    def testGetAcceleration(self):
        from proteus import Domain  
        from proteus.mprans import SpatialTools as st
        from proteus.mprans import BodyDynamics as bd   
        domain = Domain.PlanarStraightLineGraphDomain()
        pos = np.array([1.0, 1.0, 0.0])
        rot = np.array([ [1.,0.,0.], [0.,1.,0.], [0.,0.,1.] ])  
        dim=(0.3, 0.385)
        caisson = st.Rectangle(domain, dim=dim, coords=[pos[0], pos[1]])
        caisson2D = bd.RigidBody(shape=caisson)
        st.assembleDomain(domain)      
        c2d = caisson2D

        F = np.array([100.0, -200.0, 10.0])
        mass = 150.0
        acceleration = F/mass

        c2d.F = F
        c2d.mass = mass
        c2d.acceleration = c2d.getAcceleration()
        npt.assert_equal(c2d.acceleration, acceleration)

    def testGetAngularAcceleration(self):
        from proteus import Domain  
        from proteus.mprans import SpatialTools as st
        from proteus.mprans import BodyDynamics as bd   
        domain = Domain.PlanarStraightLineGraphDomain()
        pos = np.array([1.0, 1.0, 0.0])
        rot = np.array([ [1.,0.,0.], [0.,1.,0.], [0.,0.,1.] ])    
        dim=(0.3, 0.385)
        caisson = st.Rectangle(domain, dim=dim, coords=[pos[0], pos[1]])
        caisson2D = bd.RigidBody(shape=caisson)
        st.assembleDomain(domain)      
        c2d = caisson2D

        F = np.array([100.0, -200.0, 10.0])
        mass = 150.0
        acceleration = F/mass

        c2d.F = F
        c2d.mass = mass
        c2d.acceleration = c2d.getAcceleration()
        npt.assert_equal(c2d.acceleration, acceleration)

    def testGetDisplacement(self):
        from proteus import Domain  
        from proteus import SpatialTools as st    
        from proteus.mprans import SpatialTools as mst
        from proteus.mprans import BodyDynamics as bd   

        # parameters
        domain2D = Domain.PlanarStraightLineGraphDomain()
        pos = np.array([1.0, 1.0, 0.0])
        dim = (0.3, 0.385)
        caisson = mst.Rectangle(domain=domain2D, dim=dim, coords=[pos[0], pos[1]])
        caisson2D = bd.RigidBody(shape=caisson)
        mst.assembleDomain(domain2D)      
        c2d = caisson2D
        c2d.nd = 2
        Ix = ((dim[0]**2)/12.)
        Iy = ((dim[1]**2)/12.)
        It = Ix + Iy 
        mass = 50.0
        c2d.mass, c2d.It = mass, It
        Kx, Ky, Kz = K = [200000., 250000., 0.0]
        Cx, Cy, Cz = C = [2000., 4000., 0.0]
        c2d.Kx, c2d.Ky, c2d.Kz, c2d.Cx, c2d.Cy, c2d.Cz =  Kx, Ky, Kz, Cx, Cy, Cz
        init_barycenter, last_position, last_velocity = pos, pos*2.0, np.array([0.05, 0.07, 0.0])
        c2d.init_barycenter, c2d.last_position, c2d.last_velocity = init_barycenter, last_position, last_velocity 
        Fx, Fy, Fz = F = np.array([200., 300., 0.0])
        dt, substeps = 0.001, 50  
        dt_sub = c2d.dt_sub = float(dt/substeps)
        c2d.F, c2d.substeps, c2d.dt, c2d.acceleration = F, substeps, dt, c2d.getAcceleration()

        # runge-kutta calculation
        c2d.scheme = 'Runge_Kutta'
        u0 = ux0, uy0, uz0 = last_position - init_barycenter
        v0 = vx0, vy0, vz0 = last_velocity
        a0 = ax0, ay0, az0 = (F - C*v0 - K*u0) / mass
        for ii in range(substeps):
            ux, vx, ax = bd.runge_kutta(ux0, vx0, ax0, dt_sub, substeps, Fx, Kx, Cx, mass, False)
            uy, vy, ay = bd.runge_kutta(uy0, vy0, ay0, dt_sub, substeps, Fy, Ky, Cy, mass, False)
            uz, vz, az = bd.runge_kutta(uz0, vz0, az0, dt_sub, substeps, Fz, Kz, Cz, mass, False)
        h = hx, hy, hz = [ux-ux0, uy-uy0, uz-uz0]
        c2d.h = c2d.getDisplacement(dt)
        npt.assert_equal(c2d.h, h)

#    def testImposeSinusoidalMotion(self):
#        from proteus import Domain  
#        from proteus import SpatialTools as st    
#        from proteus.mprans import SpatialTools as mst
#        from proteus.mprans import BodyDynamics as bd   
#
#        # parameters
#        domain2D = Domain.PlanarStraightLineGraphDomain()
#        pos = np.array([1.0, 1.0, 0.0])
#        dim = (0.3, 0.385)
#        caisson = mst.Rectangle(domain=domain2D, dim=dim, coords=[pos[0], pos[1]])
#        caisson2D = bd.RigidBody(shape=caisson)
#        mst.assembleDomain(domain2D)      
#        c2d = caisson2D
#        c2d.nd = 2
#        c2d.attachModel(model='twp', ar='auxiliaryVariables')
#        t = c2d.model.stepController.t_model_last = 5.25
#        # transl and rot amplitudes, transl and rot periods
#        At, Tt, Ar, Tr = c2d.At, c2d.Tt, c2d.Ar, c2d.Tr = np.array([2.0, 2.0, 0.0]), np.array([4.0, 4.0, 0.0]), np.array([0.0, 0.0, 1.0]), np.array([0.0, 0.0, 3.0])
#        Wt, Wr = 2.*3.14/Tt, 2.*3.14/Tr 
#        for ii in range(len(Tt)):
#            if Tt[ii] == 0.0:
#                Wt[ii]=0.0
#            if Tr[ii] == 0.0:
#                Wr[ii]=0.0
#        displ = dx, dy, dz = At[0]*sin(Wt[0]*t), At[1]*sin(Wt[1]*t), At[2]*sin(Wt[2]*t)
#        rot   = rx, ry, rz = Ar[0]*sin(Wr[0]*t), Ar[1]*sin(Wr[1]*t), Ar[2]*sin(Wr[2]*t)
#        init_barycenter, last_position, rot0 = pos, pos*2.0, np.array([0.0, 0.0, 2.0])
#        c2d.init_barycenter, c2d.last_position, c2d.last_rotation_euler = init_barycenter, last_position, rot0
#        displTot = displ - (last_position - init_barycenter)
#        rotTot   = rot  - rot0
#        # from bodydynamics
#        displTot2, rotTot2 = c2d.imposeSinusoidalMotion()
#        npt.assert_equal(displTot, displTot2)
#        npt.assert_equal(rotTot, rotTot2)

    def testStep(self):
        from proteus import Domain  
        from proteus import SpatialTools as st    
        from proteus.mprans import SpatialTools as mst
        from proteus.mprans import BodyDynamics as bd   

        # parameters
        domain2D = Domain.PlanarStraightLineGraphDomain()
        pos = np.array([1.0, 1.0, 0.0])
        dim = (0.3, 0.385)
        caisson = mst.Rectangle(domain=domain2D, dim=dim, coords=[pos[0], pos[1]])
        caisson2D = bd.RigidBody(shape=caisson)
        mst.assembleDomain(domain2D)      
        c2d = caisson2D
        nd = c2d.nd = 2
        Ix = ((dim[0]**2)/12.)
        Iy = ((dim[1]**2)/12.)
        It = Ix + Iy 
        mass = 50.0
        I = mass*It
        c2d.mass, c2d.It = mass, It
        init_barycenter, last_position, last_velocity, last_rotation, last_ang_vel = pos, pos*2.0, np.array([0.05, 0.07, 0.0]), np.array([0.0, 0.0, 0.0001]), np.array([0., 0., 0.0000002])
        c2d.init_barycenter, c2d.last_position, c2d.last_velocity, c2d.last_rotation, c2d.last_ang_vel = init_barycenter, last_position, last_velocity, last_rotation, last_ang_vel 
        F = Fx, Fy, Fz = np.array([200., 300., 0.0])
        M = Mx, My, Mz = np.array([0., 0., 50.0])
        dt, substeps = 0.001, 20  
        dt_sub = c2d.dt_sub = float(dt/substeps)
        c2d.F, c2d.substeps, c2d.dt, c2d.acceleration = F, substeps, dt, c2d.getAcceleration()
        c2d.InputMotion = False
        c2d.scheme = 'Forward_Euler'
        h, tra_velocity = bd.forward_euler(last_position, last_velocity, F/mass, dt)
        rot, rot_velocity = bd.forward_euler(last_rotation, last_ang_vel, M/I, dt)
        if rot[2] < 0.0: rot[2]=-rot[2]
        # 2d calculation
        caisson.translate(h[:nd]), caisson.rotate(rot[2])
        posTra1, rot1 = caisson.barycenter, caisson.coords_system
        # 2d from bodydynamics
        caisson.translate(-h[:nd]), caisson.rotate(-rot[2])
        c2d.step(dt)
        posTra2, rot2 = c2d.position, c2d.rotation[:nd,:nd]
        npt.assert_allclose(posTra1, posTra2, atol=1e-10)
        npt.assert_allclose(rot1, rot2, atol=1e-10)

#    def testSetRecordValues(self):
#        from proteus import Domain  
#        from proteus import SpatialTools as st    
#        from proteus.mprans import SpatialTools as mst
#        from proteus.mprans import BodyDynamics as bd   
        
#        filename='test', all_values=True


class CaissonBody(unittest.TestCase):

    def testGetInertia(self):
        from proteus import Domain  
        from proteus import SpatialTools as st    
        from proteus.mprans import SpatialTools as mst
        from proteus.mprans import BodyDynamics as bd   
        # parameters
        domain2D = Domain.PlanarStraightLineGraphDomain()
        pos = np.array([1.0, 1.0, 0.0])
        dim = (0.3, 0.385)
        dt, substeps = 0.001, 20  
        caisson = mst.Rectangle(domain=domain2D, dim=dim, coords=[pos[0], pos[1]])
        caisson2D = bd.CaissonBody(shape=caisson, substeps=substeps)
        mst.assembleDomain(domain2D)      
        c2d = caisson2D
        c2d.nd = 2
        mass = c2d.mass = 50.0
        pivot = np.array([(caisson.vertices[0][0]+caisson.vertices[1][0])*0.5, caisson.vertices[0][1], 0.0])
        d = dx, dy, dz = pivot - caisson.barycenter
        # inertia transformation with different pivot 
        Ix = ((dim[0]**2)/12.) + dy**2
        Iy = ((dim[1]**2)/12.) + dx**2
        It = Ix + Iy 
        I = It*mass
        I1 = c2d.getInertia(vec=(0.,0.,1.), pivot=pivot)
        npt.assert_equal(I, I1)


    def testFrictionModule(self):
        from proteus import Domain  
        from proteus import SpatialTools as st    
        from proteus.mprans import SpatialTools as mst
        from proteus.mprans import BodyDynamics as bd   

        #########################
        # 1st case              #
        #########################
        # When sliding is true and caisson already experiences plastic displacements

        # parameters
        domain2D = Domain.PlanarStraightLineGraphDomain()
        pos = np.array([1.0, 1.0, 0.0])
        dim = (0.3, 0.385)
        dt, substeps = 0.001, 20  
        caisson = mst.Rectangle(domain=domain2D, dim=dim, coords=[pos[0], pos[1]])
        caisson2D = bd.CaissonBody(shape=caisson, substeps=substeps)
        mst.assembleDomain(domain2D)      
        c2d = caisson2D
        c2d.nd = 2
        dt_sub = c2d.dt_sub = float(dt/substeps)
        Ix = ((dim[0]**2)/12.)
        Iy = ((dim[1]**2)/12.)
        It = Ix + Iy 
        mass = 50.0
        c2d.mass, c2d.It = mass, It
        K = Kx, Ky, Kz = [200000., 250000., 0.0]
        C = Cx, Cy, Cz = [2000., 4000., 0.0]
        c2d.Kx, c2d.Ky, c2d.Kz, c2d.Cx, c2d.Cy, c2d.Cz =  Kx, Ky, Kz, Cx, Cy, Cz
        c2d.substeps, c2d.dt, c2d.acceleration = substeps, dt, c2d.getAcceleration()
        m_static, m_dynamic = 0.6, 0.4
        c2d.setFriction(friction=True, m_static=m_static, m_dynamic=m_dynamic, tolerance=0.000001, grainSize=0.02)
        # sliding == True. dynamic sign and dynamic friction
        F = Fx, Fy, Fz = np.array([200., -300., 0.0])
        disp = np.array([0.001, 0.001, 0.0])
        init_barycenter, last_position, last_velocity = pos, pos+disp, np.array([0.05, 0.07, 0.0])
        c2d.last_uxEl = pos[0]+disp[0]-init_barycenter[0]           
        c2d.F, c2d.init_barycenter, c2d.last_position, c2d.last_velocity = F, init_barycenter, last_position, last_velocity  
        eps = 10.**-30
        sign_static, sign_dynamic = Fx/abs(Fx+eps), last_velocity[0]/abs(last_velocity[0]+eps)
        c2d.sliding, sign, m = True, sign_dynamic, m_dynamic
        # vertical calculation and frictional force        
        uy0, vy0 = (last_position[1] - init_barycenter[1]), last_velocity[1]
        ay0 = (Fy - Cy*vy0 - Ky*uy0) / mass
        uy, vy, ay = bd.runge_kutta(uy0, vy0, ay0, dt_sub, substeps, Fy, Ky, Cy, mass, False)
        Rx = -Kx*(last_position[0]-init_barycenter[0])
        Ry = -Ky*uy
        Ftan = -sign*abs(Ry)*m
        if Ftan == 0.0: 
            Ftan = -sign*abs(Fy)*m
        # runge-kutta calculation
        c2d.scheme = 'Runge_Kutta'
        ux0, uz0 = last_position[0] - init_barycenter[0], last_position[2] - init_barycenter[2]
        vx0, vz0 = last_velocity[0], last_velocity[2]
        Fh = Fx+Ftan
        Kx, Cx = 0.0, 0.0
        ax0, az0 = (Fh - Cx*vx0 - Kx*ux0) / mass , (Fz - Cz*vz0 - Kz*uz0) / mass
        ux, vx, ax = bd.runge_kutta(ux0, vx0, ax0, dt_sub, substeps, Fh, Kx, Cx, mass, True)
        uz, vz, az = bd.runge_kutta(uz0, vz0, az0, dt_sub, substeps, Fz, Kz, Cz, mass, False)
        # bodydynamics calculation
        c2d.friction_module(dt) 
        # tests
        npt.assert_equal(c2d.ux, ux)
        EL1, PL1 = 0.0, 1.0 # elastic and plastic motion parameters
        npt.assert_equal(c2d.EL, EL1)
        npt.assert_equal(c2d.PL, PL1)
         

        #########################
        # 2nd case              #
        #########################
        # When sliding is false but the caisson starts to experience sliding motion

        # parameters
        domain2D = Domain.PlanarStraightLineGraphDomain()
        pos = np.array([1.0, 1.0, 0.0])
        dim = (0.3, 0.385)
        dt, substeps = 0.001, 20  
        caisson = mst.Rectangle(domain=domain2D, dim=dim, coords=[pos[0], pos[1]])
        caisson2D = bd.CaissonBody(shape=caisson, substeps=substeps)
        mst.assembleDomain(domain2D)      
        c2d = caisson2D
        c2d.nd = 2
        dt_sub = c2d.dt_sub = float(dt/substeps)
        Ix = ((dim[0]**2)/12.)
        Iy = ((dim[1]**2)/12.)
        It = Ix + Iy 
        mass = 50.0
        c2d.mass, c2d.It = mass, It
        K = Kx, Ky, Kz = [200000., 250000., 0.0]
        C = Cx, Cy, Cz = [2000., 4000., 0.0]
        c2d.Kx, c2d.Ky, c2d.Kz, c2d.Cx, c2d.Cy, c2d.Cz =  Kx, Ky, Kz, Cx, Cy, Cz
        c2d.substeps, c2d.dt, c2d.acceleration = substeps, dt, c2d.getAcceleration()
        m_static, m_dynamic = 0.6, 0.4
        c2d.setFriction(friction=True, m_static=m_static, m_dynamic=m_dynamic, tolerance=0.000001, grainSize=0.02)
        # sliding == False. static sign and static friction
        F = Fx, Fy, Fz = np.array([200., -300., 0.0])
        disp = np.array([0.001, 0.001, 0.0])
        init_barycenter, last_position, last_velocity = pos, pos+disp, np.array([0.05, 0.07, 0.0])
        c2d.last_uxEl = pos[0]+disp[0]-init_barycenter[0]           
        c2d.F, c2d.init_barycenter, c2d.last_position, c2d.last_velocity = F, init_barycenter, last_position, last_velocity  
        eps = 10.**-30
        sign_static, sign_dynamic = Fx/abs(Fx+eps), last_velocity[0]/abs(last_velocity[0]+eps)
        c2d.sliding, sign, m = False, sign_static, m_static
        # vertical calculation and frictional force        
        uy0, vy0 = (last_position[1] - init_barycenter[1]), last_velocity[1]
        ay0 = (Fy - Cy*vy0 - Ky*uy0) / mass
        uy, vy, ay = bd.runge_kutta(uy0, vy0, ay0, dt_sub, substeps, Fy, Ky, Cy, mass, False)
        Rx = -Kx*(last_position[0]-init_barycenter[0])
        Ry = -Ky*uy
        Ftan = -sign*abs(Ry)*m
        if Ftan == 0.0: 
            Ftan = -sign*abs(Fy)*m
        # runge-kutta calculation
        c2d.scheme = 'Runge_Kutta'
        ux0, uz0 = last_position[0] - init_barycenter[0], last_position[2] - init_barycenter[2]
        vx0, vz0 = last_velocity[0], last_velocity[2]
        Fh = Fx+Ftan
        Kx, Cx = 0.0, 0.0
        ax0, az0 = (Fh - Cx*vx0 - Kx*ux0) / mass , (Fz - Cz*vz0 - Kz*uz0) / mass
        ux, vx, ax = bd.runge_kutta(ux0, vx0, ax0, dt_sub, substeps, Fh, Kx, Cx, mass, True)
        uz, vz, az = bd.runge_kutta(uz0, vz0, az0, dt_sub, substeps, Fz, Kz, Cz, mass, False)
        # bodydynamics calculation
        c2d.friction_module(dt)
        # tests
        npt.assert_equal(c2d.ux, ux)
        EL1, PL1 = 0.0, 1.0 # elastic and plastic motion parameters
        npt.assert_equal(c2d.EL, EL1)
        npt.assert_equal(c2d.PL, PL1)
         

        #########################
        # 3rd case              #
        #########################
        # When caisson experiences vibration motion

        # parameters
        domain2D = Domain.PlanarStraightLineGraphDomain()
        pos = np.array([1.0, 1.0, 0.0])
        dim = (0.3, 0.385)
        dt, substeps = 0.001, 20  
        caisson = mst.Rectangle(domain=domain2D, dim=dim, coords=[pos[0], pos[1]])
        caisson2D = bd.CaissonBody(shape=caisson, substeps=substeps)
        mst.assembleDomain(domain2D)      
        c2d = caisson2D
        c2d.nd = 2
        dt_sub = c2d.dt_sub = float(dt/substeps)
        Ix = ((dim[0]**2)/12.)
        Iy = ((dim[1]**2)/12.)
        It = Ix + Iy 
        mass = 50.0
        c2d.mass, c2d.It = mass, It
        K = Kx, Ky, Kz = [200000., 250000., 0.0]
        C = Cx, Cy, Cz = [2000., 4000., 0.0]
        c2d.Kx, c2d.Ky, c2d.Kz, c2d.Cx, c2d.Cy, c2d.Cz =  Kx, Ky, Kz, Cx, Cy, Cz
        c2d.substeps, c2d.dt, c2d.acceleration = substeps, dt, c2d.getAcceleration()
        m_static, m_dynamic = 0.6, 0.4
        c2d.setFriction(friction=True, m_static=m_static, m_dynamic=m_dynamic, tolerance=0.000001, grainSize=0.02)
        # sliding == False. static sign and static friction. Kx and Cx different from 0!
        F = Fx, Fy, Fz = np.array([200., -300., 0.0])
        disp = np.array([0.00001, 0.00001, 0.0])
        init_barycenter, last_position, last_velocity = pos, pos+disp, np.array([0.05, 0.07, 0.0])
        c2d.last_uxEl = pos[0]+disp[0]-init_barycenter[0]           
        c2d.F, c2d.init_barycenter, c2d.last_position, c2d.last_velocity = F, init_barycenter, last_position, last_velocity  
        eps = 10.**-30
        sign_static, sign_dynamic = Fx/abs(Fx+eps), last_velocity[0]/abs(last_velocity[0]+eps)
        c2d.sliding, sign, m = False, sign_static, m_static
        # vertical calculation and frictional force        
        uy0, vy0 = (last_position[1] - init_barycenter[1]), last_velocity[1]
        ay0 = (Fy - Cy*vy0 - Ky*uy0) / mass
        uy, vy, ay = bd.runge_kutta(uy0, vy0, ay0, dt_sub, substeps, Fy, Ky, Cy, mass, False)
        Rx = -Kx*(last_position[0]-init_barycenter[0])
        Ry = -Ky*uy
        Ftan = -sign*abs(Ry)*m
        if Ftan == 0.0: 
            Ftan = -sign*abs(Fy)*m
        # runge-kutta calculation
        c2d.scheme = 'Runge_Kutta'
        ux0, uz0 = last_position[0] - init_barycenter[0], last_position[2] - init_barycenter[2]
        vx0, vz0 = last_velocity[0], last_velocity[2]
        Fh = Fx        
        ax0, az0 = (Fh - Cx*vx0 - Kx*ux0) / mass , (Fz - Cz*vz0 - Kz*uz0) / mass
        ux, vx, ax = bd.runge_kutta(ux0, vx0, ax0, dt_sub, substeps, Fh, Kx, Cx, mass, True)
        uz, vz, az = bd.runge_kutta(uz0, vz0, az0, dt_sub, substeps, Fz, Kz, Cz, mass, False)
        # bodydynamics calculation
        c2d.friction_module(dt)
        # tests
        npt.assert_equal(c2d.ux, ux)
        EL1, PL1 = 1.0, 0.0 # elastic and plastic motion parameters
        npt.assert_equal(c2d.EL, EL1)
        npt.assert_equal(c2d.PL, PL1)
      
    def testOverturningModule(self):
        from proteus import Domain  
        from proteus import SpatialTools as st    
        from proteus.mprans import SpatialTools as mst
        from proteus.mprans import BodyDynamics as bd   

        # parameters
        domain2D = Domain.PlanarStraightLineGraphDomain()
        pos = np.array([1.0, 1.0, 0.0])
        dim = (0.3, 0.385)
        dt, substeps = 0.001, 20  
        caisson = mst.Rectangle(domain=domain2D, dim=dim, coords=[pos[0], pos[1]])
        caisson2D = bd.CaissonBody(shape=caisson, substeps=substeps)
        mst.assembleDomain(domain2D)      
        c2d = caisson2D
        c2d.nd = 2
        dt_sub = c2d.dt_sub = float(dt/substeps)
        Krot = 200000.
        Crot = 2000.
        c2d.Krot, c2d.Crot =  Krot, Crot
        c2d.substeps, c2d.dt, c2d.ang_acc = substeps, dt, c2d.getAngularAcceleration()
        c2d.setFriction(friction=True, m_static=0.0, m_dynamic=0.0, tolerance=0.000001, grainSize=0.02)
        rotation, last_rotation = c2d.Shape.coords_system, c2d.Shape.coords_system
        ang_vel, last_ang_vel = np.array([0., 0., 0.]), np.array([0., 0., 0.])
        F = Fx, Fy, Fz = np.array([200., -300., 0.0])
        M =  np.array([0., 0., 10.0])
        init_barycenter, last_position = pos, pos     
        c2d.F, c2d.M, c2d.init_barycenter, c2d.last_position = F, M, init_barycenter, last_position  
        eps = 10.**-30
        mass = c2d.mass = 50.0
        # moment and inertia transformations
        pivot = c2d.pivot_friction = np.array([(caisson.vertices[0][0]+caisson.vertices[1][0])*0.5, caisson.vertices[0][1], 0.0])
        barycenter = caisson.barycenter
        distance = dx, dy, dz = pivot - barycenter
        Mpivot = np.array([0., 0., dx*Fy-dy*Fx]) # Moment calculated in pivot
        Mp = M - Mpivot # Moment transformation through the new pivot
        Ix = ((dim[0]**2)/12.) + dy**2
        Iy = ((dim[1]**2)/12.) + dx**2
        It = Ix + Iy    
        I = It*mass
        c2d.springs = True
        # runge-kutta calculation
        c2d.scheme = 'Runge_Kutta'
        rz0, vrz0 = atan2(last_rotation[0,1], last_rotation[0,0]), last_ang_vel[2]
        arz0 = (Mp[2] - Crot*vrz0 - Krot*rz0) / I
        rz, vrz, arz = bd.runge_kutta(rz0, vrz0, arz0, dt_sub, substeps, Mp[2], Krot, Crot, I, False)
        # calculation with bodydynamics mopdule
        c2d.cV_init = c2d.cV = c2d.cV_last = caisson.vertices
        c2d.inertia = c2d.getInertia(pivot=pivot)
        c2d.overturning_module(dt)    
        # test
        ang_disp = rz - rz0
        npt.assert_equal(c2d.ang_disp[2], ang_disp)
        npt.assert_equal(c2d.ang_vel[2], vrz)
        npt.assert_equal(c2d.ang_acc[2], arz)





if __name__ == '__main__':
    unittest.main(verbosity=2)


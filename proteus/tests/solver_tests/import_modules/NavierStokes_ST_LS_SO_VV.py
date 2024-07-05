from proteus.TransportCoefficients import TwophaseNavierStokes_ST_LS_SO as TPNSE_ST_LS_SO
from proteus import ctransportCoefficients

class NavierStokes_ST_LS_SO_VV(TPNSE_ST_LS_SO):
    '''This coefficients class is specially designed to allow a time
       varying viscosity parameter.  The purpose of this is to allow us
       to use a pseudo-continuity approach to solve high Reynolds number NSE.'''
    def evaluate(self,t,c):
        self.nu_0 = 10**(1-t) ;  self.nu_1 = 10**(1-t)
        if self.nu_0 == self.nu_1:
            self.nu = self.nu_0
        import math
        #self.rho_0 = 1000.0*self.rho_1
        #self.nu_0 = 0.1*self.nu_1
        if c[('u',0)].shape == self.q_phi.shape:
            phi = self.q_phi
            n   = self.q_n
            kappa = self.q_kappa
            #slopeAngle=0.1*math.pi/2.0#math.pi/4.0
            #surfaceNormal = [-sin(slopeAngle),cos(slopeAngle)]
            #waterLevel=0.5
            #for eN in range(phi.shape[0]):
            #   for k in range(phi.shape[1]):
            #       phi[eN,k] = (c['x'][eN,k,0] - 0.5)*surfaceNormal[0]+(c['x'][eN,k,1] - waterLevel)*surfaceNormal[1]
        elif c[('u',0)].shape == self.ebqe_phi.shape:
            phi   = self.ebqe_phi
            n     = self.ebqe_n
            kappa = self.ebqe_kappa
        else:
            phi   = self.ebq_phi
            n     = self.ebq_n
            kappa = self.ebq_kappa
        #mwf debug
        #waterLevelBase = 0.529
        #for i in range(len(phi.flat)):
            #if abs(phi.flat[i]) > 0.0:
            #    assert abs(phi.flat[i] - (c['x'].flat[3*i+1] - waterLevelBase)) <= 1.0e-5, "Problem with phi t=%s phi.shape=%s i=%s phi=%s y=%s wl=%s " % (t,phi.shape,i,phi.flat[i],c['x'].flat[3*i+1],waterLevelBase)
            #phi.flat[i] = c['x'].flat[3*i+1] - waterLevelBase#self.waterLevel
        #self.sd=False
        if self.nd==2:
            if self.sd:
                ctransportCoefficients.TwophaseNavierStokes_ST_LS_SO_2D_Evaluate_sd(self.eps_density,
                                                                                    self.eps_viscosity,
                                                                                    self.sigma,
                                                                                    self.rho_0,
                                                                                    self.nu_0,
                                                                                    self.rho_1,
                                                                                    self.nu_1,
                                                                                    self.g,
                                                                                    phi,
                                                                                    n,
                                                                                    kappa,
                                                                                    c[('u',0)],
                                                                                    c[('grad(u)',0)],
                                                                                    c[('u',1)],
                                                                                    c[('u',2)],
                                                                                    c[('m',1)],
                                                                                    c[('dm',1,1)],
                                                                                    c[('m',2)],
                                                                                    c[('dm',2,2)],
                                                                                    c[('f',0)],
                                                                                    c[('df',0,1)],
                                                                                    c[('df',0,2)],
                                                                                    c[('f',1)],
                                                                                    c[('df',1,1)],
                                                                                    c[('df',1,2)],
                                                                                    c[('f',2)],
                                                                                    c[('df',2,1)],
                                                                                    c[('df',2,2)],
                                                                                    c[('a',1,1)],
                                                                                    c[('a',2,2)],
                                                                                    c[('a',1,2)],
                                                                                    c[('a',2,1)],
                                                                                    c[('r',1)],
                                                                                    c[('r',2)],
                                                                                    c[('H',1)],
                                                                                    c[('dH',1,0)],
                                                                                    c[('H',2)],
                                                                                    c[('dH',2,0)])
            else:
                ctransportCoefficeints.TwophaseNavierStokes_ST_LS_SO_2D_Evaluate(self.eps_density,
                                                                                 self.eps_viscosity,
                                                                                 self.sigma,
                                                                                 self.rho_0,
                                                                                 self.nu_0,
                                                                                 self.rho_1,
                                                                                 self.nu_1,
                                                                                 self.g,
                                                                                 phi,
                                                                                 n,
                                                                                 kappa,
                                                                                 c[('u',0)],
                                                                                 c[('grad(u)',0)],
                                                                                 c[('u',1)],
                                                                                 c[('u',2)],
                                                                                 c[('m',1)],
                                                                                 c[('dm',1,1)],
                                                                                 c[('m',2)],
                                                                                 c[('dm',2,2)],
                                                                                 c[('f',0)],
                                                                                 c[('df',0,1)],
                                                                                 c[('df',0,2)],
                                                                                 c[('f',1)],
                                                                                 c[('df',1,1)],
                                                                                 c[('df',1,2)],
                                                                                 c[('f',2)],
                                                                                 c[('df',2,1)],
                                                                                 c[('df',2,2)],
                                                                                 c[('a',1,1)],
                                                                                 c[('a',2,2)],
                                                                                 c[('a',1,2)],
                                                                                 c[('a',2,1)],
                                                                                 c[('r',1)],
                                                                                 c[('r',2)],
                                                                                 c[('H',1)],
                                                                                 c[('dH',1,0)],
                                                                                 c[('H',2)],
                                                                                 c[('dH',2,0)])
            if self.stokes:
                c[('f',1)].flat[:] = 0.0
                c[('df',1,1)].flat[:] = 0.0
                c[('df',1,2)].flat[:] = 0.0
                c[('f',2)].flat[:] = 0.0
                c[('df',2,1)].flat[:] = 0.0
                c[('df',2,2)].flat[:] = 0.0
        elif  self.nd==3:
            if self.sd:
                ctransportCoefficients.TwophaseNavierStokes_ST_LS_SO_3D_Evaluate_sd(self.eps_density,
                                                                                    self.eps_viscosity,
                                                                                    self.sigma,
                                                                                    self.rho_0,
                                                                                    self.nu_0,
                                                                                    self.rho_1,
                                                                                    self.nu_1,
                                                                                    self.g,
                                                                                    phi,
                                                                                    n,
                                                                                    kappa,
                                                                                    c[('u',0)],
                                                                                    c[('grad(u)',0)],
                                                                                    c[('u',1)],
                                                                                    c[('u',2)],
                                                                                    c[('u',3)],
                                                                                    c[('m',1)],
                                                                                    c[('dm',1,1)],
                                                                                    c[('m',2)],
                                                                                    c[('dm',2,2)],
                                                                                    c[('m',3)],
                                                                                    c[('dm',3,3)],
                                                                                    c[('f',0)],
                                                                                    c[('df',0,1)],
                                                                                    c[('df',0,2)],
                                                                                    c[('df',0,3)],
                                                                                    c[('f',1)],
                                                                                    c[('df',1,1)],
                                                                                    c[('df',1,2)],
                                                                                    c[('df',1,3)],
                                                                                    c[('f',2)],
                                                                                    c[('df',2,1)],
                                                                                    c[('df',2,2)],
                                                                                    c[('df',2,3)],
                                                                                    c[('f',3)],
                                                                                    c[('df',3,1)],
                                                                                    c[('df',3,2)],
                                                                                    c[('df',3,3)],
                                                                                    c[('a',1,1)],
                                                                                    c[('a',2,2)],
                                                                                    c[('a',3,3)],
                                                                                    c[('a',1,2)],
                                                                                    c[('a',1,3)],
                                                                                    c[('a',2,1)],
                                                                                    c[('a',2,3)],
                                                                                    c[('a',3,1)],
                                                                                    c[('a',3,2)],
                                                                                    c[('r',1)],
                                                                                    c[('r',2)],
                                                                                    c[('r',3)],
                                                                                    c[('H',1)],
                                                                                    c[('dH',1,0)],
                                                                                    c[('H',2)],
                                                                                    c[('dH',2,0)],
                                                                                    c[('H',3)],
                                                                                    c[('dH',3,0)])
            else:
                ctransportCoefficients.TwophaseNavierStokes_ST_LS_SO_3D_Evaluate(self.eps_density,
                                                                                 self.eps_viscosity,
                                                                                 self.sigma,
                                                                                 self.rho_0,
                                                                                 self.nu_0,
                                                                                 self.rho_1,
                                                                                 self.nu_1,
                                                                                 self.g,
                                                                                 phi,
                                                                                 n,
                                                                                 kappa,
                                                                                 c[('u',0)],
                                                                                 c[('grad(u)',0)],
                                                                                 c[('u',1)],
                                                                                 c[('u',2)],
                                                                                 c[('u',3)],
                                                                                 c[('m',1)],
                                                                                 c[('dm',1,1)],
                                                                                 c[('m',2)],
                                                                                 c[('dm',2,2)],
                                                                                 c[('m',3)],
                                                                                 c[('dm',3,3)],
                                                                                 c[('f',0)],
                                                                                 c[('df',0,1)],
                                                                                 c[('df',0,2)],
                                                                                 c[('df',0,3)],
                                                                                 c[('f',1)],
                                                                                 c[('df',1,1)],
                                                                                 c[('df',1,2)],
                                                                                 c[('df',1,3)],
                                                                                 c[('f',2)],
                                                                                 c[('df',2,1)],
                                                                                 c[('df',2,2)],
                                                                                 c[('df',2,3)],
                                                                                 c[('f',3)],
                                                                                 c[('df',3,1)],
                                                                                 c[('df',3,2)],
                                                                                 c[('df',3,3)],
                                                                                 c[('a',1,1)],
                                                                                 c[('a',2,2)],
                                                                                 c[('a',3,3)],
                                                                                 c[('a',1,2)],
                                                                                 c[('a',1,3)],
                                                                                 c[('a',2,1)],
                                                                                 c[('a',2,3)],
                                                                                 c[('a',3,1)],
                                                                                 c[('a',3,2)],
                                                                                 c[('r',1)],
                                                                                 c[('r',2)],
                                                                                 c[('r',3)],
                                                                                 c[('H',1)],
                                                                                 c[('dH',1,0)],
                                                                                 c[('H',2)],
                                                                                 c[('dH',2,0)],
                                                                                 c[('H',3)],
                                                                                 c[('dH',3,0)])
            if self.stokes:
                c[('f',1)].flat[:] = 0.0
                c[('df',1,1)].flat[:] = 0.0
                c[('df',1,2)].flat[:] = 0.0
                c[('df',1,3)].flat[:] = 0.0
                c[('f',2)].flat[:] = 0.0
                c[('df',2,1)].flat[:] = 0.0
                c[('df',2,2)].flat[:] = 0.0
                c[('df',2,3)].flat[:] = 0.0
                c[('f',3)].flat[:] = 0.0
                c[('df',3,1)].flat[:] = 0.0
                c[('df',3,2)].flat[:] = 0.0
                c[('df',3,3)].flat[:] = 0.0
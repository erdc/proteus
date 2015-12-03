from proteus import Comm, Profiling
import numpy as np
import numpy.testing as npt
import unittest

comm = Comm.init()
Profiling.procID = comm.rank()

Profiling.logEvent("Testing SedClosure")
class GlobalVariables():
    def __init__(self):
        from proteus.mprans.SedClosure import HsuSedStress
        self.C4e = 1.
        self.C3e = 1.2
        self.eR = 0.8
        self.aDarcy = 1.
        self.bForch = 1.
        self.grain = 0.1
        self.packFraction = 0.2
        self.packMargin = 0.01
        self.sigmaC = 1.1
        self.maxFraction = 0.635
        self.frFraction = 0.57
        self.sedSt = HsuSedStress( self.aDarcy, self.bForch, self.grain, self.packFraction,  self.packMargin, self.sigmaC, self.C3e, self.C4e, self.eR)
    

class TestHsu(unittest.TestCase):
    def testGradularDrag1(self):
        gl = GlobalVariables()
        uf = np.array([5.,4.],"d")
        us = np.array([1.,1.],"d") 
        umag_temp = (uf - us)*(uf - us)
        umag = np.sqrt(sum(umag_temp))
        # Testing for sedF >> gl.packFraction
        nu = 1.
        sedF = 0.5
        drag = gl.aDarcy * nu* sedF /((1.-sedF)*gl.grain**2) +  gl.bForch * umag / gl.grain
        # For sedF > pacFraction - > drag = a * nu* sedF /((1-sedF)*gl.grain^2) + beta * umag / gl.grain
        drag2 = gl.sedSt.betaCoeff(sedF, uf, us, nu)
        self.assertTrue(drag == drag2)
    def testGradularDrag2(self):
        gl=GlobalVariables()
        uf = np.array([5.,4.],"d")
        us = np.array([1.,1.],"d") 
        umag_temp = (uf - us)*(uf - us)
        umag = np.sqrt(sum(umag_temp))
        # Testing for sedF << gl.packFraction and Rep = (1. - sed) * umag*nu/gl.grain = 0.9 * 5. * 0.1 / 1. = 0.45
        nu = 1.
        sedF = 0.1
        Rep = (1.- sedF)*umag*gl.grain / nu    
        drag =  ( 24. * (1.+0.15*Rep**(0.687))/Rep) * 0.75 * umag * (1. -sedF)**(-1.65) / gl.grain # Chen and Hsu 2014
        drag2 = gl.sedSt.betaCoeff(sedF, uf, us, nu)
        self.assertTrue(round(drag,10) == round(drag2,10))

    def testGradularDrag3(self):
        gl=GlobalVariables()
        # Constant params
        uf = np.array([5.,4.],"d")
        us = np.array([1.,1.],"d") 
        umag_temp = (uf - us)*(uf - us)
        umag = np.sqrt(sum(umag_temp))
        # Testing for sedF << gl.packFraction and Rep > 1000
        sedF = 0.1
        nu = 1e-4
        Rep = (1.- sedF) * umag * gl.grain / nu    
        drag =  ( 0.44 * 0.75 * umag * (1. -sedF)**(-1.65) )/ gl.grain # Chen and Hsu 2014
        drag2 = gl.sedSt.betaCoeff(sedF, uf, us, nu)
        self.assertTrue(round(drag,10) == round(drag2,10))
    def testGradularDrag4(self):
        gl=GlobalVariables()
        # Constant params
        uf = np.array([5.,4.],"d")
        us = np.array([1.,1.],"d") 
        umag_temp = (uf - us)*(uf - us)
        umag = np.sqrt(sum(umag_temp))
        # Testing for sedF =  gl.packFraction +0.5 packmargin and Rep > 1000
        sedF = 0.205
        nu = 1e-4
        Rep = (1.- sedF) * umag * gl.grain / nu
        draga = gl.aDarcy * nu* sedF /((1.-sedF)*gl.grain**2) +  gl.bForch * umag / gl.grain
        dragb =  ( 0.44 * 0.75 * umag * (1. -sedF)**(-1.65) )/ gl.grain # Chen and Hsu 2014
        w = 0.5 + (sedF - gl.packFraction) / (2. * gl.packMargin)
        drag = w*draga + (1.-w) * dragb
        drag2 = gl.sedSt.betaCoeff(sedF, uf, us, nu)
        self.assertTrue(round(drag,10) == round(drag2,10))
    def testgs0(self):
        gl=GlobalVariables()
        f = 8
        sedF = 0.2
        gs0 = gl.sedSt.gs0(sedF)
        self.assertTrue(gs0 == 0.5*(2-sedF)/(1-sedF)**3)
        sedF = 0.55
        gs0 = gl.sedSt.gs0(sedF)
        self.assertTrue(round(gs0,f) ==round( 0.5*(2-0.49)/(1-0.49)**3 * (0.64-0.49)/(0.64-sedF),f))
        sedF = 0.65
        gs0 = gl.sedSt.gs0(sedF)
        self.assertTrue(round(gs0,f) == round(0.5*(2-0.49)/(1-0.49)**3 * (0.64-0.49)/(0.64-0.635),f))
    def testTkeSed(self):
        gl=GlobalVariables()
        import random
        f = 10
        uf = np.array([5.,4.],"d")
        us = np.array([1.,1.],"d") 
        gradC=np.array([0.1,0.1])
        rhoS = 2000
        rhoF = 1000
        nu = 1e-4
        nuT = 1e-2
        sedF = 0.3
# Setting 0 t_c
        theta_n = random.random() + 1e-30
        kappa_n = random.random() + 1e-30
        epsilon_n = random.random() + 1e-30


        beta = gl.sedSt.betaCoeff(sedF, uf, us, nu)      
        t_p =rhoS/ beta
        l_c = np.sqrt(np.pi)*gl.grain / (24.*sedF * gl.sedSt.gs0(sedF))
        t_cl = min(l_c/np.sqrt(theta_n),0.165*kappa_n/epsilon_n)
        aa = 1/( 1. + t_p/t_cl)

       
        es1 = 2.*beta * rhoS*(1-aa)*sedF*kappa_n/((1-sedF)*rhoF)

        UgradC = np.dot((uf - us),gradC)

        es2  = beta * rhoF * nuT * UgradC / ((1-sedF)*rhoF)

        kappa_sed = gl.sedSt.kappa_sed(sedF,rhoF,rhoS,uf,us,gradC,nu,theta_n,kappa_n,epsilon_n,nuT)
        self.assertTrue(round(kappa_sed,f) ==round( -es1+es2,f))
        
    def testEpsSed(self):
        gl=GlobalVariables()
        import random
        f = 8
        uf = np.array([5.,4.],"d")
        us = np.array([1.,1.],"d") 
        gradC=np.array([0.1,0.1])
        rhoS = 2000
        rhoF = 1000
        nu = 1e-4
        nuT = 1e-2
        sedF = 0.3
# Setting 0 t_c
        theta_n = random.random() + 1e-30
        kappa_n = random.random() + 1e-30
        epsilon_n = random.random() + 1e-30


        beta = gl.sedSt.betaCoeff(sedF, uf, us, nu)      
        t_p =rhoS/ beta
        l_c = np.sqrt(np.pi)*gl.grain / (24.*sedF * gl.sedSt.gs0(sedF))
        t_cl = min(l_c/np.sqrt(theta_n),0.165*kappa_n/epsilon_n)
        aa = 1/( 1. + t_p/t_cl)

       
        es1 = 2.*beta * rhoS*(1-aa)*sedF*kappa_n/((1-sedF)*rhoF)

        UgradC = np.dot((uf - us),gradC)

        es2  = beta * rhoF * nuT * UgradC / ((1-sedF)*rhoF)

        eps_sed = gl.sedSt.eps_sed(sedF,rhoF,rhoS,uf,us,gradC,nu,theta_n,kappa_n,epsilon_n,nuT)

        self.assertTrue(round(eps_sed,f) ==round( -gl.C3e*es1*epsilon_n/kappa_n+gl.C4e*es2*epsilon_n/kappa_n,f))

    def testPsc(self):
        gl=GlobalVariables()
        import random
        f = 10        
        # Setting 0 t_c
        sedF = 0.3
        theta = random.random() + 1e-30
        rhoS = 2000
        psc = gl.sedSt.psc(sedF,rhoS,theta)
        self.assertTrue(round(psc,f) ==round(rhoS*sedF*(1. + 2*(1.+gl.eR)*sedF*gl.sedSt.gs0(sedF))*theta,f))

    def testPscTerm(self):
        gl=GlobalVariables()
        import random
        f = 10        
        # Setting 0 t_c
        sedF = 0.3
        theta = random.random() + 1e-30
        dudx = random.random() + 1e-30
        dvdy = random.random() + 1e-30
        dwdz = random.random() + 1e-30
        divU = dudx + dvdy + dwdz
        rhoS = 2000
        test = gl.sedSt.psc_term(sedF,rhoS,theta,dudx,dvdy,dwdz)
        self.assertTrue(round(test,f) == round(-2.*gl.sedSt.psc(sedF,rhoS,theta)*divU/(3.*rhoS*sedF),f))


    def testdpsc_term_dtheta(self):
        gl=GlobalVariables()
        import random
        f = 10        
        # Setting 0 t_c
        sedF = 0.3
        theta = random.random() + 1e-30
        dudx = random.random() + 1e-30
        dvdy = random.random() + 1e-30
        dwdz = random.random() + 1e-30
        divU = dudx + dvdy + dwdz
        rhoS = 2000
        test = gl.sedSt.dpsc_term_dtheta(sedF,rhoS,dudx,dvdy,dwdz)
        self.assertTrue(round(test,f) == round(-2.*gl.sedSt.psc(sedF,rhoS,theta)*divU/(3.*rhoS*sedF)/theta,f))
    def testMu_sc(self):
        gl=GlobalVariables()
        import random
        sqrt = np.sqrt
        f = 10        
        # Setting 0 t_c
        sedF = 0.3
        theta = random.random() + 1e-30
        rhoS = 2000
        test = gl.sedSt.mu_sc(sedF,rhoS,theta)
        g0 = gl.sedSt.gs0(sedF)
        valid = rhoS * gl.grain * sqrt(theta) * ( 0.8 *sedF**2 * g0 * (1. + gl.eR) / sqrt(np.pi) + (1./15) *sedF**2 * g0 * (1. + gl.eR) * sqrt(np.pi) +  (1./6.) *sedF * sqrt(np.pi) +  (5./48.) * sqrt(np.pi)/((1+gl.eR)*g0))        
        self.assertTrue(round(test,f) == round(valid,f))

    def testL_sc(self):
        gl=GlobalVariables()
        import random
        sqrt = np.sqrt
        f = 10        
        # Setting 0 t_c
        sedF = 0.3
        theta = random.random() + 1e-30
        rhoS = 2000
        test = gl.sedSt.l_sc(sedF,rhoS,theta)
        g0 = gl.sedSt.gs0(sedF)
        valid = (4./3.)*sedF**2 * rhoS * gl.grain *g0*(1+gl.eR)* (sqrt(theta) / sqrt(np.pi) )
        
        self.assertTrue(round(test,f) == round(valid,f))

    def test_tsc_term(self):
        gl=GlobalVariables()
        import random
        sqrt = np.sqrt
        f = 10        
        # Setting 0 t_c
        sedF = 0.3
        theta = random.random() + 1e-30
        dudx = random.random() + 1e-30
        dudy = random.random() + 1e-30
        dudz = random.random() + 1e-30
        dvdx = random.random() + 1e-30
        dvdy = random.random() + 1e-30
        dvdz = random.random() + 1e-30
        dwdx = random.random() + 1e-30
        dwdy = random.random() + 1e-30
        dwdz = random.random() + 1e-30
        rhoS = 2000
        test = gl.sedSt.tausc_term_theta(sedF,rhoS,theta,dudx,dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz)
        divU = dudx + dvdy + dwdz
        mu = gl.sedSt.mu_sc(sedF, rhoS, theta)
        l = gl.sedSt.l_sc(sedF, rhoS, theta)
        s_tensor =  np.array([ [ 2.*dudx ,  dudy+dvdx,  dudz+dwdx],
                                [ dudy+dvdx, 2.*dvdy,    dvdz+dwdy],
                                [ dudz+dwdx, dvdz+dwdy,  2.* dwdz]])
        t_tensor = mu * s_tensor  + (l - (2./3.) * mu) * divU * np.array([ [ 1 , 0 , 0],
                                                   [ 0,  1,  0],
                                                   [ 0,  0,  1] ])
        product = s_tensor * t_tensor
        valid = 0.
        for i in product:
            for j in i:
                valid+=j / (3.*rhoS*sedF)
                
        self.assertTrue(round(test,f) == round(valid,f))


    def testgamma_s(self):
        gl=GlobalVariables()
        import random
        sqrt = np.sqrt
        f = 10        
        # Setting 0 t_c
        sedF = 0.3
        theta = random.random() + 1e-30
        theta_np1 = random.random() + 1e-30
        rhoS = 2000
        g0 = gl.sedSt.gs0(sedF)
        dudx = random.random() + 1e-30
        dvdy = random.random() + 1e-30
        dwdz = random.random() + 1e-30
        test = gl.sedSt.gamma_s(sedF,rhoS,theta , theta_np1, dudx, dvdy, dwdz)
        divU = dudx + dvdy + dwdz
        valid =  - 3. *  (1. - gl.eR**2) * sedF**2 * rhoS * g0 * theta_np1 * ( (sqrt(theta)/sqrt(np.pi)) * (4./gl.grain) - divU) * (2./(3. * rhoS * sedF))
        self.assertTrue(round(test,f) == round(valid,f))

    def testdgamma_s_dtheta(self):
        gl=GlobalVariables()
        import random
        sqrt = np.sqrt
        f = 10        
        # Setting 0 t_c
        sedF = 0.3
        theta = random.random() + 1e-30
        rhoS = 2000
        g0 = gl.sedSt.gs0(sedF)
        dudx = random.random() + 1e-30
        dvdy = random.random() + 1e-30
        dwdz = random.random() + 1e-30

        test = gl.sedSt.dgamma_s_dtheta(sedF,rhoS,theta , dudx, dvdy, dwdz)
        divU = dudx + dvdy + dwdz
        
        valid =  - 3. *  (1. - gl.eR**2) * sedF**2 * rhoS * g0  * ( (sqrt(theta)/sqrt(np.pi)) * (4./gl.grain) - divU) * (2./(3. * rhoS * sedF))
        self.assertTrue(round(test,f) == round(valid,f))

    def testJint1(self):
        gl=GlobalVariables()
        import random
        f = 8
        uf = np.array([5.,4.],"d")
        us = np.array([1.,1.],"d") 
        gradC=np.array([0.1,0.1])
        rhoS = 2000
        rhoF = 1000
        nu = 1e-4
        nuT = 1e-2
        sedF = 0.3
# Setting 0 t_c
        theta_n = random.random() + 1e-30
        kappa_n = random.random() + 1e-30
        epsilon_n = random.random() + 1e-30


        beta = gl.sedSt.betaCoeff(sedF, uf, us, nu)      
        t_p =rhoS/ beta
        l_c = np.sqrt(np.pi)*gl.grain / (24.*sedF * gl.sedSt.gs0(sedF))
        t_cl = min(l_c/np.sqrt(theta_n),0.165*kappa_n/epsilon_n)
        aa = 1/( 1. + t_p/t_cl)

       
        es1 = 2.*beta * rhoS*(1-aa)*sedF*kappa_n/((1-sedF)*rhoF)

        UgradC = np.dot((uf - us),gradC)

        es2  = beta * rhoF * nuT * UgradC / ((1-sedF)*rhoF)

        test = gl.sedSt.jint1(sedF,uf,us,rhoS, kappa_n,epsilon_n, theta_n, nu)

        self.assertTrue(round(test,f) ==round(2*aa*beta*sedF*kappa_n*(2./(3.*sedF*rhoS)),f ))


    def testJint2(self):
        gl=GlobalVariables()
        import random
        f = 8
        uf = np.array([5.,4.],"d")
        us = np.array([1.,1.],"d") 
        gradC=np.array([0.1,0.1])
        rhoS = 2000
        rhoF = 1000
        nu = 1e-4
        nuT = 1e-2
        sedF = 0.3
# Setting 0 t_c
        theta_n = random.random() + 1e-30

        beta = gl.sedSt.betaCoeff(sedF, uf, us, nu)      

       
        test = gl.sedSt.jint2(sedF,uf,us,rhoS, theta_n, nu)

        self.assertTrue(round(test,f) ==round(-3*beta*sedF*theta_n*(2./(3.*sedF*rhoS)),f ))

    def testJint2dTheta(self):
        gl=GlobalVariables()
        import random
        f = 8
        uf = np.array([5.,4.],"d")
        us = np.array([1.,1.],"d") 
        gradC=np.array([0.1,0.1])
        rhoS = 2000
        rhoF = 1000
        nu = 1e-4
        nuT = 1e-2
        sedF = 0.3
# Setting 0 t_c
        theta_n = random.random() + 1e-30

        beta = gl.sedSt.betaCoeff(sedF, uf, us, nu)      

       
        test = gl.sedSt.djint2_dtheta(sedF,uf,us,rhoS, nu)

        self.assertTrue(round(test,f) ==round(-3*beta*sedF*(2./(3.*sedF*rhoS)),f ))


    def testK_diff(self):
        gl=GlobalVariables()
        import random
        sqrt = np.sqrt
        f = 10        
        # Setting 0 t_c
        sedF = 0.3
        theta = random.random() + 1e-30
        rhoS = 2000
        test = gl.sedSt.k_diff(sedF,rhoS,theta)
        g0 = gl.sedSt.gs0(sedF)
        valid = rhoS * gl.grain * sqrt(theta) * ( 2. *sedF**2 * g0 * (1. + gl.eR) / sqrt(np.pi) + (9./16) *sedF**2 * g0 * (1. + gl.eR) * sqrt(np.pi) +  (15./16.) *sedF * sqrt(np.pi) +  (25./64.) * sqrt(np.pi)/((1+gl.eR)*g0))
        
        self.assertTrue(round(test,f) == round(valid,f))



    def testMint(self):
        gl=GlobalVariables()
        uf = np.array([5.,4.],"d")
        us = np.array([1.,1.],"d") 
        gradc = np.array([0.1,0.1],"d") 
        sedF = 0.205
        nu = 1e-4
        nuT = 1e-2
        beta = gl.sedSt.betaCoeff(sedF, uf, us, nu)      
        mint = gl.sedSt.mInt(sedF, uf, us, uf, us , nu, nuT, gradc)
        self.assertTrue(round(mint.all(),10) == round((-sedF*beta*(uf-us) - sedF * beta * gradc * nuT / gl.sigmaC).all() , 10))
    def testdMintdUf(self):
        gl=GlobalVariables()
        uf = np.array([5.,4.],"d")
        us = np.array([1.,1.],"d") 
        gradc = np.array([0.1,0.1],"d") 
        sedF = 0.205
        nu = 1e-4
        nuT = 1e-2
        beta = gl.sedSt.betaCoeff(sedF, uf, us, nu)      
        mint = gl.sedSt.dmInt_duFluid(sedF, uf, us , nu)
        self.assertTrue(round(mint,10) == round(  - sedF*beta , 10))
    def testdMintdUs(self):
        gl=GlobalVariables()
        uf = np.array([5.,4.],"d")
        us = np.array([1.,1.],"d") 
        gradc = np.array([0.1,0.1],"d") 
        sedF = 0.205
        nu = 1e-4
        nuT = 1e-2
        beta = gl.sedSt.betaCoeff(sedF, uf, us, nu)      
        mint = gl.sedSt.dmInt_duSolid(sedF, uf, us , nu)
        self.assertTrue(round(mint,10) == round(  sedF*beta , 10))



if __name__ == '__main__':
    unittest.main(verbosity=2)

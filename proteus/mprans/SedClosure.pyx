import numpy
import cython
cimport numpy
from cpython cimport array
# pull in C++ definitions and declare interface
cdef extern from "mprans/SedClosure.h" namespace "proteus":

    cdef cppclass cppHsuSedStress2D:
        double aDarcy_
        double betaForch_
        double grain_
        double packFraction_
        double packMargin_
        double frFraction_
        double maxFraction_
        double sigmaC_
        double C3e_
        double C4e_
        double eR_
        double fContact_
        double mContact_
        double nContact_
        double angFriction_
        double small_
        double notSoLarge_
        double large_
        cppHsuSedStress2D(
            double aDarcy,  # darcy parameter for drag term. Default value from Ergun (1952) is 150
            double betaForch,  # forchheimer parameter for drag term. Default value from Ergun (1952) is 1.75
            double grain,  # Grain size, default assumed as d50
            double packFraction,  # Critical volume fraction for switching the drag relation 0.2 by default, see Chen and Hsu 2014
            double packMargin,
            double maxFraction,
            double frFraction,
            double sigmaC,
            double C3e,
            double C4e,
            double eR,
            double fContact,
            double mContact,
            double nContact,
            double angFriction

        )
        double betaCoeff(
            double sedF,  # Sediment fraction
            double rhoFluid,
            double * uFluid,  # Fluid velocity
            double * uSolid,  # Sediment velocity
            double nu  # Kinematic viscosity
        )
        double gs0(
            double sedF  # Sediment fraction
        )
        double kappa_sed1(double sedF,
                          double rhoFluid,
                          double rhoSolid,
                          double * uFluid,
                          double * uSolid,
                          double * gradC,
                          double nu,
                          double theta_n,
                          double kappa_n,
                          double kappa_np1,
                          double epsilon_n,
                          double nuT_n)
        double kappa_sed2(double sedF,
                          double rhoFluid,
                          double rhoSolid,
                          double * uFluid,
                          double * uSolid,
                          double * gradC,
                          double nu,
                          double theta_n,
                          double kappa_n,
                          double epsilon_n,
                          double nuT_n)
        double dkappa_sed1_dk(double sedF,
                              double rhoFluid,
                              double rhoSolid,
                              double * uFluid,
                              double * uSolid,
                              double * gradC,
                              double nu,
                              double theta_n,
                              double kappa_n,
                              double epsilon_n,
                              double nuT_n)
        double eps_sed(double sedF,
                       double rhoFluid,
                       double rhoSolid,
                       double * uFluid,
                       double * uSolid,
                       double * gradC,
                       double nu,
                       double theta_n,
                       double kappa_n,
                       double epsilon_n,
                       double epsilon_np1,
                       double nuT_n)
        double deps_sed_deps(double sedF,
                             double rhoFluid,
                             double rhoSolid,
                             double * uFluid,
                             double * uSolid,
                             double * gradC,
                             double nu,
                             double theta_n,
                             double kappa_n,
                             double epsilon_n,
                             double nuT_n)

        double psc(
            double sedF,
            double rhoSolid,
            double theta_n)
        double psc_term(
            double sedF,
            double rhoSolid,
            double theta_np1,
            double du_dx,
            double dv_dy,
            double dw_dz)

        double dpsc_term_dtheta(double sedF,
                                double rhoSolid,
                                double du_dx,
                                double dv_dy,
                                double dw_dz)

        double mu_sc(double sedF,
                     double rhoSolid,
                     double theta)
        double l_sc(double sedF,
                    double rhoSolid,
                    double theta)

        double tausc_term_theta(
            double sedF,
            double rhoSolid,
            double theta_n,
            double du_dx,
            double du_dy,
            double du_dz,
            double dv_dx,
            double dv_dy,
            double dv_dz,
            double dw_dx,
            double dw_dy,
            double dw_dz)

        double  gamma_s(double sedF,
                        double rhoSolid,
                        double theta_n,
                        double theta_np1,
                        double du_dx,
                        double dv_dy,
                        double dw_dz)

        double  dgamma_s_dtheta(double sedF,
                                double rhoSolid,
                                double theta_np1,
                                double du_dx,
                                double dv_dy,
                                double dw_dz)
        double  jint1(double sedF,
                      double rhoFluid,
                      double rhoSolid,
                      double * uFluid,
                      double * uSolid,
                      double kappa,
                      double epsilon,
                      double theta_n,
                      double nu)

        double  jint2(double sedF,
                      double rhoFluid,
                      double rhoSolid,
                      double * uFluid,
                      double * uSolid,
                      double theta,
                      double nu)

        double  djint2_dtheta(double sedF,
                              double rhoFluid,
                              double rhoSolid,
                              double * uFluid,
                              double * uSolid,
                              double nu)

        double k_diff(double sedF, double rhoSolid,  double theta)

        double p_friction(double sedF)
        double mu_fr(double sedF,
                     double du_dx,
                     double du_dy,
                     double du_dz,
                     double dv_dx,
                     double dv_dy,
                     double dv_dz,
                     double dw_dx,
                     double dw_dy,
                     double dw_dz)

        void  mIntFluid(double * mint2,
                        double sedF,         # Sediment fraction
                        double rhoFluid,
                        double * uFluid_n,  # Fluid velocity
                        double * uSolid_n,  # Sediment velocity
                        double * uFluid_np1,  # Fluid
                        double nu,  # Kinematic viscosity
                        double nuT,
                        double * gradc
                        )
        void  mIntSolid(double * mint2,
                        double sedF,         # Sediment fraction
                        double rhoFluid,
                        double * uFluid_n,  # Fluid velocity
                        double * uSolid_n,  # Sediment velocity
                        double * uFluid_np1,  # Fluid
                        double nu,  # Kinematic viscosity
                        double nuT,
                        double * gradc
                        )

        void  mIntgradC(double * mint2,
                        double sedF,         # Sediment fraction
                        double rhoFluid,
                        double * uFluid_n,  # Fluid velocity
                        double * uSolid_n,  # Sediment velocity
                        double nu,  # Kinematic viscosity
                        double nuT,
                        double * gradc
                        )

        double dmInt_duFluid(
            double sedF,  # Sediment fraction
            double rhoFluid,
            double * uFluid_n,  # Fluid velocity
            double * uSolid_n,  # Sediment velocity
            double nu  # Kinematic viscosity
        )
        double dmInt_duSolid(
            double sedF,  # Sediment fraction
            double rhoFluid,
            double * uFluid_n,  # Fluid velocity
            double * uSolid_n,  # Sediment velocity
            double nu  # Kinematic viscosity
        )
        double p_s(double sedF,
                   double rhoSolid,
                   double theta,
                   double du_dx,
                   double du_dy,
                   double du_dz,
                   double dv_dx,
                   double dv_dy,
                   double dv_dz,
                   double dw_dx,
                   double dw_dy,
                   double dw_dz)

# define the way we want to present to Python
cdef class HsuSedStress:
    cdef  cppHsuSedStress2D * thisptr

    def __cinit__(self, aDarcy, betaForch, grain, packFraction, packMargin, maxFraction, frFraction, sigmaC, C3e, C4e, eR, fContact, mContact, nContact, angFriction):
        """ Class for caclulating sediment / fluid momentum transfer, see Chen and Hsu, CACR 14-08, A Multidimensional TwoPhase Eulerian Model for Sediment Transport TwoPhaseEulerSedFoam (Version 1.0)
        http://www.coastal.udel.edu/~thsu/simulation_data_files/CACR-14-08.pdf
        param: aDarcy: Darcy parameter for drag term [-]. Default value from Ergun (1952) is 150
        param: betaForch: Forchheimer parameter for drag term [-]. Default value from Ergun (1952) is 1.75
        param: grain: Grain size, default assumed as d50 [L]
        param: packFraction : Critical sediment fraction [-] for switching the drag relation 0.2 by default, see Chen and Hsu 2014, equation (7)
        param: packMargin : [-] For packFraction \pm packMargin where the two braches in equation (7) are blended with linear weighting. Currently no information on the default value of this """
        self.thisptr = new cppHsuSedStress2D(aDarcy, betaForch, grain, packFraction, packMargin, maxFraction, frFraction, sigmaC, C3e, C4e, eR, fContact,  mContact, nContact, angFriction)

    @property
    def aDarcy(self):
        return self.thisptr.aDarcy_

    @property
    def betaForch(self):
        return self.thisptr.betaForch_

    @property
    def grain(self):
        return self.thisptr.grain_

    @property
    def packFraction(self):
        return self.thisptr.packFraction_

    @property
    def packMargin(self):
        return self.thisptr.packMargin_

    @property
    def maxFraction(self):
        return self.thisptr.maxFraction_

    @property
    def frFraction(self):
        return self.thisptr.frFraction_

    @property
    def sigmaC(self):
        return self.thisptr.sigmaC_

    @property
    def C3e(self):
        return self.thisptr.C3e_

    @property
    def C4e(self):
        return self.thisptr.C4e_

    @property
    def eR(self):
        return self.thisptr.eR_

    @property
    def fContact(self):
        return self.thisptr.eR_

    @property
    def mContact(self):
        return self.thisptr.mContact_

    @property
    def nContact(self):
        return self.thisptr.nContact_

    @property
    def angFriction(self):
        return self.thisptr.angFriction_

    def __dealloc__(self):
        del self.thisptr

    def betaCoeff(self,
                  sedF,
                  rhoFluid,
                  uFluid,
                  uSolid,
                  nu):
        """ Function for calculating equation (7) from Chen and Hsu, CACR 14-08, A Multidimensional TwoPhase Eulerian Model for Sediment Transport TwoPhaseEulerSedFoam (Version 1.0)
        http://www.coastal.udel.edu/~thsu/simulation_data_files/CACR-14-08.pdf
        param: sedF: Sediment fraction [-]
        param: uFluid: Fluid velocity vector [L/T]
        param: uSolid: Solid velocity vector [L/T]
        param: nu  : Fluid kinematic viscosity [L^2/T]
        """

        cython.declare(UF=double[2])
        cython.declare(US=double[2])
        for ii in range(2):
            UF[ii] = uFluid[ii]
            US[ii] = uSolid[ii]
        beta = self.thisptr.betaCoeff(sedF, rhoFluid, UF, US, nu)
        return beta

    def gs0(self, sedF):
        """ Radial distribution function for collision closure,  equation (2.31) from  Hsu et al 2004 'On two-phase sediment transport:
        sheet flow of massive particles', Proc. Royal Soc. Lond A 460, pp 2223-2250
        http://www.coastal.udel.edu/~thsu/simulation_data_files/CACR-14-08.pdf
        param: sedF: Sediment fraction [-]
        """
        return self.thisptr.gs0(sedF)

    def kappa_sed1(self,
                   sedF,
                   rhoFluid,
                   rhoSolid,
                   numpy.ndarray uFluid,
                   numpy.ndarray uSolid,
                   numpy.ndarray gradC,
                   nu,
                   theta_n,
                   kappa_n,
                   kappa_np1,
                   epsilon_n,
                   nuT_n):
        return self.thisptr.kappa_sed1(sedF,
                                       rhoFluid,
                                       rhoSolid,
                                       < double * > uFluid.data,
                                       < double * >  uSolid.data,
                                       < double * >  gradC.data,
                                       nu,
                                       theta_n,
                                       kappa_n,
                                       kappa_np1,
                                       epsilon_n,
                                       nuT_n)

    def dkappa_sed1_dk(self,
                       sedF,
                       rhoFluid,
                       rhoSolid,
                       numpy.ndarray uFluid,
                       numpy.ndarray uSolid,
                       numpy.ndarray gradC,
                       nu,
                       theta_n,
                       kappa_n,
                       epsilon_n,
                       nuT_n):
        return self.thisptr.dkappa_sed1_dk(sedF,
                                           rhoFluid,
                                           rhoSolid,
                                           < double * > uFluid.data,
                                           < double * >  uSolid.data,
                                           < double * >  gradC.data,
                                           nu,
                                           theta_n,
                                           kappa_n,
                                           epsilon_n,
                                           nuT_n)

    def kappa_sed2(self,
                   sedF,
                   rhoFluid,
                   rhoSolid,
                   numpy.ndarray uFluid,
                   numpy.ndarray uSolid,
                   numpy.ndarray gradC,
                   nu,
                   theta_n,
                   kappa_n,
                   epsilon_n,
                   nuT_n):
        return self.thisptr.kappa_sed2(sedF,
                                       rhoFluid,
                                       rhoSolid,
                                       < double * > uFluid.data,
                                       < double * >  uSolid.data,
                                       < double * >  gradC.data,
                                       nu,
                                       theta_n,
                                       kappa_n,
                                       epsilon_n,
                                       nuT_n)

    def eps_sed(self,
                sedF,
                rhoFluid,
                rhoSolid,
                numpy.ndarray uFluid,
                numpy.ndarray uSolid,
                numpy.ndarray gradC,
                nu,
                theta_n,
                kappa_n,
                epsilon_n,
                epsilon_np1,
                nuT_n):
        return self.thisptr.eps_sed(sedF,
                                    rhoFluid,
                                    rhoSolid,
                                    < double * > uFluid.data,
                                    < double * >  uSolid.data,
                                    < double * >  gradC.data,
                                    nu,
                                    theta_n,
                                    kappa_n,
                                    epsilon_n,
                                    epsilon_np1,
                                    nuT_n)

    def deps_sed_deps(self,
                      sedF,
                      rhoFluid,
                      rhoSolid,
                      numpy.ndarray uFluid,
                      numpy.ndarray uSolid,
                      numpy.ndarray gradC,
                      nu,
                      theta_n,
                      kappa_n,
                      epsilon_n,
                      nuT_n):
        return self.thisptr.deps_sed_deps(sedF,
                                          rhoFluid,
                                          rhoSolid,
                                          < double * > uFluid.data,
                                          < double * >  uSolid.data,
                                          < double * >  gradC.data,
                                          nu,
                                          theta_n,
                                          kappa_n,
                                          epsilon_n,
                                          nuT_n)

    def psc(self,
            sedF,
            rhoSolid,
            theta):
        return self.thisptr.psc(sedF, rhoSolid, theta)

    def psc_term(self,
                 sedF,
                 rhoSolid,
                 theta_np1,
                 du_dx,
                 dv_dy,
                 dw_dz):

        return self.thisptr.psc_term(
            sedF,
            rhoSolid,
            theta_np1,
            du_dx,
            dv_dy,
            dw_dz)

    def dpsc_term_dtheta(self,
                         sedF,
                         rhoSolid,
                         du_dx,
                         dv_dy,
                         dw_dz):

        return self.thisptr.dpsc_term_dtheta(
            sedF,
            rhoSolid,
            du_dx,
            dv_dy,
            dw_dz)

    def mu_sc(self,
              sedF,
              rhoSolid,
              theta):
        return self.thisptr.mu_sc(sedF, rhoSolid, theta)

    def mu_fr(self, sedF, du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz):
        return self.thisptr.mu_fr(sedF, du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz)

    def l_sc(self,
             sedF,
             rhoSolid,
             theta):
        return self.thisptr.l_sc(sedF, rhoSolid, theta)

    def tausc_term_theta(self,
                         sedF,
                         rhoSolid,
                         theta_n,
                         du_dx,
                         du_dy,
                         du_dz,
                         dv_dx,
                         dv_dy,
                         dv_dz,
                         dw_dx,
                         dw_dy,
                         dw_dz):
        return self.thisptr.tausc_term_theta(sedF, rhoSolid, theta_n, du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz)

    def gamma_s(self, sedF,
                rhoSolid,
                theta_n,
                theta_np1,
                du_dx,
                dv_dy,
                dw_dz):
        return self.thisptr.gamma_s(sedF, rhoSolid, theta_n, theta_np1,  du_dx, dv_dy, dw_dz)

    def dgamma_s_dtheta(self, sedF,
                        rhoSolid,
                        theta_n,
                        du_dx,
                        dv_dy,
                        dw_dz):
        return self.thisptr.dgamma_s_dtheta(sedF, rhoSolid, theta_n,  du_dx, dv_dy, dw_dz)

    def jint1(self, sedF, rhoFluid, rhoSolid, numpy.ndarray uFluid, numpy.ndarray uSolid, kappa, epsilon, theta_n, nu):
        return self.thisptr.jint1(sedF, rhoFluid, rhoSolid, < double * > uFluid.data, < double * > uSolid.data,  kappa, epsilon, theta_n, nu)

    def jint2(self, sedF, rhoFluid, rhoSolid, numpy.ndarray uFluid, numpy.ndarray uSolid,  theta, nu):
        return self.thisptr.jint2(sedF, rhoFluid, rhoSolid, < double * > uFluid.data, < double * > uSolid.data,  theta, nu)

    def djint2_dtheta(self, sedF, rhoFluid, rhoSolid, numpy.ndarray uFluid, numpy.ndarray uSolid, nu):
        return self.thisptr.djint2_dtheta(sedF, rhoFluid, rhoSolid, < double * > uFluid.data, < double * > uSolid.data,  nu)

    def k_diff(self,  sedF, rhoSolid, theta):
        return self.thisptr.k_diff(sedF, rhoSolid, theta)

    def p_friction(self, sedF):
        return self.thisptr.p_friction(sedF)

    def mIntFluid(self,
                  sedF,
                  rhoFluid,
                  numpy.ndarray uFluid_n,
                  numpy.ndarray uSolid_n,
                  numpy.ndarray uFluid_np1,
                  nu,
                  nuT,
                  numpy.ndarray gradc
                  ):

        cython.declare(xx=cython.double[2])
        for ii in range(2):
            xx[ii] = 0.

        self.thisptr.mIntFluid(xx, sedF, rhoFluid,

                               < double * > uFluid_n.data,
                               < double * > uSolid_n.data,
                               < double * > uFluid_np1.data,
                               nu,
                               nuT,
                               < double * > gradc.data)
        mint = numpy.zeros(2,)
        for ii in range(2):
            mint[ii] = xx[ii]
        return mint

    def mIntSolid(self,
                  sedF,
                  rhoFluid,
                  numpy.ndarray uFluid_n,
                  numpy.ndarray uSolid_n,
                  numpy.ndarray uSolid_np1,
                  nu,
                  nuT,
                  numpy.ndarray gradc
                  ):

        cython.declare(xx=cython.double[2])
        for ii in range(2):
            xx[ii] = 0.

        self.thisptr.mIntSolid(xx, sedF, rhoFluid,
                               < double * > uFluid_n.data,
                               < double * > uSolid_n.data,
                               < double * > uSolid_np1.data,
                               nu,
                               nuT,
                               < double * > gradc.data)
        mint = numpy.zeros(2,)
        for ii in range(2):
            mint[ii] = xx[ii]
        return mint

    def mIntgradC(self,
                  sedF,
                  rhoFluid,
                  numpy.ndarray uFluid_n,
                  numpy.ndarray uSolid_n,
                  nu,
                  nuT,
                  numpy.ndarray gradc
                  ):

        cython.declare(xx=cython.double[2])
        for ii in range(2):
            xx[ii] = 0.

        self.thisptr.mIntgradC(xx, sedF, rhoFluid,
                               < double * > uFluid_n.data,
                               < double * > uSolid_n.data,
                               nu,
                               nuT,
                               < double * > gradc.data)
        mint = numpy.zeros(2,)
        for ii in range(2):
            mint[ii] = xx[ii]
        return mint

    def dmInt_duFluid(self,
                      sedF,  # Sediment fraction
                      rhoFluid,
                      numpy.ndarray uFluid_n,  # Fluid velocity
                      numpy.ndarray uSolid_n,  # Sediment velocity
                      nu):  # Kinematic viscosity

        return self.thisptr.dmInt_duFluid(sedF, rhoFluid,
                                          < double * > uFluid_n.data,
                                          < double * > uSolid_n.data,
                                          nu)

    def dmInt_duSolid(self,
                      sedF,  # Sediment fraction
                      rhoFluid,
                      numpy.ndarray uFluid_n,  # Fluid velocity
                      numpy.ndarray uSolid_n,  # Sediment velocity
                      nu):  # Kinematic viscosity

        return self.thisptr.dmInt_duSolid(sedF, rhoFluid,
                                          < double * > uFluid_n.data,
                                          < double * > uSolid_n.data,
                                          nu)

    def p_s(self, sedF,  rhoSolid, theta, du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz):
        return self.thisptr.p_s(sedF,  rhoSolid, theta, du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz)

import numpy
cimport numpy
from cpython cimport array
#pull in C++ definitions and declare interface
cdef extern from "mprans/SedClosure.h" namespace "proteus":
    cdef cppclass cppHsuSedStress:
        cppHsuSedStress(
		 double aDarcy, # darcy parameter for drag term. Default value from Ergun (1952) is 150
		 double betaForch, # forchheimer parameter for drag term. Default value from Ergun (1952) is 1.75
		 double grain, # Grain size, default assumed as d50
		 double packFraction, #Critical volume fraction for switching the drag relation 0.2 by default, see Chen and Hsu 2014
		 double packMargin, #
                 double sigmaC,
                 double C3e,
                 double C4e
            
		 )
        double betaCoeff(
                            double sedF, # Sediment fraction
                            double* uFluid, #Fluid velocity
                            double* uSolid, #Sediment velocity
                            double nu #Kinematic viscosity
                           )
        double gs0(
                             double sedF # Sediment fraction
            )
        double kappa_sed(double sedF, 
                       double rhoFluid,
                       double rhoSolid,
                       double* uFluid, 
                       double* uSolid,
                       double* gradC,
                       double nu, 
                       double theta_n,
                       double kappa_n,
                       double epsilon_n,
                       double nuT_n)                          
        double eps_sed(double sedF, 
                       double rhoFluid,
                       double rhoSolid,
                       double* uFluid, 
                       double* uSolid,
                       double* gradC,
                       double nu, 
                       double theta_n,
                       double kappa_n,
                       double epsilon_n,
                       double nuT_n)                          


        double*  mInt(
                            double sedF,         # Sediment fraction
                            double* uFluid_np1,  #Fluid velocity
                            double* uSolid_np1,  #Sediment velocity
                            double* uFluid_n,    #Fluid velocity
                            double* uSolid_n,    #Sediment velocity
                            double nu,           #Kinematic viscosity
                            double nuT,
                            double* gradc
                           )
        double dmInt_duFluid(
                            double sedF, # Sediment fraction
                            double* uFluid_n, #Fluid velocity
                            double* uSolid_n, #Sediment velocity
                            double nu #Kinematic viscosity
                           )
        double dmInt_duSolid(
                            double sedF, # Sediment fraction
                            double* uFluid_n, #Fluid velocity
                            double* uSolid_n, #Sediment velocity
                            double nu #Kinematic viscosity
                           )

#define the way we want to present to Python
cdef class HsuSedStress:
    cdef  cppHsuSedStress* thisptr
    def __cinit__(self, aDarcy, betaForch, grain, packFraction,packMargin, sigmaC, C3e, C4e ):
        """ Class for caclulating sediment / fluid momentum transfer, see Chen and Hsu, CACR 14-08, A Multidimensional TwoPhase Eulerian Model for Sediment Transport TwoPhaseEulerSedFoam (Version 1.0) 
        http://www.coastal.udel.edu/~thsu/simulation_data_files/CACR-14-08.pdf
        param: aDarcy: Darcy parameter for drag term [-]. Default value from Ergun (1952) is 150
        param: betaForch: Forchheimer parameter for drag term [-]. Default value from Ergun (1952) is 1.75
        param: grain: Grain size, default assumed as d50 [L]
        param: packFraction : Critical sediment fraction [-] for switching the drag relation 0.2 by default, see Chen and Hsu 2014, equation (7)
        param: packMargin : [-] For packFraction \pm packMargin where the two braches in equation (7) are blended with linear weighting. Currently no information on the default value of this """
        self.thisptr = new cppHsuSedStress( aDarcy, betaForch, grain, packFraction, packMargin, sigmaC, C3e, C4e)
    def __dealloc__(self):
        del self.thisptr
    def betaCoeff(self, 
                     sedF,  
                     numpy.ndarray uFluid, 
                     numpy.ndarray uSolid, 
                     nu):
        """ Function for calculating equation (7) from Chen and Hsu, CACR 14-08, A Multidimensional TwoPhase Eulerian Model for Sediment Transport TwoPhaseEulerSedFoam (Version 1.0) 
        http://www.coastal.udel.edu/~thsu/simulation_data_files/CACR-14-08.pdf
        param: sedF: Sediment fraction [-]
        param: uFluid: Fluid velocity vector [L/T]
        param: uSolid: Solid velocity vector [L/T]
        param: nu  : Fluid kinematic viscosity [L^2/T]
        """
        return self.thisptr.betaCoeff(sedF, 
                                  < double * > uFluid.data,
                                  < double * > uSolid.data, 
                                  nu)

    def gs0(self,sedF):
        """ Radial distribution function for collision closure,  equation (2.31) from  Hsu et al 2004 'On two-phase sediment transport:
        sheet flow of massive particles', Proc. Royal Soc. Lond A 460, pp 2223-2250 
        http://www.coastal.udel.edu/~thsu/simulation_data_files/CACR-14-08.pdf
        param: sedF: Sediment fraction [-]        
        """
        return self.thisptr.gs0(sedF) 
    def kappa_sed(self,
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
        return self.thisptr.kappa_sed(sedF,  
                    rhoFluid,
                    rhoSolid,
                  < double *> uFluid.data, 
                  < double *>  uSolid.data, 
                  < double *>  gradC.data, 
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
                nuT_n):
        return self.thisptr.eps_sed(sedF,  
                    rhoFluid,
                    rhoSolid,
                  < double *> uFluid.data, 
                  < double *>  uSolid.data, 
                  < double *>  gradC.data, 
                  nu,
                  theta_n,
                  kappa_n,
                epsilon_n,
                nuT_n)


    def  mInt(self, 
                     sedF,  
                     numpy.ndarray uFluid_np1, 
                     numpy.ndarray uSolid_np1, 
                     numpy.ndarray uFluid_n, 
                     numpy.ndarray uSolid_n,                                                      
                     nu, 
                     nuT,
                     numpy.ndarray gradc
   ):
        
        mint = numpy.zeros(len(uFluid_n),"d")
        cdef double* carr = self.thisptr.mInt(sedF, 
                                  < double * > uFluid_np1.data,
                                  < double * > uSolid_np1.data, 
                                  < double * > uFluid_n.data,
                                  < double * > uSolid_n.data, 
                                         nu,
                                         nuT, 
                                  < double * > gradc.data)
        for ii in range(len(mint)):
            mint[ii] = carr[ii]
        return mint
                                 

    def dmInt_duFluid(self,
                            sedF, # Sediment fraction
                            numpy.ndarray uFluid_n, #Fluid velocity
                            numpy.ndarray uSolid_n, #Sediment velocity
                            nu): #Kinematic viscosity
                           
         return self.thisptr.dmInt_duFluid(sedF, 
                                  < double * > uFluid_n.data,
                                  < double * > uSolid_n.data, 
                                         nu)

        

    def dmInt_duSolid(self,
                            sedF, # Sediment fraction
                            numpy.ndarray uFluid_n, #Fluid velocity
                            numpy.ndarray uSolid_n, #Sediment velocity
                            nu): #Kinematic viscosity
                           

         return self.thisptr.dmInt_duSolid(sedF, 
                                  < double * > uFluid_n.data,
                                  < double * > uSolid_n.data, 
                                         nu)



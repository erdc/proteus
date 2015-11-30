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
                 double sigmaC
		 )
        double granularDrag(
                            double sedF, # Sediment fraction
                            double* uFluid, #Fluid velocity
                            double* uSolid, #Sediment velocity
                            double nu #Kinematic viscosity
                           )
        double turbSusp(
                            double sedF, # Sediment fraction
                            double* uFluid, #Fluid velocity
                            double* uSolid, #Sediment velocity
                            double nu, #Kinematic viscosity
                            double nuT
                           )

#define the way we want to present to Python
cdef class HsuSedStress:
    cdef  cppHsuSedStress* thisptr
    def __cinit__(self, aDarcy, betaForch, grain, packFraction,packMargin, sigmaC ):
        """ Class for caclulating sediment / fluid momentum transfer, see Chen and Hsu, CACR 14-08, A Multidimensional TwoPhase Eulerian Model for Sediment Transport TwoPhaseEulerSedFoam (Version 1.0) 
        http://www.coastal.udel.edu/~thsu/simulation_data_files/CACR-14-08.pdf
        param: aDarcy: Darcy parameter for drag term [-]. Default value from Ergun (1952) is 150
        param: betaForch: Forchheimer parameter for drag term [-]. Default value from Ergun (1952) is 1.75
        param: grain: Grain size, default assumed as d50 [L]
        param: packFraction : Critical sediment fraction [-] for switching the drag relation 0.2 by default, see Chen and Hsu 2014, equation (7)
        param: packMargin : [-] For packFraction \pm packMargin where the two braches in equation (7) are blended with linear weighting. Currently no information on the default value of this """
        self.thisptr = new cppHsuSedStress( aDarcy, betaForch, grain, packFraction, packMargin, sigmaC)
    def __dealloc__(self):
        del self.thisptr
    def granularDrag(self, 
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
        return self.thisptr.granularDrag(sedF, 
                                  < double * > uFluid.data,
                                  < double * > uSolid.data, 
                                  nu)
    def turbSusp(self, 
                     sedF,  
                     numpy.ndarray uFluid, 
                     numpy.ndarray uSolid, 
                     nu, 
                     nuT):
        """ Function for calculating equation (7) from Chen and Hsu, CACR 14-08, A Multidimensional TwoPhase Eulerian Model for Sediment Transport TwoPhaseEulerSedFoam (Version 1.0) 
        http://www.coastal.udel.edu/~thsu/simulation_data_files/CACR-14-08.pdf
        param: sedF: Sediment fraction [-]
        param: uFluid: Fluid velocity vector [L/T]
        param: uSolid: Solid velocity vector [L/T]
        param: nu  : Fluid kinematic viscosity [L^2/T]
        """
        return self.thisptr.turbSusp(sedF, 
                                  < double * > uFluid.data,
                                  < double * > uSolid.data, 
                                         nu,
                                         nuT)

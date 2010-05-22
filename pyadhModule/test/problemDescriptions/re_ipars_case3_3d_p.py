from pyadh import *
from pyadh.default_p import *

nd = 3

iparsInput = InputTranslators.Ipars("iparsinput_case3",layered=True)
L = iparsInput.L
#polyfile = iparsInput.polyfile
polyfile = "iparsinput_case2_wwells"

analyticalSolution = None

viscosity     = 8.9e-4  #kg/(m*s)
density       = 998.2   #kg/m^3
gravity       = 9.8     #m/s^2
beta          = density*gravity*4.524e-10
m_per_s_by_m_per_d = 1.1574074e-5
ktd = [9.0,6.0,1.0,4.0,1.0]
kvd = [2.7,1.8,0.3,1.2,0.3]
ktm2 = [k*0.987e-12 for k in ktd]
kvm2 = [k*0.987e-12 for k in kvd]
#ktd = [9.0,6.0,1.0,4.0,1.0]
ktd = [1.0,4.0,1.0,6.0,9.0]
#kvd = [2.7,1.8,0.3,1.2,0.3]
kvd = [0.3,1.2,0.3,1.8,2.7]
lengthScale=1.0
timeScale=1.0
kt_dls  = [(timeScale*density*gravity*k/(viscosity*lengthScale))/m_per_s_by_m_per_d for k in ktm2]
kv_dls  = [(timeScale*density*gravity*k/(viscosity*lengthScale))/m_per_s_by_m_per_d for k in kvm2]
dimensionless_density  = 1.0
dimensionless_gravity  = numpy.array([0.0,
                                      0.0,
                                      -1.0])
#psiD_psi = [2.0,2.45,6.0,3.0,6.0]
psiD_psi = [6.0,3.0,6.0,2.45,2.0]
psiD_m = [(p*math.sqrt(0.4/18.0)*6894.7573)/(density*gravity) for p in psiD_psi]
coefficients = ConservativeHeadRichardsJLeverettAni(phi_block=numpy.array([0.4,0.4,0.4,0.4,0.4]),
                                                    psiD_block=numpy.array(psiD_m),
                                                    ns_block=numpy.array([1.0/1.5,1.0/1.5,1.0/1.5,1.0/1.5,1.0/1.5,1.0/1.5]),
                                                    nk_block=numpy.array([1.5,1.5,1.5,1.5,1.5,1.5]),
                                                    S_wirr_block=numpy.array([0.02,0.02,0.02,0.02,0.02]),
                                                    S_nwr_block=numpy.array([0.0,0.0,0.0,0.0,0.0]),
                                                    kr0x_block=numpy.array(kt_dls),
                                                    kr0y_block=numpy.array(kt_dls),
                                                    kr0z_block=numpy.array(kv_dls),
                                                    gravity=dimensionless_gravity,
                                                    density=dimensionless_density,
                                                    beta=beta)
mvg_alpha     = 5.47    #1/m
mvg_n         = 4.264
mvg_m         = 1.0 - 1.0/mvg_n
bc_lambda     = mvg_n-1.0
bc_lambda     = 2.0
bc_pd         = 1.0/mvg_alpha
bc_pd         = 0.01
dimensionless_pd0    = bc_pd/lengthScale
k0 = kt_dls[0]
dimensionless_pd    = [dimensionless_pd0*math.sqrt(k0/k) for k in kt_dls]
psiD=0.0
coefficients = ConservativeHeadRichardsJLeverettAni(phi_block=numpy.array([0.4,0.4,0.4,0.4,0.4]),
                                                    psiD_block=numpy.array(dimensionless_pd),
                                                    ns_block=numpy.array([bc_lambda,bc_lambda,bc_lambda,bc_lambda,bc_lambda]),
                                                    nk_block=numpy.array([1.5,1.5,1.5,1.5,1.5]),
                                                    S_wirr_block=numpy.array([0.2,0.2,0.2,0.2,0.2]),
                                                    S_nwr_block=numpy.array([0.0,0.0,0.0,0.0,0.0,0.0]),
                                                    kr0x_block=numpy.array(kt_dls),
                                                    kr0y_block=numpy.array(kt_dls),
                                                    kr0z_block=numpy.array(kv_dls),
                                                    gravity=dimensionless_gravity,
                                                    density=dimensionless_density,
                                                    beta=beta)
#not used
presureWells=[(1.68,4.99),
              (1.68,11.33),
              (1.68,15.03),
              (1.68,19.67),
              (66.38,4.99),
              (66.38,11.33),
              (66.38,15.03),
              (66.38,19.67)]

pressureWellRadius=0.1

#injection wells, used for flux BCs
injectionWells=[(67.28,12.33),
                (67.78,12.33),
                (68.28,12.33),
                (68.78,12.33)]

injectionWellRadius=0.05*0.3048
wellArea = math.pi*injectionWellRadius**2

fluxBoundaryConditions = {}

def inInjectionWell(x):
    if x[2] >= L[2] - 1.0e-8:
        for w in injectionWells:
            if math.sqrt((x[0]-w[0]*0.3048)**2 + (x[1]-w[1]*0.3048)**2) < injectionWellRadius:
                return True
    return False

def getDBC(x,flag):
    if flag in [iparsInput.boundaryFlags['left'],
                iparsInput.boundaryFlags['right']]:
        return lambda x,t: -psiD + x[2]*dimensionless_gravity[2]
    else:
        return None

dirichletConditions = {0:getDBC}

fluxBoundaryConditions = {}

class IC:
    def uOfXT(self,x,t):
        return -psiD + x[2]*dimensionless_gravity[2]

initialConditions = {0:IC()}

fluxBoundaryConditions = {0:'noFlow'}

def waterFlowRate(x,t):
    if t < 2.0:
        return  -(5.0/0.5**2)*0.3048 #3.0 ft^3/day over 0.5 ft^2 area -> m/day
    else:
        return 0.0

def getAFBC(x,flag):
    if flag in [iparsInput.boundaryFlags['left'],
                iparsInput.boundaryFlags['right']]:
        return None
    else:
        return lambda x,t: 0.0

def getDFBC(x,flag):
    if flag in [iparsInput.boundaryFlags['left'],iparsInput.boundaryFlags['right']]:
        return None
#     elif flag == iparsInput.boundaryFlags['top']:
#         if inInjectionWell(x):
#             return waterFlowRate
#         else:
#             return lambda x,t: 0.0           
    elif flag == iparsInput.boundaryFlags['top']:
        return lambda x,t: 0.0           
    elif flag == 7:
        return waterFlowRate
    else:
        return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC}

diffusiveFluxBoundaryConditions = {0:{0:getDFBC}}

T = 5.0

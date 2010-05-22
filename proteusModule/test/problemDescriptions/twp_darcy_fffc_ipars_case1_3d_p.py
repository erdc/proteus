from pyadh import *
from pyadh.default_p import *

nd = 3
iparsInput = InputTranslators.Ipars("iparsinput_case1")
L = iparsInput.L
polyfile = iparsInput.polyfile
print iparsInput.values.keys()    
coefficients = TwophaseDarcyCoefficients.TwophaseDarcy_fc_ff(Ksw=18.0,
                                   rhon=0.002,
                                   rhow=0.433,
                                   g=iparsInput.values['DOWN'],
                                   bc_lambda=1.5,
                                   bc_pd = 2.0,
                                   omega=0.4,
                                   Sw_min=0.2,
                                   Sw_max=1.0,
                                   mun=0.015,
                                   muw=1.0,
                                   psk_model='BCB')

presureWells=[(1.68,4.99),
              (1.68,11.33),
              (1.68,15.03),
              (1.68,19.67),
              (66.38,4.99),
              (66.38,11.33),
              (66.38,15.03),
              (66.38,19.67)]
pressureWellRadius=0.1

#injection wells
injectionWells=[(33.78,12.33),
                (34.28,12.33)]
injectionWellRadius=0.05
injectionWellRadius=5.0
fluxBoundaryConditions = {}

def inInjectionWell(x):
    for w in injectionWells:
        if math.sqrt((x[0]-w[0])**2 + (x[1]-w[1])**2) < injectionWellRadius:
            return True
    return False

def getDBC_sw(x,flag):
    if inInjectionWell(x):
        print "Injection Well BC, Sw=1",x
        return lambda x,t: 1.0 #fully saturated on inflow
    if x[0] in [0.0,L[0]]:
        print "Open Air BC Sw=0.2"
        return lambda x,t: 0.2
def getDBC_psiw(x,flag):
    if x[0] in [0.0,L[0]]:
        print "Open Air BC pw = 14.7",x
        return lambda x,t: 14.7 - x[2]

dirichletConditions = {0:getDBC_sw,1:getDBC_psiw}

class IC_sw:
    def uOfXT(self,x,t):
        return 0.2

class IC_psiw:
    def uOfXT(self,x,t):
        return x[2]

initialConditions = {0:IC_sw(),
                     1:IC_psiw()}


def waterFlowRate(x,t):
    if t < 2.0:
        return -3.0#ft/day
    else:
        return 0.0

def getAFBC_pressureEqn(x,flag):
    return lambda x,t: 0.0

def getDFBC_satEqn(x,flag):
    if inInjectionWell(x):
        return waterFlowRate

def getDFBC_pressureEqn(x,flag):
    if inInjectionWell(x):
        return waterFlowRate

def getDFBC_capillary(x,flag):
    return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {1:getAFBC_pressureEqn}

diffusiveFluxBoundaryConditions = {0:{1:getDFBC_pressureEqn},
                                   1:{1:getDFBC_pressureEqn,
                                      0:getDFBC_capillary}}

T = iparsInput.values['TIMEEND'][0]

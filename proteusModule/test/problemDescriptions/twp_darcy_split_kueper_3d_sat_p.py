from pyadh import *
from pyadh.default_p import *
from pyadh.TwophaseDarcyCoefficients import *
from kueper3d_modelParams import *
name = "twp_darcy_split_kueper_sat_cgv"

useHet = True

phase = 'saturation'
coefficients = TwophaseDarcy_split_saturation(g=g, 
                                              rhon=rhon,
                                              rhow=rhow,
                                              mun    = mun,
                                              muw    = muw,
                                              Ksw=Ksw_types[source_type],
                                              psk_model=model,
                                              vg_alpha = mvg_alpha_types[source_type],
                                              vg_m  = mvg_m_types[source_type],
                                              bc_pd  = bc_pd_types[source_type], 
                                              bc_lambda = bc_lambda_types[source_type],
                                              omega  = omega_types[source_type],
                                              Sw_max = sw_max_types[source_type],
                                              Sw_min = sw_min_types[source_type])
if useHet:
    coefficients.setMaterialTypes(Ksw_types = Ksw_types,
                                  omega_types  = omega_types,
                                  Sw_max_types = sw_max_types,
                                  Sw_min_types = sw_min_types,
                                  bc_lambda_types = bc_lambda_types,
                                  bc_pd_types = bc_pd_types,
                                  vg_alpha_types = mvg_alpha_types,
                                  vg_m_types = mvg_m_types)
        
hour = 3600.0 #[s]
T = 1.0e-2/timeScale


dirichletConditions = {0:getDBC_sw}

initialConditions  = {0:sw_IC()}

fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}


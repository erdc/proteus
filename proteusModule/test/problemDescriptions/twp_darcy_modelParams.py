# Set the model parameters to be given to the impes model for two phase flow. 

#choose in problem specific class now
#model = 'VGM'#'simp'#'VGM'       # Choose from {simp,BCB,BCM,VGB,VGM}
#for now have to have PorouMedia set here
PorousMedia = 'clay-2'#'sand'#'clay-2' 

# Model Parameters
rhow = 997.0        # density wetting      (Water 997 kg/m^3)
rhon = 1.205        # density nonwetting   (Air 1.205 kg/m^3)
muw  = 1.002e-3     # viscosity nonwetting (Water := 1.002e-3 kg/m s)
mun  = 1.81e-5      # viscosity wetting    (Air := 1.81e-5 kg/ m s)

gmag  = 9.8         # magnitude of gravity             (Earth 9.8 m/s^2) 
m_per_d_to_m_per_s = 1.1574074e-5

#water density parameters
beta_w = rhow*gmag*4.524e-10#rhow*gmag*4.524e-10#mwf should be rhow*gmag*4.524e-10
density_w_exponential = {'model':'Exponential',
                         'rho_0':rhow,
                         'psi_0':0.0,
                         'beta':beta_w}

#air density parameters
Temp = 288.15       # Temperature
R    = 8.314        # Ideal Gas Constant   (N/mol k)
M    = 2.896e-2     # Molar Mass           (Air := 2.896e-2)
p_o  = 1.013e5      # (N/m^2)
W    = 2.9e-4       # molar weight  [kg/mol]
m_to_Newton_per_m2 =rhow*gmag#mwf should be rhow*gmag

density_n_ideal = {'model':'IdealGas',
                   'T':Temp,
                   'W':W,
                   'R':R,
                   'headToPressure':m_to_Newton_per_m2,
                   'rho_0':rhon,
                   'psi_0':0.0}

a    = 1.0e-6

porousMediumDatabase = {}
for medium in ['sand','loam','clay','clay-2','clay-3']:
    porousMediumDatabase[medium] = {}
#
porousMediumDatabase['sand']['sw_min']    = 0.308
porousMediumDatabase['sand']['sw_max']    = 1.0
porousMediumDatabase['sand']['omega']     = 0.301
porousMediumDatabase['sand']['Ks']        = 5.040
porousMediumDatabase['sand']['mvg_alpha'] = 5.470
porousMediumDatabase['sand']['mvg_n']     = 4.264
porousMediumDatabase['sand']['bc_pd']     = 0.120
porousMediumDatabase['sand']['bc_lambda'] = 3.264
#
porousMediumDatabase['loam']['sw_min']    = 0.308
porousMediumDatabase['loam']['sw_max']    = 1.0
porousMediumDatabase['loam']['omega']     = 0.301
porousMediumDatabase['loam']['Ks']        = 5.040
porousMediumDatabase['loam']['mvg_alpha'] = 5.470
porousMediumDatabase['loam']['mvg_n']     = 4.264
porousMediumDatabase['loam']['bc_pd']     = 0.120
porousMediumDatabase['loam']['bc_lambda'] = 3.264
#
porousMediumDatabase['clay']['sw_min']    = 0.232
porousMediumDatabase['clay']['sw_max']    = 1.0
porousMediumDatabase['clay']['omega']     = 0.410
porousMediumDatabase['clay']['Ks']        = 0.062
porousMediumDatabase['clay']['mvg_alpha'] = 1.900
porousMediumDatabase['clay']['mvg_n']     = 1.310
porousMediumDatabase['clay']['bc_pd']     = 0.044
porousMediumDatabase['clay']['bc_lambda'] = 0.310
#
porousMediumDatabase['clay-2']['sw_min']    = 0.277
porousMediumDatabase['clay-2']['sw_max']    = 1.0
porousMediumDatabase['clay-2']['omega']     = 0.368
porousMediumDatabase['clay-2']['Ks']        = 7.970
porousMediumDatabase['clay-2']['mvg_alpha'] = 3.350
porousMediumDatabase['clay-2']['mvg_n']     = 2.00
porousMediumDatabase['clay-2']['bc_pd']     = 1.0/porousMediumDatabase['clay-2']['mvg_alpha']
porousMediumDatabase['clay-2']['bc_lambda'] = porousMediumDatabase['clay-2']['mvg_n'] - 1.0
#
porousMediumDatabase['clay-3']['sw_min']    = 0.277
porousMediumDatabase['clay-3']['sw_max']    = 1.0
porousMediumDatabase['clay-3']['omega']     = 0.268
porousMediumDatabase['clay-3']['Ks']        = 7.970e-2
porousMediumDatabase['clay-3']['mvg_alpha'] = 3.350
porousMediumDatabase['clay-3']['mvg_n']     = 2.00
porousMediumDatabase['clay-3']['bc_pd']     = 1.0/porousMediumDatabase['clay-2']['mvg_alpha']
porousMediumDatabase['clay-3']['bc_lambda'] = porousMediumDatabase['clay-2']['mvg_n'] - 1.0

for medium in porousMediumDatabase.keys():
    porousMediumDatabase[medium]['mvg_m']   = 1.0-1.0/porousMediumDatabase[medium]['mvg_n']
    # Hydraulic Conductivity  in meters per second
    porousMediumDatabase[medium]['Ksw']     = porousMediumDatabase[medium]['Ks']*m_per_d_to_m_per_s   

#now set default values based on PorousMedia key
assert porousMediumDatabase.has_key(PorousMedia), 'PorousMedia=%s not found in database keys= %s ' %(PorousMedia,
                                                                                                     porousMediumDatabase.keys())
sw_min    = porousMediumDatabase[PorousMedia]['sw_min']
sw_max    = porousMediumDatabase[PorousMedia]['sw_max']
omega     = porousMediumDatabase[PorousMedia]['omega']
Ks        = porousMediumDatabase[PorousMedia]['Ks']
mvg_alpha = porousMediumDatabase[PorousMedia]['mvg_alpha']
mvg_n     = porousMediumDatabase[PorousMedia]['mvg_n']
bc_pd     = porousMediumDatabase[PorousMedia]['bc_pd']
bc_lambda = porousMediumDatabase[PorousMedia]['bc_lambda']
mvg_m     = porousMediumDatabase[PorousMedia]['mvg_m']
Ksw       = porousMediumDatabase[PorousMedia]['Ksw']  
# if (PorousMedia == 'sand'):
#     sw_min    = 0.308
#     sw_max    = 1.0
#     omega     = 0.301
#     Ks        = 5.040
#     mvg_alpha = 5.470
#     mvg_n     = 4.264
#     bc_pd     = 0.120
#     bc_lambda = 3.264
# if (PorousMedia == 'loam'):
#     sw_min    = 0.181
#     sw_max    = 1.0
#     omega     = 0.430
#     Ks        = 0.250
#     mvg_alpha = 3.600
#     mvg_n     = 1.560
#     bc_pd     = 0.050
#     bc_lambda = 0.560
# if (PorousMedia == 'clay'):
#     sw_min    = 0.232
#     sw_max    = 1.0
#     omega     = 0.410
#     Ks        = 0.062
#     mvg_alpha = 1.900
#     mvg_n     = 1.310
#     bc_pd     = 0.044
#     bc_lambda = 0.310
# if (PorousMedia == 'clay-2'):
#     sw_min    = 0.277
#     sw_max    = 1.0
#     omega     = 0.368
#     Ks        = 7.970
#     mvg_alpha = 3.350
#     mvg_n     = 2.00
#     bc_pd     = 1.0/mvg_alpha
#     bc_lambda = mvg_n - 1.0

# mvg_m     = 1.0-1.0/mvg_n
# Ksw = Ks*m_per_d_to_m_per_s   # Hydrolic Conductivity  related to medium permeability  (3.5e-10 to 3.5e-1 for Water)


### need to add 

# grab a sand from Forsyth example
#     Type 3: k_i     = 4.898e-12 m^2
#             K_s     = 4.143
#             \theta_s= 0.3250
#             \theta_r= 0.3250*0.2643
#             \alpha  = 3.455     1/m
#             n       = 5.0
#     Type 2: k_i     = 5.55e-12 m^2
#             K_s     = 4.69     m/d
#             \theta_s= 0.3510
#             \theta_r= 0.3510*0.2806
#             \alpha  = 3.63     1/m
#             n       = 1.632

#helmig sanbox example 1 (ch 6)
#     Type 1: k_i     = 6.64e-11 m^2
#             K_s     = 7.0e-4 m/s
#             \theta_s= 0.40
#             \theta_r= 0.40*0.09
#             \lambda = 2.7
#             p_d     = 755 [Pa] -> psi_d= 0.0772 [m]
#             \alpha  = 1/p_d = 13.0 [1/m]
#             n       = lambda+1 = 3.7
#     Type 2: k_i     = 7.15e-12 m^2
#             K_s     = 7.0e-5 m/s
#             \theta_s= 0.39
#             \theta_r= 0.40*0.12
#             \lambda = 2.0
#             p_d     = 2060 [Pa] -> psi_d= 0.211 [m]
#             \alpha  = 1/p_d = 4.75
#             n       = lambda+1 = 3.

#for convenience , to get bcs straight, need to replace with call to psk class?
def seVGM(psic,alVG,nVG,mVG):
    if psic <= 0: return 1.0
    tmp1 = pow(alVG*psic,nVG)
    tmp2 = pow(1.+tmp1,-mVG)
    return min(max(tmp2,0.0),1.0)
def pcVGM(se,alVG,nVG,mVG):
    if se >= 1.: return 0.0
    tmp1 = pow(se,-1./mVG)
    tmp2 = pow(tmp1-1.0,1.0/nVG)/alVG
    return tmp2
def seBCB(psic,pdBC,lamBC):
    if psic <= pdBC: return 1.0
    tmp1 = pow(pdBC/psic,lamBC)
    return min(max(tmp1,0.0),1.0)
def pcBCB(se,pdBC,lamBC):
    if se >= 1.0: return 0.0
    tmp1 = pow(se,-1.0/lamBC)
    return pdBC*tmp1

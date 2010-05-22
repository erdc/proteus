# Set the model parameters to be given to the impes model for two phase flow. 

model = 'BCB'       # Choose from {simp,BCB,BCM,VGB,VGM}
PorousMedia = 'sand' 

L=(1.0,1.0,1.0)
nd=2

# Model Parameters
rhow = 997.0        # density wetting      (Water 997 kg/m^3)
rhon = 1.205        # density nonwetting   (Air 1.205 kg/m^3)
muw  = 1.002e-3     # viscosity nonwetting (Water := 1.002e-3 kg/m s)
mun  = 1.81e-5      # viscosity wetting    (Air := 1.81e-5 kg/ m s)

Temp = 288.15       # Temperature
R    = 8.314        # Ideal Gas Constant   (N/mol k)
M    = 2.896e-2     # Molar Mass           (Air := 2.896e-2)
p_o  = 1.013e5      # (N/m^2)

g    = [0.0,9.8] #-9.8]   # gravity              (Earth 9.8 m/s^2) 

FudgeFactor = 0.01
T = 100.0

if (PorousMedia == 'sand'):
    sw_min    = 0.308
    sw_max    = 1.0
    omega     = 0.301
    Ks        = 5.040
    mvg_alpha = 5.470
    mvg_n     = 4.264
    bc_pd     = 0.120
    bc_lambda = 3.264
if (PorousMedia == 'loam'):
    sw_min    = 0.181
    sw_max    = 1.0
    omega     = 0.430
    Ks        = 0.250
    mvg_alpha = 3.600
    mvg_n     = 1.560
    bc_pd     = 0.050
    bc_lambda = 0.560
if (PorousMedia == 'clay'):
    sw_min    = 0.232
    sw_max    = 1.0
    omega     = 0.410
    Ks        = 0.062
    mvg_alpha = 1.900
    mvg_n     = 1.310
    bc_pd     = 0.044
    bc_lambda = 0.310
    
Ksw = Ks*1.1574074e-5   # Hydrolic Conductivity  related to medium permeability  (3.5e-10 to 3.5e-1 for Water)
mvg_m = 1.0-1.0/mvg_n
Q_m_per_d = 0.1 
q    = Q_m_per_d*1.1574074e-5
q = 0.0

a    = 1.0e-6


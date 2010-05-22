from pyadh import *
from pyadh.default_p import *

model = 'BCB'                # Choose from {simp,BCB,BCM,VGB,VGM}
PorousMedia = 'sand' 

# Model Parameters

rhow = 997.0                 # density wetting      (Water 997 kg/m^3)
rhon = 1.205                 # density nonwetting   (Air 1.205 kg/m^3)
muw  = 1.002e-3              # viscosity nonwetting (Water := 1.002e-3 kg/m s)
mun  = 1.81e-5               # viscosity wetting    (Air := 1.81e-5 kg/ m s)
gravity = 9.8
g    = [0.0,-gravity,0.0]    # gravity              (Earth 9.8 m/s^2) 


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

thetaS        = 0.301   
thetaR        = 0.093   
m_per_s_by_m_per_d = 1.1574074e-5
mvg_m = 1.0-1.0/mvg_n
Ksw = Ks*m_per_s_by_m_per_d  # Hydrolic Conductivity  related to medium permeability  (3.5e-10 to 3.5e-1 for Water)

Q_m_per_d = 0.1 
q    = Q_m_per_d*1.1574074e-5
a    = 1.0e-6
lengthScale   = 1.0 
timeScale     = 1.0
pondingPressure=-0.01 #0.1

FudgeFactor = 1e-1
useHet = True#False
slit_is_top = False#True
open_bottom = True

T = 100.0
DT=1.0e-1
nDTout = int(T/DT)
nLevels=4

polyfile = "kueper"
nd=2
top = 0.5
right = 0.7

# Set the model parameters to be given to the impes model for two phase flow. 

model = 'simp'#'VGM'       # Choose from {simp,BCB,BCM,VGB,VGM}

capillaryDiffusionScaling = 0.0 # to turn off capillary pressure and test BL

nd = 1
L = (1.0,1.0,1.0)
# Model Parameters
rhow = 1.0
rhon = rhow
muw  = 1.00    
mun  = 2.00    
omega= 1.0
g    = [0.0]        # gravity              (Earth 9.8 m/s^2) 
q    = 1.0

Se_top = 1.0
Se_bottom = 0.0
T = 0.5

Ksw = 1.0
sw_min = 0.0
sw_max = 1.0

#go ahead and set for keeping uniform quadrature orders?
darcy_simp_quad_order = 5


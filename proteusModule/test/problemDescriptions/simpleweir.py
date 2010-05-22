nd = 2
DT = 1.0e-2
T=9999*DT
nDTout = int(T/DT)
nLevels = 5

polyfile = "simpleWeir"

epsFact = 3.0
#epsFact = 0.0
#water
rho_0=998.2
nu_0=1.004e-6
#air
rho_1=1.205
nu_1=1.500e-5
#gravity
g=[0.0,-9.8]
inflow = 1.0
waterLevel = 0.5
weirHeight = 0.4
weirStart = 9.9
weirEnd = 10.1
flumeEnd = 20.0
flumeTop = 2.0

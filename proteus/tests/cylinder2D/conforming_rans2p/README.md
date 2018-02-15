## Introduction
This folder is to show how to use RANS2P2D module to solve NSE on the flow over the cylinder

## TODO

+ P2 does not work. check
+ 


## Note

+ RANS2P_1.py has force term implemented
+ RANS2P_2.h prints yyForce line
+ Use 
    dt_system_fixed = dt_fixed 
and 
    systemStepExact=False;
+ Use 
    triangleOptions="pAq30.0Dena"

## Data
.h5 is obtained for the case of p1, dt fixed: 0.005+0.0025, DX=0.04, nPoints_cyl=10.
# How to run it?
+ python setup.py build_ext -i
+ parun -v -l 5 poisson_so.py -b L2_batch.py


############################################################################################################
06/27/2018
# he is not right in computing error since the width and length are 2 for sbm
# add square domain with unstructured mesh

############################################################################################################
06/26/2018
# add disk domain with unstructured mesh

############################################################################################################

# Error of p1
refine = 5
t= 1; error_u_L2[0][0]= 0.0352341;
[       2]
t= 1; error_u_H1[0][0]= 0.267208;

refine = 6
t= 1; error_u_L2[0][0]= 0.0144469;
[       3]
t= 1; error_u_H1[0][0]= 0.178115;

refine=7
t= 1; error_u_L2[0][0]= 0.00758924;
[      12]
t= 1; error_u_H1[0][0]= 0.132393;

refine = 8
t= 1; error_u_L2[0][0]= 0.00338397;
[      32]
t= 1; error_u_H1[0][0]= 0.0901307;

# Error of p2
refine=5
t= 1; error_u_L2[0][0]= 0.0327477;
[       3]
t= 1; error_u_H1[0][0]= 0.282337;

refine=6
t= 1; error_u_L2[0][0]= 0.012275;
[       5]
t= 1; error_u_H1[0][0]= 0.17666;

refine=7
t= 1; error_u_L2[0][0]= 0.00687616;
[      18]
t= 1; error_u_H1[0][0]= 0.134029;

refine=8
t= 1; error_u_L2[0][0]= 0.00307069;
[      75]
t= 1; error_u_H1[0][0]= 0.0903693;


## conforming mesh
# Error of P1
refine = 4
t= 1; error_u_L2[0][0]= 0.00242525;
[       0]
t= 1; error_u_H1[0][0]= 0.0442464;

refine = 5
t= 1; error_u_L2[0][0]= 0.000547616;
[       1]
t= 1; error_u_H1[0][0]= 0.0221354;
refine = 6
t= 1; error_u_L2[0][0]= 0.000136134;
[       1]
t= 1; error_u_H1[0][0]= 0.0111423;
refine = 7
t= 1; error_u_L2[0][0]= 3.55815e-05;
[       5]
t= 1; error_u_H1[0][0]= 0.00555925;

# Error of P2: it gives so small error since the true solution is polynomial of order 2.
refine = 5
t= 1; error_u_L2[0][0]= 3.3255e-08;
[      17]
t= 1; error_u_H1[0][0]= 1.33891e-06;
refine = 6
t= 1; error_u_L2[0][0]= 3.03935e-08;
[       9]
t= 1; error_u_H1[0][0]= 1.95269e-06;
refine = 7
t= 1; error_u_L2[0][0]= 3.38378e-08;
[      21]
t= 1; error_u_H1[0][0]= 2.83192e-06;

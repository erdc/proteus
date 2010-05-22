We solve heterogeneous "Laplace" equation in 2D on three different meshes and 
with two different "extreme" examples:

1) Ex1 - Obstacles_D1: Square [0,1]x[0,1] domain with Dirichelet boundary 
   conditions 1 and 0 at x=0.0 and x=1.0;
2) Ex1.1 - Obstacles_D2: Refinement of example 1.
3) Ex2 - Clusters_D1: Square domain with low permeability clusters.

The permeability ratio between the two zones is always 1:10^{-6}.

The results are shown graphically in terms of streamlines.
Streamlines are calculated from the constant elemental velocity field
by starting from 100-2 (98) uniformly distributed points on the inflow 
boundary and drawing a line parallel to the velocity vector starting 
from the entering point of each element.

For both LN (Larson-Niklasson) and MH (RT0-P0 with hybridization) 
schemes, the elemental velocity fields are interpolated using the 
RT0 basis functions.

The following files can be found within each of the above folders:

{File Name}_LN.pdf - Streamlines obtained from the LN post-processing technique.
{File Name}_MH.pdf - Streamlines obtained from the MH technique.
{File Name}_LN_AllData.txt - All the relevant input and output 
data used in the LN scheme. In particular, we include:

1 - Node coordinates.
2 - Element incidences and properties.
3 - Incidences and reference edge vectors.
4 - Material property definition.
5 - Dirichelet boundary conditions.
6 - P1 Galerkin nodal potential values.
7 - P1 Galerkin element fluxes.
8 - LN-Corrected integrated fluxes across mesh edges.
9 - LN-Corrected element fluxes.


As you can see the main difference between LN and MH streamlines is 
close to the corners of the low permeability areas, where in the case of LN
some streamlines enter the low permeability zones. From the two 
refinements in Example 1, you can see that the problem persists, even 
though it is obviously constrained in single elements, and thus the 
convergence should not be affected. But local conservation is 
questionable!!



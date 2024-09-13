##################
# *** 2D HEX *** #
##################
echo "          ****************************"
echo "          ********** 2D HEX **********"
echo "          ****************************"
parun -l1 -v poisson_p.py poisson_n.py -b L2_batch.py -C "nd=2 refinement=0 useHex=True"
parun -l1 -v poisson_p.py poisson_n.py -b L2_batch.py -C "nd=2 refinement=1 useHex=True"

######################
# *** 2D SIMPLEX *** #
######################
echo "          ********************************"
echo "          ********** 2D SIMPLEX **********"
echo "          *********************************"
parun -l1 -v poisson_p.py poisson_n.py -b L2_batch.py -C "nd=2 refinement=0 useHex=False genMesh=True"
parun -l1 -v poisson_p.py poisson_n.py -b L2_batch.py -C "nd=2 refinement=1 useHex=False genMesh=True"

##################
# *** 3D HEX *** #
##################
echo "          ****************************"
echo "          ********** 3D HEX **********"
echo "          ****************************"
parun -l1 -v poisson_p.py poisson_n.py -b L2_batch.py -C "nd=3 refinement=0 useHex=True"
parun -l1 -v poisson_p.py poisson_n.py -b L2_batch.py -C "nd=3 refinement=1 useHex=True"

######################
# *** 2D SIMPLEX *** #
######################
echo "          ********************************"
echo "          ********** 3D SIMPLEX **********"
echo "          ********************************"
parun -l1 -v poisson_p.py poisson_n.py -b L2_batch.py -C "nd=3 refinement=0 useHex=False genMesh=True"
parun -l1 -v poisson_p.py poisson_n.py -b L2_batch.py -C "nd=3 refinement=1 useHex=False genMesh=True"

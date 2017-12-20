import os

hk = [1.0, 0.5, 0.3]
ref= [3, 4, 5]

for i in range(3):
    for j in range(3):
        print(' '.join(['parun -v -l 5', 'ls_vortex_2d_so.py', '-C', " 'hk=%0.1f lRefinement=%d applyCorrection=True applyRedistancing=True'" %(hk[i],ref[j])]))
        os.system(' '.join(['parun -v -l 5', 'ls_vortex_2d_so.py', '-C', " 'hk=%0.1f lRefinement=%d applyCorrection=True applyRedistancing=True'" %(hk[i],ref[j])]))
        print(' '.join(['parun -v -l 5', 'ls_vortex_2d_so.py', '-C', " 'hk=%0.1f lRefinement=%d applyCorrection=True applyRedistancing=False'" %(hk[i],ref[j])]))
        os.system(' '.join(['parun -v -l 5', 'ls_vortex_2d_so.py', '-C', " 'hk=%0.1f lRefinement=%d applyCorrection=True applyRedistancing=False'" %(hk[i],ref[j])]))
        print(' '.join(['parun -v -l 5', 'ls_vortex_2d_so.py', '-C', " 'hk=%0.1f lRefinement=%d applyCorrection=False applyRedistancing=False'" %(hk[i],ref[j])]))
        os.system(' '.join(['parun -v -l 5', 'ls_vortex_2d_so.py', '-C', " 'hk=%0.1f lRefinement=%d applyCorrection=False applyRedistancing=False'" %(hk[i],ref[j])]))
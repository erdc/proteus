#this is python, so you can loop over the particles
# for i in range(nParticles): ...["q:('phis',{i})".format(i)]
#import pdb
#pdb.set_trace()
#note: this has to be attached to correct model
#navier-stokes model has index 4 in this model
#from so file:
#    pnList = [("vos_p",               "vos_n"),#0
              # ("vof_p",               "vof_n"),#1
              # ("ls_p",                "ls_n"),#2
              # ("redist_p",            "redist_n"),#3
              # ("ls_consrv_p",         "ls_consrv_n"),#4
              # ("threep_navier_stokes_sed_p", "threep_navier_stokes_sed_n"),#5
              # ("twp_navier_stokes_p", "twp_navier_stokes_n"),#6
              # ("pressureincrement_p", "pressureincrement_n"),#7
              # ("pressure_p", "pressure_n"),#8
              # ("pressureInitial_p", "pressureInitial_n")]#9

simFlagsList[4]['storeQuantities'] = ["q:('phis',0)"]
start
end
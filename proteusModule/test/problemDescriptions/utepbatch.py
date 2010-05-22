#if you want to view the solution
#if simFlags['plotTimes'] != None:
#    #specify quantities that you would like to look at
#    #q: means visualize at element quadrature points
#    #*,0 means component 0 of system
#    if simFlags.has_key(('plotQuantities')): 
#        simFlags['plotQuantities'].append("q:('u',0)") #solution
#    else:
#        simFlags['plotQuantities'] = ["q:('u',0)"]
#    simFlags['plotQuantities'].append("q:('velocity',0)")#gw velocity field
#    simFlags['plotQuantities'].append("q:('r',0)")#well terms
#    simFlags['plotQuantities'].append("q:('visPerm',0)")#permeability field
#
#where to store the reults
simFlags['dataDir']='./results'
#file name in which to store results
simFlags['dataFile']='utep_example.dat'
#when to store results
simFlags['storeTimes'] = ['Last']
#what to store, u_dof means solution degrees of freeom, can also store quadrature point info and
#mesh data structure
simFlags['storeQuantities'] = ['u_dof']#['mesh','u_dof']

start
quit

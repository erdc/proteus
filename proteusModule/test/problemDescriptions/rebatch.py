#error calculation section
simFlags['errorQuantities'] = ['u','velocity'] #compute error solution
simFlags['errorTimes'] = ['Last'] #compute error solution
simFlags['errorTypes'] = ['localMassBalance','globalMassBalance'] #compute error in local mass bal
simFlags['errorTypes'].append('numericalSolution')
simFlags['errorNorms'] = ['L2','LI','TV']


simFlags['echo'] = True

#visualization section
simFlags['plotTimes'] = ['Last']
if simFlags['plotTimes'] != None:
    if opts.viewer == 'gnuplot':
        simFlags['plotOptions']    = ['setGnuplotGridSize'] #manually set gnuplot grids
    if opts.viewer == 'matlab':
        simFlags['plotOptions']    = ['usePDEtoolbox']
    if 'plotQuantities' not in simFlags:
        simFlags['plotQuantities'] = []
    simFlags['plotQuantities'].append('velocity')
    simFlags['plotQuantities'].append("q:('m',0)")
    #simFlags['plotQuantities'].append("q:('phi',0)")
    #simFlags['plotQuantities'].append("ebq_global:('velocity',0)")
    #simFlags['plotQuantities'].append("q:('mt',0)")
#leave blank for now

#disk storage section
simFlags['dataFile']       = simFlags['simulationName']+'_results.dat'
simFlags['dataDir']        = os.getcwd()+'/results'
simFlags['storeTimes']     = ['Last']        #write to file at last step
simFlags['storeQuantities']= ['simulationData','errorData'] #include errorData for mass bal
simFlags['storeQuantities'].append('u_dof') #save soln component 0 at element quadrature points
simFlags['storeQuantities'].append("q:('u',0)") #save soln component 0 at element quadrature points
simFlags['storeQuantities'].append("q:('m',0)") #save mass component 0 at element quadrature points
simFlags['storeQuantities'].append("q:('velocity',0)") #save velocity component 0 at element quadrature points
start
quit

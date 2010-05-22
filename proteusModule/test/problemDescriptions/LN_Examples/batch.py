#
for simFlags in simFlagsList[0:1]:
    if not simFlags.has_key('storeQuantities'):
        simFlags['storeQuantities']= []
    simFlags['storeQuantities'].append("q:('velocity',0)")
    #simFlags['storeQuantities'].append("q:('a',0,0)")
    for key in ['errorQuantities','errorTypes','errorNorms','errorTimes']:
        if not simFlags.has_key(key):
            simFlags[key] = []
    simFlags['errorQuantities']=['u','velocity'] #compute error in solution value
    simFlags['errorTypes']= ['numericalSolution'] #compute error in soln and glob. mass bal
    simFlags['errorTypes'].append('localMassBalance') #compute local mass balance errors
    simFlags['errorNorms']= ['L2'] #compute L2 norm in space
    simFlags['errorTimes']= ['Last']
    simFlags['echo'] = True
    #
    simFlags['dataFile']       = simFlags['simulationName']+'_results.dat'
    simFlags['dataDir']        = os.getcwd()+'/results'
    simFlags['storeTimes']     = ['Last']        #write to file at last step
    simFlags['storeQuantities'].append('errorData')
    simFlags['storeQuantities'].append('simulationData')#save error calculations and basic mesh info 

start
quit

for simFlags in simFlagsList:
    simFlags['errorQuantities']=['u','velocity']
    simFlags['errorTypes']= ['numericalSolution'] #compute error in soln and glob. mass bal
    simFlags['errorNorms']= ['L2','L1','LI','H1'] #compute L2 norm in space
    simFlags['errorTimes']= ['Last'] #All, Last
    simFlags['echo']=True
    simFlags['dataFile']       = simFlags['simulationName']+'_results'
    simFlags['dataDir']        = os.getcwd()+'/results'
    simFlags['storeQuantities']= ['simulationData','errorData'] #include errorData for mass bal
    simFlags['storeTimes']     = ['Last']
#
start
quit
simFlags['errorQuantities']=['u']
simFlags['errorTypes']= ['numericalSolution'] #compute error in soln and glob. mass bal
#simFlags['errorNorms']= ['L2'] #compute L2 norm in space
simFlags['errorNorms']= ['L2','L2_L2'] # compute L2_L2 norm in space
simFlags['errorTimes']= ['All']   # not just last
simFlags['echo']=True
simFlags['dataFile']       = simFlags['simulationName']+'_results.dat'
simFlags['dataDir']        = os.getcwd()+'/results'
simFlags['storeQuantities']= ['simulationData','errorData'] #include errorData for mass bal
simFlags['storeTimes']     = ['Last']
simFlags['plotQuantities'].append('u_exact')
#simFlags['plotQuantities'].append('velocity')
start
quit

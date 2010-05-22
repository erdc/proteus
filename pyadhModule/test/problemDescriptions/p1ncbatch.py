#error calculation section
simFlags['errorQuantities']=['u','velocity'] #compute error in solution value
        
simFlags['errorTypes']= ['numericalSolution'] #compute error in soln and glob. mass bal
simFlags['errorTypes'].append('localMassBalance') #compute local mass balance errors
simFlags['errorNorms']= ['L2'] #compute L2 norm in space
simFlags['errorTimes']= ['Last']

#visualization section
if simFlags.has_key('plotTimes') and simFlags['plotTimes'] != None:
    simFlags['plotQuantities'].append('velocity')
    simFlags['plotQuantities'].append('u_exact')
    simFlags['plotQuantities'].append('velocity_exact')
    #and example of how to look at source term simFlags['plotQuantities'].append("q:('r',1)")
#leave blank for now

#disk storage section
simFlags['dataFile']       = simFlags['simulationName']+'_results.dat'
simFlags['dataDir']        = os.getcwd()+'/results'
simFlags['storeTimes']     = ['Last']        #write to file at last step
simFlags['storeQuantities']= ['errorData','simulationData']#save error calculations and basic mesh info 
simFlags['storeQuantities'].append('u_dof') #save dofs as well
simFlags['storeQuantities'].append("q:('u',0)") #save soln component 0 at element quadrature points
simFlags['storeQuantities'].append("q:('u',1)") #save soln component 1 at element quadrature points
simFlags['storeQuantities'].append("q:('velocity',0)")#velocity component 0 at element quad points
simFlags['storeQuantities'].append("q:('velocity',1)")#velocity component 1 at element quad points 
simFlags['storeQuantities'].append("ebq:('velocity',0)")#velocity component 0 at element boundary quad points
simFlags['storeQuantities'].append("ebq:('velocity',1)")#velocity component 1 at element boundary quad points
simFlags['storeQuantities'].append("ebq_global:('velocity',0)")#velocity component 0 at element boundary quad points
simFlags['storeQuantities'].append("ebq_global:('velocity',1)")#velocity component 1 at element boundary quad points
simFlags['storeQuantities'].append("ebqe:('velocity',0)")#velocity component 0 at element boundary quad points
simFlags['storeQuantities'].append("ebqe:('velocity',1)")#velocity component 1 at element boundary quad points
#simFlags['storeQuantities'].append("q:('velocity_dofs',0)") #save velocity dofs for component 0
#simFlags['storeQuantities'].append("q:('velocity_dofs',1)") #save velocity dofs for component 1
#how to run section
simFlags['echo']=True

start
quit
#refine grid
n.nLevels += 1
simFlags['dataFile']      = simFlags['simulationName']+'_1_results.dat'

simFlags['errorNorms'].append('L1') #add L1 
#simFlags['plotQuantities'] = [None]
#simFlags['echo']=True
#simFlags['storeQuantities']= ['errorData','simulationData'] #reset so don't save dofs

start
quit
#refine grid
n.nLevels += 1
simFlags['dataFile']      = simFlags['simulationName']+'_2_results.dat'

start

quit 

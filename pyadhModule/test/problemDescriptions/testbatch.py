#error calculation section
simFlags['errorQuantities']=['u','velocity'] #compute error in solution value
        
simFlags['errorTypes']= ['numericalSolution','localMassBalance'] #compute error in soln and glob. mass bal
simFlags['errorNorms']= ['L2']#,'L1','TV','W11semi'] #compute which norms in space
simFlags['errorTimes']= ['All']
#whether or not to compute space-time norms
if n.timeIntegration != TimeIntegration.NoIntegration and 'L2' in simFlags['errorNorms']:
    simFlags['errorNorms'].append('L2_LI') 
    simFlags['errorNorms'].append('L2_L2')
if n.timeIntegration != TimeIntegration.NoIntegration and 'L1' in simFlags['errorNorms']:
    simFlags['errorNorms'].append('L1_LI')
    simFlags['errorNorms'].append('L1_L1')

#visualization section
#only plot final solution?
simFlags['plotTimes'] = ['Init','Last']
if opts.viewer == 'gnuplot' and simFlags['plotTimes'] != None:
    simFlags['plotOptions']    = ['setGnuplotGridSize'] #manually set gnuplot grids
    simFlags['plotQuantities'].append('u_exact')
    simFlags['plotQuantities'].append('velocity')
    #simFlags['plotQuantities'].append("q:('m',0)")
    simFlags['plotQuantities'].append("ebq_global:('velocity',0)")
    simFlags['plotQuantities'].append("ebq:('velocity',0)")
    simFlags['plotQuantities'].append("ebq_global:('velocityAverage',0)")
    
#leave blank for now

#disk storage section
simFlags['dataFile']       = simFlags['simulationName']+'_results.dat'
simFlags['dataDir']        = os.getcwd()+'/results'
simFlags['storeTimes']     = ['Last']        #write to file at last step
#simFlags['storeTimes']     = ['All']        #write to file every output step
simFlags['storeQuantities']= ['errorData','simulationData']#save error calculations and basic mesh info 
#simFlags['storeQuantities'].append('u_dof') #save dofs as well
#simFlags['storeQuantities'].append("q:('u',0)") #save soln component 0 at element quadrature points
#simFlags['storeQuantities'].append('mesh') #save soln component 0 at element quadrature points
#how to run section
simFlags['echo']=True

start
quit

#refine grid
n.nLevels += 1
simFlags['dataFile']      = simFlags['simulationName']+'_1_results.dat'

#simFlags['plotQuantities'] = [None]
#simFlags['echo']=True
#simFlags['storeQuantities']= ['errorData','simulationData'] #reset so don't save dofs

start
quit 

#refine again
n.nLevels += 1
simFlags['dataFile']      = simFlags['simulationName']+'_2_results.dat'

start


#refine again
n.nLevels += 1
simFlags['dataFile']      = simFlags['simulationName']+'_3_results.dat'

start

quit


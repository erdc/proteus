simFlags['errorQuantities']=['u'] #compute error in solution value
simFlags['errorTypes']= ['numericalSolution'] #compute error in soln and glob. mass bal
simFlags['errorTimes']= ['Last']
simFlags['errorNorms']= ['L2'] #compute L2 norm in space
 
simFlags['echo']=True

start
quit

from proteus import Context
from proteus import Comm
comm = Comm.get()
ctx = Context.get()


# simulation flags for error analysis
#
# simFlagsList is initialized in proteus.iproteus
#

simFlagsList[0]['errorQuantities']=['u']
simFlagsList[0]['errorTypes']= ['numericalSolution'] #compute error in soln and glob. mass bal
simFlagsList[0]['errorNorms']= ['L2'] #compute L2 norm in space or H0 or ...
simFlagsList[0]['errorTimes']= ['Last'] #'All', 'Last'
simFlagsList[0]['echo']=True
#simFlagsList[0]['echoRelativeErrors']=True

#
start
quit
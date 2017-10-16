
from proteus import Context
from proteus import Comm
comm = Comm.get()
ctx = Context.get()


# simulation flags for error analysis
#
# simFlagsList is initialized in proteus.iproteus
#

simFlagsList[0]['errorQuantities'] = ['u']
# compute error in soln and glob. mass bal
simFlagsList[0]['errorTypes'] = ['numericalSolution']
# compute L2 norm in space or H0 or ...
simFlagsList[0]['errorNorms'] = ['L1', 'L2']
simFlagsList[0]['errorTimes'] = ['Last']  # 'All', 'Last'
simFlagsList[0]['echo'] = True
simFlagsList[0]['storeTimes'] = []
simFlagsList[0]['storeQuantities'] = ['mesh', 'errorData']
simFlagsList[0]['dataDir'] = '.'
simFlagsList[0]['dataFile'] = ctx.datafile
#
start
quit

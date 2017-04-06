from proteus import Domain
import proteus.MeshTools
from proteus.MeshAdaptPUMI import MeshAdaptPUMI

nd = 2 #number of dimensions in the problem
parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.element #type of partitioning if parallel
nLayersOfOverlapForParallel = 0 #amount of ghosting if parallel
boundaries=['left','right','bottom','top'] #boundary tag dictionary
boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])

domain = Domain.PUMIDomain(dim=nd) #initialize the domain
domain.faceList=[[11],[13],[14],[12]] #model entities associated wtih boundary tags
adaptMesh = True #adapt the mesh?
adaptMesh_nSteps = 5 #amount of time steps before checking error?

hMax = 0.08
hMin = 0.00625
adaptMesh_numIter = 2 #number of iterations for mesh adaptation
errorType="ERM" #only just ERM at the moment
logSwitch="off" #on or off
target_error = 10.0 
target_element_count = 8000

domain.PUMIMesh=MeshAdaptPUMI.MeshAdaptPUMI(hmax=hMax, 
                                            hmin=hMin, 
                                            numIter=adaptMesh_numIter,
                                            sfConfig=errorType,
                                            logType=logSwitch,
                                            targetError=target_error,
                                            targetElementCount=target_element_count)

domain.PUMIMesh.loadModelAndMesh("Dambreak.null","Dambreak.msh")
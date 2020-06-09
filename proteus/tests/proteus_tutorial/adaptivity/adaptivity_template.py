#initialize PUMI domain
domain = Domain.PUMIDomain(dim=2) #initialize the domain

#read the geometry and mesh
#domain.AdaptManager.PUMIAdapter.loadModelAndMesh(b"Reconstructed.dmg", b"Reconstructed.smb")
domain.AdaptManager.PUMIAdapter.loadModelAndMesh(b"Reconstructed.dmg", b"4-Proc/.smb")

#mesh is partitioned by elements and not by nodes
from proteus import MeshTools
parallelPartitioningType = MeshTools.MeshParallelPartitioningTypes.element
domain.MeshOptions.setParallelPartitioningType('element')

#input for adaptation describing the staggered solve, needed for proper post-adapt progression
domain.AdaptManager.modelDict = {'flow':0,'phase':2,'correction':[3,4]}

#input to determine what strategy for adaptation is used
domain.AdaptManager.sizeInputs = [b'interface',b'error_vms']

#do you want to adapt?
domain.AdaptManager.adapt = 1

#max,min mesh edge length
domain.AdaptManager.hmax = he*2.0
domain.AdaptManager.hmin= he/2.0

#mesh edge length for near interface
domain.AdaptManager.hphi= he/2.0

#number of predictive timesteps for interface-based adapt
domain.AdaptManager.numAdaptSteps= 10

#target element error for error estimation-based adapt
domain.AdaptManager.targetError= 2.0

#maximum ratio between adjacent edge lengths
domain.AdaptManager.gradingFactor= 1.5

#number of iterations in adapt routine
domain.AdaptManager.numIterations= 5

#see before and after meshes and other info for adaptation
domain.AdaptManager.logging= 0

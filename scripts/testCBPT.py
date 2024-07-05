#! /usr/bin/env python
from proteus import *
import numpy
def generateTrajectories(mesh,v,x_source,t_0,t_f):
    t_traj = [0.0]; x_traj = [x_source]; e_traj = [0]
    t = t_0
    while t < t_f:
        x_0 = x_traj[-1]; t_0 = t_traj[-1]; e_0 = e_traj[-1]
        x_1 = mesh.nodeArray[mesh.elementNodesArray[e_0,1],0]
        x_12 = 0.5*(x_0+x_1)
        t_12= t_0 + (x_12-x_0)/v #include midpoint to test search alg
        t_1 = t_0 + (x_1-x_0)/v
        e_12= e_0
        e_1 = mesh.elementNeighborsArray[e_0,0]
        if e_1 == -1: e_1 = e_0
        t_traj.append(t_12); x_traj.append(x_12); e_traj.append(e_12)
        t_traj.append(t_1); x_traj.append(x_1); e_traj.append(e_1)
        t = t_1
    #
    t_traj = numpy.array(t_traj);  e_traj = numpy.array(e_traj,dtype='i')
    n_traj = t_traj.shape[0]
    x_tmp = numpy.zeros((n_traj,3),'d')
    x_tmp[:,0] = x_traj[:]; x_traj=x_tmp
    return x_traj,t_traj,e_traj
def generateMassSource(t_0,t_f,timeflag=0):
    if timeflag == 1:
        #mass source = 1 over interval intermediate interval
        massSource_t = numpy.array([t_0-1,t_0,0.1*(t_0+t_f),0.3*(t_0+t_f),0.4*(t_0+t_f),t_f,t_f+1])
        massSource_m = numpy.array([0.0,0.0,1.0,1.0,0.0,0.0,0.0])
    else:
        #constant source value is one
        massSource_t = numpy.array([t_0-1,t_0,0.5*(t_0+t_f),t_f,t_f+1])
        massSource_m = numpy.array([1.0,1.0,1.0,1.0,1.0])

    return massSource_t,massSource_m

def generatePhysicalCoefficients(mesh,nflags,testflag=0):
    decay = numpy.zeros((mesh.nElements_global,nflags),'d'); retardation = numpy.ones((mesh.nElements_global,nflags),'d')
    if testflag == 1:
        #no decay, constant retardation
        decay.fill(0.0)
        retardation.fill(1.5)
    elif testflag == 2:
        #decay, no retardation
        decay.fill(-0.1)
        retardation.fill(1.)
    else:
        pass
    return decay,retardation

def test_onesource(nx=11,Lx=1.0,ndt=20,q=1.0,theta=1.0,timeflag=0,testflag=0,plotSolution=True,wait4plot=True):
    """
    simple test example for CBPT in 1d, solves
    (\theta R c)_t + (q c)_x = -\theta \lambda c
    """
    #build the mesh on [0,Lx] with nx nodes
    mesh = MeshTools.EdgeMesh()
    mesh.generateEdgeMeshFromRectangularGrid(nx,Lx)

    #assume uniform porosity and velocity to make tracking simple
    v = q/theta;
    #source located at inflow, track over [t_0,t_f] with one particle
    x_source = 0.0; t_0 = 0.0; t_f=1.0; n_p = 1;

    #locations, times, and elements for trajectories
    x_traj,t_traj,e_traj = generateTrajectories(mesh,v,x_source,t_0,t_f)
    massSource_t,massSource_m = generateMassSource(t_0,t_f,timeflag=timeflag)

    n_traj = t_traj.shape[0]
    traj_offsets = numpy.array([0,n_traj],dtype='i')

    #only one type of particle
    nflags = 1
    particleFlags = numpy.zeros((nflags,),'i')
    #\lambda, R = (1 + \rho_b K_d/\theta)
    decay,retardation = generatePhysicalCoefficients(mesh,nflags,testflag=testflag)

    ct = numpy.zeros((mesh.nElements_global,),'d') #total concentration, c_t = R c_a

    from proteus.cellam import accumulateSourceContribution, integratePiecewiseLinearMassSource

    if plotSolution:
        from matplotlib import pylab
    #try to approximate for time tau
    taulist = numpy.linspace(t_0,t_f,ndt+1)
    xc = mesh.elementBarycentersArray[:,0]
    elementVolumeArray = numpy.array([abs(mesh.nodeArray[mesh.elementNodesArray[eN,1],0]-mesh.nodeArray[mesh.elementNodesArray[eN,0],0]) for eN in range(mesh.nElements_global)])
    #normalize by aqueous phase volume
    aqueousVolume = theta*elementVolumeArray

    for tau in taulist:
        ct.fill(0.0)
        accumulateSourceContribution(tau,traj_offsets,x_traj,t_traj,e_traj,massSource_t,massSource_m,decay,retardation,particleFlags,ct)
        ct /= aqueousVolume
        caq = ct/retardation[:,0]
        if plotSolution:
            ax = pylab.plot(xc,caq,'r*-',xc,ct,'bo-',hold=False);  pylab.legend(('c_aq','c_t')); pylab.title('tau= %s ' % tau)
            if wait4plot:
                input('hit return to continue')
        else:
            print('tau= %s c= %s ' % (tau,c))
    #
#

test_onesource()

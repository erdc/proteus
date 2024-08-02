"""
Collect higher level tools for running simulation, processing results, etc

.. inheritance-diagram:: proteus.SimTools
   :parts: 1
"""
from . import Norms
import numpy
import numpy as np
from . import FemTools
from .Profiling import logEvent

from proteus import Comm
comm = Comm.get()

#dummy classes for computing exact solution norms with error functions
class zeroFunction(object):
    def uOfX(self,x):
        return 0.0
    def uOfXT(self,x,t):
        return 0.0
class zeroVectorFunction(object):
    def __init__(self,shape):
        self.shape=shape
    def uOfX(self,x):
        return numpy.zeros(self.shape,'d')
    def uOfXT(self,x,t):
        return numpy.zeros(self.shape,'d')

class SimulationProcessor(object):
    """
    Collect some functionality for doing something with simulation results like
    calculating error, saving it to disk, etc.

    TO DO

    be able to append data to file correctly

    """
    import math
    defaultFlags = {'simulationName':None,       #label for simulation
                    'simulationNameProc':None,   #label for simulation and processor id
                    'dataFile':'results.dat',    #file for holding results
                    'dataDir' :'.',              #where file is located
                    'appendResults':False,       #append to existing data files?
                    'echo':False,                #print to screen
                    'echoRelativeErrors':False,  #print to screen also relative errors
                    'components':[0],            #list of components to monitor
                    'errorQuantities':[None],    #quantities in which to estimate error
                    'errorNorms':[None],         #norms to use for error calc
                    'errorTimes':[None],         #time levels at which to estimate error (All,Last, ...),
                    'errorTypes':[None],         #types of error to calculate (u-u_h,mass, etc)
                    'plotQuantities':[None],     #quantities to plot
                    'plotTimes':[None],          #time levels to plot (All,Last,tList ...
                    'plotOptions':{'ensight':{'on':False},
                                   'gnuplot':{'on':False},
                                   'matlab':{'on':False},
                                   'vtk':{'on':False}},         #options controlling way to plot
                    'storeQuantities':[None],    #quantities to write to disk
                    'storeTimes':[None]}         #time levels to store times (All, Last, tList ...)
    notSetVal = -12345.0
    #possible errors to measure
    ErrorQuantities = ['m','u','grad(u)','velocity']
    #possible error norms
    ErrorNorms      = ['L1','L2','LI','TV','H1','H1semi','W11','W11semi',
                       'L1_L1','L1_LI','L2_L2','L2_LI','TV_L1','TV_LI']
    #
    ErrorTimes      = ['All','Last','tList']
    #type of error to check (e.g., mass balance,
    ErrorTypes      = ['globalMassBalance','globalHeavisideMassBalance',
                       'localMassBalance',
                       'numericalSolution']
    #
    PlotTimes       = ['All','Last','Init','tList']
    #
    PlotQuantities  = ['u','u_exact','velocity','velocity_exact',"q:('%s',%d)","q:('%s',%d,%d)"]
    #
    PlotOptions     = {'ensight':{'on':False},
                       'gnuplot':{'on':False,'setGnuplotGridSize':True},
                       'matlab':{'on':False,'usePDEtoolbox':False},
                       'vtk':{'on':False,'pause':False,'hardcopy':False}}
    #
    StoreTimes      = ['All','Last','Init','tList']
    #
    StoreQuantities  = ['u','u_dof','errorData','simulationData','mesh',"q:('%s',%d)",'multilevelModel']
    #mwf temporary step before setting up fill archiving
    #data members that need to be stored to reconstruct a solution:
    SolutionDataStorageData = ['nSpace_global','name','dof','dofMap.l2g','dim_dof',
                               'isVector',
                               'CGDOFMap.lagrangesNodesArray','CGDOFMap.l2g']
    def __init__(self,flags=None,nLevels=1,pFile=None,nFile=None,
                 analyticalSolution={}):
        """
        save flags for what to compute, how to compute it,
        and labels for storing things, etc
        """
        self.analyticalSolution = {}
        if analyticalSolution is not None:
            self.analyticalSolution = analyticalSolution
        self.timeValues = []
        self.plotOffSet = None
        self.stepPlotCalled = {}
        self.stepPlotCalled['exact']=False; self.stepPlotCalled['elementQuantities']=False
        self.stepPlotCalled['ensight']=False; self.stepPlotCalled['ensightElementQuantities']=False
        self.plotWindowStart= {}
        self.nLevels    = nLevels
        self.flags = {}#force a deep copy?
        for key,val in SimulationProcessor.defaultFlags.items():
            self.flags[key] = val
        #mwf for postprocessing nodal values of coefficients etc
        self.nodalQuadratureInfo  = None
        #store p and n files now
        self.pFile = pFile; self.nFile = nFile
        if flags is not None:
            for key in list(self.flags.keys()):
                if key in list(flags.keys()):
                    self.flags[key]=flags[key]
                #end key found
            #end for all keys
        #end input flags given
        #need to check flags given are in allowed ranges
        self.errorData = {}
        for ci in self.flags['components']:
            self.errorData[ci] = {}
            for il in range(self.nLevels):
                self.errorData[ci][il] = {}
        #end ci
        if 'globalMassBalance' in self.flags['errorTypes']:
            for ci in self.flags['components']:
                for il in range(self.nLevels):
                    self.errorData[ci][il]['globalMass0']    = SimulationProcessor.notSetVal
                    self.errorData[ci][il]['globalMassF']    = [SimulationProcessor.notSetVal]
        if 'globalHeavisideMassBalance' in self.flags['errorTypes']:
            for ci in self.flags['components']:
                for il in range(self.nLevels):
                    self.errorData[ci][il]['globalHeavisideMass0']    = SimulationProcessor.notSetVal
                    self.errorData[ci][il]['globalHeavisideMassF']    = [SimulationProcessor.notSetVal]
        if 'localMassBalance' in self.flags['errorTypes']:
            self.conservationResidual = {}
            self.elementResidual = {}
            for ci in self.flags['components']:
                for il in range(self.nLevels):
                    self.errorData[ci][il]['localMassBalance'] = [SimulationProcessor.notSetVal]
                    self.conservationResidual[il] = None
                    self.elementResidual[il] = None
        #
        if 'numericalSolution' in self.flags['errorTypes']:
            for equant in self.flags['errorQuantities']:
                for enorm in self.flags['errorNorms']:
                    ktmp = 'error_'+equant+'_'+enorm
                    xtmp = 'exact_'+equant+'_'+enorm
                    #put in way to control cross component quantities
                    for ci in self.flags['components']:
                        for il in range(self.nLevels):
                            self.errorData[ci][il][ktmp] = []
                            self.errorData[ci][il][xtmp] = []
                        #for il
                    #for ci
                #end for enorm
            #end for
        #end if

        #save info about the spatial mesh?
        self.simulationData = {}
        self.simulationData['spatialMesh'] = {}
        #save values for levels and time steps
        for il in range(self.nLevels):
            self.simulationData['spatialMesh'][il] = {}
            self.simulationData['spatialMesh'][il]['nNodes_global']= []
            self.simulationData['spatialMesh'][il]['h'] = []
            self.simulationData['spatialMesh'][il]['hMin'] = []
        #end il
        self.solutionData = {}
        if ('u_dof' in self.flags['storeQuantities']):
            for ci in self.flags['components']:
                self.solutionData[ci] = {}
                for il in range(self.nLevels):
                    self.solutionData[ci][il] = {}
                    if 'u_dof' in self.flags['storeQuantities']:
                        self.solutionData[ci][il]['u_dof']= None
                        self.solutionData[ci][il]['l2g']  = None
                #il
            #ci
        #u_dof
        #mwf hack, decide if storing heavy data here, about to move fully to archiving
        self.storeHeavyData = False
        self.dataStorage = None
        #mwf debug
        #print "SimTools before dataStorage decision dataStorage=%s " % (self.dataStorage)
        if self.flags['storeTimes'] != [None]:
            import os
            import shelve
            absfile = os.path.join(self.flags['dataDir'],self.flags['dataFile'])
            if self.flags['appendResults'] == True:
                if not os.path.exists(absfile):
                    #figure out which exception to raise
                    assert 0, "SimTools append=True but couldn't find storage file=%s" % absfile
            else:
                if not os.path.exists(self.flags['dataDir']) and comm.isMaster():
                    os.makedirs(self.flags['dataDir'])
                if os.path.exists(absfile):
                    logEvent("Warning SimTools storing data removing old data in %s " % absfile)
                    #try:
                    os.remove(absfile)
                    #except:
                    #pass
            #
            self.dataStorage = shelve.open(absfile)
            #mwf debug
            #print "SimTools opening dataStorage file=%s dataStorage=%s " % (absfile,self.dataStorage)
            assert self.dataStorage is not None, "dataStorage is None storeTimes=%s absfile=%s " % (self.flags['storeTimes'],
                                                                                                absfile)

        #end storing something

    #end init

    def preprocess(self,mlvt,tsim):
        from . import Viewers
        """
        calculate desired quantities before simulation starts
        input :
          p    --- problem definition
          n    --- numerics definition
          mlvt --- multilevel vector transport object that has quantities for measuring
          tsim --- simulation time

          TO DO:
            local mass balance
            solution plotting of initial condition?
        """
        self.timeValues.append(tsim)
        p = self.pFile; n = self.nFile
        if 'globalMassBalance' in self.flags['errorTypes']:
            for il,m in enumerate(mlvt.levelModelList):
                for ci in range(p.coefficients.nc):
                    if ci in self.flags['components'] and ('m',ci) in m.q:
                        self.errorData[ci][il]['globalMass0'] = Norms.globalScalarDomainIntegral(m.q['abs(det(J))'],
                                                                                             m.elementQuadratureWeights[('m',ci)],
                                                                                             m.q[('m',ci)])
                        if self.flags['echo']:
                            logEvent("""t= %g globalMass0[%d][%d] = %g """ % (tsim,ci,il,
                                                                           self.errorData[ci][il]['globalMass0']))
                        #end if
                    #end if
                #end ci
            #end il
        #end calc globalMassBalance
        if 'globalHeavisideMassBalance' in self.flags['errorTypes']:
            for il,m in enumerate(mlvt.levelModelList):
                for ci in range(p.coefficients.nc):
                    if ci in self.flags['components'] and ('m',ci) in m.q:
                        hm = numpy.where(m.q[('m',ci)] >= 0.0,1.0,0.0)
                        self.errorData[ci][il]['globalHeavisideMass0'] = Norms.globalScalarDomainIntegral(m.q['abs(det(J))'],
                                                                                                  m.elementQuadratureWeights[('m',ci)],
                                                                                                  hm)
                        if self.flags['echo']:
                            logEvent("""t= %g globalHeavisideMass0[%d][%d] = %g """ % (tsim,ci,il,
                                                                                       self.errorData[ci][il]['globalHeavisideMass0']))
                        #end if
                    #end if
                #end ci
            #end il
        #end calc globalHeavisideMassBalance
        #
        #
        #now include velocity as element quantity if velocity flag is given in plotQuanties
        if 'velocity' in self.flags['plotQuantities']:
            has_q_velocity = False ; has_ebq_global_velocity = False
            for ci in self.flags['components']:
                pcikey = "q:('velocity',%s)" % ci
                if ('velocity',ci) in mlvt.levelModelList[-1].q:
                    has_q_velocity = True
                    if pcikey not in self.flags['plotQuantities']:
                        self.flags['plotQuantities'].append(pcikey)
                pcikey = "ebq:('velocity',%s)" % ci
                if ('velocity',ci) in mlvt.levelModelList[-1].ebq_global:
                    has_ebq_global_velocity = True
                    if pcikey not in self.flags['plotQuantities']:
                        self.flags['plotQuantities'].append(pcikey)
                #ebq_global
            #q
        #end velocity key fix
        ### set options for various output types ...
        for plotter in list(SimulationProcessor.defaultFlags['plotOptions'].keys()):
            if plotter not in self.flags['plotOptions']:
                self.flags['plotOptions'][plotter] = {'on':False}
        #mwf debug
        #import pdb
        #pdb.set_trace()
        #set on flags based on viewer
        if 'viewerType' in dir(Viewers):
            for plotter in list(self.flags['plotOptions'].keys()):
                if plotter == Viewers.viewerType:
                    self.flags['plotOptions'][plotter]['on']=True
        if (self.flags['plotOptions']['gnuplot']['on'] and
            'setGnuplotGridSize' not in self.flags['plotOptions']['gnuplot']):
            self.flags['plotOptions']['gnuplot']['setGnuplotGridSize'] = True
        if (self.flags['plotOptions']['matlab']['on'] and
            'usePDEtoolbox' not in self.flags['plotOptions']['matlab']):
            self.flags['plotOptions']['matlab']['usePDEtoolbox'] = False
        if self.flags['plotOptions']['vtk']['on']:
            if 'pause' not in self.flags['plotOptions']['vtk']:
                self.flags['plotOptions']['vtk']['pause']=False
            if 'hardcopy' not in self.flags['plotOptions']['vtk']:
                self.flags['plotOptions']['vtk']['hardcopy']=False
        #try to setup ensight header files correctly
        if self.flags['plotOptions']['ensight']['on']:
            if 'caseFileName' not in self.flags['plotOptions']['ensight']:
                self.flags['plotOptions']['ensight']['caseFileName'] = self.flags['simulationNameProc']
            mlvt.levelModelList[-1].u[0].femSpace.writeMeshEnsight(self.flags['plotOptions']['ensight']['caseFileName'],
                                                              self.flags['plotOptions']['ensight']['caseFileName'])
            self.ensightTimeSeries = []

            ensight_q_header_written = False; ensight_ebq_global_header_written = False;
            #NOTE: looks like ensight only allows one set of measured nodes
            #      so we enforce that either q: entries are plotted or ebq_global
            #element quadrature entries
            self.plottingQuadratureValuesForEnsight={'elements':False,'elementBoundaries':False}

            for quant in self.flags['plotQuantities']:
                recType = quant.split(':')
                if len(recType) > 1 and recType[0] == 'q': #found element quadrature quantity
                    stval = eval(recType[1])
                    if stval in mlvt.levelModelList[-1].q and not ensight_q_header_written:
                        self.writeEnsightMeshForElementQuantities(self.flags['plotOptions']['ensight']['caseFileName'],mlvt)
                        ensight_q_header_written = True
                        self.plottingQuadratureValuesForEnsight['elements']=True
                    #
                #q
            #quant
            if not self.plottingQuadratureValuesForEnsight['elements']:
                for quant in self.flags['plotQuantities']:
                    recType = quant.split(':')
                    if len(recType) > 1 and recType[0] == 'ebq_global': #found element boundary quadrature (global) quantity
                        stval = eval(recType[1])
                        if stval in mlvt.levelModelList[-1].ebq_global and not ensight_ebq_global_header_written:
                            self.writeEnsightMeshForElementBoundaryQuantities(self.flags['plotOptions']['ensight']['caseFileName'],mlvt)
                            ensight_ebq_global_header_written = True
                            self.plottingQuadratureValuesForEnsight['elementBoundaries']=True
                        #
                    #ebq_global
            #quantities
            #go ahead and write some or all of preamble?
            case_filename = self.flags['plotOptions']['ensight']['caseFileName']
            caseOut=open(case_filename+'.case','a')
            caseOut.write('VARIABLE\n')
            caseOut.close()
            mFinest = mlvt.levelModelList[-1]
            #solution values
            for ci in self.flags['components']:
                mFinest.u[ci].femSpace.writeFunctionHeaderEnsight(mFinest.u[ci],case_filename,append=False,
                                                                  firstVariable=False)
            #velocity dofs
            if mFinest.coefficients.vectorComponents is not None:
                if len(mFinest.coefficients.vectorComponents) == 2:
                    vcomp = [mFinest.coefficients.vectorComponents[0],
                             mFinest.coefficients.vectorComponents[1]]
                    mFinest.u[vcomp[0]].femSpace.writeE2VectorFunctionHeaderEnsight(mFinest.u[vcomp[0]],
                                                                                    mFinest.u[vcomp[1]],
                                                                                    case_filename,
                                                                                    nOutput=mFinest.u[vcomp[0]].femSpace.nOutput-1,
                                                                                    append=False,
                                                                                    firstVariable=False)
                elif len(mFinest.coefficients.vectorComponents) == 3:
                    vcomp = [mFinest.coefficients.vectorComponents[0],
                             mFinest.coefficients.vectorComponents[1],
                             mFinest.coefficients.vectorComponents[2]]
                    mFinest.u[vcomp[0]].femSpace.writeE3VectorFunctionHeaderEnsight(mFinest.u[vcomp[0]],
                                                                                    mFinest.u[vcomp[1]],
                                                                                    mFinest.u[vcomp[2]],
                                                                                    case_filename,
                                                                                    nOutput=mFinest.u[vcomp[0]].femSpace.nOutput-1,
                                                                                    append=False,
                                                                                    firstVariable=False)
            #
            for quant in self.flags['plotQuantities']:
                recType = quant.split(':')
                if len(recType) > 1 and recType[0] == 'q': #found element quadrature quantity
                    stval = eval(recType[1])
                    if (stval in mlvt.levelModelList[-1].q and
                        len(mlvt.levelModelList[-1].q[stval].shape) == 2): #found quantity and it's a scalar
                        self.writeScalarElementFunctionHeaderEnsight(stval,case_filename,append=False,firstVariable=False)
                    elif (stval in mlvt.levelModelList[-1].q and
                          len(mlvt.levelModelList[-1].q[stval].shape) == 3): #found quantity and it's a vector
                        self.writeVectorElementFunctionHeaderEnsight(stval,case_filename,append=False,firstVariable=False)
                    #vec
                #in q dict
            #element boundary (global) quadrature entries
            if not self.plottingQuadratureValuesForEnsight['elements']:
                for quant in self.flags['plotQuantities']:
                    recType = quant.split(':')
                    if len(recType) > 1 and recType[0] == 'ebq_global': #found element quadrature quantity
                        stval = eval(recType[1])
                        if (stval in mlvt.levelModelList[-1].ebq_global and
                            len(mlvt.levelModelList[-1].ebq_global[stval].shape) == 2): #found quantity and it's a scalar
                            self.writeScalarElementFunctionHeaderEnsight(stval,case_filename,append=False,firstVariable=False)
                        elif (stval in mlvt.levelModelList[-1].q and
                              len(mlvt.levelModelList[-1].q[stval].shape) == 3): #found quantity and it's a vector
                            self.writeVectorElementFunctionHeaderEnsight(stval,case_filename,append=False,firstVariable=False)
                        #vec
                    #in ebq_global dict
                #quant
            #end plot quantities

        #end ensight configuration
#cek moving to Viewers.V_base
#         if (('Init' in self.flags['plotTimes'] or 'All' in self.flags['plotTimes']) and
#             'u' in self.flags['plotQuantities'] and 'viewerType' in dir(Viewers)):#
#             #and  p.initialConditions is not None ):
#             dgrid = (n.nn-1)*(2**n.nLevels) #default should be 50
#             #mwf debug
#             #import pdb
#             #pdb.set_trace()

#             if self.plotOffSet is None:
#                self.plotOffSet = Viewers.windowNumber #keep from orphaning windows?
#             #don't reset window number
#             pause = False
#             if self.flags['plotOptions']['vtk']['on']:
#                 pause = self.flags['plotOptions']['vtk']['pause']
#             windowNumberTmp= mlvt.levelModelList[-1].viewSolution(plotOffSet=None,titleModifier=': Initial Condition',
#                                                                   dgridnx=dgrid,dgridny=dgrid,pause=pause)
#             #
#             self.stepPlotEnsight(mlvt,tsim)

#             #should create new windows if plotted here
#             self.stepPlotExact(mlvt,tsim)
#             self.stepPlotElementQuantities(mlvt,tsim)
#             self.stepPlotElementQuantitiesEnsight(mlvt,tsim)
#             if self.flags['plotOptions']['ensight']['on']:
#                 self.ensightTimeSeries.append(tsim)
#cek
        #
        if (('Init' in self.flags['storeTimes'] or 'All' in self.flags['storeTimes']) and
            p.initialConditions is not None):
            if 'u' in self.flags['storeQuantities']:
                mlvt.levelModelList[-1].saveSolution()
            self.stepStoreQuantities(mlvt,tsim)
        #end if
        if 'mesh' in self.flags['storeQuantities'] and 'mesh' not in self.dataStorage:
            #write out mesh information that is needed by at least ensight?
            meshDict = {}
            pm = mlvt.levelModelList[-1].mesh
            meshDict['nSpace_global']=mlvt.levelModelList[-1].nSpace_global
            mmemb = ['nNodes_global','nElements_global','nodeArray','elementNodesArray','nElementBoundaries_global',
                     'elementBoundaryNodesArray']
            for mm in mmemb:
                meshDict[mm] = getattr(pm,mm)
            self.dataStorage['mesh']=meshDict
            #mwf hack
            #self.dataStorage['wholeMesh']=mlvt.levelModelList[-1].mesh
    #end preproc

    def processTimeLevel(self,mlvt,tsim=None,plotOffSet=None):
        """calculate desired quantities after each macro time step

        Parameters
        ----------

        mlvt : multilevel vector transport that holds the quantities to measure
        tsim : simulation time
        
        """
        
#        input :
#          p    --- problem definition
#          n    --- numerics definition
#
        p = self.pFile; n = self.nFile
        if tsim is None:
            mlvt.levelModelList[-1].timeIntegration.t
        self.timeValues.append(tsim)
        if plotOffSet is not None:
            self.plotOffSet = plotOffSet
        if 'All' in self.flags['errorTimes'] or tsim in self.flags['errorTimes']:
            self.stepProcessError(mlvt,tsim)
        if 'All' in self.flags['storeTimes'] or tsim in self.flags['storeTimes']:
            mlvt.levelModelList[-1].saveSolution()
            self.stepStoreQuantities(mlvt,tsim)
        return self.plotOffSet
    def postprocess(self,mlvt,tsim):
        """
        calculate desired quantities after simulation ends

        Parameters
        ----------

          mlvt : multilevel vector transport that holds the quantities to measure
          tsim : simulation time

        """
        p = self.pFile; n = self.nFile

        if 'Last' in self.flags['errorTimes'] or tsim in self.flags['errorTimes']:
            self.stepProcessError(mlvt,tsim)

        #now look at global mass balance (could do similar thing each step too)
        if 'globalMassBalance' in self.flags['errorTypes']:
            for il,m in enumerate(mlvt.levelModelList):
                for ci in range(p.coefficients.nc):
                    if ci in self.flags['components'] and ('m',ci) in m.q:
                        self.errorData[ci][il]['globalMassF'].append(Norms.globalScalarDomainIntegral(m.q['abs(det(J))'][0:m.mesh.subdomainMesh.nElements_owned],
                                                                                                      m.elementQuadratureWeights[('m',ci)],
                                                                                                      m.q[('m',ci)][0:m.mesh.subdomainMesh.nElements_owned]))
                        if self.flags['echo']:
                            logEvent("""t= %g globalMassF[%d][%d] = %g globalMassDiff[%d][%d]= %g""" % \
                                         (tsim,ci,il,self.errorData[ci][il]['globalMassF'][-1],ci,il,
                                          self.errorData[ci][il]['globalMassF'][-1]-self.errorData[ci][il]['globalMass0']))
                        #end if
                    #end if
                #end for ci
            #end for il
        #end calc globalMassBalance
        if 'globalHeavisideMassBalance' in self.flags['errorTypes']:
            for il,m in enumerate(mlvt.levelModelList):
                for ci in range(p.coefficients.nc):
                    if ci in self.flags['components'] and ('m',ci) in m.q:
                        hm = numpy.where(m.q[('m',ci)][0:m.mesh.subdomainMesh.nElements_owned] >= 0.0,1.0,0.0)
                        self.errorData[ci][il]['globalHeavisideMassF'].append(Norms.globalScalarDomainIntegral(m.q['abs(det(J))'][0:m.mesh.subdomainMesh.nElements_owned],
                                                                                                       m.elementQuadratureWeights[('m',ci)],
                                                                                                       hm))
                        if self.flags['echo']:
                            logEvent("""t= %g globalHeavisideMassF[%d][%d] = %g globalHeavisideMassDiff[%d][%d]= %g""" % \
                                         (tsim,ci,il,self.errorData[ci][il]['globalHeavisideMassF'][-1],ci,il,
                                          abs(self.errorData[ci][il]['globalHeavisideMassF'][-1]-self.errorData[ci][il]['globalHeavisideMass0'])))
                        #end if
                    #end if
                #end for ci
            #end for il
        #end calc globalMassBalance

        #compute space-time norms, ...
        if 'numericalSolution' in self.flags['errorTypes']:
            for il,m in enumerate(mlvt.levelModelList):
                for ci in range(p.coefficients.nc):
                    if ci in self.flags['components']:
                        for snorm in ['L2','L1','LI','H1','H1semi','W11','W11semi']:
                            kerr = 'error_'+'u'+'_'+snorm
                            kexa = 'exact_'+'u'+'_'+snorm
                            calcNorm = kerr in self.errorData[ci][il]
                            if calcNorm and 'L2_'+snorm in self.flags['errorNorms']:
                                errTL2 = 0.0
                                eLast  = 0.0
                                exaTL2 = 0.0
                                exaLast= 0.0
                                for it in range(1,len(self.timeValues)):
                                    dt = self.timeValues[it]-self.timeValues[it-1] #timeValues has t0
                                    errTL2 += 0.5*dt*(self.errorData[ci][il][kerr][it-1]**2 +
                                                      eLast**2)
                                    eLast = self.errorData[ci][il][kerr][it-1]
                                    exaTL2 += 0.5*dt*(self.errorData[ci][il][kexa][it-1]**2 +
                                                      exaLast**2)
                                    exaLast = self.errorData[ci][il][kexa][it-1]
                                #end it
                                kerrtL2 = kerr+'_L2'
                                kexatL2 = kexa+'_L2'
                                self.errorData[ci][il][kerrtL2] = self.math.sqrt(errTL2)
                                self.errorData[ci][il][kexatL2] = self.math.sqrt(exaTL2)
                                if self.flags['echo']:
                                    print("""t= %g; %s[%d][%d]= %g;""" % (tsim,kerrtL2,ci,il,errTL2))
                                #end if
                            if calcNorm and 'L1_'+snorm in self.flags['errorNorms']:
                                errTL1 = 0.0
                                eLast  = 0.0
                                #mwf debug
                                #print """postproc timeValues= %s """ % self.timeValues
                                exaTL1 = 0.0
                                for it in range(1,len(self.timeValues)):
                                    dt = self.timeValues[it]-self.timeValues[it-1] #timeValues has t0
                                    errTL1+= 0.5*dt*(self.errorData[ci][il][kerr][it-1] +
                                                     eLast)
                                    eLast = self.errorData[ci][il][kerr][it-1]
                                    #mwf debug
                                    #print """postproc timeValues[%d]=%g err[%d]=%g errTL1=%g """ \
                                    #      % (it,self.timeValues[it],it-1,eLast,errTL1)
                                    exaTL1+= 0.5*dt*(self.errorData[ci][il][kexa][it-1] +
                                                     exaLast)
                                    exaLast = self.errorData[ci][il][kexa][it-1]
                                #end it
                                kerrtL1 = kerr+'_L1'
                                kexatL1 = kexa+'_L1'
                                self.errorData[ci][il][kerrtL1] = errL1TV
                                self.errorData[ci][il][kexatL1] = exaL1TV
                                if self.flags['echo']:
                                    print("""t= %g; %s[%d][%d]= %g;""" % (tsim,kerrtL1,ci,il,errL1TV))
                                #end if
                            #if calcL1+snorm
                            if calcNorm and 'LI_'+snorm in self.flags['errorNorms']:
                                errTLI = 0.0
                                exaTLI = 0.0
                                for it in range(1,len(self.timeValues)):
                                    errTLI = max(errTLI,self.errorData[ci][il][kerr][it-1])
                                    exaTLI = max(exaTLI,self.errorData[ci][il][kexa][it-1])
                                #it
                                kerrtLI = kerr+'_LI' ;
                                kexatLI = kexa+'_LI' ;
                                self.errorData[ci][il][kerrtLI] = errTLI
                                self.errorData[ci][il][kexatLI] = exaTLI
                                if self.flags['echo']:
                                    print("""t= %g; %s[%d][%d]= %g;""" % (tsim,kerrtLI,ci,il,errTLI))
                                #end if
                            #calcLI norm
                        #end space norms
                    #end if calc
                #end for ci
            #end for il
        if self.flags['plotOptions']['ensight']['on']:
            mlvt.levelModelList[-1].u[0].femSpace.endTimeSeriesEnsight(self.ensightTimeSeries,
                                                                  self.flags['plotOptions']['ensight']['caseFileName'],
                                                                  self.flags['plotOptions']['ensight']['caseFileName'])

        #
        #assume that 'u' taken care of in step process???
        if 'Last' in self.flags['storeTimes']:
            self.stepStoreQuantities(mlvt,tsim)
            #now in stepStoreQuantities?
            #if 'u_dof' in self.flags['storeQuantities']:
            #    for ci in self.flags['components']:
            #        for il,m in enumerate(mlvt.levelModelList):
            #            self.solutionData[ci][il]['u_dof'] = m.u[ci].dof
            #            self.solutionData[ci][il]['l2g'] = m.u[ci].femSpace.dofMap.l2g
            #        #for
            #    #for
            ##if
            #store everything basically? ...
            #if 'multilevelModel' in self.flags['storeQuantities']:
            #    #need to make
            #    self.dataStorage[('multilevelModel',tsim)]=mlvt
#
        #if
        self.saveToDisk()


    def saveToDisk(self):
        """
        TO DO
        save error to disk
        make sure can append if necessary?
        """
        if self.flags['storeTimes'] != [None]:
            assert self.dataStorage is not None, "dataStorage None storeTimes= %s " % self.flags['storeTimes']
            if 'simulationData' in self.flags['storeQuantities']:
                self.dataStorage['timeValues']    = self.timeValues
                self.dataStorage['simulationData']= self.simulationData
                self.dataStorage['flags']         = self.flags
            if 'errorData' in self.flags['storeQuantities']:
                self.dataStorage['errorData']=self.errorData
            #what else to store?
            if ('u_dof' in self.flags['storeQuantities']):
                self.dataStorage['solutionData'] = self.solutionData
            #end solutionData
            self.dataStorage.close()
            #
        #don't store anything
    #end def

    def stepProcessError(self,mlvt,tsim):
        """ calculate desired error quantities for a single step

        Parameters
        ----------

          mlvt : multilevel vector transport that holds the quantities to measure
          tsim : simulation time
        
        """
#        TO DO:
#          synchronize Norms L*error*AF[,2] functions used to calculate error
#          setup to work in parallel
        p = self.pFile; n = self.nFile
        for il,m in enumerate(mlvt.levelModelList):
            self.simulationData['spatialMesh'][il]['nNodes_global'].append(m.mesh.nNodes_global)
            self.simulationData['spatialMesh'][il]['h'].append(m.mesh.h)
            self.simulationData['spatialMesh'][il]['hMin'].append(m.mesh.hMin)
        #end il
        #assumes this is a correct time to compute error
        hasAnalyticalSolution = {}
        hasAnalyticalSolutionVelocity = {}
        for ci in range(p.coefficients.nc):
            hasAnalyticalSolution[ci] = (ci in self.analyticalSolution  and
                                         self.analyticalSolution[ci] is not None)
            hasAnalyticalSolutionVelocity[ci] = ('analyticalSolutionVelocity' in dir(p) and
                                                 p.analyticalSolutionVelocity is not None and
                                                 ci in p.analyticalSolutionVelocity and
                                                 p.analyticalSolutionVelocity[ci] is not None)
        #ci
        class gradWrapper(object):
            def __init__(self,ex):
                self.ex = ex
            def  uOfX(self,X):
                return self.ex.duOfX(X)
            def uOfXT(self,X,T):
                return self.ex.duOfXT(X,T)
        #grad wrapper
        import math
        if 'numericalSolution' in self.flags['errorTypes']:
            for il,m in enumerate(mlvt.levelModelList):
                #first see if need to project solution for calculations at all
                needProj = False
                uproj = None; uprojGrad = None
                velproj = {}
                for ci in range(p.coefficients.nc):
                    if ci in self.flags['components'] and not hasAnalyticalSolution[ci]:
                        needProj = True
                if needProj:
                    uproj,uprojGrad = projectToFinestLevel(mlvt,il,tsim)#assumes conforming unless used MultilevelTransfer...NC
                #mwf hack now allow evaluation of postprocessed velocities on fine grid
                #import testStuff
                for ci in range(p.coefficients.nc):
                    if (ci in self.flags['components']and
                        not hasAnalyticalSolutionVelocity[ci] and
                        n.conservativeFlux is not None and 'velocity' in self.flags['errorQuantities']):
                        #mwf debug
                        logEvent("SimTools proj velocity for error calling projectVelocityToFinestLevelNC")
                        velproj[ci] = projectVelocityToFinestLevelNC(mlvt,il,ci)

                # CALCULATE THE L2 ERROR IN PRESSURE 
                if 'p' in self.flags['errorQuantities'] and hasattr(m,'analyticalPressureSolution'):
                    assert hasattr(m,'analyticalPressureSolution'), "analyticalPressureSolution must be provided"
                    # COMPUTE MEAN VALUE OF PRESSURE
                    pressureAnalyticalSolution = m.analyticalPressureSolution[0]
                    x = m.q['x'][0:m.mesh.subdomainMesh.nElements_owned]                        
                    abs_det_J = m.q['abs(det(J))'][0:m.mesh.subdomainMesh.nElements_owned]
                    quad_weight = list(m.elementQuadratureWeights.values())[0]
                    pressureNumericalSolution = m.q['p'][0:m.mesh.subdomainMesh.nElements_owned]
                    # compute mean values
                    mean_value_exact_p = 0.0
                    mean_value_numerical_p = 0.0
                    for eN in range (x.shape[0]):
                        for k in range(x.shape[1]):
                            mean_value_exact_p += pressureAnalyticalSolution.uOfXT(x[eN,k],tsim)*quad_weight[k]*abs_det_J[eN,k]
                            mean_value_numerical_p += pressureNumericalSolution[eN,k]*quad_weight[k]*abs_det_J[eN,k]
                    # remove mean value of numerical solution and add mean value of exact solution 
                    pressureNumericalSolution += mean_value_exact_p - mean_value_numerical_p
                    err = Norms.L2errorSFEMvsAF2(pressureAnalyticalSolution, 
                                                 x, 
                                                 abs_det_J, 
                                                 quad_weight,
                                                 pressureNumericalSolution, 
                                                 T=tsim)
                    kerr = 'error_'+'p'+'_'+'L2'
                    if self.flags['echo']:
                        logEvent("""\nt= %g; %s= %g;""" % (tsim,kerr,err),level=0)
                # END OF COMPUTING THE L2 ERROR OF THE PRESSURE 

                for ci in range(p.coefficients.nc):
                    if ci in self.flags['components']:
                        if not hasAnalyticalSolution[ci]:
                            udense     = mlvt.levelModelList[-1].q[('u',ci)]
                            gradu_dense= mlvt.levelModelList[-1].q[('grad(u)',ci)]
                        if not hasAnalyticalSolutionVelocity[ci] and 'velocity' in self.flags['errorQuantities']:
                            veldense = mlvt.levelModelList[-1].q[('velocity',ci)]
                        mFine = mlvt.levelModelList[-1]

                        calcL2u = ('L2' in self.flags['errorNorms'] and
                                   'u' in self.flags['errorQuantities'])
                        if calcL2u:
                            err = -12345.0
                            exa = 1.0
                            if hasAnalyticalSolution[ci]:
                                err = Norms.L2errorSFEMvsAF2(self.analyticalSolution[ci],
                                                             m.q['x'][0:m.mesh.subdomainMesh.nElements_owned],
                                                             m.q['dV'][0:m.mesh.subdomainMesh.nElements_owned],
                                                             np.ones_like(list(m.elementQuadratureWeights.values())[0]),
                                                             m.q[('u',ci)][0:m.mesh.subdomainMesh.nElements_owned],
                                                             T=tsim)
                                exa = Norms.L2errorSFEMvsAF2(zeroFunction(),
                                                             m.q['x'][0:m.mesh.subdomainMesh.nElements_owned],
                                                             m.q['dV'][0:m.mesh.subdomainMesh.nElements_owned],
                                                             np.ones_like(list(m.elementQuadratureWeights.values())[0]),
                                                             m.q[('u',ci)][0:m.mesh.subdomainMesh.nElements_owned],
                                                             T=tsim)
                            else:
                                err = Norms.L2errorSFEM(mFine.q[('dV_u',ci)][0:mFine.mesh.subdomainMesh.nElements_owned],udense[0:mFine.mesh.subdomainMesh.nElements_owned],uproj[ci][0:mFine.mesh.subdomainMesh.nElements_owned])
                                exa = Norms.L2errorSFEM(mFine.q[('dV_u',ci)][0:mFine.mesh.subdomainMesh.nElements_owned],udense[0:mFine.mesh.subdomainMesh.nElements_owned],
                                                        numpy.zeros(udense[0:mFine.mesh.subdomainMesh.nElements_owned].shape,'d'))
                            kerr = 'error_'+'u'+'_'+'L2'
                            kexa = 'exact_'+'u'+'_'+'L2'
                            self.errorData[ci][il][kerr].append(err)
                            self.errorData[ci][il][kexa].append(exa)
                            if self.flags['echo']:
                                if self.flags['echoRelativeErrors']:
                                    logEvent("""\nt= %g; %s[%d][%d]= %g; relative_error= %g;""" % (tsim,kerr,ci,il,err,err/(exa+1E-15)),level=0)
                                else:
                                    logEvent("""\nt= %g; %s[%d][%d]= %g;""" % (tsim,kerr,ci,il,err),level=0)
                                #end if
                        #if calcL2u
                        calcL1u = ('L1' in self.flags['errorNorms'] and
                                   'u' in self.flags['errorQuantities'])
                        if calcL1u:
                            err = -12345.0
                            exa = 1.0
                            if hasAnalyticalSolution[ci]:
                                err = Norms.L1errorSFEMvsAF2(self.analyticalSolution[ci],
                                                             m.q['x'][0:m.mesh.subdomainMesh.nElements_owned],
                                                             m.q['dV'][0:m.mesh.subdomainMesh.nElements_owned],
                                                             np.ones_like(list(m.elementQuadratureWeights.values())[0]),
                                                             m.q[('u',ci)][0:m.mesh.subdomainMesh.nElements_owned],
                                                             T=tsim)
                                exa = Norms.L1errorSFEMvsAF2(zeroFunction(),
                                                             m.q['x'][0:m.mesh.subdomainMesh.nElements_owned],
                                                             m.q['dV'][0:m.mesh.subdomainMesh.nElements_owned],
                                                             np.ones_like(list(m.elementQuadratureWeights.values())[0]),
                                                             m.q[('u',ci)][0:m.mesh.subdomainMesh.nElements_owned],
                                                             T=tsim)
                            else:
                                err = Norms.L1errorSFEM(mFine.q[('dV_u',ci)][0:mFine.mesh.subdomainMesh.nElements_owned],udense[0:mFine.mesh.subdomainMesh.nElements_owned],uproj[ci][0:mFine.mesh.subdomainMesh.nElements_owned])
                                exa = Norms.L1errorSFEM(mFine.q[('dV_u',ci)][0:mFine.mesh.subdomainMesh.nElements_owned],udense[0:mFine.mesh.subdomainMesh.nElements_owned],
                                                        numpy.zeros(udense[0:mFine.mesh.subdomainMesh.nElements_owned].shape,'d'))
                            kerr = 'error_'+'u'+'_'+'L1'
                            kexa = 'exact_'+'u'+'_'+'L1'
                            self.errorData[ci][il][kerr].append(err)
                            self.errorData[ci][il][kexa].append(exa)
                            if self.flags['echo']:
                                if self.flags['echoRelativeErrors']:
                                    logEvent("""\nt= %g; %s[%d][%d]= %g; relative_error= %g;""" % (tsim,kerr,ci,il,err,err/(exa+1E-15)),level=0)
                                else:
                                    logEvent("""\nt= %g; %s[%d][%d]= %g;""" % (tsim,kerr,ci,il,err),level=0)                                
                            #end if
                        #if calcL1u
                        calcLIu = ('LI' in self.flags['errorNorms'] and
                                   'u' in self.flags['errorQuantities'])
                        if calcLIu:
                            err = -12345.0
                            exa = 0.0

                            if hasAnalyticalSolution[ci]:
                                err = Norms.LIerrorSFEMvsAF(self.analyticalSolution[ci],
                                                            m.q['x'][0:m.mesh.subdomainMesh.nElements_owned],
                                                            m.q[('u',ci)][0:m.mesh.subdomainMesh.nElements_owned],
                                                            T=tsim)
                                exa = Norms.LIerrorSFEMvsAF(zeroFunction(),
                                                            m.q['x'][0:m.mesh.subdomainMesh.nElements_owned],
                                                            m.q[('u',ci)][0:m.mesh.subdomainMesh.nElements_owned],
                                                            T=tsim)
                            else:
                                err = max(numpy.absolute(uproj[ci][0:mFine.mesh.subdomainMesh.nElements_owned].flat[:]-udense[0:mFine.mesh.subdomainMesh.nElements_owned].flat[:]))
                                exa = max(numpy.absolute(udense[0:mFine.mesh.subdomainMesh.nElements_owned].flat))

                            kerr = 'error_'+'u'+'_'+'LI'
                            kexa = 'exact_'+'u'+'_'+'LI'
                            self.errorData[ci][il][kerr].append(err)
                            self.errorData[ci][il][kexa].append(exa)
                            if self.flags['echo']:
                                if self.flags['echoRelativeErrors']:
                                    logEvent("""\nt= %g; %s[%d][%d]= %g; relative_error= %g;""" % (tsim,kerr,ci,il,err,err/(exa+1E-15)),level=0)
                                else:
                                    logEvent("""\nt= %g; %s[%d][%d]= %g;""" % (tsim,kerr,ci,il,err),level=0)
                            #end if
                            #
                        #calcLIu
                        calcH1u = ('H1' in self.flags['errorNorms'] and
                                   'u' in self.flags['errorQuantities'])

                        if calcH1u:
                            err = -12345.0; err0 = -12345.0; err1 = -12345.0
                            exa = 1.0; exa0 = 1.0; exa1 = 1.0
                            if hasAnalyticalSolution[ci]:
                                err0 = Norms.L2errorSFEMvsAF2(self.analyticalSolution[ci],
                                                              m.q['x'][0:m.mesh.subdomainMesh.nElements_owned],
                                                              m.q['abs(det(J))'][0:m.mesh.subdomainMesh.nElements_owned],
                                                              list(m.elementQuadratureWeights.values())[0],
                                                              m.q[('u',ci)][0:m.mesh.subdomainMesh.nElements_owned],
                                                              T=tsim)

                                exa0 = Norms.L2errorSFEMvsAF2(zeroFunction(),
                                                              m.q['x'][0:m.mesh.subdomainMesh.nElements_owned],
                                                              m.q['abs(det(J))'][0:m.mesh.subdomainMesh.nElements_owned],
                                                              list(m.elementQuadratureWeights.values())[0],
                                                              m.q[('u',ci)][0:m.mesh.subdomainMesh.nElements_owned],
                                                              T=tsim)
                            else:

                                err0 = Norms.L2errorSFEM(mlvt.levelModelList[-1].q[('dV_u',ci)][0:mFine.mesh.subdomainMesh.nElements_owned],udense[0:mFine.mesh.subdomainMesh.nElements_owned],uproj[ci][0:mFine.mesh.subdomainMesh.nElements_owned])
                                exa0 = Norms.L2errorSFEM(mlvt.levelModelList[-1].q[('dV_u',ci)][0:mFine.mesh.subdomainMesh.nElements_owned],udense[0:mFine.mesh.subdomainMesh.nElements_owned],
                                                         numpy.zeros(udense[0:mFine.mesh.subdomainMesh.nElements_owned].shape,'d'))
                            #now gradients
                            if hasAnalyticalSolution[ci]:
                                err1 = Norms.L2errorVFEMvsAF(gradWrapper(p.analyticalSolution[ci]),
                                                             m.q['x'][0:m.mesh.subdomainMesh.nElements_owned],
                                                             m.q[('dV_u',ci)][0:m.mesh.subdomainMesh.nElements_owned],
                                                             m.q[('grad(u)',ci)][0:m.mesh.subdomainMesh.nElements_owned],
                                                             T=tsim)
                                exa1 = Norms.L2errorVFEMvsAF(gradWrapper(p.analyticalSolution[ci]),
                                                             m.q['x'][0:m.mesh.subdomainMesh.nElements_owned],
                                                             m.q[('dV_u',ci)][0:m.mesh.subdomainMesh.nElements_owned],
                                                             numpy.zeros(m.q[('grad(u)',ci)][0:m.mesh.subdomainMesh.nElements_owned].shape,
                                                                           'd'),
                                                             T=tsim)
                            else:
                                err1 = Norms.L2errorVFEM(mlvt.levelModelList[-1].q[('dV_u',ci)][0:mFine.mesh.subdomainMesh.nElements_owned],
                                                         gradu_dense[0:mFine.mesh.subdomainMesh.nElements_owned],uprojGrad[ci][0:mFine.mesh.subdomainMesh.nElements_owned])
                                exa1 = Norms.L2errorVFEM(mlvt.levelModelList[-1].q[('dV_u',ci)][0:mFine.mesh.subdomainMesh.nElements_owned],
                                                        gradu_dense[0:mFine.mesh.subdomainMesh.nElements_owned],
                                                        numpy.zeros(gradu_dense[0:mFine.mesh.subdomainMesh.nElements_owned].shape,
                                                                      'd'))

                            err = math.sqrt(err0**2 + err1**2)
                            exa = math.sqrt(exa0**2 + exa1**2)
                            kerr = 'error_'+'u'+'_'+'H1'
                            kexa = 'exact_'+'u'+'_'+'H1'
                            self.errorData[ci][il][kerr].append(err)
                            self.errorData[ci][il][kexa].append(exa)
                            if self.flags['echo']:
                                if self.flags['echoRelativeErrors']:
                                    logEvent("""\nt= %g; %s[%d][%d]= %g; relative_error= %g;""" % (tsim,kerr,ci,il,err,err/(exa+1E-15)),level=0)
                                else:
                                    logEvent("""\nt= %g; %s[%d][%d]= %g;""" % (tsim,kerr,ci,il,err),level=0)
                            #end if
                        #if calcH1u
                        calcH1semiU = ('H1semi' in self.flags['errorNorms'] and
                                       'u' in self.flags['errorQuantities'])

                        if calcH1semiU:
                            err = -12345.0;
                            exa = 1.0;
                            if hasAnalyticalSolution[ci]:
                                err = Norms.L2errorVFEMvsAF(gradWrapper(p.analyticalSolution[ci]),
                                                            m.q['x'][0:m.mesh.subdomainMesh.nElements_owned],
                                                            m.q[('dV_u',ci)][0:m.mesh.subdomainMesh.nElements_owned],
                                                            m.q[('grad(u)',ci)][0:m.mesh.subdomainMesh.nElements_owned],
                                                            T=tsim)
                                exa = Norms.L2errorVFEMvsAF(gradWrapper(p.analyticalSolution[ci]),
                                                            m.q['x'][0:m.mesh.subdomainMesh.nElements_owned],
                                                            m.q[('dV_u',ci)][0:m.mesh.subdomainMesh.nElements_owned],
                                                            numpy.zeros(m.q[('grad(u)',ci)][0:m.mesh.subdomainMesh.nElements_owned].shape,
                                                                          'd'),
                                                            T=tsim)
                            else:
                                err = Norms.L2errorVFEM(mlvt.levelModelList[-1].q[('dV_u',ci)][0:mFine.mesh.subdomainMesh.nElements_owned],
                                                        gradu_dense[0:mFine.mesh.subdomainMesh.nElements_owned],uprojGrad[ci][0:mFine.mesh.subdomainMesh.nElements_owned])
                                exa = Norms.L2errorVFEM(mlvt.levelModelList[-1].q[('dV_u',ci)][0:mFine.mesh.subdomainMesh.nElements_owned],
                                                        gradu_dense[0:mFine.mesh.subdomainMesh.nElements_owned],
                                                        numpy.zeros(gradu_dense[0:mFine.mesh.subdomainMesh.nElements_owned].shape,
                                                                      'd'))

                            kerr = 'error_'+'u'+'_'+'H1semi'
                            kexa = 'exact_'+'u'+'_'+'H1semi'
                            self.errorData[ci][il][kerr].append(err)
                            self.errorData[ci][il][kexa].append(exa)
                            if self.flags['echo']:
                                if self.flags['echoRelativeErrors']:
                                    logEvent("""\nt= %g; %s[%d][%d]= %g; relative_error= %g;""" % (tsim,kerr,ci,il,err,err/(exa+1E-15)),level=0)
                                else:
                                    logEvent("""\nt= %g; %s[%d][%d]= %g;""" % (tsim,kerr,ci,il,err),level=0)
                            #end if
                        #if calcH1semiu
                        calcW11u = ('W11' in self.flags['errorNorms'] and
                                    'u' in self.flags['errorQuantities'])

                        if calcW11u:
                            err = -12345.0; err0 = -12345.0; err1 = -12345.0
                            exa = 1.0; exa0 = 1.0; exa1 = 1.0
                            if hasAnalyticalSolution[ci]:
                                err0 = Norms.L1errorSFEMvsAF2(self.analyticalSolution[ci],
                                                              m.q['x'][0:m.mesh.subdomainMesh.nElements_owned],
                                                              m.q['abs(det(J))'][0:m.mesh.subdomainMesh.nElements_owned],
                                                              list(m.elementQuadratureWeights.values())[0],
                                                              m.q[('u',ci)][0:m.mesh.subdomainMesh.nElements_owned],
                                                              T=tsim)

                                exa0 = Norms.L1errorSFEMvsAF2(zeroFunction(),
                                                              m.q['x'][0:m.mesh.subdomainMesh.nElements_owned],
                                                              m.q['abs(det(J))'][0:m.mesh.subdomainMesh.nElements_owned],
                                                              list(m.elementQuadratureWeights.values())[0],
                                                              m.q[('u',ci)][0:m.mesh.subdomainMesh.nElements_owned],
                                                              T=tsim)
                            else:
                                err0 = Norms.L1errorSFEM(mlvt.levelModelList[-1].q[('dV_u',ci)][0:mFine.mesh.subdomainMesh.nElements_owned],udense[0:mFine.mesh.subdomainMesh.nElements_owned],uproj[ci][0:mFine.mesh.subdomainMesh.nElements_owned])
                                exa0 = Norms.L1errorSFEM(mlvt.levelModelList[-1].q[('dV_u',ci)][0:mFine.mesh.subdomainMesh.nElements_owned],udense[0:mFine.mesh.subdomainMesh.nElements_owned],
                                                         numpy.zeros(udense[0:mFine.mesh.subdomainMesh.nElements_owned].shape,'d'))
                            #now gradients
                            if hasAnalyticalSolution[ci]:
                                err1 = Norms.L1errorVFEMvsAF(gradWrapper(p.analyticalSolution[ci]),
                                                             m.q['x'][0:m.mesh.subdomainMesh.nElements_owned],
                                                             m.q[('dV_u',ci)][0:m.mesh.subdomainMesh.nElements_owned],
                                                             m.q[('grad(u)',ci)][0:m.mesh.subdomainMesh.nElements_owned],
                                                             T=tsim)
                                exa1 = Norms.L1errorVFEMvsAF(gradWrapper(p.analyticalSolution[ci]),
                                                             m.q['x'][0:m.mesh.subdomainMesh.nElements_owned],
                                                             m.q[('dV_u',ci)][0:m.mesh.subdomainMesh.nElements_owned],
                                                             numpy.zeros(m.q[('grad(u)',ci)][0:m.mesh.subdomainMesh.nElements_owned].shape,
                                                                           'd'),
                                                             T=tsim)
                            else:
                                err1 = Norms.L1errorVFEM(mlvt.levelModelList[-1].q[('dV_u',ci)][0:mFine.mesh.subdomainMesh.nElements_owned],
                                                         gradu_dense[0:mFine.mesh.subdomainMesh.nElements_owned],uprojGrad[ci][0:mFine.mesh.subdomainMesh.nElements_owned])
                                exa1 = Norms.L1errorVFEM(mlvt.levelModelList[-1].q[('dV_u',ci)][0:mFine.mesh.subdomainMesh.nElements_owned],
                                                        gradu_dense[0:mFine.mesh.subdomainMesh.nElements_owned],
                                                        numpy.zeros(gradu_dense[0:mFine.mesh.subdomainMesh.nElements_owned].shape,
                                                                      'd'))

                            err = err0 + err1
                            exa = exa0 + exa1
                            kerr = 'error_'+'u'+'_'+'W11'
                            kexa = 'exact_'+'u'+'_'+'W11'
                            self.errorData[ci][il][kerr].append(err)
                            self.errorData[ci][il][kexa].append(exa)
                            if self.flags['echo']:
                                if self.flags['echoRelativeErrors']:
                                    logEvent("""\nt= %g; %s[%d][%d]= %g; relative_error= %g;""" % (tsim,kerr,ci,il,err,err/(exa+1E-15)),level=0)
                                else:
                                    logEvent("""\nt= %g; %s[%d][%d]= %g;""" % (tsim,kerr,ci,il,err),level=0)
                            #end if
                        #if calcW11u
                        calcW11semiU = ('W11semi' in self.flags['errorNorms'] and
                                        'u' in self.flags['errorQuantities'])
                        if calcW11semiU:
                            err = -12345.0;
                            exa = 1.0;
                            if hasAnalyticalSolution[ci]:
                                err = Norms.L1errorVFEMvsAF(gradWrapper(p.analyticalSolution[ci]),
                                                            m.q['x'][0:m.mesh.subdomainMesh.nElements_owned],
                                                            m.q[('dV_u',ci)][0:m.mesh.subdomainMesh.nElements_owned],
                                                            m.q[('grad(u)',ci)][0:m.mesh.subdomainMesh.nElements_owned],
                                                            T=tsim)
                                exa = Norms.L1errorVFEMvsAF(gradWrapper(p.analyticalSolution[ci]),
                                                            m.q['x'][0:m.mesh.subdomainMesh.nElements_owned],
                                                            m.q[('dV_u',ci)][0:m.mesh.subdomainMesh.nElements_owned],
                                                            numpy.zeros(m.q[('grad(u)',ci)][0:m.mesh.subdomainMesh.nElements_owned].shape,
                                                                          'd'),
                                                            T=tsim)
                            else:
                                err = Norms.L1errorVFEM(mlvt.levelModelList[-1].q[('dV_u',ci)][0:mFine.mesh.subdomainMesh.nElements_owned],
                                                        gradu_dense[0:mFine.mesh.subdomainMesh.nElements_owned],uprojGrad[ci][0:mFine.mesh.subdomainMesh.nElements_owned])
                                exa = Norms.L1errorVFEM(mlvt.levelModelList[-1].q[('dV_u',ci)][0:mFine.mesh.subdomainMesh.nElements_owned],
                                                        gradu_dense[0:mFine.mesh.subdomainMesh.nElements_owned],
                                                        numpy.zeros(gradu_dense[0:mFine.mesh.subdomainMesh.nElements_owned].shape,
                                                                      'd'))

                            kerr = 'error_'+'u'+'_'+'W11semi'
                            kexa = 'exact_'+'u'+'_'+'W11semi'
                            self.errorData[ci][il][kerr].append(err)
                            self.errorData[ci][il][kexa].append(exa)
                            if self.flags['echo']:
                                if self.flags['echoRelativeErrors']:
                                    logEvent("""\nt= %g; %s[%d][%d]= %g; relative_error= %g;""" % (tsim,kerr,ci,il,err,err/(exa+1E-15)),level=0)
                                else:
                                    logEvent("""\nt= %g; %s[%d][%d]= %g;""" % (tsim,kerr,ci,il,err),level=0)
                            #end if
                        #if calcH1semiu
                        calcTVu = ('TV' in self.flags['errorNorms'] and
                                   'u' in self.flags['errorQuantities'])
                        if calcTVu:
                            err = -12345.0
                            exa = 1.0#need to compute using interpolation conditions?
                            #one version of discrete seminorm
                            #err = Norms.TVseminormSFEM(m.u[ci].dof,m.u[ci].femSpace.dofMap.l2g)
                            #what about just using W11seminorm
                            shtmp = (m.q[('grad(u)',ci)][0:m.mesh.subdomainMesh.nElements_owned].shape[-1])
                            err =  Norms.L1errorVFEMvsAF(zeroVectorFunction(shtmp),
                                                         m.q['x'][0:m.mesh.subdomainMesh.nElements_owned],
                                                         m.q[('dV_u',ci)][0:m.mesh.subdomainMesh.nElements_owned],
                                                         m.q[('grad(u)',ci)][0:m.mesh.subdomainMesh.nElements_owned])

                            #for now use L1 norm of exact solution to normalize?
                            if hasAnalyticalSolution[ci]:
                                exa = Norms.L1errorVFEMvsAF(gradWrapper(p.analyticalSolution[ci]),
                                                            m.q['x'][0:m.mesh.subdomainMesh.nElements_owned],
                                                            m.q[('dV_u',ci)][0:m.mesh.subdomainMesh.nElements_owned],
                                                            numpy.zeros(m.q[('grad(u)',ci)][0:m.mesh.subdomainMesh.nElements_owned].shape,
                                                                          'd'),
                                                            T=tsim)
                            else:
                                exa = Norms.L1errorVFEM(mlvt.levelModelList[-1].q[('dV_u',ci)][0:mFine.mesh.subdomainMesh.nElements_owned],
                                                        gradu_dense[0:mFine.mesh.subdomainMesh.nElements_owned],
                                                        numpy.zeros(gradu_dense[0:mFine.mesh.subdomainMesh.nElements_owned].shape,
                                                                      'd'))


                            kerr = 'error_'+'u'+'_'+'TV'
                            kexa = 'exact_'+'u'+'_'+'TV'
                            self.errorData[ci][il][kerr].append(err)
                            self.errorData[ci][il][kexa].append(exa)
                            if self.flags['echo']:
                                if self.flags['echoRelativeErrors']:
                                    logEvent("""\nt= %g; %s[%d][%d]= %g; relative_error= %g;""" % (tsim,kerr,ci,il,err,err/(exa+1E-15)),level=0)
                                else:
                                    logEvent("""\nt= %g; %s[%d][%d]= %g;""" % (tsim,kerr,ci,il,err),level=0)
                            #end if
                        #end calcTV
                        ############### velocity specific calculations ###############
                        calcL2vel = ('L2' in self.flags['errorNorms'] and
                                     'velocity' in self.flags['errorQuantities'])
                        if calcL2vel:
                            err = -12345.0
                            exa = 1.0
                            if hasAnalyticalSolutionVelocity[ci]:
                                err = Norms.L2errorVFEMvsAF(p.analyticalSolutionVelocity[ci],
                                                            m.q['x'][0:m.mesh.subdomainMesh.nElements_owned],
                                                            m.q[('dV_u',ci)][0:m.mesh.subdomainMesh.nElements_owned],
                                                            m.q[('velocity',ci)][0:m.mesh.subdomainMesh.nElements_owned],
                                                            T=tsim)
                                exa = Norms.L2errorVFEMvsAF(p.analyticalSolutionVelocity[ci],
                                                            m.q['x'][0:m.mesh.subdomainMesh.nElements_owned],
                                                            m.q[('dV_u',ci)][0:m.mesh.subdomainMesh.nElements_owned],
                                                            numpy.zeros(m.q[('velocity',ci)][0:m.mesh.subdomainMesh.nElements_owned].shape,
                                                                          'd'),
                                                            T=tsim)
                            else:
                                #now try to project velocity to finer grids?
                                #mwf debug
                                #print """SimTools calcL2vel ci=%d veldense.shape=%s velproj[ci].shape=%s """ % (ci,veldense.shape,velproj[ci].shape)

                                err = Norms.L2errorVFEM(mlvt.levelModelList[-1].q[('dV_u',ci)][0:mFine.mesh.subdomainMesh.nElements_owned],
                                                        veldense[0:mFine.mesh.subdomainMesh.nElements_owned],velproj[ci][0:mFine.mesh.subdomainMesh.nElements_owned])
                                exa = Norms.L2errorVFEM(mlvt.levelModelList[-1].q[('dV_u',ci)][0:mFine.mesh.subdomainMesh.nElements_owned],
                                                        veldense[0:mFine.mesh.subdomainMesh.nElements_owned],
                                                        numpy.zeros(veldense[0:mFine.mesh.subdomainMesh.nElements_owned].shape,
                                                                      'd'))
                            kerr = 'error_'+'velocity'+'_'+'L2'
                            kexa = 'exact_'+'velocity'+'_'+'L2'
                            self.errorData[ci][il][kerr].append(err)
                            self.errorData[ci][il][kexa].append(exa)
                            if self.flags['echo']:
                                if self.flags['echoRelativeErrors']:
                                    logEvent("""\nt= %g; %s[%d][%d]= %g; relative_error= %g;""" % (tsim,kerr,ci,il,err,err/(exa+1E-15)),level=0)
                                else:
                                    logEvent("""\nt= %g; %s[%d][%d]= %g;""" % (tsim,kerr,ci,il,err),level=0)
                            #end if
                        #if calcL2vel
                        calcL1vel = ('L1' in self.flags['errorNorms'] and
                                     'velocity' in self.flags['errorQuantities'])
                        if calcL1vel:
                            err = -12345.0
                            exa = 1.0
                            if hasAnalyticalSolutionVelocity[ci]:
                                err = Norms.L1errorVFEMvsAF(p.analyticalSolutionVelocity[ci],
                                                            m.q['x'][0:m.mesh.subdomainMesh.nElements_owned],
                                                            m.q[('dV_u',ci)][0:m.mesh.subdomainMesh.nElements_owned],
                                                            m.q[('velocity',ci)][0:m.mesh.subdomainMesh.nElements_owned],
                                                            T=tsim)
                                exa = Norms.L1errorVFEMvsAF(p.analyticalSolutionVelocity[ci],
                                                            m.q['x'][0:m.mesh.subdomainMesh.nElements_owned],
                                                            m.q[('dV_u',ci)][0:m.mesh.subdomainMesh.nElements_owned],
                                                            numpy.zeros(m.q[('velocity',ci)][0:m.mesh.subdomainMesh.nElements_owned].shape,
                                                                          'd'),
                                                            T=tsim)
                            else:
                                #now try to project velocity to finer grids?
                                err = Norms.L1errorVFEM(mlvt.levelModelList[-1].q[('dV_u',ci)][0:mFine.mesh.subdomainMesh.nElements_owned],
                                                        veldense,velproj[ci][0:mFine.mesh.subdomainMesh.nElements_owned])
                                exa = Norms.L1errorVFEM(mlvt.levelModelList[-1].q[('dV_u',ci)][0:mFine.mesh.subdomainMesh.nElements_owned],
                                                        veldense[0:mFine.mesh.subdomainMesh.nElements_owned],
                                                        numpy.zeros(veldense[0:mFine.mesh.subdomainMesh.nElements_owned].shape,
                                                                      'd'))

                            kerr = 'error_'+'velocity'+'_'+'L1'
                            kexa = 'exact_'+'velocity'+'_'+'L1'
                            self.errorData[ci][il][kerr].append(err)
                            self.errorData[ci][il][kexa].append(exa)
                            if self.flags['echo']:
                                if self.flags['echoRelativeErrors']:
                                    logEvent("""\nt= %g; %s[%d][%d]= %g; relative_error= %g;""" % (tsim,kerr,ci,il,err,err/(exa+1E-15)),level=0)
                                else:
                                    logEvent("""\nt= %g; %s[%d][%d]= %g;""" % (tsim,kerr,ci,il,err),level=0)
                            #end if
                        #if calcL2vel

                        calcLIvel = ('LI' in self.flags['errorNorms'] and
                                     'velocity' in self.flags['errorQuantities'])
                        if calcLIvel:
                            err = -12345.0
                            exa = 1.0
                            if hasAnalyticalSolutionVelocity[ci]:
                                err = Norms.LIerrorVFEMvsAF(p.analyticalSolutionVelocity[ci],
                                                            m.q['x'][0:m.mesh.subdomainMesh.nElements_owned],
                                                            m.q[('dV_u',ci)][0:m.mesh.subdomainMesh.nElements_owned],
                                                            m.q[('velocity',ci)][0:m.mesh.subdomainMesh.nElements_owned],
                                                            T=tsim)
                                exa = Norms.LIerrorVFEMvsAF(p.analyticalSolutionVelocity[ci],
                                                            m.q['x'][0:m.mesh.subdomainMesh.nElements_owned],
                                                            m.q[('dV_u',ci)][0:m.mesh.subdomainMesh.nElements_owned],
                                                            numpy.zeros(m.q[('velocity',ci)][0:m.mesh.subdomainMesh.nElements_owned].shape,
                                                                          'd'),
                                                            T=tsim)

                            else:
                                #now try to project velocity to finer grids?
                                err = max(numpy.absolute(veldense.flat[:]-velproj[ci].flat[:]))
                                exa = max(numpy.absolute(veldense.flat))

                            kerr = 'error_'+'velocity'+'_'+'LI'
                            kexa = 'exact_'+'velocity'+'_'+'LI'
                            self.errorData[ci][il][kerr].append(err)
                            self.errorData[ci][il][kexa].append(exa)
                            if self.flags['echo']:
                                if self.flags['echoRelativeErrors']:
                                    logEvent("""\nt= %g; %s[%d][%d]= %g; relative_error= %g;""" % (tsim,kerr,ci,il,err,err/(exa+1E-15)),level=0)
                                else:
                                    logEvent("""\nt= %g; %s[%d][%d]= %g;""" % (tsim,kerr,ci,il,err),level=0)
                            #end if
                        #if calcLIvel

                    #end if calc
                #end for ci
            #end for il
        #end if numerical solution
        if 'localMassBalance' in self.flags['errorTypes']:
            from . import cfemIntegrals
            for ci in self.flags['components']:
                for il,m in enumerate(mlvt.levelModelList):
                    #
                    if self.conservationResidual[il] is None:
                        self.conservationResidual[il] = numpy.zeros((m.mesh.nElements_global,),'d')
                    else:
                        self.conservationResidual[il].flat[:] = 0.0
                    if self.elementResidual[il] is None:
                        self.elementResidual[il] = numpy.array(m.elementResidual[ci],'d')
                    else:
                        self.elementResidual[il].flat[:] = m.elementResidual[ci].flat[:]
                    if n.conservativeFlux is None or ci not in list(n.conservativeFlux.keys()) or 'dg' in n.conservativeFlux[ci]:#have to adjust residual appropriately for different methods
                        pass
                    else:
                        flux = -1.0*m.ebq_global[('totalFlux',ci)]
                        cfemIntegrals.updateExteriorElementBoundaryFlux(m.mesh.exteriorElementBoundariesArray,
                                                                        m.mesh.elementBoundaryElementsArray,
                                                                        m.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                        flux,
                                                                        m.ebq[('w*dS_u',ci)],
                                                                        self.elementResidual[il])
                    #removing boundary flux from
                    if n.conservativeFlux is None or ci not in list(n.conservativeFlux.keys()) or 'dg' in n.conservativeFlux[ci]:
                        cfemIntegrals.calculateConservationResidualDG(self.elementResidual[il],self.conservationResidual[il])
                    else:
                        cfemIntegrals.calculateConservationResidual(m.ebq['n'],
                                                                    m.ebq[('dS_u',ci)],
                                                                    self.elementResidual[il],
                                                                    m.ebq[('velocity',ci)],
                                                                    self.conservationResidual[il])
                    maxConsError = max(numpy.absolute(self.conservationResidual[il].flat[0:m.mesh.subdomainMesh.nElements_owned]))
                    self.errorData[ci][il]['localMassBalance'].append(maxConsError)
                    if self.flags['echo']:
                        logEvent("""\nt= %g; max_localMassBalanceError[%d][%d]= %g; """ % (tsim,ci,il,maxConsError),level=0)
                #il
            #ci
        #local mass balance
        #now look at global mass balance (could do similar thing each step too)
        if 'globalMassBalance' in self.flags['errorTypes']:
            for il,m in enumerate(mlvt.levelModelList):
                for ci in range(p.coefficients.nc):
                    if ci in self.flags['components'] and ('m',ci) in m.q:
                        globalMass = Norms.globalScalarDomainIntegral(m.q['abs(det(J))'],
                                                             m.elementQuadratureWeights[('m',ci)],
                                                             m.q[('m',ci)])
                        globErr = abs(globalMass-self.errorData[ci][il]['globalMass0'])
                        self.errorData[ci][il]['globalMassF'].append(globalMass)
                        if self.flags['echo']:
                            logEvent("""t= %g globalMassF[%d][%d] = %g globalMassDiff[%d][%d]= %g""" % \
                                         (tsim,ci,il,globalMass,ci,il,globErr),level=0)
                        #end if
                    #end if
                #end for ci
            #end for il
        #end calc globalMassBalance
        if 'globalHeavisideMassBalance' in self.flags['errorTypes']:
            for il,m in enumerate(mlvt.levelModelList):
                for ci in range(p.coefficients.nc):
                    if ci in self.flags['components'] and ('m',ci) in m.q:
                        hm = numpy.where(m.q[('m',ci)] >= 0.0,1.0,0.0)
                        globalMass = Norms.globalScalarDomainIntegral(m.q['abs(det(J))'],
                                                              m.elementQuadratureWeights[('m',ci)],
                                                              hm)
                        globErr = abs(globalMass-self.errorData[ci][il]['globalHeavisideMass0'])
                        self.errorData[ci][il]['globalHeavisideMassF'].append(globalMass)
                        if self.flags['echo']:
                            logEvent("""t= %g globalHeavisideMassF[%d][%d] = %g globalHeavisideMassDiff[%d][%d]= %g""" % \
                                         (tsim,ci,il,globalMass,ci,il,globErr),level=0)
                        #end if
                    #end if
                #end for ci
            #end for il
        #end calc globalMassBalance
    #end def

    def getScalarElementStorageKeys(self,mlvt,tsim):
        """
        simple utility to pull out keys for things in element quadrature dictionary that
        need to be stored
        """
        scalarElementStorageKeys = []
        for quant in [a for a in self.flags['storeQuantities'] if a is not None]:
            recType = quant.split(':')
            if len(recType) > 1 and recType[0] == 'q': #found element quadrature quantity
                stval = eval(recType[1])
                if (stval in mlvt.levelModelList[-1].q and
                    len(mlvt.levelModelList[-1].q[stval].shape) == 2): #found quantity and it's a scalar
                    scalarElementStorageKeys.append(stval)
        return scalarElementStorageKeys
    def getVectorElementStorageKeys(self,mlvt,tsim):
        """
        simple utility to pull out keys for things in element quadrature dictionary that
        need to be stored
        """
        vectorElementStorageKeys = []
        for quant in [a for a in self.flags['storeQuantities'] if a is not None]:
            recType = quant.split(':')
            if len(recType) > 1 and recType[0] == 'q': #found element quadrature quantity
                stval = eval(recType[1])
                if (stval in mlvt.levelModelList[-1].q and
                    len(mlvt.levelModelList[-1].q[stval].shape) == 3): #found quantity and it's a vector
                    vectorElementStorageKeys.append(stval)
        return vectorElementStorageKeys
    def getTensorElementStorageKeys(self,mlvt,tsim):
        """
        simple utility to pull out keys for things in element quadrature dictionary that
        need to be stored
        """
        tensorElementStorageKeys = []
        for quant in [a for a in self.flags['storeQuantities'] if a is not None]:
            recType = quant.split(':')
            if len(recType) > 1 and recType[0] == 'q': #found element quadrature quantity
                stval = eval(recType[1])
                if (stval in mlvt.levelModelList[-1].q and
                    len(mlvt.levelModelList[-1].q[stval].shape) == 4): #found quantity and it's a tensor
                    tensorElementStorageKeys.append(stval)
        return tensorElementStorageKeys
    def getScalarElementBoundaryStorageKeys(self,mlvt,tsim):
        """
        simple utility to pull out keys for things in element  boundary quadrature dictionary that
        need to be stored
        """
        scalarElementStorageKeys = []
        for quant in [a for a in self.flags['storeQuantities'] if a is not None]:
            recType = quant.split(':')
            if len(recType) > 1 and recType[0] == 'ebq_global': #found element quadrature quantity
                stval = eval(recType[1])
                if (stval in mlvt.levelModelList[-1].ebq_global and
                    len(mlvt.levelModelList[-1].ebq_global[stval].shape) == 2): #found quantity and it's a scalar
                    scalarElementStorageKeys.append(stval)
        return scalarElementStorageKeys
    def getVectorElementBoundaryStorageKeys(self,mlvt,tsim):
        """
        simple utility to pull out keys for things in element  boundary quadrature dictionary that
        need to be stored
        """
        vectorElementStorageKeys = []
        for quant in [a for a in self.flags['storeQuantities'] if a is not None]:
            recType = quant.split(':')
            if len(recType) > 1 and recType[0] == 'ebq_global': #found element quadrature quantity
                stval = eval(recType[1])
                if (stval in mlvt.levelModelList[-1].ebq_global and
                    len(mlvt.levelModelList[-1].ebq_global[stval].shape) == 3): #found quantity and it's a vector
                    vectorElementStorageKeys.append(stval)
        return vectorElementStorageKeys
    def getTensorElementBoundaryStorageKeys(self,mlvt,tsim):
        """
        simple utility to pull out keys for things in element  boundary quadrature dictionary that
        need to be stored
        """
        tensorElementStorageKeys = []
        for quant in [a for a in self.flags['storeQuantities'] if a is not None]:
            recType = quant.split(':')
            if len(recType) > 1 and recType[0] == 'ebq_global': #found element quadrature quantity
                stval = eval(recType[1])
                if (stval in mlvt.levelModelList[-1].ebq_global and
                    len(mlvt.levelModelList[-1].ebq_global[stval].shape) == 4): #found quantity and it's a tensor
                    tensorElementStorageKeys.append(stval)
        return tensorElementStorageKeys
    def getScalarExteriorElementBoundaryStorageKeys(self,mlvt,tsim):
        """
        simple utility to pull out keys for things in exterior element boundary  quadrature dictionary that
        need to be stored
        """
        scalarElementStorageKeys = []
        for quant in [a for a in self.flags['storeQuantities'] if a is not None]:
            recType = quant.split(':')
            if len(recType) > 1 and recType[0] == 'ebqe': #found element quadrature quantity
                stval = eval(recType[1])
                if (stval in mlvt.levelModelList[-1].ebqe and
                    len(mlvt.levelModelList[-1].ebqe[stval].shape) == 2): #found quantity and it's a scalar
                    scalarElementStorageKeys.append(stval)
        return scalarElementStorageKeys
    def getVectorExteriorElementBoundaryStorageKeys(self,mlvt,tsim):
        """
        simple utility to pull out keys for things in exterior element boundary quadrature dictionary that
        need to be stored
        """
        vectorElementStorageKeys = []
        for quant in [a for a in self.flags['storeQuantities'] if a is not None]:
            recType = quant.split(':')
            if len(recType) > 1 and recType[0] == 'ebqe': #found element quadrature quantity
                stval = eval(recType[1])
                if (stval in mlvt.levelModelList[-1].ebqe and
                    len(mlvt.levelModelList[-1].ebqe[stval].shape) == 3): #found quantity and it's a vector
                    vectorElementStorageKeys.append(stval)
        return vectorElementStorageKeys
    def getTensorExteriorElementBoundaryStorageKeys(self,mlvt,tsim):
        """
        simple utility to pull out keys for things in exterior element  boundary quadrature dictionary that
        need to be stored
        """
        tensorElementStorageKeys = []
        for quant in [a for a in self.flags['storeQuantities'] if a is not None]:
            recType = quant.split(':')
            if len(recType) > 1 and recType[0] == 'ebqe': #found element quadrature quantity
                stval = eval(recType[1])
                if (stval in mlvt.levelModelList[-1].ebqe and
                    len(mlvt.levelModelList[-1].ebqe[stval].shape) == 4): #found quantity and it's a tensor
                    tensorElementStorageKeys.append(stval)
        return tensorElementStorageKeys
    def stepPlotElementQuantitiesEnsight(self,mlvt,tsim):
        """
        sort through desired quantities in quadrature dictionaries like m, dm, to plot
          p    --- problem definition
          n    --- numerics definition
          mlvt --- multilevel vector transport that holds the quantities to measure
          tsim --- simulation time
        assumes this is the correct time to plot
        and plotOffSet is set correctly
        """
        if self.flags['plotOptions']['ensight']['on'] == False:
            return False
        p = self.pFile; n = self.nFile

        plottedSomething = False
        for quant in self.flags['plotQuantities']:
            recType = quant.split(':')
            if len(recType) > 1 and recType[0] == 'q': #found element quadrature quantity
                stval = eval(recType[1])
                if (stval in mlvt.levelModelList[-1].q and
                    len(mlvt.levelModelList[-1].q[stval].shape) == 2): #found quantity and it's a scalar
                    plottedSomething = True
                    self.plotScalarElementQuantityEnsight(stval,mlvt,tsim)
                elif (stval in mlvt.levelModelList[-1].q and
                      len(mlvt.levelModelList[-1].q[stval].shape) == 3): #found quantity and it's a vector
                    plottedSomething = True
                    self.plotVectorElementQuantityEnsight(stval,mlvt,tsim)
            elif len(recType) > 1 and recType[0] == 'ebq_global': #found global element boundary quantity
                stval = eval(recType[1])
                if (stval in mlvt.levelModelList[-1].ebq_global and
                    len(mlvt.levelModelList[-1].ebq_global[stval].shape) == 2): #found quantity and it's a scalar
                    plottedSomething = True
                    self.plotScalarGlobalElementBoundaryQuantityEnsight(stval,mlvt,tsim)
                elif (stval in mlvt.levelModelList[-1].ebq_global and
                      len(mlvt.levelModelList[-1].ebq_global[stval].shape) == 3): #found quantity and its a vector
                    self.plotVectorGlobalElementBoundaryQuantityEnsight(stval,mlvt,tsim)
                    plottedSomething = False
                    #ebq_global vector
                #has key
            #ebq_global
        #quantities
    def stepStoreQuantities(self,mlvt,tsim):
        """
        shelve quantities for a given time instance, if self.storeHeavyData == True
        soon to be deprecated and will use Archiver tools instead


        Quadrature dictionary quantities are shelved in a dictionary
        whose key is the corresponding quadrature dictionary name.

        The stored dictionary has fields 'x' 't' and 'vals' which hold
        the quadrature points (assumed static for now) as well as
        lists of the desired quantities at the requested time
        levels. The quantities to be stored are specified in
        flags['storeQuantities'] in the format

        "q:(%s,%d)",ebq:(%s,%d) etc e.g., "q:('u',0)" will store component 0 solution values from q

        dataStorage['q']['x']= [[0.,0.,0.],[...],...] element quadrature points for 1st time called
        dataStorage['q']['t']= [0.0,0.1,...,1.0]
        dataStorage['q']['vals'] = [subset of q at t=0.0, subset of q at t=0.1, ...]

        if multilevelModel is set, stores basically everything

        Input
          p    --- problem definition
          n    --- numerics definition
          mlvt --- multilevel vector transport that holds the quantities to measure
          tsim --- simulation time
        assumes this is the correct time to store

        TODO: add option for storage directory
        """
        assert self.dataStorage is not None, "dataStorage None storeTimes= %s " % self.flags['storeTimes']
        if self.storeHeavyData == False:
            return
        p = self.pFile; n = self.nFile
        #mwf debug
        #print """SimTools entering stepStore t=%s dataFile=%s """ % (tsim,self.flags['dataFile'])
        m = mlvt.levelModelList[-1]
        q = {}; ebq = {}; ebq_global = {};
        solutionData = {};

        for quant in self.flags['storeQuantities']:
            recType = quant.split(':')
            if len(recType) > 1:
                quadDict = recType[0]
                stval    = eval(recType[1])
                if recType[0] == 'q' and stval in m.q:
                    q[stval]=m.q[stval]
                    #if not q.has_key('x'):
                    #    q['x']= m.q['x']
                elif recType[0] == 'ebq' and stval in m.ebq:
                    ebq[stval] = m.ebq[stval]
                    #if not ebq.has_key('x'):
                    #    ebq['x']=m.ebq['x']
                elif recType[0] == 'ebq_global' and stval in m.ebq_global:
                    ebq_global = m.ebq_global[stval]
                    #if not ebq_global.has_key('x'):
                    #    ebq_global['x']=m.ebq_global['x']
            elif recType == 'u_dof':
                for ci in self.flags['components']:
                    solutionData[ci]={}
                    for il,im in enumerate(mlvt.levelModelList):
                        solutionData[ci][il] = {}
                        solutionData[ci][il]['u_dof'] = im.u[ci].dof
                        solutionData[ci][il]['l2g'] = im.u[ci].femSpace.dofMap.l2g
            #found quad dict entry for storage
        #end quantities
        for d in ['q','ebq','ebq_global','solutionData']:
            dval = eval(d)
            #mwf debug
            #print """SimTools stepStore tsim=%s d=%s len(dval)=%d """ % (tsim,d,len(dval))
            if len(dval) > 0:
                if d in self.dataStorage:
                    dtmp = self.dataStorage[d]
                    dtmp['t'].append(tsim)
                    dtmp['vals'].append(dval)
                    self.dataStorage[d] = dtmp
                    #mwf debug
                    #print """SimTools data Has Key %s """ % d
                else:
                    #mwf debug
                    #print """SimTools data did not Have Key %s """ %d
                    self.dataStorage[d]={}
                    dtmp = {'t':[tsim],'vals':[dval],'x':getattr(m,d)['x']}
                    self.dataStorage[d] = dtmp
                #not already stored
            #end something to store
        #end loop through things to store
    #def
    def computeNodalQuadratureInfo(self,mlvt,t):
        """
        if need values of quantities at mesh nodes and don't have them already, use this
        only compute values on finest mesh for now
        """
        from . import Quadrature
        self.nodalQuadratureInfo = {}
        vt = mlvt.levelModelList[-1]
        nd = vt.nSpace_global; nq = nd+1 ; ne = vt.mesh.nElements_global
        quad = Quadrature.SimplexLobattoQuadrature(nd,1)

        self.nodalQuadratureInfo['elementQuadraturePoints'] = numpy.array(quad.points,'d')
        self.nodalQuadratureInfo['x'] = numpy.zeros((ne,nq,3),'d')
        self.nodalQuadratureInfo['J'] = numpy.zeros((ne,nq,nd,nd),'d')
        self.nodalQuadratureInfo['inverse(J)'] = numpy.zeros((ne,nq,nd,nd),'d')
        self.nodalQuadratureInfo['det(J)'] = numpy.zeros((ne,nq),'d')

        vt.u[0].femSpace.elementMaps.getValues(self.nodalQuadratureInfo['elementQuadraturePoints'],
                                               self.nodalQuadratureInfo['x'])
        vt.u[0].femSpace.elementMaps.getJacobianValues(self.nodalQuadratureInfo['elementQuadraturePoints'],
                                                       self.nodalQuadratureInfo['J'],
                                                       self.nodalQuadratureInfo['inverse(J)'],
                                                       self.nodalQuadratureInfo['det(J)'])
        self.nodalQuadratureInfo['abs(det(J))']=numpy.absolute(self.nodalQuadratureInfo['det(J)'])
        for cj in self.flags['components']:
            self.nodalQuadratureInfo[('v',cj)] = numpy.zeros((ne,nq,nd+1),'d')
            self.nodalQuadratureInfo[('grad(v)',cj)] = numpy.zeros((ne,nq,nd+1,nd),'d')
            vt.u[cj].femSpace.getBasisValues(self.nodalQuadratureInfo['elementQuadraturePoints'],
                                             self.nodalQuadratureInfo[('v',cj)])
            vt.u[cj].femSpace.getBasisGradientValues(self.nodalQuadratureInfo['elementQuadraturePoints'],
                                                     self.nodalQuadratureInfo['inverse(J)'],
                                                     self.nodalQuadratureInfo[('grad(v)',cj)])
        #cj

        #wasteful
        for key in list(vt.q.keys()):
            if key not in list(self.nodalQuadratureInfo.keys()):
                tmp = list(vt.q[key].shape)
                if len(tmp) > 1:
                    tmp[1] = nd+1
                self.nodalQuadratureInfo[key] = numpy.zeros(tuple(tmp),'d')
            #key not in already
        #keys
        #mwf debug
        #for key in self.nodalQuadratureInfo.keys():
        #    print """SimTools nodalQuadrature Dict key= %s, shape=%s \n""" % (key,self.nodalQuadratureInfo[key].shape)
        #print """SimTools nodalQuadrature x = %s """ % self.nodalQuadratureInfo['x']

        #really wasteful, but need to get spatial quantities at new points
        import copy
        self.nodalQuadratureInfo['coefficients'] = copy.deepcopy(vt.coefficients)
        self.nodalQuadratureInfo['coefficients'].initializeElementQuadrature(t,self.nodalQuadratureInfo)
    #
    def stepPlotEnsight(self,mlvt,tsim):
        """
        plot solution and solution 'velocity' for ensight/paraview
          p    --- problem definition
          n    --- numerics definition
          mlvt --- multilevel vector transport that holds the quantities to measure
          tsim --- simulation time
        assumes this is the correct time to plot
        and plotOffSet is set correctly

        assumes that initial case,sos, geo files set somewhere else
        """
        p = self.pFile; n = self.nFile

        if self.flags['plotOptions']['ensight']['on'] == False:
            return False
        mFinest = mlvt.levelModelList[-1]
        for ci in self.flags['components']:
            mFinest.u[ci].femSpace.writeFunctionEnsight(mFinest.u[ci],
                                                        self.flags['plotOptions']['ensight']['caseFileName'],
                                                        append=True,
                                                        firstVariable=False)
        #ci
        if mFinest.coefficients.vectorComponents is not None:
            if len(mFinest.coefficients.vectorComponents) == 2:
                vcomp = [mFinest.coefficients.vectorComponents[0],
                         mFinest.coefficients.vectorComponents[1]]
                mFinest.u[vcomp[0]].femSpace.writeE2VectorFunctionEnsight(mFinest.u[vcomp[0]],
                                                                          mFinest.u[vcomp[1]],
                                                                          self.flags['plotOptions']['ensight']['caseFileName'],
                                                                          nOutput=mFinest.u[vcomp[0]].femSpace.nOutput-1,
                                                                          append=True,
                                                                          firstVariable=False)
            elif len(mFinest.coefficients.vectorComponents) == 3:
                vcomp = [mFinest.coefficients.vectorComponents[0],
                         mFinest.coefficients.vectorComponents[1],
                         mFinest.coefficients.vectorComponents[2]]
                mFinest.u[vcomp[0]].femSpace.writeE3VectorFunctionEnsight(mFinest.u[vcomp[0]],
                                                                          mFinest.u[vcomp[1]],
                                                                          mFinest.u[vcomp[2]],
                                                                          self.flags['plotOptions']['ensight']['caseFileName'],
                                                                          nOutput=mFinest.u[vcomp[0]].femSpace.nOutput-1,
                                                                          append=True,
                                                                          firstVariable=False)
            #3d
        #vector components
        self.stepPlotCalled['ensight'] = True
        return False
    #stepPlotEnsight
    def writeEnsightMeshForElementQuantities(self,filename,mlvt,tsim=0.0):
        caseOut=open(filename+'.case','a')
        caseOut.write('measured: '+filename+'q.geo\n')
        caseOut.close()
        meshOut=open(filename+'q.geo','w')
        meshOut.write('Element quadrature\n')
        meshOut.write('particle coordinates\n')
        meshOut.write('%8i\n' % (mlvt.levelModelList[-1].mesh.nElements_global*mlvt.levelModelList[-1].nQuadraturePoints_element,))
        pN=1
        for eN in range(mlvt.levelModelList[-1].mesh.nElements_global):
            for k in range(mlvt.levelModelList[-1].nQuadraturePoints_element):
                meshOut.write('%8i%12.5E%12.5E%12.5E\n' % (pN,
                                                           mlvt.levelModelList[-1].q['x'][eN,k,0],
                                                           mlvt.levelModelList[-1].q['x'][eN,k,1],
                                                           mlvt.levelModelList[-1].q['x'][eN,k,2]))
                pN+=1
        meshOut.close()

    def writeEnsightMeshForElementBoundaryQuantities(self,filename,mlvt,tsim=0.0):
        caseOut=open(filename+'.case','a')
        caseOut.write('measured: '+filename+'ebq.geo\n')
        caseOut.close()
        meshOut=open(filename+'ebq.geo','w')
        meshOut.write('Element boundary quadrature\n')
        meshOut.write('particle coordinates\n')
        meshOut.write('%8i\n' % (mlvt.levelModelList[-1].mesh.nElementBoundaries_global*mlvt.levelModelList[-1].nElementBoundaryQuadraturePoints_elementBoundary,))
        pN=1
        for ebN in range(mlvt.levelModelList[-1].mesh.nElementBoundaries_global):
            for k in range(mlvt.levelModelList[-1].nElementBoundaryQuadraturePoints_elementBoundary):
                meshOut.write('%8i%12.5E%12.5E%12.5E\n' % (pN,
                                                           mlvt.levelModelList[-1].ebq_global['x'][ebN,k,0],
                                                           mlvt.levelModelList[-1].ebq_global['x'][ebN,k,1],
                                                           mlvt.levelModelList[-1].ebq_global['x'][ebN,k,2]))
                pN+=1
        meshOut.close()

    def writeScalarElementFunctionHeaderEnsight(self,ckey,filename,append=False,firstVariable=True,case_filename=None):
        if case_filename is None:
            case_filename = filename
        if not append:
            caseOut=open(case_filename+'.case','a')
            if firstVariable == True:
                caseOut.write('VARIABLE\n')
            line="scalar per measured node: %s%s %s%s%s.scl****\n" % (ckey[0],ckey[1],filename,
                                                                      ckey[0],ckey[1]);
            caseOut.write(line)
            caseOut.close()
        #
    def writeVectorElementFunctionHeaderEnsight(self,ckey,filename,append=False,firstVariable=True,case_filename=None):
        if case_filename is None:
            case_filename = filename
        if not append:
            caseOut=open(case_filename+'.case','a')
            if firstVariable == True:
                caseOut.write('VARIABLE\n')
            line= "vector per measured node: %s%s %s%s%s.vec****\n" % (ckey[0],ckey[1],filename,
                                                                       ckey[0],ckey[1]);
            caseOut.write(line)
            caseOut.close()
        #

    def plotScalarElementQuantityEnsight(self,ckey,mlvt,tsim):
        """
        Ensight plotting routine to look at scalar quantity stored in element quad dictionary q
          ckey --- what should be plotted
          p    --- problem definition
          n    --- numerics definition
          mlvt --- multilevel vector transport that holds the quantities to measure
          tsim --- simulation time
        assumes this is the correct time to plot
        and plotOffSet is set correctly

        assumes that initial case,sos, geo files set somewhere else

        """
        p = self.pFile; n = self.nFile

        if self.flags['plotOptions']['ensight']['on'] == False:
            return False
        case_filename = self.flags['plotOptions']['ensight']['caseFileName']
        filename_ci = "%s%s%s.scl%4.4i" %(case_filename,ckey[0],ckey[1],mlvt.levelModelList[-1].u[ckey[1]].femSpace.nOutput-1)
        uOut = open(filename_ci,'w')
        uOut.write("%s%s\n" % (ckey[0],ckey[1]))
        n=0
        #slow
        for eN in range(mlvt.levelModelList[-1].mesh.nElements_global):
            for k in range(mlvt.levelModelList[-1].nQuadraturePoints_element):
                #has to be 6 per line
                uOut.write('%12.5e' % mlvt.levelModelList[-1].q[ckey][eN,k])
                n+=1
                if n%6==0:
                    uOut.write('\n')
        uOut.write('\n')
        uOut.close()

        self.stepPlotCalled['ensightElementQuantities']= True
    def plotVectorElementQuantityEnsight(self,ckey,mlvt,tsim,scaleOutput=None):
        """
        Ensight plotting routine to look at vector quantity stored in element quad dictionary q
          ckey --- what should be plotted
          p    --- problem definition
          n    --- numerics definition
          mlvt --- multilevel vector transport that holds the quantities to measure
          tsim --- simulation time
        assumes this is the correct time to plot
        and plotOffSet is set correctly

        assumes that initial case,sos, geo files set somewhere else
        TODO : check format for 3d
        """
        p = self.pFile; n = self.nFile

        if self.flags['plotOptions']['ensight']['on'] == False:
            return False
        case_filename = self.flags['plotOptions']['ensight']['caseFileName']
        filename_ci = "%s%s%s.vec%4.4i" %(case_filename,ckey[0],ckey[1],mlvt.levelModelList[-1].u[ckey[1]].femSpace.nOutput-1)
        uOut = open(filename_ci,'w')
        uOut.write("%s%s\n" % (ckey[0],ckey[1]))
        n=0
        vmax=1. #no scaling by default
        if scaleOutput == 'maxComponent':
            vmax =max(mlvt.levelModelList[-1].q[ckey].flat[:])+1.0e-8
            print("WARNING SimTools Ensight_q_%s: Scaling velocity for output by %s" % (ckey,vmax))
        for eN in range(mlvt.levelModelList[-1].mesh.nElements_global):
            for k in range(mlvt.levelModelList[-1].nQuadraturePoints_element):
                for i in range(mlvt.levelModelList[-1].q[ckey].shape[-1]):
                    uOut.write('%12.5e' % (mlvt.levelModelList[-1].q[ckey][eN,k,i]/vmax))
                for i in range(mlvt.levelModelList[-1].q[ckey].shape[-1],3):
                    uOut.write('%12.5e' % (0.0))
                if n%2==1:
                    uOut.write('\n')
                n+=1
        if n%2==1:
            uOut.write('\n')
        uOut.close()

        self.stepPlotCalled['ensightElementQuantities']= True

    def plotScalarGlobalElementBoundaryQuantityEnsight(self,ckey,mlvt,tsim):
        """
        Ensight plotting routine to look at scalar quantity stored in element quad dictionary ebq_global
          ckey --- what should be plotted
          p    --- problem definition
          n    --- numerics definition
          mlvt --- multilevel vector transport that holds the quantities to measure
          tsim --- simulation time
        assumes this is the correct time to plot
        and plotOffSet is set correctly

        assumes that initial case,sos, geo files set somewhere else

        """
        p = self.pFile; n = self.nFile

        if self.flags['plotOptions']['ensight']['on'] == False:
            return False
        if (self.plottingQuadratureValuesForEnsight['elements'] == True or
            self.plottingQuadratureValuesForEnsight['elementBoundaries'] == False):
            return False
        case_filename = self.flags['plotOptions']['ensight']['caseFileName']
        filename_ci = "%s%s%s.scl%4.4i" %(case_filename,ckey[0],ckey[1],mlvt.levelModelList[-1].u[ckey[1]].femSpace.nOutput-1)
        uOut = open(filename_ci,'w')
        uOut.write("%s%s\n" % (ckey[0],ckey[1]))
        n=0
        for ebN in range(mlvt.levelModelList[-1].mesh.nElementBoundaries_global):
            for k in range(mlvt.levelModelList[-1].nElementBoundaryQuadraturePoints_elementBoundary):
                #has to be 6 per line
                uOut.write('%12.5e' % mlvt.levelModelList[-1].ebq_global[ckey][ebN,k])
                n+=1
                if n%6==0:
                    uOut.write('\n')
        uOut.write('\n')
        uOut.close()

        self.stepPlotCalled['ensightElementQuantities']= True
    def plotVectorGlobalElementBoundaryQuantityEnsight(self,ckey,mlvt,tsim,scaleOutput=None):
        """
        Ensight plotting routine to look at vector quantity stored in element quad dictionary ebq_global
          ckey --- what should be plotted
          p    --- problem definition
          n    --- numerics definition
          mlvt --- multilevel vector transport that holds the quantities to measure
          tsim --- simulation time
        assumes this is the correct time to plot
        and plotOffSet is set correctly

        assumes that initial case,sos, geo files set somewhere else
        TODO : check format for 3d
        """
        if self.flags['plotOptions']['ensight']['on'] == False:
            return False
        if (self.plottingQuadratureValuesForEnsight['elements'] == True or
            self.plottingQuadratureValuesForEnsight['elementBoundaries'] == False):
            return False
        p = self.pFile; n = self.nFile

        case_filename = self.flags['plotOptions']['ensight']['caseFileName']
        filename_ci = "%s%s%s.vec%4.4i" %(case_filename,ckey[0],ckey[1],mlvt.levelModelList[-1].u[ckey[1]].femSpace.nOutput-1)
        uOut = open(filename_ci,'w')
        uOut.write("%s%s\n" % (ckey[0],ckey[1]))
        n=0
        vmax=1. #no scaling by default
        if scaleOutput == 'maxComponent':
            vmax =max(mlvt.levelModelList[-1].q[ckey].flat[:])+1.0e-8
            print("WARNING SimTools Ensight_ebq_global_%s: Scaling velocity for output by %s" % (ckey,vmax))
        for ebN in range(mlvt.levelModelList[-1].mesh.nElementBoundaries_global):
            for k in range(mlvt.levelModelList[-1].nElementBoundaryQuadraturePoints_elementBoundary):
                for i in range(mlvt.levelModelList[-1].ebq_global[ckey].shape[-1]):
                    uOut.write('%12.5e' % (mlvt.levelModelList[-1].ebq_global[ckey][ebN,k,i]/vmax))
                for i in range(mlvt.levelModelList[-1].ebq_global[ckey].shape[-1],3):
                    uOut.write('%12.5e' % (0.0))
                if n%2==1:
                    uOut.write('\n')
                n+=1
        if n%2==1:
            uOut.write('\n')
        uOut.close()

        self.stepPlotCalled['ensightElementQuantities']= True

#end SimulationProcessor

########################################################################
#project solutions to fine grid for computing error
########################################################################

def projectToFinestLevel(mlTransport,level,tsim=0.0,verbose=0):
    """use multilevel transport prolongation to get fine grid information
    starting at level.

    returns quadrature dictionary of projected values on fine grid

    """
#TODO
#       appears broken (error values not consistent) 1/13/10
#       set uproj to be size of level down to mfine rather than full hiearachy
    import numpy
    nLevels = len(mlTransport.uList)
    assert 0 <= level and level <= nLevels, "projectToFinestLevel range= [0,%d]" % nLevels-1

    coefficients = mlTransport.levelModelList[-1].coefficients
    uqprojFine = [numpy.zeros(mlTransport.levelModelList[-1].q[('u',ci)].shape,'d')
                  for ci in range(coefficients.nc)]
    graduqProjFine = [numpy.zeros(mlTransport.levelModelList[-1].q[('grad(u)',ci)].shape,'d')
                      for ci in range(coefficients.nc)]
    uproj  = [[FemTools.FiniteElementFunction(m.u[ci].femSpace) for m in mlTransport.levelModelList]
              for ci in range(coefficients.nc)]
    mFine = mlTransport.levelModelList[-1]
    for ci in range(coefficients.nc):
        for l in range(nLevels):
            uproj[ci][l].dof[:] = 0.0
        m = mlTransport.levelModelList[level]
        uproj[ci][level].dof[:] = m.u[ci].dof
        if level < nLevels-1:
            for lf in range(level,nLevels-1):
                mlTransport.meshTransfers.prolong_bcListDict[ci][lf+1].matvec(uproj[ci][lf].dof,
                                                                              uproj[ci][lf+1].dof)
                #load Dirichlet conditions in
                for dofN,g in mlTransport.levelModelList[lf+1].dirichletConditions[ci].DOFBoundaryConditionsDict.items():
                    uproj[ci][lf+1].dof[dofN] = g(mlTransport.levelModelList[lf+1].dirichletConditions[ci].DOFBoundaryPointDict[dofN],tsim)
                #dirichlet conditions
            #lf up to fine
            uproj[ci][-1].getValues(mFine.q['v',ci],uqprojFine[ci])
            uproj[ci][-1].getGradientValues(mFine.q['grad(v)',ci],graduqProjFine[ci])
            #mwf debug
            #from proteusGraphical import vtkViewers
            #import pdb
            #pdb.set_trace()
        else: #already on fine
            uqprojFine[ci].flat[:] = mFine.q[('u',ci)].flat[:]
            graduqProjFine[ci].flat[:] = mFine.q[('grad(u)',ci)].flat[:]
        #else
    #ci
    if verbose > 2:
        from . import Viewers
        if 'viewerType' in dir(Viewers) and Viewers.viewerType == 'gnuplot' and mFine.nSpace_global == 2:
            for ci in range(coefficients.nc):
                for eN in range(uqprojFine[ci].shape[0]):
                    for k in range(uqprojFine[ci].shape[1]):
                        Viewers.datFile.write("%12.5e %12.5e %12.5e \n" % (mFine.q['x'][eN,k,0],
                                                                           mFine.q['x'][eN,k,1],
                                                                           uqprojFine[ci][eN,k]))
                    #
                #eN
                Viewers.datFile.write("\n \n#end uqproj ci=%d level=%d \n" % (ci,level))
                ggrid = (3-1)*(2**nLevels)
                title="uqproj ci=%d level=%d " % (ci,level)
                cmd = "set dgrid3d %d,%d,16; set contour base; set term x11 %i; splot \'%s\' index %i with lines title \"%s\" \n" % (ggrid,ggrid,Viewers.windowNumber,
                                                                                                                                 Viewers.datFilename,
                                                                                                                                 Viewers.plotNumber,
                                                                                                                                 title)
                Viewers.cmdFile.write(cmd)
                Viewers.viewerPipe.write(cmd)
                Viewers.newPlot()
                Viewers.newWindow()
                input('press return to continue')
            #end ci
        #end viewer typ
    return uqprojFine,graduqProjFine
########################################################################
#try computing solution values directly on fine grid using coarse grid
#for nonconforming approximations
########################################################################
def generateParentInfo(mlMesh):
    """
    get array P[l,e] = e_c, where element e_c is the parent of element e
        P[0,:] = -1
    """
    #mwf now use new interface
    P = mlMesh.elementParentsArrayList
    import numpy
    P = {}
    nLevels = len(mlMesh.meshList)
    P[0] = numpy.ones((mlMesh.meshList[0].nElements_global,),'i')
    P[0][:] = -1
    for l in range(1,nLevels):
        P[l] = mlMesh.elementParentsArrayList[l]
    return P

def getIntegrationPointsOnCoarseGrid(xf,lf,P,lc):
    """
    given array of points N^f_e x n_q x 3 on level lf
    generate dictionary on coarse
    grid that's N^c_e x n_q^c(e_c) x 3 and holds the integration points
    assigned to the correct coarse grid element. If using uniform refinement
    should have the same number of integration points per coarse grid element
    but the number depends on the level of refinement eg  n_q^4(lf-lc).
    In general the number won't be the same for nonuniform refinement
    """
    nEf = len(P[lf])
    assert nEf == xf.shape[0], "fine grid mismatch lf=%d nEf=%d xf.shape=%s " % (lf,nEf,xf.shape)
    nq  = xf.shape[1]
    nEc = len(P[lc])
    ldiff = lf-lc
    #mwf debug
    #print """testStuff generateCoarseGridX lf=%d nEf=%d xf.shape=%s lc=%d nEc=%d ldiff=%d """ % (lf,nEf,xf.shape,lc,nEc,
    #                                                                                             ldiff)

    xc = {}
    for ec in range(nEc):
        xc[ec] = {}

    for ef in range(nEf):
        ec = ef
        for l in range(ldiff):
            ep = P[lf-l][ec]
            #mwf debug
            #print "genCoarseGridX ef=%d l=%d ec=%d ep=%d " % (ef,l,ec,ep)
            ec = ep

        for iq in range(nq):
            xc[ec][(ef,iq)] = xf[ef,iq,:]
        #iq
    #ef
    #mwf debug
    #print """getIntPointsOnCoarseGrid lf=%s lc=%s \n xf=%s""" % (lf,lc,xf)
    #for ec in range(nEc):
    #    print """xc[%s] len=%d x=%s """ % (ec,len(xc[ec]),xc[ec])


    return xc
def projectVelocityToFinestLevelNC(mlTransport,level,ci=0,tsim=0.0,verbose=0):
    """
    use brute force evaluation to get coarse grid quantities on fine grid
    starting at level.

    returns quadrature dictionary of projected values on fine grid
    """
    import numpy
    nLevels = len(mlTransport.levelModelList)

    assert 0 <= level and level < nLevels, "projectVelocityToFinestLevelNC range= [0,%d]" % nLevels-1

    mFine  = mlTransport.levelModelList[-1]
    mCoarse= mlTransport.levelModelList[level]
    if mCoarse.velocityPostProcessor is None:
        return None

    P = generateParentInfo(mlTransport.mlMeshSave)
    nEf = mFine.q['x'].shape[0]; nqf = mFine.q['x'].shape[1]
    lf = nLevels-1 #index into list
    ldiff = lf-level
    if level == lf:
        velciprojFine = mFine.q[('velocity',ci)]
        xArray        = mFine.q['x']
    else:
        #mwf debug
        #import pdb
        #pdb.set_trace()
        xc= getIntegrationPointsOnCoarseGrid(mFine.q['x'],lf,P,level)
        nEc = len(xc)
        #will get extra padding for some coarse cells with nonuniform refinement
        #but these points should get skipped in error calculation
        nqc = max([len(xc[i]) for i in range(nEc)])
        xArray  = numpy.zeros((nEc,nqc,3),'d')
        for ec in range(nEc):
            iqc = 0
            for k,x in xc[ec].items():
                #mwf debug
                #print "ec=%d iqc=%d x=%s " % (ec,iqc,x)
                xArray[ec,iqc,:] = x
                iqc += 1
            #end x
            #now padd xArray using last value if number of integration points is less than max
            #over domain
            if iqc < nqc:
                for iqp in range(iqc,nqc):
                    xArray[ec,iqp,:]=xArray[ec,iqc-1,:]
        #end ec
        velciprojFine = numpy.zeros(mFine.q[('velocity',ci)].shape,'d')
        if mCoarse.velocityPostProcessor.postProcessingTypes[ci] == 'point-eval':
            #assume constant solution/potential gradient over coarse grid
            print("WARNING projectVelocityToFinestLevelNC type= point-eval assuming constant potential on coarse grid")

            for ef in range(nEf):
                ec = ef
                for l in range(ldiff):
                    ep = P[lf-l][ec]
                    ec = ep
                for iq in range(nqf):
                    if ('a',ci,ci) in mFine.q:
                        velciprojFine[ef,iq,:] = -numpy.dot(mFine.q[('a',ci,ci)][ef,iq,:,:],
                                                              mCoarse.q[('grad(phi)',ci)][ec,0,:])
                    else:
                        velciprojFine[ef,iq,:] = 0.0
                    if ('f',ci) in mFine.q:
                        velciprojFine[ef,iq,:]  += mFine.q[('f',ci)][ef,iq,:]
                #iq
            #ef
        else:
            velci0 = mCoarse.velocityPostProcessor.evaluateElementVelocityField(xArray,ci)
            for ec in range(nEc):
                iqc = 0
                for k,x in xc[ec].items():
                    ef = k[0]; iqf = k[1]
                    velciprojFine[ef,iqf,:] = velci0[ec,iqc,:]
                    iqc += 1
                #x
            #ec

        #postprocessing type
    #else on level
    if verbose > 2:
        print("""velocityProjNC \n xArray=%s velciprojFine= %s \n""" % (xArray,
                                                                        velciprojFine))
    if verbose > 1:
        from . import Viewers
        if 'viewerType' in dir(Viewers) and Viewers.viewerType == 'gnuplot' and mFine.nSpace_global == 2:
            max_u=max(numpy.absolute(numpy.take(velciprojFine,[0],2).flat))
            max_v=max(numpy.absolute(numpy.take(velciprojFine,[1],2).flat))
            L = min((max(mFine.mesh.nodeArray[:,0]),max(mFine.mesh.nodeArray[:,1])))
            scale =10.0*max((max_u,max_v,1.0e-16))/L
            for eN in range(mFine.q['x'].shape[0]):
                for iq in range(mFine.q['x'].shape[1]):
                    x = mFine.q['x'][eN,iq,:]
                    v = velciprojFine[eN,iq,:]
                    Viewers.datFile.write("%12.5e %12.5e %12.5e %12.5e \n" % (x[0],x[1],
                                                                              v[0]/scale,
                                                                              v[1]/scale))
            Viewers.datFile.write("\n \n#end velciproj ci=%d level=%d" % (ci,level))
            title = "velciproj ci=%d level=%d " % (ci,level)
            cmd = "set term x11 %i; plot \'%s\' index %i with vectors title \"%s\" \n" % (Viewers.windowNumber,
                                                                                                  Viewers.datFilename,
                                                                                                  Viewers.plotNumber,
                                                                                                  title)
            Viewers.cmdFile.write(cmd)
            Viewers.viewerPipe.write(cmd)
            Viewers.newPlot()
            Viewers.newWindow()
            input('press return to continue')

            #now try just coarse grid velocity
            max_u=max(numpy.absolute(numpy.take(mCoarse.q[('velocity',ci)],[0],2).flat))
            max_v=max(numpy.absolute(numpy.take(mCoarse.q[('velocity',ci)],[1],2).flat))
            L = min((max(mFine.mesh.nodeArray[:,0]),max(mFine.mesh.nodeArray[:,1])))
            scale =10.0*max((max_u,max_v,1.0e-16))/L
            for eN in range(mCoarse.q['x'].shape[0]):
                for iq in range(mCoarse.q['x'].shape[1]):
                    x = mCoarse.q['x'][eN,iq,:]
                    v = mCoarse.q[('velocity',ci)][eN,iq,:]
                    Viewers.datFile.write("%12.5e %12.5e %12.5e %12.5e \n" % (x[0],x[1],
                                                                              v[0]/scale,
                                                                              v[1]/scale))
            Viewers.datFile.write("\n \n#end coarse velocity ci=%d level=%d" % (ci,level))
            title = "coarse velocity ci=%d level=%d " % (ci,level)
            cmd = "set term x11 %i; plot \'%s\' index %i with vectors title \"%s\" \n" % (Viewers.windowNumber,
                                                                                          Viewers.datFilename,
                                                                                          Viewers.plotNumber,
                                                                                          title)
            Viewers.cmdFile.write(cmd)
            Viewers.viewerPipe.write(cmd)
            Viewers.newPlot()
            Viewers.newWindow()
            input('press return to continue')

        #end gnuplot
    #end verbose
    return velciprojFine
"""
Some simple modules for doing runtime visualization

.. todo::

   Added support for image cache and some global data structure that provides a map into the image cache
   Clean up Viewers.py

.. inheritance-diagram:: proteus.Viewers
   :parts: 1
"""
import subprocess
import numpy

cmdFile=None
datFile=None
datFilename=None
viewerPipe=None
viewerType=None
plotNumber=None
windowNumber=None
meshDataStructuresWritten=None


def gnuplotOn(problemName):
    global viewerPipe,cmdFile,datFilename,datFile,viewerType,plotNumber,windowNumber
    viewerPipe=subprocess.Popen('gnuplot',shell=True,bufsize=1,stdin=subprocess.PIPE).stdin
    cmdFile=open(problemName+'_gnuplot.cmd','w',1)
    datFilename = problemName+'_gnuplot.dat'
    datFile=open(datFilename,'w',1)
    viewerType='gnuplot'
    plotNumber=0
    windowNumber=0

def matlabOn(problemName,silent=True):
    global viewerPipe,cmdFile,datFilename,datFile,viewerType,plotNumber,windowNumber
    #mwf add for handling mesh data structures
    #mwf debug
    if silent == True:
        viewerPipe=open('/dev/null','w',1)
    else:
        #        viewerPipe=subprocess.Popen('matlab -nosplash -nodesktop -nojvm',shell=True,bufsize=1,
        #                                    stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE).stdin
        dummyFile = open('/dev/null','w')
        viewerPipe=subprocess.Popen('matlab -nosplash -nodesktop -nojvm',shell=True,bufsize=1,
                                    stdin=subprocess.PIPE,stdout=dummyFile,stderr=dummyFile).stdin
    cmdFile=open(problemName+'.m','w',1)
    viewerType='matlab'
    plotNumber=0
    windowNumber=0
    meshDataStructuresWritten = False

def vtkOn(problemName,silent=True):
    global viewerPipe,cmdFile,datFilename,datFile,viewerType,plotNumber,windowNumber
    #viewerPipe=subprocess.Popen('gnuplot',shell=True,bufsize=1,stdin=subprocess.PIPE).stdin
    viewerPipe=open('/dev/null','w',1)
    cmdFile=open(problemName+'_vtk_dummy.cmd','w',1)
    datFilename = problemName+'_vtk_dummy.dat'
    datFile=open(datFilename,'w',1)
    import os
    from proteusGraphical import vtkViewers
    #vtkViewers.g.ImageFolderPath = os.path.abspath('.')+"/results/"
    vtkViewers.g.ImageFolderPath = os.path.abspath('.')+"/tmp/"+problemName+"/"
    if not os.path.exists(vtkViewers.g.ImageFolderPath):
        os.makedirs(vtkViewers.g.ImageFolderPath)
    viewerType='vtk'
    plotNumber=0
    windowNumber=0
    return vtkViewers

def viewerOn(problemName,viewer):
    if viewer == 'gnuplot':
        return gnuplotOn(problemName)
    if viewer == 'matlab':
        return matlabOn(problemName)
    if viewer == 'vtk':
        return vtkOn(problemName)

def newPlot():
    global plotNumber
    #print "plot number",plotNumber
    plotNumber +=1

def newWindow():
    global windowNumber
    #print "window",windowNumber
    windowNumber +=1

class V_base(object):
    def __init__(self,
                 p=None,
                 n=None,
                 s=None):
        if p is None:
            from . import default_p as p
        if n is None:
            from . import default_n as n
        if s is None:
            from . import default_s as s
        global cmdFile,datFile,datFilename,viewerPipe,viewerType,plotNumber,windowNumber,meshDataStructuresWritten
        self.cmdFile=cmdFile
        self.datFile=datFile
        self.datFilename = datFilename
        self.viewerPipe=viewerPipe
        self.viewerType=viewerType
        self.meshDataStructuresWritten=meshDataStructuresWritten
        #check s for correctness somewhere
        self.p=p
        self.n=n
        self.s=s
        if n.nnx is not None:
            self.dgridx = (n.nnx-1)*(2**n.nLevels)
        else:
            self.dgridx = 1.0
        if n.nny is not None:
            self.dgridy = (n.nny-1)*(2**n.nLevels)
        else:
            self.dgridy = 1.0
        if n.nnz is not None:
            self.dgridz = (n.nnz-1)*(2**n.nLevels)
        else:
            self.dgridz = 1.0
        self.plotOffSet = None
        self.stepPlotCalled = {}
        self.stepPlotCalled['exact']=False; self.stepPlotCalled['elementQuantities']=False
        self.plotWindowStart= {}
        if self.s.viewComponents == 'All':
            self.s.viewComponents = list(range(self.p.coefficients.nc))
    def windowNumber(self):
        global windowNumber
        return windowNumber
    def plotNumber(self):
        global plotNumber
        return plotNumber
    def preprocess(self,mlvt,tsim):
        if (('Init' in self.s.viewTimes or
             'All' in self.s.viewTimes or
             tsim in self.s.viewTimes)
            and
            'u' in self.s.viewQuantities):
            if self.plotOffSet is None:
                self.plotOffSet = self.windowNumber()
            self.windowNumberTmp= mlvt.levelModelList[-1].viewSolution(plotOffSet=self.plotOffSet,
                                                                       titleModifier='',
                                                                       dgridnx=self.dgridx,
                                                                       dgridny=self.dgridy,
                                                                       pause=self.s.viewerPause)
            #should create new windows if plotted here
            self.stepPlotExact(mlvt,tsim)
            self.stepPlotElementQuantities(mlvt,tsim)
    def processTimeLevel(self,mlvt,tsim=None,plotOffSet=None):
        if ('All' in self.s.viewTimes or
            tsim in self.s.viewTimes):
            self.stepProcessPlot(mlvt,tsim)
    def postprocess(self,mlvt,tsim):
        if ('All' in self.s.viewTimes or
            'Last' in self.s.viewTimes or
            tsim in self.s.viewTimes):
            self.stepProcessPlot(mlvt,tsim)
    def stepProcessPlot(self,mlvt,tsim):
        """plot desired quantities for a single step
        
        Parameters
        ----------
           mlvt : multilevel vector transport that holds the quantities to measure
           tsim : simulation time

        assumes this is the correct time to plot

        """
        import pdb
        nplots = 0
        if 'u' in self.s.viewQuantities:
            self.windowNumberSave = self.windowNumber()
            mlvt.levelModelList[-1].viewSolution(plotOffSet=self.plotOffSet,titleModifier='',
                                                 dgridnx=self.dgridx,dgridny=self.dgridy,pause=self.s.viewerPause)
            if self.plotOffSet is None:
                self.plotOffSet = self.windowNumberSave
        #end if

        self.stepPlotExact(mlvt,tsim)
        self.stepPlotElementQuantities(mlvt,tsim)
    def stepPlotExact(self,mlvt,tsim):
        """plot 'exact' value of desired quantities for a single step

        Parameters
        ----------
           mlvt :  multilevel vector transport that holds the quantities to measure
           tsim : simulation time

        assumes this is the correct time to plot
        and plotOffSet is set correctly

        """
#        TO DO: Fix scaling for exact vector components to match Transport
        #mwf taking a lot of time on jade
        if ('u_exact' not in self.s.viewQuantities) and ('velocity_exact' not in self.s.viewQuantities):
            return
        global windowNumber
        try:
            from proteusGraphical import vtkViewers
        except:
            return
        vt = mlvt.levelModelList[-1]
        self.windowNumberSave = self.windowNumber()
        #try not to orphan exact plots
        if self.stepPlotCalled['exact'] == True:
            windowNumber = self.plotWindowStart['exact']
        matlabNodalPointsWritten = False #keep track of data structures written for matlab
        for ci in range(self.p.coefficients.nc):
            if (ci in self.s.viewComponents):
                plotExact= 'u_exact' in self.s.viewQuantities and \
                           self.p.analyticalSolution is not None and \
                           ci in self.p.analyticalSolution  and \
                           self.p.analyticalSolution[ci] is not None
                if plotExact:
                #copy the code from VectorTransport.viewSolution as much as possibe
                    if self.viewerType == 'gnuplot':
                        title=vt.coefficients.variableNames[ci]+'_exact: t=%12.5e' % tsim
                        if vt.nSpace_global == 1:
                            xandu = [(vt.mesh.nodeArray[nN,0],self.p.analyticalSolution[ci].uOfXT(vt.mesh.nodeArray[nN],tsim))
                                     for nN in range(vt.mesh.nNodes_global)]
                            xandu.sort()
                            for xu in xandu:
                                self.datFile.write("%12.5e %12.5e \n" % (xu[0],xu[1]))
                            self.datFile.write("\n \n")
                            cmd = "set term x11 %i; plot \'%s\' index %i with linespoints title \"%s\" \n" % (self.windowNumber(),
                                                                                                              self.datFilename,
                                                                                                              self.plotNumber(),
                                                                                                              title)
                            self.cmdFile.write(cmd)
                            self.viewerPipe.write(cmd)
                            newPlot()
                            newWindow()
                        #end if 1d
                        elif vt.nSpace_global == 2:
                            for x in vt.mesh.nodeArray[:,:]:
                                uex = self.p.analyticalSolution[ci].uOfXT(x,tsim)
                                self.datFile.write("%12.5e %12.5e %12.5e \n" % (x[0],x[1],uex))
                            self.datFile.write("\n \n")
                            cmd = "set dgrid3d %d,%d,16; set contour base; set term x11 %i; splot \'%s\' index %i with lines title \"%s\" \n" % (self.dgridx,
                                                                                                                                                 self.dgridy,
                                                                                                                                                 self.windowNumber(),
                                                                                                                                                 self.datFilename,
                                                                                                                                                 self.plotNumber(),
                                                                                                                                                 title)
                            self.cmdFile.write(cmd)
                            self.viewerPipe.write(cmd)
                            newPlot()
                            newWindow()
                        #end 2d
                        elif vt.nSpace_global == 3:
                            (slice_x,slice_y,slice_z) = vt.mesh.nodeArray[vt.mesh.nodeArray.shape[0]//2,:]
                            for x in vt.mesh.nodeArray[:,:]:
                                uex = self.p.analyticalSolution[ci].uOfXT(x,tsim)
                                if x[0] == slice_x:
                                    self.datFile.write("%12.5e %12.5e %12.5e\n" % (x[1],x[2],uex))
                            self.datFile.write("\n \n")
                            cmd = "set dgrid3d; set contour base; set term x11 %i; splot \'%s\' index %i with lines title \"%s-x-slice\" \n" % (self.windowNumber(),
                                                                                                                                        self.datFilename,
                                                                                                                                        self.plotNumber(),title)
                            self.cmdFile.write(cmd)
                            self.viewerPipe.write(cmd)
                            newPlot()
                            newWindow()
                            for x in vt.mesh.nodeArray[:,:]:
                                uex = self.p.analyticalSolution[ci].uOfXT(x,tsim)
                                if x[1] == slice_y:
                                    self.datFile.write("%12.5e %12.5e %12.5e\n" % (x[0],x[2],uex))
                            self.datFile.write("\n \n")
                            cmd = "set dgrid3d; set contour base; set term x11 %i; splot \'%s\' index %i with lines title \"%s-y-slice\" \n" % (self.windowNumber(),
                                                                                                                                        self.datFilename,
                                                                                                                                        self.plotNumber(),title)
                            self.cmdFile.write(cmd)
                            self.viewerPipe.write(cmd)
                            newPlot()
                            newWindow()
                            for x in vt.mesh.nodeArray[:,:]:
                                uex = self.p.analyticalSolution[ci].uOfXT(x,tsim)
                                if x[2] == slice_z:
                                    self.datFile.write("%12.5e %12.5e %12.5e\n" % (x[0],x[1],uex))
                            self.datFile.write("\n \n")
                            cmd = "set dgrid3d; set contour base; set term x11 %i; splot \'%s\' index %i with lines title \"%s-z-slice\" \n" % (self.windowNumber(),
                                                                                                                                        self.datFilename,
                                                                                                                                        self.plotNumber(),title)
                            self.cmdFile.write(cmd)
                            self.viewerPipe.write(cmd)
                            newPlot()
                            newWindow()
                        #end 3d
                    #end gnuplot
                    elif self.viewerType == 'matlab':
                        #assume matlab data structures will be written elsewhere
                        title=vt.coefficients.variableNames[ci]+'-exact: t=%12.5e' % tsim
                        name =vt.coefficients.variableNames[ci]

                        writer = MatlabWriter(nxgrid=50,nygrid=50,nzgrid=50)
                        nplotted = writer.viewScalarAnalyticalFunction(self.cmdFile,vt.nSpace_global,
                                                                       self.p.analyticalSolution[ci].uOfXT,tsim,
                                                                       vt.mesh.nodeArray,vt.mesh.elementNodesArray,
                                                                       name=name,storeMeshData=not self.meshDataStructuresWritten,
                                                                       figureNumber =self.windowNumber()+1,title=title)

                        windowNumber += nplotted

                    elif self.viewerType == 'vtk':
                        title=vt.coefficients.variableNames[ci]+'_exact'
                        if vt.nSpace_global == 1:
                            xvals = []; yvals = []
                            for x in vt.mesh.nodeArray:
                                uex = self.p.analyticalSolution[ci].uOfXT(x,tsim)
                                xvals.append(x[0]); yvals.append(uex)
                            #
                            vtkViewers.viewScalar_1D(xvals,yvals,"x",vt.coefficients.variableNames[ci]+'_exact',title,
                                                        self.windowNumber(),
                                                        Pause=self.viewerPause,sortPoints=True)

                            newPlot()
                            newWindow()
                        #1d
                    #vtk
                #end plotExact
                plotExactVel = ('velocity_exact' in self.s.viewQuantities and
                                'p.analyticalSolutionVelocity' in dir(p) and
                                self.p.p.analyticalSolutionVelocity is not None and
                                ('velocity',ci) in vt.q)
                if plotExactVel:
                    import math
                    if self.viewerType == 'gnuplot':
                        title=vt.coefficients.variableNames[ci]+'velocity_exact: t=%12.5e' % tsim
                        #to scale need exact solution values everywhere first
                        v = numpy.zeros(vt.q[('velocity',ci)].shape,'d')
                        if vt.nSpace_global == 1:
                            max_u = 0.0;
                            for eN in range(vt.mesh.nElements_global):
                                for k in range(vt.nQuadraturePoints_element):
                                    xtmp = vt.q['x'][eN,k,:];
                                    v[eN,k,:] = self.p.p.analyticalSolutionVelocity[ci].uOfXT(xtmp,tsim)
                                    max_u=max(abs(v[eN,k,0]),max_u)
                            scale = 10.*max_u
                            if abs(scale) < 1.0e-12:
                                scale = 1.0
                            for eN in range(vt.mesh.nElements_global):
                                for k in range(vt.nQuadraturePoints_element):
                                    xtmp = vt.q['x'][eN,k,:];
                                    vtmp = v[eN,k,:]
                            self.datFile.write("%12.5e %12.5e \n" % (xtmp[0],vtmp[0]/scale))
                            cmd = "set term x11 %i; plot \'%s\' index %i with linespoints title \"%s\" \n" % (self.windowNumber(),
                                                                                                              self.datFilename,
                                                                                                              self.plotNumber(),
                                                                                                              title)
                            self.cmdFile.write(cmd)
                            self.viewerPipe.write(cmd)
                            newPlot()
                            newWindow()
                        elif vt.nSpace_global == 2:
                            max_u = 0.0; max_v =0.0;
                            for eN in range(vt.mesh.nElements_global):
                                for k in range(vt.nQuadraturePoints_element):
                                    xtmp = vt.q['x'][eN,k,:];
                                    v[eN,k,:] = self.p.p.analyticalSolutionVelocity[ci].uOfXT(xtmp,tsim)
                                    max_u=max(max_u,abs(v[eN,k,0]))
                                    max_v=max(max_u,abs(v[eN,k,1]))
                            scale = 10.0*math.sqrt(max_u**2 + max_v**2)
                            if abs(scale) < 1.e-12:
                                scale = 1.0
                            for eN in range(vt.mesh.nElements_global):
                                for k in range(vt.nQuadraturePoints_element):
                                    xtmp = vt.q['x'][eN,k,:];
                                    vtmp = v[eN,k,:]
                                    self.datFile.write("%12.5e %12.5e %12.5e %12.5e \n" % (xtmp[0],xtmp[1],
                                                                                                  vtmp[0]/scale,vtmp[1]/scale))
                            self.datFile.write("\n \n")
                            cmd = "set term x11 %i; plot \'%s\' index %i with vectors title \"%s\" \n" % (self.windowNumber(),
                                                                                                          self.datFilename,
                                                                                                          self.plotNumber(),
                                                                                                          title)

                            self.cmdFile.write(cmd)
                            self.viewerPipe.write(cmd)
                            newPlot()
                            newWindow()
                        elif vt.nSpace_global == 3:
                            max_u = 0.0; max_v =0.0; max_w = 0.0;
                            (slice_x,slice_y,slice_z) = vt.mesh.nodeArray[vt.mesh.nodeArray.shape[0]/2,:]
                            for eN in range(vt.mesh.nElements_global):
                                for k in range(vt.nQuadraturePoints_element):
                                    xtmp = vt.q['x'][eN,k,:];
                                    v[eN,k,:] = self.p.p.analyticalSolutionVelocity[ci].uOfXT(xtmp,tsim)
                                    max_u=max(max_u,abs(v[eN,k,0]))
                                    max_v=max(max_u,abs(v[eN,k,1]))
                                    max_w=max(max_w,abs(v[eN,k,2]))
                            scale = 10.0*math.sqrt(max_u**2 + max_v**2 + max_w**2)
                            if abs(scale) < 1.e-12:
                                scale = 1.0
                            for eN in range(vt.mesh.nElements_global):
                                for k in range(vt.nQuadraturePoints_element):
                                    xtmp = vt.q['x'][eN,k,:];
                                    vtmp = v[eN,k,:]
                                    if abs(xtmp[0]- slice_x) < vt.mesh.h:
                                        self.datFile.write("%12.5e %12.5e %12.5e %12.5e \n" % (xtmp[1],xtmp[2],
                                                                                                  vtmp[1]/scale,vtmp[2]/scale))
                            self.datFile.write("\n \n")
                            cmd = "set term x11 %i; plot \'%s\' index %i with vectors title \"%s\" \n" % (self.windowNumber(),
                                                                                                          self.datFilename,
                                                                                                          self.plotNumber(),
                                                                                                          title+' x-slice')
                            self.cmdFile.write(cmd)
                            self.viewerPipe.write(cmd)
                            newPlot()
                            newWindow()
                            #yslice
                            for eN in range(vt.mesh.nElements_global):
                                for k in range(vt.nQuadraturePoints_element):
                                    xtmp = vt.q['x'][eN,k,:];
                                    vtmp = v[eN,k,:]
                                    if abs(xtmp[1]- slice_y) < vt.mesh.h:
                                        self.datFile.write("%12.5e %12.5e %12.5e %12.5e \n" % (xtmp[0],xtmp[2],
                                                                                                  vtmp[0]/scale,vtmp[2]/scale))
                            self.datFile.write("\n \n")
                            cmd = "set term x11 %i; plot \'%s\' index %i with vectors title \"%s\" \n" % (self.windowNumber(),
                                                                                                          self.datFilename,
                                                                                                          self.plotNumber(),
                                                                                                          title+' y-slice')
                            self.cmdFile.write(cmd)
                            self.viewerPipe.write(cmd)
                            newPlot()
                            newWindow()
                            #zslice
                            for eN in range(vt.mesh.nElements_global):
                                for k in range(vt.nQuadraturePoints_element):
                                    xtmp = vt.q['x'][eN,k,:];
                                    vtmp = v[eN,k,:]
                                    if abs(xtmp[2]- slice_z) < vt.mesh.h:
                                        self.datFile.write("%12.5e %12.5e %12.5e %12.5e \n" % (xtmp[0],xtmp[1],
                                                                                                  vtmp[0]/scale,vtmp[1]/scale))
                            self.datFile.write("\n \n")
                            cmd = "set term x11 %i; plot \'%s\' index %i with vectors title \"%s\" \n" % (self.windowNumber(),
                                                                                                          self.datFilename,
                                                                                                          self.plotNumber(),
                                                                                                          title+' z-slice')
                            self.cmdFile.write(cmd)
                            self.viewerPipe.write(cmd)
                            newPlot()
                            newWindow()

                        #end 3d
                    #gnuplot
                    elif self.viewerType == 'matlab':
                        title=vt.coefficients.variableNames[ci]+'velocity-exact: t=%12.5e' % tsim
                        name =vt.coefficients.variableNames[ci]+'velocity'

                        writer = MatlabWriter(nxgrid=50,nygrid=50,nzgrid=50)
                        nplotted = writer.viewVectorAnalyticalFunction(self.cmdFile,vt.nSpace_global,
                                                                       self.p.p.analyticalSolutionVelocity[ci].uOfXT,tsim,
                                                                       vt.mesh.nodeArray,vt.mesh.elementNodesArray,
                                                                       name=name,storeMeshData=not self.meshDataStructuresWritten,
                                                                       figureNumber =self.windowNumber()+1,title=title)

                        windowNumber += nplotted
                    #need vtk option

            #end components
        #end ci
        #vector components
        if vt.coefficients.vectorComponents is not None:
            title = 'velocity_exact : t=%12.5e' % tsim
            if vt.nSpace_global == 2:
                uci = vt.coefficients.vectorComponents[0]; vci = vt.coefficients.vectorComponents[1]
                plotVector = (uci in self.s.viewComponents and vci in self.s.viewComponents and
                              self.p.analyticalSolution is not None and
                              uci in self.p.analyticalSolution and vci in self.p.analyticalSolution and
                              self.p.analyticalSolution[uci] is not None and self.p.analyticalSolution[vci] is not None)
                if plotVector and self.viewerType == 'gnuplot':
                    for x in vt.mesh.nodeArray[:,:]:
                        uex = self.p.analyticalSolution[uci].uOfXT(x,tsim)
                        vex = self.p.analyticalSolution[vci].uOfXT(x,tsim)
                        self.datFile.write("%12.5e %12.5e %12.5e %12.5e \n" % (x[0],x[1],uex,vex))
                    self.datFile.write("\n \n")
                    cmd = "set term x11 %i; plot \'%s\' index %i with vectors title \"%s\" \n" % (self.windowNumber(),
                                                                                                  self.datFilename,
                                                                                                  self.plotNumber(),
                                                                                                  title)
                    self.cmdFile.write(cmd)
                    self.viewerPipe.write(cmd)
                    newPlot()
                    newWindow()
            elif vt.nSpace_global == 3:
                (slice_x,slice_y,slice_z) = vt.mesh.nodeArray[vt.mesh.nodeArray.shape[0]/2,:]
                uci = vt.coefficients.vectorComponents[0]; vci = vt.coefficients.vectorComponents[1]
                wci = vt.coefficients.vectorComponents[2]
                plotVector = (uci in self.s.viewComponents and vci in self.s.viewComponents and
                              wci in self.s.viewComponents and self.p.analyticalSolution is not None and
                              self.p.analyticalSolution is not None and
                              uci in self.p.analyticalSolution and vci in self.p.analyticalSolution and
                              wci in self.p.analyticalSolution and
                              self.p.analyticalSolution[uci] is not None and self.p.analyticalSolution[vci] is not None and
                              self.p.analyticalSolution[wci] is not None)

                if plotVector and self.viewerType == 'gnuplot':
                    for x in vt.mesh.nodeArray[:,:]:
                        uex = self.p.analyticalSolution[uci].uOfXT(x,tsim)
                        vex = self.p.analyticalSolution[vci].uOfXT(x,tsim)
                        wex = self.p.analyticalSolution[wci].uOfXT(x,tsim)
                        if x[0] == slice_x:
                            self.datFile.write("%12.5e %12.5e %12.5e %12.5e \n" % (x[1],x[2],vex,wex))

                    self.datFile.write("\n \n")
                    cmd = "set term x11 %i; plot \'%s\' index %i with vectors title \"%s\" \n" % (self.windowNumber(),
                                                                                                  self.datFilename,
                                                                                                  self.plotNumber(),
                                                                                                  title+' x-slice')
                    self.cmdFile.write(cmd)
                    self.viewerPipe.write(cmd)
                    newPlot()
                    newWindow()
                    for x in vt.mesh.nodeArray[:,:]:
                        uex = self.p.analyticalSolution[uci].uOfXT(x,tsim)
                        vex = self.p.analyticalSolution[vci].uOfXT(x,tsim)
                        wex = self.p.analyticalSolution[wci].uOfXT(x,tsim)
                        if x[1] == slice_y:
                            self.datFile.write("%12.5e %12.5e %12.5e %12.5e \n" % (x[0],x[2],uex,wex))

                    self.datFile.write("\n \n")
                    cmd = "set term x11 %i; plot \'%s\' index %i with vectors title \"%s\" \n" % (self.windowNumber(),
                                                                                                  self.datFilename,
                                                                                                  self.plotNumber(),
                                                                                                  title+' y-slice')
                    self.cmdFile.write(cmd)
                    self.viewerPipe.write(cmd)
                    newPlot()
                    newWindow()
                    for x in vt.mesh.nodeArray[:,:]:
                        uex = self.p.analyticalSolution[uci].uOfXT(x,tsim)
                        vex = self.p.analyticalSolution[vci].uOfXT(x,tsim)
                        wex = self.p.analyticalSolution[wci].uOfXT(x,tsim)
                        if x[2] == slice_z:
                            self.datFile.write("%12.5e %12.5e %12.5e %12.5e \n" % (x[0],x[1],uex,vex))

                    self.datFile.write("\n \n")
                    cmd = "set term x11 %i; plot \'%s\' index %i with vectors title \"%s\" \n" % (self.windowNumber(),
                                                                                                  self.datFilename,
                                                                                                  self.plotNumber(),
                                                                                                  title+' z-slice')
                    self.cmdFile.write(cmd)
                    self.viewerPipe.write(cmd)
                    newPlot()
                    newWindow()
                #end plot vector
            #end 3d
        #end vector components
        if self.stepPlotCalled['exact'] == False:
            self.plotWindowStart['exact'] = self.windowNumberSave
        self.stepPlotCalled['exact'] = True
    #end def
    def stepPlotElementQuantities(self,mlvt,tsim):
        """
        sort through desired quantities in quadrature dictionaries like m, dm, to plot
          p    --- problem definition
          n    --- numerics definition
          mlvt --- multilevel vector transport that holds the quantities to measure
          tsim --- simulation time
        assumes this is the correct time to plot
        and plotOffSet is set correctly
        """
        global windowNumber
        self.windowNumberSave = self.windowNumber()
        plottedSomething = False
        if self.stepPlotCalled['elementQuantities'] == True:
            windowNumber = self.plotWindowStart['elementQuantities']
        for quant in self.s.viewQuantities:
            recType = quant.split(':')
            if len(recType) > 1 and recType[0] == 'q': #found element quadrature quantity
                stval = eval(recType[1])
                if (stval in mlvt.levelModelList[-1].q and
                    len(mlvt.levelModelList[-1].q[stval].shape) == 2): #found quantity and it's a scalar
                    self.plotScalarElementQuantity(stval,mlvt,tsim)
                    plottedSomething = True
                elif (stval in mlvt.levelModelList[-1].q and
                      len(mlvt.levelModelList[-1].q[stval].shape) == 3): #found quantity and it's a vector
                    self.plotVectorElementQuantity(stval,mlvt,tsim)
                    plottedSomething = True
            elif len(recType) > 1 and recType[0] == 'ebq_global': #found global element boundary quantity
                stval = eval(recType[1])
                if stval in mlvt.levelModelList[-1].ebq_global:
                    if len(mlvt.levelModelList[-1].ebq_global[stval].shape) == 3: #found quantity and its a vector
                        self.plotVectorGlobalElementBoundaryQuantity(stval,mlvt,tsim)
                        plottedSomething = True
        if self.stepPlotCalled['elementQuantities'] == False:
            self.plotWindowStart['elementQuantities'] = self.windowNumberSave
        if plottedSomething:
            self.stepPlotCalled['elementQuantities'] = True
    #
    def plotScalarElementQuantity(self,ckey,mlvt,tsim):
        """plotting routine to look at scalar quantity stored in element quad dictionary q

        Parameters
        -----------
          ckey : what should be plotted
          mlvt : multilevel vector transport that holds the quantities to measure
          tsim : simulation time
        assumes this is the correct time to plot
        and plotOffSet is set correctly

        """
        p = self.p; n = self.n

        from proteusGraphical import vtkViewers
        vt = mlvt.levelModelList[-1]
        title = """q[%s]""" % (ckey,)
        assert ckey in vt.q
        if self.viewerType == 'gnuplot':
            if vt.nSpace_global == 1:
                npoints  = vt.q['x'].shape[0]*vt.q['x'].shape[1]
                xandu = [(vt.q['x'].flat[i*3+0],vt.q[ckey].flat[i]) for i in range(npoints)]
                xandu.sort()
                for xu in xandu:
                    self.datFile.write("%12.5e %12.5e \n" % (xu[0],xu[1]))
                self.datFile.write("\n \n")
                cmd = "set term x11 %i; plot \'%s\' index %i with linespoints title \"%s\" \n" % (self.windowNumber(),
                                                                                                  self.datFilename,
                                                                                                  self.plotNumber(),
                                                                                                  title)
                self.cmdFile.write(cmd)
                self.viewerPipe.write(cmd)
                newPlot()
                newWindow()
            elif vt.nSpace_global == 2:
                for eN in range(vt.q['x'].shape[0]):
                    for k in range(vt.q['x'].shape[1]):
                        self.datFile.write("%12.5e %12.5e %12.5e \n" % (vt.q['x'][eN,k,0],vt.q['x'][eN,k,1],vt.q[ckey][eN,k]))
                self.datFile.write("\n \n")
                ggrid = 50;
                cmd = "set dgrid3d %d,%d,16; set contour base; set term x11 %i; splot \'%s\' index %i with lines title \"%s\" \n" % (self.dgridx,
                                                                                                                                     self.dgridy,
                                                                                                                                     self.windowNumber(),
                                                                                                                                     self.datFilename,
                                                                                                                                     self.plotNumber(),
                                                                                                                                     title)

                self.cmdFile.write(cmd)
                self.viewerPipe.write(cmd)
                newPlot()
                newWindow()
            #end 2d
            elif vt.nSpace_global == 3:
                (slice_x,slice_y,slice_z) = vt.mesh.nodeArray[vt.mesh.nodeArray.shape[0]/2,:]
                for eN in range(vt.q['x'].shape[0]):
                    for k in range(vt.q['x'].shape[1]):
                        if vt.q['x'][eN,k,0] == slice_x:
                            self.datFile.write("%12.5e %12.5e %12.5e \n" % (vt.q['x'][eN,k,1],
                                                                               vt.q['x'][eN,k,2],vt.q[ckey][eN,k]))

                self.datFile.write("\n \n")
                cmd = "set dgrid3d; set contour base; set term x11 %i; splot \'%s\' index %i with lines title \"%s-x-slice\" \n" % (self.windowNumber(),
                                                                                                                                            self.datFilename,
                                                                                                                                            self.plotNumber(),
                                                                                                                                            title)
                self.cmdFile.write(cmd)
                self.viewerPipe.write(cmd)
                newPlot()
                newWindow()
                #
                for eN in range(vt.q['x'].shape[0]):
                    for k in range(vt.q['x'].shape[1]):
                        if vt.q['x'][eN,k,1] == slice_y:
                            self.datFile.write("%12.5e %12.5e %12.5e \n" % (vt.q['x'][eN,k,0],
                                                                               vt.q['x'][eN,k,2],vt.q[ckey][eN,k]))

                self.datFile.write("\n \n")
                cmd = "set dgrid3d; set contour base; set term x11 %i; splot \'%s\' index %i with lines title \"%s-y-slice\" \n" % (self.windowNumber(),
                                                                                                                                    self.datFilename,
                                                                                                                                    self.plotNumber(),
                                                                                                                                    title)
                self.cmdFile.write(cmd)
                self.viewerPipe.write(cmd)
                newPlot()
                newWindow()
                #
                for eN in range(vt.q['x'].shape[0]):
                    for k in range(vt.q['x'].shape[1]):
                        if vt.q['x'][eN,k,2] == slice_z:
                            self.datFile.write("%12.5e %12.5e %12.5e \n" % (vt.q['x'][eN,k,0],
                                                                               vt.q['x'][eN,k,1],vt.q[ckey][eN,k]))

                self.datFile.write("\n \n")
                cmd = "set dgrid3d; set contour base; set term x11 %i; splot \'%s\' index %i with lines title \"%s-z-slice\" \n" % (self.windowNumber(),
                                                                                                                                    self.datFilename,
                                                                                                                                    self.plotNumber(),
                                                                                                                                    title)
                self.cmdFile.write(cmd)
                self.viewerPipe.write(cmd)
                newPlot()
                newWindow()

            #3d
        #gnuplot
        elif self.viewerType == 'matlab':
            name = ckey[0];
            for i in range(len(ckey)-1):
                name += "_%s" % ckey[1+i]

            title = "%s t = %g " % (name,tsim)
            #does not handle window number counting internally
            writer = MatlabWriter(nxgrid=50,nygrid=50,nzgrid=50)
            nplotted = writer.viewScalarPointData(self.cmdFile,vt.nSpace_global,vt.q,ckey,name=name,
                                                  storeMeshData=not self.meshDataStructuresWritten,
                                                  useLocal = True,
                                                  figureNumber =self.windowNumber()+1,title=title)
            windowNumber += nplotted
        elif self.viewerType == 'vtk':
            title = """q[%s]""" % (ckey,)

            if vt.nSpace_global == 1:
                npoints  = vt.q['x'].shape[0]*vt.q['x'].shape[1]
                xvals = numpy.array([vt.q['x'].flat[i*3+0] for i in range(npoints)])
                yvals = numpy.array([vt.q[ckey].flat[i] for i in range(npoints)])
#                yvals = [vt.q[ckey].flat[:]
                vtkViewers.viewScalar_1D(xvals,yvals,"x",ckey[0],title,self.windowNumber(),
                                            Pause=self.s.viewerPause,sortPoints=True)
                newPlot()
                newWindow()
            elif vt.nSpace_global == 2:
                vtkViewers.viewScalar_pointSet_2D(vt.q['x'], vt.q[ckey], title, self.windowNumber(), True, self.s.viewerPause, False)
                newPlot()
                newWindow()
            elif vt.nSpace_global == 3:
                vtkViewers.viewScalar_pointSet_3D(vt.q['x'], vt.q[ckey], title, self.windowNumber(), self.s.viewerPause, False)
                newPlot()
                newWindow()

    #def
    def plotVectorGlobalElementBoundaryQuantity(self,ckey,mlvt,tsim):
        """plotting routine to look at vector quantity stored in global elementBoundary quad dictionary ebq_global

        Parameters
        -----------
          ckey : what should be plotted
          mlvt : multilevel vector transport that holds the quantities to measure
          tsim : simulation time

        assumes this is the correct time to plot
        and plotOffSet is set correctly

        """
        from proteusGraphical import vtkViewers
        p = self.p; n = self.n

        vt = mlvt.levelModelList[-1]
        title = """ebq_global[%s] t= %s""" % (ckey,tsim)
        assert ckey in vt.ebq_global
        if self.viewerType == 'gnuplot':
            if vt.nSpace_global == 1:
                max_u=max(numpy.absolute(numpy.take(vt.ebq_global[ckey],[0],2).flat))
                L = max(vt.mesh.nodeArray[:,0])
                scale = 10.*max_u/L
                if abs(scale) < 1.0e-12:
                    scale = 1.0
                npoints  = vt.ebq_global['x'].shape[0]*vt.ebq_global['x'].shape[1]
                xandu = [(vt.ebq_global['x'].flat[i*3+0],vt.ebq_global[ckey].flat[i]) for i in range(npoints)]
                xandu.sort()
                for xu in xandu:
                    self.datFile.write("%12.5e %12.5e \n" % (xu[0],xu[1]/scale))
                self.datFile.write("\n \n")
                cmd = "set term x11 %i; plot \'%s\' index %i with linespoints title \"%s\" \n" % (self.windowNumber(),
                                                                                                  self.datFilename,
                                                                                                  self.plotNumber(),
                                                                                                  title+" max= "+repr(max_u))
                self.cmdFile.write(cmd)
                self.viewerPipe.write(cmd)
                newPlot()
                newWindow()

            elif vt.nSpace_global == 2:
                max_u=max(numpy.absolute(numpy.take(vt.ebq_global[ckey],[0],2).flat))
                max_v=max(numpy.absolute(numpy.take(vt.ebq_global[ckey],[1],2).flat))
                L = min((max(vt.mesh.nodeArray[:,0]),max(vt.mesh.nodeArray[:,1])))
                scale =10.0*max((max_u,max_v,1.0e-16))/L
                if abs(scale) < 1.e-12:
                    scale = 1.0
                for ebN in range(vt.ebq_global[ckey].shape[0]):
                    for k in range(vt.ebq_global[ckey].shape[1]):
                        xtmp =vt.ebq_global['x'][ebN,k,:]; vtmp = vt.ebq_global[ckey][ebN,k,:]
                        self.datFile.write("%12.5e %12.5e %12.5e %12.5e \n" % (xtmp[0],xtmp[1],
                                                                                  vtmp[0]/scale,vtmp[1]/scale))
                self.datFile.write("\n \n")
                cmd = "set term x11 %i; plot \'%s\' index %i with vectors title \"%s\" \n" % (self.windowNumber(),
                                                                                              self.datFilename,
                                                                                              self.plotNumber(),
                                                                                              title+" max=(%s,%s) " % (max_u,max_v))
                self.cmdFile.write(cmd)
                self.viewerPipe.write(cmd)
                newPlot()
                newWindow()
                #mwf debug
                #raw_input('simTools coef press return to continue\n')
            #end 2d
            elif vt.nSpace_global == 3:
                (slice_x,slice_y,slice_z) = vt.mesh.nodeArray[vt.mesh.nodeArray.shape[0]//2,:]
                max_u=max(numpy.absolute(numpy.take(vt.ebq_global[ckey],[0],2).flat))
                max_v=max(numpy.absolute(numpy.take(vt.ebq_global[ckey],[1],2).flat))
                max_w=max(numpy.absolute(numpy.take(vt.ebq_global[ckey],[2],2).flat))
                L = min((max(vt.mesh.nodeArray[:,0]),max(vt.mesh.nodeArray[:,1]),
                         max(vt.mesh.nodeArray[:,1])))
                scale = 10.0*max((max_u,max_v,max_w,1.e-16))/L
                if abs(scale) < 1.e-12:
                    scale = 1.0
                #x slice
                for ebN in range(vt.ebq_global[ckey].shape[0]):
                    for k in range(vt.ebq_global[ckey].shape[1]):
                        xtmp = vt.ebq_global['x'][ebN,k,:]; vtmp = vt.ebq_global[ckey][ebN,k,:]
                        if abs(xtmp[0]-slice_x) < vt.mesh.h:
                            self.datFile.write("%12.5e %12.5e %12.5e %12.5e \n" % (xtmp[1],xtmp[2],
                                                                                      vtmp[1]/scale,vtmp[2]/scale))

                self.datFile.write("\n \n")
                cmd = "set term x11 %i; plot \'%s\' index %i with vectors title \"%s\" \n" % (self.windowNumber(),
                                                                                              self.datFilename,
                                                                                              self.plotNumber(),
                                                                                              title+" max=(%s,%s,%s) " % (max_u,max_v,max_w)+" : x-slice")


                self.cmdFile.write(cmd)
                self.viewerPipe.write(cmd)
                newPlot()
                newWindow()
                #y slice
                for ebN in range(vt.ebq_global[ckey].shape[0]):
                    for k in range(vt.ebq_global[ckey].shape[1]):
                        xtmp = vt.ebq_global['x'][ebN,k,:]; vtmp = vt.ebq_global[ckey][ebN,k,:]
                        if abs(xtmp[0]-slice_y) < vt.mesh.h:
                            self.datFile.write("%12.5e %12.5e %12.5e %12.5e \n" % (xtmp[1],xtmp[2],
                                                                                      vtmp[1]/scale,vtmp[2]/scale))

                self.datFile.write("\n \n")
                cmd = "set term x11 %i; plot \'%s\' index %i with vectors title \"%s\" \n" % (self.windowNumber(),
                                                                                              self.datFilename,
                                                                                              self.plotNumber(),
                                                                                              title+" max=(%s,%s,%s) " % (max_u,max_v,max_w)+" : y-slice")
                self.cmdFile.write(cmd)
                self.viewerPipe.write(cmd)
                newPlot()
                newWindow()

                #z slice
                for ebN in range(vt.ebq_global[ckey].shape[0]):
                    for k in range(vt.ebq_global[ckey].shape[1]):
                        xtmp = vt.ebq_global['x'][ebN,k,:]; vtmp = vt.ebq_global[ckey][ebN,k,:]
                        if abs(xtmp[0]-slice_z) < vt.mesh.h:
                            self.datFile.write("%12.5e %12.5e %12.5e %12.5e \n" % (xtmp[1],xtmp[2],
                                                                                      vtmp[1]/scale,vtmp[2]/scale))

                self.datFile.write("\n \n")
                cmd = "set term x11 %i; plot \'%s\' index %i with vectors title \"%s\" \n" % (self.windowNumber(),
                                                                                              self.datFilename,
                                                                                              self.plotNumber(),
                                                                                              title+" max=(%s,%s,%s) " % (max_u,max_v,max_w)+" : z-slice")


                self.cmdFile.write(cmd)
                self.viewerPipe.write(cmd)
                newPlot()
                newWindow()
            #3d
        #gnuplot
        elif self.viewerType == 'vtk':
            title = """ebq_global[%s]""" % (ckey,)
            if vt.nSpace_global == 1:
                max_u=max(numpy.absolute(numpy.take(vt.ebq_global[ckey],[0],2).flat))
                L = max(vt.mesh.nodeArray[:,0])
                scale = 10.*max_u/L
                if abs(scale) < 1.0e-12:
                    scale = 1.0
                npoints  = vt.ebq_global['x'].shape[0]*vt.ebq_global['x'].shape[1]
                xvals = [vt.ebq_global['x'].flat[i*3+0] for i in range(npoints)]
                yvals = [vt.ebq_global[ckey].flat[i]/scale for i in range(npoints)]
                vtkViewers.viewScalar_1D(xvals,yvals,"x",ckey[0],title,self.windowNumber(),
                                            Pause=self.s.viewerPause,sortPoints=True)
                newPlot()
                newWindow()

            elif vt.nSpace_global == 2:
                max_u=max(numpy.absolute(numpy.take(vt.ebq_global[ckey],[0],2).flat))
                max_v=max(numpy.absolute(numpy.take(vt.ebq_global[ckey],[1],2).flat))
                L = min((max(vt.mesh.nodeArray[:,0]),max(vt.mesh.nodeArray[:,1])))
                scale =10.0*max((max_u,max_v,1.0e-16))/L
                if abs(scale) < 1.e-12:
                    scale = 1.0
                npoints  = vt.ebq_global['x'].shape[0]*vt.ebq_global['x'].shape[1]
#                x = [vt.ebq_global['x'].flat[i*3+0] for i in range(npoints)]
#                y = [vt.ebq_global['x'].flat[i*3+1] for i in range(npoints)]
#                z = [vt.ebq_global['x'].flat[i*3+2] for i in range(npoints)]
                xvals= [vt.ebq_global[ckey].flat[i*2+0]/scale for i in range(npoints)]
                yvals= [vt.ebq_global[ckey].flat[i*2+1]/scale for i in range(npoints)]
                nodes = vt.ebq_global['x'].flat[:]
                vtkViewers.viewVector_pointSet_2D(nodes,xvals,yvals,None,title,self.windowNumber(),
                                                         arrows=True,streamlines=False,
                                                         Pause=self.s.viewerPause)
#                vtkDisplay2DVectorMeshGeneric(x,y,z,xvals,yvals,None,title,self.windowNumber(),
#                                                         arrows=True,streamlines=False,
#                                                         Pause=self.flags['plotOptions']['vtk']['pause'])
                newPlot()
                newWindow()
            elif vt.nSpace_global == 3:
                max_u=max(numpy.absolute(numpy.take(vt.ebq_global[ckey],[0],2).flat))
                max_v=max(numpy.absolute(numpy.take(vt.ebq_global[ckey],[1],2).flat))
                max_w=max(numpy.absolute(numpy.take(vt.ebq_global[ckey],[2],2).flat))
                L = min((max(vt.mesh.nodeArray[:,0]),max(vt.mesh.nodeArray[:,1]),max(vt.mesh.nodeArray[:,2])))
                scale =10.0*max((max_u,max_v,max_w,1.0e-16))/L
                if abs(scale) < 1.e-12:
                    scale = 1.0
                npoints  = vt.ebq_global['x'].shape[0]*vt.ebq_global['x'].shape[1]
#                x = [vt.ebq_global['x'].flat[i*3+0] for i in range(npoints)]
#                y = [vt.ebq_global['x'].flat[i*3+1] for i in range(npoints)]
#                z = [vt.ebq_global['x'].flat[i*3+2] for i in range(npoints)]
                nodes = vt.ebq_global['x'].flat[:]
                xvals= [vt.ebq_global[ckey].flat[i*3+0]/scale for i in range(npoints)]
                yvals= [vt.ebq_global[ckey].flat[i*3+1]/scale for i in range(npoints)]
                zvals= [vt.ebq_global[ckey].flat[i*3+2]/scale for i in range(npoints)]
                vtkViewers.viewVector_pointSet_3D(nodes,xvals,yvals,zvals,title,self.windowNumber(),
                                                         arrows=True,streamlines=False,
                                                         Pause=self.s.viewerPause)
#                vtkDisplay3DVectorMeshGeneric(x,y,z,xvals,yvals,zvals,title,self.windowNumber(),
#                                                         arrows=True,streamlines=False,
#                                                         Pause=self.flags['plotOptions']['vtk']['pause'])
                newPlot()
                newWindow()
     #def
    def plotVectorElementQuantity(self,ckey,mlvt,tsim,nVectorPlotPointsPerElement=1):
        """
        plotting routine to look at vector quantity stored in global element quad dictionary q
         input :
          ckey --- what should be plotted
          p    --- problem definition
          n    --- numerics definition
          mlvt --- multilevel vector transport that holds the quantities to measure
          tsim --- simulation time
        assumes this is the correct time to plot
        and plotOffSet is set correctly

        """
        from proteusGraphical import vtkViewers
        p = self.p; n = self.n

        vt = mlvt.levelModelList[-1]
        title = """q[%s] t= %s""" % (ckey,tsim)
        assert ckey in vt.q
        if self.viewerType == 'gnuplot':
            if vt.nSpace_global == 1:
                max_u=max(numpy.absolute(numpy.take(vt.q[ckey],[0],2).flat))
                L = max(vt.mesh.nodeArray[:,0])
                scale = 10.*max_u/L
                if abs(scale) < 1.0e-12:
                    scale = 1.0
                npoints  = vt.q['x'].shape[0]*vt.q['x'].shape[1]
                xandu = [(vt.q['x'].flat[i*3+0],vt.q[ckey].flat[i]) for i in range(npoints)]
                xandu.sort()
                for xu in xandu:
                    self.datFile.write("%12.5e %12.5e \n" % (xu[0],xu[1]/scale))
                self.datFile.write("\n \n")
                ptitle = title+" max= %g" % max_u
                cmd = "set term x11 %i; plot \'%s\' index %i with linespoints title \"%s\" \n" % (self.windowNumber(),
                                                                                                  self.datFilename,
                                                                                                  self.plotNumber(),
                                                                                                  ptitle)
                self.cmdFile.write(cmd)
                self.viewerPipe.write(cmd)
                newPlot()
                newWindow()
            elif vt.nSpace_global == 2:
                max_u=max(numpy.absolute(numpy.take(vt.q[ckey],[0],2).flat))
                max_v=max(numpy.absolute(numpy.take(vt.q[ckey],[1],2).flat))
                L = min((max(vt.mesh.nodeArray[:,0]),max(vt.mesh.nodeArray[:,1])))
                scale =10.0*max((max_u,max_v,1.0e-16))/L
                if abs(scale) < 1.e-12:
                    scale = 1.0
                for eN in range(vt.mesh.nElements_global):
                    #mwf what about just 1 point per element for k in range(vt.nQuadraturePoints_element):
                    for k in range(min(nVectorPlotPointsPerElement,vt.nQuadraturePoints_element)):
                        xtmp = vt.q['x'][eN,k,:]; vtmp = vt.q[ckey][eN,k,:]
                        self.datFile.write("%12.5e %12.5e %12.5e %12.5e \n" % (xtmp[0],xtmp[1],
                                                                                  vtmp[0]/scale,vtmp[1]/scale))
                self.datFile.write("\n \n")
                ptitle = title + "max=(%s,%s)" % (max_u,max_v)
                cmd = "set term x11 %i; plot \'%s\' index %i with vectors title \"%s\" \n" % (self.windowNumber(),
                                                                                              self.datFilename,
                                                                                              self.plotNumber(),
                                                                                              ptitle)

                self.cmdFile.write(cmd)
                self.viewerPipe.write(cmd)
                newPlot()
                newWindow()
            elif vt.nSpace_global == 3:
                (slice_x,slice_y,slice_z) = vt.mesh.nodeArray[vt.mesh.nodeArray.shape[0]//2,:]
                max_u=max(numpy.absolute(numpy.take(vt.q[ckey],[0],2).flat))
                max_v=max(numpy.absolute(numpy.take(vt.q[ckey],[1],2).flat))
                max_w=max(numpy.absolute(numpy.take(vt.q[ckey],[2],2).flat))
                L = min((max(vt.mesh.nodeArray[:,0]),max(vt.mesh.nodeArray[:,1]),
                         max(vt.mesh.nodeArray[:,1])))
                scale = 10.0*max((max_u,max_v,max_w,1.e-16))/L
                if abs(scale) < 1.e-12:
                    scale = 1.0
                for eN in range(vt.mesh.nElements_global):
                    #mwf now try one point per element for k in range(vt.nQuadraturePoints_element):
                    for k in range(min(nVectorPlotPointsPerElement,vt.nQuadraturePoints_element)):
                        xtmp = vt.q['x'][eN,k,:]; vtmp = vt.q[ckey][eN,k,:]
                        if abs(xtmp[0]-slice_x) < vt.mesh.h:
                            self.datFile.write("%12.5e %12.5e %12.5e %12.5e \n" % (xtmp[1],xtmp[2],
                                                                                      vtmp[1]/scale,vtmp[2]/scale))
                self.datFile.write("\n \n")
                ptitle = title + " max=(%s,%s,%s) " % (max_u,max_v,max_w)+" : x-slice"
                cmd = "set term x11 %i; plot \'%s\' index %i with vectors title \"%s\" \n" % (self.windowNumber(),
                                                                                              self.datFilename,
                                                                                              self.plotNumber(),
                                                                                              ptitle)
                self.cmdFile.write(cmd)
                self.viewerPipe.write(cmd)
                newPlot()
                newWindow()
                #y slice
                for eN in range(vt.mesh.nElements_global):
                    #mwf now try one point per element for k in range(vt.nQuadraturePoints_element):
                    for k in range(min(nVectorPlotPointsPerElement,vt.nQuadraturePoints_element)):
                        xtmp = vt.q['x'][eN,k,:]; vtmp = vt.q[ckey][eN,k,:]
                        if abs(xtmp[1]-slice_y) < vt.mesh.h:
                            self.datFile.write("%12.5e %12.5e %12.5e %12.5e \n" % (xtmp[0],xtmp[2],
                                                                                      vtmp[0]/scale,vtmp[2]/scale))
                self.datFile.write("\n \n")
                ptitle = title + " max=(%s,%s,%s) " % (max_u,max_v,max_w)+" : y-slice"
                cmd = "set term x11 %i; plot \'%s\' index %i with vectors title \"%s\" \n" % (self.windowNumber(),
                                                                                              self.datFilename,
                                                                                              self.plotNumber(),
                                                                                              ptitle)
                self.cmdFile.write(cmd)
                self.viewerPipe.write(cmd)
                newPlot()
                newWindow()
                #z slice
                for eN in range(vt.mesh.nElements_global):
                    #mwf now try one point per element for k in range(vt.nQuadraturePoints_element):
                    for k in range(min(nVectorPlotPointsPerElement,vt.nQuadraturePoints_element)):
                        xtmp = vt.q['x'][eN,k,:]; vtmp = vt.q[ckey][eN,k,:]
                        if abs(xtmp[2]-slice_z) < vt.mesh.h:
                            self.datFile.write("%12.5e %12.5e %12.5e %12.5e \n" % (xtmp[0],xtmp[1],
                                                                                      vtmp[0]/scale,vtmp[1]/scale))
                self.datFile.write("\n \n")
                ptitle = title + " max=(%s,%s,%s) " % (max_u,max_v,max_w)+" : z-slice"
                cmd = "set term x11 %i; plot \'%s\' index %i with vectors title \"%s\" \n" % (self.windowNumber(),
                                                                                              self.datFilename,
                                                                                              self.plotNumber(),
                                                                                              ptitle)
                self.cmdFile.write(cmd)
                self.viewerPipe.write(cmd)
                newPlot()
                newWindow()
            #end 3d
        #gnuplot
        elif self.viewerType == 'matlab':
            name = ckey[0];
            for i in range(len(ckey)-1):
                name += "_%s" % ckey[1+i]

            title = "%s t = %g " % (name,tsim)
            #does not handle window number counting internally
            writer = MatlabWriter(nxgrid=50,nygrid=50,nzgrid=50)
            nplotted = writer.viewVectorPointData(self.cmdFile,vt.nSpace_global,vt.q,ckey,name=name,
                                                  storeMeshData=not self.meshDataStructuresWritten,
                                                  useLocal=False,#not implemented yed
                                                  figureNumber =self.windowNumber()+1,title=title)
            windowNumber += nplotted

        elif self.viewerType == 'vtk':
            title = """q[%s]""" % (ckey,)
            if vt.nSpace_global == 1:
                max_u=max(numpy.absolute(numpy.take(vt.q[ckey],[0],2).flat))
                L = max(vt.mesh.nodeArray[:,0])
                scale = 1.0
                if abs(scale) < 1.0e-12:
                    scale = 1.0
                npoints  = vt.q['x'].shape[0]*vt.q['x'].shape[1]
                xvals = [vt.q['x'].flat[i*3+0] for i in range(npoints)]
                yvals = [vt.q[ckey].flat[i]/scale for i in range(npoints)]
                vtkViewers.viewVector_1D(xvals,yvals,"x",ckey[0],title,self.windowNumber(),
                                            Pause=self.s.viewerPause)
                newPlot()
                newWindow()
            #1d
            elif vt.nSpace_global == 2:
                vtkViewers.viewVector_pointSet_2D(vt.q['x'],vt.q[ckey],title)
                newPlot()
                newWindow()
            elif vt.nSpace_global == 3:
                vtkViewers.viewVector_pointSet_3D(vt.q['x'],vt.q[ckey],title,self.windowNumber(),
                                                  Pause=self.s.viewerPause)
                newPlot()
                newWindow()

     #def


class MatlabWriter(object):
    """
    collect functionality for generating visualation data and commands in matlab
    TODO:
      C0P2 in 3d
      DG monomials
    """
    def __init__(self,nxgrid=50,nygrid=50,nzgrid=10,verbose=0):
        self.verbose = 0
        self.ngrid=[nxgrid,nygrid,nzgrid] #default grid size if converting to regular mesh
    def storePointMeshData(self,cmdFile,x,name):
        """
        write out spatial locations for generic point data
        """
        cmdFile.write("%s_x_q = [ ... \n" % name)
        for eN in range(x.shape[0]):
            for k in range(x.shape[1]):
                cmdFile.write("%g %g %g \n" % (x[eN,k,0],x[eN,k,1],x[eN,k,2]))
        cmdFile.write("];")
        #

    def viewScalarPointData(self,cmdFile,nSpace,q,ckey,name=None,
                            storeMeshData=True,useLocal=True,
                            figureNumber=1,title=None):
        """
        wrapper for visualling element quadrature points, can try to use a local representation
        or build a global one depending on useLocal.

        If useLocal and  nPoints_elemnet < nSpace+1 calls global routine

        """
        if not useLocal:
            return self.viewGlobalScalarPointData(cmdFile,nSpace,q,ckey,name=name,
                                                  storeMeshData=storeMeshData,
                                                  figureNumber=figureNumber,title=title)
        nPoints_element = q['x'].shape[1]
        if nPoints_element <= nSpace:
            print("""
Warning! viewScalarPointData nPoints_element=%s < %s too small for useLocal, using global interp""" % (nPoints_element,
                                                                                                       nSpace+1))
            return self.viewGlobalScalarPointData(cmdFile,nSpace,q,ckey,name=name,
                                                  storeMeshData=storeMeshData,
                                                  figureNumber=figureNumber,title=title)
        return self.viewLocalScalarPointData(cmdFile,nSpace,q,ckey,name=name,
                                             storeMeshData=storeMeshData,
                                             figureNumber=figureNumber,title=title)

    def viewVectorPointData(self,cmdFile,nSpace,q,ckey,name=None,
                            storeMeshData=True,useLocal=True,
                            figureNumber=1,title=None):
        """
        wrapper for visualling element quadrature points, can try to use a local representation
        or build a global one depending on useLocal.
        TODO: implement local view for vectors

        If useLocal and  nPoints_elemnet < nSpace+1 calls global routine

        """
        if not useLocal:
            return self.viewGlobalVectorPointData(cmdFile,nSpace,q,ckey,name=name,
                                                  storeMeshData=storeMeshData,
                                                  figureNumber=figureNumber,title=title)
        else:
            print("viewLocalVectorPointData not implemented, using global!")
            return self.viewGlobalVectorPointData(cmdFile,nSpace,q,ckey,name=name,
                                                  storeMeshData=storeMeshData,
                                                  figureNumber=figureNumber,title=title)
#         nPoints_element = q['x'].shape[1]
#         if nPoints_element <= nSpace:
#             print """
# Warning! viewScalarPointData nPoints_element=%s < %s too small for useLocal, using global interp""" % (nPoints_element,
#                                                                                                        nSpace+1)
#             return self.viewGlobalScalarPointData(cmdFile,nSpace,q,ckey,name=name,
#                                                   storeMeshData=storeMeshData,
#                                                   figureNumber=figureNumber,title=title)
#         return self.viewLocalScalarPointData(cmdFile,nSpace,q,ckey,name=name,
#                                              storeMeshData=storeMeshData,
#                                              figureNumber=figureNumber,title=title)

    def viewGlobalScalarPointData(self,cmdFile,nSpace,q,ckey,name=None,
                                  storeMeshData=True,figureNumber=1,title=None):
        """
        input scalar variable and coordinates stored in dictionary
           q['x'], q[ckey]
        respectively, generate global continuous interpolant
        should work for q, ebq_global, and ebqe quadrature dictionaries

        uses delaunay triangulation in 2d and 3d

        scalar data is stored in
            name_q

        if storeMeshData = True, writes out
            name_x_q    -- point data
            tri_name_q  -- Delaunay representation (2d,3d)

        returns number of figures actually plotted
        """
        nplotted = 0
        ###simple visualization commands (%s --> name)
        #1d
        cmd1dData = """
[x_tmp,i_tmp] = sort(%s_x_q(:,1)); %s_x_q = %s_x_q(i_tmp); %s_q = %s_q(i_tmp);
"""
        cmd1dView = """
figure(%i) ; plot(%s_x_q(:,1),%s_q) ;  title('%s');
"""
        #2d
        cmd2dData = """
tri_%s_q = delaunay(%s_x_q(:,1),%s_x_q(:,2));
"""
        cmd2dView = """
figure(%i) ; trimesh(tri_%s_q,%s_x_q(:,1),%s_x_q(:,2),%s_q); title('%s');
%%also could be used
%%trisurf(tri_name_q,name_x_q(:,1),name_x_q(:,2),name_q);
"""
        cmd2dGrid = """
nx = %i; ny = %i;
XYZ = [min(%s_x_q(:,1)) max(%s_x_q(:,1)) ; min(%s_x_q(:,2)) max(%s_x_q(:,2))];
dx = (XYZ(1,2)-XYZ(1,1))/nx; dy = (XYZ(2,2)-XYZ(2,1))/ny;
[%s_qxg,%s_qyg] = meshgrid(XYZ(1,1):dx:XYZ(1,2),XYZ(2,1):dy:XYZ(2,2));
%s_qg = griddata(%s_x_q(:,1),%s_x_q(:,2),%s_q,%s_qxg,%s_qyg);
"""
        #3d
        cmd3dData = """
tri_%s_q = delaunay3(%s_x_q(:,1),%s_x_q(:,2),%s_x_q(:,3));
"""
        cmd3dView = """
%%Warning not very good right now
%%figure(%i) ; tetramesh(tri_%s_q,%s_x_q);  title('%s');
"""
        cmd3dGrid = """
nx = %i; ny = %i; nz = %i ;
XYZ = [min(%s_x_q(:,1)) max(%s_x_q(:,1)) ; min(%s_x_q(:,2)) max(%s_x_q(:,2)) ; min(%s_x_q(:,3)) max(%s_x_q(:,3))];
dx = (XYZ(1,2)-XYZ(1,1))/nx; dy = (XYZ(2,2)-XYZ(2,1))/ny; dz = (XYZ(3,2)-XYZ(3,1))/nz;
[%s_qxg,%s_qyg,%s_qzg] = meshgrid(XYZ(1,1):dx:XYZ(1,2),XYZ(2,1):dy:XYZ(2,2),XYZ(3,1):dz:XYZ(3,2));
%s_qg = griddata3(%s_x_q(:,1),%s_x_q(:,2),%s_x_q(:,3),%s_q,%s_qxg,%s_qyg,%s_qzg);
"""

        ###
        if name is None:
            name = ckey[0];
            for i in range(len(ckey)-1):
                name += "_%s" % ckey[1+i]
        if title is None:
            title = name
        assert ckey in q, " ckey = %s missing from q " % ckey
        assert len(q[ckey].shape) == 2, " q[%s].shape= %s should be ( , ) " % (ckey,q[ckey].shape)

        if storeMeshData:
            self.storePointMeshData(cmdFile,q['x'],name)
            #cmdFile.write("%s_x_q = [ ... \n" % name)
            #for eN in range(q['x'].shape[0]):
            #    for k in range(q['x'].shape[1]):
            #        cmdFile.write("%g %g %g \n" % (q['x'][eN,k,0],q['x'][eN,k,1],q['x'][eN,k,2]))
            #cmdFile.write("];")
        #

        cmdFile.write("%s_q = [ ... \n" % name)
        for eN in range(q[ckey].shape[0]): #ebq_global or ebqe would work
            for k in range(q[ckey].shape[1]):
                cmdFile.write("%g \n" % q[ckey][eN,k])
        cmdFile.write("];");

        if nSpace == 1:
            cmd = cmd1dData % (name,name,name,name,name)
            cmdFile.write(cmd)
            cmd = cmd1dView % (figureNumber,name,name,title)
            cmdFile.write(cmd)
            nplotted = 1
        elif nSpace == 2:
            cmd = cmd2dData % (name,name,name)
            cmdFile.write(cmd)
            cmd = cmd2dView % (figureNumber,name,name,name,name,title)
            cmdFile.write(cmd)
            cmd = cmd2dGrid % (self.ngrid[0],self.ngrid[1],
                               name,name,name,name,
                               name,name,
                               name,name,name,name,name,name)

            cmdFile.write(cmd)
            nplotted = 1
        else:
            cmd = cmd3dData % (name,name,name,name)
            cmdFile.write(cmd)
            cmd = cmd3dView % (figureNumber,name,name,title)
            cmdFile.write(cmd)
            cmd = cmd3dGrid % (self.ngrid[0],self.ngrid[1],self.ngrid[2],
                               name,name,name,name,name,name,
                               name,name,name,
                               name,name,name,name,name,name,name,name)
            cmdFile.write(cmd)
            nplotted = 0 #tetramesh no good right now

        #
        return nplotted

    def viewGlobalVectorPointData(self,cmdFile,nSpace,q,ckey,name=None,
                                  storeMeshData=True,figureNumber=1,title=None):
        """
        input vector variable and coordinates stored in dictionary
           q['x'], q[ckey]
        respectively, generate global continuous interpolant
        should work for q, ebq_global, and ebqe quadrature dictionaries

        uses delaunay triangulation in 2d and 3d

        vector data is stored in
            name_q

        if storeMeshData = True, writes out
            name_x_q    -- point data
            tri_name_q  -- Delaunay representation (2d,3d)

        returns number of figures actually generated
        """
        nplotted = 0
        ###simple visualization commands (%s --> name)
        #1d
        cmd1dData = """
[x_tmp,i_tmp] = sort(%s_x_q(:,1)); %s_x_q = %s_x_q(i_tmp); %s_q = %s_q(i_tmp);
"""
        cmd1dView = """
figure(%i) ; plot(%s_x_q(:,1),%s_q) ;  title('%s');
"""
        #2d
        cmd2dData = """
tri_%s_q = delaunay(%s_x_q(:,1),%s_x_q(:,2));
"""
        cmd2dView = """
%s_q_mag = (%s_q(:,1).^2 + %s_q(:,2).^2).^(0.5);
figure(%i) ; quiver(%s_x_q(:,1),%s_x_q(:,2),%s_q(:,1),%s_q(:,2));title('%s');
%%could also use
%% trimesh(tri_name_q,name_x_q(:,1),name_x_q(:,2),name_q_mag);
"""
        cmd2dGrid = """
nx = %i; ny = %i;
XYZ = [min(%s_x_q(:,1)) max(%s_x_q(:,1)) ; min(%s_x_q(:,2)) max(%s_x_q(:,2))];
dx = (XYZ(1,2)-XYZ(1,1))/nx; dy = (XYZ(2,2)-XYZ(2,1))/ny;
[%s_qxg,%s_qyg] = meshgrid(XYZ(1,1):dx:XYZ(1,2),XYZ(2,1):dy:XYZ(2,2));
%s_qg_x = griddata(%s_x_q(:,1),%s_x_q(:,2),%s_q(:,1),%s_qxg,%s_qyg);
%s_qg_y = griddata(%s_x_q(:,1),%s_x_q(:,2),%s_q(:,2),%s_qxg,%s_qyg);
"""
        #3d
        cmd3dData = """
tri_%s_q = delaunay3(%s_x_q(:,1),%s_x_q(:,2),%s_x_q(:,3));
"""
        cmd3dView = """
%s_q_mag = (%s_q(:,1).^2 + %s_q(:,2).^2 + %s_q(:,3).^2).^(0.5);
figure(%i) ; quiver3(%s_x_q(:,1),%s_x_q(:,2),%s_x_q(:,3),%s_q(:,1),%s_q(:,2),%s_q(:,3));title('%s');
%%Warning not very good right now
%% tetramesh(tri_name_q,name_x_q);
"""
        cmd3dGrid = """
nx = %i; ny = %i; nz = %i ;
XYZ = [min(%s_x_q(:,1)) max(%s_x_q(:,1)) ; min(%s_x_q(:,2)) max(%s_x_q(:,2)) ; min(%s_x_q(:,3)) max(%s_x_q(:,3))];
dx = (XYZ(1,2)-XYZ(1,1))/nx; dy = (XYZ(2,2)-XYZ(2,1))/ny; dz = (XYZ(3,2)-XYZ(3,1))/nz;
[%s_qxg,%s_qyg,%s_qzg] = meshgrid(XYZ(1,1):dx:XYZ(1,2),XYZ(2,1):dy:XYZ(2,2),XYZ(3,1):dz:XYZ(3,2));
%s_qg_x = griddata3(%s_x_q(:,1),%s_x_q(:,2),%s_x_q(:,3),%s_q(:,1),%s_qxg,%s_qyg,%s_qzg);
%s_qg_y = griddata3(%s_x_q(:,1),%s_x_q(:,2),%s_x_q(:,3),%s_q(:,2),%s_qxg,%s_qyg,%s_qzg);
%s_qg_z = griddata3(%s_x_q(:,1),%s_x_q(:,2),%s_x_q(:,3),%s_q(:,3),%s_qxg,%s_qyg,%s_qzg);
"""

        ###
        if name is None:
            name = ckey[0];
            for i in range(len(ckey)-1):
                name += "_%s" % ckey[1+i]
        if title is None:
            title = name
        assert ckey in q, " ckey = %s missing from q " % ckey
        assert len(q[ckey].shape) == 3, " q[%s].shape= %s should be ( , , ) " % (ckey,q[ckey].shape)

        if storeMeshData:
            self.storePointMeshData(cmdFile,q['x'],name)

        #

        cmdFile.write("%s_q = [ ... \n" % name)
        for eN in range(q[ckey].shape[0]): #ebq_global or ebqe would work
            for k in range(q[ckey].shape[1]):
                for I in range(q[ckey].shape[2]):
                    cmdFile.write(" %g " % q[ckey][eN,k,I])
                cmdFile.write(" ; \n ")
        cmdFile.write("];");

        if nSpace == 1:
            cmd = cmd1dData % (name,name,name,name,name)
            cmdFile.write(cmd)
            cmd = cmd1dView % (figureNumber,name,name,title)
            cmdFile.write(cmd)
            nplotted = 1
        elif nSpace == 2:
            cmd = cmd2dData % (name,name,name)
            cmdFile.write(cmd)
            cmd = cmd2dView % (name,name,name,
                               figureNumber,name,name,name,name,title)
            cmdFile.write(cmd)
            cmd = cmd2dGrid % (self.ngrid[0],self.ngrid[1],
                               name,name,name,name,
                               name,name,
                               name,name,name,name,name,name,
                               name,name,name,name,name,name)

            cmdFile.write(cmd)
            nplotted = 1
        else:
            cmd = cmd3dData % (name,name,name,name)
            cmdFile.write(cmd)
            cmd = cmd3dView % (name,name,name,name,
                               figureNumber,name,name,name,name,name,name,title)
            cmdFile.write(cmd)
            cmd = cmd3dGrid % (self.ngrid[0],self.ngrid[1],self.ngrid[2],
                               name,name,name,name,name,name,
                               name,name,name,
                               name,name,name,name,name,name,name,name,
                               name,name,name,name,name,name,name,name,
                               name,name,name,name,name,name,name,name)
            cmdFile.write(cmd)
            nplotted = 1

        #
        return nplotted

    def viewLocalScalarPointData(self,cmdFile,nSpace,q,ckey,name=None,
                                 storeMeshData=True,figureNumber=1,title=None):
        """
        input scalar variable and coordinates stored in dictionary
           q['x'], q[ckey]
        respectively, generate local  interpolant that is element-wise continuous
        should work for q, ebq_global, and ebqe quadrature dictionaries

        uses delaunay triangulation in 2d and 3d on each element (stored in cell array)

        scalar data is stored in
            name_q

        if storeMeshData = True, writes out
            name_x_q    -- point data
            tri_name_q  -- Delaunay representation (2d,3d)

        returns number of figures actually plotted
        """
        #1d
        cmd1dData = """
nElements_global = %i; nPoints_element = %i;
tri_%s_q = [];
for eN = 1:nElements_global
  [x_tmp,i_tmp] = sort(%s_x_q(nPoints_element*(eN-1)+1:nPoints_element*eN,1));
  tmp = [];
  for j = 1:nPoints_element-1
    tmp = [tmp ; i_tmp(j) i_tmp(j+1)];
  end
  tri_%s_q = [tri_%s_q ; nPoints_element*(eN-1) + tmp];
end
%s_q_loc = %s_x_q; %s_q_loc(:,2) = %s_q;
"""
        cmd1dView = """
figure(%i) ; patch('vertices',%s_q_loc,'faces',tri_%s_q,'FaceColor','none','EdgeColor','black'); title('%s');
"""
        #2d
        cmd2dData = """
nElements_global = %i; nPoints_element = %i;
tri_%s_q = [];
for eN = 1:nElements_global
  tmp = delaunay(%s_x_q(nPoints_element*(eN-1)+1:nPoints_element*eN,1),...
                 %s_x_q(nPoints_element*(eN-1)+1:nPoints_element*eN,2));
  tri_%s_q = [tri_%s_q ; nPoints_element*(eN-1)+tmp];
end
%s_q_loc = %s_x_q; %s_q_loc(:,3) = %s_q;

"""
        cmd2dView = """
figure(%i) ; patch('vertices',%s_q_loc,'faces',tri_%s_q,'FaceVertexCData',%s_q,'FaceColor','interp','EdgeColor','none'); title('%s');
"""
        cmd2dGrid = """
nx = %i; ny = %i;
XYZ = [min(%s_x_q(:,1)) max(%s_x_q(:,1)) ; min(%s_x_q(:,2)) max(%s_x_q(:,2))];
dx = (XYZ(1,2)-XYZ(1,1))/nx; dy = (XYZ(2,2)-XYZ(2,1))/ny;
[%s_qxg,%s_qyg] = meshgrid(XYZ(1,1):dx:XYZ(1,2),XYZ(2,1):dy:XYZ(2,2));
%s_qg = griddata(%s_x_q(:,1),%s_x_q(:,2),%s_q,%s_qxg,%s_qyg);
"""
        #3d
        cmd3dData = """
nElements_global = %i; nPoints_element = %i;
tri_%s_q = [];
for eN = 1:nElements_global
  tmp = delaunay3(%s_x_q(nPoints_element*(eN-1)+1:nPoints_element*eN,1),...
                  %s_x_q(nPoints_element*(eN-1)+1:nPoints_element*eN,2),...
                  %s_x_q(nPoints_element*(eN-1)+1:nPoints_element*eN,3));
  tri_%s_q = [tri_%s_q ; nPoints_element*(eN-1)+tmp];
end
%s_q_loc = %s_x_q; %s_q_loc(:,3) = %s_q;

"""
        cmd3dView = """
%%good luck
figure(%i); patch('vertices',%s_x_q,'faces',tri_%s_q,'FaceVertexCData',%s_q,'FaceColor','none','EdgeColor','interp'); title('%s')
"""
        cmd3dGrid = """
nx = %i; ny = %i; nz = %i ;
XYZ = [min(%s_x_q(:,1)) max(%s_x_q(:,1)) ; min(%s_x_q(:,2)) max(%s_x_q(:,2)) ; min(%s_x_q(:,3)) max(%s_x_q(:,3))];
dx = (XYZ(1,2)-XYZ(1,1))/nx; dy = (XYZ(2,2)-XYZ(2,1))/ny; dz = (XYZ(3,2)-XYZ(3,1))/nz;
[%s_qxg,%s_qyg,%s_qzg] = meshgrid(XYZ(1,1):dx:XYZ(1,2),XYZ(2,1):dy:XYZ(2,2),XYZ(3,1):dz:XYZ(3,2));
%s_qg = griddata3(%s_x_q(:,1),%s_x_q(:,2),%s_x_q(:,3),%s_q,%s_qxg,%s_qyg,%s_qzg);
"""
        nplotted = 0
        if name is None:
            name = ckey[0];
            for i in range(len(ckey)-1):
                name += "_%s" % ckey[1+i]
        if title is None:
            title = name
        assert ckey in q, " ckey = %s missing from q " % ckey
        assert len(q[ckey].shape) == 2, " q[%s].shape= %s should be ( , ) " % (ckey,q[ckey].shape)

        nElements_global = q['x'].shape[0]; nPoints_element = q['x'].shape[1];
        if storeMeshData:
            self.storePointMeshData(cmdFile,q['x'],name)

        cmdFile.write("%s_q = [ ... \n" % name)
        for eN in range(q[ckey].shape[0]): #ebq_global or ebqe would work
            for k in range(q[ckey].shape[1]):
                cmdFile.write("%g \n" % q[ckey][eN,k])
        cmdFile.write("];");

        if nSpace == 1:
            cmd = cmd1dData % (nElements_global,nPoints_element,
                               name,
                               name,
                               name,name,
                               name,name,name,name)

            cmdFile.write(cmd)
            cmd = cmd1dView % (figureNumber,name,name,title)
            cmdFile.write(cmd)
            nplotted = 1

        elif nSpace == 2:
            cmd = cmd2dData % (nElements_global,nPoints_element,
                               name,
                               name,name,
                               name,name,
                               name,name,name,name)
            cmdFile.write(cmd)
            cmd = cmd2dView % (figureNumber,name,name,name,title)
            cmdFile.write(cmd)
            cmd = cmd2dGrid % (self.ngrid[0],self.ngrid[1],
                               name,name,name,name,
                               name,name,
                               name,name,name,name,name,name)

            cmdFile.write(cmd)
            nplotted = 1

        else:
            cmd = cmd3dData % (nElements_global,nPoints_element,
                               name,
                               name,name,name,
                               name,name,
                               name,name,name,name)
            cmdFile.write(cmd)
            cmd = cmd3dView % (figureNumber,name,name,name,title)
            cmdFile.write(cmd)
            cmd = cmd3dGrid % (self.ngrid[0],self.ngrid[1],self.ngrid[2],
                               name,name,name,name,name,name,
                               name,name,name,
                               name,name,name,name,name,name,name,name)
            cmdFile.write(cmd)
            nplotted = 1
        return nplotted
    #
    def viewScalarAnalyticalFunction(self,cmdFile,nSpace,f,t,x,elementNodesConnectivity=None,
                                     name='exact',storeMeshData=True,figureNumber=1,title=None):
        """
        input scalar analytical function f(x,t) and array of points

        respectively, generate global continuous interpolant

        uses delaunay triangulation in 2d and 3d if element - node triangulation not already defined

        scalar data is stored in
            name_ex

        if storeMeshData = True, writes out
            name_x_ex   -- point data
            tri_name_ex -- element-node representation (2d,3d)

        returns number of figures actually plotted
        """
        nplotted = 0
        ###simple visualization commands (%s --> name)
        #1d
        cmd1dData = """
[x_tmp,i_tmp] = sort(%s_x_ex(:,1)); %s_x_ex = %s_x_ex(i_tmp); %s_ex = %s_ex(i_tmp);
"""
        cmd1dView = """
figure(%i) ; plot(%s_x_ex(:,1),%s_ex) ;  title('%s');
"""
        #2d
        #if does not have element-node connectivity already
        cmd2dData = """
tri_%s_ex = delaunay(%s_x_ex(:,1),%s_x_ex(:,2));
"""
        cmd2dView = """
figure(%i) ; trimesh(tri_%s_ex,%s_x_ex(:,1),%s_x_ex(:,2),%s_ex); title('%s');
%%also could be used
%%trisurf(tri_name_ex,name_x_ex(:,1),name_x_ex(:,2),name_ex);
"""
        cmd2dGrid = """
nx = %i; ny = %i;
XYZ = [min(%s_x_ex(:,1)) max(%s_x_ex(:,1)) ; min(%s_x_ex(:,2)) max(%s_x_ex(:,2))];
dx = (XYZ(1,2)-XYZ(1,1))/nx; dy = (XYZ(2,2)-XYZ(2,1))/ny;
[%s_qxg,%s_qyg] = meshgrid(XYZ(1,1):dx:XYZ(1,2),XYZ(2,1):dy:XYZ(2,2));
%s_qg = griddata(%s_x_ex(:,1),%s_x_ex(:,2),%s_ex,%s_qxg,%s_qyg);
"""
        #3d
        #if does not have element-node connectivity already
        cmd3dData = """
tri_%s_ex = delaunay3(%s_x_ex(:,1),%s_x_ex(:,2),%s_x_ex(:,3));
"""
        cmd3dView = """
%%Warning not very good right now
%%figure(%i) ; tetramesh(tri_%s_ex,%s_x_ex);  title('%s');
"""
        cmd3dGrid = """
nx = %i; ny = %i; nz = %i ;
XYZ = [min(%s_x_ex(:,1)) max(%s_x_ex(:,1)) ; min(%s_x_ex(:,2)) max(%s_x_ex(:,2)) ; min(%s_x_ex(:,3)) max(%s_x_ex(:,3))];
dx = (XYZ(1,2)-XYZ(1,1))/nx; dy = (XYZ(2,2)-XYZ(2,1))/ny; dz = (XYZ(3,2)-XYZ(3,1))/nz;
[%s_qxg,%s_qyg,%s_qzg] = meshgrid(XYZ(1,1):dx:XYZ(1,2),XYZ(2,1):dy:XYZ(2,2),XYZ(3,1):dz:XYZ(3,2));
%s_qg = griddata3(%s_x_ex(:,1),%s_x_ex(:,2),%s_x_ex(:,3),%s_ex,%s_qxg,%s_qyg,%s_qzg);
"""

        ###
        if title is None:
            title = name

        nPoints = 1
        for i in range(len(x.shape)-1):
            nPoints *= x.shape[i]
        if storeMeshData:
            cmdFile.write("%s_x_ex = [ ... \n" % name)
            for i in range(nPoints):
                cmdFile.write("%g %g %g \n" % (x.flat[3*i+0],x.flat[3*i+1],x.flat[3*i+2]))
            #
            cmdFile.write("];")
        #
        if elementNodesConnectivity is not None:
            assert len(elementNodesConnectivity.shape) == 2
            cmdFile.write("tri_%s_ex = [ ... \n" % name)
            for eN in range(elementNodesConnectivity.shape[0]):
                for nN in range(elementNodesConnectivity.shape[1]):
                    cmdFile.write(" %i " % (elementNodesConnectivity[eN,nN]+1))
                cmdFile.write("; \n ")
            #
            cmdFile.write("];")
        #
        cmdFile.write("%s_ex = [ ... \n" % name)
        for i in range(nPoints):
            uex = f(x.flat[3*i:3*(i+1)],t)
            cmdFile.write("%g \n" % uex)
        cmdFile.write("];");


        if nSpace == 1:
            cmd = cmd1dData % (name,name,name,name,name)
            cmdFile.write(cmd)
            cmd = cmd1dView % (figureNumber,name,name,title)
            cmdFile.write(cmd)
            nplotted = 1
        elif nSpace == 2:
            if elementNodesConnectivity is None:
                cmd = cmd2dData % (name,name,name)
                cmdFile.write(cmd)
            cmd = cmd2dView % (figureNumber,name,name,name,name,title)
            cmdFile.write(cmd)
            cmd = cmd2dGrid % (self.ngrid[0],self.ngrid[1],
                               name,name,name,name,
                               name,name,
                               name,name,name,name,name,name)

            cmdFile.write(cmd)
            nplotted = 1
        else:
            if elementNodesConnectivity is None:
                cmd = cmd3dData % (name,name,name,name)
                cmdFile.write(cmd)
            cmd = cmd3dView % (figureNumber,name,name,title)
            cmdFile.write(cmd)
            cmd = cmd3dGrid % (self.ngrid[0],self.ngrid[1],self.ngrid[2],
                               name,name,name,name,name,name,
                               name,name,name,
                               name,name,name,name,name,name,name,name)
            cmdFile.write(cmd)
            nplotted = 0 #tetramesh no good right now

        #
        return nplotted

    def viewVectorAnalyticalFunction(self,cmdFile,nSpace,f,t,x,elementNodesConnectivity=None,
                                     name='exact',storeMeshData=True,figureNumber=1,title=None):
        """
        input vector analytical function f(x,t) and array of points

        respectively, generate global continuous interpolant

        uses delaunay triangulation in 2d and 3d if element - node triangulation not already defined

        vector data is stored in
            name_ex

        if storeMeshData = True, writes out
            name_x_ex   -- point data
            tri_name_ex -- element-node representation (2d,3d)

        returns number of figures actually plotted
        """
        nplotted = 0
        ###simple visualization commands (%s --> name)
        #1d
        cmd1dData = """
[x_tmp,i_tmp] = sort(%s_x_ex(:,1)); %s_x_ex = %s_x_ex(i_tmp); %s_ex = %s_ex(i_tmp);
"""
        cmd1dView = """
figure(%i) ; plot(%s_x_ex(:,1),%s_ex) ;  title('%s');
"""
        #2d
        #if does not have element-node connectivity already
        cmd2dData = """
tri_%s_ex = delaunay(%s_x_ex(:,1),%s_x_ex(:,2));
"""
        cmd2dView = """
%s_ex_mag = (%s_ex(:,1).^2 + %s_ex(:,2).^2).^(0.5);
figure(%i) ; quiver(%s_x_ex(:,1),%s_x_ex(:,2),%s_ex(:,1),%s_ex(:,2));title('%s');
%%could also use
%% trimesh(tri_name_ex,name_x_ex(:,1),name_x_ex(:,2),name_ex_mag);
"""
        cmd2dGrid = """
nx = %i; ny = %i;
XYZ = [min(%s_x_ex(:,1)) max(%s_x_ex(:,1)) ; min(%s_x_ex(:,2)) max(%s_x_ex(:,2))];
dx = (XYZ(1,2)-XYZ(1,1))/nx; dy = (XYZ(2,2)-XYZ(2,1))/ny;
[%s_qxg,%s_qyg] = meshgrid(XYZ(1,1):dx:XYZ(1,2),XYZ(2,1):dy:XYZ(2,2));
%s_qg_x = griddata(%s_x_ex(:,1),%s_x_ex(:,2),%s_ex(:,1),%s_qxg,%s_qyg);
%s_qg_y = griddata(%s_x_ex(:,1),%s_x_ex(:,2),%s_ex(:,2),%s_qxg,%s_qyg);
"""
        #3d
        #if does not have element-node connectivity already
        cmd3dData = """
tri_%s_ex = delaunay3(%s_x_ex(:,1),%s_x_ex(:,2),%s_x_ex(:,3));
"""
        cmd3dView = """
%s_ex_mag = (%s_ex(:,1).^2 + %s_ex(:,2).^2 + %s_ex(:,3).^2).^(0.5);
figure(%i) ; quiver3(%s_x_ex(:,1),%s_x_ex(:,2),%s_x_ex(:,3),%s_ex(:,1),%s_ex(:,2),%s_ex(:,3));title('%s');
%%Warning not very good right now
%% tetramesh(tri_name_ex,name_x_ex);
"""
        cmd3dGrid = """
nx = %i; ny = %i; nz = %i ;
XYZ = [min(%s_x_ex(:,1)) max(%s_x_ex(:,1)) ; min(%s_x_ex(:,2)) max(%s_x_ex(:,2)) ; min(%s_x_ex(:,3)) max(%s_x_ex(:,3))];
dx = (XYZ(1,2)-XYZ(1,1))/nx; dy = (XYZ(2,2)-XYZ(2,1))/ny; dz = (XYZ(3,2)-XYZ(3,1))/nz;
[%s_qxg,%s_qyg,%s_qzg] = meshgrid(XYZ(1,1):dx:XYZ(1,2),XYZ(2,1):dy:XYZ(2,2),XYZ(3,1):dz:XYZ(3,2));
%s_qg_x = griddata3(%s_x_ex(:,1),%s_x_ex(:,2),%s_x_ex(:,3),%s_ex(:,1),%s_qxg,%s_qyg,%s_qzg);
%s_qg_y = griddata3(%s_x_ex(:,1),%s_x_ex(:,2),%s_x_ex(:,3),%s_ex(:,2),%s_qxg,%s_qyg,%s_qzg);
%s_qg_z = griddata3(%s_x_ex(:,1),%s_x_ex(:,2),%s_x_ex(:,3),%s_ex(:,3),%s_qxg,%s_qyg,%s_qzg);
"""


        ###
        if title is None:
            title = name

        nPoints = 1
        for i in range(len(x.shape)-1):
            nPoints *= x.shape[i]

        if storeMeshData:
            cmdFile.write("%s_x_ex = [ ... \n" % name)
            for i in range(nPoints):
                cmdFile.write("%g %g %g \n" % (x.flat[3*i+0],x.flat[3*i+1],x.flat[3*i+2]))
            #
            cmdFile.write("];")
        #
        if elementNodesConnectivity is not None:
            assert len(elementNodesConnectivity.shape) == 2
            cmdFile.write("tri_%s_ex = [ ... \n" % name)
            for eN in range(elementNodesConnectivity.shape[0]):
                for nN in range(elementNodesConnectivity.shape[1]):
                    cmdFile.write(" %i " % (elementNodesConnectivity[eN,nN]+1))
                cmdFile.write("; \n ")
            #
            cmdFile.write("];")
        #
        cmdFile.write("%s_ex = [ ... \n" % name)
        for i in range(nPoints):
            uex = f(x.flat[3*i:3*(i+1)],t)
            for I in range(nSpace):
                cmdFile.write(" %g " % uex[I])
            cmdFile.write("; \n")
        cmdFile.write("];");


        if nSpace == 1:
            cmd = cmd1dData % (name,name,name,name,name)
            cmdFile.write(cmd)
            cmd = cmd1dView % (figureNumber,name,name,title)
            cmdFile.write(cmd)
            nplotted = 1
        elif nSpace == 2:
            if elementNodesConnectivity is None:
                cmd = cmd2dData % (name,name,name)
                cmdFile.write(cmd)
            cmd = cmd2dView % (name,name,name,
                               figureNumber,name,name,name,name,title)
            cmdFile.write(cmd)
            cmd = cmd2dGrid % (self.ngrid[0],self.ngrid[1],
                               name,name,name,name,
                               name,name,
                               name,name,name,name,name,name,
                               name,name,name,name,name,name)
            cmdFile.write(cmd)
            nplotted = 1
        else:
            if elementNodesConnectivity is None:
                cmd = cmd3dData % (name,name,name,name)
                cmdFile.write(cmd)
            cmd = cmd3dView % (name,name,name,name,
                               figureNumber,name,name,name,name,name,name,title)
            cmdFile.write(cmd)
            cmd = cmd3dGrid % (self.ngrid[0],self.ngrid[1],self.ngrid[2],
                               name,name,name,name,name,name,
                               name,name,name,
                               name,name,name,name,name,name,name,name,
                               name,name,name,name,name,name,name,name,
                               name,name,name,name,name,name,name,name)
            cmdFile.write(cmd)
            nplotted = 1 #

        #
        return nplotted

    def viewScalar_LagrangeC0P1(self,cmdFile,nSpace,nodeArray,elementNodesArray,u_dof,
                                name="u",storeMeshData=True,figureNumber=1,title=None):
        # """given C0 P1 function with nodal Lagrange representation generate a
        # matlab representation that is as faithful as possible to the
        # actual finite element function structure

        # C0-P1 output: mesh vertices and element connectivity degrees of freedom at vertices

        # scalar data is stored in name

        # if storeMeshData = True, writes out
        
        # name_x      -- mesh vertices
        
        # tri_name_v  -- element-node representation

        # returns number of figures actually plotted
        
        # """
        nplotted = 0
        ###simple visualization commands (%s --> name)
        #1d
        cmd1dData = """
[x_tmp,i_tmp] = sort(%s_x(:,1)); %s_x = %s_x(i_tmp); %s = %s(i_tmp);
"""
        cmd1dView = """
figure(%i) ; plot(%s_x(:,1),%s) ;  title('%s');
"""
        #2d
        cmd2dData = """ \n """

        cmd2dView = """
figure(%i) ; trimesh(tri_%s,%s_x(:,1),%s_x(:,2),%s); title('%s');
%%also could be used
%%trisurf(tri_name,name_x(:,1),name_x(:,2),name);
%%tmp = name_x; tmp(:,3) = name;
%%patch('vertices',tmp,'faces',tri_name,'FaceVertexCData',name,'FaceColor','interp','EdgeColor','none')
"""

        cmd2dGrid = """
nx = %i; ny = %i;
XYZ = [min(%s_x(:,1)) max(%s_x(:,1)) ; min(%s_x(:,2)) max(%s_x(:,2))];
dx = (XYZ(1,2)-XYZ(1,1))/nx; dy = (XYZ(2,2)-XYZ(2,1))/ny;
[%s_xg,%s_yg] = meshgrid(XYZ(1,1):dx:XYZ(1,2),XYZ(2,1):dy:XYZ(2,2));
%s_g = griddata(%s_x(:,1),%s_x(:,2),%s,%s_xg,%s_yg);
"""
        #3d
        cmd3dData = """ \n """

        cmd3dView = """
%%figure out reasonable default using tetramesh or patch
"""
        cmd3dGrid = """
nx = %i; ny = %i; nz = %i ;
XYZ = [min(%s_x(:,1)) max(%s_x(:,1)) ; min(%s_x(:,2)) max(%s_x(:,2)) ; min(%s_x(:,3)) max(%s_x(:,3))];
dx = (XYZ(1,2)-XYZ(1,1))/nx; dy = (XYZ(2,2)-XYZ(2,1))/ny; dz = (XYZ(3,2)-XYZ(3,1))/nz;
[%s_xg,%s_yg,%s_zg] = meshgrid(XYZ(1,1):dx:XYZ(1,2),XYZ(2,1):dy:XYZ(2,2),XYZ(3,1):dz:XYZ(3,2));
%s_g = griddata3(%s_x(:,1),%s_x(:,2),%s_x(:,3),%s,%s_xg,%s_yg,%s_zg);
"""
        if title is None:
            title = name

        nNodes_global = nodeArray.shape[0]; nElements_global = elementNodesArray.shape[0]
        nNodes_element= elementNodesArray.shape[1]
        assert nNodes_element == nSpace+1, "affine simplicial geometry only"

        if storeMeshData:
            cmdFile.write("%s_x = [ ... \n" % name)
            for nN in range(nNodes_global):
                cmdFile.write("%g %g %g \n" % (nodeArray[nN,0],nodeArray[nN,1],nodeArray[nN,2]))
            #
            cmdFile.write("];")
        #
        cmdFile.write("tri_%s = [ ... \n" % name)
        for eN in range(nElements_global):
            for nN in range(nNodes_element):
                cmdFile.write(" %i " % (elementNodesArray[eN,nN]+1))
            cmdFile.write("; \n ")
        #
        cmdFile.write("];")
        #

        cmdFile.write("%s = [ ... \n" % name)
        for i in range(u_dof.shape[0]):
            cmdFile.write("%g \n" % u_dof[i])
        cmdFile.write("];");

        if nSpace == 1:
            cmd = cmd1dData % (name,name,name,name,name)
            cmdFile.write(cmd)
            cmd = cmd1dView % (figureNumber,name,name,title)
            cmdFile.write(cmd)
            nplotted = 1

        elif nSpace == 2:
            #nothing to be done with mesh-data representation here
            cmd = cmd2dView % (figureNumber,name,name,name,name,title)
            cmdFile.write(cmd)
            cmd = cmd2dGrid % (self.ngrid[0],self.ngrid[1],
                               name,name,name,name,
                               name,name,
                               name,name,name,name,name,name)

            cmdFile.write(cmd)
            nplotted = 1
        else:
            #nothing to be done with mesh-data representation here
            #cmd = cmd3dView % (figureNumber,name,name,title)
            #cmdFile.write(cmd)
            cmd = cmd3dGrid % (self.ngrid[0],self.ngrid[1],self.ngrid[2],
                               name,name,name,name,name,name,
                               name,name,name,
                               name,name,name,name,name,name,name,name)
            cmdFile.write(cmd)
            nplotted = 0 #tetramesh no good right now

        #
        return nplotted
    #
    def viewVector_LagrangeC0P1(self,cmdFile,nSpace,nodeArray,elementNodesArray,
                                u_dof,v_dof=None,w_dof=None,
                                name="velocity",storeMeshData=True,figureNumber=1,title=None):
        # """
        # given a vector valued C0 P1 function with nodal Lagrange representation
        # generate a matlab representation that is as faithful as possible to
        # the actual finite element function structure

        # Assumes 1 <= number of components <= nSpace
        # C0-P1 output:
        #    mesh vertices and element connectivity
        #    degrees of freedom at vertices

        # vector data is stored in
        #     name which is [nNodes,nCoord]

        # if storeMeshData = True, writes out
        #     name_x      -- mesh vertices
        #     tri_name_v  -- element-node representation

        # returns number of figures actually plotted

        # """
        nplotted = 0
        ###simple visualization commands (%s --> name)
        #1d
        cmd1dData = """
[x_tmp,i_tmp] = sort(%s_x(:,1)); %s_x = %s_x(i_tmp); %s = %s(i_tmp);
"""
        cmd1dView = """
figure(%i) ; plot(%s_x(:,1),%s) ;  title('%s');
"""
        #2d
        cmd2dData = """ \n """

        cmd2dView = """
%s_mag = (%s(:,1).^2 + %s(:,2).^2).^(0.5);
figure(%i) ; quiver(%s_x(:,1),%s_x(:,2),%s(:,1),%s(:,2));title('%s');
%%could also use
%% trimesh(tri_name,name_x(:,1),name_x(:,2),name_mag);
"""

        cmd2dGrid = """
nx = %i; ny = %i;
XYZ = [min(%s_x(:,1)) max(%s_x(:,1)) ; min(%s_x(:,2)) max(%s_x(:,2))];
dx = (XYZ(1,2)-XYZ(1,1))/nx; dy = (XYZ(2,2)-XYZ(2,1))/ny;
[%s_xg,%s_yg] = meshgrid(XYZ(1,1):dx:XYZ(1,2),XYZ(2,1):dy:XYZ(2,2));
%s_gxcoord = griddata(%s_x(:,1),%s_x(:,2),%s(:,1),%s_xg,%s_yg);
%s_gycoord = griddata(%s_x(:,1),%s_x(:,2),%s(:,2),%s_xg,%s_yg);
"""
        #3d
        cmd3dData = """ \n """

        cmd3dView = """
%s_mag = (%s_x(:,1).^2 + %s_x(:,2).^2 + %s_x(:,3).^2).^(0.5);
figure(%i) ; quiver3(%s_x(:,1),%s_x(:,2),%s_x(:,3),%s(:,1),%s(:,2),%s(:,3));title('%s');
%%figure out reasonable default using tetramesh or patch
"""
        cmd3dGrid = """
nx = %i; ny = %i; nz = %i ;
XYZ = [min(%s_x(:,1)) max(%s_x(:,1)) ; min(%s_x(:,2)) max(%s_x(:,2)) ; min(%s_x(:,3)) max(%s_x(:,3))];
dx = (XYZ(1,2)-XYZ(1,1))/nx; dy = (XYZ(2,2)-XYZ(2,1))/ny; dz = (XYZ(3,2)-XYZ(3,1))/nz;
[%s_xg,%s_yg,%s_zg] = meshgrid(XYZ(1,1):dx:XYZ(1,2),XYZ(2,1):dy:XYZ(2,2),XYZ(3,1):dz:XYZ(3,2));
%s_gxcoord = griddata3(%s_x(:,1),%s_x(:,2),%s_x(:,3),%s(:,1),%s_xg,%s_yg,%s_zg);
%s_gycoord = griddata3(%s_x(:,1),%s_x(:,2),%s_x(:,3),%s(:,2),%s_xg,%s_yg,%s_zg);
%s_gzcoord = griddata3(%s_x(:,1),%s_x(:,2),%s_x(:,3),%s(:,3),%s_xg,%s_yg,%s_zg);
"""
        if title is None:
            title = name
        #
        nNodes_global = nodeArray.shape[0]; nElements_global = elementNodesArray.shape[0]
        nNodes_element= elementNodesArray.shape[1]
        assert nNodes_element == nSpace+1, "affine simplicial geometry only"
        nCoords = 1
        if v_dof is not None:
            nCoords += 1
            assert v_dof.shape == u_dof.shape
        if w_dof is not None:
            nCoords += 1
            assert w_dof.shape == u_dof.shape
        assert (1 <= nCoords and nCoords <= nSpace), "nCoords= %s nSpace= %s wrong " % (nCoords,nSpace)


        if storeMeshData:
            cmdFile.write("%s_x = [ ... \n" % name)
            for nN in range(nNodes_global):
                cmdFile.write("%g %g %g \n" % (nodeArray[nN,0],nodeArray[nN,1],nodeArray[nN,2]))
            #
            cmdFile.write("];")
        #
        cmdFile.write("tri_%s = [ ... \n" % name)
        for eN in range(nElements_global):
            for nN in range(nNodes_element):
                cmdFile.write(" %i " % (elementNodesArray[eN,nN]+1))
            cmdFile.write("; \n ")
        #
        cmdFile.write("];")
        #

        cmdFile.write("%s = [ ... \n" % name)
        for i in range(u_dof.shape[0]):
            cmdFile.write(" %g " % u_dof[i])
            if nCoords > 1: cmdFile.write(" %g " % v_dof[i])
            if nCoords > 2: cmdFile.write(" %g " % w_dof[i])
            cmdFile.write(" \n ")

        cmdFile.write("];");

        if nSpace == 1:
            cmd = cmd1dData % (name,name,name,name,name)
            cmdFile.write(cmd)
            cmd = cmd1dView % (figureNumber,name,name,title)
            cmdFile.write(cmd)
            nplotted = 1

        elif nSpace == 2:
            #nothing to be done with mesh-data representation here
            cmd = cmd2dView % (name,name,name,
                               figureNumber,name,name,name,name,title)
            cmdFile.write(cmd)
            cmd = cmd2dGrid % (self.ngrid[0],self.ngrid[1],
                               name,name,name,name,
                               name,name,
                               name,name,name,name,name,name,
                               name,name,name,name,name,name)

            cmdFile.write(cmd)
            nplotted = 1
        else:
            #nothing to be done with mesh-data representation here
            cmd = cmd3dView % (name,name,name,name,
                               figureNumber,name,name,name,name,name,name,title)
            cmdFile.write(cmd)
            cmd = cmd3dGrid % (self.ngrid[0],self.ngrid[1],self.ngrid[2],
                               name,name,name,name,name,name,
                               name,name,name,
                               name,name,name,name,name,name,name,name,
                               name,name,name,name,name,name,name,name,
                               name,name,name,name,name,name,name,name)
            cmdFile.write(cmd)
            nplotted = 1 #

        #
        return nplotted
    #
    def viewScalar_LagrangeC0P2(self,cmdFile,nSpace,lagrangeNodesArray,elementNodesArray,
                                l2g,u_dof,
                                name="u",storeMeshData=True,figureNumber=1,title=None):
        # """TODO 3D
        # given C0 P2 function with nodal Lagrange representation
        # generate a matlab representation that is as faithful as possible to
        # the actual finite element function structure

        # C0-P2 output: matlab mesh vertices stored according [geometric
        # vertices, midNodesVertices] and element connectivity for
        # refined mesh degrees of freedom at all vertices 'mid-edge'
        # vertices stored in midNodesArray

        # scalar data is stored in
        #     name

        # if storeMeshData = True, writes out
        #     name_x      -- matlab mesh vertices
        #     tri_name_v  -- matlab element-node representation

        # returns number of figures actually plotted

        # """
        nplotted = 0
        ###simple visualization commands (%s --> name)
        #1d
        cmd1dData = """
[x_tmp,i_tmp] = sort(%s_x(:,1)); %s_x = %s_x(i_tmp); %s = %s(i_tmp);
"""
        cmd1dView = """
figure(%i) ; plot(%s_x(:,1),%s) ;  title('%s');
"""
        #2d
        cmd2dData = """ \n """

        cmd2dView = """
figure(%i) ; trimesh(tri_%s,%s_x(:,1),%s_x(:,2),%s); title('%s');
%%also could be used
%%trisurf(tri_name,name_x(:,1),name_x(:,2),name);
%%tmp = name_x; tmp(:,3) = name;
%%patch('vertices',tmp,'faces',tri_name,'FaceVertexCData',name,'FaceColor','interp','EdgeColor','none')
"""

        cmd2dGrid = """
nx = %i; ny = %i;
XYZ = [min(%s_x(:,1)) max(%s_x(:,1)) ; min(%s_x(:,2)) max(%s_x(:,2))];
dx = (XYZ(1,2)-XYZ(1,1))/nx; dy = (XYZ(2,2)-XYZ(2,1))/ny;
[%s_xg,%s_yg] = meshgrid(XYZ(1,1):dx:XYZ(1,2),XYZ(2,1):dy:XYZ(2,2));
%s_g = griddata(%s_x(:,1),%s_x(:,2),%s,%s_xg,%s_yg);
"""
        #3d
        cmd3dData = """ \n """

        cmd3dView = """
%%figure out reasonable default using tetramesh or patch
"""
        cmd3dGrid = """
nx = %i; ny = %i; nz = %i ;
XYZ = [min(%s_x(:,1)) max(%s_x(:,1)) ; min(%s_x(:,2)) max(%s_x(:,2)) ; min(%s_x(:,3)) max(%s_x(:,3))];
dx = (XYZ(1,2)-XYZ(1,1))/nx; dy = (XYZ(2,2)-XYZ(2,1))/ny; dz = (XYZ(3,2)-XYZ(3,1))/nz;
[%s_xg,%s_yg,%s_zg] = meshgrid(XYZ(1,1):dx:XYZ(1,2),XYZ(2,1):dy:XYZ(2,2),XYZ(3,1):dz:XYZ(3,2));
%s_g = griddata3(%s_x(:,1),%s_x(:,2),%s_x(:,3),%s,%s_xg,%s_yg,%s_zg);
"""
        if title is None:
            title = name


        nLagrangeNodes_global = lagrangeNodesArray.shape[0]; nElements_global = elementNodesArray.shape[0]
        nNodes_element= elementNodesArray.shape[1]
        assert nNodes_element == nSpace+1, "affine simplicial geometry only"
        assert u_dof.shape[0] == nLagrangeNodes_global, "u_dof len=%s but nLagrangeNodes_global=%s " % (u_dof.shape[0],
                                                                                                        nLagrangeNodes_global)
        nMidNodes_element = int((nSpace+2)*(nSpace+1)/2) - nNodes_element
        #store geometric nodes first, then "quadratic" ones
        if storeMeshData:
            cmdFile.write("%s_x = [ ... \n" % name)
            for nN in range(nLagrangeNodes_global):
                cmdFile.write("%g %g %g \n" % (lagrangeNodesArray[nN,0],lagrangeNodesArray[nN,1],lagrangeNodesArray[nN,2]))
            #

            cmdFile.write("];")
        #

        if nSpace == 1:
            cmdFile.write("tri_%s = [ ... \n" % name)
            for eN in range(nElements_global):
                #not necessary going to get positive jacobians
                cmdFile.write(" %i " % (l2g[eN,0]+1))
                cmdFile.write(" %i " % (l2g[eN,2]+1)) #use l2g?
                cmdFile.write(" \n ")
                cmdFile.write(" %i " % (l2g[eN,2]+1))#use l2g?
                cmdFile.write(" %i " % (l2g[eN,1]+1))
                cmdFile.write(" \n ")
            cmdFile.write("];")
        elif nSpace == 2:
            cmdFile.write("tri_%s = [ ... \n" % name)
            for eN in range(nElements_global):
                #four triangles per element, see FemTools Quadratic space for numbering convention
                cmdFile.write(" %i " % (l2g[eN,0]+1))
                cmdFile.write(" %i " % (l2g[eN,nNodes_element+0]+1))
                cmdFile.write(" %i " % (l2g[eN,nNodes_element+2]+1))
                cmdFile.write(" \n ")
                cmdFile.write(" %i " % (l2g[eN,nNodes_element+0]+1))
                cmdFile.write(" %i " % (l2g[eN,1]+1))
                cmdFile.write(" %i " % (l2g[eN,nNodes_element+1]+1))
                cmdFile.write(" \n ")
                cmdFile.write(" %i " % (l2g[eN,nNodes_element+1]+1))
                cmdFile.write(" %i " % (l2g[eN,2]+1))
                cmdFile.write(" %i " % (l2g[eN,nNodes_element+2]+1))
                cmdFile.write(" \n ")
                cmdFile.write(" %i " % (l2g[eN,nNodes_element+0]+1))
                cmdFile.write(" %i " % (l2g[eN,nNodes_element+1]+1))
                cmdFile.write(" %i " % (l2g[eN,nNodes_element+2]+1))
                cmdFile.write(" \n ")
            cmdFile.write("];")
        else:
            print("""3d LagrangeC0P2 not implemented yet""")
            return 0
        #
        #assumes l2g layout consistent with matlab one
        cmdFile.write("%s = [ ... \n" % name)
        for i in range(u_dof.shape[0]):
            cmdFile.write("%g \n" % u_dof[i])
        cmdFile.write("];");

        if nSpace == 1:
            cmd = cmd1dData % (name,name,name,name,name)
            cmdFile.write(cmd)
            cmd = cmd1dView % (figureNumber,name,name,title)
            cmdFile.write(cmd)
            nplotted = 1

        elif nSpace == 2:
            #nothing to be done with mesh-data representation here
            cmd = cmd2dView % (figureNumber,name,name,name,name,title)
            cmdFile.write(cmd)
            cmd = cmd2dGrid % (self.ngrid[0],self.ngrid[1],
                               name,name,name,name,
                               name,name,
                               name,name,name,name,name,name)

            cmdFile.write(cmd)
            nplotted = 1
        else:
            #nothing to be done with mesh-data representation here
            #cmd = cmd3dView % (figureNumber,name,name,title)
            #cmdFile.write(cmd)
            cmd = cmd3dGrid % (self.ngrid[0],self.ngrid[1],self.ngrid[2],
                               name,name,name,name,name,name,
                               name,name,name,
                               name,name,name,name,name,name,name,name)
            cmdFile.write(cmd)
            nplotted = 0 #tetramesh no good right now

        #
        return nplotted
    #
    def viewVector_LagrangeC0P2(self,cmdFile,nSpace,lagrangeNodesArray,elementNodesArray,
                                l2g,u_dof,v_dof=None,w_dof=None,
                                name="velocity",storeMeshData=True,figureNumber=1,title=None):
        # """
        # TODO: 3D

        # given a vector valued C0 P2 function with nodal Lagrange representation
        # generate a matlab representation that is as faithful as possible to
        # the actual finite element function structure

        # Assumes 1 <= number of components <= nSpace
        # C0-P2 output:
        #    matlab mesh vertices stored according [geometric vertices, midNodesVertices]
        #    and element connectivity for refined mesh
        #    degrees of freedom at all vertices
        #    'mid-edge' vertices stored in midNodesArray

        # vector data is stored in
        #     name which is [nNodes,nCoord]

        # if storeMeshData = True, writes out
        #     name_x      -- mesh vertices
        #     tri_name_v  -- element-node representation

        # returns number of figures actually plotted

        # """
        nplotted = 0
        ###simple visualization commands (%s --> name)
        #1d
        cmd1dData = """
[x_tmp,i_tmp] = sort(%s_x(:,1)); %s_x = %s_x(i_tmp); %s = %s(i_tmp);
"""
        cmd1dView = """
figure(%i) ; plot(%s_x(:,1),%s) ;  title('%s');
"""
        #2d
        cmd2dData = """ \n """

        cmd2dView = """
%s_mag = (%s(:,1).^2 + %s(:,2).^2).^(0.5);
figure(%i) ; quiver(%s_x(:,1),%s_x(:,2),%s(:,1),%s(:,2));title('%s');
%%could also use
%% trimesh(tri_name,name_x(:,1),name_x(:,2),name_mag);
"""

        cmd2dGrid = """
nx = %i; ny = %i;
XYZ = [min(%s_x(:,1)) max(%s_x(:,1)) ; min(%s_x(:,2)) max(%s_x(:,2))];
dx = (XYZ(1,2)-XYZ(1,1))/nx; dy = (XYZ(2,2)-XYZ(2,1))/ny;
[%s_xg,%s_yg] = meshgrid(XYZ(1,1):dx:XYZ(1,2),XYZ(2,1):dy:XYZ(2,2));
%s_gxcoord = griddata(%s_x(:,1),%s_x(:,2),%s(:,1),%s_xg,%s_yg);
%s_gycoord = griddata(%s_x(:,1),%s_x(:,2),%s(:,2),%s_xg,%s_yg);
"""
        #3d
        cmd3dData = """ \n """

        cmd3dView = """
%s_mag = (%s_x(:,1).^2 + %s_x(:,2).^2 + %s_x(:,3).^2).^(0.5);
figure(%i) ; quiver3(%s_x(:,1),%s_x(:,2),%s_x(:,3),%s(:,1),%s(:,2),%s(:,3));title('%s');
%%figure out reasonable default using tetramesh or patch
"""
        cmd3dGrid = """
nx = %i; ny = %i; nz = %i ;
XYZ = [min(%s_x(:,1)) max(%s_x(:,1)) ; min(%s_x(:,2)) max(%s_x(:,2)) ; min(%s_x(:,3)) max(%s_x(:,3))];
dx = (XYZ(1,2)-XYZ(1,1))/nx; dy = (XYZ(2,2)-XYZ(2,1))/ny; dz = (XYZ(3,2)-XYZ(3,1))/nz;
[%s_xg,%s_yg,%s_zg] = meshgrid(XYZ(1,1):dx:XYZ(1,2),XYZ(2,1):dy:XYZ(2,2),XYZ(3,1):dz:XYZ(3,2));
%s_gxcoord = griddata3(%s_x(:,1),%s_x(:,2),%s_x(:,3),%s(:,1),%s_xg,%s_yg,%s_zg);
%s_gycoord = griddata3(%s_x(:,1),%s_x(:,2),%s_x(:,3),%s(:,2),%s_xg,%s_yg,%s_zg);
%s_gzcoord = griddata3(%s_x(:,1),%s_x(:,2),%s_x(:,3),%s(:,3),%s_xg,%s_yg,%s_zg);
"""
        if title is None:
            title = name
        #
        nLagrangeNodes_global = lagrangeNodesArray.shape[0]; nElements_global = elementNodesArray.shape[0]
        nNodes_element= elementNodesArray.shape[1]
        assert nNodes_element == nSpace+1, "affine simplicial geometry only"
        assert u_dof.shape[0] == nLagrangeNodes_global, "u_dof len=%s but nLagrangeNodes_global=%s " % (u_dof.shape[0],
                                                                                                        nLagrangeNodes_global)

        nMidNodes_element = int((nSpace+2)*(nSpace+1)/2) - nNodes_element

        nCoords = 1
        if v_dof is not None:
            nCoords += 1
            assert v_dof.shape == u_dof.shape
        if w_dof is not None:
            nCoords += 1
            assert w_dof.shape == u_dof.shape
        assert (1 <= nCoords and nCoords <= nSpace), "nCoords= %s nSpace= %s wrong " % (nCoords,nSpace)


        #store geometric nodes first, then "quadratic" ones
        if storeMeshData:
            cmdFile.write("%s_x = [ ... \n" % name)
            for nN in range(nLagrangeNodes_global):
                cmdFile.write("%g %g %g \n" % (lagrangeNodesArray[nN,0],lagrangeNodesArray[nN,1],lagrangeNodesArray[nN,2]))
            #

            cmdFile.write("];")
        #
        if nSpace == 1:
            cmdFile.write("tri_%s = [ ... \n" % name)
            for eN in range(nElements_global):
                #not necessary going to get positive jacobians
                cmdFile.write(" %i " % (l2g[eN,0]+1))
                cmdFile.write(" %i " % (l2g[eN,2]+ 1)) #use l2g?
                cmdFile.write(" \n ")
                cmdFile.write(" %i " % (l2g[eN,2] + 1))#use l2g?
                cmdFile.write(" %i " % (l2g[eN,1]+1))
                cmdFile.write(" \n ")
            cmdFile.write("];")
        elif nSpace == 2:
            cmdFile.write("tri_%s = [ ... \n" % name)
            for eN in range(nElements_global):
                #four triangles per element, see FemTools Quadratic space for numbering convention
                cmdFile.write(" %i " % (l2g[eN,0]+1))
                cmdFile.write(" %i " % (l2g[eN,nNodes_element+0]+1))
                cmdFile.write(" %i " % (l2g[eN,nNodes_element+2]+1))
                cmdFile.write(" \n ")
                cmdFile.write(" %i " % (l2g[eN,nNodes_element+0]+1))
                cmdFile.write(" %i " % (l2g[eN,1]+1))
                cmdFile.write(" %i " % (l2g[eN,nNodes_element+1]+1))
                cmdFile.write(" \n ")
                cmdFile.write(" %i " % (l2g[eN,nNodes_element+1]+1))
                cmdFile.write(" %i " % (l2g[eN,2]+1))
                cmdFile.write(" %i " % (l2g[eN,nNodes_element+2]+1))
                cmdFile.write(" \n ")
                cmdFile.write(" %i " % (l2g[eN,nNodes_element+0]+1))
                cmdFile.write(" %i " % (l2g[eN,nNodes_element+1]+1))
                cmdFile.write(" %i " % (l2g[eN,nNodes_element+2]+1))
                cmdFile.write(" \n ")
            cmdFile.write("];")
        else:
            print("""3d LagrangeC0P2 not implemented yet""")
            return 0
        #

        #assumes l2g layout consistent with matlab one
        cmdFile.write("%s = [ ... \n" % name)
        for i in range(u_dof.shape[0]):
            cmdFile.write(" %g " % u_dof[i])
            if nCoords > 1: cmdFile.write(" %g " % v_dof[i])
            if nCoords > 2: cmdFile.write(" %g " % w_dof[i])
            cmdFile.write(" \n ")

        cmdFile.write("];");

        if nSpace == 1:
            cmd = cmd1dData % (name,name,name,name,name)
            cmdFile.write(cmd)
            cmd = cmd1dView % (figureNumber,name,name,title)
            cmdFile.write(cmd)
            nplotted = 1

        elif nSpace == 2:
            #nothing to be done with mesh-data representation here
            cmd = cmd2dView % (name,name,name,
                               figureNumber,name,name,name,name,title)
            cmdFile.write(cmd)
            cmd = cmd2dGrid % (self.ngrid[0],self.ngrid[1],
                               name,name,name,name,
                               name,name,
                               name,name,name,name,name,name,
                               name,name,name,name,name,name)

            cmdFile.write(cmd)
            nplotted = 1
        else:
            #nothing to be done with mesh-data representation here
            cmd = cmd3dView % (name,name,name,name,
                               figureNumber,name,name,name,name,name,name,title)
            cmdFile.write(cmd)
            cmd = cmd3dGrid % (self.ngrid[0],self.ngrid[1],self.ngrid[2],
                               name,name,name,name,name,name,
                               name,name,name,
                               name,name,name,name,name,name,name,name,
                               name,name,name,name,name,name,name,name,
                               name,name,name,name,name,name,name,name)
            cmdFile.write(cmd)
            nplotted = 1 #

        #
        return nplotted
    #
    def viewScalar_DGP0(self,cmdFile,nSpace,nodeArray,elementNodesArray,l2g,u_dof,
                        name="u",storeMeshData=True,figureNumber=1,title=None):
        # """
        # given DG P0 function
        # generate a matlab representation that is as faithful as possible to
        # the actual finite element function structure

        # Assumes local dof associated with local node numbering
        # DG-P0 output:
        #    element-wise list of vertices and local elemntwise-connectivity matrices
        #    degrees of freedom at element vertices (as if a DG-P1 function)

        # scalar data is stored in
        #     name

        # if storeMeshData = True, writes out
        #     name_x      -- mesh vertices
        #     tri_name    -- element-node representation

        # returns number of figures actually plotted

        # """
        nplotted = 0
        ###simple visualization commands (%s --> name)
        #1d
        cmd1dData = """
%s_dg = %s_x; %s_dg(:,2) = %s;
nNodes = max(max(tri_%s_cg)); %s_cg = zeros(nNodes,1); nn_cg = zeros(nNodes,1);
nElements = size(tri_%s_cg,1); nDOF_element =  size(tri_%s_cg,2);
for eN = 1:nElements
  for i = 1:nDOF_element
     I = tri_%s_cg(eN,i);
     %s_cg(I) = %s_cg(I) + %s(tri_%s(eN,i));
     nn_cg(I) = nn_cg(I) + 1.;
  end
end
%s_cg = %s_cg ./ nn_cg;

"""
        cmd1dView = """
figure(%i); patch('vertices',%s_dg(:,1:2),'faces',tri_%s); title('%s');
"""
        #2d
        cmd2dData = """
%s_dg = %s_x; %s_dg(:,3) = %s;
nNodes = max(max(tri_%s_cg)); %s_cg = zeros(nNodes,1); nn_cg = zeros(nNodes,1);
nElements = size(tri_%s_cg,1); nDOF_element =  size(tri_%s_cg,2);
for eN = 1:nElements
  for i = 1:nDOF_element
     I = tri_%s_cg(eN,i);
     %s_cg(I) = %s_cg(I) + %s(tri_%s(eN,i));
     nn_cg(I) = nn_cg(I) + 1.;
  end
end
%s_cg = %s_cg ./ nn_cg;

"""

        cmd2dView = """
figure(%i); patch('vertices',%s_dg,'faces',tri_%s,'FaceVertexCData',%s_dg,'FaceColor','interp','EdgeColor','none'); title('%s');
"""

        cmd2dGrid = """
nx = %i; ny = %i;
XYZ = [min(%s_x(:,1)) max(%s_x(:,1)) ; min(%s_x(:,2)) max(%s_x(:,2))];
dx = (XYZ(1,2)-XYZ(1,1))/nx; dy = (XYZ(2,2)-XYZ(2,1))/ny;
[%s_xg,%s_yg] = meshgrid(XYZ(1,1):dx:XYZ(1,2),XYZ(2,1):dy:XYZ(2,2));
%%note, uses average of duplicate values
%s_g = griddata(%s_x(:,1),%s_x(:,2),%s,%s_xg,%s_yg);
"""
        #3d
        cmd3dData = """
nNodes = max(max(tri_%s_cg)); %s_cg = zeros(nNodes,1); nn_cg = zeros(nNodes,1);
nElements = size(tri_%s_cg,1); nDOF_element =  size(tri_%s_cg,2);
for eN = 1:nElements
  for i = 1:nDOF_element
     I = tri_%s_cg(eN,i);
     %s_cg(I) = %s_cg(I) + %s(tri_%s(eN,i));
     nn_cg(I) = nn_cg(I) + 1.;
  end
end
%s_cg = %s_cg ./ nn_cg;
"""

        cmd3dView = """
figure(%i); patch('vertices',%s_x,'faces',tri_%s,'FaceVertexCData',%s,'FaceColor','interp','EdgeColor','none'); title('%s');
"""
        cmd3dGrid = """
nx = %i; ny = %i; nz = %i ;
XYZ = [min(%s_x(:,1)) max(%s_x(:,1)) ; min(%s_x(:,2)) max(%s_x(:,2)) ; min(%s_x(:,3)) max(%s_x(:,3))];
dx = (XYZ(1,2)-XYZ(1,1))/nx; dy = (XYZ(2,2)-XYZ(2,1))/ny; dz = (XYZ(3,2)-XYZ(3,1))/nz;
[%s_xg,%s_yg,%s_zg] = meshgrid(XYZ(1,1):dx:XYZ(1,2),XYZ(2,1):dy:XYZ(2,2),XYZ(3,1):dz:XYZ(3,2));
%%note, uses average of duplicate values
%s_g = griddata3(%s_x(:,1),%s_x(:,2),%s_x(:,3),%s,%s_xg,%s_yg,%s_zg);
"""
        if title is None:
            title = name

        nNodes_global = nodeArray.shape[0]; nElements_global = l2g.shape[0]
        nNodes_element= elementNodesArray.shape[1]
        nDOF_element = l2g.shape[1]
        if storeMeshData:
            cmdFile.write("%s_x = [ ... \n" % name)
            for eN in range(nElements_global):
                for nN_local in range(nNodes_element):
                    nN = elementNodesArray[eN,nN_local]
                    cmdFile.write("%g %g %g \n" % (nodeArray[nN,0],nodeArray[nN,1],nodeArray[nN,2]))
            #
            cmdFile.write("];")
            cmdFile.write("%s_x_cg = [ ... \n" % name)
            for nN in range(nNodes_global):
                cmdFile.write("%g %g %g \n" % (nodeArray[nN,0],nodeArray[nN,1],nodeArray[nN,2]))
            #
            cmdFile.write("];")
        #
        cmdFile.write("tri_%s = [ ... \n" % name)
        for eN in range(nElements_global):
            for nN in range(nNodes_element):#assumes laid out line element nodes
                cmdFile.write(" %i " % (eN*nNodes_element+nN+1))
            cmdFile.write("; \n ")
        #
        cmdFile.write("];")
        #
        #provide CG connectivity if wanted nodal averages for some reason
        cmdFile.write("tri_%s_cg = [ ... \n" % name)
        for eN in range(nElements_global):
            for i in range(nDOF_element):#assumes laid out line element nodes
                cmdFile.write(" %i " % (elementNodesArray[eN,i]+1))
            cmdFile.write("; \n ")
        #
        cmdFile.write("];")
        #

        cmdFile.write("%s = [ ... \n" % name)
        for eN in range(nElements_global):
            for nN in range(nNodes_element):
                cmdFile.write("%g \n" % u_dof[l2g[eN,0]])
        cmdFile.write("];");

        if nSpace == 1:
            cmd = cmd1dData % (name,name,name,name,
                               name,name,
                               name,name,
                               name,
                               name,name,name,name,
                               name,name)
            cmdFile.write(cmd)
            cmd = cmd1dView % (figureNumber,name,name,title)
            cmdFile.write(cmd)
            nplotted = 1

        elif nSpace == 2:
            cmd = cmd2dData % (name,name,name,name,
                               name,name,
                               name,name,
                               name,
                               name,name,name,name,
                               name,name)
            cmdFile.write(cmd)
            cmd = cmd2dView % (figureNumber,name,name,name,title)
            cmdFile.write(cmd)
            cmd = cmd2dGrid % (self.ngrid[0],self.ngrid[1],
                               name,name,name,name,
                               name,name,
                               name,name,name,name,name,name)

            cmdFile.write(cmd)
            nplotted = 1
        else:
            cmd = cmd3dData % (name,name,
                               name,name,
                               name,
                               name,name,name,name,
                               name,name)
            cmdFile.write(cmd)
            cmd = cmd3dView % (figureNumber,name,name,name,title)
            cmdFile.write(cmd)
            cmd = cmd3dGrid % (self.ngrid[0],self.ngrid[1],self.ngrid[2],
                               name,name,name,name,name,name,
                               name,name,name,
                               name,name,name,name,name,name,name,name)
            cmdFile.write(cmd)
            nplotted = 1

        #
        return nplotted
    #

    def viewScalar_LagrangeDGP1(self,cmdFile,nSpace,nodeArray,elementNodesArray,l2g,u_dof,
                                name="u",storeMeshData=True,figureNumber=1,title=None):
        # """
        # given DG P1 function with nodal Lagrange representation
        # generate a matlab representation that is as faithful as possible to
        # the actual finite element function structure

        # Assumes local dof associated with local node numbering
        # DG-P1 output:
        #    element-wise list of vertices and local elemntwise-connectivity matrices
        #    degrees of freedom at vertices

        # scalar data is stored in
        #     name

        # if storeMeshData = True, writes out
        #     name_x      -- mesh vertices
        #     tri_name    -- element-node representation

        # returns number of figures actually plotted

        # """
        nplotted = 0
        ###simple visualization commands (%s --> name)
        #1d
        cmd1dData = """
%s_dg = %s_x; %s_dg(:,2) = %s;
nNodes = max(max(tri_%s_cg)); %s_cg = zeros(nNodes,1); nn_cg = zeros(nNodes,1);
nElements = size(tri_%s_cg,1); nDOF_element =  size(tri_%s_cg,2);
for eN = 1:nElements
  for i = 1:nDOF_element
     I = tri_%s_cg(eN,i);
     %s_cg(I) = %s_cg(I) + %s(tri_%s(eN,i));
     nn_cg(I) = nn_cg(I) + 1.;
  end
end
%s_cg = %s_cg ./ nn_cg;

"""
        cmd1dView = """
figure(%i); patch('vertices',%s_dg(:,1:2),'faces',tri_%s); title('%s');
"""
        #2d
        cmd2dData = """
%s_dg = %s_x; %s_dg(:,3) = %s;
nNodes = max(max(tri_%s_cg)); %s_cg = zeros(nNodes,1); nn_cg = zeros(nNodes,1);
nElements = size(tri_%s_cg,1); nDOF_element =  size(tri_%s_cg,2);
for eN = 1:nElements
  for i = 1:nDOF_element
     I = tri_%s_cg(eN,i);
     %s_cg(I) = %s_cg(I) + %s(tri_%s(eN,i));
     nn_cg(I) = nn_cg(I) + 1.;
  end
end
%s_cg = %s_cg ./ nn_cg;

"""

        cmd2dView = """
figure(%i); patch('vertices',%s_dg,'faces',tri_%s,'FaceVertexCData',%s_dg,'FaceColor','interp','EdgeColor','none'); title('%s');
"""

        cmd2dGrid = """
nx = %i; ny = %i;
XYZ = [min(%s_x(:,1)) max(%s_x(:,1)) ; min(%s_x(:,2)) max(%s_x(:,2))];
dx = (XYZ(1,2)-XYZ(1,1))/nx; dy = (XYZ(2,2)-XYZ(2,1))/ny;
[%s_xg,%s_yg] = meshgrid(XYZ(1,1):dx:XYZ(1,2),XYZ(2,1):dy:XYZ(2,2));
%%note, uses average of duplicate values
%s_g = griddata(%s_x(:,1),%s_x(:,2),%s,%s_xg,%s_yg);
"""
        #3d
        cmd3dData = """
nNodes = max(max(tri_%s_cg)); %s_cg = zeros(nNodes,1); nn_cg = zeros(nNodes,1);
nElements = size(tri_%s_cg,1); nDOF_element =  size(tri_%s_cg,2);
for eN = 1:nElements
  for i = 1:nDOF_element
     I = tri_%s_cg(eN,i);
     %s_cg(I) = %s_cg(I) + %s(tri_%s(eN,i));
     nn_cg(I) = nn_cg(I) + 1.;
  end
end
%s_cg = %s_cg ./ nn_cg;
"""

        cmd3dView = """
figure(%i); patch('vertices',%s_x,'faces',tri_%s,'FaceVertexCData',%s,'FaceColor','interp','EdgeColor','none'); title('%s');
"""
        cmd3dGrid = """
nx = %i; ny = %i; nz = %i ;
XYZ = [min(%s_x(:,1)) max(%s_x(:,1)) ; min(%s_x(:,2)) max(%s_x(:,2)) ; min(%s_x(:,3)) max(%s_x(:,3))];
dx = (XYZ(1,2)-XYZ(1,1))/nx; dy = (XYZ(2,2)-XYZ(2,1))/ny; dz = (XYZ(3,2)-XYZ(3,1))/nz;
[%s_xg,%s_yg,%s_zg] = meshgrid(XYZ(1,1):dx:XYZ(1,2),XYZ(2,1):dy:XYZ(2,2),XYZ(3,1):dz:XYZ(3,2));
%%note, uses average of duplicate values
%s_g = griddata3(%s_x(:,1),%s_x(:,2),%s_x(:,3),%s,%s_xg,%s_yg,%s_zg);
"""
        if title is None:
            title = name

        nNodes_global = nodeArray.shape[0]; nElements_global = l2g.shape[0]
        nNodes_element= elementNodesArray.shape[1]
        nDOF_element = l2g.shape[1]
        if storeMeshData:
            cmdFile.write("%s_x = [ ... \n" % name)
            for eN in range(nElements_global):
                for nN_local in range(nNodes_element):
                    nN = elementNodesArray[eN,nN_local]
                    cmdFile.write("%g %g %g \n" % (nodeArray[nN,0],nodeArray[nN,1],nodeArray[nN,2]))
            #
            cmdFile.write("];")
            cmdFile.write("%s_x_cg = [ ... \n" % name)
            for nN in range(nNodes_global):
                cmdFile.write("%g %g %g \n" % (nodeArray[nN,0],nodeArray[nN,1],nodeArray[nN,2]))
            #
            cmdFile.write("];")
        #
        cmdFile.write("tri_%s = [ ... \n" % name)
        for eN in range(nElements_global):
            for i in range(nDOF_element):#assumes laid out line element nodes
                cmdFile.write(" %i " % (l2g[eN,i]+1))
            cmdFile.write("; \n ")
        #
        cmdFile.write("];")
        #
        #provide CG connectivity
        cmdFile.write("tri_%s_cg = [ ... \n" % name)
        for eN in range(nElements_global):
            for i in range(nDOF_element):#assumes laid out line element nodes
                cmdFile.write(" %i " % (elementNodesArray[eN,i]+1))
            cmdFile.write("; \n ")
        #
        cmdFile.write("];")
        #

        cmdFile.write("%s = [ ... \n" % name)
        for i in range(u_dof.shape[0]):
            cmdFile.write("%g \n" % u_dof[i])
        cmdFile.write("];");

        if nSpace == 1:
            cmd = cmd1dData % (name,name,name,name,
                               name,name,
                               name,name,
                               name,
                               name,name,name,name,
                               name,name)
            cmdFile.write(cmd)
            cmd = cmd1dView % (figureNumber,name,name,title)
            cmdFile.write(cmd)
            nplotted = 1

        elif nSpace == 2:
            cmd = cmd2dData % (name,name,name,name,
                               name,name,
                               name,name,
                               name,
                               name,name,name,name,
                               name,name)
            cmdFile.write(cmd)
            cmd = cmd2dView % (figureNumber,name,name,name,title)
            cmdFile.write(cmd)
            cmd = cmd2dGrid % (self.ngrid[0],self.ngrid[1],
                               name,name,name,name,
                               name,name,
                               name,name,name,name,name,name)

            cmdFile.write(cmd)
            nplotted = 1
        else:
            cmd = cmd3dData % (name,name,
                               name,name,
                               name,
                               name,name,name,name,
                               name,name)
            cmdFile.write(cmd)
            cmd = cmd3dView % (figureNumber,name,name,name,title)
            cmdFile.write(cmd)
            cmd = cmd3dGrid % (self.ngrid[0],self.ngrid[1],self.ngrid[2],
                               name,name,name,name,name,name,
                               name,name,name,
                               name,name,name,name,name,name,name,name)
            cmdFile.write(cmd)
            nplotted = 1

        #
        return nplotted
    #
    def viewScalar_LagrangeDGP2(self,cmdFile,nSpace,nodeArray,elementNodesArray,midNodesArray,
                                l2g,cg_l2g,u_dof,
                                name="u",storeMeshData=True,figureNumber=1,title=None):
        # """
        # TODO: 3d
        # given DG P1 function with nodal Lagrange representation
        # generate a matlab representation that is as faithful as possible to
        # the actual finite element function structure

        # Assumes local dof associated with local node numbering
        # and then local edge numbering

        # uses cg numbering for accessing mid-edge vertices
        # DG-P2 output:
        #    element-wise list of vertices and mid-edge vertices along
        #    with local elemntwise-connectivity matrices
        #    degrees of freedom at vertices and mid-edge vertices

        # scalar data is stored in
        #     name

        # if storeMeshData = True, writes out
        #     name_x      -- mesh vertices and mid-edge vertices
        #     tri_name    -- element-node representation

        # returns number of figures actually plotted

        # """
        nplotted = 0
        ###simple visualization commands (%s --> name)
        #1d
        cmd1dData = """
%s_dg = %s_x; %s_dg(:,2) = %s;

"""
        cmd1dView = """
figure(%i); patch('vertices',%s_dg(:,1:2),'faces',tri_%s); title('%s');
"""
        #2d
        cmd2dData = """
%s_dg = %s_x; %s_dg(:,3) = %s;

"""

        cmd2dView = """
figure(%i); patch('vertices',%s_dg,'faces',tri_%s,'FaceVertexCData',%s_dg,'FaceColor','interp','EdgeColor','none'); title('%s');
"""

        cmd2dGrid = """
nx = %i; ny = %i;
XYZ = [min(%s_x(:,1)) max(%s_x(:,1)) ; min(%s_x(:,2)) max(%s_x(:,2))];
dx = (XYZ(1,2)-XYZ(1,1))/nx; dy = (XYZ(2,2)-XYZ(2,1))/ny;
[%s_xg,%s_yg] = meshgrid(XYZ(1,1):dx:XYZ(1,2),XYZ(2,1):dy:XYZ(2,2));
%%note, uses average of duplicate values
%s_g = griddata(%s_x(:,1),%s_x(:,2),%s,%s_xg,%s_yg);
"""
        #3d
        cmd3dData = """
"""

        cmd3dView = """
figure(%i); patch('vertices',%s_dg,'faces',tri_%s,'FaceVertexCData',%s,'FaceColor','interp','EdgeColor','none'); title('%s');
"""
        cmd3dGrid = """
nx = %i; ny = %i; nz = %i ;
XYZ = [min(%s_x(:,1)) max(%s_x(:,1)) ; min(%s_x(:,2)) max(%s_x(:,2)) ; min(%s_x(:,3)) max(%s_x(:,3))];
dx = (XYZ(1,2)-XYZ(1,1))/nx; dy = (XYZ(2,2)-XYZ(2,1))/ny; dz = (XYZ(3,2)-XYZ(3,1))/nz;
[%s_xg,%s_yg,%s_zg] = meshgrid(XYZ(1,1):dx:XYZ(1,2),XYZ(2,1):dy:XYZ(2,2),XYZ(3,1):dz:XYZ(3,2));
%%note, uses average of duplicate values
%s_g = griddata3(%s_x(:,1),%s_x(:,2),%s_x(:,3),%s,%s_xg,%s_yg,%s_zg);
"""
        if title is None:
            title = name

        nNodes_global = nodeArray.shape[0]; nElements_global = l2g.shape[0]
        nNodes_element= elementNodesArray.shape[1]
        assert nNodes_element == nSpace+1, "affine simplicial geometry only"

        nMidNodes_global     = midNodesArray.shape[0]
        nMidNodes_element = int((nSpace+2)*(nSpace+1)/2) - nNodes_element
        nDOF_element = l2g.shape[1]
        assert nDOF_element == nMidNodes_element + nNodes_element

        #assumes local ordering vertices then midedge nodes
        if storeMeshData:
            cmdFile.write("%s_x = [ ... \n" % name)
            for eN in range(nElements_global):
                for nN_local in range(nNodes_element):
                    nN = elementNodesArray[eN,nN_local]
                    cmdFile.write("%g %g %g \n" % (nodeArray[nN,0],nodeArray[nN,1],nodeArray[nN,2]))
                for nM_local in range(nMidNodes_element):
                    nM = cg_l2g[eN,nNodes_element+nM_local] - nNodes_global
                    cmdFile.write("%g %g %g \n" % (midNodesArray[nM,0],midNodesArray[nM,1],midNodesArray[nM,2]))
            #
            cmdFile.write("];")
        #
        cmdFile.write("tri_%s = [ ... \n" % name)
        if nSpace == 1:
            for eN in range(nElements_global):
                cmdFile.write(" %i " % (l2g[eN,0]+1))
                cmdFile.write(" %i " % (l2g[eN,nNodes_element+0]+1))
                cmdFile.write(" \n ")
                cmdFile.write(" %i " % (l2g[eN,nNodes_element+0]+1))
                cmdFile.write(" %i \n" % (l2g[eN,1]+1))
            cmdFile.write("; \n ")
        elif nSpace == 2:
            for eN in range(nElements_global):
                #four triangles per element, see FemTools Quadratic space for numbering convention
                cmdFile.write(" %i " % (l2g[eN,0]+1))
                cmdFile.write(" %i " % (l2g[eN,nNodes_element+0]+1))
                cmdFile.write(" %i " % (l2g[eN,nNodes_element+2]+1))
                cmdFile.write(" \n ")
                cmdFile.write(" %i " % (l2g[eN,nNodes_element+0]+1))
                cmdFile.write(" %i " % (l2g[eN,1]+1))
                cmdFile.write(" %i " % (l2g[eN,nNodes_element+1]+1))
                cmdFile.write(" \n ")
                cmdFile.write(" %i " % (l2g[eN,nNodes_element+1]+1))
                cmdFile.write(" %i " % (l2g[eN,2]+1))
                cmdFile.write(" %i " % (l2g[eN,nNodes_element+2]+1))
                cmdFile.write(" \n ")
                cmdFile.write(" %i " % (l2g[eN,nNodes_element+0]+1))
                cmdFile.write(" %i " % (l2g[eN,nNodes_element+1]+1))
                cmdFile.write(" %i " % (l2g[eN,nNodes_element+2]+1))
                cmdFile.write(" \n ")
            cmdFile.write("; \n ")
        else:
            print("""3d LagrangeDGP2 not implemented yet""")
            return 0

        #
        cmdFile.write("];")
        #

        cmdFile.write("%s = [ ... \n" % name)
        for i in range(u_dof.shape[0]):
            cmdFile.write("%g \n" % u_dof[i])
        cmdFile.write("];");

        if nSpace == 1:
            cmd = cmd1dData % (name,name,name,name)
            cmdFile.write(cmd)
            cmd = cmd1dView % (figureNumber,name,name,title)
            cmdFile.write(cmd)
            nplotted = 1

        elif nSpace == 2:
            cmd = cmd2dData % (name,name,name,name)
            cmdFile.write(cmd)
            cmd = cmd2dView % (figureNumber,name,name,name,title)
            cmdFile.write(cmd)
            cmd = cmd2dGrid % (self.ngrid[0],self.ngrid[1],
                               name,name,name,name,
                               name,name,
                               name,name,name,name,name,name)

            cmdFile.write(cmd)
            nplotted = 1
        else:
            cmd = cmd3dView % (figureNumber,name,name,name,title)
            cmdFile.write(cmd)
            cmd = cmd3dGrid % (self.ngrid[0],self.ngrid[1],self.ngrid[2],
                               name,name,name,name,name,name,
                               name,name,name,
                               name,name,name,name,name,name,name,name)
            cmdFile.write(cmd)
            nplotted = 1

        #
        return nplotted
    #
    def viewScalar_CrouzeixRaviartP1(self,cmdFile,nSpace,nodeArray,elementNodesArray,l2g,u_dof,
                                     name="u",storeMeshData=True,figureNumber=1,title=None):
        # """
        # given FEM function with local Crouzeix-Raviart representation
        # generate a matlab representation that is as faithful as possible to
        # the actual finite element function structure

        # Assumes local dof associated with local node numbering
        # CR output:
        #    Just treat as a DGP1 function
        #    element-wise list of vertices and local elemntwise-connectivity matrices
        #    degrees of freedom at vertices instead of face barycenters

        # scalar data is stored in
        #     name

        # if storeMeshData = True, writes out
        #     name_x      -- mesh vertices
        #     tri_name    -- element-node representation

        # returns number of figures actually plotted

        # """
        nplotted = 0
        ###simple visualization commands (%s --> name)
        #1d
        cmd1dData = """
%%actully just C0P1 in 1d
%s_cr = %s_x; %s_cr(:,2) = %s;
"""
        cmd1dView = """
figure(%i); patch('vertices',%s_cr(:,1:2),'faces',tri_%s); title('%s');
"""
        #2d
        cmd2dData = """
%s_cr = %s_x; %s_cr(:,3) = %s;
nNodes = max(max(tri_%s_cg)); %s_cg = zeros(nNodes,1); nn_cg = zeros(nNodes,1);
nElements = size(tri_%s_cg,1); nDOF_element =  size(tri_%s_cg,2);
for eN = 1:nElements
  for i = 1:nDOF_element
     I = tri_%s_cg(eN,i);
     %s_cg(I) = %s_cg(I) + %s(tri_%s(eN,i));
     nn_cg(I) = nn_cg(I) + 1.;
  end
end
%s_cg = %s_cg ./ nn_cg;
"""

        cmd2dView = """
figure(%i); patch('vertices',%s_cr,'faces',tri_%s,'FaceVertexCData',%s_cr,'FaceColor','interp','EdgeColor','none'); title('%s');
%%for average cg interpolant can also use
%%tmp = name_x_cg; tmp(:,3) = name_cg;
%%patch('vertices',tmp,'faces',tri_name_cg,'FaceVertexCData',name_cg,'FaceColor','interp','EdgeColor','none');
"""

        cmd2dGrid = """
nx = %i; ny = %i;
XYZ = [min(%s_x(:,1)) max(%s_x(:,1)) ; min(%s_x(:,2)) max(%s_x(:,2))];
dx = (XYZ(1,2)-XYZ(1,1))/nx; dy = (XYZ(2,2)-XYZ(2,1))/ny;
[%s_xg,%s_yg] = meshgrid(XYZ(1,1):dx:XYZ(1,2),XYZ(2,1):dy:XYZ(2,2));
%%note, uses average of duplicate values
%s_g = griddata(%s_x(:,1),%s_x(:,2),%s,%s_xg,%s_yg);
"""
        #3d
        cmd3dData = """
nNodes = max(max(tri_%s_cg)); %s_cg = zeros(nNodes,1); nn_cg = zeros(nNodes,1);
nElements = size(tri_%s_cg,1); nDOF_element =  size(tri_%s_cg,2);
for eN = 1:nElements
  for i = 1:nDOF_element
     I = tri_%s_cg(eN,i);
     %s_cg(I) = %s_cg(I) + %s(tri_%s(eN,i));
     nn_cg(I) = nn_cg(I) + 1.;
  end
end
%s_cg = %s_cg ./ nn_cg;
"""

        cmd3dView = """
figure(%i); patch('vertices',%s_cr,'faces',tri_%s,'FaceVertexCData',%s,'FaceColor','interp','EdgeColor','none'); title('%s');
"""
        cmd3dGrid = """
nx = %i; ny = %i; nz = %i ;
XYZ = [min(%s_x(:,1)) max(%s_x(:,1)) ; min(%s_x(:,2)) max(%s_x(:,2)) ; min(%s_x(:,3)) max(%s_x(:,3))];
dx = (XYZ(1,2)-XYZ(1,1))/nx; dy = (XYZ(2,2)-XYZ(2,1))/ny; dz = (XYZ(3,2)-XYZ(3,1))/nz;
[%s_xg,%s_yg,%s_zg] = meshgrid(XYZ(1,1):dx:XYZ(1,2),XYZ(2,1):dy:XYZ(2,2),XYZ(3,1):dz:XYZ(3,2));
%%note, uses average of duplicate values
%s_g = griddata3(%s_x(:,1),%s_x(:,2),%s_x(:,3),%s,%s_xg,%s_yg,%s_zg);
"""
        import numpy
        if title is None:
            title = name

        nNodes_global = nodeArray.shape[0]; nElements_global = l2g.shape[0]
        nNodes_element= elementNodesArray.shape[1]
        nDOF_element = l2g.shape[1]
        assert nDOF_element == nNodes_element
        #treat like a dgp1 function with Lagrange basis
        if storeMeshData:
            cmdFile.write("%s_x = [ ... \n" % name)
            for eN in range(nElements_global):
                for nN_local in range(nNodes_element):
                    nN = elementNodesArray[eN,nN_local]
                    cmdFile.write("%g %g %g \n" % (nodeArray[nN,0],nodeArray[nN,1],nodeArray[nN,2]))
            #
            cmdFile.write("];")
            cmdFile.write("%s_x_cg = [ ... \n" % name)
            for nN in range(nNodes_global):
                cmdFile.write("%g %g %g \n" % (nodeArray[nN,0],nodeArray[nN,1],nodeArray[nN,2]))
            #
            cmdFile.write("];")
        #
        cmdFile.write("tri_%s = [ ... \n" % name)
        for eN in range(nElements_global):
            for i in range(nDOF_element):#assumes laid out line element nodes
                cmdFile.write(" %i " % (eN*nNodes_element+i+1))
            cmdFile.write("; \n ")
        #
        cmdFile.write("];")
        #
        #provide CG connectivity
        cmdFile.write("tri_%s_cg = [ ... \n" % name)
        for eN in range(nElements_global):
            for i in range(nDOF_element):#assumes laid out line element nodes
                cmdFile.write(" %i " % (elementNodesArray[eN,i]+1))
            cmdFile.write("; \n ")
        #
        cmdFile.write("];")
        #

        #dofs
        cmdFile.write("%s = [ ... \n" % name)
        #u(vertex_i) = u^i(1-nSpace) + \sum_{j \ne i} u^j where u^j is local dof j

        if nSpace == 1:
            for eN in range(nElements_global):
                #dof associated with face id, so opposite usual C0P1 ordering here
                cmdFile.write(" %g \n" % u_dof[l2g[eN,1]])
                cmdFile.write(" %g \n" % u_dof[l2g[eN,0]])
            #
        elif nSpace == 2:
            nodeVal = numpy.zeros(3,'d')
            for eN in range(nElements_global):
                nodeVal *= 0.0
                #assume vertex associated with face across from it
                nodeVal[0] = u_dof[l2g[eN,1]]
                nodeVal[0]+= u_dof[l2g[eN,2]]
                nodeVal[0]-= u_dof[l2g[eN,0]]

                nodeVal[1] = u_dof[l2g[eN,0]]
                nodeVal[1]+= u_dof[l2g[eN,2]]
                nodeVal[1]-= u_dof[l2g[eN,1]]

                nodeVal[2] = u_dof[l2g[eN,0]]
                nodeVal[2]+= u_dof[l2g[eN,1]]
                nodeVal[2]-= u_dof[l2g[eN,2]]

                cmdFile.write(" %g \n %g \n %g \n " % (nodeVal[0],nodeVal[1],nodeVal[2]))
            #
        else:
            for eN in range(nElements_global):
                for i in range(nNodes_element):
                    nodalVal = u_dof[l2g[eN,i]]*(1.0-float(nSpace)) + sum([u_dof[l2g[eN,j]] for j in range(nDOF_element) if j != i])
                    cmdFile.write(" %g \n " % nodeVal)
                #
            #
        #
        cmdFile.write("];");

        if nSpace == 1:
            cmd = cmd1dData % (name,name,name,name)
            cmdFile.write(cmd)
            cmd = cmd1dView % (figureNumber,name,name,title)
            cmdFile.write(cmd)
            nplotted = 1

        elif nSpace == 2:
            cmd = cmd2dData % (name,name,name,name,
                               name,name,
                               name,name,
                               name,
                               name,name,name,name,
                               name,name)
            cmdFile.write(cmd)
            cmd = cmd2dView % (figureNumber,name,name,name,title)
            cmdFile.write(cmd)
            cmd = cmd2dGrid % (self.ngrid[0],self.ngrid[1],
                               name,name,name,name,
                               name,name,
                               name,name,name,name,name,name)

            cmdFile.write(cmd)
            nplotted = 1
        else:
            cmd = cmd3dData % (name,name,
                               name,name,
                               name,
                               name,name,name,name,
                               name,name)
            cmdFile.write(cmd)
            cmd = cmd3dView % (figureNumber,name,name,name,title)
            cmdFile.write(cmd)
            cmd = cmd3dGrid % (self.ngrid[0],self.ngrid[1],self.ngrid[2],
                               name,name,name,name,name,name,
                               name,name,name,
                               name,name,name,name,name,name,name,name)
            cmdFile.write(cmd)
            nplotted = 1

        #
        return nplotted
    #
    def viewScalar_MonomialDGPK(self,cmdFile,nSpace,nodeArray,elementNodesArray,
                                interpolationPoints,u_interpolationPoints,
                                name="u",storeMeshData=True,figureNumber=1,title=None):
        # """
        # given DG PK function with monomial basis
        # generate a matlab representation that is as faithful as possible to
        # the actual finite element function structure

        # input is u values at interpolation points (not dof)
        # Generates an array of triangulations of interpolation points on each element
        # DG-PK output:
        #    element-wise list of interpolation points and local elementwise-connectivity matrices
        #    degrees of freedom at interpolation points

        # scalar data is stored in
        #     name

        # if storeMeshData = True, writes out
        #     name_x      -- mesh vertices
        #     tri_name    -- element-node representation

        # returns number of figures actually plotted

        # """
        #1d
        cmd1dData = """
nElements_global = %i; nPoints_element = %i;
tri_%s = [];
for eN = 1:nElements_global
  [x_tmp,i_tmp] = sort(%s_x(nPoints_element*(eN-1)+1:nPoints_element*eN,1));
  tmp = [];
  for j = 1:nPoints_element-1
    tmp = [tmp ; i_tmp(j) i_tmp(j+1)];
  end
  tri_%s = [tri_%s ; nPoints_element*(eN-1) + tmp];
end
%s_loc = %s_x; %s_loc(:,2) = %s;
"""
        cmd1dView = """
figure(%i) ; patch('vertices',%s_loc,'faces',tri_%s,'FaceColor','none','EdgeColor','black'); title('%s');
"""
        #2d
        cmd2dData = """
nElements_global = %i; nPoints_element = %i;
tri_%s = [];
for eN = 1:nElements_global
  tmp = delaunay(%s_x(nPoints_element*(eN-1)+1:nPoints_element*eN,1),...
                 %s_x(nPoints_element*(eN-1)+1:nPoints_element*eN,2));
  tri_%s = [tri_%s ; nPoints_element*(eN-1)+tmp];
end
%s_loc = %s_x; %s_loc(:,3) = %s;

"""
        cmd2dView = """
figure(%i) ; patch('vertices',%s_loc,'faces',tri_%s,'FaceVertexCData',%s,'FaceColor','interp','EdgeColor','none'); title('%s');
"""
        cmd2dGrid = """
nx = %i; ny = %i;
XYZ = [min(%s_x(:,1)) max(%s_x(:,1)) ; min(%s_x(:,2)) max(%s_x(:,2))];
dx = (XYZ(1,2)-XYZ(1,1))/nx; dy = (XYZ(2,2)-XYZ(2,1))/ny;
[%s_xg,%s_yg] = meshgrid(XYZ(1,1):dx:XYZ(1,2),XYZ(2,1):dy:XYZ(2,2));
%s_g = griddata(%s_x(:,1),%s_x(:,2),%s,%s_xg,%s_yg);
"""
        #3d
        cmd3dData = """
nElements_global = %i; nPoints_element = %i;
tri_%s = [];
for eN = 1:nElements_global
  tmp = delaunay3(%s_x(nPoints_element*(eN-1)+1:nPoints_element*eN,1),...
                  %s_x(nPoints_element*(eN-1)+1:nPoints_element*eN,2),...
                  %s_x(nPoints_element*(eN-1)+1:nPoints_element*eN,3));
  tri_%s = [tri_%s ; nPoints_element*(eN-1)+tmp];
end
%s_loc = %s_x; %s_loc(:,3) = %s;

"""
        cmd3dView = """
%%good luck
figure(%i); patch('vertices',%s_x,'faces',tri_%s,'FaceVertexCData',%s,'FaceColor','none','EdgeColor','interp'); title('%s')
"""
        cmd3dGrid = """
nx = %i; ny = %i; nz = %i ;
XYZ = [min(%s_x(:,1)) max(%s_x(:,1)) ; min(%s_x(:,2)) max(%s_x(:,2)) ; min(%s_x(:,3)) max(%s_x(:,3))];
dx = (XYZ(1,2)-XYZ(1,1))/nx; dy = (XYZ(2,2)-XYZ(2,1))/ny; dz = (XYZ(3,2)-XYZ(3,1))/nz;
[%s_xg,%s_yg,%s_zg] = meshgrid(XYZ(1,1):dx:XYZ(1,2),XYZ(2,1):dy:XYZ(2,2),XYZ(3,1):dz:XYZ(3,2));
%s_g = griddata3(%s_x(:,1),%s_x(:,2),%s_x(:,3),%s,%s_xg,%s_yg,%s_zg);
"""
        nplotted = 0
        if title is None:
            title = name

        nElements_global = interpolationPoints.shape[0]; nPoints_element = interpolationPoints.shape[1];
        #if constants use DGP0
        if nPoints_element == 1:
            return self.viewScalar_DGP0(cmdFile,nSpace,nodeArray,elementNodesArray,l2g,u_dof,
                                        name=name,storeMeshData=storeMeshData,
                                        figureNumber=figureNumber,title=title)
        assert nPoints_element > nSpace
        if storeMeshData:
            cmdFile.write("%s_x = [ ... \n" % name)
            for eN in range(nElements_global):
                for k in range(nPoints_element):
                    cmdFile.write("%g %g %g \n" % (interpolationPoints[eN,k,0],
                                                   interpolationPoints[eN,k,1],
                                                   interpolationPoints[eN,k,2]))
            cmdFile.write("];")

        cmdFile.write("%s = [ ... \n" % name)
        for eN in range(nElements_global):
            for k in range(nPoints_element):
                cmdFile.write("%g \n" % u_interpolationPoints[eN,k])
        cmdFile.write("];");

        if nSpace == 1:
            cmd = cmd1dData % (nElements_global,nPoints_element,
                               name,
                               name,
                               name,name,
                               name,name,name,name)

            cmdFile.write(cmd)
            cmd = cmd1dView % (figureNumber,name,name,title)
            cmdFile.write(cmd)
            nplotted = 1

        elif nSpace == 2:
            cmd = cmd2dData % (nElements_global,nPoints_element,
                               name,
                               name,name,
                               name,name,
                               name,name,name,name)
            cmdFile.write(cmd)
            cmd = cmd2dView % (figureNumber,name,name,name,title)
            cmdFile.write(cmd)
            cmd = cmd2dGrid % (self.ngrid[0],self.ngrid[1],
                               name,name,name,name,
                               name,name,
                               name,name,name,name,name,name)

            cmdFile.write(cmd)
            nplotted = 1

        else:
            cmd = cmd3dData % (nElements_global,nPoints_element,
                               name,
                               name,name,name,
                               name,name,
                               name,name,name,name)
            cmdFile.write(cmd)
            cmd = cmd3dView % (figureNumber,name,name,name,title)
            cmdFile.write(cmd)
            cmd = cmd3dGrid % (self.ngrid[0],self.ngrid[1],self.ngrid[2],
                               name,name,name,name,name,name,
                               name,name,name,
                               name,name,name,name,name,name,name,name)
            cmdFile.write(cmd)
            nplotted = 1
        return nplotted
    #

import vtk
from vtk import *
from proteus.MeshTools import *
from proteus.FemTools import *
from proteus import Quadrature
import sys

hasQt=True
useCoPro=False

if hasQt:
    from vtk.qt4.QVTKRenderWindowInteractor import *
    from PyQt4 import QtGui,Qt,QtCore


if useCoPro:
    def initializeCoProcessor(comm,controller):
        import paraview
        from paraview import servermanager
        if comm.rank() == 0:
            print "initializing coprocessor"
        paraview.options.batch = True
        paraview.options.symmetric = True
        pm = servermanager.vtkProcessModule.GetProcessModule()
        globalController = pm.GetGlobalController()
        globalController.SetGlobalController(controller)
#        if globalController == None:
#            globalController = vtk.vtkMPIController()
#            globalController.Initialize()
#            globalController.SetGlobalController(globalController)
#            newGlobalController = globalController.PartitionController(0, comm.size()-comm.rank()-1)
#            newGlobalController.SetGlobalController(newGlobalController)
#        elif globalController.IsA("vtkDummyController") == True:
#            globalController = vtk.vtkMPIController()
#            globalController.Initialize()
#            globalController.SetGlobalController(globalController)
#            newGlobalController = globalController.PartitionController(0, comm.rank())
#            newGlobalController.SetGlobalController(newGlobalController)
#        else:
#            print globalController.GetLocalProcessId(), ' we have a global controller from the process module. the mpcontroller is ', globalController.GetGlobalController(), pm.GetReferenceCount()
#end coprocessing init
from proteus import flcbdfWrappers
#cek todo
#put back in hard copy capability
#improve 3D views (cut and slice planes with plane widget)
#do viewMesh and viewDomain
#clean up interfaces, connect with, domain-mesh-model-variable names, move some redundant code into functions
#see if vtk memory leak complaints need addressing
#see if something is wrong with quadratics
#see if something is wrong with 3d pointSet
#
useMainWindow = False#True
#
#Utilities and high level functions
#
class Window:
    def __init__(self,name,title):
        import proteus.Comm
        comm = proteus.Comm.get()
        self.hardCopies=0
        #mwf vtk on laptop needs update
        #import pdb
        #pdb.set_trace()
        skipComm = False
        self.comm = comm
        self.compManager = vtkCompositeRenderManager()
        if not skipComm:
            self.communicator = vtkMPICommunicator()
            self.controller = vtkMPIController()
            self.controller.SetCommunicator(self.communicator.GetWorldCommunicator())
            self.compManager.SetController(self.controller)
            if useCoPro:
                self.copro = initializeCoProcessor(self.comm,self.controller)
        self.myProcId = comm.rank()
        self.numProcs = comm.size()
        self.isMaster = comm.isMaster()
        self.background=(1,1,1)
        #self.background=(0,0,0)
        self.textColor=(0,0,0)
        self.name=name
        #vtk object dictionary
        self.vod={}
        if comm.rank() == 0:
            if hasQt:
                if useMainWindow:
                    #Qt widgets
                    self.frameWidget = QtGui.QFrame(g.mainWindow)
                    self.hbox = QtGui.QHBoxLayout()
                    self.iren = QVTKRenderWindowInteractor(self.frameWidget)
                    if comm.size() > 1:
                        self.iren.Disable() 
                    else:
                        self.iren.SetInteractorStyle(vtkInteractorStyleTrackballCamera())
                    #self.renWin = self.iren.GetRenderWindow()
                    self.renWin = self.compManager.MakeRenderWindow()
                    self.renWin.SetWindowName(name)
                    self.iren.SetRenderWindow(self.renWin)
                    self.iren.Initialize()
                    self.hbox.addWidget(self.iren)
                    self.frameWidget.setLayout(self.hbox)
                    g.tabWidget.addTab(self.frameWidget,title)
                    g.tabWidget.setCurrentWidget(self.frameWidget)
                    screen = QtGui.QDesktopWidget().screenGeometry()
                    size = g.mainWindow.geometry()
                    (x,y) = (screen.width()-comm.size()*size.width())/2, (screen.height()-size.height())/2
                    g.mainWindow.move(x+comm.rank()*size.width(),y)
                    g.mainWindow.show()
                else:
                    #self.iren = vtkRenderWindowInteractor()#QVTKRenderWindowInteractor()
                    self.iren = QVTKRenderWindowInteractor()
                    if comm.size() > 1:
                        self.iren.Disable() 
                    else:
                        self.iren.SetInteractorStyle(vtkInteractorStyleTrackballCamera())
                    #self.renWin = self.iren.GetRenderWindow()
                    self.renWin = self.compManager.MakeRenderWindow()
                    self.iren.SetRenderWindow(self.renWin)
                    self.renWin.SetWindowName(name)
                    self.iren.Initialize()
                    self.iren.show()
            else:
                self.ren = vtkRenderer()
                self.iren = vtkRenderWindowInteractor()
                self.iren.SetInteractorStyle(vtkInteractorStyleTrackballCamera())
                self.renWin = self.compManager.MakeRenderWindow()
                self.renWin.AddRenderer(self.ren)
                self.iren.SetRenderWindow(self.renWin)
                self.renWin.SetWindowName(name)
                self.iren.Initialize()
        else:
            self.renWin = self.compManager.MakeRenderWindow()
            if not hasQt:
                self.ren = vtkRenderer()
                self.renWin.AddRenderer(self.ren)
            self.renWin.OffScreenRenderingOn()

windowDict={}
imageDict={}
        
class vtkGlobals:
    def __init__(self):
        import proteus.Comm
        comm = proteus.Comm.get()
        if comm.rank() == 0:
            if hasQt:
                if useMainWindow:
                    self.app = QtGui.QApplication(sys.argv)
                    self.mainWindow = QtGui.QMainWindow()
                    self.mainWindow.setWindowTitle('Proteus')
                    self.tabWidget = QtGui.QTabWidget(self.mainWindow)
                    self.frameWidgetDict={}
                    self.hboxDict={}
                    self.mainWindow.setCentralWidget(self.tabWidget)
                else:
                    self.app = QtGui.QApplication(sys.argv)

class ModelDataStorage:
    def __init__(self):
        self.CellArray = []
        self.ModelData = None
        self.VectorData = None
        self.VectorDataArrayPointer = None
        
# Module global variables
global g

g = vtkGlobals()

def createRenderers(viewTypes,window):   
    nViewports = len(viewTypes)
    nx = int(math.sqrt(nViewports))
    dx = 1.0/float(nx)
    ny = nViewports/nx + nViewports%nx
    dy = 1.0/float(ny)
    nV=0
    for i in range(nx):
        for j in range(ny):
            ren = window.compManager.MakeRenderer()
            ren.SetViewport(i*dx,
                            1.0 - (j+1)*dy,
                            (i+1)*dx,
                            1.0 - j*dy)
            ren.SetBackground(window.background)
            if nV == nViewports-1:
                ren.SetViewport(i*dx,
                                0.0,
                                1.0,
                                1.0 - j*dy)
            window.vod['ren_'+viewTypes[nV]] = ren
            window.renWin.AddRenderer(ren)
            nV+=1
            if nV == nViewports:
                break
    if window.numProcs > 1:
        window.compManager.SetRenderWindow(window.renWin)
        window.compManager.InitializePieces()

def createLegendActor(lut):
    legend = vtk.vtkScalarBarActor()
    legend.SetLookupTable(lut)
    legend.SetOrientationToHorizontal()
    legend.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
    legend.GetPositionCoordinate().SetValue(0.1,0.01)
    legend.SetWidth(0.8)
    legend.SetHeight(0.1)
    legend.SetNumberOfLabels(5)
    # Set up the legend text
    tp = legend.GetTitleTextProperty()
    tp.SetFontSize(12)
    tp.SetFontFamilyToArial()
    tp.SetColor(0,0,0)
    tp = legend.GetLabelTextProperty()
    tp.SetFontSize(12)
    tp.SetFontFamilyToArial()
    tp.SetColor(0,0,0)
    return legend

def ViewMesh(mesh,title="Mesh",viewMaterialTypes=True,hardcopy=False):
    ##\todo improve legend
    ##      add probe filter for querying material types
    ##      1d mesh visualization
    #mwf debug
    #import pdb
    #pdb.set_trace()
    #put in separate routine for 1d
    if mesh.nNodes_element == 2:
        return
    import cvtkviewers
    global g,windowDict
    windowName = title
    window = Window(windowName,title)
    windowDict[window.name] = window
    window.vod['vtkMesh'] = cvtkviewers.getUnstructuredGridFromMesh(mesh.nodeArray,
                                                                    mesh.elementNodesArray)
    window.vod['vtkMesh'].Update()
    #crude viewing for  now
    createRenderers(['wireframe'],window)
    ren = window.vod['ren_wireframe']
    #ren = vtk.vtkRenderer()
    #renWin = vtk.vtkRenderWindow()
    #renWin.SetWindowName(title)
    #ren.SetBackground(1, 1, 1)
    #renWin.SetSize(500, 500)
    #renWin.AddRenderer(ren)
    #iren = QVTKRenderWindowInteractor()
    #iren.SetInteractorStyle(vtkInteractorStyleTrackballCamera())
    #iren.SetRenderWindow(renWin)
    if viewMaterialTypes and mesh.elementMaterialTypes != None:
        window.vod['matData'] = vtk.vtkIntArray()
        window.vod['matData'].SetNumberOfComponents(1)
        window.vod['matData'].SetNumberOfValues(mesh.nElements_global)
        saveArray=1
        window.vod['matData'].SetVoidArray(mesh.elementMaterialTypes,mesh.nElements_global,saveArray)
        #copy in manually ... 
        #for eN in range(mesh.nElements_global):
        #    window.vod['matData'].InsertValue(eN,mesh.elementMaterialTypes[eN])
        window.vod['vtkMesh'].GetCellData().SetScalars(window.vod['matData'])
        window.vod['vtkMesh'].Update()

        mrange = window.vod['vtkMesh'].GetCellData().GetScalars().GetRange()
        window.vod['meshMapper'] = vtk.vtkDataSetMapper()
        window.vod['meshMapper'].SetInputData(window.vod['vtkMesh'])
        window.vod['meshMapper'].SetScalarRange(mrange)
        window.vod['meshActor'] = vtk.vtkActor()
        window.vod['meshActor'].SetMapper(window.vod['meshMapper'])
        window.vod['meshActor'].GetProperty().SetRepresentationToWireframe()
    
        window.vod['meshActor'].GetProperty().SetDiffuseColor(1,0,0)
        ren.AddActor(window.vod['meshActor'])
        
        # Create the LookupTable for the Mesh
        window.vod['lut'] = vtk.vtkLookupTable()
        window.vod['lut'].SetTableRange(mrange[0], mrange[1])
        window.vod['lut'].Build()
        #
        window.vod['legend'] = vtk.vtkScalarBarActor()
        window.vod['legend'].SetLookupTable(window.vod['lut'])
        window.vod['legend'].SetTitle("Material Flags")
        window.vod['legend'].SetOrientationToVertical()
        window.vod['legend'].GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
        window.vod['legend'].GetPositionCoordinate().SetValue(0.9,0.25)
        window.vod['legend'].SetWidth(0.1)
        window.vod['legend'].SetHeight(0.75)
        window.vod['legend'].SetNumberOfLabels(int(mrange[1]-mrange[0]+1))
        # Set up the window.vod['legend'] text
        window.vod['tp'] = window.vod['legend'].GetTitleTextProperty()
        window.vod['tp'].SetFontSize(12)
        window.vod['tp'].SetFontFamilyToArial()
        window.vod['tp'].SetColor(0,0,0)
        window.vod['tp'] = window.vod['legend'].GetLabelTextProperty()
        window.vod['tp'].SetFontSize(12)
        window.vod['tp'].SetFontFamilyToArial()   
        window.vod['tp'].SetColor(0,0,0)

        ren.AddActor(window.vod['legend'])
    else:
        #view mesh this way
        window.vod['meshMapper'] = vtk.vtkDataSetMapper()
        window.vod['meshMapper'].SetInputData(window.vod['vtkMesh'])
        window.vod['meshActor'] = vtk.vtkActor()
        window.vod['meshActor'].SetMapper(window.vod['meshMapper'])
        window.vod['meshActor'].GetProperty().SetRepresentationToWireframe()
    
        window.vod['meshActor'].GetProperty().SetDiffuseColor(1,0,0)
        ren.AddActor(window.vod['meshActor'])
    if window.myProcId == 0:
        if window.numProcs > 1:
            window.compManager.ResetAllCameras()
        window.renWin.Render()
        if hasQt:
            g.app.processEvents()
        window.compManager.StopServices()
    else:
        window.compManager.StartServices()
    
def ViewBoundaryMesh(mesh,title="Boundary Mesh",viewBoundaryMaterialTypes=True,hardcopy=False):
    ##\todo improve legend
    ##      add probe filter for querying material types
    ##      1d mesh visualization
    #mwf debug
    #import pdb
    #pdb.set_trace()
    #put in separate routine for 1d
    if mesh.nNodes_element == 2:
        return
    import cvtkviewers
    global g,windowDict
    #windowName = mesh.domain.name+mesh.name+variableName
    windowName = title
    t = 0.0
    window = Window(windowName,title)
    windowDict[window.name] = window
    window.vod['pdata'] = cvtkviewers.getPolyDataBoundaryMesh(mesh.nodeArray,mesh.elementBoundaryNodesArray,
                                                              mesh.elementBoundaryMaterialTypes)
    window.vod['pdata'].Update()
    brange = [0,1]
    #
    if  viewBoundaryMaterialTypes and mesh.elementBoundaryMaterialTypes != None:
        brange = window.vod['pdata'].GetCellData().GetScalars().GetRange()
    #crude viewing for  now
    createRenderers(['colormapped'],window)
    ren = window.vod['ren_colormapped']
    #view mesh this way
    window.vod['meshMapper'] = vtk.vtkPolyDataMapper()
    window.vod['meshMapper'].SetInputData(window.vod['pdata'])
    window.vod['meshMapper'].SetScalarRange(brange)
    window.vod['meshActor'] = vtk.vtkActor()
    window.vod['meshActor'].SetMapper(window.vod['meshMapper'])
    ren.AddActor(window.vod['meshActor'])
    if viewBoundaryMaterialTypes and mesh.elementBoundaryMaterialTypes != None and brange[1] > brange[0]:
        # Create the LookupTable for the Mesh
        window.vod['lut'] = vtk.vtkLookupTable()
        window.vod['lut'].SetTableRange(brange[0], brange[1])
        window.vod['lut'].ForceBuild()
        #
        window.vod['legend'] = vtk.vtkScalarBarActor()
        window.vod['legend'].SetLookupTable(window.vod['lut'])
        window.vod['legend'].SetTitle("Boundary Flags")
        window.vod['legend'].SetOrientationToVertical()
        window.vod['legend'].GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
        window.vod['legend'].GetPositionCoordinate().SetValue(0.9,0.25)
        window.vod['legend'].SetWidth(0.1)
        window.vod['legend'].SetHeight(0.75)
        window.vod['legend'].SetNumberOfLabels(int(brange[1]-brange[0]+1))
        # Set up the legend text
        window.vod['tp'] = window.vod['legend'].GetTitleTextProperty()
        window.vod['tp'].SetFontSize(12)
        window.vod['tp'].SetFontFamilyToArial()
        window.vod['tp'].SetColor(0,0,0)
        window.vod['tp'] = window.vod['legend'].GetLabelTextProperty()
        window.vod['tp'].SetFontSize(12)
        window.vod['tp'].SetFontFamilyToArial()   
        window.vod['tp'].SetColor(0,0,0)
        ren.AddActor(window.vod['legend'])
    #
    if window.myProcId == 0:
        if window.numProcs > 1:
            window.compManager.ResetAllCameras()
        window.renWin.Render()
        if hasQt:
            g.app.processEvents()
        window.compManager.StopServices()
    else:
        window.compManager.StartServices()
#
# 1D Scalars
#
def viewScalar_1D(xVals, yVals, xTitle, yTitle, title, winNum, Pause = True, sortPoints=True, Hardcopy = False):
    """
    Display an XY Plot
    
    Arguments:
        xVals:  Array of X values for the XY plot.
        yVals:  Array of Y values for the XY plot.
        xTitle: String.  Title of X axis.
        yTitle: String.  Title of Y axis.
        title:  String.  Title of the chart.
        winNum:  Integer.  Window number.  Used like the window number argument for gnuplot windows.  If window number is the same as one that
        Pause:  Boolean.  Specify whether to open VTK in interactive mode (if True) or just display mode.
    mwf: added sort option because nodes not in ascending order anymore
    """
    import cvtkviewers
    global windowDict
    windowName = "1D Mesh"+title
    viewTypes=['plot']
    t=0
    if not windowDict.has_key(windowName):
        window=Window(windowName,title)
        windowDict[window.name] = window
        createRenderers(viewTypes,window)
        window.xSorted = numpy.zeros(xVals.shape,'d')
        window.ySorted = numpy.zeros(xVals.shape,'d')
        #sort later
        window.vod['xArray'] = cvtkviewers.prepareScalarValueArray(window.xSorted)
        window.vod['yArray'] = cvtkviewers.prepareScalarValueArray(window.ySorted)
        # Create the Data Object for the xyPlot class to interpret
        window.vod['fieldData'] = vtk.vtkFieldData()
        window.vod['fieldData'].AllocateArrays(2)
        window.vod['fieldData'].AddArray(window.vod['xArray'])
        window.vod['fieldData'].AddArray(window.vod['yArray'])
        window.vod['dataSet'] = vtk.vtkDataObject()
        window.vod['dataSet'].SetFieldData(window.vod['fieldData'])   
        # Create an xy-plot
        # The x-values we are plotting are the underlying point data values.
        window.vod['xyPlotActor'] = vtk.vtkXYPlotActor()
        window.vod['xyPlotActor'].AddDataObjectInput(window.vod['dataSet']) 
        window.vod['xyPlotActor'].SetXValuesToValue()
        window.vod['xyPlotActor'].SetDataObjectXComponent(0, 0)
        window.vod['xyPlotActor'].SetDataObjectYComponent(0, 1)
        window.vod['xyPlotActor'].GetPositionCoordinate().SetValue(0.05, 0.05, 0)
        window.vod['xyPlotActor'].GetPosition2Coordinate().SetValue(0.9, 0.9, 0) #relative to Position    
        window.vod['xyPlotActor'].SetNumberOfXLabels(6)
        window.vod['xyPlotActor'].SetTitle(title)
        window.vod['xyPlotActor'].SetXTitle(xTitle)
        window.vod['xyPlotActor'].SetYTitle(yTitle)
        window.vod['xyPlotActor'].PlotPointsOn()
        window.vod['xyPlotActor'].GetProperty().SetColor(0, 0, 1)
        window.vod['xyPlotActor'].GetProperty().SetPointSize(3)
        # Set text prop color (same color for backward compat with test)
        # Assign same object to all text props
        tprop = window.vod['xyPlotActor'].GetTitleTextProperty()
        tprop.SetColor(0,0,0)
        window.vod['xyPlotActor'].SetAxisTitleTextProperty(tprop)
        window.vod['xyPlotActor'].SetAxisLabelTextProperty(tprop)   
        window.vod['ren_plot'].AddActor2D(window.vod['xyPlotActor'])
        xMax = max(xVals)
        yMax = max(yVals)
        xMin = min(xVals)
        yMin = min(yVals)
        Ly = yMax-yMin
        window.xMax = xMax
        window.yMax = yMax + 0.05*Ly
        window.xMin = xMin
        window.yMin = yMin - 0.05*Ly
        window.vod['xyPlotActor'].SetXRange(window.xMin, window.xMax)
        window.vod['xyPlotActor'].SetYRange(window.yMin, window.yMax)
    else:
        window = windowDict[windowName]
    # Find the minimum and maximum data values to set the axis min/maxes from
    xMax = max(xVals)
    yMax = max(yVals)
    xMin = min(xVals)
    yMin = min(yVals)
    if (xMax > window.xMax or 
        xMin < window.xMin):
        window.xMax = xMax
        window.xMin = xMin
        window.vod['xyPlotActor'].SetXRange(window.xMin, window.xMax)
    if (yMax > window.yMax or
        yMin < window.yMin):
        Ly = yMax-yMin
        window.yMax = yMax + 0.05*Ly
        window.yMin = yMin - 0.05*Ly
        window.vod['xyPlotActor'].SetYRange(window.yMin, window.yMax)
    isort = xVals.argsort()
    window.xSorted[:]=xVals[isort]
    window.ySorted[:]=yVals[isort]
    window.vod['xArray'].Modified()
    window.vod['yArray'].Modified()
    if window.myProcId == 0:
        if window.numProcs > 1:
            window.compManager.ResetAllCameras()
        window.renWin.Render()
        if hasQt:
            g.app.processEvents()
        window.compManager.StopServices()
    else:
        window.compManager.StartServices()

#
#2D Scalars
#
#def viewScalar_tri3_2D(mesh, scalars, title, winNum, viewTypes=['colorMapped','contour','warp'],IsoSurface = True, Pause = True, Hardcopy = False):
def viewScalar_tri3_2D(mesh, scalars, title, winNum, viewTypes=['colorMapped','hardcopy'],IsoSurface = True, Pause = True, Hardcopy = False, Adapted = False):
    #
    #build an unstructured grid data set and pass to viewScalar_2DMesh
    #
    #Domain name (from mesh)
    #Mesh name (supplied)
    #Model name (supplied)
    #Variable name (supplied)
    #Coefficient name
    #
    #=> Unique Window key
    #views types => unique viewport key
    #
    #
    import cvtkviewers
    global windowDict
    #windowName = mesh.domain.name+mesh.name+variableName
    mesh.name = "Triangular mesh"
    windowName = "Triangular Mesh"+title
    t = 0.0
    if not windowDict.has_key(windowName):
        windowCreated=True
        window = Window(windowName,title)
        windowDict[window.name] = window
        #data set
        window.vod['dataSet'] = cvtkviewers.getUnstructuredGridFromMesh(mesh.nodeArray, 
                                                                        mesh.elementNodesArray)
        window.vod['scalars'] = cvtkviewers.prepareScalarValueArray(scalars)
        window.vod['dataSet'].GetPointData().SetScalars(window.vod['scalars'])
        #mapper
        window.vod['dataSetMapper'] = vtk.vtkDataSetMapper()
        window.vod['dataSetMapper'].SetInputData(window.vod['dataSet'])
        window.vod['dataSetMapper'].SetScalarVisibility(1)
        window.vod['lut'] = vtk.vtkLookupTable()
        range = window.vod['scalars'].GetRange()
        range = (flcbdfWrappers.globalMin(range[0]),flcbdfWrappers.globalMax(range[1]))
        #window.vod['lut'].SetTableRange(window.vod['scalars'].GetRange())
        window.vod['lut'].SetTableRange(range)
        window.vod['lut'].SetHueRange(0.66667,0.0)
        window.vod['dataSetMapper'].SetLookupTable(window.vod['lut'])
        window.vod['dataSetMapper'].UseLookupTableScalarRangeOn()
    elif Adapted:
        windowCreated=False
        window = windowDict[windowName]
        #data set
        window.vod['dataSet'] = cvtkviewers.getUnstructuredGridFromMesh(mesh.nodeArray, 
                                                                        mesh.elementNodesArray)
        window.vod['scalars'] = cvtkviewers.prepareScalarValueArray(scalars)
        window.vod['dataSet'].GetPointData().SetScalars(window.vod['scalars'])
        #mapper
        window.vod['dataSetMapper'].SetInputData(window.vod['dataSet'])
        window.vod['scalars'].Modified()
        window.vod['lut'].SetTableRange(window.vod['scalars'].GetRange())
    else:
        windowCreated=False
        window = windowDict[windowName]
        window.vod['scalars'].Modified()
        #window.vod['lut'].SetTableRange(window.vod['scalars'].GetRange())
        range = window.vod['scalars'].GetRange()
        range = (flcbdfWrappers.globalMin(range[0]),flcbdfWrappers.globalMax(range[1]))
        #window.vod['lut'].SetTableRange(window.vod['scalars'].GetRange())
        window.vod['lut'].SetTableRange(range)
    #window and viewports 
    viewScalar_2D(window,
                  windowCreated,
                  viewTypes,
                  Adapted = Adapted)

def viewScalar_tri6_2D(mesh, dofMap, scalars, title, winNum, viewTypes=['colorMapped'],IsoSurface = True, Pause = True, Hardcopy = False):
    #
    #build an unstructured grid data set and pass to viewScalar_2DMesh
    #
    #Domain name (from mesh)
    #Mesh name (supplied)
    #Model name (supplied)
    #Variable name (supplied)
    #Coefficient name
    #
    #=> Unique Window key
    #views types => unique viewport key
    #
    #
    import cvtkviewers
    global windowDict
    #windowName = mesh.domain.name+mesh.name+variableName
    mesh.name = "Quadratic Triangular mesh"
    windowName = "Quadratic Triangular Mesh"+title
    t = 0.0
    if not windowDict.has_key(windowName):
        windowCreated=True
        window = Window(windowName,title)
        windowDict[window.name] = window
        #data set
        #mwf change to match new def of lagrangeNodesArray (all points)
        window.nodeArray = dofMap.lagrangeNodesArray#numpy.append(mesh.nodeArray, dofMap.lagrangeNodesArray)
        window.vod['dataSet'] = cvtkviewers.getUnstructuredGridFromQuadraticTriangleMesh(mesh.nNodes_global, 
                                                                                         len(dofMap.lagrangeNodesArray),#mesh.nNodes_global + len(dofMap.lagrangeNodesArray),
                                                                                         mesh.nElements_global, 
                                                                                         window.nodeArray, 
                                                                                         dofMap.l2g, 
                                                                                         mesh.edgeNodesArray)
        window.vod['scalars'] = cvtkviewers.prepareScalarValueArray(scalars)
        window.vod['dataSet'].GetPointData().SetScalars(window.vod['scalars'])
        #mapper
        window.vod['dataSetMapper'] = vtk.vtkDataSetMapper()
        window.vod['dataSetMapper'].SetInputData(window.vod['dataSet'])
        window.vod['dataSetMapper'].SetScalarVisibility(1)
        window.vod['lut'] = vtk.vtkLookupTable()
        window.vod['lut'].SetTableRange(window.vod['scalars'].GetRange())
        window.vod['lut'].SetHueRange(0.66667,0.0)
        window.vod['dataSetMapper'].SetLookupTable(window.vod['lut'])
        window.vod['dataSetMapper'].UseLookupTableScalarRangeOn()
    else:
        windowCreated=False
        window = windowDict[windowName]
        window.vod['scalars'].Modified()
        window.vod['lut'].SetTableRange(window.vod['scalars'].GetRange())
    #window and viewports 
    viewScalar_2D(window,
                  windowCreated,
                  viewTypes)

def viewScalar_pointSet_2D(nodes, scalars, title, winNum,IsoSurface = True, Pause = True, Hardcopy = False,viewTypes=['colorMapped']):#,'contour','warp']):
    """
    Creates a mesh from the x,y,z point locations passed in and calls vtkDisplay2DScalarMesh to display it
    
    Arguments:
        nodes:       Double[]   Array of doubles in X, Y, Z, X, Y, Z format for node locations
        scalars:     Double[].  Array of scalar values, one value for each node
        title:       String.  Title of the plot.
        winNum:      Integer.  Window number.  Used like the window number argument for gnuplot windows.  
                     If window number is the same as one that already exists, it will update the existing window.  
                     If not, a new window is created.
        IsoSurface:  Boolean.  Specify whether or not to display the isosurface.  Optional.  Default is true.
        Pause:       Boolean.  Specify whether to open VTK in interactive mode (if True) or just display mode.
    """
    import cvtkviewers
    global windowDict
    #windowName = mesh.domain.name+mesh.name+variableName
    windowName = "Quadrature Mesh"+title
    t = 0.0
    if not windowDict.has_key(windowName):
        windowCreated=True
        window = Window(windowName,title)
        windowDict[windowName] = window
        window.vod['scalars'] = cvtkviewers.prepareScalarValueArray(scalars)
        window.vod['points'] = cvtkviewers.prepareVTKPoints3(nodes)
	# Triangulate the points
        window.vod['polyData'] = vtk.vtkPolyData()
        window.vod['polyData'].SetPoints(window.vod['points'])
	window.vod['delny'] = vtk.vtkDelaunay2D()
	window.vod['delny'].SetInputData(window.vod['polyData'])
	window.vod['delny'].SetTolerance(0.001)
        window.vod['polyData'] = window.vod['delny'].GetOutput()
        window.vod['polyData'].Update()
        #form the mesh
        window.vod['cells']= window.vod['polyData'].GetPolys()
        window.vod['dataSet'] = vtk.vtkUnstructuredGrid()
        window.vod['dataSet'].SetCells(5,window.vod['cells'])
        window.vod['dataSet'].SetPoints(window.vod['points'])
	window.vod['dataSet'].GetPointData().SetScalars(window.vod['scalars'])
        #mapper
        window.vod['dataSetMapper'] = vtk.vtkDataSetMapper()
        window.vod['dataSetMapper'].SetInputData(window.vod['dataSet'])
        window.vod['dataSetMapper'].SetScalarVisibility(1)
        window.vod['lut'] = vtk.vtkLookupTable()
        window.vod['lut'].SetTableRange(window.vod['scalars'].GetRange())
        window.vod['lut'].SetHueRange(0.66667,0.0)
        window.vod['dataSetMapper'].SetLookupTable(window.vod['lut'])
        window.vod['dataSetMapper'].UseLookupTableScalarRangeOn()
    else:
        windowCreated=False
        window=windowDict[windowName]
        window.vod['scalars'].Modified()
        window.vod['lut'].SetTableRange(window.vod['scalars'].GetRange())
    viewScalar_2D(window,windowCreated,viewTypes)

def viewScalar_2D(window,windowCreated,viewTypes,Adapted=False,Hardcopy=True):
    global g,windowDict
    if windowCreated:
        createRenderers(viewTypes,window)
        #this should be a loop over viewTypes
        if 'colorMapped' in viewTypes:
            window.vod['gridActor_colorMapped'] = vtk.vtkActor()
            window.vod['gridActor_colorMapped'].SetMapper(window.vod['dataSetMapper'])
            ren = window.vod['ren_colorMapped']
            if window.isMaster:
                window.vod['axesActor_colorMapped'] = vtk.vtkCubeAxesActor2D()
                window.vod['axesActor_colorMapped'].SetInputData(window.vod['dataSet'])
                window.vod['axesActor_colorMapped'].SetCamera(window.vod['ren_colorMapped'].GetActiveCamera())
                window.vod['axesActor_colorMapped'].SetFlyModeToClosestTriad()
                window.vod['axesActor_colorMapped'].SetZAxisVisibility(0)
                tp=window.vod['axesActor_colorMapped'].GetAxisTitleTextProperty()
                tp.SetFontSize(14)
                tp.SetFontFamilyToArial()
                tp.SetColor(0,0,0)
                tp=window.vod['axesActor_colorMapped'].GetAxisLabelTextProperty()
                tp.SetFontSize(14)
                tp.SetFontFamilyToArial()
                tp.SetColor(0,0,0)
                window.vod['legendActor_colorMapped']= createLegendActor(window.vod['lut'])
                window.vod['legendActor_colorMapped'].SetLookupTable(window.vod['lut'])
                ren.AddActor(window.vod['legendActor_colorMapped'])
                ren.AddActor2D(window.vod['axesActor_colorMapped'])
            ren.AddActor(window.vod['gridActor_colorMapped'])
            ren.ResetCamera()
        if 'contour' in viewTypes:
            window.vod['algo_contour'] = vtk.vtkContourGrid()
            window.vod['algo_contour'].SetInputData(window.vod['dataSet'])
            window.vod['algo_contour'].GenerateValues(10,window.vod['dataSet'].GetPointData().GetScalars().GetRange())            
            window.vod['contourMapper'] = vtk.vtkPolyDataMapper()
            window.vod['contourMapper'].SetInputConnection(window.vod['algo_contour'].GetOutputPort())
            #window.vod['contourMapper'].SetScalarRange(window.vod['dataSet'].GetPointData().GetScalars().GetRange())
            window.vod['contourMapper'].SetLookupTable(window.vod['lut'])
            window.vod['contourMapper'].UseLookupTableScalarRangeOn()
            window.vod['actorcontour'] = vtk.vtkActor()
            window.vod['actorcontour'].SetMapper(window.vod['contourMapper'])
            window.vod['axesActor_contour'] = vtk.vtkCubeAxesActor2D()
            window.vod['axesActor_contour'].SetInputData(window.vod['dataSet'])
            window.vod['axesActor_contour'].SetCamera(window.vod['ren_contour'].GetActiveCamera())
            window.vod['axesActor_contour'].SetFlyModeToClosestTriad()
            window.vod['axesActor_contour'].SetZAxisVisibility(0)
            tp=window.vod['axesActor_contour'].GetAxisTitleTextProperty()
            tp.SetFontSize(14)
            tp.SetFontFamilyToArial()
            tp.SetColor(0,0,0)
            tp=window.vod['axesActor_contour'].GetAxisLabelTextProperty()
            tp.SetFontSize(14)
            tp.SetFontFamilyToArial()
            tp.SetColor(0,0,0)
            window.vod['legendActor_contour']= createLegendActor(window.vod['lut'])
            window.vod['legendActor_contour'].SetLookupTable(window.vod['lut'])
            ren = window.vod['ren_contour']
            ren.AddActor(window.vod['actorcontour'])
            ren.AddActor(window.vod['legendActor_contour'])
            ren.AddActor2D(window.vod['axesActor_contour'])
            ren.ResetCamera()
        if 'warp' in viewTypes:
            window.vod['algo_geo'] = vtk.vtkGeometryFilter()
            window.vod['algo_geo'].SetInputData(window.vod['dataSet'])
            window.vod['algo_warp'] = vtk.vtkWarpScalar()
            window.vod['algo_warp'].SetInputConnection(window.vod['algo_geo'].GetOutputPort())
            window.vod['warpMapper'] = vtk.vtkPolyDataMapper()
            window.vod['warpMapper'].SetInputConnection(window.vod['algo_warp'].GetOutputPort())
            window.vod['warpMapper'].SetLookupTable(window.vod['lut'])
            window.vod['warpMapper'].UseLookupTableScalarRangeOn()
            window.vod['warpActor'] = vtk.vtkActor()
            window.vod['warpActor'].SetMapper(window.vod['warpMapper'])
            window.vod['axesActor_warp'] = vtk.vtkCubeAxesActor2D()
            window.vod['axesActor_warp'].SetInputData(window.vod['dataSet'])
            window.vod['axesActor_warp'].SetCamera(window.vod['ren_warp'].GetActiveCamera())
            window.vod['axesActor_warp'].SetFlyModeToClosestTriad()
            window.vod['axesActor_warp'].SetZAxisVisibility(0)
            tp=window.vod['axesActor_warp'].GetAxisTitleTextProperty()
            tp.SetFontSize(14)
            tp.SetFontFamilyToArial()
            tp.SetColor(0,0,0)
            tp=window.vod['axesActor_warp'].GetAxisLabelTextProperty()
            tp.SetFontSize(14)
            tp.SetFontFamilyToArial()
            tp.SetColor(0,0,0)
            window.vod['legendActor_warp']= createLegendActor(window.vod['lut'])
            window.vod['legendActor_warp'].SetLookupTable(window.vod['lut'])
            ren = window.vod['ren_warp']
            ren.AddActor(window.vod['warpActor'])
            ren.AddActor(window.vod['legendActor_warp'])
            ren.AddActor2D(window.vod['axesActor_warp'])
            ren.ResetCamera()
        if Hardcopy:
            if window.myProcId == 0:
                window.hardCopies = 0
                # exp = vtk.vtkGL2PSExporter()
                # exp.SetRenderWindow(window.renWin)
                # #exp.SetFileFormatToTeX()
                # #exp.OcclusionCullOff()
                # exp.SetFilePrefix(window.name)
                # # Turn off compression so PIL can read file.
                # exp.CompressOff() 
                # exp.SetSortToBSP()
                # exp.DrawBackgroundOff()
                # exp.SimpleLineOffsetOn()
                # exp.BestRootOn()
                # exp.TextOff()
                # #
                # exp.Write()
    elif Adapted:
        if 'colorMapped' in viewTypes:
            window.vod['gridActor_colorMapped'].SetMapper(window.vod['dataSetMapper'])
            window.vod['axesActor_colorMapped'].SetInputData(window.vod['dataSet'])
            window.vod['legendActor_colorMapped'].SetLookupTable(window.vod['lut'])
            ren = window.vod['ren_colorMapped']
            ren.ResetCamera()
        if 'contour' in viewTypes:
            window.vod['algo_contour'].SetInputData(window.vod['dataSet'])
            window.vod['algo_contour'].GenerateValues(10,window.vod['dataSet'].GetPointData().GetScalars().GetRange())            
            window.vod['contourMapper'].SetInputConnection(window.vod['algo_contour'].GetOutputPort())
            window.vod['contourMapper'].SetLookupTable(window.vod['lut'])
            window.vod['axesActor_contour'].SetInputData(window.vod['dataSet'])
            ren = window.vod['ren_contour']
            ren.ResetCamera()
        if 'warp' in viewTypes:
            window.vod['algo_geo'].SetInputData(window.vod['dataSet'])
            window.vod['algo_warp'].SetInputConnection(window.vod['algo_geo'].GetOutputPort())
            window.vod['warpMapper'].SetInputConnection(window.vod['algo_warp'].GetOutputPort())
            window.vod['warpMapper'].SetLookupTable(window.vod['lut'])
            window.vod['axesActor_warp'].SetInputData(window.vod['dataSet'])
            ren = window.vod['ren_warp']
            ren.ResetCamera()

    if 'contour' in viewTypes:
        window.vod['algo_contour'].GenerateValues(10,window.vod['dataSet'].GetPointData().GetScalars().GetRange())
    if window.myProcId == 0:
        if window.numProcs > 1:
            window.compManager.ResetAllCameras()
        window.renWin.Render()
        if Hardcopy:
            window.vod['w2if'] = vtk.vtkWindowToImageFilter()
            window.vod['writer'] = vtk.vtkPNGWriter()
            window.vod['w2if'].SetInput(window.renWin)
            window.vod['writer'].SetInputConnection(window.vod['w2if'].GetOutputPort())   
            filename = window.name+`window.hardCopies`+'.png'
            window.vod['writer'].SetFileName(filename)
            window.hardCopies+=1
            window.vod['writer'].Write()        
            window.png = open(filename).read()
        if hasQt:
            #g.app.sendPostedEvents()
            #g.app.flush()
            #g.app.processEvents(QtCore.QEventLoop.ExcludeUserInputEvents)
            #g.app.processEvents(QtCore.QEventLoop.ExcludeSocketNotifiers)
            g.app.processEvents()
        window.compManager.StopServices()
        #
    else:
        window.compManager.StartServices()
    if useCoPro:
        #coprocess
        from vtk import vtkCoProcessorPython
        import samplecoprocessingscript as cpscript
        datadescription = \
            vtkCoProcessorPython.vtkCPDataDescription()
        #to do --- pass down real value of time and integer step number
        cptime=0.0
        cpstep=0
        #if step == lastStep-1:
        #    datadescription.SetForceOuputOn()#or something like that
        datadescription.SetTimeData(cptime, cpstep)
        datadescription.AddInput("input")
        cpscript.RequestDataDescription(datadescription)
        inputdescription = \
            datadescription.GetInputDescriptionByName("input")
        if inputdescription.GetIfGridIsNecessary() == False:
            return
        inputdescription.SetGrid(windowDict[window.name].vod['dataSet'])
        cpscript.DoCoProcessing(datadescription)

#
#2D Vectors
#
def viewVector_tri3_2D(mesh, u,v, title,IsoSurface = True, Pause = True, Hardcopy = False,
                       viewTypes=['colorMapped'],Adapted=False):#'arrows','streamlines','colorMapped','contour','warp']):
    #
    #build an unstructured grid data set and pass to viewScalar_2DMesh
    #
    #Domain name (from mesh)
    #Mesh name (supplied)
    #Model name (supplied)
    #Variable name (supplied)
    #Coefficient name
    #
    #=> Unique Window key
    #views types => unique viewport key
    #
    #
    import cvtkviewers
    global windowDict
    #windowName = mesh.domain.name+mesh.name+variableName
    mesh.name = "Triangular mesh"
    windowName = "Triangular Mesh"+title
    t = 0.0
    if not windowDict.has_key(windowName):
        windowCreated=True
        window = Window(windowName,title)
        windowDict[window.name] = window
        #data set
        window.vod['dataSet'] = cvtkviewers.getUnstructuredGridFromMesh(mesh.nodeArray, 
                                                                        mesh.elementNodesArray)
        window.w=numpy.zeros(v.shape,'d')
        window.vectors=numpy.column_stack((u,v,window.w)).flatten()
        window.vod['vectors'] = cvtkviewers.prepareVectorValueArray(window.vectors)
        window.vod['dataSet'].GetPointData().SetVectors(window.vod['vectors'])
        #mapper
        window.vod['normFilter'] = vtk.vtkVectorNorm()
        window.vod['normFilter'].SetInputData(window.vod['dataSet'])
        window.vod['normFilter'].SetAttributeModeToUsePointData()
        window.vod['normFilter'].Update()
        window.vod['dataSetMapper'] = vtk.vtkDataSetMapper()
        #window.vod['dataSetMapper'].SetInputConnection(window.vod['normFilter'].GetOutputPort())
        window.vod['dataSetMapper'].SetInputData(window.vod['normFilter'].GetOutput())
        window.vod['dataSetMapper'].SetScalarVisibility(1)
        window.vod['lut'] = vtk.vtkLookupTable()
        #window.vod['lut'].SetTableRange(window.vod['normFilter'].GetOutput().GetPointData().GetScalars().GetRange())
        range = window.vod['normFilter'].GetOutput().GetPointData().GetScalars().GetRange()
        range = (flcbdfWrappers.globalMin(range[0]),flcbdfWrappers.globalMax(range[1]))
        window.vod['lut'].SetTableRange(range)
        window.vod['lut'].SetHueRange(0.66667,0.0)
        window.vod['dataSetMapper'].SetLookupTable(window.vod['lut'])
        window.vod['dataSetMapper'].UseLookupTableScalarRangeOn()
    elif Adapted:
        windowCreated=True
        window = windowDict[windowName]
        window.w = numpy.zeros(v.shape,'d')
        window.vectors = numpy.column_stack((u,v,window.w)).flatten()
        window.vod['vectors'] = cvtkviewers.prepareVectorValueArray(window.vectors)
        window.vod['dataSet'].GetPointData().SetVectors(window.vod['vectors'])
        #mapper
        window.vod['normFilter'].SetInputData(window.vod['dataSet'])
        window.vod['normFilter'].Update()
        window.vod['dataSetMapper'].SetInputConnection(window.vod['normFilter'].GetOutputPort())
        window.vod['lut'].SetTableRange(window.vod['normFilter'].GetOutput().GetPointData().GetScalars().GetRange())
        window.vod['dataSetMapper'].SetLookupTable(window.vod['lut'])
    else:
        windowCreated=False
        window = windowDict[windowName]
        window.vectors[:] = numpy.column_stack((u,v,window.w)).flatten()
        window.vod['vectors'].Modified()
        #window.vod['lut'].SetTableRange(window.vod['normFilter'].GetOutput().GetPointData().GetScalars().GetRange())
        range = window.vod['normFilter'].GetOutput().GetPointData().GetScalars().GetRange()
        range = (flcbdfWrappers.globalMin(range[0]),flcbdfWrappers.globalMax(range[1]))
        window.vod['lut'].SetTableRange(range)
    #window and viewports 
    viewVector_2D(window,
                  windowCreated,
                  viewTypes)

def viewVector_tri6_2D(mesh, dofMap,u,v, title,IsoSurface = True, Pause = True, Hardcopy = False,
                       viewTypes=['streamlines']):#'arrows','streamlines','colorMapped','contour','warp']):
    #
    #build an unstructured grid data set and pass to viewScalar_2DMesh
    #
    #Domain name (from mesh)
    #Mesh name (supplied)
    #Model name (supplied)
    #Variable name (supplied)
    #Coefficient name
    #
    #=> Unique Window key
    #views types => unique viewport key
    #
    #
    import cvtkviewers
    global windowDict
    #windowName = mesh.domain.name+mesh.name+variableName
    mesh.name = "Quadratic Triangular mesh"
    windowName = "Quadratic Triangular Mesh"+title
    t = 0.0
    if not windowDict.has_key(windowName):
        windowCreated=True
        window = Window(windowName,title)
        windowDict[window.name] = window
        #data set
        #mwf change to match new def of lagrangeNodesArray (all points)
        window.nodeArray = dofMap.lagrangeNodesArray#numpy.append(mesh.nodeArray, dofMap.lagrangeNodesArray)
        window.vod['dataSet'] = cvtkviewers.getUnstructuredGridFromQuadraticTriangleMesh(mesh.nNodes_global, 
                                                                                         len(dofMap.lagrangeNodesArray),#mesh.nNodes_global + len(dofMap.lagrangeNodesArray),
                                                                                         mesh.nElements_global, 
                                                                                         window.nodeArray, 
                                                                                         dofMap.l2g, 
                                                                                         mesh.edgeNodesArray)
        window.w=numpy.zeros(v.shape,'d')
        window.vectors=numpy.column_stack((u,v,window.w)).flatten()
        window.vod['vectors'] = cvtkviewers.prepareVectorValueArray(window.vectors)
        window.vod['dataSet'].GetPointData().SetVectors(window.vod['vectors'])
        #mapper
        window.vod['normFilter'] = vtk.vtkVectorNorm()
        window.vod['normFilter'].SetInputData(window.vod['dataSet'])
        window.vod['normFilter'].SetAttributeModeToUsePointData()
        window.vod['normFilter'].Update()
        window.vod['dataSetMapper'] = vtk.vtkDataSetMapper()
        window.vod['dataSetMapper'].SetInputConnection(window.vod['normFilter'].GetOutputPort())
        window.vod['dataSetMapper'].SetScalarVisibility(1)
        window.vod['lut'] = vtk.vtkLookupTable()
        window.vod['lut'].SetTableRange(window.vod['normFilter'].GetOutput().GetPointData().GetScalars().GetRange())
        window.vod['lut'].SetHueRange(0.66667,0.0)
        window.vod['dataSetMapper'].SetLookupTable(window.vod['lut'])
        window.vod['dataSetMapper'].UseLookupTableScalarRangeOn()
    else:
        windowCreated=False
        window = windowDict[windowName]
        window.vectors[:] = numpy.column_stack((u,v,window.w)).flatten()
        window.vod['vectors'].Modified()
        window.vod['lut'].SetTableRange(window.vod['normFilter'].GetOutput().GetPointData().GetScalars().GetRange())
    #window and viewports 
    viewVector_2D(window,
                  windowCreated,
                  viewTypes)

def viewVector_pointSet_2D(nodes, vectors, title,IsoSurface = True, Pause = True, Hardcopy = False,
                           viewTypes=['streamlines']):#,'arrows','streamlines','contour','warp']):
    #
    #build an unstructured grid data set and pass to viewScalar_2DMesh
    #
    #Domain name (from mesh)
    #Mesh name (supplied)
    #Model name (supplied)
    #Variable name (supplied)
    #Coefficient name
    #
    #=> Unique Window key
    #views types => unique viewport key
    #
    #
    import cvtkviewers
    global windowDict
    #windowName = mesh.domain.name+mesh.name+variableName
    #mesh.name = "Triangular mesh"
    windowName = "Quadrature Mesh"+title
    t = 0.0
    if not windowDict.has_key(windowName):
        windowCreated=True
        window = Window(windowName,title)
        windowDict[window.name] = window
        window.w=numpy.zeros(vectors.shape[:-1],'d')
        window.vectors=numpy.column_stack((vectors.flat[::2],vectors.flat[1::2],window.w.flat)).flatten()
        window.vod['vectors'] = cvtkviewers.prepareVectorValueArray(window.vectors)
        window.vod['points'] = cvtkviewers.prepareVTKPoints3(nodes)
	# Triangulate the points
        window.vod['polyData'] = vtk.vtkPolyData()
        window.vod['polyData'].SetPoints(window.vod['points'])
	window.vod['delny'] = vtk.vtkDelaunay2D()
	window.vod['delny'].SetInputData(window.vod['polyData'])
	window.vod['delny'].SetTolerance(0.001)
        window.vod['polyData'] = window.vod['delny'].GetOutput()
        window.vod['polyData'].Update()
        #form the mesh
        window.vod['cells']= window.vod['polyData'].GetPolys()
        window.vod['dataSet'] = vtk.vtkUnstructuredGrid()
        window.vod['dataSet'].SetCells(5,window.vod['cells'])
        window.vod['dataSet'].SetPoints(window.vod['points'])
	window.vod['dataSet'].GetPointData().SetVectors(window.vod['vectors'])
        #mapper
        window.vod['normFilter'] = vtk.vtkVectorNorm()
        window.vod['normFilter'].SetInputData(window.vod['dataSet'])
        window.vod['normFilter'].SetAttributeModeToUsePointData()
        window.vod['normFilter'].Update()
        window.vod['dataSetMapper'] = vtk.vtkDataSetMapper()
        window.vod['dataSetMapper'].SetInputConnection(window.vod['normFilter'].GetOutputPort())
        window.vod['dataSetMapper'].SetScalarVisibility(1)
        window.vod['lut'] = vtk.vtkLookupTable()
        window.vod['lut'].SetTableRange(window.vod['normFilter'].GetOutput().GetPointData().GetScalars().GetRange())
        window.vod['lut'].SetHueRange(0.66667,0.0)
        window.vod['dataSetMapper'].SetLookupTable(window.vod['lut'])
        window.vod['dataSetMapper'].UseLookupTableScalarRangeOn()
    else:
        windowCreated=False
        window = windowDict[windowName]
        window.vectors.flat[:]=numpy.column_stack((vectors.flat[::2],vectors.flat[1::2],window.w.flat)).flatten()
        window.vod['vectors'].Modified()
        window.vod['lut'].SetTableRange(window.vod['normFilter'].GetOutput().GetPointData().GetScalars().GetRange())
    #window and viewports 
    viewVector_2D(window,
                  windowCreated,
                  viewTypes)

def viewVector_2D(window,
                  windowCreated,
                  viewTypes):
    global g,windowDict
    import pdb
    if windowCreated:
        createRenderers(viewTypes,window)
        #this should be a loop over viewTypes
        if 'arrows' in viewTypes:
            window.vod['glyph'] = vtk.vtkGlyph3D()
            window.vod['glyph'].SetInputData(window.vod['dataSet'])
            window.vod['glyph'].SetScaleModeToScaleByVector()
            window.vod['arrow'] = vtk.vtkArrowSource()
            window.vod['glyph'].SetSource(window.vod['arrow'].GetOutput())
            window.vod['glyphMapper'] = vtk.vtkPolyDataMapper()
            window.vod['glyphMapper'].SetInput(window.vod['glyph'].GetOutput())
            #window.vod['glyphMapper'].SetLookupTable(window.vod['lut'])
            #window.vod['glyphMapper'].UseLookupTableScalarRangeOn()
            window.vod['glyphActor'] = vtk.vtkActor()
            window.vod['glyphActor'].SetMapper(window.vod['glyphMapper'])
            window.vod['glyphActor'].GetProperty().SetColor(0,1,0)
            ren = window.vod['ren_arrows']
            ren.AddActor(window.vod['glyphActor'])
            ren.ResetCamera()
        if 'streamlines' in viewTypes:
            #
            seeds = vtk.vtkPolyData()
            window.vod['seeds'] = seeds
            #
            lineWidget = vtk.vtkLineWidget()
            window.vod['lineWidget'] = lineWidget
            lineWidget.SetInputData(window.vod['dataSet'])
            lineWidget.SetAlignToYAxis()
            lineWidget.SetResolution(50)
            lineWidget.PlaceWidget()
            lineWidget.GetPolyData(seeds)
            lineWidget.ClampToBoundsOn()
            #
            streamer = vtk.vtkStreamTracer()
            window.vod['streamLineFilter'] = streamer
            streamer.SetInputData(window.vod['dataSet'])
            streamer.SetIntegratorTypeToRungeKutta45()
            streamer.SetIntegrationDirectionToBoth()
            streamer.SetSource(seeds)
            #
            window.vod['rk4'] = vtk.vtkRungeKutta4()
            streamer.SetIntegrator(window.vod['rk4'])
            #
            window.vod['streamLineMapper']=vtk.vtkPolyDataMapper()
            window.vod['streamLineMapper'].SetInputConnection(streamer.GetOutputPort())
            #
            window.vod['streamLineActor'] = vtk.vtkActor()
            streamline = window.vod['streamLineActor']
            window.vod['streamLineActor'].SetMapper(window.vod['streamLineMapper'])
            window.vod['streamLineActor'].VisibilityOn()
            #if window.isMaster:
            #    lineWidget.SetInteractor(window.iren)
            window.vod['gridActor_colorMapped'] = vtk.vtkActor()
            window.vod['gridActor_colorMapped'].SetMapper(window.vod['dataSetMapper'])
            window.vod['axesActor_colorMapped'] = vtk.vtkCubeAxesActor2D()
            window.vod['axesActor_colorMapped'].SetInputData(window.vod['dataSet'])
            window.vod['axesActor_colorMapped'].SetCamera(window.vod['ren_streamlines'].GetActiveCamera())
            window.vod['axesActor_colorMapped'].SetFlyModeToClosestTriad()
            window.vod['axesActor_colorMapped'].SetZAxisVisibility(0)
            tp=window.vod['axesActor_colorMapped'].GetAxisTitleTextProperty()
            tp.SetFontSize(14)
            tp.SetFontFamilyToArial()
            tp.SetColor(0,0,0)
            tp=window.vod['axesActor_colorMapped'].GetAxisLabelTextProperty()
            tp.SetFontSize(14)
            tp.SetFontFamilyToArial()
            tp.SetColor(0,0,0)
            window.vod['legendActor_colorMapped']= createLegendActor(window.vod['lut'])
            window.vod['legendActor_colorMapped'].SetLookupTable(window.vod['lut'])
            ren = window.vod['ren_streamlines']
            ren.AddActor(window.vod['gridActor_colorMapped'])
            ren.AddActor(window.vod['legendActor_colorMapped'])
            ren.AddActor2D(window.vod['axesActor_colorMapped'])
            ren.AddActor(window.vod['streamLineActor'])
            ren.ResetCamera()
        if 'colorMapped' in viewTypes:
            window.vod['gridActor_colorMapped'] = vtk.vtkActor()
            window.vod['gridActor_colorMapped'].SetMapper(window.vod['dataSetMapper'])
            ren = window.vod['ren_colorMapped']
            if window.isMaster:
                window.vod['axesActor_colorMapped'] = vtk.vtkCubeAxesActor2D()
                window.vod['axesActor_colorMapped'].SetInputData(window.vod['dataSet'])
                window.vod['axesActor_colorMapped'].SetCamera(window.vod['ren_colorMapped'].GetActiveCamera())
                window.vod['axesActor_colorMapped'].SetFlyModeToClosestTriad()
                window.vod['axesActor_colorMapped'].SetZAxisVisibility(0)
                tp=window.vod['axesActor_colorMapped'].GetAxisTitleTextProperty()
                tp.SetFontSize(14)
                tp.SetFontFamilyToArial()
                tp.SetColor(0,0,0)
                tp=window.vod['axesActor_colorMapped'].GetAxisLabelTextProperty()
                tp.SetFontSize(14)
                tp.SetFontFamilyToArial()
                tp.SetColor(0,0,0)
                window.vod['legendActor_colorMapped']= createLegendActor(window.vod['lut'])
                window.vod['legendActor_colorMapped'].SetLookupTable(window.vod['lut'])
                ren.AddActor(window.vod['legendActor_colorMapped'])
                ren.AddActor2D(window.vod['axesActor_colorMapped'])
            ren.AddActor(window.vod['gridActor_colorMapped'])
            ren.ResetCamera()
        if 'contour' in viewTypes:
            window.vod['algo_contour'] = vtk.vtkContourGrid()
            window.vod['algo_contour'].SetInput(window.vod['normFilter'].GetOutput())
            window.vod['algo_contour'].GenerateValues(10,window.vod['normFilter'].GetOutput().GetPointData().GetScalars().GetRange())            
            window.vod['contourMapper'] = vtk.vtkPolyDataMapper()
            window.vod['contourMapper'].SetInputConnection(window.vod['algo_contour'].GetOutputPort())
            window.vod['contourMapper'].SetLookupTable(window.vod['lut'])
            window.vod['contourMapper'].UseLookupTableScalarRangeOn()
            window.vod['actorcontour'] = vtk.vtkActor()
            window.vod['actorcontour'].SetMapper(window.vod['contourMapper'])
            window.vod['axesActor_contour'] = vtk.vtkCubeAxesActor2D()
            window.vod['axesActor_contour'].SetInputData(window.vod['dataSet'])
            window.vod['axesActor_contour'].SetCamera(window.vod['ren_contour'].GetActiveCamera())
            window.vod['axesActor_contour'].SetFlyModeToClosestTriad()
            window.vod['axesActor_contour'].SetZAxisVisibility(0)
            tp=window.vod['axesActor_contour'].GetAxisTitleTextProperty()
            tp.SetFontSize(14)
            tp.SetFontFamilyToArial()
            tp.SetColor(0,0,0)
            tp=window.vod['axesActor_contour'].GetAxisLabelTextProperty()
            tp.SetFontSize(14)
            tp.SetFontFamilyToArial()
            tp.SetColor(0,0,0)
            window.vod['legendActor_contour']= createLegendActor(window.vod['lut'])
            window.vod['legendActor_contour'].SetLookupTable(window.vod['lut'])
            ren = window.vod['ren_contour']
            ren.AddActor(window.vod['actorcontour'])
            ren.AddActor(window.vod['legendActor_contour'])
            ren.AddActor2D(window.vod['axesActor_contour'])
            ren.ResetCamera()
        if 'warp' in viewTypes:
            window.vod['algo_geo'] = vtk.vtkGeometryFilter()
            window.vod['algo_geo'].SetInput(window.vod['normFilter'].GetOutput())
            window.vod['algo_warp'] = vtk.vtkWarpScalar()
            window.vod['algo_warp'].SetInputConnection(window.vod['algo_geo'].GetOutputPort())
            window.vod['warpMapper'] = vtk.vtkPolyDataMapper()
            window.vod['warpMapper'].SetInputConnection(window.vod['algo_warp'].GetOutputPort())
            window.vod['warpMapper'].SetLookupTable(window.vod['lut'])
            window.vod['warpMapper'].UseLookupTableScalarRangeOn()
            window.vod['warpActor'] = vtk.vtkActor()
            window.vod['warpActor'].SetMapper(window.vod['warpMapper'])
            window.vod['axesActor_warp'] = vtk.vtkCubeAxesActor2D()
            window.vod['axesActor_warp'].SetInputData(window.vod['dataSet'])
            window.vod['axesActor_warp'].SetCamera(window.vod['ren_warp'].GetActiveCamera())
            window.vod['axesActor_warp'].SetFlyModeToClosestTriad()
            window.vod['axesActor_warp'].SetZAxisVisibility(0)
            tp=window.vod['axesActor_warp'].GetAxisTitleTextProperty()
            tp.SetFontSize(14)
            tp.SetFontFamilyToArial()
            tp.SetColor(0,0,0)
            tp=window.vod['axesActor_warp'].GetAxisLabelTextProperty()
            tp.SetFontSize(14)
            tp.SetFontFamilyToArial()
            tp.SetColor(0,0,0)
            window.vod['legendActor_warp']= createLegendActor(window.vod['lut'])
            window.vod['legendActor_warp'].SetLookupTable(window.vod['lut'])
            ren = window.vod['ren_warp']
            ren.AddActor(window.vod['warpActor'])
            ren.AddActor(window.vod['legendActor_warp'])
            ren.AddActor2D(window.vod['axesActor_warp'])
            ren.ResetCamera()
    if 'contour' in viewTypes:
        window.vod['algo_contour'].GenerateValues(10,window.vod['normFilter'].GetOutput().GetPointData().GetScalars().GetRange())
    if 'streamlines' in viewTypes:
        window.vod['lineWidget'].GetPolyData(window.vod['seeds'])
    if window.myProcId == 0:
        if window.numProcs > 1:
            window.compManager.ResetAllCameras()
        window.renWin.Render()
        if hasQt:
            #g.app.sendPostedEvents()
            #g.app.flush()
            #g.app.processEvents(QtCore.QEventLoop.ExcludeUserInputEvents)
            g.app.processEvents()
        window.compManager.StopServices()
    else:
        window.compManager.StartServices()
#
#3D Scalars
#
#def viewScalar_tet4_3D(mesh, scalars, title, winNum, viewTypes=['colorMapped','contour','warp'],IsoSurface = True, Pause = True, Hardcopy = False):
def viewScalar_tet4_3D(mesh, scalars, title, winNum, viewTypes=['colorMapped'],IsoSurface = True, Pause = True, Hardcopy = False, Adapted = False):
    #
    #build an unstructured grid data set and pass to viewScalar_3DMesh
    #
    #Domain name (from mesh)
    #Mesh name (supplied)
    #Model name (supplied)
    #Variable name (supplied)
    #Coefficient name
    #
    #=> Unique Window key
    #views types => unique viewport key
    #
    #
    import cvtkviewers
    global windowDict
    #windowName = mesh.domain.name+mesh.name+variableName
    mesh.name = "Tetrahedral Mesh"
    windowName = "Tetrahedral Mesh"+title
    t = 0.0
    if not windowDict.has_key(windowName):
        windowCreated=True
        window = Window(windowName,title)
        windowDict[window.name] = window
        window.vod['dataSet'] = cvtkviewers.getUnstructuredGridFromMesh(mesh.nodeArray, 
                                                                        mesh.elementNodesArray)
        window.vod['scalars'] = cvtkviewers.prepareScalarValueArray(scalars)
        window.vod['dataSet'].GetPointData().SetScalars(window.vod['scalars'])
        #mapper
        window.vod['dataSetMapper'] = vtk.vtkDataSetMapper()
        window.vod['dataSetMapper'].SetInputData(window.vod['dataSet'])
        window.vod['dataSetMapper'].SetScalarVisibility(1)
        window.vod['lut'] = vtk.vtkLookupTable()
        window.vod['lut'].SetTableRange(window.vod['scalars'].GetRange())
        window.vod['lut'].SetHueRange(0.66667,0.0)
        window.vod['dataSetMapper'].SetLookupTable(window.vod['lut'])
        window.vod['dataSetMapper'].UseLookupTableScalarRangeOn()
    elif Adapted:
        windowCreated=False
        window = windowDict[windowName]
        #data set
        window.vod['dataSet'] = cvtkviewers.getUnstructuredGridFromMesh(mesh.nodeArray, 
                                                                        mesh.elementNodesArray)
        window.vod['scalars'] = cvtkviewers.prepareScalarValueArray(scalars)
        window.vod['dataSet'].GetPointData().SetScalars(window.vod['scalars'])
        #mapper
        window.vod['dataSetMapper'].SetInputData(window.vod['dataSet'])
        window.vod['scalars'].Modified()
        window.vod['lut'].SetTableRange(window.vod['scalars'].GetRange())
    else:
        windowCreated=False
        window = windowDict[windowName]
        window.vod['scalars'].Modified()
        window.vod['lut'].SetTableRange(window.vod['scalars'].GetRange())
    #window and viewports 
    viewScalar_3D(window,
                  windowCreated,
                  viewTypes,
                  Adapted = Adapted)

def viewScalar_tet10_3D(mesh, dofMap, scalars, title, winNum, viewTypes=['colorMapped'],IsoSurface = True, Pause = True, Hardcopy = False):
    #
    #build an unstructured grid data set and pass to viewScalar_3DMesh
    #
    #Domain name (from mesh)
    #Mesh name (supplied)
    #Model name (supplied)
    #Variable name (supplied)
    #Coefficient name
    #
    #=> Unique Window key
    #views types => unique viewport key
    #
    #
    import cvtkviewers
    global windowDict
    #windowName = mesh.domain.name+mesh.name+variableName
    mesh.name = "Tetrahedral Mesh"
    windowName = "Tetrahedral Mesh"+title
    t = 0.0
    if not windowDict.has_key(windowName):
        windowCreated=True
        window = Window(windowName,title)
        windowDict[window.name] = window
        #mwf change to match new def of lagrangeNodesArray (all points)
        window.nodeArray = dofMap.lagrangeNodesArray#numpy.append(mesh.nodeArray, dofMap.lagrangeNodesArray)
        window.vod['dataSet'] = cvtkviewers.getUnstructuredGridFromQuadraticTetMesh(mesh.nNodes_global, 
                                                                                    len(dofMap.lagrangeNodesArray),#mesh.nNodes_global + len(dofMap.lagrangeNodesArray),
                                                                                    mesh.nElements_global, 
                                                                                    window.nodeArray, 
                                                                                    dofMap.l2g, 
                                                                                    mesh.edgeNodesArray)
        window.vod['scalars'] = cvtkviewers.prepareScalarValueArray(scalars)
        window.vod['dataSet'].GetPointData().SetScalars(window.vod['scalars'])
        #mapper
        window.vod['dataSetMapper'] = vtk.vtkDataSetMapper()
        window.vod['dataSetMapper'].SetInputData(window.vod['dataSet'])
        window.vod['dataSetMapper'].SetScalarVisibility(1)
        window.vod['lut'] = vtk.vtkLookupTable()
        window.vod['lut'].SetTableRange(window.vod['scalars'].GetRange())
        window.vod['lut'].SetHueRange(0.66667,0.0)
        window.vod['dataSetMapper'].SetLookupTable(window.vod['lut'])
        window.vod['dataSetMapper'].UseLookupTableScalarRangeOn()
    else:
        windowCreated=False
        window = windowDict[windowName]
        window.vod['scalars'].Modified()
        window.vod['lut'].SetTableRange(window.vod['scalars'].GetRange())
    #window and viewports 
    viewScalar_3D(window,
                  windowCreated,
                  viewTypes)

def viewScalar_pointSet_3D(nodes, scalars, title, winNum,IsoSurface = True, Pause = True, Hardcopy = False,viewTypes=['colorMapped']):#,'contour','warp']):
    import cvtkviewers
    import pdb,gc
    global windowDict
    #windowName = mesh.domain.name+mesh.name+variableName
    windowName = "Quadrature Mesh"+title
    t = 0.0
    if not windowDict.has_key(windowName):
        windowCreated=True
        window = Window(windowName,title)
        windowDict[windowName] = window
        window.vod['scalars'] = cvtkviewers.prepareScalarValueArray(scalars)
        window.vod['points'] = cvtkviewers.prepareVTKPoints3(nodes)
	# Triangulate the points
        window.vod['polyData'] = vtk.vtkPolyData()
        window.vod['polyData'].SetPoints(window.vod['points'])
	window.vod['delny'] = vtk.vtkDelaunay3D()
	window.vod['delny'].SetInputData(window.vod['polyData'])
	window.vod['delny'].SetTolerance(0.001)
        #form the mesh
        window.vod['dataSet'] = window.vod['delny'].GetOutput()
        window.vod['dataSet'].Update()
        #window.vod['cells']= window.vod['polyData'].GetPolys()
        #window.vod['dataSet'] = vtk.vtkUnstructuredGrid()
        #window.vod['dataSet'].SetCells(5,window.vod['cells'])
        #window.vod['dataSet'].SetPoints(window.vod['points'])
	window.vod['dataSet'].GetPointData().SetScalars(window.vod['scalars'])
        #mapper
        window.vod['dataSetMapper'] = vtk.vtkDataSetMapper()
        window.vod['dataSetMapper'].SetInputData(window.vod['dataSet'])
        window.vod['dataSetMapper'].SetScalarVisibility(1)
        window.vod['lut'] = vtk.vtkLookupTable()
        window.vod['lut'].SetTableRange(window.vod['scalars'].GetRange())
        window.vod['lut'].SetHueRange(0.66667,0.0)
        window.vod['dataSetMapper'].SetLookupTable(window.vod['lut'])
        window.vod['dataSetMapper'].UseLookupTableScalarRangeOn()
    else:
        windowCreated=False
        window=windowDict[windowName]
        window.vod['scalars'].Modified()
        window.vod['lut'].SetTableRange(window.vod['scalars'].GetRange())
    viewScalar_3D(window,windowCreated,viewTypes)

def viewScalar_3D(window,windowCreated,viewTypes,Adapted=False):
    global g,windowDict
    if windowCreated:
        createRenderers(viewTypes,window)
        #this should be a loop over viewTypes
        if 'colorMapped' in viewTypes:
            window.vod['gridActor_colorMapped'] = vtk.vtkActor()
            window.vod['gridActor_colorMapped'].SetMapper(window.vod['dataSetMapper'])
            window.vod['axesActor_colorMapped'] = vtk.vtkCubeAxesActor2D()
            window.vod['axesActor_colorMapped'].SetInputData(window.vod['dataSet'])
            window.vod['axesActor_colorMapped'].SetCamera(window.vod['ren_colorMapped'].GetActiveCamera())
            window.vod['axesActor_colorMapped'].SetFlyModeToClosestTriad()
            tp=window.vod['axesActor_colorMapped'].GetAxisTitleTextProperty()
            tp.SetFontSize(14)
            tp.SetFontFamilyToArial()
            tp.SetColor(0,0,0)
            tp=window.vod['axesActor_colorMapped'].GetAxisLabelTextProperty()
            tp.SetFontSize(14)
            tp.SetFontFamilyToArial()
            tp.SetColor(0,0,0)
            window.vod['legendActor_colorMapped']= createLegendActor(window.vod['lut'])
            window.vod['legendActor_colorMapped'].SetLookupTable(window.vod['lut'])
            window.vod['legendActor_colorMapped'].SetTitle('')
            ren = window.vod['ren_colorMapped']
            ren.AddActor(window.vod['gridActor_colorMapped'])
            ren.AddActor(window.vod['legendActor_colorMapped'])
            ren.AddActor2D(window.vod['axesActor_colorMapped'])
            ren.ResetCamera()
        if 'contour' in viewTypes:
            window.vod['algo_contour'] = vtk.vtkContourGrid()
            window.vod['algo_contour'].SetInputData(window.vod['dataSet'])
            window.vod['algo_contour'].GenerateValues(10,window.vod['dataSet'].GetPointData().GetScalars().GetRange())            
            window.vod['contourMapper'] = vtk.vtkPolyDataMapper()
            window.vod['contourMapper'].SetInputConnection(window.vod['algo_contour'].GetOutputPort())
            #window.vod['contourMapper'].SetScalarRange(window.vod['dataSet'].GetPointData().GetScalars().GetRange())
            window.vod['contourMapper'].SetLookupTable(window.vod['lut'])
            window.vod['contourMapper'].UseLookupTableScalarRangeOn()
            window.vod['actorcontour'] = vtk.vtkActor()
            window.vod['actorcontour'].SetMapper(window.vod['contourMapper'])
            window.vod['axesActor_contour'] = vtk.vtkCubeAxesActor2D()
            window.vod['axesActor_contour'].SetInputData(window.vod['dataSet'])
            window.vod['axesActor_contour'].SetCamera(window.vod['ren_contour'].GetActiveCamera())
            window.vod['axesActor_contour'].SetFlyModeToClosestTriad()
            tp=window.vod['axesActor_contour'].GetAxisTitleTextProperty()
            tp.SetFontSize(14)
            tp.SetFontFamilyToArial()
            tp.SetColor(0,0,0)
            tp=window.vod['axesActor_contour'].GetAxisLabelTextProperty()
            tp.SetFontSize(14)
            tp.SetFontFamilyToArial()
            tp.SetColor(0,0,0)
            window.vod['legendActor_contour']= createLegendActor(window.vod['lut'])
            window.vod['legendActor_contour'].SetLookupTable(window.vod['lut'])
            ren = window.vod['ren_contour']
            ren.AddActor(window.vod['actorcontour'])
            ren.AddActor(window.vod['legendActor_contour'])
            ren.AddActor2D(window.vod['axesActor_contour'])
            ren.ResetCamera()
        if 'warp' in viewTypes:
            window.vod['algo_geo'] = vtk.vtkGeometryFilter()
            window.vod['algo_geo'].SetInputData(window.vod['dataSet'])
            window.vod['algo_warp'] = vtk.vtkWarpScalar()
            window.vod['algo_warp'].SetInputConnection(window.vod['algo_geo'].GetOutputPort())
            window.vod['warpMapper'] = vtk.vtkPolyDataMapper()
            window.vod['warpMapper'].SetInputConnection(window.vod['algo_warp'].GetOutputPort())
            window.vod['warpMapper'].SetLookupTable(window.vod['lut'])
            window.vod['warpMapper'].UseLookupTableScalarRangeOn()
            window.vod['warpActor'] = vtk.vtkActor()
            window.vod['warpActor'].SetMapper(window.vod['warpMapper'])
            window.vod['axesActor_warp'] = vtk.vtkCubeAxesActor2D()
            window.vod['axesActor_warp'].SetInputData(window.vod['dataSet'])
            window.vod['axesActor_warp'].SetCamera(window.vod['ren_warp'].GetActiveCamera())
            window.vod['axesActor_warp'].SetFlyModeToClosestTriad()
            tp=window.vod['axesActor_warp'].GetAxisTitleTextProperty()
            tp.SetFontSize(14)
            tp.SetFontFamilyToArial()
            tp.SetColor(0,0,0)
            tp=window.vod['axesActor_warp'].GetAxisLabelTextProperty()
            tp.SetFontSize(14)
            tp.SetFontFamilyToArial()
            tp.SetColor(0,0,0)
            window.vod['legendActor_warp']= createLegendActor(window.vod['lut'])
            window.vod['legendActor_warp'].SetLookupTable(window.vod['lut'])
            ren = window.vod['ren_warp']
            ren.AddActor(window.vod['warpActor'])
            ren.AddActor(window.vod['legendActor_warp'])
            ren.AddActor2D(window.vod['axesActor_warp'])
            ren.ResetCamera()
    elif Adapted:
        if 'colorMapped' in viewTypes:
            window.vod['gridActor_colorMapped'].SetMapper(window.vod['dataSetMapper'])
            window.vod['axesActor_colorMapped'].SetInputData(window.vod['dataSet'])
            window.vod['legendActor_colorMapped'].SetLookupTable(window.vod['lut'])
            ren = window.vod['ren_colorMapped']
            ren.ResetCamera()
        if 'contour' in viewTypes:
            window.vod['algo_contour'].SetInputData(window.vod['dataSet'])
            window.vod['algo_contour'].GenerateValues(10,window.vod['dataSet'].GetPointData().GetScalars().GetRange())            
            window.vod['contourMapper'].SetInputConnection(window.vod['algo_contour'].GetOutputPort())
            window.vod['contourMapper'].SetLookupTable(window.vod['lut'])
            window.vod['axesActor_contour'].SetInputData(window.vod['dataSet'])
            ren = window.vod['ren_contour']
            ren.ResetCamera()
        if 'warp' in viewTypes:
            window.vod['algo_geo'].SetInputData(window.vod['dataSet'])
            window.vod['algo_warp'].SetInputConnection(window.vod['algo_geo'].GetOutputPort())
            window.vod['warpMapper'].SetInputConnection(window.vod['algo_warp'].GetOutputPort())
            window.vod['warpMapper'].SetLookupTable(window.vod['lut'])
            window.vod['axesActor_warp'].SetInputData(window.vod['dataSet'])
            ren = window.vod['ren_warp']
            ren.ResetCamera()
        
    if 'contour' in viewTypes:
        window.vod['algo_contour'].GenerateValues(10,window.vod['dataSet'].GetPointData().GetScalars().GetRange())
    if window.myProcId == 0:
        if window.numProcs > 1:
            window.compManager.ResetAllCameras()
        window.renWin.Render()
        if hasQt:
            g.app.processEvents()
        window.compManager.StopServices()
    else:
        window.compManager.StartServices()

#
#3D Vectors
#
def viewVector_tet4_3D(mesh, u,v, w,title,IsoSurface = True, Pause = True, Hardcopy = False,
                       viewTypes=['streamlines','arrows','streamlines','colorMapped','contour','warp']):
    #
    #build an unstructured grid data set and pass to viewScalar_2DMesh
    #
    #Domain name (from mesh)
    #Mesh name (supplied)
    #Model name (supplied)
    #Variable name (supplied)
    #Coefficient name
    #
    #=> Unique Window key
    #views types => unique viewport key
    #
    #
    import cvtkviewers
    global windowDict
    #windowName = mesh.domain.name+mesh.name+variableName
    mesh.name = "Tetrahedral mesh"
    windowName = "Tetrahedral Mesh"+title
    t = 0.0
    if not windowDict.has_key(windowName):
        windowCreated=True
        window = Window(windowName,title)
        windowDict[window.name] = window
        #data set
        window.vod['dataSet'] = cvtkviewers.getUnstructuredGridFromMesh(mesh.nodeArray, 
                                                                        mesh.elementNodesArray)
        window.vectors=numpy.column_stack((u,v,w)).flatten()
        window.vod['vectors'] = cvtkviewers.prepareVectorValueArray(window.vectors)
        window.vod['dataSet'].GetPointData().SetVectors(window.vod['vectors'])
        #mapper
        window.vod['normFilter'] = vtk.vtkVectorNorm()
        window.vod['normFilter'].SetInputData(window.vod['dataSet'])
        window.vod['normFilter'].SetAttributeModeToUsePointData()
        window.vod['normFilter'].Update()
        window.vod['dataSetMapper'] = vtk.vtkDataSetMapper()
        window.vod['dataSetMapper'].SetInputConnection(window.vod['normFilter'].GetOutputPort())
        window.vod['dataSetMapper'].SetScalarVisibility(1)
        window.vod['lut'] = vtk.vtkLookupTable()
        window.vod['lut'].SetTableRange(window.vod['normFilter'].GetOutput().GetPointData().GetScalars().GetRange())
        window.vod['lut'].SetHueRange(0.66667,0.0)
        window.vod['dataSetMapper'].SetLookupTable(window.vod['lut'])
        window.vod['dataSetMapper'].UseLookupTableScalarRangeOn()
    else:
        windowCreated=False
        window = windowDict[windowName]
        window.vectors[:] = numpy.column_stack((u,v,w)).flatten()
        window.vod['vectors'].Modified()
        window.vod['lut'].SetTableRange(window.vod['normFilter'].GetOutput().GetPointData().GetScalars().GetRange())
    #window and viewports 
    viewVector_3D(window,
                  windowCreated,
                  viewTypes)

def viewVector_tet10_3D(mesh, dofMap, u,v, w,title,IsoSurface = True, Pause = True, Hardcopy = False,
                       viewTypes=['streamlines','arrows','streamlines','colorMapped','contour','warp']):
    #
    #build an unstructured grid data set and pass to viewScalar_2DMesh
    #
    #Domain name (from mesh)
    #Mesh name (supplied)
    #Model name (supplied)
    #Variable name (supplied)
    #Coefficient name
    #
    #=> Unique Window key
    #views types => unique viewport key
    #
    #
    import cvtkviewers
    global windowDict
    #windowName = mesh.domain.name+mesh.name+variableName
    mesh.name = "Tetrahedral mesh"
    windowName = "Tetrahedral Mesh"+title
    t = 0.0
    if not windowDict.has_key(windowName):
        windowCreated=True
        window = Window(windowName,title)
        windowDict[window.name] = window
        #data set
        #mwf change to match new def of lagrangeNodesArray (all points)
        window.nodeArray = dofMap.lagrangeNodesArray#numpy.append(mesh.nodeArray, dofMap.lagrangeNodesArray)
        window.vod['dataSet'] = cvtkviewers.getUnstructuredGridFromQuadraticTetMesh(mesh.nNodes_global, 
                                                                                    len(dofMap.lagrangeNodesArray),#mesh.nNodes_global + len(dofMap.lagrangeNodesArray),
                                                                                    mesh.nElements_global, 
                                                                                    window.nodeArray, 
                                                                                    dofMap.l2g, 
                                                                                    mesh.edgeNodesArray)
        window.vectors=numpy.column_stack((u,v,w)).flatten()
        window.vod['vectors'] = cvtkviewers.prepareVectorValueArray(window.vectors)
        window.vod['dataSet'].GetPointData().SetVectors(window.vod['vectors'])
        #mapper
        window.vod['normFilter'] = vtk.vtkVectorNorm()
        window.vod['normFilter'].SetInputData(window.vod['dataSet'])
        window.vod['normFilter'].SetAttributeModeToUsePointData()
        window.vod['normFilter'].Update()
        window.vod['dataSetMapper'] = vtk.vtkDataSetMapper()
        window.vod['dataSetMapper'].SetInputConnection(window.vod['normFilter'].GetOutputPort())
        window.vod['dataSetMapper'].SetScalarVisibility(1)
        window.vod['lut'] = vtk.vtkLookupTable()
        window.vod['lut'].SetTableRange(window.vod['normFilter'].GetOutput().GetPointData().GetScalars().GetRange())
        window.vod['lut'].SetHueRange(0.66667,0.0)
        window.vod['dataSetMapper'].SetLookupTable(window.vod['lut'])
        window.vod['dataSetMapper'].UseLookupTableScalarRangeOn()
    else:
        windowCreated=False
        window = windowDict[windowName]
        window.vectors[:] = numpy.column_stack((u,v,w)).flatten()
        window.vod['vectors'].Modified()
        window.vod['lut'].SetTableRange(window.vod['normFilter'].GetOutput().GetPointData().GetScalars().GetRange())
    #window and viewports 
    viewVector_3D(window,
                  windowCreated,
                  viewTypes)

def viewVector_pointSet_3D(nodes, vectors, title,IsoSurface = True, Pause = True, Hardcopy = False,
                           viewTypes=['colorMapped']):#,'arrows','streamlines','contour','warp']):
    #
    #build an unstructured grid data set and pass to viewScalar_2DMesh
    #
    #Domain name (from mesh)
    #Mesh name (supplied)
    #Model name (supplied)
    #Variable name (supplied)
    #Coefficient name
    #
    #=> Unique Window key
    #views types => unique viewport key
    #
    #
    import cvtkviewers
    global windowDict
    #windowName = mesh.domain.name+mesh.name+variableName
    windowName = "Quadrature Mesh"+title
    t = 0.0
    if not windowDict.has_key(windowName):
        windowCreated=True
        window = Window(windowName,title)
        windowDict[window.name] = window
        window.vod['vectors'] = cvtkviewers.prepareVectorValueArray(vectors)
        window.vod['points'] = cvtkviewers.prepareVTKPoints3(nodes)
	# Triangulate the points
        window.vod['polyData'] = vtk.vtkPolyData()
        window.vod['polyData'].SetPoints(window.vod['points'])
        window.vod['polyData'].Update()
	window.vod['delny'] = vtk.vtkDelaunay3D()
	window.vod['delny'].SetInputData(window.vod['polyData'])
	window.vod['delny'].SetTolerance(0.0)
        #form the mesh
        window.vod['dataSet'] = window.vod['delny'].GetOutput()
        window.vod['dataSet'].Update()
        window.vod['dataSet'].SetPoints(window.vod['points'])
	window.vod['dataSet'].GetPointData().SetVectors(window.vod['vectors'])
        #mapper
        window.vod['normFilter'] = vtk.vtkVectorNorm()
        window.vod['normFilter'].SetInputData(window.vod['dataSet'])
        window.vod['normFilter'].SetAttributeModeToUsePointData()
        window.vod['normFilter'].Update()
        window.vod['dataSetMapper'] = vtk.vtkDataSetMapper()
        window.vod['dataSetMapper'].SetInputConnection(window.vod['normFilter'].GetOutputPort())
        window.vod['dataSetMapper'].SetScalarVisibility(1)
        window.vod['lut'] = vtk.vtkLookupTable()
        window.vod['lut'].SetTableRange(window.vod['normFilter'].GetOutput().GetPointData().GetScalars().GetRange())
        window.vod['lut'].SetHueRange(0.66667,0.0)
        window.vod['dataSetMapper'].SetLookupTable(window.vod['lut'])
        window.vod['dataSetMapper'].UseLookupTableScalarRangeOn()
    else:
        windowCreated=False
        window = windowDict[windowName]
        window.vod['vectors'].Modified()
        window.vod['lut'].SetTableRange(window.vod['normFilter'].GetOutput().GetPointData().GetScalars().GetRange())
    #window and viewports 
    viewVector_3D(window,
                  windowCreated,
                  viewTypes)

def viewVector_3D(window,
                  windowCreated,
                  viewTypes):
    global g,windowDict
    import pdb
    if windowCreated:
        createRenderers(viewTypes,window)
        #this should be a loop over viewTypes
        if 'arrows' in viewTypes:
            window.vod['glyph'] = vtk.vtkGlyph3D()
            window.vod['glyph'].SetInputData(window.vod['dataSet'])
            window.vod['glyph'].SetScaleModeToScaleByVector()
            window.vod['arrow'] = vtk.vtkArrowSource()
            window.vod['glyph'].SetSource(window.vod['arrow'].GetOutput())
            window.vod['glyphMapper'] = vtk.vtkPolyDataMapper()
            window.vod['glyphMapper'].SetInput(window.vod['glyph'].GetOutput())
            #window.vod['glyphMapper'].SetLookupTable(window.vod['lut'])
            #window.vod['glyphMapper'].UseLookupTableScalarRangeOn()
            window.vod['glyphActor'] = vtk.vtkActor()
            window.vod['glyphActor'].SetMapper(window.vod['glyphMapper'])
            window.vod['glyphActor'].GetProperty().SetColor(0,1,0)
            ren = window.vod['ren_arrows']
            ren.AddActor(window.vod['glyphActor'])
            ren.ResetCamera()
        if 'streamlines' in viewTypes:
            #
            seeds = vtk.vtkPolyData()
            window.vod['seeds'] = seeds
            #
            planeWidget = vtk.vtkPlaneWidget()
            window.vod['planeWidget'] = planeWidget
            planeWidget.SetInputData(window.vod['dataSet'])
            planeWidget.NormalToXAxisOn()
            planeWidget.SetResolution(10)
            planeWidget.PlaceWidget()
            planeWidget.GetPolyData(seeds)
            #
            streamer = vtk.vtkStreamTracer()
            window.vod['streamLineFilter'] = streamer
            streamer.SetInputData(window.vod['dataSet'])
            streamer.SetIntegratorTypeToRungeKutta45()
            streamer.SetIntegrationDirectionToBoth()
            streamer.SetSource(seeds)
            #
            window.vod['rk4'] = vtk.vtkRungeKutta4()
            streamer.SetIntegrator(window.vod['rk4'])
            #
            tubeFilter = vtk.vtkTubeFilter()
            window.vod['tubeFilter'] = tubeFilter
            tubeFilter.SetInputConnection(streamer.GetOutputPort())
            tubeFilter.SetRadius(0.02)
            window.vod['streamLineMapper']=vtk.vtkPolyDataMapper()
            #window.vod['streamLineMapper'].SetInputConnection(streamer.GetOutputPort())
            window.vod['streamLineMapper'].SetInputConnection(tubeFilter.GetOutputPort())
            #
            window.vod['streamLineActor'] = vtk.vtkActor()
            streamline = window.vod['streamLineActor']
            window.vod['streamLineActor'].SetMapper(window.vod['streamLineMapper'])
            window.vod['streamLineActor'].VisibilityOn()
            #if window.isMaster:
            #    import pdb
            #    pdb.set_trace()
            #    planeWidget.SetInteractor(window.iren)
#             window.vod['gridActor_colorMapped'] = vtk.vtkActor()
#             window.vod['gridActor_colorMapped'].SetMapper(window.vod['dataSetMapper'])
            window.vod['axesActor_colorMapped'] = vtk.vtkCubeAxesActor2D()
            window.vod['axesActor_colorMapped'].SetInputData(window.vod['dataSet'])
            window.vod['axesActor_colorMapped'].SetCamera(window.vod['ren_streamlines'].GetActiveCamera())
            window.vod['axesActor_colorMapped'].SetFlyModeToClosestTriad()
            tp=window.vod['axesActor_colorMapped'].GetAxisTitleTextProperty()
            tp.SetFontSize(14)
            tp.SetFontFamilyToArial()
            tp.SetColor(0,0,0)
            tp=window.vod['axesActor_colorMapped'].GetAxisLabelTextProperty()
            tp.SetFontSize(14)
            tp.SetFontFamilyToArial()
            tp.SetColor(0,0,0)
#             window.vod['legendActor_colorMapped']= createLegendActor(window.vod['lut'])
#             window.vod['legendActor_colorMapped'].SetLookupTable(window.vod['lut'])
            ren = window.vod['ren_streamlines']
            #ren.AddActor(window.vod['gridActor_colorMapped'])
            #ren.AddActor(window.vod['legendActor_colorMapped'])
            ren.AddActor2D(window.vod['axesActor_colorMapped'])
            ren.AddActor(window.vod['streamLineActor'])
            ren.ResetCamera()
        if 'colorMapped' in viewTypes:
            window.vod['gridActor_colorMapped'] = vtk.vtkActor()
            window.vod['gridActor_colorMapped'].SetMapper(window.vod['dataSetMapper'])
            window.vod['axesActor_colorMapped'] = vtk.vtkCubeAxesActor2D()
            window.vod['axesActor_colorMapped'].SetInputData(window.vod['dataSet'])
            window.vod['axesActor_colorMapped'].SetCamera(window.vod['ren_colorMapped'].GetActiveCamera())
            window.vod['axesActor_colorMapped'].SetFlyModeToClosestTriad()
            tp=window.vod['axesActor_colorMapped'].GetAxisTitleTextProperty()
            tp.SetFontSize(14)
            tp.SetFontFamilyToArial()
            tp.SetColor(0,0,0)
            tp=window.vod['axesActor_colorMapped'].GetAxisLabelTextProperty()
            tp.SetFontSize(14)
            tp.SetFontFamilyToArial()
            tp.SetColor(0,0,0)
            window.vod['legendActor_colorMapped']= createLegendActor(window.vod['lut'])
            window.vod['legendActor_colorMapped'].SetLookupTable(window.vod['lut'])
            ren = window.vod['ren_colorMapped']
            ren.AddActor(window.vod['gridActor_colorMapped'])
            ren.AddActor(window.vod['legendActor_colorMapped'])
            ren.AddActor2D(window.vod['axesActor_colorMapped'])
            ren.ResetCamera()
        if 'contour' in viewTypes:
            window.vod['algo_contour'] = vtk.vtkContourGrid()
            window.vod['algo_contour'].SetInput(window.vod['normFilter'].GetOutput())
            window.vod['algo_contour'].GenerateValues(10,window.vod['normFilter'].GetOutput().GetPointData().GetScalars().GetRange())            
            window.vod['contourMapper'] = vtk.vtkPolyDataMapper()
            window.vod['contourMapper'].SetInputConnection(window.vod['algo_contour'].GetOutputPort())
            window.vod['contourMapper'].SetLookupTable(window.vod['lut'])
            window.vod['contourMapper'].UseLookupTableScalarRangeOn()
            window.vod['actorcontour'] = vtk.vtkActor()
            window.vod['actorcontour'].SetMapper(window.vod['contourMapper'])
            window.vod['axesActor_contour'] = vtk.vtkCubeAxesActor2D()
            window.vod['axesActor_contour'].SetInputData(window.vod['dataSet'])
            window.vod['axesActor_contour'].SetCamera(window.vod['ren_contour'].GetActiveCamera())
            window.vod['axesActor_contour'].SetFlyModeToClosestTriad()
            tp=window.vod['axesActor_contour'].GetAxisTitleTextProperty()
            tp.SetFontSize(14)
            tp.SetFontFamilyToArial()
            tp.SetColor(0,0,0)
            tp=window.vod['axesActor_contour'].GetAxisLabelTextProperty()
            tp.SetFontSize(14)
            tp.SetFontFamilyToArial()
            tp.SetColor(0,0,0)
            window.vod['legendActor_contour']= createLegendActor(window.vod['lut'])
            window.vod['legendActor_contour'].SetLookupTable(window.vod['lut'])
            ren = window.vod['ren_contour']
            ren.AddActor(window.vod['actorcontour'])
            ren.AddActor(window.vod['legendActor_contour'])
            ren.AddActor2D(window.vod['axesActor_contour'])
            ren.ResetCamera()
        if 'warp' in viewTypes:
            window.vod['algo_geo'] = vtk.vtkGeometryFilter()
            window.vod['algo_geo'].SetInput(window.vod['normFilter'].GetOutput())
            window.vod['algo_warp'] = vtk.vtkWarpScalar()
            window.vod['algo_warp'].SetInputConnection(window.vod['algo_geo'].GetOutputPort())
            window.vod['warpMapper'] = vtk.vtkPolyDataMapper()
            window.vod['warpMapper'].SetInputConnection(window.vod['algo_warp'].GetOutputPort())
            window.vod['warpMapper'].SetLookupTable(window.vod['lut'])
            window.vod['warpMapper'].UseLookupTableScalarRangeOn()
            window.vod['warpActor'] = vtk.vtkActor()
            window.vod['warpActor'].SetMapper(window.vod['warpMapper'])
            window.vod['axesActor_warp'] = vtk.vtkCubeAxesActor2D()
            window.vod['axesActor_warp'].SetInputData(window.vod['dataSet'])
            window.vod['axesActor_warp'].SetCamera(window.vod['ren_warp'].GetActiveCamera())
            window.vod['axesActor_warp'].SetFlyModeToClosestTriad()
            tp=window.vod['axesActor_warp'].GetAxisTitleTextProperty()
            tp.SetFontSize(14)
            tp.SetFontFamilyToArial()
            tp.SetColor(0,0,0)
            tp=window.vod['axesActor_warp'].GetAxisLabelTextProperty()
            tp.SetFontSize(14)
            tp.SetFontFamilyToArial()
            tp.SetColor(0,0,0)
            window.vod['legendActor_warp']= createLegendActor(window.vod['lut'])
            window.vod['legendActor_warp'].SetLookupTable(window.vod['lut'])
            ren = window.vod['ren_warp']
            ren.AddActor(window.vod['warpActor'])
            ren.AddActor(window.vod['legendActor_warp'])
            ren.AddActor2D(window.vod['axesActor_warp'])
            ren.ResetCamera()
    if 'contour' in viewTypes:
        window.vod['algo_contour'].GenerateValues(10,window.vod['normFilter'].GetOutput().GetPointData().GetScalars().GetRange())
    if 'streamlines' in viewTypes:
        window.vod['planeWidget'].GetPolyData(window.vod['seeds'])
    if window.myProcId == 0:
        if window.numProcs > 1:
            window.compManager.ResetAllCameras()
        window.renWin.Render()
        if hasQt:
            g.app.processEvents()
        window.compManager.StopServices()
    else:
        window.compManager.StartServices()

# def vtkDisplay3DScalarMesh(mesh, scalars, title, winNum, Pause = True, Hardcopy = False):
#     """
#     Display and contour a tetrahedral unstructured grid with scalar data
    
#     Arguments:
#         mesh:            vtkPolyData
#         scalars:         double[]. Scalar values, one value for each node
#         cutPlaneOrigin:  double[3].  A point that exists on the initial cut plane.
#         cutPlaneNormal:  double[3].  Normal vector of the initial cut plane.
#         title:           String.  Title of the plot.
#         winNum:          Integer.  Window number.  Used like the window number argument for gnuplot windows.  
#                          If window number is the same as one that already exists, it will update the existing window.  
#                          If not, a new window is created.
#         IsoSurface:      Boolean.  Specify whether or not to display the isosurface.  Optional.  Default is true.
#         Pause:           Boolean.  Specify whether to open VTK in interactive mode (if True) or just display mode.
#     """

#     global plane, cutter, cutWarper, bounds, textActor, cutWarpActor, dist, isoActor

#     firstTime = True
#     if(WindowExists(winNum)):
# 	    firstTime = False

#     (created,renWin) = GetOrCreateUpdatingWindow(winNum,title)
#     rens = renWin.GetRenderers()
#     ren = rens.GetFirstRenderer()
#     renWin.GetInteractor().ReInitialize()
    
#     # Create the mapper and actor for the main mesh
#     aTetraMapper = vtk.vtkDataSetMapper()
#     aTetraMapper.SetInputData(mesh)
#     aTetraActor = vtk.vtkActor()
#     aTetraActor.SetMapper(aTetraMapper)
#     aTetraActor.GetProperty().SetOpacity(0.1)
#     ren.AddActor(aTetraActor)

#     data = GetMeshData(winNum)
#     data.gridMapper = aTetraMapper

#     # Get the data mins and maxes
#     bounds = mesh.GetBounds()
#     dataRange = scalars.GetRange()
#     xMag = abs(bounds[0]-bounds[1])
#     yMag = abs(bounds[2]-bounds[3])
#     zDim = abs(bounds[4]-bounds[5])
#     zMag = abs(dataRange[0] - dataRange[1])
#     minMag = min([xMag, yMag, zMag])

#     # Create the LookupTable for the Mesh
#     lut = vtk.vtkLookupTable()
#     lut.SetTableRange(dataRange[0], dataRange[1])
#     lut.ForceBuild()
    
#     # Define the plane that will be used to slice the data
#     plane = vtk.vtkPlane()
#     # Origin and Normal will be passed into the funciton eventually
#     #mwf orig
#     #plane.SetOrigin(g.CutPlaneOrigin[0], g.CutPlaneOrigin[1], g.CutPlaneOrigin[2])
#     plane.SetOrigin(g.CutPlaneOrigin[0]*(bounds[1]-bounds[0])+bounds[0],
#                     g.CutPlaneOrigin[1]*(bounds[3]-bounds[2])+bounds[2],
#                     g.CutPlaneOrigin[2]*(bounds[5]-bounds[4])+bounds[4])
#     plane.SetNormal(g.CutPlaneNormal[0], g.CutPlaneNormal[1], g.CutPlaneNormal[2])
    
#     # Define the cutter that will make our slice
#     cutter = vtk.vtkCutter()
#     cutter.SetInputData(mesh)
#     cutter.SetCutFunction(plane)
#     cutMapper = vtk.vtkPolyDataMapper()
#     cutMapper.SetInputData(cutter.GetOutput())
#     cutActor = vtk.vtkActor()
#     cutActor.SetMapper(cutMapper)
#     cutActor.GetProperty().SetOpacity(0.3)
#     ren.AddActor(cutActor)
    
#     # Create the contours for the slice
#     contItems = DisplayContours(cutter.GetOutput(), lut, ren)
#     data.ContourMapper = contItems[1]
#     data.ContourFilter = contItems[2]

    
#     # Define the iso-surface for the slice
#     scale = CalculateIsosurfaceScaleFactor(bounds, dataRange)
#     cutWarper = vtk.vtkWarpScalar()
#     isoActor = DisplayIsosurface(cutter.GetOutput(), plane.GetNormal(), scale, dataRange, ren, cutWarper)
#     data.IsoMapper = isoActor.GetMapper()
#     # Move the isoActor so it is some distance from the cutting plane
#     # Try using halfspace as the distance
#     dist = (dataRange[0]) * scale
#     o = plane.GetOrigin()
#     n = plane.GetNormal()
#     isoActor.SetPosition(0-n[0]*dist,0-n[1]*dist,0-n[2]*dist)
#     #isoActor.SetPosition(o[0]-n[0]*dist,o[1]-n[1]*dist,o[2]-n[2]*dist)

    
#     # Display the plane Data in the bottom corner
#     planeText = "Plane Data:\n   pX = %10.1e \n   pY = %10.1e \n   pZ = %10.1e \n   nX = %10.1e \n   nY = %10.1e \n   nZ = %10.1e" % (plane.GetOrigin()[0],plane.GetOrigin()[1],plane.GetOrigin()[2],plane.GetNormal()[0],plane.GetNormal()[1],plane.GetNormal()[2])
#     # Create a scaled text actor. 
#     # Set the text, font, justification, and properties (bold, italics, etc.).
#     textActor = vtk.vtkTextActor()
#     #DMH: VTK requested that we remove ScaledTextOn() and call SetTextScaleModelToProp() instead.
#     #mwf in case have an older version (say from macports or fink)
#     try:
#         textActor.SetTextScaleModeToProp()
#     except AttributeError:
#         textActor.ScaledTextOn() 
#     textActor.SetInputData(planeText)
#     textActor.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
#     textActor.GetPositionCoordinate().SetValue(0.8,0.02)
#     textActor.SetWidth(0.2)
#     textActor.SetHeight(0.2)
#     # Adjust text properties
#     #tprop = textActor.GetTextProperty()
#     #tprop.SetFontSize(12)
#     #tprop.SetFontFamilyToArial()
#     ren.AddActor(textActor)

#     if (firstTime):
# 	# The callback function to update on widget interaction
# 	def myCallback(obj, event):
# 		obj.GetPlane(plane)
# 		o = plane.GetOrigin()
# 		n = plane.GetNormal()
# 		g.CutPlaneOrigin = o
# 		g.CutPlaneNormal = n
# 		planeText = "Plane Data:\n   pX = %10.1e \n   pY = %10.1e \n   pZ = %10.1e \n   nX = %10.1e \n   nY = %10.1e \n   nZ = %10.1e" % (plane.GetOrigin()[0],plane.GetOrigin()[1],plane.GetOrigin()[2],plane.GetNormal()[0],plane.GetNormal()[1],plane.GetNormal()[2])
# 		textActor.SetInputData(planeText)
# 		cutWarper.SetNormal(n)
# 		isoActor.SetPosition(0-n[0]*dist,0-n[1]*dist,0-n[2]*dist)
		
# 	# Associate the line widget with the interactor
# 	planeWidget = vtk.vtkImplicitPlaneWidget()
# 	renWin.planeWidget = planeWidget
# 	planeWidget.SetOrigin(g.CutPlaneOrigin)
# 	planeWidget.SetNormal(g.CutPlaneNormal)
# 	planeWidget.SetInteractor(renWin.GetInteractor())
# 	planeWidget.SetPlaceFactor(1.25)
# 	planeWidget.SetInputData(mesh)
# 	planeWidget.PlaceWidget()
# 	planeWidget.AddObserver("InteractionEvent", myCallback)
   
#     #
#     #mwf debug see if can fix scaling problems
#     # Create the Legend for the plot
#     #mwf debug doesn't cause problem
#     data.legend = DisplayGenericLegend(ren, lut, title)
    
#     # Create the XYZ axis for the plot
#     DisplayAxis(ren, xMag)
   
#     # Put the data on the screen
#     RenderWindow(ren, renWin, Pause)
#     if created:
#         ren.ResetCamera()
#     if Hardcopy == True:
#         GenerateHardcopy(renWin, title, format= g.DefaultImageFormat, filepath = g.ImageFolderPath,
#                          is3D = True, isVectorPlot = False, forceWriter = 'WriteImage')

##################################################
#start adding some tools for particles

def viewScalar_pointCloud_2D(nodes, scalars, title, winNum, nodeArray = None, elementNodesArray = None, vectors = None,
                             Pause = True, Hardcopy = False,viewTypes=['spheres']):#,'contour','warp']):
    """
    Display group of points and a value given the x,y,z point locations passed in 
    
    Arguments:
        nodes:       Double[]   Array of doubles in X, Y, Z, X, Y, Z format for node locations
        scalars:     Double[].  Array of scalar values, one value for each node
        title:       String.  Title of the plot.
        winNum:      Integer.  Window number.  Used like the window number argument for gnuplot windows.  
                     If window number is the same as one that already exists, it will update the existing window.  
                     If not, a new window is created.
        Pause:       Boolean.  Specify whether to open VTK in interactive mode (if True) or just display mode.
    """
    import cvtkviewers
    global windowDict
    #windowName = mesh.domain.name+mesh.name+variableName
    windowName = "Quadrature Mesh"+title
    t = 0.0
    if not windowDict.has_key(windowName):
        windowCreated=True
        window = Window(windowName,title)
        windowDict[windowName] = window
        window.vod['scalars'] = cvtkviewers.prepareScalarValueArray(scalars)
        window.vod['points'] = cvtkviewers.prepareVTKPoints3(nodes)
        # Create a cloud with trivial connectivity
        window.vod['polyData'] = vtk.vtkPolyData()
        window.vod['polyData'].SetPoints(window.vod['points'])
        #for now use a stupid python loop
        window.vod['cells'] = vtk.vtkCellArray()
        npts = window.vod['points'].GetNumberOfPoints()
        window.vod['cells'].Allocate(vtk.VTK_VERTEX,5)
        for i in range(npts):
            window.vod['cells'].InsertNextCell(vtk.VTK_VERTEX)
            window.vod['cells'].InsertCellPoint(i)
        #add in the connectivity
        window.vod['polyData'].SetVerts(window.vod['cells'])
        window.vod['dataSet']  = window.vod['polyData'] #all that is needed for now
        #and the property
        window.vod['dataSet'].GetPointData().SetScalars(window.vod['scalars'])
        #mapper
        window.vod['dataSetMapper'] = vtk.vtkDataSetMapper()
        window.vod['dataSetMapper'].SetInputData(window.vod['dataSet'])
        window.vod['dataSetMapper'].SetScalarVisibility(0)
        #background mesh
        if nodeArray == None or elementNodesArray == None:
            window.vod['dataSet_background'] = None
        else:
            window.vod['dataSet_background'] = cvtkviewers.getUnstructuredGridFromMesh(nodeArray,
                                                                                       elementNodesArray)
            window.vod['dataSet_background'].Update()
            window.vod['dataSetMapper_background']=vtk.vtkDataSetMapper()
            window.vod['dataSetMapper_background'].SetInputData(window.vod['dataSet_background'])

        if vectors != None:
            window.w=numpy.zeros(vectors.shape[:-1],'d')
            window.vectors=numpy.column_stack((vectors.flat[::2],vectors.flat[1::2],window.w.flat)).flatten()
            window.vod['vectors'] = cvtkviewers.prepareVectorValueArray(window.vectors)
            window.vod['dataSet'].GetPointData().SetVectors(window.vod['vectors'])
#         window.vod['lut'] = vtk.vtkLookupTable()
#         window.vod['lut'].SetTableRange(window.vod['scalars'].GetRange())
#         window.vod['lut'].SetHueRange(0.66667,0.0)
#         window.vod['dataSetMapper'].SetLookupTable(window.vod['lut'])
#         window.vod['dataSetMapper'].UseLookupTableScalarRangeOn()
        
    else:
        windowCreated=False
        window=windowDict[windowName]
        #assume for now number of points stays the same
        #window.vod['points'] = cvtkviewers.prepareVTKPoints3(nodes)
        #window.vod['polyData'].SetPoints(window.vod['points'])
        window.vod['polyData'].Modified()
        window.vod['scalars'].Modified()
        
        window.vod['polyData'].Update()
        window.vod['glyph'].Update()
#         window.vod['lut'].SetTableRange(window.vod['scalars'].GetRange())
    viewParticles_2D(window,windowCreated,viewTypes)
def viewParticles_2D(window,
                  windowCreated,
                  viewTypes):
    """
    TODO:
          need scaling for spheres
    """
    global g,windowDict
    import pdb
    if windowCreated:
        createRenderers(viewTypes,window)
        if window.vod.has_key('dataSet_background') and window.vod['dataSet_background'] !=  None:
            window.vod['actor_background'] = vtk.vtkActor()
            window.vod['actor_background'].SetMapper(window.vod['dataSetMapper_background'])
            window.vod['actor_background'].GetProperty().SetRepresentationToWireframe()
        if 'spheres' in viewTypes:
            window.vod['glyph'] = vtk.vtkGlyph3D()
            window.vod['glyph'].SetInputData(window.vod['dataSet'])
            window.vod['glyph'].SetColorModeToColorByScalar()
            window.vod['glyph'].SetScaleModeToDataScalingOff()
            #window.vod['glyph'].SetScaleFactor(0.1)
            window.vod['sphere'] = vtk.vtkSphereSource()
            window.vod['sphere'].SetRadius(0.005)
            window.vod['glyph'].SetSource(window.vod['sphere'].GetOutput())
            window.vod['glyphMapper'] = vtk.vtkPolyDataMapper()
            window.vod['glyphMapper'].SetInput(window.vod['glyph'].GetOutput())
             #window.vod['glyphMapper'].SetLookupTable(window.vod['lut'])
             #window.vod['glyphMapper'].UseLookupTableScalarRangeOn()
            window.vod['glyphActor'] = vtk.vtkActor()
            window.vod['glyphActor'].SetMapper(window.vod['glyphMapper'])
            window.vod['glyphActor'].GetProperty().SetColor(0,1,0)
            ren = window.vod['ren_spheres']
            ren.AddActor(window.vod['glyphActor'])
            if window.vod.has_key('actor_background'):
                ren.AddActor(window.vod['actor_background'])
            ren.ResetCamera()
        if 'streamlines' in viewTypes and window.vod.has_key('vectors'):
            #
            #mwf debug
            #pdb.set_trace()
            streamer = vtk.vtkStreamTracer()
            window.vod['streamLineFilter'] = streamer
            streamer.SetInputData(window.vod['dataSet'])
            streamer.SetIntegratorTypeToRungeKutta45()
            streamer.SetIntegrationDirectionToForward()
            #need to pass in a time argument
            streamer.SetMaximumPropagation(0,1.0)
            #supposed to use SetSourceConnection?
            streamer.SetSource(window.vod['dataSet'])
            #
            window.vod['rk4'] = vtk.vtkRungeKutta4()
            streamer.SetIntegrator(window.vod['rk4'])
            #
            window.vod['streamLineMapper']=vtk.vtkPolyDataMapper()
            window.vod['streamLineMapper'].SetInputConnection(streamer.GetOutputPort())
            #
            window.vod['streamLineActor'] = vtk.vtkActor()
            streamline = window.vod['streamLineActor']
            window.vod['streamLineActor'].SetMapper(window.vod['streamLineMapper'])
            window.vod['streamLineActor'].VisibilityOn()
#             window.vod['gridActor_colorMapped'] = vtk.vtkActor()
#             window.vod['gridActor_colorMapped'].SetMapper(window.vod['dataSetMapper'])
#             window.vod['axesActor_colorMapped'] = vtk.vtkCubeAxesActor2D()
#             window.vod['axesActor_colorMapped'].SetInputData(window.vod['dataSet'])
#             window.vod['axesActor_colorMapped'].SetCamera(window.vod['ren_streamlines'].GetActiveCamera())
#             window.vod['axesActor_colorMapped'].SetFlyModeToClosestTriad()
#             window.vod['axesActor_colorMapped'].SetZAxisVisibility(0)
#             tp=window.vod['axesActor_colorMapped'].GetAxisTitleTextProperty()
#             tp.SetFontSize(14)
#             tp.SetFontFamilyToArial()
#             tp.SetColor(0,0,0)
#             tp=window.vod['axesActor_colorMapped'].GetAxisLabelTextProperty()
#             tp.SetFontSize(14)
#             tp.SetFontFamilyToArial()
#             tp.SetColor(0,0,0)
#             window.vod['legendActor_colorMapped']= createLegendActor(window.vod['lut'])
#             window.vod['legendActor_colorMapped'].SetLookupTable(window.vod['lut'])
            ren = window.vod['ren_streamlines']
#             ren.AddActor(window.vod['gridActor_colorMapped'])
#             ren.AddActor(window.vod['legendActor_colorMapped'])
#             ren.AddActor2D(window.vod['axesActor_colorMapped'])
            ren.AddActor(window.vod['streamLineActor'])
            ren.ResetCamera()
#         if 'colorMapped' in viewTypes:
#             window.vod['gridActor_colorMapped'] = vtk.vtkActor()
#             window.vod['gridActor_colorMapped'].SetMapper(window.vod['dataSetMapper'])
#             window.vod['axesActor_colorMapped'] = vtk.vtkCubeAxesActor2D()
#             window.vod['axesActor_colorMapped'].SetInputData(window.vod['dataSet'])
#             window.vod['axesActor_colorMapped'].SetCamera(window.vod['ren_colorMapped'].GetActiveCamera())
#             window.vod['axesActor_colorMapped'].SetFlyModeToClosestTriad()
#             window.vod['axesActor_colorMapped'].SetZAxisVisibility(0)
#             tp=window.vod['axesActor_colorMapped'].GetAxisTitleTextProperty()
#             tp.SetFontSize(14)
#             tp.SetFontFamilyToArial()
#             tp.SetColor(0,0,0)
#             tp=window.vod['axesActor_colorMapped'].GetAxisLabelTextProperty()
#             tp.SetFontSize(14)
#             tp.SetFontFamilyToArial()
#             tp.SetColor(0,0,0)
#             window.vod['legendActor_colorMapped']= createLegendActor(window.vod['lut'])
#             window.vod['legendActor_colorMapped'].SetLookupTable(window.vod['lut'])
#             ren = window.vod['ren_colorMapped']
#             ren.AddActor(window.vod['gridActor_colorMapped'])
#             ren.AddActor(window.vod['legendActor_colorMapped'])
#             ren.AddActor2D(window.vod['axesActor_colorMapped'])
#             ren.ResetCamera()
#         if 'contour' in viewTypes:
#             window.vod['algo_contour'] = vtk.vtkContourGrid()
#             window.vod['algo_contour'].SetInput(window.vod['normFilter'].GetOutput())
#             window.vod['algo_contour'].GenerateValues(10,window.vod['normFilter'].GetOutput().GetPointData().GetScalars().GetRange())            
#             window.vod['contourMapper'] = vtk.vtkPolyDataMapper()
#             window.vod['contourMapper'].SetInputConnection(window.vod['algo_contour'].GetOutputPort())
#             window.vod['contourMapper'].SetLookupTable(window.vod['lut'])
#             window.vod['contourMapper'].UseLookupTableScalarRangeOn()
#             window.vod['actorcontour'] = vtk.vtkActor()
#             window.vod['actorcontour'].SetMapper(window.vod['contourMapper'])
#             window.vod['axesActor_contour'] = vtk.vtkCubeAxesActor2D()
#             window.vod['axesActor_contour'].SetInputData(window.vod['dataSet'])
#             window.vod['axesActor_contour'].SetCamera(window.vod['ren_contour'].GetActiveCamera())
#             window.vod['axesActor_contour'].SetFlyModeToClosestTriad()
#             window.vod['axesActor_contour'].SetZAxisVisibility(0)
#             tp=window.vod['axesActor_contour'].GetAxisTitleTextProperty()
#             tp.SetFontSize(14)
#             tp.SetFontFamilyToArial()
#             tp.SetColor(0,0,0)
#             tp=window.vod['axesActor_contour'].GetAxisLabelTextProperty()
#             tp.SetFontSize(14)
#             tp.SetFontFamilyToArial()
#             tp.SetColor(0,0,0)
#             window.vod['legendActor_contour']= createLegendActor(window.vod['lut'])
#             window.vod['legendActor_contour'].SetLookupTable(window.vod['lut'])
#             ren = window.vod['ren_contour']
#             ren.AddActor(window.vod['actorcontour'])
#             ren.AddActor(window.vod['legendActor_contour'])
#             ren.AddActor2D(window.vod['axesActor_contour'])
#             ren.ResetCamera()
#         if 'warp' in viewTypes:
#             window.vod['algo_geo'] = vtk.vtkGeometryFilter()
#             window.vod['algo_geo'].SetInput(window.vod['normFilter'].GetOutput())
#             window.vod['algo_warp'] = vtk.vtkWarpScalar()
#             window.vod['algo_warp'].SetInputConnection(window.vod['algo_geo'].GetOutputPort())
#             window.vod['warpMapper'] = vtk.vtkPolyDataMapper()
#             window.vod['warpMapper'].SetInputConnection(window.vod['algo_warp'].GetOutputPort())
#             window.vod['warpMapper'].SetLookupTable(window.vod['lut'])
#             window.vod['warpMapper'].UseLookupTableScalarRangeOn()
#             window.vod['warpActor'] = vtk.vtkActor()
#             window.vod['warpActor'].SetMapper(window.vod['warpMapper'])
#             window.vod['axesActor_warp'] = vtk.vtkCubeAxesActor2D()
#             window.vod['axesActor_warp'].SetInput(window.vod['dataSet'])
#             window.vod['axesActor_warp'].SetCamera(window.vod['ren_warp'].GetActiveCamera())
#             window.vod['axesActor_warp'].SetFlyModeToClosestTriad()
#             window.vod['axesActor_warp'].SetZAxisVisibility(0)
#             tp=window.vod['axesActor_warp'].GetAxisTitleTextProperty()
#             tp.SetFontSize(14)
#             tp.SetFontFamilyToArial()
#             tp.SetColor(0,0,0)
#             tp=window.vod['axesActor_warp'].GetAxisLabelTextProperty()
#             tp.SetFontSize(14)
#             tp.SetFontFamilyToArial()
#             tp.SetColor(0,0,0)
#             window.vod['legendActor_warp']= createLegendActor(window.vod['lut'])
#             window.vod['legendActor_warp'].SetLookupTable(window.vod['lut'])
#             ren = window.vod['ren_warp']
#             ren.AddActor(window.vod['warpActor'])
#             ren.AddActor(window.vod['legendActor_warp'])
#             ren.AddActor2D(window.vod['axesActor_warp'])
#             ren.ResetCamera()
#     if 'contour' in viewTypes:
#         window.vod['algo_contour'].GenerateValues(10,window.vod['normFilter'].GetOutput().GetPointData().GetScalars().GetRange())
#     if 'streamlines' in viewTypes:
#         window.vod['lineWidget'].GetPolyData(window.vod['seeds'])
    #mwf debug
    #pdb.set_trace()
    window.renWin.Render()
    g.app.processEvents()

def viewScalar_pointCloud_3D(nodes, scalars, title, winNum, nodeArray = None, elementNodesArray = None, vectors = None,
                             Pause = True, Hardcopy = False,viewTypes=['spheres']):#,'contour','warp']):
    """
    Display group of points and a value given the x,y,z point locations passed in 
    
    Arguments:
        nodes:       Double[]   Array of doubles in X, Y, Z, X, Y, Z format for node locations
        scalars:     Double[].  Array of scalar values, one value for each node
        title:       String.  Title of the plot.
        winNum:      Integer.  Window number.  Used like the window number argument for gnuplot windows.  
                     If window number is the same as one that already exists, it will update the existing window.  
                     If not, a new window is created.
        Pause:       Boolean.  Specify whether to open VTK in interactive mode (if True) or just display mode.
    """
    import cvtkviewers
    global windowDict
    #windowName = mesh.domain.name+mesh.name+variableName
    windowName = "Quadrature Mesh"+title
    t = 0.0
    if not windowDict.has_key(windowName):
        windowCreated=True
        window = Window(windowName,title)
        windowDict[windowName] = window
        window.vod['scalars'] = cvtkviewers.prepareScalarValueArray(scalars)
        window.vod['points'] = cvtkviewers.prepareVTKPoints3(nodes)
        # Create a cloud with trivial connectivity
        window.vod['polyData'] = vtk.vtkPolyData()
        window.vod['polyData'].SetPoints(window.vod['points'])
        #for now use a stupid python loop
        window.vod['cells'] = vtk.vtkCellArray()
        npts = window.vod['points'].GetNumberOfPoints()
        window.vod['cells'].Allocate(vtk.VTK_VERTEX,5)
        for i in range(npts):
            window.vod['cells'].InsertNextCell(vtk.VTK_VERTEX)
            window.vod['cells'].InsertCellPoint(i)
        #add in the connectivity
        window.vod['polyData'].SetVerts(window.vod['cells'])
        window.vod['dataSet']  = window.vod['polyData'] #all that is needed for now
        #and the property
        window.vod['dataSet'].GetPointData().SetScalars(window.vod['scalars'])
        #mapper
        window.vod['dataSetMapper'] = vtk.vtkDataSetMapper()
        window.vod['dataSetMapper'].SetInputData(window.vod['dataSet'])
        window.vod['dataSetMapper'].SetScalarVisibility(0)
        #background mesh
        if nodeArray == None or elementNodesArray == None:
            window.vod['dataSet_background'] = None
        else:
            window.vod['dataSet_background'] = cvtkviewers.getUnstructuredGridFromMesh(nodeArray,
                                                                                       elementNodesArray)
            window.vod['dataSet_background'].Update()
            window.vod['dataSetMapper_background']=vtk.vtkDataSetMapper()
            window.vod['dataSetMapper_background'].SetInputData(window.vod['dataSet_background'])

        if vectors != None:
            window.vectors=numpy.column_stack((vectors.flat[::3],vectors.flat[1::3],vectors.flat[2::3])).flatten()
            window.vod['vectors'] = cvtkviewers.prepareVectorValueArray(window.vectors)
            window.vod['dataSet'].GetPointData().SetVectors(window.vod['vectors'])
#         window.vod['lut'] = vtk.vtkLookupTable()
#         window.vod['lut'].SetTableRange(window.vod['scalars'].GetRange())
#         window.vod['lut'].SetHueRange(0.66667,0.0)
#         window.vod['dataSetMapper'].SetLookupTable(window.vod['lut'])
#         window.vod['dataSetMapper'].UseLookupTableScalarRangeOn()
        
    else:
        windowCreated=False
        window=windowDict[windowName]
        #assume for now number of points stays the same
        #window.vod['points'] = cvtkviewers.prepareVTKPoints3(nodes)
        #window.vod['polyData'].SetPoints(window.vod['points'])
        window.vod['polyData'].Modified()
        window.vod['scalars'].Modified()
        
        window.vod['polyData'].Update()
        window.vod['glyph'].Update()
#         window.vod['lut'].SetTableRange(window.vod['scalars'].GetRange())
    viewParticles_3D(window,windowCreated,viewTypes)

def viewParticles_3D(window,
                     windowCreated,
                     viewTypes):
    """
    TODO:
          need scaling for spheres
    """
    global g,windowDict
    import pdb
    if windowCreated:
        createRenderers(viewTypes,window)
        if window.vod.has_key('dataSet_background') and window.vod['dataSet_background'] !=  None:
            window.vod['actor_background'] = vtk.vtkActor()
            window.vod['actor_background'].SetMapper(window.vod['dataSetMapper_background'])
            window.vod['actor_background'].GetProperty().SetRepresentationToWireframe()
        if 'spheres' in viewTypes or True:
            window.vod['glyph'] = vtk.vtkGlyph3D()
            window.vod['glyph'].SetInputData(window.vod['dataSet'])
            window.vod['glyph'].SetColorModeToColorByScalar()
            window.vod['glyph'].SetScaleModeToDataScalingOff()
            #window.vod['glyph'].SetScaleFactor(0.1)
            window.vod['sphere'] = vtk.vtkSphereSource()
            window.vod['sphere'].SetRadius(0.005)
            window.vod['glyph'].SetSource(window.vod['sphere'].GetOutput())
            window.vod['glyphMapper'] = vtk.vtkPolyDataMapper()
            window.vod['glyphMapper'].SetInput(window.vod['glyph'].GetOutput())
             #window.vod['glyphMapper'].SetLookupTable(window.vod['lut'])
             #window.vod['glyphMapper'].UseLookupTableScalarRangeOn()
            window.vod['glyphActor'] = vtk.vtkActor()
            window.vod['glyphActor'].SetMapper(window.vod['glyphMapper'])
            window.vod['glyphActor'].GetProperty().SetColor(0,1,0)
            ren = window.vod['ren_spheres']
            ren.AddActor(window.vod['glyphActor'])
            if window.vod.has_key('actor_background'):
                ren.AddActor(window.vod['actor_background'])
            ren.ResetCamera()
        if 'streamlines' in viewTypes and window.vod.has_key('vectors'):
            #
            streamer = vtk.vtkStreamTracer()
            window.vod['streamLineFilter'] = streamer
            streamer.SetInputData(window.vod['dataSet'])
            streamer.SetIntegratorTypeToRungeKutta45()
            streamer.SetIntegrationDirectionToForward()
            #supposed to use SetSourceConnection?
            streamer.SetSource(window.vod['dataSet'])
            #
            window.vod['rk4'] = vtk.vtkRungeKutta4()
            streamer.SetIntegrator(window.vod['rk4'])
            #
            window.vod['streamLineMapper']=vtk.vtkPolyDataMapper()
            window.vod['streamLineMapper'].SetInputConnection(streamer.GetOutputPort())
            #
            window.vod['streamLineActor'] = vtk.vtkActor()
            streamline = window.vod['streamLineActor']
            window.vod['streamLineActor'].SetMapper(window.vod['streamLineMapper'])
            window.vod['streamLineActor'].VisibilityOn()
#             window.vod['gridActor_colorMapped'] = vtk.vtkActor()
#             window.vod['gridActor_colorMapped'].SetMapper(window.vod['dataSetMapper'])
#             window.vod['axesActor_colorMapped'] = vtk.vtkCubeAxesActor2D()
#             window.vod['axesActor_colorMapped'].SetInputData(window.vod['dataSet'])
#             window.vod['axesActor_colorMapped'].SetCamera(window.vod['ren_streamlines'].GetActiveCamera())
#             window.vod['axesActor_colorMapped'].SetFlyModeToClosestTriad()
#             window.vod['axesActor_colorMapped'].SetZAxisVisibility(0)
#             tp=window.vod['axesActor_colorMapped'].GetAxisTitleTextProperty()
#             tp.SetFontSize(14)
#             tp.SetFontFamilyToArial()
#             tp.SetColor(0,0,0)
#             tp=window.vod['axesActor_colorMapped'].GetAxisLabelTextProperty()
#             tp.SetFontSize(14)
#             tp.SetFontFamilyToArial()
#             tp.SetColor(0,0,0)
#             window.vod['legendActor_colorMapped']= createLegendActor(window.vod['lut'])
#             window.vod['legendActor_colorMapped'].SetLookupTable(window.vod['lut'])
            ren = window.vod['ren_streamlines']
#             ren.AddActor(window.vod['gridActor_colorMapped'])
#             ren.AddActor(window.vod['legendActor_colorMapped'])
#             ren.AddActor2D(window.vod['axesActor_colorMapped'])
            ren.AddActor(window.vod['streamLineActor'])
            ren.ResetCamera()
    #mwf debug
    #pdb.set_trace()
    window.renWin.Render()
    g.app.processEvents()

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XDMF Reader'
periodic_solitary_wavesxmf = XDMFReader(FileNames=['/home/cekees/proteus/proteus/tests/periodic/hp025-new-weight/periodic_solitary_waves.xmf'])
periodic_solitary_wavesxmf.PointArrayStatus = ['nodeMaterialTypes', 'p', 'phi', 'phiCorr', 'phi_analytical', 'phid', 'u', 'v', 'velocity', 'vof', 'vof_analytical']
periodic_solitary_wavesxmf.CellArrayStatus = ['elementMaterialTypes']

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# Properties modified on periodic_solitary_wavesxmf
periodic_solitary_wavesxmf.GridStatus = ['Grid_2', 'Grid_5', 'Grid_8', 'Grid_11', 'Grid_14', 'Grid_17', 'Grid_20', 'Grid_23', 'Grid_26', 'Grid_29', 'Grid_32', 'Grid_35', 'Grid_38', 'Grid_41', 'Grid_44', 'Grid_47', 'Grid_50', 'Grid_53', 'Grid_56', 'Grid_59', 'Grid_62', 'Grid_65', 'Grid_68', 'Grid_71', 'Grid_74', 'Grid_77', 'Grid_80', 'Grid_83', 'Grid_86', 'Grid_89', 'Grid_92', 'Grid_95', 'Grid_98', 'Grid_101', 'Grid_104', 'Grid_107', 'Grid_110', 'Grid_113', 'Grid_116', 'Grid_119', 'Grid_122', 'Grid_125', 'Grid_128', 'Grid_131', 'Grid_134', 'Grid_137', 'Grid_140', 'Grid_143', 'Grid_146', 'Grid_149', 'Grid_152', 'Grid_155', 'Grid_158', 'Grid_161', 'Grid_164', 'Grid_167', 'Grid_170', 'Grid_173', 'Grid_176', 'Grid_179', 'Grid_182', 'Grid_185']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [2113, 1315]

# get color transfer function/color map for 'nodeMaterialTypes'
nodeMaterialTypesLUT = GetColorTransferFunction('nodeMaterialTypes')
nodeMaterialTypesLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 2.0, 0.865003, 0.865003, 0.865003, 4.0, 0.705882, 0.0156863, 0.14902]
nodeMaterialTypesLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'nodeMaterialTypes'
nodeMaterialTypesPWF = GetOpacityTransferFunction('nodeMaterialTypes')
nodeMaterialTypesPWF.Points = [0.0, 0.0, 0.5, 0.0, 4.0, 1.0, 0.5, 0.0]
nodeMaterialTypesPWF.ScalarRangeInitialized = 1

# show data in view
periodic_solitary_wavesxmfDisplay = Show(periodic_solitary_wavesxmf, renderView1)
# trace defaults for the display properties.
periodic_solitary_wavesxmfDisplay.Representation = 'Surface'
periodic_solitary_wavesxmfDisplay.ColorArrayName = ['POINTS', 'nodeMaterialTypes']
periodic_solitary_wavesxmfDisplay.LookupTable = nodeMaterialTypesLUT
periodic_solitary_wavesxmfDisplay.OSPRayScaleArray = 'nodeMaterialTypes'
periodic_solitary_wavesxmfDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
periodic_solitary_wavesxmfDisplay.SelectOrientationVectors = 'velocity'
periodic_solitary_wavesxmfDisplay.ScaleFactor = 1.5
periodic_solitary_wavesxmfDisplay.SelectScaleArray = 'nodeMaterialTypes'
periodic_solitary_wavesxmfDisplay.GlyphType = 'Arrow'
periodic_solitary_wavesxmfDisplay.GlyphTableIndexArray = 'nodeMaterialTypes'
periodic_solitary_wavesxmfDisplay.DataAxesGrid = 'GridAxesRepresentation'
periodic_solitary_wavesxmfDisplay.PolarAxes = 'PolarAxesRepresentation'
periodic_solitary_wavesxmfDisplay.ScalarOpacityFunction = nodeMaterialTypesPWF
periodic_solitary_wavesxmfDisplay.ScalarOpacityUnitDistance = 0.41837568517341406
periodic_solitary_wavesxmfDisplay.GaussianRadius = 0.75
periodic_solitary_wavesxmfDisplay.SetScaleArray = ['POINTS', 'nodeMaterialTypes']
periodic_solitary_wavesxmfDisplay.ScaleTransferFunction = 'PiecewiseFunction'
periodic_solitary_wavesxmfDisplay.OpacityArray = ['POINTS', 'nodeMaterialTypes']
periodic_solitary_wavesxmfDisplay.OpacityTransferFunction = 'PiecewiseFunction'
#periodic_solitary_wavesxmfDisplay.SelectInputVectors = ['POINTS', 'velocity']
#periodic_solitary_wavesxmfDisplay.WriteLog = ''

# reset view to fit data
renderView1.ResetCamera()

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [7.5, 1.0, 10000.0]
renderView1.CameraFocalPoint = [7.5, 1.0, 0.0]

# show color bar/color legend
periodic_solitary_wavesxmfDisplay.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Contour'
contour1 = Contour(Input=periodic_solitary_wavesxmf)
contour1.ContourBy = ['POINTS', 'nodeMaterialTypes']
contour1.Isosurfaces = [2.0]
contour1.PointMergeMethod = 'Uniform Binning'

# Properties modified on contour1
contour1.ContourBy = ['POINTS', 'phi']
contour1.Isosurfaces = [0.0]

# get color transfer function/color map for 'elementMaterialTypes'
elementMaterialTypesLUT = GetColorTransferFunction('elementMaterialTypes')
elementMaterialTypesLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 5.878906683738906e-39, 0.865003, 0.865003, 0.865003, 1.1757813367477812e-38, 0.705882, 0.0156863, 0.14902]
elementMaterialTypesLUT.ScalarRangeInitialized = 1.0

# show data in view
contour1Display = Show(contour1, renderView1)
# trace defaults for the display properties.
contour1Display.Representation = 'Surface'
contour1Display.ColorArrayName = ['CELLS', 'elementMaterialTypes']
contour1Display.LookupTable = elementMaterialTypesLUT
contour1Display.OSPRayScaleArray = 'elementMaterialTypes'
contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour1Display.SelectOrientationVectors = 'velocity'
contour1Display.ScaleFactor = 1.5
contour1Display.SelectScaleArray = 'elementMaterialTypes'
contour1Display.GlyphType = 'Arrow'
contour1Display.GlyphTableIndexArray = 'elementMaterialTypes'
contour1Display.DataAxesGrid = 'GridAxesRepresentation'
contour1Display.PolarAxes = 'PolarAxesRepresentation'
contour1Display.GaussianRadius = 0.75
contour1Display.SetScaleArray = ['POINTS', 'nodeMaterialTypes']
contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
contour1Display.OpacityArray = ['POINTS', 'nodeMaterialTypes']
contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
#contour1Display.SelectInputVectors = ['POINTS', 'velocity']
#contour1Display.WriteLog = ''

# hide data in view
Hide(periodic_solitary_wavesxmf, renderView1)

# show color bar/color legend
contour1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(elementMaterialTypesLUT, renderView1)

# Properties modified on contour1Display
contour1Display.LineWidth = 4.0

# create a new 'Calculator'
calculator1 = Calculator(Input=contour1)
calculator1.Function = ''

# Properties modified on calculator1
calculator1.ResultArrayName = 'wave_profile'
calculator1.Function = 'coordsY'

# get color transfer function/color map for 'wave_profile'
wave_profileLUT = GetColorTransferFunction('wave_profile')
wave_profileLUT.RGBPoints = [1.000295525832002, 0.231373, 0.298039, 0.752941, 1.2251358584516647, 0.865003, 0.865003, 0.865003, 1.4499761910713274, 0.705882, 0.0156863, 0.14902]
wave_profileLUT.ScalarRangeInitialized = 1.0

# show data in view
calculator1Display = Show(calculator1, renderView1)
# trace defaults for the display properties.
calculator1Display.Representation = 'Surface'
calculator1Display.ColorArrayName = ['POINTS', 'wave_profile']
calculator1Display.LookupTable = wave_profileLUT
calculator1Display.OSPRayScaleArray = 'wave_profile'
calculator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
calculator1Display.SelectOrientationVectors = 'velocity'
calculator1Display.ScaleFactor = 1.5
calculator1Display.SelectScaleArray = 'wave_profile'
calculator1Display.GlyphType = 'Arrow'
calculator1Display.GlyphTableIndexArray = 'wave_profile'
calculator1Display.DataAxesGrid = 'GridAxesRepresentation'
calculator1Display.PolarAxes = 'PolarAxesRepresentation'
calculator1Display.GaussianRadius = 0.75
calculator1Display.SetScaleArray = ['POINTS', 'wave_profile']
calculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
calculator1Display.OpacityArray = ['POINTS', 'wave_profile']
calculator1Display.OpacityTransferFunction = 'PiecewiseFunction'
#calculator1Display.SelectInputVectors = ['POINTS', 'velocity']
#calculator1Display.WriteLog = ''

# hide data in view
Hide(contour1, renderView1)

# show color bar/color legend
calculator1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Python Calculator'
pythonCalculator1 = PythonCalculator(Input=calculator1)
pythonCalculator1.Expression = ''

# Properties modified on pythonCalculator1
pythonCalculator1.Expression = 'max(wave_profile)'
pythonCalculator1.ArrayName = 'max_wave_profile'

# show data in view
pythonCalculator1Display = Show(pythonCalculator1, renderView1)
# trace defaults for the display properties.
pythonCalculator1Display.Representation = 'Surface'
pythonCalculator1Display.ColorArrayName = ['POINTS', 'wave_profile']
pythonCalculator1Display.LookupTable = wave_profileLUT
pythonCalculator1Display.OSPRayScaleArray = 'wave_profile'
pythonCalculator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
pythonCalculator1Display.SelectOrientationVectors = 'velocity'
pythonCalculator1Display.ScaleFactor = 1.5
pythonCalculator1Display.SelectScaleArray = 'wave_profile'
pythonCalculator1Display.GlyphType = 'Arrow'
pythonCalculator1Display.GlyphTableIndexArray = 'wave_profile'
pythonCalculator1Display.DataAxesGrid = 'GridAxesRepresentation'
pythonCalculator1Display.PolarAxes = 'PolarAxesRepresentation'
pythonCalculator1Display.GaussianRadius = 0.75
pythonCalculator1Display.SetScaleArray = ['POINTS', 'wave_profile']
pythonCalculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
pythonCalculator1Display.OpacityArray = ['POINTS', 'wave_profile']
pythonCalculator1Display.OpacityTransferFunction = 'PiecewiseFunction'
#pythonCalculator1Display.SelectInputVectors = ['POINTS', 'velocity']
#pythonCalculator1Display.WriteLog = ''

# hide data in view
Hide(calculator1, renderView1)

# show color bar/color legend
pythonCalculator1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(periodic_solitary_wavesxmf)

# create a new 'Contour'
contour2 = Contour(Input=periodic_solitary_wavesxmf)
contour2.ContourBy = ['POINTS', 'nodeMaterialTypes']
contour2.Isosurfaces = [2.0]
contour2.PointMergeMethod = 'Uniform Binning'

# Properties modified on contour2
contour2.ContourBy = ['POINTS', 'phi_analytical']
contour2.Isosurfaces = [0.0]

# show data in view
contour2Display = Show(contour2, renderView1)
# trace defaults for the display properties.
contour2Display.Representation = 'Surface'
contour2Display.ColorArrayName = ['CELLS', 'elementMaterialTypes']
contour2Display.LookupTable = elementMaterialTypesLUT
contour2Display.OSPRayScaleArray = 'elementMaterialTypes'
contour2Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour2Display.SelectOrientationVectors = 'velocity'
contour2Display.ScaleFactor = 1.5
contour2Display.SelectScaleArray = 'elementMaterialTypes'
contour2Display.GlyphType = 'Arrow'
contour2Display.GlyphTableIndexArray = 'elementMaterialTypes'
contour2Display.DataAxesGrid = 'GridAxesRepresentation'
contour2Display.PolarAxes = 'PolarAxesRepresentation'
contour2Display.GaussianRadius = 0.75
contour2Display.SetScaleArray = ['POINTS', 'nodeMaterialTypes']
contour2Display.ScaleTransferFunction = 'PiecewiseFunction'
contour2Display.OpacityArray = ['POINTS', 'nodeMaterialTypes']
contour2Display.OpacityTransferFunction = 'PiecewiseFunction'
#contour2Display.SelectInputVectors = ['POINTS', 'velocity']
#contour2Display.WriteLog = ''

# hide data in view
Hide(periodic_solitary_wavesxmf, renderView1)

# show color bar/color legend
contour2Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(elementMaterialTypesLUT, renderView1)

# change solid color
contour2Display.DiffuseColor = [0.0, 0.0, 0.0]

# Properties modified on contour2Display
contour2Display.LineWidth = 4.0

# create a new 'Calculator'
calculator2 = Calculator(Input=contour2)
calculator2.Function = ''

# Properties modified on calculator2
calculator2.ResultArrayName = 'wave_profile_analytical'
calculator2.Function = 'coordsY'

# get color transfer function/color map for 'wave_profile_analytical'
wave_profile_analyticalLUT = GetColorTransferFunction('wave_profile_analytical')
wave_profile_analyticalLUT.RGBPoints = [1.000295525832002, 0.231373, 0.298039, 0.752941, 1.2251358584516647, 0.865003, 0.865003, 0.865003, 1.4499761910713274, 0.705882, 0.0156863, 0.14902]
wave_profile_analyticalLUT.ScalarRangeInitialized = 1.0

# show data in view
calculator2Display = Show(calculator2, renderView1)
# trace defaults for the display properties.
calculator2Display.Representation = 'Surface'
calculator2Display.ColorArrayName = ['POINTS', 'wave_profile_analytical']
calculator2Display.LookupTable = wave_profile_analyticalLUT
calculator2Display.OSPRayScaleArray = 'wave_profile_analytical'
calculator2Display.OSPRayScaleFunction = 'PiecewiseFunction'
calculator2Display.SelectOrientationVectors = 'velocity'
calculator2Display.ScaleFactor = 1.5
calculator2Display.SelectScaleArray = 'wave_profile_analytical'
calculator2Display.GlyphType = 'Arrow'
calculator2Display.GlyphTableIndexArray = 'wave_profile_analytical'
calculator2Display.DataAxesGrid = 'GridAxesRepresentation'
calculator2Display.PolarAxes = 'PolarAxesRepresentation'
calculator2Display.GaussianRadius = 0.75
calculator2Display.SetScaleArray = ['POINTS', 'wave_profile_analytical']
calculator2Display.ScaleTransferFunction = 'PiecewiseFunction'
calculator2Display.OpacityArray = ['POINTS', 'wave_profile_analytical']
calculator2Display.OpacityTransferFunction = 'PiecewiseFunction'
#calculator2Display.SelectInputVectors = ['POINTS', 'velocity']
#calculator2Display.WriteLog = ''

# hide data in view
Hide(contour2, renderView1)

# show color bar/color legend
calculator2Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Python Calculator'
pythonCalculator2 = PythonCalculator(Input=calculator2)
pythonCalculator2.Expression = ''

# Properties modified on pythonCalculator2
pythonCalculator2.Expression = 'max(wave_profile_analytical)'
pythonCalculator2.ArrayName = 'max(wave_profile_analytical)'

# show data in view
pythonCalculator2Display = Show(pythonCalculator2, renderView1)
# trace defaults for the display properties.
pythonCalculator2Display.Representation = 'Surface'
pythonCalculator2Display.ColorArrayName = ['POINTS', 'wave_profile_analytical']
pythonCalculator2Display.LookupTable = wave_profile_analyticalLUT
pythonCalculator2Display.OSPRayScaleArray = 'wave_profile_analytical'
pythonCalculator2Display.OSPRayScaleFunction = 'PiecewiseFunction'
pythonCalculator2Display.SelectOrientationVectors = 'velocity'
pythonCalculator2Display.ScaleFactor = 1.5
pythonCalculator2Display.SelectScaleArray = 'wave_profile_analytical'
pythonCalculator2Display.GlyphType = 'Arrow'
pythonCalculator2Display.GlyphTableIndexArray = 'wave_profile_analytical'
pythonCalculator2Display.DataAxesGrid = 'GridAxesRepresentation'
pythonCalculator2Display.PolarAxes = 'PolarAxesRepresentation'
pythonCalculator2Display.GaussianRadius = 0.75
pythonCalculator2Display.SetScaleArray = ['POINTS', 'wave_profile_analytical']
pythonCalculator2Display.ScaleTransferFunction = 'PiecewiseFunction'
pythonCalculator2Display.OpacityArray = ['POINTS', 'wave_profile_analytical']
pythonCalculator2Display.OpacityTransferFunction = 'PiecewiseFunction'
#pythonCalculator2Display.SelectInputVectors = ['POINTS', 'velocity']
#pythonCalculator2Display.WriteLog = ''

# hide data in view
Hide(calculator2, renderView1)

# show color bar/color legend
pythonCalculator2Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on pythonCalculator2
pythonCalculator2.ArrayName = 'max_wave_profile_analytical'

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(calculator2)

# Properties modified on calculator2
calculator2.ResultArrayName = 'wave_profile_analytical'

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on calculator2Display
calculator2Display.SelectScaleArray = 'None'

# Properties modified on calculator2Display
calculator2Display.GlyphTableIndexArray = 'None'

# Properties modified on calculator2Display
calculator2Display.SetScaleArray = ['POINTS', 'nodeMaterialTypes']

# Properties modified on calculator2Display
calculator2Display.OpacityArray = ['POINTS', 'nodeMaterialTypes']

# Properties modified on calculator2Display
calculator2Display.OSPRayScaleArray = 'nodeMaterialTypes'

# set active source
SetActiveSource(pythonCalculator2)

# set active source
SetActiveSource(contour2)

# set active source
SetActiveSource(contour2)

# show data in view
contour2Display = Show(contour2, renderView1)

# set active source
SetActiveSource(contour1)

# show data in view
contour1Display = Show(contour1, renderView1)

# hide data in view
Hide(pythonCalculator2, renderView1)

# hide data in view
Hide(pythonCalculator1, renderView1)

# set active source
SetActiveSource(periodic_solitary_wavesxmf)

# show data in view
periodic_solitary_wavesxmfDisplay = Show(periodic_solitary_wavesxmf, renderView1)

# show color bar/color legend
periodic_solitary_wavesxmfDisplay.SetScalarBarVisibility(renderView1, True)

# set scalar coloring
ColorBy(periodic_solitary_wavesxmfDisplay, ('POINTS', 'velocity', 'Magnitude'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(nodeMaterialTypesLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
periodic_solitary_wavesxmfDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
periodic_solitary_wavesxmfDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'velocity'
velocityLUT = GetColorTransferFunction('velocity')
velocityLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 0.918608115132005, 0.865003, 0.865003, 0.865003, 1.83721623026401, 0.705882, 0.0156863, 0.14902]
velocityLUT.ScalarRangeInitialized = 1.0

# change representation type
periodic_solitary_wavesxmfDisplay.SetRepresentationType('Surface LIC')

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
velocityLUT.ApplyPreset('Rainbow Desaturated', True)

# get color legend/bar for velocityLUT in view renderView1
velocityLUTColorBar = GetScalarBar(velocityLUT, renderView1)
velocityLUTColorBar.Title = 'velocity'
velocityLUTColorBar.ComponentTitle = 'Magnitude'

# change scalar bar placement
velocityLUTColorBar.Orientation = 'Horizontal'
velocityLUTColorBar.WindowLocation = 'AnyLocation'
velocityLUTColorBar.Position = [0.5357122574538571, 0.6381749049429658]
velocityLUTColorBar.ScalarBarLength = 0.3300000000000003

# create a new 'Annotate Time'
annotateTime1 = AnnotateTime()

# Properties modified on annotateTime1
annotateTime1.Format = 'Time: %2.2f s'

# show data in view
annotateTime1Display = Show(annotateTime1, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0

# change scalar bar placement
velocityLUTColorBar.Position = [0.4758343016948945, 0.8404562737642586]

# set active source
SetActiveSource(periodic_solitary_wavesxmf)

# Properties modified on velocityLUTColorBar
velocityLUTColorBar.TitleFontSize = 22
velocityLUTColorBar.LabelFontSize = 22

# Properties modified on velocityLUTColorBar
velocityLUTColorBar.TitleFontSize = 25
velocityLUTColorBar.LabelFontSize = 25

# Properties modified on velocityLUTColorBar
velocityLUTColorBar.TitleFontSize = 30
velocityLUTColorBar.LabelFontSize = 30

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [7.5, 1.0, 10000.0]
renderView1.CameraFocalPoint = [7.5, 1.0, 0.0]
renderView1.CameraParallelScale = 4.271020289569919

# save animation
#SaveAnimation('/home/cekees/proteus/proteus/tests/periodic/hp025-new-weight/soliton.ogv', renderView1, ImageResolution=[2622, 1315],
#    FrameWindow=[0, 61])

time=[]
max_numerical=[]
max_analytical=[]
wave_ele=[]
for i in range(61):
    # create a new 'Calculator'
    calculator1 = Calculator(Input=contour1)
    calculator1.Function = ''

    # Properties modified on calculator1
    calculator1.ResultArrayName = 'wave_profile'
    calculator1.Function = 'coordsY'

    # get color transfer function/color map for 'wave_profile'
    wave_profileLUT = GetColorTransferFunction('wave_profile')
    wave_profileLUT.RGBPoints = [1.000295525832002, 0.231373, 0.298039, 0.752941, 1.2251358584516647, 0.865003, 0.865003, 0.865003, 1.4499761910713274, 0.705882, 0.0156863, 0.14902]
    wave_profileLUT.ScalarRangeInitialized = 1.0

    # show data in view
    calculator1Display = Show(calculator1, renderView1)
    # trace defaults for the display properties.
    calculator1Display.Representation = 'Surface'
    calculator1Display.ColorArrayName = ['POINTS', 'wave_profile']
    calculator1Display.LookupTable = wave_profileLUT
    calculator1Display.OSPRayScaleArray = 'wave_profile'
    calculator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    calculator1Display.SelectOrientationVectors = 'velocity'
    calculator1Display.ScaleFactor = 1.5
    calculator1Display.SelectScaleArray = 'wave_profile'
    calculator1Display.GlyphType = 'Arrow'
    calculator1Display.GlyphTableIndexArray = 'wave_profile'
    calculator1Display.DataAxesGrid = 'GridAxesRepresentation'
    calculator1Display.PolarAxes = 'PolarAxesRepresentation'
    calculator1Display.GaussianRadius = 0.75
    calculator1Display.SetScaleArray = ['POINTS', 'wave_profile']
    calculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
    calculator1Display.OpacityArray = ['POINTS', 'wave_profile']
    calculator1Display.OpacityTransferFunction = 'PiecewiseFunction'
    #calculator1Display.SelectInputVectors = ['POINTS', 'velocity']
    #calculator1Display.WriteLog = ''

    # hide data in view
    Hide(contour1, renderView1)

    # show color bar/color legend
    calculator1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # set active source
    SetActiveSource(calculator1)

    # create a new 'Python Calculator'
    pythonCalculator1 = PythonCalculator(Input=calculator1)
    pythonCalculator1.Expression = ''

    # Properties modified on pythonCalculator1
    pythonCalculator1.Expression = 'max(wave_profile)'
    pythonCalculator1.ArrayName = 'max_wave_profile'

    # show data in view
    pythonCalculator1Display = Show(pythonCalculator1, renderView1)
    # trace defaults for the display properties.
    pythonCalculator1Display.Representation = 'Surface'
    pythonCalculator1Display.ColorArrayName = ['POINTS', 'wave_profile']
    pythonCalculator1Display.LookupTable = wave_profileLUT
    pythonCalculator1Display.OSPRayScaleArray = 'wave_profile'
    pythonCalculator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    pythonCalculator1Display.SelectOrientationVectors = 'velocity'
    pythonCalculator1Display.ScaleFactor = 1.5
    pythonCalculator1Display.SelectScaleArray = 'wave_profile'
    pythonCalculator1Display.GlyphType = 'Arrow'
    pythonCalculator1Display.GlyphTableIndexArray = 'wave_profile'
    pythonCalculator1Display.DataAxesGrid = 'GridAxesRepresentation'
    pythonCalculator1Display.PolarAxes = 'PolarAxesRepresentation'
    pythonCalculator1Display.GaussianRadius = 0.75
    pythonCalculator1Display.SetScaleArray = ['POINTS', 'wave_profile']
    pythonCalculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
    pythonCalculator1Display.OpacityArray = ['POINTS', 'wave_profile']
    pythonCalculator1Display.OpacityTransferFunction = 'PiecewiseFunction'
    #pythonCalculator1Display.SelectInputVectors = ['POINTS', 'velocity']
    #pythonCalculator1Display.WriteLog = ''

    # hide data in view
    Hide(calculator1, renderView1)

    # show color bar/color legend
    pythonCalculator1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # set active source
    SetActiveSource(periodic_solitary_wavesxmf)

    # create a new 'Calculator'
    calculator2 = Calculator(Input=calculator1)
    calculator2.Function = ''

    # set active source
    SetActiveSource(calculator1)

    # create a new 'Python Calculator'
    pythonCalculator1 = PythonCalculator(Input=calculator1)
    pythonCalculator1.Expression = ''

    # Properties modified on pythonCalculator1
    pythonCalculator1.Expression = 'max(wave_profile)'
    pythonCalculator1.ArrayName = 'max_wave_profile'

    # show data in view
    pythonCalculator1Display = Show(pythonCalculator1, renderView1)
    # trace defaults for the display properties.
    pythonCalculator1Display.Representation = 'Surface'
    pythonCalculator1Display.ColorArrayName = ['POINTS', 'wave_profile']
    pythonCalculator1Display.LookupTable = wave_profileLUT
    pythonCalculator1Display.OSPRayScaleArray = 'wave_profile'
    pythonCalculator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    pythonCalculator1Display.SelectOrientationVectors = 'velocity'
    pythonCalculator1Display.ScaleFactor = 1.5
    pythonCalculator1Display.SelectScaleArray = 'wave_profile'
    pythonCalculator1Display.GlyphType = 'Arrow'
    pythonCalculator1Display.GlyphTableIndexArray = 'wave_profile'
    pythonCalculator1Display.DataAxesGrid = 'GridAxesRepresentation'
    pythonCalculator1Display.PolarAxes = 'PolarAxesRepresentation'
    pythonCalculator1Display.GaussianRadius = 0.75
    pythonCalculator1Display.SetScaleArray = ['POINTS', 'wave_profile']
    pythonCalculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
    pythonCalculator1Display.OpacityArray = ['POINTS', 'wave_profile']
    pythonCalculator1Display.OpacityTransferFunction = 'PiecewiseFunction'
    #pythonCalculator1Display.SelectInputVectors = ['POINTS', 'velocity']
    #pythonCalculator1Display.WriteLog = ''

    # hide data in view
    Hide(calculator1, renderView1)

    # show color bar/color legend
    pythonCalculator1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # set active source
    SetActiveSource(periodic_solitary_wavesxmf)

    # create a new 'Calculator'
    calculator2 = Calculator(Input=contour2)
    calculator2.Function = ''

    # Properties modified on calculator2
    calculator2.ResultArrayName = 'wave_profile_analytical'
    calculator2.Function = 'coordsY'

    # get color transfer function/color map for 'wave_profile_analytical'
    wave_profile_analyticalLUT = GetColorTransferFunction('wave_profile_analytical')
    wave_profile_analyticalLUT.RGBPoints = [1.000295525832002, 0.231373, 0.298039, 0.752941, 1.2251358584516647, 0.865003, 0.865003, 0.865003, 1.4499761910713274, 0.705882, 0.0156863, 0.14902]
    wave_profile_analyticalLUT.ScalarRangeInitialized = 1.0

    # show data in view
    calculator2Display = Show(calculator2, renderView1)
    # trace defaults for the display properties.
    calculator2Display.Representation = 'Surface'
    calculator2Display.ColorArrayName = ['POINTS', 'wave_profile_analytical']
    calculator2Display.LookupTable = wave_profile_analyticalLUT
    calculator2Display.OSPRayScaleArray = 'wave_profile_analytical'
    calculator2Display.OSPRayScaleFunction = 'PiecewiseFunction'
    calculator2Display.SelectOrientationVectors = 'velocity'
    calculator2Display.ScaleFactor = 1.5
    calculator2Display.SelectScaleArray = 'wave_profile_analytical'
    calculator2Display.GlyphType = 'Arrow'
    calculator2Display.GlyphTableIndexArray = 'wave_profile_analytical'
    calculator2Display.DataAxesGrid = 'GridAxesRepresentation'
    calculator2Display.PolarAxes = 'PolarAxesRepresentation'
    calculator2Display.GaussianRadius = 0.75
    calculator2Display.SetScaleArray = ['POINTS', 'wave_profile_analytical']
    calculator2Display.ScaleTransferFunction = 'PiecewiseFunction'
    calculator2Display.OpacityArray = ['POINTS', 'wave_profile_analytical']
    calculator2Display.OpacityTransferFunction = 'PiecewiseFunction'
    #calculator2Display.SelectInputVectors = ['POINTS', 'velocity']
    #calculator2Display.WriteLog = ''

    # hide data in view
    Hide(contour2, renderView1)

    # show color bar/color legend
    calculator2Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # create a new 'Python Calculator'
    pythonCalculator2 = PythonCalculator(Input=calculator2)
    pythonCalculator2.Expression = ''

    # Properties modified on pythonCalculator2
    pythonCalculator2.Expression = 'max(wave_profile_analytical)'
    pythonCalculator2.ArrayName = 'max(wave_profile_analytical)'

    # show data in view
    pythonCalculator2Display = Show(pythonCalculator2, renderView1)
    # trace defaults for the display properties.
    pythonCalculator2Display.Representation = 'Surface'
    pythonCalculator2Display.ColorArrayName = ['POINTS', 'wave_profile_analytical']
    pythonCalculator2Display.LookupTable = wave_profile_analyticalLUT
    pythonCalculator2Display.OSPRayScaleArray = 'wave_profile_analytical'
    pythonCalculator2Display.OSPRayScaleFunction = 'PiecewiseFunction'
    pythonCalculator2Display.SelectOrientationVectors = 'velocity'
    pythonCalculator2Display.ScaleFactor = 1.5
    pythonCalculator2Display.SelectScaleArray = 'wave_profile_analytical'
    pythonCalculator2Display.GlyphType = 'Arrow'
    pythonCalculator2Display.GlyphTableIndexArray = 'wave_profile_analytical'
    pythonCalculator2Display.DataAxesGrid = 'GridAxesRepresentation'
    pythonCalculator2Display.PolarAxes = 'PolarAxesRepresentation'
    pythonCalculator2Display.GaussianRadius = 0.75
    pythonCalculator2Display.SetScaleArray = ['POINTS', 'wave_profile_analytical']
    pythonCalculator2Display.ScaleTransferFunction = 'PiecewiseFunction'
    pythonCalculator2Display.OpacityArray = ['POINTS', 'wave_profile_analytical']
    pythonCalculator2Display.OpacityTransferFunction = 'PiecewiseFunction'
    #pythonCalculator2Display.SelectInputVectors = ['POINTS', 'velocity']
    #pythonCalculator2Display.WriteLog = ''

    # hide data in view
    Hide(calculator2, renderView1)

    # show color bar/color legend
    pythonCalculator2Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()
    
    # Properties modified on pythonCalculator2
    pythonCalculator2.ArrayName = 'max_wave_profile_analytical'
    
    # update the view to ensure updated data information
    renderView1.Update()
    
    # set active source
    SetActiveSource(calculator2)
    
    # Properties modified on calculator2
    calculator2.ResultArrayName = 'wave_profile_analytical'
    
    # update the view to ensure updated data information
    renderView1.Update()
    data=servermanager.Fetch(pythonCalculator1)
    max_numerical.append(data.GetPointData().GetArray('max_wave_profile').GetValue(0))
    data2=servermanager.Fetch(pythonCalculator2)
    max_analytical.append(data2.GetPointData().GetArray('max_wave_profile_analytical').GetValue(0))
    time.append(periodic_solitary_wavesxmf.TimestepValues[i])
    animationScene1.GoToNext()
    wave_ele.append([time[-1],max_numerical[-1],max_analytical[-1]])
import numpy as np
np.savetxt("wave_ele.txt",np.array(wave_ele))
import os
import numpy as np
from visit import *

domain_folders=['PostProc']



## Immutables ##

alphaFolder='alpha'
velocityFolder='velocity'
interfaceFolder='interface'
pressureFolder='pressure'
volumefracFolder='volumeFrac'

ActiveWindow1=1
ActiveWindow2=2

## Clean VisIt
def Clean():
  DeleteAllPlots()
  ClearAllWindows()

for domain_folder in domain_folders:
  DeleteAllPlots()
  
  folderPath=f'/home/exy214/Documents/cavitation/data/jetting_ws_2025/runs/{domain_folder}/domain/'
  alphawritePath=f'{folderPath}{alphaFolder}'
  velocitywritePath=f'{folderPath}{velocityFolder}'
  interfacewritePath=f'{folderPath}{interfaceFolder}'
  pressurewritePath=f'{folderPath}{pressureFolder}'
  volumefracwritePath=f'{folderPath}{volumefracFolder}'

  os.makedirs(interfacewritePath,exist_ok=True)
  os.makedirs(pressurewritePath,exist_ok=True)
  os.makedirs(alphawritePath,exist_ok=True)
  os.makedirs(velocitywritePath,exist_ok=True)
  os.makedirs(volumefracwritePath,exist_ok=True)

  DefineScalarExpression('velocityY',"velocity[1]")

  OpenDatabase(f'{folderPath}domain_0.*.xdmf database', 0)
  print(f'Opened Database for {domain_folder}')
# Add Plot and Operators
  AddPlot("Pseudocolor", "diffuse volume fraction 1", 1, 1)
  AddOperator("Reflect", 1)
  AddOperator("Slice", 1)

# Reflect Attributes
  ReflectAtts = ReflectAttributes()
  ReflectAtts.octant = ReflectAtts.PXPYPZ  # Reflection type
  ReflectAtts.useXBoundary = 1
  ReflectAtts.specifiedX = 0
  ReflectAtts.useYBoundary = 0
  ReflectAtts.useZBoundary = 1
  ReflectAtts.specifiedZ = 0
  ReflectAtts.reflections = (1, 1, 0, 0, 0, 0, 0, 0)
  SetOperatorOptions(ReflectAtts, 0, 1)

# Slice Attributes
  SliceAtts = SliceAttributes()
  SliceAtts.originType = SliceAtts.Intercept
  SliceAtts.originPoint = (0, 0, 0)
  SliceAtts.normal = (0, 0, 1)
  SliceAtts.axisType = SliceAtts.ZAxis
  SliceAtts.project2d = 1
  SetOperatorOptions(SliceAtts, 1, 1)


  DrawPlots()
  nt=TimeSliderGetNStates()
  f1=open(f'{folderPath}Centroid.txt',"w+")

  for ts in range(0,nt):
    TimeSliderSetState(ts)
    time=Query("Time")
    t=GetQueryOutputValue()
    area=Query("Centroid")
    a=GetQueryOutputValue()
    f1.write("%15.12e %15.12e \n" % (t, a[1]))
    f1.flush()
    print(f'Centroid File Written, at time index: {ts}')
  f1.close()
  Clean()

  ## Volume Loop ##
  # Add Plot and Operators
  print(f'Volume post processing for {domain_folder}')
  AddPlot("Pseudocolor", "diffuse volume fraction 1", 1, 1)
  AddOperator("Slice", 1)
#Slice Attributes
  SliceAtts = SliceAttributes()
  SliceAtts.originType = SliceAtts.Intercept
  SliceAtts.originPoint = (0, 0, 0)
  SliceAtts.normal = (0, 0, 1)
  SliceAtts.axisType = SliceAtts.ZAxis
  SliceAtts.project2d = 1
  SetOperatorOptions(SliceAtts, 1, 1)

  AddOperator("Threshold", 1)
  ThresholdAtts = ThresholdAttributes()
  ThresholdAtts.outputMeshType = 0
  ThresholdAtts.boundsInputType = 0
  ThresholdAtts.listedVarNames = ("default")
  ThresholdAtts.zonePortions = (1)
  ThresholdAtts.lowerBounds = (0.1)
  ThresholdAtts.upperBounds = (1e+37)
  ThresholdAtts.defaultVarName = "diffuse volume fraction 1"
  ThresholdAtts.defaultVarIsScalar = 1
  ThresholdAtts.boundsRange = ("0.1:1e+37")
  SetOperatorOptions(ThresholdAtts, 2, 1)


#Revolve Attribtue
  AddOperator("Revolve", 1)
  RevolveAtts = RevolveAttributes()
  RevolveAtts.meshType = RevolveAtts.ZR  # Auto, XY, RZ, ZR
  RevolveAtts.autoAxis = 1
  RevolveAtts.axis = (1, 0, 0)
  RevolveAtts.startAngle = 0
  RevolveAtts.stopAngle = 360
  RevolveAtts.steps = 120
  SetOperatorOptions(RevolveAtts, 2, 1)

  DrawPlots()
  nt=TimeSliderGetNStates()
  print("Opening file")
  f1=open(f"{folderPath}Volume.txt","w+")

  for ts in range(0,nt):
    TimeSliderSetState(ts)
    time=Query("Time")
    t=GetQueryOutputValue()
    area=Query("Weighted Variable Sum")
    a=GetQueryOutputValue()
    f1.write("%15.12e %15.12e \n" % (t, a))
    f1.flush()
    print(f"File written at timestep: {ts}")
  f1.close()
  
  Clean()

  ## Interface loop ##
  # Add Plot and Operators
  print(f'Interface post processing loop for {domain_folder}')
  AddPlot("Pseudocolor", "diffuse volume fraction 1", 1, 1)
  AddOperator("Slice", 1)
  AddOperator("Isosurface", 1)

  # Slice Attributes
  SliceAtts = SliceAttributes()
  SliceAtts.originType = SliceAtts.Intercept
  SliceAtts.originPoint = (0, 0, 0)
  SliceAtts.normal = (0, 0, 1)
  SliceAtts.axisType = SliceAtts.ZAxis
  SliceAtts.project2d = 1
  SetOperatorOptions(SliceAtts, 1, 1)

  IsosurfaceAtts = IsosurfaceAttributes()
  IsosurfaceAtts.contourNLevels = 10
  IsosurfaceAtts.contourValue = (0.5)
  IsosurfaceAtts.contourPercent = ()
  IsosurfaceAtts.contourMethod = IsosurfaceAtts.Value  # Level, Value, Percent
  IsosurfaceAtts.minFlag = 0
  IsosurfaceAtts.min = 0
  IsosurfaceAtts.maxFlag = 0
  IsosurfaceAtts.max = 1
  IsosurfaceAtts.scaling = IsosurfaceAtts.Linear  # Linear, Log
  IsosurfaceAtts.variable = "default"
  SetOperatorOptions(IsosurfaceAtts, 1, 1)

  # End spontaneous state
  DrawPlots()

  nt=TimeSliderGetNStates()
  SaveWindowAtts = SaveWindowAttributes()
  SaveWindowAtts.outputToCurrentDirectory=0
  SaveWindowAtts.outputDirectory = interfacewritePath
  SaveWindowAtts.fileName = "interface"
  SaveWindowAtts.family = 1
  SaveWindowAtts.format = SaveWindowAtts.CURVE

  for ts in range(0,nt):
    TimeSliderSetState(ts)
    SaveWindowAtts.fileName=f"Interface"
    SetSaveWindowAttributes(SaveWindowAtts)
    SaveWindow()
  
  Clean()
  print(f'Centerline extraction for {domain_folder}')
  ## Pressure Loop ##
  AddPlot("Pseudocolor", "pressure", 1, 1)
  DrawPlots()

  SetActiveWindow(ActiveWindow1)
  nt=TimeSliderGetNStates()
  for ts in range(0,nt):
      SetActiveWindow(ActiveWindow1)
      TimeSliderSetState(ts)
      Lineout([0,6e-3,0],[0,10e-3,0])
      SetActiveWindow(ActiveWindow2)
      SetActivePlots(0)
      Data=GetPlotInformation()["Curve"]
      Pressure=Data[1::2]
      Position=Data[::2]
      print(f"Writing file at timestep: {ts}")
      np.savetxt(f"{pressurewritePath}/pressure{ts:04d}.txt", Pressure)
      if ts == 0:
        np.savetxt(f"{pressurewritePath}/position{ts:04d}.txt", Position)
  
  Clean()

  ## Alpha Loop ##
  AddPlot("Pseudocolor", "diffuse volume fraction 1", 1, 1)
  DrawPlots()
  ActiveWindow1=2
  ActiveWindow2=3
  SetActiveWindow(ActiveWindow1)
  nt=TimeSliderGetNStates()
  for ts in range(0,nt):
      SetActiveWindow(ActiveWindow1)
      TimeSliderSetState(ts)
      Lineout([0,6e-3,0],[0,10e-3,0])
      SetActiveWindow(ActiveWindow2)
      SetActivePlots(0)
      Data=GetPlotInformation()["Curve"]
      Alpha=Data[1::2]
      Position=Data[::2]
      print(f'Writing Files to folder at timestep {ts}')
      np.savetxt(f"{alphawritePath}/alpha{ts:04d}.txt", Alpha)
      if ts == 0:
        np.savetxt(f"{alphawritePath}/position{ts:04d}.txt", Position)

  Clean()

  ## Velocity Loop ##
  AddPlot("Pseudocolor", "velocityY", 1, 1)
  DrawPlots()
  ActiveWindow1=3
  ActiveWindow2=4
  SetActiveWindow(ActiveWindow1)
  nt=TimeSliderGetNStates()
  for ts in range(0,nt):
      SetActiveWindow(ActiveWindow1)
      TimeSliderSetState(ts)
      Lineout([0,4e-3,0],[0,12e-3,0],500)
      SetActiveWindow(ActiveWindow2)
      SetActivePlots(0)
      Data=GetPlotInformation()["Curve"]
      Velocity=Data[1::2]
      Position=Data[::2]
      print(f'Writing file at timestep: {ts}')
      np.savetxt(f"{velocitywritePath}/velocity{ts:04d}.txt", Velocity)
      if ts == 0:
        np.savetxt(f"{velocitywritePath}/position{ts:04d}.txt", Position)
  
  Clean()

  ## Pressure/Schlieren Imaging Loop ##

  # Add Plot and Operators
  AddPlot("Pseudocolor", "pressure", 1, 0)
  AddOperator("Slice", 0)
  SetActivePlots(0)

  ReflectAtts = ReflectAttributes()
  ReflectAtts.octant = ReflectAtts.PXPYPZ  # PXPYPZ, NXPYPZ, PXNYPZ, NXNYPZ, PXPYNZ, NXPYNZ, PXNYNZ, NXNYNZ
  ReflectAtts.useXBoundary = 1
  ReflectAtts.specifiedX = 0
  ReflectAtts.useYBoundary = 1
  ReflectAtts.specifiedY = 0
  ReflectAtts.useZBoundary = 0
  ReflectAtts.specifiedZ = 0
  ReflectAtts.reflections = (1, 0, 0, 0, 0, 0, 0, 0)
  ReflectAtts.planePoint = (0, 0, 0)
  ReflectAtts.planeNormal = (0, 0, 0)
  ReflectAtts.reflectType = ReflectAtts.Axis  # Plane, Axis
  SetOperatorOptions(ReflectAtts, 1, 0)

  # Slice Attributes
  SliceAtts = SliceAttributes()
  SliceAtts.originType = SliceAtts.Intercept
  SliceAtts.originPoint = (0, 0, 0)
  SliceAtts.normal = (0, 0, 1)
  SliceAtts.axisType = SliceAtts.ZAxis
  SliceAtts.project2d = 1
  SetOperatorOptions(SliceAtts, 1, 0)

  #Pseudoclolor Attribtutes
  PseudocolorAtts=PseudocolorAttributes()
  PseudocolorAtts.colorTableName="YlGnBu"
  SetPlotOptions(PseudocolorAtts)

  AddPlot("Pseudocolor","schlieren",1,0)

  AddOperator("Slice",0)
  AddOperator("Reflect",0)
  SetActivePlots(1)

  # Slice Attributes
  SliceAtts = SliceAttributes()
  SliceAtts.originType = SliceAtts.Intercept
  SliceAtts.originPoint = (0, 0, 0)
  SliceAtts.normal = (0, 0, 1)
  SliceAtts.axisType = SliceAtts.ZAxis
  SliceAtts.project2d = 1
  SetOperatorOptions(SliceAtts, 0, 0)

  ReflectAtts = ReflectAttributes()
  ReflectAtts.octant = ReflectAtts.PXPYPZ  # PXPYPZ, NXPYPZ, PXNYPZ, NXNYPZ, PXPYNZ, NXPYNZ, PXNYNZ, NXNYNZ
  ReflectAtts.useXBoundary = 1
  ReflectAtts.specifiedX = 0
  ReflectAtts.useYBoundary = 1
  ReflectAtts.specifiedY = 0
  ReflectAtts.useZBoundary = 0
  ReflectAtts.specifiedZ = 0
  ReflectAtts.reflections = (0, 1, 0, 0, 0, 0, 0, 0)
  ReflectAtts.planePoint = (0, 0, 0)
  ReflectAtts.planeNormal = (0, 0, 0)
  ReflectAtts.reflectType = ReflectAtts.Axis  # Plane, Axis
  SetOperatorOptions(ReflectAtts, 1, 0)

  #Pseudoclolor Attribtutes
  PseudocolorAtts=PseudocolorAttributes()
  PseudocolorAtts.min=0
  PseudocolorAtts.max=600000
  PseudocolorAtts.colorTableName="xray"
  SetPlotOptions(PseudocolorAtts)

  '''
  AnnotationAtts=AnnotationAttributes()
  AnnotationAtts.axes2D.visible=0
  AnnotationAtts.axes3D.visible=0
  AnnotationAtts.axesArray.visible=0

  AnnotationAtts.timeInfoFlag=0
  AnnotationAtts.databaseInfoFlag=0
  AnnotationAtts.userInfoFlag=0
  AnnotationAtts.legendInfoFlag=0
  SetAnnotationAttributes(AnnotationAtts)
  '''
  DrawPlots()

  # Begin spontaneous state
  View2DAtts = View2DAttributes()
  View2DAtts.windowCoords = (-0.000815831, 0.000815831, 0.0071878, 0.0088122)
  View2DAtts.viewportCoords = (0.2, 0.95, 0.15, 0.95)
  View2DAtts.fullFrameActivationMode = View2DAtts.Auto  # On, Off, Auto
  View2DAtts.fullFrameAutoThreshold = 100
  View2DAtts.xScale = View2DAtts.LINEAR  # LINEAR, LOG
  View2DAtts.yScale = View2DAtts.LINEAR  # LINEAR, LOG
  View2DAtts.windowValid = 1
  SetView2D(View2DAtts)

  # End spontaneous state
  nt=TimeSliderGetNStates()
  SaveWindowAtts = SaveWindowAttributes()
  SaveWindowAtts.outputToCurrentDirectory = 0
  SaveWindowAtts.outputDirectory = f'{volumefracwritePath}'
  SaveWindowAtts.format = SaveWindowAtts.PNG
  SaveWindowAtts.height = 1024
  SaveWindowAtts.width= 1024

  for ts in range(0,nt):
    TimeSliderSetState(ts)
    SaveWindowAtts.fileName=f"volumeFrac"
    SetSaveWindowAttributes(SaveWindowAtts)
    print(f'Writing image at time step: {ts}')
    SaveWindow()

  Clean()
exit()

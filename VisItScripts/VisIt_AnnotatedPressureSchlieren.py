from visit import *
import numpy as np
import os

domain_folders=['PostProc']
for domain_folder in domain_folders:
   
  # Create folder and write paths, then create the write folders
  folderPath=f'/home/exy214/Documents/cavitation/data/jetting_ws_2025/runs/{domain_folder}/domain/'
  writeFolder=f'volumeFrac_{domain_folder}'
  writePath=f'{folderPath}{writeFolder}'
  os.makedirs(writePath,exist_ok=True)


  DeleteAllPlots()
  OpenDatabase(f'{folderPath}domain_0.*.xdmf database', 0)

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
  SaveWindowAtts.outputDirectory = f'{writePath}'
  SaveWindowAtts.format = SaveWindowAtts.PNG
  SaveWindowAtts.height = 1024
  SaveWindowAtts.width= 1024

  for ts in range(0,nt):
    TimeSliderSetState(ts)
    SaveWindowAtts.fileName=f"volumeFrac"
    SetSaveWindowAttributes(SaveWindowAtts)
    print(f'Writing image at time step: {ts}')
    SaveWindow()
exit()

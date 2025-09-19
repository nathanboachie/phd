from visit import *
import os

domain_folders=['PostProc']

for domain_folder in domain_folders:
  DeleteAllPlots()

  writeFolder=f'Interface_{domain_folder}'
  folderPath=f'/home/exy214/Documents/cavitation/data/jetting_ws_2025/runs/{domain_folder}/domain/'
  writePath=f'{folderPath}{writeFolder}'
  os.makedirs(writePath,exist_ok=True)

  OpenDatabase(f"{folderPath}domain_0.*.xdmf database", 0)

  # Add Plot and Operators
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
  SaveWindowAtts.outputDirectory = writePath
  SaveWindowAtts.fileName = "interface"
  SaveWindowAtts.family = 1
  SaveWindowAtts.format = SaveWindowAtts.CURVE 

  for ts in range(0,nt):
    TimeSliderSetState(ts)
    SaveWindowAtts.fileName=f"Interface"
    SetSaveWindowAttributes(SaveWindowAtts)
    SaveWindow()
exit()

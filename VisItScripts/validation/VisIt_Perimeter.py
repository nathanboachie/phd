from visit import *

domain_folders=['PostProc']

for domain_folder in domain_folders:
  DeleteAllPlots()
  folderPath=f'/home/exy214/Documents/cavitation/data/jetting_ws_2025/runs/{domain_folder}/domain/'
  OpenDatabase(f"{folderPath}domain_0.*.xdmf database", 0)
  print('Opened Database')
# Add Plot and Operators
  AddPlot("Pseudocolor", "diffuse volume fraction 1", 1, 1)
  AddOperator("Reflect", 1)
  AddOperator("Slice", 1)
  AddOperator("Isosurface",1)

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
  DrawPlots()
  nt=TimeSliderGetNStates()
  f1=open(f"{folderPath}Perimeter.txt","w+")

  for ts in range(0,nt):
    TimeSliderSetState(ts)
    time=Query("Time")
    t=GetQueryOutputValue()
    area=Query("Total Length")
    a=GetQueryOutputValue()
    print(f'Printing length at timestep {nt}')
    f1.write("%15.12e %15.12e \n" % (t, a))
    f1.flush()
  f1.close()
exit()

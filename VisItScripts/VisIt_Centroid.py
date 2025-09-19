from visit import *

domain_folders=['PostProc']

for domain_folder in domain_folders:
  DeleteAllPlots()
  folderPath=f'/home/exy214/Documents/cavitation/data/jetting_ws_2025/runs/{domain_folder}/domain/'

  OpenDatabase(f'{folderPath}domain_0.*.xdmf database', 0)
  print('Opened Database')
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
    print(f'File Written, at time index: {ts}')
  f1.close()
exit()

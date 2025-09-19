from visit import *

domain_folders=['PostProc']

for domain_folder in domain_folders:
  DeleteAllPlots()
  folderPath=f'/home/exy214/Documents/cavitation/data/jetting_ws_2025/runs/{domain_folder}/domain/'

  OpenDatabase(f'{folderPath}domain_0.*.xdmf database', 0)
  print("Database loaded")
# Add Plot and Operators
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
exit()

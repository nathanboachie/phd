from visit import *
import numpy as np
import os

domain_folders=['PostProc']

ActiveWindow1=1
ActiveWindow2=2
for domain_folder in domain_folders:

  DeleteAllPlots()

  writeFolder=f'pressure'
  folderPath=f'/home/exy214/Documents/cavitation/data/jetting_ws_2025/runs/{domain_folder}/domain/'
  writePath=f'{folderPath}{writeFolder}'
  os.makedirs(writePath,exist_ok=True)

  OpenDatabase(f"{folderPath}domain_0.*.xdmf database", 0)
  AddPlot("Pseudocolor", "pressure", 1, 1)
  DrawPlots()

  SetActiveWindow(ActiveWindow1)
  nt=TimeSliderGetNStates()

  for ts in range(0,nt):
      SetActiveWindow(ActiveWindow1)
      TimeSliderSetState(ts)
      Lineout([0,4e-3,0],[0,12e-3,0],500)
      SetActiveWindow(ActiveWindow2)
      SetActivePlots(0)
      Data=GetPlotInformation()["Curve"]
      Pressure=Data[1::2]
      Position=Data[::2]
      print(f"Writing file at timestep: {ts}")
      np.savetxt(f"{writePath}/pressure{ts:04d}.txt", Pressure)
      if ts == 0:
        np.savetxt(f"{writePath}/position{ts:04d}.txt", Position)
  ActiveWindow1+=1
  ActiveWindow2+=1
exit()

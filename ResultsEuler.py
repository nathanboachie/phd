import numpy as np
import os
import matplotlib.pyplot as plt



def AppendArrays(arr1,val1,arr2,val2,arr3,val3):
    if(val1 == None):
        return
    arr1.append(val1)
    arr2.append(val2)
    arr3.append(val3)
    return

plt.rcParams.update({
    'font.size': 9,
    'font.family': 'serif',
    'axes.labelsize': 9,
    'axes.titlesize': 9,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 8,
    'lines.linewidth': 1,
    'lines.markersize': 4,
    'axes.linewidth': 0.8,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.major.size': 4,
    'ytick.major.size': 4,
    'savefig.dpi': 600
})

def PlotSpecs(xlabel,ylabel):
        fig, ax = plt.subplots(figsize=(3.5,2.5))
        # Show ticks on all 4 sides
        ax.tick_params(
        axis='both',
        which='both',
        direction='in',
        top=True,
        right=True
        )
        # --- Labels and Limits ---
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        return fig,ax

def JettingImpactPlot(tJetting,minVolumes,pressures,dX,Cl):
        fig,ax = PlotSpecs(xlabel=r'$P (MPa)$',ylabel=r'$\Delta T_{jet}/ \Delta T_{Collapse}$')
        tShockImpact=dX/Cl
        dParamJet=(tJetting-tShockImpact)/(minVolumes-tShockImpact)
        ax.plot(pressures,dParamJet,'-',color='black')
        fig.savefig(f"JettingImpactPlot.png")
        plt.close(fig)

def ShockVolumePlot(minVolumes,pressures,R0,Cl):
        fig,ax = PlotSpecs(xlabel=r'$P (MPa)$',ylabel=r'$T_{diffract}/T_{Collapse}$')
        tShockDiffract=2*R0/Cl
        dParamJet=tShockDiffract/minVolumes
        ax.plot(pressures,dParamJet,'-',color='black')
        fig.savefig(f"ShockVolumePlot.png")
        plt.close(fig)



class simulationRun:
    def __init__(self,filepath,pressure,tfinal,name,particle_speed):
        #Initial Stuff
        self.filepath=filepath
        self.pressure=pressure
        self.tfinal=tfinal
        self.name=name
        self.R0=0.000338
        self.Cl=1462
        self.Rho_l=998
        self.ParticleSpeed=particle_speed
        # For Area and Sphericity splots
        self.AreaArr=np.loadtxt(os.path.join(filepath, "Area.txt"))
        self.VolumeArr=np.loadtxt(os.path.join(filepath, "Volume.txt"))
        self.SurfaceAreaArr=np.loadtxt(os.path.join(filepath, "SurfaceArea.txt"))
        self.SphericityArr=np.column_stack((self.VolumeArr[:,0],(np.pi)**(1/3)*(6*self.VolumeArr[:,1])**(2/3)/(self.SurfaceAreaArr[:,1])))
        
        #For Pressure Matrix
        self.PressureTemplate = os.path.join(filepath, "pressure", "pressure{:04d}.txt")
        self.NumRuns = len(self.AreaArr)
        self.PressureList = [np.loadtxt(self.PressureTemplate.format(i)) for i in range(self.NumRuns)]
        self.PressureMatrix=np.stack(self.PressureList,axis=1)

        #For Velocity Matrix
        self.VelocityTemplate = os.path.join(filepath, "velocity", "velocity{:04d}.txt")
        self.NumRuns = len(self.AreaArr)
        self.VelocityList = [np.loadtxt(self.VelocityTemplate.format(i)) for i in range(self.NumRuns)]
        self.VelocityMatrix=np.stack(self.VelocityList,axis=1)

        #For Diffuse Volume Fraction Matrix
        self.AlphaTemplate = os.path.join(filepath, "alpha", "alpha{:04d}.txt")
        self.NumRuns = len(self.AreaArr)
        self.AlphaList = [np.loadtxt(self.AlphaTemplate.format(i)) for i in range(self.NumRuns)]
        self.AlphaMatrix=np.stack(self.AlphaList,axis=1)

        #For Position Matrix
        self.PositionMatrix=np.loadtxt(os.path.join(filepath,"pressure/position0000.txt"))+0.022 # Starts measuring from 0.022

        #For Keller Miksis array
        self.KellerMiksis=np.loadtxt(os.path.join(filepath,"KM.txt"))
        self.KM_U=self.KellerMiksis[:,2]
        self.KM_T=self.KellerMiksis[:,0]
        self.KM_R=self.KellerMiksis[:,1]
            
        # Time Matrix
        self.tArray=np.linspace(0,tfinal,self.NumRuns)

        # Interface Values
        self.yInterfaceUpstream = np.zeros(((self.NumRuns,1)))
        self.uInterfaceUpstream = np.zeros(((self.NumRuns,1)))
        self.pInterfaceUpstream = np.zeros(((self.NumRuns,1)))
        self.aInterfaceUpstream = np.zeros(((self.NumRuns,1)))
        self.yInterfaceDownstream = np.zeros(((self.NumRuns,1)))
        self.uInterfaceDownstream = np.zeros(((self.NumRuns,1)))
        self.pInterfaceDownstream = np.zeros(((self.NumRuns,1)))
        self.aInterfaceDownstream = np.zeros(((self.NumRuns,1)))


        # Trimmed Interface Values
        '''
        Why? This is a trim matrix that ends at the point of jetting
        The data looks weird including information past jetting
        '''
        self.TrimyInterfaceUpstream = np.zeros(((self.NumRuns,1)))
        self.TrimuInterfaceUpstream = np.zeros(((self.NumRuns,1)))
        self.TrimpInterfaceUpstream = np.zeros(((self.NumRuns,1)))
        self.TrimaInterfaceUpstream = np.zeros(((self.NumRuns,1)))
        self.TrimyInterfaceDownstream = np.zeros(((self.NumRuns,1)))
        self.TrimuInterfaceDownstream = np.zeros(((self.NumRuns,1)))
        self.TrimpInterfaceDownstream = np.zeros(((self.NumRuns,1)))
        self.TrimaInterfaceDownstream = np.zeros(((self.NumRuns,1)))
        self.TrimtArray = np.zeros(((self.NumRuns,1)))
        self.TrimPressureMatrix=np.zeros((self.NumRuns,len(self.PositionMatrix)))
        self.TrimVelocityMatrix=np.zeros((self.NumRuns,len(self.PositionMatrix)))
        self.TrimAlphaMatrix=np.zeros((self.NumRuns,len(self.PositionMatrix)))
        
        # Values
        self.y0=0.025
        self.tJetting=0

    
    '''
    Interface, look at interface, if there are no points w/volume fraction above 0.5 it means there is no interface along centerline
    The data at this point is all from the interface
    Find the most upstream point and downstream point and extract the values at these locations, which have a volume fraction above a point
    The interface value is then one-dimensional
    Repeats for downstream point
    '''
    def Interface(self):
        for i in range(self.NumRuns):
            Iarr=np.where(self.AlphaMatrix[:,i]>0.9)
            if(len(Iarr[0])==1):
                self.tJetting=self.tArray[i] 
                print(f'Jetting Achieved @ {self.tArray[i]}')
                break
            else:
                YGasPosition=self.PositionMatrix[Iarr]
                YMax=np.max(YGasPosition)
                I=int(np.where(self.PositionMatrix==YMax)[0][0])
                self.yInterfaceDownstream[i]=self.PositionMatrix[I]
                self.uInterfaceDownstream[i]=self.VelocityMatrix[:,i][I]
                self.pInterfaceDownstream[i]=self.PressureMatrix[:,i][I]
                self.aInterfaceDownstream[i]=self.AlphaMatrix[:,i][I]

                YMin=np.min(YGasPosition)
                I=int(np.where(self.PositionMatrix==YMin)[0][0])
                self.yInterfaceUpstream[i]=self.PositionMatrix[I]
                self.uInterfaceUpstream[i]=self.VelocityMatrix[:,i][I]
                self.pInterfaceUpstream[i]=self.PressureMatrix[:,i][I]
                self.aInterfaceUpstream[i]=self.AlphaMatrix[:,i][I]
    '''
    Returns collapse index based on the point where the differnce in upstream and downstream velocities begins to change to a certain degree
    '''
    def CollapseStartTime(self):
        diffUpstream=np.diff(self.TrimuInterfaceUpstream.flatten())
        UpstreamCollapse=np.min(np.where(diffUpstream > 0.1))
        return UpstreamCollapse
    '''
    Trimming arryas based upon length of velocity interface where nothing is being output
    '''
    def TrimArrays(self):
        self.TrimuInterfaceUpstream=np.trim_zeros(self.uInterfaceUpstream,'b')
        self.TrimuInterfaceDownstream=np.trim_zeros(self.uInterfaceDownstream,'b')
        self.TrimyInterfaceUpstream=self.yInterfaceUpstream[0:len(self.TrimuInterfaceDownstream)]
        self.TrimyInterfaceDownstream=self.yInterfaceDownstream[0:len(self.TrimuInterfaceDownstream)]
        self.TrimaInterfaceUpstream=self.aInterfaceUpstream[0:len(self.TrimuInterfaceDownstream)]
        self.TrimaInterfaceDownstream=self.aInterfaceDownstream[0:len(self.TrimuInterfaceDownstream)]
        self.TrimtArray=self.tArray[0:len(self.TrimuInterfaceDownstream)]
        self.TrimPressureMatrix=self.PressureMatrix[:,:len(self.TrimuInterfaceDownstream)]
        self.TrimVelocityMatrix=self.VelocityMatrix[:,:len(self.TrimuInterfaceDownstream)]
        self.TrimAlphaMatrix=self.AlphaMatrix[:,:len(self.TrimuInterfaceDownstream)]

    def BokmanPressure(self):
        PjetMPa=(self.Rho_l*(0.915*0.05*self.Cl*0.5)**2 + 1e5)/1e6
        print(f'Bokman jetting pressure has to be greater than {PjetMPa} MPa')



    def MinVolumeTime(self):
        MinVolumeIdx=np.argmin(self.VolumeArr[:,1])
        print(MinVolumeIdx)
        tMinVolume=self.tArray[MinVolumeIdx]
        print(tMinVolume)
        return tMinVolume

    def BokmanCollapse(self):
        tCollapse=self.MinVolumeTime()
        tShockDiffract=2*self.R0/self.Cl
        if(tShockDiffract/tCollapse > 0.05):
            print(f"Bokman predicts collapse for {self.pressure}")
        else:
            print(f"Bokman et al. do not predict collapse for {self.pressure}")

    def PlotSpecs(self,xlabel,ylabel):
        fig, ax = plt.subplots(figsize=(3.5,2.5))
        # Show ticks on all 4 sides
        ax.tick_params(
        axis='both',
        which='both',
        direction='in',
        top=True,
        right=True
        )
        # --- Labels and Limits ---
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        return fig, ax
    def VelocityPlot(self):
        fig,ax = self.PlotSpecs(xlabel=r'$t/(R_{0}/C_{l})$',ylabel=r'$u$')
        tCollapse=self.TrimtArray[self.CollapseStartTime()]
        ax.plot((self.TrimtArray-tCollapse)/(self.R0/self.Cl),-1*self.TrimuInterfaceDownstream,'--',color='black',label='Downstream')
        ax.plot((self.TrimtArray-tCollapse)/(self.R0/self.Cl),self.TrimuInterfaceUpstream,'-',color='black',label='Upstream')
        ax.plot(self.KM_T/(self.R0/self.Cl),-1*self.KM_U,linestyle='dashdot',color='black')
        ax.plot(self.KM_T/(self.R0/self.Cl),self.KM_U,linestyle='dashdot',color='black',label='Keller Miksis')
        ax.legend()
        if(self.tJetting!=0):
            ax.set_xlim(0,(self.tJetting-tCollapse)/(self.R0/self.Cl))
        else:
            ax.set_xlim(0,None) 
        fig.savefig(os.path.join(self.filepath, f"VelocityPlot_{self.name}.png"))
        plt.close(fig)
    def PositionPlot(self):
        fig,ax = self.PlotSpecs(xlabel=r'$t/\tau_{c}$',ylabel=r'$\frac{(y-y_{0})}{R_{0}}$')
        ax.plot(self.TrimtArray/(self.R0/self.Cl),(self.TrimyInterfaceDownstream-self.y0)/self.R0,'--',color='black',label='Downstream')
        ax.plot(self.TrimtArray/(self.R0/self.Cl),(self.TrimyInterfaceUpstream-self.y0)/self.R0,'-',color='black',label='Upstream')
        ax.legend()
        fig.savefig(os.path.join(self.filepath, f"PositionPlot_{self.name}.png"))
        plt.close(fig)
    def VelocityContour(self):
        fig,ax = self.PlotSpecs(xlabel=r'$\frac{(y-y_{0})}{R_{0}}$',ylabel=r'$t/\tau_{c}$')
        ax.contourf(self.TrimtArray/(self.R0/self.Cl),(self.PositionMatrix.squeeze()-self.y0)/self.R0,self.TrimVelocityMatrix/(self.Cl**2),cmap='Greys')
        ax.plot(self.TrimtArray/(self.R0/self.Cl),(self.TrimyInterfaceDownstream-self.y0)/self.R0,'--',color='black')
        ax.plot(self.TrimtArray/(self.R0/self.Cl),(self.TrimyInterfaceUpstream-self.y0)/self.R0,'-',color='black')
        fig.savefig(os.path.join(self.filepath, f"VelocityContour_{self.name}.png"))
        plt.close(fig)
    def PressureContour(self):
        fig,ax = self.PlotSpecs(xlabel=r'$\frac{(y-y_{0})}{R_{0}}$',ylabel=r'$t/\tau_{c}$')
        ax.contourf(self.TrimtArray/(self.R0/self.Cl),(self.PositionMatrix.squeeze()-self.y0)/self.R0,self.TrimPressureMatrix,cmap='Greys')
        ax.plot(self.TrimtArray/(self.R0/self.Cl),(self.TrimyInterfaceDownstream-self.y0)/self.R0,'--',color='black')
        ax.plot(self.TrimtArray/(self.R0/self.Cl),(self.TrimyInterfaceUpstream-self.y0)/self.R0,'-',color='black')
        fig.savefig(os.path.join(self.filepath, f"PressureContour_{self.name}.png"))
        plt.close(fig)
    def VolumePlot(self):
        fig,ax = self.PlotSpecs(xlabel=r'$t/\tau_{c}$',ylabel=r'$V/V_{0}$')
        ax.plot(self.VolumeArr[:,0]/(self.R0/self.Cl),self.VolumeArr[:,1]/self.VolumeArr[0,1],'-',color='black')
        if(self.tJetting != 0):
            ax.axvline(x=self.tJetting/(self.R0/self.Cl),linestyle='--',color='black')
        fig.savefig(os.path.join(self.filepath, f"VolumePlot_{self.name}.png"))
        plt.close(fig)
    def SphericityPlot(self):
        fig,ax = self.PlotSpecs(xlabel=r'$t/\tau_{c}$',ylabel=r'$Sphericity$')
        ax.plot(self.SphericityArr[:,0]/(self.R0/self.Cl),self.SphericityArr[:,1],'-',color='black')
        fig.savefig(os.path.join(self.filepath, f"SphericityPlot_{self.name}.png"))
        plt.close(fig)
    def JettingParameters(self):
        if(self.tJetting==0):
            return None,None
        MinVolumeTime=self.MinVolumeTime()
        return self.tJetting, MinVolumeTime



    def Run(self):
        self.Interface()
        self.TrimArrays()
        self.BokmanCollapse()
        self.BokmanPressure()
        print('Plotting Volume')
        self.VolumePlot()
        print('Plotting Sphericity')
        self.SphericityPlot()
        print('Plotting Interface Velocity')
        self.VelocityPlot()
        print('Plotting Interface Position')
        self.PositionPlot()
        print('Plotting Velocity Contour')
        self.VelocityContour()
        print('Plotting Pressure Contour')
        self.PressureContour()
        tJ,mV=self.JettingParameters()
        return tJ,mV
jettingTimes=[]*4
minVolumes=[]*4
pressures=[]*4
WS100e6=simulationRun('/home/exy214/Documents/cavitation/data/jetting_ws_2025/WS_100e6_Euler','100e6',5e-6,'WS100e6',65.4931)
tJ,mV=WS100e6.Run()
AppendArrays(jettingTimes,tJ,minVolumes,mV,pressures,100e6)
WS10e6=simulationRun('/home/exy214/Documents/cavitation/data/jetting_ws_2025/WS_10e6_Euler','10e6',8e-6,'WS10e6',6.7532)
tJ,mV=WS10e6.Run()
AppendArrays(jettingTimes,tJ,minVolumes,mV,pressures,10e6)
WS1e6=simulationRun('/home/exy214/Documents/cavitation/data/jetting_ws_2025/WS_1e6_Euler','1e6',2e-5,'WS1e6',0.6165)
tJ,mV=WS1e6.Run()
AppendArrays(jettingTimes,tJ,minVolumes,mV,pressures,1e6)
WS500e5=simulationRun('/home/exy214/Documents/cavitation/data/jetting_ws_2025/WS_500e5_Euler','500e5',5e-5,'WS500e5',0.2741)
tJ,mv=WS500e5.Run()
AppendArrays(jettingTimes,tJ,minVolumes,mV,pressures,500e5)
jettingTimes=np.array(jettingTimes)
minVolumes=np.array(minVolumes)
pressures=np.array(pressures)
JettingImpactPlot(jettingTimes,minVolumes,pressures,0.004662,1462)
ShockVolumePlot(minVolumes,pressures,0.000338,1462)

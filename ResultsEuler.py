import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from KM import KellerMiksis



def AppendArrays(arr1,val1,arr2,val2,arr3,val3,arr4,val4,arr5,val5):
    """   AppendArrays(arr[i],val[i]) for i in range(N) is used to append a value to an array
    Parameters:
        arrI: Array 'I' you wish to append to 
        valI: Value 'I' you wish to add to 
    """
    if(val1 == None):
        return
    arr1.append(val1)
    arr2.append(val2)
    arr3.append(val3)
    arr4.append(val4)
    arr5.append(val5)
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
    """ PlotSpecs(xlabel,ylabel) --> Creates figure and axes objects with set plot specification appropriate
                                     for this level of plotting
    Parameters:
        xlabel: r-String of the latex definition of your x label
        ylabel: r-String of the latex defintion of your y label
    Return:
        fig: Figure object
        ax: Axes object
    """
    fig, ax = plt.subplots(figsize=(3,3))
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
    """ JettingImpactPlot(tJetting,minVolumes,pressures,dX,Cl) --> Creates a plot of shock pressure ratios vs deltaTjetting/deltaTcollapse against
                                                                   Basically plots shock pressure vs when the re-entrnant jet forms vs when the bubble collapses
    Parameters:
        tJetting: Array of times when re-entrant jet formed
        minVolumes: Array of times when bubble collapsed (reached minimum volume)
        pressures: Array of shock overpressures 
        dX: Distance until proximal side of bubble (most upstream point)
        Cl: Speed of sound in water
    """
    fig,ax = PlotSpecs(xlabel=r'$p_{s}/p_{0}$',ylabel=r'$\Delta T_{j}/ \Delta T_{c}$')
    tShockImpact=dX/Cl
    dParamJet=(tJetting-tShockImpact)/(minVolumes-tShockImpact)
    ax.plot(pressures/1e5,dParamJet,color='black',linestyle=None,marker='o',markersize=3)
    fig.savefig(f"JettingImpactPlot.png",bbox_inches='tight')
    plt.close(fig)

def ShockVolumePlot(minVolumes,pressures,R0,Cl):
    """ ShockVolumePlot(minVolumes,pressures,R0,Cl) --> Creates a plot of shock pressure ratios vs the time taken for the shock to diffract / minimum volume (collapse)
    Parameters:
        minVolumes: Array of times when bubble collapsed (reached minimum volume)
        pressures: Array of shock overpressures
        R0: Initial bubble radius
        Cl: Speed of sound in water
    """
    fig,ax = PlotSpecs(xlabel=r'$p_{s}/p_{0}$',ylabel=r'$T_{d}/T_{c}$')
    tShockDiffract=2*R0/Cl
    dParamJet=tShockDiffract/minVolumes
    ax.plot(pressures/1e5,dParamJet,color='black',linestyle=None,marker='o',markersize=3)
    ax.axhline(y=0.05,linestyle='--',color='Black')
    fig.tight_layout()
    fig.savefig(f"ShockVolumePlot.png",bbox_inches='tight')
    plt.close(fig)

def maxJetVelocity(maxJetting,pressures,Cl):
    """ maxJetVelocity(maxJetting,pressures,Cl) --> Creates a plot of shock pressures vs the maximum jet velocity
    Parameters:
        maxJetting: Maximum jet velocity
        pressures: Array of shock overpressures
        Cl: Speed of sound in water (potentially unused)
    """
    fig,ax = PlotSpecs(xlabel=r'$p_{s}/p_{0}$',ylabel=r'$u (m/s)$')
    ax.plot(pressures/1e5,maxJetting,color='black',linestyle=None,marker='o',markersize=3)
    fig.savefig(f"MaxJettingVelocity.png",bbox_inches='tight')
    plt.close(fig)

class simulationRun:
    def __init__(self,filepath,pressure,tfinal,particle_speed,name):
        #Initial Stuff
        self.filepath=filepath
        self.pressure=pressure
        self.tfinal=tfinal
        self.name=name
        self.R0=0.000338
        self.Cl=1462
        self.Rho_l=998
        self.P0=1e5
        self.up=particle_speed
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
        self.KMT,self.RRU=KellerMiksis(float(self.pressure),self.R0)
        self.KMUDO=self.RRU[:,1]
        self.KMUUP=-1*self.KMUDO

            
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

        self.tRayleigh=0.915*self.R0*np.sqrt(self.Rho_l/(float(self.pressure)-self.P0))
        
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
            if(len(Iarr[0]) < 2):
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

    '''
    Returns collapse index based on the point where the differnce in upstream and downstream velocities begins to change to a certain degree
    '''
    def CollapseStartTime(self):
        diffUpstream=np.diff(self.TrimuInterfaceUpstream.flatten())
        UpstreamCollapse=np.min(np.where(diffUpstream > 0.1))
        return UpstreamCollapse
    
    def VectorDiff(self,vec,diff=1e-5):
        VecDiff=np.abs(np.diff(vec))
        CollapseIdx=np.min(np.where(VecDiff > diff))
        return CollapseIdx
    

    def MinVolumeTime(self):
        MinVolumeIdx=np.argmin(self.VolumeArr[:,1])
        tMinVolume=self.tArray[MinVolumeIdx]
        return tMinVolume


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
        fig,ax = self.PlotSpecs(xlabel=r'$t/t_{\tau}$',ylabel=r'$u$')
        tCollapse=self.TrimtArray[self.CollapseStartTime()]
        ax.plot((self.TrimtArray-tCollapse)/self.tRayleigh, self.TrimuInterfaceDownstream,'--',color='black',label='Downstream')
        ax.plot((self.TrimtArray-tCollapse)/self.tRayleigh,self.TrimuInterfaceUpstream,'-',color='black',label='Upstream')
        ax.plot(self.KMT/self.tRayleigh,self.KMUUP,linestyle='dashdot',color='black')
        ax.plot(self.KMT/self.tRayleigh,self.KMUDO,linestyle='dashdot',color='black',label='Keller-Miksis')
        ax.set_xlim(0,1.05)
        fig.tight_layout()
        plt.show()
        fig.savefig(os.path.join(self.filepath, f"VelocityPlot_{self.name}.png"),bbox_inches='tight')
        plt.close(fig)
    def PositionPlot(self):
        fig,ax = self.PlotSpecs(xlabel=r'$t/\tau_{c}$',ylabel=r'$\frac{(y-y_{0})}{R_{0}}$')
        ax.plot(self.TrimtArray/(self.R0/self.Cl),(self.TrimyInterfaceDownstream-self.y0)/self.R0,'--',color='black',label='Downstream')
        ax.plot(self.TrimtArray/(self.R0/self.Cl),(self.TrimyInterfaceUpstream-self.y0)/self.R0,'-',color='black',label='Upstream')
        ax.legend()
        fig.savefig(os.path.join(self.filepath, f"PositionPlot_{self.name}.png"))
        plt.close(fig)
    def VelocityContour(self):
        fig,ax = self.PlotSpecs(xlabel=r'$t/t_{\tau}$',ylabel=r'$\frac{(y-y_{0})}{R_{0}}$')
        tCollapse=self.TrimtArray[self.CollapseStartTime()]
        ax.contourf((self.TrimtArray-tCollapse)/self.tRayleigh,(self.PositionMatrix.squeeze()-self.y0)/self.R0,self.TrimVelocityMatrix/(self.Cl**2),cmap='Greys')
        ax.plot((self.TrimtArray-tCollapse)/self.tRayleigh,(self.TrimyInterfaceDownstream-self.y0)/self.R0,'--',color='black')
        ax.plot((self.TrimtArray-tCollapse)/self.tRayleigh,(self.TrimyInterfaceUpstream-self.y0)/self.R0,'-',color='black')
        ax.set_xlim(0,None)
        fig.tight_layout()
        fig.savefig(os.path.join(self.filepath, f"VelocityContour_{self.name}.png"),bbox_inches='tight')
        plt.close(fig)
    def PressureContour(self):
        fig,ax = self.PlotSpecs(xlabel=r'$\frac{(y-y_{0})}{R_{0}}$',ylabel=r'$t/\tau_{c}$')
        ax.contourf(self.TrimtArray/(self.R0/self.Cl),(self.PositionMatrix.squeeze()-self.y0)/self.R0,self.TrimPressureMatrix,cmap='Greys')
        ax.plot(self.TrimtArray/(self.R0/self.Cl),(self.TrimyInterfaceDownstream-self.y0)/self.R0,'--',color='black')
        ax.plot(self.TrimtArray/(self.R0/self.Cl),(self.TrimyInterfaceUpstream-self.y0)/self.R0,'-',color='black')
        fig.savefig(os.path.join(self.filepath, f"PressureContour_{self.name}.png"))
        plt.close(fig)
    def VolumePlot(self):
        fig,ax = self.PlotSpecs(xlabel=r'$t/t_{\tau}$',ylabel=r'$V_{B}/V_{0}$')
        vB=self.VolumeArr[:,1]/self.VolumeArr[0,1]
        tShift=self.VolumeArr[:,0][self.VectorDiff(vB,1e-6)]
        ax.plot((self.VolumeArr[:,0]-tShift)/self.tRayleigh,vB,'-',color='black')
        if(self.tJetting != 0):
            ax.axvline(x=(self.tJetting-tShift)/self.tRayleigh,linestyle='--',color='black')
        fig.tight_layout()
        ax.set_xlim(0,None)
        fig.savefig(os.path.join(self.filepath, f"VolumePlot_{self.name}.png"),bbox_inches='tight')
        plt.close(fig)
    def SphericityPlot(self):
        fig,ax = self.PlotSpecs(xlabel=r'$t/\tau_{c}$',ylabel=r'$Sphericity$')
        ax.plot(self.SphericityArr[:,0]/(self.R0/self.Cl),self.SphericityArr[:,1],'-',color='black')
        fig.savefig(os.path.join(self.filepath, f"SphericityPlot_{self.name}.png"))
        plt.close(fig)
    def JettingParameters(self):
        if(self.tJetting==0):
            return None,None,None,None
        MinVolumeTime=self.MinVolumeTime()
        JettingVolume=np.max(self.uInterfaceUpstream)
        tMove=self.TrimtArray[self.CollapseStartTime()]
        return self.tJetting, MinVolumeTime,JettingVolume,tMove



    def Run(self):
        self.Interface()
        self.TrimArrays()

        print('Plotting Volume')
        self.VolumePlot()

        print('Plotting Interface Velocity')
        self.VelocityPlot()

        #print('Plotting Interface Position')
        #self.PositionPlot()

        print('Plotting Velocity Contour')
        self.VelocityContour()

        #print('Plotting Pressure Contour')
        #self.PressureContour()

        tJ,mV,jV,tM=self.JettingParameters()
        return tJ,mV,jV,tM

jettingTimes=[]*4
minVolumes=[]*4
maxJetting=[]*4
tMove=[]*4
pressures=[]*4

WS100e6=simulationRun('/home/exy214/Documents/cavitation/data/jetting_ws_2025/WS_100e6_Euler','100e6',5e-6,65.4931,'WS100e6')
tJ,mV,jV,tM=WS100e6.Run()
AppendArrays(jettingTimes,tJ,minVolumes,mV,maxJetting,jV,pressures,100e6,tMove,tM)

WS10e6=simulationRun('/home/exy214/Documents/cavitation/data/jetting_ws_2025/WS_10e6_Euler','10e6',8e-6,6.7532,'WS10e6')
tJ,mV,jV,tM=WS10e6.Run()
AppendArrays(jettingTimes,tJ,minVolumes,mV,maxJetting,jV,pressures,10e6,tMove,tM)

WS1e6=simulationRun('/home/exy214/Documents/cavitation/data/jetting_ws_2025/WS_1e6_Euler','1e6',2e-5,0.616480,'WS1e6')
tJ,mV,jV,tM=WS1e6.Run()
AppendArrays(jettingTimes,tJ,minVolumes,mV,maxJetting,jV,pressures,1e6,tMove,tM)

WS500e5=simulationRun('/home/exy214/Documents/cavitation/data/jetting_ws_2025/WS_500e5_Euler','0.5e6',5e-5,0.274504,'WS500e5')
tJ,mv,jV,tM=WS500e5.Run()
AppendArrays(jettingTimes,tJ,minVolumes,mV,maxJetting,jV,pressures,500e5,tMove,tM)


jettingTimes=np.array(jettingTimes)
minVolumes=np.array(minVolumes)
pressures=np.array(pressures)
maxJetting=np.array(maxJetting)
tMove=np.array(tMove)

JettingImpactPlot(jettingTimes,minVolumes,pressures,0.004662,1462)
ShockVolumePlot(minVolumes,pressures,0.000338,1462)
maxJetVelocity(maxJetting,pressures,1462)

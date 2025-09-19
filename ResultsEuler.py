import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from shockInfo import timeShift
from BubbleCurvature import curvature
from BubbleCurvature import FlipYAxAndJoin


# Some Immutable Stuff

## Plotting Colours
colorBlue = "#0072B2"
colorOrange = "#E69F00"
colorBlack = "#000000"
colorGrey = "#999999"
Colours = [colorBlue, colorOrange, colorBlack, colorGrey]
numSims=4
## Immutable Block Over

plt.rc('text',usetex=True)
plt.rcParams.update({
    'font.size': 14,
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
    'savefig.dpi': 600,
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
    fig, ax = plt.subplots(figsize=(4,3))
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

def JettingImpactPlot(Simulations):
    """ JettingImpactPlot(Simulations) --> Creates a plot of shock pressure ratios vs deltaTjetting/deltaTcollapse against
                                                                   Basically plots shock pressure vs when the re-entrnant jet forms vs when the bubble collapses
    Parameters:
        tJetting: Array of times when re-entrant jet formed
        minVolumes: Array of times when bubble collapsed (reached minimum volume)
        pressures: Array of shock overpressures 
        dX: Distance until proximal side of bubble (most upstream point)
        Cl: Speed of sound in water
    """
    fig,ax = PlotSpecs(xlabel=r'$p_{s}/p_{0}$',ylabel=r'$\Delta T_{j}/ \Delta T_{c}$')
    for idx,sim in enumerate(Simulations):
        tShockImpact=sim.timeShift
        minVolumes=sim.tMinVolume
        tJetting=sim.tJetting
        dParamJet=(tJetting-tShockImpact)/(minVolumes-tShockImpact)
        ax.semilogx(float(sim.pressure)/1e5,dParamJet,color=colorBlue,linestyle=None,marker='o',markersize=3)
    fig.savefig(f"paper-images/JettingImpactPlot.png")
    plt.close(fig)

def ShockVolumePlot(Simulations):
    """ ShockVolumePlot(minVolumes,pressures,R0,Cl) --> Creates a plot of shock pressure ratios vs the time taken for the shock to diffract / minimum volume (collapse)
    Parameters:
        minVolumes: Array of times when bubble collapsed (reached minimum volume)
        pressures: Array of shock overpressures
        R0: Initial bubble radius
        Cl: Speed of sound in water
    """
    fig,ax = PlotSpecs(xlabel=r'$p_{s}/p_{0}$',ylabel=r'$T_{d}/T_{c}$')
    for idx,sim in enumerate(Simulations):
        tShockDiffract=2*sim.R0/sim.Speed
        dParamJet=tShockDiffract/sim.tMinVolume
        tRayleigh=0.915*sim.R0*np.sqrt(998/(float(sim.pressure)-1e5))
        tRayleigh=tShockDiffract/tRayleigh
        ax.semilogx(float(sim.pressure)/1e5,dParamJet,color=Colours[0],linestyle=None,marker='o',markersize=3)
        ax.semilogx(float(sim.pressure)/1e5,tRayleigh,color=Colours[1],linestyle=None,marker='o',markersize=3)
        ax.axhline(y=0.05,linestyle='--',color=colorBlack)
    fig.tight_layout()
    fig.savefig(f"paper-images/ShockVolumePlot.png")
    plt.close(fig)

def VolumePlot(Simulations):
        """
            VolumePlot() --> Plots of volume over time, with the point at which the re-entrant jet forms, if at all
        """
        fig,ax = PlotSpecs(xlabel=r'$\tau$',ylabel=r'$V/V_{0}$')
        for idx,sim in enumerate(Simulations):
            vB=sim.VolumeArr[:,1]/sim.VolumeArr[0,1]
            ax.plot((sim.VolumeArr[:,0]-sim.timeShift)/sim.tRayleigh,vB,'-',color=Colours[idx])
            if(sim.tJetting != 0):
                ax.axvline(x=(sim.tJetting-sim.timeShift)/sim.tRayleigh,linestyle='--',color=Colours[idx])
        fig.tight_layout()
        ax.set_xlim(0,1.5)
        #ax.legend(loc='upper right')
        fig.savefig(f"paper-images/Volume.png")
        plt.close(fig)


class simulationRun:
    """
        class simulationRun --> A post-processing simulation object for each shock induced bubble collapse simulation
    Important Functions
        interface() --> Does all the manual labour regarding interface extraction
        Run() --> Calls class methods which are actually useful
    """
    def __init__(self,filepath,pressure,tfinal,letter,name):
        #Initial Stuff
        self.filepath=filepath
        self.pressure=pressure
        self.tfinal=tfinal
        self.name=name
        self.letter=letter

        # Parameters
        self.R0=5e-4 # Initial Bubble Radius (m)
        self.Cl=1462 # Speed of sound in water (m/s)
        self.Rho_l=998 # Density of water (kg/m^3)
        self.P0=1e5 # Pressure of bubble at t=0
        self.yMeasure=4e-3 # Position that the data begins measuring from 
        
        # Volume Array
        self.VolumeArr=np.loadtxt(os.path.join(filepath, "Volume.txt"))
        
        #For Pressure Matrix
        self.PressureTemplate = os.path.join(filepath, "pressure", "pressure{:04d}.txt")
        self.NumRuns = len(self.VolumeArr)
        self.PressureList = [np.loadtxt(self.PressureTemplate.format(i)) for i in range(self.NumRuns)]
        self.PressureMatrix=np.stack(self.PressureList,axis=1)

        #For Velocity Matrix
        self.VelocityTemplate = os.path.join(filepath, "velocity", "velocity{:04d}.txt")
        self.NumRuns = len(self.VolumeArr)
        self.VelocityList = [np.loadtxt(self.VelocityTemplate.format(i)) for i in range(self.NumRuns)]
        self.VelocityMatrix=np.stack(self.VelocityList,axis=1)

        #For Diffuse Volume Fraction Matrix
        self.AlphaTemplate = os.path.join(filepath, "alpha", "alpha{:04d}.txt")
        self.NumRuns = len(self.VolumeArr)
        self.AlphaList = [np.loadtxt(self.AlphaTemplate.format(i)) for i in range(self.NumRuns)]
        self.AlphaMatrix=np.stack(self.AlphaList,axis=1)

        #For Interface Positional Matrix

        self.InterfaceTemplate = os.path.join(filepath, "interface", "Interface{:04d}.curve")
        self.NumRuns = len(self.VolumeArr)
        self.InterfaceList = [np.loadtxt(self.InterfaceTemplate.format(i)) for i in range(self.NumRuns)]
        self.InterfaceMatrix=np.stack(self.AlphaList,axis=1)

        #For Position Matrix
        self.PositionMatrix=np.loadtxt(os.path.join(filepath,"pressure/position0000.txt"))+self.yMeasure

        # Time Matrix
        self.startTime=4e-6
        self.tArray=np.linspace(self.startTime,tfinal,self.NumRuns)
        self.tArray=np.insert(self.tArray,0,0)

        # Interface Values
        self.yInterfaceUpstream = np.zeros(((self.NumRuns,1)))
        self.uInterfaceUpstream = np.zeros(((self.NumRuns,1)))
        self.pInterfaceUpstream = np.zeros(((self.NumRuns,1)))
        self.aInterfaceUpstream = np.zeros(((self.NumRuns,1)))
        self.yInterfaceDownstream = np.zeros(((self.NumRuns,1)))
        self.uInterfaceDownstream = np.zeros(((self.NumRuns,1)))
        self.pInterfaceDownstream = np.zeros(((self.NumRuns,1)))
        self.aInterfaceDownstream = np.zeros(((self.NumRuns,1)))
        self.Curvature=np.zeros(((self.NumRuns,1)))


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

        self.tRayleigh=0.915*self.R0*np.sqrt(self.Rho_l/(float(self.pressure)-self.P0)) #Time taken for bubble collape according to Rayleigh model
        
        # Values
        self.y0=8e-3 # Assuming the start position of bubble --> Actually could be wrong unfortunaetly 
        self.dx=(8e-3)-(5e-4) #Distance from bottom boundary to bubble

        # Water params
        self.gmaW=2.955
        self.pinfW=7.22e8
        self.rhoW=998
        self.v1W=0

        # Jetting Parameters
        self.tJetting=0 #Time at which jetting occurs (initialising)
        self.tMinVolume=0 # Time at which minimum volume is reached (collape time)



        self.writePath='/home/exy214/Documents/cavitation/data/jetting_ws_2025/paper-images/' # Self explanatory
        self.Speed,self.timeShift=timeShift(float(self.pressure),self.P0,self.gmaW,self.pinfW,self.rhoW,self.v1W,self.dx) # Time taken for shockwave to hit bubble
    def Interface(self):
        """
            Interface(self) --> Defines interfacial properties of position, velocity and volume fraction along centerline at most upstream and downstream points
        """
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
    
    def TrimArrays(self):
        """
            TrimArrays(self) --> Trim arrays based upon the number of zeros trailing the upstream velocity interface vector
        """
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

    def MinVolumeTime(self):
        """
            MinVolumeTime() --> Index at which the minimum volume is reached
        """
        MinVolumeIdx=np.argmin(self.VolumeArr[:,1])
        tMinVolume=self.tArray[MinVolumeIdx]
        return tMinVolume
    
    def BubbleCurvature(self):
        for i in range(self.NumRuns):
            Interface=self.InterfaceMatrix[:,i]
            idxMUP=np.argmin(Interface[:,0])  # Left most point
            Interface=Interface[idxMUP:idxMUP+int(0.05*len(Interface))] #Find next 10 points and add to array (nearest most upstream point)
            x=Interface[:,0]
            y=Interface[:,1]
            x,y=FlipYAxAndJoin(x,y)

            x=savgol_filter(x,window_length=int(0.15*len(x)),polyorder=2)
            y=savgol_filter(y,window_length=int(0.15*len(y)),polyorder=2)

            # Bubble curvature
            bubbleCur=curvature(x,y)
            self.Curvature[i]=bubbleCur[0]

    def PlotSpecs(self,xlabel,ylabel):
        """
            Defined Above
        """
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
        """
            VelocityPlot() --> Plots the interfacial velocity for most upstream and downstream points
        """
        fig,ax = self.PlotSpecs(xlabel=r'$\tau$',ylabel=r'$u$')
        ax.plot((self.TrimtArray-self.timeShift)/self.tRayleigh, self.TrimuInterfaceDownstream,'--',color='black',label='Downstream')
        ax.plot((self.TrimtArray-self.timeShift)/self.tRayleigh,self.TrimuInterfaceUpstream,'-',color='black',label='Upstream')
        if(self.tJetting!=0):
            ax.set_xlim(0,(self.tJetting-self.timeShift)/self.tRayleigh)
        else:
            ax.set_xlim(0,1.75)
        ax.text(0.025,0.95, self.letter, transform=ax.transAxes, fontsize=12, va='top', ha='left')
        fig.tight_layout()
        fig.savefig(os.path.join(self.writePath, f"VelocityPlot_{self.name}.png"))
        plt.close(fig)

    def VelocityContour(self):
        """
            VelocityContour() --> Plots the velocity contour for most upstream and downstream points overlayed with position plot
        """
        fig,ax = self.PlotSpecs(xlabel=r'$\tau$',ylabel=r'$\frac{(z-z_{0})}{R_{0}}$')
        cntr = ax.contourf((self.TrimtArray-self.timeShift)/self.tRayleigh,(self.PositionMatrix.squeeze()-self.y0)/self.R0,self.TrimVelocityMatrix,cmap='RdYlBu_r',levels=256)
        ax.plot((self.TrimtArray-self.timeShift)/self.tRayleigh,(self.TrimyInterfaceDownstream-self.y0)/self.R0,'--',color='black')
        ax.plot((self.TrimtArray-self.timeShift)/self.tRayleigh,(self.TrimyInterfaceUpstream-self.y0)/self.R0,'-',color='black')
        cbar=fig.colorbar(cntr) 
        if(self.tJetting!=0):
            ax.set_xlim(0,(self.tJetting-self.timeShift)/self.tRayleigh)
        else:
            ax.set_xlim(0,1.75)
        ax.set_ylim(-2,2)
        ax.text(0.025,0.95, self.letter, transform=ax.transAxes, fontsize=14, va='top', ha='left')
        fig.tight_layout()
        fig.savefig(os.path.join(self.writePath, f"VelocityContour_{self.name}.png"))
        plt.close(fig)
    def PressureContour(self):
        """
            PressureContour() --> Plots the pressure contour for the most upstream and downstream points overlayed with position plot
        """
        fig,ax = self.PlotSpecs(xlabel=r'$\tau$',ylabel=r'$\frac{(z-z_{0})}{R_{0}}$')
        cntr = ax.contourf((self.TrimtArray-self.timeShift)/self.tRayleigh,(self.PositionMatrix.squeeze()-self.y0)/self.R0,self.TrimPressureMatrix,cmap='RdYlBu_r',levels=256)
        ax.plot((self.TrimtArray-self.timeShift)/self.tRayleigh,(self.TrimyInterfaceDownstream-self.y0)/self.R0,'--',color='black')
        ax.plot((self.TrimtArray-self.timeShift)/self.tRayleigh,(self.TrimyInterfaceUpstream-self.y0)/self.R0,'-',color='black')
        ax.text(0.025,0.95, self.letter, transform=ax.transAxes, fontsize=12, va='top', ha='left')
        cbar=fig.colorbar(cntr)
        if(self.tJetting!=0):
            ax.set_xlim(0,(self.tJetting-self.timeShift)/self.tRayleigh)
        else:
            ax.set_xlim(0,1.75)
        ax.set_ylim(-2,2)
        fig.tight_layout()
        fig.savefig(os.path.join(self.writePath, f"PressureContour_{self.name}.png"))
        plt.close(fig)
    
    def JettingParameters(self):
        """
            JettingParameters() --> Returns paramters of minimum volume
        Return:
            tJetting: Time at which re-entrant jet forms
            MinVolumeTime: Time at which the minimum bubble volume is reached
        """
        if(self.tJetting==0):
            return None,None,None
        self.tMinVolume=self.MinVolumeTime()

    def Run(self):
        """
            Run() --> Runs class methods
        """
        self.Interface()
        self.TrimArrays()

        print('Plotting Interface Velocity')
        self.VelocityPlot()

        print('Plotting Velocity Contour')
        self.VelocityContour()

        print('Plotting Pressure Contour')
        self.PressureContour()

        print('Getting Minimum Volume time')
        self.JettingParameters()




## Simulation Objects
Simulations=[]*numSims #Empty simulations array

#Running objects and appending to the arrays
PostProc=simulationRun('/home/exy214/Documents/cavitation/data/jetting_ws_2025/PostProc','100e6',7.5e-6,'(a)','PostProc')
PostProc.Run()
Simulations.append(PostProc)

# Plotting
JettingImpactPlot(Simulations)
ShockVolumePlot(Simulations)
VolumePlot(Simulations)

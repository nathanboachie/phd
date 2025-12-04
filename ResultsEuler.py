import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
from scipy.signal import resample
from scipy.signal import savgol_filter
from shockInfo import timeShift
from BubbleCurvature import curvature
from BubbleCurvature import FlipYAxAndJoin
from PIL import Image


# Some Immutable Stuff
## Plotting Colours
colorBlue = "#0072B2"
colorOrange = "#E69F00"
colorBlack = "#000000"
colorGrey = "#999999"
colorPurple = "#CC79A7"
Colours = [colorBlue, colorOrange, colorBlack, colorGrey]
numSims=4
## Immutable Block Over


def Jumps(arr):
    for i in range(len(arr)-1):
        if np.abs(arr[i+1]-arr[i])>0.5:
            print(f'Jump happens at index {i}')

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
    fig,ax = PlotSpecs(xlabel=r'$p_{s}/p_{0}$',ylabel=r'$log(t_{j}/t_{c})$')
    Params=[]
    Pressures=[]
    for idx,sim in enumerate(Simulations):
        if(sim.tJetting == 0):
            continue
        tShockImpact=sim.timeShift
        minVolumes=sim.tMinVolume
        tJetting=sim.tJetting
        dParamJet=(tJetting-tShockImpact)/(minVolumes-tShockImpact)
        ax.semilogx((float(sim.pressure)/1e5),(dParamJet),color=colorBlue,linestyle='-',marker='o',markersize=3)
        print(f'For Simulation {sim.name} is: {float(sim.pressure)/1e5},{dParamJet}')
    ax.axhline(y=1.0,color=colorBlue,linestyle='dashed')
    ax.set_xlim([4.5,1200])
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
    fig,ax = PlotSpecs(xlabel=r'$p_{s}/p_{0}$',ylabel=r'$t_{d}/t_{c}$')

    rcPressure=np.logspace(np.log10(0.4e6), np.log10(100e6), 200)
    R0=5e-4
    RhoW=998
    P0=1e5
    rcSpeed=np.zeros(len(rcPressure))
    dx=(8e-3)-(5e-4) #Distance from bottom boundary to bubble
    for i in range(len(rcPressure)):
        rcSpeed[i],ts=timeShift(rcPressure[i],P0,2.955,7.22e8,RhoW,0,dx)
    rcTd=(2*R0/rcSpeed)
    rcTc=(0.915*R0*np.sqrt(RhoW/((rcPressure-P0))))
    rcdParamJet=rcTd/rcTc

    for idx,sim in enumerate(Simulations):
        tShockDiffract=2*sim.R0/sim.Speed 
        dParamJet=tShockDiffract/(sim.tMinVolume-sim.timeShift)
        tRayleigh=0.915*sim.R0*np.sqrt(998/(float(sim.pressure)-1e5))

        if (sim.tJetting !=0):
            ax.semilogx(float(sim.pressure)/1e5,dParamJet,color=Colours[0],linestyle='-',marker='o',markersize=3)
            ax.semilogx(rcPressure/1e5,rcdParamJet,color=Colours[2],linestyle='-',marker=None,markersize=3)
        else:
            ax.semilogx(float(sim.pressure)/1e5,dParamJet,color=colorPurple,linestyle=None,marker='o',markersize=3)
        ax.axhline(y=0.05,linestyle='--',color=colorBlack)
    fig.tight_layout()
    fig.savefig(f"paper-images/ShockVolumePlot.png")
    plt.close(fig)

def ProlongationPlot(Simulations):
    """ Prolongation(minVolumes,pressures,R0,Cl) --> Creates a plot of prolongation factor
    Parameters:
        minVolumes: Array of times when bubble collapsed (reached minimum volume)
        pressures: Array of shock overpressures
        R0: Initial bubble radius
        Cl: Speed of sound in water
    """
    fig,ax = PlotSpecs(xlabel=r'$p_{s}/p_{0}$',ylabel=r'$t_{c}/\tau$')
    for idx,sim in enumerate(Simulations):
        tRayleigh=0.915*sim.R0*np.sqrt(998/(float(sim.pressure)-1e5))
        prg=(sim.tMinVolume-sim.timeShift)/tRayleigh
        if (sim.tJetting !=0):
            ax.semilogx(float(sim.pressure)/1e5,prg,color=Colours[0],linestyle='-',marker='o',markersize=3)
        else:
            ax.semilogx(float(sim.pressure)/1e5,prg,color=colorPurple,linestyle=None,marker='o',markersize=3)
    fig.tight_layout()
    fig.savefig(f"paper-images/ProlongationPlot.png")
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
                ax.axvline(x=(sim.tMinVolume-sim.timeShift)/sim.tRayleigh,linestyle='-',color=Colours[idx])
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
    def __init__(self,filepath,writePath,R0,y0,yMeasure,pressure,startTime,tfinal,tVisJetting,sigma,velShock,EOS,letter,name):
        #Initial Stuff
        self.filepath=filepath
        self.pressure=pressure
        self.tfinal=tfinal
        self.startTime=startTime
        self.tVisJetting=tVisJetting #Time when bubbl 'jets' visually
        self.name=name
        self.letter=letter

        # Parameters
        self.R0=R0 # Initial Bubble Radius (m)
        self.Cl=1462 # Speed of sound in water (m/s)
        self.Rho_l=998 # Density of water (kg/m^3)
        self.sigma=sigma
        self.P0 = 1e5 + 2*(sigma)/(self.R0)
        self.yMeasure=yMeasure # Position that the data begins measuring from 
        self.jetTol=0.025 # Percentage of a value such that jetting can be said to have occured
        
        # Volume Array
        self.VolumeArr=np.loadtxt(os.path.join(filepath, "Volume.txt"))
        
        #For Pressure Matrix
        self.PressureTemplate = os.path.join(filepath, "pressure", "pressure{:04d}.txt")
        self.NumRuns = len(self.VolumeArr)
        self.PressureList = [np.loadtxt(self.PressureTemplate.format(i)) for i in range(self.NumRuns)]
        self.PressureMatrix=np.array(self.PressureList,dtype=object)

        #For Velocity Matrix
        self.VelocityTemplate = os.path.join(filepath, "velocity", "velocity{:04d}.txt")
        self.NumRuns = len(self.VolumeArr)
        self.VelocityList = [np.loadtxt(self.VelocityTemplate.format(i)) for i in range(self.NumRuns)]
        self.VelocityMatrix=np.array(self.VelocityList,dtype=object)

        #For Diffuse Volume Fraction Matrix
        self.AlphaTemplate = os.path.join(filepath, "alpha", "alpha{:04d}.txt")
        self.NumRuns = len(self.VolumeArr)
        self.AlphaList = [np.loadtxt(self.AlphaTemplate.format(i)) for i in range(self.NumRuns)]
        self.AlphaMatrix=np.array(self.AlphaList,dtype=object)

        #For Interface Positional Matrix
        self.InterfaceTemplate = os.path.join(filepath, "interface", "Interface{:04d}.curve")
        self.NumRuns = len(self.VolumeArr)
        self.InterfaceList = [np.loadtxt(self.InterfaceTemplate.format(i)) for i in range(self.NumRuns)]
        self.InterfaceMatrix=np.array(self.InterfaceList,dtype=object)

        #For Position Matrix
        self.posMeasure=self.yMeasure
        self.PositionTemplate = os.path.join(filepath, "pressure", "pressureposition{:04d}.txt")
        self.NumRuns = len(self.VolumeArr)
        self.PositionList = [np.loadtxt(self.PositionTemplate.format(i))+self.posMeasure for i in range(self.NumRuns)]
        self.PositionMatrix=np.array(self.PositionList,dtype=object)

        # Time Matrix
        self.tArray=self.VolumeArr[:,0]

        # Interface Values
        self.yInterfaceUpstream = np.zeros(((self.NumRuns,1)))
        self.uInterfaceUpstream = np.zeros(((self.NumRuns,1)))
        self.pInterfaceUpstream = np.zeros(((self.NumRuns,1)))
        self.aInterfaceUpstream = np.zeros(((self.NumRuns,1)))
        self.yInterfaceDownstream = np.zeros(((self.NumRuns,1)))
        self.uInterfaceDownstream = np.zeros(((self.NumRuns,1)))
        self.pInterfaceDownstream = np.zeros(((self.NumRuns,1)))
        self.aInterfaceDownstream = np.zeros(((self.NumRuns,1)))
        self.dudtInterfaceDownstream= np.zeros(((self.NumRuns,1)))
        self.dudtInterfaceUpstream= np.zeros(((self.NumRuns,1)))
        self.Curvature=np.zeros(((self.NumRuns,1)))

        # Water Interface
        self.WaterpInterfaceUpstream=np.zeros(((self.NumRuns,1)))
        self.WaterpInterfaceDownstream=np.zeros(((self.NumRuns,1)))


        #Sampling/Interpolating
        self.nSamples=200
        self.InterpPresM=np.zeros((self.NumRuns,self.nSamples))
        self.InterpVelM=np.zeros((self.NumRuns,self.nSamples))
        self.InterpPosM=np.zeros((self.NumRuns,self.nSamples))
        self.InterpTimeM=np.zeros((self.NumRuns,self.nSamples))

        self.tRayleigh=0.915*self.R0*np.sqrt(self.Rho_l/(float(self.pressure)-self.P0)) #Time taken for bubble collape according to Rayleigh model
        
        # Values
        self.y0=y0 # Assuming the start position of bubble --> Actually could be wrong unfortunaetly 
        self.dx=(self.y0)-(self.R0) #Distance from bottom boundary to bubble

        # Water params
        self.gmaW=EOS[0]
        self.pinfW=EOS[1]
        self.rhoW=EOS[2]
        self.v1W=velShock

        # Jetting Parameters
        self.tJetting=0 #Time at which jetting occurs (initialising)
        self.tMinVolume=0 # Time at which minimum volume is reached (collape time)
        # Trimmed Interface Values
        '''
        Why? This is a trim matrix that ends at the point of jetting
        The data looks weird including information past jetting
        '''
        self.MaxCenterlinePressure = np.zeros(((self.NumRuns,1)))
        self.TrimyInterfaceUpstream = np.zeros(((self.NumRuns,1)))
        self.TrimuInterfaceUpstream = np.zeros(((self.NumRuns,1)))
        self.TrimpInterfaceUpstream = np.zeros(((self.NumRuns,1)))
        self.TrimaInterfaceUpstream = np.zeros(((self.NumRuns,1)))
        self.TrimyInterfaceDownstream = np.zeros(((self.NumRuns,1)))
        self.TrimuInterfaceDownstream = np.zeros(((self.NumRuns,1)))
        self.TrimpInterfaceDownstream = np.zeros(((self.NumRuns,1)))
        self.TrimaInterfaceDownstream = np.zeros(((self.NumRuns,1)))
        self.TrimdudtInterfaceDownstream = np.zeros(((self.NumRuns,1)))
        self.TrimdudtInterfaceUpstream = np.zeros(((self.NumRuns,1)))
        self.TrimtArray = np.zeros(((self.NumRuns,1)))
        self.TrimPressureMatrix=np.zeros((self.NumRuns,len(self.PositionMatrix)))
        self.TrimVelocityMatrix=np.zeros((self.NumRuns,len(self.PositionMatrix)))
        self.TrimAlphaMatrix=np.zeros((self.NumRuns,len(self.PositionMatrix)))

        # Water Trim
        self.WaterTrimpInterfaceDownstream = np.zeros(((self.NumRuns,1)))
        self.WaterTrimpInterfaceUpstream = np.zeros(((self.NumRuns,1)))

        



        self.writePath=writePath # Self explanatory
        self.Speed,self.timeShift=timeShift(float(self.pressure),self.P0,self.gmaW,self.pinfW,self.rhoW,self.v1W,self.dx) # Time taken for shockwave to hit bubble
        print(f'For {self.name} Shock reaches the bubble at {self.timeShift} s')

    def maxPressure(self):
        for i in range(self.NumRuns):
            PressureArr=self.PressureMatrix[i]
            maxP=np.max(PressureArr)
            self.MaxCenterlinePressure[i]=maxP

    def Interface(self):
        """
            Interface(self) --> Defines interfacial properties of position, velocity and volume fraction along centerline at most upstream and downstream points
        """
        for i in range(self.NumRuns):
            alphaI=self.AlphaList[i]
            pressureI=self.PressureList[i]
            velocityI=self.VelocityList[i]
            positionI=self.PositionList[i]
            Iarr=np.where(alphaI>0.95)[0]
            if(len(Iarr) < 2):
                print(f'No more than one cell on the centerline')
                break
            else:
                YGasPosition=positionI[Iarr]
                YMax=np.max(YGasPosition)
                I=int(np.where(positionI==YMax)[0][0])
                self.yInterfaceDownstream[i]=positionI[I]
                self.uInterfaceDownstream[i]=velocityI[I]
                self.pInterfaceDownstream[i]=pressureI[I]
                self.aInterfaceDownstream[i]=alphaI[I]

                YMin=np.min(YGasPosition)
                I=int(np.where(positionI==YMin)[0][0])
                self.yInterfaceUpstream[i]=positionI[I]
                self.uInterfaceUpstream[i]=velocityI[I]
                self.pInterfaceUpstream[i]=pressureI[I]
                self.aInterfaceUpstream[i]=alphaI[I]
        self.dudtInterfaceDownstream=np.gradient(np.asarray(self.uInterfaceDownstream).squeeze(),np.asarray(self.tArray).squeeze())
        self.dudtInterfaceUpstream=np.gradient(np.asarray(self.uInterfaceUpstream).squeeze(),np.asarray(self.tArray).squeeze())
    
    def WaterInterface(self):
        """
            Interface(self) --> Defines interfacial properties of position, velocity and volume fraction along centerline at most upstream and downstream points
        """
        for i in range(self.NumRuns):
            alphaI=self.AlphaList[i]
            pressureI=self.PressureList[i]
            velocityI=self.VelocityList[i]
            positionI=self.PositionList[i]
            Iarr=np.where(alphaI>1e-9)[0]
            CellShift=2
            if(len(Iarr) < 2):
                break
            else:
               Id=Iarr[np.argmax(positionI[Iarr])]
               Iu=Iarr[np.argmin(positionI[Iarr])]

               Idshift=min(Id+CellShift,len(positionI)-1)
               Iushift=max(Iu-CellShift,0)

               self.WaterpInterfaceDownstream[i]=pressureI[Idshift]
               self.WaterpInterfaceUpstream[i]=pressureI[Iushift]
    def Interpolate(self):

        xmin=min(pos[0] for pos in self.PositionMatrix)
        xmax=max(pos[-1] for pos in self.PositionMatrix)

        spatialGrid=np.linspace(xmin,xmax,self.nSamples)

        pressureI=[np.interp(spatialGrid,self.PositionMatrix[i],self.PressureMatrix[i]) for i in range(self.NumRuns)]
        velocityI=[np.interp(spatialGrid,self.PositionMatrix[i],self.VelocityMatrix[i]) for i in range(self.NumRuns)]
        positionI=[np.interp(spatialGrid,self.PositionMatrix[i],self.PositionMatrix[i]) for i in range(self.NumRuns)]

        self.InterpPresM=np.vstack(pressureI)
        self.InterpVelM=np.vstack(velocityI)
        self.InterpPosM=np.vstack(positionI)
        self.InterpTimeM = np.tile(self.tArray[:, np.newaxis], (1, 200))
    def TrimArrays(self):
        """
            TrimArrays(self) --> Trim arrays based upon the number of zeros trailing the upstream velocity interface vector
        """
        self.TrimuInterfaceUpstream=np.trim_zeros(self.uInterfaceUpstream,'b')
        self.TrimuInterfaceDownstream=np.trim_zeros(self.uInterfaceDownstream,'b')
        self.TrimyInterfaceUpstream=self.yInterfaceUpstream[0:len(self.TrimuInterfaceDownstream)]
        self.TrimyInterfaceDownstream=self.yInterfaceDownstream[0:len(self.TrimuInterfaceDownstream)]
        self.TrimpInterfaceDownstream=self.pInterfaceDownstream[0:len(self.TrimuInterfaceDownstream)]
        self.TrimpInterfaceUpstream=self.pInterfaceUpstream[0:len(self.TrimuInterfaceDownstream)]
        self.TrimaInterfaceUpstream=self.aInterfaceUpstream[0:len(self.TrimuInterfaceDownstream)]
        self.TrimaInterfaceDownstream=self.aInterfaceDownstream[0:len(self.TrimuInterfaceDownstream)]
        self.TrimtArray=self.tArray[0:len(self.TrimuInterfaceDownstream)]

        self.TrimdudtInterfaceDownstream=self.dudtInterfaceDownstream[0:len(self.TrimuInterfaceDownstream)]
        self.TrimdudtInterfaceUpstream=self.dudtInterfaceUpstream[0:len(self.TrimuInterfaceDownstream)]
        
        self.WaterTrimpInterfaceUpstream=self.WaterpInterfaceUpstream[0:len(self.TrimuInterfaceDownstream)]
        self.WaterTrimpInterfaceDownstream=self.WaterpInterfaceDownstream[0:len(self.TrimuInterfaceDownstream)]

        self.TrimInterpPresM=self.InterpPresM[:len(self.TrimuInterfaceDownstream),:]
        self.TrimInterpPosM=self.InterpPosM[:len(self.TrimuInterfaceDownstream),:]
        self.TrimInterpVelM=self.InterpVelM[:len(self.TrimuInterfaceDownstream),:]
        self.TrimInterpTimeM=self.InterpTimeM[:len(self.TrimuInterfaceDownstream),:]

    def MinVolumeTime(self):
        """
            MinVolumeTime() --> Index at which the minimum volume is reached
        """
        MinVolumeIdx=np.argmin(self.VolumeArr[:,1])
        tMinVolume=self.tArray[MinVolumeIdx]

        return tMinVolume
        
    def JettingParameters(self):
        """
            JettingParameters() --> Returns paramters of minimum volume
        Return:
            tJetting: Time at which re-entrant jet forms
            MinVolumeTime: Time at which the minimum bubble volume is reached
        """
        self.tDiffract=self.timeShift+2*self.R0/self.Speed
        print(f'For {self.name} the time which the shock diffracts around the bubble, to the most downstream point is: {self.tDiffract}')
        self.tMinVolume=self.MinVolumeTime()
        print(f'For {self.name} the time when bubble collapses is: {self.tMinVolume} s')
        print(f'For {self.name} collapse time according to the Rayleigh model is: {self.tRayleigh + (self.startTime)} s')
        
        BokmandudtDownstream=self.dudtInterfaceDownstream/(self.P0/self.rhoW*self.R0) 
        BokmandudtUpstream=self.dudtInterfaceUpstream/(self.P0/self.rhoW*self.R0) 

        #self.tAccShift=self.TrimtArray[np.argwhere(BokmandudtUpstream)[0]]
        #print(f'Point at which bubbles upstream acceleration becomes negative is: {self.tAccShift}')

        InterfaceDx=np.abs((self.TrimyInterfaceUpstream-self.TrimyInterfaceDownstream))
        JettingIdx=np.argwhere(InterfaceDx<self.jetTol*self.R0)
        if(JettingIdx.size<=0 and self.tVisJetting==0):
            self.tJetting=0
            print(f'For {self.name} the bubble doesnt jet')
        else:
            self.tJetting=self.tVisJetting
            print(f'For {self.name} the bubble jets visually at {self.tJetting} s')
                
        # Maximum Pressure #
        self.MaxPressure=np.max(self.TrimInterpPresM)
    
    def BubbleCurvature(self):
        for i in range(self.NumRuns):
            Interface=self.InterfaceMatrix[i]
            idxMUP=np.argmin(Interface[:,0])  # Left most point
            LenInterfaceArr=51
            Interface=Interface[idxMUP:idxMUP+LenInterfaceArr] #Find next 10 points and add to array (nearest most upstream point)
            x=Interface[:,0]
            y=Interface[:,1]
            x,y=FlipYAxAndJoin(x,y)
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

    def MaxPressurePlot(self):
        """
        """
        fig,ax = self.PlotSpecs(xlabel=r'$t/\tau$',ylabel=r'$P_{max}/P_{s}$')
        ax.semilogy((self.tArray-self.timeShift)/self.tRayleigh,self.MaxCenterlinePressure/float(self.pressure),'-',color='black',label='Downstream')
        if(self.tJetting!=0):
            ax.axvline(x=(self.tJetting-self.timeShift)/self.tRayleigh,linestyle='-',color='black')
        ax.set_xlim(0,None)
        ax.text(0.025,0.95, self.letter, transform=ax.transAxes, fontsize=12, va='top', ha='left')
        fig.tight_layout()
        fig.savefig(os.path.join(self.writePath, f"MaxPressurePlot_{self.name}.png"))
        plt.close(fig)
    def VolumePlot(self):
        """
            VelocityPlot() --> Plots the interfacial velocity for most upstream and downstream points
        """
        fig,ax = self.PlotSpecs(xlabel=r'$t$',ylabel=r'$v$')
        ax.plot((self.tArray-self.timeShift)/self.tRayleigh,self.VolumeArr[:,1]/self.VolumeArr[0,1],'-',color='black',label='Downstream')
        ax.axvline(x=(self.tMinVolume-self.timeShift)/self.tRayleigh,linestyle='--',color='black')
        if(self.tJetting!=0):
            ax.axvline(x=(self.tJetting-self.timeShift)/self.tRayleigh,linestyle='--',color='blue')
        ax.set_xlim(0,None)
        ax.text(0.025,0.95, self.letter, transform=ax.transAxes, fontsize=12, va='top', ha='left')
        fig.tight_layout()
        fig.savefig(os.path.join(self.writePath, f"VolumePlot_{self.name}.png"))
        plt.close(fig)

    def PositionPlot(self):
        """
            PositionPlot() --> Plots the interfacial position for most upstream and downstream points
        """
        fig,ax = self.PlotSpecs(xlabel=r'$t/\tau$',ylabel=r'$(y-y0)/r0$')
        ax.plot((self.TrimtArray-self.timeShift)/self.tRayleigh, (self.TrimyInterfaceDownstream-self.y0)/self.R0,'--',color='black',label='Downstream')
        ax.plot((self.TrimtArray-self.timeShift)/self.tRayleigh,(self.TrimyInterfaceUpstream-self.y0)/self.R0,'-',color='black',label='Upstream')
        ax.plot((self.TrimtArray-self.timeShift)/self.tRayleigh,(self.TrimyInterfaceUpstream-self.TrimyInterfaceDownstream),'-',color='blue',label='Upstream')
        ax.set_xlim(0,None)
        ax.text(0.025,0.95, self.letter, transform=ax.transAxes, fontsize=12, va='top', ha='left')
        fig.tight_layout()
        fig.savefig(os.path.join(self.writePath, f"PositionPlot_{self.name}.png"))
        plt.close(fig)
    
    def PressurePlot(self):
        """
            PressurePlot() --> Plots the interfacial position for most upstream and downstream points
        """
        fig,ax = self.PlotSpecs(xlabel=r'$t/t_{c}$',ylabel=r'$(P_{u}-P_{d})/P_{s}$')
        ax.plot((self.TrimtArray-self.timeShift)/self.tRayleigh, self.TrimpInterfaceDownstream/float(self.pressure),'--',color='black',label='Downstream')
        ax.plot((self.TrimtArray-self.timeShift)/self.tRayleigh,self.TrimpInterfaceUpstream/float(self.pressure),'-',color='black',label='Upstream')
        #ax.plot((self.TrimtArray-self.timeShift)/(self.tMinVolume-self.timeShift),(self.WaterTrimpInterfaceUpstream-self.WaterTrimpInterfaceDownstream)/float(self.pressure),'-',color='black',label='Upstream')
        ax.set_xlim(0,None)
        fig.tight_layout()
        fig.savefig(os.path.join(self.writePath, f"PressurePlot_{self.name}.png"))
        plt.close(fig)

    def DiffPositionPlot(self):
        """
            PositionPlot() --> Plots the interfacial position for most upstream and downstream points
        """
        fig,ax = self.PlotSpecs(xlabel=r'$t/\tau$',ylabel=r'$(y-y0)/r0$')
        InterfaceDx=np.abs((self.TrimyInterfaceUpstream-self.TrimyInterfaceDownstream))
        ax.plot((self.TrimtArray-self.timeShift)/self.tRayleigh,InterfaceDx,'-',color='black',label='Upstream')
        if(self.tJetting!=0):
            ax.axvline(x=(self.tJetting-self.timeShift)/self.tRayleigh,linestyle='-',color='black')
            ax.axvline(x=(self.tVisJetting-self.timeShift)/self.tRayleigh,linestyle='--',color='blue')
        ax.set_xlim(0,None)
        ax.text(0.025,0.95, self.letter, transform=ax.transAxes, fontsize=12, va='top', ha='left')
        fig.tight_layout()
        fig.savefig(os.path.join(self.writePath, f"DiffPositionPlot_{self.name}.png"))
        plt.close(fig)


    def VelocityPlot(self):
        """
            VelocityPlot() --> Plots the interfacial velocity for most upstream and downstream points
        """
        fig,ax = self.PlotSpecs(xlabel=r'$t/\tau$',ylabel=r'$u$')
        ax.plot((self.TrimtArray-self.timeShift)/self.tRayleigh, self.TrimuInterfaceDownstream,'--',color='black',label='Downstream')
        ax.plot((self.TrimtArray-self.timeShift)/self.tRayleigh,self.TrimuInterfaceUpstream,'-',color='black',label='Upstream')
        ax.set_xlim(0,None)
        ax.text(0.025,0.95, self.letter, transform=ax.transAxes, fontsize=12, va='top', ha='left')
        fig.tight_layout()
        fig.savefig(os.path.join(self.writePath, f"VelocityPlot_{self.name}.png"))
        plt.close(fig)
    
    def AccelerationPlot(self):
        """
            AccelerationPlot() --> Plots the interfacial acceleration for most upstream and downstream points
        """
        fig,ax = self.PlotSpecs(xlabel=r'$t/t_{c}$',ylabel=r'$a_{j}/(p_{b,0}/\rho r_{0})$')
        #ax.plot((self.TrimtArray-self.timeShift)/(self.tMinVolume-self.timeShift),self.TrimdudtInterfaceDownstream/(self.P0/self.rhoW*self.R0),'--',color='black',label='Downstream')
        ax.plot(self.TrimtArray,self.TrimdudtInterfaceUpstream/(self.P0/self.rhoW*self.R0),'-',color='black',label='Upstream')
        ax.set_xlim(0,None)
        ax.text(0.025,0.95, self.letter, transform=ax.transAxes, fontsize=12, va='top', ha='left')
        fig.tight_layout()
        fig.savefig(os.path.join(self.writePath, f"AccelerationPlot_{self.name}.png"))
        plt.close(fig)

    def VelocityContour(self):
        """
            VelocityContour() --> Plots the velocity contour for most upstream and downstream points overlayed with position plot
        """
        fig,ax = self.PlotSpecs(xlabel=r'$\tau$',ylabel=r'$\frac{(z-z_{0})}{R_{0}}$')
        cntr = ax.contourf((self.TrimInterpTimeM-self.timeShift)/self.tRayleigh,(self.TrimInterpPosM-self.y0)/self.R0,self.TrimInterpVelM,cmap='RdYlBu_r',levels=256)
        ax.plot((self.TrimtArray-self.timeShift)/self.tRayleigh,(self.TrimyInterfaceDownstream-self.y0)/self.R0,'--',color='black')
        ax.plot((self.TrimtArray-self.timeShift)/self.tRayleigh,(self.TrimyInterfaceUpstream-self.y0)/self.R0,'-',color='black')
        if(self.tJetting!=0):
            ax.set_xlim(0,(self.tVisJetting-self.timeShift)/self.tRayleigh)
        else:
            ax.set_xlim(0,None)
        ax.set_ylim(-2,2)
        ax.text(0.025,0.95, self.letter, transform=ax.transAxes, fontsize=14, va='top', ha='left')
        cbar=fig.colorbar(cntr)
        fig.tight_layout()
        fig.savefig(os.path.join(self.writePath, f"VelocityContour_{self.name}.png"))
        plt.close(fig)
    def PressureContour(self):
        """
            PressureContour() --> Plots the pressure contour for the most upstream and downstream points overlayed with position plot
        """
        fig,ax = self.PlotSpecs(xlabel=r'$t/t_{c}$',ylabel=r'$\frac{(z-z_{0})}{R_{0}}$')
        minColor=0.5*self.MaxPressure
        maxColor=self.MaxPressure
        levels=np.linspace(minColor,maxColor,256)
        cntr = ax.contourf((self.TrimInterpTimeM-self.timeShift)/(self.tMinVolume-self.timeShift),(self.TrimInterpPosM-self.y0)/self.R0,
                           self.TrimInterpPresM,levels=levels,cmap='RdYlBu_r',extend='min')
        ax.plot((self.TrimtArray-self.timeShift)/(self.tMinVolume-self.timeShift),(self.TrimyInterfaceDownstream-self.y0)/self.R0,'--',color='black')
        ax.plot((self.TrimtArray-self.timeShift)/(self.tMinVolume-self.timeShift),(self.TrimyInterfaceUpstream-self.y0)/self.R0,'-',color='black')
        #ax.text(0.025,0.95, self.letter, transform=ax.transAxes, fontsize=12, va='top', ha='left')
        cbar=fig.colorbar(cntr)
        ax.set_xlim(0,None)
        ax.set_ylim(-1.5,1.5)

        fig.tight_layout()
        fig.savefig(os.path.join(self.writePath, f"PressureContour_{self.name}.png"))
        plt.close(fig)

    def CurvaturePlot(self):
        """
            Curvature Plot() --> Plots of curvature over time, with the point at which the re-entrant jet forms, if at all
        """
        fig,ax = PlotSpecs(xlabel=r'$t$',ylabel=r'$\kappa$')
        ax.plot(self.tArray,self.Curvature,'-',color='black')
        fig.tight_layout()
        fig.savefig(f"paper-images/Curvature_{self.name}.png")
        plt.close(fig)

    def Run(self):
        """
            Run() --> Runs class methods
        """
        ### Post - Processning - Shouldn't be switched off ###
        self.maxPressure()
        self.Interface()
        #self.WaterInterface()
        self.Interpolate()
        self.TrimArrays()
       
        #### Parameters and Curvature ###
        print('Getting Minimum Volume time')
        self.JettingParameters()

        print('Plotting curvature at upstream pole')
        self.BubbleCurvature()
        self.CurvaturePlot()
        
        ####### Plotting ##############
        print('Plotting maximum centerline pressure')
        self.MaxPressurePlot()
        
        #print('Plotting Volume')
        #self.VolumePlot()
       
        print('Plotting Pressure Difference')
        self.PressurePlot()

        print('Plotting Interface Velocity')
        self.VelocityPlot()

        print('Plotting Interface Acceleration')
        self.AccelerationPlot()
        
        print('Plotting Velocity Contour')
        self.VelocityContour()
            
        print('Plotting Pressure Contour')
        self.PressureContour()

## Simulation Objects
numSims=15
Simulations=[]*numSims #Empty simulations array
MainSimulations=[]*4

writePath='/media/vyntra/Disk_1/nathan/jetting-sims/paper-images'
#Running objects and appending to the arrays


EOS_Colonius=np.array([6.68,4050e5,998])
EOS_Paper=np.array([2.955,7.22e8,998])


Shock0P5=simulationRun('/media/vyntra/Disk_1/nathan/jetting-sims/Shock0P5',writePath,5e-4,8e-3,6e-3,'0.5e6',4e-6,4.5e-5,0,0,0.2741,EOS_Paper,'(a)','Shock0P5')
Shock0P5.Run() # Doesn't Jet
Simulations.append(Shock0P5)
MainSimulations.append(Shock0P5)

Shock0P7=simulationRun('/media/vyntra/Disk_1/nathan/jetting-sims/Shock0P7',writePath,5e-4,8e-3,6e-3,'0.7e6',4e-6,4.5e-5,4.30069e-05,0,0.3426,EOS_Paper,'(a)','Shock0P7')
Shock0P7.Run()
Simulations.append(Shock0P7)

Shock0P9=simulationRun('/media/vyntra/Disk_1/nathan/jetting-sims/Shock0P9',writePath,5e-4,8e-3,6e-3,'0.9e6',4e-6,3.5e-5,2.63582e-05,0,0.5480,EOS_Paper,'(a)','Shock0P9')
Shock0P9.Run()
Simulations.append(Shock0P9)

Shock1P0=simulationRun('/media/vyntra/Disk_1/nathan/jetting-sims/Shock1P0',writePath,5e-4,8e-3,6e-3,'1e6',4e-6,3e-5,2.39552e-05,0,0.6165,EOS_Paper,'(b)','Shock1P0')
Shock1P0.Run()
Simulations.append(Shock1P0)
MainSimulations.append(Shock1P0)

Shock2P0=simulationRun('/media/vyntra/Disk_1/nathan/jetting-sims/Shock2P0',writePath,5e-4,8e-3,6e-3,'2e6',4e-6,3e-5,1.64505e-05,0,1.301,EOS_Paper,'(b)','Shock2P0')
Shock2P0.Run()
Simulations.append(Shock2P0)

Shock4P0=simulationRun('/media/vyntra/Disk_1/nathan/jetting-sims/Shock4P0',writePath,5e-4,8e-3,6e-3,'4e6',4e-6,2.5e-5,1.29027e-05,0,2.6678,EOS_Paper,'(b)','Shock4P0')
Shock4P0.Run()
Simulations.append(Shock4P0)

Shock8P0=simulationRun('/media/vyntra/Disk_1/nathan/jetting-sims/Shock8P0',writePath,5e-4,8e-3,6e-3,'8e6',4e-6,2e-5,1.05658e-05,0,5.3939,EOS_Paper,'(b)','Shock8P0')
Shock8P0.Run()
Simulations.append(Shock8P0)

Shock10P0=simulationRun('/media/vyntra/Disk_1/nathan/jetting-sims/Shock10P0',writePath,5e-4,8e-3,6e-3,'10e6',4e-6,1.5e-5,1.00013e-05,0,6.7532,EOS_Paper,'(c)','Shock10P0')
Shock10P0.Run()
Simulations.append(Shock10P0)
MainSimulations.append(Shock10P0)

Shock100P0=simulationRun('/media/vyntra/Disk_1/nathan/jetting-sims/Shock100P0',writePath,5e-4,8e-3,6e-3,'100e6',4e-6,7.5e-6,6.53978e-06,0,65.4931,EOS_Paper,'(d)','Shock100P0')
Shock100P0.Run()
Simulations.append(Shock100P0)
MainSimulations.append(Shock100P0)

'''
ColoniusM=simulationRun('/media/vyntra/Disk_1/nathan/jetting-sims/ColoniusM',writePath,5e-4,8e-3,6e-3,'35.3e6',0,10e-6,7.06953e-06,0,20.9037,EOS_Colonius,'(d)','ColoniusM')
ColoniusM.Run()
Simulations.append(ColoniusM)




Shock1P0ST50R=simulationRun('/media/vyntra/Disk_1/nathan/jetting-sims/Shock1P0ST50R',writePath,25e-6,8e-4,6e-4,'1e6',2.7e-7,1.6e-6,1.5365e-06,0.073,0.6165,EOS_Paper,'(b)','Shock1P0ST50R')
Shock1P0ST50R.Run()
'''


JettingImpactPlot(Simulations)
ShockVolumePlot(Simulations)
ProlongationPlot(Simulations)
VolumePlot(MainSimulations)

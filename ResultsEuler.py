import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter


def AppendArrays(arr1,val1,arr2,val2,arr3,val3,arr4,val4):
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
    return

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
    fig,ax = PlotSpecs(xlabel=r'$\log(p_{s}/p_{0})$',ylabel=r'$\Delta T_{j}/ \Delta T_{c}$')
    tShockImpact=dX/Cl
    dParamJet=(tJetting-tShockImpact)/(minVolumes-tShockImpact)
    ax.plot(np.log(pressures/1e5),dParamJet,color='black',linestyle=None,marker='o',markersize=3)
    fig.savefig(f"paper-images/JettingImpactPlot.png")
    plt.close(fig)

def ShockVolumePlot(minVolumes,pressures,R0,Cl):
    """ ShockVolumePlot(minVolumes,pressures,R0,Cl) --> Creates a plot of shock pressure ratios vs the time taken for the shock to diffract / minimum volume (collapse)
    Parameters:
        minVolumes: Array of times when bubble collapsed (reached minimum volume)
        pressures: Array of shock overpressures
        R0: Initial bubble radius
        Cl: Speed of sound in water
    """
    fig,ax = PlotSpecs(xlabel=r'$\log(p_{s}/p_{0})$',ylabel=r'$T_{d}/T_{c}$')
    tShockDiffract=2*R0/Cl
    dParamJet=tShockDiffract/minVolumes
    rho=998
    tRayleigh=0.915*R0*np.sqrt(998/(pressures-1e5))
    tRayleigh=tShockDiffract/tRayleigh
    ax.plot(np.log(pressures/1e5),dParamJet,color='black',linestyle=None,marker='o',markersize=3)
    ax.plot(np.log(pressures/1e5),tRayleigh,color='black',linestyle='dashed',marker='o',markersize=3)
    ax.axhline(y=0.05,linestyle='--',color='Black')
    fig.tight_layout()
    fig.savefig(f"paper-images/ShockVolumePlot.png")
    plt.close(fig)

def VolumePlot(Simulations,Colours,Labels,timeShift):
        """
            VolumePlot() --> Plots of volume over time, with the point at which the re-entrant jet forms, if at all
        """
        fig,ax = PlotSpecs(xlabel=r'$t/\tau$',ylabel=r'$V/V_{0}$')
        for idx,sim in enumerate(Simulations):
            vB=sim.VolumeArr[:,1]/sim.VolumeArr[0,1]
            ax.plot((sim.VolumeArr[:,0]-timeShift)/sim.tRayleigh,vB,'-',color=Colours[idx], label=Labels[idx])
            if(sim.tJetting != 0):
                ax.axvline(x=(sim.tJetting-timeShift)/sim.tRayleigh,linestyle='--',color=Colours[idx])
        fig.tight_layout()
        ax.set_xlim(0,1.5)
        ax.legend(loc='upper right')
        fig.savefig(f"paper-images/Volume.png")
        plt.close(fig)


class simulationRun:
    """
        class simulationRun --> A post-processing simulation object for each shock induced bubble collapse simulation
    Important Functions
        interface() --> Does all the manual labour regarding interface extraction
        Run() --> Calls class methods which are actually useful
    """
    def __init__(self,filepath,pressure,tfinal,particle_speed,letter,name):
        #Initial Stuff
        self.filepath=filepath
        self.pressure=pressure
        self.tfinal=tfinal
        self.name=name
        self.letter=letter

        # Paramaters
        self.R0=0.000338 # Initial Bubble Radius (m)
        self.Cl=1462 # Speed of sound in water (m/s)
        self.Rho_l=998 # Density of water (kg/m^3)
        self.P0=1e5 # Pressure of bubble at t=0
        self.up=particle_speed # Speed of particles behind shockwave front
        
        # For Area and Sphericity plots
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

        #For Diffuse Volume Fraction Matrix
        '''
        self.InterfaceTemplate = os.path.join(filepath, "interface", "alpha{:04d}.txt")
        self.NumRuns = len(self.AreaArr)
        self.InterfaceList = [np.loadtxt(self.InterfaceTemplate.format(i)) for i in range(self.NumRuns)]
        self.InterfaceMatrix=np.stack(self.AlphaList,axis=1)
        '''

        #For Position Matrix
        self.PositionMatrix=np.loadtxt(os.path.join(filepath,"pressure/position0000.txt"))+0.022 # Starts measuring from 0.022

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

        self.tRayleigh=0.915*self.R0*np.sqrt(self.Rho_l/(float(self.pressure)-self.P0)) # Time taken for bubble collape according to Rayleigh model
        
        # Values
        self.y0=0.025 #Point which I'm uncertain of
        self.tJetting=0 #Time at which jetting occurs (initialising)

        self.writePath='/home/exy214/Documents/cavitation/data/jetting_ws_2025/paper-images/'
        self.timeShift=2.6e-6
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
    
    # Curvature Calculation
    def curvature(self,x,y):

        dx=np.gradient(x)
        dy=np.gradient(y)

        ddx=np.gradient(dx)
        ddy=np.gradient(dy)

        c1=np.abs(dx*ddy-dy*ddx)
        c2=(dx*dx+dy*dy)**1.5
        cur=c1/c2
        return cur

    # Flip in y axis
    def FlipYAxAndJoin(self,x,y):
        tempX=-1*x[1:][::-1]
        tempY=y[1::][::-1]
        x=np.concatenate((tempX,x))
        y=np.concatenate((tempY,y))
        return x,y
    '''
    def BubbleCurvature(self):
        for i in range(self.NumRuns):
            Interface=self.InterfaceMatrix[:,i]
            idxMUP=np.argmin(Interface[:,0])  # Left most point
            Interface=Interface[idxMUP:idxMUP+int(0.05*len(Interface))] #Find next 10 points and add to array (nearest most upstream point)
            x=Interface[:,0]
            y=Interface[:,1]
            x,y=self.FlipYAxAndJoin(x,y)

            x=savgol_filter(x,window_length=int(0.15*len(x)),polyorder=2)
            y=savgol_filter(y,window_length=int(0.15*len(y)),polyorder=2)

            # Bubble curvature
            bubbleCur=self.curvature(x,y)
            self.Curvature[i]=bubbleCur[0]
    '''

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
        fig,ax = self.PlotSpecs(xlabel=r'$t/t_{\tau}$',ylabel=r'$u$')
        ax.plot((self.TrimtArray-self.timeShift)/self.tRayleigh, self.TrimuInterfaceDownstream,'--',color='black',label='Downstream')
        ax.plot((self.TrimtArray-self.timeShift)/self.tRayleigh,self.TrimuInterfaceUpstream,'-',color='black',label='Upstream')
        #ax.plot((self.TrimtArray-tCollapse)/self.tRayleigh,self.Curvature*np.max(self.TrimuInterfaceUpstream),'*',color='black',label='Upstream')
        ax.set_xlim(0,1.05)
        ax.text(0.025,0.95, self.letter, transform=ax.transAxes, fontsize=12, va='top', ha='left')
        fig.tight_layout()
        fig.savefig(os.path.join(self.writePath, f"VelocityPlot_{self.name}.png"))
        plt.close(fig)
    def PositionPlot(self):
        """
            PositionPlot() --> Plots the interfacial position for most upstream and downstream points
        """
        fig,ax = self.PlotSpecs(xlabel=r'$t/\tau_{c}$',ylabel=r'$\frac{(y-y_{0})}{R_{0}}$')
        ax.plot((self.TrimtArray-self.timeShift)/self.tRayleigh,self.TrimyInterfaceDownstream,'--',color='black',label='Downstream')
        ax.plot((self.TrimtArray-self.timeShift)/self.tRayleigh,self.TrimyInterfaceUpstream,'-',color='black',label='Upstream')
        if(self.tJetting != 0):
            ax.axvline(x=(self.tJetting-self.timeShift)/self.tRayleigh,linestyle='--',color='black')
        ax.set_xlim(0,None)
        ax.text(0.025,0.95, self.letter, transform=ax.transAxes, fontsize=12, va='top', ha='left')
        ax.legend()
        fig.tight_layout()
        fig.savefig(os.path.join(self.writePath, f"PositionPlot_{self.name}.png"))
        plt.close(fig)
    def VelocityContour(self):
        """
            VelocityContour() --> Plots the velocity contour for most upstream and downstream points overlayed with position plot
        """
        fig,ax = self.PlotSpecs(xlabel=r'$t/t_{\tau}$',ylabel=r'$\frac{(y-y_{0})}{R_{0}}$')
        ax.contourf((self.TrimtArray-self.timeShift)/self.tRayleigh,(self.PositionMatrix.squeeze()-self.y0)/self.R0,self.TrimVelocityMatrix,cmap='jet')
        ax.plot((self.TrimtArray-self.timeShift)/self.tRayleigh,(self.TrimyInterfaceDownstream-self.y0)/self.R0,'--',color='black')
        ax.plot((self.TrimtArray-self.timeShift)/self.tRayleigh,(self.TrimyInterfaceUpstream-self.y0)/self.R0,'-',color='black')
        ax.set_xlim(0,None)
        ax.set_ylim(-2,2)
        ax.text(0.025,0.95, self.letter, transform=ax.transAxes, fontsize=14, va='top', ha='left')
        fig.tight_layout()
        fig.savefig(os.path.join(self.writePath, f"VelocityContour_{self.name}.png"))
        plt.close(fig)
    def PressureContour(self):
        """
            PressureContour() --> Plots the pressure contour for the most upstream and downstream points overlayed with position plot
        """
        fig,ax = self.PlotSpecs(xlabel=r'$\frac{(y-y_{0})}{R_{0}}$',ylabel=r'$t/\tau_{c}$')
        ax.contourf(self.TrimtArray/(self.R0/self.Cl),(self.PositionMatrix.squeeze()-self.y0)/self.R0,self.TrimPressureMatrix,cmap='Greys')
        ax.plot(self.TrimtArray/(self.R0/self.Cl),(self.TrimyInterfaceDownstream-self.y0)/self.R0,'--',color='black')
        ax.plot(self.TrimtArray/(self.R0/self.Cl),(self.TrimyInterfaceUpstream-self.y0)/self.R0,'-',color='black')
        ax.text(0.025,0.95, self.letter, transform=ax.transAxes, fontsize=12, va='top', ha='left')
        fig.savefig(os.path.join(self.writePath, f"PressureContour_{self.name}.png"))
        plt.close(fig)
    def JettingParameters(self):
        """
            JettingParameters() --> Returns paramters of jetting, minimum volume, jetting maximum jet velocity and response time
        Return:
            tJetting: Time at which re-entrant jet forms
            MinVolumeTime: Time at which the minimum bubble volume is reached
            JettingVolume: Maximum jet velocity, name needs changing...
            tMove: Response time, time at which the bubble begins to collapse. Useless. 
        """
        if(self.tJetting==0):
            return None,None,None
        MinVolumeTime=self.MinVolumeTime()
        JettingVolume=np.max(self.uInterfaceUpstream)
        return self.tJetting, MinVolumeTime,JettingVolume



    def Run(self):
        """
            Run() --> Runs class methods
        Return:
            tJ: Time at which re-entrant jet forms
            mV: Time at which the minimum bubble volume is reached
            jV: Maximum jet velocity
            tMove: Time at which the bubble collapse begins
        """
        self.Interface()
        self.TrimArrays()
        print('Plotting Interface Velocity')
        self.VelocityPlot()

        #print('Plotting Interface Position')
        #self.PositionPlot()

        #print('Plotting Velocity Contour')
        self.VelocityContour()

        #print('Plotting Pressure Contour')
        #self.PressureContour()

        tJ,mV,jV=self.JettingParameters()
        return tJ,mV,jV

## Plotting Colours
colorBlue = "#0072B2"
colorOrange = "#E69F00"
colorBlack = "#000000"
colorGrey = "#999999"

Colours = [colorBlue, colorOrange, colorBlack, colorGrey]

Labels = [r'$p_{s}/p_{0} = 1000$',r'$p_{s}/p_{0} = 100$', r'$p_{s}/p_{0} = 10$', r'$p_{s}/p_{0} = 5$']

## Empty Arrays
jettingTimes=[]*4
minVolumes=[]*4
maxJetting=[]*4
tMove=[]*4
pressures=[]*4

## Simulation Objects

Simulations=[]*4

#Running objects and appending to the arrays
WS100e6=simulationRun('/home/exy214/Documents/cavitation/data/jetting_ws_2025/WS_100e6_Euler','100e6',5e-6,65.4931,'(a)','WS100e6')
tJ,mV,jV=WS100e6.Run()
AppendArrays(jettingTimes,tJ,minVolumes,mV,maxJetting,jV,pressures,100e6)
Simulations.append(WS100e6)

WS10e6=simulationRun('/home/exy214/Documents/cavitation/data/jetting_ws_2025/WS_10e6_Euler','10e6',8e-6,6.7532,'(b)','WS10e6')
tJ,mV,jV=WS10e6.Run()
AppendArrays(jettingTimes,tJ,minVolumes,mV,maxJetting,jV,pressures,10e6)
Simulations.append(WS10e6)

WS1e6=simulationRun('/home/exy214/Documents/cavitation/data/jetting_ws_2025/WS_1e6_Euler','1e6',2e-5,0.616480,'(c)','WS1e6')
tJ,mV,jV=WS1e6.Run()
AppendArrays(jettingTimes,tJ,minVolumes,mV,maxJetting,jV,pressures,1e6)
Simulations.append(WS1e6)

WS500e5=simulationRun('/home/exy214/Documents/cavitation/data/jetting_ws_2025/WS_500e5_Euler','0.5e6',5e-5,0.274504,'(d)','WS500e5')
tJ,mv,jV=WS500e5.Run()
AppendArrays(jettingTimes,tJ,minVolumes,mV,maxJetting,jV,pressures,500e5)
Simulations.append(WS500e5)



jettingTimes=np.array(jettingTimes)
minVolumes=np.array(minVolumes)
pressures=np.array(pressures)
maxJetting=np.array(maxJetting)
tMove=np.array(tMove)

# Plotting
JettingImpactPlot(jettingTimes,minVolumes,pressures,0.004662,1462)
ShockVolumePlot(minVolumes,pressures,0.000338,1462)
VolumePlot(Simulations,Colours,Labels,3.15e-6)

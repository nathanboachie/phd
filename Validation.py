import numpy as np 
import os
import matplotlib.pyplot as plt

'''
    @brief Set plotting parameters
'''
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
    'savefig.dpi': 600
})

blue = "#0072B2"   # Okabe-Ito blue
orange = "#E69F00" # Okabe-Ito orange
black = "#000000"

## Global Parameters ##
filepath='/home/exy214/Documents/cavitation/data/jetting_ws_2025/'

'''
    @brief Define figure and axes objects
    @param xlabel String of xlabel
    @param ylabel String of y label
    @return Returns a figure and axes object
'''
def PlotSpecs(xlabel,ylabel) :
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

        return fig, ax

'''
    @brief Validation plot for each bubble resolution
    @param dict1 Dictinonary containing arrays for a bubble resolution
    @param dict2 See dict1
    @param X0 Initial location of bubble centroid
    @param R0 Initial bubble radius
    @param Cl Speed of sound within liquid
    @param tShift Time at which shock reaches bubble proximal side
'''
def ValidationPlot(dict1,dict2,X0,R0,Cl,tShift,Comparison):
    fig,ax=PlotSpecs(xlabel=r'$t/(R_{0}/c_{L})$',ylabel=None)
    Cl=1687
    
    ax.plot((dict1["Volume"][:,0]-tShift)/(R0/Cl),dict1["Volume"][:,1]/dict1["Volume"][0,1],'-',marker='o',markevery=4,markersize=2,color='black',label=dict1["Name"])
    ax.plot((dict2["Volume"][:,0]-tShift)/(R0/Cl),dict2["Volume"][:,1]/dict2["Volume"][0,1],'-',marker='^',markevery=4,markersize=2,color='blue',label=dict2["Name"])

    ax.plot((dict1["Centroid"][:,0]-tShift)/(R0/Cl),-1*(dict1["Centroid"][:,1]-X0)/R0,linestyle='--',marker='o',markevery=4,markersize=2,color='black')
    ax.plot((dict2["Centroid"][:,0]-tShift)/(R0/Cl),-1*(dict2["Centroid"][:,1]-X0)/R0,linestyle='--',marker='^',markevery=4,markersize=2,color='blue')


    #ax.plot((dict1["NS"][:,0]-tShift)/(R0/Cl),dict1["NS"][:,1],linestyle='dashdot',marker='o',markevery=4,markersize=2,color='black',label=dict1["Name"])
    #ax.plot((dict2["NS"][:,0]-tShift)/(R0/Cl),dict2["NS"][:,1],linestyle='dashdot',marker='^',markevery=4,markersize=2,color='black',label=dict2["Name"])
    
    ax.set_xlim(0,10)
    ax.set_ylim(-0.6,1.05)

    if Comparison=='Colonius': 

        ValidationVolume=np.loadtxt(f"{filepath}Colonius-Data/Volume.txt")
        ValidationNonSphericity=np.loadtxt(f"{filepath}Colonius-Data/NonSphericity.txt",delimiter=',')
        ValidationCentroid=np.loadtxt(f"{filepath}Colonius-Data/Centroid.txt",delimiter=',')

        ax.plot(ValidationVolume[:,0],ValidationVolume[:,1],'-',marker='*',markevery=4,markersize=2,color='orange',label="Johnsen \\& Colonius (2009)")
        ax.plot(ValidationCentroid[:,0],ValidationCentroid[:,1],linestyle='--',marker='*',markevery=4,markersize=2,color='orange')
        #ax.plot(ValidationNonSphericity[:,0],ValidationNonSphericity[:,1],linestyle='dashdot',marker='*',markevery=4,markersize=2,color='orange',label="Johnsen & Colonius (2009)")
        print('Printing Colonius Parameter Validation')
        ax.text(0.025,0.95, '(a)', transform=ax.transAxes, fontsize=14, va='top', ha='left')
        ax.legend(loc='upper right')
        fig.savefig(os.path.join("paper-images/ParameterValidationColonius.png"))

    elif Comparison =='FronTier':
        
        ValidationVolume=np.loadtxt(f"{filepath}FronTier-Data/Volume.txt")
        #ValidationNonSphericity=np.loadtxt(f"{filepath}FronTier-Data/NonSphericity.txt",delimiter=',')
        ValidationCentroid=np.loadtxt(f"{filepath}FronTier-Data/Centroid.txt")


        ax.plot((ValidationVolume[:,0]-tShift)/(R0/Cl),ValidationVolume[:,1]/ValidationVolume[0,1],'-',marker='*',markevery=4,markersize=2,color='orange',label="Bempedelis \\& Ventikos (2021)")
        ax.plot((ValidationCentroid[:,0]-tShift)/(R0/Cl),-1*(ValidationCentroid[:,1]-X0)/R0,linestyle='--',marker='*',markevery=4,markersize=2,color='orange')
        #ax.plot(ValidationNonSphericity[:,0],ValidationNonSphericity[:,1],linestyle='dashdot',marker='*',markevery=4,markersize=2,color='blue',label="Johnsen")
        print('Printing FronTier Parameter Validation')
        ax.text(0.025,0.95, '(a)', transform=ax.transAxes, fontsize=14, va='top', ha='left')
        ax.legend(loc='upper right')
        fig.savefig(os.path.join("paper-images/ParameterValidationFronTier.png"))
    else:
        print(f'There is no data available using the {Comparison} simulation')


def VelocityValidationPlot(dict1,dict2,R0,Cl,tSimShift,tJohnsen,Comparison):
    fig,ax=PlotSpecs(xlabel=r'$t/(R_{0}/c_{L})$',ylabel=r'$u/c_{L}$')
       
    ax.plot((dict1["Time"]-tSimShift)/(R0/Cl),dict1["UUP"]/Cl,linestyle='-',marker='o',markevery=4,markersize=2,color='black',label=dict1["Name"])
    ax.plot((dict1["Time"]-tSimShift)/(R0/Cl),dict1["UDO"]/Cl,linestyle='--',marker='o',markevery=4,markersize=2,color='black')
    
    ax.plot((dict2["Time"]-tSimShift)/(R0/Cl),dict2["UUP"]/Cl,linestyle='-',marker='^',markevery=4,markersize=2,color='blue',label=dict2["Name"])
    ax.plot((dict2["Time"]-tSimShift)/(R0/Cl),dict2["UDO"]/Cl,linestyle='--',marker='^',markevery=4,markersize=2,color='blue')

    ax.set_xlim(0,10)
    ax.set_ylim(-0.7,1.05)

    if Comparison == 'Colonius':
    
        JUUP=np.loadtxt(f"{filepath}Colonius-Data/VelocityUpstream.txt",delimiter=',')
        JUDO=np.loadtxt(f"{filepath}Colonius-Data/VelocityDownstream.txt",delimiter=',')

        ax.plot(JUUP[:,0]-tJohnsen,-1*JUUP[:,1],linestyle='-',marker='*',markevery=4,markersize=2,color='orange',label="Johnsen \\& Colonius (2009)")
        ax.plot(JUDO[:,0]-tJohnsen,-1*JUDO[:,1],linestyle='--',marker='*',markevery=4,markersize=2,color='orange')

        print('Printing Colonius Velocity Validation')
        ax.text(0.025,0.95, '(b)', transform=ax.transAxes, fontsize=14, va='top', ha='left')
        ax.legend(loc='upper right')
        fig.savefig(os.path.join("paper-images/VelocityValidationColonius.png"))
    
    if Comparison == 'FronTier':
        
        FP=np.loadtxt(f"{filepath}FronTier-Data/Velocity.txt")
        FT=FP[:,0]
        UP=FP[:,1]
        DP=FP[:,2]

        FUUP=np.gradient(UP,FT)
        FUDO=np.gradient(DP,FT)

        ax.plot((FT-tSimShift)/(R0/Cl),FUUP/Cl,linestyle='-',marker='o',markevery=4,markersize=2,color='orange',label='Bempedelis \\& Ventikos (2021)')
        ax.plot((FT-tSimShift)/(R0/Cl),FUDO/Cl,linestyle='-',marker='o',markevery=4,markersize=2,color='orange')
       
        print('Printing FronTier Velocity Validation')
        ax.text(0.025,0.95, '(b)', transform=ax.transAxes, fontsize=14, va='top', ha='left')
        ax.legend(loc='upper right')
        fig.savefig(os.path.join("paper-images/VelocityValidationFronTier.png"))


class simulationRun:
    def __init__(self,filepath,pressure,tfinal,name):
        #Initial Stuff
        self.filepath=filepath
        self.pressure=pressure
        self.tfinal=tfinal
        self.name=name
        self.R0=0.0005
        self.Rho_l=998
        # For Area and Sphericity splots
        self.AreaArr=np.loadtxt(os.path.join(filepath, "Area.txt"))
        self.CentroidArr=np.loadtxt(os.path.join(filepath,"Centroid.txt"))
        self.PerimeterArr=np.loadtxt(os.path.join(filepath,"Perimeter.txt"))
        self.VolumeArr=np.loadtxt(os.path.join(filepath, "Volume.txt"))

        # Derived Values

        #Average Radius
        self.RadiusArr=np.empty_like(self.AreaArr)
        self.RadiusArr[:,0]=self.AreaArr[:,0]
        self.RadiusArr[:,1]=self.R0*np.sqrt(self.AreaArr[:,1]/self.AreaArr[0,1])

        # Non Sphericity
        self.NSArr=np.empty_like(self.AreaArr)
        self.NSArr[:,0]=self.AreaArr[:,0]
        self.NSArr[:,1]=self.AreaArr[:,1]/(self.RadiusArr[:,1]*self.PerimeterArr[:,1])
       
        #Validation Plot Info
        self.X0=0.5*16e-3
        self.R0=5e-4
        self.Cl=1646.66
        self.tShiftBottom=5e-6
        self.V0 = (4/3)*np.pi*self.R0**3
        
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
        self.PositionMatrix=np.loadtxt(os.path.join(filepath,"velocity/position0000.txt"))+4e-3
    
        # Time Matrix
        self.tArray=np.linspace(0,tfinal,self.NumRuns)

        # Interface Values
        #self.yInterfaceUpstream = np.zeros(((self.NumRuns,1)))
        self.uInterfaceUpstream = np.zeros(((self.NumRuns,1)))
        #self.pInterfaceUpstream = np.zeros(((self.NumRuns,1)))
        #self.yInterfaceDownstream = np.zeros(((self.NumRuns,1)))
        self.uInterfaceDownstream = np.zeros(((self.NumRuns,1)))
        #self.pInterfaceDownstream = np.zeros(((self.NumRuns,1)))
        
        '''
        # Trimmed Interface Values
        
        Why? This is a trim matrix that ends at the point of jetting
        The data looks weird including information past jetting
        '''

        #self.TrimyInterfaceUpstream = np.zeros(((self.NumRuns,1)))
        self.TrimuInterfaceUpstream = np.zeros(((self.NumRuns,1)))
        #self.TrimpInterfaceUpstream = np.zeros(((self.NumRuns,1)))
        self.TrimaInterfaceUpstream = np.zeros(((self.NumRuns,1)))
        #self.TrimyInterfaceDownstream = np.zeros(((self.NumRuns,1)))
        self.TrimuInterfaceDownstream = np.zeros(((self.NumRuns,1)))
        #self.TrimpInterfaceDownstream = np.zeros(((self.NumRuns,1)))
        self.TrimaInterfaceDownstream = np.zeros(((self.NumRuns,1)))
        self.TrimtArray = np.zeros(((self.NumRuns,1)))
        #self.TrimPressureMatrix=np.zeros((self.NumRuns,len(self.PositionMatrix)))
        self.TrimVelocityMatrix=np.zeros((self.NumRuns,len(self.PositionMatrix)))
        self.TrimAlphaMatrix=np.zeros((self.NumRuns,len(self.PositionMatrix)))

    '''
    Interface, look at interface, if there are no points w/volume fraction above 0.5 it means there is no interface along centerline
    The data at this point is all from the interface
    Find the most upstream point and downstream point and extract the values at these locations, which have a volume fraction above a point
    The interface value is then one-dimensional
    Repeats for downstream point
    '''

    def Interface(self):
        for i in range(self.NumRuns):
            Iarr=np.where(self.AlphaMatrix[:,i]>0.5)
            if(len(Iarr[0]) < 2):
                break
            else:
                YGasPosition=self.PositionMatrix[Iarr]
                YMax=np.max(YGasPosition)
                I=int(np.where(self.PositionMatrix==YMax)[0][0])
                self.uInterfaceDownstream[i]=self.VelocityMatrix[:,i][I]

                YMin=np.min(YGasPosition)
                I=int(np.where(self.PositionMatrix==YMin)[0][0])
                self.uInterfaceUpstream[i]=self.VelocityMatrix[:,i][I]

    def JetTime(self):
        for i in range(self.NumRuns):
            Iarr=np.where(self.AlphaMatrix[:,i]>0.5)
            if(len(Iarr[0]) < 1):
                self.tJetting=self.tArray[i] 
                print(f'Jetting Achieved @ {self.tArray[i]}')
                break
    '''
    Trimming arryas based upon length of velocity interface where nothing is being output
    '''

    def TrimArrays(self):
        self.TrimuInterfaceUpstream=np.trim_zeros(self.uInterfaceUpstream,'b')
        self.TrimuInterfaceDownstream=np.trim_zeros(self.uInterfaceDownstream,'b')
        self.TrimtArray=self.tArray[0:len(self.TrimuInterfaceDownstream)]
        self.TrimVelocityMatrix=self.VelocityMatrix[:,:len(self.TrimuInterfaceDownstream)]
        self.TrimAlphaMatrix=self.AlphaMatrix[:,:len(self.TrimuInterfaceDownstream)]
   
    def Run(self):
        self.Interface()
        self.TrimArrays()
        self.JetTime()
        return self.tJetting

ColoniusL=simulationRun('/home/exy214/Documents/cavitation/data/jetting_ws_2025/ColoniusL','35.3e6',8e-6,'ColoniusL')
ColoniusL.Run()
ColoniusM=simulationRun('/home/exy214/Documents/cavitation/data/jetting_ws_2025/ColoniusM','35.3e6',8e-6,'ColoniusM')
tJet=ColoniusM.Run()


## Plotting betweeen simulations
x0=0.0075
X0=0.5*16e-3
R0 = 5e-4
Cl = 1687


## Dictionaries 
ColoniusLDict={
        "Volume": ColoniusL.VolumeArr,
        "Centroid": ColoniusL.CentroidArr,
        "NS": ColoniusL.NSArr,
        "UUP":ColoniusL.TrimuInterfaceUpstream,
        "UDO":ColoniusL.TrimuInterfaceDownstream,
        "Time": ColoniusL.TrimtArray,
        "Name": "128 ppbr"
        }

ColoniusMDict={
        "Volume": ColoniusM.VolumeArr,
        "Centroid": ColoniusM.CentroidArr,
        "NS": ColoniusM.NSArr,
        "UUP":ColoniusM.TrimuInterfaceUpstream,
        "UDO":ColoniusM.TrimuInterfaceDownstream,
        "Time":ColoniusM.TrimtArray,
        "Name": "256 ppbr"
        }

tShock=x0/Cl
print(tShock)
tCollapseJohnsen=1.4

ValidationPlot(ColoniusLDict,ColoniusMDict,X0,R0,Cl,tShock,'FronTier')
VelocityValidationPlot(ColoniusLDict,ColoniusMDict,R0,Cl,tShock,tCollapseJohnsen,'FronTier')

import numpy as np
import os
import matplotlib.pyplot as plt

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

        return fig, ax

def ValidationPlot(dict1,dict2,X0,R0,Cl,tCollapse=4e-6):
    fig,ax=PlotSpecs(xlabel=r'$t/(R_{0}/C_{l})$',ylabel=None)
    V0=(4.0/3.0)*np.pi*R0**3
    Cl=1642
    
    ax.plot((dict1["Volume"][:,0]-tCollapse)/(R0/Cl),dict1["Volume"][:,1]/V0,'-',marker='o',markevery=8,markersize=2,color='black',label=dict1["Name"])
    ax.plot((dict2["Volume"][:,0]-tCollapse)/(R0/Cl),dict2["Volume"][:,1]/V0,'-',marker='^',markevery=8,markersize=2,color='black',label=dict2["Name"])

    ax.plot((dict1["Centroid"][:,0]-tCollapse)/(R0/Cl),-1*(dict1["Centroid"][:,1]-X0)/R0,linestyle='--',marker='o',markevery=8,markersize=2,color='black',label=dict1["Name"])
    ax.plot((dict2["Centroid"][:,0]-tCollapse)/(R0/Cl),-1*(dict2["Centroid"][:,1]-X0)/R0,linestyle='--',marker='^',markevery=8,markersize=2,color='black',label=dict2["Name"])


    ax.plot((dict1["NS"][:,0]-tCollapse)/(R0/Cl),dict1["NS"][:,1],linestyle='dashdot',marker='o',markevery=8,markersize=2,color='black',label=dict1["Name"])
    ax.plot((dict2["NS"][:,0]-tCollapse)/(R0/Cl),dict2["NS"][:,1],linestyle='dashdot',marker='^',markevery=8,markersize=2,color='black',label=dict2["Name"])


    ax.set_xlim(0,None)
    fig.savefig(os.path.join("res_param_validation.png"),bbox_inches='tight')


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
        self.tShiftBottom=3e-6
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
        self.PositionMatrix=np.loadtxt(os.path.join(filepath,"velocity/position0000.txt"))
    
        self.ValidationVolume=np.loadtxt("/home/exy214/Documents/cavitation/data/jetting_ws_2025/Colonius-Data/Volume.txt",delimiter=',')
        self.ValidationNonSphericity=np.loadtxt("/home/exy214/Documents/cavitation/data/jetting_ws_2025/Colonius-Data/NonSphericity.txt",delimiter=',')
        self.ValidationCentroid=np.loadtxt("/home/exy214/Documents/cavitation/data/jetting_ws_2025/Colonius-Data/Centroid.txt",delimiter=',')
        self.VelocityUpstream=np.loadtxt("/home/exy214/Documents/cavitation/data/jetting_ws_2025/Colonius-Data/VelocityUpstream.txt",delimiter=',')
        self.VelocityDownstream=np.loadtxt("/home/exy214/Documents/cavitation/data/jetting_ws_2025/Colonius-Data/VelocityDownstream.txt",delimiter=',')
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
            Iarr=np.where(self.AlphaMatrix[:,i]>0.9)
            if(len(Iarr[0])==1):
                self.tJetting=self.tArray[i] 
                print(f'Jetting Achieved @ {self.tArray[i]}')
                break
            else:
                YGasPosition=self.PositionMatrix[Iarr]
                YMax=np.max(YGasPosition)
                I=int(np.where(self.PositionMatrix==YMax)[0][0])
                self.uInterfaceDownstream[i]=self.VelocityMatrix[:,i][I]

                YMin=np.min(YGasPosition)
                I=int(np.where(self.PositionMatrix==YMin)[0][0])
                self.uInterfaceUpstream[i]=self.VelocityMatrix[:,i][I]

    '''
    Returns collapse index based on the point where the differnce in upstream and downstream velocities begins to change to a certain degree
    '''

    def CollapseTime(self):
        diffUpstream=np.diff(self.TrimuInterfaceUpstream.flatten())
        UpstreamCollapse=np.min(np.where(diffUpstream > 0.1))
        return UpstreamCollapse

    '''
    Trimming arryas based upon length of velocity interface where nothing is being output
    '''

    def TrimArrays(self):
        self.TrimuInterfaceUpstream=np.trim_zeros(self.uInterfaceUpstream,'b')
        self.TrimuInterfaceDownstream=np.trim_zeros(self.uInterfaceDownstream,'b')
        self.TrimtArray=self.tArray[0:len(self.TrimuInterfaceDownstream)]
        self.TrimVelocityMatrix=self.VelocityMatrix[:,:len(self.TrimuInterfaceDownstream)]
        self.TrimAlphaMatrix=self.AlphaMatrix[:,:len(self.TrimuInterfaceDownstream)]

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

    '''
    def VelocityPlot(self):
        fig,ax = self.PlotSpecs(xlabel=r'$t/(R_{0}/C_{l})$',ylabel=r'$u$')
        #tCollapse=self.TrimtArray[self.CollapseTime()]
        #ax.plot((self.TrimtArray-tCollapse)/(self.R0/self.Cl),self.TrimuInterfaceDownstream,'--',color='black',label='Downstream')
        #ax.plot((self.TrimtArray-tCollapse)/(self.R0/self.Cl),self.TrimuInterfaceUpstream,'-',color='black',label='Upstream')
        #ax.plot(self.KM_T/(self.R0/self.Cl),-1*self.KM_U,linestyle='dashdot',color='black',label='Keller Miksis')
        #ax.plot(self.KM_T/(self.R0/self.Cl),self.KM_U,linestyle='dashdot',color='black',label='Keller Miksis')
        #ax.legend()
        if(self.tJetting!=0):
            #ax.set_xlim(0,(self.tJetting-tCollapse)/(self.R0/self.Cl))
            print('.')
        else:
            #ax.set_xlim(0,None) 
            print('.')
        plt.show()
        #fig.savefig(os.path.join(self.filepath, f"VelocityPlot_{self.name}.png"))
        plt.close(fig)
    '''

    def ValidationPlot(self):
        fig,ax=self.PlotSpecs(xlabel=r'$t/(R_{0}/C_{l})$',ylabel=None)
        ax.plot(self.ValidationVolume[:,0],self.ValidationVolume[:,1],'-',color='black')
        ax.plot((self.VolumeArr[:,0]-self.tShiftBottom)/(self.R0/self.Cl),self.VolumeArr[:,1]/self.V0,marker='o',markevery=8,markersize=3,color='black')

        ax.plot(self.ValidationCentroid[:,0],self.ValidationCentroid[:,1],'--',color='black')
        ax.plot((self.CentroidArr[:,0]-self.tShiftBottom)/(self.R0/self.Cl),-1*((self.CentroidArr[:,1]-self.X0)/self.R0),marker='o',markevery=8,markersize=3,color='black')
        

        ax.plot(self.ValidationNonSphericity[:,0],self.ValidationNonSphericity[:,1],linestyle='dashdot',color='black')
        #ax.plot((self.NSArr[:,0]-self.tShiftBottom)/(self.R0/self.Cl),self.NSArr[:,1],marker='o',markevery=8,markersize=3,color='black')

        ax.set_xlim(0,10)
        fig.tight_layout()
        fig.savefig(os.path.join("sim_param_validation.png"),bbox_inches='tight')


    def VelocityValidationPlot(self):
        fig,ax=self.PlotSpecs(xlabel=r'$t/(R_{0}/C_{l})$',ylabel=r'$u/c_{L}$')
        
        ax.plot((self.VelocityUpstream[:,0]),self.VelocityUpstream[:,1],linestyle='-',color='black')
        ax.plot((self.VelocityDownstream[:,0]),self.VelocityDownstream[:,1],linestyle='--',color='black')

        ax.plot((self.TrimtArray-self.tShiftBottom)/(self.R0/self.Cl),-1*self.TrimuInterfaceUpstream/self.Cl,linestyle='None',marker='o',markevery=8,markersize=3,color='black')
        ax.plot((self.TrimtArray-self.tShiftBottom)/(self.R0/self.Cl),-1*self.TrimuInterfaceDownstream/self.Cl,linestyle='None',marker='o',markevery=8,markersize=3,color='black')

        fig.savefig(os.path.join("sim_velocity_validation.png"),bbox_inches='tight')
   
    def Run(self):
        self.Interface()
        self.TrimArrays()
        #print('Plotting Interface Velocity')
        #self.VelocityPlot()
        #self.ValidationPlot()
        self.VelocityValidationPlot() 

ColoniusLL=simulationRun('/home/exy214/Documents/cavitation/data/jetting_ws_2025/ColoniusLL','35.3e6',8e-6,'ColoniusLL')
ColoniusLL.Run()
#ColoniusL=simulationRun('/home/exy214/Documents/cavitation/data/jetting_ws_2025/ColoniusL','35.3e6',8e-6,'ColoniusL')
#ColoniusL.Run()
#ColoniusM=simulationRun('/home/exy214/Documents/cavitation/data/jetting_ws_2025/ColoniusM','35.3e6',8e-6,'ColoniusM')
#ColoniusM.Run()

## Plotting betweeen simulations
X0 = 0.5*16e-3
R0 = 5e-4
Cl = 1642

ColoniusLLDict={
        "Volume": ColoniusLL.VolumeArr,
        "Centroid": ColoniusLL.CentroidArr,
        "NS": ColoniusLL.NSArr,
        "Name": "ColoniusLL"
        }
'''
## Dictionaries 
ColoniusLDict={
        "Volume": ColoniusL.VolumeArr,
        "Centroid": ColoniusL.CentroidArr,
        "NS": ColoniusL.NSArr,
        "Name": "ColoniusL"
        }

ColoniusMDict={
        "Volume": ColoniusM.VolumeArr,
        "Centroid": ColoniusM.CentroidArr,
        "NS": ColoniusM.NSArr,
        "Name": "ColoniusM"
        }


ValidationPlot(ColoniusLDict,ColoniusMDict,X0,R0,Cl)
'''

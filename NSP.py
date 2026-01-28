import os
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from ShockProperties import ShockProperties
from shockInfo import timeShift
import numpy as np

## Plotting Colours
colorBlue = "#0072B2"
colorOrange = "#E69F00"
colorBlack = "#000000"
colorGrey = "#999999"
colorPurple = "#CC79A7"
Colours = [colorBlue, colorOrange, colorBlack, colorGrey]

plt.rc('text',usetex=True)

plt.rcParams.update({
    'font.size': 14,
    'font.family': 'serif',
    'axes.labelsize': 9,
    'axes.titlesize': 9,
    'legend.fontsize': 8,
    'lines.linewidth': 1,
    'lines.markersize': 4,
    'axes.linewidth': 0.8,

    # Disable all ticks and labels
    'xtick.bottom': False,
    'xtick.top': False,
    'ytick.left': False,
    'ytick.right': False,
    'xtick.labelbottom': False,
    'ytick.labelleft': False,

    'savefig.dpi': 600,
})



class PerturbationBubble:
    """
        class simulationRun --> A post-processing simulation object for each shock induced bubble collapse simulation
    Important Functions
        interface() --> Does all the manual labour regarding interface extraction
        Run() --> Calls class methods which are actually useful
    """
    def __init__(self,filepath,writePath,R0,y0,yMeasure,pressure,startTime,tfinal,tVisJetting,sigma,EOS,letter,name):
        
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

        self.CentroidArr=np.loadtxt(os.path.join(filepath, "Centroid.txt"))
        
        #For Interface Positional Matrix
        self.InterfaceTemplate = os.path.join(filepath, "interface", "Interface{:04d}.curve")
        self.NumRuns = len(self.CentroidArr)
        self.InterfaceList = [np.loadtxt(self.InterfaceTemplate.format(i)) for i in range(self.NumRuns)]
        self.InterfaceMatrix=np.array(self.InterfaceList,dtype=object)

        # Time Matrix
        self.tArray=self.CentroidArr[:,0]
        self.tRayleigh=0.915*self.R0*np.sqrt(self.Rho_l/(float(self.pressure)-self.P0)) #Time taken for bubble collape according to Rayleigh model

        # Values
        self.y0=y0 # Assuming the start position of bubble --> Actually could be wrong unfortunaetly
        self.dx=(self.y0)-(self.R0) #Distance from bottom boundary to bubble

        # Water params
        self.gmaW=EOS[0]
        self.pinfW=EOS[1]
        self.rhoW=EOS[2]
        _,self.velW,_,_=ShockProperties(float(pressure),0)
        
        # Jetting Parameters
        self.tJetting=0 #Time at which jetting occurs (initialising)
        self.tMinVolume=0 # Time at which minimum volume is reached (collape time)
        
        self.writePath=writePath # Self explanatory
        self.Speed,self.timeShift=timeShift(float(self.pressure),self.P0,self.gmaW,self.pinfW,self.rhoW,self.velW,self.dx) # Time taken for shockwave to hit bubble

        self.Type=list(self.name)[0]
        if(self.Type =='R'):
            self.timeShift=0
        print(f'For {self.name} Shock reaches the bubble at {self.timeShift} s')

        self.tDiffract=self.timeShift+((2*R0)/(self.Speed))

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

    def PlotInterface(self,Interface):
        fig,ax=self.PlotSpecs(xlabel=r'$x$',ylabel=r'$y$')
        xCoord=Interface[:,0]
        yCoord=Interface[:,1]
        ax.plot(xCoord,yCoord,'.',linestyle='solid',markersize=2,markevery=1)
        #ax.set_xlim(-self.R0+self.y0,self.R0+self.y0)
        #ax.set_ylim(-self.R0+self.y0,self.R0+self.y0)
        fig.savefig(os.path.join(self.writePath,f'InitialDiffractPosition_{self.name}'))
          
    def InitialPositions(self):
        for i in range(len(self.tArray)):
            if self.tArray[i] > self.tDiffract:
                self.tDiffractIDX=i
                print(f'Time step at which shock has diffracted around bubble is nt = {i}')
                break
        self.InitialCoords=self.InterfaceMatrix[self.tDiffractIDX]
        ref=self.InitialCoords.copy()
        ref[:,0]*=-1
        ref=ref[::-1]
        ref=ref[1:-1]
        self.InitialCoords=np.vstack([self.InitialCoords,ref])
        self.tPerturbation=self.tArray[self.tDiffractIDX:-1]
    
    def LegPoly4(self,x):
        p0=np.ones_like(x)
        p1=x
        p2=0.5*((3*x**2)-1)
        p3=0.5*((5*x**3)-3*x)
        p4=0.125*((35*x**4)-(30*x**2)+3)
        return p0,p1,p2,p3,p4

    #### Bubble dynamics with the KM equation
    def RayleighModel(self,t,initialAmps,gma=1.4):
        pinf=float(self.pressure)
        rho=self.Rho_l
        ucd=2.0*(pinf-self.P0)/(self.Rho_l*self.Speed)
        def RayleighAmplitude(y,t):
            R, R_dot, a0, a0_dot,a1, a1_dot, a2, a2_dot = y
            pb=self.P0*(self.R0/R)**(3*gma)
            R_ddot = (pb-pinf)/(rho*R)-1.5*(R_dot**2)/R
            a0_ddot= -3*(R_dot/R)*a0_dot+(0-1.0)*(R_ddot/R)*a0
            a1_ddot=-3*(R_dot/R)*a1_dot
            a2_ddot=-3*(R_dot/R)*a2_dot+(2-1.0)*(R_ddot/R)*a2
            return [R_dot, R_ddot,a0_dot,a0_ddot,a1_dot,a1_ddot,a2_dot,a2_ddot]

        y0=[self.R0,ucd,initialAmps[0],initialAmps[1],initialAmps[2],initialAmps[3],initialAmps[4],initialAmps[5]]
        sol=odeint(RayleighAmplitude,y0,t)
        # Solve the KM equation
        R=sol[:,0]
        a0=sol[:,2]
        a1=sol[:,4]
        a2=sol[:,6]
        return R,a0,a1,a2
    
    def NonSphericalPerturbation(self):
        
        ### Generate Radial Array
        tRayleigh=self.tPerturbation
        xcoords=self.InitialCoords[:,0]
        ycoords=self.InitialCoords[:,1]
        
        centroid=self.CentroidArr[self.tDiffractIDX,1]
        ycoords=ycoords-centroid
        ### Generate Initial Radial Positions
        r0=np.sqrt(xcoords**2+ycoords**2)
        theta=np.arctan2(xcoords,ycoords)
        R0=np.mean(r0)
        dr=r0-R0

        lp0,lp1,lp2,_,_=self.LegPoly4(np.cos(theta))
        A=np.column_stack([lp0,lp1,lp2])

        ## Fitting radius(theta,0) = a0P0+a1P1+a2P2
        coeffs,*_=np.linalg.lstsq(A,dr,rcond=None)
        a0_0,a1_0,a2_0=coeffs

        print(f'Coefficients: a0(0), a1(0) and a2(0) = {a0_0,a1_0,a2_0}')
        
        V0=(2*self.velW)*lp1
        Rdot0=0.0
        dV=V0-Rdot0
        coeffs_dot,*_=np.linalg.lstsq(A,dV,rcond=None)
        a0dot_0,a1dot_0,a2dot_0=coeffs_dot
        self.coeffs=[a0_0,a0dot_0,a1_0,a1dot_0,a2_0,a2dot_0]
        
        print(f'Coefficients: a0_dot(0), a1_dot(0) and a2_dot(0) = {a0dot_0,a1dot_0,a2dot_0}')

        ## Necessary as Rayleigh Model is not valid past t > tc
        nt=len(tRayleigh)
        nx=len(xcoords)
        radMatrix=np.zeros((nt,nx))
        nspMatrix=np.zeros((nt,nx))

        rad,a0,a1,a2=self.RayleighModel(tRayleigh,self.coeffs)

        for i in range(len(r0)):
            nspMatrix[:,i]=a0*lp0[i]+a1*lp1[i]+a2*lp2[i]
            radMatrix[:,i]=rad
        RS_Matrix=nspMatrix+radMatrix
        xTransient=RS_Matrix*np.sin(theta)
        yTransient=centroid+RS_Matrix*np.cos(theta)
        return xTransient,yTransient,tRayleigh
    
    def PlotTransient(self,x,y,t):
        coordsSim=self.InterfaceMatrix[self.tDiffractIDX:-1]
        timeSteps=[0,40,80,115,130,160]
        fig,ax=self.PlotSpecs(xlabel='$x$',ylabel='$y$')
        cmap=plt.cm.Greys
        colours=cmap(np.linspace(0.75,0.25,len(timeSteps)))

        for c,i in zip(colours,timeSteps):
            ax.plot(x[i,:],y[i,:],linestyle='solid',color=c,linewidth=2)
            ax.set_xlim(-1*self.R0,self.R0) 
            ax.set_aspect('equal',adjustable='box') 
            fig.tight_layout()
        fig.tight_layout()
        plt.savefig(os.path.join(self.writePath,f'NSP.png'),bbox_inches='tight')
        plt.close(fig)
    def Run(self):
        self.InitialPositions()
        x,y,t=self.NonSphericalPerturbation()
        self.PlotTransient(x,y,t)

EOS_Paper=np.array([2.955,7.22e8,998])
writePath='/media/vyntra/Disk_1/nathan/jetting-sims/nsp'

Shock1P0L5=PerturbationBubble('/media/vyntra/Disk_1/nathan/jetting-sims/Shock1P0L5/domain/',writePath,5e-4,8e-3,6e-3,'1e6',4e-6,3e-5,2.42538e-05,0,EOS_Paper,'a','Shock1P0')
Shock1P0L5.Run()

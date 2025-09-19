import math

def timeShift(shockP,initP,gamma,pinf,rho,v1,dx):
    """ timeShift(pressure) --> Calculates time taken for shockwave to reach bubble
    Parameters:
        shockP: Overpressure of shock wave
        initP: Initial Pressure of fluid
        Gamma: Ratio of specific heats
        Pinf: Background pressure for stiffened gas EOS
        rho: Density of fluid
        v1: Initial fluid speed ahead of shockwave
        dx: Distance from bottom of domain to proximal bubble interface
    Return: 
        time: Time taken for shock to propagate from bottom boundary to proximal bubble interface
    """
    
    # Initialise Pressure
    p1=initP+pinf
    c1=math.sqrt(gamma*p1/rho) # Speed of sound in fluid ahead of shock
    p2=shockP+pinf

    # Shock calculations in F.O.R of shock wave
    pr=p2/p1 # Prssure ratio
    M1sqr = (pr-1)*(gamma+1)/(2*gamma) + 1 # Probably Hugonoiot jump conditions...not actually sure
    M1=math.sqrt(M1sqr)
    vs=v1+M1*c1
    Speed=vs
    print(f'Shock speed in this case is {Speed} m/s')
    time=dx/Speed
    print(f'Time shift is {time/(1e-6)} microseconds')
    return Speed,time







    

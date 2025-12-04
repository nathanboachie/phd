import math


# Set chosen material - 0: Air, 1: Water, 2: Solid (HMX-Type expl.)
def ShockProperties(shockP,printStates,g=2.955,pinf=7.22e8,rhoW=998):
    
    rho1=rhoW
    
    v1_glob = 0.0  # Initially quiescent flow
    p1 = 1e5
    p1 += pinf

    c1 = math.sqrt(g * p1 / rho1)

    # Give strength of shock
    p2=shockP
    p2 += pinf

    # Shock calculations in the reference frame of the shockwave
    pr_rat = p2 / p1
    M1_sqr = (pr_rat - 1) * (g + 1) / (2 * g) + 1
    M1 = math.sqrt(M1_sqr)

    vs = v1_glob + M1 * c1
    v1_loc = vs - v1_glob

    pr_rat = 1 + (2 * g) / (g + 1) * (M1_sqr - 1)

    v_rat = 1 - 2 / (g + 1) * (1 - 1 / M1_sqr)
    v2_loc = v_rat * v1_loc
    v2_glob = vs - v2_loc

    rho_rat = 1 / v_rat
    rho2 = rho_rat * rho1

    c2 = math.sqrt(g * p2 / rho2)
    M2 = v2_loc / c2
    
    if printStates==1:
        # Output
        print("\nAMBIENT STATE")
        print("Pressure:", p1 - pinf)
        print("Velocity:", v1_glob)
        print("Density:", rho1)
        print("\nSHOCKED STATE")
        print("Pressure:", p2 - pinf)
        print("Velocity:", v2_glob)
        print("Density:", rho2)
        print("\nShock-wave speed")
        print(vs)
        print("\nPre- and Post-Shock speed of sound")
        print(c1, c2)
        print("\nPre- and Post-Shock Local Mach number")
        print(M1, M2)
    shockPressure=p2-pinf
    shockVelocity=v2_glob
    shockDensity=rho2
    shockSpeed=vs
    return shockPressure,shockVelocity,shockDensity,shockSpeed

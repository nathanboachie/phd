import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.integrate import odeint
from scipy.optimize import newton


def KellerMiksis(pinf,R0):
    ## Global Parameters

    # Liquid
    rho=998
    c0=np.sqrt(2.955*(pinf+7.22e8)/rho)

    # Gas bubble
    p0=1e5
    gamma=1.4
    sigma=0 #0.073
    mu=0 #1e-3
    pb0=p0+2*sigma/R0
    rhob=1.16

    #Time
    tc_Rayleigh=0.915*R0*np.sqrt(rho/(pinf-pb0))
    nt=10000
    t=np.linspace(0,3*tc_Rayleigh,nt)

    #### Bubble dynamics with the KM equation

    def KM(y,t):
        R, R_dot = y
        P = pb0 * (R0/R)**(3*gamma) - 2*sigma/R - 4*mu*R_dot/R
        dpb_dt = (- pb0*(R0/R)**(3*gamma)*3*gamma/R*R_dot
                    + 2*sigma/R**2*R_dot
                    + 4*mu*R_dot**2/R**2
                    )

        R_ddot = ( (P - pinf)/rho*(1+R_dot/c0)
                    - 3/2*R_dot**2*(1-R_dot/(3*c0))
                    + R*dpb_dt/(rho*c0)
                    ) / ((1-R_dot/c0+4*mu/(rho*c0*R))*R)

        return [R_dot, R_ddot]

    def solve_KM():
        # Choice of initial condition. Not clearly defined for initial interface disequilibrium.
        y0 = [R0, -(pinf-pb0)/(rho*c0)]
        y0 = [R0, 0]

        solution = odeint(KM, y0, t)

        return t, solution

    # Solve the KM equation
    t_sol, sol = solve_KM()

    return t_sol,sol

'''
# Collapse time
tc_index = np.argmin(sol[:,0])
tc_KM = t_sol[tc_index]
print(tc_KM,sol[tc_index,0]/R0)
exp_mat = np.vstack([t_sol, sol[:,0], sol[:,1]]).T
np.savetxt('./KM.txt', exp_mat)
# Plot results
fig, ax = plt.subplots(1,2, figsize=(10,3))

ax[0].plot(t_sol, sol[:, 0]/sol[0,0],'k')
ax[0].set_xlabel('t[s]')
ax[0].set_ylabel('R/R_0')
ax[0].set_xlim([t_sol[0],t_sol[-1]])

ax[0].axvline(tc_Rayleigh)

ax[1].plot(t_sol, sol[:, 1],'k')
ax[1].set_xlabel('t[s]')
ax[1].set_ylabel('rdot m/s')
ax[1].set_xlim([t_sol[0],t_sol[-1]])

plt.tight_layout()
plt.show()
'''

# Export results

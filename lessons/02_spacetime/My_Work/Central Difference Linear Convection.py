#Christopher Bell
#Notebook 02-02 Space & Time
#Stability and the CFL Condition

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

def linconv(nx):
    """Solves the linear convection eqn
    
    Solves d/dt u + c d/x u = 0
    *wavespeed, c, is 1
    *domain is x in [0,2]
    *20 timesteps, w/ delta t computed using the CFL 0.5
    *intial data is hat function
    
    Produces Plot of results
    
    Parameters
    --------
    
    nx : integer
        number of internal grid points
        
    Returns
    -------
    None :  none
    """
    
def centraldiff (unm1, u, linconv, dt):
    """Returns the solution at n+1 using a central difference Eularian Method
    
    Parameters
    ----------
    unm1 : array of float
        solution at time step n-1
    u : array of float
        solution at time step n
    f : function
        function to compute RHS of SOE's
    dt : float
        time-increment
        
    Returns
    -------
    u_n_plus_1 : array of float
        solution at time step  n+1
    """
    
    return unm1 + 2*dt*linconv(u)
    
def rk2_step(u, linconv, dt):
    """Returns the solution at the next time step using 2nd order Runge-Kutta
    
    Parameter:
    ----------
    u : array of float
        solution at the previous time-step.
    f : function
        function to compute the RHS of the system of equations.
    dt : float
        time-increment.
    
    Returns:
    --------
    u_n_plus_1 : array of float
        solution at the next time step.
    """
    u_star = u + 0.5*dt*linconv(u)
    return u + dt*linconv(u_star)

nx = 50.0   
dx = 2.0/(nx-1)
nt = 20
sigma = 0.5
c = 1.0
dt = (sigma * dx)/c
    
u = np.ones(nx)
#u[.5/dx : 1/dx+1] = 2
un = np.ones(nx)
x = 0.0

f = linconv(nx)
   
for i in range(len(u)):
    x=x+dx
    if x<=0.9 or x>=1.1:
        u[i] = 0
    elif x<=1.0 and x>0.9:
        u[i]=10*(x-0.9)
    elif x<=1.1 or x>1:
        u[i]=-10*(1.1-x)

#first step
u[1] = rk2_step(u[0],f,dt)
            
for n in range(1,nt-1):
    un = u.copy()
    u[n+1] = centraldiff(u[n-1], u[n], f, dt)
    
    plt.plot(np.linspace(0,2,nx), u, color='#003366', ls='--',lw=3)
    #plt.ylim(0,2.5)
    plt.show()


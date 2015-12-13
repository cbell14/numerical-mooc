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
    dx = 2.0/(nx-1)
    nt = 20
    sigma = 0.5
    c = 1.0
    dt = (sigma * dx)/c
    
    u = np.ones(nx)
    u[.5/dx : 1/dx+1] = 2
    un = np.ones(nx)
    
    for n in range(nt):
        un = u.copy()
        u[1:] = un[1:] - ((c*(dt/dx))*(un[1:]-un[0:-1]))
    
    plt.plot(np.linspace(0,2,nx), u, color='#003366', ls='--',lw=3)
    plt.ylim(0,2.5)
    plt.show()

